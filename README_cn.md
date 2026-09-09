# PASApipeline-turbo（中文说明）

**[English README](README.md)**

一个面向性能与稳健性优化的 [PASApipeline](https://github.com/PASApipeline/PASApipeline) 分支，为大型真核基因组上的"按染色体分片 + SQLite 驱动"注释更新流程而开发。

基线：上游 master @ `cc1f8d7`（v2.5.3 + 2 个提交）。许可与版权归 Broad Institute 所有（BSD 3-Clause，见 [LICENSE](LICENSE)）。

## 为什么有这个分支

我们经常在大型植物基因组上用 PASA 更新 EVM 基因模型，采用的是按染色体分片的方案：每条染色体一个 SQLite 库、一条独立 PASA 进程链，转录本直接从 GTF 导入（StringTie/GeMoSeq/长读长装配），不再做重头比对。这个规模下原版流程暴露出多个性能瓶颈和一些正确性/稳健性问题；本分支在保持与标准 `Launch_PASA_pipeline.pl` 工作流行为兼容的前提下解决它们。

## 主要改动

### 性能（第一轮）

- **SQLite pragmas**（`PerlLib/DB_connect.pm`）：`journal_mode=MEMORY`、`synchronous=OFF`、`temp_store=MEMORY`、每连接约 1 GB 页缓存。PASA 库每次运行都重建，持久化 fsync 纯属浪费。
- **蛋白比对**（`PerlLib/fasta.ph`）：相同序列直接短路；结果按 MD5 对缓存；`fasta` 二进制每进程只解析一次。仅此一项就消除了 `cDNA_annotation_comparer.dbi` 中每条约染色体 10^5–10^6 次外部 `fasta` 进程调用——注释比较步骤原来的最大开销。
- **单染色体库上的真并行**：`validate_alignments_in_db.dbi` 现在把待验证比对分块到多线程（原来是每 contig 一线程——单染色体库上等于串行）；`assemble_clusters.dbi` 把簇列表分块到多线程并写 `*.partN.assemblies` 分块文件（与加载器的 glob 兼容）。两者都按基因组大小限制并发，避免大染色体撑爆内存（每个线程私有持有一份染色体序列）。
- **不再每簇起进程**：`SingleLinkageClusterer.pm` 用纯 Perl 并查集替代调用 `slclust` 子进程；`pasa` 二进制路径每进程只解析一次；临时文件用 pid+线程号+计数器命名（并行下不撞名）。
- **文件句柄复用**（`PerlLib/Fasta_retriever.pm`），并移除验证阶段多余的 `cdbfasta` 建索引；有 `samtools faidx` 索引就直接用。
- **导入路径批量 SQL**（`scripts/import_spliced_alignments.dbi`、`PerlLib/Ath1_cdnas.pm`）：预处理语句按连接复用，`last_insert_id()` 替代额外的 SELECT 往返，`avg_per_id` 并入 INSERT 省掉后续 UPDATE。
- `assign_clusters_by_stringent_alignment_overlap.dbi` 的**批量跨度预取**（两条 GROUP BY 代替 2N 次点查）、进度打印节流、以及末尾单次提交（中断运行可干净回滚）。
- `subcluster_builder.dbi` 的**排序成对循环 + 提前终止**（避免超大簇上的 O(k²) 爆炸）。
- 热路径调试打印收敛到 verbose 开关后（`assemble_clusters.dbi`、`cDNA_annotation_comparer.dbi`、`subcluster_builder.dbi`），每条染色体少写几百 MB 日志 IO。

### 性能（第二轮：装配 + 注释比较核心）

- **批量加载比对对象**（`PerlLib/Ath1_cdnas.pm`，用于 `assemble_clusters.dbi`、`subcluster_builder.dbi`、`cDNA_annotation_comparer.dbi`）：`create_alignment_objs_bulk()` / `get_alignment_objs_via_align_accs()` 每 500 条比对一次 SQL，取代原来的每条 2–3 次查询。`create_alignment_obj()` 通过共用的"行→对象"构建器保持行为完全一致（经 mock-DB 等价测试验证）。
- **1–2 条比对的纯 Perl 装配**（`PerlLib/CDNA/PASA_alignment_assembler.pm`）：`pasa_cpp_assemblies()` 对单条/成对装配在进程内计算，逐行镜像 `pasa_cpp` 的 `canMerge()`/`mergeAlignments()`（跨度重叠 → 同向判定 → 按位置 lockstep 检查剪接边界，fuzzlength=20 → 剪接感知合并）。这消掉了每对 2 次外部 `pasa` 进程调用——那是 comparer 里"（基因模型 × 转录本）"成对兼容性检查和装配步骤中单例簇的最大开销。设 `PASA_NO_PERL_PAIR_ASSEMBLY=1` 可回退到二进制。≥3 条比对才需要 `pasa` 二进制（惰性解析）。
- **亚簇比对对象只建一次**（`cDNA_annotation_comparer.dbi`）：FL-mode 和 nonFL-mode 两遍遍历共享缓存的比对对象（复用时还原 spliced orientation 到初建状态）。缓存设上限（约 30 万对象）防止超大染色体撑内存。
- **按基因缓存外显子重叠检查**：`model_overlaps_exon_segment()` 不再为每个候选重叠重新查库+解冻基因模型；缓存的外显子坐标在每次 `annotation_updates` 写入/失效处精确失效（`get_update_id`、反义/合并/拆分状态变更）。
- **每个基因对象只算一次序列**：`compare_updated_proteins()` 在蛋白/cDNA/CDS 序列已算过时跳过 `create_all_sequence_types()`（同一个注释模型要和很多装配比较）。
- **其它**：已并入的 EST 比对缓存起来而不是每次拼接重取；`stitch_nonFL_alignments_into_FL_alignments` 里 FL 推导的基因对象每个 FL-cDNA 只算一次，而不是每（EST × FL）对一次；亚簇列表跨两遍缓存；`DB_connect::get_last_insert_id()` 用驱动原生调用省掉一次 SELECT 往返；`Fasta_retriever::get_seq()` 边读边逐行去空白，不再对整条染色体序列跑第二遍正则。

### 性能（第三轮：注释比较，NYTProf 实测驱动）

- **ORF 扫描**（`PerlLib/Longest_orf.pm`）：`get_orfs()` 原来对每条序列做全（起始 × 终止）位置对扫描；现在对每个起始密码子二分查找其后的第一个终止子（选择结果完全一致，8000 条随机序列验证）。
- **翻译结果记忆化**（`PerlLib/Nuc_translator.pm`）：`translate_sequence()` / `get_protein()` 按精确序列记忆化（有容量上限，换遗传密码表时自动清空）——同一 CDS 序列在一个位点的众多装配中反复出现。Chr6 实测中翻译调用量降约 11 倍（198 万 → 17.6 万），耗时 68s → 约 10s。
- **不再重复算 ORF**：`validate_FLcdna_inferred_geneObjs` 原来对同一比对重复计算两次完全相同的 ORF；现在复用第一次结果。
- **批量 + 预处理写库**（`cDNA_annotation_comparer.dbi`、`PerlLib/DB_connect.pm`）：`status_link` / `annotation_link` 改成 150 行一条的多行 INSERT；comparer 的热点 SELECT/INSERT 走按连接缓存的预处理语句；单 contig 的比较写入包进一个事务（多 contig 库上不启用——见"稳健性"）。

### 真实数据实测

*Fragaria*（草莓属）基因组 Chr6 染色体，用 187,797 条 StringTie+GeMoSeq 证据转录本更新（WSL2，8 线程）：

| 步骤 | 原版 | turbo | 加速比 |
|---|---|---|---|
| 装配（`assemble_clusters.dbi`） | 103 s | 39 s | 2.7× |
| 注释比较（`cDNA_annotation_comparer.dbi`） | 318 s | 约 200 s | 1.6× |

输出与原版代码逐字节一致（装配 GFF3、更新 GFF3、数据库全量内容），经多轮 A/B 重复运行验证，外加用真实 `pasa` 二进制对纯 Perl 成对装配器做的 18,000 例模糊测试（零不一致）。

### 按染色体分片支持与唯一稳定标识符

- `assembly_db_loader.dbi`：装配 ID 带库（染色体）标签——`asmbl_Chr1_1`——分染色体产生的装配文件合并时不撞名。
- `cDNA_annotation_comparer.dbi`：新基因/新模型 ID 确定且带染色体标签——`novel_gene_15_Chr1`、`novel_model_27_Chr1`（多 contig 库追加 contig：`novel_gene_2_others_scaf99`）。编号用按 contig 计数器；不再有任何 `time()` 令牌，同样输入多次运行得到同样 ID。
- `dump_valid_annot_updates.dbi`：可变剪接模型后缀用确定性的递增编号（不再用时间戳后缀）。

### 稳健性

- `scripts/Pasa_init.pm`：本树自带的 `PerlLib` 现在优先于 `$PASAHOME/PerlLib`，conda 提供的 PASAHOME 不会再遮蔽本安装的模块。
- 二进制查找在 `which` 失败时回退到本树自带的 `bin/pasa`（以及 `$PASAHOME/bin/fasta`）。
- `subcluster_loader.dbi` 加载前先清空亚簇表，重跑幂等。
- **多 contig SQLite 安全**：当库里有多个 contig 时（比如装未定位 scaffold 的 "others" 桶，或任何标准单库 PASA 运行），注释比较器不再开跨全程的长事务——那会锁住其它 contig 工作线程并报 `database is locked`。单 contig 库（单工作线程）保留事务及其写批处理加速。
- **全仓库行尾统一为 LF** 并用 `.gitattributes` 固定（`* text=auto eol=lf`）；检出文件里的 CRLF 在 Linux 上会破坏 shebang 执行（`/usr/bin/env: 'perl\r': No such file or directory`）。
- **多次运行结果确定性**：簇 ID、装配 ID、亚簇 ID 和 gff3 输出顺序不再受 Perl 每进程哈希随机化影响——`assign_clusters_by_stringent_alignment_overlap.dbi` 按排序顺序遍历分组，`SingleLinkageClusterer` 以确定顺序返回成员和簇，`import_spliced_alignments.dbi` 按排序后的 contig 顺序分配 `align_id`，`PASA_alignment_assembler.pm` 里装配方向多数投票的打平裁决也是确定的。同样输入现在每次运行都得到逐字节相同的输出（上游 PASA 每次运行都会随机重排簇/装配编号）。

## 兼容性

通过 `Launch_PASA_pipeline.pl` 的标准（单库）运行照常工作；唯一可见区别是生成的标识符（装配、新基因）现在带库标签。比对验证阈值和更新逻辑不变。

设 `PASA_NO_PERL_PAIR_ASSEMBLY=1` 可让成对装配回退到外部 `pasa` 二进制（纯 Perl 快路径是默认行为）。

## 测试 / 回归脚手架

本分支用专门搭建的脚手架验证（合成 PASA SQLite 库、18,000 例成对装配器与编译版 `pasa` 二进制的模糊对比、以及完整真实数据 A/B 对照运行）。如果你改动这里的代码，发布前请重跑 A/B 对照。

## 引用

请引用上游 PASA 的论文（Haas 等，*Nucleic Acids Res.* 2003，PMID: 12829561；Haas 等，*BMC Bioinformatics* 2008，PMID: 18673596）；如果你使用了这些修改，也请引用本仓库。

## 上游 README

见 [wiki](https://github.com/PASApipeline/PASApipeline/wiki) 标签页文档。
