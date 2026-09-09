# PASApipeline-turbo

**[中文说明见 README_cn.md](README_cn.md)**

A performance- and robustness-oriented fork of [PASApipeline](https://github.com/PASApipeline/PASApipeline), developed for chromosome-sharded, SQLite-backed annotation updates on large eukaryotic genomes.

Base: upstream master @ `cc1f8d7` (v2.5.3 + 2 commits). License and copyright remain with the Broad Institute (BSD 3-Clause, see [LICENSE](LICENSE)).

## Why this fork

We routinely update EVM gene models with PASA on large plant genomes, using a chromosome-sharded setup: one SQLite database per chromosome, one PASA process chain per chromosome, with transcripts imported directly from GTF (StringTie/GeMoSeq/long-read assemblies) instead of de novo spliced alignment. At this scale the stock pipeline shows several bottlenecks and a few correctness/robustness issues; this fork addresses them while staying behavior-compatible with the standard `Launch_PASA_pipeline.pl` workflow.

## Main changes

### Performance

- **SQLite pragmas** (`PerlLib/DB_connect.pm`): `journal_mode=MEMORY`, `synchronous=OFF`, `temp_store=MEMORY`, ~1 GB page cache per connection. PASA databases are rebuilt every run, so durability fsyncs are wasted work.
- **Protein comparison** (`PerlLib/fasta.ph`): identical sequences short-circuit; results cached by MD5 pair; the `fasta` binary is resolved once per process. This removes ~10^5–10^6 external `fasta` process spawns per chromosome in `cDNA_annotation_comparer.dbi`, the dominant cost of the annotation-comparison step.
- **Real parallelism on single-chromosome databases**: `validate_alignments_in_db.dbi` now chunks pending alignments across threads (was: one thread per contig, i.e. serial on a one-chromosome db); `assemble_clusters.dbi` chunks cluster lists across threads and writes per-chunk `*.partN.assemblies` files (compatible with the loader's glob). Both cap concurrency by genome size so big chromosomes do not exhaust RAM (each thread holds a private copy of the chromosome sequence).
- **No more per-cluster process spawning**: `SingleLinkageClusterer.pm` uses a pure-Perl union-find instead of shelling out to `slclust`; the `pasa` binary path is resolved once per process instead of per cluster; temp files use pid+thread-id+counter names (collision-proof under parallelism).
- **Filehandle reuse** (`PerlLib/Fasta_retriever.pm`) and removal of unused `cdbfasta` index building in validation; `samtools faidx` indices are used when present.
- **Batch SQL in the import path** (`scripts/import_spliced_alignments.dbi`, `PerlLib/Ath1_cdnas.pm`): prepared statements are reused per connection, `last_insert_id()` replaces an extra SELECT round-trip, and `avg_per_id` is folded into the INSERT instead of a follow-up UPDATE.
- **Bulk span prefetch** in `assign_clusters_by_stringent_alignment_overlap.dbi` (two GROUP BY queries instead of 2N point queries), throttled progress printing, and a single final commit so an interrupted run rolls back cleanly.
- **Sorted pairwise loop with early break** in `subcluster_builder.dbi` (avoids O(k^2) blowups on giant clusters).
- Hot-path debug prints are now guarded behind verbose flags in `assemble_clusters.dbi`, `cDNA_annotation_comparer.dbi` and `subcluster_builder.dbi` (removes hundreds of MB of log I/O per chromosome).

### Performance, round 2 (assembly + annotation-comparison cores)

- **Bulk alignment loading** (`PerlLib/Ath1_cdnas.pm`, used by `assemble_clusters.dbi`, `subcluster_builder.dbi` and `cDNA_annotation_comparer.dbi`): `create_alignment_objs_bulk()` / `get_alignment_objs_via_align_accs()` build all alignment objects of a cluster/subcluster with one SQL query per 500 alignments instead of two/three queries per alignment. `create_alignment_obj()` keeps its exact behavior through a shared row-to-object builder (verified by mock-DB equivalence tests).
- **Pure-Perl assembly for 1–2 alignments** (`PerlLib/CDNA/PASA_alignment_assembler.pm`): `pasa_cpp_assemblies()` now computes singleton and pairwise assemblies in-process, mirroring `canMerge()`/`mergeAlignments()` of `pasa_cpp` line by line (span overlap → equal fixed orientation → lockstep segment walk with structural splice-junction flags and fuzzlength=20 → splice-aware merge). This removes two external `pasa` process spawns per pair, the dominant cost of the comparer's per-(gene model × transcript) compatibility checks and of singleton clusters in the assembly step. Set `PASA_NO_PERL_PAIR_ASSEMBLY=1` to fall back to the binary. The `pasa` binary is only required for 3+ alignments now (resolved lazily).
- **Build subcluster alignment objects once, not twice** (`cDNA_annotation_comparer.dbi`): the FL-mode and nonFL-mode passes share cached alignment objects per subcluster (spliced orientations are restored on reuse so pass two sees freshly-built state). Cache is capped (~300k objects) to bound memory on very large chromosomes.
- **Exon-overlap checks cached per gene**: `model_overlaps_exon_segment()` no longer re-queries and re-thaws gene models for every candidate overlap; cached exon coordinates are invalidated exactly where `annotation_updates` rows are stored or (in)validated (`get_update_id`, antisense/merge/split status changes).
- **Sequence computation once per gene object**: `compare_updated_proteins()` skips `create_all_sequence_types()` when the object's protein/cDNA/CDS sequences are already computed (the same annotated model is compared against many assemblies).
- **Other**: previously-incorporated EST alignments are cached instead of re-fetched per stitching attempt; FL-inferred gene objects in `stitch_nonFL_alignments_into_FL_alignments` are computed once per FL-cDNA instead of per (EST × FL) pair; subcluster listings are cached across the two comparer passes; `DB_connect::get_last_insert_id()` uses the driver-native call instead of an extra SELECT round-trip; `Fasta_retriever::get_seq()` strips whitespace per line while reading instead of a second full pass over chromosome-length strings.

### Performance, round 3 (annotation comparer, NYTProf-profiled)

- **ORF scanning** (`PerlLib/Longest_orf.pm`): `get_orfs()` used to scan all (start × stop) position pairs per sequence; it now binary-searches the first stop past each start codon (identical picks, verified on 8000 random sequences).
- **Translation memoization** (`PerlLib/Nuc_translator.pm`): `translate_sequence()` / `get_protein()` are memoized per exact sequence (size-capped, cleared on genetic-code change) — the same CDS sequences recur across the many assemblies of a locus. Translation count dropped ~11× (1.98M → 176k) in the Chr6 profile, 68s → ~10s.
- **No double ORF computation**: `validate_FLcdna_inferred_geneObjs` recomputed the identical ORF for the same alignment twice; it now reuses the first result.
- **Batched + prepared DB writes** (`cDNA_annotation_comparer.dbi`, `PerlLib/DB_connect.pm`): `status_link` / `annotation_link` are written as 150-row multi-INSERTs; the comparer's hot SELECTs/INSERTs go through per-connection prepared-statement caching; each single-contig comparison's writes are wrapped in one transaction (disabled on multi-contig databases — see Robustness).

### Per-chromosome sharding support and unique, stable identifiers

- `assembly_db_loader.dbi`: assembly IDs carry the db (chromosome) tag — `asmbl_Chr1_1` — so per-chromosome assembly files merge without collisions.
- `cDNA_annotation_comparer.dbi`: novel gene/model IDs are deterministic and chromosome-tagged — `novel_gene_15_Chr1`, `novel_model_27_Chr1` (multi-contig databases append the contig: `novel_gene_2_others_scaf99`). Numbering uses per-contig counters; no `time()` tokens anywhere, so identical inputs give identical IDs across re-runs.
- `dump_valid_annot_updates.dbi`: alt-splice model suffixes are deterministic incrementers (no time-based suffix).

### Measured on real data

Chromosome Chr6 of a *Fragaria* genome, updated with 187,797 StringTie+GeMoSeq evidence transcripts (WSL2, 8 threads):

| step | stock | turbo | speedup |
|---|---|---|---|
| assembly (`assemble_clusters.dbi`) | 103 s | 39 s | 2.7× |
| annotation comparison (`cDNA_annotation_comparer.dbi`) | 318 s | ~200 s | 1.6× |

Outputs are byte-identical to the stock code (assemblies GFF3, updated GFF3, full database content), verified by repeated A/B runs plus an 18,000-case fuzz test of the pure-Perl pair assembler against the real `pasa` binary.

### Robustness

- `scripts/Pasa_init.pm`: the tree's own `PerlLib` now takes precedence over `$PASAHOME/PerlLib`, so a conda-provided PASAHOME no longer shadows this installation's modules.
- Binary discovery falls back to the tree's bundled `bin/pasa` (and `$PASAHOME/bin/fasta`) when `which` fails in batch-job environments.
- `subcluster_loader.dbi` purges subcluster tables before loading, making re-runs idempotent.
- **Multi-contig SQLite safety**: the annotation comparer no longer opens a long-lived transaction when the database holds more than one contig (e.g. an "others" bucket of unplaced scaffolds, or any standard single-database PASA run) — doing so locked out the other contig worker threads with `database is locked`. Single-contig databases (one worker thread) keep the transaction and its write-batching speedup.
- **Line endings normalized to LF repository-wide** and pinned via `.gitattributes` (`* text=auto eol=lf`); CRLF in checked-out scripts breaks shebang execution on Linux (`/usr/bin/env: 'perl\r': No such file or directory`).
- **Run-to-run determinism**: cluster IDs, assembly IDs, subcluster IDs and gff3 output order no longer depend on Perl's per-process hash randomization — `assign_clusters_by_stringent_alignment_overlap.dbi` iterates groups in sorted order, `SingleLinkageClusterer` returns members and clusters in deterministic order, `import_spliced_alignments.dbi` assigns `align_id`s in sorted contig order, and the assembly-orientation majority vote in `PASA_alignment_assembler.pm` breaks ties deterministically. Identical inputs now give byte-identical outputs across runs (upstream PASA renumbers clusters/assemblies randomly on every run).

## Compatibility

Standard (single-database) runs through `Launch_PASA_pipeline.pl` keep working; the only visible difference is that generated identifiers (assemblies, novel genes) now carry the database tag. Alignment validation thresholds and update logic are unchanged.

Set `PASA_NO_PERL_PAIR_ASSEMBLY=1` to fall back to the external `pasa` binary for pairwise assemblies (the pure-Perl fast path is the default).

## Testing / regression harness

The fork was validated with a purpose-built harness (synthetic PASA SQLite fixtures, an 18k-case fuzz comparison of the in-process pair assembler against the compiled `pasa` binary, and a full real-data A/B run). If you modify this code, re-run an A/B comparison before shipping.

## Citation

Please cite the upstream PASA publications (Haas et al., *Nucleic Acids Res.* 2003, PMID: 12829561; Haas et al., *BMC Bioinformatics* 2008, PMID: 18673596) and, if you use these modifications, this repository.

## Upstream README

See the [wiki](https://github.com/PASApipeline/PASApipeline/wiki) tab for documentation.

remaining to include in documentation:

1.  Defaults to SQLite. Users can override this by setting the DBI_DRIVER environment variable (see the DBI perldoc) to 'mysql' (if defaulting to SQLite is too bold for the next release, the default could be set the other way around, and SQLite chosen with DBI_DRIVER=SQLite).

2. The user instead sets a DATABASE parameter to either the name of the MySQL database, or the absolute pathname of the SQLite database.

3. When using SQLite, the $PASAHOME/pasa_conf/conf.txt file is optional (defaulting to ${PASAHOME}/pasa_conf/pasa.TEMPLATE for hooks) when SQLite is used. If it exists, it still will be used, and users can still override the pathname with PASACONF environment variable.
