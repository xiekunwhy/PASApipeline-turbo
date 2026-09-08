# PASApipeline-turbo

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

### Per-chromosome sharding support and unique, stable identifiers

- `assembly_db_loader.dbi`: assembly IDs carry the db (chromosome) tag — `asmbl_Chr1_1` — so per-chromosome assembly files merge without collisions.
- `cDNA_annotation_comparer.dbi`: novel gene/model IDs are deterministic and chromosome-tagged — `novel_gene_15_Chr1`, `novel_model_27_Chr1` (multi-contig databases append the contig: `novel_gene_2_others_scaf99`). Numbering uses per-contig counters; no `time()` tokens anywhere, so identical inputs give identical IDs across re-runs.
- `dump_valid_annot_updates.dbi`: alt-splice model suffixes are deterministic incrementers (no time-based suffix).

### Robustness

- `scripts/Pasa_init.pm`: the tree's own `PerlLib` now takes precedence over `$PASAHOME/PerlLib`, so a conda-provided PASAHOME no longer shadows this installation's modules.
- Binary discovery falls back to the tree's bundled `bin/pasa` (and `$PASAHOME/bin/fasta`) when `which` fails in batch-job environments.
- `subcluster_loader.dbi` purges subcluster tables before loading, making re-runs idempotent.

## Compatibility

Standard (single-database) runs through `Launch_PASA_pipeline.pl` keep working; the only visible difference is that generated identifiers (assemblies, novel genes) now carry the database tag. Alignment validation thresholds and update logic are unchanged.

## Citation

Please cite the upstream PASA publications (Haas et al., *Nucleic Acids Res.* 2003, PMID: 12829561; Haas et al., *BMC Bioinformatics* 2008, PMID: 18673596) and, if you use these modifications, this repository.

## Upstream README

See the [wiki](https://github.com/PASApipeline/PASApipeline/wiki) tab for documentation.

remaining to include in documentation:

1.  Defaults to SQLite. Users can override this by setting the DBI_DRIVER environment variable (see the DBI perldoc) to 'mysql' (if defaulting to SQLite is too bold for the next release, the default could be set the other way around, and SQLite chosen with DBI_DRIVER=SQLite).

2. The user instead sets a DATABASE parameter to either the name of the MySQL database, or the absolute pathname of the SQLite database.

3. When using SQLite, the $PASAHOME/pasa_conf/conf.txt file is optional (defaulting to ${PASAHOME}/pasa_conf/pasa.TEMPLATE for hooks) when SQLite is used. If it exists, it still will be used, and users can still override the pathname with PASACONF environment variable.
