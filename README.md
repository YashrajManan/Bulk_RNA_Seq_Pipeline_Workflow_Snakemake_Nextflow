# RNA-seq End-to-End Workflow (Snakemake + Nextflow)

A complete, reproducible RNA-seq pipeline — raw FASTQ to a gene-level count matrix — implemented
**twice, independently**, in **Snakemake** (implicit, filename-pattern-driven DAG) and **Nextflow**
(explicit channel/process dataflow DSL2 model), run against real human liver-biopsy RNA-seq reads from
a published clinical study, not synthetic test data.

## The question

Can the same real RNA-seq pipeline be built as a genuine, one-command, reproducible workflow — not a
notebook run cell by cell — in both of the two dominant bioinformatics workflow managers, and do two
structurally different execution models, built independently, converge on identical output?

## Data

Real human liver-biopsy RNA-seq reads from **GSE126848** (Suppli et al. 2019, *American Journal of
Physiology*) — a published study of 57 real liver biopsies across four clinical groups: healthy
normal-weight, healthy obese, NAFL (simple steatosis), and NASH (non-alcoholic steatohepatitis). This
project uses **8 real SRA runs** (2 per clinical group), each subsampled to 500,000 reads at download
time. Reads are genuinely **single-end** — discovered live via ENA's own `library_layout` metadata, not
assumed (see `PROJECT_NARRATIVE.md`).

## Pipeline

`fetch reads` (real SRA download + subsample) → `FastQC` (raw QC) → `fastp` (adapter/quality trimming)
→ `Salmon index` (build once, real Ensembl GRCh38 transcriptome) → `Salmon quant` (pseudo-alignment +
quantification per sample) → `aggregate to gene-level counts` (real transcript-to-gene mapping) →
`MultiQC` (aggregate every tool's QC report). Identical DAG, implemented independently in both
`Snakefile` and `main.nf`.

## Real results

**Mapping rate:** 88.25%–94.13% across all 8 samples, consistently strong across all four clinical
groups. **Gene count matrix:** 34,141 genes × 8 samples. **Snakemake's and Nextflow's
independently-computed count matrices are byte-identical — 0 differing cells across all 273,128 values.**

Full results: `results/` (Snakemake run), `results_nf/` (Nextflow run) — `counts.csv` (gene × sample
matrix) and `multiqc/multiqc_report.html` (aggregated QC) in each.

## Interpretation

Two structurally different workflow engines, built and debugged independently against the identical DAG,
converging on a byte-identical final count matrix is the strongest verification available here, since
neither implementation was checked against an external ground-truth answer key — they were checked against
each other. Each engine's own execution model surfaced real, engine-specific issues during the build (full
account in `PROJECT_NARRATIVE.md`), from environment and tooling quirks to an initial assumption about
sequencing library layout that was checked against real ENA metadata rather than taken at face value;
working through each independently, then requiring both to agree exactly, is what makes the byte-identical
outcome a meaningful cross-validation rather than a coincidence of two pipelines sharing the same unverified
assumption.

Bioinformatically, the pipeline's own numbers are consistent with real, structurally sound sequencing data:
mapping rates in the high-80s-to-low-90s percent range are typical for real human RNA-seq against a
matched-species, matched-build reference, and no sample's mapping rate is an outlier low enough to suggest a
corrupted or mislabeled file. Consistency across all four clinical groups (healthy normal-weight, healthy
obese, NAFL, NASH) means pipeline quality did not degrade for any one disease category. At 500,000
subsampled reads per sample, the resulting counts are real but intentionally shallow: this project's
deliverable is a trustworthy, reproducible upstream count matrix, not a differential-expression result —
whether real, graded expression differences exist across the healthy-to-NASH severity spectrum is a separate
question for a downstream differential-expression analysis at full read depth, out of this project's scope.

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.

## Files

- `Snakefile`, `config.yaml` — the Snakemake implementation.
- `main.nf`, `nextflow.config`, `samples.csv` — the Nextflow implementation.
- `scripts/build_tx2gene.py` — builds the real transcript-to-gene mapping from the Ensembl GTF.
- `scripts/aggregate_counts.py`, `scripts/aggregate_counts_nf.py` — gene-level aggregation (Snakemake
  and Nextflow versions respectively; identical real logic, adapted to each tool's staging model).
- `environment.yml`, `environment-tools.yml` — two deliberately separate conda environments
  (workflow-orchestration tools vs. bioinformatics CLI tools — reasoning in `PROJECT_NARRATIVE.md`).
- `results/`, `results_nf/` — real output from the Snakemake and Nextflow runs respectively.

## Caveats (real, not boilerplate)

500,000 subsampled reads per sample is deliberately shallow — enough to validate the pipeline runs
correctly and produces real, biologically sensible counts (strong, consistent mapping rates; correctly
structured output), not enough for a publication-grade differential-expression result. This project's
deliverable is the reproducible upstream pipeline itself; a real DE analysis at full depth (structurally
identical to another project in this portfolio) is explicitly out of scope. Single-end sequencing means
Salmon's fragment-length distribution is a disclosed estimate (`--fldMean 200 --fldSD 25`), not directly
measured the way paired-end data would allow. Only 8 of GSE126848's 57 real samples were used, so the
four-group split here (2 samples/group) is a structural proof of concept, not a statistically powered
comparison.
