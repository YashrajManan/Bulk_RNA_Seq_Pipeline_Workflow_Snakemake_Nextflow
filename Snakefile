# =============================================================================
# #1B -- Complete End-to-End RNA-seq Workflow -- SNAKEMAKE implementation
#
# Snakemake infers the DAG IMPLICITLY from filename patterns: each `rule`
# declares `input:`/`output:` paths, and Snakemake works out execution order
# by matching one rule's output pattern to another rule's input pattern --
# you never write "run step 2 after step 1" yourself.
# =============================================================================

configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
N_READS = config["subsample_reads"]
THREADS = config["threads"]

rule all:
    input:
        "results/counts.csv",
        "results/multiqc/multiqc_report.html"

# =============================================================================
# RULE: fetch_reads
# WHAT: for each real SRA run accession in config["samples"], download and
#       subsample reads to config["subsample_reads"] spots.
# WHY: real reads from the real GSE126848 liver-biopsy study, kept small at
#      DOWNLOAD time (not downloaded-then-subsampled) so the pipeline stays
#      fast.
# HOW: `prefetch {accession}` (SRA-tools) downloads the .sra file, then
#      `fastq-dump --split-files -X {n_reads}` extracts and subsamples it.
#      (fasterq-dump, the newer tool, hit a real, documented, unresolved
#      sra-tools bug -- github.com/ncbi/sra-tools issues #467/#598/#900,
#      "Failed to call external services" -- so this pipeline uses the older
#      fastq-dump instead, which doesn't use that cloud mechanism.)
#      This study's real reads are SINGLE-END (confirmed via ENA's
#      library_layout metadata) -- `--split-files` on single-end data writes
#      one file per sample, `{accession}_1.fastq`.
# =============================================================================
rule fetch_reads:
    output:
        r1="data/raw/{sample}_1.fastq"
    params:
        n=N_READS
    shell:
        """
        prefetch {wildcards.sample} -O data/raw
        fastq-dump --split-files -X {params.n} -O data/raw {wildcards.sample}
        """

# =============================================================================
# RULE: fastqc_raw
# WHAT: run FastQC on each sample's raw (pre-trim) single-end FASTQ file.
# WHY: the first real QC checkpoint -- confirms the downloaded/subsampled
#      reads are usable before spending time trimming/quantifying them.
# =============================================================================
rule fastqc_raw:
    input:
        r1="data/raw/{sample}_1.fastq"
    output:
        directory("results/fastqc_raw/{sample}")
    shell:
        "mkdir -p {output} && fastqc {input.r1} -o {output}"

# =============================================================================
# RULE: trim_reads
# WHAT: adapter + quality trim each sample's raw single-end reads with fastp.
# WHY: real sequencers make more errors near read ends -- trimming before
#      quantification stops those artifacts from biasing Salmon's mapping.
# HOW: single-end fastp mode (no `-I`/`-O` mate flags, since this study's
#      real reads are single-end). fastp's own JSON/HTML report doubles as a
#      second QC checkpoint, aggregated by MultiQC later.
# =============================================================================
rule trim_reads:
    input:
        r1="data/raw/{sample}_1.fastq"
    output:
        r1="data/trimmed/{sample}_1.trimmed.fastq",
        json="results/fastp/{sample}.fastp.json",
        html="results/fastp/{sample}.fastp.html"
    threads: THREADS
    shell:
        "fastp -i {input.r1} -o {output.r1} \
            --json {output.json} --html {output.html} --thread {threads}"

# =============================================================================
# RULE: download_reference + salmon_index (run ONCE, not per sample)
# WHAT: download the real reference transcriptome and build a Salmon index.
# WHY: Salmon pseudo-aligns reads against a TRANSCRIPTOME index, not the raw
#      genome -- this is the one-time setup step every sample's quantification
#      step below depends on. k=31 is Salmon's own recommended default k-mer
#      size for typical (>=75bp) reads.
# =============================================================================
rule download_reference:
    output:
        "data/reference/transcriptome.fa.gz"
    params:
        url=config["reference_transcriptome_url"]
    shell:
        "wget -O {output} {params.url}"

rule salmon_index:
    input:
        "data/reference/transcriptome.fa.gz"
    output:
        directory("data/reference/salmon_index")
    threads: THREADS
    shell:
        "salmon index -t {input} -i {output} -k 31 -p {threads}"

# =============================================================================
# RULE: salmon_quant
# WHAT: pseudo-align + quantify each sample's TRIMMED single-end reads
#       against the shared Salmon index.
# WHY: transcript-level abundance estimates per sample, the input to
#      gene-level count aggregation next.
# HOW: single-end Salmon mode (`-r` instead of `-1/-2`; `-l A` auto-detects
#      library strandedness). Single-end reads carry no directly observable
#      fragment length, so `--fldMean 200 --fldSD 25` supply Salmon's own
#      documented default estimate for typical Illumina single-end libraries.
# =============================================================================
rule salmon_quant:
    input:
        r1="data/trimmed/{sample}_1.trimmed.fastq",
        index="data/reference/salmon_index"
    output:
        directory("results/salmon/{sample}")
    threads: THREADS
    shell:
        "salmon quant -i {input.index} -l A -r {input.r1} \
             -p {threads} -o {output} --validateMappings \
             --fldMean 200 --fldSD 25"

# =============================================================================
# RULE: aggregate_counts
# WHAT: collapse every sample's Salmon transcript-level `quant.sf` output to
#       ONE gene x sample count matrix, using a real transcript-to-gene
#       mapping built from the same Ensembl release as the transcriptome.
# WHY: most downstream analyses (e.g. a DE tool) want gene-level counts, not
#      per-transcript abundances.
# =============================================================================
rule aggregate_counts:
    input:
        expand("results/salmon/{sample}", sample=SAMPLES)
    output:
        "results/counts.csv"
    script:
        "scripts/aggregate_counts.py"

# =============================================================================
# RULE: multiqc
# WHAT: aggregate every tool's per-sample logs (FastQC, fastp, Salmon) into
#       ONE HTML report.
# WHY: the real, standard way a bioinformatics pipeline reports its own QC
#      status to a human reviewer -- one report, not N scattered files.
# =============================================================================
rule multiqc:
    input:
        expand("results/fastqc_raw/{sample}", sample=SAMPLES),
        expand("results/fastp/{sample}.fastp.json", sample=SAMPLES),
        expand("results/salmon/{sample}", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        "multiqc results/ -o results/multiqc"
