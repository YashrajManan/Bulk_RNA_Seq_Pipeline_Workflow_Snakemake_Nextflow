#!/usr/bin/env nextflow
// =============================================================================
// #1B -- Complete End-to-End RNA-seq Workflow -- NEXTFLOW implementation
//
// Nextflow's DAG is EXPLICIT: each `process` declares input/output channels,
// and you WIRE processes together yourself by piping one process's output
// channel into the next process's call -- a genuinely different mental model
// from Snakemake's implicit filename-pattern matching. Same real pipeline,
// same real tools as the Snakemake twin -- independently implemented, not a
// mechanical translation.
// =============================================================================

nextflow.enable.dsl=2

params.samples_csv = "samples.csv"   // sample,group
params.n_reads = 500000
params.transcriptome_url = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.outdir = "results"
// tx2gene.tsv is built once (scripts/build_tx2gene.py) before running this pipeline.
params.tx2gene = "data/reference/tx2gene.tsv"

// =============================================================================
// PROCESS: fetch_reads
// WHAT/WHY: real reads, subsampled at download time. fasterq-dump (the newer
// tool) hits a real, documented, unresolved sra-tools bug (github.com/
// ncbi/sra-tools issues #467/#598/#900, "Failed to call external services")
// -- use the older `fastq-dump --split-files -X` instead. This study's real
// reads are SINGLE-END (confirmed via ENA's library_layout metadata), so
// there is only ONE output FASTQ per sample, not a mate pair.
// =============================================================================
process fetch_reads {
    tag "$sample"
    publishDir "${params.outdir}/raw", mode: 'copy'

    input:
    val sample

    output:
    tuple val(sample), path("${sample}_1.fastq")

    script:
    """
    prefetch ${sample}
    fastq-dump --split-files -X ${params.n_reads} ${sample}
    """
}

// =============================================================================
// PROCESS: fastqc_raw -- single-end QC, same real purpose as the Snakemake twin.
// =============================================================================
process fastqc_raw {
    tag "$sample"
    publishDir "${params.outdir}/fastqc_raw/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(r1)

    output:
    path "*"

    script:
    "fastqc ${r1}"
}

// =============================================================================
// PROCESS: trim_reads -- single-end fastp trimming: no `-I`/`-O` mate flags,
// since there is no second read file.
// =============================================================================
process trim_reads {
    tag "$sample"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(sample), path(r1)

    output:
    tuple val(sample), path("${sample}_1.trimmed.fastq"), emit: trimmed
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    "fastp -i ${r1} -o ${sample}_1.trimmed.fastq \
       --json ${sample}.fastp.json --html ${sample}.fastp.html"
}

// =============================================================================
// PROCESS: download_reference + salmon_index (run ONCE) -- same real one-time
// index build as the Snakemake twin.
// =============================================================================
process download_reference {
    output:
    path "transcriptome.fa.gz"

    script:
    "wget -O transcriptome.fa.gz ${params.transcriptome_url}"
}

process salmon_index {
    input:
    path transcriptome

    output:
    path "salmon_index"

    script:
    "salmon index -t ${transcriptome} -i salmon_index -k 31"
}

// =============================================================================
// PROCESS: salmon_quant -- single-end mode: `-r {r1}` instead of `-1/-2`,
// plus `--fldMean`/`--fldSD` since single-end reads carry no directly
// observable fragment length (200/25 is Salmon's own documented default
// recommendation). Salmon always names its per-sample output file
// `quant.sf`, so this writes into a directory NAMED after the sample
// (`-o ${sample}`) -- otherwise every sample's task would produce an
// identically-named file, which collides when `aggregate_counts` collects
// all 8 samples' outputs into one input list (Nextflow has to stage them
// all into a single work directory).
// =============================================================================
process salmon_quant {
    tag "$sample"
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample), path(r1)
    path index

    output:
    tuple val(sample), path("${sample}")

    script:
    "salmon quant -i ${index} -l A -r ${r1} -o ${sample} --validateMappings \
       --fldMean 200 --fldSD 25"
}

// =============================================================================
// PROCESS: aggregate_counts -- same real gene-level aggregation logic as the
// Snakemake twin's Python script, adapted for Nextflow's isolated per-task
// work directories: the script is referenced via Nextflow's built-in
// `${projectDir}` (resolves to wherever main.nf itself lives, regardless of
// the task's own isolated directory -- a plain relative path would not
// resolve there), and tx2gene.tsv is declared as a genuine `path` input so
// Nextflow stages a real copy into the task directory under its own
// basename (see aggregate_counts_nf.py's `tx2gene_path = "tx2gene.tsv"`).
// =============================================================================
process aggregate_counts {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quant_files
    path tx2gene

    output:
    path "counts.csv"

    script:
    "python3 \"${projectDir}/scripts/aggregate_counts_nf.py\""
}

// =============================================================================
// PROCESS: multiqc
// =============================================================================
process multiqc {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    "multiqc ."
}

// =============================================================================
// WORKFLOW block -- this is where Nextflow's channel-wiring actually happens,
// the explicit-DAG counterpart to Snakemake's implicit filename matching.
// =============================================================================
workflow {
    samples_ch = Channel.fromPath(params.samples_csv).splitCsv(header:true).map { row -> row.sample }

    raw_reads = fetch_reads(samples_ch)
    fastqc_raw(raw_reads)
    trimmed = trim_reads(raw_reads)

    ref = download_reference()
    index = salmon_index(ref)

    quant = salmon_quant(trimmed.trimmed, index)
    tx2gene_ch = Channel.fromPath(params.tx2gene)
    aggregate_counts(quant.map { it[1] }.collect(), tx2gene_ch)

    qc_files = fastqc_raw.out.mix(trimmed.json, quant.map { it[1] }).collect()
    multiqc(qc_files)
}
