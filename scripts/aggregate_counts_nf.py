# =============================================================================
# aggregate_counts_nf.py -- called by the Nextflow `aggregate_counts` process.
# Same real logic as the Snakemake twin's scripts/aggregate_counts.py,
# adapted to read from Nextflow's working-directory-staged quant.sf files.
#
# WHAT/WHY: identical to the Snakemake twin -- see aggregate_counts.py.
# HOW: glob every staged quant.sf in the current work directory (Nextflow
# stages each input `path` into the process's own working directory), infer
# the sample name from its parent directory, then the same tx2gene groupby
# logic as the Snakemake twin. tx2gene.tsv is referenced by its staged
# basename ("tx2gene.tsv") -- see main.nf's aggregate_counts process, which
# declares it as a genuine `path` input so Nextflow stages a real copy into
# this task's own directory (a Nextflow task runs in its own isolated work
# directory, so a path relative to the PROJECT folder would not resolve here).
# =============================================================================
import glob
import pandas as pd

tx2gene_path = "tx2gene.tsv"
tx2gene = pd.read_csv(tx2gene_path, sep="\t", header=None, names=["transcript_id", "gene_id"])

quant_files = glob.glob("*/quant.sf")  # each sample's quant.sf staged under its own subdir

sample_counts = {}
for qf in quant_files:
    sample = qf.split("/")[0]
    quant = pd.read_csv(qf, sep="\t")
    quant["transcript_id"] = quant["Name"].str.split(".").str[0]
    merged = quant.merge(tx2gene, on="transcript_id", how="left")
    gene_counts = merged.groupby("gene_id")["NumReads"].sum()
    sample_counts[sample] = gene_counts

counts_df = pd.DataFrame(sample_counts).fillna(0).round(0).astype(int)
counts_df.to_csv("counts.csv")

print(f"wrote counts.csv: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
