# =============================================================================
# aggregate_counts.py -- called by the Snakemake `aggregate_counts` rule
# (Snakemake's `script:` directive injects a `snakemake` object with
# `snakemake.input`, `snakemake.output`, etc.)
#
# WHAT: read every sample's Salmon quant.sf, map transcript IDs to gene IDs,
#       sum to gene-level counts per sample, and write one wide CSV.
# HOW: pandas for the per-sample table ops; a tx2gene mapping built once
#      from the GTF referenced in config.yaml (parsed from the GTF's
#      `transcript_id "..."; gene_id "...";` attribute fields).
# =============================================================================
import pandas as pd

samples = [s.split("/")[-1] for s in snakemake.input]
tx2gene_path = "data/reference/tx2gene.tsv"

tx2gene = pd.read_csv(tx2gene_path, sep="\t", header=None, names=["transcript_id", "gene_id"])

sample_counts = {}
for sample in samples:
    quant_path = f"results/salmon/{sample}/quant.sf"
    quant = pd.read_csv(quant_path, sep="\t")
    quant["transcript_id"] = quant["Name"].str.split(".").str[0]  # strip Ensembl version suffix
    merged = quant.merge(tx2gene, on="transcript_id", how="left")
    gene_counts = merged.groupby("gene_id")["NumReads"].sum()
    sample_counts[sample] = gene_counts

counts_df = pd.DataFrame(sample_counts).fillna(0).round(0).astype(int)
counts_df.to_csv(snakemake.output[0])

print(f"wrote {snakemake.output[0]}: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
