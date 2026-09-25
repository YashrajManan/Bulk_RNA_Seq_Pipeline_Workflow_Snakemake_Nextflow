# =============================================================================
# GIVEN -- generic GTF transcript-to-gene mapping builder. Run once, manually,
# after downloading the real GTF from config["tx2gene_url"]:
#   python3 scripts/build_tx2gene.py data/reference/annotation.gtf.gz data/reference/tx2gene.tsv
#
# This is TRUE boilerplate: it knows nothing project-specific --
# it just parses any standard GTF's `transcript_id "..."; gene_id "...";`
# attribute fields off `transcript` feature rows into a 2-column TSV. Safe to
# reuse for any future project that needs a tx2gene mapping.
# =============================================================================
import gzip
import re
import sys

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def build_tx2gene(gtf_path, out_path):
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    gene_re = re.compile(r'gene_id "([^"]+)"')
    n = 0
    with open_maybe_gz(gtf_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            tx_match = tx_re.search(attrs)
            gene_match = gene_re.search(attrs)
            if tx_match and gene_match:
                fout.write(f"{tx_match.group(1)}\t{gene_match.group(1)}\n")
                n += 1
    print(f"wrote {out_path}: {n} transcript -> gene mappings")

if __name__ == "__main__":
    build_tx2gene(sys.argv[1], sys.argv[2])
