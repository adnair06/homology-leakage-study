import shutil, subprocess
import pandas as pd

from homology.settings import OUT, MIN_SEQ_ID, COVERAGE, COV_MODE, THREADS
from homology.report import log

M8_COLS = ["query", "target", "fident", "alnlen", "qcov", "tcov", "evalue", "bits"]

def write_fasta(df, path):
    with open(path, "w") as f:
        for idx, seq in zip(df.index, df["seq"]):
            f.write(f">{idx}\n{seq}\n")

def mmseqs(*args):
    cmd = ["mmseqs", *[str(arg) for arg in args], "--threads", str(THREADS)]
    subprocess.run(cmd, check = True, capture_output = True)
    shutil.rmtree(OUT / "tmp", ignore_errors = True)

def search(query, target, m8):
    mmseqs("easy-search", query, target, m8, OUT / "tmp", "--format-output", ",".join(M8_COLS))

def best_hits(path):
    hits = pd.read_csv(path, sep = "\t", header = None, names = M8_COLS)
    hits = hits.sort_values("bits", ascending = False, kind = "stable")
    return hits.drop_duplicates(subset = "query", keep = "first").set_index("query")

def cluster():
    df = pd.read_parquet(OUT / "clean.parquet")
    write_fasta(df, OUT / "all.fasta")
    mmseqs("easy-cluster", OUT / "all.fasta", OUT / "cluster", OUT / "tmp",
           "--min-seq-id", MIN_SEQ_ID, "-c", COVERAGE, "--cov-mode", COV_MODE)

    pairs = pd.read_csv(OUT / "cluster_cluster.tsv", sep = "\t", header = None, names = ["rep", "member"])
    df["cluster"] = df.index.map(pairs.set_index("member")["rep"]).astype(int)
    sizes = df["cluster"].value_counts()
    df["cluster_size"] = df["cluster"].map(sizes)
    multi = df[df["cluster_size"] >= 2].groupby("cluster")["label"].mean()
    pure = ((multi == 0) | (multi == 1)).mean()

    log(f"\nMMseqs2 clustering (identity >= {MIN_SEQ_ID:.0%}, coverage >= {COVERAGE:.0%}, cov-mode {COV_MODE}) "
        f"gives {len(sizes):,} clusters.")
    log(f"{(df['cluster_size'] >= 2).mean():.1%} of proteins share a cluster with another protein, and "
        f"{pure:.1%} of multi-protein clusters hold only one class, so relatives almost always share a label.")

    df.to_parquet(OUT / "clustered.parquet")
