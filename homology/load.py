import pandas as pd

from homology.settings import DATA, OUT, MIN_LEN, MAX_LEN, AMINO_ACIDS
from homology.report import log

def drop(df, mask, reason, removed):
    lost = df[mask]
    removed.append((len(lost), int(lost["label"].sum()), int((lost["label"] == 0).sum()), reason))
    return df[~mask]

def load():
    frames = []
    for name, expected in [("enzymes", 1), ("non_enzymes", 0)]:
        frame = pd.read_csv(DATA / f"{name}.tsv.gz", sep = "\t", dtype = str, keep_default_na = False)
        frame["expected"] = expected
        frames.append(frame)
    df = pd.concat(frames, ignore_index = True)
    df["seq"] = df["Sequence"].str.replace(r"\s+", "", regex = True).str.upper()
    df["length"] = df["seq"].str.len()
    df["label"] = (df["EC number"].str.strip() != "").astype(int)
    total = len(df)

    removed = []
    df = drop(df, df["label"] != df["expected"], "in the enzyme file but without an EC number", removed)
    df = drop(df, df["seq"].str.contains(f"[^{AMINO_ACIDS}]"), "contain a non-standard residue (X, B, Z, U, O)", removed)
    df = drop(df, df["length"] < MIN_LEN, f"shorter than {MIN_LEN} residues", removed)
    df = drop(df, df["length"] > MAX_LEN, f"longer than {MAX_LEN:,} residues", removed)
    labels = df.groupby("seq")["label"].transform("nunique")
    df = drop(df, labels > 1, "identical sequences with conflicting labels (all copies)", removed)
    df = drop(df, df["seq"].duplicated(keep = "first"), "extra copies of identical sequences", removed)

    log(f"Loaded {total:,} reviewed UniProt proteins. Removed (enzymes / non-enzymes):")
    for count, enzymes, others, reason in removed:
        log(f"  {count:>6,}  ({enzymes:,} / {others:,})  {reason}")
    log(f"Kept {len(df):,} proteins: {df['label'].sum():,} enzymes ({df['label'].mean():.1%}) "
        f"and {(df['label'] == 0).sum():,} non-enzymes.")

    df[["Entry", "seq", "length", "label"]].reset_index(drop = True).to_parquet(OUT / "clean.parquet")
