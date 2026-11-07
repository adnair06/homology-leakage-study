import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold, train_test_split

from homology.settings import OUT, SEED, N_FOLDS
from homology.report import log

SPLITS = {"split_random": "random", "split_cluster": "clustered"}

def split():
    df = pd.read_parquet(OUT / "clustered.parquet")

    folds = StratifiedGroupKFold(n_splits = N_FOLDS, shuffle = True, random_state = SEED)
    _, positions = next(folds.split(df, df["label"], groups = df["cluster"]))
    cluster_test = df.index[positions]
    _, random_test = train_test_split(df.index, test_size = len(cluster_test),
                                      stratify = df["label"], random_state = SEED)

    df["split_cluster"] = "train"
    df.loc[cluster_test, "split_cluster"] = "test"
    df["split_random"] = "train"
    df.loc[random_test, "split_random"] = "test"

    rows = []
    for name in SPLITS:
        train, test = df[df[name] == "train"], df[df[name] == "test"]
        shared = set(train["cluster"]) & set(test["cluster"])
        rows.append({
            "split": name, "n_train": len(train), "n_test": len(test),
            "pos_train": train["label"].mean(), "pos_test": test["label"].mean(),
            "test_clusters": test["cluster"].nunique(),
            "shared_clusters": len(shared),
            "test_in_shared_cluster": test["cluster"].isin(shared).mean(),
            "test_len_median": test["length"].median(),
        })
    summary = pd.DataFrame(rows).set_index("split")

    rand, clust = summary.loc["split_random"], summary.loc["split_cluster"]
    log(f"\nBoth test sets hold {int(rand['n_test']):,} proteins ({rand['n_test'] / len(df):.0%}), "
        f"{rand['pos_test']:.1%} of them enzymes.")
    log(f"Test proteins that share a cluster with a training protein: {rand['test_in_shared_cluster']:.1%} "
        f"in the random split, {clust['test_in_shared_cluster']:.1%} in the clustered split.")

    summary.to_csv(OUT / "split_summary.csv")
    df.to_parquet(OUT / "splits.parquet")
