import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import accuracy_score, matthews_corrcoef, roc_auc_score

from homology.settings import OUT, BANDS, MIN_BAND_N
from homology.report import log
from homology.split import SPLITS
from homology.features import MODELS

STYLES = {
    ("split_random", "kmer"): ("#2a78d6", "-", "random, 3-mer"),
    ("split_random", "comp"): ("#2a78d6", "--", "random, composition"),
    ("split_cluster", "kmer"): ("#eb6834", "-", "clustered, 3-mer"),
    ("split_cluster", "comp"): ("#eb6834", "--", "clustered, composition"),
}

def band_names(edges):
    return [f"{a:g}-{b:g}" for a, b in zip(edges[:-1], edges[1:])]

def band_labels(edges):
    return [f"{a * 100:.0f}-{b * 100:.0f}%" for a, b in zip(edges[:-1], edges[1:])]

def band_cut(fident, edges):
    names = band_names(edges)
    return pd.cut(fident, edges, labels = names, include_lowest = True), names

def by_band():
    df = pd.read_parquet(OUT / "splits.parquet")
    preds = np.load(OUT / "preds_logistic.npz")
    best = pd.read_parquet(OUT / "best_hits.parquet")

    rows = []
    for split in SPLITS:
        test_idx = df.index[df[split] == "test"]
        hits = best[best["split"] == split].set_index("idx").reindex(test_idx)
        band, names = band_cut(hits["fident"], BANDS)
        band = band.cat.add_categories("no hit").fillna("no hit")
        for model in list(MODELS) + ["lookup"]:
            frame = pd.DataFrame({"y": df.loc[test_idx, "label"].to_numpy(), "band": band.to_numpy()},
                                 index = test_idx)
            if model == "lookup":
                frame["pred"] = hits["lookup_pred"].to_numpy()
                frame["prob"] = np.nan
            else:
                frame["pred"] = preds[f"{split}_{model}_pred"]
                frame["prob"] = preds[f"{split}_{model}_prob"]
            for name in names + ["no hit"]:
                group = frame[frame["band"] == name]
                ok = len(group) >= MIN_BAND_N and group["y"].nunique() == 2
                rows.append({
                    "split": split, "model": model, "band": name, "n": len(group),
                    "pos_frac": group["y"].mean() if len(group) else np.nan,
                    "auc": roc_auc_score(group["y"], group["prob"]) if ok and model != "lookup" else np.nan,
                    "acc": accuracy_score(group["y"], group["pred"]) if ok else np.nan,
                    "mcc": matthews_corrcoef(group["y"], group["pred"]) if ok else np.nan,
                })
    table = pd.DataFrame(rows)

    log("\nTest ROC-AUC by sequence identity to the nearest training protein (n test proteins):")
    for model, label in MODELS.items():
        log(f"  {label:<13}{'random':<17}clustered")
        for name, pretty in zip(["no hit"] + names, ["no hit"] + band_labels(BANDS)):
            cells = []
            for split in SPLITS:
                row = table[(table["split"] == split) & (table["model"] == model) & (table["band"] == name)].iloc[0]
                cells.append(f"{row['auc']:.3f} ({row['n']:,})")
            log(f"  {pretty:<13}{cells[0]:<17}{cells[1]}")

    table.to_csv(OUT / "results_bands.csv", index = False)

def plot():
    table = pd.read_csv(OUT / "results_bands.csv")
    order = band_names(BANDS)

    n = table[table["model"] == "comp"].pivot(index = "band", columns = "split", values = "n")
    nohit = n.loc["no hit"]
    xlabels = []
    for label, name in zip(band_labels(BANDS), order):
        xlabels.append(f"{label}\n{n.loc[name, 'split_random']:,} / {n.loc[name, 'split_cluster']:,}")

    fig, ax = plt.subplots(figsize = (10, 5.8))
    x_nohit = -1.3
    for (split, model), (color, ls, label) in STYLES.items():
        sub = table[(table["split"] == split) & (table["model"] == model)].set_index("band")
        ax.plot(range(len(order)), [sub.loc[name, "auc"] for name in order], color = color, ls = ls,
                marker = "o", lw = 2, ms = 6, label = label)
        ax.plot([x_nohit + (0.12 if split == "split_cluster" else -0.12)],
                [sub.loc["no hit", "auc"]], color = color, marker = "s" if model == "kmer" else "D",
                mfc = color if model == "kmer" else "white", ls = "none", ms = 7)

    ax.axhline(0.5, color = "grey", lw = 1, ls = ":")
    ax.text(len(order) - 0.55, 0.505, "chance", color = "grey", ha = "right", va = "bottom", fontsize = 9)
    ax.axvline(-0.65, color = "lightgrey", lw = 1)
    ax.set_xticks([x_nohit] + list(range(len(order))))
    ax.set_xticklabels([f"no hit\n{nohit['split_random']:,} / {nohit['split_cluster']:,}"] + xlabels,
                       fontsize = 9)
    ax.set_xlim(x_nohit - 0.6, len(order) - 0.5)
    ax.set_xlabel("sequence identity to nearest training protein\n"
                  "(n test proteins: random / clustered)")
    ax.set_ylabel("ROC-AUC")
    ax.grid(axis = "y", color = "#eeeeee")
    ax.spines[["top", "right"]].set_visible(False)
    handles, _ = ax.get_legend_handles_labels()
    handles += [Line2D([], [], color = "#52514e", marker = "s", ls = "none", ms = 7, label = "no hit, 3-mer"),
                Line2D([], [], color = "#52514e", marker = "D", mfc = "white", ls = "none", ms = 7,
                       label = "no hit, composition")]
    ax.legend(handles = handles, loc = "lower left", bbox_to_anchor = (0.27, 0.02), ncol = 3,
              fontsize = 9, frameon = False)
    plt.tight_layout()
    plt.savefig(OUT / "auc_by_identity.png", dpi = 150)
