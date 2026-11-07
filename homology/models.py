import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef,
                             average_precision_score, roc_auc_score)
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from homology.settings import OUT, MIN_SEQ_ID, COVERAGE, SEED, CV_FOLDS, N_JOBS, C_COMP, C_KMER, BANDS
from homology.report import log
from homology.cluster import M8_COLS, write_fasta, search, best_hits
from homology.split import SPLITS
from homology.features import MODELS
from homology.bands import band_cut, band_labels

def scores(y, pred):
    return {"accuracy": accuracy_score(y, pred),
            "balanced_acc": balanced_accuracy_score(y, pred),
            "f1": f1_score(y, pred, zero_division = 0),
            "mcc": matthews_corrcoef(y, pred)}

def lookup():
    df = pd.read_parquet(OUT / "splits.parquet")
    labels = df["label"]
    folder = OUT / "search"
    folder.mkdir(exist_ok = True)

    results, consistency, best_all = [], [], []
    for split in SPLITS:
        name = split.replace("split_", "")
        train_idx = df.index[df[split] == "train"]
        test_idx = df.index[df[split] == "test"]
        m8 = folder / f"{name}.m8"
        write_fasta(df.loc[train_idx], folder / f"{name}_train.fasta")
        write_fasta(df.loc[test_idx], folder / f"{name}_test.fasta")
        search(folder / f"{name}_test.fasta", folder / f"{name}_train.fasta", m8)

        hits = pd.read_csv(m8, sep = "\t", header = None, names = M8_COLS)
        best = best_hits(m8)

        majority = int(labels[train_idx].mean() >= 0.5)
        neighbour = best["target"].map(labels).reindex(test_idx)
        hit = neighbour.notna()
        pred = neighbour.fillna(majority).astype(int)
        y = labels[test_idx]

        result = {"split": split, "n_test": len(test_idx), "hit_rate": hit.mean(),
                  "fallback_class": majority}
        result.update(scores(y, pred))

        strict = (hits["fident"] >= MIN_SEQ_ID) & (hits["qcov"] >= COVERAGE) & (hits["tcov"] >= COVERAGE)
        nearest = best.reindex(test_idx)
        nearest_strict = (nearest["fident"] >= MIN_SEQ_ID) & (nearest["qcov"] >= COVERAGE) & (nearest["tcov"] >= COVERAGE)
        result["homolog_any_hit"] = hits.loc[strict, "query"].nunique() / len(test_idx)
        result["homolog_best_hit"] = nearest_strict.mean()
        results.append(result)

        pairs = best.loc[best.index.intersection(test_idx)].copy()
        pairs["agree"] = pairs.index.map(labels).to_numpy() == pairs["target"].map(labels).to_numpy()
        pairs["band"], _ = band_cut(pairs["fident"], BANDS)
        agreement = pairs.groupby("band", observed = False)["agree"].agg(["mean", "size"])
        agreement["split"] = split
        consistency.append(agreement.reset_index())

        kept = nearest[["target", "fident", "qcov", "tcov", "bits"]].copy()
        kept["split"] = split
        kept["nb_label"] = neighbour
        kept["lookup_pred"] = pred
        best_all.append(kept.rename_axis("idx").reset_index())

    table = pd.DataFrame(results).set_index("split")
    agreement = pd.concat(consistency, ignore_index = True)
    agreement.columns = ["band", "label_agreement", "n", "split"]

    log("\nNearest-neighbour lookup: each test protein takes the label of its most similar training protein")
    log("(MMseqs2 search; proteins without a hit get the majority class).")
    for split, label in SPLITS.items():
        row = table.loc[split]
        log(f"  {label:<10} MCC {row['mcc']:.3f}, accuracy {row['accuracy']:.3f}, "
            f"{row['hit_rate']:.1%} of test proteins have a hit")
    log(f"Test proteins with a training hit passing the clustering thresholds: "
        f"{table.loc['split_random', 'homolog_any_hit']:.1%} random, {table.loc['split_cluster', 'homolog_any_hit']:.1%} clustered")
    log("(clustering is greedy, so a few relatives still land in different clusters).")

    log("\nLabel agreement between a test protein and its nearest training protein, by sequence identity (n):")
    log(f"  {'identity':<10}{'random':<17}clustered")
    for band, label in zip(agreement["band"].unique(), band_labels(BANDS)):
        cells = []
        for split in SPLITS:
            row = agreement[(agreement["band"] == band) & (agreement["split"] == split)].iloc[0]
            cells.append(f"{row['label_agreement']:.3f} ({row['n']:,})")
        log(f"  {label:<10}{cells[0]:<17}{cells[1]}")

    table.to_csv(OUT / "results_lookup.csv")
    agreement.to_csv(OUT / "label_consistency.csv", index = False)
    pd.concat(best_all, ignore_index = True).to_parquet(OUT / "best_hits.parquet")

def logistic():
    df = pd.read_parquet(OUT / "splits.parquet")
    comp = pd.read_parquet(OUT / "features_comp.parquet")
    matrices = {"comp": comp.to_numpy(), "kmer": sparse.load_npz(OUT / "features_kmer.npz")}
    y = df["label"].to_numpy()
    groups = df["cluster"].to_numpy()

    results, cv_rows, preds = [], [], {}
    for split in SPLITS:
        train = (df[split] == "train").to_numpy()
        test = ~train
        if split == "split_cluster":
            cv = StratifiedGroupKFold(n_splits = CV_FOLDS, shuffle = True, random_state = SEED)
            fit_args = {"groups": groups[train]}
        else:
            cv = StratifiedKFold(n_splits = CV_FOLDS, shuffle = True, random_state = SEED)
            fit_args = {}

        for model in MODELS:
            if model == "comp":
                estimator = make_pipeline(StandardScaler(), LogisticRegression(max_iter = 5000))
                param = "logisticregression__C"
                grid = C_COMP
            else:
                estimator = LogisticRegression(solver = "liblinear", max_iter = 5000)
                param = "C"
                grid = C_KMER

            grid_search = GridSearchCV(estimator, {param: grid}, scoring = "roc_auc", cv = cv,
                                       n_jobs = N_JOBS, refit = True)
            grid_search.fit(matrices[model][train], y[train], **fit_args)

            best = grid_search.best_params_[param]
            for value, mean, sd in zip(grid, grid_search.cv_results_["mean_test_score"],
                                       grid_search.cv_results_["std_test_score"]):
                cv_rows.append({"split": split, "model": model, "C": value, "cv_auc": mean, "cv_sd": sd})

            prob = grid_search.best_estimator_.predict_proba(matrices[model][test])[:, 1]
            pred = (prob >= 0.5).astype(int)
            result = {"split": split, "model": model, "C": best, "cv_auc": grid_search.best_score_}
            result.update(scores(y[test], pred))
            result.update({"roc_auc": roc_auc_score(y[test], prob),
                           "pr_auc": average_precision_score(y[test], prob)})
            results.append(result)
            preds[f"{split}_{model}_prob"] = prob
            preds[f"{split}_{model}_pred"] = pred

    diag = []
    for split in SPLITS:
        train = (df[split] == "train").to_numpy()
        clf = LogisticRegression(solver = "liblinear", max_iter = 5000, C = 1.0)
        clf.fit(matrices["kmer"][train], y[train])
        prob = clf.predict_proba(matrices["kmer"][train])[:, 1]
        diag.append({"split": split, "n_iter": int(clf.n_iter_.max()), "prob_min": prob.min(),
                     "prob_max": prob.max(), "n_pos_pred": int((prob >= 0.5).sum()),
                     "n_pos_true": int(y[train].sum()), "n": int(train.sum()),
                     "train_auc": roc_auc_score(y[train], prob)})

    table = pd.DataFrame(results).set_index(["split", "model"])
    diag = pd.DataFrame(diag).set_index("split")

    log(f"\nLogistic regression, test ROC-AUC. C is chosen by {CV_FOLDS}-fold cross-validation on the training set,")
    log("grouped by cluster for the clustered split so that tuning never sees relatives of the validation proteins.")
    for model, label in MODELS.items():
        rand, clust = table.loc[("split_random", model)], table.loc[("split_cluster", model)]
        log(f"  {label:<12} random {rand['roc_auc']:.3f} (C = {rand['C']:g}), "
            f"clustered {clust['roc_auc']:.3f} (C = {clust['C']:g})")
    log(f"At the default C = 1 the 3-mer model underfits: it labels only {diag.loc['split_random', 'n_pos_pred']:,} (random) "
        f"and {diag.loc['split_cluster', 'n_pos_pred']:,} (clustered)")
    log(f"of the {diag.loc['split_random', 'n_pos_true']:,} training enzymes as enzymes. 3-mer frequencies are small, "
        f"so the model needs weak regularisation (large C).")

    table.to_csv(OUT / "results_logistic.csv")
    pd.DataFrame(cv_rows).to_csv(OUT / "cv_results.csv", index = False)
    diag.reset_index().to_csv(OUT / "diag_default_C.csv", index = False)
    np.savez(OUT / "preds_logistic.npz", **preds)
