# Writeup: how much does train/test homology inflate protein classifiers?

## Question

Proteins come in families of related sequences. If a random split puts close relatives on
both sides, a model can score well by recognising the relative rather than learning anything
that transfers to new families. This project asks how large that inflation is, how it depends
on the type of model, and how it changes with how similar a test protein is to its nearest
training protein.

## Data

Reviewed UniProt (Swiss-Prot) proteins, downloaded 30 Oct 2026 (queries in the README):

- **Enzymes**: 50,160 proteins with experimentally evidenced catalytic activity.
- **Non-enzymes**: 65,708 proteins with no EC number, a function annotation, and protein-level existence evidence.

Label = 1 if the EC number field is non-empty. 1,759 rows from the enzyme file have catalytic
activity but no EC number, so their rule label (0) disagrees with their file. They were dropped
as ambiguous.

Filtering, from 115,868 proteins to **101,304** (44.0% enzymes):

| removed | rows |
|---|---|
| rule label disagrees with source file | 1,759 |
| any non-standard residue (X, B, Z, U, O) | 877 |
| length < 50 | 5,159 |
| length > 1500 | 3,157 |
| exact duplicate sequence with conflicting labels (all copies) | 6 |
| exact duplicate sequence, consistent label (extra copies) | 3,606 |

Sequences are only kept or dropped whole, never edited.

## Method

1. **Clustering.** MMseqs2 `easy-cluster` at 30% identity and 80% coverage of both sequences
   gives 34,043 clusters. 81.6% of proteins have at least one relative in their cluster. The
   clusters are nearly label-pure: 97.9% of multi-member clusters contain only one label.
   Homology is therefore a strong shortcut for this task.
2. **Splits.** The clustered test set is one fold of a `StratifiedGroupKFold` over clusters.
   The random test set is a stratified random sample of exactly the same size (20,261), with
   the same positive fraction (0.440). Any gap between the two comes only from how proteins
   were assigned, not from test size or class balance. Under the random split, 78.9% of test
   proteins share a cluster with training; under the clustered split, none do.
3. **Models**, each trained and tested on both splits:
   - **Lookup baseline**: predict the label of the best MMseqs2 hit in training, or the training majority class if there is no hit. This measures the shortcut directly.
   - **Composition LR**: logistic regression on the 20 amino-acid frequencies (standardised).
   - **3-mer LR**: logistic regression on the 8,000 normalised 3-mer frequencies.

   C is tuned by 3-fold inner cross-validation on the training set only. On the clustered split
   that inner CV is itself grouped by cluster, so tuning does not leak homology either.
4. **Identity bands.** Each test protein is placed in a band by its identity to its best
   training hit (0-30%, 30-40%, 40-50%, 50-70%, 70-90%, 90-100%, or no hit). Performance is
   then reported per band.

## Results

### Whole test set

| model | metric | random | clustered | gap |
|---|---|---|---|---|
| lookup | MCC | 0.953 | 0.799 | 0.154 |
| lookup | hit rate | 0.930 | 0.793 | |
| 3-mer LR | ROC-AUC | 0.886 | 0.816 | 0.070 |
| composition LR | ROC-AUC | 0.759 | 0.743 | 0.016 |

The ordering is what you would expect if homology drives the gap. The lookup baseline is pure
homology and loses the most. 3-mer features carry family information and lose a moderate amount.
Composition is a coarse summary that is mostly blind to family and loses very little.

### By identity to the nearest training protein (3-mer LR, ROC-AUC)

| identity band | random (n) | random AUC | clustered (n) | clustered AUC |
|---|---|---|---|---|
| no hit | 1,410 | 0.814 | 4,201 | 0.835 |
| 0-30% | 894 | 0.769 | 4,258 | 0.809 |
| 30-40% | 2,046 | 0.808 | 5,574 | 0.832 |
| 40-50% | 2,023 | 0.833 | 2,839 | 0.826 |
| 50-70% | 3,851 | 0.869 | 2,445 | 0.765 |
| 70-90% | 4,899 | 0.918 | 718 | 0.880 |
| 90-100% | 5,138 | 0.944 | 226 | 0.823 |

Under the random split, AUC rises steadily with identity. Half of the random test set sits
above 70% identity, against 4.7% of the clustered test set. On dissimilar proteins
(below 50%) the two splits score similarly, with the clustered split slightly higher. The overall gap therefore comes from *who is in the
test set*: the random split is dominated by near-copies of training proteins.

The clustered split still has proteins above 50% identity. MMseqs2 clustering is greedy and
coverage-based, so some similar pairs end up in different clusters. Those bands are small,
and the clustered 0.5-0.7 dip (0.765) is based on fewer proteins. It should not be read as a
real effect without repeating the split with other seeds.

### Label consistency check

For each test protein and its best training hit, how often do the two share a label? On the
random split this rises from 0.943 (below 30% identity) to 0.997 (above 90%). Close relatives
almost always share a label here, which is what makes homology a usable shortcut.

## Earlier experiment: DeepSol solubility (archive/deepsol)

The same design was first run on the DeepSol solubility dataset (69,420 proteins after
cleaning). There, the random and clustered splits scored the same: 3-mer AUC 0.666 vs 0.667,
composition 0.621 vs 0.616. That looked like "no homology inflation", but the label
consistency check showed why. In the random split's >90%-identity band, only **9.9%** of best-hit
pairs shared a label (1,491 pairs), so near-identical sequences mostly had *opposite*
solubility labels. Close relatives were actively misleading, the expected inflation cancelled
out, and the random split had about three times as many of these pairs as the clustered split
(1,491 vs 445). The enzyme task, with consistent labels, was run as a positive control. It
shows that the same method does detect inflation when the labels allow it.

## Limitations

- One clustered split and one random split (seed 42). The gaps are large relative to the
  test size, but there are no confidence intervals. Repeating with several seeds, or
  bootstrapping the test set, would put error bars on them.
- Only linear models on simple features. Deep models (for example protein language models)
  may rely on homology more or less than these do.
- Identity bands use the single best MMseqs2 hit. Proteins with no hit may still have remote
  relatives below MMseqs2's sensitivity.
- The non-enzyme class is defined by the absence of an EC number plus a function annotation.
  Some proteins in it may be uncharacterised enzymes.
