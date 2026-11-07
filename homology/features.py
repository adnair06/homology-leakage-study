from itertools import product
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

from homology.settings import OUT, AMINO_ACIDS

MODELS = {"comp": "composition", "kmer": "3-mer"}

def features():
    df = pd.read_parquet(OUT / "splits.parquet")

    counts = CountVectorizer(analyzer = "char", ngram_range = (1, 1), lowercase = False,
                             vocabulary = list(AMINO_ACIDS)).transform(df["seq"])
    comp = pd.DataFrame(normalize(counts, norm = "l1").toarray(), columns = list(AMINO_ACIDS), index = df.index)

    kmers = ["".join(p) for p in product(AMINO_ACIDS, repeat = 3)]
    counts = CountVectorizer(analyzer = "char", ngram_range = (3, 3), lowercase = False,
                             vocabulary = kmers).transform(df["seq"])
    kmer = normalize(counts.astype(np.float64), norm = "l1").tocsr()

    comp.to_parquet(OUT / "features_comp.parquet")
    sparse.save_npz(OUT / "features_kmer.npz", kmer)
