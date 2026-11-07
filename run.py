from homology.settings import OUT
from homology.load import load
from homology.cluster import cluster
from homology.split import split
from homology.features import features
from homology.models import lookup, logistic
from homology.bands import by_band, plot

if __name__ == "__main__":
    OUT.mkdir(exist_ok = True)
    (OUT / "log.txt").unlink(missing_ok = True)
    load()
    cluster()
    split()
    features()
    lookup()
    logistic()
    by_band()
    plot()
