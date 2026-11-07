from homology.settings import OUT

def log(line = ""):
    print(line, flush = True)
    with open(OUT / "log.txt", "a") as f:
        f.write(line + "\n")
