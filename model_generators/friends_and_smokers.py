import numpy as np
from sys import argv

if len(argv) < 2 or "-h" in argv:
    print("This script creates file fr_sm.wfomcs that encodes the friends and smokers experiment.\n" \
    "argv[1] = domain size. Domain size % 3 must be 0.")
    exit(0)

domain_size = int(argv[1])

with open("fr_sm.wfomcs", "w") as f:
    f.write("\\forall X: (~friends(X, X)) &\n")
    f.write("\\forall X: (\\forall Y: (friends(X, Y) -> friends(Y, X))) &\n")
    f.write("\\forall X: (\\forall Y: ((friends(X, Y) & smokes(X)) -> smokes(Y))) \n\n")

    f.write("V = {")
    for i in range(domain_size - 1):
        f.write(f"v{i}, ")
    f.write("v" + f"{domain_size - 1}" + "}\n\n")

    for i in range(0, domain_size, 3):
        f.write(f"friends(v{i}, v{i + 1})\n")
        f.write(f"friends(v{i + 1}, v{i})\n")
        f.write(f"friends(v{i + 1}, v{i + 2})\n")
        f.write(f"friends(v{i + 2}, v{i + 1})\n")
        f.write(f"friends(v{i}, v{i + 2})\n")
        f.write(f"friends(v{i + 2}, v{i})\n")

