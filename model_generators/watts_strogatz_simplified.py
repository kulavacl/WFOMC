import numpy as np
from sys import argv

if len(argv) < 2 or "-h" in argv:
    print("This script creates file fr_sm.wfomcs that encodes the friends and smokers experiment.\n" \
    "argv[1] = domain size. Domain size % 3 must be 0.\n" \
    "optional: -w, set weight, default = 0.8")
    exit(0)

domain_size = int(argv[1])

weight = 0.8
if "-w" in argv:
    weight = float(argv[argv.index("-w") + 1])

with open("watts_strogatz_simplified.wfomcs", "w") as f:
    f.write("\\forall X: (~edge(X, X)) &\n")
    f.write("\\forall X: (\\forall Y: (edge(X, Y) -> edge(Y, X))) &\n")
    f.write("\\forall X: (\\forall Y: (weight(X, Y) <-> ((smokes(X) & edge(X, Y))-> smokes(Y))))\n\n")

    f.write("V = {")
    for i in range(domain_size - 1):
        f.write(f"v{i}, ")
    f.write("v" + f"{domain_size - 1}" + "}\n\n")

    f.write(f"{np.exp(weight)} 1 weight\n\n")

    f.write(f"|edge| = {domain_size * 2 + (domain_size + (1 if domain_size  % 2 == 1 else 0))}\n\n")

    
    for i in range(0, domain_size - 1):
        f.write(f"edge(v{i}, v{i + 1})\n")
    f.write(f"edge(v{domain_size - 1}, v{0})\n")
    
    #uncomment the following for the cliques model
    """
    for i in range(0, domain_size, 3):
        f.write(f"edge(v{i}, v{i + 1})\n")
        f.write(f"edge(v{i + 1}, v{i})\n")
        f.write(f"edge(v{i + 1}, v{i + 2})\n")
        f.write(f"edge(v{i + 2}, v{i + 1})\n")
        f.write(f"edge(v{i}, v{i + 2})\n")
        f.write(f"edge(v{i + 2}, v{i})\n")
    """