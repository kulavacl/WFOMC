import numpy as np
from sys import argv

if len(argv) < 2 or "-h" in argv:
    print("This script creates file fr_sm.wfomcs that encodes the friends and smokers experiment.\n" \
    "argv[1] = domain size. Domain size % 3 must be 0.\n" \
    "optional: -beta [f], set beta, default 0.5\n" \
    "optional: -w, set weight, default = 0.8")
    exit(0)

beta = 0.5
if "-beta" in argv:
    beta = float(argv[argv.index("-beta") + 1])

weight = 0.8
if "-w" in argv:
    weight = float(argv[argv.index("-w") + 1])

domain_size = int(argv[1])

with open("watts_strogatz.wfomcs", "w") as f:
    f.write("\\forall X: (~wired_edge(X, X)) &\n")
    f.write("\\forall X: (~edge(X, X)) &\n")
    f.write("\\forall X: (~fake_edge(X, X)) &\n")
    f.write("\\forall X: (~beta(X, X)) &\n")
    f.write("\\forall X: (weight(X, X)) &\n")
    f.write("\\forall X: (~neg_beta(X, X)) &\n")
    f.write("\\forall X: (\\forall Y: (~wired_edge(X, Y) | ~wired_edge(Y, X))) &\n")
    f.write("\\forall X: (\\forall Y: (beta(X, Y) <-> (wired_edge(X, Y) & ~fake_edge(X, Y))))&\n")
    f.write("\\forall X: (\\forall Y: (neg_beta(X, Y) <-> (wired_edge(X, Y) & (fake_edge(X, Y)))))&\n")
    f.write("\\forall X: (\\exists_{=1} Y: (wired_edge(X, Y)))&\n")
    f.write("\\forall X: (\\forall Y: (edge(X, Y) -> edge(Y, X))) &\n")

    f.write("\\forall X: (\\forall Y: (weight(X, Y) <-> ((smokes(X) & edge(X, Y)) -> smokes(Y))))")

    f.write("V = {")
    for i in range(domain_size - 1):
        f.write(f"v{i}, ")
    f.write("v" + f"{domain_size - 1}" + "}\n\n")

    f.write(f"{np.exp(beta)} 1 beta \n")
    f.write(f"{np.exp(1 - beta)} 1 neg_beta \n\n")
    f.write(f"{np.exp(weight)} 1 weight \n\n")

    #for i in range(0, domain_size - 1):
    #    f.write(f"fake_edge(v{i}, v{i + 1})\n")
    #f.write(f"fake_edge(v{domain_size - 1}, v{0})\n")

    #uncomment the following for the cliques watts strogatz
    """
    for i in range(0, domain_size, 3):
        f.write(f"fake_edge(v{i}, v{i + 1})\n")
        f.write(f"fake_edge(v{i + 1}, v{i})\n")
        f.write(f"fake_edge(v{i + 1}, v{i + 2})\n")
        f.write(f"fake_edge(v{i + 2}, v{i + 1})\n")
        f.write(f"fake_edge(v{i}, v{i + 2})\n")
        f.write(f"fake_edge(v{i + 2}, v{i})\n")
    """
    
    for i in range(0, domain_size - 1):
        f.write(f"fake_edge(v{i}, v{i + 1})\n")
    f.write(f"fake_edge(v{domain_size - 1}, v{0})\n")
    #f.write(f"fake_edge(v{domain_size - 1}, v{1})\n")
    
    f.write("\n[fake_edge]")
