import matplotlib as mpl
from glob import glob
import pandas as pd
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse

RED = "#EB6231"
BLUE = "#5D80B4"

my_dpi = 256
fontsize = 14
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)
plt.ylabel("Cov($X_i , X_j$)", fontsize=fontsize)
plt.xlabel(r"$\frac{V_G}{N_e} t$", fontsize=fontsize)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    cov = []
    v_genet = []
    for filepath in glob(args.folder + "/*.tsv"):
        df = pd.read_csv(filepath, sep="\t")
        t_max = max(df["Generation"])
        m1 = float(df[(df["Population"] == 1) & (df["Generation"] == t_max)]["Phenotype mean"])
        m2 = float(df[(df["Population"] == 2) & (df["Generation"] == t_max)]["Phenotype mean"])
        cov.append(np.mean(m1 * m2 - (m1 + m2) ** 2 / 4))
        v_genet.append(np.mean(df[df["Population"] != 0]["Phenotype var"]) * t_max / (100 * 2))

    plt.scatter(v_genet, cov, color=RED)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")

    print('Plot completed')
