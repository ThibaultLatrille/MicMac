import matplotlib as mpl
import pandas as pd

mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse

my_dpi = 256
fontsize = 14
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input", help="Input file")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    groups = [v for k, v in pd.read_csv(args.input, sep="\t").groupby("Lineage")]
    for ddf in groups:
        plt.plot(ddf["Generation"], ddf["Phenotype mean"])
    plt.ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
    plt.xlabel("time ($t$)", fontsize=fontsize)
    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.clf()

    for ddf in groups:
        plt.plot(ddf["Generation"], ddf["Phenotype var"])
    plt.ylabel("Phenotype variance ($V_G$)", fontsize=fontsize)
    plt.xlabel("time ($t$)", fontsize=fontsize)
    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output.replace(".pdf", "_Vg.pdf"), format="pdf")
    plt.savefig(args.output.replace(".pdf", "_Vg.png"), format="png")
    plt.clf()

    for ddf in groups:
        plt.plot(ddf["Generation"], ddf["Population size"])
    plt.ylabel("Population size ($N_e$)", fontsize=fontsize)
    plt.xlabel("time ($t$)", fontsize=fontsize)
    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output.replace(".pdf", "_Ne.pdf"), format="pdf")
    plt.savefig(args.output.replace(".pdf", "_Ne.png"), format="png")
    plt.clf()
