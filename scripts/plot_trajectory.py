from os.path import basename
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import argparse

mpl.use('Agg')

my_dpi = 128
fontsize = 14
fontsize_legend = 12

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input", help="Input file")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    groups = [v for k, v in pd.read_csv(args.input, sep="\t").groupby("Lineage")]
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='all', figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    fig.suptitle(basename(args.input).replace(".tsv", '').replace("_", ' '))

    for ddf in groups:
        ax1.plot(ddf["Generation"], ddf["Phenotype mean"])
        ax2.plot(ddf["Generation"], ddf["Population size"])
        ax3.plot(ddf["Generation"], ddf["Fitness mean"])
        ax4.plot(ddf["Generation"], ddf["Genotype var"])
        ax5.plot(ddf["Generation"], ddf["Phenotype var"])
        ax6.plot(ddf["Generation"], ddf["Heritability"])

    ax1.set_ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
    ax2.set_ylabel("Population size ($N_e$)", fontsize=fontsize)
    ax3.set_ylabel("Fitness mean ($\\bar{w}$)", fontsize=fontsize)
    ax4.set_ylabel("Genotype variance ($V_G$)", fontsize=fontsize)
    ax5.set_ylabel("Phenotype variance ($V_P$)", fontsize=fontsize)
    ax6.set_ylabel("Heritability ($h^2$)", fontsize=fontsize)

    ax4.set_xlabel("time ($t$)", fontsize=fontsize)
    ax5.set_xlabel("time ($t$)", fontsize=fontsize)
    ax6.set_xlabel("time ($t$)", fontsize=fontsize)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.clf()
