from os.path import basename
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

mpl.use('Agg')

my_dpi = 128
fontsize = 14
fontsize_legend = 12


def plot_trajectory(ax, x, y, mean=False):
    ax.plot(x, y)
    if mean:
        y_mean = np.mean(y)
        ax.plot((np.min(x), np.max(x)), (y_mean, y_mean), linewidth=3, color='black')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input", help="Input file")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    groups = [v for k, v in pd.read_csv(args.input, sep="\t").groupby("Lineage")]
    fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='all', linewidth=0.5,
                                                                     figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    fig.suptitle(basename(args.input).replace(".tsv", '').replace("_", ' '))

    for ddf in groups:
        plot_trajectory(ax1, ddf["Generation"], ddf["Phenotype mean"])
        plot_trajectory(ax2, ddf["Generation"], ddf["Population size"])
        plot_trajectory(ax3, ddf["Generation"], ddf["Fitness mean"])
        plot_trajectory(ax4, ddf["Generation"], ddf["Mutational var"], True)
        plot_trajectory(ax5, ddf["Generation"], ddf["Genotype var"], True)
        plot_trajectory(ax6, ddf["Generation"], ddf["Phenotype var"], True)
        plot_trajectory(ax7, ddf["Generation"], ddf["Heritability"], True)
        plot_trajectory(ax8, ddf["Generation"], ddf["Genotype var"] / ddf["Phenotype var"], True)

    ax1.set_ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
    ax2.set_ylabel("Population size ($N_e$)", fontsize=fontsize)
    ax3.set_ylabel("Fitness mean ($\\bar{w}$)", fontsize=fontsize)
    ax4.set_ylabel("Mutational variance ($V_M$)", fontsize=fontsize)
    ax5.set_ylabel("Genotype variance ($V_G$)", fontsize=fontsize)
    ax6.set_ylabel("Phenotype variance ($V_P$)", fontsize=fontsize)
    ax7.set_ylabel("Heritability parent-offsprings ($h^2$)", fontsize=fontsize)
    ax8.set_ylabel("Heritability ratio ($V_G / V_P$)", fontsize=fontsize)

    ax5.set_xlabel("time ($t$)", fontsize=fontsize)
    ax6.set_xlabel("time ($t$)", fontsize=fontsize)
    ax7.set_xlabel("time ($t$)", fontsize=fontsize)
    ax8.set_xlabel("time ($t$)", fontsize=fontsize)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.clf()
