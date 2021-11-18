from os.path import basename
import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

mpl.use('Agg')

my_dpi = 256
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
    fig, ax1 = plt.subplots(1, 1, sharex='all', linewidth=0.5, figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
    for ddf in groups:
        plot_trajectory(ax1, ddf["Generation"], ddf["Phenotype mean"])
    ax1.set_ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
    ax1.set_xlabel("Generations", fontsize=fontsize)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")
    plt.clf()

    fig, ax1 = plt.subplots(1, 1, sharex='all', linewidth=0.5, figsize=(1920 / my_dpi, 1440 / my_dpi), dpi=my_dpi)
    for ddf in groups:
        std = 1.96 * np.std(ddf["Phenotype var"])
        ax1.fill_between(ddf["Generation"], ddf["Phenotype mean"] - std, ddf["Phenotype mean"] + std, alpha=0.2)
        plot_trajectory(ax1, ddf["Generation"], ddf["Phenotype mean"])
    ax1.set_ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
    ax1.set_xlabel("Generations", fontsize=fontsize)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output.replace(".pdf", ".fill.pdf"), format="pdf")
    plt.savefig(args.output.replace(".pdf", ".fill.png"), format="png")
    plt.clf()

