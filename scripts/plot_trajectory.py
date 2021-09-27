import matplotlib as mpl
import pandas as pd

mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse

RED = "#EB6231"
BLUE = "#5D80B4"
GREEN = "#8FB03E"

my_dpi = 256
fontsize = 14
fontsize_legend = 12
plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
plt.subplot(1, 1, 1)
plt.ylabel("Mean phenotype ($\\bar{X}$)", fontsize=fontsize)
plt.xlabel("time ($t$)", fontsize=fontsize)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input", help="Input file")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    ddf = df[df["Population"] == 0]
    plt.plot(ddf["Generation"], ddf["Phenotype mean"], color=GREEN)
    ddf = df[df["Population"] == 1]
    plt.plot(ddf["Generation"], ddf["Phenotype mean"], color=BLUE)
    ddf = df[df["Population"] == 2]
    plt.plot(ddf["Generation"], ddf["Phenotype mean"], color=RED)

    plt.xticks(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")
    plt.savefig(args.output.replace(".pdf", ".png"), format="png")

    print('Plot completed')
