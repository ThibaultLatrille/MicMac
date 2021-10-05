import os.path
import argparse
import itertools
from glob import glob
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')

colors = ["#5D80B4", "#E29D26", "#8FB03E", "#EB6231", "#857BA1"]

my_dpi = 128
fontsize = 14
fontsize_legend = 12


def scatter_plot(x, y, x_label, y_label, output):
    fig = plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    max_x = max([max(i) for i in x.values()])
    max_y = max([max(i) for i in y.values()])

    idf = np.linspace(0, max_x, 30)
    ax.plot(idf, idf, '-', linewidth=0.5, color="black")
    for id_m, m in enumerate(x.keys()):
        color = colors[id_m % len(colors)]
        ax.scatter(x[m], y[m], color=color)
        results = sm.OLS(y[m], x[m]).fit()
        a = results.params[0]
        linear = a * idf
        reg = '{0} - slope of {1:.2g} ($r^2$={2:.2g})'.format(m, a, results.rsquared)
        ax.plot(idf, linear, '-', linestyle="--", label=reg, color=color)
        ax_histx.hist(x[m], bins=25, color=color, alpha=0.3)
        ax_histy.hist(y[m], bins=25, orientation='horizontal', color=color, alpha=0.3)
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.set_xlim((0, max_x * 1.05))
    ax.set_ylim((0, max_y * 1.05))
    ax.legend(fontsize=fontsize_legend)
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.clf()
    print(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    d_var_between, d_var_scaled_within, d_var_within, d_var_mutational, d_var_e_mutational = {}, {}, {}, {}, {}
    for model_path in glob(args.folder + "/*"):
        if not os.path.isdir(model_path):
            continue
        model = os.path.basename(model_path).replace("_", " ").capitalize()

        var_between, var_within, var_scaled_within, var_mutational, var_expected_mutational = [], [], [], [], []
        for filepath in glob(model_path + "/*.tsv"):
            print(filepath)
            gr = [v for k, v in pd.read_csv(filepath, sep="\t").groupby("Lineage") if k != 0]
            gr_extant = [v[v["Generation"] == max(v["Generation"])] for v in gr]

            for i, j in itertools.combinations(range(len(gr)), 2):
                df_extant = pd.concat([gr_extant[i], gr_extant[j]])
                var_between.append(np.var(df_extant["Phenotype mean"], ddof=1))

                ddf = pd.concat([gr[i], gr[j]])
                delta_t = max(ddf["Generation"]) - min(ddf["Generation"])
                var_scaled_within.append(np.mean(ddf["Phenotype var"]) * delta_t / np.mean(ddf["Population size"]))

            gr_last_gens = [v[v["Generation"] > (max(v["Generation"]) - 50)] for v in gr]
            for i, df_extant in enumerate(gr_extant):
                var_within.append(df_extant["Phenotype var"])
                var_mutational.append(np.mean(gr_last_gens[i]["Mutational var"]) * df_extant["Population size"] * 2)
                var_expected_mutational.append(df_extant["Mutational expected var"] * df_extant["Population size"] * 2)

        d_var_within[model] = np.array(var_within)
        d_var_scaled_within[model] = np.array(var_scaled_within)
        d_var_between[model] = np.array(var_between)
        d_var_mutational[model] = np.array(var_mutational)
        d_var_e_mutational[model] = np.array(var_expected_mutational)

    scatter_plot(d_var_scaled_within, d_var_between,
                 "Variance within $\\left( \\frac{\\bar{V_G}}{\\bar{N_e}} t \\right)$",
                 "Variance between (Var$[\\bar{X}]$)", args.output)

    scatter_plot(d_var_mutational, d_var_within, "$2 N_e V_M$", "$V_G$",
                 args.output.replace(".pdf", "_Vg_Vm.pdf"))

    scatter_plot(d_var_e_mutational, d_var_within, "Expected $2 N_e V_M$", "$V_G$",
                 args.output.replace(".pdf", "_Vg_Vme.pdf"))
