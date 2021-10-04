import os.path
import argparse
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
    max_x = max([max(i) for i in x.values()])
    max_y = max([max(i) for i in y.values()])
    idf = np.linspace(0, max_x, 30)
    plt.plot(idf, idf, '-', linewidth=0.5, color="black")
    for id_m, m in enumerate(x.keys()):
        color = colors[id_m % len(colors)]
        plt.scatter(x[m], y[m], color=color)
        results = sm.OLS(y[m], x[m]).fit()
        a = results.params[0]
        linear = a * idf
        reg = '{0} - slope of {1:.2g} ($r^2$={2:.2g})'.format(m, a, results.rsquared)
        plt.plot(idf, linear, '-', linestyle="--", label=reg, color=color)
    plt.xlabel(x_label, fontsize=fontsize)
    plt.ylabel(y_label, fontsize=fontsize)
    plt.xlim((0, max_x * 1.05))
    plt.ylim((0, max_y * 1.05))
    plt.xticks(fontsize=fontsize_legend)
    plt.legend(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.clf()


if __name__ == '__main__':
    population_size = 100
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()

    plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    plt.subplot(1, 1, 1)

    d_var_between, d_var_scaled_within, d_var_within, d_var_mutational, d_var_e_mutational = {}, {}, {}, {}, {}
    for model_path in glob(args.folder + "/*"):
        if not os.path.isdir(model_path):
            continue
        model = os.path.basename(model_path).replace("_", " ").capitalize()

        var_between, var_within, var_mutational, var_expected_mutational, rep_t_max_array = [], [], [], [], []
        for filepath in glob(model_path + "/*.tsv"):
            df = pd.read_csv(filepath, sep="\t")
            t_max = max(df["Generation"])
            populations = list(sorted([i for i in set(df["Population"]) if i != 0]))
            df_pops = [df[(df["Population"] == pop) & (df["Generation"] == t_max)] for pop in populations]
            var_between.append(np.var([ddf["Phenotype mean"] for ddf in df_pops]))
            var_within.append([ddf["Phenotype var"] for ddf in df_pops])
            var_mutational.append(np.mean([ddf["Mutational var"] for ddf in df_pops]))
            var_expected_mutational.append(np.mean([ddf["Mutational expected var"] for ddf in df_pops]))
            rep_t_max_array.append(t_max - max(df[df["Population"] == 0]["Generation"]))

        t_max = set(rep_t_max_array)
        assert len(t_max) == 1
        t_max = t_max.pop()
        d_var_within[model] = np.array(var_within).mean(axis=1)
        d_var_scaled_within[model] = np.array(var_within).mean(axis=1) * t_max / population_size
        d_var_between[model] = np.array(var_between)
        d_var_mutational[model] = np.array(var_mutational) * population_size * 2
        d_var_e_mutational[model] = np.array(var_expected_mutational) * population_size * 2

    scatter_plot(d_var_scaled_within, d_var_between, "Variance within $\\left( \\frac{V_G}{N_e} t \\right)$",
                 "Variance between (Var$[\\bar{X}]$)", args.output)

    scatter_plot(d_var_mutational, d_var_within, "$2 N_e V_M$", "$V_G$",
                 args.output.replace(".pdf", "_Vg_Vm.pdf"))

    scatter_plot(d_var_e_mutational, d_var_within, "Expected $2 N_e V_M$", "$V_G$",
                 args.output.replace(".pdf", "_Vg_Vme.pdf"))
