from os.path import basename, isdir
import threading
import argparse
import itertools
from glob import glob
import pandas as pd
import numpy as np
from scipy.stats import hmean, gmean
import statsmodels.api as sm
from ete3 import Tree
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')
hist_filled = {'bins': 'auto', 'density': True, 'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'bins': 'auto', 'density': True, 'histtype': 'step'}
my_dpi = 128
fontsize = 14
fontsize_legend = 12


def colors(x):
    cs = ["#EB6231", "#8FB03E", "#E29D26", "#5D80B4", "#857BA1"]
    return [cs[idx % len(cs)] for idx in range(len(x))]


def get_node_names(root, target1, target2):
    n1 = next(root.iter_search_nodes(name=target1))
    n2 = next(root.iter_search_nodes(name=target2))
    ancestor = root.get_common_ancestor(n1, n2)
    nodes = []
    for n in [n2, n1]:
        current = n
        while current != ancestor:
            nodes.append(current.name)
            current = current.up
    return nodes


def scatter_plot(x, y, x_label, y_label, output, scale, nbr_genes, title=""):
    color_models = colors(x)
    fig = plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    set_nb_points = set([len(nb) for nb in x.values()])
    assert len(set_nb_points) == 1
    nb_points = set_nb_points.pop()
    if scale == "pairs":
        fig.suptitle("{0} genes for each pairs of species ({1} points)".format(nbr_genes, nb_points) + title)
    elif scale == "species":
        fig.suptitle("{0} genes for each species ({1} points)".format(nbr_genes, nb_points) + title)

    gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    max_x = max([np.max(x_i) for x_i in x.values()])
    max_y = max([np.max(x_i) for x_i in y.values()])

    idf = np.linspace(0, max_x, 30)
    ax.plot(idf, idf, '-', linewidth=0.5, color="black")
    for id_m, m in enumerate(x):
        ax.scatter(x[m], y[m], color=color_models[id_m], alpha=0.2)
        results = sm.OLS(y[m], x[m]).fit()
        a = results.params[0]
        linear = a * idf
        reg = '{0} - slope of {1:.2g} ($r^2$={2:.2g})'.format(m.replace("_", " ").capitalize(), a, results.rsquared)
        ax.plot(idf, linear, '-', linestyle="--", label=reg, color=color_models[id_m])

    ax_histx.hist(x.values(), color=color_models, **hist_filled)
    ax_histx.hist(x.values(), color=color_models, **hist_step)
    ax_histy.hist(y.values(), color=color_models, orientation='horizontal', **hist_filled)
    ax_histy.hist(y.values(), color=color_models, orientation='horizontal', **hist_step)
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.set_xlim((0, max_x * 1.05))
    ax.set_ylim((0, max_y * 1.05))
    ax.legend(fontsize=fontsize_legend)
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.clf()
    print(output)


def hist_plot(x, x_label, output, nbr_genes, title=""):
    output = output.replace("scatter_plot", "histogram")
    color_models = colors(x)
    fig = plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    fig.suptitle("{0} genes".format(nbr_genes) + title)
    ax = fig.add_subplot(1, 1, 1)

    hist, _, _ = ax.hist(x.values(), color=color_models, **hist_filled)
    ax.hist(x.values(), color=color_models, **hist_step)
    max_y = 1.2 * max([max(h) for h in hist])

    for id_m, m in enumerate(x):
        x_mean = np.mean(x[m])
        x_median = np.median(x[m])
        ax.plot((x_mean, x_mean), (0, max_y), linewidth=3, color=color_models[id_m],
                label=(m.replace("_", " ").capitalize() + ' (mean {0:.2g}; median {0:.2g})'.format(x_mean, x_median)))
        ax.plot((x_median, x_median), (0, max_y), linewidth=2, linestyle='--', color=color_models[id_m])

    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel("Density", fontsize=fontsize)
    ax.set_xlim((0, 1.05 * max([np.max(x_i) for x_i in x.values()])))
    ax.set_ylim((0, max_y))
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

    Vi, Vg_pair, Vg_scaled_pair, Vg_mean_pair, Vg_harm_pair = {}, {}, {}, {}, {}
    ratio_Vi_Vg_scaled, ratio_Vi_Vg_harm, ratio_Vi_Vg_mean = {}, {}, {}
    Vg, Vm, Vm_theo = {}, {}, {}

    models_path = {basename(p): p for p in glob(args.folder + "/*") if isdir(p)}
    models = [i for i in ["moving_optimum", "directional", "stabilizing", "neutral"] if i in models_path]

    replicates = {m: glob(p + "/*.tsv") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    rep_nhx = {m: set([open(f + ".nhx", 'r').read() for f in g]) for m, g in replicates.items()}
    nhx = set([t_set.pop() for m, t_set in rep_nhx.items() if len(t_set) == 1])
    assert len(nhx) == 1
    tree = Tree(nhx.pop(), format=3)
    nb_leaves = len(tree.get_leaf_names())

    for model, v in itertools.product(models, [Vi, Vg_pair, Vg_scaled_pair, Vg_mean_pair, Vg_harm_pair]):
        v[model] = np.zeros((len(replicates[model]), nb_leaves * (nb_leaves - 1) // 2))

    for model, v in itertools.product(models, [ratio_Vi_Vg_scaled, ratio_Vi_Vg_harm, ratio_Vi_Vg_mean]):
        v[model] = np.nan

    for model, v in itertools.product(models, [Vg, Vm, Vm_theo]):
        v[model] = np.zeros((len(replicates[model]), nb_leaves))


    def threading_replicates(m):
        for f, filepath in enumerate(replicates[m]):
            print(filepath)
            gr = {k: gb for k, gb in pd.read_csv(filepath, sep="\t").groupby("Lineage")}
            gr_leaves = {k: gb for k, gb in gr.items() if k in tree.get_leaf_names()}
            gr_extant = {k: gb[gb["Generation"] == max(gb["Generation"])] for k, gb in gr_leaves.items()}
            gr_last_gens = {k: gb[gb["Generation"] > (max(gb["Generation"]) - 50)] for k, gb in gr_leaves.items()}

            idx = 0
            for i, j in itertools.combinations(gr_leaves, 2):
                dfe = pd.concat([gr_extant[i], gr_extant[j]])
                Vi[m][f, idx] = np.var(dfe["Phenotype mean"], ddof=1)
                Vg_pair[m][f, idx] = np.mean(dfe["Phenotype var"])
                time = tree.get_distance(i, j)
                Vg_scaled_pair[m][f, idx] = np.mean(dfe["Phenotype var"]) * time / (2 * np.mean(dfe["Population size"]))

                dfp = pd.concat([gr[b] for b in get_node_names(tree, i, j)])
                Vg_mean_pair[m][f, idx] = np.mean(dfp["Phenotype var"]) * time / (2 * np.mean(dfp["Population size"]))
                Vg_harm_pair[m][f, idx] = np.mean(dfp["Phenotype var"]) * time / (2 * hmean(dfp["Population size"]))
                idx += 1

            for idx, (i, dfe) in enumerate(gr_extant.items()):
                Vg[m][f, idx] = float(dfe["Phenotype var"])
                Vm[m][f, idx] = float(np.mean(gr_last_gens[i]["Mutational var"]) * dfe["Population size"] * 2)
                Vm_theo[m][f, idx] = float(dfe["Mutational expected var"] * dfe["Population size"] * 2)

        ratio_Vi_Vg_scaled[m] = np.mean(Vi[m] / Vg_scaled_pair[m], axis=1).flatten()
        ratio_Vi_Vg_mean[m] = np.mean(Vi[m] / Vg_mean_pair[m], axis=1).flatten()
        ratio_Vi_Vg_harm[m] = np.mean(Vi[m] / Vg_harm_pair[m], axis=1).flatten()


    threads = list()
    for model in models:
        t = threading.Thread(target=threading_replicates, args=(model,))
        threads.append(t)
        t.start()
    for index, thread in enumerate(threads):
        thread.join()

    replace = args.output.replace
    nb_genes = set([len(v) for v in ratio_Vi_Vg_scaled.values()]).pop()
    if "neutral" in models:
        ratio_Vi_Vg_neutral = {}
        for model in [m for m in models if m != "neutral"]:
            r = (Vi[model] * Vg_pair["neutral"]) / (Vg_pair[model] * Vi["neutral"])
            ratio_Vi_Vg_neutral[model] = gmean(r, axis=1).flatten()
        x_str = r"Variance ratio $\left(\frac{Var[\bar{X^S}]}{Var[\bar{X^N}]}\frac{\bar{V_G^N}}{\bar{V_G^S}}\right)$"
        hist_plot(ratio_Vi_Vg_neutral, x_str, replace(".pdf", "_SN.pdf"), nb_genes)

    for model, d in itertools.product(models, [Vi, Vg_scaled_pair, Vg_mean_pair, Vg_harm_pair, Vg, Vm, Vm_theo]):
        d[model] = d[model].flatten()

    x_str = r"Variance within $\left( \frac{\bar{V_G}}{\bar{N_e}} t \right)$"
    y_str = r"Variance between (Var$[\bar{X}]$)"
    scatter_plot(Vg_scaled_pair, Vi, x_str, y_str, args.output, "pairs", nb_genes, ' - Contemporary data')
    scatter_plot(Vg_harm_pair, Vi, x_str, y_str, replace(".pdf", "_mean.pdf"), "pairs", nb_genes,
                 ' - Mean along ancestry')
    scatter_plot(Vg_mean_pair, Vi, x_str, y_str, replace(".pdf", "_harm.pdf"), "pairs", nb_genes,
                 ' - Harmonic mean along ancestry')
    scatter_plot(Vm, Vg, "$2 N_e V_M$", "$V_G$", replace(".pdf", "_Vg_Vm.pdf"), "species", nb_genes)
    scatter_plot(Vm_theo, Vg, "Expected $2 N_e V_M$", "$V_G$", replace(".pdf", "_Vg_Vme.pdf"), "species", nb_genes)

    x_str = r"Variance ratio $\left( \frac{Var[\bar{X}]}{\bar{V_G}} \frac{\bar{N_e}}{t} \right)$"
    hist_plot(ratio_Vi_Vg_scaled, x_str, args.output, nb_genes, ' - Contemporary data')
    hist_plot(ratio_Vi_Vg_mean, x_str, replace(".pdf", "_mean.pdf"), nb_genes, ' - Mean along ancestry')
    hist_plot(ratio_Vi_Vg_harm, x_str, replace(".pdf", "_harm.pdf"), nb_genes, ' - Harmonic mean along ancestry')
