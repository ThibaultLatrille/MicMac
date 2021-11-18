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
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
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


def scatter_plot(x, y, x_label, y_label, output, nbr_genes, title="", histy_log=False, loglog=False):
    color_models = colors(x)
    fig = plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    set_nb_points = set([len(nb) for nb in x.values()])
    assert len(set_nb_points) == 1
    nb_points = set_nb_points.pop()
    assert nb_points == nbr_genes
    fig.suptitle("{0} genes".format(nbr_genes) + title)

    gs = fig.add_gridspec(2, 2, width_ratios=(7, 2), height_ratios=(2, 7),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    min_x = min([np.min(x_i) for x_i in x.values()])
    max_x = max([np.max(x_i) for x_i in x.values()])
    min_y = min([np.min(y_i) for y_i in y.values()])
    max_y = max([np.max(y_i) for y_i in y.values()])

    idf = np.linspace(min_x, max_x, 30)
    ax.plot(idf, idf, '-', linewidth=0.5, color="black")
    for id_m, m in enumerate(x):
        ax.scatter(x[m], y[m], color=color_models[id_m], alpha=0.1)
        results = sm.OLS(y[m], x[m]).fit()
        a = results.params[0]
        linear = a * idf
        reg = '{0} - slope of {1:.2g} ($r^2$={2:.2g})'.format(m.replace("_", " ").capitalize(), a, results.rsquared)
        ax.plot(idf, linear, '-', linestyle="--", label=reg, color=color_models[id_m])

    bins_x = np.geomspace(max(1e-6, min_x), max_x, 100) if loglog else 100
    bins_y = np.geomspace(max(1e-6, min_y), max_y, 100) if loglog else 100
    ax_histx.hist(x.values(), bins=bins_x, color=color_models, **hist_filled)
    ax_histx.hist(x.values(), bins=bins_x, color=color_models, **hist_step)
    ax_histy.hist(y.values(), bins=bins_y, color=color_models, log=histy_log, orientation='horizontal', **hist_filled)
    ax_histy.hist(y.values(), bins=bins_y, color=color_models, log=histy_log, orientation='horizontal', **hist_step)
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
    else:
        ax.set_xlim((0.95 * min_x, max_x * 1.05))
        ax.set_ylim((0.95 * min_y, max_y * 1.05))
    ax.legend(fontsize=fontsize_legend)
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.close('all')
    print(output)


def hist_plot(x, x_label, output, nbr_genes, title=""):
    output = output.replace("scatter_plot", "histogram")
    color_models = colors(x)
    fig = plt.figure(figsize=(1920 / my_dpi, 880 / my_dpi), dpi=my_dpi)
    fig.suptitle("{0} genes".format(nbr_genes) + title)
    ax = fig.add_subplot(1, 1, 1)
    min_x = max(1e-6, np.min([min(i) for i in x.values()]))
    max_x = np.max([max(i) for i in x.values()])
    logbins = np.geomspace(min_x, max_x, 100)
    hist, _, _ = ax.hist(x.values(), bins=logbins, color=color_models, **hist_filled)
    hist, _, _ = ax.hist(x.values(), bins=logbins, color=color_models, **hist_step)
    max_y = 1.2 * max([max(h) for h in hist])

    for id_m, m in enumerate(x):
        x_mean = np.mean(x[m])
        x_median = np.median(x[m])
        ax.plot((x_mean, x_mean), (0, max_y), linewidth=3, color=color_models[id_m],
                label=(m.replace("_", " ").capitalize() + ' (mean {0:.2g}; median {1:.2g})'.format(x_mean, x_median)))
        ax.plot((x_median, x_median), (0, max_y), linewidth=2, linestyle='--', color=color_models[id_m])

    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_xlim((0.95 * min_x, 1.05 * max_x))
    ax.set_xscale("log")
    ax.set_ylim((0, max_y))
    ax.set_ylabel("Density", fontsize=fontsize)
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

    models_path = {basename(p): p for p in glob(args.folder + "/*") if isdir(p)}
    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: glob(p + "/*.tsv") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    rep_nhx = {m: set([open(f + ".nhx", 'r').read() for f in g]) for m, g in replicates.items()}
    nhx = set([t_set.pop() for m, t_set in rep_nhx.items() if len(t_set) == 1])
    assert len(nhx) == 1
    tree = Tree(nhx.pop(), format=3)
    nb_leaves = len(tree.get_leaf_names())
    nb_pairs = nb_leaves * (nb_leaves - 1) // 2
    distances = {(i, j): tree.get_distance(i, j) for i, j in itertools.combinations(sorted(tree.get_leaf_names()), 2)}

    Vg, Vm, Vm_theo, Vi_scaled, Vg_scaled_pair, Vg_mean_pair, Vg_harm_pair = {}, {}, {}, {}, {}, {}, {}
    for model, v in itertools.product(models, [Vg, Vm, Vm_theo, Vi_scaled, Vg_scaled_pair, Vg_mean_pair, Vg_harm_pair]):
        v[model] = np.zeros(len(replicates[model]))

    Vi, Vg_pair = {}, {}
    for model, v in itertools.product(models, [Vi, Vg_pair]):
        v[model] = np.zeros((len(replicates[model]), nb_pairs))


    def threading_replicates(m):
        for f, filepath in enumerate(replicates[m]):
            print(filepath)
            gr = {k: gb for k, gb in pd.read_csv(filepath, sep="\t").groupby("Lineage")}
            gr_leaves = {k: gb for k, gb in gr.items() if k in tree.get_leaf_names()}
            gr_extant = {k: gb[gb["Generation"] == max(gb["Generation"])] for k, gb in gr_leaves.items()}
            gr_last_gens = {k: gb[gb["Generation"] > (max(gb["Generation"]) - 50)] for k, gb in gr_leaves.items()}

            idx = 0
            for i, j in itertools.combinations(sorted(tree.get_leaf_names()), 2):
                dfe = pd.concat([gr_extant[i], gr_extant[j]])
                v_inter = np.var(dfe["Phenotype mean"], ddof=1)
                v_intra = np.mean(dfe["Phenotype var"])

                Vi[m][f, idx] = v_inter
                Vg_pair[m][f, idx] = v_intra
                idx += 1

                time = distances[(i, j)]
                Vi_scaled[m][f] += v_inter / (time * nb_pairs)
                Vg_scaled_pair[m][f] += v_intra / (2 * np.mean(dfe["Population size"]) * nb_pairs)

                dfp = pd.concat([gr[b] for b in get_node_names(tree, i, j)])
                v_intra_dfp = np.mean(dfp["Phenotype var"]) / nb_pairs
                Vg_mean_pair[m][f] += v_intra_dfp / (2 * np.mean(dfp["Population size"]))
                Vg_harm_pair[m][f] += v_intra_dfp / (2 * hmean(dfp["Population size"]))

            for i, dfe in gr_extant.items():
                Vg[m][f] += float(dfe["Phenotype var"] / (2 * dfe["Population size"] * nb_leaves))
                Vm[m][f] += float(np.mean(gr_last_gens[i]["Mutational var"])) / nb_leaves
                Vm_theo[m][f] += float(dfe["Mutational expected var"]) / nb_leaves

    threads = list()
    for model in models:
        t = threading.Thread(target=threading_replicates, args=(model,))
        threads.append(t)
        t.start()
    for index, thread in enumerate(threads):
        thread.join()

    replace = args.output.replace

    nb_genes = set([len(v) for v in Vg.values()]).pop()
    if "neutral" in models:
        ratio_Vi_Vg_neutral = {}
        for model in [m for m in models if m != "neutral"]:
            r = (Vi[model] * Vg_pair["neutral"]) / (Vg_pair[model] * Vi["neutral"])
            ratio_Vi_Vg_neutral[model] = gmean(r, axis=1).flatten()
        x_str = r"Variance ratio $\left(\frac{Var[\overline{X}] / V_G}{Var[\overline{X}^0] / V_G^0} \right)$"
        hist_plot(ratio_Vi_Vg_neutral, x_str, replace(".pdf", "_SN.pdf"), nb_genes)

    x_str = r"Variance within $\left( \frac{V_P}{2N} \right)$"
    y_str = r"Variance between $\left( \frac{Var[\overline{X}]}{t} \right)$"
    scatter_plot(Vg_scaled_pair, Vi_scaled, x_str, y_str, args.output, nb_genes,
                 title=' - Contemporary data', histy_log=True)
    scatter_plot(Vg_scaled_pair, Vi_scaled, x_str, y_str, replace(".pdf", "_loglog.pdf"), nb_genes,
                 title=' - Contemporary data', loglog=True)

    scatter_plot(Vg_harm_pair, Vi_scaled, x_str, y_str, replace(".pdf", "_mean.pdf"), nb_genes,
                 title=' - Mean along ancestry', histy_log=True)
    scatter_plot(Vg_harm_pair, Vi_scaled, x_str, y_str, replace(".pdf", "_loglog_mean.pdf"), nb_genes,
                 title=' - Mean along ancestry', loglog=True)

    scatter_plot(Vg_mean_pair, Vi_scaled, x_str, y_str, replace(".pdf", "_harm.pdf"), nb_genes,
                 title=' - Harmonic mean along ancestry', histy_log=True)
    scatter_plot(Vg_mean_pair, Vi_scaled, x_str, y_str, replace(".pdf", "_loglog_harm.pdf"), nb_genes,
                 title=' - Harmonic mean along ancestry', loglog=True)

    x_str = r"Variance mutational ($V_M$)"
    y_str = r"Variance within $\left( \frac{V_P}{2N} \right)$"
    scatter_plot(Vm, Vg, x_str, y_str, replace(".pdf", "_Vg_Vm.pdf"), nb_genes)
    scatter_plot(Vm_theo, Vg, x_str, y_str, replace(".pdf", "_Vg_Vme.pdf"), nb_genes)

    x_str = r"Variance ratio $\left( \frac{Var[\overline{X}] 2 N}{V_P t} \right)$"
    ratio_Vi_Vg_scaled = {m: Vi_scaled[m] / Vg_scaled_pair[m] for m in models}
    ratio_Vi_Vg_mean = {m: Vi_scaled[m] / Vg_mean_pair[m] for m in models}
    ratio_Vi_Vg_harm = {m: Vi_scaled[m] / Vg_harm_pair[m] for m in models}
    hist_plot(ratio_Vi_Vg_scaled, x_str, args.output, nb_genes, ' - Contemporary data')
    hist_plot(ratio_Vi_Vg_mean, x_str, replace(".pdf", "_mean.pdf"), nb_genes, ' - Mean along ancestry')
    hist_plot(ratio_Vi_Vg_harm, x_str, replace(".pdf", "_harm.pdf"), nb_genes, ' - Harmonic mean along ancestry')
