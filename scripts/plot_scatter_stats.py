import argparse
import itertools
from glob import glob
from gzip import open as gzopen
from collections import defaultdict
from os.path import basename, isdir
import numpy as np
from ete3 import Tree
import statsmodels.api as sm
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
        ax.plot(idf, linear, linestyle="--", label=reg, color=color_models[id_m])

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


def main(folder, output):
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p)}
    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: glob(p + "/*.nhx.gz") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            newick = gzopen(filepath).read().decode()
            tree = Tree(newick, format=1)
            combinations, weights = defaultdict(list), defaultdict(list)
            for i, j in itertools.combinations(tree.get_leaves(), 2):
                means = [float(getattr(n, "Phenotype_mean")) for n in [i, j]]
                var_x = np.var(means, ddof=1)

                ancestor = tree.get_common_ancestor(i, j)
                if not ancestor.is_root():
                    continue
                d = 0.0
                for n in [i, j]:
                    current = n
                    while current != ancestor:
                        d += float(getattr(n, "d"))
                        current = current.up
                combinations['between'].append(var_x / (2.0 * d))
                weights['d'].append(d)
                weights['1/d'].append(1.0 / d)
            for n in tree.get_leaves():
                vp = float(getattr(n, "Phenotype_var"))
                combinations['within_tajima'].append(vp / float(getattr(n, "Theta_Tajima")))
                combinations['within_watterson'].append(vp / float(getattr(n, "Theta_Watterson")))
                combinations["within_fay_wu"].append(vp / float(getattr(n, "Theta_Fay_Wu")))
                combinations["within_sampled"].append(vp / float(getattr(n, "Theta")))
                r = float(getattr(n, "Mutational_var")) / (2 * float(getattr(n, "Mutational_rate")))
                combinations["mutational"].append(r)

            for k, v in combinations.items():
                var_dict[k][m].append(np.mean(v))
            var_dict['between_d'][m].append(np.average(combinations['between'], weights=weights['d']))
            var_dict['between_1/d'][m].append(np.average(combinations['between'], weights=weights['1/d']))
    nb_genes = set([len(v) for v in var_dict["between"].values()]).pop()
    mut_str = r"Variance mutational $\left( \frac{V_M}{2u} \right)$"
    within_str = r"Variance within $\left( \frac{V_P}{\theta} \right)$"
    between_str = r"Variance between $\left( \frac{V_B}{2d} \right)$"
    scatter_plot(var_dict["within_sampled"], var_dict["between"], within_str, between_str,
                 output, nb_genes, title=' - Sample theta', histy_log=True)
    scatter_plot(var_dict["within_sampled"], var_dict["between_d"], within_str, between_str,
                 output.replace(".pdf", "_d.pdf"), nb_genes, title=' - w=q', histy_log=True)
    scatter_plot(var_dict["within_sampled"], var_dict["between_1/d"], within_str, between_str,
                 output.replace(".pdf", "_d-1.pdf"), nb_genes, title=' - w=1/d', histy_log=True)
    scatter_plot(var_dict["within_tajima"], var_dict["between"], within_str, between_str,
                 output.replace(".pdf", "_tajima.pdf"), nb_genes, title=' - Tajima', histy_log=True)
    scatter_plot(var_dict["within_watterson"], var_dict["between"], within_str, between_str,
                 output.replace(".pdf", "_watterson.pdf"), nb_genes, title=' - Watterson', histy_log=True)
    scatter_plot(var_dict["within_fay_wu"], var_dict["between"], within_str, between_str,
                 output.replace(".pdf", "_fay_Wu.pdf"), nb_genes, title=' - Fay Wu', histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["within_sampled"], mut_str, within_str,
                 output.replace(".pdf", "_mut_within.pdf"), nb_genes, title=' - Vm versus Vp', histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["between"], mut_str, between_str,
                 output.replace(".pdf", "_mut_between.pdf"), nb_genes, title=' - Vm versus Vd', histy_log=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
