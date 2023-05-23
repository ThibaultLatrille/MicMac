import numpy as np
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 18


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
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.close('all')
    print(output)


def hist_plot(x, x_label, output, nbr_genes, title="", xscale="log"):
    output = output.replace("scatter_plot", "histogram")
    color_models = colors(x)
    fig = plt.figure(figsize=(1920 / my_dpi, 480 / my_dpi), dpi=my_dpi)
    fig.suptitle("{0} genes".format(nbr_genes) + title)
    ax = fig.add_subplot(1, 1, 1)
    x = {k: [i for i in v if np.isfinite(i)] for k, v in x.items()}
    min_x = max(1e-6, np.min([min(i) for i in x.values()]))
    max_x = np.max([max(i) for i in x.values()])
    if xscale == "log":
        bins = np.geomspace(min_x, max_x, 100)
    else:
        min_x = np.min([min(i) for i in x.values()])
        bins = np.linspace(min_x, max_x, 100)
    hist, _, _ = ax.hist(x.values(), bins=bins, color=color_models, **hist_filled)
    hist, _, _ = ax.hist(x.values(), bins=bins, color=color_models, **hist_step)
    max_y = 1.2 * (max([max(h) for h in hist]) if len(x) > 1 else max(hist))

    if xscale == "log":
        for id_m, m in enumerate(x):
            x_mean = np.mean(x[m])
            ax.plot((x_mean, x_mean), (0, max_y), linewidth=3, color=color_models[id_m],
                    label=f'{m.replace("_", " ").capitalize()} (mean {x_mean:.2g})')

    ax.set_xlabel(x_label, fontsize=fontsize)
    if xscale == "log":
        ax.set_xlim((0.95 * min_x, 1.05 * max_x))
        ax.set_xscale("log")
    else:
        ax.set_xlim((-0.01, 1.01))
        ax.set_yscale("log")
    # Change x and y ticks font size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    ax.set_ylim((0, max_y))
    ax.set_ylabel("Density", fontsize=fontsize)
    ax.legend(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.savefig(output.replace(".pdf", ".png"), format="png")
    plt.clf()
    print(output)
