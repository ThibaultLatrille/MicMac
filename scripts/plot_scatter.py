from os.path import basename, isdir
import threading
import argparse
import itertools
from gzip import open as gzopen
from glob import glob
import pandas as pd
from scipy.stats import hmean, gmean
from ete3 import Tree
from libraries import *


def main(folder, output):
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p)}
    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: glob(f"{p}/*.tsv.gz") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    rep_nhx = {m: [Tree(gzopen(f.replace('.tsv.gz', '.nhx.gz')).read().decode(), format=1) for f in g] for m, g in
               replicates.items()}
    # Get first tree to get the number of genes
    for m, g in rep_nhx.items():
        assert len(set(["-".join(sorted(t.get_leaf_names())) for t in g])) == 1
        assert len(set([t.get_distance(t.get_tree_root()) for t in g])) == 1

    tree = rep_nhx[models[0]][0]
    nb_leaves = len(tree.get_leaf_names())
    nb_pairs = nb_leaves * (nb_leaves - 1) // 2
    print(f"Number of leaves: {nb_leaves}")
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
            print(f"Number of group: {len(gr)}")
            print(' '.join(sorted(gr.keys())))
            gr_leaves = {k: gb for k, gb in gr.items() if k in tree.get_leaf_names()}
            print(f"Number of leaves group: {len(gr_leaves)}")
            print(' '.join(sorted(gr_leaves.keys())))
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

    replace = output.replace

    nb_genes = set([len(v) for v in Vg.values()]).pop()
    non_neutral_models = [m for m in models if m != "neutral"]
    if "neutral" in models and len(non_neutral_models) > 0:
        ratio_Vi_Vg_neutral = {}
        for model in [m for m in models if m != "neutral"]:
            r = (Vi[model] * Vg_pair["neutral"]) / (Vg_pair[model] * Vi["neutral"])
            ratio_Vi_Vg_neutral[model] = gmean(r, axis=1).flatten()
        x_str = r"Variance ratio $\left(\frac{Var[\overline{X}] / V_G}{Var[\overline{X}^0] / V_G^0} \right)$"
        hist_plot(ratio_Vi_Vg_neutral, x_str, replace(".pdf", "_SN.pdf"), nb_genes)

    x_str = r"Variance within $\left( \frac{V_P}{2N} \right)$"
    y_str = r"Variance between $\left( \frac{Var[\overline{X}]}{t} \right)$"
    scatter_plot(Vg_scaled_pair, Vi_scaled, x_str, y_str, output, nb_genes,
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
    hist_plot(ratio_Vi_Vg_scaled, x_str, output, nb_genes, ' - Contemporary data')
    hist_plot(ratio_Vi_Vg_mean, x_str, replace(".pdf", "_mean.pdf"), nb_genes, ' - Mean along ancestry')
    hist_plot(ratio_Vi_Vg_harm, x_str, replace(".pdf", "_harm.pdf"), nb_genes, ' - Harmonic mean along ancestry')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
