import os
import argparse
import itertools
import pandas as pd
from glob import glob
from gzip import open as gzopen
from collections import defaultdict
from os.path import basename, isdir
from ete3 import Tree
from libraries import *


def main(folder, output):
    os.makedirs(os.path.basename(output), exist_ok=True)
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
                r = float(getattr(n, "Expected_Mutational_var")) / (2 * float(getattr(n, "Mutational_rate")))
                combinations["mutational"].append(r)

            var_dict['replicate'][m].append(os.path.basename(filepath).replace(".nhx.gz", "").split('_')[-1])
            for k, v in combinations.items():
                var_dict[k][m].append(np.mean(v))
            var_dict['between_d'][m].append(np.average(combinations['between'], weights=weights['d']))
            var_dict['between_1/d'][m].append(np.average(combinations['between'], weights=weights['1/d']))
    # Output the dictionary as a dataframe

    dict_of_df = {k: pd.DataFrame(v) for k, v in var_dict.items()}
    df = pd.concat(dict_of_df, axis=1)
    # to break out the lists into columns
    df.to_csv(output, sep="\t", index=False)
    nb_genes = set([len(v) for v in var_dict["between"].values()]).pop()
    mut_str = r"Variance mutational $\left( \frac{V_M}{2u} \right)$"
    within_str = r"Variance within $\left( \frac{V_P}{\theta} \right)$"
    between_str = r"Variance between $\left( \frac{V_B}{2d} \right)$"
    scatter_plot(var_dict["within_sampled"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "-H.pdf"), nb_genes, title=' - Sample theta', histy_log=True)
    scatter_plot(var_dict["within_sampled"], var_dict["between_d"], within_str, between_str,
                 output.replace(".tsv.gz", "_Hwq.pdf"), nb_genes, title=' - w=q', histy_log=True)
    scatter_plot(var_dict["within_sampled"], var_dict["between_1/d"], within_str, between_str,
                 output.replace(".tsv.gz", "_Hwd-1.pdf"), nb_genes, title=' - w=1/d', histy_log=True)
    scatter_plot(var_dict["within_tajima"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_tajima.pdf"), nb_genes, title=' - Tajima', histy_log=True)
    scatter_plot(var_dict["within_watterson"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_watterson.pdf"), nb_genes, title=' - Watterson', histy_log=True)
    scatter_plot(var_dict["within_fay_wu"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_fay_Wu.pdf"), nb_genes, title=' - Fay Wu', histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["within_sampled"], mut_str, within_str,
                 output.replace(".tsv.gz", "_mut_within.pdf"), nb_genes, title=' - Vm versus Vp', histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["between"], mut_str, between_str,
                 output.replace(".tsv.gz", "_mut_between.pdf"), nb_genes, title=' - Vm versus Vd', histy_log=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
