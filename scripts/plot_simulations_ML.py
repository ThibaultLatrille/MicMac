import os
import argparse
import numpy as np
import pandas as pd
from glob import glob
from collections import defaultdict
from os.path import basename, isdir
from neutrality_index import brownian_fitting, open_tree
from libraries import hist_plot, scatter_plot


def main(folder, output):
    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))
    assert "neutral" in models

    replicates = {m: glob(p + "/*.nhx.gz") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    neutral_tree = open_tree(
        max(replicates["neutral"], key=lambda x: int(os.path.basename(x).replace(".", "_").split("_")[1])))
    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            tree = open_tree(filepath)
            assert tree.get_leaf_names() == neutral_tree.get_leaf_names()
            for node, node_neutral in zip(tree.traverse("preorder"), neutral_tree.traverse("preorder")):
                assert node.name == node_neutral.name
                node.dist = float(getattr(node_neutral, "d"))

            leaves_value = defaultdict(list)
            for n, n_neutral in zip(tree.get_leaves(), neutral_tree.get_leaves()):
                vg = float(getattr(n, "Phenotype_var")) * abs(float(getattr(n, "Heritability")))
                leaves_value['within_tajima'].append(vg / float(getattr(n_neutral, "Theta_Tajima")))
                leaves_value['within_watterson'].append(vg / float(getattr(n_neutral, "Theta_Watterson")))
                leaves_value["within_fay_wu"].append(vg / float(getattr(n_neutral, "Theta_Fay_Wu")))
                leaves_value["within_sampled"].append(vg / float(getattr(n_neutral, "Theta")))
                r = float(getattr(n, "Expected_Mutational_var")) / (2 * float(getattr(n_neutral, "Mutational_rate")))
                leaves_value["mutational"].append(r)

            var_dict['replicate'][m].append(os.path.basename(filepath).replace(".nhx.gz", "").split('_')[-1])
            for k, v in leaves_value.items():
                var_dict[k][m].append(np.mean(v))

            z, var_between = brownian_fitting(tree, trait="Phenotype_mean")
            var_dict['between'][m].append(var_between)
    # Output the dictionary as a dataframe

    dict_of_df = {k: pd.DataFrame(v) for k, v in var_dict.items()}
    df = pd.concat(dict_of_df, axis=1)
    # to break out the lists into columns
    df.to_csv(output, sep="\t", index=False)
    nb_genes = set([len(v) for v in var_dict["between"].values()]).pop()
    x_str = r"Neutrality index $\left( \frac{\sigma^2_{Phy}}{\sigma^2_{Pop}} \right)$"
    ratio_scaled = {m: np.array(var_dict["between"][m]) / np.array(var_dict["within_sampled"][m]) for m in models}
    hist_plot(ratio_scaled, x_str, output.replace(".tsv.gz", ".hist.pdf"), nb_genes)

    mut_str = r"Variance mutational $\left( \frac{V_M}{2u} \right)$"
    within_str = r"Variance within $\left( \frac{V_G}{\theta} \right)$"
    between_str = r"Variance between $\left( \frac{Var[\overline{X}_t]}{4d} \right)$"
    scatter_plot(var_dict["within_sampled"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "-H.pdf"), nb_genes, title=' - Sample theta', histy_log=True)
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
