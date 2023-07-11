import os
import argparse
import numpy as np
import pandas as pd
from natsort import natsorted
from glob import glob
from collections import defaultdict
from os.path import basename, isdir
from neutrality_index import brownian_fitting, open_tree
from libraries import hist_plot, scatter_plot


def main(folder, neutral_tree_path, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models_path = {basename(p): p for p in glob(f"{folder}/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))
    assert "neutral" in models

    replicates = {m: natsorted(glob(f"{p}/*.nhx.gz")) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    neutral_tree = open_tree(neutral_tree_path)
    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            tree = open_tree(filepath)
            assert tree.get_leaf_names() == neutral_tree.get_leaf_names()
            for node, node_neutral in zip(tree.traverse("preorder"), neutral_tree.traverse("preorder")):
                assert node.name == node_neutral.name
                if node.is_root():
                    node.dist = 0
                else:
                    node.dist = float(getattr(node_neutral, "d"))

            leaves_value = defaultdict(list)
            for n, n_neutral in zip(tree.get_leaves(), neutral_tree.get_leaves()):
                var_pheno = float(getattr(n, "Phenotype_var")) * float(getattr(n, "BranchMean_Heritability"))
                leaves_value["within"].append(var_pheno / float(getattr(n_neutral, "Theta")))
                leaves_value['within_tajima'].append(var_pheno / float(getattr(n_neutral, "Theta_Tajima")))
                leaves_value['within_watterson'].append(var_pheno / float(getattr(n_neutral, "Theta_Watterson")))
                leaves_value["within_fay_wu"].append(var_pheno / float(getattr(n_neutral, "Theta_Fay_Wu")))

                var_geno = float(getattr(n, "Genotype_var"))
                leaves_value["within_geno"].append(var_geno / float(getattr(n_neutral, "Theta")))

                r = float(getattr(n, "Expected_Mutational_var")) / (2 * float(getattr(n_neutral, "Mutational_rate")))
                leaves_value["mutational"].append(r)

            var_dict['replicate'][m].append(os.path.basename(filepath).replace(".nhx.gz", "").split('_')[-1])
            for k, v in leaves_value.items():
                var_dict[k][m].append(np.mean(v))

            pheno_anc, var_between_pheno = brownian_fitting(tree, trait="Phenotype_mean")
            var_dict['between'][m].append(var_between_pheno)
            geno_anc, var_between_geno = brownian_fitting(tree, trait="Genotype_mean")
            var_dict['between_geno'][m].append(var_between_geno)

    # Output the dictionary as a dataframe
    dict_of_df = {k: pd.DataFrame(v) for k, v in var_dict.items()}
    df = pd.concat(dict_of_df, axis=1)
    # to break out the lists into columns
    df.to_csv(output, sep="\t", index=False)
    x_str = r"Neutrality index $\left( \rho \right)$"
    ratio_scaled = {m: np.array(var_dict["between"][m]) / np.array(var_dict["within"][m]) for m in models}
    hist_plot(ratio_scaled, x_str, output.replace(".tsv.gz", ".hist.pheno.pdf"))

    ratio_scaled_geno = {m: np.array(var_dict["between_geno"][m]) / np.array(var_dict["within_geno"][m]) for m in models}
    hist_plot(ratio_scaled_geno, x_str, output.replace(".tsv.gz", ".hist.geno.pdf"))

    mut_str = r"Variance mutational $\left( \sigma^2_{M} \right)$"
    within_str = r"Variance within $\left( \sigma^2_{W} \right)$"
    between_str = r"Variance between $\left( \sigma^2_{B} \right)$"
    scatter_plot(var_dict["within"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "-H.pheno.pdf"), histy_log=True)
    scatter_plot(var_dict["within_tajima"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_tajima.pheno.pdf"), histy_log=True)
    scatter_plot(var_dict["within_watterson"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_watterson.pheno.pdf"), histy_log=True)
    scatter_plot(var_dict["within_fay_wu"], var_dict["between"], within_str, between_str,
                 output.replace(".tsv.gz", "_fay_Wu.pheno.pdf"), histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["within"], mut_str, within_str,
                 output.replace(".tsv.gz", "_mut_within.pheno.pdf"), histy_log=True)
    scatter_plot(var_dict["mutational"], var_dict["between"], mut_str, between_str,
                 output.replace(".tsv.gz", "_mut_between.pheno.pdf"), histy_log=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-t', '--neutral_tree', required=True, type=str, dest="neutral_tree", help="Input neutral tree")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.neutral_tree, args.output)
