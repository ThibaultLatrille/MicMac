import numpy as np
from natsort import natsorted
from glob import glob
from os.path import basename, isdir
from neutrality_index import *
from libraries import *
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


def compute_distance(tree: Tree, trait: str, dico_nuc_d: dict, dico_pheno_d: dict) -> (float, float):
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        leaves = node.get_leaves()
        x = [float(getattr(leaf, trait)) for i, leaf in enumerate(leaves)]
        n = len(leaves)
        list_nuc, list_pheno = [], []
        for i, j in itertools.product(range(n), range(n)):
            if i >= j:
                continue
            d = tree.get_distance(leaves[i], leaves[j]) * 4
            list_nuc.append(d)
            list_pheno.append(np.var([x[i], x[j]]))
        dico_nuc_d[node.name].append(np.mean(list_nuc))
        dico_pheno_d[node.name].append(np.mean(list_pheno))


def saturation_fit(x, a, b):
    return a * (x / (b + x))


def linear_fit(x, a, b):
    return a * x + b


# define function to calculate r-squared
def polyfit(x_array, y_array, func, label, ax, color):
    popt, pcov = curve_fit(func, x_array, y_array, bounds=(0, np.inf))
    y_pred = func(x_array, *popt)
    r2 = r2_score(y_array, y_pred)
    if label == "saturation":
        reg = f"y = {popt[0]:.2g}x / ({popt[1]:.2g} + x) ($r^2$={r2:.2g})"
    else:
        reg = f"y = {popt[0]:.2g}x + {popt[1]:.2g} ($r^2$={r2:.2g})"
    x_line = np.linspace(np.min(x_array), np.max(x_array), 100)
    ax.plot(x_line, func(x_line, *popt), linestyle="--", label=reg, color=color)


def distance_plot(x_dico, y_dico, output):
    x_mean = {k: np.mean(v) for k, v in x_dico.items()}
    y_mean = {k: np.mean(v) for k, v in y_dico.items()}
    y_lower = {k: np.percentile(v, 10) for k, v in y_dico.items()}
    y_upper = {k: np.percentile(v, 90) for k, v in y_dico.items()}
    fig = plt.figure(figsize=(1280 / my_dpi, 640 / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)
    x_array = np.array(list(x_mean.values()))
    y_array = np.array(list(y_mean.values()))
    polyfit(x_array, y_array, saturation_fit, "saturation", ax, "blue")
    polyfit(x_array, y_array, linear_fit, "linear", ax, "red")
    ax.scatter(x_array, y_array, s=12, color="black")
    for k in x_mean.keys():
        ax.errorbar(x_mean[k], y_mean[k], color="black", alpha=0.05, linewidth=0.5,
                    yerr=[[y_mean[k] - y_lower[k]], [y_upper[k] - y_mean[k]]])
    ax.set_xlabel('Nucleotide distance', weight='bold', fontsize=fontsize)
    ax.set_ylabel('Phenotypic distance', weight='bold', fontsize=fontsize)
    ax.legend(fontsize=fontsize_legend)
    # Change x and y ticks font size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.savefig(replace_last(output, ".pdf", ".png"), format="png")
    plt.close('all')
    print(output)


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
        nuc_distance = defaultdict(list)
        pheno_distance = defaultdict(list)
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

            compute_distance(tree, trait="Phenotype_mean", dico_nuc_d=nuc_distance, dico_pheno_d=pheno_distance)
            assert len(nuc_distance) == len(pheno_distance), "Error in the computation of the distances."

        distance_plot(nuc_distance, pheno_distance, replace_last(output, ".tsv.gz", f".distance.{m}.pdf"))
    # Output the dictionary as a dataframe
    dict_of_df = {k: pd.DataFrame(v) for k, v in var_dict.items()}
    df = pd.concat(dict_of_df, axis=1)
    # to break out the lists into columns
    df.to_csv(output, sep="\t", index=False)
    x_str = r"Neutrality index $\left( \rho \right)$"
    ratio_scaled = {m: np.array(var_dict["between"][m]) / np.array(var_dict["within"][m]) for m in models}
    hist_plot(ratio_scaled, x_str, replace_last(output, ".tsv.gz", ".hist.pheno.pdf"))

    ratio_scaled_geno = {m: np.array(var_dict["between_geno"][m]) / np.array(var_dict["within_geno"][m]) for m in
                         models}
    hist_plot(ratio_scaled_geno, x_str, replace_last(output, ".tsv.gz", ".hist.geno.pdf"))

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
