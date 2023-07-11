import argparse
from collections import defaultdict
import pandas as pd
from neutrality_index import open_tree, prune_tree


def main(input_path, input_neutral_path, tree_path, traits_path, var_within_path):
    tree = open_tree(input_path, format_ete3=1)

    dico_traits = defaultdict(list)
    for n in tree.get_leaves():
        dico_traits["TaxonName"].append(n.name)
        pheno_mean = float(getattr(n, "Phenotype_mean"))
        dico_traits["Phenotype_mean"].append(pheno_mean if pheno_mean != 0.0 else "NaN")
        geno_mean = float(getattr(n, "Genotype_mean"))
        dico_traits["Genotype_mean"].append(geno_mean)

    df_traits = pd.DataFrame(dico_traits)
    df_traits.to_csv(traits_path, sep=("\t" if traits_path.endswith(".tsv") else ","), index=False)

    neutral_tree = open_tree(input_neutral_path, format_ete3=1)
    assert neutral_tree.get_leaf_names() == tree.get_leaf_names()
    for n in neutral_tree.traverse():
        n.dist = 0. if n.is_root() else float(getattr(n, "d"))

    neutral_tree = prune_tree(neutral_tree)
    # Write the topology to a newick file
    neutral_tree.write(outfile=tree_path, format=3)

    var_within = defaultdict(list)
    for n, n_neutral in zip(tree.get_leaves(), neutral_tree.get_leaves()):
        assert n.name == n_neutral.name
        var_within["TaxonName"].append(n.name)
        var_within["Nucleotide_diversity"].append(float(getattr(n_neutral, "Theta")))
        var_within["Phenotype_variance"].append(float(getattr(n, "Phenotype_var")))
        var_within["Phenotype_heritability"].append(float(getattr(n, "BranchMean_Heritability")))
        var_within["Genotype_variance"].append(float(getattr(n, "Genotype_var")))
        var_within["Genotype_heritability"].append(1.0)
    df = pd.DataFrame(var_within)
    df.to_csv(var_within_path, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tree file", required=True)
    parser.add_argument("--neutral_tree", help="Input neutral tree file", required=True)
    parser.add_argument("--tree", help="Output tree file", required=True)
    parser.add_argument("--traitsfile", help="Output traits file", required=True)
    parser.add_argument("--var_within", help="Output var_within file", required=True)
    args = parser.parse_args()
    main(args.input, args.neutral_tree, args.tree, args.traitsfile, args.var_within)
