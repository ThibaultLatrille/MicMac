import os
import argparse
import itertools
from gzip import open as gzopen
from collections import defaultdict
import pandas as pd
import numpy as np
from ete3 import Tree


def open_tree(tree_path: str, format_ete3: int = 1) -> Tree:
    if tree_path.endswith(".gz"):
        newick = gzopen(tree_path).read().decode()
        return Tree(newick, format=format_ete3)
    else:
        return Tree(tree_path, format=format_ete3)


def prune_tree(input_tree: Tree, list_taxa: list = None) -> Tree:
    # Prune tree
    tree = input_tree.copy()
    if list_taxa is not None:
        tree.prune(list_taxa)
    # Add polytomies if branch length are 0
    remove_nodes = set([n for n in tree.traverse() if (n.dist == 0.0 and not n.is_root())])
    for n in remove_nodes:
        n.delete()
    assert len(set([n for n in tree.traverse()]).intersection(remove_nodes)) == 0
    for n in tree.traverse():
        if not n.is_root():
            assert n.dist > 0.0
    return tree


def brownian_fitting(tree: Tree, trait: str) -> (float, float):
    leaves = tree.get_leaves()
    n = len(leaves)
    C = np.zeros((n, n))
    x = np.zeros(n)
    for i, leaf in enumerate(leaves):
        x[i] = float(getattr(leaf, trait))

    for node in tree.traverse("levelorder"):
        if node.is_root():
            node.t = 0
        else:
            # Don't forget the 4 factor (sigma^2_phy = Var(x)/4t)
            node.t = 4 * node.dist + node.up.t

    for i, j in itertools.product(range(n), range(n)):
        if i == j:
            C[i, j] = leaves[i].t
        elif i < j:
            ancestor = leaves[i].get_common_ancestor(leaves[j])
            C[i, j] = ancestor.t
            C[j, i] = ancestor.t

    invC = np.linalg.inv(C)
    ones = np.ones(n)
    v = np.dot(ones.T, invC)
    anc_z = float(np.dot(v, x)) / float(np.dot(v, ones))
    assert np.isfinite(anc_z)

    d = (x - anc_z * ones)
    var = np.dot(d.T, np.dot(invC, d)) / (n - 1)
    return anc_z, var


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(input_traits: str, input_tree: str, input_variance_pop: str, output_tsv: str):
    for path in [input_traits, input_tree, input_variance_pop]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_tsv]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Open tree, branch length should be in substitution per site
    tree = open_tree(input_tree)
    set_taxa = set(tree.get_leaf_names())

    # Open mean trait file
    # The column names should be in the format "*_mean" (where * is the trait name)
    # The first column should be "TaxonName"
    # "*_mean" is mandatory, but can contain missing values (ideally not)
    trait_df = pd.read_csv(input_traits, sep="\t")
    traits = [replace_last(col, "_mean", "") for col in trait_df.columns if col.endswith("_mean")]
    assert len(traits) > 0, f"No column ending with '_mean' found in {input_traits}."
    assert "TaxonName" == trait_df.columns[0], f"Column 'TaxonName' is not the first column in {input_traits}."
    trait_set_taxa = set(trait_df["TaxonName"])
    for taxa in trait_set_taxa:
        assert taxa in set_taxa, f"Taxon {taxa} in file {input_traits} not found in {input_tree}."

    # Open variance file containing the population variance and heritability
    # The column names should be in the format "*_variance" and "*_heritability" (where * is the trait name)
    # The trait names should be the same as in the mean trait file
    # The first column should be "TaxonName"
    # The second column should be "pS", and not contain any missing value
    # "*_variance" is mandatory, but can contain missing values (ideally not)
    # "*_heritability" is optional, if not found, it will be set to 1.0
    var_pop_df = pd.read_csv(input_variance_pop, sep="\t")
    assert "TaxonName" == var_pop_df.columns[0], f"Column 'TaxonName' is not the first column in {input_variance_pop}."
    assert "pS" == var_pop_df.columns[1], f"Column 'pS' is not the second column in {input_variance_pop}."
    for trait in traits:
        assert f"{trait}_variance" in var_pop_df.columns, f"Column {trait}_variance not found in {input_variance_pop}."
    for taxa in var_pop_df["TaxonName"]:
        assert taxa in set_taxa, f"Taxon {taxa} in file {input_variance_pop} not found in {input_tree}."
        assert taxa in trait_set_taxa, f"Taxon {taxa} in file {input_variance_pop} not found in {input_traits}."
    print("The traits found are:")
    print("\t" + "\n\t".join(traits))
    output_dict = defaultdict(list)
    for trait in traits:
        print(f"\nProcessing phenotype {trait}.")
        if f"{trait}_heritability" not in var_pop_df.columns:
            print(f"Warning: column {trait}_heritability not found in {input_variance_pop}.")
            print("Assuming heritability = 1.0.")
            var_pop_df[f"{trait}_heritability"] = 1.0
        notna = np.isfinite(var_pop_df[f"{trait}_variance"])
        # Computing the genetic variance (geno = h² * pheno)
        genetic_variance_array = var_pop_df[f"{trait}_heritability"][notna] * var_pop_df[f"{trait}_variance"][notna]
        variance_pop_array = genetic_variance_array / (var_pop_df["pS"][notna])
        print(f'Found {len(genetic_variance_array)} species with a variance for {trait}.')
        variance_pop = np.mean(variance_pop_array)

        sp_mean_pheno = {sp: v for sp, v in zip(trait_df["TaxonName"], trait_df[f"{trait}_mean"]) if np.isfinite(v)}
        keep_leaf = [leaf.name for leaf in tree.get_leaves() if leaf.name in sp_mean_pheno]
        print(f'Found {len(keep_leaf)} species with a mean {trait}.')
        pruned_tree = prune_tree(tree, keep_leaf)
        assert len(pruned_tree.get_leaves()) == len(keep_leaf), "Error in the pruning of the tree."
        for leaf in pruned_tree.get_leaves():
            assert leaf.name in sp_mean_pheno, f"Leaf {leaf.name} not found in {input_traits}."
            setattr(leaf, trait, sp_mean_pheno[leaf.name])

        # Fitting the brownian model to obtain the phylogenetic variance
        anc_phenotype_mean, variance_phy = brownian_fitting(pruned_tree, trait=trait)
        ratio = variance_phy / variance_pop
        print(f"variance_phy = {variance_phy}")
        print(f"variance_pop = {variance_pop}.")
        print(f"ratio = {ratio}.")
        output_dict["trait"].append(trait)
        output_dict["nbr_taxa_pop"].append(len(genetic_variance_array))
        output_dict["nbr_taxa_phy"].append(len(keep_leaf))
        output_dict["variance_pop"].append(variance_pop)
        output_dict["variance_phy"].append(variance_phy)
        output_dict["ratio"].append(ratio)
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(output_tsv, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--traits", help="Input trait file", required=True)
    parser.add_argument("--variance_pop", help="Input variance_pop file", required=True)
    parser.add_argument("--tree", help="Input tree file", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()
    main(args.traits, args.tree, args.variance_pop, args.output)
