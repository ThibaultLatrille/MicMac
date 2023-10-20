#!/usr/bin/env python3
import os
import argparse
from gzip import open as gzopen
import pandas as pd
from ete3 import Tree
from Bio import Phylo


def open_tree(tree_path: str, format_ete3: int = 1) -> Tree:
    if tree_path.endswith(".gz"):
        newick = gzopen(tree_path).read().decode()
        return Tree(newick, format=format_ete3)
    else:
        return Tree(tree_path, format=format_ete3)


def prune_tree(input_tree: Tree, list_taxa: list = None) -> Tree:
    tree = input_tree.copy()

    # Prune tree
    if list_taxa is not None:
        tree.prune(list_taxa, preserve_branch_length=True)
        assert len(tree.get_leaves()) == len(list_taxa), f"Pruning failed: {len(tree.get_leaves())} != {len(list_taxa)}"

    # Add polytomies if branch length are 0
    remove_nodes = set([n for n in tree.traverse() if (n.dist == 0.0 and not n.is_root())])
    for n in remove_nodes:
        n.delete()
    assert len(
        set([n for n in tree.traverse()]).intersection(remove_nodes)) == 0, "Some polytomies could not be removed"
    for n in tree.traverse():
        if not n.is_root():
            assert n.dist > 0.0, f"Branch length is 0.0 for node {n.name}"
    return tree


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(input_tree: str, input_traits: str, output_tree: str, output_traits: str):
    for path in [input_traits, input_tree]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_traits, output_tree]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Convert newick to nexus including root
    Phylo.convert(input_tree, "newick", output_tree, "nexus")
    traits = pd.read_csv(input_traits, sep="\t")
    with open(output_traits, "w") as f:
        f.write('#NEXUS\n\n')
        f.write("Begin data;\n")
        f.write(f"Dimensions ntax={len(traits)} nchar=2;\n")
        f.write("Format datatype=Continuous missing=? gap=-;\n")
        f.write("Matrix\n")
        for row in traits.itertuples():
            f.write(f"{row.TaxonName}\t{row.Phenotype_mean}\t{row.Genotype_mean}\n")
        f.write(";\n")
        f.write("End;\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tree", help="Input tree file", required=True)
    parser.add_argument("--input_traits", help="Input traits file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    args = parser.parse_args()
    main(args.input_tree, args.input_traits, args.output_tree, args.output_traits)
