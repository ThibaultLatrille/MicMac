import numpy as np
import pandas as pd
from ete3 import Tree


def total_tree_length(t: Tree):
    return sum([n.dist for n in t.traverse() if not n.is_root()])


def distances_from_root_to_leaves(t: Tree):
    assert t.is_root()
    return [n.get_distance(t) for n in t.get_leaves()]


if __name__ == '__main__':
    var_within_path = "MamGenTimescale.m.var_within.tsv"
    dS_tree_path = "MamGenTimescale.m.tree"
    time_tree_path = "MamGenTimescale.m.species.nwk"

    var_within = pd.read_csv(var_within_path, sep="\t")
    pS = var_within["Nucleotide_diversity"].mean()
    print(f"pS = {pS}")
    dS_tree = Tree(dS_tree_path, format=1)
    dS_total_tree_length = total_tree_length(dS_tree)
    print(f"Tree length (dS) = {dS_total_tree_length}")
    time_tree = Tree(time_tree_path, format=1)
    time_total_tree_length = total_tree_length(time_tree)
    print(f"Tree length (My) = {time_total_tree_length}")
    ages_list = distances_from_root_to_leaves(time_tree)
    root_age = np.mean(ages_list)
    print(f"Root age (My) = {root_age}")
    root_dS = dS_total_tree_length * root_age / time_total_tree_length
    print(f"Root age (dS) = {root_dS}")

    pop_size = 50
    nbr_loci = 10000
    nbr_sites_var = nbr_loci * pS
    print(f"nbr_sites_var = {nbr_sites_var}")
    # pS = 4 * N * u
    # u = pS / (4 * N)
    u = pS / (4 * pop_size)
    print(f"u = {u}")

    nbr_generations = root_dS / u
    print(f"nbr_generations = {nbr_generations}")

