import argparse
from collections import defaultdict
import pandas as pd
from libraries import open_tree


def main(input_path, input_neutral_path, tree_path, traits_path, fossils_path):
    tree = open_tree(input_path, format=1)

    dico_traits = defaultdict(list)
    for n in tree.get_leaves():
        dico_traits["TaxonName"].append(n.name)
        pheno_mean = float(getattr(n, "Phenotype_mean"))
        dico_traits["Phenotype_mean"].append(pheno_mean if pheno_mean != 0.0 else "NaN")
        geno_mean = float(getattr(n, "Genotype_mean"))
        dico_traits["Genotype_mean"].append(geno_mean)

    df_traits = pd.DataFrame(dico_traits)
    df_traits.to_csv(traits_path, sep=("\t" if fossils_path.endswith(".tsv") else ","), index=False)

    dico_fossils = defaultdict(list)
    neutral_tree = open_tree(input_neutral_path, format=1)
    assert neutral_tree.get_leaf_names() == tree.get_leaf_names()
    for n in neutral_tree.traverse():
        n.dist = 0. if n.is_root() else float(getattr(n, "d"))

    # Add polytomies if branch length are 0
    neutral_tree.prune([n.name for n in neutral_tree.traverse() if (n.dist != 0.0 or n.is_root())])
    # Write the topology to a newick file
    neutral_tree.write(outfile=tree_path, format=3)

    assert neutral_tree.is_root()
    farthest, eccentricity = neutral_tree.get_farthest_leaf()
    for n in neutral_tree.traverse():
        dico_fossils["NodeName"].append(n.name)
        age_from_root = n.get_distance(neutral_tree)
        age = max(eccentricity - age_from_root, 0.0)
        if n.is_root():
            assert age == eccentricity
        dico_fossils["Age"].append(age)
        dico_fossils["LowerBound"].append(age)
        dico_fossils["UpperBound"].append(age)
    df_fossils = pd.DataFrame(dico_fossils)
    df_fossils.to_csv(fossils_path, sep=("\t" if fossils_path.endswith(".tsv") else ","), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tree file", required=True)
    parser.add_argument("--neutral", help="Input neutral tree file", required=True)
    parser.add_argument("--tree", help="Output tree file", required=True)
    parser.add_argument("--traits", help="Output traits file", required=True)
    parser.add_argument("--fossils", help="Output fossils file", required=True)
    args = parser.parse_args()
    main(args.input, args.neutral, args.tree, args.traits, args.fossils)
