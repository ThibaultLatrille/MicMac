import argparse
from collections import defaultdict
from gzip import open as gzopen
import numpy as np
import pandas as pd
from ete3 import Tree


def main(input_path, tree_path, traits_path, fossils_path):
    newick = gzopen(input_path).read().decode()
    tree = Tree(newick, format=1)

    # Write the topology to a newick file
    tree.write(outfile=tree_path, format=6)

    dico_traits = defaultdict(list)
    for n in tree.get_leaves():
        dico_traits["TaxonName"].append(n.name)
        pheno_mean = float(getattr(n, "Phenotype_mean"))
        dico_traits["Phenotype_mean"].append(pheno_mean if pheno_mean != 0.0 else "NaN")
        geno_mean = float(getattr(n, "Genotype_mean"))
        dico_traits["Genotype_mean"].append(geno_mean if geno_mean != 0.0 else "NaN")
    df_traits = pd.DataFrame(dico_traits)
    df_traits.to_csv(traits_path, sep=("\t" if fossils_path.endswith(".tsv") else ","), index=False)

    dico_fossils = defaultdict(list)
    for n in tree.traverse():
        n.dist = 0. if n.is_root() else float(getattr(n, "d"))

    assert tree.is_root()
    farthest, eccentricity = tree.get_farthest_leaf()
    for n in tree.traverse():
        dico_fossils["NodeName"].append(n.name)
        age_from_root = n.get_distance(tree)
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
    parser.add_argument("--tree", help="Output tree file", required=True)
    parser.add_argument("--traits", help="Output traits file", required=True)
    parser.add_argument("--fossils", help="Output fossils file", required=True)
    args = parser.parse_args()
    main(args.input, args.tree, args.traits, args.fossils)