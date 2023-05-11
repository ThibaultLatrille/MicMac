import os
import argparse
from gzip import open as gzopen
from collections import defaultdict
import pandas as pd
import numpy as np
from ete3 import Tree


def open_tree(tree_path: str, format_ete3: int = 1) -> Tree:
    if tree_path.endswith(".gz"):
        newick = gzopen(tree_path).read().decode()
        t = Tree(newick, format=format_ete3)
    else:
        t = Tree(tree_path, format=format_ete3)
    for node in t.traverse("levelorder"):
        if node.is_leaf():
            node.name = "_".join(node.name.replace("_E", "").split("_")[:-1])
    return t


def main(path_input_traits, path_input_pS, path_input_dS, path_output_tree, path_output_traits,
         path_output_variance_pop, log_transform, sex):
    for path in [path_input_traits, path_input_pS, path_input_dS]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [path_output_tree, path_output_traits, path_output_variance_pop]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    dS_tree = open_tree(path_input_dS, format_ete3=1)
    print(f"The tree has {len(dS_tree.get_leaves())} leaves")
    read_pS = open(path_input_pS, "r").readlines()
    if read_pS[0].strip() == "#TRAITS":
        columns = read_pS[1].strip().split(" ")
        n_taxa = int(columns[0])
        n_traits = int(columns[1])
        columns = ["specie"] + columns[2:]
        assert len(columns) == n_traits + 1
        pS_df = pd.read_csv(path_input_pS, sep=" ", skiprows=2, names=columns)
        # Replace -1 by NaN
        pS_df = pS_df.replace(-1, np.nan)
        # Remove "_E" from the specie name
        pS_df["specie"] = pS_df["specie"].apply(lambda x: x.replace("_E", ""))
        assert len(pS_df) == n_taxa
    else:
        pS_df = pd.read_csv(path_input_pS, sep=",")
    print(f"The pS dataframe has {len(pS_df)} rows before filtering.")
    pS_df = pS_df[pS_df["specie"].isin([n.name for n in dS_tree.get_leaves()])]
    pS_col = "pS" if "pS" in pS_df.columns else "pS_uniq"
    pS_df = pS_df[np.isfinite(pS_df[pS_col]) & (pS_df[pS_col] > 0)]
    print(f"The pS dataframe has {len(pS_df)} rows after filtering taxon name and available pS.")
    set_dS_tree = set([n.name for n in dS_tree.get_leaves()])
    set_pS_df = set(pS_df["specie"])
    inter_set = set_dS_tree.intersection(set_pS_df)
    print(f"The intersection between the tree and the pS dataframe has {len(inter_set)} elements.")

    df_traits = pd.read_csv(path_input_traits)
    assert "Genus_Species" in df_traits.columns
    print(f"The trait dataframe has {len(df_traits)} rows before filtering taxon name.")
    df_traits = df_traits[df_traits["Genus_Species"].isin(inter_set)]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering taxon name.")
    # Keep only the male or the female
    df_traits = df_traits[df_traits["Sex"] == sex]
    print(f"The trait dataframe has {len(df_traits)} rows after keeping {'male' if sex == 'm' else 'female'}.")
    # Keep only the adult
    df_traits = df_traits[df_traits["Age Class"] == "Adult"]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering out not adult.")
    df_traits.to_csv(path_output_traits.replace("traits.tsv", "dataframe.tsv"), sep="\t", index=False)
    inter_set = inter_set.intersection(set(df_traits["Genus_Species"]))
    print(f"The intersection of the three sets has {len(inter_set)} taxa.")

    dS_tree.prune(inter_set)
    assert len(dS_tree.get_leaves()) == len(inter_set)
    dico_variance_pop, dico_traits = defaultdict(list), defaultdict(list)
    for leaf in dS_tree.get_leaves():
        pS = float(pS_df[pS_df["specie"] == leaf.name][pS_col])
        dico_variance_pop["TaxonName"].append(leaf.name)
        dico_variance_pop[f"pS"].append(pS)
        dico_traits["TaxonName"].append(leaf.name)

    df_traits["Ratio of brain mass to body mass"] = df_traits["Brain mass (g)"] / df_traits["Body mass (g)"]
    for trait in ["Body mass (g)", "Brain Volume (ml)", "Brain mass (g)", "Ratio of brain mass to body mass"]:
        name = trait.replace(" ", "_").replace("(", "").replace(")", "")
        print(f"\nPhenotype considered is {trait}")
        trait_df = df_traits[np.isfinite(df_traits[trait])].copy()
        print(f"The trait dataframe has {len(trait_df)} rows after filtering not finite values.")
        if log_transform:
            # Log transform the trait
            trait_df[trait] = np.log(trait_df[trait])
            print(f"The trait is log transformed.")

        # Filter out the species with a unique row in the trait dataframe
        var_df = trait_df[trait_df["Sample size (Brain)"] == 1].copy()
        print(f"The trait dataframe has {len(var_df)} rows after filtering for Nsize != 1.")
        var_df = var_df.groupby("Genus_Species").filter(lambda x: len(x) > 1)
        print(f"The trait dataframe has {len(var_df)} individuals after filtering out species with no variance.")
        dico_count = {k: len(v) for k, v in var_df.groupby("Genus_Species")}
        print(f"The number of individuals per species is {dico_count}")
        for leaf in dS_tree.get_leaves():
            if leaf.name in dico_count:
                phenotype_var = np.var(var_df[var_df["Genus_Species"] == leaf.name][trait], ddof=1)
            else:
                phenotype_var = np.nan
            dico_variance_pop[f"{name}_variance"].append(phenotype_var)
        assert sum(np.isfinite(dico_variance_pop[f"{name}_variance"])) == len(dico_count)
        print(f"{len(dico_count)} species with variance computed.")

        for leaf in dS_tree.get_leaves():
            filtered = trait_df[trait_df["Genus_Species"] == leaf.name]
            filtered = filtered[np.isfinite(filtered[trait])]
            if len(filtered) == 0:
                phenotype_mean = np.nan
            else:
                phenotype_mean = np.average(filtered[trait], weights=filtered["Sample size (Brain)"])
            dico_traits[f"{name}_mean"].append(phenotype_mean)
        print(f"{sum(np.isfinite(dico_traits[f'{name}_mean']))} species with mean computed.")

    df_traits = pd.DataFrame(dico_traits)
    df_traits = df_traits[np.isfinite(df_traits.drop(["TaxonName"], axis=1)).any(axis=1)]
    df_traits.to_csv(path_output_traits, sep="\t", index=False, na_rep="NaN")

    # Remove the species with only nan values across the traits
    df_variance_pop = pd.DataFrame(dico_variance_pop)
    df_variance_pop = df_variance_pop[df_variance_pop["TaxonName"].isin(df_traits["TaxonName"])]
    df_variance_pop = df_variance_pop[np.isfinite(df_variance_pop.drop(["TaxonName", "pS"], axis=1)).any(axis=1)]
    # Write NaN for the species with no variance
    df_variance_pop.to_csv(path_output_variance_pop, sep="\t", index=False, na_rep="NaN")
    dS_tree.prune(set(df_traits["TaxonName"]))
    set_taxa = set(dS_tree.get_leaf_names())

    for taxa in df_variance_pop["TaxonName"]:
        assert taxa in set_taxa, f"{taxa} not in the tree"
    print(f"The final tree has {len(dS_tree.get_leaves())} taxa.")
    # Add polytomies if branch length are 0
    remove_nodes = set([n for n in dS_tree.traverse() if (n.dist == 0.0 and not n.is_root())])
    for n in remove_nodes:
        n.delete()
    assert len(set([n for n in dS_tree.traverse()]).intersection(remove_nodes)) == 0
    node_i = 0
    for n in dS_tree.traverse():
        if not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    # Write the topology to a newick file
    dS_tree.write(outfile=path_output_tree, format=3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_traits", help="Input trait file", required=True)
    parser.add_argument("--input_pS", help="Input pS file", required=True)
    parser.add_argument("--input_dS", help="Input dS tree file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    parser.add_argument('--output_variance_pop', help="Output variance_pop file", required=True)
    parser.add_argument('--log_transform', help="log transformed valued", default="false", type=str, required=False)
    parser.add_argument('--sex', help="Sex (m or f)", default="m", type=str, required=False)
    args = parser.parse_args()
    args.log_transform = args.log_transform.lower() == "true"
    assert args.sex in ["m", "f"]
    main(args.input_traits, args.input_pS, args.input_dS, args.output_tree, args.output_traits,
         args.output_variance_pop, args.log_transform, args.sex)
