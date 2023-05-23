import os
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from ete3 import Tree
from neutrality_index import open_tree, prune_tree


def rename_tree(tree: Tree) -> Tree:
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            if node.name.endswith("_E"):
                node.name = "_".join(node.name.replace("_E", "").split("_")[:-1])
            else:
                node.name = node.name.replace(" ", "_")
    return tree


def name_internal_nodes(tree: Tree) -> Tree:
    node_i = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    return tree


def main(path_input_traits, path_input_pS, path_input_dS, path_output_tree, path_output_traits,
         path_output_variance_pop, log_transform, sex):
    for path in [path_input_traits, path_input_pS, path_input_dS]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [path_output_tree, path_output_traits, path_output_variance_pop]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    tree = open_tree(path_input_dS, format_ete3=1)
    tree = rename_tree(tree)
    set_taxa_names = set(tree.get_leaf_names())
    print(f"The tree has {len(set_taxa_names)} leaves")

    # Open the pS/hererozygosity file
    read_pS = open(path_input_pS, "r").readlines()
    if read_pS[0].strip() == "#TRAITS":
        columns = read_pS[1].strip().split(" ")
        n_taxa = int(columns[0])
        n_traits = int(columns[1])
        columns = ["species"] + columns[2:]
        assert len(columns) == n_traits + 1
        pS_df = pd.read_csv(path_input_pS, sep=" ", skiprows=2, names=columns)
        # Replace -1 by NaN
        pS_df = pS_df.replace(-1, np.nan)
        # Remove "_E" from the specie name
        pS_df["species"] = pS_df["species"].apply(lambda x: x.replace("_E", ""))
        assert len(pS_df) == n_taxa
    else:
        pS_df = pd.read_csv(path_input_pS, sep=",")

    # Remove all spaces in the specie name
    pS_df["species"] = pS_df["species"].apply(lambda x: x.replace(" ", "_"))
    print(f"The pS dataframe has {len(pS_df)} rows before filtering.")
    pS_df = pS_df[pS_df["species"].isin(set_taxa_names)]
    pS_col_list = ["heterozygosity", "pS", 'Hzoo']
    for pS_col in pS_col_list:
        if pS_col in pS_df.columns:
            print(f"Keeping only species with available '{pS_col}'.")
            pS_df = pS_df[np.isfinite(pS_df[pS_col]) & (pS_df[pS_col] > 0)]
            break
    print(f"The pS dataframe has {len(pS_df)} rows after filtering taxon name and available pS.")
    assert len(pS_df) >= 5, "Not enough species with pS. Exiting."

    df_traits = pd.read_csv(path_input_traits)
    assert "Genus_Species" in df_traits.columns
    print(f"The trait dataframe has {len(df_traits)} rows before filtering taxa.")
    df_traits = df_traits[df_traits["Genus_Species"].isin(set_taxa_names)]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering taxa also in the tree.")
    if sex in ["m", "f"]:
        # Keep only the male or the female
        df_traits = df_traits[df_traits["Sex"] == sex]
        print(f"The trait dataframe has {len(df_traits)} rows after keeping {'male' if sex == 'm' else 'female'}s.")
    else:
        print("Keeping all individuals (males and females).")
    # Keep only the adult
    df_traits = df_traits[df_traits["Age Class"] == "Adult"]
    print(f"The trait dataframe has {len(df_traits)} rows after keeping adults.")

    # Convert the brain volume in brain mass
    # density of brain tissue is 1.036 g/ml
    bm, bv = "Brain mass (g)", "Brain Volume (ml)"
    density_brain_tissue = 1.036
    df_traits[bm] = df_traits.apply(lambda r: r[bm] if np.isfinite(r[bm]) else r[bv] * density_brain_tissue, axis=1)
    ss_brain, ss_body = "Sample size (Brain)", "Sample size (Body)"
    df_traits[ss_body] = df_traits.apply(lambda r: r[ss_body] if np.isfinite(r[ss_body]) else r[ss_brain], axis=1)
    df_traits.to_csv(path_output_traits.replace("traits.tsv", "dataframe.tsv"), sep="\t", index=False)

    # Filter the tree and create dictionaries for variance and mean of traits
    set_taxa_names = set_taxa_names.intersection(set(df_traits["Genus_Species"]))
    print(f"The intersection of the tree and trait dataframe has {len(set_taxa_names)} taxa.")
    tree = prune_tree(tree, list(set_taxa_names))
    assert len(tree.get_leaves()) == len(set_taxa_names)
    dico_variance_pop, dico_traits = defaultdict(list), defaultdict(list)
    for taxa_name in set_taxa_names:
        leaf_pS_df = pS_df[pS_df["species"] == taxa_name]
        if len(leaf_pS_df) == 0:
            pS = np.nan
        else:
            assert len(leaf_pS_df) == 1
            pS = float(leaf_pS_df[pS_col])
        dico_variance_pop["TaxonName"].append(taxa_name)
        dico_variance_pop[f"pS"].append(pS)
        dico_traits["TaxonName"].append(taxa_name)

    for trait in ["Body mass (g)", "Brain mass (g)"]:
        trait_name = trait.replace(" ", "_").replace("(", "").replace(")", "")
        print(f"\nPhenotype considered is {trait}")
        trait_df = df_traits.copy()
        trait_df = trait_df[np.isfinite(trait_df[trait])]
        print(f"The trait dataframe has {len(trait_df)} rows after filtering not finite values.")
        if log_transform:
            # Log transform the trait
            trait_df[trait] = np.log(trait_df[trait])
            print(f"The trait is log transformed.")

        # Filter out the species with a unique row in the trait dataframe
        var_df = trait_df.copy()
        var_df = var_df[var_df["Sample size (Brain)"] == 1]
        print(f"The trait dataframe has {len(var_df)} rows after filtering for Nsize != 1.")
        var_grouped = {k: v for k, v in var_df.groupby("Genus_Species") if len(v) > 1}
        print(f"The trait dataframe has {len(var_grouped)} taxa after keeping taxon with more than 1 individuals.")
        for taxa_name in set_taxa_names:
            if taxa_name in var_grouped:
                phenotype_var = np.var(var_grouped[taxa_name][trait], ddof=1)
            else:
                phenotype_var = np.nan
            dico_variance_pop[f"{trait_name}_variance"].append(phenotype_var)
        assert len(set_taxa_names) == len(dico_variance_pop[f"{trait_name}_variance"])
        print(f"{len(var_grouped)} species with variance computed.")

        mean_df = trait_df.copy()
        mean_grouped = {k: v for k, v in mean_df.groupby("Genus_Species")}
        for taxa_name in set_taxa_names:
            if taxa_name in mean_grouped:
                filtered = mean_grouped[taxa_name]
                phenotype_mean = np.average(filtered[trait], weights=filtered["Sample size (Body)"])
            else:
                phenotype_mean = np.nan
            dico_traits[f"{trait_name}_mean"].append(phenotype_mean)
        print(f"{sum(np.isfinite(dico_traits[f'{trait_name}_mean']))} species with mean computed.")

    df_traits = pd.DataFrame(dico_traits)
    df_traits = df_traits[np.isfinite(df_traits.drop(["TaxonName"], axis=1)).any(axis=1)]
    df_traits.to_csv(path_output_traits, sep="\t", index=False, na_rep="NaN")
    set_taxa_names_traits = set(df_traits["TaxonName"])

    # Remove the species with only nan values across the traits
    df_variance_pop = pd.DataFrame(dico_variance_pop)
    df_variance_pop = df_variance_pop[df_variance_pop["TaxonName"].isin(set(pS_df["species"]))]
    df_variance_pop = df_variance_pop[df_variance_pop["TaxonName"].isin(set_taxa_names_traits)]
    df_variance_pop = df_variance_pop[np.isfinite(df_variance_pop.drop(["TaxonName", "pS"], axis=1)).any(axis=1)]
    # Write NaN for the species with no variance
    df_variance_pop.to_csv(path_output_variance_pop, sep="\t", index=False, na_rep="NaN")

    # Prune the tree and write it
    set_taxa_names = set(df_traits["TaxonName"])
    tree = prune_tree(tree, list(set_taxa_names))
    for taxa in df_variance_pop["TaxonName"]:
        assert taxa in set_taxa_names, f"{taxa} not in the tree"
    print(f"The final tree has {len(tree.get_leaves())} taxa.")
    tree_length = sum([node.dist for node in tree.traverse()])
    print(f"The tree length is {tree_length}.")
    tree.write(outfile=path_output_tree, format=3)


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
    main(args.input_traits, args.input_pS, args.input_dS, args.output_tree, args.output_traits,
         args.output_variance_pop, args.log_transform, args.sex)
