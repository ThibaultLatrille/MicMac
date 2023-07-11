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
            elif "PD_" in node.name[:3]:
                node.name = "_".join(node.name.split("_")[2:])
                assert node.name.count("_") == 1
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


def main(path_input_traits, path_input_nuc_div, path_input_dS, path_output_tree, path_output_traits,
         path_output_var_within, log_transform, sex):
    for path in [path_input_traits, path_input_nuc_div, path_input_dS]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [path_output_tree, path_output_traits, path_output_var_within]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    t = open_tree(path_input_dS, format_ete3=1)
    tree = name_internal_nodes(rename_tree(t))
    set_taxa_names = set(tree.get_leaf_names())
    print(f"The tree has {len(set_taxa_names)} leaves")

    # Open the nucleotide diversity file
    read_nuc_div = open(path_input_nuc_div, "r").readlines()
    if read_nuc_div[0].strip() == "#TRAITS":
        columns = read_nuc_div[1].strip().split(" ")
        n_taxa = int(columns[0])
        n_traits = int(columns[1])
        columns = ["species"] + columns[2:]
        assert len(columns) == n_traits + 1
        nuc_div_df = pd.read_csv(path_input_nuc_div, sep=" ", skiprows=2, names=columns)
        # Replace -1 by NaN
        nuc_div_df = nuc_div_df.replace(-1, np.nan)
        # Remove "_E" from the specie name
        nuc_div_df["species"] = nuc_div_df["species"].apply(lambda x: x.replace("_E", ""))
        assert len(nuc_div_df) == n_taxa
    else:
        nuc_div_df = pd.read_csv(path_input_nuc_div, sep=",")
        if "SPECIES_BINOMIAL" in nuc_div_df.columns:
            nuc_div_df["species"] = nuc_div_df["SPECIES_BINOMIAL"]

    # Remove all spaces in the specie name
    nuc_div_df["species"] = nuc_div_df["species"].apply(lambda x: x.replace(" ", "_"))
    print(f"The dataframe of nucleotide diversity has {len(nuc_div_df)} rows before filtering.")
    nuc_div_df = nuc_div_df[nuc_div_df["species"].isin(set_taxa_names)]
    nuc_div_col = [p for p in ["MEDIAN_HETEROZYGOSITY", "heterozygosity", "pS", 'Hzoo'] if p in nuc_div_df.columns][0]
    print(f"Keeping only species with available '{nuc_div_col}'.")
    nuc_div_df = nuc_div_df[np.isfinite(nuc_div_df[nuc_div_col]) & (nuc_div_df[nuc_div_col] > 0)]
    print(f"The dataframe has {len(nuc_div_df)} rows after filtering taxon name and available nucleotide diversity.")
    assert len(nuc_div_df) >= 5, "Not enough species with nucleotide diversity. Exiting."

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
    dico_var_within, dico_traits = defaultdict(list), defaultdict(list)
    for taxa_name in set_taxa_names:
        leaf_nuc_div_df = nuc_div_df[nuc_div_df["species"] == taxa_name]
        if len(leaf_nuc_div_df) == 0:
            nuc_div = np.nan
        else:
            assert len(leaf_nuc_div_df) == 1
            nuc_div = float(leaf_nuc_div_df[nuc_div_col])
        dico_var_within["TaxonName"].append(taxa_name)
        dico_var_within[f"Nucleotide_diversity"].append(nuc_div)
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
        var_df = var_df[var_df[ss_brain] == 1]
        print(f"The trait dataframe has {len(var_df)} rows after filtering for Nsize != 1.")
        var_grouped = {k: v for k, v in var_df.groupby("Genus_Species") if len(v) > 1}
        print(f"The trait dataframe has {len(var_grouped)} taxa after keeping taxon with more than 1 individuals.")
        for taxa_name in set_taxa_names:
            if taxa_name in var_grouped:
                phenotype_var = np.var(var_grouped[taxa_name][trait], ddof=1)
                phenotype_mean = np.average(var_grouped[taxa_name][trait])
            else:
                phenotype_var = np.nan
                phenotype_mean = np.nan
            dico_var_within[f"{trait_name}_variance"].append(phenotype_var)
            # dico_var_within[f"{trait_name}_heritability_lower"].append(0.1)
            # dico_var_within[f"{trait_name}_heritability_upper"].append(0.4)
            dico_traits[f"{trait_name}_mean"].append(phenotype_mean)
        assert len(set_taxa_names) == len(dico_var_within[f"{trait_name}_variance"])
        print(f"{len(var_grouped)} species with variance computed.")

    df_traits = pd.DataFrame(dico_traits)
    df_traits = df_traits[np.isfinite(df_traits.drop(["TaxonName"], axis=1)).any(axis=1)]
    df_traits.to_csv(path_output_traits, sep="\t", index=False, na_rep="NaN")
    set_taxa_names_traits = set(df_traits["TaxonName"])

    # Remove the species with only nan values across the traits
    df_var_within = pd.DataFrame(dico_var_within)
    df_var_within = df_var_within[df_var_within["TaxonName"].isin(set(nuc_div_df["species"]))]
    df_var_within = df_var_within[df_var_within["TaxonName"].isin(set_taxa_names_traits)]
    df_var_within = df_var_within[
        np.isfinite(df_var_within.drop(["TaxonName", "Nucleotide_diversity"], axis=1)).any(axis=1)]
    # Write NaN for the species with no variance
    df_var_within.to_csv(path_output_var_within, sep="\t", index=False, na_rep="NaN")

    # Prune the tree and write it
    set_taxa_names = set(df_traits["TaxonName"])
    tree = prune_tree(tree, list(set_taxa_names))
    for taxa in df_var_within["TaxonName"]:
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
    parser.add_argument('--output_var_within', help="Output var_within file", required=True)
    parser.add_argument('--log_transform', help="log transformed valued", default="false", type=str, required=False)
    parser.add_argument('--sex', help="Sex (m or f)", default="m", type=str, required=False)
    args = parser.parse_args()
    args.log_transform = args.log_transform.lower() == "true"
    main(args.input_traits, args.input_pS, args.input_dS, args.output_tree, args.output_traits,
         args.output_var_within, args.log_transform, args.sex)
