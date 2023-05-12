import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict


def open_trace_file(trace_file: str, burn_in: int) -> dict[str: np.ndarray]:
    df = pd.read_csv(trace_file, sep="\t")
    var_phy_dict = {}
    for col in df.columns:
        if col.startswith("Var_"):
            assert col.endswith("_mean")
            var_phy_dict[replace_last(col[4:], "_mean", "")] = df[col].values[burn_in:] / 4.0
    n_taxa = int(df["n_taxa"][0])
    return var_phy_dict, n_taxa


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(inference_path: str, variance_pop_path: str, burn_in: int, output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    traces_var_phy_dict, n_taxa = open_trace_file(inference_path, burn_in=burn_in)
    var_pop_df = pd.read_csv(variance_pop_path, sep="\t")

    output_dict = defaultdict(list)
    for trait, var_phy_trace in traces_var_phy_dict.items():
        variance_phy = np.mean(var_phy_trace)

        assert f"{trait}_variance" in var_pop_df.columns
        notna = np.isfinite(var_pop_df[f"{trait}_variance"])
        if f"{trait}_heritability" not in var_pop_df.columns:
            print(f"Warning: column {f'{trait}_heritability'} not found in {variance_pop_path}.")
            print("Assuming heritability = 1.0.")
            var_pop_df[f"{trait}_heritability"] = 1.0
        # Computing the genetic variance (geno = h² * pheno)
        genetic_variance_array = var_pop_df[f"{trait}_heritability"][notna] * var_pop_df[f"{trait}_variance"][notna]
        variance_pop_array = genetic_variance_array / (var_pop_df["pS"][notna])
        print(f'Found {len(genetic_variance_array)} species with a variance for {trait}.')
        variance_pop = np.mean(variance_pop_array)
        assert variance_pop > 0

        ratio = variance_phy / variance_pop
        ratio_pv = np.mean(var_phy_trace > variance_pop)
        print(f"variance_phy = {variance_phy}")
        print(f"variance_pop = {variance_pop}.")
        print(f"ratio = {ratio}.")
        output_dict["trait"].append(trait)
        output_dict["nbr_taxa_pop"].append(n_taxa)
        output_dict["nbr_taxa_phy"].append(n_taxa)
        output_dict["variance_pop"].append(variance_pop)
        output_dict["variance_phy"].append(variance_phy)
        output_dict["ratio"].append(ratio)
        output_dict["ratio_pv"].append(ratio_pv)
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(output_path, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inference', required=True, type=str, dest="inference", help="Input trace file path")
    parser.add_argument('--variance_pop', required=True, type=str, dest="variance_pop",
                        help="Input variance pop file path")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output file path")
    parser.add_argument('--burn_in', required=False, type=int, dest="burn_in", default=100, help="Burn in")
    args = parser.parse_args()
    main(args.inference, args.variance_pop, args.burn_in, args.output)
