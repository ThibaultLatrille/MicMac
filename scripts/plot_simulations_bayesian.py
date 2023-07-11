import os
import argparse
from glob import glob
from natsort import natsorted
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries import hist_plot, scatter_plot


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: natsorted(glob(f"{p}/*.ratio.tsv")) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            df = pd.read_csv(filepath, sep="\t")
            var_dict["model"][m].append(m)
            for trait in ["Phenotype", "Genotype"]:
                df_trait = df[df["trait"] == trait]
                assert len(df_trait) == 1
                for col in df_trait.columns:
                    var_dict[f"{trait}_{col}"][m].append(df_trait[col].values[0])

    df = pd.DataFrame(var_dict)
    df.to_csv(output, sep="\t", index=False)
    rename = lambda x: output.replace(".tsv.gz", x)
    for trait in ["Phenotype", "Genotype"]:
        hist_plot(var_dict[f"{trait}_ratio"], "Neutrality index ($\\rho$)", rename(f".{trait}.pdf"))
        hist_plot(var_dict[f"{trait}_pp_ratio_greater_1"], "$\\mathbb{P} ( \\rho > 1 )$",
                  rename(f".{trait}.pvalues.pdf"), xscale="linear")
        scatter_plot(var_dict[f"{trait}_var_within"], var_dict[f"{trait}_var_between"],
                     r"Variance within $\left( \sigma^2_{W} \right)$",
                     r"Variance between $\left( \sigma^2_{B} \right)$",
                     rename(f".{trait}.scatter.pdf"), histy_log=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
