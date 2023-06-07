import os
import argparse
from glob import glob
from natsort import natsorted
import numpy as np
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries import hist_plot, scatter_plot


def open_covar_file(covar_file):
    """
    The file has the following format:
    entries are in the following order:
    Phenotype_mean
    Genotype_mean

    covariances

    62.266  61.58
    61.58 61.114

    correlation coefficients

    1  0.998
    0.998      1

    posterior probs

    -      1
    1      -
    """
    with open(covar_file) as f:
        f.readline()
        header = []
        for line in f:
            if line.strip() == "":
                break
            header.append(line.strip())
        matrix = [[0. for _ in range(len(header))] for _ in range(len(header))]
        assert f.readline().strip() == "covariances"
        assert f.readline().strip() == ""
        row = 0
        for line in f:
            if line.strip() == "":
                break
            matrix[row] = [float(x.strip()) for x in line.strip().split(" ") if x.strip() != ""]
            if matrix[row][row] < 0:
                matrix[row][row] = 0
                exit(1)
            assert len(matrix[row]) == len(header)
            row += 1
    return header, matrix


def open_trace_file(trace_file, burn_in):
    df = pd.read_csv(trace_file, sep="\t")
    return df["Var_Phenotype_mean"].values[burn_in:] / 4.0


def main(folder, tsv_path, burn_in, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    # Read the tsv file with multiindex columns
    df = pd.read_csv(tsv_path, sep="\t", header=None)
    df.columns = pd.MultiIndex.from_arrays([df.iloc[0], df.iloc[1]])
    df = df.iloc[2:]  # trim off the first 2 rows after we put them in the MultiIndex
    dict_tree = dict()
    for name, model in df.columns:
        if "replicate" == name:
            continue
        reps = df[('replicate', model)].values
        for rep, val in zip(reps, df[(name, model)].values):
            dict_tree[(name, model, rep)] = float(val)

    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: natsorted(glob(f"{p}/*.trace.gz")) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            rep = os.path.basename(filepath).replace(".trace.gz", "").split("_")[-1]
            var_phy_trace = open_trace_file(filepath, burn_in=burn_in)
            var_phy = np.mean(var_phy_trace)
            var_pop = dict_tree[("within_sampled", m, rep)]
            if var_pop <= 0:
                continue
            var_mut = dict_tree[("mutational", m, rep)]
            assert var_phy > 0
            assert var_mut > 0
            var_dict["phy"][m].append(var_phy)
            var_dict["pop"][m].append(var_pop)
            var_dict["mut"][m].append(var_mut)
            var_dict["phy_pop"][m].append(var_phy / var_pop)
            var_dict["phy_pop_pv"][m].append(np.mean(var_phy_trace > var_pop))
            var_dict["phy_mut"][m].append(var_phy / var_mut)

    df = pd.DataFrame(var_dict)
    df.to_csv(output, sep="\t", index=False)
    rename = lambda x: output.replace(".tsv.gz", x)
    hist_plot(var_dict["phy_pop"], "$\\rho$", rename(".pdf"))
    hist_plot(var_dict["phy_pop_pv"], "$P ( \\rho > 1 )$", rename(".pvalues.pdf"), xscale="linear")

    mut_str = r"Variance mutational $\left( \sigma^2_{M} \right)$"
    within_str = r"Variance within $\left( \sigma^2_{W} \right)$"
    between_str = r"Variance between $\left( \sigma^2_{B} \right)$"
    scatter_plot(var_dict["mut"], var_dict["phy"], mut_str, between_str, rename(".mutphy.pdf"), histy_log=True)
    scatter_plot(var_dict["mut"], var_dict["pop"], mut_str, within_str, rename(".mutpop.pdf"), histy_log=True)
    scatter_plot(var_dict["pop"], var_dict["phy"], within_str, between_str, rename(".popphy.pdf"), histy_log=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=True, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    parser.add_argument('-b', '--burn_in', required=False, type=int, dest="burn_in", default=100, help="Burn in")
    args = parser.parse_args()
    main(args.folder, args.tsv, args.burn_in, args.output)
