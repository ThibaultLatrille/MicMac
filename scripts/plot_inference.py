import os
import argparse
from glob import glob
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries import *


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
            assert matrix[row][row] >= 0
            assert len(matrix[row]) == len(header)
            row += 1
    return header, matrix


def main(folder, tsv_path, output):
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

    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p)}
    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: glob(p + "/*.cov") for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    header_var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, filepath in enumerate(replicates[m]):
            rep = os.path.basename(filepath).replace(".cov", "").split("_")[-1]
            header, matrix = open_covar_file(filepath)
            for i, h in enumerate(header):
                sigma_phy = matrix[i][i]
                sigma_pop = dict_tree[("within_sampled", m, rep)]
                assert sigma_phy >= 0
                assert sigma_pop >= 0
                header_var_dict[h][m].append(sigma_phy / sigma_pop)

    df = pd.DataFrame(header_var_dict)
    df.to_csv(output, sep="\t", index=False)
    for h, var_dict in header_var_dict.items():
        nb_genes = set([len(v) for v in var_dict.values()]).pop()
        hist_plot(var_dict, "$\\frac{\\sigma}{V_P / \\Theta}$", output.replace(".tsv.gz", f".{h}.pdf"), nb_genes,
                  'Var inter divided by Var intra')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tsv', required=True, type=str, dest="tsv", help="Input tsv file")
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.tsv, args.output)
