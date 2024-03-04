import os
import argparse
from glob import glob

import numpy as np
from natsort import natsorted
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries import hist_plot, scatter_plot


def open_log(path):
    df = pd.read_csv(path, sep="\t")
    df["t_half"] = np.log(2) / df["alpha"]
    return df


def posterior(df, col="alpha", burnin=0.5):
    return np.nanmean(df[col][int(len(df) * burnin):])


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    model_prefs = {"moving_optimum": 0, "directional": 1, "stabilizing": 2, "neutral": 3}
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: natsorted([i for i in glob(f"{p}/*") if isdir(i)]) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    simple_OU_dict = defaultdict(lambda: defaultdict(list))
    relaxed_OU_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, folderpath in enumerate(replicates[m]):
            simple_df = open_log(f"{folderpath}/simple_OU_RJ.log.gz")
            simple_OU_dict["model"][m].append(m)
            for col in ["alpha", "is_BM", "is_OU", "sigma2", "theta", "t_half"]:
                simple_OU_dict[col][m].append(posterior(simple_df, col=col, burnin=0.5))

            relaxed_df = open_log(f"{folderpath}/relaxed_OU_RJ.log.gz")
            relaxed_OU_dict["model"][m].append(m)
            for col in ["alpha", "num_theta_changes", "sigma2", "theta_root", "t_half"]:
                relaxed_OU_dict[col][m].append(posterior(relaxed_df, col=col, burnin=0.5))

    rename = lambda x: output.replace(".tsv.gz", x)

    hist_plot(simple_OU_dict["is_OU"], "p[OU]", rename(f".simple_OU.is_OU.pdf"), xscale="uniform")
    hist_plot(simple_OU_dict["is_BM"], "p[BM]", rename(f".simple_OU.is_BM.pdf"), xscale="uniform")
    hist_plot(simple_OU_dict["t_half"], "t 1/2", rename(f".simple_OU.t_half.pdf"), xscale="log")
    hist_plot(simple_OU_dict["alpha"], "alpha", rename(f".simple_OU.alpha.pdf"), xscale="log")
    hist_plot(simple_OU_dict["sigma2"], "sigma2", rename(f".simple_OU.sigma2.pdf"), xscale="log")
    hist_plot(simple_OU_dict["theta"], "theta", rename(f".simple_OU.theta.pdf"), xscale="linear")

    hist_plot(relaxed_OU_dict["t_half"], "t 1/2", rename(f".relaxed_OU.t_half.pdf"), xscale="log")
    hist_plot(relaxed_OU_dict["alpha"], "alpha", rename(f".relaxed_OU.alpha.pdf"), xscale="log")
    hist_plot(relaxed_OU_dict["sigma2"], "sigma2", rename(f".relaxed_OU.sigma2.pdf"), xscale="log")
    hist_plot(relaxed_OU_dict["theta_root"], "theta_root", rename(f".relaxed_OU.theta_root.pdf"), xscale="linear")
    hist_plot(relaxed_OU_dict["num_theta_changes"], "num_theta_changes", rename(f".relaxed_OU.num_theta_changes.pdf"),
              xscale="linear")

    out_dict = defaultdict(list)
    for m in models:
        for col, values in simple_OU_dict.items():
            out_dict[f"simple_OU.{col}"].extend(values[m])
        for col, values in relaxed_OU_dict.items():
            out_dict[f"relaxed_OU.{col}"].extend(values[m])
        print(f"{m}: {np.mean(simple_OU_dict['is_OU'][m])}")
    df = pd.DataFrame(out_dict)
    df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
