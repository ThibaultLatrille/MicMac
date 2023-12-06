import os
import argparse
from glob import glob

import numpy as np
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

    replicates = {m: natsorted([i for i in glob(f"{p}/*") if isdir(i)]) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    var_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, folderpath in enumerate(replicates[m]):
            filepath = f"{folderpath}/simple_OU_RJ.log"
            df = pd.read_csv(filepath, sep="\t")
            var_dict["model"][m].append(m)
            for col in ["alpha", "is_BM", "is_OU", "sigma2", "theta"]:
                var_dict[col][m].append(np.mean(df[col][len(df) // 2:]))
            var_dict["t_half"][m].append(np.log(2) / np.mean(df["alpha"][len(df) // 2:]))

    rename = lambda x: output.replace(".tsv.gz", x)
    hist_plot(var_dict["is_OU"], "p[OU]", rename(f".is_OU.pdf"), xscale="uniform")
    hist_plot(var_dict["is_BM"], "p[BM]", rename(f".is_BM.pdf"), xscale="uniform")
    hist_plot(var_dict["t_half"], "t 1/2", rename(f".t_half.pdf"), xscale="log")
    hist_plot(var_dict["alpha"], "alpha", rename(f".alpha.pdf"), xscale="log")
    hist_plot(var_dict["sigma2"], "sigma2", rename(f".sigma2.pdf"), xscale="log")
    hist_plot(var_dict["theta"], "theta", rename(f".theta.pdf"), xscale="linear")

    out_dict = defaultdict(list)
    for m in models:
        for col, values in var_dict.items():
            out_dict[col].extend(values[m])
        print(f"{m}: {np.mean(var_dict['is_OU'][m])}")
    df = pd.DataFrame(out_dict)
    df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
