#!python3
import os
from glob import glob
import numpy as np
import pandas as pd


def minipage(size, file, name=""):
    out = "\\begin{minipage}{" + str(size) + "\\linewidth} \n"
    out += "\\centering " + name + "\n"
    out += "\\includegraphics[width=\\linewidth, page=1]{" + file + "} \n"
    out += "\\end{minipage}\n"
    return out


def write(exp, o):
    o.write("\\begin{center}\n")
    for method in ["neutral", "moving_optimum", "stabilizing"]:
        distance_path = f"{exp}/results/inference_ML.distance.{method}.pdf"
        assert os.path.exists(distance_path), f"Missing {distance_path}"
        o.write(minipage(0.32, distance_path, method.replace("_", " ").title()))

    for param in ["is_OU", "t_half", "sigma2"]:
        hist_path = f"{exp}/results/inference_RevBayes.{param}.pdf"
        assert os.path.exists(hist_path), f"Missing {hist_path}"
        o.write(minipage(0.32, hist_path, param.replace("_", " ").title()))
    o.write("\\end{center}\n")


def alpha(exp):
    tsv_path = f"{exp}/results/inference_RevBayes.tsv.gz"
    if not os.path.exists(tsv_path):
        print(f"Missing {tsv_path}")
        return
    tsv = pd.read_csv(tsv_path, sep="\t")
    tsv["t_half"] = np.log(2) / tsv["alpha"]
    groups = tsv.groupby("model")
    for model, group in groups:
        print(group["is_BM"])
        for col in ["t_half", "is_BM"]:
            print(f"\t{model}\t{col}\t{group[col].mean()}")


def main():
    tex_include = "docs/include-figures.tex"
    tex_target = "docs/supp-mat.tex"
    experiments = glob("experiments/*_u")
    o = open(tex_include, 'w')
    o.write("\\section{Experiments} \n")
    for exp in experiments:
        subsection = os.path.basename(exp).replace("_", " ").title()
        o.write("\\subsection{" + subsection + "}\n")
        write(exp, o)
    o.close()
    output_dir = os.path.dirname(tex_target)
    tex_to_pdf = f"pdflatex -synctex=1 -interaction=nonstopmode -output-directory={output_dir} {tex_target}"
    os.system(tex_to_pdf)
    os.system(tex_to_pdf)

    for exp in experiments:
        print(exp)
        alpha(exp)


if __name__ == '__main__':
    main()
