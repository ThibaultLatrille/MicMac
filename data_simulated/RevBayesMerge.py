#!python3
import os
from glob import glob


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
        if not os.path.exists(distance_path):
            continue
        o.write(minipage(0.32, distance_path, method.replace("_", " ").title()))

    for param in ["is_OU", "t_half", "sigma2"]:
        hist_path = f"{exp}/results/inference_RevBayes.simple_OU.{param}.pdf"
        if not os.path.exists(hist_path):
            continue
        o.write(minipage(0.32, hist_path, param.replace("_", " ").title()))

    for param in ["num_theta_changes", "t_half", "sigma2"]:
        hist_path = f"{exp}/results/inference_RevBayes.relaxed_OU.{param}.pdf"
        if not os.path.exists(hist_path):
            continue
        o.write(minipage(0.32, hist_path, param.replace("_", " ").title()))

    o.write("\\end{center}\n")


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


if __name__ == '__main__':
    main()
