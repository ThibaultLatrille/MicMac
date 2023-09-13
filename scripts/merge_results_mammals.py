import os
import argparse
import pandas as pd


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(tsv_ML_list: str, tsv_Bayes_list: str, output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    list_df = []
    for method, tsv_list in [("ML", tsv_ML_list), ("Bayesian", tsv_Bayes_list)]:
        for path in tsv_list:
            name_split = os.path.basename(path).replace(".ratio.tsv", "").split("_")
            df = pd.read_csv(path, sep='\t')
            df["method"] = method
            for i, p in enumerate(["dataset", "sex", "logT", "h2"]):
                df[p] = name_split[i]
            list_df.append(df)
    df_out = pd.concat(list_df)
    # Sort by trait, dataset, sex, logT and then method
    df_out = df_out.sort_values(by=["dataset", "trait", "logT", "h2", "sex", "method"],
                                ascending=[True, True, False, False, False, False])
    # Fill missing values with NaN
    df_out = df_out.fillna("NaN")
    df_out.to_csv(output, sep="\t", index=False, float_format="%.3f")

    # Filter out the ML results
    df_out = df_out[df_out["method"] == "Bayesian"]
    dico_label = {"dataset": "dataset", "trait": "trait", "h2": "$h^2$", "sex": "sex",
                  "nbr_taxa_between": "$n$", "ratio": "$\\rho$", "ratio_CI95": "$\\rho_{CI95}$",
                  "pp_ratio_greater_1": "$\\mathbb{P}(\\rho > 1)$"}
    columns = [i for i in dico_label if i in df_out]
    df_out = df_out[columns]
    df_out["dataset"] = df_out["dataset"].apply(lambda x: "Mammals" if x == "NeutralPhyloP" else x)
    df_out["trait"] = df_out["trait"].apply(lambda x: x.replace("_g", "").replace("_", " "))
    df_out["ratio_CI95"] = df_out["ratio_CI95"].apply(lambda x: "-".join([f"{float(i):.3f}" for i in x.split("--")]))
    df_out = df_out.rename(columns={i: (dico_label[i] if i in dico_label else i) for i in columns})
    df_out.to_latex(replace_last(output, ".tsv", ".tex"), index=False, float_format="%.3f",
                    escape=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tsv_ML", nargs="+", help="Input tsv ML file", required=True)
    parser.add_argument("--tsv_Bayes", nargs="+", help="Input tsv BayesCode file", required=True)
    parser.add_argument("--output", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.tsv_ML, args.tsv_Bayes, args.output)
