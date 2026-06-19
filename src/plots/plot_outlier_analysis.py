#!/usr/bin/env python

import pandas as pd
from src.data_processing.outlier_analysis import all_cancer_run_dbscan


def main():
    scatter = pd.read_csv("/results/scatter_df.csv")
    tcga_meta = pd.read_csv("/ref/tcga_metadata.csv")

    loss_cph_summary, loss_cox_fig, loss_sample_list = all_cancer_run_dbscan(
        scatter, tcga_meta, direction="negative", eps=0.7, min_samples=2,
        scaling_method="z-score", strat=True, prefix="loss"
    )

    gain_cph_summary, gain_cox_fig, gain_sample_list = all_cancer_run_dbscan(
        scatter, tcga_meta, direction="positive", eps=0.7, min_samples=2,
        scaling_method="z-score", strat=True, prefix="gain"
    )


if __name__ == "__main__":
    main()
