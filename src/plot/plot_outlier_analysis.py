#!/usr/bin/env python

import pandas as pd
import src.data_processing.outlier_analysis as outlier_analysis
from src.data_processing.outlier_analysis import all_cancer_run_dbscan


def main():
    scatter = pd.read_csv("/results/scatter_df.csv")

    tcga_meta = pd.read_csv("/ref/tag_metadata.csv")


    outlier_analysis.scatter = scatter
    outlier_analysis.tcga_meta = tcga_meta

    filtered_loss = scatter[scatter["direction"] == "negative"].copy()
    filtered_loss["patient_id"] = filtered_loss["Participants"]

    loss_prefix = "loss"  
    loss_cph_summary, loss_cox_fig, loss_sample_list = all_cancer_run_dbscan(
        filtered_loss,
        eps=0.7,
        min_samples=2,
        scaling_method="z-score",
        strat=True,
        prefix=loss_prefix,
    )

    filtered_gain = scatter[scatter["direction"] == "positive"].copy()
    filtered_gain["patient_id"] = filtered_gain["Participants"]

    gain_prefix = "gain"
    gain_cph_summary, gain_cox_fig, gain_sample_list = all_cancer_run_dbscan(
        filtered_gain,
        eps=0.7,
        min_samples=2,
        scaling_method="z-score",
        strat=True,
        prefix=gain_prefix,
    )


if __name__ == "__main__":
    main()
