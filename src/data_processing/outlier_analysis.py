def all_cancer_run_dbscan(scatter, # scatterplot dataframe
                          tcga_meta, # clinical metadata dataframe from recount
                          direction='positive',
                          eps=0.7,
                          min_samples=2,
                          scaling_method='z-score',
                          strat=True,
                          prefix=''):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.cluster import DBSCAN
    from sklearn.preprocessing import StandardScaler, MinMaxScaler
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test


    if scaling_method not in ['none', 'min-max', 'z-score']:
        raise ValueError("scaling_method must be one of 'none', 'min-max', or 'z-score'")

    if scaling_method == 'min-max':
        scaler = MinMaxScaler()
    elif scaling_method == 'z-score':
        scaler = StandardScaler()
    else:
        scaler = None  

    # Helper for tumor stage (used later in Cox model)
    def simplify_stage(stage):
        if pd.isnull(stage):
            return np.nan
        s = stage.lower()
        ret = 'Unknown'
        if 'stage i' in s:
            ret = 'Stage I'
        if 'stage ii' in s:
            ret = 'Stage II'
        if 'stage iii' in s:
            ret = 'Stage III'
        if 'stage iv' in s:
            ret = 'Stage IV'
        return ret

    input_df = scatter[scatter['direction'] == direction].copy()
    input_df['patient_id'] = input_df['Participants']

    all_cancer_results = []

    # Run DBSCAN within each cancer type
    for cancer_type in input_df['cancer_type'].unique():
        cancer_results = input_df[input_df['cancer_type'] == cancer_type].copy()

        required_columns = ['patient_id', 'alternative_se_sig_pd', 'lcc_edges']
        missing_cols = [col for col in required_columns if col not in cancer_results.columns]
        if missing_cols:
            raise ValueError(f"Missing columns in input_df: {missing_cols}")

        cancer_results = cancer_results[['patient_id', 'alternative_se_sig_pd', 'lcc_edges']].copy()
        cancer_results['cancer_type'] = cancer_type

        # Clean + impute
        cancer_results[['alternative_se_sig_pd', 'lcc_edges']] = (
            cancer_results[['alternative_se_sig_pd', 'lcc_edges']]
            .replace([np.inf, -np.inf], np.nan)
        )
        cancer_results[['alternative_se_sig_pd', 'lcc_edges']] = (
            cancer_results[['alternative_se_sig_pd', 'lcc_edges']]
            .fillna(cancer_results[['alternative_se_sig_pd', 'lcc_edges']].mean())
        )

        # Scaling
        if scaler is not None:
            scaled_features = scaler.fit_transform(
                cancer_results[['alternative_se_sig_pd', 'lcc_edges']]
            )
            cancer_results['scaled_se_sig_pd'] = scaled_features[:, 0]
            cancer_results['scaled_lcc_edges'] = scaled_features[:, 1]
            clustering_features = ['scaled_se_sig_pd', 'scaled_lcc_edges']
        else:
            clustering_features = ['alternative_se_sig_pd', 'lcc_edges']

        # DBSCAN
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        cancer_results['dbscan_outlier'] = dbscan.fit_predict(cancer_results[clustering_features])
        cancer_results['is_outlier'] = cancer_results['dbscan_outlier'] == -1

        # Merge with clinical meta 
        merged_df = cancer_results.merge(
            tcga_meta,
            left_on='patient_id',
            right_on='submitter_id',
            how='left'
        )

        merged_df.replace('--', np.nan, inplace=True)
        merged_df['days_to_death'] = pd.to_numeric(merged_df['days_to_death'], errors='coerce')
        merged_df['days_to_last_follow_up'] = pd.to_numeric(
            merged_df['days_to_last_follow_up'],
            errors='coerce'
        )
        merged_df['vital_status'] = (
            merged_df['vital_status']
            .fillna('Alive')
            .replace({'not reported': 'Alive'})
        )

        def get_time(row):
            if isinstance(row['vital_status'], str) and row['vital_status'].lower() == 'dead':
                return row['days_to_death']
            else:
                return row['days_to_last_follow_up']

        merged_df['event'] = merged_df['vital_status'].apply(
            lambda x: 1 if isinstance(x, str) and x.lower() == 'dead' else 0
        )
        merged_df['survival_time'] = merged_df.apply(get_time, axis=1)

        survival_df = merged_df.dropna(subset=['survival_time', 'event'])
        all_cancer_results.append(survival_df)


    combined_df = pd.concat(all_cancer_results, ignore_index=True)
    outliers = combined_df[combined_df['is_outlier']].copy()
    non_outliers = combined_df[~combined_df['is_outlier']].copy()

    # Kaplan–Meier curves
    kmf_outlier = KaplanMeierFitter()
    kmf_non_outlier = KaplanMeierFitter()

    km_fig, km_ax = plt.subplots(figsize=(6, 6))
    kmf_outlier.fit(
        outliers['survival_time'],
        event_observed=outliers['event'],
        label='Outliers'
    )
    kmf_non_outlier.fit(
        non_outliers['survival_time'],
        event_observed=non_outliers['event'],
        label='Non-Outliers'
    )

    kmf_outlier.plot_survival_function(ci_show=True, ax=km_ax)
    kmf_non_outlier.plot_survival_function(ci_show=True, ax=km_ax)

    km_ax.set_title("Kaplan–Meier: Outliers vs. Non-Outliers (All Cancers)")
    km_ax.set_xlabel("Time (days)")
    km_ax.set_ylabel("Survival probability")
    km_ax.legend()
    for spine in ['top', 'bottom', 'left', 'right']:
        km_ax.spines[spine].set_color('black')
    plt.tight_layout()

    # Optional save
    if prefix:
        km_fig.savefig(f"{prefix}_kaplan_meier.pdf", bbox_inches='tight')

    # Log-rank test
    combined_results_val = logrank_test(
        outliers['survival_time'],
        non_outliers['survival_time'],
        event_observed_A=outliers['event'],
        event_observed_B=non_outliers['event']
    )
    print("Log-rank test p-value (outliers vs non-outliers):",
          combined_results_val.p_value)

    # Cox proportional hazards model
    combined_df['simplified_stage'] = combined_df['tumor_stage'].apply(simplify_stage)

    cox_df = combined_df.dropna(
        subset=['survival_time', 'event', 'gender', 'simplified_stage',
                'age_at_diagnosis', 'is_outlier']
    ).copy()

    cox_df['age_at_diagnosis'] = pd.to_numeric(
        cox_df['age_at_diagnosis'],
        errors='coerce'
    )
    cox_df = cox_df.dropna(subset=['age_at_diagnosis'])

    # Remove "Unknown" stage from Cox
    cox_df = cox_df[cox_df['simplified_stage'] != 'Unknown']

    cox_df['simplified_stage'] = pd.Categorical(
        cox_df['simplified_stage'],
        categories=['Stage I', 'Stage II', 'Stage III', 'Stage IV'],
        ordered=True
    )

    # Clean column names for lifelines
    cox_df.columns = cox_df.columns.str.replace(' ', '_')

    cph = CoxPHFitter()
    formula = 'age_at_diagnosis + C(gender) + C(is_outlier) + C(simplified_stage)'

    if strat:
        cph.fit(
            cox_df,
            duration_col='survival_time',
            event_col='event',
            formula=formula,
            strata=['cancer_type']
        )
    else:
        cph.fit(
            cox_df,
            duration_col='survival_time',
            event_col='event',
            formula=formula
        )

    cph_summary = cph.summary

    cox_fig, cox_ax = plt.subplots(figsize=(6, 4))
    cph.plot(ax=cox_ax)
    cox_ax.set_title('Cox model coefficients')
    for spine in ['top', 'bottom', 'left', 'right']:
        cox_ax.spines[spine].set_color('black')
    plt.tight_layout()

    if prefix:
        cox_fig.savefig(f"{prefix}_cox.pdf", bbox_inches='tight')

    # Scatterplots with outliers annotated (marker-style)
    #   'cancer_type', 'direction', 'alternative_se_sig_pd', 'lcc_edges', 'Participants'
    sns.set_style("whitegrid")

    cancer_types = scatter['cancer_type'].unique()
    # Custom color dict from your original code (kept)
    color_dict = {
        "OV": "#1f4f9c", "MESO": "#379f7c", "CHOL": "#9e9ac8", "SKCM": "#006d2c",
        "BRCA": "#de2d26", "STAD": "#54278f", "LIHC": "#6a51a3", "THCA": "#bcbddc",
        "GBM": "#02818a", "LUAD": "#56b87a", "READ": "#ae017e", "UCEC": "#fdae6b",
        "BLCA": "#fdb462", "LUSC": "#379f7c", "KIRP": "#cb181d", "COAD": "#d94a9f",
        "LGG": "#42a7b1", "ESCA": "#f768a1", "PRAD": "#fee391", "KIRC": "#f03b20",
        "LAML": "#980043", "ACC": "#2c7f2b", "CESC": "#fc9272", "UCS": "#fdae6b",
        "TGCT": "#d95f0e", "PAAD": "#969696", "KICH": "#cb464d"
    }

    outlier_ids = set(outliers['patient_id'])

    scatter_fig, axes = plt.subplots(2, 1, figsize=(15, 12), sharex=False)

    for i, direction in enumerate(['positive', 'negative']):
        ax = axes[i]
        ax.set_title(f'Outliers within cancer in the {direction} direction')

        sns.scatterplot(
            data=scatter[
                (scatter['direction'] == direction) &
                (~scatter['Participants'].isin(outlier_ids))
            ],
            x='alternative_se_sig_pd',
            y='lcc_edges',
            hue='cancer_type',
            palette=color_dict,
            s=50,
            alpha=0.3,
            legend=False,
            ax=ax
        )

        # Outliers: highlighted and effectively "annotated" via marker
        sns.scatterplot(
            data=scatter[
                (scatter['direction'] == direction) &
                (scatter['Participants'].isin(outlier_ids))
            ],
            x='alternative_se_sig_pd',
            y='lcc_edges',
            hue='cancer_type',
            palette=color_dict,
            marker='X',
            s=80,
            edgecolor='black',
            legend=False,
            ax=ax
        )

        # Cancer-type mean points
        avg_points = (
            scatter[scatter['direction'] == direction]
            .groupby('cancer_type')[['alternative_se_sig_pd', 'lcc_edges']]
            .mean()
            .reset_index()
        )
        sns.scatterplot(
            data=avg_points,
            x='alternative_se_sig_pd',
            y='lcc_edges',
            hue='cancer_type',
            palette=color_dict,
            marker='o',
            s=70,
            edgecolor='black',
            legend=False,
            ax=ax
        )

        ax.set_xlabel('Number of significant alternative splicing events')
        if i == 0:
            ax.set_ylabel('Edges in largest connected component')
        else:
            ax.set_ylabel('')

        ax.grid(False)
        for spine in ['top', 'bottom', 'left', 'right']:
            ax.spines[spine].set_color('black')

    plt.tight_layout()

    if prefix:
        scatter_fig.savefig(f"{prefix}_outlier_scatter.pdf", bbox_inches='tight')

    sample_list = outliers['patient_id'].tolist()

    return cph_summary, cox_fig, sample_list
