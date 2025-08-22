import pandas as pd
import pyarrow.parquet as pq
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind, ttest_rel
from scipy.stats import t
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go
from gprofiler import GProfiler
import plotly.express as px
import streamlit as st
import scipy.stats as stats
from pathlib import Path
import sys
import matplotlib.ticker as mtick
from matplotlib.patches import Ellipse
from scipy.stats import gaussian_kde


#BASIC AND DATA FUNCTIONS
def bool_to_str(value):
    if value is True:
        return "TRUE"
    elif value is False:
        return "FALSE"
    else:
        return ""


def number_to_str(value):
    if isinstance(value, (int, float)):
        return str(value)
    else:
        return ""


def make_columns_unique(df):
    cols = pd.Series(df.columns)
    for dup in cols[cols.duplicated()].unique():
        dup_idx = cols[cols == dup].index.tolist()
        for i, idx in enumerate(dup_idx[1:], start=1):
            cols[idx] = f"{cols[idx]}_{i}"
    df.columns = cols
    return df


def sanitize_dataframe(df):
    df = df.copy()
    df = make_columns_unique(df)
    for col in df.columns:
        if isinstance(df[col], pd.Series):
            if not pd.api.types.is_numeric_dtype(df[col]):
                df[col] = df[col].astype(str)
        else:
            df[col] = df[col].astype(str)
    return df


def resource_path(relative_path: str) -> str:
    if hasattr(sys, "_MEIPASS"):
        base_path = Path(sys._MEIPASS)
    else:
        base_path = Path(__file__).parent
    return str(base_path / relative_path)


def rename_cols(df):
    rename_map = {
        "Run": "File.Name",
        "PG.ProteinNames": "ProteinNames",
        "PG.ProteinGroups": "ProteinNames",
        "Protein.IDs": "ProteinNames",
        "PG.Genes": "GeneNames",
        "Protein.Names": "ProteinNames",
        "Protein IDs": "ProteinNames",
        "Genes": "GeneNames",
        "Gene names": "GeneNames",
        "PEP.StrippedSequence": "Stripped.Sequence",
        "Sequence": "Stripped.Sequence",
        "EG.ModifiedSequence": "Modified.Sequence",
        "R.FileName": "File.Name",
        "EG.iRTPredicted": "Predicted.RT",
        "EG.iRTEmpirical": "RT"
    }

    total_quantity_cols = [col for col in df.columns if re.match(r"EG\.TotalQuantity", col)]
    if total_quantity_cols:
        rename_map[total_quantity_cols[0]] = "Precursor.Quantity"

    df = df.rename(columns=rename_map)

    if {"GeneNames", "ProteinNames", "Raw file", "Intensity"}.issubset(df.columns):
        df['Intensity'] = pd.to_numeric(df['Intensity'], errors='coerce').fillna(0)
        agg_df = df.groupby(['GeneNames', 'ProteinNames', 'Raw file'], dropna=False)['Intensity'].sum().reset_index()
        df = agg_df.pivot(index=['GeneNames', 'ProteinNames'], columns='Raw file', values='Intensity').reset_index()
        df.columns.name = None

    if "Reverse" in df.columns:
        df = df[~df["Reverse"].str.contains(r"\+", na=False)]
    if "Only identified by site" in df.columns:
        df = df[~df["Only identified by site"].str.contains(r"\+", na=False)]

    return df


@st.cache_data
def read_data(file):
    if hasattr(file, 'name'):
        filename = file.name
    else:
        filename = file

    ext = os.path.splitext(filename)[1].lower().strip(".")

    if ext not in ["csv", "tsv", "txt", "xlsx", "parquet"]:
        raise ValueError("Invalid file type")

    if ext == "csv":
        df = pd.read_csv(file)
    elif ext == "tsv" or ext == "txt":
        df = pd.read_csv(file, sep="\t")
    elif ext == "xlsx":
        df = pd.read_excel(file)
    elif ext == "parquet":
        table = pq.read_table(file)
        df = table.to_pandas()
    else:
        raise ValueError("Unsupported file type")

    df = rename_cols(df)

    return df


def extract_id_or_number(x):
    matches = re.findall(r'\d+', x)
    return matches[-1] if matches else x


@st.cache_data
def log2_transform_data(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns].copy()
    data_filtered.replace(0, np.nan, inplace=True)
    log2_data = np.log2(data_filtered)

    remaining_columns = [col for col in data.columns if col not in annotated_columns]
    combined_data = pd.concat([log2_data, data[remaining_columns]], axis=1)

    return combined_data


@st.cache_data
def inverse_log2_transform_data(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    annotated_columns = meta['sample'].tolist()
    log2_data = data[annotated_columns].copy()
    original_data = 2 ** log2_data

    remaining_columns = [col for col in data.columns if col not in annotated_columns]
    combined_data = pd.concat([original_data, data[remaining_columns]], axis=1)

    return combined_data


@st.cache_data
def filter_data(data, meta, num, filterops="per group"):
    meta = meta.copy()
    meta['id'] = meta['sample'].apply(extract_id_or_number)

    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns].copy()

    conditions = meta['condition'].unique()

    mask = pd.DataFrame(False, index=data_filtered.index, columns=conditions)
    for cond in conditions:
        cols = meta.loc[meta['condition'] == cond, 'sample']
        mask[cond] = (data_filtered[cols].notna().sum(axis=1) >= num)

    if filterops == "per group":
        rows_to_keep = mask.all(axis=1)
    else:
        rows_to_keep = mask.any(axis=1)

    return data.loc[rows_to_keep]


@st.cache_data
def qqnorm_plot(data, meta):
    sample_columns = meta['sample'].tolist()
    data_values = data[sample_columns].replace(0, np.nan)

    values_vector = data_values.values.flatten()
    values_vector = values_vector[~np.isnan(values_vector)]
    values_vector = values_vector[values_vector > 0]

    if len(values_vector) > 100_000:
        np.random.seed(187)
        values_vector = np.random.choice(values_vector, 100_000, replace=False)

    fig, ax = plt.subplots()
    stats.probplot(values_vector, plot=ax)
    ax.set_title("QQ Plot")
    ax.set_xlabel("Theoretical Quantiles")
    ax.set_ylabel("Sample Quantiles")
    return fig


@st.cache_data
def first_digit_distribution(data, meta):
    sample_columns = meta["sample"].tolist()
    data_values = data[sample_columns]
    data_values = data_values.replace(0, np.nan)
    values_vector = data_values.values.flatten()
    values_vector = values_vector[~np.isnan(values_vector)]
    values_vector = values_vector[values_vector > 0]
    first_digits = [int(str(int(v))[0]) for v in values_vector if v > 0]
    first_digits = [d for d in first_digits if d in range(1, 10)]
    digit_freq = pd.Series(first_digits).value_counts(normalize=True).sort_index()
    benford = {d: np.log10(1 + 1/d) for d in range(1, 10)}
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(digit_freq.index, digit_freq.values, color="darkgreen", label="Observed")
    for idx, val in zip(digit_freq.index, digit_freq.values):
        ax.text(idx, val + 0.01, f"{val:.1%}", ha="center", va="bottom", fontsize=9)
    benford_x = list(benford.keys())
    benford_y = list(benford.values())
    ax.plot(benford_x, benford_y, color="red", linewidth=2, marker="o", label="Benford")
    ax.set_xlabel("First Digit")
    ax.set_ylabel("Relative Frequency")
    ax.set_xticks(range(1, 10))
    ax.set_yticks(np.linspace(0, max(max(digit_freq.values), max(benford_y)) + 0.1, 6))
    ax.set_yticklabels([f"{y:.0%}" for y in ax.get_yticks()])
    ax.legend()
    ax.set_title("First Digit Distribution vs Benford’s Law")
    plt.tight_layout()
    return fig


@st.cache_data
def data_pattern_structure(data, meta):
    sample_columns = meta["sample"].tolist()
    data_values = data[sample_columns]
    all_values = data_values.values.flatten()
    all_values = all_values[~pd.isna(all_values)]
    value_freq = pd.Series(all_values).value_counts()
    freq_of_freq = value_freq.value_counts().sort_index()
    freq_of_freq_percent = 100 * freq_of_freq / freq_of_freq.sum()
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(freq_of_freq_percent.index.astype(str), freq_of_freq_percent.values,
                  color="skyblue")
    for bar, val in zip(bars, freq_of_freq_percent.values):
        ax.text(bar.get_x() + bar.get_width()/2, val, f"{val:.3f}%", ha="center",
                va="bottom", fontsize=8, color="blue")
    ax.set_xlabel("Duplicate Number Occurrences")
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, max(freq_of_freq_percent.values) * 1.2)
    plt.tight_layout()
    return fig, freq_of_freq_percent


def impute_values(data: pd.DataFrame, meta: pd.DataFrame, q=0.01, adj_std=1, ret=0, sample_wise=False, seed=69):
    np.random.seed(seed)
    data = data.replace(0, np.nan).copy()
    annotated_columns = meta["sample"].tolist()
    data_filtered = data[annotated_columns]
    combined_data = data_filtered.values.flatten()
    missing_value_mask = data_filtered.isna()
    data_without_na = data_filtered.loc[~data_filtered.isna().any(axis=1)]
    data_with_na = data_filtered.loc[data_filtered.isna().any(axis=1)]
    combined_data_without_na = data_without_na.values.flatten()
    combined_data_with_na = data_with_na.values.flatten()

    if ret == 1:
        fig, ax = plt.subplots()
        ax.hist([combined_data_without_na, combined_data_with_na],
                bins=50, density=True, alpha=0.6, histtype='stepfilled',
                label=["Without Missing Values", "With Missing Values"])
        ax.legend()
        ax.set_title("Overall distribution of data with and without missing values")
        ax.set_xlabel("log2 Intensity")
        ax.set_ylabel("Density")
        return fig

    mean_val = np.nanmean(combined_data)
    sd_val = np.nanstd(combined_data)
    quantile_val = np.nanquantile(combined_data, q)

    if ret == 2:
        fig, ax = plt.subplots()
        ax.hist(combined_data[~np.isnan(combined_data)],
                bins=50, density=True, alpha=0.6, color="blue", histtype='stepfilled')
        x = np.linspace(np.nanmin(combined_data), np.nanmax(combined_data), 200)
        ax.plot(x, (1/(sd_val*np.sqrt(2*np.pi))) * np.exp(-0.5*((x-mean_val)/sd_val)**2),
                color="red", label="Normal fit")
        ax.set_title("Overall data distribution and norm fit")
        ax.set_xlabel("log2 Intensity")
        ax.set_ylabel("Density")
        ax.legend()
        return fig

    if sample_wise:
        for col in annotated_columns:
            col_vals = data_filtered[col]
            col_sd = np.nanstd(col_vals)
            num_na = col_vals.isna().sum()
            col_quantile_val = np.nanquantile(col_vals, q)
            if num_na > 0:
                data_filtered.loc[col_vals.isna(), col] = np.random.normal(col_quantile_val, col_sd, num_na)
    else:
        for col in annotated_columns:
            col_vals = data_filtered[col]
            num_na = col_vals.isna().sum()
            if num_na > 0:
                data_filtered.loc[col_vals.isna(), col] = np.random.normal(quantile_val, sd_val*adj_std, num_na)

    imputed_data = data_filtered[missing_value_mask]
    non_imputed_data = data_filtered[~missing_value_mask]
    imputed_combined = imputed_data.values.flatten()
    non_imputed_combined = non_imputed_data.values.flatten()

    if ret == 3:
        fig, ax = plt.subplots()
        ax.hist([non_imputed_combined, imputed_combined],
                bins=50, density=True, alpha=0.6, histtype='stepfilled',
                label=["Non-Imputed Values", "Imputed Values"])
        ax.legend()
        ax.set_title("Distribution of data after imputation")
        ax.set_xlabel("log2 Intensity")
        ax.set_ylabel("Density")
        return fig

    data_imputed = data.copy()
    data_imputed[annotated_columns] = data_filtered
    return data_imputed


#QC FUNCTIONS
@st.cache_data
def coverage_plot(data, meta, id=True, header=True, legend=True, plot_colors=None, width=20, height=10, dpi=300):
    data = data.replace(0, np.nan)
    conditions = meta['condition'].unique()
    meta['sample'] = meta['sample'].astype(str)
    meta['id'] = meta['sample'].apply(extract_id_or_number)

    if id:
        meta['new_sample'] = meta.groupby('condition').cumcount() + 1
        meta['new_sample'] = meta.apply(
            lambda row: f"{row['condition']}_{row['new_sample']}\n ({row['id']})", axis=1
        )
    else:
        meta['new_sample'] = meta.groupby('condition').cumcount() + 1
        meta['new_sample'] = meta.apply(
            lambda row: f"{row['condition']}_{row['new_sample']}", axis=1
        )

    rename_dict = dict(zip(meta['sample'], meta['new_sample']))
    data_filtered = data.rename(columns=rename_dict)
    annotated_columns = meta['new_sample'].tolist()
    data_filtered = data_filtered[annotated_columns]
    data_filtered = data_filtered.notna().astype(int)
    data_long = data_filtered.reset_index().melt(id_vars=data_filtered.index.name or None,
                                                var_name='Sample', value_name='Value')
    data_annotated = data_long.merge(meta[['new_sample', 'condition']],
                                    left_on='Sample', right_on='new_sample')
    summary = data_annotated.groupby(['Sample', 'condition'], as_index=False)['Value'].sum()
    summary['Sample'] = pd.Categorical(summary['Sample'], categories=annotated_columns, ordered=True)
    summary['condition'] = pd.Categorical(summary['condition'], categories=conditions, ordered=True)
    summary = summary.sort_values(['condition', 'Sample'])

    fig, ax = plt.subplots(figsize=(width/2.54, height/2.54), dpi=dpi)
    colors = plot_colors if plot_colors is not None else plt.cm.tab10.colors
    color_map = {cond: colors[i % len(colors)] for i, cond in enumerate(conditions)}

    for cond in conditions:
        cond_data = summary[summary['condition'] == cond]
        ax.bar(cond_data['Sample'], cond_data['Value'], label=cond, color=color_map[cond])

    ax.set_xlabel('Condition')

    if 'ProteinNames' in data.columns:
        ax.set_ylabel('Number of proteins')
        ax.axhline(y=len(data['ProteinNames']), color='red', linestyle='--')
    elif 'PTM_Collapse_key' in data.columns:
        ax.set_ylabel('Number of phosphosites')
        ax.axhline(y=len(data['PTM_Collapse_key']), color='red', linestyle='--')
    else:
        ax.set_ylabel('Number')

    if header:
        if 'ProteinNames' in data.columns:
            fig.suptitle('Proteins per sample', fontsize=14)
        elif 'PTM_Collapse_key' in data.columns:
            fig.suptitle('Phosphosites per sample', fontsize=14)
        else:
            fig.suptitle('')
    else:
        fig.suptitle('')

    ax.tick_params(axis='x', rotation=90, labelsize=6)

    if legend:
        ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(right=0.75, top=0.88)
    else:
        if ax.legend_:
            ax.legend_.remove()
        plt.subplots_adjust(right=0.95, top=0.88)
    return fig


@st.cache_data
def coverage_plot_summary(data, meta, header=True, plot_colors=None, width=20, height=10, dpi=300):
    data = data.replace(0, np.nan)
    meta = meta.copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["id"] = meta["sample"].str.extract(r'(\d+|[A-Za-z]+)', expand=False)

    meta["new_sample"] = [
        f"{cond}_{i+1}\n({sid})"
        for cond in meta["condition"].unique()
        for i, sid in enumerate(meta.loc[meta["condition"] == cond, "id"])
    ]

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data = data.rename(columns=rename_dict)
    annotated_cols = meta["new_sample"].tolist()
    data_filtered = data[annotated_cols]
    data_binary = data_filtered.notna().astype(int)
    melted = data_binary.melt(var_name="Sample", value_name="Value")
    data_annotated = melted.merge(meta, left_on="Sample", right_on="new_sample")

    sample_summary = (
        data_annotated.groupby("Sample")
        .agg(Value=("Value", "sum"), condition=("condition", "first"))
        .reset_index()
    )

    condition_summary = (
        sample_summary.groupby("condition")
        .agg(mean_value=("Value", "mean"), sd_value=("Value", "std"))
        .reset_index()
    )

    conditions = condition_summary["condition"].tolist()
    if plot_colors is None:
        plot_colors = plt.cm.tab10.colors
    if len(plot_colors) < len(conditions):
        plot_colors = (plot_colors * (len(conditions) // len(plot_colors) + 1))[:len(conditions)]
    color_map = dict(zip(conditions, plot_colors))

    fig, ax = plt.subplots(figsize=(width/2.54, height/2.54), dpi=dpi)
    bar_positions = np.arange(len(conditions))
    bar_means = condition_summary["mean_value"].values
    bar_sds = condition_summary["sd_value"].values

    ax.bar(
        bar_positions,
        bar_means,
        yerr=bar_sds,
        capsize=5,
        color=[color_map[c] for c in conditions],
        edgecolor="black"
    )

    rng = np.random.default_rng()
    for i, cond in enumerate(conditions):
        cond_values = sample_summary.loc[sample_summary["condition"] == cond, "Value"].values
        jitter = rng.uniform(-0.2, 0.2, size=len(cond_values))
        ax.scatter(
            np.full(len(cond_values), bar_positions[i]) + jitter,
            cond_values,
            color="black",
            alpha=0.7,
            s=20,
            zorder=3
        )

    plot_title, y_label = "", "Number"
    red_line_value = None
    if header:
        if "ProteinNames" in data.columns:
            plot_title = "Proteins per sample"
            y_label = "Number of proteins"
            red_line_value = len(data["ProteinNames"])
            ax.axhline(red_line_value, linestyle="--", color="red")
        elif "PTM_Collapse_key" in data.columns:
            plot_title = "Phosphosites per sample"
            y_label = "Number of phosphosites"
            red_line_value = len(data["PTM_Collapse_key"])
            ax.axhline(red_line_value, linestyle="--", color="red")

    ax.set_title(plot_title)
    ax.set_ylabel(y_label)
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(conditions, rotation=90, ha="center")

    ymax = max(
        max(bar_means + np.nan_to_num(bar_sds, nan=0)),
        red_line_value if red_line_value is not None else 0
    )
    ax.set_ylim(0, ymax * 1.1)

    plt.tight_layout()
    return fig


@st.cache_data
def coverage_plot_pep(data, meta, id=True, header=True, legend=True,
                      plot_colors=None, width=10, height=6, dpi=100):
    data = data.copy()
    meta = meta.copy()

    meta['sample'] = meta['sample'].astype(str)
    meta['id'] = meta['sample'].apply(extract_id_or_number)

    if id:
        meta['new_sample'] = meta.groupby(meta['condition'].astype(str)).cumcount() + 1
        meta['new_sample'] = (
            meta['condition'].astype(str) + "_" +
            meta['new_sample'].astype(str) + "\n(" +
            meta['id'].astype(str) + ")"
        )
    else:
        meta['new_sample'] = meta.groupby(meta['condition'].astype(str)).cumcount() + 1
        meta['new_sample'] = meta['condition'].astype(str) + "_" + meta['new_sample'].astype(str)

    if "File.Name" in data.columns:
        df_wide = data.pivot_table(
            index="File.Name",
            columns="Stripped.Sequence",
            values="Precursor.Quantity",
            aggfunc="max",
            fill_value=0
        ).T.reset_index()
        df_wide.rename(columns={"index": "Stripped.Sequence"}, inplace=True)
        df_wide.columns = [extract_id_or_number(c) if c != "Stripped.Sequence" else c for c in df_wide.columns]
        rename_dict = dict(zip(meta['id'], meta['new_sample']))
        df_wide = df_wide.rename(columns=rename_dict)
        annotated_cols = [c for c in meta['new_sample'] if c in df_wide.columns]
        data_filtered = df_wide[["Stripped.Sequence"] + annotated_cols]
    else:
        rename_dict = dict(zip(meta['sample'], meta['new_sample']))
        data_filtered = data.rename(columns=rename_dict)
        annotated_cols = [c for c in meta['new_sample'] if c in data_filtered.columns]
        data_filtered = data_filtered[annotated_cols]

    for col in annotated_cols:
        data_filtered[col] = np.where(data_filtered[col] != 0, 1, 0)

    plot_data = data_filtered.melt(id_vars="Stripped.Sequence", var_name="Sample", value_name="Value")
    summary = plot_data.groupby(['Sample']).agg({'Value': 'sum'}).reset_index()
    summary["Sample"] = pd.Categorical(summary["Sample"], categories=meta["new_sample"].astype(str), ordered=True)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    colors = plot_colors if plot_colors is not None else plt.cm.tab10.colors
    cond_order = meta['condition'].astype(str).unique()
    color_map = {cond: colors[i % len(colors)] for i, cond in enumerate(cond_order)}
    bar_colors = [color_map[meta.loc[meta['new_sample'] == str(s), 'condition'].values[0]] for s in summary['Sample']]

    ax.bar(summary['Sample'], summary['Value'], color=bar_colors)
    ax.axhline(y=data_filtered.shape[0], color='red', linestyle='--')
    ax.set_ylabel("Number of peptides")
    if header:
        ax.set_title("Peptides per sample")

    ax.set_xticks(range(len(meta["new_sample"])))
    ax.set_xticklabels(meta["new_sample"].astype(str), rotation=90)

    if legend:
        handles = [plt.Rectangle((0,0),1,1, color=color_map[cond]) for cond in cond_order]
        ax.legend(handles, cond_order, title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    return fig


@st.cache_data
def missing_value_plot(data, meta, bin=0, header=True, text=True, text_size=3.88, width=7.87, height=3.94, dpi=300):
    annotated_columns = meta['sample']
    data_filtered = data[annotated_columns].copy()
    data_filtered = data_filtered.replace(0, np.nan)
    na_count = data_filtered.isna().sum(axis=1)

    if bin > 0:
        na_count = na_count.apply(lambda x: f">{bin}" if x > bin else str(x))
        levels_vec = [str(i) for i in range(bin + 1)] + [f">{bin}"]
    else:
        na_count = na_count.astype(str)
        max_count = int(na_count[na_count != 'nan'].astype(float).max()) if not na_count.empty else 0
        levels_vec = [str(i) for i in range(max_count + 1)]

    miss_vals = na_count.value_counts().reset_index()
    miss_vals.columns = ['na_count', 'Freq']
    miss_vals = miss_vals[miss_vals['na_count'] != str(len(annotated_columns))]
    miss_vals['na_count'] = pd.Categorical(miss_vals['na_count'], categories=levels_vec, ordered=True)
    miss_vals = miss_vals.sort_values('na_count')

    if "ProteinNames" in data.columns:
        plot_title = "Missing Value Plot - Protein Level"
    elif "PTM_Collapse_key" in data.columns:
        plot_title = "Missing Value Plot - Phosphosite Level"
    else:
        plot_title = "Missing value plot"

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    ax.bar(miss_vals['na_count'], miss_vals['Freq'], color='blue')
    ax.set_xlabel("Number of Missing Values")
    ax.set_ylabel("Frequency")

    if text:
        for idx, row in miss_vals.iterrows():
            ax.text(row['na_count'], row['Freq'] + max(miss_vals['Freq']) * 0.01, str(row['Freq']),
                    ha='center', fontsize=text_size)

    if header:
        ax.set_title(plot_title)

    ax.tick_params(axis='x', rotation=0)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    return fig


@st.cache_data
def missing_value_plot_prec(data, meta, bin=0, header=True, text=True,
                            text_size=8, width=10, height=6, dpi=100):
    df_wide = data.pivot_table(
        index='File.Name',
        columns='Precursor.Id',
        values='Precursor.Quantity',
        aggfunc='first'
    )

    df_wide = df_wide.T
    df_wide.columns = df_wide.columns.map(lambda x: extract_id_or_number(x))

    meta['sample'] = meta['sample'].map(lambda x: extract_id_or_number(x))

    annotated_columns = meta['sample'].tolist()
    data_filtered = df_wide.loc[:, df_wide.columns.isin(annotated_columns)]

    na_count = data_filtered.isna().sum(axis=1)

    if bin > 0:
        na_count = na_count.apply(lambda x: f">{bin}" if x > bin else str(x))
    else:
        na_count = na_count.astype(str)

    miss_vals = na_count.value_counts().sort_index()

    if bin > 0:
        levels_vec = [str(i) for i in range(bin + 1)] + [f">{bin}"]
    else:
        levels_vec = sorted(miss_vals.index, key=lambda x: int(x))

    miss_vals = miss_vals.reindex(levels_vec, fill_value=0)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    ax.bar(miss_vals.index, miss_vals.values, color='blue')

    if header:
        ax.set_title("Missing Value Plot - Precursor Level")
    ax.set_xlabel("Number of Missing Values")
    ax.set_ylabel("Frequency")

    if text:
        for i, val in enumerate(miss_vals.values):
            if val > 0:
                ax.text(i, val + 0.5, str(val), ha='center', va='bottom', fontsize=text_size)

    plt.tight_layout()
    return fig


@st.cache_data
def missing_value_plot_pep(data, meta, bin=0, header=True, text=True,
                           text_size=8, width=10, height=6, dpi=100):
    if "File.Name" in data.columns:
        df_wide = data.pivot_table(
            index='File.Name',
            columns='Stripped.Sequence',
            values='Precursor.Quantity',
            aggfunc='max'
        )
        df_wide = df_wide.T
        df_wide.columns = df_wide.columns.map(lambda x: extract_id_or_number(x))
    else:
        df_wide = data.copy()

    meta['sample'] = meta['sample'].map(lambda x: extract_id_or_number(x))
    annotated_columns = meta['sample'].tolist()
    data_filtered = df_wide.loc[:, df_wide.columns.isin(annotated_columns)]

    na_count = data_filtered.isna().sum(axis=1)

    if bin > 0:
        na_count = na_count.apply(lambda x: f">{bin}" if x > bin else str(x))
    else:
        na_count = na_count.astype(str)

    miss_vals = na_count.value_counts().sort_index()

    if bin > 0:
        levels_vec = [str(i) for i in range(bin + 1)] + [f">{bin}"]
    else:
        levels_vec = sorted(miss_vals.index, key=lambda x: int(x))

    miss_vals = miss_vals.reindex(levels_vec, fill_value=0)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    ax.bar(miss_vals.index, miss_vals.values, color='blue')

    if header:
        ax.set_title("Missing Value Plot - Peptide Level")
    ax.set_xlabel("Number of Missing Values")
    ax.set_ylabel("Frequency")

    if text:
        for i, val in enumerate(miss_vals.values):
            if val > 0:
                ax.text(i, val + 0.5, str(val), ha='center', va='bottom', fontsize=text_size)

    plt.tight_layout()
    return fig


def histo_int(data, meta, plot_colors=None, header=True, legend=True,
              width=20, height=10, dpi=300):
    figsize_inches = (width / 2.54, height / 2.54)
    fig, ax = plt.subplots(figsize=figsize_inches, dpi=dpi)

    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns]

    mean_intensities = []
    for condition in meta['condition'].unique():
        columns = meta.loc[meta['condition'] == condition, 'sample'].tolist()
        mean_intensity = data_filtered[columns].mean(axis=1, skipna=True)
        mean_intensity = mean_intensity[np.isfinite(mean_intensity)]
        mean_intensities.append((condition, mean_intensity))

    n_conditions = len(mean_intensities)

    if plot_colors is None:
        colors = plt.cm.tab10.colors
        if n_conditions > len(colors):
            colors = plt.cm.get_cmap('tab20').colors
    else:
        if len(plot_colors) < n_conditions:
            plot_colors = (plot_colors * (n_conditions // len(plot_colors) + 1))[:n_conditions]
        colors = plot_colors

    all_values = np.concatenate([vals for _, vals in mean_intensities])
    xmin, xmax = np.min(all_values), np.max(all_values)
    x_grid = np.linspace(xmin, xmax, 1000)

    for idx, (condition, vals) in enumerate(mean_intensities):
        kde = gaussian_kde(vals)
        density = kde(x_grid)
        ax.plot(x_grid, density, label=condition, color=colors[idx])

    ax.set_xlabel('log2 Intensity')
    ax.set_ylabel('Density')

    if header:
        ax.set_title('Distribution of measured intensity (log2)')

    if legend:
        ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
        plt.tight_layout()

    return fig


@st.cache_data
def boxplot_int(data, meta, outliers=False, header=True, plot_colors=None,
                width_cm=20, height_cm=10, dpi=100):
    n_conditions = meta['condition'].nunique()

    if plot_colors is not None:
        if len(plot_colors) < n_conditions:
            plot_colors = (plot_colors * (n_conditions // len(plot_colors) + 1))[:n_conditions]
    else:
        base_colors = plt.cm.tab10.colors
        if n_conditions > len(base_colors):
            base_colors = plt.cm.get_cmap('tab20').colors
        plot_colors = base_colors[:n_conditions]

    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns]

    mean_intensities = []
    conditions_order = meta['condition'].unique()
    for condition in conditions_order:
        columns = meta.loc[meta['condition'] == condition, 'sample'].tolist()
        mean_intensity = data_filtered[columns].mean(axis=1, skipna=True)
        mean_intensity = mean_intensity[np.isfinite(mean_intensity)]
        mean_intensities.append(mean_intensity)

    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54), dpi=dpi)

    boxprops = dict(linewidth=1.5, color='black')
    whiskerprops = dict(linewidth=1.5, color='black')
    capprops = dict(linewidth=1.5, color='black')
    medianprops = dict(linewidth=2.5, color='firebrick')
    flierprops = dict(marker='o', markerfacecolor='red', markersize=5, linestyle='none') if outliers else dict(marker='')

    bplot = ax.boxplot(
        mean_intensities,
        patch_artist=True,
        showfliers=outliers,
        boxprops=boxprops,
        whiskerprops=whiskerprops,
        capprops=capprops,
        medianprops=medianprops,
        flierprops=flierprops,
        positions=np.arange(1, n_conditions + 1)
    )

    for patch, color in zip(bplot['boxes'], plot_colors):
        patch.set_facecolor(color)

    ax.set_xticks(np.arange(1, n_conditions + 1))
    ax.set_xticklabels(conditions_order)
    ax.set_xlabel("Condition")
    ax.set_ylabel("log2 Intensity")

    if header:
        ax.set_title("Measured protein intensity values (log2)")

    if ax.get_legend() is not None:
        ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    return fig


@st.cache_data
def boxplot_int_single(
    data, meta, outliers=False, id=True, header=True, legend=True, plot_colors=None,
    width_cm=20, height_cm=10, dpi=100
):
    meta = meta.copy()
    meta["id"] = meta["sample"].apply(extract_id_or_number)

    if id:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = meta.apply(
            lambda row: f"{row['condition']}_{row['new_sample']} ({row['id']})", axis=1
        )
    else:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = meta.apply(
            lambda row: f"{row['condition']}_{row['new_sample']}", axis=1
        )

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data = data.rename(columns=rename_dict)

    intensities = data.melt(var_name="sample", value_name="intensity")
    intensities = intensities.merge(
        meta[["new_sample", "condition"]],
        left_on="sample", right_on="new_sample"
    )
    intensities = intensities.replace([np.inf, -np.inf], np.nan).dropna()

    samples = meta["new_sample"].tolist()

    grouped_data = [
        intensities.loc[intensities["sample"] == s, "intensity"].values
        for s in samples
    ]

    unique_conditions = meta["condition"].unique()
    if plot_colors is None:
        cmap = plt.cm.get_cmap("tab10", len(unique_conditions))
        condition_colors = dict(zip(unique_conditions, [cmap(i) for i in range(len(unique_conditions))]))
    else:
        condition_colors = dict(zip(unique_conditions, plot_colors))

    sample_colors = [condition_colors[c] for c in meta["condition"]]

    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54), dpi=dpi)
    bp = ax.boxplot(grouped_data, patch_artist=True, showfliers=outliers)

    for patch, color in zip(bp["boxes"], sample_colors):
        patch.set_facecolor(color)

    ax.set_xticks(range(1, len(samples) + 1))
    ax.set_xticklabels(samples, rotation=90, ha="center")
    ax.set_ylabel("log2 Intensity")

    if header:
        ax.set_title("Measured protein intensity values (log2)")
    else:
        ax.set_title("")

    if legend:
        from matplotlib.patches import Patch
        legend_handles = [
            Patch(facecolor=condition_colors[c], edgecolor="black", label=c)
            for c in unique_conditions
        ]
        ax.legend(handles=legend_handles, title="Condition",
                  loc="center left", bbox_to_anchor=(1, 0.5))
        plt.subplots_adjust(right=0.8)
    else:
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        plt.subplots_adjust(right=0.95)

    plt.tight_layout()
    return fig


@st.cache_data
def cov_plot(data, meta, outliers=False, header=True, legend=True, plot_colors=None,
             width_cm=20, height_cm=10, dpi=100):
    conditions = meta['condition'].unique()
    valid_conditions = []
    cv_data = []

    for cond in conditions:
        cols = meta.loc[meta['condition'] == cond, 'sample'].tolist()
        if len(cols) >= 2:
            subset = data[cols]
            means = subset.mean(axis=1, skipna=True)
            sds = subset.std(axis=1, skipna=True)
            cv = (sds / means) * 100
            cv = cv[np.isfinite(cv)]
            if len(cv) > 0:
                valid_conditions.append(cond)
                cv_data.append(cv)

    n_conditions = len(valid_conditions)

    if plot_colors is not None:
        if len(plot_colors) < n_conditions:
            plot_colors = (plot_colors * (n_conditions // len(plot_colors) + 1))[:n_conditions]
    else:
        base_colors = plt.cm.tab10.colors
        if n_conditions > len(base_colors):
            base_colors = plt.cm.get_cmap('tab20').colors
        plot_colors = base_colors[:n_conditions]

    fig, ax = plt.subplots(figsize=(width_cm / 2.54, height_cm / 2.54), dpi=dpi)

    boxprops = dict(linewidth=1.5, color='black')
    whiskerprops = dict(linewidth=1.5, color='black')
    capprops = dict(linewidth=1.5, color='black')
    medianprops = dict(linewidth=2.5, color='firebrick')
    flierprops = dict(marker='o', markerfacecolor='red', markersize=5, linestyle='none') if outliers else dict(marker='')

    positions = np.arange(1, n_conditions + 1)

    bplot = ax.boxplot(
        cv_data,
        patch_artist=True,
        showfliers=outliers,
        boxprops=boxprops,
        whiskerprops=whiskerprops,
        capprops=capprops,
        medianprops=medianprops,
        flierprops=flierprops,
        positions=positions
    )

    for patch, color in zip(bplot['boxes'], plot_colors):
        patch.set_facecolor(color)

    ax.set_xticks(positions)
    ax.set_xticklabels(valid_conditions)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Coefficient of Variation (%)")

    if header:
        ax.set_title("Coefficient of Variation")

    if legend:
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=plot_colors[i], edgecolor='black', label=cond) for i, cond in enumerate(valid_conditions)]
        ax.legend(handles=legend_handles, title="Condition", loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    plt.tight_layout()
    return fig


@st.cache_data
def pca_plot(data, meta, header=True, legend=True, dot_size=3, width_cm=20, height_cm=10, dpi=100, plot_colors=None):
    meta['condition'] = pd.Categorical(meta['condition'], categories=meta['condition'].unique(), ordered=True)
    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns].dropna()

    transposed_expr = data_filtered.T
    zero_variance_columns = transposed_expr.var(axis=0) == 0
    transposed_expr = transposed_expr.loc[:, ~zero_variance_columns]

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(transposed_expr)
    explained_variance = pca.explained_variance_ratio_ * 100

    pca_scores = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_scores['sample'] = transposed_expr.index
    pca_scores = pca_scores.merge(meta, on='sample')
    pca_scores['condition'] = pd.Categorical(pca_scores['condition'], categories=meta['condition'].unique(), ordered=True)

    figsize = (width_cm / 2.54, height_cm / 2.54)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    colors = plot_colors if plot_colors is not None else plt.cm.tab10.colors
    conditions = pca_scores['condition'].unique()
    color_map = {cond: colors[i % len(colors)] for i, cond in enumerate(conditions)}

    for cond in conditions:
        subset = pca_scores[pca_scores['condition'] == cond]
        ax.scatter(subset['PC1'], subset['PC2'], label=cond, s=dot_size*10, alpha=0.7, color=color_map[cond])

    ax.set_xlabel(f"Principal Component 1 - {explained_variance[0]:.2f}% variance")
    ax.set_ylabel(f"Principal Component 2 - {explained_variance[1]:.2f}% variance")

    if header:
        ax.set_title("PCA Plot")

    if legend:
        ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        legend_obj = ax.get_legend()
        if legend_obj is not None:
            legend_obj.remove()

    plt.tight_layout()
    return fig


def pca_plot_interactive(data, meta, plot_colors=None, header=True, legend=True):
    annotated_columns = meta["sample"].tolist()
    data_filtered = data.loc[:, annotated_columns]
    data_filtered = data_filtered.dropna()

    transposed_expr = data_filtered.T
    variances = transposed_expr.var(axis=0)
    transposed_expr = transposed_expr.loc[:, variances > 0]

    try:
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(transposed_expr)
        explained_variance = pca.explained_variance_ratio_ * 100

        pca_scores = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
        pca_scores["sample"] = transposed_expr.index
        pca_scores = pca_scores.merge(meta, on="sample")

        x_label = f"Principal Component 1 - {explained_variance[0]:.2f}% variance"
        y_label = f"Principal Component 2 - {explained_variance[1]:.2f}% variance"

        fig = px.scatter(
            pca_scores,
            x="PC1",
            y="PC2",
            color="condition",
            color_discrete_sequence=plot_colors,
            hover_data={"sample": True, "condition": True},
            opacity=0.8
        )

        fig.update_traces(marker=dict(size=8))

        layout_kwargs = dict(
            xaxis=dict(title=x_label, zeroline=False),
            yaxis=dict(title=y_label, zeroline=False),
            showlegend=legend
        )
        if header:
            layout_kwargs["title"] = "PCA Plot"

        fig.update_layout(**layout_kwargs)
        return fig

    except Exception as e:
        print(f"An error occurred during PCA computation: {e}")
        return None


@st.cache_data
def abundance_plot(data, meta, workflow="Protein", plot_colors=None, width_cm=20, height_cm=10, dpi=300,
                   legend=True, header=True):
    data = data.replace(0, np.nan)
    unique_conditions = meta['condition'].unique()

    if plot_colors is not None:
        if len(plot_colors) < len(unique_conditions):
            reps = int(np.ceil(len(unique_conditions) / len(plot_colors)))
            plot_colors = (plot_colors * reps)[:len(unique_conditions)]
    else:
        plot_colors = plt.cm.get_cmap('tab10').colors[:len(unique_conditions)]

    annotated_columns = meta['sample'].tolist()

    if workflow == "Protein":
        key_col = "ProteinNames"
    elif workflow == "Phosphosite":
        key_col = "PTM_Collapse_key"
    else:
        raise ValueError("workflow must be 'Protein' or 'Phosphosite'")

    data_filtered = data[[key_col] + annotated_columns].copy()

    mean_intensities = pd.DataFrame()
    mean_intensities[key_col] = data_filtered[key_col]

    for condition in unique_conditions:
        columns = meta.loc[meta['condition'] == condition, 'sample'].tolist()
        condition_data = data_filtered[[key_col] + columns].copy()
        means = condition_data[columns].mean(axis=1, skipna=True)
        means_log = np.log10(means + 1)
        mean_intensities[condition] = means_log

    if mean_intensities.shape[1] > 2:
        mean_intensities = mean_intensities[mean_intensities.iloc[:, 1:].isna().sum(axis=1) < mean_intensities.shape[1] - 2]
    else:
        mean_intensities = mean_intensities.dropna(subset=[mean_intensities.columns[1]])

    long_intensities = mean_intensities.melt(id_vars=[key_col], var_name="Condition", value_name="log10Intensity")
    long_intensities['Rank'] = long_intensities.groupby('Condition')['log10Intensity'].rank(ascending=False, method='first')
    long_intensities['Condition'] = pd.Categorical(long_intensities['Condition'], categories=unique_conditions, ordered=True)

    width_in = width_cm / 2.54
    height_in = height_cm / 2.54

    fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=dpi)

    for cond, color in zip(unique_conditions, plot_colors):
        subset = long_intensities[long_intensities['Condition'] == cond]
        ax.scatter(subset['Rank'], subset['log10Intensity'], label=cond, color=color, s=5)

    ax.set_xlabel(f"{workflow} Rank")
    ax.set_ylabel(f"log10 {workflow} Intensity")

    if header:
        ax.set_title("Abundance plot - all conditions")
    else:
        ax.set_title("")

    if legend:
        ax.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
        fig.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        if ax.get_legend() is not None:
            ax.legend_.remove()
        fig.tight_layout()

    return fig


@st.cache_data
def interactive_abundance_plot(data, meta, condition, workflow="Protein", search=None):
    data = data.replace(0, np.nan)
    annotated_columns = meta.loc[meta["condition"] == condition, "sample"].tolist()

    if workflow == "Protein":
        data_filtered = data[["ProteinNames"] + annotated_columns]
        feature_col = "ProteinNames"
    elif workflow == "Phosphosite":
        data_filtered = data[["PTM_Collapse_key"] + annotated_columns]
        feature_col = "PTM_Collapse_key"
    else:
        raise ValueError("workflow must be 'Protein' or 'Phosphosite'")

    condition_means = data_filtered[annotated_columns].mean(axis=1, skipna=True)
    condition_means = np.log10(condition_means + 1)

    mean_intensities = pd.DataFrame({
        "Feature": data_filtered[feature_col],
        "log10Intensity": condition_means
    }).dropna()

    mean_intensities["Feature"] = mean_intensities["Feature"].astype(str).str.split(";").str[0]

    mean_intensities = mean_intensities.sort_values("log10Intensity", ascending=False).reset_index(drop=True)
    mean_intensities["Rank"] = mean_intensities.index + 1

    if search is None:
        search = []
    mean_intensities["Color"] = np.where(mean_intensities["Feature"].isin(search), "red", "blue")

    mean_intensities["TextLabel"] = mean_intensities.apply(
        lambda row: row["Feature"] if row["Feature"] in search else "",
        axis=1
    )

    fig = px.scatter(
        mean_intensities,
        x="Rank",
        y="log10Intensity",
        color="Color",
        text="TextLabel",
        hover_data={"Feature": True, "Rank": True, "log10Intensity": True, "Color": False, "TextLabel": False},
        color_discrete_map={"red": "red", "blue": "blue"},
        title=f"Abundance plot - {condition}"
    )

    fig.update_traces(marker=dict(size=6), textposition="top center")

    fig.update_layout(
        xaxis_title=f"{workflow} Rank",
        yaxis_title=f"log10 {workflow} Intensity",
        showlegend=False
    )

    return fig


@st.cache_data
def corr_plot(data, meta, method="Matrix", id=True, full_range=False, width=10, height=8, dpi=100):
    meta = meta.copy()
    meta['sample'] = meta['sample'].astype(str)
    meta['id'] = meta['sample'].apply(extract_id_or_number)

    if id:
        meta['new_sample'] = meta.groupby('condition').cumcount() + 1
        meta['new_sample'] = meta.apply(lambda row: f"{row['condition']}_{row['new_sample']} ({row['id']})", axis=1)
    else:
        meta['new_sample'] = meta.groupby('condition').cumcount() + 1
        meta['new_sample'] = meta.apply(lambda row: f"{row['condition']}_{row['new_sample']}", axis=1)

    rename_map = dict(zip(meta['sample'], meta['new_sample']))
    data = data.rename(columns=rename_map)

    annotated_columns = meta['new_sample'].tolist()
    data_filtered = data[annotated_columns]

    correlation_matrix = data_filtered.corr(method='pearson')
    distance_matrix = 1 - correlation_matrix
    linkage_matrix = linkage(squareform(distance_matrix), method='complete')
    ordered_indices = leaves_list(linkage_matrix)
    ordered_corr = correlation_matrix.iloc[ordered_indices, ordered_indices]

    vmin, vmax = (-1, 1) if full_range else (ordered_corr.min().min(), ordered_corr.max().max())

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    if method == "Matrix":
        cax = ax.matshow(ordered_corr, cmap='coolwarm', vmin=vmin, vmax=vmax)
        plt.xticks(range(len(ordered_corr.columns)), ordered_corr.columns, rotation=90, ha='center', va='bottom')
        plt.yticks(range(len(ordered_corr.index)), ordered_corr.index)
        fig.colorbar(cax)

    elif method == "Ellipse":
        n = len(ordered_corr)
        ax.set_xlim(0, n)
        ax.set_ylim(0, n)
        ax.set_xticks(np.arange(n) + 0.5)
        ax.set_yticks(np.arange(n) + 0.5)
        ax.set_xticklabels(ordered_corr.columns, rotation=90, ha='center', va='top')
        ax.set_yticklabels(ordered_corr.index)
        ax.invert_yaxis()

        for i in range(n):
            for j in range(n):
                r = ordered_corr.iloc[i, j]
                width_ellipse = 0.9
                height_ellipse = 0.9 * (1 - abs(r))
                angle = 45 if r > 0 else -45
                ellipse = Ellipse(
                    (j + 0.5, i + 0.5),
                    width=width_ellipse,
                    height=height_ellipse,
                    angle=angle,
                    facecolor="red" if r > 0 else "blue",
                    alpha=0.6
                )
                ax.add_patch(ellipse)

        sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=-1, vmax=1))
        fig.colorbar(sm, ax=ax)

    else:
        raise ValueError("method must be 'matrix' or 'ellipse'")

    plt.tight_layout()
    return fig


#STATISTIC FUNCTIONS
@st.cache_data
def volcano_plot(volcano_df, condition1, condition2, in_pval=0.05, in_log2fc=1, uncorrected=False):
    color_mapping = {'Downregulated': 'blue', 'Not significant': 'gray', 'Upregulated': 'red'}
    y_hline = -np.log10(in_pval)
    max_y = max(y_hline, volcano_df['neg_log10_adj_pval'].max())

    shapes = [
        dict(type="line", x0=in_log2fc, x1=in_log2fc, y0=0, y1=max_y, line=dict(color="black", dash="dash")),
        dict(type="line", x0=-in_log2fc, x1=-in_log2fc, y0=0, y1=max_y, line=dict(color="black", dash="dash")),
        dict(type="line", x0=volcano_df['log2FC'].min(), x1=volcano_df['log2FC'].max(),
             y0=y_hline, y1=y_hline, line=dict(color="black", dash="dash"))
    ]

    label_col = 'ProteinNames' if 'ProteinNames' in volcano_df.columns else 'PTM_Collapse_key'

    fig = go.Figure()
    for sig, color in color_mapping.items():
        subset = volcano_df[volcano_df['significance'] == sig]
        fig.add_trace(go.Scatter(
            x=subset['log2FC'],
            y=subset['neg_log10_adj_pval'],
            mode='markers',
            marker=dict(color=color, size=5),
            text=subset[label_col],
            name=sig
        ))

    fig.update_layout(
        title=f"{condition1} vs {condition2}",
        xaxis_title=f"log2 fold change ({condition2} - {condition1})",
        yaxis_title='-log10 adj. p-value' if not uncorrected else '-log10 p-value',
        shapes=shapes,
        hovermode='closest'
    )
    return fig


@st.cache_data
def volcano_data_f(data, meta, condition1, condition2, in_pval=0.05, in_log2fc=1,
                   workflow="Protein", paired="Unpaired", uncorrected=False):
    label_col = "ProteinNames" if workflow == "Protein" else "PTM_Collapse_key"
    annotated_columns1 = meta.loc[meta['condition'] == condition1, 'sample'].tolist()
    annotated_columns2 = meta.loc[meta['condition'] == condition2, 'sample'].tolist()
    cols_to_use = [label_col] + annotated_columns1 + annotated_columns2

    df = data[cols_to_use].replace(0, np.nan).copy()
    arr_x = df[annotated_columns1].to_numpy(dtype=float)
    arr_y = df[annotated_columns2].to_numpy(dtype=float)

    n1_real = np.sum(~np.isnan(arr_x), axis=1)
    n2_real = np.sum(~np.isnan(arr_y), axis=1)
    valid_mask = (n1_real >= 2) & (n2_real >= 2)

    df = df.loc[valid_mask]
    arr_x = arr_x[valid_mask]
    arr_y = arr_y[valid_mask]
    n1 = n1_real[valid_mask]
    n2 = n2_real[valid_mask]

    mean1 = np.nanmean(arr_x, axis=1)
    mean2 = np.nanmean(arr_y, axis=1)
    log2fc = mean2 - mean1

    if paired == "Unpaired":
        var1 = np.nanvar(arr_x, axis=1, ddof=1)
        var2 = np.nanvar(arr_y, axis=1, ddof=1)
        pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
        se = np.sqrt(pooled_var * (1 / n1 + 1 / n2))
        t_stat = log2fc / se
        dfree = n1 + n2 - 2
        pvals = 2 * t.sf(np.abs(t_stat), dfree)
    else:
        diffs = arr_y - arr_x
        min_len = np.minimum(n1, n2)
        mean_diff = np.nanmean(diffs, axis=1)
        se_diff = np.nanstd(diffs, axis=1, ddof=1) / np.sqrt(min_len)
        t_stat = mean_diff / se_diff
        pvals = 2 * t.sf(np.abs(t_stat), min_len - 1)

    adj_pvals = pvals if uncorrected else multipletests(pvals, method='fdr_bh')[1]

    significance = np.full(len(log2fc), "Not significant", dtype=object)
    significance[(adj_pvals < in_pval) & (log2fc > in_log2fc)] = "Upregulated"
    significance[(adj_pvals < in_pval) & (log2fc < -in_log2fc)] = "Downregulated"

    return pd.DataFrame({
        label_col: df[label_col].values,
        "log2FC": log2fc,
        "pval": pvals,
        "adj_pval": adj_pvals,
        "significance": significance,
        "neg_log10_adj_pval": -np.log10(adj_pvals)
    })


gp = GProfiler(return_dataframe=True)

@st.cache_data
def different_genes(data, meta, condition1, condition2, in_pval=0.05, in_log2fc=1,
                    workflow="Gene", paired="Unpaired", uncorrected=False):
    annotated_columns1 = meta.loc[meta['condition'] == condition1, 'sample'].tolist()
    annotated_columns2 = meta.loc[meta['condition'] == condition2, 'sample'].tolist()

    key_col = 'GeneNames'

    cols_to_use = [key_col] + annotated_columns1 + annotated_columns2
    data_filtered = data[cols_to_use].copy()
    threshold = len(annotated_columns1) + len(annotated_columns2) - 2
    data_filtered = data_filtered[data_filtered.iloc[:, 1:].isna().sum(axis=1) < threshold]

    log2fc_list = []
    pvals_list = []

    for idx, row in data_filtered.iterrows():
        values1 = row[annotated_columns1].dropna().astype(float)
        values2 = row[annotated_columns2].dropna().astype(float)

        if paired == "Unpaired":
            if len(values1) < 2 or len(values2) < 2:
                log2fc_list.append(np.nan)
                pvals_list.append(np.nan)
                continue
            log2fc_list.append(values2.mean() - values1.mean())
            pvals_list.append(ttest_ind(values2, values1, equal_var=True).pvalue)

        elif paired == "Paired":
            paired_values = pd.concat([values1, values2], axis=1).dropna()
            if paired_values.shape[0] < 1:
                log2fc_list.append(np.nan)
                pvals_list.append(np.nan)
                continue
            log2fc_list.append(paired_values.iloc[:, 1].mean() - paired_values.iloc[:, 0].mean())
            if paired_values.shape[0] < 2:
                pvals_list.append(np.nan)
            else:
                pvals_list.append(ttest_rel(paired_values.iloc[:, 1], paired_values.iloc[:, 0]).pvalue)

    mask = ~np.isnan(pvals_list)
    log2fc_array = np.array(log2fc_list)[mask]
    pvals_array = np.array(pvals_list)[mask]
    data_filtered = data_filtered.iloc[mask, :]

    if uncorrected:
        adj_pvals = pvals_array
    else:
        adj_pvals = multipletests(pvals_array, method='fdr_bh')[1]

    result_data = pd.DataFrame({
        'Gene': data_filtered[key_col],
        'log2FC': log2fc_array,
        'pval': pvals_array,
        'adj_pval': adj_pvals
    })

    up_genes = result_data.loc[
        (result_data['adj_pval'] < in_pval) & (result_data['log2FC'] > in_log2fc), 'Gene'
    ].apply(lambda x: x.split(";")[0]).tolist()

    down_genes = result_data.loc[
        (result_data['adj_pval'] < in_pval) & (result_data['log2FC'] < -in_log2fc), 'Gene'
    ].apply(lambda x: x.split(";")[0]).tolist()

    return {'Upregulated': up_genes, 'Downregulated': down_genes}


@st.cache_data
def enrichment_analysis(gene_list, top_n=10, min_num=20, max_num=300):
    gene_list = [g for g in gene_list if g]
    if not gene_list:
        return None

    res = gp.profile(
        organism='hsapiens',
        query=gene_list,
        sources=['GO:BP', 'GO:CC', 'GO:MF'],
        user_threshold=0.05,
        significance_threshold_method='g_SCS'
    )

    if res.empty:
        return None

    res = res[(res['term_size'] >= min_num) & (res['term_size'] <= max_num)].copy()
    res = res.sort_values('p_value').head(top_n)

    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(
        x=res['intersection_size'] / res['term_size'] * 100,
        y=res['name'],
        s=res['intersection_size']*10,
        c=-np.log10(res['p_value']),
        cmap='RdBu_r'
    )

    ax.set_xlabel('Hits (%)')
    ax.set_ylabel('GO Term')
    ax.set_title('Gene Set Enrichment Analysis')

    cbar = fig.colorbar(scatter, ax=ax, pad=0.15)
    cbar.set_label('-log10(p-value)')

    plt.tight_layout()
    return fig


@st.cache_data
def volcano_plot_sim(data, meta, condition1, condition2, in_pval=0.05, in_log2fc=1,
                     workflow="Protein", mod_var=1, mod_n=0):
    annotated_columns1 = meta.loc[meta['condition'] == condition1, 'sample'].tolist()
    annotated_columns2 = meta.loc[meta['condition'] == condition2, 'sample'].tolist()
    label_col = "ProteinNames" if workflow == "Protein" else "PTM_Collapse_key"
    cols_to_use = [label_col] + annotated_columns1 + annotated_columns2
    data_filtered = data[cols_to_use].copy()
    threshold = len(annotated_columns1) + len(annotated_columns2) - 2
    mask = data_filtered.iloc[:, 1:].isna().sum(axis=1) < threshold
    data_filtered = data_filtered.loc[mask]
    x_data = data_filtered[annotated_columns1].astype(float)
    y_data = data_filtered[annotated_columns2].astype(float)
    n1_real = x_data.notna().sum(axis=1)
    n2_real = y_data.notna().sum(axis=1)
    valid_mask = (n1_real >= 2) & (n2_real >= 2)
    x_data = x_data[valid_mask]
    y_data = y_data[valid_mask]
    data_filtered = data_filtered.loc[valid_mask]
    mean1 = x_data.mean(axis=1)
    mean2 = y_data.mean(axis=1)
    var1 = x_data.var(axis=1, ddof=1)
    var2 = y_data.var(axis=1, ddof=1)
    n1 = n1_real[valid_mask] if mod_n in (0, 1) else mod_n
    n2 = n2_real[valid_mask] if mod_n in (0, 1) else mod_n
    pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    pooled_var *= mod_var
    se = np.sqrt(pooled_var * (1 / n1 + 1 / n2))
    t_stat = (mean2 - mean1) / se
    df = n1 + n2 - 2
    pvals = 2 * t.sf(np.abs(t_stat), df)
    adj_pvals = multipletests(pvals, method='fdr_bh')[1]
    log2fc = mean2 - mean1
    significance = np.where((adj_pvals < in_pval) & (log2fc > in_log2fc), "Upregulated",
                            np.where((adj_pvals < in_pval) & (log2fc < -in_log2fc), "Downregulated",
                                     "Not significant"))
    volcano_data = pd.DataFrame({
        'Feature': data_filtered[label_col],
        'log2FC': log2fc,
        'pval': pvals,
        'adj_pval': adj_pvals,
        'significance': significance
    })
    volcano_data['neg_log10_pval'] = -np.log10(volcano_data['adj_pval'])
    color_mapping = {"Downregulated": "blue", "Not significant": "gray", "Upregulated": "red"}
    fig = px.scatter(
        volcano_data,
        x='log2FC',
        y='neg_log10_pval',
        color='significance',
        hover_data=['Feature'],
        color_discrete_map=color_mapping,
        title=f"{condition1} vs. {condition2}",
    )
    vline_height = max(volcano_data['neg_log10_pval'].max(), -np.log10(in_pval))
    fig.add_shape(type="line", x0=in_log2fc, x1=in_log2fc, y0=0, y1=vline_height,
                  line=dict(dash='dash', color='black'))
    fig.add_shape(type="line", x0=-in_log2fc, x1=-in_log2fc, y0=0, y1=vline_height,
                  line=dict(dash='dash', color='black'))
    fig.add_shape(type="line", x0=volcano_data['log2FC'].min(), x1=volcano_data['log2FC'].max(),
                  y0=-np.log10(in_pval), y1=-np.log10(in_pval),
                  line=dict(dash='dash', color='black'))
    fig.update_layout(
        xaxis_title=f"log2 fold change ({condition2} - {condition1})",
        yaxis_title="-log10 adj. p-value",
        hovermode='closest'
    )
    return fig


#PEPTIDE
@st.cache_data
def rt_vs_pred_rt_plot(data, method="Hexbin Plot", add_line=False, bins=1000, header=True, width=8, height=6, dpi=100):
    x = data['Predicted.RT']
    y = data['RT']

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    if method == "Scatter Plot":
        if len(data) > 100000:
            data = data.sample(n=100000, random_state=69)
        ax.scatter(data['Predicted.RT'], data['RT'], alpha=0.2, s=1)
    elif method == "Density Plot":
        from matplotlib.colors import LogNorm
        hb = ax.hist2d(data['Predicted.RT'], data['RT'], bins=200, norm=LogNorm())
        plt.colorbar(hb[3], ax=ax)
    elif method == "Hexbin Plot":
        hb = ax.hexbin(data['Predicted.RT'], data['RT'], gridsize=bins, cmap='Blues')
        plt.colorbar(hb, ax=ax)

    if add_line:
        ax.plot([x.min(), x.max()], [x.min(), x.max()], ls='--', color='red')

    ax.set_xlabel("Predicted RT")
    ax.set_ylabel("Actual RT")
    if header:
        ax.set_title("Retention Time Plot")

    return fig



@st.cache_data
def modification_plot(data2, meta, id=True, header=True, legend=True, width=10, height=6, dpi=100):
    data2 = data2.copy()
    meta = meta.copy()

    data2_names = data2['File.Name'].unique()
    for name in data2_names:
        meta['sample'] = np.where(meta['sample'] == name, name, meta['sample'])

    data2_wide = (
        data2.pivot_table(
            index="File.Name",
            columns="Modified.Sequence",
            values="Precursor.Quantity",
            aggfunc="max"
        )
        .reset_index()
    )

    data2_wide = data2_wide.set_index("File.Name").T.reset_index()
    data2_wide.rename(columns={"index": "Modified.Sequence"}, inplace=True)

    new_cols = [extract_id_or_number(c) for c in data2_wide.columns]
    data2_wide.columns = new_cols

    meta["sample"] = meta["sample"].astype(str)
    meta["id"] = meta["sample"].apply(extract_id_or_number)

    if id:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = (
            meta["condition"].astype(str).str.strip() + "_" +
            meta["new_sample"].astype(str).str.strip() +
            "\n (" + meta["id"].astype(str).str.strip() + ")"
        )
    else:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = (
            meta["condition"].astype(str).str.strip() + "_" +
            meta["new_sample"].astype(str).str.strip()
        )

    rename_dict = dict(zip(meta["id"], meta["new_sample"]))
    data2_wide = data2_wide.rename(columns=rename_dict)

    annotated_cols = [c for c in meta["new_sample"] if c in data2_wide.columns]
    data_filtered = data2_wide[["Modified.Sequence"] + annotated_cols]

    def get_mod_count(df, pattern):
        df_mod = df[df["Modified.Sequence"].str.contains(pattern, na=False)].copy()
        for col in df_mod.columns[1:]:
            df_mod[col] = df_mod[col].notna().astype(int)
        return df_mod.iloc[:, 1:].sum(axis=0)

    col_sums_carb = get_mod_count(data_filtered, "UniMod:4|Carbamidomethyl")
    col_sums_oxi  = get_mod_count(data_filtered, "UniMod:35|Oxidation")
    col_sums_ace  = get_mod_count(data_filtered, "UniMod:1|Acetyl")

    plot_data = pd.DataFrame({
        "Sample": list(col_sums_carb.index) * 3,
        "Count": list(col_sums_carb.values) + list(col_sums_oxi.values) + list(col_sums_ace.values),
        "Modification": (["Carbamylation"] * len(col_sums_carb)) +
                        (["Oxidation"] * len(col_sums_oxi)) +
                        (["Acetylation"] * len(col_sums_ace))
    })

    plot_data = plot_data[plot_data["Count"] > 0]
    plot_data["Sample"] = pd.Categorical(plot_data["Sample"], categories=meta["new_sample"], ordered=True)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    modifications = plot_data["Modification"].unique()
    x = np.arange(len(plot_data["Sample"].unique()))
    bar_width = 0.25
    offsets = np.linspace(-bar_width, bar_width, len(modifications))

    for i, mod in enumerate(modifications):
        subset = plot_data[plot_data["Modification"] == mod]
        ax.bar(x + offsets[i],
               subset["Count"],
               width=bar_width,
               label=mod)

    ax.set_xticks(x)
    ax.set_xticklabels(plot_data["Sample"].unique(), rotation=90)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Number of modified peptides")
    if header:
        ax.set_title("Modifications per sample")

    if legend:
        ax.legend(title="Modification", bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    plt.tight_layout()
    return fig


@st.cache_data
def missed_cleavage_plot(data2, meta, id=True, text=True, text_size=8, header=True, width=10, height=6, dpi=100):
    data2 = data2.copy()
    meta = meta.copy()

    data2_names = data2['File.Name'].unique()
    for name in data2_names:
        meta['sample'] = np.where(meta['sample'] == name, name, meta['sample'])

    df_wide = data2.pivot_table(
        index="File.Name",
        columns="Stripped.Sequence",
        values="Precursor.Quantity",
        aggfunc="max"
    ).T.reset_index()

    df_wide.rename(columns={"index": "Stripped.Sequence"}, inplace=True)
    df_wide.columns = [extract_id_or_number(c) if c != "Stripped.Sequence" else c for c in df_wide.columns]

    meta["sample"] = meta["sample"].astype(str)
    meta["id"] = meta["sample"].apply(extract_id_or_number)

    if id:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = (
            meta["condition"].astype(str) + "_" +
            meta["new_sample"].astype(str) + "\n(" +
            meta["id"].astype(str) + ")"
        )
    else:
        meta["new_sample"] = meta.groupby("condition").cumcount() + 1
        meta["new_sample"] = meta["condition"].astype(str) + "_" + meta["new_sample"].astype(str)

    rename_dict = dict(zip(meta["id"], meta["new_sample"]))
    df_wide = df_wide.rename(columns=rename_dict)

    annotated_cols = [c for c in meta["new_sample"] if c in df_wide.columns]
    data_filtered = df_wide[["Stripped.Sequence"] + annotated_cols]

    def count_all_RK(sequence):
        count_R = sequence.count("R")
        count_K = sequence.count("K")
        return max(0, count_R + count_K - 1)

    rk_counts = [count_all_RK(seq) for seq in data_filtered["Stripped.Sequence"]]

    for col in data_filtered.columns[1:]:
        data_filtered[col] = np.where(data_filtered[col].notna(), rk_counts, np.nan)

    plot_data = data_filtered.melt(id_vars="Stripped.Sequence", var_name="Sample", value_name="Count")
    plot_data = plot_data.dropna(subset=["Count"])
    plot_data_grouped = (
        plot_data.groupby(["Sample", "Count"])
        .size()
        .reset_index(name="Occurrences")
    )
    plot_data_grouped["Total"] = plot_data_grouped.groupby("Sample")["Occurrences"].transform("sum")
    plot_data_grouped["Percentage"] = plot_data_grouped["Occurrences"] / plot_data_grouped["Total"] * 100

    samples_order = meta["new_sample"].tolist()
    plot_data_grouped["Sample"] = pd.Categorical(plot_data_grouped["Sample"], categories=samples_order, ordered=True)
    plot_data_grouped["Count"] = plot_data_grouped["Count"].astype(int)

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    counts_sorted = sorted(plot_data_grouped["Count"].unique())
    bottom = np.zeros(len(samples_order))

    for c in counts_sorted:
        subset = (
            plot_data_grouped[plot_data_grouped["Count"] == c]
            .set_index("Sample")
            .reindex(samples_order)["Percentage"]
            .fillna(0)
        )
        ax.bar(samples_order, subset, bottom=bottom, label=str(c))
        bottom += subset

    ax.set_ylabel("Percentage of peptides with missed cleavage (%)")
    ax.set_ylim(0, max(bottom) * 1.05)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    if text:
        single_mc = (
            plot_data_grouped[plot_data_grouped["Count"] == 1]
            .set_index("Sample")
            .reindex(samples_order)["Percentage"]
            .fillna(0)
        )
        for i, val in enumerate(single_mc):
            if val > 0:
                ax.text(i, bottom[i] * 1, f"{val:.1f}%", ha="center", va="bottom", color="black", fontsize=text_size)

    if header:
        ax.set_title("Missed cleavages per sample")

    ax.set_xticks(range(len(samples_order)))
    ax.set_xticklabels(samples_order, rotation=90)

    ax.legend(title="Number", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    return fig


#SINGLE PROTEIN
@st.cache_data
def compare_prot_line(data, meta, conditions, inputs, id=True, header=True, legend=True,
                      workflow="Protein", plot_colors=None, width=10, height=6, dpi=100):
    meta = meta.copy()
    data = data.copy()
    meta["sample"] = meta["sample"].astype(str)

    if id:
        meta = (meta.groupby("condition", group_keys=False)
                    .apply(lambda g: g.assign(new_sample=[f"{g['condition'].iloc[i]}_{i+1}\n({g['sample'].iloc[i]})"
                                                          for i in range(len(g))]))
                    .reset_index(drop=True))
    else:
        meta = (meta.groupby("condition", group_keys=False)
                    .apply(lambda g: g.assign(new_sample=[f"{g['condition'].iloc[i]}_{i+1}"
                                                          for i in range(len(g))]))
                    .reset_index(drop=True))

    relevant_samples = meta["new_sample"].tolist()
    rename_vector = dict(zip(meta["sample"], meta["new_sample"]))
    data.rename(columns=rename_vector, inplace=True)

    if workflow == "Protein":
        data_filtered = data[data["ProteinNames"].isin(inputs)][["ProteinNames"] + relevant_samples]
        data_melted = data_filtered.melt(id_vars="ProteinNames", var_name="Sample", value_name="Value")
    elif workflow == "Phosphosite":
        data_filtered = data[data["PTM_Collapse_key"].isin(inputs)][["PTM_Collapse_key"] + relevant_samples]
        data_melted = data_filtered.melt(id_vars="PTM_Collapse_key", var_name="Sample", value_name="Value")
    else:
        raise ValueError("workflow must be 'Protein' or 'Phosphosite'")

    data_melted = data_melted.merge(meta[["new_sample", "condition"]],
                                    left_on="Sample", right_on="new_sample", how="left").dropna()

    ordered_samples = []
    for cond in conditions:
        ordered_samples.extend(meta.loc[meta["condition"] == cond, "new_sample"].tolist())

    sample_to_x = {s: i for i, s in enumerate(ordered_samples)}
    data_melted["x"] = data_melted["Sample"].map(sample_to_x)

    if workflow == "Protein":
        data_melted["ProteinNames"] = data_melted["ProteinNames"].apply(lambda x: x.split(";")[0])

    if plot_colors is None:
        import matplotlib.pyplot as plt
        plot_colors = plt.cm.tab10.colors

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)

    if workflow == "Protein":
        for i, prot in enumerate(data_melted["ProteinNames"].unique()):
            prot_data = data_melted[data_melted["ProteinNames"] == prot].sort_values("x")
            ax.plot(prot_data["x"], prot_data["Value"], marker="o", label=prot,
                    color=plot_colors[i % len(plot_colors)])
    elif workflow == "Phosphosite":
        for i, site in enumerate(data_melted["PTM_Collapse_key"].unique()):
            site_data = data_melted[data_melted["PTM_Collapse_key"] == site].sort_values("x")
            ax.plot(site_data["x"], site_data["Value"], marker="o", label=site,
                    color=plot_colors[i % len(plot_colors)])

    ax.set_xticks(range(len(ordered_samples)))
    ax.set_xticklabels(ordered_samples, rotation=90)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Log2 intensity")

    if header:
        ax.set_title(f"{workflow} Expression Across Samples")

    if legend:
        ax.legend(title=workflow, bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.grid(True)
    plt.tight_layout()
    fig.subplots_adjust(right=0.75)

    return fig


@st.cache_data
def boxplot_int_single_prot(data, meta, protein, outliers=False, header=True, legend=True,
                       plot_colors=None, width=8, height=6, dpi=100):
    meta = meta.copy()
    data = data.copy()
    data = data[data.index == protein]
    if data.empty:
        raise ValueError(f"No data found for protein: {protein}")

    meta = meta.groupby("condition", group_keys=False).apply(
        lambda g: g.assign(new_sample=[f"{str(g['condition'].iloc[i])}_{i + 1}" for i in range(len(g))])
    ).reset_index(drop=True)

    rename_vector = dict(zip(meta["sample"], meta["new_sample"]))
    data.rename(columns=rename_vector, inplace=True)
    annotated_columns = meta["new_sample"].tolist()
    data_filtered = data.loc[:, annotated_columns]

    if plot_colors is None:
        plot_colors = plt.cm.tab10.colors

    intensities = []
    labels = []
    colors = []
    plotted_conditions = []
    condition_colors = {cond: plot_colors[i % len(plot_colors)] for i, cond in enumerate(meta['condition'].unique())}

    for condition in meta["condition"].unique():
        cols = meta.loc[meta["condition"] == condition, "new_sample"]
        if len(cols) < 3:
            continue
        vals = data_filtered[cols].values.flatten()
        vals = vals[~pd.isna(vals)]
        if len(vals) == 0:
            continue
        intensities.append(vals)
        labels.append(condition)
        colors.append(condition_colors[condition])
        plotted_conditions.append(condition)

    if not intensities:
        fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
        ax.text(0.5, 0.5, "No conditions with ≥3 samples", ha="center", va="center")
        ax.axis("off")
        return fig

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    bplots = ax.boxplot(intensities, patch_artist=True, showfliers=outliers)
    for patch, color in zip(bplots['boxes'], colors):
        patch.set_facecolor(color)

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlabel("Condition")
    ax.set_ylabel("log2 Intensity")
    if header:
        display_protein = protein.split(";")[0]
        ax.set_title(f"Measured {display_protein} intensity values (log2)")
    if legend:
        handles = [plt.Line2D([0], [0], color=condition_colors[c], lw=4) for c in plotted_conditions]
        ax.legend(handles, plotted_conditions, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    return fig


#Phospho-specific
@st.cache_data
def phossite_coverage_plot(data, meta, id=False, header=True, legend=True, width=20, height=10, dpi=300):
    conditions = meta["condition"].unique()
    meta["sample"] = meta["sample"].astype(str)
    meta["id"] = meta["sample"].apply(lambda x: extract_id_or_number(x))

    if id:
        meta["new_sample"] = meta.groupby("condition").cumcount().astype(str)
        meta["new_sample"] = meta["condition"].astype(str) + "_" + meta["new_sample"] + "\n(" + meta["id"].astype(str) + ")"
    else:
        meta["new_sample"] = meta.groupby("condition").cumcount().astype(str)
        meta["new_sample"] = meta["condition"].astype(str) + "_" + meta["new_sample"]

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data = data.rename(columns=rename_dict)
    annotated_columns = meta["new_sample"].tolist()

    data_classI = data.loc[data["PTM_localization"] >= 0.75, annotated_columns]
    data_notclassI = data.loc[data["PTM_localization"] < 0.75, annotated_columns]

    count_classI = data_classI.notna().sum()
    count_notclassI = data_notclassI.notna().sum()

    plot_data = pd.DataFrame({
        "sample": annotated_columns * 2,
        "PTM_localization": ["Class I"] * len(annotated_columns) + ["Not Class I"] * len(annotated_columns),
        "count": list(count_classI) + list(count_notclassI)
    })

    plot_data["sample"] = pd.Categorical(plot_data["sample"], categories=meta["new_sample"], ordered=True)

    plot_title = "Phosphosites per sample" if header else ""
    y_label = "Number of phosphosites" if header else "Number"

    fig, ax = plt.subplots(figsize=(width, height), dpi=dpi)
    classI_data = plot_data[plot_data["PTM_localization"] == "Class I"]
    notclassI_data = plot_data[plot_data["PTM_localization"] == "Not Class I"]

    x = range(len(annotated_columns))
    bars1 = ax.bar(x, classI_data["count"], color="blue", label="Class I")
    bars2 = ax.bar(x, notclassI_data["count"], bottom=classI_data["count"], color="orange", label="Not Class I")

    ax.set_xticks(x)
    ax.set_xticklabels(meta["new_sample"], rotation=90, fontsize=5)
    ax.set_xlabel("Sample")
    ax.set_ylabel(y_label)
    ax.set_title(plot_title)

    if legend:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    return fig


@st.cache_data
def phossite_coverage_plot_summary(data, meta, header=True, legend=True, width=20, height=10, dpi=300):
    meta = meta.copy()
    meta["sample"] = meta["sample"].astype(str)

    meta["new_sample"] = meta.groupby("condition").cumcount().add(1).astype(str)
    meta["new_sample"] = meta.apply(lambda row: f'{row.condition}_{row["new_sample"]}', axis=1)

    rename_vector = dict(zip(meta["sample"], meta["new_sample"]))
    data = data.copy()
    data.columns = [rename_vector.get(col, col) for col in data.columns]
    annotated_columns = meta["new_sample"]

    data_classI = data[data["PTM_localization"] >= 0.75][annotated_columns]
    data_notclassI = data[data["PTM_localization"] < 0.75][annotated_columns]

    count_classI = data_classI.notna().sum()
    count_notclassI = data_notclassI.notna().sum()

    plot_data = pd.DataFrame({
        "sample": list(annotated_columns) * 2,
        "PTM_localization": ["Class I"] * len(annotated_columns) + ["Not Class I"] * len(annotated_columns),
        "count": list(count_classI) + list(count_notclassI)
    })
    plot_data["condition"] = plot_data["sample"].str.extract(r'(^.*?)_')[0]

    summarized_plot_data = plot_data.groupby(["condition", "PTM_localization"]).agg(
        mean_count=("count", "mean"),
        sd_count=("count", "std")
    ).reset_index()

    plot_title = "Phosphosites per condition" if header else ""
    y_label = "Number of phosphosites" if header else "Number"

    fig, ax = plt.subplots(figsize=(width / 2.54, height / 2.54), dpi=dpi)

    conditions = summarized_plot_data["condition"].unique()
    locs = np.arange(len(conditions))
    width_bar = 0.35

    classI_means = summarized_plot_data[summarized_plot_data["PTM_localization"] == "Class I"]["mean_count"].values
    classI_sd = summarized_plot_data[summarized_plot_data["PTM_localization"] == "Class I"]["sd_count"].values
    notclassI_means = summarized_plot_data[summarized_plot_data["PTM_localization"] == "Not Class I"]["mean_count"].values
    notclassI_sd = summarized_plot_data[summarized_plot_data["PTM_localization"] == "Not Class I"]["sd_count"].values

    bars1 = ax.bar(locs - width_bar/2, classI_means, width=width_bar, yerr=classI_sd, capsize=3, color="blue", label="Class I")
    bars2 = ax.bar(locs + width_bar/2, notclassI_means, width=width_bar, yerr=notclassI_sd, capsize=3, color="orange", label="Not Class I")

    for i, condition in enumerate(conditions):
        y_classI = plot_data[(plot_data["condition"] == condition) & (plot_data["PTM_localization"] == "Class I")]["count"].values
        y_notclassI = plot_data[(plot_data["condition"] == condition) & (plot_data["PTM_localization"] == "Not Class I")]["count"].values
        ax.scatter(np.random.normal(i - width_bar/2, 0.05, len(y_classI)), y_classI, color='black', alpha=0.7)
        ax.scatter(np.random.normal(i + width_bar/2, 0.05, len(y_notclassI)), y_notclassI, color='black', alpha=0.7)

    ax.set_xticks(locs)
    ax.set_xticklabels(conditions, rotation=45, ha='right')
    ax.set_ylabel(y_label)
    ax.set_title(plot_title)

    if legend:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    return fig


@st.cache_data
def simple_phos_site_plot(data, filter_value=0, width_cm=15, height_cm=10, dpi=100):
    data = data[data["PTM_localization"] >= filter_value].copy()
    phosprot = data["Protein_group"].nunique()
    data["Seq"] = data["UPD_seq"].astype(str).str.replace(r"[^\w]", "", regex=True).str.upper()
    phospep = data["Seq"].nunique()
    phossite = data["PTM_Collapse_key"].nunique()
    df = pd.DataFrame({
        "term": ["Phosphoproteins", "Phosphopeptides", "Phosphosites"],
        "count": [phosprot, phospep, phossite]
    }).sort_values("count", ascending=False)
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54
    fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=dpi)
    ax.barh(df["term"], df["count"], color="skyblue")
    for i, v in enumerate(df["count"]):
        ax.text(v + max(df["count"]) * 0.01, i, str(v), color="black", fontweight="bold", va="center")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.grid(False)
    plt.tight_layout()
    return fig