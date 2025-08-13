import pandas as pd
import pyarrow.parquet as pq
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import PCA

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
    import re
    m = re.search(r'\d+', x)
    return m.group(0) if m else x


def log2_transform_data(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    annotated_columns = meta['sample'].tolist()
    data_filtered = data[annotated_columns].copy()
    data_filtered.replace(0, np.nan, inplace=True)
    log2_data = np.log2(data_filtered)

    remaining_columns = [col for col in data.columns if col not in annotated_columns]
    combined_data = pd.concat([log2_data, data[remaining_columns]], axis=1)

    return combined_data


def inverse_log2_transform_data(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    annotated_columns = meta['sample'].tolist()
    log2_data = data[annotated_columns].copy()
    original_data = 2 ** log2_data

    remaining_columns = [col for col in data.columns if col not in annotated_columns]
    combined_data = pd.concat([original_data, data[remaining_columns]], axis=1)

    return combined_data


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


def histo_int(data, meta, plot_colors=None, header=True, legend=True, ax=None):
    import numpy as np
    from scipy.stats import gaussian_kde
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

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


def boxplot_int(data, meta, outliers=False, header=True, legend=True, plot_colors=None,
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

    if legend:
        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor=plot_colors[i], edgecolor='black', label=cond) for i, cond in enumerate(conditions_order)]
        ax.legend(handles=legend_handles, title="Condition",
                  loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    return fig


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
        ax.legend(handles=legend_handles, title="Condition")
    else:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    plt.tight_layout()
    return fig


def pca_plot(data, meta, header=True, legend=True, dot_size=3, width_cm=20, height_cm=10, dpi=100):
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

    colors = plt.cm.tab10.colors
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
