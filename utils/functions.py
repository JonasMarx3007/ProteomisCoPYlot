import pandas as pd
import pyarrow.parquet as pq
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind, ttest_rel
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go

#BASIC AND DATA FUNCTIONS
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


#QC FUNCTIONS
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


def corr_plot(data, meta, method=False, id=True, full_range=False, width=10, height=8, dpi=100):
    meta['sample'] = meta['sample'].astype(str)
    meta['id'] = meta['sample'].apply(extract_id_or_number)
    meta = meta.copy()
    if id:
        meta['new_sample'] = meta.groupby('condition').cumcount().add(1).astype(str)
        meta['new_sample'] = meta.apply(lambda row: f"{row['condition']}_{row.name+1}\n ({row['id']})", axis=1)
    else:
        meta['new_sample'] = meta.groupby('condition').cumcount().add(1).astype(str)
        meta['new_sample'] = meta.apply(lambda row: f"{row['condition']}_{row.name+1}", axis=1)
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
    cax = ax.matshow(ordered_corr, cmap='coolwarm', vmin=vmin, vmax=vmax)
    plt.xticks(range(len(ordered_corr.columns)), ordered_corr.columns, rotation=45, ha='left')
    plt.yticks(range(len(ordered_corr.index)), ordered_corr.index)
    fig.colorbar(cax)
    plt.tight_layout()
    return fig


#STATISTIC FUNCTIONS
def volcano_plot(data, meta, condition1, condition2, in_pval=0.05, in_log2fc=1,
                 workflow="Protein", paired="Unpaired", uncorrected=False):
    annotated_columns1 = meta.loc[meta['condition'] == condition1, 'sample'].tolist()
    annotated_columns2 = meta.loc[meta['condition'] == condition2, 'sample'].tolist()

    label_col = "ProteinNames" if workflow == "Protein" else "PTM_Collapse_key"
    cols_to_use = [label_col] + annotated_columns1 + annotated_columns2
    data_filtered = data[cols_to_use].copy()
    data_filtered = data_filtered[data_filtered.iloc[:, 1:].notna().sum(axis=1) >= 2]

    log2fc_list = []
    pval_list = []

    for _, row in data_filtered.iterrows():
        vals1 = row[annotated_columns1].astype(float).dropna()
        vals2 = row[annotated_columns2].astype(float).dropna()

        if paired == "Unpaired":
            if len(vals1) < 2 or len(vals2) < 2:
                log2fc_list.append(np.nan)
                pval_list.append(np.nan)
                continue
            log2fc_list.append(vals2.mean() - vals1.mean())
            pval_list.append(ttest_ind(vals2, vals1, equal_var=True).pvalue)
        else:  # Paired
            paired_vals = pd.DataFrame({'v1': vals1, 'v2': vals2}).dropna()
            if len(paired_vals) < 2:
                log2fc_list.append(np.nan)
                pval_list.append(np.nan)
                continue
            log2fc_list.append(paired_vals['v2'].mean() - paired_vals['v1'].mean())
            pval_list.append(ttest_rel(paired_vals['v2'], paired_vals['v1']).pvalue)

    volcano_data = data_filtered.copy()
    volcano_data['log2FC'] = log2fc_list
    volcano_data['pval'] = pval_list
    volcano_data = volcano_data.dropna(subset=['pval'])

    volcano_data['adj_pval'] = (volcano_data['pval'] if uncorrected
                                else multipletests(volcano_data['pval'], method='fdr_bh')[1])
    volcano_data['neg_log10_pval'] = -np.log10(volcano_data['adj_pval'])
    volcano_data['significance'] = np.where(
        (volcano_data['adj_pval'] < in_pval) & (volcano_data['log2FC'] > in_log2fc), 'Upregulated',
        np.where((volcano_data['adj_pval'] < in_pval) & (volcano_data['log2FC'] < -in_log2fc), 'Downregulated',
                 'Not significant')
    )

    color_mapping = {'Downregulated': 'blue', 'Not significant': 'gray', 'Upregulated': 'red'}
    y_hline = -np.log10(in_pval)
    max_y = max(y_hline, volcano_data['neg_log10_pval'].max())

    shapes = [
        dict(type="line", x0=in_log2fc, x1=in_log2fc, y0=0, y1=max_y, line=dict(color="black", dash="dash")),
        dict(type="line", x0=-in_log2fc, x1=-in_log2fc, y0=0, y1=max_y, line=dict(color="black", dash="dash")),
        dict(type="line", x0=volcano_data['log2FC'].min(), x1=volcano_data['log2FC'].max(), y0=y_hline, y1=y_hline,
             line=dict(color="black", dash="dash"))
    ]

    fig = go.Figure()
    for sig in color_mapping:
        subset = volcano_data[volcano_data['significance'] == sig]
        fig.add_trace(go.Scatter(
            x=subset['log2FC'],
            y=subset['neg_log10_pval'],
            mode='markers',
            marker=dict(color=color_mapping[sig], size=5),
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


def volcano_data_f(data, meta, condition1, condition2, in_pval=0.05, in_log2fc=1,
                   workflow="Protein", paired="Unpaired", uncorrected=False):
    data = data.replace(0, np.nan)

    annotated_columns1 = meta.loc[meta['condition'] == condition1, 'sample'].tolist()
    annotated_columns2 = meta.loc[meta['condition'] == condition2, 'sample'].tolist()

    if workflow == "Protein":
        data_filtered = data[['ProteinNames'] + annotated_columns1 + annotated_columns2].copy()
    elif workflow == "Phosphosite":
        data_filtered = data[['PTM_Collapse_key'] + annotated_columns1 + annotated_columns2].copy()

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

    if workflow == "Protein":
        volcano_df = pd.DataFrame({
            'Protein': data_filtered['ProteinNames'],
            'log2FC': log2fc_array,
            'pval': pvals_array,
            'adj_pval': adj_pvals
        })
    elif workflow == "Phosphosite":
        volcano_df = pd.DataFrame({
            'Phossite': data_filtered['PTM_Collapse_key'],
            'log2FC': log2fc_array,
            'pval': pvals_array,
            'adj_pval': adj_pvals
        })

    volcano_df['significance'] = np.where(
        (volcano_df['adj_pval'] < in_pval) & (volcano_df['log2FC'] > in_log2fc), "Upregulated",
        np.where((volcano_df['adj_pval'] < in_pval) & (volcano_df['log2FC'] < -in_log2fc), "Downregulated",
                 "Not significant")
    )

    volcano_df['neg_log10_adj_pval'] = -np.log10(volcano_df['adj_pval'])

    return volcano_df