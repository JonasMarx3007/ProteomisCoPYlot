import streamlit as st
from utils.functions import compare_prot_line, boxplot_int_single

#MAIN
def single_protein_ui():
    sp_tabs = st.tabs([
        "Protein Lineplot",
        "Protein Boxplot"
    ])

    with sp_tabs[0]:
        protein_line_ui()
    with sp_tabs[1]:
        protein_box_ui()


#SUB
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


def protein_line_ui():
    plot_colors = st.session_state.get("plot_colors", None)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.checkbox("Toggle IDs", value=False, key="toggle_id18")
        st.checkbox("Toggle Header", value=True, key="toggle_header18")
        st.checkbox("Toggle Legend", value=True, key="toggle_legend18")

        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level18")
        protein_col = "ProteinNames" if level == "Protein" else "PTM_Collapse_key"

        if "log2_data" in st.session_state:
            log2_data = st.session_state.log2_data
            protein_options = sorted(log2_data[protein_col].dropna().unique().tolist())
        else:
            protein_options = []

        selected_proteins = st.multiselect("Select Proteins:", options=protein_options, key="protein18")

        if "meta" in st.session_state:
            meta = st.session_state.meta
            condition_options = sorted(meta["condition"].unique().tolist())
        else:
            meta = None
            condition_options = []

        selected_conditions = st.multiselect(
            "Select Conditions:",
            options=condition_options,
            default=condition_options,
            key="conditions18"
        )

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="plotWidth18")
        st.number_input("Height (cm):", value=10, key="plotHeight18")
        st.number_input("DPI:", value=300, key="plotDPI18")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat18")
        st.download_button("Download Plot", data="", file_name="protein_lineplot.png", key="downloadProteinLine")

    with col2:
        if "log2_data" in st.session_state and "meta" in st.session_state:
            if selected_proteins and selected_conditions:
                fig = compare_prot_line(
                    data=st.session_state.log2_data,
                    meta=st.session_state.meta,
                    conditions=selected_conditions,
                    inputs=selected_proteins,
                    id=st.session_state.get("toggle_id18", True),
                    header=st.session_state.get("toggle_header18", True),
                    legend=st.session_state.get("toggle_legend18", True),
                    workflow=level,
                    plot_colors=plot_colors,
                    width=st.session_state.get("plotWidth18", 20) / 2.54,
                    height=st.session_state.get("plotHeight18", 10) / 2.54,
                    dpi=st.session_state.get("plotDPI18", 300)
                )
                st.pyplot(fig)
            else:
                st.info("Select at least one protein and one condition.")
        else:
            st.info("Log2_data or meta not available. UI controls are active.")


def protein_box_ui():
    plot_colors = st.session_state.get("plot_colors", None)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.checkbox("Toggle Header", value=True, key="toggle_header17")
        st.checkbox("Toggle Legend", value=True, key="toggle_legend17")
        st.checkbox("Show Outliers", value=False, key="toggle_outliers17")
        st.markdown("---")

        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level17")
        protein_col = "ProteinNames" if level == "Protein" else "PTM_Collapse_key"

        if "log2_data" in st.session_state:
            log2_data = st.session_state.log2_data
            protein_options = sorted(log2_data[protein_col].dropna().unique().tolist())
        else:
            protein_options = []

        selected_protein = st.selectbox("Select Protein:", options=protein_options, key="protein17")

        if "meta" in st.session_state:
            meta = st.session_state.meta
            condition_options = sorted(meta["condition"].unique().tolist())
        else:
            meta = None
            condition_options = []

        selected_conditions = st.multiselect(
            "Select Conditions:",
            options=condition_options,
            default=condition_options,
            key="conditions17"
        )

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="plotWidth17")
        st.number_input("Height (cm):", value=10, key="plotHeight17")
        st.number_input("DPI:", value=300, key="plotDPI17")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat17")
        st.download_button("Download Plot", data="", file_name="protein_boxplot.png", key="downloadProteinBox")

    with col2:
        if "log2_data" in st.session_state and "meta" in st.session_state:
            if selected_protein and selected_conditions:
                meta_filtered = meta[meta["condition"].isin(selected_conditions)]
                data_filtered = log2_data[log2_data[protein_col] == selected_protein]
                data_filtered = data_filtered.set_index(protein_col)
                selected_samples = meta_filtered["sample"].tolist()
                data_filtered = data_filtered[selected_samples]

                fig = boxplot_int_single(
                    data_filtered,
                    meta_filtered,
                    protein=selected_protein,
                    outliers=st.session_state.get("toggle_outliers17", False),
                    header=st.session_state.get("toggle_header17", True),
                    legend=st.session_state.get("toggle_legend17", True),
                    plot_colors=plot_colors,
                    width=st.session_state.get("plotWidth17", 20) / 2.54,
                    height=st.session_state.get("plotHeight17", 10) / 2.54,
                    dpi=st.session_state.get("plotDPI17", 300)
                )
                st.pyplot(fig)
            else:
                st.info("Select a protein and at least one condition.")
        else:
            st.info("Log2_data or meta not available. UI controls are active.")
