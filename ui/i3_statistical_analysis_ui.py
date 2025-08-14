import streamlit as st
from utils.functions import volcano_plot, volcano_plot_sim, volcano_data_f, enrichment_analysis, different_genes

#MAIN
def statistical_analysis_ui():
    stat_tabs = st.tabs([
        "Volcano Plot",
        "GSEA",
        "Simulation"
    ])

    with stat_tabs[0]:
        volcano_plot_ui()
    with stat_tabs[1]:
        gsea_ui()
    with stat_tabs[2]:
        simulation_ui()


#SUB
def volcano_plot_ui():
    st.header("Volcano Plot")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Thresholds")
        pval_threshold = st.number_input("P-value threshold:", value=0.05, step=0.01, key="in_pval10")
        log2fc_threshold = st.number_input("log₂ FC threshold:", value=1.0, step=0.1, key="in_log2fc10")
        uncorrected = st.checkbox("Use uncorrected p-values", value=False, key="uncorrected10")
        st.markdown("---")
        test_type = st.selectbox("Test Type:", options=["Unpaired", "Paired"], key="paired10")

        available_levels = []
        if "log2_data" in st.session_state and "meta" in st.session_state:
            available_levels.append("Protein")
        if "log2_data3" in st.session_state and "meta2" in st.session_state:
            available_levels.append("Phosphosite")

        if not available_levels:
            st.warning("No data loaded. Please load Protein or Phosphosite data.")
            return

        level = st.selectbox("Level:", options=available_levels, key="level10")

    with col2:
        st.subheader("Conditions")

        if level == "Protein":
            data_key, meta_key = "log2_data", "meta"
        else:
            data_key, meta_key = "log2_data3", "meta2"

        conditions_list = []
        if data_key in st.session_state and meta_key in st.session_state:
            meta_df = st.session_state[meta_key]
            if not meta_df.empty:
                counts = meta_df['condition'].value_counts()
                conditions_list = sorted(counts[counts >= 2].index.tolist())

        condition1 = st.selectbox("Condition 1:", options=conditions_list, key="condition1_10")
        condition2 = st.selectbox("Condition 2:", options=conditions_list, key="condition2_10")

        st.markdown("---")
        st.subheader("Volcano Plot")

        if condition1 and condition2 and condition1 != condition2:
            volcano_df = volcano_data_f(
                st.session_state[data_key],
                st.session_state[meta_key],
                condition1=condition1,
                condition2=condition2,
                in_pval=pval_threshold,
                in_log2fc=log2fc_threshold,
                workflow=level,
                paired=test_type,
                uncorrected=uncorrected
            )

            fig = volcano_plot(
                volcano_df,
                condition1=condition1,
                condition2=condition2,
                in_pval=pval_threshold,
                in_log2fc=log2fc_threshold,
                uncorrected=uncorrected
            )

            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(volcano_df)
        else:
            st.info("Please select two different conditions to generate the volcano plot.")


def gsea_ui():
    st.header("GSEA")

    if "log2_data" not in st.session_state or "meta" not in st.session_state:
        st.warning("Please load data and generate a volcano plot first.")
        return

    if "GeneNames" not in st.session_state.log2_data.columns:
        st.warning("The loaded data does not contain a 'GeneNames' column, which is required for GSEA.")
        return

    col1, col2 = st.columns([1, 2])

    with col1:
        top_n = st.number_input("Number of top gene sets to show:", min_value=1, value=10, step=1, key="top_n11")
        st.markdown("---")
        st.subheader("Term Size Filter")
        min_term = st.number_input("Min term size:", value=20, key="min_term11")
        max_term = st.number_input("Max term size:", value=300, key="max_term11")
        st.markdown("---")
        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width (cm):", value=20, key="plot_width11")
        plot_height = st.number_input("Height (cm):", value=10, key="plot_height11")
        dpi = st.number_input("DPI:", value=300, key="dpi11")
        st.markdown("---")
        st.download_button("Download Plot (UP)", data=b"", key="downloadUpPlot11")
        st.download_button("Download Plot (DOWN)", data=b"", key="downloadDownPlot11")

    with col2:
        st.subheader("Upregulated Gene Sets")
        up_genes_placeholder = st.empty()
        st.markdown("---")
        st.subheader("Downregulated Gene Sets")
        down_genes_placeholder = st.empty()
        st.markdown("---")
        list_col1, list_col2 = st.columns(2)
        with list_col1:
            st.subheader("Upregulated Gene List")
            up_list_placeholder = st.empty()
        with list_col2:
            st.subheader("Downregulated Gene List")
            down_list_placeholder = st.empty()

    condition1 = st.session_state.get("condition1_10", None)
    condition2 = st.session_state.get("condition2_10", None)
    in_pval = st.session_state.get("in_pval10", 0.05)
    in_log2fc = st.session_state.get("in_log2fc10", 1.0)
    use_uncorrected = st.session_state.get("uncorrected10", False)

    if condition1 and condition2 and condition1 != condition2:
        diff_res = different_genes(
            st.session_state.log2_data,
            st.session_state.meta,
            condition1,
            condition2,
            in_pval=in_pval,
            in_log2fc=in_log2fc,
            uncorrected=use_uncorrected
        )
        up_genes = diff_res["Upregulated"]
        down_genes = diff_res["Downregulated"]

        up_list_placeholder.write(up_genes)
        down_list_placeholder.write(down_genes)

        if up_genes:
            up_fig = enrichment_analysis(up_genes, top_n=top_n, min_num=min_term, max_num=max_term)
            if up_fig:
                up_genes_placeholder.pyplot(up_fig, dpi=dpi)
        if down_genes:
            down_fig = enrichment_analysis(down_genes, top_n=top_n, min_num=min_term, max_num=max_term)
            if down_fig:
                down_genes_placeholder.pyplot(down_fig, dpi=dpi)
    else:
        st.info("Please select two different conditions in the volcano plot to run GSEA.")


def simulation_ui():
    st.header("Simulation")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Thresholds")
        pval_threshold = st.number_input("P-value threshold:", value=0.05, step=0.01, key="in_pval10.5")
        log2fc_threshold = st.number_input("log₂ FC threshold:", value=1.0, step=0.1, key="in_log2fc10.5")

        st.markdown("---")

        available_levels = []
        if "log2_data" in st.session_state and "meta" in st.session_state:
            available_levels.append("Protein")
        if "log2_data3" in st.session_state and "meta2" in st.session_state:
            available_levels.append("Phosphosite")

        if not available_levels:
            st.warning("No data available for simulation.")
            return

        level = st.selectbox("Level:", options=available_levels, key="level10.5")

        st.markdown("---")

        mod_var = st.slider("Variance multiplier:", min_value=0.25, max_value=4.0, value=1.0, step=0.05, key="mod_var10.5")
        mod_n = st.slider("Sample size override (n):", min_value=1, max_value=20, value=1, step=1, key="mod_n10.5")

    with col2:
        st.subheader("Conditions")

        if level == "Protein":
            meta_df = st.session_state.meta if "meta" in st.session_state else None
        else:
            meta_df = st.session_state.meta2 if "meta2" in st.session_state else None

        conditions_list = []
        if meta_df is not None and not meta_df.empty:
            counts = meta_df['condition'].value_counts()
            conditions_list = sorted(counts[counts >= 2].index.tolist())

        condition1 = st.selectbox("Condition 1:", options=conditions_list, key="condition1_10.5")
        condition2 = st.selectbox("Condition 2:", options=conditions_list, key="condition2_10.5")

        st.markdown("---")
        st.subheader("Simulation Volcano Plot")

        data_key = "log2_data" if level == "Protein" else "log2_data3"
        meta_key = "meta" if level == "Protein" else "meta2"

        if condition1 and condition2 and condition1 != condition2:
            if data_key in st.session_state and meta_key in st.session_state:
                fig = volcano_plot_sim(
                    data=st.session_state[data_key],
                    meta=st.session_state[meta_key],
                    condition1=condition1,
                    condition2=condition2,
                    in_pval=pval_threshold,
                    in_log2fc=log2fc_threshold,
                    workflow=level,
                    mod_var=mod_var,
                    mod_n=mod_n
                )
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("No data loaded for the selected level.")
        else:
            st.info("Please select two different conditions to generate the volcano plot.")
