import streamlit as st

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
        pval_threshold = st.number_input("P‑value threshold:", value=0.05, step=0.01, key="in_pval10")
        log2fc_threshold = st.number_input("log₂ FC threshold:", value=1.0, step=0.1, key="in_log2fc10")
        uncorrected = st.checkbox("Use uncorrected p‑values", value=False, key="uncorrected10")

        st.markdown("---")

        test_type = st.selectbox("Test Type:", options=["Unpaired", "Paired"], key="paired10")
        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level10")

        st.markdown("---")

        st.download_button("Download Volcano Data", data=b"", key="download10")
        st.download_button("Download All Volcano Data", data=b"", key="download10a")

        st.markdown("---")

        text_position = st.selectbox("Position:", options={"Above": "up", "Below": "down"}, index=0, key="textPosition10")
        annotation_text = st.text_area("Annotation Text:", value="", height=100, key="text10")

        add_col, del_col = st.columns(2)
        with add_col:
            st.button("Add", key="addText10")
        with del_col:
            st.button("Delete", key="deleteText10")

        st.markdown("---")

        st.button("Add Plot to Report", key="addVolc")
        st.button("Remove Plot from Report", key="removeVolc")

    with col2:
        st.subheader("Conditions")
        condition1 = st.selectbox("Condition 1:", options=[], key="condition1_10")
        condition2 = st.selectbox("Condition 2:", options=[], key="condition2_10")

        st.markdown("---")

        st.subheader("Volcano Plot")
        st.empty()

        st.markdown("---")

        st.subheader("Data Table")
        st.empty()


def gsea_ui():
    st.header("GSEA")

    col1, col2 = st.columns([1, 2])

    with col1:
        top_n = st.number_input("Number of top gene sets to show:", min_value=1, value=10, step=1, key="top_n11")

        st.markdown("---")

        st.subheader("Term Size Filter")
        min_term = st.number_input("Min term size:", value=20, key="filter11min")
        max_term = st.number_input("Max term size:", value=300, key="filter11max")

        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth11")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight11")
        dpi = st.number_input("DPI:", value=300, key="plotDPI11")

        st.markdown("---")

        st.download_button("Download Plot (UP)", data=b"", key="downloadUpPlot")
        st.download_button("Download Plot (DOWN)", data=b"", key="downloadDownPlot")

    with col2:
        st.subheader("Upregulated Gene Sets")
        st.empty()

        st.markdown("---")

        st.subheader("Downregulated Gene Sets")
        st.empty()

        st.markdown("---")

        list_col1, list_col2 = st.columns(2)
        with list_col1:
            st.subheader("Upregulated Gene List")
            st.empty()
        with list_col2:
            st.subheader("Downregulated Gene List")
            st.empty()


def simulation_ui():
    st.header("Simulation")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Thresholds")
        pval_threshold = st.number_input("P‑value threshold:", value=0.05, step=0.01, key="in_pval10.5")
        log2fc_threshold = st.number_input("log₂ FC threshold:", value=1.0, step=0.1, key="in_log2fc10.5")

        st.markdown("---")

        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level10.5")

        st.markdown("---")

        st.slider("Variance multiplier:", min_value=0.25, max_value=4.0, value=1.0, step=0.05, key="mod_var10.5")
        st.slider("Sample size override (n):", min_value=1, max_value=20, value=1, step=1, key="mod_n10.5")

    with col2:
        st.subheader("Conditions")
        condition1 = st.selectbox("Condition 1:", options=[], key="condition1_10.5")
        condition2 = st.selectbox("Condition 2:", options=[], key="condition2_10.5")

        st.markdown("---")

        st.subheader("Simulation Volcano Plot")
        st.empty()