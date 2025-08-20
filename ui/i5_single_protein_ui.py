import streamlit as st
from utils.functions import compare_prot_line, boxplot_int_single
from io import BytesIO

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
