import streamlit as st
from utils.functions import *

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
    log2_data = st.session_state.log2_data
    meta = st.session_state.meta

    col1, col2 = st.columns([1, 2])

    with col1:
        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level18")
        if level == "Protein":
            protein_col = "ProteinNames"
        else:
            protein_col = "PTM_Collapse_key"

        protein_options = sorted(log2_data[protein_col].dropna().unique().tolist())
        selected_proteins = st.multiselect("Select Proteins:", options=protein_options, key="protein18")

        condition_options = sorted(meta["condition"].unique().tolist())
        selected_conditions = st.multiselect("Select Conditions:",
                                             options=condition_options,
                                             default=condition_options,
                                             key="conditions18")

    with col2:
        if selected_proteins and selected_conditions:
            fig = compare_prot_line(
                data=log2_data,
                meta=meta,
                conditions=selected_conditions,
                inputs=selected_proteins,
                id=True,
                workflow=level
            )
            st.pyplot(fig)
        else:
            st.info("Select at least one protein and one condition.")


def protein_box_ui():
    log2_data = st.session_state.log2_data
    meta = st.session_state.meta
    plot_colors = st.session_state.get("plot_colors", None)

    col1, col2 = st.columns([1, 2])

    with col1:
        level = st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level17")
        protein_col = "ProteinNames" if level == "Protein" else "PTM_Collapse_key"

        protein_options = sorted(log2_data[protein_col].dropna().unique().tolist())
        selected_protein = st.selectbox("Select Protein:", options=protein_options, key="protein17")

        condition_options = sorted(meta["condition"].unique().tolist())
        selected_conditions = st.multiselect(
            "Select Conditions:", options=condition_options, default=condition_options, key="conditions17"
        )

    with col2:
        if selected_protein and selected_conditions:
            meta_filtered = meta[meta["condition"].isin(selected_conditions)]
            data_filtered = log2_data[log2_data[protein_col] == selected_protein]
            selected_samples = meta_filtered["sample"].tolist()
            data_filtered = data_filtered[selected_samples]

            fig = boxplot_int_single(
                data_filtered,
                meta_filtered,
                protein="ZZ",
                outliers=False,
                id=True,
                header=True,
                legend=True,
                plot_colors=plot_colors
            )
            st.pyplot(fig)
        else:
            st.info("Select a protein and at least one condition.")
