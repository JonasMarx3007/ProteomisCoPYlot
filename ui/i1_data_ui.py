import streamlit as st
from utils.functions import *

#MAIN
def data_ui():
    data_tabs = st.tabs([
        "Data Upload",
        "Data Annotation",
        "Impute Data",
        "Distribution",
        "Verification"
    ])

    with data_tabs[0]:
        data_upload_ui()
    with data_tabs[1]:
        data_annotation_ui()
    with data_tabs[2]:
        impute_data_ui()
    with data_tabs[3]:
        distribution_ui()
    with data_tabs[4]:
        verification_ui()


#SUB
def data_upload_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        st.header("Upload Files")
        st.file_uploader("Protein Level Data", type=["csv", "tsv", "txt", "xlsx", "parquet"], key="file")
        st.file_uploader("Phospho Data", type=["csv", "tsv", "txt", "xlsx", "parquet"], key="file3")
        st.file_uploader("Full Report", type=["csv", "tsv", "txt", "xlsx", "parquet"], key="file2")

        st.markdown("---")
        st.header("Collapse Options")
        st.radio("Select Collapse Option:", ["Collapsed", "Not Collapsed"], index=0, key="collapse_option")
        st.number_input("Collapse Cutoff:", value=0, key="collapse_cutoff")

    with col2:
        st.subheader("Protein Level Data Preview")
        if st.session_state.get("file"):
            st.session_state["data"] = read_data(st.session_state["file"])
            st.dataframe(st.session_state["data"])
        else:
            st.dataframe([], key="table1")

        st.markdown("---")
        st.subheader("Phospho Data Preview")
        if st.session_state.get("file3"):
            st.session_state["data3"] = read_data(st.session_state["file3"])
            st.dataframe(st.session_state["data3"])
        else:
            st.dataframe([], key="table3")

        st.markdown("---")
        st.subheader("Full Report Preview")
        if st.session_state.get("file2"):
            st.session_state["data2"] = read_data(st.session_state["file2"])
            st.dataframe(st.session_state["data2"])
        else:
            st.dataframe([], key="table2")


def data_annotation_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.header("Normal Data Annotation")
        meta_file = st.file_uploader(
            "Upload Metadata (Protein Group)",
            type=["csv", "xlsx", "txt", "tsv"],
            key="upload_meta"
        )
        if meta_file is not None:
            try:
                st.session_state["meta"] = read_data(meta_file)
            except Exception as e:
                st.error(f"Failed to read metadata: {e}")

        st.write("Is the data log2 transformed?")
        col_a, col_b = st.columns(2)

        if col_a.button("Yes", key="log2_yes"):
            if "data" in st.session_state:
                st.session_state["org_data"] = inverse_log2_transform_data(st.session_state["data"], st.session_state.get("meta", pd.DataFrame()))
                st.session_state["log2_data"] = st.session_state["data"]
            else:
                st.error("No data loaded yet.")

        if col_b.button("No", key="log2_no"):
            if "data" in st.session_state:
                st.session_state["org_data"] = st.session_state["data"]
                st.session_state["log2_data"] = log2_transform_data(st.session_state["data"], st.session_state.get("meta", pd.DataFrame()))
            else:
                st.error("No data loaded yet.")

        st.markdown("---")
        st.number_input("Filter: At least", value=3, min_value=1, key="filter_num")
        st.selectbox("Value(s)", ["per group", "in at least one group"], key="filterop1")
        st.button("Apply Filter", key="apply_filter")

        st.markdown("---")
        st.header("Phospho Data Annotation")
        meta_file2 = st.file_uploader(
            "Upload Metadata (Phospho)",
            type=["csv", "xlsx", "txt", "tsv"],
            key="upload_meta2"
        )
        if meta_file2 is not None:
            try:
                st.session_state["meta2"] = read_data(meta_file2)
            except Exception as e:
                st.error(f"Failed to read phospho metadata: {e}")

        st.write("Is the data log2 transformed?")
        col_c, col_d = st.columns(2)

        if col_c.button("Yes", key="log2_yes2"):
            if "data3" in st.session_state:
                st.session_state["org_data3"] = inverse_log2_transform_data(st.session_state["data3"], st.session_state.get("meta2", pd.DataFrame()))
                st.session_state["log2_data3"] = st.session_state["data3"]
            else:
                st.error("No phospho data loaded yet.")

        if col_d.button("No", key="log2_no2"):
            if "data3" in st.session_state:
                st.session_state["org_data3"] = st.session_state["data3"]
                st.session_state["log2_data3"] = log2_transform_data(st.session_state["data3"], st.session_state.get("meta2", pd.DataFrame()))
            else:
                st.error("No phospho data loaded yet.")

        st.markdown("---")
        st.number_input("Filter: At least", value=3, min_value=1, key="filter_num2")
        st.selectbox("Value(s)", ["per group", "in at least one group"], key="filterop2")
        st.button("Apply Filter", key="apply_filter2")

        st.markdown("---")
        st.header("Color Scheme")
        st.selectbox(
            "Choose a Color Palette:",
            ["Default", "Mario Document Input", "Default16", "Warm/Cold", "Black/Grey", "Yue7"],
            key="color_palette"
        )
        st.button("Reload All Plots", key="reloadButton")

    with col2:
        st.subheader("Normal Data Condition Setup")
        st.number_input("Number of Conditions:", value=1, min_value=1, key="num_conditions")
        st.empty()
        st.markdown("---")
        st.subheader("Annotated Normal Data")
        st.dataframe(st.session_state.get("log2_data", pd.DataFrame()), key="displayed_data")
        st.markdown("---")
        st.subheader("Phospho Data Condition Setup")
        st.number_input("Number of Conditions:", value=1, min_value=1, key="num_conditions2")
        st.empty()
        st.markdown("---")
        st.subheader("Annotated Phospho Data")
        st.dataframe(st.session_state.get("log2_data3", pd.DataFrame()), key="displayed_data2")


def impute_data_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.button("Impute Data Values", key="ImputeEVE")
        st.markdown("---")
        st.header("Imputation Settings")
        st.number_input("q-Value:", value=0.01, key="qn1")
        st.number_input("Adjust Standard Deviation:", value=1, key="adj_stdn1")
        st.number_input("Random Seed:", value=1337, key="seedn1")
        st.selectbox("Level:", ["Protein", "Phosphosite"], key="leveln1")
        st.markdown("---")
        st.download_button("Download Imputed Data", data="", file_name="imputed.csv", key="imputed_data_down")
    with col2:
        st.subheader("Imputation Plots")
        st.empty()
        st.subheader("Imputed Data Table")
        st.dataframe([], key="imputed_data_tab")


def distribution_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.selectbox("Level:", ["Protein", "Phosphosite"], key="level4.5")
    with col2:
        st.subheader("QQ Norm Plot")
        st.empty()


def verification_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.selectbox("Level:", ["Protein", "Phosphosite"], key="level3.5")
    with col2:
        st.subheader("Verification Plots")
        st.empty()