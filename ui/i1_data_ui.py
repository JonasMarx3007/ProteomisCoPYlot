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
        meta_file = st.file_uploader("Upload Metadata (Protein Group)", type=["csv", "xlsx", "txt", "tsv"], key="upload_meta")
        if meta_file:
            try:
                st.session_state["meta"] = read_data(meta_file)
            except Exception as e:
                st.error(f"Failed to read metadata: {e}")

        st.header("Phospho Data Annotation")
        meta_file2 = st.file_uploader("Upload Metadata (Phospho)", type=["csv", "xlsx", "txt", "tsv"], key="upload_meta2")
        if meta_file2:
            try:
                st.session_state["meta2"] = read_data(meta_file2)
            except Exception as e:
                st.error(f"Failed to read phospho metadata: {e}")

        st.markdown("---")
        st.header("Color Scheme")
        st.selectbox("Choose a Color Palette:", ["Default", "Mario Document Input", "Default16", "Warm/Cold", "Black/Grey", "Yue7"], key="color_palette")
        st.button("Reload All Plots", key="reloadButton")

    with col2:
        st.subheader("Normal Data Condition Setup")
        all_cols = list(st.session_state["data"].columns) if "data" in st.session_state else []
        num_conditions = st.number_input("Number of Conditions:", value=1, min_value=1, key="num_conditions")
        condition_meta_rows = []
        used_cols = []
        for i in range(int(num_conditions)):
            condition_name = st.text_input(f"Condition {i+1} Name", key=f"cond_name_{i}")
            remaining_cols = [c for c in all_cols if c not in used_cols]
            selected_cols = st.multiselect(f"Columns for {condition_name}", remaining_cols, key=f"cond_cols_{i}")
            used_cols += selected_cols
            for col in selected_cols:
                condition_meta_rows.append({"sample": col, "condition": condition_name})
        if st.button("Generate Meta for Normal Data"):
            if condition_meta_rows:
                st.session_state["meta"] = pd.DataFrame(condition_meta_rows)

        if "meta" in st.session_state:
            st.write("Is the data log2 transformed?")
            col_a, col_b = st.columns(2)
            if col_a.button("Yes", key="log2_yes"):
                if "data" in st.session_state:
                    st.session_state["org_data"] = inverse_log2_transform_data(st.session_state["data"], st.session_state["meta"])
                    st.session_state["log2_data"] = st.session_state["data"]
                else:
                    st.error("No data loaded yet.")
            if col_b.button("No", key="log2_no"):
                if "data" in st.session_state:
                    st.session_state["org_data"] = st.session_state["data"]
                    st.session_state["log2_data"] = log2_transform_data(st.session_state["data"], st.session_state["meta"])
                else:
                    st.error("No data loaded yet.")

            filter_num = st.number_input("Filter: At least", value=3, min_value=1, key="filter_num")
            filterop = st.selectbox("Value(s)", ["per group", "in at least one group"], key="filterop1")
            if st.button("Apply Filter", key="apply_filter"):
                if "log2_data" in st.session_state:
                    before_count = st.session_state["log2_data"].shape[0]
                    st.session_state["filtered_log2_data"] = filter_data(st.session_state["log2_data"], st.session_state["meta"], filter_num, filterop)
                    after_count = st.session_state["filtered_log2_data"].shape[0]
                    st.success(f"Filtered Normal Data: {before_count} → {after_count} proteins")
                else:
                    st.error("No log2 data available for filtering.")

            st.subheader("Annotated Normal Data (Log2)")
            st.dataframe(st.session_state.get("filtered_log2_data", st.session_state.get("log2_data", pd.DataFrame())), key="displayed_data")

        st.subheader("Phospho Data Condition Setup")
        all_cols2 = list(st.session_state["data3"].columns) if "data3" in st.session_state else []
        num_conditions2 = st.number_input("Number of Conditions:", value=1, min_value=1, key="num_conditions2")
        condition_meta_rows2 = []
        used_cols2 = []
        for i in range(int(num_conditions2)):
            condition_name2 = st.text_input(f"Phospho Condition {i+1} Name", key=f"cond_name2_{i}")
            remaining_cols2 = [c for c in all_cols2 if c not in used_cols2]
            selected_cols2 = st.multiselect(f"Columns for {condition_name2}", remaining_cols2, key=f"cond_cols2_{i}")
            used_cols2 += selected_cols2
            for col in selected_cols2:
                condition_meta_rows2.append({"sample": col, "condition": condition_name2})
        if st.button("Generate Meta for Phospho Data"):
            if condition_meta_rows2:
                st.session_state["meta2"] = pd.DataFrame(condition_meta_rows2)

        if "meta2" in st.session_state:
            st.write("Is the data log2 transformed?")
            col_c, col_d = st.columns(2)
            if col_c.button("Yes", key="log2_yes2"):
                if "data3" in st.session_state:
                    st.session_state["org_data3"] = inverse_log2_transform_data(st.session_state["data3"], st.session_state["meta2"])
                    st.session_state["log2_data3"] = st.session_state["data3"]
                else:
                    st.error("No phospho data loaded yet.")
            if col_d.button("No", key="log2_no2"):
                if "data3" in st.session_state:
                    st.session_state["org_data3"] = st.session_state["data3"]
                    st.session_state["log2_data3"] = log2_transform_data(st.session_state["data3"], st.session_state["meta2"])
                else:
                    st.error("No phospho data loaded yet.")

            filter_num2 = st.number_input("Filter: At least", value=3, min_value=1, key="filter_num2")
            filterop2 = st.selectbox("Value(s)", ["per group", "in at least one group"], key="filterop2")
            if st.button("Apply Filter", key="apply_filter2"):
                if "log2_data3" in st.session_state:
                    before_count = st.session_state["log2_data3"].shape[0]
                    st.session_state["filtered_log2_data3"] = filter_data(st.session_state["log2_data3"], st.session_state["meta2"], filter_num2, filterop2)
                    after_count = st.session_state["filtered_log2_data3"].shape[0]
                    st.success(f"Filtered Phospho Data: {before_count} → {after_count} proteins")
                else:
                    st.error("No log2 phospho data available for filtering.")

            st.subheader("Annotated Phospho Data (Log2)")
            st.dataframe(st.session_state.get("filtered_log2_data3", st.session_state.get("log2_data3", pd.DataFrame())), key="displayed_data2")


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