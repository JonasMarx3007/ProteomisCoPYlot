import streamlit as st
import pandas as pd
import os
import io
from utils.functions import read_data, inverse_log2_transform_data, log2_transform_data, filter_data, qqnorm_plot, first_digit_distribution, data_pattern_structure, sanitize_dataframe, impute_values

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
            st.session_state["data"] = sanitize_dataframe(read_data(st.session_state["file"]))
            st.dataframe(st.session_state["data"])
        else:
            st.dataframe([], key="table1")

        st.markdown("---")
        st.subheader("Phospho Data Preview")
        if st.session_state.get("file3"):
            st.session_state["data3"] = sanitize_dataframe(read_data(st.session_state["file3"]))
            st.dataframe(st.session_state["data3"])
        else:
            st.dataframe([], key="table3")

        st.markdown("---")
        st.subheader("Full Report Preview")
        if st.session_state.get("file2"):
            st.session_state["data2"] = sanitize_dataframe(read_data(st.session_state["file2"]))
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

        st.header("Color Scheme")
        color_palette_name = st.selectbox(
            "Choose a Color Palette:",
            ["Default", "Yue7"],
            key="color_palette"
        )

        palette_dict = {
            "Default": None,
            "Yue7": ["#2D5F85", "#5184B2", "#AAD4F8", "#F2F5FA", "#F1A7B5", "#D55276", "#AB3A54"]
        }

        st.session_state["selected_colors"] = palette_dict.get(color_palette_name, palette_dict["Default"])

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

            filter_num = st.number_input("Filter: At least", value=3, min_value=0, key="filter_num")
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
        st.header("Imputation Settings")
        q_val = st.number_input("q-Value:", value=0.01, key="qn1")
        adj_std = st.number_input("Adjust Standard Deviation:", value=1, key="adj_stdn1")
        seed = st.number_input("Random Seed:", value=1337, key="seedn1")
        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="leveln1")
        st.markdown("---")
        impute_btn = st.button("Impute Data Values", key="ImputeEVE")
        st.markdown("---")
        download_placeholder = st.empty()
    with col2:
        st.subheader("Imputation Plots")
        plot1 = st.empty()
        plot2 = st.empty()
        plot3 = st.empty()
        st.subheader("Imputed Data Table")
        table_placeholder = st.empty()

    if "imputed_log2_data" not in st.session_state:
        st.session_state.imputed_log2_data = None
    if "imputed_log2_data3" not in st.session_state:
        st.session_state.imputed_log2_data3 = None

    if level == "Protein":
        base_data = st.session_state.get("filtered_log2_data", st.session_state.get("log2_data"))
        meta = st.session_state.get("meta")
        if base_data is not None and meta is not None:
            plot1.pyplot(impute_values(base_data, meta, ret=1, q=q_val, adj_std=adj_std, seed=seed))
            plot2.pyplot(impute_values(base_data, meta, ret=2, q=q_val, adj_std=adj_std, seed=seed))
            plot3.pyplot(impute_values(base_data, meta, ret=3, q=q_val, adj_std=adj_std, seed=seed))

    elif level == "Phosphosite":
        base_data = st.session_state.get("filtered_log2_data3", st.session_state.get("log2_data3"))
        meta = st.session_state.get("meta2")
        if base_data is not None and meta is not None:
            plot1.pyplot(impute_values(base_data, meta, ret=1, q=q_val, adj_std=adj_std, seed=seed))
            plot2.pyplot(impute_values(base_data, meta, ret=2, q=q_val, adj_std=adj_std, seed=seed))
            plot3.pyplot(impute_values(base_data, meta, ret=3, q=q_val, adj_std=adj_std, seed=seed))

    if impute_btn:
        if level == "Protein" and base_data is not None and meta is not None:
            st.session_state.imputed_log2_data = impute_values(base_data, meta, ret=0, q=q_val, adj_std=adj_std, seed=seed)

        elif level == "Phosphosite" and base_data is not None and meta is not None:
            st.session_state.imputed_log2_data3 = impute_values(base_data, meta, ret=0, q=q_val, adj_std=adj_std, seed=seed)

    if level == "Protein" and st.session_state.imputed_log2_data is not None:
        table_placeholder.dataframe(st.session_state.imputed_log2_data)
        csv_buf = io.StringIO()
        st.session_state.imputed_log2_data.to_csv(csv_buf, index=False)
        download_placeholder.download_button(
            "Download Imputed Data",
            data=csv_buf.getvalue(),
            file_name="imputed_log2_data.csv"
        )

    elif level == "Phosphosite" and st.session_state.imputed_log2_data3 is not None:
        table_placeholder.dataframe(st.session_state.imputed_log2_data3)
        csv_buf = io.StringIO()
        st.session_state.imputed_log2_data3.to_csv(csv_buf, index=False)
        download_placeholder.download_button(
            "Download Imputed Phospho Data",
            data=csv_buf.getvalue(),
            file_name="imputed_log2_data_phospho.csv"
        )


def distribution_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="level4.5")
    with col2:
        st.subheader("QQ Norm Plot")

        if level == "Protein":
            log_data_key = "log2_data"
        else:
            log_data_key = "log2_data3"

        if log_data_key in st.session_state and st.session_state[log_data_key] is not None:
            fig = qqnorm_plot(st.session_state[log_data_key], st.session_state["meta"])
            st.pyplot(fig)
        else:
            st.info(f"{log_data_key} is not defined yet.")

        assets_dir = os.path.join(os.path.dirname(__file__), "..", "assets")
        img1 = os.path.join(assets_dir, "qqnorm.jpg")
        img2 = os.path.join(assets_dir, "qqnorm_txt.jpg")

        st.image([img1, img2])


def verification_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="level3.5")
    with col2:
        st.subheader("Verification Plots")

        if level == "Protein":
            org_data_key = "org_data"
        else:
            org_data_key = "org_data3"

        if org_data_key in st.session_state and st.session_state[org_data_key] is not None:
            data = st.session_state[org_data_key]
            meta = st.session_state["meta"]

            fig1 = first_digit_distribution(data, meta)
            st.pyplot(fig1)

            fig2, _ = data_pattern_structure(data, meta)
            st.pyplot(fig2)
        else:
            st.info(f"{org_data_key} is not defined yet.")
