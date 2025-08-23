import streamlit as st
import pandas as pd
import io
import csv
from utils.functions import bool_to_str, number_to_str, extract_number, format_value
from utils.report_function import report_function
import streamlit.components.v1 as components

#MAIN
def summary_ui():
    summary_tabs = st.tabs([
        "Tables",
        "Log",
        "Report"
    ])

    with summary_tabs[0]:
        meta_ui()
    with summary_tabs[1]:
        log_ui()
    with summary_tabs[2]:
        report_ui()

#SUB
def meta_ui():
    if "meta" in st.session_state and isinstance(st.session_state["meta"], pd.DataFrame):
        st.subheader("Normal Data Metadata")
        st.dataframe(st.session_state["meta"], use_container_width=True)

    if "meta2" in st.session_state and isinstance(st.session_state["meta2"], pd.DataFrame):
        st.subheader("Phospho Data Metadata")
        st.dataframe(st.session_state["meta2"], use_container_width=True)

    if ("meta" not in st.session_state) and ("meta2" not in st.session_state):
        st.info("No metadata available yet. Please upload or generate metadata.")


def log_ui():
    st.header("System Variables Log")
    log_dict = {"Var": [], "Select": []}
    exclude_substrings = [
        "level", "plotFormat", "plotDPI", "data", "meta", "paired", "pval", "cond", "fc", "download",
        "Impute", "collapse", "n1", "corrected", "add", "delete", "Width", "Height", "DPI", "apply_filter",
        "log2", "term", "10.5", "11", "feature", "textPositions", "title", "author", "filter_num", "filterop",
        "selected_colors"
    ]
    for key in st.session_state:
        if not any(sub in key for sub in exclude_substrings):
            log_dict["Var"].append(key)
            log_dict["Select"].append(format_value(st.session_state[key]))
    log_df = pd.DataFrame(log_dict)
    file_df = log_df[log_df["Var"].str.contains("file", case=False)]
    other_df = log_df[~log_df["Var"].str.contains("file", case=False)]
    file_df = file_df.iloc[sorted(range(len(file_df)), key=lambda i: extract_number(file_df.iloc[i]["Var"]))]
    other_df = other_df.iloc[sorted(range(len(other_df)), key=lambda i: extract_number(other_df.iloc[i]["Var"]))]
    log_df = pd.concat([file_df, other_df], ignore_index=True)
    st.dataframe(log_df.style.hide(axis="index"), use_container_width=True)
    buf = io.StringIO()
    log_df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)
    buf.seek(0)
    st.download_button(
        "Download Log CSV",
        data=buf.getvalue().encode("utf-8"),
        file_name="system_variables_log.csv",
        mime="text/csv",
        key="download_log_csv"
    )

def report_ui():
    st.header("Report Generator")

    title_input = st.text_input("Report Title:", key="report_title")
    author_input = st.text_input("Author:", key="report_author")

    if st.button("Generate HTML Report"):
        html_content = report_function()

        st.download_button(
            label="Download HTML Report",
            data=html_content,
            file_name="CoPYlot_report.html",
            mime="text/html"
        )

        components.html(html_content, height=800, scrolling=True)