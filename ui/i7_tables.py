import streamlit as st
import pandas as pd
import io
import csv
from utils.functions import bool_to_str, number_to_str

#MAIN
def tables_ui():
    tables_tabs = st.tabs([
        "Meta",
        "Log",
        "Auto Log"
    ])

    with tables_tabs[0]:
        meta_ui()
    with tables_tabs[1]:
        log_ui()
    with tables_tabs[2]:
        log_ui_auto()

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

    toggle_id3 = bool_to_str(st.session_state.get("toggle_id3", None))
    toggle_header3 = bool_to_str(st.session_state.get("toggle_header3", None))
    toggle_legend3 = bool_to_str(st.session_state.get("toggle_legend3", None))
    toggle_header4 = bool_to_str(st.session_state.get("toggle_header4", None))
    bin_val4 = number_to_str(st.session_state.get("missValBin4", None))
    toggle_text4 = bool_to_str(st.session_state.get("toggle_text4", None))
    text_size4 = number_to_str(st.session_state.get("text_size4", None))
    toggle_header5 = bool_to_str(st.session_state.get("toggle_header5", None))
    toggle_legend5 = bool_to_str(st.session_state.get("toggle_legend5", None))
    toggle_legend6 = bool_to_str(st.session_state.get("toggle_legend6", None))
    toggle_id6 = bool_to_str(st.session_state.get("show_id6", None))
    outliers6 = bool_to_str(st.session_state.get("outliers6", None))
    mode6 = st.session_state.get("mode6", None)
    toggle_header6 = bool_to_str(st.session_state.get("toggle_header6", None))
    outliers7 = bool_to_str(st.session_state.get("outliers7", None))
    toggle_header7 = bool_to_str(st.session_state.get("header7", None))
    toggle_legend7 = bool_to_str(st.session_state.get("legend7", None))
    toggle_header8 = bool_to_str(st.session_state.get("header8", None))
    toggle_legend8 = bool_to_str(st.session_state.get("legend8", None))
    level12 = st.session_state.get("method12", None)
    toggle_id12 = bool_to_str(st.session_state.get("toggleId12", None))
    style14 = st.session_state.get("style14", None)
    line14 = bool_to_str(st.session_state.get("line14", None))
    bins14 = number_to_str(st.session_state.get("bins14", None))
    toggle_header14 = bool_to_str(st.session_state.get("header14", None))
    toggle_id15 = bool_to_str(st.session_state.get("toggle_id15", None))
    toggle_header15 = bool_to_str(st.session_state.get("header15", None))
    toggle_legend15 = bool_to_str(st.session_state.get("legend15", None))
    toggle_id16 = bool_to_str(st.session_state.get("id_toggle16", None))
    toggle_text16 = bool_to_str(st.session_state.get("toggle_text16", None))
    text_size16 = number_to_str(st.session_state.get("text_size16", None))
    toggle_header16 = bool_to_str(st.session_state.get("header16", None))

    log_df = pd.DataFrame({
        "Var": ["CoveragePlotID", "CoveragePlotHeader", "CoveragePlotLegend",
                "MissingValuePlotHeader", "MissingValuePlotBin", "MissingValueText", "MissingValueTextSize",
                "HistogramIntHeader", "HistogramIntLegend",
                "BoxplotIntLegend", "BoxplotIntID", "BoxplotIntOut","BoxplotIntMean", "BoxplotIntHeader",
                "CovPlotOut", "CovPlotHeader", "CovPlotLegend",
                "PCAHeader", "PCALegend",
                "CorrPlotDisplay", "CorrPlotID",
                "RTPlotStyle", "RTLine", "HexbinsRT", "RTHeader",
                "ModPlotID", "ModPlotHeader", "ModPlotLegend",
                "MissedCleavID", "MissedCleavText", "MissedCleavTextSize", "MissedCleavHeader"],
        "Select": [toggle_id3, toggle_header3, toggle_legend3,
                   toggle_header4, bin_val4, toggle_text4, text_size4,
                   toggle_header5, toggle_legend5,
                   toggle_legend6, toggle_id6, outliers6, mode6, toggle_header6,
                   outliers7, toggle_header7, toggle_legend7,
                   toggle_header8, toggle_legend8,
                   level12, toggle_id12,
                   style14, line14, bins14, toggle_header14,
                   toggle_id15, toggle_header15, toggle_legend15,
                   toggle_id16, toggle_text16, text_size16, toggle_header16]
    })

    st.dataframe(log_df.style.hide(axis="index"), use_container_width=True)

    buf = io.StringIO()
    log_df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)
    buf.seek(0)

    st.download_button(
        "Download Log CSV",
        data=buf.getvalue().encode("utf-8"),
        file_name="system_variables_log.csv",
        mime="text/csv"
    )


def log_ui_auto():
    st.header("System Variables Log")

    def format_value(val):
        if isinstance(val, bool) or val is None:
            return bool_to_str(val)
        elif isinstance(val, (int, float)):
            return number_to_str(val)
        else:
            return str(val)

    # Build DataFrame dynamically
    log_dict = {"Var": [], "Select": []}
    for key in st.session_state:
        log_dict["Var"].append(key)
        log_dict["Select"].append(format_value(st.session_state[key]))

    log_df = pd.DataFrame(log_dict)

    st.dataframe(log_df.style.hide(axis="index"), use_container_width=True)

    buf = io.StringIO()
    log_df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)
    buf.seek(0)

    st.download_button(
        "Download Log CSV",
        data=buf.getvalue().encode("utf-8"),
        file_name="system_variables_log.csv",
        mime="text/csv"
    )