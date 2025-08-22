import streamlit as st
import pandas as pd
import io
import csv
from utils.functions import bool_to_str, number_to_str

#MAIN
def tables_ui():
    tables_tabs = st.tabs([
        "Meta",
        "Log"
    ])

    with tables_tabs[0]:
        meta_ui()
    with tables_tabs[1]:
        log_ui()


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
    toggle_legend6 = bool_to_str(st.session_state.get("legend6", None))
    toggle_header6 = bool_to_str(st.session_state.get("header6", None))

    log_df = pd.DataFrame({
        "Var": ["CoveragePlotID", "CoveragePlotHeader", "CoveragePlotLegend",
                "MissingValuePlotHeader", "MissingValuePlotBin", "MissingValueText", "MissingValueTextSize",
                "HistogramIntHeader", "HistogramIntLegend",
                "BoxplotIntLegend", "BoxplotIntID", "BoxplotIntOut","BoxplotIntMean", "BoxplotIntHeader", "BoxplotIntHeader",
                "CovPlotOut", "CovPlotHeader", "CovPlotLegend",
                "PCAHeader", "PCALegend", "CorrPlotDisplay", "CorrPlotID",
                "RTPlotStyle", "RTLine", "HexbinsRT", "RTHeader",
                "ModPlotID", "ModPlotHeader",
                "MissedCleavID", "MissedCleavText", "MissedCleavTextSize", "MissedCleavHeader"],
        "Select": [toggle_id3, toggle_header3, toggle_legend3,
                   toggle_header4, bin_val4, toggle_text4, text_size4,
                   toggle_header5, toggle_legend5,
                   toggle_legend6, ]
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