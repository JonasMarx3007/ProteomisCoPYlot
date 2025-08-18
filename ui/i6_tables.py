import streamlit as st
import pandas as pd

#MAIN
def tables_ui():
    tables_tabs = st.tabs([
        "Tables",
    ])

    with tables_tabs[0]:
        meta_ui()


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
