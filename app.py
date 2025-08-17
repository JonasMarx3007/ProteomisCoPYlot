import streamlit as st
import sys
import os
from ui.i0_main_tabs import render_main_tabs

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def main():
    favicon_path = resource_path("CopilotIcoV1.ico")
    st.set_page_config(page_title="Proteomics Copilot", layout="wide", page_icon=favicon_path)
    st.title("Proteomics Copilot")
    render_main_tabs()

if __name__ == "__main__":
    main()