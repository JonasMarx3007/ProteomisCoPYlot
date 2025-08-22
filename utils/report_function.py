import io
import base64
from datetime import datetime
from utils.functions import *


import streamlit as st
import io
import base64
from datetime import datetime

def report_function():
    title = st.session_state.get("report_title", "").strip() or "Untitled Report"
    author = st.session_state.get("report_author", "").strip() or "Unknown Author"
    date_str = datetime.now().strftime("%Y-%m-%d")
    version = "Proteomics CoPYlot V0.01"

    html_parts = [
        f"<html><head><title>{title}</title></head><body>",
        f"<h1 style='font-size:2.5em'>{title}</h1>",
        f"<p>Author: {author}</p>",
        f"<p>Date: {date_str}</p>",
        f"<p>{version}</p>"
    ]

    #Prot - Coverage
    if "data" in st.session_state and "meta" in st.session_state:
        fig_3 = coverage_plot(
            data=st.session_state.get("data"),
            meta=st.session_state.get("meta"),
            id=st.session_state.get("id3", False),
            header=st.session_state.get("header3", True),
            legend=st.session_state.get("legend3", True),
            width=st.session_state.get("plotWidth3", 20),
            height=st.session_state.get("plotHeight3", 10),
            dpi=st.session_state.get("plotDPI3", 300),
            plot_colors=st.session_state.get("selected_colors", None)
        )
        fig_3.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig_3.savefig(buf, format="png", dpi=st.session_state.get("plotDPI3", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append("<h1>1 Protein Level</h1>")
        html_parts.append("<h2>1.1 Coverage Plot</h2>")
        html_parts.append("<p>This is some filler text describing the coverage plot.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - Missing Value
    if "data" in st.session_state and "meta" in st.session_state:
        fig4 = missing_value_plot(
            data=st.session_state.get("data"),
            meta=st.session_state.get("meta"),
            bin=st.session_state.get("bins4", 0),
            header=st.session_state.get("header4", True),
            text=st.session_state.get("plotText4", True),
            text_size=st.session_state.get("text_size4", 8),
            width=st.session_state.get("plotWidth4", 20) / 2.54,
            height=st.session_state.get("plotHeight4", 10) / 2.54,
            dpi=st.session_state.get("plotDPI4", 300)
        )
        fig4.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig4.savefig(buf, format="png", dpi=st.session_state.get("plotDPI4", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>1.2 Missing Value Plot</h2>")
        html_parts.append("<p>This plot shows the distribution of missing values.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - HistoInt
    if "log2_data" in st.session_state and "meta" in st.session_state:

        figsize5 = (st.session_state.get("plotWidth5", 20) / 2.54, st.session_state.get("plotHeight5", 10) / 2.54)
        fig, ax = plt.subplots(figsize=figsize5, dpi=st.session_state.get("plotDPI5", 300))
        fig5 = histo_int(st.session_state.get("log2_data"), st.session_state.get("meta"), header=st.session_state.get("header5", True),
                         legend=st.session_state.get("legend5", True), ax=ax, plot_colors=st.session_state.get("selected_colors", None))

        fig5.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig5.savefig(buf, format="png", dpi=st.session_state.get("plotDPI5", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append("<h2>1.3 Histogram Intensity Plot</h2>")
        html_parts.append("<p>Something about Histogram Intensity.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    html_parts.append("</body></html>")
    return "\n".join(html_parts)