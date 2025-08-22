import streamlit as st
import io
import base64
from datetime import datetime
from utils.functions import *


def report_function():
    selected_colors = st.session_state.get("selected_colors", None)
    id3 = st.session_state.get("id3", False)
    header3 = st.session_state.get("header3", True)
    legend3 = st.session_state.get("legend3", True)
    type3 = st.session_state.get("type3", "Normal")
    plotWidth3 = st.session_state.get("plotWidth3", 20)
    plotHeight3 = st.session_state.get("plotHeight3", 10)
    plotDPI3 = st.session_state.get("plotDPI3", 300)

    data = st.session_state.get("data", None)
    meta = st.session_state.get("meta", None)

    if type3 == "Normal":
        fig3 = coverage_plot(
            data=data,
            meta=meta,
            id=id3,
            header=header3,
            legend=legend3,
            width=plotWidth3,
            height=plotHeight3,
            dpi=plotDPI3,
            plot_colors=selected_colors
        )
    elif type3 == "Summary":
        fig3 = coverage_plot_summary(
            data=data,
            meta=meta,
            id=id3,
            header=header3,
            legend=legend3,
            plot_colors=selected_colors,
            width=plotWidth3,
            height=plotHeight3,
            dpi=plotDPI3
        )

    fig3.set_size_inches(12, 6)

    buf = io.BytesIO()
    fig3.savefig(buf, format="png", dpi=plotDPI3, bbox_inches="tight")
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")

    title = st.session_state.get("report_title", "").strip()
    author = st.session_state.get("report_author", "").strip()

    if not title:
        title = "Untitled Report"
    if not author:
        author = "Unknown Author"

    date_str = datetime.now().strftime("%Y-%m-%d")
    version = "Proteomics CoPYlot V0.01"

    html_content = f"""
    <html>
        <head><title>{title}</title></head>
        <body>
            <h1 style="font-size:2.5em">{title}</h1>
            <p>Author: {author}<p>
            <p>Date: {date_str}</p>
            <p>{version}</p>

            <h1>1 Protein Level</h1>
            <h2>1.1 Coverage Plot</h2>
            <p>This is some filler text describing the plot.</p>
            <img src="data:image/png;base64,{img_base64}" alt="Coverage Plot" style="width:100%;height:auto;">
        </body>
    </html>
    """
    return html_content