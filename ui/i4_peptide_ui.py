import streamlit as st
from utils.functions import rt_vs_pred_rt_plot, modification_plot, missed_cleavage_plot, calculate_coverage, vis_coverage, fasta_to_dataframe
import io

#MAIN
def peptide_level_ui():
    peptide_tabs = st.tabs([
        "RT Plot",
        "Modification Plot",
        "Missed Cleavage Plot",
        "Sequence Coverage"
    ])

    with peptide_tabs[0]:
        rt_plot_ui()
    with peptide_tabs[1]:
        modification_plot_ui()
    with peptide_tabs[2]:
        missed_cleavage_plot_ui()
    with peptide_tabs[3]:
        sequence_coverage_ui()

#SUB
def rt_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        add_line = st.checkbox("Add Line", key="line14")
        header = st.checkbox("Toggle Header", key="header14")
        st.markdown("---")

        type = st.selectbox("Type:", ["Scatter Plot", "Hexbin Plot", "Density Plot"], key="type14")
        bins = st.number_input("Bins:", value=1000, key="bins14")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plotWidth = st.number_input("Width:", value=8, key="plotWidth14")
        plotHeight = st.number_input("Height:", value=6, key="plotHeight14")
        plotDPI = st.number_input("DPI:", value=100, key="plotDPI14")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat14")

        download_placeholder = st.empty()

    with col2:
        if "data2" in st.session_state and st.session_state["data2"] is not None:
            fig = rt_vs_pred_rt_plot(
                st.session_state["data2"],
                method=type,
                bins=bins,
                add_line=add_line,
                header=header,
                width=plotWidth,
                height=plotHeight,
                dpi=plotDPI
            )
            st.pyplot(fig)

            buf = io.BytesIO()
            fig.savefig(buf, format=plotFormat, dpi=plotDPI)
            buf.seek(0)

            download_placeholder.download_button(
                "Download Plot",
                data=buf,
                file_name=f"rt_plot.{plotFormat}",
                mime=f"image/{plotFormat}"
            )
        else:
            st.info("No data available for RT Plot.")


def modification_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        id = st.checkbox("Toggle ID", value=False, key="id15")
        header = st.checkbox("Toggle Header", key="header15")
        legend = st.checkbox("Toggle Legend", value=True, key="legend15")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plotWidth = st.number_input("Width:", value=10, key="plotWidth15")
        plotHeight = st.number_input("Height:", value=6, key="plotHeight15")
        plotDPI = st.number_input("DPI:", value=100, key="plotDPI15")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat15")

        download_placeholder = st.empty()

    with col2:
        if ("data2" in st.session_state and st.session_state["data2"] is not None) and \
                ("meta" in st.session_state and st.session_state["meta"] is not None):

            fig = modification_plot(
                st.session_state["data2"],
                st.session_state["meta"],
                id=id,
                header=header,
                legend=legend,
                width=plotWidth,
                height=plotHeight,
                dpi=plotDPI
            )
            st.pyplot(fig)

            buf = io.BytesIO()
            fig.savefig(buf, format=plotFormat, dpi=plotDPI)
            buf.seek(0)

            download_placeholder.download_button(
                "Download Plot",
                data=buf,
                file_name=f"modification_plot.{plotFormat}",
                mime=f"image/{plotFormat}"
            )
        else:
            st.info("Both data and meta information are required for the Modification Plot.")


def missed_cleavage_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        header = st.checkbox("Toggle Header", key="header16")
        id = st.checkbox("Toggle ID", value=False, key="id16")
        plotText = st.checkbox("Toggle Text", key="plotText16")
        text_size = st.number_input("Text Size:", value=8, key="text_size16")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plotWidth = st.number_input("Width:", value=10, key="plotWidth16")
        plotHeight = st.number_input("Height:", value=6, key="plotHeight16")
        plotDPI= st.number_input("DPI:", value=100, key="plotDPI16")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat16")

        download_placeholder = st.empty()

    with col2:
        if ("data2" in st.session_state and st.session_state["data2"] is not None) and \
                ("meta" in st.session_state and st.session_state["meta"] is not None):

            fig = missed_cleavage_plot(
                st.session_state["data2"],
                st.session_state["meta"],
                text=plotText,
                text_size=text_size,
                header=header,
                id=id,
                width=plotWidth,
                height=plotHeight,
                dpi=plotDPI
            )
            st.pyplot(fig)

            buf = io.BytesIO()
            fig.savefig(buf, format=plotFormat, dpi=plotDPI)
            buf.seek(0)

            download_placeholder.download_button(
                "Download Plot",
                data=buf,
                file_name=f"missed_cleavage_plot.{plotFormat}",
                mime=f"image/{plotFormat}"
            )
        else:
            st.info("Both data and meta information are required for the Missed Cleavage Plot.")


def sequence_coverage_ui():
    if "data2" not in st.session_state:
        st.warning("No data2 found in session state")
        return

    data2 = st.session_state["data2"]
    left_col, right_col = st.columns([1, 2])

    with left_col:
        species = st.selectbox("Select database", ["Mouse", "Human"])
        if st.button("Load Database"):
            if species == "Mouse":
                st.session_state["db"] = fasta_to_dataframe("data/db/UP000000589_10090.fasta")
            else:
                st.session_state["db"] = fasta_to_dataframe("data/db/UP000005640_9606.fasta")

        protein = st.selectbox("Select a protein", sorted(data2["ProteinNames"].unique()))
        chunk_size = st.number_input("Chunk size", min_value=10, max_value=500, value=100, step=10)

    if "db" not in st.session_state or st.session_state["db"] is None:
        with right_col:
            st.info("Load a database first")
        return

    db = st.session_state["db"]
    coverage = calculate_coverage(data2, db, protein)

    with right_col:
        st.write(f"Coverage for {protein}: {coverage}%")
        coverage_text = vis_coverage(data2, db, protein, chunk_size=chunk_size)
        st.text(coverage_text)