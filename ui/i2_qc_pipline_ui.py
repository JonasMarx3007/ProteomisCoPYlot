import streamlit as st
from utils.functions import *
import io

#MAIN
def qc_pipeline_ui():
    qc_tabs = st.tabs([
        "Coverage Plot",
        "Missing Value Plot",
        "Histogram Intensity Plot",
        "Boxplot Intensity Plot",
        "Cov Plot",
        "Principal Component Analysis",
        "Abundance Plot",
        "Correlation Plot"
    ])

    with qc_tabs[0]:
        coverage_plot_ui()
    with qc_tabs[1]:
        missing_value_plot_ui()
    with qc_tabs[2]:
        histogram_intensity_ui()
    with qc_tabs[3]:
        boxplot_intensity_ui()
    with qc_tabs[4]:
        cov_plot_ui()
    with qc_tabs[5]:
        principal_component_analysis_ui()
    with qc_tabs[6]:
        abundance_plot_ui()
    with qc_tabs[7]:
        correlation_plot_ui()


#SUB
def coverage_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        st.checkbox("Toggle IDs", value=False, key="toggle_id3")
        st.checkbox("Toggle Header", value=True, key="toggle_header3")
        st.checkbox("Toggle Legend", value=True, key="toggle_legend3")

        st.markdown("---")
        st.selectbox("Level:", ["Protein", "Peptide", "Phosphosite"], key="level3")
        st.selectbox("Type:", ["Normal", "Summary"], key="type3")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="plotWidth3")
        st.number_input("Height (cm):", value=10, key="plotHeight3")
        st.number_input("DPI:", value=300, key="plotDPI3")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat3")
        st.download_button("Download Plot", data="", file_name="coverage_plot.png", key="downloadCoveragePlot")
        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition3")
        st.text_area("Annotation Text:", key="text3")
        st.button("Add", key="addText3")
        st.button("Delete", key="deleteText3")

    with col2:
        if "data" in st.session_state and "meta" in st.session_state:
            try:
                fig = coverage_plot(
                    data=st.session_state["data"],
                    meta=st.session_state["meta"],
                    id=st.session_state.get("toggle_id3", True),
                    header=st.session_state.get("toggle_header3", True),
                    legend=st.session_state.get("toggle_legend3", True),
                    width=st.session_state.get("plotWidth3", 20),
                    height=st.session_state.get("plotHeight3", 10),
                    dpi=st.session_state.get("plotDPI3", 300)
                )
                st.pyplot(fig)
            except Exception as e:
                st.error(f"Error generating coverage plot: {e}")
        else:
            st.info("Please upload both data and metadata to see the plot.")


def missing_value_plot_ui():
    col1, col2 = st.columns([1, 3])

    with col1:
        header_toggle = st.checkbox("Toggle Header", value=True, key="toggle_header4")
        st.markdown("---")
        level = st.selectbox("Level:", ["Protein", "Peptide", "Phosphosite", "Precursor"], key="level4")
        bin_val = st.number_input("Bin missing values (optional):", value=0, min_value=0, key="missValBin4")
        st.markdown("---")
        text_toggle = st.checkbox("Toggle Text", value=True, key="toggle_text4")
        text_size = st.number_input("Text Size:", value=3.88, key="text_size4")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width = st.number_input("Width (cm):", value=20, key="plotWidth4")
        height = st.number_input("Height (cm):", value=10, key="plotHeight4")
        dpi = st.number_input("DPI:", value=300, key="plotDPI4")
        file_format = st.selectbox("Download Format:", ["png", "jpg", "svg", "pdf"], key="missValPlotFormat")
        st.download_button("Download Plot", data="", file_name=f"missval_plot.{file_format}", key="downloadMissValPlot")
        st.markdown("---")
        position = st.selectbox("Position:", ["Above", "Below"], key="textPosition4")
        annotation_text = st.text_area("Annotation Text:", key="text4")
        st.button("Add", key="addText4")
        st.button("Delete", key="deleteText4")

    with col2:
        if "data" in st.session_state and "meta" in st.session_state:
            try:
                fig = missing_value_plot(
                    data=st.session_state["data"],
                    meta=st.session_state["meta"],
                    bin=bin_val,
                    header=header_toggle,
                    text=text_toggle,
                    text_size=text_size,
                    width=width / 2.54,
                    height=height / 2.54,
                    dpi=dpi
                )
                st.pyplot(fig)
            except Exception as e:
                st.error(f"Error generating missing value plot: {e}")
        else:
            st.info("Please upload both data and metadata to see the plot.")


def histogram_intensity_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        header = st.checkbox("Show Header", value=True, key="toggle_header5")
        legend = st.checkbox("Show Legend", value=True, key="toggle_legend5")
        st.markdown("---")
        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="level5")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width_cm = st.number_input("Width (cm):", value=20, key="plotWidth5")
        height_cm = st.number_input("Height (cm):", value=10, key="plotHeight5")
        dpi = st.number_input("DPI:", value=300, key="plotDPI5")
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat5")

        st.download_button("Download Plot", data=b"", file_name=f"hist_intensity.{file_format}",
                           key="downloadHistIntPlot")

        st.markdown("---")
        position = st.selectbox("Position:", ["Above", "Below"], key="textPosition5")
        annotation_text = st.text_area("Annotation Text:", key="text5")
        st.button("Add", key="addText5")
        st.button("Delete", key="deleteText5")

    with col2:
        if "log2_data" in st.session_state and "meta" in st.session_state:
            try:
                data = st.session_state["log2_data"]
                meta = st.session_state["meta"]
                figsize = (width_cm / 2.54, height_cm / 2.54)

                fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
                histo_int(data, meta, header=header, legend=legend, ax=ax)
                st.pyplot(fig)
                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating histogram intensity plot: {e}")
        else:
            st.info("Please define both log2_data and metadata to see the plot.")


def boxplot_intensity_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        outliers = st.checkbox("Show Outliers", value=False)
        header = st.checkbox("Show Header", value=True)
        legend = st.checkbox("Show Legend", value=True)
        width_cm = st.number_input("Width (cm):", value=20)
        height_cm = st.number_input("Height (cm):", value=10)
        dpi = st.number_input("DPI:", value=300)
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"])
        st.download_button("Download Plot", data=b"", file_name=f"boxplot_intensity.{file_format}")

        st.markdown("---")
        position = st.selectbox("Position:", ["Above", "Below"])
        annotation_text = st.text_area("Annotation Text:")
        st.button("Add")
        st.button("Delete")


    with col2:
        if "log2_data" in st.session_state and "meta" in st.session_state:
            try:
                fig = boxplot_int(
                    data=st.session_state["log2_data"],
                    meta=st.session_state["meta"],
                    outliers=outliers,
                    header=header,
                    legend=legend,
                    width_cm=width_cm,
                    height_cm=height_cm,
                    dpi=dpi,
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("Please load 'log2_data' and 'meta' in session state.")


def cov_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        outliers = st.checkbox("Show Outliers", value=False, key="toggle_outliersCO")
        header = st.checkbox("Show Header", value=True, key="toggle_header7")
        legend = st.checkbox("Show Legend", value=True, key="toggle_legend7")

        st.markdown("---")
        st.selectbox("Level:", ["Protein", "Phosphosite"], key="level7")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        width_cm = st.number_input("Width (cm):", value=20, key="plotWidth7")
        height_cm = st.number_input("Height (cm):", value=10, key="plotHeight7")
        dpi = st.number_input("DPI:", value=300, key="plotDPI7")
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat7")

        st.download_button(
            "Download Plot",
            data=b"",
            file_name=f"cov_plot.{file_format}",
            key="downloadCovPlot"
        )

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition7")
        st.text_area("Annotation Text:", key="text7")
        st.button("Add", key="addText7")
        st.button("Delete", key="deleteText7")

    with col2:
        if "org_data" in st.session_state and "meta" in st.session_state:
            try:
                fig = cov_plot(
                    data=st.session_state["org_data"],
                    meta=st.session_state["meta"],
                    outliers=outliers,
                    header=header,
                    legend=legend,
                    width_cm=width_cm,
                    height_cm=height_cm,
                    dpi=dpi
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating Coefficient of Variation plot: {e}")
        else:
            st.info("Please load 'log2_data' and 'meta' in session state.")


def principal_component_analysis_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        legend = st.checkbox("Show Legend", value=True, key="toggle_legend8")
        header = st.checkbox("Show Header", value=True, key="toggle_header8")
        st.markdown("---")
        st.selectbox("Level:", ["Protein", "Phosphosite"], key="level8")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width_cm = st.number_input("Width (cm):", value=20, key="plotWidth8")
        height_cm = st.number_input("Height (cm):", value=10, key="plotHeight8")
        dpi = st.number_input("DPI:", value=300, key="plotDPI8")
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat8")
        st.download_button("Download Plot", data="", file_name=f"pca_plot.{file_format}", key="downloadPCAPlot")
        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition8")
        st.text_area("Annotation Text:", key="text8")
        st.button("Add", key="addText8")
        st.button("Delete", key="deleteText8")

    with col2:
        if "org_data" in st.session_state and "meta" in st.session_state:
            try:
                fig = pca_plot(
                    data=st.session_state["org_data"],
                    meta=st.session_state["meta"],
                    header=header,
                    legend=legend,
                    width_cm=width_cm,
                    height_cm=height_cm,
                    dpi=dpi
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating PCA plot: {e}")
        else:
            st.info("Please load 'log2_data' and 'meta' in session state.")


def abundance_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        legend = st.checkbox("Show Legend", value=True, key="legend9")
        header = st.checkbox("Show Header", value=True, key="header9")
        st.markdown("---")
        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="level9")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width_cm = st.number_input("Width (cm):", value=20, key="plotWidth9")
        height_cm = st.number_input("Height (cm):", value=10, key="plotHeight9")
        dpi = st.number_input("DPI:", value=300, key="plotDPI9")
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat9")
        st.download_button("Download Plot", data=b"", file_name=f"abundance_plot.{file_format}", key="downloadAbPlot9")
        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition9")
        st.text_area("Enter text:", key="text9", height=150)
        st.button("Add", key="addText9")
        st.button("Delete", key="deleteText9")

    with col2:
        condition = st.selectbox("Choose Condition:", ["All Conditions"], key="condition9")
        protein_choices = st.session_state.get("protein_choices", [])
        selected_proteins = st.multiselect("Select Proteins:", options=protein_choices, key="protein9")

        if "org_data" in st.session_state and "meta" in st.session_state:
            try:
                fig = abundance_plot(
                    data=st.session_state["org_data"],
                    meta=st.session_state["meta"],
                    workflow=level,
                    width_cm=width_cm,
                    height_cm=height_cm,
                    dpi=dpi,
                    legend=legend,
                    header=header
                )
                st.pyplot(fig)
            except Exception as e:
                st.error(f"Error generating abundance plot: {e}")
        else:
            st.info("Please load 'org_data' and 'meta' in session state.")


def correlation_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        if st.button("Change Display", key="Change12"):
            st.session_state["change_display12"] = not st.session_state.get("change_display12", False)
        if st.button("Toggle ID", key="toggle_id12"):
            st.session_state["toggle_id12"] = not st.session_state.get("toggle_id12", False)

        st.markdown("---")

        level = st.selectbox("Level:", ["Protein", "Phosphosite"], key="level12")

        st.markdown("---")

        text_position = st.selectbox(
            "Position:", options={"Above": "up", "Below": "down"}, index=0, key="textPosition12"
        )
        annotation_text = st.text_area(
            "Annotation Text:", value="", height=100, key="text12"
        )

        st.button("Add", key="addText12")
        st.button("Delete", key="deleteText12")

    with col2:
        st.header("Correlation Plot (Pearson)")
        plot_placeholder = st.empty()
