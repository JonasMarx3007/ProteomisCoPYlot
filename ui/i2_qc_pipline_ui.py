import streamlit as st
import matplotlib.pyplot as plt
import io
from utils.functions import coverage_plot, missing_value_plot, histo_int, boxplot_int, cov_plot, pca_plot, abundance_plot, corr_plot, coverage_plot_pep, missing_value_plot_prec, missing_value_plot_pep

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

    available_levels = []
    if "data" in st.session_state:
        available_levels.append("Protein")
    if "data2" in st.session_state:
        available_levels.append("Peptide")
    if "data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        toggle_id = st.checkbox("Toggle IDs", value=False, key="toggle_id3")
        toggle_header = st.checkbox("Toggle Header", value=True, key="toggle_header3")
        toggle_legend = st.checkbox("Toggle Legend", value=True, key="toggle_legend3")

        st.markdown("---")

        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level3")
        else:
            st.info("No data available for plotting.")
            return

        plot_type = st.selectbox("Type:", ["Normal", "Summary"], key="type3")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        width = st.number_input("Width (cm):", value=20, key="plotWidth3")
        height = st.number_input("Height (cm):", value=10, key="plotHeight3")
        dpi = st.number_input("DPI:", value=300, key="plotDPI3")

        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat3")

        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition3")
        st.text_area("Annotation Text:", key="text3")
        st.button("Add", key="addText3")
        st.button("Delete", key="deleteText3")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("data", None)
            meta_to_use = st.session_state.get("meta", None)
            plot_func = coverage_plot
        elif level == "Peptide":
            data_to_use = st.session_state.get("data2", None)
            meta_to_use = st.session_state.get("meta", None)
            plot_func = coverage_plot_pep
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("data3", None)
            meta_to_use = st.session_state.get("meta2", None)
            plot_func = coverage_plot

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = plot_func(
                    data=data_to_use,
                    meta=meta_to_use,
                    id=toggle_id,
                    header=toggle_header,
                    legend=toggle_legend,
                    width=width,
                    height=height,
                    dpi=dpi
                )
                st.pyplot(fig)

                import io
                buf = io.BytesIO()
                fig.savefig(buf, format=file_format, dpi=dpi)
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"coverage_plot.{file_format}",
                    mime=f"image/{file_format}"
                )

            except Exception as e:
                st.error(f"Error generating coverage plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def missing_value_plot_ui():
    col1, col2 = st.columns([1, 3])

    available_levels = []
    if "data" in st.session_state:
        available_levels.append("Protein")
    if "data2" in st.session_state:
        available_levels.extend(["Peptide", "Precursor"])
    if "data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        header_toggle = st.checkbox("Toggle Header", value=True, key="toggle_header4")
        st.markdown("---")
        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level4")
        else:
            st.info("No data available for plotting.")
            return

        bin_val = st.number_input("Bin missing values (optional):", value=0, min_value=0, key="missValBin4")
        st.markdown("---")
        text_toggle = st.checkbox("Toggle Text", value=True, key="toggle_text4")
        text_size = st.number_input("Text Size:", value=8, key="text_size4")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width = st.number_input("Width (cm):", value=20, key="plotWidth4")
        height = st.number_input("Height (cm):", value=10, key="plotHeight4")
        dpi = st.number_input("DPI:", value=300, key="plotDPI4")
        file_format = st.selectbox("Download Format:", ["png", "jpg", "svg", "pdf"], key="missValPlotFormat")

        if "data" in st.session_state or "data2" in st.session_state or "data3" in st.session_state:
            try:
                if level == "Protein":
                    data_to_use = st.session_state.get("data", None)
                    meta_to_use = st.session_state.get("meta", None)
                    plot_func = missing_value_plot
                elif level == "Peptide":
                    data_to_use = st.session_state.get("data2", None)
                    meta_to_use = st.session_state.get("meta", None)
                    plot_func = missing_value_plot_pep
                elif level == "Precursor":
                    data_to_use = st.session_state.get("data2", None)
                    meta_to_use = st.session_state.get("meta", None)
                    plot_func = missing_value_plot_prec
                elif level == "Phosphosite":
                    data_to_use = st.session_state.get("data3", None)
                    meta_to_use = st.session_state.get("meta2", None)
                    plot_func = missing_value_plot

                fig = plot_func(
                    data=data_to_use,
                    meta=meta_to_use,
                    bin=bin_val,
                    header=header_toggle,
                    text=text_toggle,
                    text_size=text_size,
                    width=width / 2.54,
                    height=height / 2.54,
                    dpi=dpi
                )

                buf = io.BytesIO()
                fig.savefig(buf, format=file_format, dpi=dpi)
                buf.seek(0)
                st.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"missval_plot.{file_format}",
                    mime=f"image/{file_format}"
                )

            except Exception as e:
                st.error(f"Error generating missing value plot: {e}")

        st.markdown("---")
        position = st.selectbox("Position:", ["Above", "Below"], key="textPosition4")
        annotation_text = st.text_area("Annotation Text:", key="text4")
        st.button("Add", key="addText4")
        st.button("Delete", key="deleteText4")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("data", None)
            meta_to_use = st.session_state.get("meta", None)
            plot_func = missing_value_plot
        elif level == "Peptide":
            data_to_use = st.session_state.get("data2", None)
            meta_to_use = st.session_state.get("meta", None)
            plot_func = missing_value_plot_pep
        elif level == "Precursor":
            data_to_use = st.session_state.get("data2", None)
            meta_to_use = st.session_state.get("meta", None)
            plot_func = missing_value_plot_prec
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("data3", None)
            meta_to_use = st.session_state.get("meta2", None)
            plot_func = missing_value_plot

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = plot_func(
                    data=data_to_use,
                    meta=meta_to_use,
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
            st.info("Please upload the corresponding data and metadata to see the plot.")


def histogram_intensity_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "log2_data" in st.session_state:
        available_levels.append("Protein")
    if "log2_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        header = st.checkbox("Show Header", value=True, key="toggle_header5")
        legend = st.checkbox("Show Legend", value=True, key="toggle_legend5")
        st.markdown("---")
        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level5")
        else:
            st.info("No data available for plotting.")
            return
        st.markdown("---")
        st.header("Plot Size & Resolution")
        width_cm = st.number_input("Width (cm):", value=20, key="plotWidth5")
        height_cm = st.number_input("Height (cm):", value=10, key="plotHeight5")
        dpi = st.number_input("DPI:", value=300, key="plotDPI5")
        file_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat5")

        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition5")
        st.text_area("Annotation Text:", key="text5")
        st.button("Add", key="addText5")
        st.button("Delete", key="deleteText5")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("log2_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("log2_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                figsize = (width_cm / 2.54, height_cm / 2.54)
                fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
                histo_int(data_to_use, meta_to_use, header=header, legend=legend, ax=ax)
                st.pyplot(fig)

                import io
                buf = io.BytesIO()
                fig.savefig(buf, format=file_format, dpi=dpi)
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"hist_intensity.{file_format}",
                    mime=f"image/{file_format}"
                )

                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating histogram intensity plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def boxplot_intensity_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "log2_data" in st.session_state:
        available_levels.append("Protein")
    if "log2_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        if not available_levels:
            st.info("No data available for plotting.")
            return

        st.checkbox("Show Outliers", value=False, key="outliers6")
        st.checkbox("Show Header", value=True, key="header6")
        st.checkbox("Show Legend", value=True, key="legend6")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level6")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="width6")
        st.number_input("Height (cm):", value=10, key="height6")
        st.number_input("DPI:", value=300, key="dpi6")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="format6")
        st.download_button("Download Plot", data=b"", file_name=f"boxplot_intensity.png", key="downloadBoxplot")

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="position6")
        st.text_area("Annotation Text:", key="text6")
        st.button("Add", key="add6")
        st.button("Delete", key="delete6")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("log2_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("log2_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = boxplot_int(
                    data=data_to_use,
                    meta=meta_to_use,
                    outliers=st.session_state.get("outliers6", False),
                    header=st.session_state.get("header6", True),
                    legend=st.session_state.get("legend6", True),
                    width_cm=st.session_state.get("width6", 20),
                    height_cm=st.session_state.get("height6", 10),
                    dpi=st.session_state.get("dpi6", 300),
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def cov_plot_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "org_data" in st.session_state:
        available_levels.append("Protein")
    if "org_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        if not available_levels:
            st.info("No data available for plotting.")
            return

        st.checkbox("Show Outliers", value=False, key="outliers7")
        st.checkbox("Show Header", value=True, key="header7")
        st.checkbox("Show Legend", value=True, key="legend7")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level7")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="width7")
        st.number_input("Height (cm):", value=10, key="height7")
        st.number_input("DPI:", value=300, key="dpi7")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="format7")
        st.download_button(
            "Download Plot",
            data=b"",
            file_name=f"cov_plot.png",
            key="downloadCovPlot"
        )

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="position7")
        st.text_area("Annotation Text:", key="text7")
        st.button("Add", key="add7")
        st.button("Delete", key="delete7")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("org_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("org_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = cov_plot(
                    data=data_to_use,
                    meta=meta_to_use,
                    outliers=st.session_state.get("outliers7", False),
                    header=st.session_state.get("header7", True),
                    legend=st.session_state.get("legend7", True),
                    width_cm=st.session_state.get("width7", 20),
                    height_cm=st.session_state.get("height7", 10),
                    dpi=st.session_state.get("dpi7", 300),
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating Coefficient of Variation plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def principal_component_analysis_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "log2_data" in st.session_state:
        available_levels.append("Protein")
    if "log2_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        if not available_levels:
            st.info("No data available for PCA plot.")
            return

        st.checkbox("Show Legend", value=True, key="legend8")
        st.checkbox("Show Header", value=True, key="header8")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level8")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="width8")
        st.number_input("Height (cm):", value=10, key="height8")
        st.number_input("DPI:", value=300, key="dpi8")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="format8")
        st.download_button(
            "Download Plot",
            data=b"",
            file_name=f"pca_plot.png",
            key="downloadPCAPlot"
        )

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="position8")
        st.text_area("Annotation Text:", key="text8")
        st.button("Add", key="add8")
        st.button("Delete", key="delete8")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("log2_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("log2_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = pca_plot(
                    data=data_to_use,
                    meta=meta_to_use,
                    header=st.session_state.get("header8", True),
                    legend=st.session_state.get("legend8", True),
                    width_cm=st.session_state.get("width8", 20),
                    height_cm=st.session_state.get("height8", 10),
                    dpi=st.session_state.get("dpi8", 300),
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating PCA plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def abundance_plot_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "org_data" in st.session_state:
        available_levels.append("Protein")
    if "org_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        if not available_levels:
            st.info("No data available for abundance plot.")
            return

        st.checkbox("Show Legend", value=True, key="legend9")
        st.checkbox("Show Header", value=True, key="header9")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level9")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        st.number_input("Width (cm):", value=20, key="width9")
        st.number_input("Height (cm):", value=10, key="height9")
        st.number_input("DPI:", value=300, key="dpi9")
        st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="format9")
        st.download_button(
            "Download Plot",
            data=b"",
            file_name=f"abundance_plot.{st.session_state.get('format9','png')}",
            key="downloadAbPlot9"
        )

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="position9")
        st.text_area("Enter Text:", key="text9", height=150)
        st.button("Add", key="addText9")
        st.button("Delete", key="deleteText9")

    with col2:
        condition = st.selectbox("Choose Condition:", ["All Conditions"], key="condition9")
        protein_choices = st.session_state.get("protein_choices", [])
        selected_proteins = st.multiselect("Select Proteins:", options=protein_choices, key="protein9")

        if level == "Protein":
            data_to_use = st.session_state.get("org_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("org_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = abundance_plot(
                    data=data_to_use,
                    meta=meta_to_use,
                    workflow=level,
                    width_cm=st.session_state.get("width9", 20),
                    height_cm=st.session_state.get("height9", 10),
                    dpi=st.session_state.get("dpi9", 300),
                    legend=st.session_state.get("legend9", True),
                    header=st.session_state.get("header9", True)
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating abundance plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")


def correlation_plot_ui():
    col1, col2 = st.columns([1, 2])

    available_levels = []
    if "org_data" in st.session_state:
        available_levels.append("Protein")
    if "org_data3" in st.session_state:
        available_levels.append("Phosphosite")

    with col1:
        if not available_levels:
            st.info("No data available for correlation plot.")
            return

        change_display = st.checkbox(
            "Change Display",
            value=st.session_state.get("change_display12", False),
            key="changeDisplay12"
        )

        toggle_id = st.checkbox(
            "Toggle ID",
            value=st.session_state.get("toggle_id12", True),
            key="toggleId12"
        )

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level12")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        width = st.number_input("Width (cm):", min_value=4, max_value=20, value=10, key="width12")
        height = st.number_input("Height (cm):", min_value=4, max_value=20, value=8, key="height12")
        dpi = st.number_input("DPI:", min_value=50, max_value=300, value=100, key="dpi12")

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition12")
        st.text_area("Annotation Text:", value="", height=100, key="text12")
        st.button("Add", key="addText12")
        st.button("Delete", key="deleteText12")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("org_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("org_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = corr_plot(
                    data_to_use,
                    meta_to_use,
                    method=False,
                    id=toggle_id,  # use the variable from the checkbox directly
                    full_range=False,
                    width=width,
                    height=height,
                    dpi=dpi
                )
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating correlation plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")
