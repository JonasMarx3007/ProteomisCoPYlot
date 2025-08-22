import streamlit as st
import matplotlib.pyplot as plt
import io
from utils.functions import coverage_plot, missing_value_plot, histo_int, boxplot_int, cov_plot, pca_plot, abundance_plot, corr_plot, coverage_plot_pep, missing_value_plot_prec, missing_value_plot_pep, coverage_plot_summary, interactive_abundance_plot, pca_plot_interactive, boxplot_int_single

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
        id = st.checkbox("Toggle IDs", value=False, key="id3")
        header = st.checkbox("Toggle Header", value=True, key="header3")
        legend = st.checkbox("Toggle Legend", value=True, key="legend3")

        st.markdown("---")

        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level3")
        else:
            st.info("No data available for plotting.")
            return

        if level == "Peptide":
            plot_type = "Normal"
        else:
            plot_type = st.selectbox("Type:", ["Normal", "Summary"], key="type3")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth3")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight3")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI3")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat3")

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
        elif level == "Peptide":
            data_to_use = st.session_state.get("data2", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("data3", None)
            meta_to_use = st.session_state.get("meta2", None)
        else:
            data_to_use, meta_to_use = None, None

        if data_to_use is not None and meta_to_use is not None:
            try:
                if plot_type == "Normal":
                    if level == "Peptide":
                        fig = coverage_plot_pep(
                            data=data_to_use,
                            meta=meta_to_use,
                            id=id,
                            header=header,
                            legend=legend,
                            width=plotWidth,
                            height=plotHeight,
                            dpi=plotDPI,
                            plot_colors=st.session_state["selected_colors"]
                        )
                    else:
                        fig = coverage_plot(
                            data=data_to_use,
                            meta=meta_to_use,
                            id=id,
                            header=header,
                            legend=legend,
                            width=plotWidth,
                            height=plotHeight,
                            dpi=plotDPI,
                            plot_colors=st.session_state["selected_colors"]
                        )
                else:
                    fig = coverage_plot_summary(
                        data=data_to_use,
                        meta=meta_to_use,
                        header=header,
                        plot_colors=st.session_state["selected_colors"],
                        width=plotWidth,
                        height=plotHeight,
                        dpi=plotDPI
                    )

                st.pyplot(fig)

                import io
                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"coverage_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
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
        header = st.checkbox("Toggle Header", value=True, key="header4")
        st.markdown("---")
        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level4")
        else:
            st.info("No data available for plotting.")
            return

        bins = st.number_input("Bin missing values (optional):", value=0, min_value=0, key="bins4")
        st.markdown("---")
        plotText = st.checkbox("Toggle Text", value=True, key="plotText4")
        text_size = st.number_input("Text Size:", value=8, key="text_size4")
        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth4")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight4")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI4")
        plotFormat = st.selectbox("Download Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat4")

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
                    bin=bins,
                    header=header,
                    text=plotText,
                    text_size=text_size,
                    width=plotWidth / 2.54,
                    height=plotHeight / 2.54,
                    dpi=plotDPI
                )

                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                buf.seek(0)
                st.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"missval_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
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
                    bin=bins,
                    header=header,
                    text=plotText,
                    text_size=text_size,
                    width=plotWidth / 2.54,
                    height=plotHeight / 2.54,
                    dpi=plotDPI
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
        header = st.checkbox("Show Header", value=True, key="header5")
        legend = st.checkbox("Show Legend", value=True, key="legend5")
        st.markdown("---")
        if available_levels:
            level = st.selectbox("Level:", available_levels, key="level5")
        else:
            st.info("No data available for plotting.")
            return
        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth5")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight5")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI5")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat5")

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
                fig = histo_int(
                    data=data_to_use,
                    meta=meta_to_use,
                    header=header,
                    legend=legend,
                    width=plotWidth,
                    height=plotHeight,
                    dpi=plotDPI,
                    plot_colors=st.session_state.get("selected_colors")
                )
                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI, bbox_inches="tight")
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"hist_intensity.{plotFormat}",
                    mime=f"image/{plotFormat}"
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

        outliers = st.checkbox("Toggle Outliers", value=False, key="outliers6")
        header = st.checkbox("Toggle Header", value=True, key="header6")
        id = st.checkbox("Toggle ID", value=True, key="id6")
        legend = st.checkbox("Toggle Legend", value=True, key="legend6")

        st.markdown("---")
        mode = st.selectbox("Plot Type:", ["Mean", "Single"], key="mode6")
        level = st.selectbox("Level:", available_levels, key="level6")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth6")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight6")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI6")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat6")
        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition6")
        st.text_area("Annotation Text:", key="text6")
        st.button("Add", key="addText6")
        st.button("Delete", key="deleteText6")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("log2_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("log2_data3", None)
            meta_to_use = st.session_state.get("meta2", None)
        else:
            data_to_use, meta_to_use = None, None

        if data_to_use is not None and meta_to_use is not None:
            try:
                if mode == "Mean":
                    fig = boxplot_int(
                        data=data_to_use,
                        meta=meta_to_use,
                        outliers=st.session_state.get("outliers6", False),
                        header=st.session_state.get("header6", True),
                        width_cm=st.session_state.get("plotWidth6", 20),
                        height_cm=st.session_state.get("plotHeight6", 10),
                        dpi=st.session_state.get("plotDPI6", 300),
                        plot_colors=st.session_state["selected_colors"]
                    )
                else:
                    fig = boxplot_int_single(
                        data=data_to_use,
                        meta=meta_to_use,
                        outliers=st.session_state.get("outliers6", False),
                        id=st.session_state.get("id6", True),
                        header=st.session_state.get("header6", True),
                        legend=st.session_state.get("legend6", True),
                        width_cm=st.session_state.get("plotWidth6", 20),
                        height_cm=st.session_state.get("plotHeight6", 10),
                        dpi=st.session_state.get("plotDPI6", 300),
                        plot_colors=st.session_state["selected_colors"]
                    )

                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"boxplot_intensity.{plotFormat}",
                    mime=f"image/{plotFormat}"
                )

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

        outliers = st.checkbox("Toggle Outliers", value=False, key="outliers7")
        header = st.checkbox("Toggle Header", value=True, key="header7")
        legend = st.checkbox("Toggle Legend", value=True, key="legend7")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level7")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth7")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight7")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI7")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat7")
        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition7")
        st.text_area("Annotation Text:", key="text7")
        st.button("Add", key="addText7")
        st.button("Delete", key="deleteText7")

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
                    width_cm=st.session_state.get("plotWidth7", 20),
                    height_cm=st.session_state.get("plotHeight7", 10),
                    dpi=st.session_state.get("plotDPI7", 300),
                    plot_colors=st.session_state["selected_colors"]
                )
                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"cov_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
                )

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

        legend = st.checkbox("Toggle Legend", value=True, key="legend8")
        header = st.checkbox("Toggle Header", value=True, key="header8")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level8")

        st.markdown("---")
        type = st.selectbox("Type:", ["Normal", "Interactive"], key="type8")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth8")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight8")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI8")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat8")
        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition8")
        st.text_area("Annotation Text:", key="text8")
        st.button("Add", key="addText8")
        st.button("Delete", key="deleteText8")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("log2_data", None)
            meta_to_use = st.session_state.get("meta", None)
        elif level == "Phosphosite":
            data_to_use = st.session_state.get("log2_data3", None)
            meta_to_use = st.session_state.get("meta2", None)

        if data_to_use is not None and meta_to_use is not None:
            try:
                if type == "Normal":
                    fig = pca_plot(
                        data=data_to_use,
                        meta=meta_to_use,
                        header=header,
                        legend=legend,
                        width_cm=plotWidth,
                        height_cm=plotHeight,
                        dpi=plotDPI,
                        plot_colors=st.session_state["selected_colors"]
                    )
                    st.pyplot(fig)

                    buf = io.BytesIO()
                    fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                    buf.seek(0)

                    download_placeholder.download_button(
                        "Download Plot",
                        data=buf,
                        file_name=f"pca_plot.{plotFormat}",
                        mime=f"image/{plotFormat}"
                    )

                    plt.close(fig)

                elif type == "Interactive":
                    fig = pca_plot_interactive(
                        data=data_to_use,
                        meta=meta_to_use,
                        plot_colors=st.session_state["selected_colors"],
                        header=header,
                        legend=legend
                    )

                    cm_to_px = 96 / 2.54
                    fig.update_layout(
                        width=int(plotWidth * cm_to_px),
                        height=int(plotHeight * cm_to_px),
                        showlegend=legend
                    )

                    st.plotly_chart(fig, use_container_width=False)

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

        legend = st.checkbox("Toggle Legend", value=True, key="legend9")
        header = st.checkbox("Toggle Header", value=True, key="header9")

        st.markdown("---")
        level = st.selectbox("Level:", available_levels, key="level9")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth9")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight9")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI9")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat9")
        download_placeholder = st.empty()

        st.markdown("---")
        st.selectbox("Position:", ["Above", "Below"], key="textPosition9")
        st.text_area("Enter Text:", key="text9", height=150)
        st.button("Add", key="addText9")
        st.button("Delete", key="deleteText9")

    with col2:
        if level == "Protein":
            data_to_use = st.session_state.get("org_data", None)
            meta_to_use = st.session_state.get("meta", None)
            feature_col = "ProteinNames"
        else:
            data_to_use = st.session_state.get("org_data3", None)
            meta_to_use = st.session_state.get("meta2", None)
            feature_col = "PTM_Collapse_key"

        if data_to_use is not None and meta_to_use is not None:
            conditions = meta_to_use["condition"].unique().tolist()
            condition = st.selectbox("Choose Condition:", ["All Conditions"] + conditions, key="condition9")

            feature_choices = data_to_use[feature_col].dropna().unique().tolist()
            selected_features = st.multiselect(f"Select {level}(s) to highlight:", options=feature_choices, key="feature9")

            try:
                if condition == "All Conditions":
                    fig = abundance_plot(
                        data=data_to_use,
                        meta=meta_to_use,
                        workflow=level,
                        width_cm=plotWidth,
                        height_cm=plotHeight,
                        dpi=plotDPI,
                        legend=legend,
                        header=header,
                        plot_colors=st.session_state["selected_colors"]
                    )
                    st.pyplot(fig)

                    buf = io.BytesIO()
                    fig.savefig(buf, format=plotFormat, dpi=plotDPI)
                    buf.seek(0)

                    download_placeholder.download_button(
                        "Download Plot",
                        data=buf,
                        file_name=f"abundance_plot.{plotFormat}",
                        mime=f"image/{plotFormat}"
                    )
                    plt.close(fig)
                else:
                    fig = interactive_abundance_plot(
                        data=data_to_use,
                        meta=meta_to_use,
                        condition=condition,
                        workflow=level,
                        search=selected_features
                    )
                    st.plotly_chart(fig, use_container_width=True)
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


        id = st.checkbox(
            "Toggle ID",
            value=st.session_state.get("id12", True),
            key="id12"
        )

        st.markdown("---")
        method = st.selectbox(
            "Correlation Plot Method:",
            options=["Matrix", "Ellipse"],
            index=0,
            key="method12"
        )
        level = st.selectbox("Level:", available_levels, key="level12")

        st.markdown("---")
        st.header("Plot Size & Resolution")
        plotWidth = st.number_input("Width (cm):", min_value=4, max_value=20, value=10, key="plotWidth12")
        plotHeight = st.number_input("Height (cm):", min_value=4, max_value=20, value=8, key="plotHeight12")
        plotDPI = st.number_input("DPI:", min_value=50, max_value=300, value=100, key="plotDPI12")
        plotFormat = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], key="plotFormat12")
        download_placeholder = st.empty()

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
        else:
            data_to_use, meta_to_use = None, None

        if data_to_use is not None and meta_to_use is not None:
            try:
                fig = corr_plot(
                    data=data_to_use,
                    meta=meta_to_use,
                    method=method,
                    id=id,
                    full_range=False,
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
                    file_name=f"correlation_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
                )

                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating correlation plot: {e}")
        else:
            st.info("Please upload the corresponding data and metadata to see the plot.")
