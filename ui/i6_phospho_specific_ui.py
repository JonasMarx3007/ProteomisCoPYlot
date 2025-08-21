import streamlit as st
import matplotlib.pyplot as plt
import io
from utils.functions import phossite_coverage_plot, simple_phos_site_plot, phossite_coverage_plot_summary

#MAIN
def phospho_specific_ui():
    phos_tabs = st.tabs([
        "Phossite Plot",
        "Phossite Coverage Plot"
    ])

    with phos_tabs[0]:
        phossite_plot_ui()
    with phos_tabs[1]:
        phossite_coverage_plot_ui()


#SUB
def phossite_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        cutoff = st.number_input("Cutoff:", value=0, key="cutoff19")
        st.markdown("---")

        plot_width = st.number_input("Width (cm):", value=15, key="plotWidth19")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight19")
        plot_dpi = st.number_input("DPI:", value=100, key="plotDPI19")
        plot_format = st.selectbox("Format:", ["png", "jpg", "svg", "pdf"], index=0)

        download_placeholder = st.empty()

        st.markdown("---")
        text_input = st.text_area("Enter text:", value="", key="text19", height=100)
        st.button("Add", key="addText19")
        st.button("Delete", key="deleteText19")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            try:
                fig = simple_phos_site_plot(
                    data,
                    filter_value=cutoff,
                    width_cm=plot_width,
                    height_cm=plot_height,
                    dpi=plot_dpi
                )
                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plot_format, dpi=plot_dpi, bbox_inches="tight")
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"phossite_plot.{plot_format}",
                    mime=f"image/{plot_format}"
                )

                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("No data available to plot.")


def phossite_coverage_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        toggle_id = st.checkbox("Show IDs", value=False, key="toggle_id21")
        toggle_header = st.checkbox("Show Header", value=True, key="toggle_header21")
        toggle_legend = st.checkbox("Show Legend", value=True, key="toggle_legend21")
        st.markdown("---")

        plot_type = st.selectbox("Type:", options=["Normal", "Summary"], key="type21")
        st.markdown("---")

        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth21")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight21")
        plot_dpi = st.number_input("DPI:", value=300, key="plotDPI21")
        plot_format = st.selectbox("File Format:", options=["png", "jpg", "svg", "pdf"], key="plotFormat21")
        download_placeholder = st.empty()
        st.markdown("---")

        text_position = st.selectbox("Position:", options=["Above", "Below"], key="textPosition21")
        text_input = st.text_area("Enter text:", value="", key="text21", height=100)
        st.button("Add", key="addText21")
        st.button("Delete", key="deleteText21")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            meta = st.session_state["meta2"]

            try:
                if plot_type == "Normal":
                    fig = phossite_coverage_plot(
                        data,
                        meta,
                        id=toggle_id,
                        header=toggle_header,
                        legend=toggle_legend,
                        width=plot_width,
                        height=plot_height,
                        dpi=plot_dpi
                    )
                elif plot_type == "Summary":
                    fig = phossite_coverage_plot_summary(
                        data,
                        meta,
                        header=toggle_header,
                        legend=toggle_legend,
                        width=plot_width,
                        height=plot_height,
                        dpi=plot_dpi
                    )

                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plot_format, dpi=plot_dpi)
                buf.seek(0)
                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"phossite_plot.{plot_format}",
                    mime=f"image/{plot_format}"
                )

                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("No data available to plot.")