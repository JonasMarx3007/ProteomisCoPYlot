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

        plotWidth = st.number_input("Width (cm):", value=15, key="plotWidth19")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight19")
        plotDPI = st.number_input("DPI:", value=100, key="plotDPI19")
        plotFormat = st.selectbox("Format:", ["png", "jpg", "svg", "pdf"], index=0, key="plotFormat19")

        download_placeholder = st.empty()

        st.markdown("---")
        text = st.text_area("Enter text:", value="", key="text19", height=100)
        st.button("Add", key="addText19")
        st.button("Delete", key="deleteText19")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            try:
                fig = simple_phos_site_plot(
                    data,
                    filter_value=cutoff,
                    width_cm=plotWidth,
                    height_cm=plotHeight,
                    dpi=plotDPI
                )
                st.pyplot(fig)

                buf = io.BytesIO()
                fig.savefig(buf, format=plotFormat, dpi=plotDPI, bbox_inches="tight")
                buf.seek(0)

                download_placeholder.download_button(
                    "Download Plot",
                    data=buf,
                    file_name=f"phossite_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
                )

                plt.close(fig)
            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("No data available to plot.")


def phossite_coverage_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        id = st.checkbox("Toggle ID", value=False, key="id21")
        header = st.checkbox("Toggle Header", value=True, key="header21")
        legend = st.checkbox("Toggle Legend", value=True, key="legend21")
        st.markdown("---")

        type = st.selectbox("Type:", options=["Normal", "Summary"], key="type21")
        st.markdown("---")

        plotWidth = st.number_input("Width (cm):", value=20, key="plotWidth21")
        plotHeight = st.number_input("Height (cm):", value=10, key="plotHeight21")
        plotDPI = st.number_input("DPI:", value=300, key="plotDPI21")
        plotFormat = st.selectbox("File Format:", options=["png", "jpg", "svg", "pdf"], key="plotFormat21")
        download_placeholder = st.empty()
        st.markdown("---")

        textPosition = st.selectbox("Position:", options=["Above", "Below"], key="textPosition21")
        text = st.text_area("Enter text:", value="", key="text21", height=100)
        st.button("Add", key="addText21")
        st.button("Delete", key="deleteText21")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            meta = st.session_state["meta2"]

            try:
                if type == "Normal":
                    fig = phossite_coverage_plot(
                        data,
                        meta,
                        id=id,
                        header=header,
                        legend=legend,
                        width=plotWidth,
                        height=plotHeight,
                        dpi=plotDPI
                    )
                elif type == "Summary":
                    fig = phossite_coverage_plot_summary(
                        data,
                        meta,
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
                    file_name=f"phossite_plot.{plotFormat}",
                    mime=f"image/{plotFormat}"
                )

                plt.close(fig)

            except Exception as e:
                st.error(f"Error generating plot: {e}")
        else:
            st.info("No data available to plot.")