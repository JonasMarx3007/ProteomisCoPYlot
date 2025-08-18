import streamlit as st
from utils.functions import phossite_coverage_plot, simple_phos_site_plot

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
        plot_width = st.number_input("Width (inch):", value=6, key="plotWidth19")
        plot_height = st.number_input("Height (inch):", value=4, key="plotHeight19")
        plot_dpi = st.number_input("DPI:", value=100, key="plotDPI19")

        text_position = st.selectbox("Position:", options={"Above": "up", "Below": "down"}, key="textPosition19")
        text_input = st.text_area("Enter text:", value="", key="text19", height=100)
        st.button("Add", key="addText19")
        st.button("Delete", key="deleteText19")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            fig = simple_phos_site_plot(
                data,
                filter_value=cutoff,
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("No data available to plot.")


def phossite_coverage_plot_ui():
    col1, col2 = st.columns([1, 2])

    with col1:
        toggle_id = st.checkbox("Show IDs", value=False, key="toggle_id21")
        toggle_header = st.checkbox("Show Header", value=True, key="toggle_header21")
        toggle_legend = st.checkbox("Show Legend", value=True, key="toggle_legend21")

        plot_type = st.selectbox("Type:", options=["Normal", "Summary"], key="type21")
        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth21")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight21")
        plot_dpi = st.number_input("DPI:", value=300, key="plotDPI21")
        plot_format = st.selectbox("File Format:", options=["png", "jpg", "svg", "pdf"], key="plotFormat21")

        st.download_button(
            "Download Plot", data=b"", file_name=f"phossite_plot.{plot_format}", key="downloadPhossitePlot"
        )

        text_position = st.selectbox("Position:", options={"Above": "up", "Below": "down"}, key="textPosition21")
        text_input = st.text_area("Enter text:", value="", key="text21", height=100)
        st.button("Add", key="addText21")
        st.button("Delete", key="deleteText21")

    with col2:
        if "log2_data3" in st.session_state and st.session_state["log2_data3"] is not None:
            data = st.session_state["log2_data3"]
            meta = st.session_state["meta2"]

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

            st.pyplot(fig)
        else:
            st.info("No data available to plot.")
