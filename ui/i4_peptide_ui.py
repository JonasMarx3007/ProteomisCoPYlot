import streamlit as st
from utils.functions import rt_plot, modification_plot, missed_cleavage_plot

#MAIN
def peptide_level_ui():
    peptide_tabs = st.tabs([
        "RT Plot",
        "Modification Plot",
        "Missed Cleavage Plot"
    ])

    with peptide_tabs[0]:
        rt_plot_ui()
    with peptide_tabs[1]:
        modification_plot_ui()
    with peptide_tabs[2]:
        missed_cleavage_plot_ui()

#SUB
def rt_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.button("Add Line", key="line14")
        st.button("Toggle Header", key="toggle_header14")
        st.markdown("---")

        plot_type = st.selectbox("Select Plot Type:", ["Scatter Plot", "Hexbin Plot", "Density Plot"], key="style14")
        bins = st.number_input("Bins:", value=1000, key="bins14")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth14")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight14")
        plot_dpi = st.number_input("DPI:", value=300, key="plotDPI14")
        plot_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], index=0, key="plotFormat14")
        st.download_button("Download Plot", data=b"", file_name=f"RTPlot.{plot_format}", key="downloadRTPlot")
        st.markdown("---")

        text_position = st.selectbox("Position:", {"Above": "up", "Below": "down"}, key="textPosition14")
        text_input = st.text_area("Annotation Text:", value="", height=100, key="text14")
        st.button("Add", key="addText14")
        st.button("Delete", key="deleteText14")

    with col2:
        if "rt_data" in st.session_state and st.session_state["rt_data"] is not None:
            fig = rt_plot(
                st.session_state["rt_data"],
                plot_type=plot_type,
                bins=bins,
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("No data available for RT Plot.")

def modification_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.button("Toggle ID", key="toggle_id15")
        st.button("Toggle Header", key="toggle_header15")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth15")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight15")
        plot_dpi = st.number_input("DPI:", value=300, key="plotDPI15")
        plot_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], index=0, key="plotFormat15")
        st.download_button("Download Plot", data=b"", file_name=f"ModPlot.{plot_format}", key="downloadModPlot")
        st.markdown("---")

        text_position = st.selectbox("Position:", {"Above": "up", "Below": "down"}, key="textPosition15")
        text_input = st.text_area("Enter text:", value="", height=200, key="text15")
        st.button("Add", key="addText15")
        st.button("Delete", key="deleteText15")

    with col2:
        if "mod_data" in st.session_state and st.session_state["mod_data"] is not None:
            fig = modification_plot(
                st.session_state["mod_data"],
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("No data available for Modification Plot.")

def missed_cleavage_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.button("Toggle ID", key="toggle_id16")
        st.button("Toggle Header", key="toggle_header16")
        st.markdown("---")

        text_size = st.number_input("Text Size:", value=3.88, key="text_size16")
        st.button("Toggle Text", key="toggle_text16")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width (cm):", value=20, key="plotWidth16")
        plot_height = st.number_input("Height (cm):", value=10, key="plotHeight16")
        plot_dpi = st.number_input("DPI:", value=300, key="plotDPI16")
        plot_format = st.selectbox("File Format:", ["png", "jpg", "svg", "pdf"], index=0, key="plotFormat16")
        st.download_button("Download Plot", data=b"", file_name=f"MCPlot.{plot_format}", key="downloadMCPlot")
        st.markdown("---")

        text_position = st.selectbox("Position:", {"Above": "up", "Below": "down"}, key="textPosition16")
        text_input = st.text_area("Enter text:", value="", height=200, key="text16")
        st.button("Add", key="addText16")
        st.button("Delete", key="deleteText16")

    with col2:
        if "mc_data" in st.session_state and st.session_state["mc_data"] is not None:
            fig = missed_cleavage_plot(
                st.session_state["mc_data"],
                text_size=text_size,
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("No data available for Missed Cleavage Plot.")
