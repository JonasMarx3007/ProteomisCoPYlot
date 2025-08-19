import streamlit as st
from utils.functions import rt_vs_pred_rt_plot, modification_plot, missed_cleavage_plot

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
        add_line = st.checkbox("Add Line", key="line14")
        header = st.checkbox("Show Header", key="header14")
        st.markdown("---")

        plot_type = st.selectbox("Select Plot Type:", ["Scatter Plot", "Hexbin Plot", "Density Plot"], key="style14")
        bins = st.number_input("Bins:", value=1000, key="bins14")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width:", value=8, key="plotWidth14")
        plot_height = st.number_input("Height:", value=6, key="plotHeight14")
        plot_dpi = st.number_input("DPI:", value=100, key="plotDPI14")

    with col2:
        if "data2" in st.session_state and st.session_state["data2"] is not None:
            fig = rt_vs_pred_rt_plot(
                st.session_state["data2"],
                method=plot_type,
                bins=bins,
                add_line=add_line,
                header=header,
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
        id_toggle = st.checkbox("Toggle ID", value=False,  key="toggle_id15")
        header = st.checkbox("Show Header", key="header15")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width:", value=10, key="plotWidth15")
        plot_height = st.number_input("Height:", value=6, key="plotHeight15")
        plot_dpi = st.number_input("DPI:", value=100, key="plotDPI15")

    with col2:
        if "data2" in st.session_state and st.session_state["data2"] is not None and \
           "meta" in st.session_state and st.session_state["meta"] is not None:
            fig = modification_plot(
                st.session_state["data2"],
                st.session_state["meta"],
                id=id_toggle,
                header=header,
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("Both data and meta information are required for the Modification Plot.")


def missed_cleavage_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        header = st.checkbox("Show Header", key="header16")
        id_toggle = st.checkbox("Toggle ID", value=False, key="id_toggle16")
        text_toggle = st.checkbox("Show Text", key="toggle_text16")
        text_size = st.number_input("Text Size:", value=8, key="text_size16")
        st.markdown("---")

        st.subheader("Plot Size & Resolution")
        plot_width = st.number_input("Width:", value=10, key="plotWidth16")
        plot_height = st.number_input("Height:", value=6, key="plotHeight16")
        plot_dpi = st.number_input("DPI:", value=100, key="plotDPI16")

    with col2:
        if "data2" in st.session_state and st.session_state["data2"] is not None and \
           "meta" in st.session_state and st.session_state["meta"] is not None:
            fig = missed_cleavage_plot(
                st.session_state["data2"],
                st.session_state["meta"],
                text=text_toggle,
                text_size=text_size,
                header=header,
                id=id_toggle,
                width=plot_width,
                height=plot_height,
                dpi=plot_dpi
            )
            st.pyplot(fig)
        else:
            st.info("Both data and meta information are required for the Missed Cleavage Plot.")