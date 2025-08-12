import streamlit as st

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
        st.number_input("Cutoff:", value=0, key="cutoff19")
        st.selectbox("Position:", options={"Above": "up", "Below": "down"}, key="textPosition19")
        st.text_area("Enter text:", value="", key="text19", height=100)
        st.button("Add", key="addText19")
        st.button("Delete", key="deleteText19")
    with col2:
        st.subheader("Phossite Plot")
        st.empty()

def phossite_coverage_plot_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.button("Toggle IDs", key="toggle_id21")
        st.button("Toggle Header", key="toggle_header21")
        st.button("Toggle Legend", key="toggle_legend21")
        st.selectbox("Type:", options=["Normal", "Summary"], key="type21")
        st.number_input("Width (cm):", value=20, key="plotWidth21")
        st.number_input("Height (cm):", value=10, key="plotHeight21")
        st.number_input("DPI:", value=300, key="plotDPI21")
        st.selectbox("File Format:", options=["png", "jpg", "svg", "pdf"], key="plotFormat21")
        st.download_button("Download Plot", data=b"", file_name="phossite_plot.png", key="downloadPhossitePlot")
        st.selectbox("Position:", options={"Above": "up", "Below": "down"}, key="textPosition21")
        st.text_area("Enter text:", value="", key="text21", height=100)
        st.button("Add", key="addText21")
        st.button("Delete", key="deleteText21")
    with col2:
        st.subheader("Phossite Coverage Plot")
        st.empty()