import streamlit as st

#MAIN
def single_protein_ui():
    sp_tabs = st.tabs([
        "Protein Boxplot",
        "Protein Lineplot"
    ])

    with sp_tabs[0]:
        protein_box_ui()
    with sp_tabs[1]:
        protein_line_ui()


#SUB
def protein_box_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.selectbox("Select Protein:", options=["Protein A", "Protein B"], key="protein17")
        st.multiselect("Select Conditions:", options=["Cond1", "Cond2"], key="conditions17")
        st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level17")
    with col2:
        st.subheader("Protein Boxplot")
        st.empty()

def protein_line_ui():
    col1, col2 = st.columns([1, 2])
    with col1:
        st.multiselect("Select Proteins:", options=["Protein A", "Protein B"], key="protein18")
        st.multiselect("Select Conditions:", options=["Cond1", "Cond2"], key="conditions18")
        st.selectbox("Level:", options=["Protein", "Phosphosite"], key="level18")
    with col2:
        st.subheader("Protein Lineplot")
        st.empty()