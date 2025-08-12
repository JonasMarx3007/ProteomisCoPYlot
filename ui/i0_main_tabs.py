import time
from ui.i1_data_ui import *
from ui.i2_qc_pipline_ui import *
from ui.i3_statistical_analysis_ui import *
from ui.i4_single_protein_ui import *
from ui.i5_phospho_specific_ui import *

def render_main_tabs():
    if st.sidebar.button("Stop App"):
        st.sidebar.write("Stopping app...")
        time.sleep(3)
        os._exit(0)

    main_tabs = st.tabs(["Data", "QC Pipeline", "Statistical Analysis", "Single Protein", "Phospho-specific"])

    with main_tabs[0]:
        data_ui()

    with main_tabs[1]:
        qc_pipeline_ui()

    with main_tabs[2]:
        statistical_analysis_ui()

    with main_tabs[3]:
        single_protein_ui()

    with main_tabs[4]:
        phospho_specific_ui()