from ui.i1_data_ui import *
from ui.i2_qc_pipline_ui import *
from ui.i3_statistical_analysis_ui import *
from ui.i4_peptide_ui import *
from ui.i5_single_protein_ui import *
from ui.i6_phospho_specific_ui import *
from ui.i7_tables import *

def render_main_tabs():
    main_tabs = st.tabs(
        ["Data", "QC Pipeline", "Statistical Analysis", "Peptide Level", "Single Protein", "Phospho-specific",
         "Tables"])

    with main_tabs[0]:
        data_ui()

    with main_tabs[1]:
        qc_pipeline_ui()

    with main_tabs[2]:
        statistical_analysis_ui()

    with main_tabs[3]:
        peptide_level_ui()

    with main_tabs[4]:
        single_protein_ui()

    with main_tabs[5]:
        phospho_specific_ui()

    with main_tabs[6]:
        tables_ui()