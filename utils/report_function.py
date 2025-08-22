import io
import base64
from utils.functions import *
from datetime import datetime

def report_function():
    title = st.session_state.get("report_title", "").strip() or "Untitled Report"
    author = st.session_state.get("report_author", "").strip() or "Unknown Author"
    date_str = datetime.now().strftime("%Y-%m-%d")
    version = "Proteomics CoPYlot V0.01"

    html_parts = [
        f"<html><head><title>{title}</title></head><body>",
        f"<h1 style='font-size:2.5em'>{title}</h1>",
        f"<p>Author: {author}</p>",
        f"<p>Date: {date_str}</p>",
        f"<p>{version}</p>"
    ]

    pep_num = 1

    if "data2" in st.session_state:
        prot_num = 2
    else:
        prot_num = 1

    if "data2" and "data3" in st.session_state:
        phos_num = 3
    elif "data2" or "data3" in st.session_state:
        phos_num = 2
    else:
        phos_num = 1

    toc = ["<br><h1>Table of Contents</h1>"]

    if "data2" in st.session_state:
        toc.append(f'<h2 style="font-weight: normal;">{pep_num} Peptide Level Plots</h2>')
        toc.append(f'<h3 style="font-weight: normal;">{pep_num}.1 Coverage Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{pep_num}.2 Missing Value Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{pep_num}.3 Retention Time Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{pep_num}.4 Modification Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{pep_num}.5 Missed Cleavage Plot</h3>')

    if "data" in st.session_state:
        toc.append(f'<h2 style="font-weight: normal;">{prot_num} Protein Level Plots</h2>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.1 Coverage Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.2 Missing Values Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.3 Histogram Intensity</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.4 Boxplot Intensity</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.5 Coefficient of Variation Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.6 PCA Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.7 Abundance Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{prot_num}.8 Correlation Plot</h3>')

    if "data3" in st.session_state:
        toc.append(f'<h2 style="font-weight: normal;">{phos_num} Phosphosite Level Plots</h2>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.1 Overview Phosphosites</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.2 Coverage Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.3 Coverage Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.4 Missed Value Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.5 Histogram Intensity</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.6 Boxplot Intensity</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.7 Coefficient of Variation Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.8 PCA Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.9 Abundance Plot</h3>')
        toc.append(f'<h3 style="font-weight: normal;">{phos_num}.10 Correlation Plot</h3>')

    toc.append("<br>")
    toc = "\n".join(toc)

    html_parts.append(toc)

    #Pep - Coverage
    if "data2" in st.session_state and "meta" in st.session_state:
        fig3pep = coverage_plot_pep(
            data=st.session_state.get("data2"),
            meta=st.session_state.get("meta"),
            id=st.session_state.get("id3", False),
            header=st.session_state.get("header3", True),
            legend=st.session_state.get("legend3", True),
            width=st.session_state.get("plotWidth3", 20),
            height=st.session_state.get("plotHeight3", 10),
            dpi=st.session_state.get("plotDPI3", 300),
            plot_colors=st.session_state["selected_colors"]
        )
        fig3pep.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig3pep.savefig(buf, format="png", dpi=st.session_state.get("plotDPI3", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h1>{pep_num} Peptide Level</h1>")
        html_parts.append(f"<h2>{pep_num}.1 Coverage Plot - Peptide Level</h2>")
        html_parts.append("<p>Numbers of peptide Identifications across all samples is a measure of non-missing measurements by replicate. Protein inference, identification and quantification are all based on reproducible peptide identification and the peptide’s ion intensity.</p>")
        html_parts.append("<p>Low peptide counts in a sample / MS-run may suggest a systematic flaw in the experiment (e.g. low protein input material, protein to peptide processing errors, technical problems during LC-MS data acquisition).</p>")
        html_parts.append("<p>Ideally, numbers of peptide identifications in all samples are high and of similar count.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Pep - Missing Value
    if "data2" in st.session_state and "meta" in st.session_state:
        fig4pep = missing_value_plot_pep(
            data=st.session_state.get("data2"),
            meta=st.session_state.get("meta"),
            bin=st.session_state.get("bins4", 0),
            header=st.session_state.get("header4", True),
            text=st.session_state.get("plotText4", True),
            text_size=st.session_state.get("text_size4", 8),
            width=st.session_state.get("plotWidth4", 20) / 2.54,
            height=st.session_state.get("plotHeight4", 10) / 2.54,
            dpi=st.session_state.get("plotDPI4", 300)
        )

        fig4pep.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig4pep.savefig(buf, format="png", dpi=st.session_state.get("plotDPI4", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{pep_num}.2 Missing Value Plot - Peptide Level</h2>")
        html_parts.append("<p>The amount of missing values can be affected by the biological condition or by technical factors during sample processing methods or MS-measurements and it can vary largely between experiments. For a healthy experiment we expect: the distribution of available measurements by replicate to be similar across replicates, especially within the same biological/experimental condition. An unusually low value in one or more replicates can be symptomatic of technical problems and should be taken into account when interpreting the final differential protein expression results after pairwise comparison.</p>")
        html_parts.append("<p>A high frequency of zero missing peptides indicates reproducible identifications of given peptides in all samples, a prerequisite for protein quantification.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Pep - RT
    if "data2" in st.session_state and "meta" in st.session_state:
        fig14 = rt_vs_pred_rt_plot(
            st.session_state.get("data2"),
            method=st.session_state.get("type14", "Scatter Plot"),
            bins=st.session_state.get("bins14", 1000),
            add_line=st.session_state.get("line14", False),
            header=st.session_state.get("header14", True),
            width=st.session_state.get("plotWidth14", 8),
            height=st.session_state.get("plotHeight14", 6),
            dpi=st.session_state.get("plotDPI14", 100)
        )

        fig14.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig14.savefig(buf, format="png", dpi=st.session_state.get("plotDPI14", 100), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{pep_num}.3 Retention Time Plot - Peptide Level</h2>")
        html_parts.append("<p>During LC-MS/MS analysis complex peptide mixtures of a given sample are separated via reversed-phase liquid chromatography (LC). The requirement for reproducible peptide and protein identification and quantification is the reproducible separation of peptides over the entire gradient length utilized for analysis in all samples of the data set.</p>")
        html_parts.append("<p>Ideally, retention times do align over the entire gradient length. Significant deviations from this pattern indicate technical problems during liquid chromatography.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Pep - Mod
    if "data2" in st.session_state and "meta" in st.session_state:
        fig15 = modification_plot(
            st.session_state.get("data2"),
            st.session_state.get("meta"),
            id=st.session_state.get("id15", False),
            header=st.session_state.get("header15", True),
            legend=st.session_state.get("legend15", True),
            width=st.session_state.get("plotWidth15", 10),
            height=st.session_state.get("plotHeight15", 6),
            dpi=st.session_state.get("plotDPI15", 100),
        )

        fig15.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig15.savefig(buf, format="png", dpi=st.session_state.get("plotDPI15", 100), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{pep_num}.4 Modification Plot - Peptide Level</h2>")
        html_parts.append("<p>Amino acid modifications can occur via post-translational modifications in the investigated biological system / condition, during sample storage or protein to peptide processing for LC-MS analysis.</p>")
        html_parts.append("<p>Typically, only expected amino acid modifications will be taken into account to be present in the data set. Numbers of identified peptide modifications are expected to be in a similar range. Outliers may indicate experimental or technical abnormalities. Usually, we expect Carbamidomethylation of Cysteine residues as a constant modification – induced during sample processing - and Methionine oxidation as a dynamic modification into account, if not otherwise stipulated.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Pep - Missed
    if "data2" in st.session_state and "meta" in st.session_state:
        fig16 = missed_cleavage_plot(
            st.session_state.get("data2"),
            st.session_state.get("meta"),
            text=st.session_state.get("plotText16", True),
            text_size=st.session_state.get("text_size16", 8),
            header=st.session_state.get("header16", True),
            id=st.session_state.get("id16", True),
            width=st.session_state.get("plotWidth16", 10),
            height=st.session_state.get("plotHeight16", 6),
            dpi=st.session_state.get("plotDPI16", 100),
        )

        fig16.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig16.savefig(buf, format="png", dpi=st.session_state.get("plotDPI16", 100), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{pep_num}.5 Missed Cleavage Plot - Peptide Level</h2>")
        html_parts.append("<p>Sample processing methods in proteomic bottom-up experiments include a proteolytic digestion step for peptide generation. Missed cleavages of a given protein by a specific protease can occur, nevertheless the frequency of missed cleavages should be low and comparable in all samples.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - Coverage
    if "data" in st.session_state and "meta" in st.session_state:
        if st.session_state.get("type3", "Normal") == "Normal":
            fig3 = coverage_plot(
                data=st.session_state.get("data"),
                meta=st.session_state.get("meta"),
                id=st.session_state.get("id3", False),
                header=st.session_state.get("header3", True),
                legend=st.session_state.get("legend3", True),
                width=st.session_state.get("plotWidth3", 20),
                height=st.session_state.get("plotHeight3", 10),
                dpi=st.session_state.get("plotDPI3", 300),
                plot_colors=st.session_state.get("selected_colors", None)
            )
        elif st.session_state.get("type3", "Normal") == "Summary":
            fig3 = coverage_plot_summary(
                data=st.session_state.get("data"),
                meta=st.session_state.get("meta"),
                header=st.session_state.get("header3", True),
                plot_colors=st.session_state["selected_colors"],
                width=st.session_state.get("plotWidth3", 20),
                height=st.session_state.get("plotHeight3", 10),
                dpi=st.session_state.get("plotDPI3", 300)
            )

        fig3.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig3.savefig(buf, format="png", dpi=st.session_state.get("plotDPI3", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h1>{prot_num} Protein Level</h1>")
        html_parts.append(f"<h2>{prot_num}.1 Coverage Plot - Protein Level</h2>")
        html_parts.append("<p>Identification of proteins is a measure of the number of non-missing measurements by replicate.</p>")
        html_parts.append("<p>Low protein counts in a sample / LC-MS/MS run may suggest a systematic flaw in the experiment (low protein input at the beginning, technical problems during sample processing, etc.) that needs to be addressed prior to further interpretation of the results in the statistical analysis section.</p>")
        html_parts.append("<p>Ideally, numbers of protein identifications in all samples of the data set are high and in a similar range.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

        # Prot - Missing Value
        if "data" in st.session_state and "meta" in st.session_state:
            fig4 = missing_value_plot(
                data=st.session_state.get("data"),
                meta=st.session_state.get("meta"),
                bin=st.session_state.get("bins4", 0),
                header=st.session_state.get("header4", True),
                text=st.session_state.get("plotText4", True),
                text_size=st.session_state.get("text_size4", 8),
                width=st.session_state.get("plotWidth4", 20) / 2.54,
                height=st.session_state.get("plotHeight4", 10) / 2.54,
                dpi=st.session_state.get("plotDPI4", 300)
            )

            fig4.set_size_inches(12, 6)
            buf = io.BytesIO()
            fig4.savefig(buf, format="png", dpi=st.session_state.get("plotDPI4", 300), bbox_inches="tight")
            buf.seek(0)
            img_base64 = base64.b64encode(buf.read()).decode("utf-8")

            html_parts.append(f"<h2>{prot_num}.2 Missing Values Plot - Protein Level</h2>")
            html_parts.append("<p>Missing values on protein level may depend on the experimental set-up, e.g. differences in experimental conditions or on technical flaws during sample collection, sample processing, or measurements.</p>")
            html_parts.append("<p>For technical replicates a very low level of missing values can be expected. A high frequency of zero missing values indicates high reproducibility of protein identifications in all samples in the data set.</p>")
            html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - HistoInt
    if "log2_data" in st.session_state and "meta" in st.session_state:
        fig5 = histo_int(
            data=st.session_state.get("log2_data"),
            meta=st.session_state.get("meta"),
            header=st.session_state.get("header5", True),
            legend=st.session_state.get("legend5", True),
            width=st.session_state.get("plotWidth5", 20),
            height=st.session_state.get("plotHeight5", 10),
            dpi=st.session_state.get("plotDPI5", 300),
            plot_colors=st.session_state.get("selected_colors", None)
        )

        fig5.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig5.savefig(buf, format="png", dpi=st.session_state.get("plotDPI5", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.3 Histogram Intensity - Protein Level</h2>")
        html_parts.append("<p>The distributions of log2-transformed means of quantitative protein values for each experimental group are plotted in the next graph.</p>")
        html_parts.append("<p>In a healthy experiment, the distributions of protein quantitative mean values of all groups are similar, depicted by the same shape of the distribution curves.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    # Prot - BoxInt
    if "log2_data" in st.session_state and "meta" in st.session_state:
        if st.session_state.get("mode6", "Mean") == "Mean":
            fig6 = boxplot_int(
                data=st.session_state.get("log2_data"),
                meta=st.session_state.get("meta"),
                outliers=st.session_state.get("outliers6", False),
                header=st.session_state.get("header6", True),
                width_cm=st.session_state.get("plotWidth6", 20),
                height_cm=st.session_state.get("plotHeight6", 10),
                dpi=st.session_state.get("plotDPI6", 300),
                plot_colors=st.session_state["selected_colors"]
            )
        elif st.session_state.get("mode6", "Mean") == "Single":
            fig6 = boxplot_int_single(
                data=st.session_state.get("log2_data"),
                meta=st.session_state.get("meta"),
                outliers=st.session_state.get("outliers6", False),
                id=st.session_state.get("id6", True),
                header=st.session_state.get("header6", True),
                legend=st.session_state.get("legend6", True),
                width_cm=st.session_state.get("plotWidth6", 20),
                height_cm=st.session_state.get("plotHeight6", 10),
                dpi=st.session_state.get("plotDPI6", 300),
                plot_colors=st.session_state["selected_colors"]
            )

        fig6.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig6.savefig(buf, format="png", dpi=st.session_state.get("plotDPI6", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.4 Boxplot Intensity - Protein Level</h2>")
        html_parts.append("<p>Protein intensity distributions - Depicted are the measured log2-raw protein intensities in a boxplot of each sample.</p>")
        html_parts.append("<p>Missing values are not considered when creating the boxplots. Zero intensity values are considered as missing values.</p>")
        html_parts.append("<p>For a healthy experiment, we expect similar log2-intensity distributions in all samples and all conditions.</p>")
        html_parts.append("<p>Significant deviations indicate unequal peptide loading during LC-MS/MS or unequal protein amounts in the starting material. Also, major contaminants in the samples can influence accurate peptide concentration measurements and, therefore, the number of protein identifications and protein quantitative values in a given sample.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - Cov
    if "org_data" in st.session_state and "meta" in st.session_state:
        fig7 = cov_plot(
            data=st.session_state.get("org_data"),
            meta=st.session_state.get("meta"),
            outliers=st.session_state.get("outliers7", False),
            header=st.session_state.get("header7", True),
            legend=st.session_state.get("legend7", True),
            width_cm=st.session_state.get("plotWidth7", 20),
            height_cm=st.session_state.get("plotHeight7", 10),
            dpi=st.session_state.get("plotDPI7", 300),
            plot_colors=st.session_state["selected_colors"]
        )

        fig7.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig7.savefig(buf, format="png", dpi=st.session_state.get("plotDPI7", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.5 Coefficient of Variation Plot - Protein Level</h2>")
        html_parts.append("<p>The CV or relative Standard Deviation is calculated by the ratio of the standard deviation to the mean (on quantitative protein level). It is used to measure the precision of a measure, in this case protein intensity.</p>")
        html_parts.append("<p>The plot below shows the distribution of the CVs by experimental conditions (groups) in boxplots where each CV is calculated by all proteins and by experimental condition. The CV is displayed as %CV, which is the percentage of dispersion of data points around the mean for all identified proteins in the data set per experimental group.</p>")
        html_parts.append("<p>For a healthy experiment, we expect: The distribution of the %CVs across all conditions to be mostly overlapping, and of similar modes. The modes of the %CVs not to be too high, ideally not above 50%, for best results in the statistical section not above 25%. If the distributions show worryingly large %CVs, the precision of LC-MS/MS measurements across all samples in the data set is low which could affect the quality of the differential expression analysis in the statistics section, due to high variance in the data set.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Protein - PCA
    if "log2_data" in st.session_state and "meta" in st.session_state:
        fig8 = pca_plot(
            data=st.session_state.get("log2_data"),
            meta=st.session_state.get("meta"),
            header=st.session_state.get("header8", True),
            legend=st.session_state.get("legend8", True),
            width_cm=st.session_state.get("plotWidth8", 20),
            height_cm=st.session_state.get("plotHeight8", 10),
            dpi=st.session_state.get("plotDPI8", 300),
            plot_colors=st.session_state["selected_colors"]
        )

        fig8.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig8.savefig(buf, format="png", dpi=st.session_state.get("plotDPI8", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.6 PCA Plot - Protein Level</h2>")
        html_parts.append("<p>The Principal Component Analysis (PCA) plot is used to visualize differences between samples that are constituted by their protein intensity profiles. PCA transforms high-dimensional data, like thousands of measured proteins or peptides intensities, into a reduced set of dimensions. The first two dimensions explain the greatest variability between the samples and they are a useful visual tool to confirm known clustering of the samples (of the same experimental group or replicates) or to identify potential problems in the data set.</p>")
        html_parts.append("<p>For a healthy experiment, we expect: 1.) Technical replicates to cluster tightly together. 2.) Biological replicates to cluster more closely together than non-replicates. 3. Clustering of samples of the same condition of interest (or experimental group) should be visible and different experimental groups should be separated from each other.</p>")
        html_parts.append("<p>If unexpected clusters occur or replicates do not cluster together, it can be due to high individual sample variability, or extra variability introduced by factors such as technical processing, high variability of protein content in the starting material (e.g. FBS from the culture medium still present), other unexplored biological differences, or unintentional sample swaps, etc).</p>")
        html_parts.append("<p>The interpretation and trust in the differential protein expression results in the statistical section should take the above mentioned considerations into account. If you think that the samples in the experiment show largely unexpected patterns, it is advisable to request support from a data analyst to help with data interpretation.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    #Prot - Abundance
    if "org_data" in st.session_state and "meta" in st.session_state:
        fig9 = abundance_plot(
            data=st.session_state.get("org_data"),
            meta=st.session_state.get("meta"),
            workflow="Protein",
            width_cm=st.session_state.get("plotWidth9", 20),
            height_cm=st.session_state.get("plotWidth9", 10),
            dpi=st.session_state.get("plotDPI9", 300),
            legend=st.session_state.get("legend9", True),
            plot_colors=st.session_state["selected_colors"]
        )

        fig9.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig9.savefig(buf, format="png", dpi=st.session_state.get("plotDPI9", 300), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.7 Abundance Plot - Protein Level</h2>")
        html_parts.append("<p>The rank plot of mean protein intensities (log10) within experimental conditions / groups illustrates the coverage of the dynamic ranges of protein intensities within the experimental conditions / groups in the data set.</p>")
        html_parts.append("<p>For a healthy experiment, we expect high similarity between the curves.</p>")
        html_parts.append("<p>If the curves look worryingly different, it could affect the quality of the differential protein expression analysis in the statistics section.</p>")
        html_parts.append("<p>Usually, the covered dynamic range of protein quantitative values should be four or more orders of magnitude. Though, the observable dynamic range is also dependent on the biological sample analyzed (e.g. body fluid or cell lysate). for example the vast dynamic range of plasma or serum proteins, estimated at 12-13 orders of magnitude, renders plasma proteome analysis based on mass spectrometry extremely challenging, as the high abundant plasma proteins hinder the identification of low abundant proteins and negatively impact the accessible dynamic range.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    if "org_data" in st.session_state and "meta" in st.session_state:
        fig12 = corr_plot(
            data=st.session_state.get("org_data"),
            meta=st.session_state.get("meta"),
            method=st.session_state.get("method12", "Matrix"),
            id=st.session_state.get("legend12",False),
            full_range=False,
            width=st.session_state.get("plotWidth12", 10),
            height=st.session_state.get("plotHeight12", 8),
            dpi=st.session_state.get("plotDPI12", 100)
        )

        fig12.set_size_inches(12, 6)
        buf = io.BytesIO()
        fig12.savefig(buf, format="png", dpi=st.session_state.get("plotDPI12", 100), bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode("utf-8")

        html_parts.append(f"<h2>{prot_num}.8 Correlation Plot - Protein Level</h2>")
        html_parts.append("<p>The correlation plot shows the Pearson’s correlation coefficient between all samples in the experiment. Hierarchical clustering was used to order the samples in the matrix.</p>")
        html_parts.append("<p>For a healthy experiment we expect: a) technical replicates have high correlations (1 = 100% identical), b) biological replicates (samples of the same group) to have higher correlations than non-replicates.</p>")
        html_parts.append("<p>Low correlation between samples may indicate many missing values in a given sample or low similarity between different experimental conditions or groups. Low correlation between samples can also be related to normalization errors due to unequal peptide amount loading during LC-MS/MS measurements or major contaminants in the samples (please, compare with plots 2.3 Histogram Intensity and 2.4 Boxplot Intensity).</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    html_parts.append("</body></html>")
    return "\n".join(html_parts)