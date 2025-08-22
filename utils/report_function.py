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

        html_parts.append(f"<h2>1.1 Coverage Plot</h2>")
        html_parts.append("<p>This plot coverage.</p>")
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

        html_parts.append(f"<h2>1.2 Missing Value Plot</h2>")
        html_parts.append("<p>This plot show missing values.</p>")
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

        html_parts.append(f"<h2>1.3 RT Plot</h2>")
        html_parts.append("<p>This plot something RT.</p>")
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

        html_parts.append(f"<h2>1.4 Mod Plot</h2>")
        html_parts.append("<p>This plot something Mods.</p>")
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

        html_parts.append(f"<h2>1.5 Missed CL Plot</h2>")
        html_parts.append("<p>This plot something Misses.</p>")
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

        html_parts.append("<h1>1 Protein Level</h1>")
        html_parts.append("<h2>1.1 Coverage Plot</h2>")
        html_parts.append("<p>This is some filler text describing the coverage plot.</p>")
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

            html_parts.append(f"<h2>1.2 Missing Value Plot</h2>")
            html_parts.append("<p>This plot shows the distribution of missing values.</p>")
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

        html_parts.append("<h2>1.3 Histogram Intensity Plot</h2>")
        html_parts.append("<p>Something about Histogram Intensity.</p>")
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

        html_parts.append("<h2>1.4 Boxplot Intensity Plot</h2>")
        html_parts.append("<p>Something about Boxplot Intensity.</p>")
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

        html_parts.append("<h2>1.5 Cov Plot</h2>")
        html_parts.append("<p>Something about Cov.</p>")
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

        html_parts.append("<h2>1.6 PCA Plot</h2>")
        html_parts.append("<p>Something about PCA.</p>")
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

        html_parts.append("<h2>1.7 Abundance Plot</h2>")
        html_parts.append("<p>Something about Abundance.</p>")
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

        html_parts.append("<h2>1.8 Heatmap</h2>")
        html_parts.append("<p>Something about a Heatmap.</p>")
        html_parts.append(f"<img src='data:image/png;base64,{img_base64}' style='width:100%;height:auto;'>")

    html_parts.append("</body></html>")
    return "\n".join(html_parts)