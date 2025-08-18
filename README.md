# Proteomics CoPYlot (Python Version)

Proteomics CoPYlot is a **Python-based tool for proteomics data analysis and visualization**, designed to simplify common workflows and provide interactive plots and light statistics for your experimental results. It can be run **locally** as a Python executable or accessed via the **web** at [proteomics-data.com](https://proteomics-data.com).

---

## Features

- Upload and analyze **Spectronaut, DIA-NN, or MaxQuant output files**  
- Generate **quality control (QC) plots**: retention times, modifications, missed cleavages, and more  
- **Interactive visualizations** with adjustable plot size, resolution, and headers  
- Export **publication-ready figures** directly from the interface  
- Lightweight **statistical summaries** for quick assessment of datasets  

---

## Installation

```bash
# Clone the repository
git clone https://github.com/JonasMarx3007/ProteomisCoPYlot.git

# Navigate to the project folder
cd ProteomisCoPYlot

# Install dependencies
pip install -r requirements.txt
```

## Example Workflow

1. Upload your proteomics dataset (TSV or compatible format)  
2. Explore proteins and phosphosites interactively in the **Streamlit interface**  
3. Generate **quality control plots** with customizable size, resolution, and headers  
4. Export figures for reports or publications  

---

## Dependencies

Python ≥ 3.8  

Required Python packages:

- streamlit  
- matplotlib  
- pandas  
- numpy  
- scikit-learn  
- statsmodels  
- plotly  
- gprofiler-official  

Install via pip:

```bash
pip install streamlit matplotlib pandas numpy scikit-learn statsmodels plotly gprofiler-official
```

## Usage
Run the locally using Streamlit:
```bash
streamlit run app.py
```
Open http://localhost:8501 in your browser to interact with the application. Or alternatively run the run_app.py.

## License
This project is licensed under the **MIT License** — see the LICENSE file for details.