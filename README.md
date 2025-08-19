Exploring miRNA Research in Colorectal Cancer: Insights from a Bibliometric Analysis

This repository contains the code and data for the study:

"Exploring miRNA Research in Colorectal Cancer: Insights by A Bibliometric Analysis"

It presents a bibliometric analysis of microRNA (miRNA) research in colorectal cancer, based on scientific publications. The workflow includes data screening, miRNA extraction and validation, trend analysis, geographical mapping, and author-level focus assessment.

Scopus Query

(TITLE (colorectal) AND TITLE (mir-*) AND KEY (colorectal AND cancer) AND KEY (mir-*)) AND (LIMIT-TO (DOCTYPE, "ar")) 
(accessed on 25 September 2024)

This query retrieves scientific articles where:
- 'colorectal' and a miRNA term (e.g., mir-21) are in the title,
- 'colorectal', 'cancer', and a miRNA term are in the keywords.

You can reproduce the initial dataset by running this query on Scopus.

Files in This Repository

- 00_scopus_1385.csv
  Scopus downloaded database in data (25 September 2024) 

- 01_Screened_dataset_miRNA_CRC_Review.csv  
  Manually reviewed and curated dataset. Includes only articles selected after title and abstract screening. Used as input for all downstream analyses.  
  This file contains publications from Scopus that were manually screened and selected based on title and abstract. It is derived from scholarly metadata used under permitted access. Full-text content is not included. The dataset and all code are open and shared to support transparency and reproducibility.

- 02_extract_miRNAs_per_year.py  
  Python script that:
  - Extracts miRNA identifiers (e.g., miR-21) from titles, abstracts, and keywords.
  - Validates them against miRBase (requires local miRBase files: merged_miRNA.dat and merged_mature.fa).
  - Generates yearly publication trends.
  - Outputs a ranked list of the most studied miRNAs and creates a heatmap of the top 10 miRNAs over time.

- 03_r_world_map.R  
  R script that visualizes the geographical distribution of research contributions.  
  Takes country-level publication data (derived from author affiliations) and generates a world map using the rworldmap package.

- 04_extract_mirna_from_10_most_productive_authors.py  
  Python script that:
  - Identifies the 10 most productive authors in the field.
  - Extracts and validates miRNAs from their article titles.
  - Builds a matrix of miRNA occurrences per author.
  - Generates a heatmap showing shared and specific research interests.

Requirements

Python (3.8+): pandas, openpyxl, matplotlib, seaborn, numpy  
R: rworldmap, classInt, RColorBrewer

Install Python packages:
pip install pandas openpyxl matplotlib seaborn numpy

Citation

If you use this work, please cite:

"Exploring miRNA Research in Colorectal Cancer: Insights by A Bibliometric Analysis"
