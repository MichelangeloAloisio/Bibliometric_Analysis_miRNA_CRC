1) Scopus database download: using the following query ": ( TITLE ( colorectal ) AND TITLE ( mir-* ) AND KEY ( colorectal AND cancer ) AND KEY ( mir-* ) ) AND ( LIMIT-TO ( DOCTYPE , "ar" ) )", the dataset was downloaded in 2 formats (.csv, .ris). The scopus output files are available:  00_CRC_miRNA_scopus.csv and 00_CRC_miRNA_scopus.ris.
2) The Manually screened dataset is available: 01_Screened_dataset_miRNA_CRC_Review.csv
3)Python Script: 02_extract_top_10_miRNAs_per_year.py This script performs two core functions: A) miRNA ID Extraction: Parses titles and author-provided keywords of each article to extract miRNA identifiers (e.g., miR-21, miR-155). Outputs a cleaned dataset with miRNA IDs linked to their respective articles. B) Publications per Year Analysis Generates a "Publications per Year" plot integrating: Annual publication counts. The top 10 most-studied miRNAs per year (based on frequency in titles/keywords). The plot visualizes both metrics simultaneously, enabling analysis of miRNA research trends alongside overall field activity. Input : 00_CRC_miRNA_scopus.csv (raw Scopus data). 02_processed_data.xlsx.
Publication trends plot (saved as publications_per_year.png).
4)
5)
6)
7) The pyhton script to  elaborate the number of publication per year is available (01_extract_mirnaID_from_scopus_dataset.py). Specifically the script is able to extract from Title and Authors keyword mirna ID and it is also useful to generate the piblication per year plot.
8) The  processed_data.xlsx produced by the python script is available in the following folder.
9) 



