import pandas as pd
from collections import Counter
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

########################################################################################################
###################### Definitions:

def extract_mir_words(text):
    processed_words = set()
    for word in text:
        if (word.startswith('mir') or word.startswith('hsa-mir')) and len(word) > 4 and word not in ['mirna', 'mirnas']:
            word_clean = word.replace('–', '-').replace('—', '-')
            suffix = re.split(r'(?:mir-|mirna-|mir\.)', word_clean)[-1]
            merged_prefix_mir = f"miR-{suffix}"
            processed_words.add(merged_prefix_mir)
    return processed_words

def miRNA_ID_validation(input_set, known_miRNA_IDs):
    validated_ID = set()
    for miRNA in input_set:
        found = False
        temp = miRNA
        while len(temp) >= 4:
            if temp in known_miRNA_IDs:
                validated_ID.add(temp)
                found = True
                break
            temp = temp[:-1]
        if found:
            continue
    return validated_ID

def miRNA_ID_validation_immature(input_set, known_miRNA_IDs):
    validated_ID = set()
    for miRNA in input_set:
        found = False
        if '-5p' in miRNA:
            temp = miRNA.split('-5p')[0]
        elif '-3p' in miRNA:
            temp = miRNA.split('-3p')[0]
        else:
            temp = miRNA
        while len(temp) >= 4:
            if temp in known_miRNA_IDs:
                validated_ID.add(temp)
                found = True
                break
            temp = temp[:-1]
        if found:
            continue
    return validated_ID

def transform_row(row):
    first_positive_index = next((i for i, value in enumerate(row) if value > 0), None)
    if first_positive_index is None:
        return row
    transformed_row = [
        np.nan if i < first_positive_index else value
        for i, value in enumerate(row)
    ]
    return transformed_row

##############################################################################################################
## Load miRNA databases

miRNA_DAT = open('00_database_miRna/merged_miRNA.dat', 'r')
miRNA_mature = open('00_database_miRna/merged_mature.fa', 'r')

all_miRNA_IDs_deduplicated = set()

for x in miRNA_DAT:
    if 'DE' in x.split(' ')[0]:
        id_part = x.strip().split(' ')[5]
        clean_id = id_part.replace('–', '-').replace('—', '-')
        all_miRNA_IDs_deduplicated.add(clean_id)

for x in miRNA_mature:
    id_part = x.strip().split(' ')[-1]
    clean_id = id_part.replace('–', '-').replace('—', '-')
    all_miRNA_IDs_deduplicated.add(clean_id)

known_miRNA_IDs = sorted(list(all_miRNA_IDs_deduplicated))

###############################################################################################################
## Process SCOPUS dataset

df = pd.read_csv("FORMATTED_FOR_BIBLIOMETRIX_DATASET_relevant_miRNA_CRC_Review.csv")
small_df = df[["Author Keywords", "Year", "Title", "Abstract"]].copy()

small_df[["Author Keywords", "Title", "Abstract"]] = small_df[["Author Keywords", "Title", "Abstract"]].apply(lambda x: x.str.lower().fillna(""))
small_df["Author Keywords"] = small_df["Author Keywords"].str.split("; ")
small_df["Title"] = small_df["Title"].str.split(" ")
small_df["Abstract"] = small_df["Abstract"].str.split(" ")

small_df["Mir_Words"] = small_df.apply(lambda row: row["Author Keywords"] + row["Title"] + row["Abstract"], axis=1)
small_df["Mir_Words"] = small_df["Mir_Words"].apply(extract_mir_words)
small_df['validated_Mir_IDs'] = small_df['Mir_Words'].apply(lambda x: miRNA_ID_validation(x, known_miRNA_IDs))
small_df["Mir_Words_Count"] = small_df["validated_Mir_IDs"].apply(len)

with pd.ExcelWriter("processed_data.xlsx", engine="openpyxl") as writer:
    small_df.to_excel(writer, sheet_name='Dataset', index=False)

melted_df = small_df.explode("validated_Mir_IDs").dropna(subset=["validated_Mir_IDs"])
yearly_counts = melted_df.groupby(["Year", "validated_Mir_IDs"]).size().unstack(fill_value=0)
yearly_counts_T = yearly_counts.T.sort_index(axis=1)
yearly_counts_T.columns = yearly_counts_T.columns.astype(str)

current_years = yearly_counts_T.columns
numeric_years = [int(year) for year in current_years]
min_year = min(numeric_years)
max_year = max(numeric_years)
all_years = list(range(min_year, max_year + 1))
existing_years = set(numeric_years)
missing_years = [year for year in all_years if year not in existing_years]

for year in missing_years:
    yearly_counts_T[str(year)] = 0

sorted_columns = [str(year) for year in all_years]
yearly_counts_T = yearly_counts_T[sorted_columns]

transformed_rows = [transform_row(row) for _, row in yearly_counts_T.iterrows()]
df_transformed = pd.DataFrame(transformed_rows, columns=sorted_columns, index=yearly_counts_T.index)

annual_columns = [col for col in df_transformed.columns if col.isdigit()]
df_transformed['Total_Mentions'] = df_transformed[annual_columns].sum(axis=1)
df_transformed['Average_Mentions'] = df_transformed[annual_columns].mean(axis=1).round(2)

sorted_df = df_transformed.sort_values('Total_Mentions', ascending=False)
sorted_df.to_csv("miRNA_NAME_per_year.csv", index=True)

top_n = 10
top_mirnas = sorted_df.head(top_n).index
top_data = sorted_df.loc[top_mirnas]

heatmap_data = pd.concat([
    top_data[annual_columns],
    top_data[['Total_Mentions']],
    top_data[['Average_Mentions']]
], axis=1)

plt.figure(figsize=(14, 8))

annual_min = top_data[annual_columns].min().min()
annual_max = top_data[annual_columns].max().max()

# Fix FutureWarning by using column names directly
heatmap_for_plot = heatmap_data.copy()
heatmap_for_plot[heatmap_for_plot.columns[-2:]] = 0  # Prevent color scale interference

# Format annotations: integers for annual + total, 2 decimals for average
annot_data = heatmap_data.copy().astype(object)
for col in annual_columns + ['Total_Mentions']:
    annot_data[col] = annot_data[col].fillna(0).astype(int).astype(str)
annot_data['Average_Mentions'] = annot_data['Average_Mentions'].apply(lambda x: f"{x:.2f}")

# Plot heatmap
ax = sns.heatmap(
    heatmap_for_plot,
    cmap="YlOrRd",
    annot=annot_data.values,
    fmt='s',  # Use string formatting for per-cell annotation
    linewidths=0.5,
    cbar_kws={"label": "Articles (Annual Data)"},
    vmin=annual_min,
    vmax=annual_max
)

# Draw border around Total & Average columns
for col_idx in range(len(annual_columns), len(heatmap_data.columns)):
    ax.add_patch(plt.Rectangle((col_idx, 0), 1, len(top_data),
                               fill=False, edgecolor='grey', lw=2))

# Set x-axis labels
column_labels = (
    [f"{col}" for col in annual_columns] +
    ["Total Articles", "Average Articles"]
)
ax.set_xticks(np.arange(len(column_labels)))
ax.set_xticklabels(column_labels, rotation=45, ha="center")

# Separator lines
num_annual = len(annual_columns)
ax.axvline(x=num_annual, color='black', linewidth=2)
ax.axvline(x=num_annual + 1, color='black', linewidth=2)

plt.xlabel("Years and Metrics", fontsize=12)
plt.ylabel("miRNA", fontsize=12)
plt.tight_layout()
plt.show()
