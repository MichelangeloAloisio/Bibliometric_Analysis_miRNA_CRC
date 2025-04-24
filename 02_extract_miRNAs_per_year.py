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
            # Normalize hyphens
            word_clean = word.replace('–', '-').replace('—', '-')
            suffix = re.split(r'(?:mir-|mirna-|mir\.)', word_clean)[-1]
            merged_prefix_mir = f"miR-{suffix}"
            processed_words.add(merged_prefix_mir)
    return processed_words

def miRNA_ID_validation_considering_mature(input_set, Tolta_microRNA_ID):
    validated_ID = set()
    for miRNA in input_set:
        found = False
        temp = miRNA
        while len(temp) >= 4:
            if temp in Tolta_microRNA_ID:
                validated_ID.add(temp)
                found = True
                break
            temp = temp[:-1]  # Delete last character
        if found:
            continue
    return validated_ID


def miRNA_ID_validation_immature(input_set, Tolta_microRNA_ID):
    validated_ID = set()
    for miRNA in input_set:
        found = False
        # Extract the prefix before '-5p' or '-3p', or use the full miRNA if neither is present
        if '-5p' in miRNA:
            temp = miRNA.split('-5p')[0]
        elif '-3p' in miRNA:
            temp = miRNA.split('-3p')[0]
        else:
            temp = miRNA  # Handle miRNAs without '-5p' or '-3p' suffix
        
        # Check progressively shorter substrings of the prefix
        while len(temp) >= 4:
            if temp in Tolta_microRNA_ID:
                validated_ID.add(temp)
                found = True
                break
            temp = temp[:-1]  # Remove the last character and try again
        
        if found:
            continue  # Skip to next miRNA if a valid ID was found
    return validated_ID



def transform_row(row):
    # Find the index of the first positive value
    first_positive_index = next((i for i, value in enumerate(row) if value > 0), None)
    
    # If no positive values are found, return the original row
    if first_positive_index is None:
        return row
    
    # Transform all previous elements into NaN
    transformed_row = [
        np.nan if i < first_positive_index else value
        for i, value in enumerate(row)
    ]
    
    return transformed_row

##############################################################################################################
## Open miRNA databases


miRNA_DAT = open('00_database_miRna/merged_miRNA.dat', 'r')
miRNA_mature = open('00_database_miRna/merged_mature.fa', 'r')

Toltal_microRNA_ID_deduplication = set()

# Extract from miRNA_DAT (normalize hyphens)
for x in miRNA_DAT:
    if 'DE' in x.split(' ')[0] and 'homo' in x.lower():
        id_part = x.strip().split(' ')[5]
        clean_id = id_part.replace('–', '-').replace('—', '-')
        Toltal_microRNA_ID_deduplication.add(clean_id)

# Extract from mature.fa (normalize hyphens)
for x in miRNA_mature:
    if 'homo' in x.lower():
        id_part = x.strip().split(' ')[-1]
        clean_id = id_part.replace('–', '-').replace('—', '-')
        Toltal_microRNA_ID_deduplication.add(clean_id)

Tolta_microRNA_ID = sorted(list(Toltal_microRNA_ID_deduplication))

###############################################################################################################
## Elaboration of SCOPUS dataset
df = pd.read_csv("FORMATTED_FOR_BIBLIOMETRIX_DATASET_relevant_miRNA_CRC_Review.csv")

small_df = df[["Author Keywords", "Year", "Title", "Abstract"]].copy()





# Convert everything to lowercase and clean empty values
small_df[["Author Keywords", "Title", "Abstract"]] = small_df[["Author Keywords", "Title", "Abstract"]].apply(lambda x: x.str.lower().fillna(""))

# Split into lists
small_df["Author Keywords"] = small_df["Author Keywords"].str.split("; ")
small_df["Title"] = small_df["Title"].str.split(" ")
small_df["Abstract"] = small_df["Abstract"].str.split(" ")


# Merge lists to extract miRNA words
small_df["Mir_Words"] = small_df.apply(lambda row: row["Author Keywords"] + row["Title"] + row["Abstract"], axis=1)

# Extract miRNA terms (normalize hyphens)
small_df["Mir_Words"] = small_df["Mir_Words"].apply(extract_mir_words)

# Validate with the database for immature mirna, this means withaout -5p or -3p at the end of the miRNA ID

small_df['validated_Mir_IDs'] = small_df['Mir_Words'].apply(
    lambda x: miRNA_ID_validation_considering_mature(x, Tolta_microRNA_ID)
)




small_df["Mir_Words_Count"] = small_df["validated_Mir_IDs"].apply(len)


#### Save Excel file for processed data

with pd.ExcelWriter("processed_data.xlsx", engine="openpyxl") as writer:
    small_df.to_excel(writer, sheet_name='Dataset', index=False)




melted_df = small_df.explode("validated_Mir_IDs").dropna(subset=["validated_Mir_IDs"])
yearly_counts = melted_df.groupby(["Year", "validated_Mir_IDs"]).size().unstack(fill_value=0)
yearly_counts_T = yearly_counts.T.sort_index(axis=1)


# Ensure column names are strings
yearly_counts_T.columns = yearly_counts_T.columns.astype(str)

# Extract column names (years)
current_years = yearly_counts_T.columns

#  Convert years to integers
numeric_years = [int(year) for year in current_years]

#  Calculate minimum and maximum years
min_year = min(numeric_years)
max_year = max(numeric_years)

#  Generate all years between min_year and max_year
all_years = list(range(min_year, max_year + 1))

# Identify missing years
existing_years = set(numeric_years)
missing_years = [year for year in all_years if year not in existing_years]

#  Add columns for missing years with value 0
for year in missing_years:
    yearly_counts_T[str(year)] = 0  # Add column with 0

# Sort columns by year in ascending order
sorted_columns = [str(year) for year in all_years]
yearly_counts_T = yearly_counts_T[sorted_columns]


# ---------------------------
# Apply the function to all rows
# ---------------------------
transformed_rows = [transform_row(row) for _, row in yearly_counts_T.iterrows()]

# ---------------------------
# Reconstruct the transformed DataFrame using sorted columns
# ---------------------------

df_transformed = pd.DataFrame(
    transformed_rows,
    columns=sorted_columns,  # Use sorted columns
    index=yearly_counts_T.index  # Use the original index
)


# ---------------------------
# CORRECT CALCULATION OF METRICS
# ---------------------------
annual_columns = [col for col in df_transformed.columns if col.isdigit()]
df_transformed['Total_Mentions'] = df_transformed[annual_columns].sum(axis=1)
df_transformed['Average_Mentions'] = df_transformed[annual_columns].mean(axis=1).round(2)  # Arrotondamento a 2 decimali

# ---------------------------
# ORDERING AND SELECTION OF THE TOP 10
# ---------------------------

sorted_df = df_transformed.sort_values('Total_Mentions', ascending=False)
sorted_df.to_csv("miRNA_NAME_per_year.csv", index=True)

top_n = 10
top_mirnas = sorted_df.head(top_n).index
top_data = sorted_df.loc[top_mirnas]


# ---------------------------
# CREATING THE DATAFRAME FOR THE HEATMAP
# ---------------------------

heatmap_data = pd.concat([
    top_data[annual_columns],
    top_data[['Total_Mentions']],
    top_data[['Average_Mentions']]
], axis=1)

# ---------------------------
# PLOT: HEATMAP WITH REQUESTED MODIFICATIONS
# ---------------------------

plt.figure(figsize=(14, 8))

# Calculate colormap limits only on the years
annual_min = top_data[annual_columns].min().min()
annual_max = top_data[annual_columns].max().max()

# Create a copy for the plot (last columns do not affect the colormap)
heatmap_for_plot = heatmap_data.copy()
heatmap_for_plot.iloc[:, -2:] = 0  # Imposto a 0 per non influenzare il colormap

ax = sns.heatmap(
    heatmap_for_plot,
    cmap="YlOrRd",
    annot=heatmap_data.values,  # Usa i valori originali per le annotazioni
    fmt='.2f',  # Formato con due decimali per la media
    linewidths=0.5,
    cbar_kws={"label": "Articles (Annual Data)"},
    vmin=annual_min,
    vmax=annual_max
)

# ---------------------------
# VISUAL MODIFICATIONS
# ---------------------------


for col_idx in range(len(annual_columns), len(heatmap_data.columns)):
    ax.add_patch(plt.Rectangle((col_idx, 0), 1, len(top_data), 
                              fill=False, edgecolor='grey', lw=2))

# 2. Center column labels
columns_labels = (
    [f"{col}" for col in annual_columns] +
    ["Total Articles", "Average Articles"]
)
ax.set_xticks(np.arange(len(columns_labels)))
ax.set_xticklabels(columns_labels, rotation=45, ha="center")

# 3. Vertical lines to separate sections
num_annual = len(annual_columns)
ax.axvline(x=num_annual, color='black', linewidth=2)  # Dopo gli anni
ax.axvline(x=num_annual + 1, color='black', linewidth=2)  # Dopo la "Total"

# ---------------------------
# FINAL FORMATTING
# ---------------------------
plt.xlabel("Years and Metrics", fontsize=12)
plt.ylabel("miRNA", fontsize=12)
plt.tight_layout()
plt.show()
