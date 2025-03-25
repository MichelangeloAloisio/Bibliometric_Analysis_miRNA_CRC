
import pandas as pd
from collections import Counter
import re
import matplotlib.pyplot as plt
import seaborn as sns



# Read the CSV file
df = pd.read_csv("FORMATTED_FOR_BIBLIOMETRIX_DATASET_relevant_miRNA_CRC_Review.csv")
print(df[["Year"]])




# Create a smaller DataFrame with relevant columns
small_df = df[["Author Keywords", "Year", "Title"]].copy()
small_df[["Author Keywords", "Title"]] = small_df[["Author Keywords", "Title"]].apply(lambda x: x.str.lower().fillna(""))

# Split keywords and titles into lists
small_df["Author Keywords"] = small_df["Author Keywords"].str.split("; ")
small_df["Title"] = small_df["Title"].str.split(" ")


def process_suffix(x):
    """
    Process miRNA suffixes to standardize naming:
    - Retains -5p and -3p suffixes
    - Removes other suffixes like -m, -wnt, -d, etc.
    """
    if "-m" in x:
        if x.startswith('hsa-mir'):
            if "-5p" in x:
                base = x.split('hsa-mir')[1].split('-5p')[0]
                return f"mir{base}-5p"
            elif "-3p" in x:
                base = x.split('hsa-mir')[1].split('-3p')[0]
                return f"mir{base}-3p"
            elif "-wnt" in x:
                base = x.split('hsa-mir')[1].split('-wnt')[0]
                return f"mir{base}"
            elif "-d" in x:
                base = x.split('hsa-mir')[1].split('-d')[0]
                return f"mir{base}"
            elif "," in x:
                base = x.split('hsa-mir')[1].split(',')[0]
                return f"mir{base}"
            elif " fam" in x:
                base = x.split('hsa-mir')[1].split(',')[0]
                return f"mir{base}"
            elif "*" in x:
                base= x.split('hsa-mir')[1].split('*')[0]
                return f"mir{base}"
            else:
                return f"mir{x.split('hsa-mir')[1]}"

        else:
            if "-5p" in x:
                return x.split("-5p")[0] + "-5p"
            elif "-3p" in x:
                return x.split("-3p")[0] + "-3p"
            elif "-wnt" in x:
                return x.split("-wnt")[0]
            elif "-d" in x:
                return x.split("-d")[0]
            elif "," in x:
                return x.split(",")[0]
            elif "*" in x: 
                return x.split("*")[0]
            elif " fam" in x:
                return x.split(" fam")[0]
            else:
                return x
    else:
        if "-5p" in x:
            return x.split("-5p")[0] + "-5p"
        elif "-3p" in x:
            return x.split("-3p")[0] + "-3p"
        elif "-wnt" in x:
            return x.split("-wnt")[0]
        elif "-d" in x:
            return x.split("-d")[0]
        elif "," in x:
            return x.split(",")[0]
        elif "*" in x:
            return x.split("*")[0]
        elif " fam" in x:
            return x.split(" fam")[0]
        else:
            return x

def extract_mir_words(text):
    """
    Extract standardized miRNA terms from text:
    - Terms must start with 'mir' or 'hsa-mir'
    - Exclude 'mirna' and 'mirnas'
    - Handle suffix processing
    """
    processed_words = set()
    for x in text:
        if (x.startswith('mir') or x.startswith('hsa-mir')) and len(x) > 4 and x not in ['mirna', 'mirnas']:
            if "/" in x:
                processed_words.add(process_suffix(x.split("/")[0]))
            else:
                processed_words.add(process_suffix(x))
    return list(processed_words)

# Combine Author Keywords and Title for miRNA extraction
small_df["Mir_Words"] = small_df.apply(
    lambda row: row["Author Keywords"] + row["Title"], axis=1
)

# Extract miRNA terms
small_df["Mir_Words"] = small_df["Mir_Words"].apply(extract_mir_words)

# Add count column
small_df["Mir_Words_Count"] = small_df["Mir_Words"].apply(len)

# Save processed data to Excel (requires openpyxl)



with pd.ExcelWriter("processed_data.xlsx", engine="openpyxl") as writer:
    small_df.to_excel(writer, sheet_name='Dataset', index=False)
    # Create annual miRNA count matrix
    melted_df = small_df.explode("Mir_Words").dropna(subset=["Mir_Words"])
    yearly_counts = melted_df.groupby(["Year", "Mir_Words"]).size().unstack(fill_value=0)
    yearly_counts_T = yearly_counts.T.sort_index(axis=1)
    yearly_counts_T.to_excel(writer, sheet_name='Annual_Counts', index=True)

####################################################################



# Assuming you have a DataFrame like 'small_df' with:
# - 'Year': publication year
# - 'Mir_Words': list of miRNAs in the paper

# 1. Calculate the number of PAPERS per year (not miRNA mentions)
total_papers = small_df["Year"].value_counts().sort_index()

# 2. Create the miRNA paper count matrix (number of papers per miRNA/year)
# Explode Mir_Words and count unique papers per miRNA/year
melted_df = (
    small_df
    .assign(Paper_ID=small_df.index)  # Create a unique paper identifier
    .explode("Mir_Words")
    .dropna(subset=["Mir_Words"])
)

# Count unique papers per miRNA/year
yearly_paper_counts = (
    melted_df
    .groupby(["Year", "Mir_Words"])["Paper_ID"]
    .nunique()  # Count distinct papers
    .unstack(fill_value=0)
    .T  # Transpose for miRNA rows and years columns
    .sort_index(axis=1)  # Sort years chronologically
)

# 3. Plot setup
plt.figure(figsize=(20, 8))
sns.set(style="whitegrid")

# Plot total papers per year
total_papers.plot(
    kind="line",
    marker="o",
    color="b",
    linewidth=2,
    markersize=8,
    label="Total Papers",
)

plt.xlabel("Year", fontsize=14)
plt.ylabel("Number of Publications", fontsize=14)
plt.xticks(rotation=45)

# Identify duplicate miRNAs (present in multiple years)
all_mirnas = []
for year in total_papers.index:
    top5 = yearly_paper_counts[year].nlargest(10)
    all_mirnas.extend(top5.index)
mirna_counts = Counter(all_mirnas)
duplicate_mirnas = {mirna for mirna, count in mirna_counts.items() if count > 1}

# Adjust y-axis to accommodate annotations
max_count = total_papers.max()
plt.ylim(top=max_count + 30)  # Add buffer space

# Add top 5 miRNAs annotations for each year with symbols
for year in total_papers.index:
    top5 = yearly_paper_counts[year].nlargest(10)
    top5 = top5[top5 > 0]
    
    if not top5.empty:
        formatted_mirnas = []
        for mirna in top5.index:
            symbol = "★" if mirna in duplicate_mirnas else ""
            formatted_mirnas.append(f"{symbol}{mirna} ")
        
        mirna_list = "\n".join(formatted_mirnas)
        
        plt.annotate(
            mirna_list,
            xy=(year, total_papers[year]),
            xytext=(year, total_papers[year] + 10),
            fontsize=9,
            color="black",
            bbox=dict(
                boxstyle="round,pad=0.5",
                fc="white",
                ec="black",
                lw=1.5
            ),
            ha='center',
            va='bottom'
        )

# Highlight peak year
peak_year = total_papers.idxmax()
peak_count = total_papers.max()
plt.annotate(
    f"      {peak_count} papers",
    xy=(peak_year, peak_count),
    xytext=(peak_year, peak_count - 15),
    fontsize=12,
    color="red",
    arrowprops=dict(facecolor="red", shrink=0.05),
    ha='center',
    va='top'
)

plt.tight_layout()
plt.savefig("temporal_distribution.tiff", dpi=600, bbox_inches='tight')
plt.show()

exit()



########### QUESTO FA LE CONTE SU TOTALE MIRNA 
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# 1. Temporal Distribution Plot
plt.figure(figsize=(20, 8))
sns.set(style="whitegrid")

total_papers = yearly_counts_T.sum(axis=0)

total_papers.plot(kind="line", marker="o", color="b", linewidth=2, markersize=8, label="Total Papers")

#plt.title("Temporal Distribution of miRNA and CRC Related Publications", fontsize=18, pad=20)
plt.xlabel("Year", fontsize=14)
plt.ylabel("Number of Publications", fontsize=14)
plt.xticks(rotation=45)

# Identify duplicate miRNAs (present in multiple years)
all_mirnas = []
for year in total_papers.index:
    top5 = yearly_counts_T[year].nlargest(10)
    all_mirnas.extend(top5.index)
mirna_counts = Counter(all_mirnas)
duplicate_mirnas = {mirna for mirna, count in mirna_counts.items() if count > 1}

# Adjust y-axis to accommodate annotations
max_count = total_papers.max()
plt.ylim(top=max_count + 30)  # Add buffer space

# Add top 5 miRNAs annotations for each year with symbols
for year in total_papers.index:
    top5 = yearly_counts_T[year].nlargest(10)
    top5 = top5[top5 > 0]
    
    if not top5.empty:
        # Create formatted list with symbols
        formatted_mirnas = []
        for mirna in top5.index:
            symbol = "★" if mirna in duplicate_mirnas else ""
            formatted_mirnas.append(f"{symbol}{mirna} ")
        
        mirna_list = "\n".join(formatted_mirnas)
        
        plt.annotate(mirna_list, 
                     xy=(year, total_papers[year]),
                     xytext=(year, total_papers[year] + 10),
                     fontsize=9, color="black",
                     bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=1.5),
                     ha='center', va='bottom')

# Highlight peak year
peak_year = total_papers.idxmax()
peak_count = total_papers.max()
plt.annotate(f"      {peak_count} papers", 
             xy=(peak_year, peak_count),
             xytext=(peak_year, peak_count - 15),
             fontsize=12, color="red",
             arrowprops=dict(facecolor="red", shrink=0.05),
             ha='center', va='top')

plt.tight_layout()
plt.savefig("temporal_distribution.tiff", dpi=600, bbox_inches='tight')
plt.show()




exit()



##################àà QUESTO FUNZIONA MA NON FA GRASSETTO
# 1. Temporal Distribution Plot
plt.figure(figsize=(16, 8))
sns.set(style="whitegrid")

total_papers = yearly_counts_T.sum(axis=0)
total_papers.plot(kind="line", marker="o", color="b", linewidth=2, markersize=8, label="Total Papers")

plt.title("Temporal Distribution of miRNA-Related Publications", fontsize=18, pad=20)
plt.xlabel("Year", fontsize=14)
plt.ylabel("Number of Publications", fontsize=14)
plt.xticks(rotation=45)

# Adjust y-axis to accommodate annotations
max_count = total_papers.max()
plt.ylim(top=max_count + 30)  # Add buffer space

# Add top 3 miRNAs annotations for each year
for year in total_papers.index:
    top3 = yearly_counts_T[year].nlargest(5)
    top3 = top3[top3 > 0]
    if not top3.empty:
        mirna_list = "\n".join(top3.index)
        plt.annotate(mirna_list, 
                     xy=(year, total_papers[year]),
                     xytext=(year, total_papers[year] + 10),
                     fontsize=10, color="black",
                     bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", lw=1.5),
                     ha='center', va='bottom')

# Highlight peak year
peak_year = total_papers.idxmax()
peak_count = total_papers.max()
plt.annotate(f"Peak: {peak_count} papers", 
             xy=(peak_year, peak_count),
             xytext=(peak_year, peak_count - 15),
             fontsize=12, color="red",
             arrowprops=dict(facecolor="red", shrink=0.05),
             ha='center', va='top')

plt.tight_layout()
plt.savefig("temporal_distribution.png", dpi=300, bbox_inches='tight')
plt.show()




exit()







# 2. Top 50 miRNA Bar Chart
top50_mirnas = yearly_counts_T.sum(axis=1).nlargest(50)
top50_data = yearly_counts_T.loc[top50_mirnas.index]

plt.figure(figsize=(14, 12))
top50_data.plot(kind="barh", stacked=True, colormap="viridis")
plt.title("Top 50 Most Frequent miRNAs by Year", fontsize=16)
plt.xlabel("Count", fontsize=14)
plt.ylabel("miRNA", fontsize=14)
plt.legend(title="Year", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
plt.tight_layout()
plt.savefig("top50_mirnas.png", dpi=300, bbox_inches='tight')
plt.show()