import pandas as pd
import re
import ast
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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

###################################################################################################################################
## Open miRNA databases


miRNA_DAT = open('/home/miki/desktop/progetti_v1/20_Serino_colon_bibliometrix/01_Analisi_Bibliometrica/01_extract_miRNA/00_database_miRna/merged_miRNA.dat', 'r')
miRNA_mature = open('/home/miki/desktop/progetti_v1/20_Serino_colon_bibliometrix/01_Analisi_Bibliometrica/01_extract_miRNA/00_database_miRna/merged_mature.fa', 'r')

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

###################################################################################################################################



# Specifica il percorso del file Excel
file_path = "Table_S8_Author_Prod_over_Time_Docs_bibliometrix_2.xlsx"

# Leggi il file Excel
df = pd.read_excel(file_path, sheet_name="Sheet1")


# Imposta la riga 0 come intestazione delle colonne
df.columns = df.iloc[0]  # Prende la riga 0 e la imposta come nomi delle colonne
df = df[1:]  # Rimuove la riga 0 dal DataFrame (ora è l'intestazione)

# Resetta gli indici (opzionale, ma utile per pulire il DataFrame)
df = df.reset_index(drop=True)
df['TI'] = df['TI'].str.lower()
#df["TI"] = df["TI"].str.split(" ")

df["TI"] = df["TI"].apply(lambda x: re.split(r'[ /]', str(x)))



# Extract miRNA terms (normalize hyphens)
df["Mir_Words"] = df["TI"].apply(extract_mir_words)

df['validated_Mir_IDs'] = df['Mir_Words'].apply(
    lambda x: miRNA_ID_validation_considering_mature(x, Tolta_microRNA_ID)
)



df['validated_Mir_IDs'] = df['validated_Mir_IDs'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)



# Genera la matrice di conteggi
counts_matrix = (
    df.explode('validated_Mir_IDs')
    .groupby(['Author', 'validated_Mir_IDs'])
    .size()
    .unstack()  # Crea una riga per autore e colonna per miRNA
)


print(counts_matrix)


# Filtra i miRNA presenti in almeno 3 Autori
boolean_condition = (counts_matrix > 0).sum(axis=0) >= 1
counts_matrix_filtered = counts_matrix.loc[:, boolean_condition]

# Converti NaN in 0
#counts_matrix_filtered = counts_matrix_filtered.fillna(0)
print(counts_matrix_filtered)


# Heatmap senza clustering
plt.figure(figsize=(20, 10))
sns.heatmap(
    counts_matrix_filtered.T,
    cmap='YlOrRd',
    annot=False,
    fmt='g',
    cbar=True,
    cbar_kws={'label': 'number of miRNA per Author'},
    linewidths=0.5,
)

plt.title("miRNA Studied by at Least 3 Authors", fontsize=14)
plt.xlabel("Authors", fontsize=12)
plt.ylabel("miRNA ID", fontsize=12)
plt.xticks(rotation=90, ha='right')  # Ruota le etichette miRNA
plt.tight_layout()

# Salva il grafico
plt.savefig("heatmap_senza_dendrogramma.png", dpi=300, bbox_inches='tight')
plt.show()