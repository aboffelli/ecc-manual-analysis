
# Import pandas for data manipulation
import pandas as pd


# Access the Snakemake inputs and outputs
patient_file = snakemake.input["patient_file"]      # Patient data file
tma_map = snakemake.input["tma_map"]                # TMA map file
patho_scores = snakemake.input["patho_score"]       # Pathologist scores file

patient_table_cores = snakemake.output["patient_cores"]   # Output: patient-core table
patient_table_unique = snakemake.output["patient_ecc"]    # Output: unique patient table



# Load the patient file as a DataFrame
patient_df = pd.read_csv(patient_file, delimiter='\t')
# print(patient_df.head(10))  # Uncomment for debugging


# Load the TMA map file, drop missing values, and ensure TMAnr is integer
tma_df = pd.read_csv(tma_map, delimiter='\t').dropna()
tma_df['TMAnr'] = tma_df['TMAnr'].astype(int)


# Prepare patient DataFrame: add PACC column, clean tma_nr
patient_df['PACC'] = None
patient_df['tma_nr'] = patient_df['tma_nr'].replace('#NULL!', 0)
patient_df['tma_nr'] = patient_df['tma_nr'].astype(float)



# Merge TMA map and patient data on TMA number
result = pd.merge(tma_df, patient_df, right_on='tma_nr', left_on='TMAnr', 
                  how='left')


# Parse the pathologist scores file and build a nested dictionary: glas_scores[glas][core] = pacc
glas_scores = {}
with open(patho_scores, "r") as scores:
    header = scores.readline()
    for line in scores:
        line_split = line.strip('\n').split('\t')
        glas, core, pacc = line_split[0:3]

        # Clean and convert PACC score
        if pacc:
            if len(pacc) > 1:
                try:
                    pacc = int(pacc[0])
                except ValueError:
                    pacc = 'NA'
            else:
                pacc = int(pacc)
        else:
            pacc = 'NA'

        glas_scores[glas] = glas_scores.get(glas, {})
        # Store the PACC score for each core
        glas_scores[glas][core] = pacc

# Assign the PACC score to each row in the merged table
result["PACC"] = None

for index, row in result.iterrows():
    tma = row['PARENT']
    core = row['ID']
    score = glas_scores.get(tma, {}).get(core, 'NA')
    result.at[index, 'PACC'] = score

# Save the patient-core table with PACC scores
result.sort_values(by='patnr').to_csv(patient_table_cores,
                                      index=False, sep='\t')
# Create a unique patient table by merging cores and scores for the same patients
df_unique = result.drop_duplicates(subset=['patnr', 'tumor'])
df_unique = df_unique.drop(['TMAnr', 'PARENT', 'ID', 'PACC'], axis=1)

# Aggregate: for each patient, get the max PACC score across all cores
result['PACC'] = pd.to_numeric(result['PACC'], errors='coerce').fillna(-1)
pacc_agg = result.groupby('patnr')['PACC'].max().reset_index()

# Merge unique patient info with aggregated PACC scores
result_unique = pd.merge(df_unique, pacc_agg, on='patnr', 
                         how='left').sort_values(by='patnr')

# Find patients in the original file not present in the unique table
differing_values = patient_df[~patient_df['patnr'].isin(result_unique['patnr'])]

# Add patients with no TMA or tumor == 3.00
no_tma = patient_df[patient_df['tma_nr'] == 0]
tumor_three = patient_df[patient_df['tumor'] == "3.00"]

to_add = (
    pd.concat([no_tma, tumor_three], ignore_index=True)
    .drop_duplicates(subset=['patnr', 'tumor'])
)

# Concatenate and deduplicate the final unique patient table
result_unique = (
    pd.concat([result_unique, to_add], ignore_index=True)
    .drop_duplicates(subset=['patnr', 'tumor'])
    .sort_values('patnr')
)

# Clean up PACC column: replace -1 and NaN with 'NA'
result_unique['PACC'] = (
    result_unique['PACC']
    .replace(-1, 'NA')
    .fillna('NA')
)

# Save the unique patient table
result_unique.to_csv(patient_table_unique, index=False, sep='\t')

