import pandas as pd

# Access the inputs and outputs
patient_file = snakemake.input["patient_file"]
tma_map = snakemake.input["tma_map"]
patho_scores = snakemake.input["patho_score"]

patient_table_cores = snakemake.output["patient_cores"]
patient_table_unique = snakemake.output["patient_ecc"]


# Load the patient file
patient_df = pd.read_csv(patient_file, delimiter='\t')
#print(patient_df.head(10))

# Load the TMA map file and make TMAnr into int
tma_df = pd.read_csv(tma_map, delimiter='\t').dropna()
tma_df['TMAnr'] = tma_df['TMAnr'].astype(int)

patient_df['PACC'] = None
patient_df['tma_nr'] = patient_df['tma_nr'].replace('#NULL!', 0)
patient_df['tma_nr'] = patient_df['tma_nr'].astype(float)


result = pd.merge(tma_df, patient_df, right_on='tma_nr', left_on='TMAnr', 
                  how='left')

# Get the PACC score for each core
glas_scores = {}
with open(patho_scores, "r") as scores:
    header = scores.readline()
    for line in scores:
        line_split=line.strip('\n').split('\t')
        glas,core,pacc= line_split[0:3]

        if pacc:
            if len(pacc)>1:
                try:
                    pacc = int(pacc[0])
                except ValueError:
                    pacc= 'NA'
            else:
                pacc = int(pacc)
        else: 
            pacc = 'NA'
        
        glas_scores[glas] = glas_scores.get(glas, {})
        
        # Considering PACC score
        glas_scores[glas][core] = pacc
        
# Put the PACC score in the table
result["PACC"] = None

for index, row in result.iterrows():
    tma = row['PARENT']
    core = row['ID']
    score = glas_scores.get(tma, {}).get(core, 'NA')
    result.at[index, 'PACC'] = score

result.sort_values(by='patnr').to_csv(patient_table_cores,
                                      index=False, sep='\t')

# Merge the cores and scores for the same patients.
df_unique = result.drop_duplicates(subset=['patnr','tumor'])
df_unique = df_unique.drop(['TMAnr', 'PARENT','ID', 'PACC'], axis=1)      

result['PACC'] = pd.to_numeric(result['PACC'], errors='coerce').fillna(-1)
pacc_agg = result.groupby('patnr')['PACC'].max().reset_index()

result_unique = pd.merge(df_unique, pacc_agg, on='patnr', 
                         how='left').sort_values(by='patnr')


differing_values = patient_df[~patient_df['patnr'].isin(result_unique['patnr'])]

no_tma = patient_df[patient_df['tma_nr'] == 0]
tumor_three = patient_df[patient_df['tumor'] == "3.00"]

to_add = (
    pd.concat([no_tma, tumor_three], ignore_index=True)
    .drop_duplicates(subset=['patnr','tumor'])
)

result_unique = (
    pd.concat([result_unique, to_add], ignore_index=True)
    .drop_duplicates(subset=['patnr','tumor'])
    .sort_values('patnr')
)

result_unique['PACC'] = (
    result_unique['PACC']
    .replace(-1, 'NA')
    .fillna('NA')
)


result_unique.to_csv(patient_table_unique, index=False, sep='\t')

