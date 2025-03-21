import pandas as pd

# Set the path for the files.
patient_file = ('/home/aboffelli/PACC/TMA/Data/Tables/'
                'SB91RT_final_20170821_PrimoLR_EjPnr_orginal.csv')
tma_map = '/home/aboffelli/PACC/TMA/Data/Tables/Koppling_PathXL_till_TMA-nr.csv'

# Load the patient file
patient_df = pd.read_csv(patient_file, delimiter='\t')
print(patient_df.head(10))

# Load the TMA map file and make TMAnr into int
tma_df = pd.read_csv(tma_map, delimiter='\t').dropna()
tma_df['TMAnr'] = tma_df['TMAnr'].astype(int)

patient_df['PACC'] = None  # Create a column for the ECC score

# Replace the #NULL! to 0 to avoid errors when transforming it into float.
patient_df['tma_nr'] = patient_df['tma_nr'].replace('#NULL!', 0)  
patient_df['tma_nr'] = patient_df['tma_nr'].astype(float)
print(patient_df.head(10))

# Merge the tma map and patient table based on the TMA number.
result = pd.merge(tma_df, patient_df, right_on='tma_nr', left_on='TMAnr', 
                  how='left')

# Get the PACC score for each core
glas_scores = {}
with open("/home/aboffelli/PACC/TMA/Results/ManualScoreTables/complete_patho_results.csv", 
          "r") as scores:
    header = scores.readline()  # read the header
    for line in scores:
        line_split=line.strip('\n').split('\t')
        glas,core,pacc= line_split[0:3]  # Save the important info in variables to easy access

        if pacc:
            if len(pacc)>1:
              """
              Since the pathologists filled the table, sometimes there is text together 
              with the score, e.g. "2?", or "1 good".
              So we try to extract the first caracter and transform it into int.
              Since we can't account for all the possibilities, if the first character is not a number,
              we treat it as an NA.
              """
                try:
                    pacc = int(pacc[0])
                except ValueError:
                    pacc= 'NA'
            else:
                pacc = int(pacc)  # If there is only one character, set that as an ECC score.
        else: 
            pacc = 'NA'  # If the variable is empty, put it as an NA.

        # Get the TMA dictonary or create a new entry if it doesn't exist.
        glas_scores[glas] = glas_scores.get(glas, {})
        
        # Add the ECC score to the dictionary
        glas_scores[glas][core] = pacc
        
# Add the ECC score in the table
result["PACC"] = None

for index, row in result.iterrows():
    tma = row['PARENT']
    core = row['ID']
    score = glas_scores.get(tma, {}).get(core, 'NA')  # Retrieve the ECC score from the dictionary.
    result.at[index, 'PACC'] = score  # add to the dataframe.

# Save the file. This file contain the TMA numbers and duplicated patients (2 cores per patient).
result.sort_values(by='patnr').to_csv('/home/aboffelli/PACC/TMA/Results/all_patients_with_TMA_PACC_new.csv',
                                      index=False, sep='\t')

# Remove duplicate entries based on patient ID ('patnr') and tumor type ('tumor')
df_unique = result.drop_duplicates(subset=['patnr','tumor'])
# Drop unnecessary columns
df_unique = df_unique.drop(['TMAnr', 'PARENT','ID', 'PACC'], axis=1)      

# Convert the ECC score to numeric values, replacing errors with -1 to indicate missing or invalid data
# We need it as a number to calculate the max.
result['PACC'] = pd.to_numeric(result['PACC'], errors='coerce').fillna(-1)
# Aggregate ECC score by taking the maximum value per patient ('patnr')
pacc_agg = result.groupby('patnr')['PACC'].max().reset_index()

# Merge the unique patient data with the aggregated ECC scores
result_unique = pd.merge(df_unique, pacc_agg, on='patnr', 
                         how='left').sort_values(by='patnr')

# Identify patients present in 'patient_df' but missing from 'result_unique'
differing_values = patient_df[~patient_df['patnr'].isin(result_unique['patnr'])]

# Select patients who do not have TMA data
no_tma = patient_df[patient_df['tma_nr'] == 0]

# Select patients with tumor type labeled as "3.00"
tumor_three = patient_df[patient_df['tumor'] == "3.00"]

# Combine both groups and remove duplicates
to_add = pd.concat([no_tma, tumor_three], ignore_index=True).drop_duplicates(subset=['patnr','tumor'])


# Append missing patients back into the dataset
result_unique = pd.concat([result_unique, to_add], ignore_index=True)

# Remove any duplicates and ensure the dataset is sorted by 'patnr'
result_unique = result_unique.drop_duplicates(subset=['patnr','tumor']).sort_values('patnr')

# Change the -1 values back to NA.
result_unique['PACC'] = result_unique['PACC'].replace(-1, 'NA').fillna('NA')

# Save the final dataframe to a file
result_unique.to_csv('/home/aboffelli/PACC/TMA/Results/all_patient_unique_with_pacc_full_score_new.csv', 
                     index=False, sep='\t')
