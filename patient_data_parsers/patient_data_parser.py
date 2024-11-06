import pandas as pd

patient_file = ('/home/aboffelli/PACC/TMA/Data/Tables/'
                'SB91RT_final_20170821_PrimoLR_EjPnr_orginal.csv')
tma_map = '/home/aboffelli/PACC/TMA/Data/Tables/Koppling_PathXL_till_TMA-nr.csv'

# Load the patient file
patient_df = pd.read_csv(patient_file, delimiter='\t')
print(patient_df.head(10))

# Load the TMA map file and make TMAnr into int
tma_df = pd.read_csv(tma_map, delimiter='\t').dropna()
tma_df['TMAnr'] = tma_df['TMAnr'].astype(int)

patient_df['PACC'] = None; print(patient_df.head(10))
patient_df['tma_nr'] = patient_df['tma_nr'].replace('#NULL!', 0); print(patient_df.head(10))
patient_df['tma_nr'] = patient_df['tma_nr'].astype(float)
print(patient_df.head(10))

result = pd.merge(tma_df, patient_df, right_on='tma_nr', left_on='TMAnr', 
                  how='left'); print(result.sort_values('patnr').head(20))

# Get the PACC score for each core
glas_scores = {}
with open("/home/aboffelli/PACC/TMA/Results/ManualScoreTables/complete_patho_results.csv", 
          "r") as scores:
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

result.sort_values(by='patnr').to_csv('/home/aboffelli/PACC/TMA/Results/all_patients_with_TMA_PACC_new.csv',
                                      index=False, sep='\t')
paired_samples = result[result['tumor'].isin(["1.00","2.00"])]

#result = result[result['tumor'] == '1.00']

# Merge the cores and scores for the same patients.
print(result)
df_unique = result.drop_duplicates(subset=['patnr','tumor']); print(df_unique) 
df_unique = df_unique.drop(['TMAnr', 'PARENT','ID', 'PACC'], axis=1)      

result['PACC'] = pd.to_numeric(result['PACC'], errors='coerce').fillna(-1)
pacc_agg = result.groupby('patnr')['PACC'].max().reset_index()

result_unique = pd.merge(df_unique, pacc_agg, on='patnr', 
                         how='left').sort_values(by='patnr'); print(result_unique)
#result_unique.to_csv('/home/aboffelli/TMA/Results/all_patient_unique_with_pacc_full_score.csv', 
#                     index=False, sep='\t')

# Merge the cores and score for the paired sample patients.
df_grouped = paired_samples.groupby('patnr')['tumor'].nunique()
patients_with_both_types = df_grouped[df_grouped == 2].index


result_paired = paired_samples[paired_samples['patnr'].isin(patients_with_both_types)]

df_unique = result_paired.drop_duplicates(subset=['patnr','tumor'])

df_unique = df_unique.drop(['TMAnr', 'PARENT','ID', 'PACC'], axis=1)

result_paired['PACC'] = pd.to_numeric(result_paired['PACC'], 
                                      errors='coerce').fillna(-1).astype('Int64')
pacc_agg = result_paired.groupby(['patnr','tumor'])['PACC'].max().reset_index()
result_paired_unique = pd.merge(df_unique, pacc_agg, on=['patnr', 'tumor'], 
                        how='left').sort_values(by='patnr')
#result_paired_unique.to_csv('/home/aboffelli/TMA/Results/paired_patient_unique_with_pacc_full_score.csv', 
 #                    index=False, sep='\t')




differing_values = patient_df[~patient_df['patnr'].isin(result_unique['patnr'])]

print(differing_values)

no_tma = patient_df[patient_df['tma_nr'] == 0]; print(no_tma)
tumor_three = patient_df[patient_df['tumor'] == "3.00"]; print(tumor_three)
to_add = pd.concat([no_tma, tumor_three], ignore_index=True).drop_duplicates(subset=['patnr','tumor'])
print(to_add)

result_unique = pd.concat([result_unique, to_add], 
                          ignore_index=True).drop_duplicates(subset=['patnr','tumor']).sort_values('patnr')
result_unique['PACC'] = result_unique['PACC'].replace(-1, 'NA').fillna('NA')
print(result_unique.head(10))

result_unique.to_csv('/home/aboffelli/PACC/TMA/Results/all_patient_unique_with_pacc_full_score_new.csv', 
                     index=False, sep='\t')

# print(patient_df)