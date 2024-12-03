# Workflow of the TMA scripts

## PACC score
Merge the manual analysis scores to the patient list, using the TMA map to 
track down which TMAs belong to each patient.

Files needed:
- SB91RT_final_20170821_PrimoLR_EjPnr_orginal.csv
- Koppling_PathXL_till_TMA-nr.csv
- PACC_TMA_CH_RM_final.xlsx

### File ***complete_patho_results.csv***
```sh
for f in *scores.csv; do 
    awk 'OFS=FS="\t" {print $1,$2,$3}' $f > ${f/.csv/_upd.csv}; 
done

awk -F '\t' 'FNR == 1 && NR != 1 {next} {print $0}' \
Scores/*upd.csv > complete_patho_results.csv
```

### File all_patients_with_TMA_PACC_new.csv
Comes from the python script ***patient_data_parser.py***

## HIF Score merging

Files needed:
- SWEBCG_pathXL_export_all_230405.xlsx
- all_patients_with_TMA_PACC_new.csv

Original file: ***SWEBCG_pathXL_export_all_230405.xlsx***
Each sheet was saved as a csv tab delimited.

### File merged.csv
Merge all the csv files into one.
```sh
awk 'FNR==1 && NR!=1 {next} {print}' \
SWEBCG_pathXL_export_all_230405*.csv > merged.csv
```

### File hif_scores.csv
Get only the HIF1A score and convert the percentages to scores.

```sh
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $2,$1,$28,$29,$30} 
$6 ~ /HIF1a/ {split($2, parent, " "); 
print parent[2],$1,$28,$29,$30}' merged.csv | \
sed -e 's/0%/0/g' -e 's/1-10/1/g' -e 's/11-50/2/g' -e 's/>50/3/g' \
-e 's/No score//g' -e 's/a/A/'> hif_scores.csv
```

### File hif_scores_merged.csv
Merge both pathologists scores, keeping the highest score.

```sh
awk 'BEGIN {FS=OFS="\t"; last_core="1A"; highest=0} 
/ID/ { print $1,$2,$4; next } { if ($2 == last_core) 
{ if ($4 > highest) {highest=$4}} else {print $1,last_core,highest; 
last_core=$2; highest=$4}} END {print $1,last_core,highest,high_int}' \
hif_scores.csv > hif_scores_merged.csv
```

### File all_patients_with_TMA_PACC_HIF.csv
Add the HIF1A score in the patient data table.
Comes from the R script ***hif1a_pacc.R***

## Gene expression data merging

Files needed:
- GEX_46k_with_gene_names.csv 
- all_patients_with_TMA_PACC_new.csv

### File gex_interesting_genes.csv
Isolate only the interesting genes out of the expression data file 
(cdk1,epas1,hif1a,mki67)
```sh
awk -F "," 'BEGIN {OFS="\t"} {print $46057,$46056,$10708,$13533}' \
GEX_46k_with_gene_names.csv | sed 's/\r//g' > gex_interesting_genes.csv
```
### Generating patient file with GEX
Merge the gene expression to the patient data table.

```sh
 awk 'BEGIN {FS=OFS="\t"} 
 NR==FNR {gex[$1]=$2"\t"$3"\t"$4"\t"$5; next} 
 {print $0, gex[$1]}' ../Data/Tables/gex_interesting_genes.csv \
 all_patients_with_TMA_PACC_HIF.csv > all_patients_with_TMA_PACC_HIF_GEX.csv
```

## QuPath detection



## Calculating nuclear area fold change

Files needed:
- _measurement.txt (measurements file from QuPath)
- merged_ecc_sizes.txt (manual analysed ECCs sizes)

### Fit the Gaussian mixture model
The script ***gaussian_mixture.R*** reads the measurements table and outputs a combined gmm table. Make sure you input all your files in the script.

### Calculate ECC fold change
Using the new file gmm_table.txt as an input in the python script ***manual_scoring_sizes.py*** together with the merged_ecc_sizes, you generate the file ecc_fold_change.txt

### Plotting the fold change
Generate a file with all measurements from qupath
```sh
awk 'BEGIN {FS=OFS="\t"; print "Parent","Image","Area"} FNR == 1 {next} {print FILENAME,$1,$6}' *A_EPCAM_measurements.txt > all_cells_area.txt
```

Then use the script ***fold_change_plots.R*** with all the three files: ***gmm_table.txt***, ***ecc_fold_change.txt***, and ***all_cells_area.txt***.
