import pandas as pd
import pyreadstat

# Load the data
df, meta = pyreadstat.read_dta("/media/aboffelli/ArthurExternalHD/BackUp/TMA/Data/Tables/GEX_46k_with_gene_names.dta")

# Select a subset of columns
subset = df[['patnr', 'tumor', 'hif1a', 'epas1']] 
print(subset)