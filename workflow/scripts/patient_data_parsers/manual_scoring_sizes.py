
# Import regular expressions module
import re


# Access the input files from Snakemake
filename = snakemake.input["merged_ecc_sizes"]      # File with ECC sizes per core
gmm_file = snakemake.input["gmm_table"]             # File with GMM means per core

# Access the output file from Snakemake
output = snakemake.output[0]


# Dictionary to store GMM means for each parent/core
gmm_dictionary = {}

gmm_pattern = re.compile(  # Compile the pattern for gmm table.
    r"""                  
    (\d+A_EPCAM)\t    # Parent - group 1
    (\d+\w)\t   # Image (core ID) - group 2
    (\d+\.\d+)\t    # mean component 1 - group 3
    (\d+\.\d+)\t    # mean component 2 - group 4
    \d+\.\d+\t    # Variance component 1
    \d+\.\d+\t    # Variance component 2
    (0\.\d+)\t  # proportion component 1 - group 5
    (0\.\d+)\t  # proportion component 2 - group 6
    (\d+\.\d+)\t  # 99th percentile - group 7
    """,
    re.VERBOSE  # VERBOSE ignores the whitespaces inside the raw string. Makes it more visual.
)

size_pattern = re.compile(  # Compile the pattern from the manual size table.
    r"""
    (Pathologist\s\d)\t    # Pathologist group 1
    (\d+A_EPCAM)\t    # Parent - group 2 
    (\d+\w)\t   # Core ID - group 3
    (.*)   # All the ECC sizes - group 4                             
    """, 
    re.VERBOSE  # VERBOSE ignores the whitespaces inside the raw string. Makes it more visual.
)



# Parse the GMM table and store the relevant mean for each parent/core
with open(gmm_file, 'r') as gmm:
    header = gmm.readline()
    # print(header)  # Uncomment for debugging
    for line in gmm:
        match = gmm_pattern.search(line)  # Get the pattern in the gmm file
        if match:
            # Split the line into variables for easy access.
            parent = match.group(1)
            core_id = match.group(2)
            mean_comp1 = float(match.group(3))
            mean_comp2 = float(match.group(4))
            prop_comp1 = float(match.group(5))
            prop_comp2 = float(match.group(6))
            perc_99 = float(match.group(7))
        if parent not in gmm_dictionary:
            gmm_dictionary[parent] = {}  # Create a dictionary key with the TMA number, and a nested dictionary as value 
        
        if prop_comp2 > 0.1:
            """
            Since the gaussian mixture does not work well in non bi-modal distributions
            we check if the proportion of the 2nd component is higher than 0.1.
            This avoids using a very high mean when the core has many big cells in a 
            non bi-modal distribution.
            """
            cancer_cell_mean = (mean_comp2, perc_99)  # Use the 2nd comp mean if the proportion is > 0.1
        else:
            cancer_cell_mean = (mean_comp1, perc_99)  # Use the 1st comp mean if the proportion is <= 0.1
        
        gmm_dictionary[parent][core_id] = cancer_cell_mean  # Add the mean and 99th percentile to the dictionary as {TMA:{Core: (mean, 99_perc)}}


with open(filename) as sizes, open(output,'w') as out:
    print('Pathologist\tParent\tCore\tCoreMean\t99thPercentile\tECCArea\tFoldChange', file=out)  # Print the header in the output file.
    for line in sizes:
        match = size_pattern.search(line)  # Match the line on the manual size file.
        if match:
            # Split the line in variables for easy access.
            patho = match.group(1)
            parent = match.group(2)
            core_id = match.group(3)
            ecc_list = match.group(4).strip().split("\t")
            core_mean, p_99 = gmm_dictionary.get(parent, {}).get(core_id, {})  # Extract the mean and 99th percentile from the gmm dict using the TMA and core number.
            
            for size in ecc_list:
                size = float(re.sub(",",".", size))  # Fix the decimal ponctuation, change from "," to ".".
                
                if core_mean:
                    fold_change = round(size/core_mean, 2)  # Calculate the fold change 
                print(f"""{patho}\t{parent}\t{core_id}\t{round(core_mean, 2)}\t\
{round(p_99, 2)}\t{size}\t{fold_change}""", file=out)  # Print to output file
    
