
# Import regular expressions module
import re


# Access the input files from Snakemake
filename = snakemake.input["merged_ecc_sizes"]      # File with ECC sizes per core
gmm_file = snakemake.input["gmm_table"]             # File with GMM means per core

# Access the output file from Snakemake
output = snakemake.output[0]


# Dictionary to store GMM means for each parent/core
gmm_dictionary = {}

# Regex pattern to parse GMM table lines
gmm_pattern = re.compile(
    r"""                  
    (\d+A_EPCAM)\t    # Parent - group 1
    (\d+\w)\t   # Image (core ID) - group 2
    (\d+\.\d+)\t    # mean component 1 - group 3
    (\d+\.\d+)\t    # mean component 2 - group 4
    \d+\.\d+\t    # Variance component 1
    \d+\.\d+\t    # Variance component 2
    (0\.\d+)\t  # proportion component 1 - group 5
    (0\.\d+)$  # proportion component 2 - group 6
    """,
    re.VERBOSE
)

# Regex pattern to parse ECC size lines
size_pattern = re.compile(
    r"""
    (Pathologist\s\d)\t    # Pathologist group 1
    (\d+A_EPCAM)\t    # Parent - group 2 
    (\d+\w)\t   # Core ID - group 3
    (.*)   # All the ECC sizes - group 4                             
    """, 
    re.VERBOSE
)



# Parse the GMM table and store the relevant mean for each parent/core
with open(gmm_file, 'r') as gmm:
    header = gmm.readline()
    # print(header)  # Uncomment for debugging
    for line in gmm:
        match = gmm_pattern.search(line)
        if match:
            parent = match.group(1)
            core_id = match.group(2)
            mean_comp1 = float(match.group(3))
            mean_comp2 = float(match.group(4))
            prop_comp1 = float(match.group(5))
            prop_comp2 = float(match.group(6))

            # Initialize dictionary for parent if needed
            if parent not in gmm_dictionary:
                gmm_dictionary[parent] = {}

            # Use mean of component 2 if its proportion is >0.1, else use component 1
            if prop_comp2 > 0.1:
                cancer_cell_mean = mean_comp2
            else:
                cancer_cell_mean = mean_comp1

            # Store the rounded mean in the dictionary
            gmm_dictionary[parent][core_id] = round(cancer_cell_mean, 2)



# Parse the ECC sizes file and write output with fold changes
with open(filename) as sizes, open(output, 'w') as out:
    # Write header to output file
    print('Pathologist\tParent\tCore\tCoreMean\tPaccArea\tFoldChange', file=out)
    for line in sizes:
        match = size_pattern.search(line)
        if match:
            patho = match.group(1)
            parent = match.group(2)
            core_id = match.group(3)
            ecc_list = match.group(4).strip().split("\t")
            # Get the mean value for this parent/core from the GMM dictionary
            core_mean = gmm_dictionary.get(parent, {}).get(core_id, {})
            for size in ecc_list:
                # Convert size to float, handling comma as decimal separator
                size = float(re.sub(",", ".", size))
                if core_mean:
                    fold_change = round(size / core_mean, 2)
                else:
                    fold_change = ''
                # Write the result to the output file
                print(f"""{patho}\t{parent}\t{core_id}\t{core_mean}\
\t{size}\t{fold_change}""", file=out)
    

