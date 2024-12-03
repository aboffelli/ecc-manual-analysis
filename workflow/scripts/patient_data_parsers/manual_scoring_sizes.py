import re

filename = '/home/aboffelli/PACC/TMA/Results/ManualScoreTables/merged_ecc_sizes.txt'
output = '/home/aboffelli/PACC/TMA/Results/ManualScoreTables/ecc_fold_change.txt'
gmm_file = '/home/aboffelli/PACC/TMA/Results/QuPathScoreTables/gmm_table.txt'


gmm_dictionary = {}
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

size_pattern = re.compile(
    r"""
    (Pathologist\s\d)\t    # Pathologist group 1
    (\d+A_EPCAM)\t    # Parent - group 2 
    (\d+\w)\t   # Core ID - group 3
    (.*)   # All the ECC sizes - group 4                             
    """, 
    re.VERBOSE
)


with open(gmm_file, 'r')as gmm:
    header = gmm.readline()
    #print(header)
    for line in gmm:
        match = gmm_pattern.search(line)
        if match:
            parent = match.group(1)
            core_id = match.group(2)
            mean_comp1 = float(match.group(3))
            mean_comp2 = float(match.group(4))
            prop_comp1 = float(match.group(5))
            prop_comp2 = float(match.group(6))

        if parent not in gmm_dictionary:
            gmm_dictionary[parent] = {}
        
        if prop_comp2 > 0.1:
            cancer_cell_mean = mean_comp2
        else:
            cancer_cell_mean = mean_comp1
        
        gmm_dictionary[parent][core_id] = round(cancer_cell_mean,2)


with open(filename) as sizes, open(output,'w') as out:
    print('Pathologist\tParent\tCore\tCoreMean\tPaccArea\tFoldChange', file=out)
    for line in sizes:
        match = size_pattern.search(line)
        if match:
            patho = match.group(1)
            parent = match.group(2)
            core_id = match.group(3)
            ecc_list = match.group(4).strip().split("\t")
            core_mean = gmm_dictionary.get(parent, {}).get(core_id, {})
            for size in ecc_list:
                size = float(re.sub(",",".", size))
                if core_mean:
                    fold_change = round(size/core_mean, 2)
                print(f"""{patho}\t{parent}\t{core_id}\t{core_mean}\
\t{size}\t{fold_change}""", file=out)
    

