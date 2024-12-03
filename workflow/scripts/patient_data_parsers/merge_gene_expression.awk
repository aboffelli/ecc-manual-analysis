#!/usr/bin/awk -f

BEGIN {
    FS=OFS="\t"
    } 

 NR==FNR {  # Process the first file
    gex[$1]=$2"\t"$3"\t"$4"\t"$5    # Create an array for the patient number
    next
    } 

 {
    print $0, gex[$1]   # Add the correct gex based on the patient number
    }