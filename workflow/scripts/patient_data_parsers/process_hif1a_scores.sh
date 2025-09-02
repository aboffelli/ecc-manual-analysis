#!/bin/bash

output_file="$1"
shift  # Remove the first argument, so $@ only has inputs

# Step 1: Merge files using awk and store in a variable (avoid intermediate file)
merged=$(mktemp)  # Create a temporary file
awk 'FNR==1 && NR!=1 {next} {print}' "$@" > "$merged"

# Step 2: Extract HIF1A scores, process inline
hif_scores=$(mktemp)
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {print $2,$1,$28,$29,$30} 
$6 ~ /HIF1a/ {split($2, parent, " "); 
print parent[2],$1,$28,$29,$30}' "$merged"  | \
sed -e 's/0%/0/g' -e 's/1-10/1/g' -e 's/11-50/2/g' -e 's/>50/3/g' \
    -e 's/No score//g' -e 's/a/A/' > "$hif_scores"

# Step 3: Merge pathologist scores, keeping the highest score
awk 'BEGIN {FS=OFS="\t"; last_core="1A"; highest4=0; highest5=0} 
/ID/ { print $1,$2,$4,$5; next } 
{ 
  if ($2 == last_core) {
    if ($4 > highest4) highest4=$4
    if ($5 > highest5) highest5=$5
  } else {
    print "Glas "$1,last_core,highest4,highest5
    last_core=$2
    highest4=$4
    highest5=$5
  }
} 
END {print "Glas "$1,last_core,highest4,highest5}' "$hif_scores" > "$output_file"