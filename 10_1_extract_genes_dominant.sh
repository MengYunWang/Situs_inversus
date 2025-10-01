#!/bin/sh


#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

###############################################
# Part I: extract the genes in the solved cases
###############################################

output_file="genes_in_solved_cases_control_genes_excluded_dominant.txt"
> "$output_file"

samples=(SIT043 SIT044 SIT051 SIT055 SIT056 SIT059 SIT060 SIT061 SIT062  SIT066)
# Loop through each sample, append all the text into one 
for sample in "${samples[@]}"; do
	cat $sample/genes_control_genes_excluded_in_$sample.dominant.txt >> "$output_file"
	
done
sort -u -o genes_in_solved_cases_control_genes_excluded_dominant.txt genes_in_solved_cases_control_genes_excluded_dominant.txt && dos2unix genes_in_solved_cases_control_genes_excluded_dominant.txt


##################################################
# Part II: extract the genes in the unsloved cases
##################################################

output_file="genes_in_unsolved_cases_control_genes_excluded_dominant.txt"
> "$output_file"

samples=(SIT045 SIT047 SIT049 SIT052 SIT054 SIT058 SIT063 SIT064 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 7, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
	cat $sample/genes_control_genes_excluded_in_$sample.dominant.txt >> "$output_file"
done
sort -u -o genes_in_unsolved_cases_control_genes_excluded_dominant.txt genes_in_unsolved_cases_control_genes_excluded_dominant.txt && dos2unix genes_in_unsolved_cases_control_genes_excluded_dominant.txt


##############################################
# Part III: extract the genes in SI with PCD
##############################################

output_file="genes_in_PCD_cases_control_genes_excluded_dominant.txt"
> "$output_file"

samples=(SIT044 SIT055)
# Loop through each sample, extract column 7, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  cat $sample/genes_control_genes_excluded_in_$sample.dominant.txt >> "$output_file"
done
sort -u -o genes_in_PCD_cases_control_genes_excluded_dominant.txt genes_in_PCD_cases_control_genes_excluded_dominant.txt && dos2unix genes_in_PCD_cases_control_genes_excluded_dominant.txt


##############################################
# Part IV: extrac the genes in SI without PCD
##############################################

output_file="genes_in_nPCD_cases_control_genes_excluded_dominant.txt"
> "$output_file"

samples=(SIT043 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 7, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  cat $sample/genes_control_genes_excluded_in_$sample.dominant.txt >> "$output_file"
done
sort -u -o genes_in_nPCD_cases_control_genes_excluded_dominant.txt genes_in_nPCD_cases_control_genes_excluded_dominant.txt && dos2unix genes_in_nPCD_cases_control_genes_excluded_dominant.txt



##############################################
# Part V: extrac the genes in controls
##############################################

output_file="genes_in_control_case_genes_excluded_dominant.txt"
> "$output_file"

samples=(SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069 SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076 SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086 SIT087 SIT088)
# Loop through each sample, extract column 7, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  
    file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_case_genes_excluded.txt"
  
  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | sort -u > $sample/genes_case_genes_excluded_in_$sample.dominant.txt
  
  cat $sample/genes_case_genes_excluded_in_$sample.dominant.txt >> "$output_file" 
  
  # rm $sample/genes_case_genes_excluded_in_$sample.txt
done
sort -u -o genes_in_control_case_genes_excluded_dominant.txt genes_in_control_case_genes_excluded_dominant.txt && dos2unix genes_in_control_case_genes_excluded_dominant.txt

##################################################
## output the unique genes in the filtered results
##################################################
# Loop through each sample, extract column 7, and output the unique genes in the filtered results
# for sample in "${samples[@]}"; do
  # file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_genes_filtered.txt"
  # # Extract the 7th column (adjust the field separator if needed)
  # awk '{print $7}' "$file"
# done | grep -v '^Gene\.refGeneWithVer$' | sort -u > genes_in_case_after_gene_filtering.txt


# #############################################
# ## output the unique genes and its frequency in the filtered results
# #############################################
# # Loop through each sample, extract the 7th column, and output the unique genes in the filtered results and its frequency
# for sample in "${samples[@]}"; do
  # file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_genes_filtered.txt"
  # if [[ -f "$file" ]]; then
    # # Extract column 7, remove the header line, and get unique genes from this sample.
    # awk '{print $7}' "$file" | grep -v '^Gene\.refGeneWithVer$' | sort -u
  # else
    # echo "Warning: File '$file' not found." >&2
  # fi
# done | sort | uniq -c > gene_freq_across_samples_after_gene_filtering.txt