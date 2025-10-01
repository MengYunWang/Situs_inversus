#!/bin/sh


#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

##############################################
# Part I: query the gene in customized gene list
##############################################

grep -v '^$' ../../code/genes_to_query.txt \
  | tr '[:lower:]' '[:upper:]' \
  | sort -u \
  > genes_to_query.txt

dos2unix genes_to_query.txt


samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 12, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_control_genes_excluded.txt"
  
  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | sort -u > $sample/genes_control_genes_excluded_in_$sample.dominant.txt
 
  # Find the common genes between cases and query genes (case-insensitive via uppercase)
  comm -12 <(tr '[:lower:]' '[:upper:]' < genes_to_query.txt | sort) <(tr '[:lower:]' '[:upper:]' < $sample/genes_control_genes_excluded_in_$sample.dominant.txt | sort) > $sample/genes_control_genes_excluded_query_in_overlap_dominant.txt
  
  dos2unix $sample/genes_control_genes_excluded_query_in_overlap_dominant.txt
  
  # keep the overlap genes in cases
  (head -n 1 "$file"; grep -F -f $sample/genes_control_genes_excluded_query_in_overlap_dominant.txt "$file" | tail -n +2) > $sample/snp_indel_combined.filtered.$sample.dominant.wCADD20_control_genes_excluded_queried.txt
	
done


##############################################
# Part II: query gene in the clinvar gene list
##############################################

awk -F'\t' 'NR>1 {print $2; print $3}' Clinvar_gene_condition_source_id.txt | sort -u > genes_to_query_clinvar.txt
dos2unix genes_to_query_clinvar.txt 

samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 12, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_control_genes_excluded.txt"

  # Find the common genes between controls and cases
  comm -12 genes_to_query_clinvar.txt $sample/genes_control_genes_excluded_in_$sample.dominant.txt > $sample/genes_control_genes_excluded_clinvar_in_overlap_dominant.txt
  
  dos2unix $sample/genes_control_genes_excluded_clinvar_in_overlap_dominant.txt
  
  # keep the overlap genes in cases
  (head -n 1 "$file"; grep -F -f $sample/genes_control_genes_excluded_clinvar_in_overlap_dominant.txt "$file") > $sample/snp_indel_combined.filtered.$sample.dominant.wCADD20_control_genes_excluded_clinvar.txt

  # keep the overlap genes in clinvar
  (head -n 1 Clinvar_gene_condition_source_id.txt; grep -F -w -f $sample/genes_control_genes_excluded_clinvar_in_overlap_dominant.txt Clinvar_gene_condition_source_id.txt) > $sample/genes_control_genes_excluded_in_clinvar_dominant.txt
  
  # rm $sample/genes_$sample.in_clinvar.txt
done


##################################################
## output the unique genes in the filtered results
##################################################
# Loop through each sample, extract column 12, and output the unique genes in the filtered results
# for sample in "${samples[@]}"; do
  # file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_control_genes_excluded.txt"
  # # Extract the 12th column (adjust the field separator if needed)
  # awk -F'\t' '{print $12}' "$file"
# done | grep -v '^Gene\.ensGene$' | sort -u > genes_in_case_after_gene_filtering.txt


# #############################################
# ## output the unique genes and its frequency in the filtered results
# #############################################
# # Loop through each sample, extract the 12th column, and output the unique genes in the filtered results and its frequency
# for sample in "${samples[@]}"; do
  # file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_control_genes_excluded.txt"
  # if [[ -f "$file" ]]; then
    # # Extract column 12, remove the header line, and get unique genes from this sample.
    # awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | sort -u
  # else
    # echo "Warning: File '$file' not found." >&2
  # fi
# done | sort | uniq -c > gene_freq_across_samples_after_gene_filtering.txt