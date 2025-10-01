#!/bin/sh


#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

##############################################
# Part I: Extract the gene list in the control
##############################################
samples=(SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069 SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076 SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086 SIT087 SIT088)
# Loop through each sample, extract column 12, and output the unique genes in controls
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"

  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file"
  
done | grep -v '^Gene\.ensGene$' | sort -u > genes_in_control.txt
dos2unix genes_in_control.txt

awk -F';' '{for(i=1;i<=NF;i++) print $i}' genes_in_control.txt | sort -u > genes_in_control_recessive.txt
dos2unix genes_in_control_recessive.txt
rm genes_in_control.txt


###################################################
# Part II: Filter out the (control) genes from case
###################################################
samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 12, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"
  
  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | sort -u > $sample/genes_in_$sample.recessive.txt
  dos2unix $sample/genes_in_$sample.recessive.txt
  
  # Find the common genes between controls and cases
  comm -12 genes_in_control_recessive.txt $sample/genes_in_$sample.recessive.txt > $sample/genes_in_overlap_recessive.txt
  dos2unix $sample/genes_in_overlap_recessive.txt
  
  # Filter out the overlap genes in cases
  out=${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_control_genes_excluded.txt
  rm $out  
  
  # append filtered body
  awk -F'\t' -v OFS='\t' -v overlap="$sample/genes_in_overlap_recessive.txt" '
    # load overlap-list
    FNR==NR { bad[$1]; next }
	{
      # split field12 by ;
      n = split($12, genes, ";")
      # if any sub-gene is in bad[], skip row
      for(i=1; i<=n; i++){
        if(genes[i] in bad) next
      }
      # otherwise print entire row
      print
    }
  ' "$sample/genes_in_overlap_recessive.txt" "$file" \
  >> "$out"
  
  
  #(head -n 1 "$file"; grep -v -F -f $sample/genes_in_overlap_recessive.txt "$file" | tail -n +2) > $sample/snp_indel_combined.filtered.$sample.recession.wCADD20_control_genes_excluded.txt
  
  # rm $sample/genes_in_overlap.txt  $sample/genes_in_$sample.txt
done

## do the mirrorr exclusion

##############################################
# Part III: Extract the gene list in the case
##############################################
samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 12, and output the unique genes in controls
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"

  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file"
  
done | grep -v '^Gene\.ensGene$' | sort -u > genes_in_case.txt

awk -F';' '{for(i=1;i<=NF;i++) print $i}' genes_in_case.txt | sort -u > genes_in_case_recessive.txt
dos2unix genes_in_case_recessive.txt
rm genes_in_case.txt


###################################################
# Part IV: Filter out the (case) genes from control
###################################################
samples=(SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069 SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076 SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086 SIT087 SIT088)
# Loop through each sample, extract column 12, and filter out the common gene and output the filtered results
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"
  
  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | sort -u > $sample/genes_in_$sample.recessive.txt
  dos2unix $sample/genes_in_$sample.recessive.txt
  
  # Find the common genes between controls and cases
  comm -12 genes_in_case_recessive.txt $sample/genes_in_$sample.recessive.txt > $sample/genes_in_overlap_recessive.txt
  dos2unix $sample/genes_in_overlap_recessive.txt
  
  # Filter out the overlap genes in cases
  out=$sample/snp_indel_combined.filtered.$sample.recession.wCADD20_case_genes_excluded.txt
  
  rm $out
  
  # append filtered body
  awk -F'\t' -v OFS='\t' -v overlap="$sample/genes_in_overlap_recessive.txt" '
    # load overlap-list
    FNR==NR { bad[$1]; next }
	{
      # split field12 by ;
      n = split($12, genes, ";")
      # if any sub-gene is in bad[], skip row
      for(i=1; i<=n; i++){
        if(genes[i] in bad) next
      }
      # otherwise print entire row
      print
    }
  ' "$sample/genes_in_overlap_recessive.txt" "$file" \
  >> "$out"
  
  #(head -n 1 "$file"; grep -v -F -f $sample/genes_in_overlap_recessive.txt "$file" | tail -n +2) > $sample/snp_indel_combined.filtered.$sample.recession.wCADD20_case_genes_excluded.txt

  # rm $sample/genes_in_overlap.txt  $sample/genes_in_$sample.txt

done


##################################################
## output the unique genes in the filtered results
##################################################
samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract column 12, and output the unique genes in the filtered results
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_control_genes_excluded.txt"
  # Extract the 12th column (adjust the field separator if needed)
  awk -F'\t' '{print $12}' "$file"
done | grep -v '^Gene\.ensGene$' | sort -u > genes_in_case_control_genes_excluded.txt

awk -F';' '{for(i=1;i<=NF;i++) print $i}' genes_in_case_control_genes_excluded.txt | sort -u > genes_in_case_control_genes_excluded_recessive.txt
dos2unix genes_in_case_control_genes_excluded_recessive.txt
rm genes_in_case_control_genes_excluded.txt


####################################################################
## output the unique genes and its frequency in the filtered results
####################################################################
samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)
# Loop through each sample, extract the 12th column, and output the unique genes in the filtered results and its frequency
for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_control_genes_excluded.txt"
  if [[ -f "$file" ]]; then
    # Extract column 12, remove the header line, and get unique genes from this sample.
    awk -F'\t' '{print $12}' "$file" | grep -v '^Gene\.ensGene$' | awk -F';' '{ for(i=1;i<=NF;i++) print $i }' | sort -u
  else
    echo "Warning: File '$file' not found." >&2
  fi
done | sort | uniq -c > genes_in_case_wfreq_control_genes_excluded_recessive.txt