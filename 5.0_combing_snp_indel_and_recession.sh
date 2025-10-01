#!/bin/sh


working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data
cd $working_dir

# combing indel and snp together

########################################################################
# combine SNP and INDEL; find double hit and homo; combing them together
########################################################################
for sampleID in $(seq -f "SIT%03g" 43 88); do

  rm -rf combined/"$sampleID" 
  mkdir -p combined/"$sampleID"
  
	{
	head -n 1 SNP/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ.txt
	tail -n +2 SNP/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ.txt
	tail -n +2 INDEL/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ.txt
	} | sort -k1,1 -k2,2n > combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt
    
# Extract the gynotype (268>GT) starting with 0* >>> hetero
  awk -F'\t' 'NR == 1 || $268 ~ /^0[\/|0-9]/' combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.txt
 
# Find variants with genes appearing more than once in hetero >> potentientially double hit
# Substep 1: Find duplicate values in the 12th column>gene names
  awk -F'\t' '{count[$12]++} END {for (val in count) if (count[val] > 1) print val}' combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.txt > combined/$sampleID/duplicate_genes.txt
# Substep 2: Extract rows from the original file based on duplicate values
  awk -F'\t' 'FNR==1 && NR!=FNR { print; next } NR==FNR { duplicates[$1]=1; next } $12 in duplicates' \
  combined/$sampleID/duplicate_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes.txt


# Find the genes in hetero gynotype with "/" appears >> at least one site hit another parent 
# substep 1: Find rows where column 268>>GT contains "/"
  awk -F'\t' 'NR == 1 ||$268 ~ /\// { print $0 }' combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes.txt > combined/$sampleID/rows_with_slash.txt
# subtep 2: Extract unique gene names from column 12
  awk -F'\t' '{ print $12 }' combined/$sampleID/rows_with_slash.txt | sort | uniq > combined/$sampleID/unique_genes.txt
# substep 3: Filter original file using the extracted gene names
  awk -F'\t' 'NR==FNR { genes[$1]=1; next } $12 in genes' combined/$sampleID/unique_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes_hit.txt

# Extract the gynotype starting not with 0* >>> homo 
  awk -F'\t' 'NR == 1 || $268 !~ /^0[\/|0-9]/' combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.homo.txt

# Combing the homo and double hit hetero  
  	{
	head -n 1 combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.txt
	tail -n +2 combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes_hit.txt
	tail -n +2 combined/$sampleID/snp_indel_combined.filtered.$sampleID.homo.txt
	} | sort -k1,1 -k2,2n > combined/$sampleID/snp_indel_combined.filtered.$sampleID.recession.txt
  
  rm combined/$sampleID/duplicate_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.hetero.genes.txt combined/$sampleID/rows_with_slash.txt combined/$sampleID/unique_genes.txt
done

###########################################
# strict model with same procedure as above
# nonframeshift discarded
###########################################

for sampleID in $(seq -f "SIT%03g" 43 88); do
	
# Exclude nonframeshift variants 14th>functional annotation
  awk -F'\t' 'NR == 1 || $14 !~ /^nonframe/' combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.txt
    
# Extract the gynotype starting with 0* >>> hetero
  awk -F'\t' 'NR == 1 || $268 ~ /^0[\/|0-9]/' combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.txt
 
# Find variants with genes appearing more than once in hetero >> potentientially double hit
# Substep 1: Find duplicate values in the 12th column
  awk -F'\t' '{count[$12]++} END {for (val in count) if (count[val] > 1) print val}' combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.txt > combined/$sampleID/duplicate_genes.txt
# Substep 2: Extract rows from the original file based on duplicate values
  awk -F'\t' 'FNR==1 && NR!=FNR { print; next } NR==FNR { duplicates[$1]=1; next } $12 in duplicates' \
  combined/$sampleID/duplicate_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes.txt


# Find the genes in hetero gynotype with "/" appears   >> at least one site hit another parent 
# substep 1: Find rows where column 268 contains "/"
  awk -F'\t' 'NR == 1 ||$268 ~ /\// { print $0 }' combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes.txt > combined/$sampleID/rows_with_slash.txt
# subtep 2: Extract unique gene names from column 12
  awk -F'\t' '{ print $12 }' combined/$sampleID/rows_with_slash.txt | sort | uniq > combined/$sampleID/unique_genes.txt
# substep 3: Filter original file using the extracted gene names
  awk -F'\t' 'NR==FNR { genes[$1]=1; next } $12 in genes' combined/$sampleID/unique_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes_hit.txt

# Extract the gynotype starting without 0*
  awk -F'\t' 'NR == 1 || $268 !~ /^0[\/|0-9]/' combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.homo.txt
  
  	{
	head -n 1 combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.txt
	tail -n +2 combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes_hit.txt
	tail -n +2 combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.homo.txt
	} | sort -k1,1 -k2,2n > combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.txt
  
  rm combined/$sampleID/duplicate_genes.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.strict.hetero.genes.txt combined/$sampleID/rows_with_slash.txt combined/$sampleID/unique_genes.txt
done



	
# exluding the uniq genes which appears only once
# Step 1: Extract duplicate values from the 12th column
 # cut -f12 combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt | sort | uniq -d > combined/$sampleID/duplicate_values_$sampleID.txt
## Step 2: Filter rows where the 12th column matches the duplicate values
  #awk -F'\t' '
  #NR==FNR { keep[$1]=1; next }  # Process the first file (duplicates) and store values in "keep" array
  #FNR==1 { print; next }         # Print the header row from the second file
  #$12 in keep { print }          # Print rows where column 12 matches a value in "keep"
  #' combined/$sampleID/duplicate_values_$sampleID.txt combined/$sampleID/snp_indel_combined.filtered.$sampleID.txt > combined/$sampleID/snp_indel_combined.filtered.$sampleID.genes.txt