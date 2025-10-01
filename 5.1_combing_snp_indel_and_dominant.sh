#!/bin/sh


working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data
cd $working_dir

# combing indel and snp together

########################################################################
# combine SNP and INDEL; find double hit and homo; combing them together
########################################################################
for sampleID in $(seq -f "SIT%03g" 43 88); do
  
	{
	head -n 1 SNP/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ_dominant.txt
	tail -n +2 SNP/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ_dominant.txt
	tail -n +2 INDEL/per_sample_all/$sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ_dominant.txt
	} | sort -k1,1 -k2,2n > combined/$sampleID/snp_indel_combined.filtered.$sampleID.dominant.txt
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