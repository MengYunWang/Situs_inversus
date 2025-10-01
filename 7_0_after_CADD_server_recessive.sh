#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir


#######################################################################
# Part I: After creating the vcf files and upload to the CADD server, then we will get the CADD score;
# Merging the CADD annotation with the original Annovar annotations
#######################################################################
sex_chrom_samples=(SIT044 SIT047 SIT048 SIT053 SIT054 SIT058 SIT059 SIT061 SIT066 SIT070 SIT073 SIT074 SIT076 SIT078 SIT082 SIT086)

#for sampleID in "${samples[@]}"; do
for sampleID in $(seq -f "SIT%03g" 43 88); do
# 1. Normalize the Chromosome Format and Sort Both Files by CHROM and POS
	cat $sampleID/GRCh38-v1.7_recessive.tsv.gz | tail -n +3 | sed 's/^/chr/' | sort -k1,1 -k2,2n > $sampleID/GRCh38-v1.7_recessive_sorted.txt # Add "chr" to GRCh38 if missing
	grep -v '^#' $sampleID/snp_indel_combined.filtered.$sampleID.recession.vcf | sort -k1,1 -k2,2n  > $sampleID/snp_indel_sorted.txt
	dos2unix $sampleID/GRCh38-v1.7_recessive_sorted.txt $sampleID/snp_indel_sorted.txt

# 2. Merge the Files
# Use awk to merge on CHROM and POS, keeping all rows from the VCF file
	awk -F'\t' '
	NR==FNR { key[$1 FS $2] = $0; next }  # Build a hash from GRCh38-v1.7_recessive_sorted.txt
	{
	if ($1 FS $2 in key)                # If CHROM and POS match in the hash
		print $0, key[$1 FS $2];          # Print the VCF row with the GRCh38 data
	else
		print $0, ".\t.";                 # If no match, add "."
	}
	' $sampleID/GRCh38-v1.7_recessive_sorted.txt $sampleID/snp_indel_sorted.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.CADD.txt
	
	# Extract the last column and add "CADD_PHRED" as the header
	(echo "CADD_PHRED"; awk -F'\t' '{print $NF}' $sampleID/snp_indel_combined.filtered.$sampleID.recession.CADD.txt) > $sampleID/CADD_extracted.txt

	# sorting the recession file first before the append
	{
	head -n1 "$sampleID/snp_indel_combined.filtered.${sampleID}.recession.txt"
	tail -n +2 "$sampleID/snp_indel_combined.filtered.${sampleID}.recession.txt" | sort -t$'\t' -k1,1 -k2,2n
	} > "$sampleID/snp_indel_combined.filtered.${sampleID}.recession.sorted.txt"

	# Append the extracted column to the target file
	paste $sampleID/snp_indel_combined.filtered.$sampleID.recession.sorted.txt $sampleID/CADD_extracted.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD.txt
	
	rm $sampleID/snp_indel_sorted.txt $sampleID/CADD_extracted.txt $sampleID/snp_indel_combined.filtered.${sampleID}.recession.CADD.txt
	
#if CADD_pred >20 keep
	# awk -F'\t' 'NR==1 || ($284 == "NA" || $284 == "." || $284 >= 20 || $284 !~ /^[0-9]/)'\
	# $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20.txt
	
	awk -F'\t' '
	NR == 1 || ($276 > 20) || 
	($276 < 20 && $199 == "P") || 
	(($276 == "NA" || $276 == ".") && $199 != "B")
	' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt
	
   # Find variants with genes appearing more than once in hetero >> potentientially double hit
	# Substep 1: Find duplicate values in the 12th column>gene names
   awk -F'\t' '{count[$12]++} END {for (val in count) if (count[val] > 1) print val}' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt > $sampleID/duplicate_genes.txt
	# Substep 2: Extract rows from the original file based on duplicate values
   awk -F'\t' 'FNR==1 && NR!=FNR { print; next } NR==FNR { duplicates[$1]=1; next } $12 in duplicates' \
  $sampleID/duplicate_genes.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt
  
   # Extract the gynotype starting not with 0* >>> homo 
   awk -F'\t' 'NR == 1 || $268 !~ /^0[\/|0-9]/' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt

	if [[ " ${sex_chrom_samples[*]} " =~ " ${sampleID} " ]]; then
	   awk -F'\t' 'NR == 1 || $1 ~ /^chr[XY]/' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt
	fi  
	
	{
	head -n 1 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt"
	[[ -f "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt" ]] && tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt"
	} | sort -t$'\t' -k1,1 -k2,2n -u \
	> "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20.txt"

   rm $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_.txt 
   rm $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt
done


#######################################################################
# Part II: After creating the vcf files and upload to the CADD server, then we will get the CADD score;
# Merging the CADD annotation with the original Annovar annotations
#######################################################################

# #for sampleID in "${samples[@]}"; do
# for sampleID in $(seq -f "SIT%03g" 43 88); do
# # 1. Normalize the Chromosome Format and Sort Both Files by CHROM and POS
	# cat $sampleID/GRCh38-strict_recession-v1.7*.tsv.gz | tail -n +3 | sed 's/^/chr/' | sort -k1,1 -k2,2n > $sampleID/GRCh38_strict_recession_sorted.txt # Add "chr" to GRCh38 if missing
	# grep -v '^#' $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.vcf | sort -k1,1 -k2,2n  > $sampleID/snp_indel_sorted_strict.txt
	# dos2unix $sampleID/GRCh38_strict_recession_sorted.txt $sampleID/snp_indel_sorted_strict.txt
# # 2. Merge the Files
# # Use awk to merge on CHROM and POS, keeping all rows from the VCF file
	# awk -F'\t' '
	# NR==FNR { key[$1 FS $2] = $0; next }  # Build a hash from GRCh38_strict_recession_sorted.txt
	# {
	# if ($1 FS $2 in key)                # If CHROM and POS match in the hash
		# print $0, key[$1 FS $2];          # Print the VCF row with the GRCh38 data
	# else
		# print $0, ".\t.";                 # If no match, add "."
	# }
	# ' $sampleID/GRCh38_strict_recession_sorted.txt $sampleID/snp_indel_sorted_strict.txt > $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.CADD.txt
	
	# # Extract the last column and add "CADD_PHRED" as the header
	# (echo "CADD_PHRED"; awk -F'\t' '{print $NF}' $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.CADD.txt) > $sampleID/CADD_extracted_strict.txt

	# # Append the extracted column to the target file
	# paste $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.txt $sampleID/CADD_extracted_strict.txt > $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD.txt
	
	# rm $sampleID/GRCh38_strict_recession_sorted.txt $sampleID/snp_indel_sorted_strict.txt $sampleID/CADD_extracted_strict.txt $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.CADD.txt
	
# #if CADD_pred >20 keep
	# # awk -F'\t' 'NR==1 || ($284 == "NA" || $284 == "." || $284 >= 20 || $284 !~ /^[0-9]/)'\
	# # $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD20.txt
	
	# awk -F'\t' '
	# NR == 1 || ($276 > 20) || 
	# ($276 < 20 && $199 == "P") || 
	# (($276 == "NA" || $276 == ".") && $199 != "B")
	# ' $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD20.txt
	
	# rm $sampleID/snp_indel_combined.filtered.$sampleID.strict.recession.wCADD.txt
# done