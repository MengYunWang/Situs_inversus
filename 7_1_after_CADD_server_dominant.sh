#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir


#######################################################################
# Part I: After creating the vcf files and upload to the CADD server, then we will get the CADD score;
# Merging the CADD annotation with the original Annovar annotations
#######################################################################

#for sampleID in "${samples[@]}"; do
for sampleID in $(seq -f "SIT%03g" 43 88); do
# 1. Normalize the Chromosome Format and Sort Both Files by CHROM and POS
	cat $sampleID/GRCh38-v1.7_dominant.tsv.gz | tail -n +3 | sed 's/^/chr/' | sort -k1,1 -k2,2n > $sampleID/GRCh38-v1.7_dominant_sorted.txt # Add "chr" to GRCh38 if missing
	grep -v '^#' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.vcf | sort -k1,1 -k2,2n  > $sampleID/snp_indel_sorted.txt
	dos2unix $sampleID/GRCh38-v1.7_dominant_sorted.txt $sampleID/snp_indel_sorted.txt

# 2. Merge the Files
# Use awk to merge on CHROM and POS, keeping all rows from the VCF file
	awk -F'\t' '
	NR==FNR { key[$1 FS $2] = $0; next }  # Build a hash from GRCh38-v1.7_dominant_sorted.txt
	{
	if ($1 FS $2 in key)                # If CHROM and POS match in the hash
		print $0, key[$1 FS $2];          # Print the VCF row with the GRCh38 data
	else
		print $0, ".\t.";                 # If no match, add "."
	}
	' $sampleID/GRCh38-v1.7_dominant_sorted.txt $sampleID/snp_indel_sorted.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.CADD.txt
	
	# Extract the last column and add "CADD_PHRED" as the header
	(echo "CADD_PHRED"; awk -F'\t' '{print $NF}' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.CADD.txt) > $sampleID/CADD_extracted.txt

	# sorting the dominant file first before the append
		{
	head -n1 "$sampleID/snp_indel_combined.filtered.${sampleID}.dominant.txt"
	tail -n +2 "$sampleID/snp_indel_combined.filtered.${sampleID}.dominant.txt" | sort -t$'\t' -k1,1 -k2,2n
	} > "$sampleID/snp_indel_combined.filtered.${sampleID}.dominant.sorted.txt"
	
	# Append the extracted column to the target file
	paste $sampleID/snp_indel_combined.filtered.$sampleID.dominant.sorted.txt $sampleID/CADD_extracted.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD.txt
	
	rm  $sampleID/snp_indel_sorted.txt $sampleID/CADD_extracted.txt $sampleID/snp_indel_combined.filtered.$sampleID.dominant.CADD.txt
	
#if CADD_pred >20 keep
	# awk -F'\t' 'NR==1 || ($284 == "NA" || $284 == "." || $284 >= 20 || $284 !~ /^[0-9]/)'\
	# $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20.txt
	
	awk -F'\t' '
	NR == 1 || ($276 > 20) || 
	($276 < 20 && $199 == "P") || 
	(($276 == "NA" || $276 == ".") && $199 != "B")
	' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20.txt
	

   rm $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD.txt
done
