#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

sex_chrom_samples=(SIT044 SIT047 SIT048 SIT053 SIT054 SIT058 SIT059 SIT061 SIT066 SIT070 SIT073 SIT074 SIT076 SIT078 SIT082 SIT086)

for sampleID in $(seq -f "SIT%03g" 43 88); do	
# 1 extract the AD values based on the gynotype values
awk -F'\t' '
NR == 1 {
    print $0 "\tSelected_Values"
    next
}
{
    split($268, alleles, /[\/|]/)   # split column 268 on / or |
    split($269, values, ",")        # split column 269 on ,
    
    selected=""
    for (i in alleles) {
        allele_index = alleles[i]
        # Adjust for 0-based index: 0 = first value
        selected = selected values[allele_index + 1] ","   # values[] is 1-based in awk

    }
    sub(/,$/, "", selected)   # remove trailing comma
    
    print $0 "\t" selected
}
' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_extracted.txt

# 2 keep rows where the AD ratio is between 1/3 and 3
awk -v sampleID="$sampleID" -F'\t' '
NR == 1 {
    header = $0
	print header > (sampleID "/snp_indel_combined.filtered." sampleID ".dominant.wCADD20_.txt")
    next
}
{
    split($NF, selected_values, ",")
    if (length(selected_values) < 2) next

    val1 = selected_values[1] + 0
    val2 = selected_values[2] + 0

	# Check both values are non-zero
    if (val1 == 0 || val2 == 0) next
	
    # Calculate ratio
    ratio = (val1 >= val2) ? val1 / val2 : val2 / val1

    # Filter ratio between 1/3 and 3
    if (ratio >= (1/3) && ratio <= 3) {
        gene = $12
        gene_count[gene]++
        lines[NR] = $0
    }
}
END {
    for (line in lines) {
        split(lines[line], fields, "\t")
        gene = fields[12]
        if (gene_count[gene] > 0) {
            # Print all columns except the last
            out = fields[1]
            for (i = 2; i < NF; i++) {
                out = out "\t" fields[i]
            }
            print out >> (sampleID "/snp_indel_combined.filtered." sampleID ".dominant.wCADD20_.txt")
        }
    }
}
' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_extracted.txt

# 3. if it is male, then drop Y and keep X with double hits; if female, then rename it
if [[ " ${sex_chrom_samples[*]} " =~ " ${sampleID} " ]]; then
  awk -F'\t' '
    NR == 1 {
      print
      next
    }
    $1 ~ /^chrX/ {
      # keep only rows where 268 has same alleles
      if ($268 ~ /^[0-9]+\/[0-9]+$/) {
        if (substr($268, 1, index($268, "/")-1) == substr($268, index($268, "/")+1)) {
          print
        }
      }
      next
    }
    $1 ~ /^chrY/ {
      # skip all chrY rows
      next
    }
    {
      # keep everything else
      print
    }
  ' $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_.txt > $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_cleaned.txt
  else
  # SampleID is NOT in the sex chromosome samples list — just copy the file
  cp $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_.txt $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_cleaned.txt
fi

	
rm $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_extracted.txt $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_.txt $sampleID/snp_indel_combined.filtered.$sampleID.dominant.wCADD20_clean.txt
	
done