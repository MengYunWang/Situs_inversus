#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

sex_chrom_samples=(SIT044 SIT047 SIT048 SIT053 SIT054 SIT058 SIT059 SIT061 SIT066 SIT070 SIT073 SIT074 SIT076 SIT078 SIT082 SIT086)

for sampleID in $(seq -f "SIT%03g" 43 88); do
# 1. extract the homo genotypes without phasing
awk -F'\t' '
NR == 1 || 
($268 ~ /^[0-9]+\/[0-9]+$/ && 
  substr($268, 1, index($268, "/")-1) == substr($268, index($268, "/")+1))
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt

# 2.0 extract the hetero genotypes, double or more hits on one gene
awk -F'\t' '
NR == 1 || 
!($268 ~ /^[0-9]+\/[0-9]+$/ && 
  substr($268, 1, index($268, "/")-1) == substr($268, index($268, "/")+1))
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_.txt

# 2.1 extract the AD values based on the gynotype values
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
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_extracted_.txt



# 2.2 keep rows where the AD ratio is between 1/3 and 3, and these genes appears more than once
awk -v sampleID="$sampleID" -F'\t' '
NR == 1 {
    header = $0
	print header > (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w0_.txt")
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
        if (gene_count[gene] > 1) {
            # Print all columns except the last
            out = fields[1]
            for (i = 2; i < NF; i++) {
                out = out "\t" fields[i]
            }
            print out >> (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w0_.txt")
        }
    }
}
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_extracted_.txt

awk -v sampleID="$sampleID" -F'\t' '
NR == 1 {
    header = $0
	print header > (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w0.txt")
    next
}
{
    gene = $12
    lines[NR] = $0
    gene_lines[NR] = gene

    # Check for a forward slash in column 268
    if ($268 ~ /\//) {
        gene_has_slash[gene] = 1
    }
}
END {
    for (line in lines) {
        gene = gene_lines[line]
        if (gene_has_slash[gene]) {
            print lines[line] >> (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w0.txt")
        }
    }
}
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w0_.txt


# 2.3 keep rows where the AD ratio is between 1/3 and 3, and genotype contains no-0
awk -v sampleID="$sampleID" -F'\t' '
NR == 1 {
    header = $0
	print header > (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w_0.txt")
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
        # Filter: column 268 must NOT contain 0
        if ($268 !~ /0/) {
            lines[NR] = $0
        }
    }
}
END {
    for (line in lines) {
        split(lines[line], fields, "\t")
        # Print all columns except the last
        out = fields[1]
        for (i = 2; i < NF; i++) {
            out = out "\t" fields[i]
        }
        print out >> (sampleID "/snp_indel_combined.filtered." sampleID ".recession.wCADD20_hetero_w_0.txt")
    }
}
' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_extracted_.txt


# 2.4 combing the hetero file together
	{
	head -n 1 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w_0.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w0.txt"
	} | sort -t$'\t' -k1,1 -k2,2n -u \
	> $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt

rm $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w_0.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w0_.txt
rm $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_w0.txt

# 3 extract the Chrom X genotype only keep those values with same values
if [[ " ${sex_chrom_samples[*]} " =~ " ${sampleID} " ]]; then
  awk -F'\t' '
    NR == 1 || 
    ($1 ~ /^chrX/ &&
      $268 ~ /^[0-9]+\/[0-9]+$/ &&
      substr($268, 1, index($268, "/")-1) == substr($268, index($268, "/")+1))
  ' $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20.txt > $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt
fi

	
# 4. Combing the homo and double hit hetero and sex
	{
	head -n 1 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero.txt"
    tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_homo.txt"
	[[ -f "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt" ]] && tail -n +2 "$sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_sex.txt"
	} | sort -t$'\t' -k1,1 -k2,2n -u \
	> $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_cleaned.txt
	
rm $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_.txt $sampleID/snp_indel_combined.filtered.$sampleID.recession.wCADD20_hetero_extracted_.txt
	
done