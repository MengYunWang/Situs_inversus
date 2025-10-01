#!/bin/sh
#$ -N Filter_SI_SNP
#$ -cwd
#$ -q big.q
#$ -S /bin/bash
#$ -M mengyun.wang@mpi.nl
#$ -m beas


#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/SNP/per_sample_all/
cd $working_dir


for sampleID in $(seq -f "SIT%03g" 43 88); do

#1. filter out the exonic variants
    awk -F'\t' 'NR==1 || $11 == "exonic" || $11 == "exonic;splicing" || $11 == "ncRNA_exonic;splicing" || $11 == "splicing"' \
	$sampleID/rare_variants_SI_all.$sampleID.hg38_multianno_with_genotypes.txt > $sampleID/$sampleID.filtered.txt

#2. excluding synonymous SNV, unknown and .	
	awk -F'\t' 'NR == 1 || ($14 != "synonymous SNV" && !(($14 == "unknown" || $14 == ".") && $11 == "exonic"))'\
	$sampleID/$sampleID.filtered.txt > $sampleID/$sampleID.filtered.wo_syn_na_unkn.txt

 #3. excluding gnomad41_exome_AF > 0.01	
	#awk -F'\t' 'NR==1 || ($84 == "NA" || $84 == "." || $84 <= 0.01 || $84 !~ /^[0-9]/)'\
	#$sampleID/$sampleID.filtered.wo_syn_na_unkn.txt > $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.txt
	
	# excluding the gnome_AF>0.005
	awk -F'\t' 'NR==1 || ($94 == "NA" || $94 == "." || $94 <= 0.005 || $94 !~ /^[0-9]/)'\
	$sampleID/$sampleID.filtered.wo_syn_na_unkn.txt > $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.txt
	
#4. if CADD_pred >20 keep; if CADD_pred <20 but AlphaMissense_pred with P keep; if CADD_pred NA or . exclude only with AlphaMissense_pred with B
	awk -F'\t' '
	NR == 1 || ($204 > 20) || 
	($204 < 20 && $199 == "P") || 
	(($204 == "NA" || $204 == ".") && $199 != "B")
	' $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.txt > $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.txt

#5. Exclude if DP<7 or GQ <20	
	awk -F'\t' 'NR == 1 || ($271 >= 20 && $270 >= 7)' $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.txt > $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.DP_GQ.txt
    
	rm $sampleID/$sampleID.filtered.txt $sampleID/$sampleID.filtered.wo_syn_na_unkn.txt
	rm $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.txt $sampleID/$sampleID.filtered.wo_syn_na_unkn.AF.alpha_cadd.txt
	
done