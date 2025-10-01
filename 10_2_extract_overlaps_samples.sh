#!/bin/sh
set -euo pipefail


#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
cd $working_dir

samples=(
  SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069
  SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076
  SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086
  SIT087 SIT088
)

# prepare output
echo -e "Sample\tRecessiveCount\tOverlapCount" > genes_in_control_counts_recessive.txt

for sample in "${samples[@]}"; do
  f1="${sample}/genes_in_${sample}.recessive.txt"
  f2="${sample}/genes_in_overlap_recessive.txt"

  # count lines (using < so we just get the number, no filename)
  c1=$(wc -l < "$f1")
  c2=$(wc -l < "$f2")

  # print a tab-separated summary
  echo -e "${sample}\t${c1}\t${c2}"
done >> genes_in_control_counts_recessive.txt


# prepare output
echo -e "Sample\tRecessiveCount\tOverlapCount" > genes_in_control_counts_dominant.txt

for sample in "${samples[@]}"; do
  f1="${sample}/genes_in_${sample}.dominant.txt"
  f2="${sample}/genes_in_overlap_dominant.txt"

  # count lines (using < so we just get the number, no filename)
  c1=$(wc -l < "$f1")
  c2=$(wc -l < "$f2")

  # print a tab-separated summary
  echo -e "${sample}\t${c1}\t${c2}"
done >> genes_in_control_counts_dominant.txt




samples=(SIT043 SIT044 SIT045 SIT047 SIT049 SIT051 SIT052 SIT054 SIT055 SIT056 SIT058 SIT059 SIT060 SIT061 SIT062 SIT063 SIT064 SIT066 SIT067 SIT068 SIT078 SIT080 SIT085)

# prepare output
echo -e "Sample\tRecessiveCount\tOverlapCount" > genes_in_cases_counts_recessive.txt

for sample in "${samples[@]}"; do
  f1="${sample}/genes_in_${sample}.recessive.txt"
  f2="${sample}/genes_in_overlap_recessive.txt"

  # count lines (using < so we just get the number, no filename)
  c1=$(wc -l < "$f1")
  c2=$(wc -l < "$f2")

  # print a tab-separated summary
  echo -e "${sample}\t${c1}\t${c2}"
done >> genes_in_cases_counts_recessive.txt


# prepare output
echo -e "Sample\tRecessiveCount\tOverlapCount" > genes_in_cases_counts_dominant.txt

for sample in "${samples[@]}"; do
  f1="${sample}/genes_in_${sample}.dominant.txt"
  f2="${sample}/genes_in_overlap_dominant.txt"

  # count lines (using < so we just get the number, no filename)
  c1=$(wc -l < "$f1")
  c2=$(wc -l < "$f2")

  # print a tab-separated summary
  echo -e "${sample}\t${c1}\t${c2}"
done >> genes_in_cases_counts_dominant.txt