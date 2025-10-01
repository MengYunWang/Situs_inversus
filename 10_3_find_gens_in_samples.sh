
#!/usr/bin/env bash
set -euo pipefail

# control sample ids
samples=(
  SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069
  SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076
  SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086
  SIT087 SIT088
)

for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20.txt"

  # Extract column 12, split on ';', then check for DNAH6
  if awk -F'\t' '{print $12}' "$file" \
       | grep -v '^Gene\.ensGene$' \
       | tr ';' '\n' \
       | grep -Fxq "DNAH6" 
  then
    echo "DNAH6 found in sample $sample" 
  fi
done




#!/usr/bin/env bash
set -euo pipefail

# control sample ids
samples=(
  SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069
  SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076
  SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086
  SIT087 SIT088
)

for sample in "${samples[@]}"; do
  file="${sample}/snp_indel_combined.filtered.${sample}.dominant.wCADD20_cleaned.txt"

  # Extract column 12, split on ';', then check for DNAH8
  if awk -F'\t' '{print $12}' "$file" \
       | grep -v '^Gene\.ensGene$' \
       | tr ';' '\n' \
       | grep -Fxq "DNAH8" 
  then
    echo "DNAH8 found in sample $sample" 
  fi
done



set -euo pipefail

# Define control sample IDs
samples=(SIT045 SIT047 SIT049 SIT052 SIT054 SIT058 SIT063 SIT064 SIT067 SIT068 SIT078 SIT080 SIT085)

# Define gene list (exact match expected)
genes=(
  CCDC66  
  CFAP77  
  CLIP2  
  DCTN1  
  DNAI7  
  EML6  
  KATNAL2  
  MAP2  
  MYH7B  
  MYO7B  
  RIMS2  
  RP1L1
)

# Loop over genes and samples
for gene in "${genes[@]}"; do
  for sample in "${samples[@]}"; do
    file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"

    if awk -F'\t' '{print $12}' "$file" \
         | grep -v '^Gene\.ensGene$' \
         | tr ';' '\n' \
         | grep -Fxq "$gene"; then
      echo "$gene found in sample $sample"
    fi
  done
done


set -euo pipefail

# Define control sample IDs
samples=(
  SIT046 SIT048 SIT050 SIT053 SIT057 SIT065 SIT069
  SIT070 SIT071 SIT072 SIT073 SIT074 SIT075 SIT076
  SIT077 SIT079 SIT081 SIT082 SIT083 SIT084 SIT086
  SIT087 SIT088
)

# Define gene list (exact match expected)
genes=(
APC2  
CFAP46  
DNAH6  
DNHD1  
HOOK1  
MID1IP1  
NAV3  
PCNT  
RADIL  
RP1  
SKA3  
SNPH  
TEKT4
)

# Loop over genes and samples
for gene in "${genes[@]}"; do
  for sample in "${samples[@]}"; do
    file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"

    if awk -F'\t' '{print $12}' "$file" \
         | grep -v '^Gene\.ensGene$' \
         | tr ';' '\n' \
         | grep -Fxq "$gene"; then
      echo "$gene found in sample $sample"
    fi
  done
done


set -euo pipefail

# Define control sample IDs
samples=(
  SIT003 SIT004 SIT005 SIT007 SIT010 SIT011 SIT013
  SIT016 SIT017 SIT018 SIT019 SIT022 SIT023 SIT026
  SIT028 SIT029 SIT030 SIT031 SIT032 SIT033 SIT034
  SIT035 SIT036 SIT037 SIT038 SIT039 SIT040 SIT041
  SIT042 SIT009
)

# Define gene list (exact match expected)
genes=(
  CCDC66  
  CFAP77  
  CLIP2  
  DCTN1  
  DNAI7  
  EML6  
  KATNAL2  
  MAP2  
  MYH7B  
  MYO7B  
  RIMS2  
  RP1L1
APC2  
CFAP46  
DNAH6  
DNHD1  
HOOK1  
MID1IP1  
NAV3  
PCNT  
RADIL  
RP1  
SKA3  
SNPH  
TEKT4
)

# Loop over genes and samples
# Loop over genes and samples
for gene in "${genes[@]}"; do
  for sample in "${samples[@]}"; do
    file="${sample}_df.xls"

    if [[ -f "$file" ]]; then
      if awk -F'\t' 'NR > 1 {split($1, a, ";"); for (i in a) if (a[i] == "'$gene'") found = 1} END {exit !found}' "$file"; then
        echo "$gene found in sample $sample"
      fi
    else
      echo "File not found: $file"
    fi
  done
done




set -euo pipefail

# Define unsolved cases
samples=(SIT045 SIT047 SIT049 SIT052 SIT054 SIT058 SIT062 SIT063 SIT064 SIT067 SIT068 SIT078 SIT080 SIT085)

# Define gene list from unsolved cases (exact match expected)
genes=(
ABCG5
AC092170.1
AC093323.1
ACAP1
ADGRL4
AL390778.1
AL591806.1
ANKAR
AP003041.2
ARPP21
ASAP3
ATP13A5
BCL11A
BLK
BRI3BP
C2CD5
CARD10
CCDC66
CD34
CFAP77
CLIP2
COL4A6
CSTF2
DCTN1
DHX34
DIS3L
DMBT1
DNAI7
DRP2
ECH1
ELP2
EML6
FCGBP
GNRH2
GRID2IP
HS6ST2
HTR2C
IGHV3-42
IRS4
KATNAL2
KIAA0586
KRTAP5-1
LOC100132874
LOC100653133
LOC100996350
LOC101928884
LOC101928917
LOC101929113
LZTR1
MAP2
MATN4
MED15
MEGF8
MUC3A
MYH7B
MYO7B
NCOA3
OR4F17
OTOGL
OTUD7A
PASK
PCDH19
PCSK4
PCSK5
PCSK9
PID1
PLA2G4E
POLR1G
PSKH2
RAB40AL
RECQL5
RFX2
RGS22
RIMS2
RP11-119F19.2
RP13-672B3.2
RP1L1
SGK223
SLC38A5
TASOR2
TFR2
TMPO
TRAF3IP2
TRIO
TRIOBP
UBE4B
UMODL1
USPL1
ZACN
ZDHHC11
ZNF302
ZNF66
)

cd /data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/combined
# Loop over genes and samples
for gene in "${genes[@]}"; do
  for sample in "${samples[@]}"; do
    file="${sample}/snp_indel_combined.filtered.${sample}.recession.wCADD20_cleaned.txt"

    if awk -F'\t' '{print $12}' "$file" \
         | grep -v '^Gene\.ensGene$' \
         | tr ';' '\n' \
         | grep -Fxq "$gene"; then
      echo "$gene found in sample $sample"
    fi
  done
done


cd /data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/batch1/recessive_model
samples=(SIT004 SIT005 SIT007 SIT010 SIT013)

for gene in "${genes[@]}"; do
  for sample in "${samples[@]}"; do
    file="${sample}_df.xls"

    if [[ -f "$file" ]]; then
      if awk -F'\t' 'NR > 1 {split($1, a, ";"); for (i in a) if (a[i] == "'$gene'") found = 1} END {exit !found}' "$file"; then
        echo "$gene found in sample $sample"
      fi
    else
      echo "File not found: $file"
    fi
  done
done