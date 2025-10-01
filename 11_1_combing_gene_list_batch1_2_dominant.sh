
#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/
cd $working_dir

set -euo pipefail


# build unsolved_batch1_2_dominant
cat batch1/dominant_model/unsolved_genes_control_gene_excluded.dominant.txt combined/genes_in_unsolved_cases_control_genes_excluded_dominant.txt \
  | sort -u > batch1_2/dominant/unsolved_batch1_2_dominant.txt
echo "→ unsolved_batch1_2_dominant.txt"


# build solved_batch1_2_dominant
cat batch1/dominant_model/solved_genes_control_gene_excluded.dominant.txt combined/genes_in_solved_cases_control_genes_excluded_dominant.txt \
  | sort -u > batch1_2/dominant/solved_batch1_2_dominant.txt
echo "→ solved_batch1_2_dominant.txt"


# build control_batch1_2_dominant
cat batch1/dominant_model/controls_genes_dominant.txt combined/genes_in_control_dominant.txt \
  | sort -u > batch1_2/dominant/control_batch1_2_dominant.txt
echo "→ control_batch1_2_dominant.txt"

# build cases_batch1_2_dominant
cat batch1/dominant_model/case_genes_dominant.txt combined/genes_in_case_dominant.txt \
  | sort -u > batch1_2/dominant/case_batch1_2_dominant.txt
echo "→ case_batch1_2_dominant.txt"

####### check if there is any overlap
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/batch1_2/dominant/
cd $working_dir

comm -12 <(sort case_batch1_2_dominant.txt) <(sort control_batch1_2_dominant.txt) \
  > overlap_batch1_2_dominant.txt
echo "→ overlap_batch1_2_dominant.txt"


# 3) subtract overlap from each original list
grep -Fvxf overlap_batch1_2_dominant.txt solved_batch1_2_dominant.txt  \
  > solved_genes_control_gene_excluded.batch1_2_dominant.txt
grep -Fvxf overlap_batch1_2_dominant.txt unsolved_batch1_2_dominant.txt \
  > unsolved_genes_control_gene_excluded.batch1_2_dominant.txt
grep -Fvxf overlap_batch1_2_dominant.txt control_batch1_2_dominant.txt \
  > control_genes_case_gene_excluded.batch1_2_dominant.txt

echo "→ solved_genes_control_gene_excluded.batch1_2_dominant.txt"
echo "→ unsolved_genes_control_gene_excluded.batch1_2_dominant.txt"
echo "→ control_genes_case_gene_excluded.batch1_2_dominant.txt"


