
#!/bin/sh

#Define directories
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data
cd $working_dir

set -euo pipefail


####### build case_genes_recessive.txt
cat batch1/recessive_model/unsolved_cases_control_genes_excluded.recessive.txt combined/genes_in_unsolved_cases_control_genes_excluded_recessive.txt \
  | sort -u > batch1_2/recessive/unsolved_batch1_2_recessive.txt
echo "→ unsolved_batch1_2_recessive.txt"


cat batch1/recessive_model/solved_cases_control_genes_excluded.recessive.txt combined/genes_in_solved_cases_control_genes_excluded_recessive.txt \
  | sort -u > batch1_2/recessive/solved_batch1_2_recessive.txt
echo "→ solved_batch1_2_recessive.txt"


cat batch1/recessive_model/left_handers_cases_control_genes_excluded.recessive.txt combined/genes_in_left_handers_control_genes_excluded_recessive.txt \
  | sort -u > batch1_2/recessive/left_handers_batch1_2_recessive.txt
echo "→ left_handers_batch1_2_recessive.txt"


cat batch1/recessive_model/case_genes_recessive.txt combined/genes_in_case_recessive.txt \
  | sort -u > batch1_2/recessive/case_batch1_2_recessive.txt
echo "→ case_batch1_2_recessive.txt"

cat batch1/recessive_model/controls_genes_recessive.txt combined/genes_in_control_recessive.txt \
  | sort -u > batch1_2/recessive/control_batch1_2_recessive.txt
echo "→ control_batch1_2_recessive.txt"

####### check if there is any overlap
working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/data/batch1_2/recessive
cd $working_dir


comm -12 <(sort case_batch1_2_recessive.txt) <(sort control_batch1_2_recessive.txt) \
  > overlap_batch1_2_recessive.txt
echo "→ overlap_batch1_2_recessive.txt"


####### subtract overlap from each original list
grep -Fvxf overlap_batch1_2_recessive.txt solved_batch1_2_recessive.txt  \
  > solved_cases_control_genes_excluded.batch1_2_recessive.txt

grep -Fvxf overlap_batch1_2_recessive.txt unsolved_batch1_2_recessive.txt \
  > unsolved_cases_control_genes_excluded.batch1_2_recessive.txt

grep -Fvxf overlap_batch1_2_recessive.txt left_handers_batch1_2_recessive.txt \
  > left_handers_cases_control_genes_excluded.batch1_2_recessive.txt
  
grep -Fvxf overlap_batch1_2_recessive.txt control_batch1_2_recessive.txt \
  > control_case_genes_excluded.batch1_2_recessive.txt

echo "→ solved_cases_control_genes_excluded.batch1_2_recessive.txt"
echo "→ unsolved_cases_control_genes_excluded.batch1_2_recessive.txt"
echo "→ control_case_genes_excluded.batch1_2_recessive.txt"