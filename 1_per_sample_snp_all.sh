#!/bin/sh


working_dir=/data/workspaces/lag/workspaces/lg-sit/working_data/MengYun/code
cd $working_dir


# Loop through sample IDs and submit the job
for sampleID in $(seq -f "SIT%03g" 43 88); do
    
	qsub -v "sampleID=$sampleID" per_sample_snp_all.sge
done



#!/bin/bash

# Loop through sample IDs and submit the job
#for sampleID in $(seq -f "SIT%03g" 43 88); do
#    export sampleID  # Export sampleID as an environment variable for qsub
#    qsub per_sample_all.sge
#done

