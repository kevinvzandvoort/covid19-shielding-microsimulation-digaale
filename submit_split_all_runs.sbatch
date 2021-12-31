#!/bin/bash
#SBATCH --job-name=covid_svk --array=1-1000 --cpus-per-task=1
/nfs/general/covid_ibm_shielding/split_run.sh $SLURM_ARRAY_TASK_ID
