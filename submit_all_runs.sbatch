#!/bin/bash
#SBATCH --job-name=covid_svk --array=270,280,282-289 --cpus-per-task=2
/nfs/general/covid_ibm_shielding/single_run.sh $SLURM_ARRAY_TASK_ID
