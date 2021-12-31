#!/bin/bash
#SBATCH --job-name=covid_svk --array=3303,3307,4431-4509,4511-4516,4518-4523,4525-4554,4556-4676,4678-4726,4728-4803,4805-4810,4812-4991,4993-5086,5088-5135,5137-5227,5229-5292 --cpus-per-task=1
/nfs/general/covid_ibm_shielding/combine_single_scen.sh $SLURM_ARRAY_TASK_ID
