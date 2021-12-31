#!/bin/bash

RUN=$1
WD="/nfs/general/covid_ibm_shielding"
FILE="/nfs/general/covid_ibm_shielding/single_run.R"

Rscript ${FILE} ${WD} ${RUN}
