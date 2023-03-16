#!/bin/bash

RUN=$1
WD="/nfs/general/covid_ibm_shielding"
FILE="/nfs/general/covid_ibm_shielding/single_run.R"
ONLY_UNMITIGATED=0
Rscript ${FILE} ${WD} ${RUN} ${ONLY_UNMITIGATED}
