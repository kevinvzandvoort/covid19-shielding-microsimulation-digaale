#!/bin/bash

RUN=$1
WD="/nfs/general/covid_ibm_shielding"
FILE="/nfs/general/covid_ibm_shielding/combine_single_scen.R"

Rscript ${FILE} ${WD} ${RUN}
