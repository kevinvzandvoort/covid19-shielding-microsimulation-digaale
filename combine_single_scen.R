pacman::p_load(qs, magrittr, data.table)

.args = if(interactive())
  c("/nfs/general/covid_ibm_shielding/", 5292) else
    commandArgs(trailingOnly = TRUE)

setwd(.args[1])
s = as.numeric(.args[2])

out_folder = "./output/"

runs = 1:1000

message(sprintf("Starting scen %s", s))

if(!"pacman" %in% rownames(installed.packages()))
  install.packages("pacman")
pacman::p_load(qs, data.table, magrittr)

model_out = list()
model_interventions = list()
for(r in runs){
  message(r)
  #if(!file.exists(sprintf("%s/byscen/model_out/out_run_%s_scen_%s_model_out.qs", out_folder, r, s))) message(sprintf("No: run %s, scen %s", r, s))
  model_out[[length(model_out)+1]] = sprintf("%s/byscen/model_out/out_run_%s_scen_%s_model_out.qs", out_folder, r, s) %>% qread %>% .[, run := r]
  model_interventions[[length(model_out)+1]] = sprintf("%s/byscen/interventions/out_run_%s_scen_%s_model_interventions.qs", out_folder, r, s) %>% qread %>% .[, run := r]
}
model_out = rbindlist(model_out)
model_interventions = rbindlist(model_interventions)

model_out %>%
  qsave(sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, s))

model_interventions %>%
  qsave(sprintf("%s/byscen/interventions/byscen/out_scen_%s_model_interventions.qs", out_folder, s))

message("Finished")
