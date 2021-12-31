pacman::p_load(qs, magrittr, data.table)

.args = if(interactive())
  c("~/workspace/covid_ibm_shielding/", 1) else
    commandArgs(trailingOnly = TRUE)

setwd(.args[1])
run = as.numeric(.args[2])

out_folder = "./output/"

dir.create(sprintf("%s/byscen", out_folder))
dir.create(sprintf("%s/byscen/model_out", out_folder))
dir.create(sprintf("%s/byscen/interventions", out_folder))

#scenario_list = readRDS("./data/clean/scenarios.RDS")
scens = 1:5292

message(sprintf("Starting run %s", run))

if(!"pacman" %in% rownames(installed.packages()))
  install.packages("pacman")
pacman::p_load(qs, data.table, magrittr)

out_run = sprintf("%s/out_run_%s.qs", out_folder, run) %>% qread
for(s in scens){
  message(sprintf("Scen %s/%s", s, length(scens)))
  out_run[["model_out"]] %>% .[scen == s] %>%
    qsave(sprintf("%s/byscen/model_out/out_run_%s_scen_%s_model_out.qs", out_folder, run, s))
  out_run[["model_interventions"]] %>% .[scen == s] %>%
    qsave(sprintf("%s/byscen/interventions/out_run_%s_scen_%s_model_interventions.qs", out_folder, run, s))
}
out_run[["population"]] %>%
  qsave(sprintf("%s/byscen/population/out_run_%s_population.qs", out_folder, run, s))

message("Finished")