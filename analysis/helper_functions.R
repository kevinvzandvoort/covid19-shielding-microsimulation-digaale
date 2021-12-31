formatItRange = function(range){
  formatSet = function(set){
    if(length(set) > 1)
      out = paste(c(set[1], set[length(set)]), collapse = "-")
    else
      out = set
    return(out)
  }
  
  range = range %>% as.numeric %>% na.omit %>% unique %>% sort
  
  if(length(range) < 2)
    return(range)
  
  current = range[1]
  current_set = current
  out = NULL
  for(i in 2:length(range)){
    if(range[i] - current == 1)
      current_set = c(current_set, range[i])
    else {
      out = c(out, formatSet(current_set))
      current_set = range[i]
    }
    current = range[i]
  }
  out = c(out, formatSet(current_set))
  paste0(out, collapse=",") 
}

prepareSlurmBatchRuns = function(analysis_name, singleFunction, iterations, main_folder, out_folder, packages, ..., cpus_per_task = 2, overwrite = FALSE){
  #' Create folder
  if(dir.exists(analysis_name) & !overwrite){
    stop(sprintf("Directory analysis %s already exists and overwrite = %s", analysis_name, overwrite))
  }
  dir.create(analysis_name)
  
  #' Save function and any additional datasets to file
  save(singleFunction, ..., file = sprintf("./%s/function_data.RData", analysis_name))
  
  #' Save iterations to be used
  iterations %>% qs::qsave(sprintf("./%s/iterations.qs", analysis_name))
  
  #' Save R-script to be used in each run (calls singleFunction)
  script_iteration = sprintf('
.args = if (interactive()) c("%s/analysis/%s", 1) else commandArgs(trailingOnly = TRUE)
setwd(.args[1])
it = as.numeric(.args[2])

message(sprintf("Starting iteration %%s", it))

main_folder = "%s"
out_folder = "%s"

if(!"pacman" %%in%% rownames(installed.packages())) install.packages("pacman")
pacman::p_load(char = c("%s"))
analysis_name = "%s"

scenario_list = readRDS(sprintf("%%s/data/clean/scenarios.RDS", main_folder))
iterations = qs::qread(sprintf("iterations.qs"))

attach(sprintf("function_data.RData"))
result = singleFunction(iterations[i == it, -"i"])
result %%>%% qs::qsave(sprintf("result_%%s.qs", it))
', main_folder, analysis_name, main_folder, out_folder, paste0(packages, collapse='","'), analysis_name)
  write(script_iteration, file = sprintf("./%s/single.R", analysis_name), append = FALSE)
  
  #' Save shell-script to be used in each run (calls R-script)
  script_shell = sprintf('#!/bin/bash
RUN=$1
WD="%s/analysis/%s/"
FILE="%s/analysis/%s/single.R"
Rscript ${FILE} ${WD} ${RUN}
', main_folder, analysis_name, main_folder, analysis_name)
  write(script_shell, file = sprintf("./%s/single.sh", analysis_name), append = FALSE)
  Sys.chmod(sprintf("./%s/single.sh", analysis_name), mode = "0775", use_umask = TRUE)
  
  #' Save batch-script to be used in each run (will call shell-script for each iteration)
  script_sbatch = sprintf('#!/bin/bash
#SBATCH --job-name=analyis_%s --array=%s --cpus-per-task=%s
%s/analysis/%s/single.sh $SLURM_ARRAY_TASK_ID
', analysis_name, iterations[, i] %>% formatItRange(), cpus_per_task, main_folder, analysis_name)
  write(script_sbatch, file = sprintf("./%s/%s.sbatch", analysis_name, analysis_name), append = FALSE)
  
  #' Use this command to submit the runs
  command = sprintf("sbatch %s/analysis/%s/%s.sbatch", main_folder, analysis_name, analysis_name)
  message(sprintf("Use the following command to submit SLURM batch run:\n%s", command))
  return(command)
}

combineSlurmResults = function(analysis_name, iterations){
  result = iterations[, i] %>% lapply(function(it){
    if(!file.exists(sprintf("%s/result_%s.qs", analysis_name, it))){
      message(sprintf("File %s does not exist", sprintf("result_%s.qs", it)))
      return(NULL)
    }
    
    sprintf("%s/result_%s.qs", analysis_name, it) %>% qs::qread()}) %>%
    rbindlist
  
  return(result)
}

###

stratified_options = c("target_R0", "prior_immunity","iv_shielded_prop", "iv_shielded_gz_size")

baseline_options = c(
  iv_shielded_unshielded_contacts = 0.2, iv_shielded_shielded_contacts = 1,
  iv_shielded_symptomatic_contacts = 1, iv_shielded_mitigate = -1, iv_shielded_enter = 1)

replBasl = function(x, baseline = baseline_options){
  baseline[names(x)] = x
  return(baseline)
}
namedOptions = function(x, func){
  setNames(x, sapply(x, func))
}

mitigation_options <- list(
  "no shielding" = list("NA" = c(iv_shielded_prop = 0)),
  "reduction contacts shielded unshielded" = c(0, 0.2, 0.4, 0.6, 0.8, 1.0) %>%
    namedOptions(function(x) paste0(100 * (1 - x), "%")) %>%
    lapply(function(x) replBasl(c(iv_shielded_unshielded_contacts = x))),
  "change contacts shielded" = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2) %>%
    namedOptions(function(x) paste0(100 * x, "%")) %>%
    lapply(function(x) replBasl(c(iv_shielded_shielded_contacts = x))),
  "reduction shielded symptomatic contacts" = c(0, 0.2, 0.4, 0.6, 0.8, 1) %>%
    namedOptions(function(x) paste0(100 * (1 - x), "%")) %>%
    lapply(function(x) replBasl( c(iv_shielded_symptomatic_contacts = x))),
  "symptomatic leave after delay" = c(0, 2, 4) %>%
    namedOptions(function(x) paste0(x, " days")) %>%
    lapply(function(x) replBasl(c(iv_shielded_mitigate = x))),
  "dismantle greenzone after case" = list("NA" = replBasl(c(iv_shielded_mitigate = 99))),
  "symptomatic individuals do not enter" = list("NA" = replBasl(c(iv_shielded_enter = 2))),
  "individuals tested before they enter" = list("NA" = replBasl(c(iv_shielded_enter = 3))))

mitigation_options_plotnames <- c(
  "No shielding",
  "Relative reduction in contacts\nbetween shielded and\nunshielded individuals",
  "Relative change in contacts\nin shielded compartment",
  "Relative reduction of\nsymptomatic contacts in\nshielded compartment",
  "Symptomatic cases leave\nshielded compartment\nafter symptom onset",
  "Shielded compartment is dismantled\nafter first symptomatic case",
  "Symptomatic individuals are\nnot shielded",
  "Individuals are tested before\nentering the shielded compartment")

mitigation_options_colours = c(
  "no shielding" = "#000000",
  "reduction contacts shielded unshielded" = "#00BF6F",
  "change contacts shielded" = "#FFB81C",
  "reduction shielded symptomatic contacts" = "#FE5000",
  "symptomatic leave after delay" = "#00AEC7",
  "dismantle greenzone after case" = "#0D5257",
  "symptomatic individuals do not enter" = "#0D5257",
  "individuals tested before they enter" = "#0D5257")

#calculate x-position for elements
plotLabelledX = function(scenario_list, mitigation_options){
  large_x = 20
  small_x = 10
  minwidth = 40
  
  y = mitigation_options %>% rev %>% lapply(function(mit_options_levels){
    mit_options_levels %>% rev %>% lapply(function(x) x %>% t %>% as.data.table) %>%
      (function(x){y = rbindlist(x); y[, option_name := names(x)]})()})
  
  for(i in 1:length(y)){
    if(i == 1) y[[i]][, x := (seq_len(.N) - 1) * small_x] else
      y[[i]][, x := max(y[[i - 1]][, x]) + large_x + (seq_len(.N) - 1) * small_x]
  }
  
  z = list()
  for(i in 1:length(y)){
    shielding = TRUE
    if("iv_shielded_prop" %in% colnames(y[[i]])){
      if(y[[i]][, iv_shielded_prop] == 0){
        shielding = FALSE 
      }
    }
    
    z[[i]] = y[[i]] %>% merge(scenario_list[(iv_shielded_prop > 0) == shielding], colnames(.) %>% subset(!. %in% c("option_name", "x")), all.x = TRUE) %>%
      .[, c("x", "scen", "option_name")] %>% .[, name := names(y)[i]]
  }
  z = z %>% rbindlist
  z[option_name == "NA", option_name := ""]
  return(z)
}
