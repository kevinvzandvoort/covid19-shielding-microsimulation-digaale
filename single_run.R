.args = if(interactive())
  c("~/workspace/covid_ibm_shielding/", 1) else
    commandArgs(trailingOnly = TRUE)

setwd(.args[1])
run = as.numeric(.args[2])

message(sprintf("Starting run %s", run))
message("Setting up R")

if(!"pacman" %in% rownames(installed.packages()))
  install.packages("pacman")

pacman::p_load(Rcpp, data.table, distcrete, magrittr, socialmixr, qs)

source("./scripts/functions.R")
source("./scripts/read_data.R")

dir.create("./output")
dir.create("./model/build")
dir.create(sprintf("./model/build/%s", Sys.info() %>% .["nodename"]))

message("Loading model")

#' compile c++ code for model using Rcpp
#' * compile with openmp flag to increase model speed
Sys.unsetenv("PKG_CPPFLAGS")
Sys.setenv(PKG_CPPFLAGS = "-fopenmp")
sourceCpp("./model/covidIBM_shielding.cpp",
          rebuild = FALSE, cacheDir = sprintf("./model/build/%s/", Sys.info() %>% .["nodename"]), verbose = T)
#sourceCpp("~/workspace/covid_ibm_shielding/model/covidIBM.cpp", rebuild = T)
#sourceCpp("./model/covidIBMshielding_quick.cpp", rebuild = T)

#' Set a seed
seed_ = run

#Number of timesteps used per day
tstep_day = 1

#' Number of days the simulation is ran for
days = 365

#' Initial number of infected people at day 1
init_infected = 1

#' From what age should individuals be shielded (6 = 60+)
shielding_age_min = age_groups[age_low == 60, age_group]

#' Prevalence level when interventions start
prevalence_start_value = 5/1000

#' Scenarios to model
scenario_list = readRDS("./data/clean/scenarios.RDS")

#' Mitigation options needed in model run
iv_shielded_mitigate = c(delay0 = 0, delay2 = 2, delay4 = 4, none= -1, allexit = 99)

#' Prepare output lists
out_run_population = out_run_interventions = list()

message("Preparing population")
set.seed(run)

#' get bootstrapped estimate of population
bootstrapped_pop = bootstrapPopulation(n_households);
bootstrapped_pop = bootstrapped_pop %>%
  setPopDistributions(tstep_day) %>%
  setClinicalFraction(age_groups) %>%
  setSusceptibility(covid_parameters$base_u)
  
#' contact matrix for this population
bootstrap_population_popsize = bootstrapped_pop[, .N, by=age_group] %>% .[order(age_group), N]
bootstrap_sample_household_total = totalAgePairs(bootstrapped_pop)
contact_matrices_bootstrap_sample_household =
  constructContactMatrix(bootstrap_sample_household_total,bootstrap_population_popsize, bootstrap_population_popsize,
                         contact_prob_household)

#' this is the empirical matrix for contants outside of the household from Digaale
contact_matrices_bootstrap_sample_outside_household = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population %>% copy %>%
                               .[, population := bootstrap_population_popsize],
                             filter = list(intra_household_contact = FALSE),
                             age.limits = digaale_survey_population$lower.age.limit,
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE,
                             n = 1, bootstrap = TRUE, quiet = TRUE) %>% .[["matrix"]] %>% t()
  
contact_matrices_bootstrap_sample_outside_household_prob =
  contact_matrices_bootstrap_sample_outside_household / bootstrap_population_popsize
  
contact_matrices_bootstrap_sample_outside_household_prob =
  divideContactProbabilities(contact_matrices_bootstrap_sample_outside_household_prob, tstep_day)

message("Starting scenario runs")
scens = scenario_list[, scen]
for(i in seq_len(length(scens))){
  message(sprintf("Run %s/%s", i, length(scens)))
  set.seed(run+1)
  
  bootstrapped_pop = bootstrapped_pop %>%
    setPriorImmunity(scenario_list[i, prior_immunity], scenario_list[i, target_R0])
  
  #' Calculate R0 (assuming a compartmental model structure) and calculate u to reach the target_R0 value  
  u = calculateSusceptibility(matrices = list(contact_matrices_bootstrap_sample_household$rate_adjusted,
                                              contact_matrices_bootstrap_sample_outside_household),
                              scenario_list[i, target_R0])
  bootstrapped_pop = setSusceptibility(bootstrapped_pop, u)
    
  model_params = list(
    n_agegroups = n_agegroups,
    population = bootstrapped_pop,
    n_households = n_households,
    n_init_infected = init_infected,
    days = days,
    tstep_day = tstep_day,
    contact_household = 1 - (1 - contact_prob_household) ^ (1 / tstep_day),
    contact_matrix = contact_matrices_bootstrap_sample_outside_household_prob,
    interventions = list(),
    test_specificity = covid_parameters$test_specificity,
    test_sensitivity_preinfectious = covid_parameters$test_sensitivity_preinfectious,
    test_sensitivity_preclinical = covid_parameters$test_sensitivity_infectious_preclinical,
    test_sensitivity_subclinical = covid_parameters$test_sensitivity_infectious_subclinical,
    test_sensitivity_clinical = covid_parameters$test_sensitivity_infectious_clinical)
  
  if(scenario_list[i, iv_shielded_prop] != 0){
    interventions = list(
      shielding = list(
        start_event = c("time", "prevalence")[2], start_value = prevalence_start_value,
        processed = 0, start_time = -1,
        shielding_group_size = scenario_list[i, iv_shielded_gz_size],
        shielding_age_min = shielding_age_min,
        shielding_proportion = scenario_list[i, iv_shielded_prop],
        shielding_enter_control = scenario_list[i, iv_shielded_enter],
        shielding_mitigate_type = names(iv_shielded_mitigate) %>%
          subset(iv_shielded_mitigate == scenario_list[i, iv_shielded_mitigate]) %>%
          (function(x) ifelse(grepl("delay", x, fixed = TRUE), "delay", x))(),
        shielding_mitigate_value = scenario_list[i, iv_shielded_mitigate],
        contacts_shielded_shielded = c(gz = scenario_list[i, iv_shielded_shielded_contacts], nongz = 0), #gz are members in same greenzone
        contacts_shielded_unshielded = scenario_list[i, iv_shielded_unshielded_contacts],
        contacts_unshielded_shielded = scenario_list[i, iv_shielded_unshielded_contacts]),
      self_isolation = list(
        start_event = c("time", "prevalence")[2], start_value = prevalence_start_value,
        processed = 0, start_time = -1,
        contacts_unshielded_infectious_clinical = c(hh = 1, nonhh = 1),
        contacts_shielded_infectious_clinical = c(gz = scenario_list[i, iv_shielded_symptomatic_contacts],
                                                  nongz = scenario_list[i, iv_shielded_symptomatic_contacts])))
  } else {
    interventions = list()
  }
  model_params[["interventions"]] = interventions
  
  output = covidIBM(model_params, seed_)
  out_population = output$population %>% as.data.table %>% .[, scen := scens[i]]
  
  if(length(interventions) == 0){
    out_intervention = data.table(processed = -1, start_time = -1, intervention = NA_character_, scen = i)
  } else {
    out_intervention = output$intervention %>% lapply(as.data.table) %>% rbindlist %>% 
      .[, c("intervention", "scen") := .(names(model_params$interventions), scens[i])]  
  }
  
  out_intervention %>% .[]
  
  out_run_population[[length(out_run_population) + 1]] = out_population
  out_run_interventions[[length(out_run_interventions) + 1]] = out_intervention
}

#save the run
message("Processing and saving output")

out_run = list(
  model_out = rbindlist(out_run_population),
  model_interventions = rbindlist(out_run_interventions),
  population = bootstrapped_pop[, -"immune"])

qs::qsave(out_run, sprintf("./output/out_run_%s.qs", run))

message("Finished")
