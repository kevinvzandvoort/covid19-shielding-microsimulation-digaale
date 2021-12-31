#' Get data from Zenodo
digaale_contact_data =
  readRDS("./data/raw/digaale_contact_data.RDS")
  #socialmixr::get_survey("https://zenodo.org/record/5226281#.YR-TzlvTVH6")
  #saveRDS(digaale_contact_data, "./data/raw/digaale_contact_data.RDS")

digaale_survey_population =
  readRDS("./data/raw/digaale_survey_population.RDS")
  #data.table::fread("https://zenodo.org/record/5226281/files/espicc_somaliland_digaale_survey_population.csv")
  #saveRDS(digaale_survey_population, "./data/raw/digaale_survey_population.RDS")
digaale_contact_data$participants[, dayofweek := ifelse(dayofweek == 6, 0, dayofweek + 1)]

digaale_contact_data$contacts[, intra_household_contact := contact_relationship == "Household member"]

#' The observed intra-household contact matrix
digaale_contact_matrix_intra_household = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                             filter = list(intra_household_contact = TRUE),
                             age.limits = digaale_survey_population$lower.age.limit,
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE, quiet = TRUE)

intra_hh_observed = digaale_contact_matrix_intra_household$matrix %>% t() %>% eigen(only.values = TRUE) %>% .[["values"]] %>% max

#population_data = tempfile()
#download.file("https://github.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019/raw/main/data/household_data_members.RDS",
#              population_data)
#population_data = population_data %>% readRDS %>% .[, c("id", "household_id", "age_group_80")]
#saveRDS(population_data, "./data/raw/population_data.RDS")
population_data = readRDS("./data/raw/population_data.RDS")

#' age groups to be used in the model
age_groups = data.table(
  age_group = c(1:9)-1,
  age_low = seq(0, 80, 10),
  age_high = c(seq(9, 79, 10), 100))

age_groups_names = age_groups[, .(name = ifelse(age_group == age_groups[, max(age_group)],
                                                paste0(age_groups[age_group == max(age_group), age_low], "+"),
                                                paste0(age_low,"-",age_high))), by=age_group]

population_data = population_data[, age_group := as.numeric(age_group_80) - 1] %>% .[, c("household_id", "age_group")]
population_popsize = population_data[, .N, by=age_group][order(age_group), N]

#' factor by which collected population estimates were corrected
population_fpc_correction = (digaale_survey_population[, population]/population_popsize)[1]

n_households = population_data[, unique(household_id) %>% length * population_fpc_correction] %>% round
n_agegroups = age_groups[, .N]

#' The expected intra-household contact matrix, assuming all household members make (at least) one contact per day
population_house_total = totalAgePairs(population_data)
contact_matrices_household_homogeneous = constructContactMatrix(population_house_total, population_popsize, population_popsize)

intra_hh_expected = contact_matrices_household_homogeneous$rate_adjusted %>% eigen(only.values = TRUE) %>% .[["values"]] %>% max

#' Probability any two household members come into contact per day
contact_prob_household = intra_hh_observed/intra_hh_expected

#' Observed extra household contact matrix
digaale_contact_matrix_extra_household = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                             filter = list(intra_household_contact = FALSE),
                             age.limits = digaale_survey_population$lower.age.limit,
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE, quiet = TRUE)

#' Transpose so contactors are columns and contactees are rows
contact_matrices_inside_household_empirical = digaale_contact_matrix_intra_household$matrix %>% t
contact_matrices_outside_household_empirical = digaale_contact_matrix_extra_household$matrix %>% t

#population_data <- fread("./data/clean/population_data.csv")
#contact_data <- fread("./data/clean/contact_data.csv")

#contact_matrices_household_homogeneous <- readRDS("./data/clean/contact_matrices_population_homogeneous_mixing.rds")
#contact_matrices_household_empirical <- readRDS("./data/clean/contact_matrices_empirical_household.rds")
#contact_matrices_outside_household_empirical <- readRDS("./data/clean/contact_matrices_empirical_outside_household.rds")

#contact_prob_household <- (
#    eigen(contact_matrices_household_empirical$rate_adjusted)$values[1]
#  ) / (
#    eigen(contact_matrices_household_homogeneous$rate_adjusted)$values[1])

prior_immunity_age_patterns = readRDS("./prevalence_by_age.RDS")

covid_parameters = list(
  #baseline susceptibility - changed to reach target R0 value
  "base_u" = 0.08,
  #duration of incubation period
  "dE" = list(distribution = "gamma",
              values = list("mu" = 2.5, "k" = 4)),
  #duration of pre-clinical period
  "dP" = list(distribution = "gamma",
              values = list("mu" = 2.5, "k" = 4)),
  #duration of clinical period
  "dC" = list(distribution = "gamma",
              values = list("mu" = 2.5, "k" = 4)),
  #duration of subclinical period
  "dS" = list(distribution = "gamma",
              values = list("mu" = 5.0, "k" = 4)),
  #clinical fraction by age; See: getClinicalFraction()
  "clinical_fraction" = list(values = list("age_y" = 19, "age_m" = 50, "age_o" = 68,
                                           "symp_y" = 0.037, "symp_m" = 0.3, "symp_o" = 0.65)),
  #relative infectiousness of those that are preclinical
  "fIp" = 1,
  #relative infectiousness of those that are clinical
  "fIc" = 1,
  #relative infectiousness of those that are subclinical
  "fIs" = 0.5,
  #https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-16-testing/
  "test_specificity" = 1,
  #assumption
  "test_sensitivity_preinfectious" = 0,
  #https://www.medrxiv.org/content/10.1101/2020.05.26.20112565v1.full.pdf
  #"test_sensitivity_infectious_clinical" = 0.91, 
  #"test_sensitivity_infectious_preclinical" = 0.91,
  #"test_sensitivity_infectious_subclinical" = 0.91)
  #https://www.bmj.com/content/374/bmj.n1676
  "test_sensitivity_infectious_clinical" = 0.842, 
  "test_sensitivity_infectious_preclinical" = 0.587,
  "test_sensitivity_infectious_subclinical" = 0.587)

