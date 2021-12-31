#' Calculate probability of clinical disease by age group
#' * following posterior estimates by Davies et al
#' * source: https://www.nature.com/articles/s41591-020-0962-9
getClinicalFraction <- function(age_groups){
  #' smoothly interpolate between points (x0, y0) and (x1, y1) using cosine interpolation.
  #' for x < x0, returns y0; for x > x1, returns y1; for x0 < x < x1, returns the cosine interpolation between y0 and y1
  interpolate_cos = function(x, x0, y0, x1, y1)
  {
    ifelse(x < x0, y0, ifelse(x > x1, y1, y0 + (y1 - y0) * (0.5 - 0.5 * cos(pi * (x - x0) / (x1 - x0)))))
  }
  
  age_groups[, mid := mean(c(age_low, age_high)), by=seq_len(nrow(age_groups))]
  
  age_y = covid_parameters[["clinical_fraction"]][["values"]][["age_y"]]
  age_m = covid_parameters[["clinical_fraction"]][["values"]][["age_m"]]
  age_o = covid_parameters[["clinical_fraction"]][["values"]][["age_o"]]
  
  # definition of "young", "middle", and "old"
  young  = interpolate_cos(age_groups[, mid], age_y, 1, age_m, 0);
  old    = interpolate_cos(age_groups[, mid], age_m, 0, age_o, 1);
  middle = 1 - young - old;
  
  symp_y = covid_parameters[["clinical_fraction"]][["values"]][["symp_y"]]
  symp_m = covid_parameters[["clinical_fraction"]][["values"]][["symp_m"]]
  symp_o = covid_parameters[["clinical_fraction"]][["values"]][["symp_o"]]
  
  return(young * symp_y + middle * symp_m + old * symp_o)
}

setClinicalFraction <- function(population_sample, age_groups){
  init_colnames <- colnames(population_sample)
  age_groups[, "y"] <- getClinicalFraction(age_groups)
  population_sample <- merge(population_sample, age_groups[, c("age_group", "y")], by="age_group")
  setorder(population_sample, participant_id)
  
  return(population_sample[, c(init_colnames, "y"), with=FALSE])
}

#' Set susceptibility to infection
setSusceptibility <- function(population_sample, u=covid_parameters$base_u){
  population_sample[, "u"] <- u
  
  return(population_sample)
}

#' Get random values from discrete gamma functions
#' * paramerised using mu and k
rgammaAltDiscrete <- function(n, mu, k, tstep=1){
  distcrete::distcrete("gamma", tstep, shape = k, scale = mu/k)$r(n)
}

#' Get bootstrapped estimate of population
#' * requires global dataset: population_data (source UNWPP)
bootstrapPopulation = function(n_households, population_data. = population_data){
  household_sample = sample(unique(population_data[, household_id]), n_households, TRUE)
  population_sample = seq_len(n_households) %>% lapply(function(x, household_sample){
    population_data[household_id == household_sample[x]] %>% copy %>%
      .[, household_id := x] %>% return}, household_sample) %>% rbindlist
  
  setorder(population_sample, household_id)
  population_sample[, participant_id := seq_len(.N)]
  
  return(population_sample[, c("participant_id", "household_id", "age_group")])
}

#' Get distribution to use
getRandDist <- function(parameter, n=1, tstep=1){
  p <- covid_parameters[[parameter]]
  if("distribution" %in% names(p)){
    
    return(switch(
      p[["distribution"]],
      "gamma" = rgammaAltDiscrete(n, p[["values"]][["mu"]], p[["values"]][["k"]], tstep),
      NULL
    ))
  } else {
    stop(sprintf("No distribution set for parameter %s", parameter))
  }
}

#' Sample random values for distributions
#' * This will be the duration an individual will be in a compartment, if they would get there
#' * dE: incubation period
#' * dP: duration of pre-clinical period
#' * dC: duration of clinical period
#' * dS: duration of subclinical period
setPopDistributions <- function(population_sample, tstep=1){
  population_sample[, "dE"] <- tstep * getRandDist("dE", nrow(population_sample), 1/tstep)
  population_sample[, "dP"] <- tstep * getRandDist("dP", nrow(population_sample), 1/tstep)
  population_sample[, "dC"] <- tstep * getRandDist("dC", nrow(population_sample), 1/tstep)
  population_sample[, "dS"] <- tstep * getRandDist("dS", nrow(population_sample), 1/tstep)
  
  return(population_sample)
}

#' Set prior immunity for each individual
#' * The probability for any individual of age a depends on the selected R0 and prior_immunity
#' * If prior_immunity == 0, function will return population table with no individual susceptible
setPriorImmunity = function(population_sample, prior_immunity, R0, age_pattern = prior_immunity_age_pattern,
                            age_group_names. = age_groups_names){
  if(prior_immunity == 0){
    population_sample[, immune := 0]
    return(population_sample)
  }
  
  age_pattern = prior_immunity_age_patterns[r0 == R0 & target_prevalence == prior_immunity]
  
  if(nrow(age_pattern) == 0)
    stop("No age pattern data for parameters requested")
  
  age_pattern %>% merge(age_groups_names, by.x="age", by.y="name") %>% .[, c("age_group", "prev")] %>%
    merge(population_sample, by="age_group") %>%
    .[, immune := rbinom(1, 1, prev), by=participant_id] %>% .[, -"prev"] %>%
    .[, c("household_id", "participant_id", "age_group", colnames(.) %>%
            subset(!. %in% c("household_id", "participant_id", "age_group"))), with=FALSE] %>%
    setorder(household_id, participant_id, age_group) %>% return()
}

#' Calculate total number of age-pairs in each household included
#' * used to create within-household contact matrix corresponding with household structures
#' * account for not being able to contact yourself
totalAgePairs = function(households){
  households = households[, c("household_id", "age_group")]
  
  household_agepairs_template = matrix(0, n_agegroups, n_agegroups)
  n_byage_template = numeric(n_agegroups)
  
  out = unique(households[, household_id]) %>% lapply(function(x){
    y = households[household_id == x]
    
    household_agetable = y[, table(age_group)]
    household_agepairs = household_agepairs_template
    household_agepairs[household_agetable %>% names %>% as.numeric + 1, ] = household_agetable %>% as.numeric
    diag(household_agepairs) = pmax(diag(household_agepairs) - 1, 0)
    
    n_byage = n_byage_template
    n_byage[household_agetable %>% names %>% as.numeric + 1] = household_agetable %>% as.numeric
    household_agepairs = household_agepairs %>% sweep(2, n_byage, "*")
    
    return(household_agepairs)})
  
  out = Reduce("+", out)
  
  return(out)
}





#' Construct contact matrix from total number of contacts by age
constructContactMatrix = function(total_contacts, sample_popsize, population_popsize, adjust_ratio=1){
  rate_raw = total_contacts
  for(j in 1:ncol(rate_raw)){
    rate_raw[, j] = total_contacts[, j] / sample_popsize[j]
  }
  
  rate_adjusted = rate_raw
  for(j in 1:ncol(rate_adjusted)){
    for(i in 1:nrow(rate_adjusted)){
      rate_adjusted[i, j] = (rate_raw[i, j] + rate_raw[j, i] * population_popsize[i] / population_popsize[j])/2
    }
  }
  rate_adjusted = rate_adjusted * adjust_ratio
  
  prob = rate_adjusted / population_popsize
  #' cannot contact oneself, so need to subtract N with 1 for i == j
  diag(prob) = diag(rate_adjusted) / (population_popsize - 1)
  
  list(total = total_contacts, rate_raw = rate_raw, rate_adjusted = rate_adjusted, prob = prob) %>% return()
}

#' Calculate R0 if this would be a compartmental model
#' * Create Next Generation Matrix
#' * Differs with IBM as it does not constrain contacts within the household, but uses
#'   the same average number of contacts as used in the IBM
calculateR0cmm <- function(matrices, covid_parameters. = covid_parameters, age_groups. = age_groups) {
  dIp = covid_parameters$dP$values$mu
  dIs = covid_parameters$dC$values$mu
  dIa = covid_parameters$dS$values$mu
  
  y = getClinicalFraction(age_groups)
  
  cm = Reduce("+", matrices)
  
  ngm = covid_parameters$base_u * t(t(cm) * (
    y * (covid_parameters$fIp * dIp + covid_parameters$fIc * dIs) + 
      (1 - y) * covid_parameters$fIs * dIa)
  )
  
  return(abs(eigen(ngm)$values[1]))
}

#' Calculate new susceptibility to reach target R0
#' * using compartmental model structure
calculateSusceptibility = function(matrices, target_R0, covid_parameters. = covid_parameters,
                                   age_groups. = age_groups){
  (covid_parameters$base_u * (target_R0 / calculateR0cmm(matrices, covid_parameters))) %>% return()
}

divideContactProbabilities = function(probability_matrix, by=1){
  if(by%%1 > 0 | by < 1){ stop("by needs to be a positive integer") }
  
  for(i in 1:nrow(probability_matrix)){
    for(j in 1:ncol(probability_matrix)){
      probability_matrix[i, j] = 1 - (1 - probability_matrix[i, j]) ^ (1 / by)
    }
  }
  
  return(probability_matrix)
}
