#' Calculate median prevalence by age to use in models with prior immunity

#' Process for each run
runs=1:1000
out_folder = "./output/"

model_output = lapply(1:2, function(i){
  out = sprintf("%s/out_run_%s.qs", out_folder, run) %>% qread
  
  #' Loop through each scenario, to ensure we get data for all age groups
  out$model_out[, .(infected=.N), by=c("age", "scen")] %>%
    .[, merge(out$population[, .N, by="age_group"], .SD, by.x="age_group", by.y="age", all.x=TRUE), by="scen"] %>%
    .[, infected := ifelse(is.na(infected), 0, infected)] %>%
    .[, c("prevalence", "run") := .(infected/N, run)] %>%
    return}) %>% rbindlist

#' Calculate average prevalence in the population in each run
model_output[, wprev := sum(infected)/sum(N), by=c("scen", "run")]

#' Restrict to prevalence of interest and calculate median value in all runs
target_prevs = c(0.25, 0.50)
prevalence_by_age = target_prevs %>% lapply(function(target_prev){
  model_output[wprev >= target_prev] %>%
    .[, .(prev=quantile(prevalence, 0.5)), by=c("age_group", "scen")] %>%
    merge(scenario_list[iv_shielded_prop == 0 & prior_immunity == 0, c("scen", "target_R0")] %>%
            setNames(c("scen", "r0")), by="scen") %>%
    .[, c("target_prevalence", "type") := .(target_prev, "direct")]}) %>% rbindlist %>%
  .[, c("age_group", "prev", "r0", "type", "target_prevalence")] %>%
  setNames(c("age", "prev", "r0", "type", "target_prevalence"))

saveRDS(prevalence_by_age, "./data/raw/prevalence_by_age.RDS")
