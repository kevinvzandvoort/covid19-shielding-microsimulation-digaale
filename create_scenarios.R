#' target R0s per simulation
target_R0 = c(1.5, 2.5, 3.5, 5)

#' Proportion of people aged >60 who are shielded
iv_shielded_prop = c(1, 0.8, 0.6, 0.4, 0.2, 0)

#' Size of each shielding compartment
iv_shielded_gz_size <- c(0, 2, 4, 8, 16)

#' Options to mitigate infections when individuals enter shielded compartment
iv_shielded_enter <- c(all = 1, nosymptomatic = 2, test = 3)

#' Ratio adjusting the contacts between shielded and unshielded contacts
iv_shielded_unshielded_contacts <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

#' Ratio adjusting the contacts between shielded contacts
iv_shielded_shielded_contacts <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

#' Ratio adjusting the contact rate of symptomatic individuals who are shielded
iv_shielded_symptomatic_contacts <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

#' Mitigation option when an individual with symptoms is detected in a shielded space
iv_shielded_mitigate <- c(delay0 = 0, delay2 = 2, delay4 = 4, none= -1, allexit = 99)

#' Immunity in the population before starting the simulation
prior_immunity = c(0, 0.25, 0.5)

base_scenario = data.table(
  iv_shielded_enter = 1,
  iv_shielded_unshielded_contacts = 0.2,
  iv_shielded_shielded_contacts = 1,
  iv_shielded_symptomatic_contacts = 1,
  iv_shielded_mitigate = -1)

variables = c("iv_shielded_enter", "iv_shielded_unshielded_contacts", "iv_shielded_shielded_contacts",
              "iv_shielded_symptomatic_contacts", "iv_shielded_mitigate")

scenario_list = variables %>% lapply(function(v){
  get(v) %>%
    subset(!. %in% base_scenario[, get(v)]) %>%
    lapply(function(x){
      base_scenario %>% copy %>% .[, {v} := x]}) %>%
    rbindlist}) %>% rbindlist

scenario_list = rbind(base_scenario, scenario_list)
scenario_list[, sub_scen := 1:.N]

#' Duplicate for all shielded_prop, shielded_gz_size, target_r0, and prior_immunity considered
scenario_list = expand.grid(iv_shielded_prop = iv_shielded_prop,
                            iv_shielded_gz_size = iv_shielded_gz_size,
                            target_R0 = target_R0,
                            prior_immunity = prior_immunity) %>% as.data.table %>%
  .[!(iv_shielded_prop == 0 & iv_shielded_gz_size > 0) & !(iv_shielded_gz_size == 0 & iv_shielded_prop > 0)] %>%
  .[, main_scen := 1:.N] %>%
  .[, cbind(.SD, scenario_list), by="main_scen"] %>%
  .[!(iv_shielded_prop == 0 & sub_scen > 1)] %>%
  .[, scen := 1:.N] %>%
  .[, c("scen", "main_scen", "sub_scen", colnames(.) %>%
          subset(!. %in% c("scen", "main_scen", "sub_scen"))), with=F]

#save scenario list
saveRDS(scenario_list, "./data/clean/scenarios.RDS")
