analysis_name = "breached_mitigation"
packages_func = c("qs", "magrittr", "data.table")
packages_other = c("ggplot2")

#' Set global plot styling options
page_outer_margin = 0.523 #inch
page_inner_margin = page_outer_margin/2 #inch
page_width = 8 #inch
page_height = 11 #inch
plot_single_col = (page_width - page_outer_margin*2 - page_inner_margin)/2 #inch
plot_double_col = page_width - page_outer_margin*2 #inch

#' Additional ggplot theme options to be used when combining plots (utilizing patchwork)
theme_multiplot = ggplot2::theme(legend.title = element_text(size=9), legend.text = element_text(size = 8),
                                 legend.box.background = element_rect(fill="#FFFFFF", colour="#000000", size = 0.2),
                                 legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
                                 legend.key.size = unit(0.3, 'cm'), axis.title = element_text(size=10),
                                 axis.text = element_text(size = 8),
                                 plot.title = element_text(size=10, face = "bold", hjust=0.62),
                                 plot.title.position = "plot", plot.tag = element_text(size=10),
                                 plot.tag.position = c(0, 0.99), strip.text = element_text(face = "plain", size=9))

pacman::p_load(char = c(packages_func, packages_other))
source("./helper_functions.R")

main_folder = c("/nfs/general/covid_ibm_shielding", "~/workspace/covid_ibm_shielding")[1]
scenario_list = readRDS(sprintf("%s/data/clean/scenarios.RDS", main_folder))
#populations = readRDS(sprintf("%s/modelled_populations.RDS", main_folder))
#populations_60p = populations[age_group >= 6, .(N=sum(N)), by="run"]

out_folder = sprintf("%s/output", main_folder)

#' iterations for singleFunction (column names should correspond with variables used in function)
target_R0_values = scenario_list[, unique(target_R0)]
target_prior_immunity_values = scenario_list[, unique(prior_immunity)]
iterations = expand.grid(r0 = target_R0_values,
                         immunity = target_prior_immunity_values) %>% as.data.table %>% .[, i := seq_len(.N)]

singleFunction = function(iteration_data){
  iteration_data %>% attach()
  
  message("Get required scenarios")
  #' Only scenarios that affect transmission into shielding compartments
  scens = scenario_list[target_R0 == r0 & prior_immunity == immunity &
                          iv_shielded_prop > 0 & iv_shielded_shielded_contacts == 1 &
                          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_unshielded_contacts == 0.2, scen]
  interventions = scens %>% lapply(function(s){
    sprintf("%s/byscen/interventions/byscen/out_scen_%s_model_interventions.qs", out_folder, s) %>% qread})  %>% rbindlist
  
  message("Limit to runs of interest")
  #' Limited to runs in scenarios where shielding was implemented
  scens_of_interest = interventions[intervention == "shielding" & processed == 1, c("scen", "run")]  
  if(nrow(scens_of_interest) == 0){
    message("No scenarios applicable")
    return(NULL)
  }
  
  scens_of_interest = scens_of_interest %>% merge(scenario_list[, c("scen", "main_scen", "sub_scen")], by="scen",
                                                  all.x=TRUE)
  
  message("Get data for each scenario")
  scens_nomitigation = scens_of_interest[sub_scen == 1, unique(scen)]
  
  model_out_nomitigation = seq_len(length(scens_nomitigation)) %>% lapply(function(s, scens){
    message(sprintf("%s/%s", s, length(scens)))
    model_out = sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, scens[s]) %>% qread %>%
      .[run %in% scens_of_interest[scen == scens[s], run]]
    #' Calculate proportion breached
    breached_group = model_out[age >= 6 & shielded_gz_id > -1, .(infected = sum(infected_shielded),
                                                                 breached = sum(infected_shielded) > 0),
                               by=c("scen", "run", "shielded_gz_id")] %>%
      .[, .(breached = sum(breached), model_max_gz = (max(shielded_gz_id) + 1)), by=c("scen", "run")] %>%
      .[, prop_breached := breached/model_max_gz] %>%
      .[, main_scen := scenario_list[scen == scens[s], main_scen]]
    
    return(breached_group)
  }, scens_nomitigation) %>% rbindlist()
  
  scens_mitigation = scens_of_interest[sub_scen > 1, unique(scen)]
  model_out_mitigation = seq_len(length(scens_mitigation)) %>% lapply(function(s, scens){
    message(sprintf("%s/%s", s, length(scens)))
    model_out = sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, scens[s]) %>% qread %>%
      .[run %in% scens_of_interest[scen == scens[s], run]]
    #' Calculate proportion breached
    breached_group = model_out[age >= 6 & shielded_gz_id > -1, .(infected = sum(infected_shielded),
                                                                 breached = sum(infected_shielded) > 0),
                               by=c("scen", "run", "shielded_gz_id")] %>%
      .[, .(breached = sum(breached), model_max_gz = (max(shielded_gz_id) + 1)), by=c("scen", "run")] %>%
      .[, prop_breached := breached/model_max_gz] %>%
      .[, main_scen := scenario_list[scen == scens[s], main_scen]]
    return(breached_group)
    
    }, scens_mitigation) %>% rbindlist
  
  model_out = model_out_mitigation[, c("run", "prop_breached", "main_scen", "scen")] %>%
    merge(model_out_nomitigation[, c("run", "prop_breached", "main_scen")], by=c("main_scen", "run"), all.x=TRUE) %>%
    .[, c("abs_diff", "rel_diff") := .(prop_breached.y - prop_breached.x, prop_breached.x/prop_breached.y)] %>%
    melt(measure.vars = c("abs_diff", "rel_diff")) %>%
    .[, .(median = median(value, na.rm = TRUE),
          mean = mean(value = TRUE, na.rm = TRUE),
          low95 = quantile(value, 0.025, na.rm = TRUE),
          high95 = quantile(value, 0.975, na.rm = TRUE),
          low50 = quantile(value, 0.25, na.rm = TRUE),
          high50 = quantile(value, 0.75, na.rm = TRUE)), by=c("variable", "scen")]
  
  return(model_out)
}

#sbatch_command = prepareSlurmBatchRuns(analysis_name, singleFunction, iterations, main_folder, out_folder, packages_func)
result = combineSlurmResults(analysis_name, iterations)

result[variable == "rel_diff"] %>% merge(scenario_list, by="scen") %>% .[iv_shielded_prop == 0.4] %>% ggplot(aes(x=iv_shielded_gz_size, y=median))+facet_grid(iv_shielded_enter~target_R0+prior_immunity)+geom_point()

#' In none of the scenarios did the mitigation options have a significant impact on the proportion of breached greenzones