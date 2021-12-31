analysis_name = "risk_bytime_full"
packages_func = c("qs", "magrittr", "data.table")
packages_other = c("ggplot2")

pacman::p_load(char = c(packages_func, packages_other))
source("./helper_functions.R")

main_folder = "/nfs/general/covid_ibm_shielding"
#main_folder = "~/workspace/covid_ibm_shielding/"
scenario_list = readRDS(sprintf("%s/data/clean/scenarios.RDS", main_folder))
populations = readRDS(sprintf("%s/modelled_populations.RDS", main_folder))

populations_60p = populations[age_group >= 6, .(N=sum(N)), by="run"]
out_folder = sprintf("%s/output", main_folder)

#' iterations for singleFunction (column names should correspond with variables used in function)
target_R0_values = scenario_list[, unique(target_R0)]
target_prior_immunity_values = scenario_list[, unique(prior_immunity)]
iv_shielded_prop_values = scenario_list[, unique(iv_shielded_prop)]
iterations = expand.grid(r0 = target_R0_values,
                         immunity = target_prior_immunity_values,
                         shielded_prop = iv_shielded_prop_values) %>% as.data.table %>% .[, i := seq_len(.N)]

singleFunction = function(iteration_data){
  iteration_data %>% attach()
  
  message("Get required scenarios")
  
  #scens = scenario_list[target_R0 == r0 & prior_immunity == immunity & iv_shielded_prop == shielded_prop &
  #                        iv_shielded_shielded_contacts == 1 & iv_shielded_symptomatic_contacts == 1 &
  #                        iv_shielded_mitigate == -1 & iv_shielded_enter == 1, scen]
  scens = scenario_list[target_R0 == r0 & prior_immunity == immunity & iv_shielded_prop == shielded_prop, scen]
  interventions = scens %>% lapply(function(s){
    sprintf("%s/byscen/interventions/byscen/out_scen_%s_model_interventions.qs", out_folder, s) %>% qread})  %>% rbindlist
  
  message("Limit to runs of interest")
  
  #' Limited to runs in scenarios where shielding was implemented
  scens_of_interest = interventions[is.na(intervention) | (intervention == "shielding" & processed == 1), c("scen", "run")]  
  if(nrow(scens_of_interest) == 0){
    message("No scenarios applicable")
    return(NULL)
  }
  
  message("Get data for each scenario")
  scens = scens_of_interest[, unique(scen)]
  model_out = seq_len(length(scens)) %>% lapply(function(s, scens){
    message(sprintf("%s/%s", s, length(scens)))
    model_out = sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, scens[s]) %>% qread %>%
      .[run %in% scens_of_interest[scen == scens[s], run]]
    
    #' Exclude runs without outbreak (at least 50 cases)
    run_with_outbreak = model_out[, .(I = sum(infected_at > 0)), by="run"] %>% .[I > 50, run]
    model_out = model_out[run %in% run_with_outbreak]
    if(nrow(model_out) == 0){
      message(sprintf("No runs with outbreak for scenario %s", scens[s]))
      return(NULL)
    }
    
    #' Calculate prevalence over time
    times = 1:200
    prevalence_by_time_by_age = times %>% lapply(function(day, data){
      data[, .(time = day, Rx = sum(recovered_at <= day & infected_at >= 0)),
           by = c("scen", "run", "age", "shielded_status")]},
      data = model_out) %>%
      rbindlist() %>%
      merge(populations, by.x=c("run", "age"), by.y=c("run", "age_group"))
    
    prevalence_by_time_60m = prevalence_by_time_by_age %>%
      .[age < 6, .(prevalence=sum(Rx)/sum(N)), by=c("run", "scen", "time", "shielded_status")] %>%
      .[, age := "<60"]
    prevalence_by_time_60p = prevalence_by_time_by_age %>%
      .[age >= 6, .(prevalence = sum(Rx)/sum(N * ifelse(shielded_status, shielded_prop, 1 - shielded_prop))),
        by=c("run", "scen", "time", "shielded_status")] %>% .[, age := ">=60"] %>%
      .[, prevalence := prevalence %>% pmin(1)] #We use estimated proportion immune, so calculated prevalence could be >1 for some (few) runs
    
    prevalence_by_time = rbind(prevalence_by_time_60m, prevalence_by_time_60p)
    
    #' Calculate mean/median/quantile across all runs in data
    summary_all = prevalence_by_time %>% .[, .(mean = mean(prevalence),
                                               median = median(prevalence),
                                               low95 = quantile(prevalence, 0.025),
                                               high95 = quantile(prevalence, 0.975),
                                               low50 = quantile(prevalence, 0.25),
                                               high50 = quantile(prevalence, 0.75)), by = c("scen", "time", "shielded_status", "age")] %>%
      melt(measure.vars = c("mean", "median", "low95", "high95", "low50", "high50"))
    
    #' Get runs that are closest to the mean/median/quantile run based on final attack rate in the (unshielded) <60yo
    prevalence_by_time_forchecks = prevalence_by_time[age == "<60" & time == 200]
    match_checks = prevalence_by_time_forchecks[, c("mean", "median", "low95", "high95", "low50", "high50") :=
                                                  .(prevalence - mean(prevalence), prevalence - median(prevalence),
                                                    prevalence - quantile(prevalence, 0.025), prevalence - quantile(prevalence, 0.975),
                                                    prevalence - quantile(prevalence, 0.25), prevalence - quantile(prevalence, 0.75))] %>%
      melt(measure.vars=c("mean", "median", "low95", "high95", "low50", "high50")) %>%
      .[, first(.SD[abs(value) == min(abs(value))]), by="variable"] %>% .[, c("variable", "run")]
    
    summary_match = prevalence_by_time %>% merge(match_checks, by="run", all.y=TRUE, allow.cartesian = T)
    
    #' Combine in one dataset
    summary_data = summary_match[, c("value", "type") := .(prevalence, "match")] %>% .[, -c("prevalence", "run")] %>%
      rbind(summary_all %>% .[, type := "all"])
    return(summary_data)}, scens) %>% rbindlist
  return(model_out)
}

#sbatch_command = prepareSlurmBatchRuns(analysis_name, singleFunction, iterations, main_folder, out_folder, packages_func, populations)
result = combineSlurmResults(analysis_name, iterations)

risk_bytime_colours = c(
  "Unmitigated: individuals aged 60+y" = "#FFB81C", "Shielding: individuals aged <60y who are not shielded" = "#0D5257",
  "Shielding: individuals aged 60+y who are not shielded" = "#00BF6F",
  "Shielding: individuals aged 60+y who are shielded" = "#00AEC7")

plotdat = scenario_list %>%
  plotLabelledX(mitigation_options[c("no shielding", "reduction contacts shielded unshielded")]) %>% #"change contacts shielded")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE, all.x = TRUE) %>%
  .[target_R0 == 2.5 & prior_immunity == 0]

plotdat_sub = plotdat %>%
  .[type == c("all", "match")[1] & iv_shielded_unshielded_contacts %in% c(0.2) & iv_shielded_prop %in% c(0, 0.2, 0.8)]

plotdat_sub = plotdat_sub[iv_shielded_prop == 0, -c("iv_shielded_prop", "iv_shielded_gz_size", "shielded_status",
                                                    "iv_shielded_unshielded_contacts")] %>%
  merge(unique(plotdat_sub[iv_shielded_prop > 0, c("iv_shielded_prop", "iv_shielded_gz_size", "time", "age", "variable", "value",
                                         "shielded_status", "iv_shielded_unshielded_contacts")]),
        by=c("age", "time", "variable"), all.y = TRUE, suffixes = c("_unmit", "_mit")) %>%
  melt(measure.vars=c("value_unmit", "value_mit"), variable.name = "scen_type") %>%
  dcast(...~variable) %>%
  .[, clrname := factor(paste0(scen_type, age, ": ", factor(shielded_status, 0:1, c("unshielded", "shielded"))),
                        c("value_unmit>=60: unshielded", "value_mit<60: unshielded", "value_mit>=60: unshielded",
                          "value_mit>=60: shielded"),
                        c("Unmitigated: individuals aged 60+y", "Shielding: individuals aged <60y who are not shielded",
                          "Shielding: individuals aged 60+y who are not shielded",
                          "Shielding: individuals aged 60+y who are shielded"))] %>% .[time <= 100]

x11(width = plot_double_col, height = plot_single_col)
png(filename = "~/workspace/covid_ibm_shielding/analysis/risk_bytime_full.png", width = plot_double_col, height = plot_single_col, units = "in", res = 300)

plotdat_sub %>% .[, iv_shielded_prop := factor(iv_shielded_prop, c(0.8, 0.2))]

ggplot(mapping=aes(x=time, y=median, ymin=low95, ymax=high95, colour=clrname, fill=clrname, group=clrname))+
  labs(x = "Days since first infectious individual",
       colour="Group", fill="Group", linetype="Group", y = "Cumulative attack rate")+
  coord_cartesian(ylim=c(0, max(plotdat_sub$high95)), xlim=c(0, 100), clip = "off")+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(factor(iv_shielded_prop, c(0.8, 0.2), c("80%", "20%")) ~ iv_shielded_gz_size, labeller = labeller(
    iv_shielded_prop = function(x) paste0(as.numeric(x)*100,"%")))+
  geom_text(data = data.table(iv_shielded_prop = 0.8, iv_shielded_gz_size = 2, label = "Shielded residence size", x = 0,
                              y = 1.35), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL, linetype=NULL,
                                             fill=NULL, group=NULL), hjust = 0, vjust = 0, size = 8/(72.27/25.4),
            show.legend = FALSE)+
  geom_text(data = data.table(iv_shielded_prop = 0.8, iv_shielded_gz_size = 16, label = "Proportion 60+y shielded",
                              x = 130, y = 1), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL,
                                                   linetype=NULL, fill=NULL, group=NULL), hjust = 0, vjust = 0,
            angle = -90, size = 8/(72.27/25.4), show.legend = FALSE)+
  geom_ribbon(data = plotdat_sub[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], alpha=0.2,
              colour=NA)+
  geom_ribbon(data = plotdat_sub[scen_type == "value_mit"], alpha=0.2, colour=NA)+
  geom_line(data = plotdat_sub[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], size=1)+
  geom_line(data = plotdat_sub[scen_type == "value_mit"], aes(linetype = age), size=1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_none())+
  scale_colour_manual(values=risk_bytime_colours)+
  scale_fill_manual(values=risk_bytime_colours)+
  scale_linetype_manual(values = c("<60" = 2, ">=60" = 1))+
  theme_bw()+
  theme(plot.margin = unit(c(1.1, 1.2, 3, 0.2), "lines"), axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5), axis.title.x = element_text(size = 8), legend.position = c(0.5, -0.32),
        strip.background = element_rect(fill="#FFFFFF"), strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(), legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8), legend.box.margin = margin(0.2, 0.2, 0.2, 0.2),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2), legend.direction = "horizontal")
  
dev.off()  

plotdat = scenario_list %>%
  plotLabelledX(mitigation_options[c("no shielding", "reduction contacts shielded unshielded")]) %>% #"change contacts shielded")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE, all.x = TRUE) %>%
  .[prior_immunity == 0]

plotdat_sub = plotdat %>%
  .[type == c("all", "match")[1] & iv_shielded_unshielded_contacts %in% c(0.2) & iv_shielded_prop %in% c(0, 0.2, 0.8)]

plotdat_sub = plotdat_sub[iv_shielded_prop == 0, -c("iv_shielded_prop", "iv_shielded_gz_size", "shielded_status",
                                                    "iv_shielded_unshielded_contacts")] %>%
  merge(unique(plotdat_sub[iv_shielded_prop > 0, c("iv_shielded_prop", "iv_shielded_gz_size", "time", "age", "variable", "value",
                                                   "shielded_status", "iv_shielded_unshielded_contacts", "target_R0")]),
        by=c("age", "time", "variable", "target_R0"), all.y = TRUE, suffixes = c("_unmit", "_mit")) %>%
  melt(measure.vars=c("value_unmit", "value_mit"), variable.name = "scen_type") %>%
  dcast(...~variable) %>%
  .[, clrname := factor(paste0(scen_type, age, ": ", factor(shielded_status, 0:1, c("unshielded", "shielded"))),
                        c("value_unmit>=60: unshielded", "value_mit<60: unshielded", "value_mit>=60: unshielded",
                          "value_mit>=60: shielded"),
                        c("Unmitigated: individuals aged 60+y", "Shielding: individuals aged <60y who are not shielded",
                          "Shielding: individuals aged 60+y who are not shielded",
                          "Shielding: individuals aged 60+y who are shielded"))] %>% .[time <= 150 & iv_shielded_gz_size %in% c(2,16)]

plotdat_sub_a = plotdat_sub[iv_shielded_prop == 0.2]
plotdat_sub_b = plotdat_sub[iv_shielded_prop == 0.8]

x11(width = plot_single_col, height = plot_double_col)

plot_r0_sens_a = ggplot(mapping=aes(x=time, y=median, ymin=low95, ymax=high95, colour=clrname, fill=clrname, group=clrname))+
  labs(x = "Days since first infectious individual",
       colour="Group", fill="Group", linetype="Group", y = "Cumulative attack rate")+
  coord_cartesian(ylim=c(0, max(plotdat_sub_a$high95)), xlim=c(0, 140), clip = "off")+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(target_R0~iv_shielded_gz_size, labeller = labeller(
    iv_shielded_prop = function(x) paste0(as.numeric(x)*100,"%")))+
  geom_text(data = data.table(target_R0 = 1.5, iv_shielded_gz_size = 2, label = "Shielded residence size", x = 0,
                              y = 1.30), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL, linetype=NULL,
                                             fill=NULL, group=NULL), hjust = 0, vjust = 0, size = 8/(72.27/25.4),
            show.legend = FALSE)+
  geom_text(data = data.table(target_R0 = 1.5, iv_shielded_gz_size = 16, label = "R[0]",
                              x = 184, y = 1), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL,
                                                   linetype=NULL, fill=NULL, group=NULL), hjust = 0, vjust = 0,
            angle = -90, size = 8/(72.27/25.4), show.legend = FALSE, parse = TRUE)+
  geom_ribbon(data = plotdat_sub_a[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], alpha=0.2,
              colour=NA)+
  geom_ribbon(data = plotdat_sub_a[scen_type == "value_mit"], alpha=0.2, colour=NA)+
  geom_line(data = plotdat_sub_a[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], size=1)+
  geom_line(data = plotdat_sub_a[scen_type == "value_mit"], aes(linetype = age), size=1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_none())+
  scale_colour_manual(values=risk_bytime_colours)+
  scale_fill_manual(values=risk_bytime_colours)+
  scale_linetype_manual(values = c("<60" = 2, ">=60" = 1))+
  theme_bw()+
  theme(plot.margin = unit(c(1.1, 1.2, 3, 0.2), "lines"), axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5), axis.title.x = element_text(size = 8), legend.position = c(0.8, -0.2),
        strip.background = element_rect(fill="#FFFFFF"), strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(), legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8), legend.box.margin = margin(0.2, 0.2, 0.2, 0.2),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2), legend.direction = "horizontal")

plot_r0_sens_b = ggplot(mapping=aes(x=time, y=median, ymin=low95, ymax=high95, colour=clrname, fill=clrname, group=clrname))+
  labs(x = "Days since first infectious individual",
       colour="Group", fill="Group", linetype="Group", y = "Cumulative attack rate")+
  coord_cartesian(ylim=c(0, max(plotdat_sub_b$high95)), xlim=c(0, 140), clip = "off")+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(target_R0~iv_shielded_gz_size, labeller = labeller(
    iv_shielded_prop = function(x) paste0(as.numeric(x)*100,"%")))+
  geom_text(data = data.table(target_R0 = 1.5, iv_shielded_gz_size = 2, label = "Shielded residence size", x = 0,
                              y = 1.30), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL, linetype=NULL,
                                             fill=NULL, group=NULL), hjust = 0, vjust = 0, size = 8/(72.27/25.4),
            show.legend = FALSE)+
  geom_text(data = data.table(target_R0 = 1.5, iv_shielded_gz_size = 16, label = "R[0]",
                              x = 184, y = 1), aes(x=x, y=y, label=label, ymin=NULL, ymax=NULL, colour=NULL,
                                                   linetype=NULL, fill=NULL, group=NULL), hjust = 0, vjust = 0,
            angle = -90, size = 8/(72.27/25.4), show.legend = FALSE, parse = TRUE)+
  geom_ribbon(data = plotdat_sub_b[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], alpha=0.2,
              colour=NA)+
  geom_ribbon(data = plotdat_sub_b[scen_type == "value_mit"], alpha=0.2, colour=NA)+
  geom_line(data = plotdat_sub_b[scen_type == "value_unmit" & age == ">=60" & shielded_status == 0], size=1)+
  geom_line(data = plotdat_sub_b[scen_type == "value_mit"], aes(linetype = age), size=1)+
  guides(color = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_none())+
  scale_colour_manual(values=risk_bytime_colours)+
  scale_fill_manual(values=risk_bytime_colours)+
  scale_linetype_manual(values = c("<60" = 2, ">=60" = 1))+
  theme_bw()+
  theme(plot.margin = unit(c(0.75, 0.75, 2, 0.2), "lines"), axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5), axis.title.x = element_text(size = 8), legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF"), strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(), legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8), legend.box.margin = margin(0.2, 0.2, 0.2, 0.2),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2), legend.direction = "horizontal")

library(patchwork)

rel_width=2000
layout = c(
  area(t = 1, b = 2, l = 1, r = 1+rel_width), area(t = 1, b = 2, l = 3+rel_width, r = 3+rel_width*2))

x11(width = plot_double_col, height = plot_double_col)
png("~/workspace/covid_ibm_shielding/analysis/risk_bytime_sens_r0.png", width = plot_double_col,
    height = plot_double_col, units = "in", res = 300)

plot_r0_sens_a+theme(legend.position="none")+plot_r0_sens_b+theme(legend.position=c(-0.199, -0.125))+plot_layout(design = layout)+plot_annotation(tag_levels = "A")&theme(plot.tag = element_text(size=8, face = "bold"), plot.tag.position = c(0,1.025))

dev.off()

#' Without shielding, 60+ have similar infection rates compared to <60
#' After shielding, the unshielded 60+ have lower attach rate compared to the shielded 60+
#' - Likely because of very low contacts after shielding large proportion of their old contacts (essentially shielded as well)
#' As shielded gz increases, risk of infection increases (same levels as unshielded)
#' - Still a protective effect in large R0 for large gz
#' - Risk of unshielded (<60yo) not substantially affected by shielding the 60+
#' - Mitigation options to enter greenzones do not significantly affect risk of breaching
#' - Risk increases as shielding is less effective (reduction of contacts shielded-unshielded is lower)
#'  * But moving people from households into greenzones already gives some protection (no household contacts)
#'  Not huge effect of prior immunity in the population