analysis_name = "cumulative_risk"
packages_func = c("qs", "magrittr", "data.table")
packages_other = c("ggplot2")

pacman::p_load(char = c(packages_func, packages_other))
source("./helper_functions.R")

main_folder = c("/nfs/general/covid_ibm_shielding", "~/workspace/covid_ibm_shielding")[2]
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
  scens = scenario_list[target_R0 == r0 & prior_immunity == immunity, scen]
  interventions = scens %>% lapply(function(s){
    sprintf("%s/byscen/interventions/byscen/out_scen_%s_model_interventions.qs", out_folder, s) %>% qread})  %>% rbindlist
  
  message("Limit to runs of interest")
  
  #' Limited to runs in scenarios where shielding was implemented
  scens_of_interest = interventions[is.na(intervention) | (intervention == "shielding" & processed == 1), c("scen", "run", "intervention")]  
  scen_unmitigated = scens_of_interest[is.na(intervention)]
  if(nrow(scen_unmitigated) == 0){
    message("No unmitigated scenarios applicable")
    return(NULL)
  }
  
  scens_of_interest = scens_of_interest[!is.na(intervention)]
  if(nrow(scens_of_interest) == 0){
    message("No non-unmitigated scenarios applicable")
    return(NULL)
  }
  
  age_groups = data.table(age=0:8)
  
  scens_unmit = scen_unmitigated[, unique(scen)]
  unmitigated_model_out = seq_len(length(scens_unmit)) %>% lapply(function(s, scens_unmit){
    message(sprintf("Unmit %s/%s", s, length(scens_unmit)))
    model_out = sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, scens_unmit[s]) %>% qread %>%
      .[run %in% scen_unmitigated[scen == scens_unmit[s], run]]
    
    model_out = model_out[, .(cases = sum(age >= 6 & infected_at >= 0)), by=c("scen", "run")]
  }, scens_unmit) %>% rbindlist
  unmitigated_model_out = unmitigated_model_out[, unmitigated_cases := cases] %>% .[, -c("cases", "scen")]
  #' Exclude any runs with 0 unmitigated cases
  unmitigated_model_out = unmitigated_model_out[unmitigated_cases > 0]
  
  message("Get data for each scenario")
  scens = scens_of_interest[, unique(scen)]
  model_out = seq_len(length(scens)) %>% lapply(function(s, scens){
    message(sprintf("%s/%s", s, length(scens)))
    model_out = sprintf("%s/byscen/model_out/byscen/out_scen_%s_model_out.qs", out_folder, scens[s]) %>% qread %>%
      .[run %in% scens_of_interest[scen == scens[s], run]]
    
    model_out = model_out[, .(cases = sum(age >= 6 & infected_at >= 0)), by=c("scen", "run")]
    model_out = model_out %>% merge(unmitigated_model_out, by=c("run"), all = TRUE)
    model_out[is.na(cases), cases := 0]

    model_out = model_out[, .(diff_abs = as.double(unmitigated_cases) - cases,
                  diff_rel = cases/unmitigated_cases), by=c("scen", "run")] %>%
      melt(measure.vars=c("diff_abs", "diff_rel")) %>%
      .[, .(mean = mean(value),
            median = median(value),
            low95 = quantile(value, 0.025),
            high95 = quantile(value, 0.975),
            low50 = quantile(value, 0.25),
            high50 = quantile(value, 0.75)), by = c("scen", "variable")]
    return(model_out)}, scens) %>% rbindlist
  return(model_out)
}

#sbatch_command = prepareSlurmBatchRuns(analysis_name, singleFunction, iterations, main_folder, out_folder, packages_func, populations, cpus_per_task = 4)
result = combineSlurmResults(analysis_name, iterations)

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == 1 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & variable == "diff_rel", sprintf("%s%% (%s-%s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 16 & iv_shielded_unshielded_contacts == 0 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & variable == "diff_rel", sprintf("%s%% (%s-%s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 16 & iv_shielded_unshielded_contacts == 0.4 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == 1 & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 8 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 8 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0.4) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 16 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0.6) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 4 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 4 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == 0, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 8 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 8 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == 0, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

#
result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 1) & target_R0 == 5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 1) & target_R0 == 5 & prior_immunity == 0.25 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 1) & target_R0 == 5 & prior_immunity == 0.5 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

#
result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 0.4) & target_R0 == 5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 0.4) & target_R0 == 5 & prior_immunity == 0.25 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == (1 - 0.4) & target_R0 == 5 & prior_immunity == 0.5 & iv_shielded_prop == 1 & iv_shielded_shielded_contacts == (1 - 0) & iv_shielded_enter == 1 & variable == "diff_rel" & iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, sprintf("%s%% (%s - %s)", round((1-median)*100), round((1-high95)*100), round((1-low95)*100))]

plotdat = scenario_list[iv_shielded_prop > 0] %>%
  plotLabelledX(mitigation_options[-1]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE, all.x = TRUE) %>%
  .[target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1]

plotdat[, plottype := "adaptation"]
plotdat[iv_shielded_unshielded_contacts == 0.2 & iv_shielded_shielded_contacts == 1 &
          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_enter == 1, plottype := "baseline"]

x11(width = plot_double_col, height = plot_double_col*0.8)
png(filename = "~/workspace/covid_ibm_shielding/analysis/cumulative_risk.png",
    width = plot_double_col, height = plot_double_col*0.8, units = "in", res = 300)

plotdat %>% .[variable == "diff_rel"] %>%
  ggplot(aes(x = x, xend = x, colour = name, fill = name))+
  geom_hline(data = NULL, yintercept = 0, colour = "#777777")+
  geom_segment(aes(y=1-high95, yend=1-low95), size=8 / (72.27/25.4))+
  #geom_segment(aes(y=median-0.00125, yend=median+0.00125), size=6)+
  geom_point(aes(y=1-median, fill=plottype, shape=plottype), size=2, colour="#000000FF", stroke=1)+ #diamond is shape 23
  #geom_label(aes(y=1-median, label=round((1-median)*100)), size=7 / (72.27/25.4), colour="#000000", fill=rgb(1, 1, 1, 0.75))+
  scale_x_continuous(breaks = plotdat[iv_shielded_gz_size == 2, sort(unique(x))],
                     labels = function(xvals) xvals %>%
                       sapply(function(xval) plotdat[iv_shielded_gz_size == 2 & x == xval, unique(option_name)]))+
  scale_colour_manual(values = mitigation_options_colours)+
  #scale_fill_manual(values = mitigation_options_colours)+
  theme_bw()+
  theme(plot.margin = unit(c(1, 0.1, 0.2, 7.5), "lines"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 7.5),
        axis.title.x = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank())+
  facet_grid(. ~ iv_shielded_gz_size, labeller = labeller())+
  coord_flip(clip = "off", ylim = c(1 - plotdat[variable == "diff_rel", max(high95)], NA),
             xlim = c(5, 285))+
  scale_y_continuous(labels = scales::percent_format(suffix=""), breaks=seq(-10, 10, 2.5)/10)+
  labs(y = "Relative reduction in cases compared to no-shielding in people 60+y (%)")+
  geom_text(data = plotdat[iv_shielded_gz_size == min(iv_shielded_gz_size),
                           .(x = median(x), label = factor(name, names(mitigation_options), mitigation_options_plotnames),
                             iv_shielded_gz_size = min(iv_shielded_gz_size)), by="name"],
            aes(x=x, label=label), y=-1.8, vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  geom_text(data = data.table(x=321, y=-0.134, iv_shielded_gz_size = 2, label = "Shielded residence size"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  #annotate(geom = "text",
  #         x = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, x],
  #         y = -0.15,
  #         label = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, name] %>%
  #           factor(names(mitigation_options), mitigation_options_plotnames),
  #         vjust = 0.5, hjust=0.5, angle = 0, size = 11 / (72.27/25.4))#+ #covert pt to mm
  geom_vline(data = data.table(y = plotdat[, x] %>% sort %>% subset(c(F, diff(.) == 20)) - 10),
             aes(xintercept = y), colour = "#AEAEAE", linetype = 2)+
  scale_fill_manual(values = c("baseline" = rgb(1,1,1,1),
                               #"adaptation" = rgb(1,1,1,0.5)))+
                               "adaptation" = rgb(1,1,1,1)))+
  scale_shape_manual(values = c("baseline" = 23,
                                "adaptation" = 21))

dev.off()

plotdat = scenario_list[iv_shielded_prop > 0] %>%
  #plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded", "change contacts shielded")]) %>%
  plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE, all.x = TRUE) %>%
  .[iv_shielded_gz_size %in% c(2, 16) & prior_immunity == 0 & iv_shielded_prop == 1]

plotdat[, plottype := "adaptation"]
plotdat[iv_shielded_unshielded_contacts == 0.2 & iv_shielded_shielded_contacts == 1 &
          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_enter == 1, plottype := "baseline"]

x11(width = plot_double_col, height = plot_single_col*0.9)
png(filename = "~/workspace/covid_ibm_shielding/analysis/cumulative_risk_sens_r0.png",
    width = plot_double_col, height = plot_single_col*0.9, units = "in", res = 300)

plotdat %>% .[variable == "diff_rel"] %>%
  ggplot(aes(x = x, xend = x, colour = name, fill = name))+
  geom_hline(data = NULL, yintercept = 0, colour = "#777777")+
  geom_segment(aes(y=1-high95, yend=1-low95), size=4)+
  #geom_segment(aes(y=median-0.00125, yend=median+0.00125), size=6)+
  #geom_point(aes(y=median), size=4, shape=21, colour="#000000", stroke=1, fill="#FFFFFF")+ #diamond is shape 23
  geom_point(aes(y=1-median, fill=plottype), size=2, shape=21, colour="#000000FF", stroke=1)+ #diamond is shape 23
  #geom_label(aes(y=1-median, label=round((1-median)*100)), size=4, colour="#000000", fill=rgb(1, 1, 1, 0.75))+
  scale_x_continuous(breaks = plotdat[iv_shielded_gz_size == 2, sort(unique(x))],
                     labels = function(xvals) xvals %>%
                       sapply(function(xval) plotdat[iv_shielded_gz_size == 2 & x == xval, unique(option_name)]))+
  scale_colour_manual(values = mitigation_options_colours)+
  #scale_fill_manual(values = mitigation_options_colours)+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 0.2, 0.2), "lines"),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5),
        axis.title.x = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank())+
  facet_grid(iv_shielded_gz_size ~ target_R0, labeller = labeller(), scales="free_y")+
  coord_flip(clip = "off", xlim = c(-5, 55), ylim=c(plotdat %>% .[variable == "diff_rel"] %>% .[, 1 - max(high95)],
                                                    plotdat %>% .[variable == "diff_rel"] %>% .[, 1 - min(low95)]))+
  scale_y_continuous(labels = scales::percent_format(suffix=""), breaks=seq(-10, 10, 2.5)/10)+
  labs(y = "Relative reduction in cases compared to no-shielding in people 60+y (%)",
       x = "Relative reduction in contacts between\nshielded and unshielded individuals")+
  geom_text(data = data.table(x=78, y=-0.5, target_R0 = 1.5, iv_shielded_gz_size = 2, label = "Shielded residence size"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4),
            colour="#000000", lineheight=0.9)+
  geom_text(data = data.table(x=55, y=1.25, target_R0 = 5, iv_shielded_gz_size = 2, label = "R[0]"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = -90, size = 7.5 / (72.27/25.4),
            colour="#000000", lineheight=0.9, parse = TRUE)+
  #geom_text(data = plotdat[target_R0 == min(target_R0),
  #                         .(x = median(x),
  #                           label = factor(name, names(mitigation_options), mitigation_options_plotnames),
  #                           target_R0 = min(target_R0)), by=c("name", "iv_shielded_gz_size")],
  #          aes(x=x, label=label), y=-1.775, vjust = 0.5, hjust=0, angle = 0, size = 11 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  #annotate(geom = "text",
  #         x = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, x],
  #         y = -0.15,
  #         label = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, name] %>%
  #           factor(names(mitigation_options), mitigation_options_plotnames),
  #         vjust = 0.5, hjust=0.5, angle = 0, size = 11 / (72.27/25.4))#+ #covert pt to mm
  geom_vline(data = data.table(y = plotdat[, x] %>% sort %>% subset(c(F, diff(.) == 20)) - 10),
             aes(xintercept = y), colour = "#AEAEAE", linetype = 2)+
  scale_fill_manual(values = c("baseline" = rgb(1,1,1,1),
                               "adaptation" = rgb(1,1,1,1)))

dev.off()

plotdat = scenario_list[iv_shielded_prop > 0] %>%
  #plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded", "change contacts shielded")]) %>%
  plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE, all.x = TRUE) %>%
  .[iv_shielded_gz_size %in% c(2, 16) & iv_shielded_prop == 1 & !(target_R0 <= 2.5 & prior_immunity == 0.5)]

plotdat[, plottype := "adaptation"]
plotdat[iv_shielded_unshielded_contacts == 0.2 & iv_shielded_shielded_contacts == 1 &
          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_enter == 1, plottype := "baseline"]

x11(width = plot_double_col, height = plot_double_col*0.8)
png(filename = "~/workspace/covid_ibm_shielding/analysis/cumulative_risk_sens_r0_immune.png",
    width = plot_double_col, height = plot_double_col*0.8, units = "in", res = 300)

plotdat %>% .[variable == "diff_rel"] %>%
  ggplot(aes(x = x, xend = x, colour = factor(iv_shielded_gz_size, c("2", "16"), c("2", "16")), fill = name))+
  geom_hline(data = NULL, yintercept = 0, colour = "#777777")+
  geom_segment(aes(y=1-high95, yend=1-low95), size=4, alpha=0.75)+
  #geom_segment(aes(y=median-0.00125, yend=median+0.00125), size=6)+
  #geom_point(aes(y=median), size=4, shape=21, colour="#000000", stroke=1, fill="#FFFFFF")+ #diamond is shape 23
  geom_point(aes(y=1-median), size=2, shape=21, colour="#000000FF", fill="#FFFFFF", stroke=1)+ #diamond is shape 23
  #geom_label(aes(y=1-median, label=round((1-median)*100)), size=4, colour="#000000", fill=rgb(1, 1, 1, 0.75))+
  scale_x_continuous(breaks = plotdat[iv_shielded_gz_size == 2, sort(unique(x))],
                     labels = function(xvals) xvals %>%
                       sapply(function(xval) plotdat[iv_shielded_gz_size == 2 & x == xval, unique(option_name)]))+
  #scale_colour_manual(values = mitigation_options_colours)+
  #scale_fill_manual(values = mitigation_options_colours)+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 1.2, 0.2), "lines"),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5),
        axis.title.x = element_text(size = 8),
        legend.position = c(0.5, -0.095),
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.box.margin = margin(0.2, 0.2, 0.2, 0.2),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        #legend.background = element_rect(fill="red"),
        #legend.box.background = element_rect(fill="yellow"),
        legend.direction = "horizontal")+
  facet_grid(prior_immunity ~ target_R0, labeller = labeller(
    prior_immunity = function(x) paste0(as.numeric(x)*100,"%")), scales="free_y")+
  coord_flip(clip = "off", ylim = c(1 - plotdat %>% .[variable == "diff_rel", max(high95)],
                                    1 - plotdat %>% .[variable == "diff_rel", min(low95)]),
             xlim = c(-2.5, 52.5))+
  scale_y_continuous(labels = scales::percent_format(suffix=""), breaks=seq(-10, 10, 2.5)/10)+
  labs(y = "Relative reduction in cases compared to no-shielding in people 60+y (%)",
       x = "Relative reduction in contacts between shielded and unshielded individuals",
       colour = "Shielded residence size")+
  geom_text(data = data.table(x=69, y=-0.80, target_R0 = 1.5, prior_immunity = 0, label = "R[0]"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4),
            colour="#000000", lineheight=0.9, parse = TRUE)+
  geom_text(data = data.table(x=55, y=1.35, target_R0 = 5, prior_immunity = 0, label = "Prior immunity"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = -90, size = 7.5 / (72.27/25.4),
            colour="#000000", lineheight=0.9)+
  #geom_text(data = plotdat[target_R0 == min(target_R0),
  #                         .(x = median(x),
  #                           label = factor(name, names(mitigation_options), mitigation_options_plotnames),
  #                           target_R0 = min(target_R0)), by=c("name", "iv_shielded_gz_size")],
  #          aes(x=x, label=label), y=-1.775, vjust = 0.5, hjust=0, angle = 0, size = 11 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  #annotate(geom = "text",
  #         x = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, x],
  #         y = -0.15,
  #         label = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, name] %>%
  #           factor(names(mitigation_options), mitigation_options_plotnames),
  #         vjust = 0.5, hjust=0.5, angle = 0, size = 11 / (72.27/25.4))#+ #covert pt to mm
  geom_vline(data = data.table(y = plotdat[, x] %>% sort %>% subset(c(F, diff(.) == 20)) - 10),
             aes(xintercept = y), colour = "#AEAEAE", linetype = 2)+
  scale_colour_manual(values = c("2" = "#00BF6F",
                                 "16" = "#0D5257"))+
  geom_rect(data = data.table(target_R0 = c(1.5, 2.5),
                              prior_immunity = 0.5,
                              x=NA_real_, y=NA_real_, iv_shielded_gz_size=NA_character_, name=NA),
            xmin=-7.5, xmax=57.5, ymin=-0.85, ymax=0.93, fill="#EEEEEE", colour="#FFFFFF00")

dev.off()

#' Substantial proportion already immune, shielded people immune as well :: shielding with 0% reduction in contacts is larger impact
#' When R0 >= 2.5, if protection is more leaky (less effective), impact is higher when more people are already immune (additional cases prevented in susceptibles) More immune, less contacts with infected, less cases, higher reduction (makes sense). If effective protection, less contact with infected regardless. More opportunity for prevention when shielding is very effective compared to when people are already immune (makes sense).