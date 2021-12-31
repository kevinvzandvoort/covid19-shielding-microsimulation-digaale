analysis_name = "breached"
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
  scens = scenario_list[target_R0 == r0 & prior_immunity == immunity &
                          iv_shielded_prop > 0 & iv_shielded_shielded_contacts == 1 &
                          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1, scen]
  interventions = scens %>% lapply(function(s){
    sprintf("%s/byscen/interventions/byscen/out_scen_%s_model_interventions.qs", out_folder, s) %>% qread})  %>% rbindlist
  
  message("Limit to runs of interest")
  #' Limited to runs in scenarios where shielding was implemented
  scens_of_interest = interventions[intervention == "shielding" & processed == 1, c("scen", "run")]  
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
    #' Calculate proportion breached
    breached_group = model_out[age >= 6 & shielded_gz_id > -1, .(infected = sum(infected_shielded),
                                                                 breached = sum(infected_shielded) > 0),
                               by=c("scen", "run", "shielded_gz_id")] %>%
      .[, .(breached = sum(breached), model_max_gz = (max(shielded_gz_id) + 1)), by=c("scen", "run")] %>%
      .[, prop_breached := breached/model_max_gz] %>% .[, .(mean = mean(prop_breached), median = median(prop_breached),
                                                            low95 = quantile(prop_breached, 0.025),
                                                            high95 = quantile(prop_breached, 0.975),
                                                            low50 = quantile(prop_breached, 0.25),
                                                            high50 = quantile(prop_breached, 0.75)), by = "scen"]
    return(breached_group)}, scens) %>% rbindlist
  return(model_out)
}

#sbatch_command = prepareSlurmBatchRuns(analysis_name, singleFunction, iterations, main_folder, out_folder, packages_func)
result = combineSlurmResults(analysis_name, iterations) #%>%

result %>% merge(scenario_list[, c("scen", "target_R0", "prior_immunity", "iv_shielded_prop", "iv_shielded_gz_size",
                        "iv_shielded_enter", "iv_shielded_unshielded_contacts")], by="scen") %>%
  .[iv_shielded_gz_size == 2 & iv_shielded_unshielded_contacts == 0.2 & target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1 & iv_shielded_enter == 1, sprintf("%s%% (%s-%s)", round(median*100), round(low95*100), round(high95*100))]

plotdat = scenario_list[iv_shielded_prop > 0 & iv_shielded_shielded_contacts == 1 & iv_shielded_symptomatic_contacts == 1 &
                iv_shielded_mitigate == -1] %>%
  plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded", "symptomatic individuals do not enter",
                                     "individuals tested before they enter")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE) %>%
  .[target_R0 == 2.5 & prior_immunity == 0 & iv_shielded_prop == 1]

plotdat[, plottype := "adaptation"]
plotdat[iv_shielded_unshielded_contacts == 0.2 & iv_shielded_shielded_contacts == 1 &
          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_enter == 1, plottype := "baseline"]

x11(width = plot_double_col, height = plot_single_col*0.8)
png("~/workspace/covid_ibm_shielding/analysis/breached.png", width = plot_double_col, height = plot_single_col*0.8,
    units = "in", res = 300)

plotdat %>%
  ggplot(aes(x = x, xend = x, colour = name, fill = name))+
  geom_segment(aes(y=low95, yend=high95), size=8 / (72.27/25.4))+
  #geom_segment(aes(y=median-0.00125, yend=median+0.00125), size=6)+
  #geom_point(aes(y=median), size=4, shape=21, colour="#000000", stroke=1, fill="#FFFFFF")+ #diamond is shape 23
  #geom_label(aes(y=median, label=round(median*100)), size=4, colour="#000000", fill=rgb(1, 1, 1, 0.75))+
  geom_point(aes(y=median, shape = plottype), size=2, fill="#FFFFFF", colour="#000000FF", stroke=1)+
  scale_x_continuous(breaks = plotdat[iv_shielded_gz_size == 2, sort(unique(x))],
    labels = function(xvals) xvals %>%
      sapply(function(xval) plotdat[iv_shielded_gz_size == 2 & x == xval, unique(option_name)]))+
    scale_colour_manual(values = mitigation_options_colours)+
    scale_fill_manual(values = mitigation_options_colours)+
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
  coord_flip(clip = "off", ylim = c(0, 1.02), xlim=c(0, plotdat[, max(x)]))+
  scale_y_continuous(labels = scales::percent_format(suffix=""), breaks=seq(0, 10, 2.5)/10)+
  labs(y = "Proportion of breached residences (%)")+
  geom_text(data = plotdat[iv_shielded_gz_size == min(iv_shielded_gz_size),
                           .(x = median(x), label = factor(name, names(mitigation_options), mitigation_options_plotnames),
                             iv_shielded_gz_size = min(iv_shielded_gz_size)), by="name"],
            aes(x=x, label=label), y=-1.7, vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  geom_text(data = data.table(x=112, y=0, iv_shielded_gz_size = 2, label = "Shielded residence size"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  #annotate(geom = "text",
  #         x = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, x],
  #         y = -0.15,
  #         label = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, name] %>%
  #           factor(names(mitigation_options), mitigation_options_plotnames),
  #         vjust = 0.5, hjust=0.5, angle = 0, size = 11 / (72.27/25.4))#+ #covert pt to mm
  geom_vline(data = data.table(y = plotdat[, x] %>% sort %>% subset(c(F, diff(.) == 20)) - 10),
             aes(xintercept = y), colour = "#AEAEAE", linetype = 2)+
  scale_shape_manual(values = c("baseline" = 23,
                                "adaptation" = 21))

dev.off()

#' Probability that a greenzone is breached is mostly affected by the greenzone size (exponential increase)
#' - Not affected by how many people are shielded (as size of greenzones governed by different parameter)
#' - Mitigation options to enter greenzones do not significantly affect risk of breaching
#' - Risk increases as shielding is less effective (reduction of contacts shielded-unshielded is lower)
#'  * But moving people from households into greenzones already gives some protection (no household contacts)

plotdat = scenario_list[iv_shielded_prop > 0 & iv_shielded_shielded_contacts == 1 & iv_shielded_symptomatic_contacts == 1 &
iv_shielded_mitigate == -1] %>%
  plotLabelledX(mitigation_options[c("reduction contacts shielded unshielded")]) %>%
  merge(scenario_list, by="scen", all.x = TRUE) %>%
  merge(result, by="scen", allow.cartesian = TRUE) %>%
  .[iv_shielded_prop == 1 & prior_immunity == 0]

plotdat[, plottype := "adaptation"]
plotdat[iv_shielded_unshielded_contacts == 0.2 & iv_shielded_shielded_contacts == 1 &
          iv_shielded_symptomatic_contacts == 1 & iv_shielded_mitigate == -1 & iv_shielded_enter == 1, plottype := "baseline"]

x11(width = plot_double_col, height = plot_double_col*0.8)
png("~/workspace/covid_ibm_shielding/analysis/breached_sens_r0.png", width = plot_double_col, height = plot_double_col*0.8,
    units = "in", res = 300)

plotdat %>%
  ggplot(aes(x = x, xend = x, colour = name, fill = name))+
  geom_segment(aes(y=low95, yend=high95), size=8 / (72.27/25.4))+
  #geom_segment(aes(y=median-0.00125, yend=median+0.00125), size=6)+
  #geom_point(aes(y=median), size=4, shape=21, colour="#000000", stroke=1, fill="#FFFFFF")+ #diamond is shape 23
  #geom_label(aes(y=median, label=round(median*100)), size=4, colour="#000000", fill=rgb(1, 1, 1, 0.75))+
  geom_point(aes(y=median, shape = plottype), size=2, fill="#FFFFFF", colour="#000000FF", stroke=1)+
  scale_x_continuous(breaks = plotdat[iv_shielded_gz_size == 2, sort(unique(x))],
                     labels = function(xvals) xvals %>%
                       sapply(function(xval) plotdat[iv_shielded_gz_size == 2 & x == xval, unique(option_name)]))+
  scale_colour_manual(values = mitigation_options_colours)+
  scale_fill_manual(values = mitigation_options_colours)+
  theme_bw()+
  theme(plot.margin = unit(c(1, 1, 0.2, 0.2), "lines"),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7.5),
        axis.title.x = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(fill="#FFFFFF"),
        strip.text = element_text(size = 8),
        panel.grid.minor = element_blank())+
  facet_grid(target_R0 ~ iv_shielded_gz_size, labeller = labeller())+
  coord_flip(clip = "off", ylim = c(0, 1.02), xlim=c(-2.5,52.5))+
  scale_y_continuous(labels = scales::percent_format(suffix=""), breaks=seq(0, 10, 2.5)/10)+
  labs(y = "Proportion of breached residences (%)", x = "Relative reduction in contacts between shielded and unshielded individuals")+
  geom_text(data = plotdat[iv_shielded_gz_size == min(iv_shielded_gz_size),
                           .(x = median(x), label = factor(name, names(mitigation_options), mitigation_options_plotnames),
                             iv_shielded_gz_size = min(iv_shielded_gz_size)), by="name"],
            aes(x=x, label=label), y=-1.7, vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  geom_text(data = data.table(x=72.5, y=0, iv_shielded_gz_size = 2, target_R0 = 1.5, label = "Shielded residence size"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = 0, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9)+
  geom_text(data = data.table(x=50, y=1.32, iv_shielded_gz_size = 16, target_R0 = 1.5, label = "R[0]"),
            aes(x=x, label=label, y=y, fill=NULL), vjust = 0.5, hjust=0, angle = -90, size = 7.5 / (72.27/25.4), colour="#000000", lineheight=0.9, parse=TRUE)+
  #annotate(geom = "text",
  #         x = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, x],
  #         y = -0.15,
  #         label = plotdat %>% .[iv_shielded_gz_size == 2, c("x", "name")] %>% .[, .(x=median(x)), by="name"] %>% .[, name] %>%
  #           factor(names(mitigation_options), mitigation_options_plotnames),
  #         vjust = 0.5, hjust=0.5, angle = 0, size = 11 / (72.27/25.4))#+ #covert pt to mm
  geom_vline(data = data.table(y = plotdat[, x] %>% sort %>% subset(c(F, diff(.) == 20)) - 10),
             aes(xintercept = y), colour = "#AEAEAE", linetype = 2)+
  scale_shape_manual(values = c("baseline" = 23,
                                "adaptation" = 21))
dev.off()

result %>% merge(scenario_list, by="scen") %>% .[iv_shielded_gz_size == 16 & iv_shielded_prop == 1 & target_R0 == 1.5 & prior_immunity == 0 & iv_shielded_unshielded_contacts == 0, sprintf("%s%% (%s - %s)", round(median*100), round(low95*100), round(high95*100))]
