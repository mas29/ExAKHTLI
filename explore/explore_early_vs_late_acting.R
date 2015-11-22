#############################################################################################################
# ---- Get curves for each timepoint as they surpass or dip below negative control confidence interval ---- #
#############################################################################################################

# @param df -- data frame of data with metrics 
# @param phenotypic_marker_name -- phenotypic marker name
# @param confidence_intervals -- confidence intervals for phenotypic marker of interest
# @param above -- separate curves by whether they dip above (TRUE) or below (FALSE) the negative controls
get_early_vs_late_acting <- function(df, confidence_intervals, phenotypic_marker_name, above_or_below) {
  
  if (above_or_below == "Above Negative Control Upperbound") { # Separated by Time Point at which Compound Surpasses Negative Control Upperbound
    ggplot(df) +
      geom_line(aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group = Compound), alpha = 0.6) +
      xlab("Time Elapsed") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name, " Curves\nSeparated by Time when Curve Surpasses Negative Control Upperbound", sep="")) +
      geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper,
                                                             fill = "red", colour = NULL), alpha = 0.6) +
      scale_fill_manual(name = "Legend",
                        values = c('red'),
                        labels = c('Negative Control\n99.9% C.I.')) +
      facet_wrap(~phenotype_value_exceeds_NC_upperbound.timepoint, scales = "fixed") +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            legend.key.size = unit(0.4, "cm"),
            legend.direction="vertical") 
  } else { # Separated by Time Point at which Compound Falls Below Negative Control Lowerbound
    ggplot(df) +
      geom_line(aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group = Compound), alpha = 0.6) +
      xlab("Time Elapsed") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name, " Curves\nSeparated by Time when Curve Falls Below Negative Control Lowerbound", sep="")) +
      geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper,
                                                             fill = "red", colour = NULL), alpha = 0.6) +
      scale_fill_manual(name = "Legend",
                        values = c('red'),
                        labels = c('Negative Control\n99.9% C.I.')) +
      facet_wrap(~phenotype_value_falls_below_NC_lowerbound.timepoint, scales = "fixed") +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            legend.key.size = unit(0.4, "cm"),
            legend.direction="vertical") 
  }
}