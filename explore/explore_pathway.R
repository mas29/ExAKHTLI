#####################################################################################
###### ----------- Plot individual sparkline for this pathway ---------------- ######
#####################################################################################


# Function to plot sparkline for pathway and phenotypic marker of interest
# param df -- data in tall format for particular phenotypic marker and pathway
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param phenotypic_marker_name -- name of phenotypic marker (e.g. "Sytox Green")
# param target -- name of target
plot_pathway_sparklines <- function(df, confidence_intervals, phenotypic_marker_name, pathway) {
  if (!(pathway == "")) {
    ggplot() +
      geom_line(data = df, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = "blue")) +
      geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, 
                                                             fill = "red", colour = NULL), alpha = 0.2) +
      scale_fill_manual(name = "Legend", values = 'red', labels = paste(phenotypic_marker_name,'\nNegative Control\n99.9% C.I.',sep="")) +
      xlab("Time Elapsed (hours)") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name," Levels for [Pathway: ",pathway,"] Over Time",sep="")) +
      scale_colour_manual(values=c("blue"), guide = FALSE) + 
      theme(panel.grid = element_blank(),
            strip.text=element_blank(),
            legend.key.height = unit(.85, "cm"),
            panel.background = element_rect(fill = "white"),
            panel.margin = unit(.085, "cm"))
  }
}

# Function to plot clusters for the phenotypic marker of interest, hilighting the pathway of interest
# param data_w_clusters -- df of data with clusters for the phenotypic marker of interest
# param pathway -- pathway of interest
# param penotypic_marker_name -- name of phenotypic marker of interest
# param confidence_intervals_each_marker -- list of confidence intervals data frames for each marker
plot_pathway_clusters <- function(data_w_clusters, pathway, phenotypic_marker_name, confidence_intervals_each_marker) {
  num_time_intervals <- length(unique(data_w_clusters$time_elapsed)) # Number of time intervals
  if (!(pathway == "")) {
    if (nrow(data_w_clusters) == num_time_intervals) { # If there's only one compound in the pathway, we only want one colour
      palette <- c("blue")
    } else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
      palette <- c("black", "blue")
    }
    df_isnt_pathway <- data_w_clusters[data_w_clusters$is.pathway == FALSE,]
    df_is_pathway <- data_w_clusters[data_w_clusters$is.pathway == TRUE,]
    ggplot() +
      geom_ribbon(data = confidence_intervals_each_marker[[phenotypic_marker_name]], 
                  mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, fill = "red", colour = NULL), alpha = 0.2) +
      geom_line(data = df_isnt_pathway, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.pathway), alpha = 0.2) +
      geom_line(data = df_is_pathway, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.pathway)) +
      scale_colour_manual(name = pathway, values=palette) + 
      scale_fill_manual(name = "Legend",
                        values = c('red'),
                        labels = c('Negative Control\n99.9% C.I.')) +
      xlab("Time Elapsed") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name, " Sparklines for Each Cluster\nHighlighting [pathway: ", pathway, "]", sep="")) +
      facet_grid(empty~cluster) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_blank())
  }
}