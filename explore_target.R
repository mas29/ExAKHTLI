#####################################################################################
###### ----------- Plot individual sparkline for this target ----------------- ######
#####################################################################################


# Function to plot sparkline for target and phenotypic marker of interest
# param df -- data in tall format for particular phenotypic marker and target
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param phenotypic_marker_name -- name of phenotypic marker (e.g. "Sytox Green")
# param target -- name of target
plot_target_sparklines <- function(df, confidence_intervals, phenotypic_marker_name, target) {
  if (!(target == "")) {
    ggplot() +
      geom_line(data = df, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = "blue")) +
      geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, 
                                                             fill = "red", colour = NULL), alpha = 0.2) +
      scale_fill_manual(name = "Legend", values = 'red', labels = paste(phenotypic_marker_name,'\nNegative Control\n99.9% C.I.',sep="")) +
      xlab("Time Elapsed (hours)") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name," Levels for [Target: ",target,"] Over Time",sep="")) +
      scale_colour_manual(values=c("blue"), guide = FALSE) + 
      theme(panel.grid = element_blank(),
            strip.text=element_blank(),
            legend.key.height = unit(.85, "cm"),
            panel.background = element_rect(fill = "white"),
            panel.margin = unit(.085, "cm"))
  }
}

##############################################################################
###### ----------- Plot target in context of pathways ----------------- ######
##############################################################################

# Function to plot sparklines for target and phenotypic marker of interest in the context of each pathway
# param df -- data in tall format for particular phenotypic marker and target
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param phenotypic_marker_name -- name of phenotypic marker (e.g. "Sytox Green")
# param target -- name of target
# param confidence_intervals_each_marker -- list of confidence intervals data frames for each marker
plot_pathway_hilight_target <- function(df, confidence_intervals, phenotypic_marker_name, target, confidence_intervals_each_marker) {
  num_time_intervals <- length(unique(df$time_elapsed)) # Number of time intervals
  if (!(target == "")) {
    if (nrow(df) == num_time_intervals) { # If there's only one compound in the target, we only want one colour
      palette <- c("blue")
    } else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
      palette <- c("black", "blue")
    }
    df_isnt_target <- df[df$is.target == FALSE,]
    df_is_target <- df[df$is.target == TRUE,]
    ggplot() +
      geom_ribbon(data = confidence_intervals_each_marker[[phenotypic_marker_name]], 
                  mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, fill = "red", colour = NULL), alpha = 0.2) +
      geom_line(data = df_isnt_target, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.target), alpha = 0.2) +
      geom_line(data = df_is_target, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.target)) +
      scale_colour_manual(name = target, values=palette) + 
      scale_fill_manual(name = "Legend",
                        values = c('red'),
                        labels = c('Negative Control\n99.9% C.I.')) +
      xlab("Time Elapsed") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name, " Sparklines for Each Pathway\nHighlighting [Target: ", target, "]", sep="")) +
      facet_wrap(~Pathway) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_blank())
  }
}

##############################################################################
###### ----------- Plot target in context of clusters ----------------- ######
##############################################################################

# Function to plot clusters for the phenotypic marker of interest, hilighting the target of interest
# param data_w_clusters -- df of data with clusters for the phenotypic marker of interest
# param target -- target of interest
# param penotypic_marker_name -- name of phenotypic marker of interest
# param confidence_intervals_each_marker -- list of confidence intervals data frames for each marker
plot_target_clusters <- function(data_w_clusters, target, phenotypic_marker_name, confidence_intervals_each_marker) {
  num_time_intervals <- length(unique(data_w_clusters$time_elapsed)) # Number of time intervals
  if (!(target == "")) {
    if (nrow(data_w_clusters) == num_time_intervals) { # If there's only one compound in the target, we only want one colour
      palette <- c("blue")
    } else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
      palette <- c("black", "blue")
    }
    df_isnt_target <- data_w_clusters[data_w_clusters$is.target == FALSE,]
    df_is_target <- data_w_clusters[data_w_clusters$is.target == TRUE,]
    ggplot() +
      geom_ribbon(data = confidence_intervals_each_marker[[phenotypic_marker_name]], 
                  mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, fill = "red", colour = NULL), alpha = 0.2) +
      geom_line(data = df_isnt_target, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.target), alpha = 0.2) +
      geom_line(data = df_is_target, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.target)) +
      scale_colour_manual(name = target, values=palette) + 
      scale_fill_manual(name = "Legend",
                        values = c('red'),
                        labels = c('Negative Control\n99.9% C.I.')) +
      xlab("Time Elapsed") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name, " Sparklines for Each Cluster\nHighlighting [Target: ", target, "]", sep="")) +
      facet_grid(empty~cluster) +
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_blank())
  }
}