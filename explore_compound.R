#####################################################################################
###### ----------- Plot individual sparkline for this compound --------------- ######
#####################################################################################


# Function to plot sparkline for compound and phenotypic marker of interest
# param df -- data in tall format for particular phenotypic marker and compound
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param compound -- compound of interest
plot_sparkline <- function(df, confidence_intervals, phenotypic_marker_name, compound) {
  if (!(compound == "")) {
    ggplot() +
      geom_line(data = df, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = "blue")) +
      geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, 
                                                             fill = "red", colour = NULL), alpha = 0.2) +
      scale_fill_manual(name = "Legend", values = 'red', labels = paste(phenotypic_marker_name,'\nNegative Control\n99.9% C.I.',sep="")) +
      xlab("Time Elapsed (hours)") +
      ylab(phenotypic_marker_name) +
      ggtitle(paste(phenotypic_marker_name," Levels for [Compound: ",compound,"] Over Time",sep="")) +
      scale_colour_manual(values=c("blue"), guide = FALSE) + 
      theme(panel.grid = element_blank(),
            strip.text=element_blank(),
            legend.key.height = unit(.85, "cm"),
            panel.background = element_rect(fill = "white"),
            panel.margin = unit(.085, "cm"))
  }
}

#####################################################################################
###### ----------- Plot sparklines for this compound's target ---------------- ######
#####################################################################################

# Function to plot the sparklines for the target of the specified compound
# param df -- data in tall format for particular phenotypic marker and compound
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param compound -- compound of interest
# param target -- target of compount of interest
plot_target <- function(df, confidence_intervals, phenotypic_marker_name, compound, target) {
  num_time_intervals <- length(unique(df$time_elapsed)) # Number of time intervals
  if (nrow(df) == num_time_intervals) { # If there's only one compound in the target, we only want one colour
    palette <- c("blue")
  } else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
    palette <- c("black", "blue")
  }
  df_isnt_compound <- df[df$is.compound == FALSE,]
  df_is_compound <- df[df$is.compound == TRUE,]
  ggplot() +
    geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper,
                                                           fill = "red", colour = NULL), alpha = 0.2) +
    geom_line(data = df_isnt_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound), alpha = 0.3) +
    geom_line(data = df_is_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound)) +
    scale_colour_manual(name = compound, values=palette) + 
    scale_fill_manual(name = "", values = c('red'), labels = c('Negative Control\n99.9% C.I.')) +
    xlab("Time Elapsed (hours)") +
    ylab(phenotypic_marker_name) +
    ggtitle(paste(phenotypic_marker_name," Levels for [Target: ",target,"] Over Time",sep="")) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.key.size = unit(0.4, "cm"),
          legend.direction="vertical")
}

#####################################################################################
###### ----------- Plot sparklines for this compound's pathway --------------- ######
#####################################################################################

# Function to plot the sparklines for the pathway of the specified compound
# param df -- data in tall format for particular phenotypic marker and compound
# param confidence_intervals -- confidence intervals for the negative control of this phenotypic marker
# param compound -- compound of interest
# param pathway -- pathway of compount of interest
plot_pathway <- function(df, confidence_intervals, phenotypic_marker_name, compound, pathway) {
  num_time_intervals <- length(unique(df$time_elapsed)) # Number of time intervals
  if (nrow(df) == num_time_intervals) { # If there's only one compound in the target, we only want one colour
    palette <- c("blue")
  }
  else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
    palette <- c("black", "blue")
  }
  df_isnt_compound <- df[df$is.compound == FALSE,]
  df_is_compound <- df[df$is.compound == TRUE,]
  ggplot(df) +
    geom_ribbon(data = confidence_intervals, mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper,
                                                           fill = "red", colour = NULL), alpha = 0.2) +
    geom_line(data = df_isnt_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound), alpha = 0.3) +
    geom_line(data = df_is_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound)) +
    scale_colour_manual(name = compound, values=palette) + 
    scale_fill_manual(name = "", values = c('red'), labels = c('Negative Control\n99.9% C.I.')) +
    xlab("Time Elapsed (hours)") +
    ylab(phenotypic_marker_name) +
    ggtitle(paste(phenotypic_marker_name," Levels for [Pathway: ",pathway,"] Over Time",sep="")) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.key.size = unit(0.4, "cm"),
          legend.direction="vertical")
}

#####################################################################################
###### --- Plot the clusters, and colour this compound's sparkline ----------- ######
#####################################################################################

# Function to get the clustering of your data
# param prelim_df -- the df as a result of the preliminary_processing() function in the GetData script
# param df -- the df in "tall" format, for only your phenotypic marker of interest
# param phenotypic_Marker_name -- name of your phenotypic marker, as is shown in the input dataset (e.g. "SG", "Con")
# param num_clusters -- how many clusters you want
get_data_w_clusters <- function(prelim_df, df, phenotypic_Marker_name, num_clusters) {
  # Parameters
  start_index <- which(colnames(prelim_df) == "0")
  end_index <- which(colnames(prelim_df) == "46")
  num_time_intervals <- length(unique(df$time_elapsed)) # Number of time intervals
  
  # Get Sytox Green raw time series data only, for all compounds (including negative controls)
  data_for_heatmap <- prelim_df[prelim_df$phenotypic_Marker == phenotypic_Marker_name, start_index:end_index]
  rownames(data_for_heatmap) <- prelim_df[prelim_df$phenotypic_Marker == phenotypic_Marker_name,]$Compound
  data_matrix <- as.matrix(data_for_heatmap)
  
  # Get distances (Euclidean) and clusters
  distMatrix <- dist(data_matrix, method="euclidean")
  hr <- hclust(distMatrix, method="average")
  
  # Cut the tree and create color vector for clusters.
  mycl <- cutree(hr, k = num_clusters) # Clusters assigned to each compound.
  mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9)
  mycolhc <- mycolhc[as.vector(mycl)] 
  
  # Plot all clusters
  compound_clusters <- as.data.frame(mycl)
  colnames(compound_clusters) <- "cluster"
  compound_clusters$cluster <- as.factor(compound_clusters$cluster)
  df_w_clusters <- merge(df, compound_clusters, by.x="Compound", by.y="row.names")
}

# Function to plot clusters for the phenotypic marker of interest, and the compound of interest
# param data_w_clusters -- df of data with clusters for the phenotypic marker of interest
# param compound -- compound of interest
# param penotypic_marker_name -- name of phenotypic marker of interest
# param confidence_intervals_each_marker -- list of confidence intervals data frames for each marker
plot_clusters <- function(data_w_clusters, compound, phenotypic_marker_name, confidence_intervals_each_marker) {
  num_time_intervals <- length(unique(data_w_clusters$time_elapsed)) # Number of time intervals
  if (nrow(data_w_clusters) == num_time_intervals) { # If there's only one compound in the target, we only want one colour
    palette <- c("blue")
  } else { # Otherwise, we want two colours, one for the compound of focus, the other for the other compounds
    palette <- c("black", "blue")
  }
  df_isnt_compound <- data_w_clusters[data_w_clusters$is.compound == FALSE,]
  df_is_compound <- data_w_clusters[data_w_clusters$is.compound == TRUE,]
  ggplot() +
    geom_ribbon(data = confidence_intervals_each_marker[[phenotypic_marker_name]], 
                mapping = aes(x = time_elapsed, ymin = phenotype_value.NC.lower, ymax = phenotype_value.NC.upper, fill = "red", colour = NULL), alpha = 0.2) +
    geom_line(data = df_isnt_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound), alpha = 0.2) +
    geom_line(data = df_is_compound, aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound, colour = is.compound)) +
    scale_colour_manual(name = compound, values=palette) + 
    scale_fill_manual(name = "Legend",
                      values = c('red'),
                      labels = c('Negative Control\n99.9% C.I.')) +
    xlab("Time Elapsed") +
    ylab(phenotypic_marker_name) +
    ggtitle(paste(phenotypic_marker_name, "Sparklines for Each Cluster", sep=" ")) +
    facet_grid(empty~cluster) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.text.x = element_blank())
}
