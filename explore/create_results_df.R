#' Create Spreadsheet of Results
#'
#' Creates a spreadsheet of results, sends it to the Output folder
#' @export
save_results <- function(phenotypic_marker_names, data_tall_no_NC_each_marker, output_dir) {

  # Function to get ordered information for a particular curve metric for each phenotypic marker
  # param df -- the data frame list (one for each phenotypic marker) of interest
  # param metric_name -- the name of the metric in df
  # param asc_or_desc -- if the data frame should be ordered ascendingly or descendingly
  get_ordered_metric_df <- function(df, metric_name, asc_or_desc, phenotypic_marker_names) {
    curr_metric.each_marker <- list()
    for (i in 1:length(phenotypic_marker_names)) { # For each phenotypic marker
      # Get df for this marker
      marker_df <- df[[phenotypic_marker_names[[i]]]]
      # Get compound and curr_metric value for this marker
      curr_metric.new_marker <- unique(marker_df[,which(colnames(marker_df) %in% c("Compound", metric_name))])
      # Order the compounds by the curr_metric value
      if (asc_or_desc == "asc") {
        curr_metric.new_marker <- curr_metric.new_marker[ order(curr_metric.new_marker[,2]), ]
      }
      if (asc_or_desc == "desc") {
        curr_metric.new_marker <- curr_metric.new_marker[ order(-curr_metric.new_marker[,2]), ]
      }
      colnames(curr_metric.new_marker)[which(colnames(curr_metric.new_marker) == metric_name)] <- paste("value",metric_name,sep=".")
      colnames(curr_metric.new_marker)[which(colnames(curr_metric.new_marker) == "Compound")] <- paste("compound",metric_name,sep=".")
      curr_metric.each_marker[[i]] <- curr_metric.new_marker
    }
    curr_metric.each_marker <- as.list(setNames(curr_metric.each_marker, phenotypic_marker_names))
    curr_metric.all_markers <- do.call("cbind", curr_metric.each_marker) # cbind for all markers

    return(curr_metric.all_markers)
  }

  # Get data for each metric

  time_x_distance.lower.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_x_distance.lower", "desc", phenotypic_marker_names)
  time_x_distance.upper.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_x_distance.upper", "desc", phenotypic_marker_names)
  phenotype_value_exceeds_NC_upperbound.timepoint.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "phenotype_value_exceeds_NC_upperbound.timepoint", "desc", phenotypic_marker_names)
  phenotype_value_falls_below_NC_lowerbound.timepoint.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "phenotype_value_falls_below_NC_lowerbound.timepoint", "desc", phenotypic_marker_names)
  time_to_max.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_to_max", "asc", phenotypic_marker_names)
  time_to_min.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_to_min", "asc", phenotypic_marker_names)
  max.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "max", "desc", phenotypic_marker_names)
  min.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "min", "asc", phenotypic_marker_names)
  AUC.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "AUC_trapezoidal_integration", "desc", phenotypic_marker_names)
  Delta_max_min.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "delta_min_max", "desc", phenotypic_marker_names)
  most_positive_slope.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "most_positive_slope", "desc", phenotypic_marker_names)
  most_negative_slope.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "most_negative_slope", "desc", phenotypic_marker_names)
  time_to_most_positive_slope.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_to_most_positive_slope", "asc", phenotypic_marker_names)
  time_to_most_negative_slope.all_markers <- get_ordered_metric_df(data_tall_no_NC_each_marker, "time_to_most_negative_slope", "asc", phenotypic_marker_names)

  # Create list of markers
  results <- list(time_x_distance.upper.all_markers=time_x_distance.upper.all_markers, 
                  time_x_distance.lower.all_markers=time_x_distance.lower.all_markers,
                  phenotype_value_exceeds_NC_upperbound.timepoint.all_markers=phenotype_value_exceeds_NC_upperbound.timepoint.all_markers,
                  phenotype_value_falls_below_NC_lowerbound.timepoint.all_markers=phenotype_value_falls_below_NC_lowerbound.timepoint.all_markers,
                  max.all_markers=max.all_markers,
                  min.all_markers=min.all_markers,
                  time_to_max.all_markers=time_to_max.all_markers,
                  time_to_min.all_markers=time_to_min.all_markers,
                  AUC.all_markers=AUC.all_markers,
                  Delta_max_min.all_markers=Delta_max_min.all_markers,
                  most_positive_slope.all_markers=most_positive_slope.all_markers,
                  most_negative_slope.all_markers=most_negative_slope.all_markers,
                  time_to_most_positive_slope.all_markers=time_to_most_positive_slope.all_markers,
                  time_to_most_negative_slope.all_markers=time_to_most_negative_slope.all_markers)
  
  # Save R object list of results
  save(results, file=paste(output_dir,"Results.RData",sep=""))

}

