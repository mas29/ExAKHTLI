# Get list of compounds
compounds <- as.character(sort(unique(data_wide$Compound)))
compound_list <- as.list(setNames(compounds, compounds))
overview_compounds <- c("All Compounds", compounds)
overview_compounds_list <- as.list(setNames(overview_compounds, overview_compounds))

# Get list of targets
targets <- as.character(sort(unique(data_wide$target)))
target_list <- as.list(setNames(targets, targets))
overview_targets <- c("All Targets", targets)
overview_targets_list <- as.list(setNames(overview_targets, overview_targets))

# Get list of pathways
pathways <- as.character(sort(unique(data_wide$Pathway)))
pathway_list <- as.list(setNames(pathways, pathways))
overview_pathways <- c("All Pathways", pathways)
overview_pathways_list <- as.list(setNames(overview_pathways, overview_pathways))

# Get list of metrics
metrics <- c("min", "mean", "max", "delta_min_max", "time_to_min", "time_to_max", "time_x_distance.upper", "time_x_distance.lower", "most_positive_slope", "most_negative_slope", "time_to_most_positive_slope", "time_to_most_negative_slope", "AUC_trapezoidal_integration", "phenotype_value_exceeds_NC_upperbound.timepoint", "phenotype_value_falls_below_NC_lowerbound.timepoint")
metric_names <- c("Minimum", "Mean", "Maximum", "Delta (max-min)", "Time To Minimum Value", "Time To Maximum Value", "Time*Distance Above N.C. C.I.", "Time*Distance Below N.C. C.I.", "Most Positive Slope", "Most Negative Slope", "Time To Most Positive Slope", "Time To Most Negative Slope", "AUC (trapezoidal integration)", "Time when Curve Surpasses N.C. C.I. Upperbound", "Time when Curve Falls below N.C. C.I. Lowerbound")
metrics <- setNames(metrics, metric_names)
metric_list <- as.list(setNames(metric_names, metric_names))

# Above or below
above_or_below <- c("Above Negative Control Upperbound", "Below Negative Control Lowerbound")
above_or_below_list <- as.list(setNames(above_or_below, above_or_below))

# Lists of metrics minima and maxima
max_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$max))) # The maximum maxima
max_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$max))) # The minimum maxima 
min_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$min))) # The maximum minima
min_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$min))) # The minimum minima
delta_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$delta_min_max))) # The maximum delta
delta_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$delta_min_max))) # The minimum delta
auc_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$AUC_trapezoidal_integration))) # The maximum auc
auc_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$AUC_trapezoidal_integration))) # The minimum auc
timexdistance.upper_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$time_x_distance.upper))) # The maximum time*distance upper
timexdistance.upper_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$time_x_distance.upper))) # The minimum time*distance upper
timexdistance.lower_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$time_x_distance.lower))) # The maximum time*distance lower
timexdistance.lower_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$time_x_distance.lower))) # The minimum time*distance lower
most.positive.slope_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$most_positive_slope))) # The maximum most positive slope
most.positive.slope_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$most_positive_slope))) # The minimum most positive slope
most.negative.slope_max_list <- lapply(data_tall_each_marker, function(x) max(ceiling(x$most_negative_slope))) # The maximum most negative slope
most.negative.slope_min_list <- lapply(data_tall_each_marker, function(x) min(floor(x$most_negative_slope))) # The minimum most negative slope

# Get list of phenotypic markers
phenotypic_marker_names_list <- as.list(setNames(phenotypic_markers, phenotypic_markers))

# Get list of NC, no NC, both
filter_nc <- c("Negative Controls Only", "Treatments Only", "All Sparklines")
filter_nc_list <- as.list(setNames(filter_nc, filter_nc))
