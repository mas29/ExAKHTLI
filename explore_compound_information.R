#####################################################################################
###### ------------------ Get information for this compound ------------------ ######
#####################################################################################

# Function to get information for the compound and phenotypic marker as specified in the Shiny app
# param input -- input from user Shiny selections
# param data_tall_each_marker -- list of data for each marker
get_compound_info <- function(input, data_tall_each_marker) {
  df <- data_tall_each_marker[[input$marker]][data_tall_each_marker[[input$marker]]$Compound == input$compound, ]
  HTML(paste(paste("Compound: ", df$Compound[1], sep=""),
             paste("Phenotypic Marker: ", input$marker, sep=""),
             paste("Catalog No.: ", df$Catalog.No.[1], sep=""),
             paste("M.w.: ", df$M.w.[1], sep=""),
             paste("Rack Number: ", df$Rack.Number[1], sep=""),
             paste("CAS Number: ", df$CAS.Number[1], sep=""),
             paste("Form: ", df$Form[1], sep=""),
             paste("Target Class: ", df$Target.class..11Mar15.[1], sep=""),
             paste("Target Type: ", df$Target.Type[1], sep=""),
             paste("Target Species: ", df$Target.Species[1], sep=""),
             paste("Molecule Type: ", df$Molecule.Type[1], sep=""),
             paste("Information: ", df$Information[1], sep=""),
             paste("Smiles: ", df$Smiles[1], sep=""),
             paste("URL: ", df$URL[1], sep=""),
             paste("Pathway: ", df$Pathway[1], sep=""),
             paste("Plate: ", df$Plate[1], sep=""),
             paste("Position: ", df$Position[1], sep=""),
             paste("Screen: ", df$Screen[1], sep=""),
             paste("Minimum Value: ", df$min[1], sep=""),
             paste("Mean Value: ", df$mean[1], sep=""),
             paste("Max Value: ", df$max[1], sep=""),
             paste("Time To Minimum Value: ", df$time_to_min[1], sep=""),
             paste("Time To Maximum Value: ", df$time_to_max[1], sep=""),
             paste("AUC (Trapezoidal Integration): ", df$AUC_trapezoidal_integration[1], sep=""),
             paste("Delta (max-min): ", df$delta_min_max[1], sep=""),
             paste("Most Positive Slope: ", df$most_positive_slope[1], sep=""),
             paste("Most Negative Slope: ", df$most_negative_slope[1], sep=""),
             paste("Time To Most Positive Slope: ", df$time_to_most_positive_slope[1], sep=""),
             paste("Time To Most Negative Slope: ", df$time_to_most_negative_slope[1], sep=""),
             paste("Time X Distance to Upper Negative Control Confidence Interval: ", df$time_x_distance.upper[1], sep=""),
             paste("Time X Distance to Lower Negative Control Confidence Interval: ", df$time_x_distance.lower[1], sep=""),
             paste("Timepoint at which Value Exceeds Negative Control Upperbound Confidence Interval: ", df$phenotype_value_exceeds_NC_upperbound.timepoint[1], sep=""),
             paste("Timepoint at which Value Falls Below Negative Control Lower Confidence Interval: ", df$phenotype_value_falls_below_NC_lowerbound.timepoint[1], sep=""),
             sep = '<br/>'))
}
