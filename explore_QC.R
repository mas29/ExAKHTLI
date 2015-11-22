################################################################################################
###### ------ Plot sparklines for phenotypic marker of interest and each plate ---------- ######
################################################################################################

# Function to plot QC of each plate for phenotypic marker of interest
# param df -- data in tall format for particular phenotypic marker
# param phenotypic_marker_name -- phenotypic marker of interest
plot_QC_by_plate <- function(df, phenotypic_marker_name) {
  ggplot(df, 
       aes(x=as.numeric(time_elapsed), y=as.numeric(phenotype_value), group=Compound)) +
  geom_line() +
  xlab("Time Elapsed") +
  ylab(phenotypic_marker_name) +
  ggtitle(paste("Quality Control\n",phenotypic_marker_name," Data Separated By Plate",sep="")) +
  facet_grid(empty ~ Plate, scales = "fixed") +
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        panel.background = element_rect(fill = "white"))
}


