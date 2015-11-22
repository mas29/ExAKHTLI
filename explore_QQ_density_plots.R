###############################################################################
###### ------ Plot QQ plot and histogram for one curve metric ---------- ######
###############################################################################

# Get density plot for one metric

# @param df -- data frame of data with metrics 
# @param metric_name -- name of metric to explore
get_single_metric_density <- function(df, metric_name) {
  # Distributions to compare
  metric_NC <- df[df$empty == "Negative Control", colnames(df) == metric_name]
  metric_Treatment <- df[df$empty == "Treatment", colnames(df) == metric_name]
  
  # Plot
  ggplot() + 
    geom_density(data = data.frame(metric_NC), aes(x = metric_NC, fill = 'Negative Control'), alpha = 0.5) + 
    geom_density(data = data.frame(metric_Treatment), aes(x = metric_Treatment, fill = 'Treatment'), alpha = 0.5) + 
    xlab(metric_name) +
    guides(fill=guide_legend(title="Legend", direction="vertical")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.key.width = unit(0.5, "cm")) 
}

# Get qqplot for one metric

# @param df -- data frame of data with metrics 
# @param metric_name -- name of metric to explore
get_single_metric_qqplot <- function(df, metric_name) {
  # Distributions to compare
  metric_NC <- df[df$empty == "Negative Control", colnames(df) == metric_name]
  metric_Treatment <- df[df$empty == "Treatment", colnames(df) == metric_name]
  
  # Calculated quantiles
  q1 <- quantile(scale(metric_NC),0:100/100)
  q2 <- quantile(scale(metric_Treatment),0:100/100)
  
  # Plot
  ggplot(data=data.frame(a=q1,b=q2)) + 
    geom_point(aes(x=a,y=b)) +
    geom_abline(intercept=0,slope=1) +
    theme_bw() +
    xlab("Negative Control") +
    ylab("Treatment") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 

}
