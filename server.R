# new shiny server Oct 24, 2015

# install any packages that need to be installed
list.of.packages <- c("shiny", "XLConnect", "plyr", "reshape", "pracma", "car", "Rmisc", "dplyr", "ggvis", "ggplot2", "grid", "tiff", "jpeg", "Rcpp", "rsconnect")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# libraries
library(shiny)
library(XLConnect)
library(plyr)
library(reshape)
library(pracma)
library(car)
library(Rmisc)
library(dplyr)
library(ggvis)
library(ggplot2)
library(grid)
library(tiff)
library(jpeg)
library(Rcpp)
# library(rsconnect)

source("parse_input_excel_files.R", local=TRUE)
source("process_data.R", local=TRUE)
source("explore_compound.R", local=TRUE)
source("explore_target.R", local=TRUE)
source("explore_pathway.R", local=TRUE)
source("explore_QC.R", local=TRUE)
source("explore_QQ_density_plots.R", local=TRUE)
source("explore_early_vs_late_acting.R", local=TRUE)
source("explore_compound_information.R", local=TRUE)
source("get_compound_images.R", local=TRUE)
source("create_results_df.R", local=TRUE)

shinyServer(function(input, output, session) {
  
  observe({
    tryCatch({
      # ---------- parse compound key, compound info, and incucyte data ----------- #
      key_filename <- input$file.key
      compound_info_filename <- input$file.selleck_info
      incucyte_data_filename <- input$file.experimental_results
      trim <- function (x) gsub("^\\s+|\\s+$", "", x)
      image_types <- trim(strsplit(input$image_types, ",")[[1]])
      if (!is.null(key_filename)) { key_data <- parse_compound_key(key_filename) }
      if (!is.null(compound_info_filename)) { input_compound_info <- parse_compound_info(compound_info_filename) }
      if (!is.null(incucyte_data_filename)) { plate_positions <- get_plate_positions(incucyte_data_filename) }
      if (!is.null(incucyte_data_filename) && !is.null(key_data)) { incucyte_data_matrices <- parse_incucyte_data(incucyte_data_filename, key_data) }
      if (!is.null(input$output_dir)) { # make a temporary www directory for the current images
        dir.create(file.path(paste(input$output_dir, "www", sep="")), showWarnings = FALSE)
      }
      
      # --------- process data, add metrics ---------- #
      # if all the files are in place
      if (!is.null(key_filename) && !is.null(compound_info_filename) && !is.null(incucyte_data_filename)) {
        dfs <- process_data(incucyte_data_matrices, input_compound_info, plate_positions)
        data_tall <- dfs[[1]] 
        data_tall_each_marker <- dfs[[2]]
        data_tall_no_NC_each_marker <- dfs[[3]]
        confidence_intervals_each_marker <- dfs[[4]]
        data_wide <- dfs[[5]]
        time_elapsed <- dfs[[6]]
        phenotypic_markers <- unique(data_wide$phenotypic_Marker)
        print("INPUT COMPLETE")
        options(error = browser)
        source("lists.R", local=TRUE)
  
        # update select inputs
        updateSelectInput(session, "overview.nc", choices = filter_nc_list)
        updateSelectInput(session, "overview.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "overview.pathway", choices = overview_pathways_list)
        updateSelectInput(session, "overview.target", choices = overview_targets_list)
        updateSelectInput(session, "overview.compound", choices = overview_compounds_list)
        updateSelectInput(session, "compound", choices = compound_list)
        updateSelectInput(session, "target", choices = target_list)
        updateSelectInput(session, "target.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "pathway", choices = pathway_list)
        updateSelectInput(session, "pathway.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "QC.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "metric.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "metric", choices = metric_names)
        updateSelectInput(session, "early.vs.late.marker", choices = phenotypic_marker_names_list)
        updateSelectInput(session, "above.or.below.NC", choices = above_or_below_list)
        updateSliderInput(session, "time.elapsed", min = head(time_elapsed, n=1), 
                          max = tail(time_elapsed, n=1), 
                          value = head(time_elapsed, n=1), 
                          step= (head(time_elapsed, n=2)[2]-head(time_elapsed, n=2)[1]))
        updateSelectInput(session, "image.type", choices = image_types)
        updateSelectInput(session, "marker", choices = phenotypic_marker_names_list)
        
        # if the input select values are now updated
        if (input$overview.marker != "Loading...") {
          
          # save the data
          save_results(phenotypic_markers, data_tall_no_NC_each_marker, input$output_dir)
          save(data_tall_each_marker, file = paste(input$output_dir, "data_tall_each_marker.R", sep=""))
          save(confidence_intervals_each_marker, file = paste(input$output_dir, "confidence_intervals_each_marker.R", sep=""))

          # Generates a sparkline for the selected compound and phenotypic marker.
          output$compound.sparklines <- renderPlot({
            
            # Get data for specified compound and marker
            data_tall.compound_and_marker <- subset(data_tall_each_marker[[input$marker]], Compound == input$compound)
            
            # Plot data
            plot_sparkline(data_tall.compound_and_marker, confidence_intervals_each_marker[[input$marker]], input$marker, input$compound)
          })
          
          # Generates sparklines for the selected target and phenotypic marker.
          output$target.sparklines <- renderPlot({
            
            # Get data for specified target and marker
            data_tall.target_and_marker <- subset(data_tall_each_marker[[input$target.marker]], target == input$target)
            
            # Plot data
            plot_target_sparklines(data_tall.target_and_marker, confidence_intervals_each_marker[[input$target.marker]], input$target.marker, input$target)
          })
          
          # Generates sparklines for the selected pathway and phenotypic marker.
          output$pathway.sparklines <- renderPlot({
            
            # Get data for specified pathway and marker
            data_tall.target_and_marker <- subset(data_tall_each_marker[[input$target.marker]], Pathway == input$pathway)
            
            # Plot data
            plot_pathway_sparklines(data_tall.target_and_marker, confidence_intervals_each_marker[[input$target.marker]], input$pathway.marker, input$pathway)
          })
          
          # Generates sparklines for target of the selected compound.
          output$compound.target <- renderPlot({
            
            # Get target for specified compound
            data_tall_compound.target <- subset(data_tall_each_marker[[input$marker]], Compound == input$compound)
            
            # Get all data with this target
            data_tall.target <- subset(data_tall_each_marker[[input$marker]], target == data_tall_compound.target$target[1])
            data_tall.target$is.compound <- FALSE
            data_tall.target[which(data_tall.target$Compound == input$compound), ]$is.compound <- TRUE
            
            # Plot data
            plot_target(data_tall.target, confidence_intervals_each_marker[[input$marker]], input$marker, input$compound, data_tall_compound.target$target[1])
          })
          
          # Generates sparklines for pathway of the selected compound.
          output$compound.pathway <- renderPlot({
            
            # Get pathway for specified compound
            data_tall_compound.pathway <- subset(data_tall_each_marker[[input$marker]], Compound == input$compound)
            pathway <- data_tall_compound.pathway$Pathway[1]
            
            # Get all data with this pathway
            data_tall.pathway <- subset(data_tall_each_marker[[input$marker]], Pathway == pathway)
            data_tall.pathway$is.compound <- FALSE
            data_tall.pathway[which(data_tall.pathway$Compound == input$compound), ]$is.compound <- TRUE
            
            # Plot data
            plot_pathway(data_tall.pathway, confidence_intervals_each_marker[[input$marker]], input$marker, input$compound, pathway)
          })
          
          # Hilights sparklines for target in pathways.
          output$target.pathway <- renderPlot({
            
            # Get all data with target hilights
            data_tall.pathways <- data_tall_each_marker[[input$target.marker]]
            data_tall.pathways$is.target <- FALSE
            data_tall.pathways[which(data_tall.pathways$target == input$target), ]$is.target <- TRUE
            
            plot_pathway_hilight_target(data_tall.pathways, confidence_intervals_each_marker[[input$target.marker]], input$target.marker, input$target, confidence_intervals_each_marker)
          })
          
          # Generates clusters, hilighting selected compound.
          output$compound.cluster <- renderPlot({
            
            # Mark the compound within the data frame
            data_tall.cluster <- get_data_w_clusters(data_wide, data_tall_each_marker[[input$marker]], input$marker, input$clusters)
            data_tall.cluster$is.compound <- FALSE
            data_tall.cluster[which(data_tall.cluster$Compound == input$compound), ]$is.compound <- TRUE
            
            # Plot data
            plot_clusters(data_tall.cluster, input$compound, input$marker, confidence_intervals_each_marker)
          })
          
          # Hilights sparklines for target in clusters.
          output$target.cluster <- renderPlot({
            
            # Get all data with target hilights
            data_tall.clusters <- get_data_w_clusters(data_wide, data_tall_each_marker[[input$target.marker]], input$target.marker, input$target.clusters)
            data_tall.clusters$is.target <- FALSE
            data_tall.clusters[which(data_tall.clusters$target == input$target), ]$is.target <- TRUE
            
            plot_target_clusters(data_tall.clusters, input$target, input$target.marker, confidence_intervals_each_marker)
          })
          
          # Hilights sparklines for pathway in clusters.
          output$pathway.cluster <- renderPlot({
            
            # Get all data with pathway hilights
            data_tall.clusters <- get_data_w_clusters(data_wide, data_tall_each_marker[[input$pathway.marker]], input$pathway.marker, input$pathway.clusters)
            data_tall.clusters$is.pathway <- FALSE
            data_tall.clusters[which(data_tall.clusters$Pathway == input$pathway), ]$is.pathway <- TRUE
            
            plot_pathway_clusters(data_tall.clusters, input$pathway, input$pathway.marker, confidence_intervals_each_marker)
          })
          
          curr_compound <- NULL
          
          output$display.image <- renderImage({
            
            # Get images for that compound into www folder of shiny app
            image_types <- setNames(image_types, image_types)
            suppressWarnings(get_images(input$compound, input$archive_dir, data_wide, image_types, time_elapsed, input$output_dir))
            
            position <- data_wide[which(data_wide$Compound == input$compound),]$Position[1] # Position of compound in plate
            plate <- data_wide[which(data_wide$Compound == input$compound),]$Plate[1] # Plate of compound
            
            image_file <- paste(input$output_dir,"www/Plate",plate, "_Position", position, "_Image", image_types[[input$image.type]],"_t_",input$time.elapsed,".jpeg",sep="")
            
            
            return(list(
              src = image_file,
              filetype = "image/jpeg",
              height = 520,
              width = 696
            ))
            
          }, deleteFile = FALSE)
          
          
          # Plot the negative controls vs treatment for specified phenotypic marker
          output$QC.by.plate <- renderPlot({
            
            plot_QC_by_plate(data_tall_each_marker[[input$QC.marker]], input$QC.marker)
            
          })
          
          # QQ plot for specified phenotypic marker and curve metric
          output$qq.plot.one.metric <- renderPlot({
            
            get_single_metric_qqplot(data_tall_each_marker[[input$metric.marker]], metrics[[input$metric]])
            
          })
          
          # QQ plot for specified phenotypic marker and curve metric
          output$density.plot.one.metric <- renderPlot({
            
            get_single_metric_density(data_tall_each_marker[[input$metric.marker]], metrics[[input$metric]])
            
          })
          
          # QQ plot for specified phenotypic marker and curve metric
          output$early.vs.late.acting <- renderPlot({
            
            get_early_vs_late_acting(data_tall_each_marker[[input$early.vs.late.marker]], confidence_intervals_each_marker[[input$early.vs.late.marker]], input$early.vs.late.marker, input$above.or.below.NC)
            
          }, height = 800, width = 1000 )
          
          ### Overview of sparklines ###
          # Subset data for phenotypic marker of interest
          selected <- reactive({
              data <- NA
              if (input$overview.nc == "Treatments Only") {
                data <- data_tall_no_NC_each_marker[[input$overview.marker]]
              }
              if (input$overview.nc == "Negative Controls Only") {
                data <- data_tall_each_marker[[input$overview.marker]]
                data <- subset(data, Pathway == "NegControl")
              }
              if (input$overview.nc == "All Sparklines") {
                data <- data_tall_each_marker[[input$overview.marker]]
              }
              if (input$overview.target != "All Targets") {
                data <- subset(data, target == input$overview.target)
              }
              if (input$overview.pathway != "All Pathways") {
                data <- subset(data, Pathway == input$overview.pathway)
              }
              if (input$overview.compound != "All Compounds") {
                data <- subset(data, Compound == input$overview.compound)
              }
              data$description = as.character(paste0(data$Compound))
              data
          })
          
          # Plot all compound sparklines for the selected phenotypic marker
          selected %>%
            group_by(Compound) %>%
            ggvis(x = ~as.numeric(time_elapsed), y = ~as.numeric(phenotype_value)) %>%
            layer_lines(stroke.hover := "red") %>%
            add_axis("x", title = "Time Elapsed") %>%
            add_axis("y", title = "Phenotype Value") %>%
            add_tooltip(function(data){
              paste0(as.character(data$Compound))
            }, "hover") %>%
            bind_shiny("allsparklines")
          
          
          # All information for the specified compound and phenotypic marker
          output$compound.additional_info <- renderUI({
            get_compound_info(input, data_tall_each_marker)
          })
        }
      }
    }, error = function(err) {
      print(err)
    })
  })
})
