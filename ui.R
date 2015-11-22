# new shiny ui Oct 24, 2015

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
library(RcppEigen)
library(tiff)
library(jpeg)
library(rsconnect)

shinyUI(navbarPage("Perspective:",
  tabPanel("Inputs",
     fileInput("file.key", "Compounds Key File", multiple = FALSE, accept = c("xlsx")),
     fileInput("file.selleck_info", "Selleck Compound Information File (March, 2015 Edition)", multiple = FALSE, accept = c("xlsx")),
     fileInput("file.experimental_results", "Experimental Results (Reconfigured) File", multiple = FALSE, accept = c("xlsx")),
     textInput("image_types", "Enter the types of images output by the IncuCyte machine, as they appear in the archive image file: ", value = "P, C1"),
     textInput("archive_dir", "Enter the image archive directory: ", value = "/Volumes/G-Drive/dc140908 c2c12 diff tun 1-5/EssenFiles/ScanData/"),
     textInput("output_dir", "Enter the directory for output: ", value = "/Users/maiasmith/Desktop/Output/")
  )
  ,
  tabPanel("Overview",
           
     # Application title
     titlePanel("Interactive Plot"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("overview.nc", 'Treatment vs Negative Control (N.C.)', choices = list("Loading..."), multiple = FALSE, selected = "All Sparklines"),
         selectizeInput("overview.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("overview.pathway", 'Pathway', choices = list("Loading..."), multiple = FALSE), 
         selectizeInput("overview.target", 'Target', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("overview.compound", 'Compound', choices = list("Loading..."), multiple = FALSE)
       ),
       mainPanel(
         # ggvisOutput("allsparklines")
         )
     )
  )
  ,
  tabPanel("Compound",
           
           # Application title
           titlePanel("Compound Perspective"),
           
           # Layout
           sidebarLayout(
             
             sidebarPanel(
               selectizeInput("compound", 'Compound of Interest (Select or Type)', choices = list("Loading..."), multiple = FALSE),
               selectizeInput("marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE)
             ),
             
             mainPanel(
               tabsetPanel(
                 tabPanel("Live Images", 
                          h3("Live images of cells responding to a compound over time."),
                          imageOutput("display.image"),
                          br(), br(), br(), br(), br(), br(), br(), br(), br(),
                          sliderInput("time.elapsed", "Time Elapsed:", 
                                      min = -10, 
                                      max = 0, 
                                      value = -1, 
                                      step= 1),
                          selectizeInput("image.type", 'Image Type', choices = list("Loading..."), multiple = FALSE)
                  ),
                 tabPanel("Information", 
                          h3("Information for the selected compound."),
                          htmlOutput("compound.additional_info")),
                 tabPanel("Sparkline", plotOutput("compound.sparklines")), 
                 tabPanel("Target", plotOutput("compound.target")),
                 tabPanel("Pathway", plotOutput("compound.pathway")),
                 tabPanel("Clusters", plotOutput("compound.cluster"),
                          sliderInput("clusters", "Number of Clusters:", min = 1, max = 25, value = 10))
               )
             )
           )
  )
  ,
  tabPanel("Target",
     
     # Application title
     titlePanel("Target Perspective"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("target", 'Target of Interest (Select or Type)', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("target.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE),
         sliderInput("target.clusters","Number of Clusters:", min = 1, max = 25, value = 10)
       ),
       
       mainPanel(
         tabsetPanel(
           tabPanel("Sparklines", plotOutput("target.sparklines")), 
           tabPanel("Pathways", plotOutput("target.pathway")),
           tabPanel("Clusters", plotOutput("target.cluster"))
         )
       )
     )
  ),
  tabPanel("Pathway",
     
     # Application title
     titlePanel("Pathway Perspective"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("pathway", 'Pathway of Interest (Select or Type)', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("pathway.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE),
         sliderInput("pathway.clusters", "Number of Clusters:", min = 1, max = 25, value = 10)
       ),
       mainPanel(
         tabsetPanel(
           tabPanel("Sparklines", plotOutput("pathway.sparklines")),
           tabPanel("Clusters", plotOutput("pathway.cluster"))
         )
       )
     )
  ),
  tabPanel("QC",
           
     # Application title
     titlePanel("Quality Control"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("QC.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE)
       ),
       mainPanel(plotOutput("QC.by.plate"))
     )
  ),
  tabPanel("Curve Metrics",
           
     # Application title
     titlePanel("QQ Plots and Density Plots of Curve Metrics"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("metric.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("metric", 'Metric', choices = list("Loading..."), multiple = FALSE)
       ),
       mainPanel(
         tabsetPanel(
           tabPanel("QQ Plot", plotOutput("qq.plot.one.metric")),
           tabPanel("Density Plot", plotOutput("density.plot.one.metric"))
         )
       )
     )
  ),
  tabPanel("Timing",
           
     # Application title
     titlePanel("Early- vs Late-Acting Compounds"),
     
     # Layout
     sidebarLayout(
       sidebarPanel(
         selectizeInput("early.vs.late.marker", 'Phenotypic Marker', choices = list("Loading..."), multiple = FALSE),
         selectizeInput("above.or.below.NC", 'With respect to negative control confidence interval:', choices = list("Loading..."),multiple = FALSE)
       ),
       mainPanel(plotOutput("early.vs.late.acting"))
     )
  )
))