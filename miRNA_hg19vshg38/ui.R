#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

mirnaList = read.table('./data/hg19_hg38_miRNAOverlapList.txt',header = F)
studyList = read.table('./data/StudyList.txt',header = F)
studyList = as.vector(studyList$V1)
studyList = append(studyList,c('All Tumor types','Pan-Cancer'),after = length(studyList))
studyList = rev(studyList)
# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  
  titlePanel("miRNA Expression - hg19 vs. hg38"),
  
  sidebarLayout
  (
    
    sidebarPanel
    (
      selectInput("mirna", "Select a miRNA",choices=mirnaList$V1,width = "75%"),
      selectInput("study",'Select a tumor type',choices = studyList, width = "75%"),
      radioButtons("tech","Color by platform/protocol: ",choices = c('BCR Protocol'='BCR_protocol','BCGSC Platform' = 'BCGSC_platform','BCGSC Protocol'='BCGSC_protocol'),selected = "BCR_protocol"),
      actionButton(inputId="submit",label = "Submit")
      
    ),
    
    mainPanel 
    (
      plotlyOutput("distPlot",height = "700px")
    )
  )
 
)
)