#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(plotly)
library(googleAuthR)

options("googleAuthR.scopes.selected" = c("https://www.googleapis.com/auth/bigquery"))
#options("googleAuthR.webapp.client_id" = "907668440978-3qbm41pompmp3m2uto4aoogklr1patbt.apps.googleusercontent.com")
#options("googleAuthR.webapp.client_secret" = "qFS-5BQPg-LbZT0V-U9Y37jo")

service_token <- gar_auth_service(json_file="./data/miRNA project-c2d2cf5124f1.json")
responseDF = data.frame()
drawPlot = function(mirID,study,colorByGroup,projectID)
{
  sqlQuery = paste("SELECT expr.AnalyteBarcode as AnalyteBarcode, expr.Study as Study,hg19.mirna_id as mirnaID, expr.hg19_RPM as hg19_RPM, expr.hg38_RPM as hg38_RPM, tech.BCGSC_protocol as BCGSC_protocol,tech.BCR_protocol as BCR_protocol,tech.BCGSC_platform as BCGSC_platform",
                   "FROM  (",
                   "SELECT LEFT(hg19.aliquot_barcode,20) as AnalyteBarcode, hg19.project_short_name as Study, hg19.mirna_id , hg19.reads_per_million_miRNA_mapped as hg19_RPM, hg38.reads_per_million_miRNA_mapped as hg38_RPM", 
                   "FROM (",
                   "SELECT *",
                   "FROM [isb-cgc-04-0010:draft_new_data.TCGA_gdc_hg19_miRNA_Expression] ) as hg19",
                   "JOIN (",
                   "SELECT *", 
                   "FROM [isb-cgc-04-0010:draft_new_data.TCGA_gdc_hg38_miRNA_Expression]",
                   "WHERE mirna_id =",shQuote(mirID),") as hg38",
                   "ON hg19.aliquot_barcode = hg38.aliquot_barcode AND hg19.mirna_id = hg38.mirna_id) as expr",
                   "JOIN (",
                   "SELECT *",
                   "FROM [isb-cgc-04-0010:reference_data.miRNA_protocol_platform_info] ) as tech",
                   "ON expr.AnalyteBarcode = tech.analyte_barcode",sep=" ")
  

  body = list(
    query=sqlQuery,
    defaultDataset.datasetId="draft_new_data",
    defaultDataset.projectId=projectID
  )
  
  f = gar_api_generator("https://www.googleapis.com/bigquery/v2",
                        "POST",
                        path_args = list(projects=projectID,queries=""))
  
  response = f(the_body=body)
  
  responseDF = as.data.frame(do.call("rbind",lapply(response$content$rows$f,FUN = t)))
  
  row.names(responseDF) = NULL
  colnames(responseDF) = c("AnalyteBarcode","Study","mirna_id","hg19_RPM","hg38_RPM","BCGSC_protocol","BCR_protocol","BCGSC_platform")
  #convert NA string to NA value
  responseDF[responseDF=='NA'] = NA
  responseDF$hg19_RPM = log2(as.numeric(as.vector(responseDF$hg19_RPM))+1)
  responseDF$hg38_RPM = log2(as.numeric(as.vector(responseDF$hg38_RPM))+1)
  #edit responseDF$Study to chop off TCGA prefix  
  responseDF$Study = substring(responseDF$Study,6)

  if(study=="Pan-Cancer")
  {
   
    spearman = cor(responseDF$hg19_RPM,responseDF$hg38_RPM,method = "spearman")
    plotDF = responseDF
    plotDF[,'colorByColumn'] = responseDF[,colorByGroup]
    p=ggplot(data=plotDF,aes(hg19_RPM,hg38_RPM,text=AnalyteBarcode,colour=factor(colorByColumn)))+
      geom_point(size=2,shape=21,stroke=2) +
      geom_abline(slope=1,intercept = 0,aes(size=2)) + 
      #geom_smooth(method = "loess",se=FALSE,colour="red") +
      labs(x="hg19 log2(RPM)",y="hg38 log2(RPM)") +
      ggtitle(paste(mirID," (rho = ",ifelse(is.numeric(spearman),round(spearman,4),"NA"),")",sep="")) + 
      theme(plot.title= element_text(size=10,face="bold",hjust = 0.3,debug = T),axis.title = element_text(size=10),axis.text = element_text(size=10)) +
      scale_color_discrete(name =colorByGroup)
    p = ggplotly(p,tooltip='text') 
    return(p)
   }  
  else if(startsWith(study,"All"))
    {
  
      spearman = by(responseDF,responseDF$Study,function(x) cor(x$hg19_RPM,x$hg38_RPM,method = "spearman"))
      spearman = as.data.frame(as.vector(spearman))
      rownames(spearman) = levels(as.factor(responseDF$Study))
      plotDF = responseDF
      plotDF[,'colorByColumn'] = responseDF[,colorByGroup]
      #edit study names to include spearman coeffs
      #responseDF$Study = vapply(responseDF$Study,function(x) paste(x,"(rho = ",round(spearman[x,1],4),")",sep=""),FUN.VALUE = character(1))
      p=ggplot(data=plotDF,aes(hg19_RPM,hg38_RPM,text=AnalyteBarcode,colour=factor(colorByColumn)))+
        geom_point(shape=21,stroke=2,size=1) +
        geom_abline(slope=1,intercept = 0,aes(size=2)) + 
        #geom_smooth(method = "loess",se=FALSE,colour="red") + 
        facet_wrap(~Study,ncol=6) +
        labs(title=mirID,x="hg19 log2(RPM)",y="hg38 log2(RPM)") + 
        theme(plot.title = element_text(size=10,face="bold",hjust = 0.5),axis.title = element_text(size=10),axis.text = element_text(size=8)) +
        scale_color_discrete(name=colorByGroup)
      p = ggplotly(p,tooltip='text')
      return(p)
    }else
    {
  
      responseDF = responseDF[responseDF$Study==study,]
      spearman = cor(responseDF$hg19_RPM,responseDF$hg38_RPM,method = "spearman")
      plotDF = responseDF
      plotDF[,'colorByColumn'] = responseDF[,colorByGroup]
      p=ggplot(data=plotDF,aes(hg19_RPM,hg38_RPM,text=AnalyteBarcode,colour=factor(colorByColumn)))+
        geom_point(shape=21,stroke=2,size=2) +
        geom_abline(slope=1,intercept = 0,aes(size=2)) + 
        #geom_smooth(method = "loess",se=FALSE,colour="red") +
        labs(x="hg19 log2(RPM)",y="hg38 log2(RPM)") +
        ggtitle(paste(mirID," (rho = ",ifelse(is.numeric(spearman),round(spearman,4),"NA"),")",sep="")) + 
        theme(plot.title = element_text(size=10,face="bold",hjust = 0.5),axis.title = element_text(size=10),axis.text = element_text(size=10)) + 
        scale_color_discrete(name =colorByGroup)
      p = ggplotly(p,tooltip='text')
      return(p)
    }
 
}


# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  
  projectID = "isb-cgc-04-0010" 
  outputPlot = eventReactive(input$submit,{
    drawPlot(mirID = tolower(input$mirna),study = input$study,colorByGroup = input$tech,projectID = projectID)
  
  })
  
  output$distPlot = renderPlotly({outputPlot()}) 
  
  
})