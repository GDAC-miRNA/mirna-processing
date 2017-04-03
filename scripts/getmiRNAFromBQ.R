#get miRNA data from bigquery for a list of samples and a list of miRNAs
#Usage:
#R getmiRNAData.R tableToQuery tumorType

library(googleAuthR)
library(bigQueryR)
library(googleCloudStorageR)
library(utils) #for downloading files from GCS using their URL
library(httr) #for content()
library(tidyr) #for spread


options("googleAuthR.scopes.selected" = c("https://www.googleapis.com/auth/bigquery","https://www.googleapis.com/auth/devstorage.full_control"))
service_token = gar_auth_service(json_file = "/Users/varsha/Documents/miRNA/scripts/ServiceAccountKey/miRNA project-c2d2cf5124f1.json")
#args = commandArgs(trailingOnly = T)
#parse args
# queryProject = ""
# queryDS = ""
# queryTable = ""
#queryTumorType = NULL
# queryColumns = "" expected in order of row IDs, columns IDs, values
directDownload = F
for(thisArg in commandArgs(trailingOnly = T))
{
  print(thisArg)
  if(startsWith(thisArg,"project"))
  {
    queryProject = unlist(strsplit(thisArg,split = '=',fixed = T))[2]
  }
  else if(startsWith(thisArg,"ds"))
  {
    queryDS = unlist(strsplit(thisArg,split='=',fixed = T))[2]
  }
  else if(startsWith(thisArg,"table"))
  {
    queryTable = unlist(strsplit(thisArg,split='=',fixed = T))[2]
  }
  else if(startsWith(thisArg,"tumor"))
  {
    queryTumorType = unlist(strsplit(thisArg,split='=',fixed = T))[2]
  }
  else if(startsWith(thisArg,"columns"))
  {
    queryColumns = unlist(strsplit(thisArg,split='=',fixed = T))[2]
  }
  else if(startsWith(thisArg,"--directDownload"))
  {
    #'T' mean directl download to a data frame, 
    #'F' implies large results and requires the query result to be stored on GCS and downloaded as smaller chunks
    directDownload =  T
  }
  
}
#TODO::assert valid values in input arguments

fullTableName = paste('[',queryProject,':',queryDS,'.',queryTable,']',sep="")

query = paste('SELECT',queryColumns,'FROM',fullTableName,ifelse(!is.null(queryTumorType),paste("WHERE Study='",queryTumorType,"'",sep=""),""),sep=" ")
print(query)  

## run query
if (directDownload) #limit to results with 100000 rows
{
  queryResult = bqr_query(projectId = queryProject,datasetId = queryDS,query,maxResults = 10000)
  print(class(queryResult))
  print(dim(queryResult))
  print(head(queryResult))
  resultMatrix = spread_(queryResult,unlist(strsplit(queryColumns,','))[2],unlist(strsplit(queryColumns,','))[3])
  print(resultMatrix[1:5,1:10])
  write.table(resultMatrix,file = "BigQueryResultMatrix.tsv",sep="\t")
  
}else
{
  #get unique row and column IDs from the expected query result
  query = paste('SELECT unique(',unlist(strsplit(queryColumns,','))[1],') FROM',fullTableName,ifelse(!is.null(queryTumorType),paste("WHERE Study='",queryTumorType,"'",sep=""),""),sep=" ")
  uniqueRowIDs = bqr_query(projectId = queryProject,datasetId = queryDS,query,maxResults = 10000)
  uniqueRowIDs = as.vector(uniqueRowIDs[,1])
  
  query = paste('SELECT unique(',unlist(strsplit(queryColumns,','))[2],') FROM',fullTableName,ifelse(!is.null(queryTumorType),paste("WHERE Study='",queryTumorType,"'",sep=""),""),sep=" ")
  uniqueColumnIDs = bqr_query(projectId = queryProject,datasetId = queryDS,query,maxResults = 10000)
  uniqueColumnIDs = as.vector(uniqueColumnIDs[,1])
  
  resultMatrix = data.frame(matrix(ncol = length(uniqueColumnIDs),nrow = length(uniqueRowIDs)))
  rownames(resultMatrix) = uniqueRowIDs
  colnames(resultMatrix) = uniqueColumnIDs
  print(dim(resultMatrix))
  print(resultMatrix[1:5,1:5])
  
  query = paste('SELECT',queryColumns,'FROM',fullTableName,ifelse(!is.null(queryTumorType),paste("WHERE Study='",queryTumorType,"'",sep=""),""),sep=" ")
  print(query) 
  
  
  job <- bqr_query_asynch(queryProject,
                       queryDS,
                       query,
                       destinationTableId = "bigResultTable")

  print('Fetching results into bigResultTable...')
  # poll the job to check its status
  # its done when job$status$state == "DONE"
  while((job = bqr_get_job(queryProject, job$jobReference$jobId))$status$state != "DONE"){}
  if(is.null(job$status$errors))
  {
    print('DONE!')
  }else
  {
    print(job$status$errors)
    #return(NULL) when making this into a function
  }
  ##once done, the query results are in "bigResultTable"

  #write to google cloud storage
  ## Create the data extract from BigQuery to Cloud Storage
  job_extract <- bqr_extract_data(queryProject,
                                  queryDS,
                                  "bigResultTable",
                                  "temp_for_bq_load")

  ## poll the extract job to check its status
  ## its done when job$status$state == "DONE"
  print('Writing to Google Cloud Storage...')
  while ((job_extract = bqr_get_job(queryProject, job_extract$jobReference$jobId))$status$state != "DONE" ) {}
  if(is.null(job_extract$status$errors))
  {
    print('DONE!')
  }else
  {
    print(job$status$errors)
    #return(NULL) when making this into a function
  }

  gcs_global_bucket("temp_for_bq_load")
  ## get object info in the default bucket
  objects <- gcs_list_objects()
  objects = objects[grep(pattern = "big-query-extract",x = objects$name),]

  ## save directly to an R object (warning, don't run out of RAM if its a big object)
  ## the download type is guessed into an appropriate R object
  thisDF = data.frame()
  for(item in objects$name)
  {
 
    print(paste('Processing',item,sep=" "))
    response <- gcs_get_object(item,meta = FALSE,parseObject = F)
    thisDF = rbind(thisDF,as.data.frame(content(response,as="parsed",type = "text/csv")))
    print(head(thisDF))
    
    # #populate resultMatrix
    # for(thisRowIndex in 1:nrow(resultMatrix))
    # {
    #   thisSampleDF = thisDF[thisDF[[1]] == rownames(resultMatrix)[thisRowIndex],]
    #   uniqueColumnNames = unique(thisSampleDF[,2])
    #   #compute mean for repeated column names (doesn't always occur)
    #   meanValues = aggregate(thisSampleDF[,3],thisSampleDF[,2],mean)
    #   resultMatrix[thisRowIndex,uniqueColumnNames[[1]]] = meanValues[match(uniqueColumnNames[[1]],meanValues[[1]]),2]
    # }
    
    #delete this object from GCS
    #gcs_delete_object(item)
  }
  #delete bigResultTable
  bqr_delete_table(queryProject,queryDS,"bigResultTable")
  #return(resultMatrix) wen making this into a function

}
