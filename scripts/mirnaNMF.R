#This program reads in mirna expression data and performs NMF clustering.

# INPUT
# 1. expression data as tsv file - miRNAs along rows, samples along columns
# 2. sample names as tsv file
# 3. miRNA names as tsv file
# 4. rank of NMF output (to be discontinued after implementing rank estimation)
# 5. number of runs of NMF


library(bigrquery)
library(NMF)
library(tidyr)

#NOTE: no BigQueries are being run currently, so billingProject is not being used for now.
billingProject = "isb-cgc-04-0010"

#command line arguments
args = commandArgs()
print(args)

#NOTE: Following commented section gets data from BigQuery. For now, we will read in a tsv instead.
#study='BRCA'
#sample type set to TP (primary tumor)
#querySql = paste("SELECT ParticipantBarcode,mirna_id,mirna_accession ,normalized_count",
#                 " FROM [isb-cgc:tcga_201607_beta.miRNA_Expression]",
#                 " WHERE SampleTypeLetterCode='TP' AND Study='",study,"'",sep = "")
#querySql = "SELECT * FROM [isb-cgc-04-0010:LargeResults.mirna_BRCA_Test]"
#resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)

#Read tsv into a data frame
dataDF = read.table(args[6],sep = "\t",header = F)
sampleNames = read.table(args[7],sep = "\t",header = F)
mirnaNames = read.table(args[8],sep = "\t",header = F)
colnames(dataDF) = as.vector(sampleNames$V1)
rownames(dataDF) = as.vector(do.call(rbind,strsplit(as.vector(mirnaNames$V1),":",fixed = T))[,8])

#dummy sample labels to try annotation tracks in NMF heatmaps
labels = rep(c("A","B","C","D"),ncol(dataDF)*c(0.1,0.2,0.65,0.05))
labels = sample(labels,ncol(dataDF))

#NOTE: mirnaCluster is a function in GDACmiRNA package,
#but if does nothing more that call in-built NMF functions as shown inline in the following code

#Perform NMF
#mirnaCluster(dataDF)

#Run NMF nrun number of times and return the best fit (use keep.all=TRUE to return a list of all runs)
# NOTE:
# 1. Parallel runs - as per documentation, default is to run in parallel using all cores
# To force parallel computation, use .opt="P" or .opt="P<n>" to specify number of cores to use (do this on command line)
# 2. Seed - specified single numeric seed; can look into other methods of generating seeds

nmfFit = nmf(dataDF,rank=as.numeric(args[9]),method = "brunet",seed=123456,nrun=as.numeric(args[10]),.opt="v")

#output to stdout
showRNG(nmfFit)
extractFeatures(nmfFit)

#heatmaps to be saved in a pdf
pdf("NMFHeatMaps.pdf")
basismap(nmfFit)
coefmap(nmfFit,annCol=labels,main="Coefficient heatmap\nWith default tracks")
coefmap(nmfFit,annCol=labels,main="Coefficient heatmap\nWithout default tracks",tracks=NULL)
consensusmap(nmfFit)
dev.off()
