# This R program computes Spearman correlation coefficients between miRNAs and their target genes.

# R packages 
library(bigrquery)

#miRNA cloud project to be billed for the BigQuery queries 
billingProject = 'isb-cgc-04-0010' 

###############################################################
#Some exploratory analyses to follow...

#1. Get number of target genes per miRNA from miRTarBase v6.1 (only Functional MTI)
querySql = paste("SELECT miRNA,count(Target_Gene ) AS N",
                 "FROM (",
                 "SELECT  miRNA ,Target_Gene ,COUNT(*)",
                 "FROM [isb-cgc-04-0010:miRTarBase_6_1.miRTarBase_MT1]", 
                 "WHERE Species_Target_Gene = 'Homo sapiens' AND Species_miRNA = 'Homo sapiens' AND Support_Type = 'Functional MTI'",
                 "GROUP BY miRNA ,Target_Gene )",
                 "GROUP BY miRNA",sep = " ")

resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)
bxp = boxplot(resultDF$N,xlab=paste("No. of miRNAs = ",nrow(resultDF), "\nNo. of miRNA-mRNA pairs = ",sum(resultDF$N)),ylab="Number of Target genes",main="miRTarBase v6.1")
text(x=c(1.25,1.25,1.25,1.25,1.25),y=bxp$stats,labels = bxp$stats,col = "darkgray")
text(x=c(0.75,0.75,0.75,0.75),y=sort(bxp$out,decreasing = T)[1:4],labels = resultDF$l_mirna_id[order(resultDF$N,decreasing = T)[1:4]],col = "darkgray")
text(x=c(1.15,1.15,1.15,1.15),y=sort(bxp$out,decreasing = T)[1:4],labels = sort(resultDF$N,decreasing = T)[1:4],col = "darkgray")

#2. Get number of target genes per miRNA from TCGA miRNA-mRNA join (only Functional MTI)
querySql = "SELECT  l_mirna_id ,count(DISTINCT Target_Gene ) as N FROM [isb-cgc-04-0010:LargeResults.miRNA_mRNA_Exp] group by l_mirna_id order by N desc"

resultDF = query_exec(querySql,project = billingProject,max_pages = Inf)
bxp = boxplot(resultDF$N,xlab=paste("No. of miRNAs = ",nrow(resultDF), "\nNo. of miRNA-mRNA pairs = ",sum(resultDF$N)),ylab="Number of Target genes",main="miRTarBase v6.1")
text(x=c(1.25,1.25,1.25,1.25,1.25),y=bxp$stats,labels = bxp$stats,col = "darkgray")
text(x=c(0.75,0.75,0.75,0.75),y=sort(bxp$out,decreasing = T)[1:4],labels = resultDF$l_mirna_id[order(resultDF$N,decreasing = T)[1:4]],col = "darkgray")
text(x=c(1.15,1.15,1.15,1.15),y=sort(bxp$out,decreasing = T)[1:4],labels = sort(resultDF$N,decreasing = T)[1:4],col = "darkgray")

#3. Get correlation data from BQ
querySql <- "SELECT * FROM [isb-cgc-04-0010:LargeResults.miRNA_mRNA_Corr_PerTumorType]"
corrDF <- query_exec(querySql, project=billingProject,max_pages = Inf)
# write.table(corrDF,file = "/Users/varsha/Documents/miRNA/TCGA_miRNA_mRNA_Corr.tsv",sep = "\t")
# corrDF = read.table(file = "/Users/varsha/Documents/miRNA/TCGA_miRNA_mRNA_Corr.tsv",header = T,sep = "\t")
par(mfrow=c(2,1))
plot(density(corrDF$pearson),main = "miRNA - mRNA Pearson Correlation")
abline(v=0,col="blue")
plot(density(corrDF$spearman),main = "miRNA - mRNA Spearman Correlation")
abline(v=0,col="blue")
par(mfrow=c(1,1))
plot(corrDF$pearson,corrDF$spearman,main="Pearson vs Spearman Correlation")


#4. Compute Coefficient of Variation
## For miRNAs
querySql = paste("SELECT l_mirna_id , mRNA_Study, count(*) as N, avg(miRNA_log2_normalized_count) as mean,stddev(miRNA_log2_normalized_count) as sd, stddev(miRNA_log2_normalized_count)/avg(miRNA_log2_normalized_count) as COV",
                 "FROM(",
                 "SELECT mRNA_SampleBarcode , mRNA_Study, l_mirna_id , miRNA_log2_normalized_count", 
                 "FROM [isb-cgc-04-0010:LargeResults.miRNA_mRNA_Exp]", 
                 "group by mRNA_SampleBarcode , mRNA_Study, l_mirna_id , miRNA_log2_normalized_count )",
                 "group by l_mirna_id , mRNA_Study",sep = " ")
miRNA_COV = query_exec(querySql,project = billingProject,max_pages = Inf)
par(mfrow=c(3,1),mar=c(4,5,4,4))
boxplot(miRNA_COV$COV,main = "Distribution of COV\n(miRNA, per Study)", cex.axis=1.2,cex.lab=1.5,ylab="Coefficient of Variation", xlab=paste("Non NA count = ",nrow(miRNA_COV[!is.na(miRNA_COV$COV),])))
plot(miRNA_COV$mean, miRNA_COV$COV, xlab="Mean log2 normalized count",ylab="COV",cex.lab=1.5,cex.axis=1.2)
lines(lowess(miRNA_COV$mean, miRNA_COV$COV),col="blue",lwd=3)
plot(miRNA_COV$mean,miRNA_COV$sd, cex.lab=1.5,cex.axis=1.2, xlab="Mean log2 normalized count",ylab="SD log2 normalized count")
lines(lowess(miRNA_COV$mean, miRNA_COV$sd),col="blue",lwd=3)

## For mRNAs
querySql = paste("SELECT Target_Gene ,mRNA_Study, count(*) as N, avg(mRNA_log2_normalized_count) as mean,stddev(mRNA_log2_normalized_count) as sd,stddev(mRNA_log2_normalized_count)/avg(mRNA_log2_normalized_count) as COV",
                 "FROM(",
                 "SELECT mRNA_SampleBarcode , mRNA_Study,Target_Gene , mRNA_log2_normalized_count", 
                 "FROM [isb-cgc-04-0010:LargeResults.miRNA_mRNA_Exp]", 
                 "group by mRNA_SampleBarcode , mRNA_Study,Target_Gene , mRNA_log2_normalized_count )",
                 "group by Target_Gene,mRNA_Study ",sep = " ")
mRNA_COV = query_exec(querySql,project = billingProject,max_pages = Inf)
par(mfrow=c(3,1),mar=c(4,5,4,4))
boxplot(mRNA_COV$COV,main = "Distribution of COV\n(mRNA, per Study)", cex.axis=1.2,cex.lab=1.5,ylab="Coefficient of Variation", xlab=paste("Non NA count = ",nrow(mRNA_COV[!is.na(mRNA_COV$COV),])))
plot(mRNA_COV$mean, mRNA_COV$COV, xlab="Mean log2 normalized count",ylab="COV",cex.axis=1.2,cex.lab=1.5)
lines(lowess(mRNA_COV$mean, mRNA_COV$COV),col="blue",lwd=3)
plot(mRNA_COV$mean,mRNA_COV$sd, cex.axis=1.2,cex.lab=1.5,xlab="Mean log2 normalized count",ylab="SD log2 normalized count")
lines(lowess(mRNA_COV$mean,mRNA_COV$sd),col="blue",lwd=3)

#5. tumor-type specific boxplots for miRNA-mRNA correlations
par(mfrow=c(1,1))
bxp = boxplot(corrDF$spearman~corrDF$mRNA_Study,ylab="Spearman Correlation",main='miRNA-mRNA Spearman Correlation\n(Tumor-type specific)',las=2)
abline(h=0,col="blue")
text(x=c(1:31),y=c(rep(-0.6,31)),labels = corrDF$N,srt=90)

#6. sample counts per miRNA-mRNA-Study group 
bxp = boxplot(corrDF$N,main="Sample count Distribution\nfor correlations")
text(x=c(1.25,1.25,1.25,1.25,1.25),y=bxp$stats,labels = bxp$stats,col = "darkgray")
#text(x=c(0.75,0.75,0.75,0.75),y=sort(bxp$out,decreasing = T)[1:4],labels = rDF$l_mirna_id[order(resultDF$N,decreasing = T)[1:4]],col = "darkgray")
text(x=c(1.15,1.15,1.15,1.15),y=sort(bxp$out,decreasing = T)[1:4],labels = sort(corrDF$N,decreasing = T)[1:4],col = "darkgray")


#7. top correlation scatter plots

expDF = read.table(file = "/Users/varsha/Documents/miRNA/TCGA_miRNA_mRNA_Exp.tsv",header = T,sep = "\t")
corrDF = read.table(file = "/Users/varsha/Documents/miRNA/TCGA_miRNA_mRNA_Corr.tsv",header = T,sep = "\t")

par(mfrow=c(2,2),mar=c(2,2,2,1))
for(i in c(1,7,10,14))
{
  temp = expDF[expDF$mRNA_Study==as.character(corrDF$mRNA_Study[i]) & expDF$l_mirna_id==as.character(corrDF$l_mirna_id[i]) & expDF$Target_Gene==as.character(corrDF$Target_Gene[i]),c('mRNA_log2_normalized_count','miRNA_log2_normalized_count')]
  plot(temp[[1]],temp[[2]],cex.axis=0.8,cex.main=0.8,main=paste(as.character(corrDF$l_mirna_id[i]),"-",as.character(corrDF$Target_Gene[i]),'-',as.character(corrDF$mRNA_Study[i]),'\n(rho=',round(cor(temp[[1]],temp[[2]],method='spearman'),2),')'))
  lines(lowess(temp[[1]],temp[[2]]),col="blue",lwd=2)
  
}

par(mfrow=c(2,2),mar=c(2,2,2,1))
for(i in c(10798,10799,10800,10801))
{
  temp = expDF[expDF$mRNA_Study==as.character(corrDF$mRNA_Study[i]) & expDF$l_mirna_id==as.character(corrDF$l_mirna_id[i]) & expDF$Target_Gene==as.character(corrDF$Target_Gene[i]),c('mRNA_log2_normalized_count','miRNA_log2_normalized_count')]
  plot(temp[[1]],temp[[2]],cex.axis=0.8,cex.main=0.8,main=paste(as.character(corrDF$l_mirna_id[i]),"-",as.character(corrDF$Target_Gene[i]),'-',as.character(corrDF$mRNA_Study[i]),'\n(rho=',round(cor(temp[[1]],temp[[2]],method='spearman'),2),')'))
  lines(lowess(temp[[1]],temp[[2]]),col="red",lwd=2)
  
}