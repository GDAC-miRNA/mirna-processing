#miRNA Biogenesis genes
#Exploring correlations between biogenesis gene alterations and miRNA expression levels

library(bigrquery)
library(ggplot2) #for bar plot
library(gridExtra) #for multiple ggplots on same figure
library(gplots) #for heatmap.2
library(reshape2)
library(effsize) #for cohen.d
billingProject = 'isb-cgc-04-0010' 

#get Biogenesis gene lengths and roles
querySql = "SELECT * FROM `isb-cgc-04-0010.reference_data.miRNABiogenesisGenes_Nykter`"
bgGenes = query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

#####################################
#read biogenesis CNVR data
querySql <- "SELECT * FROM `isb-cgc-04-0010.miRNABioGenesis.PanCan_miRNABioGenesis_CNVR`"
bgCNVR <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
#compute variation in CNVR per gene
perGeneCNVRVar = aggregate(bgCNVR$GISTIC_Calls,by=list(Gene=bgCNVR$Gene_Symbol),var)
plot(density(perGeneCNVRVar$x),main="GISTIC variance per gene")
#visualize cnvr distributions per gene
bgCNVR$Gene_Symbol<-factor(bgCNVR$Gene_Symbol, levels=perGeneCNVRVar$Gene[order(perGeneCNVRVar$x,decreasing = T)])
par(mfrow=c(1,1),mar=c(3,2,2,1))
boxplot(bgCNVR$GISTIC_Calls~bgCNVR$Gene_Symbol,las=2,main="GISTIC estimates\n(bxp sorted by variance)",cex.main=0.9,cex.axis=0.8)

#######################################
#read biogenesis GEXP data
querySql <- "SELECT * FROM `isb-cgc-04-0010.miRNABioGenesis.PanCan_miRNABioGenesis_GEXP`"
bgGEXP <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)
#compute variation in GEXP per gene
perGeneGEXPVar = aggregate(log2(bgGEXP$normalized_count+1),by=list(Gene=bgGEXP$Symbol),var)
plot(density(perGeneGEXPVar$x),main="GEXP variance per gene")
#visualize gexp distributions per gene
bgGEXP$Symbol<-factor(bgGEXP$Symbol, levels=perGeneGEXPVar$Gene[order(perGeneGEXPVar$x,decreasing = T)])
par(mfrow=c(1,1),mar=c(3,2,2,1),mgp=c(3,0.5,0))
boxplot(log2(bgGEXP$normalized_count+1)~bgGEXP$Symbol,las=2,main="Log Normalized Gene Expression\n(bxp sorted by variance)",cex.main=0.9,cex.axis=0.6)

###################################################
#cbind GEXP and GISTIC data to compute ratio GEXP/GISTIC per sample per gene
#make colnames same to simplify merging
colnames(bgCNVR)[colnames(bgCNVR)=='Gene_Symbol'] = 'Symbol'
bg_GEXP_CNVR = merge(bgGEXP,bgCNVR,by = c('SampleBarcode','Symbol'),incomparables = )

bg_GEXP_CNVR['Ratio_GEXP_CNVR'] = log2(bg_GEXP_CNVR$normalized_count+1)/bg_GEXP_CNVR$GISTIC_Calls
boxplot(bg_GEXP_CNVR$Ratio_GEXP_CNVR~bg_GEXP_CNVR$Symbol,las=2,cex.axis=0.6,main="Ratio GEXP/CNVR",yaxt="n")
axis(2, mgp=c(3, .5, 0))

#compute corr between GEXP and CNVR per gene
byCor = by(bg_GEXP_CNVR[,c('normalized_count','GISTIC_Calls')], bg_GEXP_CNVR$Symbol, function(x) {cor(log2(x$normalized_count+1), x$GISTIC_Calls,method = "spearman")})
cor_GEXP_CNVR <- as.data.frame(as.matrix(byCor))
plot(density(cor_GEXP_CNVR$V1,na.rm = T),main='Correlation between GEXP and CNVR',xlab='Spearman Correlation Coefficient',ylab="density")
#################################################

#read coding MAF data
querySql <- "SELECT * FROM `isb-cgc-04-0010.miRNABioGenesis.PanCan_miRNABioGenesis_CodingMAF`"
bgMAF <- query_exec(querySql, project=billingProject,useLegacySql = FALSE,max_pages = Inf)

perGeneMAFCount = table(bgMAF$Hugo_Symbol)
par(mfrow=c(1,1),mar=c(3,2,2,1),mgp=c(3,0.5,0))
barplot(sort(perGeneMAFCount,decreasing = T),las=2,cex.names=0.7,main="No. of coding mutations per gene")

#cbind MAF info with GEXP and CNVR
bg_GEXP_CNVR_MAF = bg_GEXP_CNVR
bg_GEXP_CNVR_MAF['Variant_Classification']=NA
lostMAF = data.frame()
for(index in 1:nrow(bgMAF))
{
print(index)
bg_GEXP_CNVR_MAF$Variant_Classification[bg_GEXP_CNVR_MAF$SampleBarcode==bgMAF$Tumor_SampleBarcode[index] & as.vector(bg_GEXP_CNVR_MAF$Symbol)==bgMAF$Hugo_Symbol[index]] = bgMAF$Variant_Classification[index]
if(length(bg_GEXP_CNVR_MAF$Variant_Classification[bg_GEXP_CNVR_MAF$SampleBarcode==bgMAF$Tumor_SampleBarcode[index] & as.vector(bg_GEXP_CNVR_MAF$Symbol)==bgMAF$Hugo_Symbol[index]])==0)
{
  lostMAF[nrow(lostMAF)+1,] = bgMAF[index,]
  #write.csv(bgMAF[index,],file = "/Users/varsha/Documents/miRNA/miRNABiogenesis/MAFLost.csv",append = T,sep = "\t")
}
}

##############################
#binarize CNVR||MAF 
bg_GEXP_CNVR_MAF['IsAltered'] = bg_GEXP_CNVR_MAF$GISTIC_Calls>1.5 | bg_GEXP_CNVR_MAF$GISTIC_Calls< (-0.5) | !is.na(bg_GEXP_CNVR_MAF$Variant_Classification)
#write this df to file so downstream steps dont have to rerun above steps again 
write.table(bg_GEXP_CNVR_MAF,file="/Users/varsha/Documents/miRNA/miRNABiogenesis/bg_GEXP_CNVR_MAF.tsv",sep = "\t")
###############################
####FOLLOWING STEPS READ FROM A SAVED DF, SO NO NEED TO RUN PRECEDING STEPS AGAIN.
##Functional alterations - association between BG gene alteration and BG gene expression levels
#read from file
bg_GEXP_CNVR_MAF = read.table("/Users/varsha/Documents/miRNA/miRNABiogenesis/bg_GEXP_CNVR_MAF.tsv",sep="\t",header = T)
bg_FuncAlt = bg_GEXP_CNVR_MAF[,c('SampleBarcode','Study.x','Symbol','normalized_count','GISTIC_Calls','IsAltered')]
bg_FuncAlt[,'log2_count'] = log2(bg_FuncAlt$normalized_count+1)
bg_FuncAlt[,'IsMutated'] = !is.na(bg_GEXP_CNVR_MAF$Variant_Classification)
bg_FuncAlt[bg_FuncAlt$GISTIC_Calls>1.5,'IsCopyNumAltered'] = "AMP"
bg_FuncAlt[bg_FuncAlt$GISTIC_Calls<(-0.5),'IsCopyNumAltered'] = "DEL"
perGeneStat = by(data = bg_FuncAlt,bg_FuncAlt$Symbol,function(x) kruskal.test(log2_count~IsAltered,data=x))
perGeneStatDF = data.frame(do.call("rbind", perGeneStat))
perGeneStatDF = perGeneStatDF[order(unlist(perGeneStatDF$p.value),decreasing = F),]

#plot relationships between cnvr and gexp, AND between mutation and gexp, AND between binarized alteration status and gexp
thisGene = "SMAD4"
thisGene_bg_FuncAlt = bg_FuncAlt[bg_FuncAlt$Symbol==thisGene,]
spearmanCorr = cor.test(x = thisGene_bg_FuncAlt$log2_count,y = thisGene_bg_FuncAlt$GISTIC_Calls,method = "spearman")
kruskal_GEXP_BinaryCN = kruskal.test(x=thisGene_bg_FuncAlt$log2_count,g = as.factor(thisGene_bg_FuncAlt$IsCopyNumAltered))
binaryCN_effsize = cohen.d(thisGene_bg_FuncAlt$log2_count,thisGene_bg_FuncAlt$IsCopyNumAltered,pooled = T,paired = F,na.rm = T,hedges.correction = T)
kruskal_GEXP_BinaryMut = kruskal.test(x=thisGene_bg_FuncAlt$log2_count,g = as.factor(thisGene_bg_FuncAlt$IsMutated))
binaryMut_effsize = cohen.d(thisGene_bg_FuncAlt$log2_count,thisGene_bg_FuncAlt$IsMutated,pooled = T,paired = F,na.rm = T,hedges.correction = T)
kruskal_GEXP_BinaryAlt = kruskal.test(x=thisGene_bg_FuncAlt$log2_count,g = as.factor(thisGene_bg_FuncAlt$IsAltered))
binaryAlt_effsize = cohen.d(thisGene_bg_FuncAlt$log2_count,thisGene_bg_FuncAlt$IsAltered,pooled = T,paired = F,na.rm = T,hedges.correction = T)


give.n <- function(x){
  return(c(y = min(x)-1, label = length(x)))
}
#cnvr vs. gexp scatter plot
cnvr_gexp_plot = ggplot(thisGene_bg_FuncAlt, aes(x=GISTIC_Calls, y=log2_count)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method = loess) +           # Add a loess smoothed fit curve with confidence region
  ylim(0,max(thisGene_bg_FuncAlt$log2_count)) + 
  annotate("text",x=0,y=1,label=paste("rho = ",round(spearmanCorr$estimate,2),sep=""))

binaryCNVR_gexp_plot = ggplot(thisGene_bg_FuncAlt, aes(IsCopyNumAltered,log2_count,fill=IsCopyNumAltered)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  guides(fill=guide_legend(title="")) + 
  ylim(0,max(thisGene_bg_FuncAlt$log2_count)) + 
  annotate("text",x=2,y=1,label=paste("Kruskal p-value = ",signif(kruskal_GEXP_BinaryCN$p.value,digits = 3),sep="")) +
  annotate("text",x=2,y=2,label=paste("Effect size = ",round(binaryCN_effsize$estimate,2)," (",binaryCN_effsize$magnitude,")",sep=""))


binaryMut_gexp_plot = ggplot(thisGene_bg_FuncAlt, aes(IsMutated,log2_count,fill=IsMutated)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  guides(fill=guide_legend(title="")) +
  ylim(0,max(thisGene_bg_FuncAlt$log2_count)) + 
  annotate("text",x=1.5,y=1,label=paste("Kruskal p-value = ",signif(kruskal_GEXP_BinaryMut$p.value,digits = 3),sep="")) +
  annotate("text",x=1.5,y=2,label=paste("Effect size = ",round(binaryMut_effsize$estimate,2)," (",binaryMut_effsize$magnitude,")",sep=""))



binaryAlt_gexp_plot = ggplot(thisGene_bg_FuncAlt, aes(IsAltered,log2_count,fill=IsAltered)) + 
  geom_boxplot() + 
  stat_summary(fun.data = give.n, geom = "text",size=2.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom") +
  guides(fill=guide_legend(title="")) +
  ylim(0,max(thisGene_bg_FuncAlt$log2_count)) + 
  annotate("text",x=1.5,y=1,label=paste("Kruskal p-value = ",signif(kruskal_GEXP_BinaryAlt$p.value,digits = 3),sep="")) +
  annotate("text",x=1.5,y=2,label=paste("Effect size = ",round(binaryAlt_effsize$estimate,2)," (",binaryAlt_effsize$magnitude,")",sep=""))

contTable = table(thisGene_bg_FuncAlt$IsCopyNumAltered,thisGene_bg_FuncAlt$IsMutated,useNA = "ifany",dnn = c('CopyNumberAltered','Mutated'))

grid.arrange(cnvr_gexp_plot,binaryCNVR_gexp_plot,binaryMut_gexp_plot,binaryAlt_gexp_plot,ncol=2,top=thisGene)
######################################
#CONTINUE HERE : draw heatmaps showing percentage of samples AMPlified, DELeted, and MUTates per BG gene per tumor type

#get total number of samples per study
uniqSampleStudyList = unique(bg_FuncAlt[,c('SampleBarcode','Study.x')])
countSamplePerStudy = as.data.frame(table(uniqSampleStudyList$Study.x))
colnames(countSamplePerStudy) = c('Study','TotalN')

#AMPlifications
ampData = bg_FuncAlt[!is.na(bg_FuncAlt$IsCopyNumAltered) & bg_FuncAlt$IsCopyNumAltered=="AMP",]
countAMPPerGeneStudy = as.matrix(table(ampData$Study.x,ampData$Symbol))
countAMPPerGeneStudy = cbind(countAMPPerGeneStudy,countSamplePerStudy$TotalN)
colnames(countAMPPerGeneStudy)[ncol(countAMPPerGeneStudy)] = 'TotalN'
percAMPPerGeneStudy = countAMPPerGeneStudy*100/countAMPPerGeneStudy[,'TotalN']
#remove TotalN column
percAMPPerGeneStudy = percAMPPerGeneStudy[,1:ncol(percAMPPerGeneStudy)-1]

#DELetions
delData = bg_FuncAlt[bg_FuncAlt$IsMutated==T,]
countDELPerGeneStudy = as.matrix(table(delData$Study.x,delData$Symbol))
countDELPerGeneStudy = cbind(countDELPerGeneStudy,countSamplePerStudy$TotalN)
colnames(countDELPerGeneStudy)[ncol(countDELPerGeneStudy)] = 'TotalN'
percDELPerGeneStudy = countDELPerGeneStudy*100/countDELPerGeneStudy[,'TotalN']
#remove TotalN column
percDELPerGeneStudy = percDELPerGeneStudy[,1:ncol(percDELPerGeneStudy)-1]

#MUTations
mutData = bg_FuncAlt[!is.na(bg_FuncAlt$IsCopyNumAltered) & bg_FuncAlt$IsCopyNumAltered=="DEL",]
countDELPerGeneStudy = as.matrix(table(delData$Study.x,delData$Symbol))
countDELPerGeneStudy = cbind(countDELPerGeneStudy,countSamplePerStudy$TotalN)
colnames(countDELPerGeneStudy)[ncol(countDELPerGeneStudy)] = 'TotalN'
percDELPerGeneStudy = countDELPerGeneStudy*100/countDELPerGeneStudy[,'TotalN']
#remove TotalN column
percDELPerGeneStudy = percDELPerGeneStudy[,1:ncol(percDELPerGeneStudy)-1]

ampData = bg_FuncAlt[!is.na(bg_FuncAlt$IsCopyNumAltered) & bg_FuncAlt$IsCopyNumAltered=="AMP",]
countAMPPerGeneStudy = as.matrix(table(ampData$Study.x,ampData$Symbol))
countAMPPerGeneStudy = cbind(countAMPPerGeneStudy,countSamplePerStudy$TotalN)
colnames(countAMPPerGeneStudy)[ncol(countAMPPerGeneStudy)] = 'TotalN'
percAMPPerGeneStudy = countAMPPerGeneStudy*100/countAMPPerGeneStudy[,'TotalN']
#remove TotalN column
percAMPPerGeneStudy = percAMPPerGeneStudy[,1:ncol(percAMPPerGeneStudy)-1]

#count number of samples altered per gene per study
countByGeneStudy = by(bg_GEXP_CNVR_MAF$'IsAltered',bg_GEXP_CNVR_MAF[,c('Study.x','Symbol')],sum)
countByGeneStudyMat = as.matrix(array(countByGeneStudy,dim(countByGeneStudy),dimnames(countByGeneStudy)))
countByGeneStudyMat= t(countByGeneStudyMat)

#melt into long form for stacked barchart
countByGeneStudyLong = melt(countByGeneStudyMat,id.vars = 'Symbol')#,measure.vars = colnames(countByGeneStudyArray),variable.name = 'Gene')
colnames(countByGeneStudyLong) = c('Symbol','Study','Count')
#aggregate counts per gene to order bars in the barplot below
countByGene = aggregate(countByGeneStudyLong$Count,by = list(Symbol=countByGeneStudyLong$Symbol),sum)
countByGeneStudyLong$Symbol = factor(countByGeneStudyLong$Symbol,levels = countByGene$Symbol[order(countByGene$x,decreasing = T)])
ggplot(countByGeneStudyLong, aes(x = Symbol , y = Count)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),plot.title = element_text(hjust = 0.5)) + 
  labs(x="miRNA Biogenesis Gene",y="Number of samples altered",title="PanCancer")
#per gene per study percentage heatmap


#divide each column of countByGeneStudyMat by TotalN count of Study
countByGeneStudyDF = as.data.frame(countByGeneStudyMat)
countByGeneStudyDF['TotalN',] = countSamplePerStudy$TotalN
countByGeneStudyDF = as.data.frame(t(countByGeneStudyDF))
#remove columns(genes) with all NAs
countByGeneStudyDF = countByGeneStudyDF[,colSums(is.na(countByGeneStudyDF))!=nrow(countByGeneStudyDF)]
countByGenePanCan = colSums(countByGeneStudyDF,na.rm = T)
countByGeneStudyDF = rbind(countByGeneStudyDF,countByGenePanCan)
rownames(countByGeneStudyDF)[nrow(countByGeneStudyDF)] = 'PanCan'
percByGeneStudyDF = countByGeneStudyDF*100/countByGeneStudyDF$TotalN

#remove TotalN column before plotting
percByGeneStudyDF = percByGeneStudyDF[,-c(ncol(percByGeneStudyDF))]

#cluster both genes as well as tumor types
par(mfrow=c(2,2),mar=c(0.5,0.5,0.5,0.5))
#lmat = rbind(c(0,3,4),c(2,1,0))
#lhei=c(1.5,4)
#lwid = c(0.5,4,2)
scaleRed<-colorRampPalette(colors=c("white","red"))(1000)
roleLevels = levels(as.factor(bgGenes$Role))
#reorder bgGenes to match order of columns of percByGeneStudyDF
bgGenes = bgGenes[match(colnames(percByGeneStudyDF),bgGenes$Gene),]

heatmap.2(as.matrix(percByGeneStudyDF),cexRow=1.1,cexCol = 1.1,cex.main=0.6,adjCol = c(1,0.5),trace = "none",main = "miRNA Biogenesis Genes\nPercentage of samples altered",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.5, ColSideColors = c("gray","blue","black")[match(bgGenes$Role,roleLevels)])
legend("topright",legend = c("Cytoplasmic Processing", "Nuclear Processing", "Transport"), col = c("gray", "blue", "black"),lty=1,lwd=5)
#one without row clustering so that PanCan can be shown in the last row of the heatmap
#heatmap.2(as.matrix(percByGeneStudyDF),Rowv=F,Colv=F,dendrogram='none',cexRow=1.1,cexCol = 1.1,adjCol = c(1,0.5), lmat = lmat,lhei = lhei,lwid = lwid,trace = "none",main = "miRNA Biogenesis Genes\nPercentage of samples altered",na.color = "gray",symbreaks = F,col=scaleRed,density.info="none",keysize = 0.4)

########################
##correlation between gene length and number of alterations observed
pancanPerc = as.data.frame(t(percByGeneStudyDF['PanCan',]))
pancanPerc[,'Gene'] = rownames(pancanPerc)
x = merge(bgGenes,pancanPerc,by="Gene")
cor.test(x$Length,x$PanCan,method = "spearman")
plot(x$Length,x$PanCan,xlab = "Gene length",ylab = "% samples altered",col="blue",lwd=2)
lo = loess(x$PanCan~x$Length)
xl <- seq(min(x$Length),max(x$Length), (max(x$Length) - min(x$Length))/1000)
lines(xl, predict(lo,xl), col='red', lwd=2)

mtext("Spearman = 0.34",3)
############################
#Get miRNA expression data (saved locally from BQ isb-cgc-04-0010.draft_new_data.TCGA_gdc_hg19_miRNA_Expression)
mirnaDF = read.table('/Users/varsha/Documents/miRNA/data/TCGA_gdc_hg19_miRNA_Expression_LongFormat.tsv',sep = "\t",header = T)
colnames(mirnaDF) = c('SampleBarcode','miRNA_ID','RPM','Study')
mirnaDFToSpread = mirnaDF[,c('SampleBarcode','miRNA_ID','RPM')]
mirnaDF_Wide = spread(mirnaDFToSpread,miRNA_ID,RPM)

#spread bg_GEXP_CNVR_MAF into a sampleXgene matrix and cbind with mirnaDF_Wide
bg_GEXP_CNVR_MAF_ToSpread = bg_GEXP_CNVR_MAF[,c('SampleBarcode','Symbol','IsAltered')]
#eliminate duplicate rows (coming from duplicates in GEXP data)
bg_GEXP_CNVR_MAF_ToSpread = unique(bg_GEXP_CNVR_MAF_ToSpread)
bg_GEXP_CNVR_MAF_Wide = spread(bg_GEXP_CNVR_MAF_ToSpread,Symbol,IsAltered)

#merge biogenesis alterations with miRNA expression data
miRNA_bgAlterations_DF = merge(mirnaDF_Wide,bg_GEXP_CNVR_MAF_Wide,by = "SampleBarcode")

#create and populate an empty matrix of miRNA X Biogenesis genes to store t-test p-values
miRNA_bg_PValues = matrix(data=NA,nrow = ncol(bg_GEXP_CNVR_MAF_Wide)-1,ncol = ncol(mirnaDF_Wide)-1)
rownames(miRNA_bg_PValues) = colnames(bg_GEXP_CNVR_MAF_Wide)[2:ncol(bg_GEXP_CNVR_MAF_Wide)]
colnames(miRNA_bg_PValues) = colnames(mirnaDF_Wide)[2:ncol(mirnaDF_Wide)]
for(rowIndex in 1:nrow(miRNA_bg_PValues)) #biogenesis genes
{
  print(rowIndex)
  res = apply(miRNA_bgAlterations_DF[2:1871],2,function(x) t.test(x~miRNA_bgAlterations_DF[,rownames(miRNA_bg_PValues)[rowIndex]])$p.value)
  miRNA_bg_PValues[rowIndex,]= res
}
write.table(miRNA_bg_PValues,file = "/Users/varsha/Documents/miRNA/miRNABiogenesis/TCGA_PanCan_miRNA_BiogenesisGenes_TTestPValues.tsv",sep = "\t")
########################

#draw boxplots of p-values per biogenesis genes
miRNA_bg_PValues = as.matrix(read.table("/Users/varsha/Documents/miRNA/miRNABiogenesis/TCGA_PanCan_miRNA_BiogenesisGenesAlterations_TTestPValues.tsv",sep = "\t",header = T))
miRNA_bg_PValues = -log10(miRNA_bg_PValues)
#order rows by median p-value
miRNA_bg_PValues = miRNA_bg_PValues[rownames(miRNA_bg_PValues)[order(rowMeans(miRNA_bg_PValues,na.rm = T),decreasing = T)],]
par(mar=c(4,4,4,4),mgp=c(2,0.5,0))
boxplot(miRNA_bg_PValues,use.cols = F,las=2,ylab="-Log10 p-values",main='miRNA Expression vs. Alteration in biogenesis genes\nT-test p-values',cex.main=0.9,cex.axis=0.95)

#get top 6 associations and visualize underlying data
miRNA_bg_PValues_Long = melt(miRNA_bg_PValues)
#order by -log10 pvalue
miRNA_bg_PValues_Long = miRNA_bg_PValues_Long[order(miRNA_bg_PValues_Long$value,decreasing = T),]
#plot top 6
par(mfrow=c(2,3),mar=c(2,2,3,2))
for(index in 1:6)
{
  boxplot(miRNA_bgAlterations_DF[,as.vector(miRNA_bg_PValues_Long$Var2)[index]]~miRNA_bgAlterations_DF[,as.vector(miRNA_bg_PValues_Long$Var1)[index]],main=paste(as.vector(miRNA_bg_PValues_Long$Var2)[index],"|",as.vector(miRNA_bg_PValues_Long$Var1)[index],"\n-log10(p-val) = ",round(as.vector(miRNA_bg_PValues_Long$value)[index],2),sep=" "),main.cex=0.9)
}


#############################
#association between biogenesis alterations and whole genome miRNA levels per sample
miRNA_WholeGenomePerSample = aggregate(mirnaDF$RPM,by=list(SampleBarcode =mirnaDF$SampleBarcode),sum)
