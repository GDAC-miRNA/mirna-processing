#stats on differences between miRNA expression calculated by BCGSC and that calculated by ISB using BQ tables and queries

df = read.table('/Users/varsha/Documents/miRNA/Diff_miRNABCGSC_miRNAGroupedByMIMAT.tsv',sep=',',header=T)

#plot(density(df$diff))
#hist(df$diff,breaks = 100)

#differences per miRNA
#boxplot(df$diff~df$mirna_accession)
#boxplot(df$diff~df$AliquotBarcode)
#qqplot(df$normalized_count,df$sumRPMM)
sampleDF = df[sample(c(1:nrow(df)),100000),]
plot(sampleDF$normalized_count,sampleDF$diff, main='Normalized count vs. Difference',xlab = 'Normalized count',ylab = 'Norm. count - Sum RPMM')
lines(lowess(sampleDF$normalized_count,sampleDF$diff),col="blue",lwd=2)

#plot diff as percentage of magnitude of expression
plot(sampleDF$normalized_count,sampleDF$diff/sampleDF$normalized_count, main='Difference as percent normalized count',xlab='Normalized count',ylab='Differrence/Normalized Count')

#plot diff. against mean expression per miRNA
subsetDF = df[,c('mirna_accession','normalized_count','diff')]
meanPerMIRNA = aggregate(subsetDF,list(subsetDF$mirna_accession),mean)
top5 = order(meanPerMIRNA$normalized_count,decreasing = T)[1:5]
plot(meanPerMIRNA$normalized_count,meanPerMIRNA$diff,col='blue',main='Mean per miRNA',xlab='Mean Normalized count',ylab='Mean Difference')
text(x=meanPerMIRNA$normalized_count[top5],y = meanPerMIRNA$diff[top5],labels = meanPerMIRNA$Group.1[top5],pos = 2,cex=0.7)


#head(meanPerMIRNA[order(meanPerMIRNA$normalized_count,decreasing = T),])
