# mirna-processing

####Initial project
miRNA-gene targeting, with validated and predicted sites

####Concept
1. TCGA data for a project. E.g. KIRC. BigQuery 'expressed' RSEM mRNAseq gene-level data, and miRNA mature strand (miR) data, using thresholds of at least X abundance in at least Y libraries (data sets), potentially considering both primary tumours and adjacent tissue normals. (Typically the RSEM and miR expression data are read into R from .txt or .tgz files.)
2. Calculate Spearman correlations between miRs and genes, using the MatrixEQTL R package, and setting p-value = 0.05 for the run.
3. Threshold the results at FDR < 0.05.
4. Filter thresholded miR-gene pairs using miRTarBase v6.0 records with 'stronger' and 'weaker' experimental evidence. 
5. Filter thresholded miR-gene pairs with TargetScan v7.1 predicted site records. As a first step, simplify the TS results by using gene symbols rather than transcript IDs. Retain site type, and weighted scores, and weighted score percentiles. Note that some 3' UTRs will have multiple binding sites for a miR, and we should add the scores, then calculate percentiles. But we could just start with taking the highest scoring site for a miR and 3' UTR. 

####Resources
**MirTarBase v6.0: one file**  
<http://mirtarbase.mbc.nctu.edu.tw/>  
<http://mirtarbase.mbc.nctu.edu.tw/php/download.php>  
We want: hsa_MTI.xlsx

**TargetScan 7.1 predicted sites: 2 files**  
<http://www.targetscan.org/cgi-bin/targetscan/data_download.vert71.cgi>  
We want two files. For filtering miR-gene anticorrelations, we'll use both:  
Conserved site context++ scores - (17.2 MB, 1.5 M rows)  
Nonconserved site context++ scores - (530.65 MB, 38.3 M rows)  
Note that we will not need all columns from these files. 

We probably want the following general files too, and they are not large, so should be easy to work with. 
miR Family - (0.17 MB)
Gene info - (0.88 MB)
UTR genome coordinates - (0.58 MB)
3P-seq tag info - (8.57 MB)

**MatrixEQTL**  
<https://CRAN.R-project.org/package=MatrixEQTL>  
Docs: <http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/>
