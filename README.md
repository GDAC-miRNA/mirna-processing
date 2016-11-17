# mirna-processing

####Initial project
miRNA-gene targeting, with validated (miRTarBase v6.0) and predicted sites (using only TargetScan v7.1 to start with). 

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
In filtering miR-gene correlations, we'll distinguish (see 'Experiments' and 'Support Type' columns) the preferred but limited 'Functional MTI' records (e.g. Luciferase reporter assay, qRT-PCR, Western blot) from the 'Functional MTI (Weak)' records (e.g. microarray, proteomics, HITS-CLIP, PAR-CLIP, ...). The stronger-evidence records are constrained and biased by which labs have done functional validation work. E.g. miR-21-5p has many 'strong' records, but miR-509-3p has none to few. We use high-scoring predicted targeting (below) to compensate for the limits of the validated records.   
Other validated data sources are miRecords and Tarbase. miRTarBase is the single best resource, if you want to use just one resource.  
miRTarBase is not necessarily up-to-date with PubMed. So for careful work we'd also check PubMed for interesting miR-gene correlations. 

**TargetScan 7.1 predicted sites: 2 files**  
<http://www.targetscan.org/cgi-bin/targetscan/data_download.vert71.cgi>  
We want two files. For filtering miR-gene anticorrelations, we'll use both:  
Conserved site context++ scores - (17.2 MB, 1.5 M rows)  
Nonconserved site context++ scores - (530.65 MB, 38.3 M rows)  
Note that we will not need all columns from these files, but the simplest thing is to load all columns. The two files can be loaded either into separate tables, or into a single table (adding a 'conserved/nonconserved' column). Correlation filtering will typically ignore which table a miR is from.  
Note that TargetScan distinguishes conserved miRNAs from conserved binging (target) sites. Given the effectiveness scores for sites, we ignore whether a site is conserved. (This is surprising, but look at a scatterplot of score vs site conservation.)   
The TargetScan **FAQs** are worth reading rather carefully (<http://www.targetscan.org/faqs.Release_7.html>). 

We probably want the following **general** TargetScan files too. They are small files, so should be easy to work with.  
miR Family - (0.17 MB, if the seed sequence matters most for targeting, here are sets of miRs that have the same seed sequence)  
Gene info - (0.88 MB)  
UTR genome coordinates - (0.58 MB)  
3P-seq tag info - (8.57 MB)  

**MatrixEQTL**  
<https://CRAN.R-project.org/package=MatrixEQTL>  
Docs: <http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/>   
We use this for all correlations. It is easy to use and very very fast. 

