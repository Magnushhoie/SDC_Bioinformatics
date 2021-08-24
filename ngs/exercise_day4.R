# Load the necessary libaries. If you do not have the packages then install the packages and then load the 
# required package.
# To install the DESeq2 package which is a bioconductor package, uncomment the following 3 lines and run the commands.
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2);
library(RColorBrewer);
library(gplots);
library(ggplot2);

description = read.table ("convertion_table.tab",  quote ="\"", sep="\t", header = T);
sample_table = data.frame(sampleName = description$description, 
                          filename= paste(description$files,".counts.txt", sep=""), 
                          condition = description$group);

##load the count data into the "dds" object
dds = DESeqDataSetFromHTSeqCount(sample_table, design = ~condition)
# Relevel the factor variable "dds$condition" to define the normal condition as the base level for
# for the statistical testing in our GLM model. So that later when we  interpret the results of fold changes we know
# what is the reference level
dds$condition <- relevel(dds$condition, ref="normal");


#Now we have loaded the count data into the "dds" object, you should see a small description of the dds object.
#dds

# You can have a look at the table with the counts as well. In rows we have the genes and in columns 
#the sample names
head(counts(dds))


# We are interested in the gene "prostate cancer associated 3" - PCA3. Because the genes in the count dataset
# are annotated using Ensembl ID we need to find the corresponding ID for PCA# gene.
# Let’s try to see the first 10 genes that we have information on:
head(rownames(dds))


# This kind of id are gene ids from Ensembl. To figure out what the id is for PCA3 you can visit the
# Ensembl website (https://www.ensembl.org/index.html) and search for PCA3. Get the ENSG identifier and 
# let’s save this (replace the id below with the id you find at Ensembl):

ensg = "ENSG00000225937"; # could you find the id?


#Solution : ENSG00000225937


# The following commands will save the plot in the file "PCA3_plotcounts.jpeg"
jpeg("PCA3.jpeg");
plotCounts(dds, gene=ensg, main="PCA3");
dev.off();
# Have a quick look in the literature. Are the counts as expected from your knowledge about PCA3?

############ Outlier detection ############
#  To check the distribution of the dataset and investigate possible outliers, use the commands below 
# to do principal component analysis and hierarchical clustering and on the data. We are using the function 
# "rlog" to convert the data in to log2 as well as stabilizing the variance. Try to get the help rlog by 
# writing "?rlog" - read the description of the function.

# principal component analysis
dds_rlog = rlog(dds)	;
dists = as.matrix(dist(t(assay(dds_rlog))));
rownames(dists) = colnames(dists) = colData(dds)$condition;
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100);

jpeg("PCAplot.jpeg");
z <- plotPCA(dds_rlog, intgroup=c("condition"));
nudge <- position_nudge(y=2);
z + geom_label(aes(label = sample_table$sampleName), position = nudge)
dev.off();
# Hierarchical clustering
jpeg("HM.jpeg");
heatmap.2(dists, trace='none', col=rev(hmcol))
dev.off();

#Question 5 : Do you see any obvious trends in the PCA plot?
#Question 6: When looking at the hierarchical clustering do you see any clusters? Do you see any possible outliers
# that would need to be excluded?
  
############  Differential expression ############

# Let us try to identify differentially expressed genes in the data. To do this we fit a negative binomial
# GLM where we compare the groups specified by the condition factor (carcinoma vs. normal) 
# (see object "sample_table" if in doubt).


dds = DESeq2::DESeq(dds);

# The command returns a fitted model object. To get hold of the differential expression test results
# from this object, you should see a brief view of what is inside of the results ("res").

res = results(dds);
write.table(as.data.frame(res),file="res.csv",sep=",")

# Try to understand what is the content of the different columns of the results dataframe. 
# What is a p-value and what is an adjusted p-value?
# Question 7: What is the log2FoldChange column?
  
# Let us order the results after the smallest adjusted p-value. You can see the most and the least 
# significant genes.

resOrdered <- res[order(res$padj),]
resOrdered

# Let us see how many genes are differentially expressed at a false discovery rate of 5%

subset(resOrdered, padj < 0.05)

# Question 8: How many genes are differentially expressed at a False Discovery Rate (FDR) of 5%?

# For you information try to search for those genes in the Ensembl web-site. Which genes are they?
  
resOrdered[ensg,]

# Question 9: What is the direction and the dimension of the foldchange and what is the adjusted pvalue ?

  
  


