## RNAseq workshop with bioconductor and DESeq2
## Instructor:  Brendan Jeffrey, Ph.D.
## Email: brendan.jeffrey@nih.gov

## DESeq2 manuscript: https://www.ncbi.nlm.nih.gov/pubmed/25516281

# package installation
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# # install packages with Bioconductor package manager
# BiocManager::install("airway", version = "3.8")
# BiocManager::install("Rsamtools", version = "3.8")
# BiocManager::install("GenomicFeatures", version = "3.8")
# BiocManager::install("DESeq2", version = "3.8")
# BiocManager::install("vsn", version = "3.8")
# BiocManager::install("hexbin", version = "3.8")
# BiocManager::install("pheatmap", version = "3.8")
# BiocManager::install("apeglm", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")
# BiocManager::install("ReportingTools", version = "3.8")
# BiocManager::install("Gviz", version = "3.8")
# BiocManager::install("ggbeeswarm", version = "3.8")
# BiocManager::install("EnrichmentBrowser", version = "3.8")


###############################################################################
# Part 1 generating the count matrix, summarized experiment, sample information
###############################################################################

#### generating summarizedExperiment object METHOD 1, using summarizeOverlaps ####
#### airway study file and sample information ####
# load airway package - see https://bioconductor.org/packages/3.8/data/experiment/html/airway.html
# "airway" manuscript - https://www.ncbi.nlm.nih.gov/pubmed/24926665
library("airway")

# get information about a package
browseVignettes("airway")

# where is this data located, and list the files in that directory
indir <- system.file("extdata", package="airway", mustWork=TRUE)
indir
list.files(indir)

# file with detailed information about each of our samples - sample_table
csvfile <- file.path(indir, "sample_table.csv")
csvfile
sampleTable <- read.csv(csvfile, row.names = 1)

# print out contents
sampleTable

filenames <- file.path(indir, paste0(sampleTable$Run, "_subset.bam"))
sampleTable$Run

file.exists(filenames)
# what else could you have used for generating the bam files?
# filenames <- file.path(indir, paste0(row.names(sampleTable), "_subset.bam"))

library("Rsamtools")
# indicate that these are bam files usiing BamFileList function from Rsamtools package
# only process 2 million reads at a time
bamfiles <- BamFileList(filenames, yieldSize=2000000)
bamfiles

# NOTE: make sure that the chromosome names of the genomic features in the annotation 
# you use are consistent with the chromosome names of the reference used for read alignment
# check chromosome names (here called seqnames)
seqinfo(bamfiles[1])

# Define gene models - read from GTF file http://www.ensembl.org/info/website/upload/gff.html
library("GenomicFeatures")

# genomic coord and gene file for airway study (GTF)
gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")

# load from this gtf file and indicate none of our sequences are circular
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
?makeTxDbFromGFF

# make a GRangesList of all exons grouped by gene
# get more information on a function
?exonsBy
ebg <- exonsBy(txdb, by="gene")
ebg

# Read counting step using summarizeOverlaps function from GenomicAlignments package
# NOTE: can use multiple cores to speed up if desired using BiocParallel package https://bioconductor.org/packages/3.8/BiocParallel

# load GenomicAlignments package
library("GenomicAlignments")

# create the SummarizedExperiment object 
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

# explore options of summarizeOverlaps
?summarizeOverlaps

# other tools that could be used:
# summarizeOverlaps, GenomicAlignments package, 'DESeqDataSet' function
# featureCounts, Rsubread package, 'DESeqDataSetFromMatrix' function
# tximport, tximport package, 'DESeqDataSetFromTximport' function
# htseq-count, HTSeq Python package (not part of R), 'DESeqDataSetFromHTSeq' function

# explore the SummarizedExperiment object
browseVignettes("SummarizedExperiment")
se

# look at assay - the count matrix
assayNames(se)
assay(se)

# rowRanges slot of se, contains information on the genomomic ranges, ie each gene
rowRanges(se)

# colData slot of se (currently empty), contains information about the samples
colData(se)

# populate the colData slot with information from sampleTable
sampleTable

# Because we used a column of  sampleTable to produce the bamfiles vector, 
# we know the columns of se are in the same order as the rows of sampleTable
colData(se) <- DataFrame(sampleTable)
colData(se)

#### generating DESeqDataSet object METHOD 2, pre-existing count matrix ####
# read in counts matrix (HTSeq generated, STAR, others)
counts <- read.table("../data/RNAseq_counts_matrix.txt", header = TRUE, row.names=1)
counts

# read in sample data
sample_table <- read.table("../data/RNAseq_experiment_design.csv", header=TRUE, sep=',')
sample_table

# experimental design
exp_design <- data.frame(row.names = sample_table$replicate, 
                         condition = sample_table$condition,
                         sample_names = sample_table$sample_name)
exp_design

# generate the DESeqDataSet using the DESeqDataSetFromMatrix function in DESeq2
# DESeqDataSet is a custom class within DESeq2, built on top of the generalize summarizedExperiment class
library("DESeq2")
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = exp_design, design = ~ condition)

# explore the DESeqDataSet, and compare to summarizedExperiment above
dds
se
colData(dds)
colData(se)

rowRanges(dds)
rowRanges(se)

#### END generating summarized experiment, DESeqDataSet ####


#### BRANCH point ####
# after generating the summarizedExperiment (DESeqDataSet) can use a variety of Bioconductor packages
# for differential gene expression. These include:
# edgeR https://bioconductor.org/packages/3.8/edgeR
# limma https://bioconductor.org/packages/3.8/limma
# EBSeq https://bioconductor.org/packages/3.8/EBSeq - isoform differential expression

###############################################################################
# Part 2 using DESeq2 for data exploration 
###############################################################################

#### starting from a prepared SummarizedExperiment ####
# similar to SummarizedExperiment above, except using all genes
# clear environment

# load airway summarized experiment data
data("airway")
airway

# rename SummarizedExperiment object
se <- airway

# explore
rowRanges(se)
colData(se)

# we want to specify that untrt (untreated) is the reference level for dex (dexamethasone)
se$dex
se$dex <- relevel(se$dex, "untrt")
se$dex

# check the millions of fragments that uniquely aligned to genes, one decimal
round( colSums(assay(se)) / 1e6, 1 )

# construct the DESeqDataSet object from the SummarizedExperiment object
# DESeqDataSet is a custom class within DESeq2, built on top of the generalized summarizedExperiment object
dds <- DESeqDataSet(se, design = ~ cell + dex)

# The simplest design formula for differential expression would be ~ condition, 
# where  condition is a column in colData(dds) that specifies which of two (or more groups) 
# the samples belong to. For the airway experiment, we will specify ~ cell + dex 
# meaning that we want to test for the effect of dexamethasone (dex) 
# controlling for the effect of different cell line (cell)
# also used when there are batch effects present in data

# pre-filter the counts data to remove those genes with zero counts
nrow(dds)

# slice of the dds values [row, column]
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

#### Exploratory analysis and visualization ####
# variance stabilization transformation and rlog
# Common statistical methods for exploratory analysis, eg clustering and principal component analysis (PCA)
# work best with data that has the same range of variance at different ranges of means (homoskedastic)
# For RNAseq data, expected variance grows with the mean
# performing PCA on raw counts, resulting plot depends mostly on genes with highest counts
# performing PCA on simple log transformed counts plus pseudo count of 1, the genes with low counts
# will contribute a lot of noise, due to taking logarithm of small counts inflates their variance

## Effects of transformations on the variance ##
# vizualizing with simulated data
library("vsn")
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)

# plot on non-transformed
meanSdPlot(cts, ranks = FALSE)

# plot on log transformed count + 1
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

# DESeq2 offers two transformation methods, VST and RLD
# VST (variance stabilizing transformation) faster, for medium to large datasets
# RLD (Regularized-logarithm transformation), better for small datasets, n < 30

# VST - transform our data
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
meanSdPlot(assay(vsd), ranks = FALSE)

# compare to raw read count
meanSdPlot(assay(dds), ranks = FALSE)

#### Sample distances ####
# A useful first step in an RNA-Seq analysis is often to assess 
# overall similarity between samples: Which samples are similarto each other, 
# which are different? Does this fit to the expectation from the experiment’s design?
# R function dist to calculate the Euclidian distance between samples
# t used to transpose the results of assay(rld), samples are now rows
sampleDists <- dist(t (assay(vsd)) )
sampleDists

# vizualize the distance matrix using pheatmap package
library("pheatmap")
library("RColorBrewer")

# sample distances to matrix, base R as.matrix function
sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix

# plot, , specify rownames and colnames for clarity
sampleDistMatrix
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
sampleDistMatrix

# default colors
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

# setting colors using RcolorBrewer
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


## visualize sample to sample distances using Principal component analysis, DESeq2 plotPCA function ##
colData(vsd)
plotPCA(vsd, intgroup = "dex")
plotPCA(vsd, intgroup = c("dex", "cell"))

# quick aside, saving figures
library("ggplot2")

plt <- plotPCA(vsd, intgroup = c("dex", "cell"))
plt
ggsave("../results/default_PCA.pdf", height=4, width=6.5)

# more control over plotting with ggplot
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

# build plot
gplt <- ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
gplt
ggsave("../results/ggplot_PCA.pdf", height=4, width=6.5)

# ggplot cheat cheet - https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf


## Heatmap of gene variance estimates ##
# extract row indexes of those genes that have largest variance estimates
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
topVarGenes

# subset of the vsd values [row, column]
mat <- assay(vsd)[ topVarGenes, ]
mat

# how does each genes transformed count deviate from the genes average across all samples
mat <- mat - rowMeans(mat)
mat
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
anno
pheatmap(mat, annotation_col = anno)
?pheatmap

###############################################################################
# Part 3 using DESeq2 for differential expression analysis
###############################################################################

## DESeq2 differential expression with DESeq function ##
dds <- DESeq(dds)
?DESeq

# rowData of dds now contains all fitted parameters
head(rowData(dds))

# extract results from dds DESeq2 results function
res <- results(dds)
res

# generating results with more specific command, more than one treatment, forcing denominator
res <- results(dds, contrast=c("dex","trt","untrt"))

# DESeq2 performs for each gene a hypothesis test to see whether evidence 
# is sufficient to decide against the null hypothesis that there is zero effect 
# of the treatment on the gene and that the observed difference between 
# treatment and control was merely caused by experimental variability

# summary of the results
summary(res)

# DESeq2 uses the Benjamini-Hochberg (BH) adjustment for FDR
# adjusting thresholds, false discovery rate threshold
res.05 <- results(dds, alpha = 0.05)
summary(res.05)

# subset results, consider a fraction of 10% false positives acceptable FDR alpha=0.1
resSig <- subset(res, padj < 0.1)
summary(resSig)

# extract up regulated genes
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigUp

# extract down regulated
resSigDown <- subset(resSig, log2FoldChange < 0)
resSigDown

## plotting results ##
# counts plot of individual genes
library("ggbeeswarm")
topGene <- rownames(resSig)[which.min(resSig$padj)]
topGene
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"), returnData = TRUE)

# generate plot
countplt <- ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + 
  geom_point(size = 3) + 
  geom_line() + 
  ggtitle(topGene)
countplt

####### MA plots? M (log ratio) and A (mean average)  ####### 
# overview of the distribution of the estimated coefficients, or comparisons of interest, across all genes
library("apeglm")

# show model coefficients
resultsNames(dds)

# moderate or shrink log2 fold changes
resMA <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(resMA, ylim = c(-5, 5))

# compare with no statistical moderation to shrink noisy log2 fold change
plotMA(res, ylim = c(-5, 5))

################################
### everyone loves a heatmap ###
# heatmap of top 50 significant genes
topSigGenes <- head(rownames(resSig[order(resSig$padj),]), 50)
topSigGenes

# extract VST transformed value for all genes for each sample
assay(vsd)
matSig <- assay(vsd)[topSigGenes, ]
matSig

# how does each genes transformed count deviate from the genes average across all samples
matSig <- matSig - rowMeans(matSig)
matSig

# annotation from vsd, remember [row, column]
anno <- as.data.frame(colData(vsd)[ , c("cell","dex")])
anno

# generate heatmap
pheatmap(matSig, annotation_col = anno)

#### Annotating and exporting results ####
library("AnnotationDbi")
library("org.Hs.eg.db")

# annotation keys available
columns(org.Hs.eg.db)

rowRanges(dds)$symbol <- mapIds(org.Hs.eg.db,
                            keys=row.names(dds),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")



# add gene symbol
resSig$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(resSig),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

# add entrez ID
resSig$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(resSig),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# google search on entrez gene to check
resSig

# export csv
resOrdered <- resSig[order(resSig$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "../results/results.csv")

## nice HTML reporting
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="../results/report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)


### additional plotting - plotting changes in genomic space ###
resGR <- results(dds, format="GRanges")
resGR$log2FoldChange <- res$log2FoldChange
resGR

# add annotation again
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")
library("Gviz")

# 1 million base pairs up and downstream from gene with smallest pvalue
topGene <- which.min(resGR$padj)

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]

# if gene has no symbol or is duplicated
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
naOrDup

# use the gene name if no symbol, else use symbol name
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
resGRsub

# create a vector specifying if the genes in this window had a low padj value
status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj), "sig", "notsig"))
status

# plot using Gviz
options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group", notsig = "grey", sig = "hotpink")

