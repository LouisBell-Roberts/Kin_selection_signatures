########Differential gene expression analysis######
###First, imports transcript-level estimates from Salmon and optionally summarizes abundances, counts, and transcript lengths to the gene-level (default) or outputs transcript-level matrices (see txOut argument)
###Second performs differential gene expression analysis using DESeq2
#Louis Bell-Roberts
#22/02/2023

library(data.table)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

#The directory where my transcript data has been quantified using Salmon which is located on my SSD hard drive
setwd("/Volumes/Pop_Gen/Differential_gene_expression/Species/Monomorium_pharaonis/quants/")

# Obtain a vector of all filenames including the path
files <- list.files(recursive = T, pattern = "quant.sf")

##Download the .gff file for that species and place it inside the "quants" directory. This file contains meta-data on the genes for that species (gene annotation)
#Get the names of the genes
gff <- list.files(pattern = ".gff")

#Make TxDb file from the GFF file
txdb <- makeTxDbFromGFF(file = gff)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi.salmon.g <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi.salmon.t <- tximport(files, type = "salmon", txOut = T)

# Gene-level: I think this is the one that we want
head(txi.salmon.g$counts)[,1:7]

# Transcript-level
head(txi.salmon.t$counts)[,1:7]


############ ############ ############ ############ ############

#########Differential gene expression analysis ############

############ ############ ############ ############ ############

#Select columns 1:7 - why?
countData = txi.salmon.g$counts[,1:7]

#Import Salmon NumReads values with tximport which contains the transcript-level counts (txIn = TRUE), abundances and average transcript lengths, then output gene-level summarization (txOut = FALSE)
txi = tximport(files, type = "salmon", tx2gene = tx2gene, txIn = TRUE, txOut = FALSE, countsFromAbundance = "no")
class(txi)
names(txi)

#Create dataframe with two columns: one for the individual's sample number and the caste that the sample comes from
id = c("DRR032045", "DRR032148", "DRR032150", "DRR032224", "DRR032226", "DRR032265", "DRR032266")
dex = c("queen", "worker", "worker", "worker", "worker", "queen", "queen")
colData = data.frame("id"=id, "dex"=dex)
colnames(countData) <- colData$id

#Create DESeqDataSet object
dds = DESeqDataSetFromMatrix(
  countData=round(countData), 
  colData = colData,
  design=~dex)

#Run the differential expression analysis
dds = DESeq(dds)

#Summarise the results
##We set the adjusted p-value cutoff (FDR) to be 0.05, hence we change the default significance cutoff used for optimizing the independent filtering alpha from 0.1 to 0.05
cbind(resultsNames(dds))
res <- results(dds, name = "dex_worker_vs_queen", alpha = 0.05)
summary(res)

#The results res object contains the follow columns
mcols(res)$description
head(res)
## Log2fold change e.g. log2(8) = 3. So if log2fold change is 3, this means it is 2^3 times greater

#Making volcano plots
##Enhanced volcano plot
###Directly from the tutorial
pCutoff = 0.05
FCcutoff = 1.0
# 
# p = EnhancedVolcano(data.frame(res), lab = NA, x = 'log2FoldChange', y = 'padj',
#                     xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
#                     pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
#                     title = "Volcano plot", subtitle = "SSA/P vs. Normal",
#                     caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
#                     legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
#                     legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

#png("DGE_VolcanoPlots.Salmon.png", width=7, height=7, units = "in", res = 300)
#print(p)
#dev.off()

###Abbreviated volcano plot
p = EnhancedVolcano(data.frame(res), lab = NA, x = 'log2FoldChange', y = 'padj', pCutoff = pCutoff, FCcutoff = FCcutoff, subtitle = bquote(italic("-ve = Q; +ve = W")))
p #About a third of the genes have Padj-value less than 0.05
# ggsave(plot = p, filename = "/Volumes/Pop_Gen/Differential_gene_expression/DGE_volcano_plots/Monomorium_pharaonis.pdf", width = 7, height = 7) #CHANGE NAME

##Using Laurie's customised plots
resultsObject = res
topT <- as.data.frame(resultsObject)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=0.1, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.1))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
#abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
abline(h=-log10(0.05), col="black", lty=4, lwd=1.0) #sets line for Padj value <0.05

#Identify whether positive or negative logfold values are queen biased or worker biased
##By looking at "Ickin1-1" we can tell that negative logfold values are queen biased
identifier <- topT
identifier$id <- rownames(identifier) # negative values are for queen-biased

#Selecting differentially expressed genes based on log2fold(1) and P-value<0.05 criteria
DEGs <- dplyr::filter(topT, (abs(log2FoldChange) > 1), padj <0.05)

#Assign genes as worker-biased or queen-biased in the data frame
# create a new column with a two-level factor
DEGs$biased_gene <- ifelse(DEGs$log2FoldChange >= 0, "worker", "queen")

#Save as csv file
write.csv(DEGs, "/Volumes/Pop_Gen/Differential_gene_expression/DGE_csv_outputs/Monomorium_pharaonis.csv", row.names = T)




##############################

####### Laurie's code ########

##############################



setwd("/Users/laurie/Dropbox/POSTDOC/PROJECT - COMPARATIVE INSECTS/differential_expression/")
library(tximport)
#install.packages('RSQLite', repos='http://cran.us.r-project.org', type="binary")
library(GenomicFeatures)
library(DESeq2)
library(EnhancedVolcano)
files=list.files(recursive=T, pattern="quant.sf")
gff = list.files(pattern="gff")
txdb <- makeTxDbFromGFF(file=gff)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
txi.salmon.g <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi.salmon.t = tximport(files, type="salmon", txOut=T)
#####
## deseq
countData = txi.salmon.g$counts[,1:7]
txi = tximport(files, type = "salmon", tx2gene = tx2gene, 
               txIn = TRUE, txOut = FALSE, countsFromAbundance = "no")
class(txi)
names(txi)
id = c("DRR032045", "DRR032148", "DRR032150", "DRR032224", "DRR032226", "DRR032265", "DRR032266")
dex = c("queen", "worker", "worker", "worker", "worker", "queen", "queen")
colData = data.frame("id"=id, "dex"=dex)
colnames(countData) <- colData$id
dds = DESeqDataSetFromMatrix(
  countData=round(countData), 
  colData = colData,
  design=~dex)
dds = DESeq(dds)
cbind(resultsNames(dds))
res <- results(dds, name = "dex_worker_vs_queen", alpha = 0.05)
summary(res)
mcols(res)$description
head(res)

pCutoff = 0.05
FCcutoff = 1.0
p = EnhancedVolcano(data.frame(res), lab = NA, x = 'log2FoldChange', y = 'padj')
p
resultsObject = res
topT <- as.data.frame(resultsObject)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)




######################

######### END #########

#######################
#Extras from below - I don't think that they're useful










samples <- list.files(path = "/Volumes/Pop_Gen/Differential_gene_expression/quants/", full.names = T)

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "/Volumes/Pop_Gen/Differential_gene_expression/quants/", "") %>% 
  str_replace(".salmon", "")


dir <- "/Volumes/Pop_Gen/Differential_gene_expression/quants/DRR032045_quant/"
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)



files <- file.path(dir, "salmon", samples$run, "quant.sf")
names(files) <- paste0("sample", 1:6)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)
