## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
##############################################

setwd("/Users/john2185/Desktop/ANSC595/2022/StudentHelp/picrust/")
list.files()

if(!dir.exists("output"))
  dir.create("output")

#Load in the metadata
metadata<-read_q2metadata("2021-31_metadata_ver.4.0_day22.txt")
rownames(metadata) <- metadata$SampleID
metadata <- metadata[,-1]
str(metadata)
#################################################################
###Differential Abundance with DESeq2
#################################################################


#Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.
#If you need help  with DESeq2 install, see this website
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")

#To use DESeq, we need no zeros in our OTU table. 
#So we will edit the table by multiplying by 2 and + 1

#First get the OTU table from physeq
physeq <- qza_to_phyloseq(
  features="pathway_abundance_day21.qza",
  metadata = "2021-31_metadata_ver.4.0_day22.txt"
)

physeq_otu_table <- data.frame(otu_table(physeq), check.names = FALSE)

OTU.clean2 <- physeq_otu_table + 1


#Now make the phyloseq object:

OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=TRUE)
#tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(metadata)


#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(OTU.physeq, meta.physeq)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.


diagdds = phyloseq_to_deseq2(physeq_deseq, ~ diet)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:

#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

#Normal CP vs. High-indigestible CP
alpha = 0.05
my_contrast = c("diet", "Normal CP", "High-indigestible CP") 

res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)

sigtab <- res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
#sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
#head(sigtab)


###Volcano Plot

with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
#x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = sigtab$log2FoldChange
names(x) <- row.names(sigtab)

x = sort(x, TRUE)

sigtab$pathway <- row.names(sigtab)
sigtab$pathway = factor(as.character(sigtab$pathway), levels=names(x))
str(sigtab$pathway)

DESeq_fig = ggplot(sigtab, aes(x = pathway, y = log2FoldChange)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  #scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
  #ylim(0,8) +
  geom_hline(yintercept = 0) +
  geom_text(x=9, y=1.5, label=my_contrast[2]) +
  geom_text(x=3, y=-.75, label=my_contrast[3]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12)) 

ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 4.5, width = 5)

