## This script goes through picrust output and then uses
## DESeq2 to find differentially abundant pathways 


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
#pathway_abundance.qza

##############################################

setwd("D:/Dropbox/Purdue Class from 2020/2022 Spring ANSC 595 Microbiome/Project Assignment/10 Picrust/day21")
list.files()

if(!dir.exists("output"))
  dir.create("output")

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

#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1

#First get the pathway table from physeq
# Note that the samples in both the metadata and the 
# picrust output need to be the same 
physeq <- qza_to_phyloseq(
  features="pathway_abundance_day21.qza",
  metadata = "2021-31_metadata_ver.4.0_day22.txt"
)

pathway_table <- data.frame(otu_table(physeq), check.names = FALSE)

pathway_table2 <- pathway_table + 1


#Now make the phyloseq object:


table.physeq = otu_table(as.matrix(pathway_table2), taxa_are_rows=TRUE)
#tax.physeq = tax_table(as.matrix(tax.clean))  
meta.physeq = sample_data(metadata.filtered)


#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(table.physeq, meta.physeq)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.

# Enter the statistical model you want to test. Eg. `~ diet`
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

sigtab = res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
#sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
#head(sigtab)




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
x = sort(x, TRUE)
sigtab$pathway = factor(as.character(sigtab$pathway), levels=names(x))

sigtab$pathway <- row.names(sigtab)

DESeq_fig = ggplot(sigtab, aes(x = pathway, y = log2FoldChange)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  #scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
  #ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))+
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12)) 

ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 6, width = 5)

