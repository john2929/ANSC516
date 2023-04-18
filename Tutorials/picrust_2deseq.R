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

setwd("~/Desktop/ANSC516/2023/from_cluster/moving-pictures-picrust2/")
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
metadata<-read_q2metadata("sample-metadata.tsv")
str(metadata)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[,-1]
colnames(metadata)[2] = "body.site"
levels(metadata$body.site)
metadata$body.site.ord = factor(metadata$body.site, c("left palm", "right palm", "gut", "tongue"))
levels(metadata$body.site.ord)
colnames(metadata)[7] <- "reported.antibiotic.usage"
colnames(metadata)[8] <- "days.since.experiment.start"
str(metadata)

#Read in the output from picrust2
picrust_table <- read_qza("pathway_abundance.qza")
picrust_table <- picrust_table$data
picrust_table <- data.frame(picrust_table, check.names = F)
picrust_table <- picrust_table + 1

#Make sure the same samples are in both the picrust table and the metadata
metadata.filtered = metadata[rownames(metadata) %in% colnames(picrust_table),]

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(picrust_table), taxa_are_rows=TRUE)

meta.physeq = sample_data(metadata.filtered)

#We then merge these into an object of class phyloseq.

physeq_picrust = phyloseq(OTU.physeq, meta.physeq)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.

# Enter the statistical model you want to test. Eg. `~ diet`
diagdds = phyloseq_to_deseq2(physeq_picrust, ~ body.site)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:

#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

alpha = 0.05
my_contrast = c("body.site", "gut", "tongue") 

res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)

sigtab = res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtab2 <- subset(sigtab, abs(log2FoldChange) > 3)
#sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
#head(sigtab)


#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.

#add in pathway names
sigtab$pathway <- row.names(sigtab)
str(sigtab)
pathway_mapfile <- read.delim("mapfiles/metacyc_pathways_info.txt", sep = "\t", header = F, quote = "", stringsAsFactors = F)
str(pathway_mapfile)
colnames(pathway_mapfile)[1] <- "abbreviation"
colnames(pathway_mapfile)[2] <- "pathway_name"

#merge `sigtab` and `pathway_mapfile`
sigtab <- inner_join(sigtab, pathway_mapfile, by = c("pathway" = "abbreviation"))

# pathway order
x = tapply(sigtab$log2FoldChange, sigtab$pathway_name, function(x) max(x))
x = sort(x, TRUE)
sigtab$pathway_name = factor(as.character(sigtab$pathway_name), levels=names(x))

DESeq_fig = ggplot(sigtab, aes(x = pathway_name, y = log2FoldChange)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  #scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
  #ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))+
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12)) 

ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 10, width = 40)

##############################################################
#We might need to reduce the number of pathways 
#we want to consider.
##############################################################

sigtab = res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtab <- subset(sigtab, abs(log2FoldChange) > 3)

#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.

#add in pathway names
sigtab$pathway <- row.names(sigtab)
str(sigtab)
pathway_mapfile <- read.delim("mapfiles/metacyc_pathways_info.txt", sep = "\t", header = F, quote = "", stringsAsFactors = F)
str(pathway_mapfile)
colnames(pathway_mapfile)[1] <- "abbreviation"
colnames(pathway_mapfile)[2] <- "pathway_name"

#merge `sigtab` and `pathway_mapfile`
sigtab <- inner_join(sigtab, pathway_mapfile, by = c("pathway" = "abbreviation"))

# pathway order
x = tapply(sigtab$log2FoldChange, sigtab$pathway_name, function(x) max(x))
x = sort(x, TRUE)
sigtab$pathway_name = factor(as.character(sigtab$pathway_name), levels=names(x))

DESeq_fig = ggplot(sigtab, aes(x = pathway_name, y = log2FoldChange)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  #scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
  #ylim(0,8) +
  geom_text(color="black", x=length(unique(sigtab$pathway_name))-5, y=max(sigtab$log2FoldChange)-1, label=my_contrast[2], show_guide = F) +
  geom_text(color="black", x=5, y=min(sigtab$log2FoldChange)+1, label=my_contrast[3], show_guide = F) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))+
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y =element_text(color = "black", size = 12)) 

ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], "-filtered.png"), DESeq_fig, height = 8, width = 15)


##############################################################
#Lets try this with a function

run_deseq2 <- function(my_factor, x, y){
  
  my_contrast <- c(my_factor, x, y)
  res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)
  
  sigtab = res[which(res$padj < alpha), ]
  sigtab <- as(sigtab, "data.frame")
  sigtab <- subset(sigtab, abs(log2FoldChange) > 3)
  #sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
  #head(sigtab)
  
  
  #Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.
  
  #add in pathway names
  sigtab$pathway <- row.names(sigtab)
  str(sigtab)
  pathway_mapfile <- read.delim("mapfiles/metacyc_pathways_info.txt", sep = "\t", header = F, quote = "", stringsAsFactors = F)
  str(pathway_mapfile)
  colnames(pathway_mapfile)[1] <- "abbreviation"
  colnames(pathway_mapfile)[2] <- "pathway_name"
  
  #merge `sigtab` and `pathway_mapfile`
  sigtab <- inner_join(sigtab, pathway_mapfile, by = c("pathway" = "abbreviation"))
  
  # pathway order
  x = tapply(sigtab$log2FoldChange, sigtab$pathway_name, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$pathway_name = factor(as.character(sigtab$pathway_name), levels=names(x))
  
  DESeq_fig = ggplot(sigtab, aes(x = pathway_name, y = log2FoldChange)) + 
    geom_point(size=3) + 
    ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
    #scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
    #ylim(0,8) +
    geom_text(color="black", x=length(unique(sigtab$pathway_name))-5, y=max(sigtab$log2FoldChange)-1, label=my_contrast[2], show_guide = F) +
    geom_text(color="black", x=5, y=min(sigtab$log2FoldChange)+1, label=my_contrast[3], show_guide = F) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.2))+
    theme(axis.text.x = element_text(color = "black", size = 10),
          axis.text.y =element_text(color = "black", size = 12)) 
  
  ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], "-filt.png"), DESeq_fig, height = 10, width = 20)
}

run_deseq2("body.site", "gut", "left palm")
run_deseq2("body.site", "gut", "right palm")
run_deseq2("body.site", "gut", "tongue") 
run_deseq2("body.site", "tongue", "left palm") 
run_deseq2("body.site", "tongue", "right palm") 
run_deseq2("body.site", "right palm", "left palm") 
