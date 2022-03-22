#Install qiime2R if you need to:
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")

#Install tidyverse if you need to:
# Install from CRAN
install.packages("tidyverse")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("tidyverse/tidyverse")

#Load the packages

library(qiime2R)

library(tidyverse)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/*.qz* .
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/core-metrics-results/*.qz* core-metrics-results/
# scp <user.name>@bell.rcac.purdue.edu:/depot/microbiome/data/ANSC595/class_materials/moving_pictures_pipeline/qiime_out/*.tsv .
##############################################


setwd("~/Desktop/ANSC595/moving-pictures/")
#setwd("/Users/john2185/Documents/")
list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("sample-metadata.tsv")
#metadata2 <- read.delim("sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#metadata2 <- metadata2[-1,]
str(metadata)
levels(metadata$`body-site`)


row.names(metadata) <- metadata[ ,1]
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

body_colors <- c("Black", "Blue", "Green", "Gray")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

my_column <- "body-site"
#my_column <- "DietTreatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "body-site"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= subject)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
  #scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

##SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Body Site")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= subject), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Body Site")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

