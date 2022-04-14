library(tidyverse)
library(qiime2R)


setwd("~/Desktop/ANSC595/2022/moving_pictures/from_bell/")

metadata<-read_q2metadata("sample-metadata.tsv")
str(metadata)
colnames(metadata)[3] = "bodysite"
levels(metadata$bodysite)

## Select and keep the columns you need to use as:
##    class, subclass (optional) and subject
metadata_filtered <- metadata %>%
  select(SampleID, bodysite)
str(metadata_filtered)

# Read in ouput from qiime2
taxonomy <- read.delim("collapse.frequency.table.txt", skip = 1)
head(taxonomy)

# Some reformatting
colnames(taxonomy)[1] <- "Taxon"
taxonomy <- taxonomy[,-(ncol(taxonomy))]
taxonomy$Feature.ID <- paste0("ASV", 1:nrow(taxonomy))

# split the taxonomy into one level per column and remove the "K__" stuff
tax.clean<-parse_taxonomy(taxonomy = taxonomy, tax_sep = ";")
tax.clean <- tax.clean[,-7]
head(tax.clean)

# Make the unclassified cells be empty so that we can edit them 
tax.clean[tax.clean=="__"] <- ""
tax.clean[is.na(tax.clean)] <- ""

ncol(tax.clean) ##This should equal 6. If not the code will need to be modified.

# Change empty cells to say "unclassified" with the last assigned taxonomic level
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("unclassified_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:ncol(tax.clean)] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("unclassified_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:ncol(tax.clean)] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("unclassified_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:ncol(tax.clean)] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("unclassified_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:ncol(tax.clean)] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("unclassified_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:ncol(tax.clean)] <- family
  } 
}

# Now finish up formatting into the format required by lefse

# Get taxonomy into one cell separated by "|"
tax.clean_lefse <- tax.clean %>%
  unite("lefse", Kingdom:Genus, sep= "|", 
        remove = FALSE)

# Merge taxonomy with relative abundance for each sample
tax.clean_lefse <- merge(tax.clean_lefse, taxonomy, by.x = 0, by.y = "Feature.ID")
tax.clean_lefse <- tax.clean_lefse[,-3:-9]


#Merge metadata with lefse table. Have to transpose the table to do this.
t.tax <- t(tax.clean_lefse)
meta_lefse <- merge(metadata_filtered, t.tax, by.x = "SampleID", by.y = 0, all = T)

#remove spaces in bodysite names in the bodysite column
meta_lefse$bodysite <- as.data.frame(gsub(" ", "_", meta_lefse$bodysite))

# final formatting to make it pretty
t.meta_lefse <- as.data.frame(t(meta_lefse))
colnames(t.meta_lefse) <- t.meta_lefse[1,]
t.meta_lefse <- t.meta_lefse[-1,]
final_lefse <- t.meta_lefse %>% select(Row.names, lefse, everything())
final_lefse <- final_lefse[,-1]
colnames(final_lefse)[1] <- "SampleID"
final_lefse[1,1] <- "bodysite"

# Save the table `final_lefse` onto your computer
write.csv(final_lefse,"lefse_taxonomy.csv", row.names = F) 

# Final Notes:
# You will need to do some manual reformatting to make this work:
# Make bodysite the first row, SampleID the second row.
# No cell can have spaces. Search for spaces and replace with underscore if needed
# DON'T FORGET: Save as a tab-separated value file in MS excel!!!
