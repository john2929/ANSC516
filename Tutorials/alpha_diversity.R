##############################################################
# title: "Alpha diversity in R - qiime2 output"
# author: "ANSC595"
# date: "March 16, 2021"
##############################################################

#The first step is very important. You need to set your working 
#directory. Just like in unix we have to `cd` to where our data is, the 
#same thing is true for R.
##############################################

#So where `~/Desktop/ANSC595/moving-pictures` is in the code below, you 
#need to enter the path to where you saved the tutorial or your data.

setwd("~/Desktop/ANSC516/2023/ANSC516-repo/ANSC516/data/moving-pictures")
list.files()

# Modified from the original online version available at 
# http://rpubs.com/dillmcfarlan/R_microbiotaSOP

# and Tutorial: Integrating QIIME2 and R for data visualization 
# and analysis using qiime2R
# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

##Goal
# The goal of this tutorial is to demonstrate basic analyses of microbiota data to determine if and how communities differ by variables of interest. In general, this pipeline can be used for any microbiota data set that has been clustered into operational taxonomic units (OTUs).
#
# This tutorial assumes some basic statistical knowledge. Please consider if your data fit the assumptions of each test (normality? equal sampling? Etc.). If you are not familiar with statistics at this level, we strongly recommend collaborating with someone who is. The incorrect use of statistics is a pervasive and serious problem in the sciences so don't become part of the problem! That said, this is an introductory tutorial and there are many, many further analyses that can be done with microbiota data. Hopefully, this is just the start for your data!

##Data
# The data used here are from the qiime2 moving pictures tutorial. 
# Please see their online tutorial for an explanation of the dataset.

##Files
# We will use the following files created using the qiime2 moving pictures tutorial.

# core-metrics-results/evenness_vector.qza (alpha diversity)
# core-metrics-results/faith_pd_vector.qza (alpha diversity)
# core-metrics-results/observed_features_vector.qza (alpha diversity)
# core-metrics-results/shannon_vector.qza (alpha diversity)
# sample-metadata.tsv (metadata)


# Data manipulation
## Load Packages

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","

# Load data

meta<-read_q2metadata("sample-metadata.tsv")
str(meta)
colnames(meta)[3] <- "body.site"
colnames(meta)[8] <- "reported.antibiotic.usage"
colnames(meta)[9] <- "days.since.experiment.start"
str(meta)

evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(meta)
#observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
#observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)



###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-1]
str(meta)


#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or usenon-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(meta$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon_entropy, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")

install.packages("ggpubr")
library("ggpubr")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

# We see that, as expected from the graphs, shannon and evenness 
# are normally distributed.


#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest are

# body site: gut, tongue, right palm, left palm
# subject: 1 and 2
# month-year: 10-2008, 1-2009, 2-2009, 3-2009, 4-2009

## Categorical variables
# Now that we know which tests can be used, let's run them. 

## Normally distributed metrics

# Since it's the closest to normalcy, we will use **Evenness** as an 
#example. First, we will test body site, which is a categorical variable 
# with more than 2 levels. Thus, we run ANOVA. If age were only two 
# levels, we could run a t-test

# Does body site impact the Evenness of the microbiota?

#Run the ANOVA and save it as an object
aov.evenness.body_site = aov(pielou_evenness ~ body.site, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.evenness.body_site)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.evenness.body_site)

# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.

levels(meta$body.site)
#Re-order the groups because the default is alphabetical order
meta$body.site.ord = factor(meta$body.site, c("left palm", "right palm", "gut", "tongue"))
levels(meta$body.site.ord)

#Plot
boxplot(pielou_evenness ~ body.site.ord, data=meta, ylab="Peilou evenness")

evenness_boxplot <- ggplot(meta, aes(body.site.ord, pielou_evenness)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/evenness_boxplot.png", evenness, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.

evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(body.site.ord) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

evenness_se <- ggplot(evenness_summary, aes(body.site.ord, mean_evenness, fill = body.site.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/evenness_se.png", evenness_se, height = 2.5, width = 3)


## **Non-normally distributed metrics**

# We will use **Faith's phylogenetic diversity** here. Since body site 
# is categorical, we use Kruskal-Wallis (non-parametric equivalent of 
# ANOVA). If we have only two levels, we would run Wilcoxon rank sum 
# test (non-parametric equivalent of t-test)

kruskal.test(faith_pd ~ body.site.ord, data=meta)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(meta$faith_pd, meta$body.site.ord, p.adjust.method="BH")

# Like evenness, we see that pd also increases with age.

#Plot
boxplot(faith_pd ~ body.site.ord, data=meta, ylab="Faith phylogenetic diversity")

# or with ggplot2

faith_pd_boxplot <- ggplot(meta, aes(body.site.ord, faith_pd)) + 
  geom_boxplot(aes(color = body.site.ord)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd.png", faith_pd_boxplot, height = 3, width = 3)

##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since days.since.experiment.start is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

# Does days.since.experiment.start impact evenness of the microbiota?

glm.evenness.time = glm(pielou_evenness ~ days.since.experiment.start, data=meta)
summary(glm.evenness.time)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(pielou_evenness ~ days.since.experiment.start, data=meta)
#Add the glm best fit line
plot(pielou_evenness ~ days.since.experiment.start, data=meta) + abline(glm.evenness.time)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.

gaussian.faith.time = glm(faith_pd ~ days.since.experiment.start, data=meta, family="gaussian")
plot(gaussian.faith.time, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.time = glm(faith_pd ~ days.since.experiment.start, data=meta, family="quasipoisson")
plot(qp.faith.time, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.time)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and days.since.experiment.start.

#Plot
plot(log(faith_pd) ~ days.since.experiment.start, data=meta, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ days.since.experiment.start, data=meta, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.time)

