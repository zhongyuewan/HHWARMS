## set WD 
#mac 
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/3-deconPostaaClus97/BOLD')
rm(list=ls())

## load lib
library(dplyr)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(ggrepel)
library(indicspecies)
library(glmmTMB) 
library(car)
library(lme4)
library(emmeans)

## the color
myColor <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

## get data
metadata <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTable <- read.csv('alphaTable.byArms.csv', header = TRUE)
SOM <- read.csv('../../4-TaxAssign/TaxAsn_shelbyOmidori.csv', header = TRUE)


# check out sharing OTUs by time 
shaARMS <- metadata[(metadata$depth=="shallow"),]$arms
midTARMS <- metadata[(metadata$depth=="middle"),]$arms
depARMS <- metadata[(metadata$depth=="deep"),]$arms

conARMS <- metadata[(metadata$treatment1=="control"),]$arms
monoARMS <- metadata[(metadata$treatment1=="mono"),]$arms
polyARMS <- metadata[(metadata$treatment1=="poly"),]$arms
mixARMS <- metadata[(metadata$treatment1=="mix"),]$arms

AcroARMS <- metadata[grepl("acropora", metadata$coral),]$arms
PavoARMS <- metadata[grepl("pavona", metadata$coral),]$arms
PlatARMS <- metadata[grepl("platygyra", metadata$coral),]$arms

# calculate datamatrix based on squaroot 
# twist the freqTable to switch among motile/sessile/ARMS
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="byARMS")
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="sessile")
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="motile")

freqTable.d2[,4:21] <- sqrt(freqTable.d2[,4:21]) 
freqTable.d3 <- freqTable.d2 %>% filter(phylum!="notAssigned")
freqTable.d4 <- freqTable.d3[,3:21]


feqTable <- as.data.frame(freqTable.d4 %>%
  group_by(phylum)%>% 
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE)
  ))

# row.names(feqTable) <- feqTable[,1]
# feqTable <- feqTable[2:19]

## make long table for analysis 
long_data <- as.data.frame(feqTable %>%
  pivot_longer(
    cols = -phylum,  # All columns except the genus column
    names_to = "arms",
    values_to = "count"))

long_data$depth <- metadata$depth
long_data$treatment <- metadata$treatment
long_data$coral <- factor(metadata$coral, 
                          levels = c('no','acropora','pavona',
                                     'platygyra','acroporaplatygyra','acroporapavonaplatygyra'))

# let's plot
plot1 <- ggplot(long_data%>% filter(phylum == "Mollusca"), aes(y=count, x = coral))+
  geom_boxplot()

# let's anova 
long_data_total <- as.data.frame(long_data%>% group_by(arms) %>% 
  summarise(total = sum(count),
            phylum="all",
            arms = unique(arms),
            depth=unique(depth),
            treatment=unique(treatment),
            coral=unique(coral)))

model1 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Arthropoda"),
                  family = gaussian(link = "identity"))


summary(model1)
Anova(model1)
emm1 <- emmeans(model1, ~ coral)
posthoc.1 <- pairs(emm1, adjust = "tukey")

model2 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Annelida"),
                  family = gaussian(link = "identity"))
summary(model2)
Anova(model2)
emm2 <- emmeans(model2, ~ coral)
posthoc.2 <- pairs(emm2, adjust = "tukey")

model3 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Bacillariophyta"),
                  family = gaussian(link = "identity"))

summary(model3)
Anova(model3)
emm3 <- emmeans(model3, ~ coral)
posthoc.3 <- pairs(emm3, adjust = "tukey")

model4 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Mollusca"),
                  family = gaussian(link = "identity"))

summary(model4)
Anova(model4)
emm4 <- emmeans(model4, ~ coral)
posthoc.4 <- pairs(emm4, adjust = "tukey")

model5 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Porifera"),
                  family = gaussian(link = "identity"))

# model5 <- glmmTMB(count ~ coral, # remove depth for motile porifera 
                  data = long_data %>% filter(phylum == "Porifera"),
                  family = gaussian(link = "identity"))

summary(model5)
Anova(model5)
emm5 <- emmeans(model5, ~ coral)
posthoc.5 <- pairs(emm5, adjust = "tukey")


model6 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Rhodophyta"),
                  family = gaussian(link = "identity"))


summary(model6)
Anova(model6)
emm6 <- emmeans(model5, ~ coral)
posthoc.6 <- pairs(emm5, adjust = "tukey")


# check model 
# Check residuals
qqnorm(residuals(model1))
qqnorm(residuals(model2))
qqnorm(residuals(model3))
qqnorm(residuals(model4))
qqnorm(residuals(model5))
qqnorm(residuals(model6))

# Shapiro-Wilk test
shapiro.test(residuals(model1))
shapiro.test(residuals(model2))
shapiro.test(residuals(model3))
shapiro.test(residuals(model4))
shapiro.test(residuals(model5))
shapiro.test(residuals(model6))

# all checked out 

