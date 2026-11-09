## set WD 
# mac
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/3-deconPostaaClus97/BOLD')
rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(vegan)
library(car)
library(lme4)
library(lmerTest)  
library(emmeans)
library(glmmTMB) 
library(DHARMa)

#### 1.0 data prep ####
metadata <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable <- read.csv('freqTableCleanByArms.csv', header = TRUE, row.names = 1)
# freqTable.d1 <- read.csv('freqTableCleanByMol.csv', header = TRUE, row.names = 1)
# freqTable <- freqTable.d1[rowSums(freqTable.d1)!=0,]
# freqTable <- read.csv('freqTableCleanBySes.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTableArms <- read.csv('alphaTable.byArms.csv', header = TRUE)
alphaTableMol <- read.csv('alphaTable.byMol.csv', header = TRUE)
alphaTableSes <- read.csv('alphaTable.bySes.csv', header = TRUE)
SOM <- read.csv('../../4-TaxAssign/BOLD/TaxAsn_shelbyOmidori.csv', header = TRUE)
countTable.d1 <- read.csv('../../9-2mmdata/2mmFeqTable.csv', header = TRUE)

# calculate with the marco data
countTable <- countTable.d1 %>% filter(X!="notSure")
richness <- specnumber(t(countTable[,2:19]))

# make data to plot 
countTable.d2 <- countTable[,-1]
countTable.d2[countTable.d2>0] <- 1
alpha2mm <- colSums(countTable.d2)
alpha2mm.ab <- colSums(countTable[,-1])


alphaTable.PLOT <- as.data.frame(matrix(0, nrow=18*5,ncol=1))
names(alphaTable.PLOT) <- "total"

alphaTable.PLOT[,1] <- as.numeric(c(alphaTableArms[1,2:19],
                                    alphaTableMol[1,2:19],alphaTableSes[1,2:19],alpha2mm, alpha2mm.ab))
alphaTable.PLOT$ARMS <- rep(names(alphaTableArms)[2:19],5)
alphaTable.PLOT$fraction <- c(rep("ARMS",18),rep("mol",18),rep("ses",18),
                              rep("2mm",18),rep('2mmab',18))
alphaTable.PLOT$depth <- c(rep(metadata$depth,5))
alphaTable.PLOT$treatment <- c(rep(metadata$treatment1,5))
alphaTable.PLOT$coral <- c(rep(metadata$coral,5))
alphaTable.PLOT$coral <- factor(alphaTable.PLOT$coral, levels=unique(alphaTable.PLOT$coral)[c(2,1,4,6,5,3)]) # reOrder

data.model <- alphaTable.PLOT

# what's the distribution of the data? 
# 2mm
hist(data.model%>% filter(fraction=='2mm') %>% .$total)
mean(data.model%>% filter(fraction=='2mm') %>% .$total) # 43.55556 
var(data.model%>% filter(fraction=='2mm') %>% .$total) # 48.84967
shapiro.test(data.model%>% filter(fraction=='2mm') %>% .$total) # p = 0.5479

# coi
hist(data.model%>% filter(fraction=='ARMS') %>% .$total)
mean(data.model%>% filter(fraction=='ARMS') %>% .$total) # 806.9444 
var(data.model%>% filter(fraction=='ARMS') %>% .$total) # 6913.82
shapiro.test(data.model%>% filter(fraction=='ARMS') %>% .$total) # p = 0.2002




#### 2.0 modelling ####
##### 2.1 GLMM #####
glmm.model_2mm <- glmmTMB(total ~ coral + (1 | depth), 
                      data = data.model%>% filter(fraction=='2mm'),
                      family = poisson(link = "log"))

glmm.model_arms <- glmmTMB(total ~ coral  + (1 | depth), 
                         data = data.model%>% filter(fraction=='ARMS'),
                         family = nbinom2(link = "log"))

##### 2.2 GLM #####
glm.model_2mm <- glmmTMB(total ~ coral, 
                        data = data.model%>% filter(fraction=='2mm'),
                        family = poisson(link = "log"))

glm.model_arms <- glmmTMB(total ~ coral, 
                            data = data.model%>% filter(fraction=='ARMS'),
                            family = nbinom2(link = "log"))


##### 2.3 LM #####
lm.model_2mm <- lm(total ~ coral, 
                           data = data.model%>% filter(fraction=='2mm'))

lm.model_arms <- lm(total ~ coral, 
                         data = data.model%>% filter(fraction=='ARMS'))


#### 3.0 model check ####
# AIC
AIC(glmm.model_2mm,glm.model_2mm,lm.model_2mm)
AIC(glmm.model_arms,glm.model_arms,lm.model_arms)

# model check 
sim_res.2mm.glmm <- simulateResiduals(glmm.model_2mm)
sim_res.2mm.glm <- simulateResiduals(glm.model_2mm)
sim_res.2mm.lm <- simulateResiduals(lm.model_2mm)

sim_res.arms.glmm <- simulateResiduals(glmm.model_arms)
sim_res.arms.glm <- simulateResiduals(glm.model_arms)
sim_res.arms.lm <- simulateResiduals(lm.model_arms)


# Check overall fit
plot(sim_res.2mm.glmm)
plot(sim_res.2mm.glm)
plot(sim_res.2mm.lm)

plot(sim_res.arms.glmm)
plot(sim_res.arms.glm)
plot(sim_res.arms.lm)


# Test for overdispersion
testOverdispersion(sim_res.2mm.glmm)
testOverdispersion(sim_res.2mm.glm)
testOverdispersion(sim_res.arms.glmm)
testOverdispersion(sim_res.arms.glm)
testOverdispersion(sim_res.arms.lm)
testOverdispersion(sim_res.2mm.lm)

testResiduals(sim_res.arms.glmm)
testResiduals(sim_res.2mm.glmm)
testResiduals(sim_res.arms.glm)
testResiduals(sim_res.2mm.glm)
testResiduals(sim_res.arms.lm)
testResiduals(sim_res.2mm.lm)

# sigularity 
performance::check_singularity(sim_res.arms.glmm)
performance::check_singularity(sim_res.2mm.glmm)

# Summary
summary(glmm.model_2mm)   
confint(glmm.model_2mm, level = 0.95)
summary(glm.model_2mm)   
confint(glm.model_2mm, level = 0.95)
summary(lm.model_2mm)
confint(lm.model_2mm, level = 0.95)

summary(glmm.model_arms)
confint(glmm.model_arms, level = 0.95)
summary(glm.model_arms)
confint(glm.model_arms, level = 0.95)
summary(lm.model_arms)
confint(lm.model_arms, level = 0.95)


Anova(lm.model_2mm)
Anova(lm.model_arms)


# Get estimated marginal means
emm.2mm <- emmeans(lm.model_2mm, ~ coral)


# Pairwise comparisons with Tukey adjustment
posthoc.2mm <- pairs(emm.2mm, adjust = "tukey")
posthoc.2mm

