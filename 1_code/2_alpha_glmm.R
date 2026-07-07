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

## get data
metadate <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable <- read.csv('freqTableCleanByArms.csv', header = TRUE, row.names = 1)
# freqTable.d1 <- read.csv('freqTableCleanByMol.csv', header = TRUE, row.names = 1)
# freqTable <- freqTable.d1[rowSums(freqTable.d1)!=0,]
# freqTable <- read.csv('freqTableCleanBySes.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTableArms <- read.csv('alphaTable.byArms.csv', header = TRUE)
alphaTableMol <- read.csv('alphaTable.byMol.csv', header = TRUE)
alphaTableSes <- read.csv('alphaTable.bySes.csv', header = TRUE)
SOM <- read.csv('../../4-TaxAssign/BOLD/TaxAsn_shelbyOmidori.csv', header = TRUE)
countTable <- read.csv('../../9-2mmdata/2mmFeqTable.csv', header = TRUE)

# calculate with the marco data
shannon <- diversity(t(countTable[,2:19]), index = "shannon") 
richness <- specnumber(t(countTable[,2:19]))
evenness <- shannon / log(richness)

# make data to plot 
countTable.d1 <- countTable[,-1]
countTable.d1[countTable.d1>0] <- 1
alpha2mm <- colSums(countTable.d1)
alpha2mm.ab <- colSums(countTable[,-1])


alphaTable.PLOT <- as.data.frame(matrix(0, nrow=18*7,ncol=1))
names(alphaTable.PLOT) <- "total"

alphaTable.PLOT[,1] <- as.numeric(c(alphaTableArms[1,2:19],
                                    alphaTableMol[1,2:19],alphaTableSes[1,2:19],alpha2mm, alpha2mm.ab,
                                    shannon,evenness))
alphaTable.PLOT$ARMS <- rep(names(alphaTableArms)[2:19],7)
alphaTable.PLOT$fraction <- c(rep("ARMS",18),rep("mol",18),rep("ses",18),
                              rep("2mm",18),rep('2mmab',18),rep("shann",18),rep('even',18))
alphaTable.PLOT$depth <- c(rep(metadate$depth,7))
alphaTable.PLOT$treatment <- c(rep(metadate$treatment1,7))
alphaTable.PLOT$coral <- c(rep(metadate$coral,7))
alphaTable.PLOT$coral <- factor(alphaTable.PLOT$coral, levels=unique(alphaTable.PLOT$coral)[c(3,2,1,4,6,5)]) # reOrder
# what if i try put pavona first 
# alphaTable.PLOT$coral <- factor(alphaTable.PLOT$coral, levels=unique(alphaTable.PLOT$coral)[c(1,2,4,6,5,3)]) # reOrder

# let's plot

plotuniOTUs <- ggplot(alphaTable.PLOT, aes(x=fraction, y=total, fill = treatment)) +
  geom_boxplot() +
  theme_classic()+
  ylab("Unique OTUs") 

plotuniOTUs.depth <- ggplot(alphaTable.PLOT, aes(x=fraction, y=total, fill = depth)) +
  geom_boxplot() +
  theme_classic()+
  ylab("Unique OTUs") 

plotuniOTUs.coral <- ggplot(alphaTable.PLOT, aes(x=fraction, y=total,fill=coral)) +
  geom_boxplot() +
  xlab("Treatments") + 
  ylab("Unique OTUs") +
  theme_classic()

plotuniOTUs.treatment <- ggplot(alphaTable.PLOT, aes(x=fraction, y=total,fill=treatment)) +
  geom_boxplot() +
  xlab("Treatments") + 
  ylab("Unique OTUs") +
  theme_classic()

# what about the stats 
# data structure 
data.model <- alphaTable.PLOT

# what about glm
# Using glmmTMB (supports negative binomial)
model_po.2mm <- glmmTMB(total ~ coral + (1 | depth), 
                      data = data.model%>% filter(fraction=='2mm'),
                      family = poisson(link = "log"))

model_po.null <- glmmTMB(total ~ 1 + (1 | depth), 
                        data = data.model%>% filter(fraction=='2mm'),
                        family = poisson(link = "log"))

model_nb.arms <- glmmTMB(total ~ coral  + (1 | depth), 
                         data = data.model%>% filter(fraction=='ARMS'),
                         family = nbinom2(link = "log"))

model_po.shann <- glmmTMB(total ~ coral + (1 | depth), 
                        data = data.model%>% filter(fraction=='shann'),
                        family = Gamma(link = "log"))

model_po.even <- glmmTMB(total ~ coral + (1 | depth), 
                          data = data.model%>% filter(fraction=='even'),
                          family = beta_family(link = "logit"))

# model check 
sim_res.2mm <- simulateResiduals(model_po.2mm)
sim_res.arms <- simulateResiduals(model_nb.arms)

# Check overall fit
plot(sim_res.2mm)
plot(sim_res.arms)

# Test for overdispersion
testOverdispersion(sim_res.2mm)
testOverdispersion(sim_res.arms)
testResiduals(sim_res.2mm)
testResiduals(sim_res.arms)


# Summary
summary(model_po.2mm)   
summary(model_po.null)
summary(model_nb.arms)
summary(model_po.shann)
summary(model_po.even)

Anova(model_po.2mm, type=2)
Anova(model_nb.arms,type=2)
Anova(model_po.shann)
Anova(model_po.even)


# Get estimated marginal means
emm.2mm <- emmeans(model_po.2mm, ~ coral)
emm.arms <- emmeans(model_nb.arms, ~ coral)
emm.shann <- emmeans(model_po.shann, ~ coral)


# Pairwise comparisons with Tukey adjustment
posthoc.2mm <- pairs(emm.2mm, adjust = "tukey")
posthoc.arms <- pairs(emm.arms, adjust = "tukey")
posthoc.shann <- pairs(emm.shann, adjust = "tukey")

# View results
posthoc.2mm
posthoc.arms
posthoc.shann


# Summary with confidence intervals
summary(posthoc, infer = c(TRUE, TRUE))


# check dispersal  
# ARMS
  data.model %>% filter(fraction == 'ARMS') %>%
  group_by(coral) %>%
  summarise(mean = mean(total), var = var(total), rio=var/mean, n = n())
# 2mm
  data.model %>% filter(fraction == '2mm') %>%
     group_by(coral) %>%
     summarise(mean = mean(total), var = var(total),rio=var/mean, n = n())

# anova? 
leveneTest(total ~ coral, data = data.model %>% filter(fraction == '2mm'))
a1 <- aov(total ~ coral, data=data.model%>% filter(fraction=='2mm'))
summary(a1)
TukeyHSD(a1)
a1.1 <- lm(total ~ coral, data=data.model%>% filter(fraction=='2mm'))
summary(a1.1)


leveneTest(total ~ coral, data = data.model %>% filter(fraction == 'ARMS'))
a2 <- aov(total ~ coral, data=data.model%>% filter(fraction=='ARMS'))
summary(a2)

leveneTest(total ~ coral, data = data.model %>% filter(fraction == 'shann'))
a3 <- aov(total ~ coral, data=data.model%>% filter(fraction=='shann'))
summary(a3)

leveneTest(total ~ coral, data = data.model %>% filter(fraction == 'even'))
a4 <- aov(total ~ coral, data=data.model%>% filter(fraction=='even'))
summary(a4)


hist(data.model%>% filter(fraction=='2mm')%>% .$total)
hist(data.model%>% filter(fraction=='ARMS')%>% .$total)
hist(data.model%>% filter(fraction=='shann')%>% .$total)
hist(data.model%>% filter(fraction=='even')%>% .$total)









