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
library(RColorBrewer)
library(DHARMa)

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
freqTable.d2.arms <- freqTable.d1%>%filter(fraction=="byARMS")
freqTable.d2.mo <- freqTable.d1%>%filter(fraction=="motile")
freqTable.d2.se <- freqTable.d1%>%filter(fraction=="sessile")


freqTable.d2.arms[,4:21] <- freqTable.d2.arms[,4:21] 
freqTable.d3.arms <- freqTable.d2.arms %>% filter(phylum!="notAssigned")
freqTable.d4.arms <- freqTable.d3.arms[,3:21]

freqTable.d2.mo[,4:21] <- freqTable.d2.mo[,4:21] 
freqTable.d3.mo <- freqTable.d2.mo %>% filter(phylum!="notAssigned")
freqTable.d4.mo <- freqTable.d3.mo[,3:21]

freqTable.d2.se[,4:21] <- freqTable.d2.se[,4:21] 
freqTable.d3.se <- freqTable.d2.se %>% filter(phylum!="notAssigned")
freqTable.d4.se <- freqTable.d3.se[,3:21]


feqTable.arms <- as.data.frame(freqTable.d4.arms %>%
  group_by(phylum)%>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)))

feqTable.mo <- as.data.frame(freqTable.d4.mo %>%
  group_by(phylum)%>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)))

feqTable.se <- as.data.frame(freqTable.d4.se %>%
  group_by(phylum)%>% 
  summarise(across(where(is.numeric), sum, na.rm = TRUE)))



## make long table for analysis 

for (i in c('.arms','.mo','.se')) {
  
  feqTable <- get(paste0('feqTable',i))

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
                  family = nbinom2(link = 'log'))
assign(paste0("model1", i), model1)


performance::check_singularity(model1)
sim_res.model1 <- simulateResiduals(model1)
testResiduals(sim_res.model1)

summary(model1)
Anova(model1)
emm1 <- emmeans(model1, ~ coral)
posthoc.1 <- pairs(emm1, adjust = "tukey")

model2 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Annelida"),
                  family=nbinom2(link = 'log'))
assign(paste0("model2", i), model2)

performance::check_singularity(model2)
sim_res.model2 <- simulateResiduals(model2)
testResiduals(sim_res.model2)
summary(model2)
Anova(model2)
emm2 <- emmeans(model2, ~ coral)
posthoc.2 <- pairs(emm2, adjust = "tukey")

model3 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Bacillariophyta"),
                  family = nbinom2(link='log'))
assign(paste0("model1", 3), model3)

performance::check_singularity(model3)
sim_res.model3 <- simulateResiduals(model3)
testResiduals(sim_res.model3)

summary(model3)
Anova(model3)
emm3 <- emmeans(model3, ~ coral)
posthoc.3 <- pairs(emm3, adjust = "tukey")
assign(paste0("model3", i), model3)


model4 <- glmmTMB(count ~ coral + (1 | depth), 
                    data = long_data %>% filter(phylum == "Cnidaria"),
                    family=nbinom2(link = "log"))
                    
assign(paste0("model4", i), model4)
summary(model4)

performance::check_singularity(model4)
sim_res.model4 <- simulateResiduals(model4)
testResiduals(sim_res.model4)
summary(model4)



Anova(model4)


emm4 <- emmeans(model4, ~ coral)
posthoc.4 <- pairs(emm4, adjust = "tukey")

model5 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Mollusca"),
                  family=nbinom2(link = "log"))

model5
summary(model5)
assign(paste0("model5", i), model5)

performance::check_singularity(model5)
sim_res.model5 <- simulateResiduals(model5)
testResiduals(sim_res.model5)
summary(model5)



Anova(model5)


emm5 <- emmeans(model5, ~ coral)
posthoc.5 <- pairs(emm4, adjust = "tukey")

model6 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Porifera"),
                  family=nbinom2(link='log'))
assign(paste0("model6", i), model6)

performance::check_singularity(model6)
sim_res.model6 <- simulateResiduals(model6)
testResiduals(sim_res.model6)

summary(model6)
Anova(model6)
emm6 <- emmeans(model6, ~ coral)
posthoc.6 <- pairs(emm6, adjust = "tukey")


model7 <- glmmTMB(count ~ coral + (1 | depth), 
                  data = long_data %>% filter(phylum == "Rhodophyta"),
                  family = nbinom2(link='log'))

performance::check_singularity(model7)
sim_res.model7 <- simulateResiduals(model7)
testResiduals(sim_res.model7)

summary(model7)
Anova(model7)
emm7 <- emmeans(model6, ~ coral)
posthoc.7 <- pairs(emm6, adjust = "tukey")
assign(paste0("model7", i), model7)

assign(paste0("posthoc.1", i), posthoc.1)
assign(paste0("posthoc.2", i), posthoc.2)
assign(paste0("posthoc.3", i), posthoc.3)
assign(paste0("posthoc.4", i), posthoc.4)
assign(paste0("posthoc.5", i), posthoc.5)
assign(paste0("posthoc.6", i), posthoc.6)
assign(paste0("posthoc.7", i), posthoc.7)

assign(paste0("sim_res.model1", i), sim_res.model1)
assign(paste0("sim_res.model2", i), sim_res.model2)
assign(paste0("sim_res.model3", i), sim_res.model3)
assign(paste0("sim_res.model4", i), sim_res.model4)
assign(paste0("sim_res.model5", i), sim_res.model5)
assign(paste0("sim_res.model6", i), sim_res.model6)
assign(paste0("sim_res.model7", i), sim_res.model7)

}

# order Arthropoda/Annelida/Bacillariophyta/Cnidaria/Mollusca/Porifera/Rhodophyta
# adjust p-value, need to manually change dataset 
p.arms <- c(Anova(model1.arms)[3],Anova(model2.arms)[3],Anova(model3.arms)[3],
            Anova(model4.arms)[3],Anova(model5.arms)[3],Anova(model6.arms)[3],Anova(model7.arms)[3])
p.motil <- c(Anova(model1.mo)[3],Anova(model2.mo)[3],Anova(model3.mo)[3],
             Anova(model4.mo)[3],Anova(model5.mo)[3],Anova(model6.mo)[3],Anova(model7.mo)[3])
p.sessile <- c(Anova(model1.se)[3],Anova(model2.se)[3],Anova(model3.se)[3],
               Anova(model4.se)[3],Anova(model5.se)[3],Anova(model6.se)[3],Anova(model7.se)[3])

p.all <- c(as.numeric(p.arms),as.numeric(p.motil),as.numeric(p.sessile))

p_adjusted <- p.adjust(p.all, method = "BH")
p.all <- cbind(p.all,p_adjusted)
# write.csv(p.all,'adj.p.csv')


# check model 
# Check residuals
# order Arthropoda/Annelida/Bacillariophyta/Cnidaria/Mollusca/Porifera/Rhodophyta
qqnorm(residuals(model1.arms))
qqnorm(residuals(model2.arms))
qqnorm(residuals(model3.arms))
qqnorm(residuals(model4.arms))
qqnorm(residuals(model5.arms))
qqnorm(residuals(model6.arms))
qqnorm(residuals(model7.arms))

# all checked out 
# order Arthropoda/Annelida/Bacillariophyta/Cnidaria/Mollusca/Porifera/Rhodophyta
# neeed to manually get the model details for AMRS/motile/sessile fractions. 
chi.arms <- c(as.numeric(Anova(model1.arms)[1]),
         as.numeric(Anova(model2.arms)[1]),
         as.numeric(Anova(model3.arms)[1]),
         as.numeric(Anova(model4.arms)[1]),
         as.numeric(Anova(model5.arms)[1]),
         as.numeric(Anova(model6.arms)[1]),
         as.numeric(Anova(model7.arms)[1]))

p.ks.arms <- c(testResiduals(sim_res.model1.arms)$uniformity$p.value,
          testResiduals(sim_res.model2.arms)$uniformity$p.value,
          testResiduals(sim_res.model3.arms)$uniformity$p.value,
          testResiduals(sim_res.model4.arms)$uniformity$p.value,
          testResiduals(sim_res.model5.arms)$uniformity$p.value,
          testResiduals(sim_res.model6.arms)$uniformity$p.value,
          testResiduals(sim_res.model7.arms)$uniformity$p.value)
p.dis <- c(testResiduals(sim_res.model1.arms)$dispersion$p.value,
           testResiduals(sim_res.model2.arms)$dispersion$p.value,
           testResiduals(sim_res.model3.arms)$dispersion$p.value,
           testResiduals(sim_res.model4.arms)$dispersion$p.value,
           testResiduals(sim_res.model5.arms)$dispersion$p.value,
           testResiduals(sim_res.model6.arms)$dispersion$p.value,
           testResiduals(sim_res.model7.arms)$dispersion$p.value)

performance::check_singularity(model1)
performance::check_singularity(model2)
performance::check_singularity(model3)
performance::check_singularity(model4)
performance::check_singularity(model5)
performance::check_singularity(model6)
performance::check_singularity(model7)

# posthoc
# order Arthropoda/Annelida/Bacillariophyta/Cnidaria/Mollusca/Porifera/Rhodophyta
posthoc.1.arms
posthoc.2.arms
posthoc.3.arms
posthoc.4.arms
posthoc.5.arms
posthoc.6.arms
posthoc.7.arms

# Cndaria .4 
p.ks.arms <- c(testResiduals(sim_res.model4.arms)$uniformity$p.value,
               testResiduals(sim_res.model4.mo)$uniformity$p.value,
               testResiduals(sim_res.model4.se)$uniformity$p.value)
p.dis <- c(testResiduals(sim_res.model4.arms)$dispersion$p.value,
           testResiduals(sim_res.model4.mo)$dispersion$p.value,
           testResiduals(sim_res.model4.se)$dispersion$p.value)
posthoc.4.arms
posthoc.4.mo
posthoc.4.se

performance::check_singularity(model4.arms)
performance::check_singularity(model4.mo)
performance::check_singularity(model4.se)

chi.arms.4 <- c(as.numeric(Anova(model4.arms)[1]),
              as.numeric(Anova(model4.mo)[1]),
              as.numeric(Anova(model4.se)[1]))
