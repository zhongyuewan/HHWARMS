## set WD 
# mac
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/3-deconPostaaClus97/BOLD')
rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(vegan)

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

# make data to plot 
alphaTable.PLOT <- as.data.frame(rbind(t(alphaTableArms)[2:19,],t(alphaTableMol)[2:19,],t(alphaTableSes)[2:19,]))
alphaTable.PLOT$ARMS <- rep(names(alphaTableArms)[2:19],3)
names(alphaTable.PLOT)[1:3] <- c("total","unique","per")
alphaTable.PLOT$fraction <- c(rep("ARMS",18),rep("mol",18),rep("ses",18))
alphaTable.PLOT$depth <- c(rep(metadate$depth,3))
alphaTable.PLOT$treatment <- c(rep(metadate$treatment1,3))
alphaTable.PLOT$coral <- c(rep(metadate$coral,3))
alphaTable.PLOT$total <- as.numeric(alphaTable.PLOT$total)
alphaTable.PLOT$unique <- as.numeric(alphaTable.PLOT$unique)
alphaTable.PLOT$per <- as.numeric(alphaTable.PLOT$per)


alphaTable.PLOT$coral <- factor(alphaTable.PLOT$coral, levels=unique(alphaTable.PLOT$coral)[c(3,2,1,4,6,5)])
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
str(alphaTable.PLOT)
model1 <- aov(total~depth, data=alphaTable.PLOT%>%filter(fraction=="ses"))
summary(model1)

model2 <- aov(total~treatment, data=alphaTable.PLOT%>%filter(fraction=="ARMS"))
summary(model2)
TukeyHSD(model2)

model3 <- lm(total~treatment, data=alphaTable.PLOT%>%filter(fraction=="ses"))
summary(model3)

# with factor 
alphaTable.PLOT$depth <- factor(alphaTable.PLOT$depth, levels = c('shallow','middle','deep'))
alphaTable.PLOT$treatment <- factor(alphaTable.PLOT$treatment, levels = c('control','mono','mix','poly'))

model1.f <- aov(total~depth, data=alphaTable.PLOT%>%filter(fraction=="ses"))
summary(model1)
summary(model1.f)

model2.f <- aov(total~treatment, data=alphaTable.PLOT%>%filter(fraction=="ARMS"))
summary(model2)
summary(model2.f)
TukeyHSD(model2)

model3.f <- lm(total~treatment, data=alphaTable.PLOT%>%filter(fraction=="ses"))
summary(model3)
summary(model3.f)

# write the alpha.table.plot out 
# be careful dont overwrite files 
write.csv(alphaTable.PLOT, #"deconPostaaClus97/byARMS/alphaTable.plot.csv")
)
          
## do one with abundance to compare the sessie 
alpha.plot.abun.d1 <- as.data.frame(matrix(0,nrow=18, ncol=2)) 
alpha.plot.abun.d1[,1] <- names(freqTable)
alpha.plot.abun.d1[,2] <- colSums(sqrt(freqTable))
names(alpha.plot.abun.d1) <- c("ARMS","sqrtRead")
alpha.plot.abun.d1$treatment <- metadate$treatment
alpha.plot.abun.d1$coral <-metadate$coral
alpha.plot.abun.d1$depth <-metadate$depth

plotuniOTUs.treatment.ab <- ggplot(alpha.plot.abun.d1, aes(x=treatment, y=sqrtRead,fill=treatment)) +
  geom_boxplot() +
  xlab("Treatments") + 
  ylab("Unique OTUs") +
  theme_classic()

model4 <- aov(sqrtRead~treatment, data=alpha.plot.abun.d1)
summary(model4) 
TukeyHSD(model4)
# so it's the same pattern with 

