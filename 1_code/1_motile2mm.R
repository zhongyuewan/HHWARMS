## set WD 
# mac
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/')
rm(list=ls())


## load lib
library(dplyr)
library(ggplot2)
library(tidyr)
library(Biostrings)
library(biomformat)

## the color
myColor <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

## get data
metadata.d1 <- read.table('sample-metadata.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')[,1:9]
freqTable <- read.csv('motileData_clean.csv', header = TRUE)

# fix metadata/freqTable
metadata.d2 <- metadata.d1 %>% filter (fraction == "sessile")
metadeta <- metadata.d2

freqTable$ckleanPhyla <- gsub(" ","", freqTable$ckleanPhyla)
freqTable$cleanGenera <- gsub(" ","", freqTable$cleanGenera)
freqTable$ARMS <-  gsub("HH","HHW", freqTable$ARMS)
allPhylum <- unique(freqTable$ckleanPhyla)
allGenera <- unique(freqTable$cleanGenera)


# make a function to add things together 
countPhylum <- function (df){
  
  result <- df %>%
    group_by(ckleanPhyla) %>%
    summarise(Count = sum(Count)) %>%
    arrange(desc(Count)) 
  
  return(result)
  
}

# make new df and loop it 
ARMSn <- unique(metadeta$arms)
temp <- as.data.frame(countPhylum(freqTable %>% filter(ARMS==ARMSn[1])))
temp <- rbind(temp,c("all",sum(temp$Count)))


temp$ARMS <- ARMSn[1]
temp$depth <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$depth), nrow(temp))
temp$treatment <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$treatment), nrow(temp))
temp$coral <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$coral), nrow(temp))

countPhyla <- temp


for (i in ARMSn[-1]) {
  
  temp <- as.data.frame(countPhylum(freqTable[freqTable$ARMS==i,]))
  temp <- rbind(temp,c("all", sum(temp$Count)))
  
  temp$ARMS <- i
  temp$depth <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$depth), nrow(temp))
  temp$treatment <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$treatment), nrow(temp))
  temp$coral <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$coral), nrow(temp))
  
  countPhyla <- rbind(countPhyla,temp)
  
}

# fix some name 
names(countPhyla) <- c("phylum","count","ARMS","depth","treatment","coral")
countPhyla$count <- as.numeric(countPhyla$count)



countPhylaLOG <- countPhyla
countPhylaLOG$count <- log(countPhyla$count+1)

# do some plot (log or no-log are same same)
countPhylaLOG$coral <- factor(countPhylaLOG$coral, levels = c("acropora","pavona","platygyra","acroporaplatygyra","acroporapavonaplatygyra","no"))
plot.all <- ggplot(countPhylaLOG %>% filter(phylum!="all"),
                   aes(y=count, x=treatment, fill=depth))+
  geom_boxplot()+
  scale_fill_manual(values = rep(myColor[-1],100))+
  theme_classic()+
  labs(title = "all") +
  theme()


model.all <- aov(count~depth, countPhyla %>% filter(phylum=="all"))
summary(model.all)
TukeyHSD(model.all)

model.all2 <- aov(count~treatment, countPhyla %>% filter(phylum=="all"))
summary(model.all2)
TukeyHSD(model.all2)




# 1 Arthropoda 
plot.artho <- ggplot(countPhyla %>% filter(phylum=="Arthropoda") ,
                     aes(y=count, x=depth))+
  geom_boxplot()+
  theme_classic()+
  labs(title = "Arthropoda") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())

model.arth <- aov(count~depth, countPhyla %>% filter(phylum=="Arthropoda"))
summary(model.arth)
TukeyHSD(model.arth)


# 2 Mollusca
plot.moll <- ggplot(countPhyla %>% filter(phylum=="Mollusca") ,
                     aes(y=count, x=depth))+
  geom_boxplot()+
  theme_classic()+
  labs(title = "Mollusca") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())

model.moll <- aov(count~treatment, countPhyla %>% filter(phylum=="Mollusca"))
summary(model.moll)
TukeyHSD(model.moll)


# 3 Mollusca
plot.annel <- ggplot(countPhyla %>% filter(phylum=="Annelida") ,
                    aes(y=count, x=depth))+
  geom_boxplot()+
  theme_classic()+
  labs(title = "Annelida") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())

model.anne <- aov(count~treatment, countPhyla %>% filter(phylum=="Annelida"))
summary(model.anne)
TukeyHSD(model.anne)



##################
# break down in genus 


# make a function to add things together 
countGenus <- function (df){
  
  result <- df %>%
    group_by(cleanGenera) %>%
    summarise(Count = sum(Count)) %>%
    arrange(desc(Count)) %>% 
    mutate(rank=row_number())
  
  return(result)
  
}

# make new df and loop it 
temp <- as.data.frame(countGenus(freqTable %>% filter(ARMS==ARMSn[1])))



temp$ARMS <- ARMSn[1]
temp$depth <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$depth), nrow(temp))
temp$treatment <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$treatment), nrow(temp))
temp$coral <- rep(unique(freqTable %>% filter(ARMS==ARMSn[1]) %>% .$coral), nrow(temp))

dataGenus <- temp


for (i in ARMSn[-1]) {
  
  temp <- as.data.frame(countGenus(freqTable[freqTable$ARMS==i,]))
  
  temp$ARMS <- i
  temp$depth <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$depth), nrow(temp))
  temp$treatment <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$treatment), nrow(temp))
  temp$coral <- rep(unique(freqTable %>% filter(ARMS==i) %>% .$coral), nrow(temp))
  
  dataGenus <- rbind(dataGenus,temp)
  
}

# fix some name 
names(dataGenus) <- c("genus","count","rank","ARMS","depth","treatment","coral")
dataGenus$count <- as.numeric(dataGenus$count)

# make feaTable 
nARMS <- length(ARMSn)
nfea <- length(unique(freqTable$cleanGenera))
namefea <- unique(freqTable$cleanGenera)

feqTable <- as.data.frame(matrix(0, ncol=nARMS,,nrow=nfea))
names(feqTable) <- ARMSn
row.names(feqTable) <- namefea
nTotalfea <- nrow(dataGenus)


for (i in 1:nTotalfea) {
  tempG <- dataGenus[i,1]
  tempA <- dataGenus[i,4]
  tempC <- dataGenus[i,2]
  feqTable[row.names(feqTable)==tempG,][[tempA]] <- tempC

}


# do some simply plots 
# by phylum -> countPhyla
# by genus -> dataGenus

countPhylaPlot <- countPhyla%>% filter(phylum!="all")%>% mutate(phylum=reorder(phylum, -count,sum))
plotP <- ggplot(countPhylaPlot, aes(x = ARMS, y = count, fill = phylum)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = myColor[-1])
plotPbyCoral <- ggplot(countPhylaPlot, aes(x = coral, y = count, fill = phylum)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = myColor[-1])

dataGenusPlot <- dataGenus%>% mutate(genus=reorder(genus, -count,sum))
plotG <- ggplot(dataGenusPlot, aes(x = ARMS, y = count, fill = genus)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = rep(myColor[-1],100))
plotGbyDepth <- ggplot(dataGenusPlot, aes(x = depth, y = count, fill = genus)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = rep(myColor[-1],100))
plotGbyTreat <- ggplot(dataGenusPlot, aes(x = treatment, y = count, fill = genus)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = rep(myColor[-1],100))
plotGbyCoral <- ggplot(dataGenusPlot, aes(x = coral, y = count, fill = genus)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values = rep(myColor[-1],100))

## some basic stats to see if anything significant 
model1 <- lm(count~depth, countPhyla%>%filter(phylum=="all"))
summary(model1)

model2 <- aov(count~treatment, countPhyla%>%filter(phylum=="all"))
summary(model2)

plot.all.coral <- ggplot(countPhyla,
                     aes(y=count, x=phylum, fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values = rep(myColor[-1],100))+
  theme_classic()+
  labs(title = "all") +
  theme()
  
plot.all.deep <- ggplot(countPhyla,
                         aes(y=count, x=phylum, fill=depth))+
  geom_boxplot()+
  scale_fill_manual(values = rep(myColor[-1],100))+
  theme_classic()+
  labs(title = "all") +
  theme()


# write data out# write data outARMS
# dont over write 
write.csv(feqTable,# "9-2mmdata/2mmFeqTable.csv")
)
