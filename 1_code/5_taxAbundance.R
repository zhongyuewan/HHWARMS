## set WD 
#mac 
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/3-deconPostaaClus97/BOLD')
rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(ggrepel)
library(tidyr)


## get data
metadate <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)



# check out sharing OTUs by time 
shaARMS <- metadate[(metadate$depth=="shallow"),]$arms
midTARMS <- metadate[(metadate$depth=="middle"),]$arms
depARMS <- metadate[(metadate$depth=="deep"),]$arms

conARMS <- metadate[(metadate$treatment1=="control"),]$arms
monoARMS <- metadate[(metadate$treatment1=="mono"),]$arms
monoARMS1 <- metadate[(metadate$coral=="acropora"),]$arms
monoARMS2 <- metadate[(metadate$coral=="pavona"),]$arms
monoARMS3 <- metadate[(metadate$coral=="platygyra"),]$arms
polyARMS <- metadate[(metadate$treatment1=="poly"),]$arms
mixARMS <- metadate[(metadate$treatment1=="mix"),]$arms

AcroARMS <- metadate[grepl("acropora", metadate$coral),]$arms
PavoARMS <- metadate[grepl("pavona", metadate$coral),]$arms
PlatARMS <- metadate[grepl("platygyra", metadate$coral),]$arms

# calculate datamatrix based on squaroot 
# twist the freqTable to switch among motile/sessile/ARMS

# freqTable.d2 <- freqTable.d1%>%filter(fraction=="motile")
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="sessile")
freqTable.d2 <- freqTable.d1%>%filter(fraction=="byARMS")
freqTable <- freqTable.d2[,3:21]
row.names(freqTable) <- freqTable.d2$id

freqTable_ranked_unique <- freqTable %>%
  count(phylum, name = "total") %>%
  mutate(rank = dense_rank(dplyr::desc(total))) %>%
  arrange(rank)

freqTableAll <- freqTable
freqTableArthropoda <- freqTable%>%filter(phylum=="Arthropoda")
freqTableBacillariophyta <- freqTable%>%filter(phylum=="Bacillariophyta")
freqTableAnnelida <- freqTable%>%filter(phylum=="Annelida")
freqTableMollusca <- freqTable%>%filter(phylum=="Mollusca")
freqTableRhodophyta <- freqTable%>%filter(phylum=="Rhodophyta")
freqTableCnidaria <- freqTable%>%filter(phylum=="Cnidaria")
freqTablePorifera <- freqTable%>%filter(phylum=="Porifera")


# get the reads 
allREAD <- row.names(freqTable) # 1842
conREAD <- row.names(freqTable[rowSums(freqTable[,conARMS])>0,]) # 621
READ1 <- row.names(freqTable[rowSums(freqTable[,monoARMS1])>0,]) # 711
READ2 <- row.names(freqTable[rowSums(freqTable[,monoARMS2])>0,]) # 812
READ3 <- row.names(freqTable[rowSums(freqTable[,monoARMS3])>0,]) # 689
READ4 <- row.names(freqTable[rowSums(freqTable[,mixARMS])>0,]) # 894
READ5 <- row.names(freqTable[rowSums(freqTable[,polyARMS])>0,]) # 735

# build a df for the sidebar
dataSB <- as.data.frame(matrix(0, nrow=40,ncol=9))
names(dataSB) <- c("treatment", "phylum", "conCOUNT", "count","share","gain","loss","pergain","perloss")
dataSB$phylum <- c(rep("All",5),rep("Arthropoda",5),rep("Annelida",5),rep("Bacillariophyta",5),
                   rep("Mollusca",5),rep("Rhodophyta",5),rep("Cnidaria",5),rep("Porifera",5)) 
dataSB$treatment <- rep(c("mono1", "mono2", "mono3", "mix", "poly"),8)
taxaRank <- dataSB$phylum[1:7]


freqTableCon <- freqTable[conREAD,]%>%
  count(phylum, name = "total") %>% 
  mutate(rank = dense_rank(dplyr::desc(total))) %>%
  arrange(rank)

conCount <- c(sum(freqTableCon$total),
              freqTableCon%>% filter(phylum=="Arthropoda")%>%.$total,
              freqTableCon%>% filter(phylum=="Annelida")%>%.$total,
              freqTableCon%>% filter(phylum=="Bacillariophyta")%>%.$total,
              freqTableCon%>% filter(phylum=="Mollusca")%>%.$total,
              freqTableCon%>% filter(phylum=="Rhodophyta")%>%.$total,
              freqTableCon%>% filter(phylum=="Cnidaria")%>%.$total,
              freqTableCon%>% filter(phylum=="Porifera")%>%.$total)
                
dataSB$conCOUNT <- c(rep(conCount[1],5),rep(conCount[2],5),rep(conCount[3],5),rep(conCount[4],5),
                     rep(conCount[5],5),rep(conCount[6],5),rep(conCount[7],5),rep(conCount[8],5))


for (k in c("All","Arthropoda","Annelida","Bacillariophyta","Cnidaria","Rhodophyta","Mollusca","Porifera")) {
  FT <- get(paste0("freqTable",k))
  READs.d1 <- row.names(FT)
  
  for (i in 1:5) {
    READs.d2 <- get(paste0("READ",i))
    READs <- READs.d1[READs.d1 %in% READs.d2]
    READsNO <- length(READs)
    shareRead <- intersect(READs,conREAD)
    shareReadNO <- length(shareRead)
    
    # total reads 
    dataSB[dataSB$phylum==k,][i,4] <- READsNO
  
    # share reads 
    dataSB[dataSB$phylum==k,][i,5] <- shareReadNO
  
    # gain / loss / pergain / perloss
    dataSB[dataSB$phylum==k,][i,6] <- READsNO - shareReadNO
    dataSB[dataSB$phylum==k,][i,7] <- dataSB[dataSB$phylum==k,][i,3] - shareReadNO
    dataSB[dataSB$phylum==k,][i,8] <- round(dataSB[dataSB$phylum==k,][i,6]/dataSB[dataSB$phylum==k,][i,3] * 100, 1)
    dataSB[dataSB$phylum==k,][i,9] <- round(dataSB[dataSB$phylum==k,][i,7]/dataSB[dataSB$phylum==k,][i,3] * 100, 1)
    
  
    }
  
  
  
  
}


# make long table 
# Pivot gain and loss into long format
dataSB[,c(7,9)] <- dataSB[,c(7,9)]*-1
data_long<- as.data.frame(pivot_longer(dataSB, cols = c(pergain, perloss), names_to = "GnL", 
             values_to = "PERchange"))

data_long$treatment <- factor(data_long$treatment, levels = c("mono1","mono2","mono3","mix","poly"))
data_long$phylum <- factor(data_long$phylum, levels = c("All","Porifera","Cnidaria","Rhodophyta",
                                                        "Mollusca","Bacillariophyta","Annelida","Arthropoda"))
# make sideway plots 
sidebarALLbyALL_GLS <- ggplot(data_long, aes(x = PERchange, y = phylum, fill = GnL)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("pergain" = "#79bfe7", "perloss" = "#D55E00")) +
  geom_text(data = subset(data_long%>%filter(GnL == "pergain")), aes(label = PERchange), 
            position = position_stack(vjust = 0), hjust=-1, color = "black", size = 3) +  # Center "gain" labels at the top
  geom_text(data = subset(data_long%>%filter(GnL == "perloss")), aes(label = PERchange), 
            position = position_stack(vjust = 0), hjust=-1,color = "yellow", size = 3) +  # Center "gain" labels at the top
  
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")+
  labs(x = NULL, y = NULL) +
  facet_wrap(~treatment)


## run chi-square 
# build a df for the chi square 
data.chi <- as.data.frame(matrix(0, nrow=8, ncol=6))
row.names(data.chi) <- c("All","Arthropoda","Annelida","Bacillariophyta","Cnidaria","Rhodophyta","Mollusca","Porifera")
names(data.chi) <- c("control","mono1","mono2","mono3","mix","poly")

data.chi[1,1] <- length(conREAD)
data.chi[2,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[2]))
data.chi[3,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[3]))
data.chi[4,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[4]))
data.chi[5,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[5]))
data.chi[6,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[6]))
data.chi[7,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[7]))
data.chi[8,1] <- nrow(freqTable[conREAD,] %>% filter(phylum==row.names(data.chi)[8]))

for (i in 1:5) {
  READs <- get(paste0("READ",i))
  data.chi[1,i+1] <- length(READs) 
  data.chi[2,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[2]))
  data.chi[3,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[3]))
  data.chi[4,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[4]))
  data.chi[5,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[5]))
  data.chi[6,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[6]))
  data.chi[7,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[7]))
  data.chi[8,i+1] <- nrow(freqTable[READs,] %>% filter(phylum==row.names(data.chi)[8]))
}


## 1 chi model on all richness 
model1 <- chisq.test(data.chi[1,])
model1 # ALL control has the lowest richness and  
model1$stdres
               

# 2 chi model on gain and loss
data.chi.2 <- data_long[1:nrow(data_long) %% 2 == 0,c(1,2,6,7)]
data.chi.2$loss <- data.chi.2$loss * -1
      
data.chi.gain <- as.data.frame(data.chi.2[,1:3] %>%      
  pivot_wider(
    names_from = treatment,    # Column to get new column names from
    values_from = gain         # Column to get values from
  ))
                            
data.chi.loss <- as.data.frame(data.chi.2[,c(1,2,4)] %>%      
                                 pivot_wider(
                                   names_from = treatment,    # Column to get new column names from
                                   values_from = loss         # Column to get values from
                                 ))

row.names(data.chi.gain) <- data.chi.gain[,1]
row.names(data.chi.loss) <- data.chi.loss[,1]

data.chi.gain.f <- data.chi.gain[,2:6]
data.chi.loss.f <- data.chi.loss[,2:6]

modelG <- chisq.test(data.chi.gain.f[1,])
modelG$stdres # mono 2 gained the least and mix gained the most 
modelL <- chisq.test(data.chi.loss.f[1,])
modelL$stdres # all lost same same 

modelg <- chisq.test(data.chi.gain.f[2:6,])
modelg$stdres # all phylum same same gain 
modell <- chisq.test(data.chi.loss.f[2:6,])
modell$stdres # all phylum same same loss  


