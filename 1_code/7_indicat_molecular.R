## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(ggrepel)
library(indicspecies)

## set WD 
#mac 
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/3-deconPostaaClus97/BOLD')
rm(list=ls())


myColor <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

## get data
metadata <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')[,-2]
# freqTable.d1 <- read.csv('freqTableCleanByMol.csv', header = TRUE, row.names = 1)
# freqTable.d1 <- read.csv('freqTableCleanBySes.csv', header = TRUE, row.names = 1)
freqTable.d1 <- read.csv('freqTableCleanByArms.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTable <- read.csv('alphaTable.byArms.csv', header = TRUE)
taxTabSheClean <- read.csv('../../4-TaxAssign/TaxAsn_shelbyOmidori.csv', header=T, row.names = 1, check.names = F)


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
freqTable.d2 <- sqrt(freqTable.d1[rowSums(freqTable.d1)>0,])
freqTable.d2$taxa <- "greenbie"

for(i in 1: nrow(freqTable.d2)){
  id <- rownames(freqTable.d2)[i]
  taxa <-ifelse(id%in%row.names(taxTabSheClean),taxTabSheClean[row.names(taxTabSheClean)==id,]$phylum,"unassigned") # change the family here
  freqTable.d2$taxa[i] <- taxa
  
}
  
freqTable.d3 <- freqTable.d2%>% filter(taxa!="unassigned" & taxa!="family" & taxa!="<NA>")
freqTable.d3[,1:18] <- sqrt(freqTable.d3[,1:18])
freqTable.d4 <- as.data.frame(freqTable.d3%>% group_by(taxa) %>% summarise(across(1:18,sum)))    
freqTable <- freqTable.d4[,2:19]
row.names(freqTable) <- freqTable.d4[,1]



##### Similarity Percentage (SIMPER) Analysis: which item is different how? 
  ## first run it by OTUs (might not work well and will run it by phylum later)
set.seed(123)
  simper_result <- simper(t(freqTable), group = metadata$depth)
  temp.d1 <- summary(simper_result) # this is good, i can find something in it now. 
  
  temp <- as.data.frame(temp.d1[1])%>% filter(is.na(.[,7])==F)
  globD <- temp[1,1]/temp[1,6]
  dataSIM <- as.data.frame(matrix(0, nrow = nrow(temp)+1, ncol=10))
  dataSIM[1:nrow(temp),3:9] <- temp
  dataSIM[1:nrow(temp),1] <- row.names(temp)
  dataSIM[nrow(temp)+1,1] <- "others"
  dataSIM[,2] <- names(temp.d1[1])
  dataSIM[nrow(temp)+1,c(8,10)] <- 1-dataSIM[nrow(temp),8]
  dataSIM[1:nrow(temp),10] <- dataSIM[1:nrow(temp),3]/globD
  
  names(dataSIM) <- c("taxa","Tpair","average","sd","ration","ava","avb","cumsum","p","perD")

set.seed(123)  
  for (i in 2:3){
    temp <- as.data.frame(temp.d1[i])%>% filter(is.na(.[,7])==F)
    globD <- temp[1,1]/temp[1,6]
    tempdf <- as.data.frame(matrix(0, nrow = nrow(temp)+1, ncol=10))
    tempdf[1:nrow(temp),3:9] <- temp
    tempdf[1:nrow(temp),1] <- row.names(temp)
    tempdf[nrow(temp)+1,1] <- "others"
    tempdf[,2] <- names(temp.d1[i])
    tempdf[nrow(temp)+1,c(8,10)] <- 1-tempdf[nrow(temp),8]
    tempdf[1:nrow(temp),10] <- tempdf[1:nrow(temp),3]/globD
    
    names(tempdf) <- c("taxa","Tpair","average","sd","ration","ava","avb","cumsum","p","perD")
    
    dataSIM <- rbind(dataSIM,tempdf)
  }

dataSIM.sig <- dataSIM%>% filter(.[,9]<0.05 & .[,1]!='others') 
dataSIM.sig.10 <- dataSIM.sig %>% filter(.[,10]>0.02)

  # make pie chart plot 
  # 1 rerank genus 
  dataSIM_plot <- as.data.frame(dataSIM %>%
    group_by(taxa) %>%               # Group by column 'a'
    summarise(total_b = sum(perD)) %>%  # Sum 'b' within groups
    arrange(desc(total_b)) %>%    # Rank descending
    left_join(dataSIM, by = "taxa"))       # Merge back with original data (optional)
  unique(dataSIM_plot$taxa)
  dataSIM_plot$taxa <- factor(dataSIM_plot$taxa, levels = c(unique(dataSIM_plot$taxa)[-4],"others"))
  
  # 2 make the pie 
  plotPie <- ggplot(dataSIM_plot, aes(x="", y=perD, fill=taxa)) +
    geom_bar(stat="identity", width=1) +
    scale_fill_manual(values=c(rep(myColor[-1],100)))+ 
    coord_polar("y", start=0)+
   # theme(legend.position = "none") +
    facet_wrap(~ interaction(Tpair)) 
  
  # anything significant? 
  dataSIM %>% filter(p>0 & p<0.05)%>%arrange(desc(Tpair))
  
  sum(dataSIM %>% filter(p>0 & p<0.05)%>% filter(Tpair=="control_poly") %>% .$perD)
  

##### # Indicator Species Analysis (ISA)##### # Indicator Species Analysissum() (ISA)
  isa <- multipatt(t(freqTable), metadata$depth, func = "IndVal", control = how(nperm = 999))
  summary(isa, alpha = .1) # Relaxed threshold
    
  
