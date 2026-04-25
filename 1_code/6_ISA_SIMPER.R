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

## the color
myColor <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

## get data
metadata <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTable <- read.csv('alphaTable.byArms.csv', header = TRUE)
SOM <- read.csv('../../4-TaxAssign/TaxAsn_shelbyOmidori.csv', header = TRUE)

# check out sharing OTUs by time 
metadata$yescoral <- "yes"
metadata[metadata$coral=="no",]$yescoral <- "no"

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
freqTable.d2 <- freqTable.d1%>%filter(fraction=="sessile")
freqTable.d2[,4:21] <- sqrt(freqTable.d2[,4:21]) 
freqTable.d3 <- freqTable.d2 %>% filter(phylum!="notAssigned")

# sqrt or not? 
freqTable.d4 <- freqTable.d3[,3:21]
freqTable.d4[,2:19] <- sqrt(freqTable.d3[,4:21])

feqTable <- as.data.frame(freqTable.d4 %>%
  group_by(phylum)%>% 
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE)
  ))

row.names(feqTable) <- feqTable[,1]
feqTable <- feqTable[2:19]


#### Similarity Percentage (SIMPER) Analysis
simper_result <- simper(t(feqTable), group = metadata$treatment)  #change between depth and treatment 
temp.d1 <- summary(simper_result) # this is good, i can find something in it now. 

temp <- as.data.frame(temp.d1[1])%>% filter(.[,6]<=0.8)
globD <- temp[1,1]/temp[1,6]
dataSIM <- as.data.frame(matrix(0, nrow = nrow(temp)+1, ncol=10))
dataSIM[1:nrow(temp),3:9] <- temp
dataSIM[1:nrow(temp),1] <- row.names(temp)
dataSIM[nrow(temp)+1,1] <- "others"
dataSIM[,2] <- names(temp.d1[1])
dataSIM[nrow(temp)+1,c(8,10)] <- 1-dataSIM[nrow(temp),8]
dataSIM[1:nrow(temp),10] <- dataSIM[1:nrow(temp),3]/globD

names(dataSIM) <- c("genus","Tpair","average","sd","ration","ava","avb","cumsum","p","perD")

for (i in 2:5){
  temp <- as.data.frame(temp.d1[i])%>% filter(.[,6]<=0.8)
  globD <- temp[1,1]/temp[1,6]
  tempdf <- as.data.frame(matrix(0, nrow = nrow(temp)+1, ncol=10))
  tempdf[1:nrow(temp),3:9] <- temp
  tempdf[1:nrow(temp),1] <- row.names(temp)
  tempdf[nrow(temp)+1,1] <- "others"
  tempdf[,2] <- names(temp.d1[i])
  tempdf[nrow(temp)+1,c(8,10)] <- 1-tempdf[nrow(temp),8]
  tempdf[1:nrow(temp),10] <- tempdf[1:nrow(temp),3]/globD
  
  names(tempdf) <- c("genus","Tpair","average","sd","ration","ava","avb","cumsum","p","perD")
  
  dataSIM <- rbind(dataSIM,tempdf)
  
  
}

dataSIM%>% filter(p<0.1 & p!= 0)




dataSIM%>% filter(p<0.1 & genus != "others")


# make pie chart plot 
# 1 rerank genus 
dataSIM_plot <- as.data.frame(dataSIM %>%
                                group_by(genus) %>%               # Group by column 'a'
                                summarise(total_b = sum(perD)) %>%  # Sum 'b' within groups
                                arrange(desc(total_b)) %>%    # Rank descending
                                left_join(dataSIM, by = "genus"))       # Merge back with original data (optional)
unique(dataSIM_plot$genus)
dataSIM_plot$genus <- factor(dataSIM_plot$genus, levels = c(unique(dataSIM_plot$genus)[-1],"others"))

# 2 make the pie 
plotPie <- ggplot(dataSIM_plot, aes(x="", y=perD, fill=genus)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c(rep(myColor[-1],100)))+ 
  facet_wrap(~ interaction(Tpair)) 


dataSIM_plot%>% filter(p>0 & p<0.05)%>%arrange(desc(Tpair))

# Indicator Species Analysis (ISA)
isa <- multipatt(t(feqTable), metadata$treatment, func = "IndVal", control = how(nperm = 999))
summary(isa, alpha = .1) # Relaxed threshold



### random drafts 
molluRows.d1 <- freqTable.d2%>% filter(phylum == "Mollusca") 
molluRows.d2 <- molluRows.d1[rowSums(molluRows.d1[,conARMS])>0,][,c("id","phylum",conARMS)] # conARMS mollus row

NconARMS <- c(monoARMS,polyARMS,mixARMS)
molluRows.d3 <- molluRows.d1[rowSums(molluRows.d1[,NconARMS])>0,][,c("id","phylum",NconARMS)] # non-conARMS mollus row

molluRows.d2$sum <- rowSums(molluRows.d2[,3:5]) 
molluRows.d3$sum <- rowSums(molluRows.d3[,3:17]) 

# one or two of the following 
molluRows.c <- molluRows.d2 %>%    arrange(desc(sum))
molluRows.nc <- molluRows.d3 %>%    arrange(desc(sum))

oysterID.c <- SOM[SOM$id %in% molluRows.c$id,] %>% filter(family=="Ostreidae") %>%.$id
oysterID.nc <- SOM[SOM$id %in% molluRows.nc$id,] %>% filter(family=="Ostreidae") %>%.$id
OysterID <- union(oysterID.c,oysterID.nc)

sum(molluRows.c[molluRows.c$id %in% oysterID.c,]$sum)/sum(molluRows.c$sum) # 62.21 % abundance accounted by oyster 
sum(molluRows.nc[molluRows.nc$id %in% oysterID.nc,]$sum)/sum(molluRows.nc$sum) # 44.8% % abundance accounted by oyster 

freqTable.d2[freqTable.d2$id %in% OysterID,]

SOM %>% filter(id %in% oysterID.c)

### 
# Cnidaria 
SOM%>%filter(phylum=="Cnidaria")
hydro.rows <- (SOM%>%filter(class=="Hydrozoa"))[,1] # 81 Hydrozoa rows 
cnida.rows <- (SOM%>%filter(phylum=="Cnidaria"))[,1] # 100 Cnidaria rows, 81% cnidaria is hydrozoa

sum(freqTable.d2[freqTable.d2$id %in% hydro.rows,4:21]) # 3331.006 Hydrozoa
sum(freqTable.d2[freqTable.d2$id %in% cnida.rows,4:21]) # 3782.646 Hydrozoa
# 88.06 abundance is hydrozoa 

sum(freqTable.d2[freqTable.d2$id %in% chor.rows,4:21]) # 3782.646 Hydrozoa

# Chordata
chor.rows <- (SOM%>%filter(phylum=="Chordata"))[,1] # 32 chordata rows 
unique((SOM%>%filter(phylum=="Chordata"))$class) # fish/ray-fish/sharks + Ascidiacea 

fish.rows <- (SOM%>%filter(class=="Actinopteri"))[,1] # 9 fish rows  
rfish.rows <- (SOM%>%filter(class=="Actinopterygii"))[,1] # 4 ray-fish rows, all sessile Gobiidae
shark.rows <- (SOM%>%filter(class=="Chondrichthyes"))[,1] # 1 shark rows 
asci.rows <- (SOM%>%filter(class=="Ascidiacea"))[,1] # 16 rows 

ALLfish.row <- c(fish.rows, rfish.rows,shark.rows)

sum(freqTable.d2[freqTable.d2$id %in% fish.rows,4:10])/7 # sha 17.37202
sum(freqTable.d2[freqTable.d2$id %in% fish.rows,11:16])/6 # mid 15.75296
sum(freqTable.d2[freqTable.d2$id %in% fish.rows,17:21])/5 # deep 71.48001

sum(freqTable.d2[freqTable.d2$id %in% rfish.rows,4:10])/7 # 29.66415
sum(freqTable.d2[freqTable.d2$id %in% rfish.rows,11:16])/6 # 39.93157
sum(freqTable.d2[freqTable.d2$id %in% rfish.rows,17:21])/5 # 57.5759

sum(freqTable.d2[freqTable.d2$id %in% shark.rows,4:10])/7 # 0
sum(freqTable.d2[freqTable.d2$id %in% shark.rows,11:16])/6 # 0
sum(freqTable.d2[freqTable.d2$id %in% shark.rows,17:21])/5 # 1.182718

sum(freqTable.d2[freqTable.d2$id %in% ALLfish.row,4:10])/7 # 47.03617
sum(freqTable.d2[freqTable.d2$id %in% ALLfish.row,11:16])/6 # 55.68452
sum(freqTable.d2[freqTable.d2$id %in% ALLfish.row,17:21])/5 # 130.2386

sum(freqTable.d2[freqTable.d2$id %in% asci.rows,4:10])/7 # 27.99615
sum(freqTable.d2[freqTable.d2$id %in% asci.rows,11:16])/6 # 40.77651
sum(freqTable.d2[freqTable.d2$id %in% asci.rows,17:21])/5 # 49.80616
# both fish and Ascidiacea are more abundance further away!!! 

# Annelida
anne.rows <- (SOM%>%filter(phylum=="Annelida"))[,1] # 218 Annelida rows 
unique((SOM%>%filter(phylum=="Annelida"))$class) # Polychaeta/Clitellata/Sipuncula
nrow(SOM%>%filter(class=="Polychaeta")) # 191 Polychaeta 87.6% Annelida rows are polychaeta 
nrow(SOM%>%filter(class=="Clitellata")) # 15 Clitellata
nrow(SOM%>%filter(class=="Sipuncula")) # 3 Sipuncula

polyc.rows <- (SOM%>%filter(class=="Polychaeta"))[,1] 
clit.rows <- (SOM%>%filter(class=="Clitellata"))[,1] 
sipu.rows <- (SOM%>%filter(class=="Sipuncula"))[,1] 

sum(freqTable.d2[freqTable.d2$id %in% polyc.rows,4:10])/7 # sha 1147.536
sum(freqTable.d2[freqTable.d2$id %in% polyc.rows,11:16])/6 # mid 1158.557
sum(freqTable.d2[freqTable.d2$id %in% polyc.rows,17:21])/5 # deep 1227.112

sum(freqTable.d2[freqTable.d2$id %in% clit.rows,4:10])/7 # 32.06948
sum(freqTable.d2[freqTable.d2$id %in% clit.rows,11:16])/6 # 25.28199
sum(freqTable.d2[freqTable.d2$id %in% clit.rows,17:21])/5 # 54.60281

sum(freqTable.d2[freqTable.d2$id %in% sipu.rows,4:10])/7 # 54.57109
sum(freqTable.d2[freqTable.d2$id %in% sipu.rows,11:16])/6 # 51.0104
sum(freqTable.d2[freqTable.d2$id %in% sipu.rows,17:21])/5 # 54.33528


##
SOM%>%filter(family=="Ostreidae")
freqTable.d2 %>% filter(id %in% fish.rows)


