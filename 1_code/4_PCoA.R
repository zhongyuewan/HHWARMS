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
library(pairwiseAdonis)
library(cowplot)


## get data
metadate <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTable <- read.csv('alphaTable.byArms.csv', header = TRUE)
countTable <- read.csv('../../9-2mmdata/2mmFeqTable.csv', header = TRUE)


# check out sharing OTUs by time 
shaARMS <- metadate[(metadate$depth=="shallow"),]$arms
midTARMS <- metadate[(metadate$depth=="middle"),]$arms
depARMS <- metadate[(metadate$depth=="deep"),]$arms

conARMS <- metadate[(metadate$treatment1=="control"),]$arms
monoARMS <- metadate[(metadate$treatment1=="mono"),]$arms
polyARMS <- metadate[(metadate$treatment1=="poly"),]$arms
mixARMS <- metadate[(metadate$treatment1=="mix"),]$arms

AcroARMS <- metadate[grepl("acropora", metadate$coral),]$arms
PavoARMS <- metadate[grepl("pavona", metadate$coral),]$arms
PlatARMS <- metadate[grepl("platygyra", metadate$coral),]$arms

# calculate datamatrix based on squaroot 
# twist the freqTable to switch among motile/sessile/ARMS
freqTable.d2 <- freqTable.d1%>%filter(fraction=="byARMS") # for plot 4
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="motile") # for plot 2
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="sessile") # for plot 3

freqTable <- freqTable.d2[,4:21]
# freqTable <- countTable[,-1] # for plot 1

# for count data


# twist the following distM to switch between bray/jaccard
# distM <- vegdist(t(freqTable), method = "bray")
distM <- vegdist(t(sqrt(freqTable)), method = "bray")
# distM <- vegdist(t(freqTable), method = "bray") # for count data
# distM <- vegdist(t(freqTable), method = "jaccard")

# plot all the baseline by square root bray curtis 

  dist_matrix <- prcomp(distM)
  evalue <- fviz_eig(dist_matrix, addlabels = TRUE) #40.7, 17.4, 10.2
  dataPlot <- data.frame(matrix(nrow=18, ncol=4))
  row.names(dataPlot) <- names(freqTable)
  colnames(dataPlot) <- c("PC1", "PC2", "treatment", "water")
  
  dataPlot$PC1 <- dist_matrix$x[,1]
  dataPlot$PC2 <- dist_matrix$x[,2]
  dataPlot$treatment <- metadate$treatment
  dataPlot$water <- metadate$depth
  dataPlot$coral <- metadate$coral
  
  baseline1 <- ggplot(dataPlot, 
         aes(x = PC1, 
             y = PC2, 
             color = water,
             shape = water)) +
    geom_point() +
    geom_point(size=8) +
    stat_ellipse(type = "norm", aes(group = water), level = 0.95, lwd = 1)+
    xlab(evalue$data[1,2]) + 
    ylab(evalue$data[2,2]) +
    theme(text = element_text(size = 35)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())+ 
    #geom_text_repel(
    #  aes(label = row.names(dataPlot)),  
    #  size = 4,                    # Adjust label size
    #  max.overlaps = Inf,          # Show all labels (even if overlapping)
    #  box.padding = 0.5,           # Spacing between labels and points
    #  show.legend = FALSE          # Hide label legend
    #)+
    theme_classic()+
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    scale_shape_manual(values = c(16, 17, 15, 18, 23, 25)) + # Different shapes for treatments
    theme(legend.position = "none")

  combined_pieplot <- plot_grid(baseline1, baseline2, 
                                baseline3,baseline4,ncol = 4)

  
# adonis 
# it's tricky to run adonis with interaction because i dont have all treatments in all depth 
# so it violates the Homogeneity assumption.... 
  
  set.seed(123)
  model1 <- adonis2(distM~treatment1, 
                    data=metadate,
                    #strata = metadate$depth,	
                    by = "terms")
  
  # controlling depth as the strata 
  # by ARMS: by treatment, p=0.006, by coral, p=0.165
  # by motile: by treatment, p=0.047, by coral, p=0.165
  # by sessile: by treatment, p=0.013, by coral, p=0.028
  
  
  pairwise_results_coral <- pairwise.adonis2(
    distM ~ treatment1, 
  #  strata = metadate$treatment1,
    data = metadate   # Add strata if needed (e.g., strata = metadate$depth)
  )
  # after 

  model2.1 <- vegan::betadisper(distM, metadate$coral)  
  anova(model2.1) # homogeneous 
  
  model2.2 <- vegan::betadisper(distM, metadate$depth)  
  anova(model2.2) # homogeneous 
  
  # Pairwise tests for 'coral'  
  set.seed(123)
  pairwise_coral <- pairwise.adonis2(
    distM ~ coral, 
    data = metadate, 
    nperm = 999
  )
  
  # Pairwise tests for 'depth'
  set.seed(123)
  pairwise_depth <- pairwise.adonis2(
    distM ~ depth, 
    data = metadate, 
    nperm = 999
  )
  
  # motile
  # shallow_vs_middle p = 0.004
  # shallow_vs_deep p = 0.002
  # middle_vs_deep p = 0.07
  
  
  # sessile 
  # shallow_vs_middle p = 0.156
  # shallow_vs_deep p = 0.002
  # middle_vs_deep p = 0.428
