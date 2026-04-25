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


## get data
metadata<- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('dna-sequences.fasta')
alphaTable <- read.csv('alphaTable.byArms.csv', header = TRUE)


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
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="motile")
# freqTable.d2 <- freqTable.d1%>%filter(fraction=="sessile")

freqTable <- freqTable.d2[,4:21]

# twist the following distM to switch between bray/jaccard
distM <- vegdist(t(sqrt(freqTable)), method = "bray")
# distM <- vegdist(t(freqTable), method = "bray") doing sqrt is the right thing 


# Run NMDS
nmds_result <- metaMDS(distM, k = 2, trymax = 100)
stress_value <- nmds_result$stress

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds_result))
nmds_scores$Sample <- rownames(nmds_scores)
# nmds_scores$Treatment <- metadata$coral
# nmds_scores$Treatment <- metadata$treatment
# nmds_scores$Treatment <- metadata$depth

# Calculate convex hulls for each treatment group
hull_data <- data.frame()
for(treatment in unique(nmds_scores$Treatment)) {
  treatment_data <- nmds_scores[nmds_scores$Treatment == treatment, ]
  hull_indices <- chull(treatment_data$NMDS1, treatment_data$NMDS2)
  hull_data <- rbind(hull_data, treatment_data[hull_indices, ])
}

# 5. CREATE THE PLOT
nmds_plot1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Add convex hulls
  geom_polygon(data = hull_data, 
               aes(fill = Treatment, color = Treatment), 
               alpha = 0.2, 
               linetype = "dashed",
               linewidth = 0.5) +
  
  # Add points
  geom_point(aes(color = Treatment, shape = Treatment), size = 3) +
  
  # Add stress value to plot
  annotate("text", 
           x = min(nmds_scores$NMDS1, na.rm = TRUE), 
           y = max(nmds_scores$NMDS2, na.rm = TRUE),
           label = paste("Stress =", round(stress_value, 3)),
           hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  
  # Customize appearance
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "NMDS Plot Based on Bray-Curtis Dissimilarity",
       x = "NMDS Axis 1",
       y = "NMDS Axis 2") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_shape_manual(values = c(16, 17, 15, 18, 23, 25)) + # Different shapes for treatments
  theme(legend.position = "none")

combined_pieplot <- plot_grid(nmds_plot, nmds_plot1,ncol = 2)
