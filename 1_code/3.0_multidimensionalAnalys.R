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
metadata <- read.csv('sample-metadata.byArms.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqTable.d1 <- read.csv('freqTablebyFraction.csv', header = TRUE, row.names = 1)
countTable <- read.csv('../../9-2mmdata/2mmFeqTable.csv', header = TRUE)


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
freqTable.d2.arms <- freqTable.d1%>%filter(fraction=="byARMS") # for ARMS
freqTable.d2.motile <- freqTable.d1%>%filter(fraction=="motile") # for motile
freqTable.d2.sessile <- freqTable.d1%>%filter(fraction=="sessile") # for sessile

freqTable.arms <- freqTable.d2.arms[,4:21] # for ARMS
freqTable.motile <- freqTable.d2.motile[,4:21] # for motile
freqTable.sessile <- freqTable.d2.sessile[,4:21] # for sessile
freqTable.count <- countTable[,-1] # for count

# for count data


# twist the following distM to switch between bray/jaccard
# distM <- vegdist(t(freqTable), method = "bray")
distM.arms <- vegdist(t(sqrt(freqTable.arms)), method = "bray")
distM.mo <- vegdist(t(sqrt(freqTable.motile)), method = "bray")
distM.si <- vegdist(t(sqrt(freqTable.sessile)), method = "bray")

distM.count <- vegdist(t(freqTable.count), method = "bray") # for count data


# plot all the baseline by square root bray curtis 

  dist_matrix1 <- prcomp(distM.count)
  dist_matrix2 <- prcomp(distM.arms)
  dist_matrix3 <- prcomp(distM.mo)
  dist_matrix4 <- prcomp(distM.si)
  
  
  
  evalue.1 <- fviz_eig(dist_matrix1, addlabels = TRUE) 
  evalue.2 <- fviz_eig(dist_matrix2, addlabels = TRUE) 
  evalue.3 <- fviz_eig(dist_matrix3, addlabels = TRUE) 
  evalue.4 <- fviz_eig(dist_matrix4, addlabels = TRUE) 
  
  
  dataPlot <- data.frame(matrix(nrow=18*4, ncol=6))
  colnames(dataPlot) <- c("ARMS","PC1", "PC2", "treatment", "water",'fraction')
  dataPlot$ARMS <- rep(names(freqTable.arms),4)
  
  dataPlot[1:18,]$PC1 <- dist_matrix1$x[,1]
  dataPlot[19:36,]$PC1 <- dist_matrix2$x[,1]
  dataPlot[37:54,]$PC1 <- dist_matrix3$x[,1]
  dataPlot[55:72,]$PC1 <- dist_matrix4$x[,1]
  
  dataPlot[1:18,]$PC2 <- dist_matrix1$x[,2]
  dataPlot[19:36,]$PC2 <- dist_matrix2$x[,2]
  dataPlot[37:54,]$PC2 <- dist_matrix3$x[,2]
  dataPlot[55:72,]$PC2 <- dist_matrix4$x[,2]
  
  
  
  dataPlot$treatment <- rep(metadata$treatment,4)
  dataPlot$water <- rep(metadata$depth,4)
  dataPlot$coral <- rep(metadata$coral,4)
  dataPlot$coral <-  factor(dataPlot$coral, levels = c('acropora', 'pavona', 'platygyra', 'acroporaplatygyra', 'acroporapavonaplatygyra', 'no'))
  dataPlot$fraction <- c(rep('count',18), rep('arms',18),rep('mo',18),rep('si',18))

#### 1. fig. for pub ####

  
  # Define distance order and labels
  water_order <- c("shallow", "middle", "deep")
  water_labels <- c("Close", "Middle", "Far")
  
  # Define coral order and labels
  coral_order <- c("no", "acropora", "pavona", "platygyra", 
                   "acroporaplatygyra", "acroporapavonaplatygyra")
  coral_labels <- c("Control", "Mono-Acropora", "Mono-Pavona", 
                    "Mono-Platygyra", "Mixed", "Polyculture")

  
  # Custom colors (colorblind-friendly, high contrast)
  custom_colors <- c(
    "Close"  = "#7DC0A6",   # teal green
    "Middle" = "#919FC7",   # blue
    "Far"    = "#ED936B"    # orange
  )
  
    
  # Extract percentage variance explained from evalue
  # evalue$data[1,2] and evalue$data[2,2] contain the percentages
  
# count data 
  # Create axis labels with percentages
  x_label <- paste0("PCo1 (", round(evalue.1$data[1, 2], digits = 1), "%)")
  y_label <- paste0("PCo2 (", round(evalue.1$data[2, 2], digits = 1), "%)")
  
  baseline1 <- dataPlot %>% filter(fraction=='count') %>%
    mutate(
      water = factor(water, levels = water_order, labels = water_labels),
      coral = factor(coral, levels = coral_order, labels = coral_labels)
    ) %>%
    ggplot(aes(x = PC1, y = PC2, color = water)) +
    # Ellipses first (background) — emphasizes distance clustering
    stat_ellipse(type = "norm", aes(group = water), 
                 level = 0.95, lwd = 1.2, alpha = 0.8) +
    # Points — shapes distinguish coral but color dominates
    geom_point(aes(shape = coral), size = 6, alpha = 0.9, stroke = 1.2) +
    # Axis labels with percentage variance explained
    xlab(x_label) +
    ylab(y_label) +
    # Custom color palette for distance (prominent)
    scale_color_manual(values = custom_colors, name = "Distance to Reef") +
    # Shape palette for coral (secondary)
    scale_shape_manual(
      values = c(4,15, 16, 17, 2, 1),
      name = "Coral Assemblage"
    ) +
    theme_classic(base_size = 18) +
    theme(
      # Remove axis tick labels (cleaner for PCoA)
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # Axis titles
      axis.title = element_text(face = "bold", size = 18),
      # Legend
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.key = element_blank(),
      # Plot margins
      plot.margin = margin(1, 1, 1, 1, "cm")
    )

# ARMS data
  # Create axis labels with percentages
  x_label <- paste0("PCo1 (", round(evalue.2$data[1, 2], digits = 1), "%)")
  y_label <- paste0("PCo2 (", round(evalue.2$data[2, 2], digits = 1), "%)")
  
  baseline2 <- dataPlot %>% filter(fraction=='arms') %>%
    mutate(
      water = factor(water, levels = water_order, labels = water_labels),
      coral = factor(coral, levels = coral_order, labels = coral_labels)
    ) %>%
    ggplot(aes(x = PC1, y = PC2, color = water)) +
    # Ellipses first (background) — emphasizes distance clustering
    stat_ellipse(type = "norm", aes(group = water), 
                 level = 0.95, lwd = 1.2, alpha = 0.8) +
    # Points — shapes distinguish coral but color dominates
    geom_point(aes(shape = coral), size = 6, alpha = 0.9, stroke = 1.2) +
    # Axis labels with percentage variance explained
    xlab(x_label) +
    ylab(y_label) +
    # Custom color palette for distance (prominent)
    scale_color_manual(values = custom_colors, name = "Distance to Reef") +
    # Shape palette for coral (secondary)
    scale_shape_manual(
      values = c(4,15, 16, 17, 2, 1),
      name = "Coral Assemblage"
    ) +
    theme_classic(base_size = 18) +
    theme(
      # Remove axis tick labels (cleaner for PCoA)
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # Axis titles
      axis.title = element_text(face = "bold", size = 18),
      # Legend
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.key = element_blank(),
      # Plot margins
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
#motile
  # Create axis labels with percentages
  x_label <- paste0("PCo1 (", round(evalue.3$data[1, 2], digits = 1), "%)")
  y_label <- paste0("PCo2 (", round(evalue.3$data[2, 2], digits = 1), "%)")
  
  baseline3 <- dataPlot %>% filter(fraction=='mo') %>%
    mutate(
      water = factor(water, levels = water_order, labels = water_labels),
      coral = factor(coral, levels = coral_order, labels = coral_labels)
    ) %>%
    ggplot(aes(x = PC1, y = PC2, color = water)) +
    # Ellipses first (background) — emphasizes distance clustering
    stat_ellipse(type = "norm", aes(group = water), 
                 level = 0.95, lwd = 1.2, alpha = 0.8) +
    # Points — shapes distinguish coral but color dominates
    geom_point(aes(shape = coral), size = 6, alpha = 0.9, stroke = 1.2) +
    # Axis labels with percentage variance explained
    xlab(x_label) +
    ylab(y_label) +
    # Custom color palette for distance (prominent)
    scale_color_manual(values = custom_colors, name = "Distance to Reef") +
    # Shape palette for coral (secondary)
    scale_shape_manual(
      values = c(4,15, 16, 17, 2, 1),
      name = "Coral Assemblage"
    ) +
    theme_classic(base_size = 18) +
    theme(
      # Remove axis tick labels (cleaner for PCoA)
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # Axis titles
      axis.title = element_text(face = "bold", size = 18),
      # Legend
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.key = element_blank(),
      # Plot margins
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
# si 
  # Create axis labels with percentages
  x_label <- paste0("PCo1 (", round(evalue.4$data[1, 2], digits = 1), "%)")
  y_label <- paste0("PCo2 (", round(evalue.4$data[2, 2], digits = 1), "%)")
  
  baseline4 <- dataPlot %>% filter(fraction=='si') %>%
    mutate(
      water = factor(water, levels = water_order, labels = water_labels),
      coral = factor(coral, levels = coral_order, labels = coral_labels)
    ) %>%
    ggplot(aes(x = PC1, y = PC2, color = water)) +
    # Ellipses first (background) — emphasizes distance clustering
    stat_ellipse(type = "norm", aes(group = water), 
                 level = 0.95, lwd = 1.2, alpha = 0.8) +
    # Points — shapes distinguish coral but color dominates
    geom_point(aes(shape = coral), size = 6, alpha = 0.9, stroke = 1.2) +
    # Axis labels with percentage variance explained
    xlab(x_label) +
    ylab(y_label) +
    # Custom color palette for distance (prominent)
    scale_color_manual(values = custom_colors, name = "Distance to Reef") +
    # Shape palette for coral (secondary)
    scale_shape_manual(
      values = c(4,15, 16, 17, 2, 1),
      name = "Coral Assemblage"
    ) +
    theme_classic(base_size = 18) +
    theme(
      # Remove axis tick labels (cleaner for PCoA)
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # Axis titles
      axis.title = element_text(face = "bold", size = 18),
      # Legend
      legend.position = "none",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.key = element_blank(),
      # Plot margins
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # combine them together 
  combined_plot <- plot_grid(
    baseline1,baseline2,baseline3,baseline4, 
    labels = c("Macrofauna", "ARMS(motile + sessile)","Motile","sessile"),
    ncol = 2,
    nrow = 2,
    align = "h"
  )

  
#### 2. adonis ####
# depth = distance 
  
  set.seed(123)
  model1.arsm <- adonis2(distM.arms~depth+coral, 
                    data=metadata,
                  strata = metadata$depth,	
                    by = "margin")
  
  set.seed(123)
  model1.mo <- adonis2(distM.mo~depth+coral, 
                    data=metadata,
                    strata = metadata$depth,	
                    by = "margin")
  
  set.seed(123)
  model1.si <- adonis2(distM.si~depth+coral, 
                    data=metadata,
                    strata = metadata$depth,	
                    by = "margin")
  
  set.seed(123)
  model1.count <- adonis2(distM.count~depth+coral, 
                    data=metadata,
                    strata = metadata$depth,	
                    by = "margin")
  
  
# dispersal check


  model2.arms <- vegan::betadisper(distM.arms, metadata$coral)  
  anova(model2.arms) # homogeneous 

  model2.mo <- vegan::betadisper(distM.mo, metadata$coral)  
  anova(model2.mo) # homogeneous 
  
  model2.si <- vegan::betadisper(distM.si, metadata$coral)  
  anova(model2.si) # homogeneous 
  
  model2.count <- vegan::betadisper(distM.count, metadata$coral)  
  anova(model2.count) # homogeneous 
  
  
    
  model3.arms <- vegan::betadisper(distM.arms, metadata$depth)  
  anova(model3.arms) # homogeneous 
  
  model3.mo <- vegan::betadisper(distM.mo, metadata$depth)  
  anova(model3.mo) # homogeneous 
  
  model3.si <- vegan::betadisper(distM.si, metadata$depth)  
  anova(model3.si) # homogeneous 
  
  model3.count <- vegan::betadisper(distM.count, metadata$depth)  
  anova(model3.count) # homogeneous 
  

#### 3. Distance-based redundancy analysis (dbRDA) ####

  dbRDA_model.a <- capscale(distM.arms ~ coral + depth,data = metadata, add = TRUE) # 0.4669
  dbRDA_model.m <- capscale(distM.mo ~ coral + depth,data = metadata, add = TRUE) # 0.4641
  dbRDA_model.s <- capscale(distM.si ~ coral + depth,data = metadata, add = TRUE) # 0.4692
  dbRDA_model.c <- capscale(distM.count ~ coral + depth,data = metadata, add = TRUE) # 0.5354

  # middle_vs_deep p 0.152