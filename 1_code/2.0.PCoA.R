#### 1.0 lib and data set WD ####
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


#### 2.0 data curation check out sharing OTUs by time #### 
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

#### 3.0 permanova ####
# calculate datamatrix based on squaroot 
# twist the freqTable to switch among motile/sessile/ARMS
freqTable.d2 <- freqTable.d1%>%filter(fraction=="byARMS") # for plot 4
freqTable.d3 <- freqTable.d1%>%filter(fraction=="motile") # for plot 2
freqTable.d4 <- freqTable.d1%>%filter(fraction=="sessile") # for plot 3

freqTable.m1 <- freqTable.d2[,4:21] # by ARMS
freqTable.m2 <- freqTable.d3[,4:21] # by motile
freqTable.m3 <- freqTable.d4[,4:21] # by sessile
freqTable.m4 <- countTable[,-1] # by count 


distM1 <- vegdist(t(sqrt(freqTable.m1)), method = "bray")
distM2 <- vegdist(t(sqrt(freqTable.m2)), method = "bray")
distM3 <- vegdist(t(sqrt(freqTable.m3)), method = "bray")
distM4 <- vegdist(t(freqTable.m4), method = "bray")

# what about coral? 

set.seed(123)
model1a <- adonis2(distM1~coral+depth, data=metadate, strata = metadate$depth, by = "margin",
                    permutations = 999)
set.seed(123)
model2a <- adonis2(distM2~coral+depth, data=metadate, strata = metadate$depth, by = "margin",
                    permutations = 999)
set.seed(123)
model3a <- adonis2(distM3~coral+depth, data=metadate, strata = metadate$depth, by = "margin",
                    permutations = 999)
set.seed(123)
model4a <- adonis2(distM4~coral+depth, data=metadate, strata = metadate$depth,  by = "margin",
                    permutations = 999)
  
model1a  # by ARMS
model2a  # by motile
model3a  # by sessil
model4a  # by count
  

# pairwise 
pairwise_1 <- pairwise.adonis2(distM1 ~ coral, #strata = metadate$depth,
  data = metadate)

#### 4.0 betadisper ####
# this is to test the community disperse
disp1 <- betadisper(distM1, metadate$coral)
disp2 <- betadisper(distM2, metadate$coral)
disp3 <- betadisper(distM3, metadate$coral)
disp4 <- betadisper(distM4, metadate$coral)

dispm1 <- anova(disp1) # ns
dispm2 <- anova(disp2) # ns
dispm3 <- anova(disp3) # ns
dispm4 <- anova(disp4) # ns

boxplot(disp1, main = "Distances to Centroid")
boxplot(disp2, main = "Distances to Centroid")
boxplot(disp3, main = "Distances to Centroid")
boxplot(disp4, main = "Distances to Centroid")
# all ns, which is awesome 


#### 5.0 dbRDA1 #### 

# ARMS dataset
dbRDA1 <- capscale(distM1 ~ coral + depth, data = metadate, add = TRUE)

# Motile dataset
dbRDA2 <- capscale(distM2 ~ coral + depth, data = metadate, add = TRUE)

# Sessile dataset
dbRDA3 <- capscale(distM3 ~ coral + depth, data = metadate, add = TRUE)

# Count dataset
dbRDA4 <- capscale(distM4 ~ coral + depth, data = metadate, add = TRUE)

summary(dbRDA1)
summary(dbRDA2)
summary(dbRDA3)
summary(dbRDA4)

#### 6.0 replace the PCoA plot #### 
# Basic plot
plot(dbRDA1)

# Nicer plot with ggplot2
# Extract scores
site_scores <- scores(dbRDA1, display = "sites")

# Create data frame for plotting
plot_data <- data.frame(
  RDA1 = site_scores[,1],
  RDA2 = site_scores[,2],
  coral = metadate$coral,
  depth = metadate$depth
)


ggplot(plot_data, aes(x = RDA1, y = RDA2, color = coral, shape = depth)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = coral), level = 0.95) +
  labs(
    title = "dbRDA - ARMS Dataset",
    x = paste0("RDA1 (", round(dbRDA1$CCA$eig[1]/sum(dbRDA1$CCA$eig)*100, 1), "%)"),
    y = paste0("RDA2 (", round(dbRDA1$CCA$eig[2]/sum(dbRDA1$CCA$eig)*100, 1), "%)")
  ) +
  theme_classic()

#### 7.0 Sample-Based Richness/Coverage #### 
# Transpose so rows = samples, columns = species
community_matrix <- t(freqTable.m1)

# Check dimensions
dim(community_matrix)  # Should be 18 x number_of_OTUs

# Rarefy to 50 individuals (or incidences)
rarefied_richness <- rarefy(community_matrix, sample = 50)

# View results
print(rarefied_richness)


# Plot rarefaction curves
prarefaction <- rarecurve(community_matrix, 
          step = 10,           # Plot every 10 individuals
          sample = min(rowSums(community_matrix)),  # Minimum sample size
          xlab = "Sample Size", 
          ylab = "Species Richness",
          col = "blue",
          label = TRUE)        # Show sample labels




# where is the elbow? 
# 1. Extract rarecurve data
rare_data <- rarecurve(community_matrix, step = 10, 
                       sample = min(rowSums(community_matrix)), # it's HHW21
                       label = TRUE)

# get all elbow
# Function to find elbow for a single rarecurve output
find_elbow <- function(rare_curve_data, threshold = 0.01) {
  # Extract richness values
  richness <- rare_curve_data
  
  # Extract sample sizes from attribute
  sample_sizes <- attr(rare_curve_data, "Subsample")
  
  # Calculate slopes
  slopes <- diff(richness) / diff(sample_sizes)
  
  # Find first point where slope < threshold
  elbow_idx <- which(slopes < threshold)[1] + 1
  
  # If no elbow found, return NA
  if (is.na(elbow_idx) || elbow_idx > length(sample_sizes)) {
    return(data.frame(
      sample_size = NA,
      richness = NA,
      slope = NA
    ))
  }
  
  # Return elbow point
  return(data.frame(
    sample_size = sample_sizes[elbow_idx],
    richness = richness[elbow_idx],
    slope = slopes[elbow_idx - 1]
  ))
}

# Apply to all 18 samples in rare_data
elbow_results <- data.frame(
  Sample = character(),
  Elbow_Sample_Size = numeric(),
  Elbow_Richness = numeric(),
  Slope_at_Elbow = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(rare_data)) {
  # Get the sample name (if available) or use index
  sample_name <- names(rare_data)[i]
  if (is.null(sample_name) || sample_name == "") {
    sample_name <- paste0("Sample_", i)
  }
  
  # Find elbow for this sample
  elbow <- find_elbow(rare_data[[i]], threshold = 0.01)
  
  # Add to results
  elbow_results <- rbind(elbow_results, data.frame(
    Sample = sample_name,
    Elbow_Sample_Size = elbow$sample_size,
    Elbow_Richness = elbow$richness,
    Slope_at_Elbow = elbow$slope
  ))
}



# View results
print(elbow_results)
elbow_results$Sample <- metadate$arms
elbow_results$sumread <- rowSums(community_matrix)

max(elbow_results$Elbow_Sample_Size) # 18071/HHW16
min(elbow_results$sumread) # 190188/HHW21




# make better looking 
# 2. Convert to tidy data frame
rare_list <- list()
for (i in 1:length(rare_data)) {
  sample_name <- metadate$arms[i]
  temp <- as.data.frame(rare_data[[i]])
  colnames(temp) <-  "richness"
  temp$sample_size <- as.numeric(sub("N", "", rownames(temp)))
  temp$sample <- sample_name
  temp$coral <- metadate$coral[match(sample_name, metadate$arms)]
  rare_list[[i]] <- temp
}
rare_df <- bind_rows(rare_list)

# mean cruve 
# rare_df.d1 <- as.data.frame(rare_df%>% group_by(coral,sample_size) %>% summarise(
#  coral = coral,
#  richness.mean = mean(richness),
#  sample_size=sample_size
# ))


# 3. Define colors
coral_colors <- c(
  "acropora" = "#E41A1C",
  "pavona" = "#377EB8",
  "platygyra" = "#4DAF4A",
  "no" = "#984EA3",
  "acroporaplatygyra" = "#FF7F00",
  "acroporapavonaplatygyra" = "#F781BF"
)

# 4. Plot
p <- ggplot(rare_df, aes(x = sample_size, y = richness, 
                         color = coral, group = sample)) +
  geom_line(size = 0.8, alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(
    values = coral_colors,
    labels = c(
      "acropora" = "Mono-Acropora",
      "pavona" = "Mono-Pavona",
      "platygyra" = "Mono-Platygyra",
      "no" = "Control",
      "acroporaplatygyra" = "Mixed",
      "acroporapavonaplatygyra" = "Polyculture"
    )
  ) +
  labs(
    x = "Sample Size (Number of Individuals)",
    y = "OTU Richness",
    color = "Coral Assemblage",
    title = "Rarefaction Curves by Coral Assemblage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

print(p)


#### 8.0 let anova again? ####
# Rarefy to the minimum sample size across all samples
min_sample <- min(rowSums(community_matrix))
rarefied_richness <- rarefy(community_matrix, sample = min_sample)

# Add to metadata
metadate$rarefied_richness <- rarefied_richness

# ANOVA
aov_rare <- aov(rarefied_richness ~ coral, data = metadate)
summary(aov_rare)

# Tukey post-hoc
TukeyHSD(aov_rare)

(alphaTable[1,c(-1,-20)] - rarefied_richness[1:18])/alphaTable[1,c(-1,-20)]


