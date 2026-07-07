#### 1.0 lib and data set WD ####
#mac 
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/')
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
metadate <- read.csv('sample-metadata.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
preRM.d1 <- read.csv('1-decontaminate/preConRM.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')[c(-1,-56),]
postRM.d1 <- read.csv('1-decontaminate/freqTable.relax.contemRM.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
finOTU.106 <- read.csv('3-deconPostaaClus97/BOLD/freqTableCleanBy100.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
finOTU.500 <- read.csv('3-deconPostaaClus97/BOLD/freqTableCleanBy500.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
finOTU.ses <- read.csv('3-deconPostaaClus97/BOLD/freqTableCleanBySes.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')
freqALL <- read.csv('3-deconPostaaClus97/BOLD/feaTableRAW.csv', row.name=1, header = TRUE, check.names = FALSE, comment.char = '')


#### 2.0 transform dataset ####
postRM.d2 <- as.data.frame(colSums(postRM.d1))
all(row.names(postRM.d2) == row.names(preRM))
names(postRM.d2) <- 'conRM'

# merge decontem with pre decontem
preRM <- preRM.d1[,c(-5,-8,-10)]
postRM <- as.data.frame(cbind(preRM, postRM.d2))
postRM[, 3:8] <- lapply(postRM[, 3:8], as.numeric)
all(postRM[,8]/postRM[,7]<1) #quick check, T all less than before 

postRM.fin <- postRM %>%  arrange(fraction,ARMS)


# merge AAtranslate and OTU to RM 
finOTU.100.fin <- as.data.frame(colSums(finOTU.106))
finOTU.500.fin <- as.data.frame(colSums(finOTU.500))
finOTU.ses.fin <- as.data.frame(colSums(finOTU.ses))
names(finOTU.100.fin) <- 'OTU'
names(finOTU.500.fin) <- 'OTU'
names(finOTU.ses.fin) <- 'OTU'
finOTU.100.fin$fraction <- '106um'
finOTU.500.fin$fraction <- '500um'
finOTU.ses.fin$fraction <- 'sessile'
finOTU.100.fin$ARMS <- row.names(finOTU.100.fin)
finOTU.500.fin$ARMS <- row.names(finOTU.500.fin)
finOTU.ses.fin$ARMS <- row.names(finOTU.ses.fin)
finOTU.fin <- rbind(finOTU.100.fin,finOTU.500.fin,finOTU.ses.fin)

# check alignment 
all(postRM.fin$ARMS == finOTU.fin$ARMS) # true 
all(postRM.fin$fraction == finOTU.fin$fraction) # true 
all.read <- cbind(postRM.fin,finOTU.fin)[,-(10:11)]

# percentage through 
all.read$perRead <- all.read[,9]/all.read[,3]


#### 7.0 Sample-Based Richness/Coverage #### 
community_matrix <- t(freqALL)

# Check dimensions
dim(community_matrix)  # Should be 54 x number_of_OTUs


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
                       sample = min(rowSums(community_matrix)), 
                       label = TRUE)

# get all elbow
# Function to find elbow for a single rarecurve output
find_elbow <- function(rare_curve_data, threshold) {
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

# Apply to all 54 samples in rare_data
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
    sample_name <- names(freqALL)[i]
  }
  
  # Find elbow for this sample
  elbow <- find_elbow(rare_data[[i]], 0.0001)
  
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
split_parts <- strsplit(elbow_results$Sample, "\\.")
elbow_results$arms <- as.numeric(sapply(split_parts, "[", 1))
elbow_results$fraction <- sapply(split_parts, "[", 2)

elbow_results.fin <- elbow_results %>% arrange(fraction,arms) 

# check alignment 
postRM.fin$ARMS  
elbow_results.fin$arms # true 
all(all.read$fraction == elbow_results.fin$fraction) # true 

# let's merge and write 
out.file <- cbind(all.read, elbow_results.fin)




