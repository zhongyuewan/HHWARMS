## set WD 
# mac
setwd('/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#2 corest/#2 3D AR/correst metabarcoding/1_data/8-coral/')
rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)

## get data
coral <- read.csv('mel.csv', header = TRUE, check.names = FALSE, comment.char = '')


Acro <- coral%>%filter(CoralGenera=='Acropora') # 126 fragments
Pav <- coral%>%filter(CoralGenera=='Pavona') # 126 fragments
Plat <- coral%>%filter(CoralGenera=='Platygyra') # 126 fragments

surA <- nrow(Acro%>% filter(status=='present'))/126 # 0.984127
surPa <- nrow(Pav%>% filter(status=='present'))/126 # 0.8492063
surPl <- nrow(Plat%>% filter(status=='present'))/126 # 0.9761905

summary(Acro%>% filter(status=='present')%>%.$netGain) # mean 11.0
sd(Acro%>% filter(status=='present')%>%.$netGain) # 4.47627

summary(Pav%>% filter(status=='present')%>%.$netGain) # mean 6.289
sd(Pav%>% filter(status=='present')%>%.$netGain) # 2.276182

summary(Plat%>% filter(status=='present')%>%.$netGain) # mean 4.814
sd(Plat%>% filter(status=='present')%>%.$netGain) # 1.736509


