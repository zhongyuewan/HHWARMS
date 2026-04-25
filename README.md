# Restoring Beyond Corals: Strong Spatial Gradients in Benthic Cryptobenthic Community Assembly


#### author list hidden for double-anonymised review 

## Abstract 

<img align="right" src="3_figures/Picture4png.png" width=450> 

Coral reef degradation has been recorded globally and increasing efforts have been put into coral restoration aiming to rebuild functional reef ecosystems. Polyculture—restoring multiple coral species together—has been identified as a promising strategy to promote coral growth and survivorship. But it remains unknown whether it’s effective at recruiting the broader cryptobenthic community—the essential but often overlooked foundation of a functional reef. Therefore, we used standardized Autonomous Reef Monitoring Structures (ARMS) and metabarcoding to compare cryptic biodiversity in a three-year restoration experiment featuring mono-, mix-, poly-culture, and unseeded control plots. We found that while proximity to a healthy reef was the strongest driver of community composition, active coral seeding significantly increased total taxonomic richness compared to controls, which trended toward oyster-dominated, lower-diversity states. However, the distinct culture types (mono-, mix-, poly-) did not yet produce substantially different communities in this early stage. Our results confirmed that active restoration is superior to passive substrate deployment for biodiversity recruitment and indicate that the hypothesized benefits of polyculture may require longer timelines to allow coral biomass and structural complexity to develop.



## Table of Contents

### Supporting Materials 
  1. [Code](1_code)
  2. [Data](2_data)
  3. [Figures](3_figures)
  4. [Tables](4_supplementaryTable)


### Sequence processing pipeline 
1. [Import & cutadap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.1_importAndCutAdapt.sh): import raw sequence data (.fastq) into Qiime artefacts (.qza) and remove PCR adaptors.
2. [Denoise-paired](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.2_denoiseAndPair.sh): remove sequences likely induced by error and merge the reverse/forward reads.
3. [Decontam](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.3_decontam.r): a process to look into the negative control and remove sequences that might have come from sample contamination.
4. [Amino Acid translation](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.4_aaTranslate.r): translate DNA sequence into amino acid and remove sequences with one of the following conditions: 1) any STOP codon, 2) >3 deletion, 3) any frameshift, 4) any insertion.
5. [Cluster all sequences](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.5_clusterReads.sh) by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.
6. [Taxonomic assignment](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.6_taxAssign.sh) with BLAST against two different libraries: 1) McIlroy et al. 2024 & 2) Medori2 (GB260).

### Data Analysis 
1. Environmental data
   - [Heatmap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.1_eData_heatmap.r) (Figure 1d, Table 1)
   - [MPA east vs west](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.2_eastVSwest.r) (Table S2)
2. Species richness by ARMS 
   - [Merge richness from all three fractions](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.3_combinFractionbyARMS.r) (Table S1)
   - [Environmental data ~ species richness](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.4_eDATAvsRichness.r) (Table 2) 
3. Community composition
   - [PCoA](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.5_PCoA.r) (Figure 2)
   - [Permutational Multivariate Analysis of Variance (adonis2)](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.6_adonis2.r) 
   - [Diverging Bar Chart & Chi-Square analysis](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.8_sidewayBar.r) (Figure 3)

     
