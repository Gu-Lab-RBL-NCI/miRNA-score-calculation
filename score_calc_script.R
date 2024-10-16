library(readxl)
library(dplyr) 
library(stringr)

setwd("~/Desktop/NGS/score_calculation/")

#import reference table with sequence information
seq_ref <- read.csv("ref_files/pre_mirna_sequences.csv")

#import reference table with 5' energy fron RNA fold and RNA eval
end_energy <- read.csv("ref_files/end_energy.csv")

#import reference table for pre and mature miRNA names
miRNA_strand <- read.csv('ref_files/miRNA_strand.csv')

#merge all reference by primary sequence
merged_ref <- merge(end_energy, miRNA_strand, all=T)

temp_5p <- merged_ref %>% 
  filter(strand == '5p') %>%
  select(c(4:6))
temp_3p <- merged_ref %>% 
  filter(strand=='3p') 
merged_ref_1 <- merge(temp_5p, temp_3p, by='pri_miRNA_unique', suffixes = c("_5p", "_3p"))

temp_5p <- merged_ref %>% 
  filter(strand == '5p')
temp_3p <- merged_ref %>% 
  filter(strand=='3p')  %>%
  select(c(4:6))

merged_ref_2 <- merge(temp_5p, temp_3p, by='pri_miRNA_unique', suffixes = c("_5p", "_3p"))
merged_ref <- rbind.data.frame(merged_ref_1, merged_ref_2) 
merged_ref <- merged_ref[is.finite(merged_ref$energy_3),]
merged_ref <- merged_ref %>%
  select(-9) %>%
  distinct()
rm(end_energy, merged_ref_1, merged_ref_2, miRNA_strand, seq_ref, temp_3p, temp_5p)

#calculate score
FN_5p <- str_sub(merged_ref$mature_seq_5p, 1, 1)
FN_3p <- str_sub(merged_ref$mature_seq_3p, 1, 1)
merged_ref <- cbind(merged_ref, FN_5p, FN_3p)

merged_ref <- merged_ref %>%
  mutate(nt_score_5p=0, 
         nt_score_3p=0)

#score calculation - values from Suzuki et al., 2015
merged_ref[merged_ref$FN_5p == 'U', 11] <- 2.63963
merged_ref[merged_ref$FN_3p == 'U', 12] <- 2.63963
merged_ref[merged_ref$FN_5p == 'A', 11] <- 1.056746
merged_ref[merged_ref$FN_3p == 'A', 12] <- 1.056746
merged_ref[merged_ref$FN_5p == 'G', 11] <- 0.37825
merged_ref[merged_ref$FN_3p == 'G', 12] <- 0.37825

merged_ref <- merged_ref %>%
  mutate(score_5p = 0.66005*energy_5/100+nt_score_5p, 
         score_3p = 0.66005*energy_3/100+nt_score_3p)
rm(FN_5p, FN_3p)

#save scores for miRNA based on precursor reference
merged_ref <- merged_ref %>%
  select(4,2,7,13,14)
write.csv(merged_ref, "miRNA_scores.csv", row.names = F)

