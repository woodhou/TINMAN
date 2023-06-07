require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Metabolon metadata for identified biochemicals.
metab.metadata <- readRDS('TINMAN_merged_feces_metadata.rds')

#' For the most part, we will exclude xenobiotics. Bacterial/fungal metabolies are still of interest, however.
endogenous <- metab.metadata %>% 
  filter(SUPER_PATHWAY != "Xenobiotics" | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

exogenous <- metab.metadata %>% 
  filter(SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY != 'Bacterial/Fungal') %>% 
  pull(CHEM_ID)

#' These data have not been imputed or scaled. 
metab <- readRDS('TINMAN_merged_feces_raw_data.rds')

#' Impute metabolites that are not xenobiotics up to half the minimum observed abundance.
metab <- metab %>% mutate_at(all_of(endogenous), ~ ifelse(is.na(.x), min(.x, na.rm = T)/2, .x))

saveRDS(metab, 
        'TINMAN_merged_metab_unscaled_20230606.rds')

#' Autoscale the imputed data (subtract the mean and divide by the SD).            
for (i in 2:ncol(metab)){
  
  mean <- metab[, i] %>% unlist() %>% mean(na.rm = T)
  sd <- metab[ , i] %>% unlist() %>% sd(na.rm = T)
  
  metab[ , i] <- (metab[, i] - mean)/sd 
  
}

saveRDS(metab,
        'TINMAN_merged_metab_20230606.rds')
