require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

# Phase one data ----------------------------------------------------------

#' Metabolon metadata for identified biochemicals.
feces.metadata <- readRDS('TINMAN_feces_metadata.rds')

#' For the most part, we will exclude xenobiotics. Bacterial/fungal metabolies are still of interest, however.
endogenous <- feces.metadata %>% 
  filter(SUPER_PATHWAY != "Xenobiotics" | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

exogenous <- feces.metadata %>% 
  filter(SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY != 'Bacterial/Fungal') %>% 
  pull(CHEM_ID)

#' These data have not been imputed or scaled. 
feces <- readRDS('TINMAN_feces_raw_data.rds')

#' Impute metabolites that are not xenobiotics up to half the minimum observed abundance.
feces <- feces %>% mutate_at(all_of(endogenous), ~ ifelse(is.na(.x), min(.x, na.rm = T)/2, .x))

saveRDS(feces, 
        'TINMAN_feces_unscaled_20221129.rds')

#' Autoscale the imputed data (subtract the mean and divide by the SD).            
for (i in 2:ncol(feces)){
  
  mean <- feces[, i] %>% unlist() %>% mean(na.rm = T)
  sd <- feces[ , i] %>% unlist() %>% sd(na.rm = T)
  
  feces[ , i] <- (feces[, i] - mean)/sd 
  
}

saveRDS(feces,
        'TINMAN_feces_20221129.rds')

# Phase two data ----------------------------------------------------------

#' New raw data.
phase.two <- readRDS('TINMAN_phase_2_feces_raw_data.rds')

#' New metabolite metadata file.
phase.two.metadata <- readRDS('TINMAN_phase_2_feces_metadata.rds')

#' For the most part, we will exclude xenobiotics. Bacterial/fungal metabolies are still of interest, however.
endogenous <- phase.two.metadata %>% 
  filter(SUPER_PATHWAY != "Xenobiotics" | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

exogenous <- phase.two.metadata %>% 
  filter(SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY != 'Bacterial/Fungal') %>% 
  pull(CHEM_ID)

#' Impute endogenous metabolites up to half the minimum observed abundance.
phase.two <- phase.two %>% mutate_at(all_of(endogenous), ~ ifelse(is.na(.x), min(.x, na.rm = T)/2, .x))

saveRDS(phase.two, 'TINMAN_phase_2_feces_unscaled_20230502.rds')

#' Autoscale the imputed data (subtract the mean and divide by the SD).            
for (i in 2:ncol(phase.two)){
  
  mean <- phase.two[, i] %>% unlist() %>% mean(na.rm = T)
  sd <- phase.two[ , i] %>% unlist() %>% sd(na.rm = T)
  
  phase.two[ , i] <- (phase.two[, i] - mean)/sd 
  
}

saveRDS(phase.two,
        'TINMAN_phase_2_feces_20230502.rds')
