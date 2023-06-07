require(tidyverse); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

metab <- readRDS('TINMAN_merged_metab_unscaled_20230606.rds')
clin <- readRDS('TINMAN_merged_clinical_data.rds')

#' Per Josie's email on 9/26/2022, for patients 53 and 61, the second sample (D, E, or F) is the baseline sample and
#' the first sample (A, B, or C) is the neutropenic sample.
#' The Metabolon group assignments do not reflect this.
clin <- readRDS('TINMAN_merged_clinical_data.rds') %>% 
  mutate(GROUP_ID = ifelse(CLIENT_IDENTIFIER %in% c('053F1','061F2'), 'NEUT_PRE',
                    ifelse(CLIENT_IDENTIFIER %in% c('053B2','061C2'), 'NEUT_POST', GROUP_ID)))

#' Group names are not the same in phase one and two.
clin <- clin %>% 
  mutate(GROUP_ID = GROUP_ID %>% toupper() %>% str_replace_all('CON', 'CTRL'),
         GROUP_NAME = ifelse(GROUP_NAME == 'neutropenic', 'Neutropenic', GROUP_NAME)) %>% 
  select(-CLIENT_IDENTIFIER)
         
#' Merge metabolite and clinical data.
tinman <- metab %>% 
  left_join(clin, by = 'PARENT_SAMPLE_NAME') %>% 
  select(PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, starts_with('C'))

saveRDS(tinman, 'TINMAN_merged_20230606.rds')
