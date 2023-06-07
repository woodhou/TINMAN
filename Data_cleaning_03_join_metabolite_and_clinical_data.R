require(tidyverse); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

# Phase one data ----------------------------------------------------------

clinical <- readRDS('TINMAN_clinical_raw_data.rds')

clinical.metadata <- readRDS('TINMAN_sample_metadata.rds')

feces <- readRDS('TINMAN_feces_20221129.rds')

#' Append clinical metadata such as BCM ID and group/neutropenia status.
feces <- feces %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, GROUP_ID, GROUP_NAME), 
            by = 'PARENT_SAMPLE_NAME')

#' Per Josie's email on 9/26/2022, for patients 53 and 61, the second sample (D, E, or F) is the baseline sample and
#' the first sample (A, B, or C) is the neutropenic sample.
#' The Metabolon group assignments do not reflect this.
feces <- feces %>% mutate(GROUP_ID = ifelse(CLIENT_IDENTIFIER %in% c('053F1','061F2'), 'NEUT_PRE',
                                     ifelse(CLIENT_IDENTIFIER %in% c('053B2','061C2'), 'NEUT_POST', GROUP_ID)),
                          TIMEPOINT = ifelse(str_detect(GROUP_ID, 'POST'), 'POST', 'PRE')) %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER:GROUP_NAME, TIMEPOINT, starts_with('C'))

saveRDS(feces, 'TINMAN_feces_20221129_2.rds')

#' A record of the samples that were assigned to each group.
assignments <- feces %>% select(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT)

write.xlsx(assignments, 
           '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/TINMAN/R outputs/sample.assignments.20221013.xlsx')

rm(list = ls()); gc()

# Phase two data ----------------------------------------------------------

new.metadata <- readRDS('TINMAN_phase_2_sample_metadata.rds')

phase.two <- readRDS('TINMAN_phase_2_feces_20230502.rds')

#' Append clinical metadata such as BCM ID and group/neutropenia status.
phase.two <- phase.two %>% 
  left_join(select(new.metadata, PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, GROUP_ID, GROUP_NAME, TIME_POINT), 
            by = 'PARENT_SAMPLE_NAME') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER:GROUP_NAME, TIME_POINT, starts_with('C'))

saveRDS(phase.two, 'TINMAN_phase_2_feces_20230502_2.rds')

#' A record of the samples that were assigned to each group.
assignments <- phase.two %>% select(CLIENT_IDENTIFIER, GROUP_NAME, TIME_POINT)

write.xlsx(assignments, 
           '//smb-main.ad.bcm.edu/genepi2/Old_genepi2/Jeremy/TINMAN/R outputs/sample.assignments.phase.two.20230502.xlsx')

rm(list = ls()); gc()
