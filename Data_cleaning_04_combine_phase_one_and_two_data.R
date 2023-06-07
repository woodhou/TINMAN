#require(tidyverse)

#' This script is deprecated. Metabolon supplied a merged file that supercedes this.
#setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Merge TINMAN phase one and two datasets.
#one <- readRDS('TINMAN_feces_20221129_2.rds')

#one.metabs <- names(one) %>% subset(., str_detect(., '^C[:digit:]{1,}'))

#two <- readRDS('TINMAN_phase_2_feces_20230502_2.rds')

#two.metabs <- names(two) %>% subset(., str_detect(., '^C[:digit:]{1,}'))

#joint.metabs <- intersect(one.metabs, two.metabs)

#combined <- bind_rows(
#  select(one, PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, all_of(joint.metabs)),
#  select(two, PARENT_SAMPLE_NAME, GROUP_ID, GROUP_NAME, all_of(joint.metabs))
#)

#saveRDS(combined, 'TINMAN_combined_stool_data_20230508.rds')
