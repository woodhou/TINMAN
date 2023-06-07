require(tidyverse); require(readxl)

# Load merged phase one and two data --------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/BAYL-03-22MD_merged/')

#' Load Metabolon data.
feces <- read_xlsx('BAYL-13-22MD CO DATA TABLES.XLSX', sheet = 'Batch-norm Data Common') 
names(feces) <- c(names(feces)[1], paste0('C',names(feces)[2:ncol(feces)]))

metab.metadata <- read_xlsx('BAYL-13-22MD CO DATA TABLES.XLSX', sheet = 'Chemical Annotation Common') %>% 
  mutate(CHEM_ID = paste0('C', CHEM_ID))

#' Load and combine clinical metadata from each phase.
clin.one.meta <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_sample_metadata.rds') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, starts_with('GROUP'))

clin.two.meta <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_sample_metadata.rds') %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, starts_with('GROUP'))

clinical.metadata <- bind_rows(clin.one.meta, clin.two.meta)

saveRDS(clinical.metadata, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_clinical_data.rds')

saveRDS(feces, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_feces_raw_data.rds')

saveRDS(metab.metadata, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_merged_feces_metadata.rds')