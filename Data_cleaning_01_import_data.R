require(tidyverse); require(readxl)

#' TODO: Recover clinical data for phase two patients, if needed.

# Phase one data ----------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/BAYL-03-22MD/')

feces <- read_xlsx('BAYL-03-22MD FECES DATA TABLES.XLSX', sheet = 'Batch-normalized Data') 
names(feces) <- c(names(feces)[1], paste0('C',names(feces)[2:ncol(feces)]))

feces.metadata <- read_xlsx('BAYL-03-22MD FECES DATA TABLES.XLSX', sheet = 'Chemical Annotation') %>% 
  mutate(CHEM_ID = paste0('C', CHEM_ID))

clinical <- read_xlsx('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/TINMAN REDCAP DATA_jms.xlsx',
                      sheet = 'cohort 1')

names(clinical) <- names(clinical) %>% tolower() %>% str_replace_all('_', '.')

clinical <- clinical %>% select(patient:anc.lowest)

clinical.metadata <- read_xlsx('BAYL-03-22MD FECES DATA TABLES.XLSX',
                               sheet = 'Sample Meta Data') 

saveRDS(clinical, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_clinical_raw_data.rds')
saveRDS(clinical.metadata, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_sample_metadata.rds')

saveRDS(feces, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_feces_raw_data.rds')
saveRDS(feces.metadata, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_feces_metadata.rds')

# Phase two data ----------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/BAYL-03-22MD_phase_two/')

feces.new <- read_xlsx('BAYL-13-22MD DATA TABLES.XLSX', sheet = 'Batch-normalized Data') 
names(feces.new) <- c(names(feces.new)[1], paste0('C',names(feces.new)[2:ncol(feces.new)]))

feces.metadata.new <- read_xlsx('BAYL-13-22MD DATA TABLES.XLSX', sheet = 'Chemical Annotation') %>% 
  mutate(CHEM_ID = paste0('C', CHEM_ID))

#clinical.new <- read_xlsx('//smb-main.ad.bcm.edu/genepi2/TINMAN/Datasets/TINMAN REDCAP DATA_jms.xlsx',
#                      sheet = 'cohort 1')

names(clinical.new) <- names(clinical.new) %>% tolower() %>% str_replace_all('_', '.')

#clinical.new <- clinical.new %>% select(patient:anc.lowest)

clinical.metadata.new <- read_xlsx('BAYL-13-22MD DATA TABLES.XLSX',
                                   sheet = 'Sample Meta Data') 

#saveRDS(clinical.new, 
#        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_clinical_raw_data.rds')
saveRDS(clinical.metadata.new, 
        '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_sample_metadata.rds')

saveRDS(feces.new, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_feces_raw_data.rds')
saveRDS(feces.metadata.new, '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_phase_2_feces_metadata.rds')
