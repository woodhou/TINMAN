require(tidyverse)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Read in autoscaled, imputed data from the new cohort.
phase.two <- readRDS('TINMAN_phase_2_feces_20230502.rds')

#' Read in autoscaled, imputed data from the old cohort.
phase.one <- readRDS('TINMAN_feces_20221012.rds')

common.metabs <- intersect( names(select(phase.one, -PARENT_SAMPLE_NAME)), names(select(phase.two, -PARENT_SAMPLE_NAME)) )

phase.one <- phase.one %>% 
  select(PARENT_SAMPLE_NAME, all_of(common.metabs))

phase.two <- phase.two %>% 
  select(PARENT_SAMPLE_NAME, all_of(common.metabs))

#' Does this list of common metabolites include all the differentially abundant ones we identified in the first phase?
control <- readRDS('TINMAN_control_unscaled_fold_changes_20221012.rds') %>% 
  filter(control.t.test.p < 0.05)
names(control) <- names(control) %>% str_remove('control.')

neutropenic <- readRDS('TINMAN_neutropenic_unscaled_fold_changes_20221012.rds') %>% 
  filter(neutropenic.t.test.p < 0.05)
names(neutropenic) <- names(neutropenic) %>% str_remove('neutropenic.') 

candidates <- union(control$metabolite, neutropenic$metabolite)

un.candidates <- subset(candidates, !candidates %in% common.metabs)

#' All eight metabolites ID'd in neutropenic subjects and 61 of 63 identified in control subjects were also detected in new dataset.
#' Exceptions are a MCFA called pelargonate (9:0) and and a lipid involved in BCAA metabolism called 2-hydroxyadipate.
candidates <- subset(candidates, candidates %in% common.metabs)
