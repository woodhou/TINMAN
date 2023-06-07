# Prep environment --------------------------------------------------------

require(tidyverse); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' A function to perform t tests and Kruskal-Wallis tests comparing metabolite abundances pre- vs. post- separately in control and neutropenic subjects.
univariate.tests <- function(x){
  
  #' Initialize an empty data frame to hold results for pre-/post- comparisons in the control and neutropenic groups.
  results <- data.frame()
  
  #' For each metabolite...
  for (i in seq_along(metabolites)){
    
    print(i)
    
    data <- select(x, GROUP_NAME, TIMEPOINT, all_of(metabolites[i]))
    names(data) <- c('GROUP_NAME','TIMEPOINT','METABOLITE')
    data <- data %>% filter(!is.na(METABOLITE))
    
    counts <- data %>% count(TIMEPOINT, .drop = F)
    
    #' If that metabolte was not detected in the index group, mark as missing.
    if (nrow(data) == 0){
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = NA,
                       mean.post = NA,
                       t.test.p = NA,
                       kruskal.test.p = NA)
      
      results <- rbind(results, result)
      
    }
    
    #' If it was detected in more than one subject from the index group at each timepoint, calculate the means at baseline and follow up and compare by univariate tests.
    else if (all(counts$n > 1) & var(data$METABOLITE) != 0){  
      
      t.test <- data %>% 
        t.test(METABOLITE ~ TIMEPOINT, data = .)
      
      kruskal.test <- data %>% 
        kruskal.test(METABOLITE ~ TIMEPOINT, data = .)
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = t.test$estimate[2],
                       mean.post = t.test$estimate[1],
                       t.test.p = t.test$p.value,
                       kruskal.test.p = kruskal.test$p.value)
      
      results <- rbind(results, result)
      
    }
    
    #' Otherwise, report means at baseline and follow up.
    else {
      
      means <- aggregate(METABOLITE ~ TIMEPOINT, data = data, mean, drop = F)
      
      result <- tibble(metabolite = metabolites[i],
                       mean.pre = means[2,2],
                       mean.post = means[1,2],
                       t.test.p = NA,
                       kruskal.test.p = NA)
      
      results <- rbind(results, result)
      
    }
    
  }
  
  return(results)
  
}

# Load data ---------------------------------------------------------------

#' Read in scaled, imputed data.
feces <- readRDS('TINMAN_feces_20221129_2.rds') %>% 
  mutate(TIMEPOINT = factor(TIMEPOINT))

#' Metabolites to test.
metabolites <- feces %>% select(-c(PARENT_SAMPLE_NAME:TIMEPOINT)) %>% names()

#' Split into control and neutropenic datasets.
control = filter(feces, GROUP_NAME == 'Control')

neutropenic = filter(feces, GROUP_NAME == 'neutropenic')

# Count the number of samples metabolites were detected, by gro --------

feces.raw.data <- readRDS('TINMAN_feces_raw_data.rds') %>% %>% %>% %>% %>% 

clinical.metadata <- readRDS('TINMAN_sample_metadata.rds')

#' Append clinical metadata such as BCM ID and group/neutropenia status.
feces.raw.data <- feces.raw.data %>% 
  left_join(select(clinical.metadata, PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, GROUP_ID, GROUP_NAME), 
            by = 'PARENT_SAMPLE_NAME')

#' Per Josie's email on 9/26/2022, for patients 53 and 61, the second sample (D, E, or F) is the baseline sample and
#' the first sample (A, B, or C) is the neutropenic sample.
#' The Metabolon group assignments do not reflect this.
feces.raw.data <- feces.raw.data %>% mutate(GROUP_ID = ifelse(CLIENT_IDENTIFIER %in% c('053F1','061F2'), 'NEUT_PRE',
                                                              ifelse(CLIENT_IDENTIFIER %in% c('053B2','061C2'), 'NEUT_POST', GROUP_ID)),
                                            TIMEPOINT = ifelse(str_detect(GROUP_ID, 'POST'), 'POST', 'PRE')) %>% 
  select(PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER:GROUP_NAME, TIMEPOINT, starts_with('C')) %>% 
  mutate(across(C30:C100022610, ~ ifelse(is.na(.x), 0, 1)))

counts <- tibble()

for (i in metabolites){
  
  data <- feces.raw.data %>% select(GROUP_ID, all_of(i))
  names(data) <- c('GROUP_ID','METABOLITE')
  
  tmp <- data %>% 
    group_by(GROUP_ID) %>% 
    tally(METABOLITE) %>% 
    mutate(METABOLITE = rep(i,4))
  
  counts <- rbind(counts, tmp)
  
}

counts <- counts %>% 
  pivot_wider(names_from = GROUP_ID, values_from = n)

names(counts) <- c('metabolite', 'control.n.post', 'control.n.pre', 'neutropenic.n.post', 'neutropenic.n.pre') 

#saveRDS(counts,
#        'TINMAN_metabolite_detection_counts_20221012.rds')

# Perform tests -----------------------------------------------------------

control.results <- univariate.tests(control) %>% 
  rename_with(~ paste0('control.', .x), !metabolite)

neutropenic.results <- univariate.tests(neutropenic) %>% 
  rename_with(~ paste0('neutropenic.', .x), !metabolite)

results <- left_join(control.results, neutropenic.results, by = 'metabolite')

results <- left_join(counts, results, by = 'metabolite')

# Calculate fold changes for metabolites detected @ both times ------------

control <- results %>% 
  select(metabolite, starts_with('control')) %>% 
  filter(control.n.post > 0, control.n.pre > 0) %>% 
  mutate(control.fold.change = control.mean.post/control.mean.pre) %>% 
  mutate(control.fold.change = ifelse(control.mean.post > control.mean.pre, abs(control.fold.change), 
                               ifelse(control.mean.post < control.mean.pre, abs(control.fold.change)*-1,         
                                      control.fold.change)))

saveRDS(control, 'TINMAN_control_fold_changes_20221012.rds')

neutropenic <- results %>% 
  select(metabolite, starts_with('neutropenic')) %>% 
  filter(neutropenic.n.post > 0, neutropenic.n.pre > 0) %>% 
  mutate(neutropenic.fold.change = neutropenic.mean.post/neutropenic.mean.pre) %>% 
  mutate(neutropenic.fold.change = ifelse(neutropenic.mean.post > neutropenic.mean.pre, abs(neutropenic.fold.change), 
                                   ifelse(neutropenic.mean.post < neutropenic.mean.pre, abs(neutropenic.fold.change)*-1,         
                                          neutropenic.fold.change)))

saveRDS(neutropenic, 'TINMAN_neutropenic_fold_changes_20221012.rds')

# Export results to Excel -------------------------------------------------

#' A convenience function for identifying metabolites that were detected in one group but not another.
id.metabolites <- function(var1, var2){
  
  output <- counts %>% 
    filter(.data[[var1]] > 0, .data[[var2]] == 0) %>% 
    left_join(metabolite.metadata, by = c('metabolite' = 'CHEM_ID'))
  
}

#' Metadata about metabolites, for annotating outputs.
metabolite.metadata <- readRDS('TINMAN_feces_metadata.rds')
metabolite.metadata <- metabolite.metadata %>% select(CHEM_ID, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY, HMDB, KEGG)

#' List of metabolites with fold changes, by group.
names(control.results) <- names(control.results) %>% str_remove('control.')
names(neutropenic.results) <- names(neutropenic.results) %>% str_remove('neutropenic.')

control.results <- control.results %>% 
  mutate(group = 'Control')

neutropenic.results <- neutropenic.results %>% 
  mutate(group = 'Neutropenic')

results <- bind_rows(control.results, neutropenic.results) %>% 
  mutate(significant.change = ifelse(t.test.p < 0.05, 1, 0))

results <- right_join(metabolite.metadata, results, by = c('CHEM_ID' = 'metabolite'))

write.xlsx(results,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
           sheetName = 'MeanAbundances',
           append = F)

#' List of metabolites detected only pre-treatment, by group.
baseline.control.only <- id.metabolites('control.n.pre','neutropenic.n.pre')
baseline.neutropenic.only <- id.metabolites('neutropenic.n.pre', 'control.n.pre')

baseline <- bind_rows(baseline.control.only, baseline.neutropenic.only)

write.xlsx(baseline,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
           sheetName = 'Baseline_Between_Groups',
           append = T)

#' List of metabolites detected only post-treatment, by group.
fu.control.only <- id.metabolites('control.n.post', 'neutropenic.n.post')
fu.neutropenic.only <- id.metabolites('neutropenic.n.post', 'control.n.post')

fu <- bind_rows(fu.control.only, fu.neutropenic.only)

write.xlsx(fu,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
           sheetName = 'FollowUp_Between_Groups',
           append = T)

#' TODO: List of metabolites detected only at one timepoint, within groups.
controls.baseline.only <- id.metabolites('control.n.pre', 'control.n.post')
controls.fu.only <- id.metabolites('control.n.post', 'control.n.pre')

controls <- bind_rows(controls.baseline.only, controls.fu.only)

write.xlsx(controls,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
           sheetName = 'Controls_Between_Timepoints',
           append = T)

neutropenic.baseline.only <- id.metabolites('neutropenic.n.pre', 'neutropenic.n.post')
neutropenic.fu.only <- id.metabolites('neutropenic.n.post', 'neutropenic.n.pre')

neutropenic <- bind_rows(neutropenic.baseline.only, neutropenic.fu.only)

write.xlsx(neutropenic,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
           sheetName = 'Neutropenic_Between_Timepoints',
           append = T)

rm(list = ls()); gc()
