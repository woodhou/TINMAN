require(tidyverse); require(eulerr)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

results <- readRDS('TINMAN_metabolite_detection_counts_20221012.rds')

# Identify metabolites detected only pre-treatment ------------------------

pre.transplant.controls <- results %>% 
  filter(control.n.pre > 0 & control.n.post == 0) %>% 
  pull(metabolite)

pre.transplant.neutropenia <- results %>% 
  filter(neutropenic.n.pre > 0 & neutropenic.n.post == 0) %>% 
  pull(metabolite)

pre.transplant.common <- intersect(pre.transplant.controls, pre.transplant.neutropenia)

# Identify metabolites detected only post-transplant ----------------------

post.transplant.controls <- results %>% 
  filter(control.n.post > 0 & control.n.pre == 0) %>% 
  pull(metabolite)

post.transplant.neutropenia <- results %>% 
  filter(neutropenic.n.post > 0 & neutropenic.n.pre == 0) %>% 
  pull(metabolite)

post.transplant.common <- intersect(post.transplant.controls, post.transplant.neutropenia)

# Euler diagrams ----------------------------------------------------------

#' Input data.
pre.euler <- c('Control' = 35, 'Neutropenic' = 32, 'Control&Neutropenic' = 8)
post.euler <- c('Control' = 7, 'Neutropenic' = 32, 'Control&Neutropenic' = 1)

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Euler_diagrams/')

#' Save plots.
svg('Euler_diagram_species_detected_pre_TX_only_20221012.svg', height = 8, width = 8)

plot(euler(pre.euler), quantities = T, main = 'Metabolites Detected Only at Baseline')

dev.off()

svg('Euler_diagram_species_detected_post_TX_only_20221012.svg', height = 8, width = 8)

plot(euler(post.euler), quantities = T, main = 'Metabolites Detected Only at Follow-up')

dev.off()

# Count species detected in each group ------------------------------------

count.metabolites <- function(x){
  
  out <- results %>% filter(.data[[x]] > 0) %>% pull(metabolite)
  
  return(out)
  
}

control.post <- count.metabolites('control.n.post')
control.pre <- count.metabolites('control.n.pre')
control.intersect <- intersect(control.post, control.pre)
control.euler <- c('Baseline' = 35, 'Follow-up' = 7, 'Baseline&Follow-up' = 1155)

neutropenic.post <- count.metabolites('neutropenic.n.post')
neutropenic.pre <- count.metabolites('neutropenic.n.pre')
neutropenic.intersect <- intersect(neutropenic.post, neutropenic.pre)
neutropenic.euler <- c('Baseline' = 32, 'Follow-up' = 32, 'Baseline&Follow-up' = 1126)

baseline.intersect <- intersect(control.pre, neutropenic.pre)
baseline.euler <- c('Control' = 1190-1146, 'Neutropenic' = 1158-1146, 'Control&Neutropenic' = 1146)

follow.up.intersect <- intersect(control.post, neutropenic.post)
follow.up.euler <- c('Control' = 1162-1131, 'Neutropenic' = 1158-1131, 'Control&Neutropenic' = 1131)

#' Save plots.
svg('Euler_diagram_control_pre_post_20221012.svg', height = 8, width = 8)

plot(venn(control.euler), quantities = T, main = 'Metabolites Detected at Baseline and Follow-up (Controls)')

dev.off()

svg('Euler_diagram_neutropenic_pre_post_20221012.svg', height = 8, width = 8)

plot(venn(neutropenic.euler), quantities = T, main = 'Metabolites Detected at Baseline and Follow-up (Neutropenic Subjects)')

dev.off()

svg('Euler_diagram_baseline_20221012.svg', height = 8, width = 8)

plot(venn(baseline.euler), quantities = T, main = 'Richness of the Control & Neutropenic Metabolomes (Baseline)')

dev.off()

svg('Euler_diagram_follow_up_20221012.svg', height = 8, width = 8)

plot(venn(follow.up.euler), quantities = T, main = 'Richness of the Control & Neutropenic Metabolomes (Follow-up)')

dev.off()
