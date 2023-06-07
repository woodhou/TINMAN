
# Prepare environment -----------------------------------------------------

require(tidyverse); require(ggrepel)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

volcano <- function(data, title, cutoff){
  
  plot <- ggplot(data = data, aes(x = log(fold.change, 2), y = -log10(t.test.p))) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 'twodash') +
    
    geom_vline(xintercept = -2, linetype = 'longdash') +
    geom_vline(xintercept = 2, linetype = 'twodash') +
    
    geom_point(aes(color = color)) +
    
    geom_text_repel(aes(label = ifelse(t.test.p < cutoff, CHEMICAL_NAME, ''))) +
    
    annotate('text', x = 11, y = 3, label = 'More abundant at follow-up', fontface = 'bold') +
    annotate('text', x = -11, y = 3, label = 'More abundant at baseline', fontface = 'bold') +
    
    labs(x = 'Log2 fold change (follow-up/baseline, unscaled values)',
         y = "-Log10 p-value (Student's t-test)",
         title = title) +
    
    scale_x_continuous(limits = c(-15,15),
                       breaks = c(-15,-10,-5,0,5,10,15)) +
    
    scale_y_continuous(limits = c(0,4)) +
    
    theme_classic() +
    
    theme(axis.title = element_text(size = 14, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 16, face = 'bold'),
          
          legend.position = 'none')
  
  return(plot)
  
}

# Prepare data ------------------------------------------------------------

#' Includes metabolite name.
metadata <- readRDS('TINMAN_feces_metadata.rds') 
metadata <- metadata %>% select(CHEM_ID, CHEMICAL_NAME)

control <- readRDS('TINMAN_control_unscaled_fold_changes_20221012.rds') %>% 
  mutate(color = ifelse(control.t.test.p < 0.05 & abs(log(control.fold.change, 2)) >= 2, "A", 
                 ifelse(control.t.test.p < 0.05 & abs(log(control.fold.change, 2)) < 2, 'B', 'C'))) %>% 
  left_join(metadata, by = c('metabolite' = 'CHEM_ID'))

names(control) <- names(control) %>% str_remove('control.')

neutropenic <- readRDS('TINMAN_neutropenic_unscaled_fold_changes_20221012.rds') %>% 
  mutate(color = ifelse(neutropenic.t.test.p < 0.05 & abs(log(neutropenic.fold.change, 2)) >= 2, "A", 
                 ifelse(neutropenic.t.test.p < 0.05 & abs(log(neutropenic.fold.change, 2)) < 2, 'B', 'C'))) %>% 
  left_join(metadata, by = c('metabolite' = 'CHEM_ID'))

names(neutropenic) <- names(neutropenic) %>% str_remove('neutropenic.') 

# Generate plots ----------------------------------------------------------

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Volcano_plots/')

svg('volcano_plot_controls_unscaled_20221017.svg', height = 8, width = 8)

volcano(control, 'CONTROLS', 0.01)

dev.off()

svg('volcano_plot_neutropenic_unscaled_20221017.svg', height = 8, width = 8)

volcano(neutropenic, 'NEUTROPENIC SUBJECTS', 0.05)

dev.off()