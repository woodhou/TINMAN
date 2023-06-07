
# Prepare environment -----------------------------------------------------

require(tidyverse); require(ggrepel)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

volcano <- function(data, title, cutoff){
  
  plot <- ggplot(data = data, aes(x = fold.change, y = -log10(t.test.p))) +
    
    geom_hline(yintercept = -log10(0.05), linetype = 'twodash') +
    
    #geom_vline(xintercept = -2, linetype = 'longdash') +
    geom_vline(xintercept = 1, linetype = 'twodash') +
    
    geom_point() +
    
    geom_text_repel(aes(label = ifelse(t.test.p < cutoff, CHEMICAL_NAME, ''))) +
    
    labs(x = 'Fold change (Follow-up/Baseline, Autoscaled Values)',
         y = "-Log10 p-value (Student's t-test)",
         title = title) +
    
    scale_x_continuous(limits = c(-9,10),
                       breaks = c(-5,1,5)) +
    
    scale_y_continuous(limits = c(0,4)) +
    
    theme_classic() +
    
    theme(axis.title = element_text(size = 14, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 16, face = 'bold'))
  
  return(plot)
  
}

# Prepare data ------------------------------------------------------------

#' Includes metabolite name.
metadata <- readRDS('TINMAN_feces_metadata.rds') 
metadata <- metadata %>% select(CHEM_ID, CHEMICAL_NAME)

control <- readRDS('TINMAN_control_fold_changes_20221012.rds') %>% 
  left_join(metadata, by = c('metabolite' = 'CHEM_ID'))

names(control) <- names(control) %>% str_remove('control.')

neutropenic <- readRDS('TINMAN_neutropenic_fold_changes_20221012.rds') %>% 
  left_join(metadata, by = c('metabolite' = 'CHEM_ID'))

names(neutropenic) <- names(neutropenic) %>% str_remove('neutropenic.')

setwd('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/Figures/Volcano_plots/')

svg('volcano_plot_controls_20221013.svg', height = 8, width = 8)

volcano(control, 'CONTROLS', 0.01)

dev.off()

svg('volcano_plot_neutropenic_20221013.svg', height = 8, width = 8)

volcano(neutropenic, 'NEUTROPENIC SUBJECTS', 0.05)

dev.off()