require(tidyverse)

tinman <- readRDS("//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_feces_20221129_2.rds")

tinman <- tinman |>
  select(-c(CLIENT_IDENTIFIER, GROUP_NAME, TIMEPOINT)) 

write.csv(tinman,
          file = '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_MetaboAnalyst_power_analysis_input_data_20230207.csv',
          row.names = F)

# Univariate analysis power, by hand --------------------------------------

require(tidyverse); require(xlsx); require(ggsci)

tinman <- readRDS("//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_feces_20221129_2.rds")

r <- 4
z.beta <- 0.84 # z score for beta=0.2
z.alpha.seq <- c(1.96, 2.81, 3.48, 4.06) # z scores for alpha=0.05, 0.005, 0.0005, 0.00005, and 0.000005)
delta.seq <- seq(0.25, 1.55, by = 0.05)
sd <- 1

power.calcs <- tibble()

total.n.calc <- function(r, z.beta, z.alpha, delta, sd){
  
  n <- ( (r+1)*(sd^2)*((z.beta+z.alpha)^2) ) / (r*(delta^2))
  
  n <- ceiling(n)*5
  
  return(n)
  
}

for (i in seq_along(delta.seq)){
  
  for (j in seq_along(z.alpha.seq)){
    
    my.n <- total.n.calc(r, z.beta, z.alpha.seq[j], delta.seq[i], sd)
    
    new.data <- tibble(delta = delta.seq[i],
                       z = z.alpha.seq[j],
                       min.sample = my.n)
    
    power.calcs <- rbind(power.calcs, new.data)
    
  }
  
}

power.calcs <- power.calcs |>
  mutate(z = factor(z, labels = c('0.05', '5E-3', '5E-4', '5E-5')))

power.plot <- ggplot(data = power.calcs, aes(x = delta, y = min.sample, color = z)) +
  
  geom_point(size = 3) +
  
  geom_line(size= 1.5) +
  
  geom_hline(yintercept = 200, linetype = 'dashed') + 
  
  scale_x_continuous(limits = c(0.25, 1.5),
                     breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) +
  
  scale_y_continuous(limits = c(0, 500),
                     breaks = c(0,100,200,300,400 )) +
  
  scale_color_npg() +
  
  labs(title = 'SAMPLE SIZE REQUIRED TO ACHIEVE 80% POWER',
       y = 'TOTAL N',
       x = 'DIFFERENCE IN AUTOSCALED MEAN ABUNDANCES',
       color = 'p-value') +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(),
        
        axis.title = element_text(face = 'bold', size = 14),
        axis.title.y = element_text(margin = margin(0,20,0,0)),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.text = element_text(face = 'bold', size = 12),
        
        legend.position = c(0.9, 0.6),
        legend.background = element_rect(color = 'black', linewidth = 1.5),
        legend.title = element_text(face = 'bold', size = 14, hjust = 0.5),
        legend.text = element_text(face = 'bold', size = 12),
        
        plot.title = element_text(face = 'bold', size = 16)
        
        )

power.plot

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/R01.aim.3.power.for.univariate.tests.20230207.svg',
    width = 8, height = 8)

print(power.plot)

dev.off()

tinman.results <- read.xlsx('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
                            sheetName = 'Mean_Abundances') |>
  as_tibble() |>
  mutate(mean.difference = abs(mean.post - mean.pre))

aggregate(mean.difference ~ significant.change, summary, data = tinman.results)
