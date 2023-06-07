require(tidyverse); require(mixOmics); require(xlsx)

setwd('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/')

#' Concentration table
tinman <- readRDS('TINMAN_feces_20221129_2.rds')

#' Metabolite metadata
metabs <- readRDS('TINMAN_feces_metadata.rds') %>% 
  filter(SUPER_PATHWAY != 'Xenobiotics' | (SUPER_PATHWAY == 'Xenobiotics' & SUB_PATHWAY == 'Bacterial/Fungal')) %>% 
  pull(CHEM_ID)

#' Remove xenobiotics, except for bacterial/fungal metabolites..
tinman <- tinman %>% dplyr::select(PARENT_SAMPLE_NAME:TIMEPOINT, all_of(metabs))

# Matrix of metabolite abundances (x) and class memberships (y) needs for sparse PLS-DA.
x <- as.matrix(tinman[,6:ncol(tinman)])
y <- as.factor(tinman$GROUP_ID)

# Tuning keepX, the number of features that will be used for each component.
#list.keepX <- c(5:10,  seq(20, 50, 10))

#set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
#tune.splsda.srbct <- tune.splsda(x, y, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
#                                 validation = 'Mfold', folds = 3, 
#                                 dist = 'max.dist', progressBar = T,
#                                 measure = "BER", test.keepX = list.keepX,
#                                 nrepeat = 50)   # we suggest nrepeat = 50

#error <- tune.splsda.srbct$error.rate
#ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
#ncomp

#select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
#select.keepX

MyResult.splsda.final <- splsda(x, y, ncomp = 2, keepX = c(50,50))

my.plot <- plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
                     ellipse = TRUE, title="sPLS-DA - final result")

#' Recover and plot loading vectors for components 1 and 2
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)
comp1.loadings <- selectVar(MyResult.splsda.final, comp=1)$value %>% 
  mutate(CHEM_ID = rownames(.))

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)
comp2.loadings <- selectVar(MyResult.splsda.final, comp=2)$value %>% 
  mutate(CHEM_ID = rownames(.))

#' Save plots
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.plot.20221212.svg',
    width = 10, height = 10)

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - final result")

dev.off()

svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.comp1.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 1)

dev.off()

#' Save plots
svg('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/plsa.da.comp2.loadings.svg',
    width = 10, height = 10)

plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean', comp = 2)

dev.off()

#' Metabolite metadata
metabs <- readRDS('TINMAN_feces_metadata.rds')

#' Append metadata to loadings
comp1.loadings <- left_join(comp1.loadings, metabs, by = 'CHEM_ID')

comp2.loadings <- left_join(comp2.loadings, metabs, by = 'CHEM_ID')

write.xlsx(comp1.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_PLS-DA_loadings_20221212.xlsx',
           sheetName = 'component1',
           row.names = F)

write.xlsx(comp2.loadings,
           '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_PLS-DA_loadings_20221212.xlsx',
           sheetName = 'component2',
           row.names = F,
           append = T)

auc.plsda <- auroc(MyResult.splsda.final)

#' Evaluate overlap between loadings and univariate results
uni <- read.xlsx('//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/TINMAN_metabolomics_analysis_20221017.xlsx',
                 sheetIndex = 2)

uni <- uni %>% filter(t.test.p < 0.05 | kruskal.test.p < 0.05)

loadings <- c(comp1.loadings$CHEM_ID, comp2.loadings$CHEM_ID)

intersect(uni$CHEM_ID, loadings)
