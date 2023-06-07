require(tidyverse); require(simstudy); require(party); require(descr)

# Simulate data -----------------------------------------------------------

tin <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/TINMAN_feces_unscaled_20221129.rds')

#' Replace missing obs with half the minimum, as in the real data.
tin <- tin %>% 
  mutate(across(-PARENT_SAMPLE_NAME, ~ ifelse(is.na(.x), min(.x, na.rm = T)/2, .x)))

#' Logarithms of variable means and SDs.
mu <- tin %>% select(-1) %>% apply(2, mean) %>% + 1 %>%  log()

sigma.squared <- tin %>% 
  select(-1) %>% 
  apply(2, var) %>% 
  ifelse(. < 1, 1.01, .) %>% 
  log()

def <- defData(varname = 'var1', dist = 'normal', formula = mu[1], variance = sigma.squared[1])

for (i in seq_along(mu)){
  
  if (i == 1){
    
    def <- defData(varname = paste0('var',i), dist = 'normal', formula = mu[i], variance = sigma.squared[i])
    
  }
  
  else {
    
    def <- defData(def, varname = paste0('var',i), dist = 'normal', formula = mu[i], variance = sigma.squared[i]) 
    
  }
  
}

dt <- genData(1000, def) %>% 
  mutate(across(starts_with('var'), ~ exp(.x)))

#' Autoscale the data.
dt <- dt %>% 
  mutate(across(starts_with('var'), ~ (.x - mean(.x))/sd(.x) ))

#' Add a binomial variable to represent neutropenia.
dt <- dt %>% 
  mutate(neutropenia = rbinom(1000, 1, 0.2)) %>% 
  select(id, neutropenia, starts_with('var')) %>% 
  as_tibble()
  
saveRDS(dt, 
        file = '//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/simulated.data.for.R01.power.calcs.rds')

# Spike in an effect ------------------------------------------------------

dt <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/simulated.data.for.R01.power.calcs.rds')

#' Choose a column at random.
my.col <- names(dt)[round(runif(1, 3, 1212))]

#' Choose an effect size at random.
effect.size <- runif(1, 0, 2)

#' Identify any columns correlated with the chosen column at r>=0.8.
c <- cor(dt[ , 3:ncol(dt)])
c <- c[,my.col] |> 
  (\(x) subset(x, x >= 0.8))() |>
  names()

#' Add the effect size to neutropenic subjects in the chosen column and any highly correlated columns.  
for (i in c){
  
  my.col <- i
  
  dt <- dt %>% 
    mutate('{my.col}' := .data[[my.col]] + neutropenia*effect.size)
  
}

# Build a random forests model --------------------------------------------

#' Sample neutropenic subjects and controls.
dt <- dt %>% 
  mutate(runif = runif(1000, 0, 1))

neut <- dt %>% 
  filter(runif < 0.2, neutropenia == 1) %>% 
  select(-runif)

control <- dt %>% 
  filter(neutropenia == 0) %>% 
  arrange(runif) %>% 
  slice(1:(nrow(neut) * 2))

data <- bind_rows(neut, control)

cif <- cforest(neutropenia ~ .,
               data = data,
               control = cforest_unbiased(ntree = 500, mtry = 10))

cif.importance <- varimp(cif)

cif.importance <- sort(cif.importance, decreasing = T)


# Put it all together -----------------------------------------------------

dt <- readRDS('//smb-main.ad.bcm.edu/genepi3/JeremySchraw/TINMAN/Datasets/simulated.data.for.R01.power.calcs.rds')

n.sims <- 7500

#' Choose a column at random.
my.col <- names(dt)[round(runif(n.sims, 3, 1212))]

#' Choose an effect size at random.
effect.size <- runif(n.sims, 0, 2)

#' Initialize an empty data frame to hold CIF results.
cif.results <- data.frame()

for (i in seq_along(my.col)){

  print(i)
  
  #' Identify any columns correlated with the chosen column at r>=0.8.
  c <- cor(dt[ , 3:ncol(dt)])
  c <- c[ , my.col[i] ] |> 
    (\(x) subset(x, x >= 0.8))() |>
    names()
  
  #' Add the effect size to neutropenic subjects in the chosen column and any highly correlated columns.  
  for (j in c){
    
    index.col <- j
    
    dt <- dt %>% 
      mutate('{index.col}' := .data[[index.col]] + neutropenia*effect.size[i])
  
  }
  
  #' Sample neutropenic subjects and controls at a 1:2 ratio.
  data <- dt %>% 
    mutate(runif = runif(1000, 0, 1))
  
  neut <- data %>% 
    filter(neutropenia == 1) %>% 
    arrange(runif) %>% 
    slice(1:36) %>% 
    select(-runif)
  
  control <- data %>% 
    filter(neutropenia == 0) %>% 
    arrange(runif) %>% 
    slice(1:(nrow(neut) * 4))
  
  data <- bind_rows(neut, control)
  
  cif <- cforest(neutropenia ~ .,
                 data = data,
                 control = cforest_unbiased(ntree = 500, mtry = 10))
  
  cif.importance <- varimp(cif) %>% sort(decreasing = T)
  
  #' Re-scale the MDA values such that the sum is 1.
  scaling.factor <- 1/sum(cif.importance)
  
  cif.importance <- cif.importance * scaling.factor
  
  new.result <- data.frame(iteration = i,
                           var = index.col,
                           effect.size = effect.size[i],
                           sigma = select(data, my.col[i]) %>% unlist() %>% sd(),
                           cif.rank = which(names(cif.importance) == my.col[i]),
                           scaled.mda = cif.importance[index.col],
                           flagged.important = cif.importance[my.col[i]] > abs(min(cif.importance)))
  
  cif.results <- rbind(cif.results, new.result)
  
  saveRDS(cif.results, 
          '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/R01.aim.3.power.for.RF.models.8.20230209.rds')

}

# Combine the results -----------------------------------------------------

dir <- '//smb-main.ad.bcm.edu/genepi/TINMAN/Metabolomics/R_outputs/'

files <- list.files(dir) %>% 
  subset(str_detect(., 'R01.aim.3.power')) %>% 
  subset(str_detect(., '.rds'))

results <- tibble()

for (i in files){
  
  new.result <- readRDS(paste0(dir, i))
  
  results <- rbind(results, new.result)
  
}

plot <- ggplot(data = results, aes(x = effect.size, y = scaled.mda, color = flagged.important)) + 
 
  geom_point() + 
  
  scale_y_continuous(limits = c(0,1))

plot


tmp <- results %>% 
  mutate(effect.size.cat = factor(round(effect.size,1)))

tab <- crosstab(tmp$effect.size.cat, tmp$flagged.important, prop.r = T)$prop.row %>% 
  as_tibble() 

names(tab) <- c('effect.size.cat', 'flagged.important', 'percent')
#%>% 
  #pivot_longer(!effect.size.cat,
   #           names_to = 'result',
    #          values_to = 'percent')

tmp <- results %>% 
  select(effect.size, flagged.important) %>% 
  mutate(effect.size = factor(round(effect.size, 1))) %>% 
  filter(effect.size != '0')

plot <- ggplot(data = tmp, aes(x = effect.size, color = flagged.important, fill = flagged.important)) + 
  
  geom_bar() +
  
  labs(x = "STANDARDIZED EFFECT SIZE",
       y = 'PERCENTAGE OF MODELS IDENTIFYING VARIABLE AS IMPORTANT')

plot

# Scratch paper -----------------------------------------------------------

tin <- tin %>% 
  mutate(across(-PARENT_SAMPLE_NAME, ~ ifelse(is.na(.x), min(.x, na.rm = T)/2, .x)))

c <- cor(tin[, c(2:30,32:60)], use = 'all.obs')

mu <- tin %>% select(-1) %>% apply(2, mean) %>% log()

sigma.squared <- tin %>% 
  select(-1) %>% 
  apply(2, var) %>% 
  ifelse(. < 1, 1.01, .) %>% 
  log()

def <- defRepeat(nVars = 1000, prefix = 'm', dist = 'gamma', formula = runif(1, 0,10), variance = runif(1,0,2)) |>
  (\(.) mutate(., formula = runif(nrow(.), 0, 10), variance = runif(nrow(.), 0, 2)))()


tmp <- tin %>% 
  as.data.frame() %>% 
  select(-1) %>% 
  apply(2, mean)

apply(tmp, 2, mean)


def <- defData(varname = "x", dist = "normal", formula = 0, variance = 1, id = "cid")
dt <- genData(1000, def)

dt <- addCorData(dt, idname = "cid", mu = c(0, 0), sigma = c(2, 0.2), rho = -0.2,
                 corstr = "cs", cnames = c("a0", "a1"))

cor(dt[,2:4])

def <- defData(varname = 'var1', dist = 'normal', formula = mu[1], variance = sigma.squared[1])

for (i in seq_along(mu)){
  
  if (i == 1){
    
    def <- defData(varname = paste0('var',i), dist = 'normal', formula = mu[i], variance = sigma.squared[i])
    
  }
  
  else {
    
    def <- defData(def, varname = paste0('var',i), dist = 'normal', formula = mu[i], variance = sigma.squared[i]) 
    
  }
  
}

dt <- genData(100, def) %>% 
  mutate(across(starts_with('var'), ~ exp(.x)))


ggplot(data = dt, aes(x=exp(var3))) +
  geom_density()




tmp <- data.frame(sim = rlnorm(100, mean(log(tin$C30)), sd(log(tin$C30)) )) 

tmp2 <- data.frame(sim = rnorm(1000, 0, 1))

ggplot(data = tmp2, aes(x=exp(sim))) +
  geom_density()

#' Generate a table that defines the properties of the variables in the simulation. 
#def <- defData(varname = "C30", dist = "gamma", formula = runif(1, 0,10), variance = runif(1,0,2))

def <- defRepeat(nVars = 1000, prefix = 'm', dist = 'gamma', formula = runif(1, 0,10), variance = runif(1,0,2)) |>
  (\(.) mutate(., formula = runif(nrow(.), 0, 10), variance = runif(nrow(.), 0, 2)))()

#def <- defRepeat(nVars = 1000, prefix = 'm', dist = 'normal', formula = runif(1, 0,10), variance = runif(1,0,10)) |>
#  (\(.) mutate(., formula = runif(nrow(.), 0, 10), variance = runif(nrow(.), 0, 6)))()

#def <- defData(def, varname = "C35", dist = "normal", formula = mean(tin$C35), variance = sd(tin$C35)^2, link = 'log')

set.seed(87261)

#' Simulate data, then add neutropenia variable, with 20% positivity rate.
dd <- genData(1000, def)
dd <- dd |>
  mutate(neutropenia = rbinom(nrow(dd), 1, 0.2))





my.var <- data.frame(my.var = rnorm(100, 0, 1))

summary(my.var)

ggplot(data = my.var, aes(x = log(my.var))) + 
  geom_density()



def <- defRepeat(nVars = 100, prefix = 'm', dist = 'normal', formula = runif(1, 0,5), variance = runif(1,0,5)) |>
  (\(.) mutate(., formula = runif(nrow(.), 0, 10), variance = runif(nrow(.), 0, 2)))()

ggplot(data = dt, aes(x=var292)) + geom_density()
