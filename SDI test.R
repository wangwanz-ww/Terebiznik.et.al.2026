### Set up ####
# Load in packages
library(readr)
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(ape)
library(httr)
library(readxl)

# GitHub working directory
HH <- "https://raw.githubusercontent.com/wangwanz-ww/Terebiznik.et.al.2026/main/"

# Load in sex ratio data
df <- read.csv(paste0(HH, "sex%20ratio%20data%20final.csv")) # sex ratios
bs <- read.csv(paste0(HH, "Terebiznik_et_al_2026_ssd_list.csv")) # SSD data
bk <- read.csv(paste0(HH, "Bokony_et_al_2019_data.csv")) # bimaturism data

# calculate size dimorphism index ####
bs$SDI <- ifelse(bs$`M.cm.maxSVL/CL` > bs$`F.cm.maxSVL/CL`,
                         (-(bs$`M.cm.maxSVL/CL` /bs$`F.cm.maxSVL/CL`) + 1),
                         ((bs$`F.cm.maxSVL/CL`/bs$`M.cm.maxSVL/CL`) - 1))

# calculate bimaturism index ####
bk$matur.ratio <- ifelse(bk$matur.male > bk$matur.female,
                         (-(bk$matur.male/bk$matur.female) + 1),
                         ((bk$matur.female/bk$matur.male) - 1))

### merge data####
# all SSD & Maturation
match_indices <- match(bk$species, bs$species)
bk$SDI <- bs$SDI[match_indices]
SMatur <- as.data.frame(bk %>% 
                          filter(!is.na(matur.ratio) & !is.na(SDI)))

# TSD turtles hatchling sex ratios & SSD 
hat <- df %>% 
  filter(Life.Stage %in% c("birth.sex.ratio"),
         sex.determination %in% c("TSD"),
         model4 %in% c("turtles"))

match_indices <- match(hat$species, bs$species)
hat$SDI <- bs$SDI[match_indices]

hat <- as.data.frame(hat %>% filter(!is.na(SDI)))

### Load in phylogeny data ####
phylo <- read.tree(paste0(HH, "ultimate%20phylogeny.phy"))
inv.phylo = inverseA(phylo, nodes="ALL", scale=TRUE)

PhyloPrec <- inv.phylo$Ainv

SMatur$phylo_id <- match(SMatur$species, rownames(PhyloPrec))
hat$phylo_id <- match(hat$species, rownames(PhyloPrec))

### Create variables for residual variance analysis ####
# To check if the numbers make sense
hat$NMale <- round(hat$N*hat$Sex.Ratio)
hat$NFemale <- round(hat$N*(1-hat$Sex.Ratio))

# make variable with unique levels for when more than 1 species is in a population
## 1. bimaturism & SSD
SMatur$pop <- SMatur$population.ID
SpPerPop <- tapply(SMatur$species, list(SMatur$population.ID),
                   function(x) length(unique(x)))
MultSp <- SMatur$pop==names(SpPerPop[SpPerPop>1])
SMatur$pop[MultSp] <- paste(SMatur$pop[MultSp], SMatur$species[MultSp], sep="_")

SMatur$phylo <- SMatur$species

## 2. TSD turtles hatchling sex ratio  & SSD 
hat$pop <- hat$population.ID
SpPerPop <- tapply(hat$species, list(hat$population.ID),
                   function(x) length(unique(x)))
MultSp <- hat$pop==names(SpPerPop[SpPerPop>1])
hat$pop[MultSp] <- paste(hat$pop[MultSp], hat$species[MultSp], sep="_")

hat$phylo <- hat$species


### Create a binary dataset ####
binary_ultimate <- function(data){
  
  # filter data so it only contains values with sex ratio (double check)
  temp <- filter(data, !is.na(Sex.Ratio))
  
  # create phylo column that is identical to species column (need for model)
  temp2 <- mutate(temp, phylo = species)
  
  # get column for number of males and females, remove those w/o sample size
  temp3 <- temp2 %>% 
    mutate(Nmales = round(Sex.Ratio * N)) %>% 
    mutate(Nfemales = round((1-Sex.Ratio)*N))
  temp4 <- filter(temp3, !is.na(N))
  
  # create empty data frame to fill in
  binary_data <- data.frame()
  
  # for each row in data; convert to binary
  for(i in 1:nrow(temp4)){
    row_fill <- temp4[i, ] # get basic information
    row_fillM <- cbind(row_fill, Sex  = c(1)) # fill in with column added if male
    row_fillF <- cbind(row_fill, Sex = c(0)) # fill in with column added if female
    
    # Get number of each sex reported
    NMales <- temp4[i, "Nmales"] # number of males (rows needed)
    NFemales <- temp4[i, "Nfemales"] # number of females(rows needed)
    
    # Create new data frame
    df <- rbind(row_fillM, row_fillF) # get data frame with both ones
    
    # Contingency in case a population is reported as 100% male or female
    if(NMales == 0|NFemales == 0){
      # If males are 0
      if(NMales == 0){  # remove male row if no males in this population reported
        df <- df[-1, ]
        df_temporary <- as.data.frame(lapply(df, rep, NFemales))
      }
      
      # If Females are 0
      if(NFemales == 0){df <- df[-2,]
      df_temporary <- as.data.frame(lapply(df, rep, NMales))
      } # remove if no females in this population reported
    }
    # Create temporary data set with replicated number of males and females
    else{df_temporary <- as.data.frame(lapply(df, rep, c(NMales, NFemales)))}
    
    # Add to total data frame
    binary_data <- rbind(binary_data, df_temporary)
    
  }
  # fix population id column name
  names(binary_data)[names(binary_data) == '?..population.ID'] <- 'population.ID'
  
  # return data
  return(binary_data)
  
}

hat_binary <- binary_ultimate(hat)

### Model of correlation between SSD vs bimaturism ####
prior_SM <- list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
  R = list (V = 1, nu = 0.002))

model_SM <- MCMCglmm(SDI ~ matur.ratio, random=~species + pop + phylo_id, rcov=~matur.ratio:units, family = "gaussian", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SM, data = SMatur,nitt=50000, burnin=3000, thin=50, verbose=F)
summary(model_SM) 

# model visualization 
new_data <- data.frame(
  matur.ratio = seq(min(SMatur$matur.ratio, na.rm = TRUE), 
                    max(SMatur$matur.ratio, na.rm = TRUE), 
                    length.out = 100)
)

X <- model.matrix(~ matur.ratio, data = new_data)
preds <- X %*% t(model_SM$Sol[, 1:2]) 

new_data$pred_mean <- rowMeans(preds)
new_data$pred_lower <- apply(preds, 1, function(x) quantile(x, 0.025))
new_data$pred_upper <- apply(preds, 1, function(x) quantile(x, 0.975))

library(ggplot2)
ggplot(new_data, aes(x = matur.ratio, y = pred_mean)) +
  geom_line(color = "black", size = 1) +  #line
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.2, fill = "grey") + #CRI
  geom_point(data = SMatur, aes(x = matur.ratio, y = SDI, shape = taxon), alpha = 0.3) + # raw data
  labs(x = "Bimaturism Index",
       y = "Size Dimorphism Index",
       shape = "Taxon") +
  theme_minimal()

stats_sum <- data.frame(
  Post_Mean = colMeans(model_SM$Sol),
  Lwr_95 = HPDinterval(model_SM$Sol, prob = 0.95)[,"lower"],
  Upr_95 = HPDinterval(model_SM$Sol, prob = 0.95)[,"upper"],
  Lwr_90 = HPDinterval(model_SM$Sol, prob = 0.90)[,"lower"],
  Upr_90 = HPDinterval(model_SM$Sol, prob = 0.90)[,"upper"]
)
print(stats_sum)

### Model of correlation between TSD turtles birth sex ratio vs SSD ####

model_sdi <- MCMCglmm(Sex ~ SDI, random=~species + pop + phylo_id, rcov=~SDI:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SM, data=hat_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
summary(model_sdi)

# model visualization
new_data <- data.frame(
  SDI = seq(min(hat$SDI, na.rm = TRUE), 
            max(hat$SDI, na.rm = TRUE), 
            length.out = 100)
)

X <- model.matrix(~ SDI, data = new_data)
logits <- X %*% t(model_sdi$Sol[, 1:2]) 
preds <- plogis(logits)

new_data$pred_mean <- rowMeans(preds)
new_data$pred_lower <- apply(preds, 1, function(x) quantile(x, 0.025))
new_data$pred_upper <- apply(preds, 1, function(x) quantile(x, 0.975))

library(ggplot2)
ggplot(new_data, aes(x = SDI, y = pred_mean)) +
  geom_line(color = "black", size = 1) +  #line
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.2, fill = "grey") + #CRI
  geom_point(data = hat, aes(x = SDI, y = Sex.Ratio), alpha = 0.3) + # raw data
  labs(x = "Size Dimorphism Index",
       y = "Sex Ratio (% Male)") +
  theme_minimal()

stats_sum <- data.frame(
  Post_Mean = colMeans(model_sdi$Sol),
  Lwr_95 = HPDinterval(model_sdi$Sol, prob = 0.95)[,"lower"],
  Upr_95 = HPDinterval(model_sdi$Sol, prob = 0.95)[,"upper"],
  Lwr_90 = HPDinterval(model_sdi$Sol, prob = 0.90)[,"lower"],
  Upr_90 = HPDinterval(model_sdi$Sol, prob = 0.90)[,"upper"]
)
print(stats_sum)

