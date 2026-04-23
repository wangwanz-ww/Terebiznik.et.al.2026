### Set up ####
# Load in packages
library(readr)
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(ape)
library(INLA)
library(readxl)

# GitHub working directory
HH <- "https://raw.githubusercontent.com/wangwanz-ww/Terebiznik.et.al.2026/main/"

## Load in sex ratio data
df <- read.csv(paste0(HH, "sex%20ratio%20data%20final.csv"))

## re-level: genetic sex determination is the reference level ####
df$sex.determination <- relevel(factor(df$sex.determination), ref = "GSD")
df$SDM.Type <- relevel(factor(df$SDM.Type), ref = "GSD")
df$model4 <- relevel(factor(df$model4), ref = "GSD")

# Create categories for evolutionary history models
df$EvoH <- case_when(df$model4 == "GSD" ~"GSD",
                     df$model4 == "lizards" ~ "TSD squamata",
                     df$model4 %in% c("crocs","turtles") ~ "TSD crocoturtle")
df$EvoH <- relevel(factor(df$EvoH), ref = "GSD")

## Load in phylogeny data ####
### 1. phylogeny for sda models
phylo <- read.tree(paste0(HH, "ultimate%20phylogeny.phy"))
inv.phylo = inverseA(phylo, nodes="ALL", scale=TRUE)
#### Create precision matrix for INLA
PhyloPrec <- inv.phylo$Ainv
df$phylo_id <- match(df$species, rownames(PhyloPrec))
A <- ape::vcv.phylo(phylo)

### 2. phylogeny for conservative models (birth + adult)
phylo_conservative <- read.tree(paste0(HH, "conservative%20BA%20phylogeny.phy"))
inv.phylo.conservative <- inverseA(phylo_conservative, nodes = "ALL", scale = TRUE)
#### Create precision matrix for INLA
PhyloPrecCon <- inv.phylo.conservative$Ainv
df$con.phylo_id <- match(df$species, rownames(PhyloPrecCon))
AA <- ape::vcv.phylo(phylo_conservative)

### 3. phylogeny for conservative models (all sex ratio)
phylo_fullcon <- read.tree(paste0(HH, "conservative%20full%20phylogeny.phy"))
inv.phylo.fullcon <- inverseA(phylo_fullcon, nodes = "ALL", scale = TRUE)

#### Create precision matrix for INLA
PhyloPrecFC <- inv.phylo.fullcon$Ainv
df$FC.phylo_id <- match(df$species, rownames(PhyloPrecFC))
AAA <- ape::vcv.phylo(phylo_fullcon)

## Create variables for residual variance analysis ####
### 1. To check if the numbers make sense
df$NMale <- round(df$N*df$Sex.Ratio)
df$NFemale <- round(df$N*(1-df$Sex.Ratio))
Diff <- df$N*df$Sex.Ratio - df$NMale

### 2. Create a variable to model extra residual variation in TSD data
df$TSDvar <- c(NA, 1)[1+(df$sex.determination=="TSD")]
df$TSDvar[!is.na(df$TSDvar)] <- cumsum(df$TSDvar[!is.na(df$TSDvar)])

### 3. make variable with unique levels for when more than 1 species is in a population
df$pop <- df$population.ID
SpPerPop <- tapply(df$species, list(df$population.ID),
                   function(x) length(unique(x)))
MultSp <- df$pop==names(SpPerPop[SpPerPop>1])
df$pop[MultSp] <- paste(df$pop[MultSp], df$species[MultSp], sep="_")

## Create subset of data for different models ####
### 1. full conservative data (birth + juvenile + adult)
conservative_full <- filter(df, Conservative == "Y")

### 2. conservative data without juvenile data (birth + adult)
conservative <- df %>% 
  filter(Conservative == "Y",
         Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"))

### 3. conservative data without TSD lizards & juvenile data (birth + adult, GSD + crocs + turtles)
conservative_CT <- df %>% 
  filter(Conservative == "Y",
         Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"),
         model4 != "lizards")

### 4. full sda (sexual dimorphic assumed data added) data (birth + juvenile + adult)
sda_full <- df

### 5. sda (sexual dimorphic assumed data added) data without juvenile data (birth + adult)
sda <- df %>% 
  filter(Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"))

### 6. sda data without TSD lizards (birth + juvenile + adult, GSD + crocs + turtles)
sda_CT <- df %>% 
  filter(model4 != "lizards")

## Create a binary dataset ####
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





## Making data into binary format####
### 1. conservative full data
conservative_binary_full <- binary_ultimate(conservative_full)
### 2. conservative B+A data 
conservative_binary <- binary_ultimate(conservative) 
### 3. conservative B+A data without TSD lizards
conservative_CT_binary <- binary_ultimate(conservative_CT)
### 4. sda full data
sdaf_binary <- binary_ultimate(sda_full)
### 5. sda B+A data 
sda_binary <- binary_ultimate(sda)
### 6. sda data without TSD lizards
sdact_binary <- binary_ultimate(sda_CT)


### Model 1: single mechanism hypothesis #######
prior_SDM <-  list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
  R = list (V = 1, nu = 0.002))

# conservative model 1.1: conservative B+A data
model_conservative1.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + con.phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_SDM, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative1.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_conservative1.1.RData")

# conservative model 2.1: conservative B+A data without TSD lizards
model_conservative2.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + con.phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_SDM, data=conservative_CT_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative2.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/BAmodel_conservative2.1.RData")

# conservative model 3.1: conservative full data (B+J+A)
model_conservative3.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + FC.phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo.fullcon$Ainv), prior=prior_SDM, data=conservative_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative3.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/FCmodel_conservative3.1.RData")


# sda model 1.1: sda B+A data 
model_sda<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SDM, data= sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda1.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.1.RData")

# sda model 2.1: sda full data
model_sda2.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SDM, data= sdaf_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda2.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda2.1.RData")

# sda model 3.1: sda data without TSD lizards
model_sda3.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~species + pop + TSDvar + phylo_id, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SDM, data= sdact_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda3.1, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda3.1.RData")

### Model 2: evolutionary history model ####
prior_phyloLH <- list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
  R = list (V = 1, nu = 0.002))

#conservative model 1.2: conservative B+A data
model_conservative1.2 <- MCMCglmm(Sex ~ Life.Stage * EvoH , random=~species + pop + TSDvar + con.phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_phyloLH, data=conservative_binary, nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative1.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_conservative1.2.1.RData")
summary(model_conservative1.2.1)

#conservative model 2.2: conservative B+A data without TSD lizards
model_conservative2.2 <- MCMCglmm(Sex ~ Life.Stage * EvoH, random=~species + pop + TSDvar + con.phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_phyloLH, data=conservative_CT_binary, nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative2.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/BAmodel_conservative2.2.RData")

# conservative model 3.2: conservative full data (B+J+A)
model_conservative3.2<- MCMCglmm(Sex ~ Life.Stage * EvoH , random=~species + pop + TSDvar + FC.phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo.fullcon$Ainv), prior=prior_phyloLH, data=conservative_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative3.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/FCmodel_conservative3.2.RData")


# sda model 1.2: sda B+A data 
model_sda1.2 <- MCMCglmm(Sex ~ Life.Stage * EvoH , random=~species + pop + TSDvar + phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_phyloLH, data=sda_binary, nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda1.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.2.1.RData")

# sda model 2.2: sda full data
model_sda2.2 <- MCMCglmm(Sex ~ Life.Stage * EvoH , random=~species + pop + TSDvar + phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_phyloLH, data=sdaf_binary, nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda2.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda2.2.1.RData")

# sda model 3.2: sda data without TSD lizards
model_sda3.2 <- MCMCglmm(Sex ~ Life.Stage * EvoH , random=~species + pop + TSDvar + phylo_id, rcov=~EvoH:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_phyloLH, data=sdact_binary, nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda3.2, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda3.2.1.RData")

### Model 3: TSD type model ####
prior_TSDtype<- list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
  R = list (V = 1, nu = 0.002))

#conservative model 1.3: conservative B+A data
model_conservative1.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + con.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_TSDtype, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative1.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_conservative1.3.RData")

#conservative model 2.3: conservative B+A data without TSD lizards
model_conservative2.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + con.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_TSDtype, data=conservative_CT_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative2.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/BAmodel_conservative2.3.RData")

#conservative model 3.3: conservative full data (B+J+A)
model_conservative3.3<- MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + FC.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.fullcon$Ainv), prior=prior_TSDtype, data=conservative_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative3.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/FCmodel_conservative3.3.RData")


# sda model 1.3: sda B+A data 
model_sda1.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda1.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.3.RData")

# sda model 2.3: sda full data
model_sda2.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sdaf_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda2.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda2.3.RData")

# sda model 3.3: sda data without TSD lizards
model_sda3.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sdact_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda3.3, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda3.3.RData")



### Model 4: Taxonomic Model ####
prior_Tax<- list(
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
           G4 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)), 
  R = list (V = 1, nu = 0.002))

# conservative model 1.4: conservative B+A data
model_conservative1.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + con.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_Tax, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative1.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_conservative1.4.RData")

# conservative model 2.4: conservative B+A data without TSD lizards
model_conservative2.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + con.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_Tax, data=conservative_CT_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative2.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/BAmodel_conservative2.4.RData")

# conservative model 3.4: conservative full data (B+J+A)
model_conservative3.4<- MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + FC.phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.fullcon$Ainv), prior=prior_Tax, data=conservative_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_conservative3.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/FCmodel_conservative3.4.RData")


# sda model 1.4: sda B+A data 
model_sda1.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_Tax, data=sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda1.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.4.RData")

# sda model 2.4: sda full data
model_sda2.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_Tax, data=sdaf_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda2.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda2.4.RData")

# sda model 3.4: sda data without TSD lizards
model_sda3.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~species + pop + TSDvar + phylo_id, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_Tax, data=sdact_binary,nitt=50000, burnin=3000, thin=50, verbose=F)
save(model_sda3.4, file = "/Users/annettew/Documents/Grad School/Mariel's paper/model_sda3.4.RData")

### Best Model ####
### Identifying Best Model ####

# conservative models
### 1. conservative B+A data
output_conservative1 <- model.sel(model_conservative1.1, model_conservative1.2, model_conservative1.3, model_conservative1.4, rank = "DIC")
output_conservative1

### 2. conservative B+A data without TSD lizards
output_conservative2 <- model.sel(model_conservative2.1, model_conservative2.2, model_conservative2.3, model_conservative2.4, rank = "DIC")
output_conservative2

### 3. conservative full data (B+J+A)
output_conservative3 <- model.sel(model_conservative3.1, model_conservative3.2, model_conservative3.3, model_conservative3.4, rank = "DIC")
output_conservative3


# sda models
### 1. sda B+A data 
output_sda1 <- model.sel(model_sda1.1, model_sda1.2, model_sda1.3, model_sda1.4,rank = "DIC")
output_sda1 
### 2. sda full data
output_sda2 <- model.sel(model_sda2.1, model_sda2.2, model_sda2.3, model_sda2.4,rank = "DIC")
output_sda2

### 3. sda data without TSD lizards
output_sda3 <- model.sel(model_sda3.1, model_sda3.2.1, model_sda3.3, model_sda3.4,rank = "DIC")
output_sda3



### Variance Tests for 4 conservative models ####

## Test for the contribution of each random effect
# PC prior for variances. 
pc.prec = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
# Create the formula: 
# The base formula has the random effects.
# The different models specify different fixed effects.
BaseForm <- NMale ~ Life.Stage +
  f(species, model="iid", hyper=pc.prec) +
  f(pop, model="iid", hyper=pc.prec) +
  f(TSDvar, model="iid", hyper=pc.prec) +
  f(FC.phylo_id, model = "generic0", Cmatrix = PhyloPrec, hyper=pc.prec)
formulae <- list(
  SM = list(form=update(BaseForm, . ~ . + Life.Stage * sex.determination),
            name="Single Mechanism Hypothesis"),
  EH = list(form=update(BaseForm, . ~ . + Life.Stage * EvoH),
            name="Supertaxa Model"),
  TSD = list(form=update(BaseForm, . ~ . + Life.Stage * SDM.Type),
             name="TSD Type Model"),
  Tax = list(form=update(BaseForm, . ~ . + Life.Stage * model4),
             name="Taxonomic Model"))

model_conservative <- lapply(formulae, function(form, data) {
  res <- inla(form$form, family="binomial", Ntrials = N, data=data,
              control.compute = list(dic=TRUE,waic=TRUE, config = TRUE))
  res$modelname <- form$name
  res
}, data = conservative_full) # change data for each subset: conservative/conservative_full/conservative_CT.

# Plots for Variance Terms
PlotVariances <- function(mod, axislab=TRUE, ...) {
  CalcInv <- function(m) inla.tmarginal(function(x) 1/x, m)
  m.var <- invisible(lapply(mod$internal.marginals.hyperpar, CalcInv))
  Vars.l <- lapply(m.var, function(m) unlist(inla.zmarginal(m, silent = TRUE)))
  Vars <- t(as.data.frame(Vars.l))
  plot(Vars[,"mean"], 1:nrow(Vars), xlim=c(0, max(Vars[,"quant0.975"])),
       yaxt="n", ylab="", xlab="", main=mod$modelname, ...)
  segments(Vars[,"quant0.025"], 1:nrow(Vars), Vars[,"quant0.975"], 1:nrow(Vars))
  segments(Vars[,"quant0.25"], 1:nrow(Vars), Vars[,"quant0.75"], 1:nrow(Vars),
           lwd=4, lend=3)
  if(axislab) axis(2, at=1:nrow(Vars), las=1,
                   labels=gsub("Log.precision.for.", "", rownames(Vars)))
}
par(mfrow=c(2,2), mar=c(4,4,4,1), oma=c(2,1,0,0))
invisible(lapply(model_conservative, PlotVariances))
mtext("Variance Term", 1, outer=TRUE)

SampTSDRatio <- function(mod, n) {
  samp <- inla.hyperpar.sample(n=n, result =mod)
  ratio <- apply(samp, 1, function(x) {
    1+x["Precision for pop"]/x["Precision for TSDvar"]
  })
  ratio
}
TSDRatios <- invisible(lapply(model_conservative, SampTSDRatio, n=1e4))
par(mfrow=c(2,2), mar=c(3,2,4,1), oma=c(2,2,0,0))
invisible(sapply(names(TSDRatios), function(name, rats) {
  ratio <- rats[[name]]
  plot(density(ratio), main=name, xlab="", ylab="")
}, rats=TSDRatios))
mtext("Variance Ratio", 1, outer=TRUE)
mtext("Density", 2, outer=TRUE)

# Rank models with other parameters
GetICs <- function(mod) {
  c(DIC=mod$dic$dic, wAIC=mod$waic$waic, MargDev = -2*mod$mlik[1])
}
ICs <- data.frame(invisible(lapply(model_conservative, GetICs)))

knitr::kable(sweep(ICs, 1, apply(ICs, 1, min), "-"), digits = 2)

knitr::kable(summary(model_conservative$SM)$fixed[,1:5])
knitr::kable(summary(model_conservative$EH)$fixed[,1:5])
knitr::kable(summary(model_conservative$TSD)$fixed[,1:5])
knitr::kable(summary(model_conservative$Tax)$fixed[,1:5])

### Variance Tests for 4 sda models ####

### Test for the contribution of each random effect 
# PC prior for variances. 
pc.prec = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))
# Create the formula. 
# The base formula has the random effects.
# The different models specify different fixed effects.
BaseForm <- NMale ~ Life.Stage +
  f(species, model="iid", hyper=pc.prec) +
  f(pop, model="iid", hyper=pc.prec) +
  f(TSDvar, model="iid", hyper=pc.prec) +
  f(phylo_id, model = "generic0", Cmatrix = PhyloPrec, hyper=pc.prec)
formulae <- list(
  SM = list(form=update(BaseForm, . ~ . + Life.Stage * sex.determination),
            name="A. Single Mechanism Model"),
  EH = list(form=update(BaseForm, . ~ . + Life.Stage * EvoH),
            name="B. Evolutionary History Model"),
  TSD = list(form=update(BaseForm, . ~ . + Life.Stage * SDM.Type),
             name="C. TSD Type Model"),
  Tax = list(form=update(BaseForm, . ~ . + Life.Stage * model4),
             name="D. Taxonomic Model"))

model_sda <- lapply(formulae, function(form, data) {
  res <- inla(form$form, family="binomial", Ntrials = N, data=data,
              control.compute = list(dic=TRUE,waic=TRUE, config = TRUE))
  res$modelname <- form$name
  res
}, data = sda_full)# change data for each subset: conservative/conservative_full/conservative_CT.

# Plots for Variance Terms
PlotVariances <- function(mod, axislab=TRUE, ...) {
  CalcInv <- function(m) inla.tmarginal(function(x) 1/x, m)
  m.var <- invisible(lapply(mod$internal.marginals.hyperpar, CalcInv))
  Vars.l <- lapply(m.var, function(m) unlist(inla.zmarginal(m, silent = TRUE)))
  Vars <- t(as.data.frame(Vars.l))
  plot(Vars[,"mean"], 1:nrow(Vars), xlim=c(0, max(Vars[,"quant0.975"])),
       yaxt="n", ylab="", xlab="", main=mod$modelname, ...)
  segments(Vars[,"quant0.025"], 1:nrow(Vars), Vars[,"quant0.975"], 1:nrow(Vars))
  segments(Vars[,"quant0.25"], 1:nrow(Vars), Vars[,"quant0.75"], 1:nrow(Vars),
           lwd=4, lend=3)
  if(axislab) axis(2, at=1:nrow(Vars), las=1,
                   labels=gsub("Log.precision.for.", "", rownames(Vars)))
}
par(mfrow=c(2,2), mar=c(4,4,4,1), oma=c(2,1,0,0))
invisible(lapply(model_sda, PlotVariances))
mtext("Variance Term", 1, outer=TRUE)

SampTSDRatio <- function(mod, n) {
  samp <- inla.hyperpar.sample(n=n, result =mod)
  ratio <- apply(samp, 1, function(x) {
    1+x["Precision for pop"]/x["Precision for TSDvar"]
  })
  ratio
}
TSDRatios <- invisible(lapply(model_sda, SampTSDRatio, n=1e4))
par(mfrow=c(2,2), mar=c(3,2,4,1), oma=c(2,2,0,0))
invisible(sapply(names(TSDRatios), function(name, rats) {
  ratio <- rats[[name]]
  plot(density(ratio), main=name, xlab="", ylab="")
}, rats=TSDRatios))
mtext("Variance Ratio", 1, outer=TRUE)
mtext("Density", 2, outer=TRUE)

# Rank models with other parameters
GetICs <- function(mod) {
  c(DIC=mod$dic$dic, wAIC=mod$waic$waic, MargDev = -2*mod$mlik[1])
}
ICs <- data.frame(invisible(lapply(model_sda, GetICs)))

knitr::kable(sweep(ICs, 1, apply(ICs, 1, min), "-"), digits = 2)
knitr::kable(summary(model_sda$SM)$fixed[,1:5])
knitr::kable(summary(model_sda$EH)$fixed[,1:5])
knitr::kable(summary(model_sda$TSD)$fixed[,1:5])
knitr::kable(summary(model_sda$Tax)$fixed[,1:5])