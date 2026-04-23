# This file will sample model data sets so we don't have to do it every single time 
# and the R files are more manageable

#### Set up ####
# Load in packages
library(readr)
library(bayestestR)
library(tidyverse)

# Load in Models
load("/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.1.RData") # SMH
load("/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.2.RData") # EvoH
load("/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.3.RData") # TSD
load("/Users/annettew/Documents/Grad School/Mariel's paper/model_sda1.4.RData") # Tax

# relevant function
inv.logit <- function(x){
  result <- exp(x)/(1+exp(x))
  return(result)
} 

## Birth + Adult####
#### Model 1: SMH #####
smh <- reshape2::melt(as.matrix(model_conservative1.1$Sol)) # change model name as you proceed with different models
smh <- smh %>% 
  spread(Var2, value)
#Intercept (AGSD): Sample from posterior estimates
# Get Values for B0, b1, b2, b3
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
for(i in 1:nrow(smh)){
  # get new values
  b0 <- smh[i,2]
  b1 <- smh[i,3]
  b2 <- smh[i,4]
  b3 <- smh[i,5]
  # save in vectors
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
}

# Sample adult GSD (intercept)
AGSD <-c()
x1 <- 0
x2 <- 0
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x1*x2
  # y = p(male)
  # B0 = adult, tsd
  # X1 = Life stage
  # x1 = 0; adult
  # x1 = 1; birth
  # X2 = SDM
  # x2 = 0; GSD
  #x2 = 1; TSD
  AGSD <- c(AGSD, y)
}

# BGSD: Sample from posterior estimates
BGSD <- c()
x1 <- 1
x2 <- 0
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x1*x2
  # y = p(male)
  # B0 = adult, tsd
  # X1 = Life stage
  # x1 = 0; adult
  # x1 = 1; birth
  # X2 = SDM
  # x2 = 0; GSD
  #x2 = 1; TSD
  BGSD <- c(BGSD, y)
}

# ATSD: Sample from posterior estimates
ATSD <- c()
x1 <- 0
x2 <- 1
for(i in 1:nrow(smh)){
  
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x1*x2
  # y = p(male)
  # B0 = adult, tsd
  # X1 = Life stage
  # x1 = 0; adult
  # x1 = 1; birth
  # X2 = SDM
  # x2 = 0; GSD
  #x2 = 1; TSD
  ATSD <- c(ATSD, y)
}

# BTSD: Sample from posterior estimates
BTSD <- c()
x1 <- 1
x2 <- 1
for(i in 1:nrow(smh)){
  
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x1*x2
  # y = p(male)
  # B0 = adult, tsd
  # X1 = Life stage
  # x1 = 0; adult
  # x1 = 1; birth
  # X2 = SDM
  # x2 = 0; GSD
  #x2 = 1; TSD
  
  BTSD <- c(BTSD, y)
}

# Combine into single data set
SMH_model <- cbind(AGSD, BGSD, ATSD, BTSD)
SMH_model <- as.data.frame(SMH_model)




# make it appropriate format with necessary columns
smh_model <- SMH_model %>%  
  mutate_all(inv.logit) %>% # transform all onto probability space 
  gather(Group) %>% 
  mutate(Life.Stage = case_when( # create new column for Life Stage
    Group == "AGSD" ~ "adult.sex.ratio",
    Group == "ATSD" ~ "adult.sex.ratio",
    Group == "BTSD" ~"birth.sex.ratio",
    Group == "BGSD" ~"birth.sex.ratio")) %>% 
  mutate(sex.determination = case_when( # create new for sex determination
    Group == "AGSD" ~ "GSD",
    Group == "ATSD" ~ "TSD",
    Group == "BTSD" ~"TSD",
    Group == "BGSD" ~"GSD")) 


# create summary data
smh_sum <- smh_model %>% 
  group_by(Group, Life.Stage, sex.determination) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  # get groups by sex determination so can draw line between them
  mutate(grp = case_when(
    sex.determination == "TSD" ~ 1,
    sex.determination == "GSD" ~ 2
  )) %>% 
  # add life stage with a bit extra for 1 side to make it fit with violin plots
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~2
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 + 0.225*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + -.225*(grp == "2"))

#### Model 2: SuperTaxa ####

sut <- reshape2::melt(as.matrix(model_conservative1.2.1$Sol))
sut <- sut %>% 
  spread(Var2, value)

# Get values for B0-15
B0 <- c() #A_GSD
B1 <- c() #B_GSD
B2 <- c() #A_CT
B3 <- c() #A_SQ
B4 <- c() #B_CT
B5 <- c() #B_SQ


for(i in 1:nrow(sut)){
  b0 <- sut[i,2]
  b1 <- sut[i,3]
  b2 <- sut[i,4]
  b3 <- sut[i,5]
  b4 <- sut[i,6]
  b5 <- sut[i,7]
  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)


}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = supertaxa (ct); x2 = 0 (GSD), x2 = 1 (ct)
# x3 = supertaxa (sq); x3 = 0 (GSD), x3 = 1 (sq)


# Intercept: A_GSD
A_GSD <- c()

x1 <- 0; x2 <- 0; x3 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_GSD <- c(A_GSD, y)
}

# B_GSD
B_GSD <- c()

x1 <- 1; x2 <- 0; x3 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_GSD <- c(B_GSD, y)
}


# A_CT
A_CT <- c()

x1 <- 0; x2 <- 1; x3 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_CT <- c(A_CT, y)
}


# A_SQ
A_SQ <- c()

x1 <- 0; x2 <- 0; x3 <- 1
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_SQ <- c(A_SQ, y)
}

# B_CT
B_CT <- c()

x1 <- 1; x2 <- 1; x3 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_CT <- c(B_CT, y)
}


# B_SQ
B_SQ <- c()

x1 <- 1; x2 <- 0; x3 <- 1
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_SQ <- c(B_SQ, y)
}


# Putting it all together
# put together into single data frame
SUT_model <- cbind(A_GSD,B_GSD,A_CT,A_SQ,B_CT,B_SQ)
SUT_model <- as.data.frame(SUT_model)

# modify and add relevant columns
sut_model <- SUT_model %>% 
  # get inv logit values to put them on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio")) %>% 
  # add taxonomic group
  mutate(EvoH = case_when( # create new column for supertaxa
    substr(Group, nchar(Group) - 1, nchar(Group)) == "SD" ~ "GSD",
    substr(Group, nchar(Group) - 1, nchar(Group)) == "CT" ~ "TSD crocoturtle",
    substr(Group, nchar(Group) - 1, nchar(Group)) == "SQ" ~ "TSD squamata")) 

sut_sum <- sut_model %>% 
  group_by(Group, Life.Stage, EvoH) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    EvoH == "GSD" ~ 0,
    EvoH == "TSD crocoturtle" ~ 1,
    EvoH == "TSD squamata" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) 



#### Model 3: TSD types ####

tsd <- reshape2::melt(as.matrix(model_conservative1.3$Sol))
tsd <- tsd %>% 
  spread(Var2, value)


# Get values for B0-8
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()


for(i in 1:nrow(tsd)){
  b0 <- tsd[i,2]
  b1 <- tsd[i,3]
  b2 <- tsd[i,4]
  b3 <- tsd[i,5]
  b4 <- tsd[i,6]
  b5 <- tsd[i,7]

  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)

}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = SDM (Ia); x2 = 0 (GSD), x2 = 1 (Ia)
# x3 = SDM (II); x3 = 0 (GSD), x3 = 1 (II)
  
# Intercept: A_GSD
A_GSD <- c()

x1 <- 0; x2 <- 0; x3 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_GSD <- c(A_GSD, y)
}

# B_GSD
B_GSD <- c()

x1 <- 1; x2 <- 0; x3 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_GSD <- c(B_GSD, y)
}

# A_Ia
A_Ia <- c()

x1 <- 0; x2 <- 1; x3 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_Ia <- c(A_Ia, y)
}

# A_II
A_II <- c()

x1 <- 0; x2 <- 0; x3 <- 1
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  A_II <- c(A_II, y)
}

# B_Ia
B_Ia <- c()

x1 <- 1; x2 <- 1; x3 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_Ia <- c(B_Ia, y)
}

# B_II
B_II <- c()

x1 <- 1; x2 <- 0; x3 <- 1
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + 
    B4[i]*x1*x2 + B5[i]*x1*x3 
  B_II <- c(B_II, y)
}


# put together into single data frame
tsd_model <- cbind(A_GSD,B_GSD,A_Ia, A_II, B_Ia,B_II)
tsd_model <- as.data.frame(tsd_model)

# modify and add relevant columns
tsd_model <- tsd_model %>% 
  # get inv logit values to putthem on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio")) %>% 
  # add sex determination
  mutate(SDM.Type = case_when( # create new column for SDM types
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 1, nchar(Group))  == "Ia" ~ "TSD Ia",
    substr(Group, nchar(Group) - 1, nchar(Group))  == "II" ~ "TSD II"))

tsd_sum <- tsd_model %>% 
  group_by(Group, Life.Stage, SDM.Type) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    SDM.Type == "TSD Ia" ~ 0,
    SDM.Type == "TSD II" ~ 1,
    SDM.Type == "GSD" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "2"))


#### Model 4: Taxonomic ####

tax <- reshape2::melt(as.matrix(model_conservative1.4$Sol))
tax <- tax %>% 
  spread(Var2, value)

# Get values for B0,B1, B2, B3 ,B4, B5, B6, B7
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()
B6 <- c()
B7 <- c()

for(i in 1:nrow(tax)){
  b0 <- tax[i,2]
  b1 <- tax[i,3]
  b2 <- tax[i,4]
  b3 <- tax[i,5]
  b4 <- tax[i,6]
  b5 <- tax[i,7]
  b6 <- tax[i,8]
  b7 <- tax[i,9]

  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)
  B6 <- c(B6, b6)
  B7 <- c(B7, b7)

}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = model4(crocs); x2 = 0 (GSD), x2 = 1 (crocs)
# x3 = model4(lizards); x3 = 0 (GSD), x3 = 1 (lizards)
# x4 = model4(turtles); x4 = 0 (GSD), x4 = 1 (turtles)

# Intercept: AdGSD
AdGSD <- c()

x1 <- 0; x2 <- 0; x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  AdGSD <- c(AdGSD, y)
}

# BrGSD
BrGSD <- c()

x1 <- 1; x2 <- 0; x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  BrGSD <- c(BrGSD, y)
}

# ACrocs
AdCrocs <- c()

x1 <- 0; x2 <- 1; x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  AdCrocs <- c(AdCrocs, y)
}

# Aliz
Adliz <- c()

x1 <- 0; x2 <- 0; x3 <- 1; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  Adliz <- c(Adliz, y)
}

# Atur
Adtur <- c()

x1 <- 0; x2 <- 0; x3 <- 0; x4 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  Adtur <- c(Adtur, y)
}

# BCrocs
BrCrocs <- c()

x1 <- 1; x2 <- 1; x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  BrCrocs <- c(BrCrocs, y)
}


# Brliz
Brliz <- c()

x1 <- 1; x2 <- 0; x3 <- 1; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  Brliz <- c(Brliz, y)
}

# Btur
Brtur <- c()

x1 <- 1; x2 <- 0; x3 <- 0; x4 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + 
    B5[i]*x1*x2 + B6[i]*x1*x3 + B7[i]*x1*x4 
  Brtur <- c(Brtur, y)
}

# Putting it all together
# put together into single data frame
tax_model <- cbind(AdGSD, BrGSD, AdCrocs, Adliz, Adtur, BrCrocs,Brliz, Brtur)
tax_model <- as.data.frame(tax_model)

# modify and add relevant columns
tax_model <- tax_model %>% 
  # get inv logit values to put them on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio")) %>% 
  # add taxonomic group
  mutate(model4 = case_when( # create new column for model4 column
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 4, nchar(Group))  == "Crocs" ~ "crocs",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "liz" ~ "lizards",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "tur" ~ "turtles"))

tax_sum <- tax_model %>% 
  group_by(Group, Life.Stage, model4) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    model4 == "GSD" ~ 0,
    model4 == "crocs" ~ 1,
    model4 == "lizards" ~ 2,
    model4 == "turtles" ~ 3
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.1*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0.1*(grp == "3"))


## Birth + Juvenile + Adult####
#### full Model 1: SMH #####
smh <- reshape2::melt(as.matrix(model_conservative3.1$Sol))
smh <- smh %>% 
  spread(Var2, value)
# Get values for B0,B1, B2, B3 ,B4, B5, B6, B7
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()

for(i in 1:nrow(smh)){
  b0 <- smh[i,2]
  b1 <- smh[i,3]
  b2 <- smh[i,4]
  b3 <- smh[i,5]
  b4 <- smh[i,6]
  b5 <- smh[i,7]

  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)

}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = Lifstage (juvenile); x2 = 0 (adult), x2 = 1 (juvenile)
# x3 = SDM; x3 = 0 (GSD), x3 = 1 (TSD)

# Intercept: A_GSD
A_GSD <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  A_GSD <- c(A_GSD, y)
}

# B_GSD
B_GSD <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  B_GSD <- c(B_GSD, y)
}

# J_GSD
J_GSD <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  J_GSD <- c(J_GSD, y)
}

# A_TSD
A_TSD <- c()

x1 <- 0; x2 <- 0; 
x3 <- 1
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  A_TSD <- c(A_TSD, y)
}

# B_TSD
B_TSD <- c()

x1 <- 1; x2 <- 0; 
x3 <- 1
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  B_TSD <- c(B_TSD, y)
}

# J_TSD
J_TSD <- c()

x1 <- 0; x2 <- 1; 
x3 <- 1
for(i in 1:nrow(smh)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x1*x3 + B5[i]*x2*x3 
  
  J_TSD <- c(J_TSD, y)
}


# Combine them all together into single data set
SMH_model <- cbind(A_GSD, B_GSD, J_GSD, A_TSD, B_TSD, J_TSD)
SMH_model <- as.data.frame(SMH_model)

# modify and add relevant columns
smh_model <- SMH_model %>% 
  # get inv logit values to putthem on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column N with number of individuals
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio",
    substr(Group, 1,1) == "J" ~ "juvenile.sex.ratio")) %>% 
  # add sex determination
  mutate(sex.determination = case_when( # create new column N with number of individuals
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "TSD" ~ "TSD"))

smh_sum <- smh_model %>% 
  group_by(Group, Life.Stage, sex.determination) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    sex.determination == "TSD" ~ 1,
    sex.determination == "GSD" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 + 0.225*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + -.225*(grp == "2"))

#### full Model 2: SuperTaxa ####

sut <- reshape2::melt(as.matrix(model_conservative3.2$Sol))
sut <- sut %>% 
  spread(Var2, value)

# Get values for B0-15
B0 <- c() #A_GSD
B1 <- c() #B_GSD
B2 <- c() #J_GSD
B3 <- c() #A_CT
B4 <- c() #A_SQ
B5 <- c() #B_CT
B6 <- c() #J_CT
B7 <- c() #B_SQ
B8 <- c() #J_SQ

for(i in 1:nrow(sut)){
  b0 <- sut[i,2]
  b1 <- sut[i,3]
  b2 <- sut[i,4]
  b3 <- sut[i,5]
  b4 <- sut[i,6]
  b5 <- sut[i,7]
  b6 <- sut[i,8]
  b7 <- sut[i,9]
  b8 <- sut[i,10]
  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)
  B6 <- c(B6, b6)
  B7 <- c(B7, b7)
  B8 <- c(B8, b8)
  
}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = Lifstage (juvenile); x2 = 0 (adult), x2 = 1 (juvenile)
# x3 = supertaxa (ct); x3 = 0 (GSD), x3 = 1 (ct)
# x4 = supertaxa (sq); x4 = 0 (GSD), x4 = 1 (sq)


# Intercept: A_GSD
A_GSD <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  A_GSD <- c(A_GSD, y)
}

# B_GSD
B_GSD <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  B_GSD <- c(B_GSD, y)
}

#J_GSD
J_GSD <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  J_GSD <- c(J_GSD, y)
}

# A_CT
A_CT <- c()

x1 <- 0; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  A_CT <- c(A_CT, y)
}


# A_SQ
A_SQ <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  A_SQ <- c(A_SQ, y)
}

# B_CT
B_CT <- c()

x1 <- 1; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  B_CT <- c(B_CT, y)
}

#J_CT
J_CT <- c()

x1 <- 0; x2 <- 1; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  J_CT <- c(J_CT, y)
}

# B_SQ
B_SQ <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  B_SQ <- c(B_SQ, y)
}

#J_SQ
J_SQ <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(sut)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  J_SQ <- c(J_SQ, y)
}


# Putting it all together
# put together into single data frame
SUT_model <- cbind(A_GSD,B_GSD,J_GSD, A_CT,A_SQ,B_CT,J_CT,B_SQ,J_SQ)
SUT_model <- as.data.frame(SUT_model)

# modify and add relevant columns
sut_model <- SUT_model %>% 
  # get inv logit values to put them on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio",
    substr(Group, 1,1) == "J" ~ "juvenile.sex.ratio")) %>% 
  # add taxonomic group
  mutate(EvoH = case_when( # create new column for supertaxa
    substr(Group, nchar(Group) - 1, nchar(Group)) == "SD" ~ "GSD",
    substr(Group, nchar(Group) - 1, nchar(Group)) == "CT" ~ "TSD crocoturtle",
    substr(Group, nchar(Group) - 1, nchar(Group)) == "SQ" ~ "TSD squamata")) 

sut_sum <- sut_model %>% 
  group_by(Group, Life.Stage, EvoH) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    EvoH == "GSD" ~ 0,
    EvoH == "TSD crocoturtle" ~ 1,
    EvoH == "TSD squamata" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) 



#### full Model 3: TSD types ####

tsd <- reshape2::melt(as.matrix(model_conservative3.3$Sol))
tsd <- tsd %>% 
  spread(Var2, value)


# Get values for B0-8
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()
B6 <- c()
B7 <- c()
B8 <- c()

for(i in 1:nrow(tsd)){
  b0 <- tsd[i,2]
  b1 <- tsd[i,3]
  b2 <- tsd[i,4]
  b3 <- tsd[i,5]
  b4 <- tsd[i,6]
  b5 <- tsd[i,7]
  b6 <- tsd[i,8]
  b7 <- tsd[i,9]
  b8 <- tsd[i,10]
  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)
  B6 <- c(B6, b6)
  B7 <- c(B7, b7)
  B8 <- c(B8, b8)
}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = Lifstage (juvenile); x2 = 0 (adult), x2 = 1 (juvenile)
# x3 = SDM (Ia); x3 = 0 (GSD), x3 = 1 (Ia)
# x4 = SDM (II); x4 = 0 (GSD), x4 = 1 (II)

# posterior means
mean(B0, na.rm = TRUE)
quantile(B0, probs = 0.025, na.rm = TRUE)
quantile(B0, probs = 0.975, na.rm = TRUE)



# Intercept: A_GSD
A_GSD <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  A_GSD <- c(A_GSD, y)
}

# B_GSD
B_GSD <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  B_GSD <- c(B_GSD, y)
}

#J_GSD
J_GSD <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  J_GSD <- c(J_GSD, y)
}

# A_Ia
A_Ia <- c()

x1 <- 0; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  A_Ia <- c(A_Ia, y)
}

# A_II
A_II <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  A_II <- c(A_II, y)
}

# B_Ia
B_Ia <- c()

x1 <- 1; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  B_Ia <- c(B_Ia, y)
}

#J_Ia
J_Ia <- c()

x1 <- 0; x2 <- 1; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  J_Ia <- c(J_Ia, y)
}

# B_II
B_II <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  B_II <- c(B_II, y)
}

#J_II
J_II <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tsd)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  
  J_II <- c(J_II, y)
}

# put together into single data frame
tsd_model <- cbind(A_GSD,B_GSD,J_GSD, A_Ia, A_II,  B_Ia,J_Ia, B_II,J_II)
tsd_model <- as.data.frame(tsd_model)

# modify and add relevant columns
tsd_model <- tsd_model %>% 
  # get inv logit values to putthem on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio",
    substr(Group, 1,1) == "J" ~ "juvenile.sex.ratio")) %>% 
  # add sex determination
  mutate(SDM.Type = case_when( # create new column for SDM types
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 1, nchar(Group))  == "Ia" ~ "TSD Ia",
    substr(Group, nchar(Group) - 1, nchar(Group))  == "II" ~ "TSD II"))

tsd_sum <- tsd_model %>% 
  group_by(Group, Life.Stage, SDM.Type) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    SDM.Type == "TSD Ia" ~ 0,
    SDM.Type == "TSD II" ~ 1,
    SDM.Type == "GSD" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "2"))


#### full Model 4: Taxonomic ####

tax <- reshape2::melt(as.matrix(model_conservative3.4$Sol))
tax <- tax %>% 
  spread(Var2, value)

# Get values for B0,B1, B2, B3 ,B4, B5, B6, B7
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()
B6 <- c()
B7 <- c()
B8 <- c()
B9 <- c()
B10 <- c()
B11 <- c()

for(i in 1:nrow(tax)){
  b0 <- tax[i,2]
  b1 <- tax[i,3]
  b2 <- tax[i,4]
  b3 <- tax[i,5]
  b4 <- tax[i,6]
  b5 <- tax[i,7]
  b6 <- tax[i,8]
  b7 <- tax[i,9]
  b8 <- tax[i,10]
  b9 <- tax[i,11]
  b10 <- tax[i,12]
  b11 <- tax[i,13]
  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)
  B6 <- c(B6, b6)
  B7 <- c(B7, b7)
  B8 <- c(B8, b8)
  B9 <- c(B9, b9)
  B10 <- c(B10, b10)
  B11 <- c(B11, b11)
}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = Lifstage (juvenile); x2 = 0 (adult), x2 = 1 (juvenile)
# x3 = model4(crocs); x3 = 0 (GSD), x3 = 1 (crocs)
# x4 = model4(lizards); x4 = 0 (GSD), x4 = 1 (lizards)
# x5 = model4(turtles); x5 = 0 (GSD), x5 = 1 (turtles)

# Intercept: AdGSD
AdGSD <- c()

x1 <- 0; x2 <- 0;
x3 <- 0; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  AdGSD <- c(AdGSD, y)
}

# BrGSD
BrGSD <- c()

x1 <- 1; x2 <- 0;
x3 <- 0; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  BrGSD <- c(BrGSD, y)
}

#JuGSD
JuGSD <- c()

x1 <- 0; x2 <- 1;
x3 <- 0; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  JuGSD <- c(JuGSD, y)
}

# ACrocs
AdCrocs <- c()

x1 <- 0; x2 <- 0;
x3 <- 1; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  AdCrocs <- c(AdCrocs, y)
}

# Aliz
Adliz <- c()

x1 <- 0; x2 <- 0;
x3 <- 0; x4 <- 1; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Adliz <- c(Adliz, y)
}

# Atur
Adtur <- c()

x1 <- 0; x2 <- 0;
x3 <- 0; x4 <- 0; x5 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Adtur <- c(Adtur, y)
}

# BCrocs
BrCrocs <- c()

x1 <- 1; x2 <- 0;
x3 <- 1; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  
  BrCrocs <- c(BrCrocs, y)
}

#JuCrocs
JuCrocs <- c()

x1 <- 0; x2 <- 1;
x3 <- 1; x4 <- 0; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  
  JuCrocs <- c(JuCrocs, y)
}

# Brliz
Brliz <- c()

x1 <- 1; x2 <- 0;
x3 <- 0; x4 <- 1; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Brliz <- c(Brliz, y)
}

#Juliz
Juliz <- c()

x1 <- 0; x2 <- 1;
x3 <- 0; x4 <- 1; x5 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Juliz <- c(Juliz, y)
}

# Btur
Brtur <- c()

x1 <- 1; x2 <- 0;
x3 <- 0; x4 <- 0; x5 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Brtur <- c(Brtur, y)
}

#Jutur
Jutur <- c()

x1 <- 0; x2 <- 1;
x3 <- 0; x4 <- 0; x5 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 + B5[i]*x5 +
    B6[i]*x1*x3 + B7[i]*x2*x3 + B8[i]*x1*x4 + B9[i]*x2*x4 + B10[i]*x1*x5 + B11[i]*x2*x5
  Jutur <- c(Jutur, y)
}

# Putting it all together
# put together into single data frame
tax_model <- cbind(AdGSD, BrGSD, JuGSD, AdCrocs, Adliz, Adtur, BrCrocs, JuCrocs, Brliz, Juliz, Brtur, Jutur)
tax_model <- as.data.frame(tax_model)

# modify and add relevant columns
tax_model <- tax_model %>% 
  # get inv logit values to put them on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio",
    substr(Group, 1,1) == "J" ~ "juvenile.sex.ratio")) %>% 
  # add taxonomic group
  mutate(model4 = case_when( # create new column for model4 column
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 4, nchar(Group))  == "Crocs" ~ "crocs",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "liz" ~ "lizards",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "tur" ~ "turtles"))

tax_sum <- tax_model %>% 
  group_by(Group, Life.Stage, model4) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    model4 == "GSD" ~ 0,
    model4 == "crocs" ~ 1,
    model4 == "lizards" ~ 2,
    model4 == "turtles" ~ 3
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.1*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0.1*(grp == "3"))



#### full noliz Model 4: Taxonomic ####

tax <- reshape2::melt(as.matrix(model_sda3.4$Sol))
tax <- tax %>% 
  spread(Var2, value)

# Get values for B0,B1, B2, B3 ,B4, B5, B6, B7
B0 <- c()
B1 <- c()
B2 <- c()
B3 <- c()
B4 <- c()
B5 <- c()
B6 <- c()
B7 <- c()
B8 <- c()

for(i in 1:nrow(tax)){
  b0 <- tax[i,2]
  b1 <- tax[i,3]
  b2 <- tax[i,4]
  b3 <- tax[i,5]
  b4 <- tax[i,6]
  b5 <- tax[i,7]
  b6 <- tax[i,8]
  b7 <- tax[i,9]
  b8 <- tax[i,10]

  
  B0 <- c(B0, b0)
  B1 <- c(B1, b1)
  B2 <- c(B2, b2)
  B3 <- c(B3, b3)
  B4 <- c(B4, b4)
  B5 <- c(B5, b5)
  B6 <- c(B6, b6)
  B7 <- c(B7, b7)
  B8 <- c(B8, b8)

}

# x1 = Lifstage (birth); x1 = 0 (adult), x1 = 1 (birth)
# x2 = Lifstage (juvenile); x2 = 0 (adult), x2 = 1 (juvenile)
# x3 = model4(crocs); x3 = 0 (GSD), x3 = 1 (crocs)
# x4 = model4(turtles); x4 = 0 (GSD), x4 = 1 (turtles)

# Intercept: AdGSD
AdGSD <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  AdGSD <- c(AdGSD, y)
}

# BrGSD
BrGSD <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  BrGSD <- c(BrGSD, y)
}

#JuGSD
JuGSD <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  JuGSD <- c(JuGSD, y)
}

# ACrocs
AdCrocs <- c()

x1 <- 0; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  AdCrocs <- c(AdCrocs, y)
}


# Atur
Adtur <- c()

x1 <- 0; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  Adtur <- c(Adtur, y)
}

# BCrocs
BrCrocs <- c()

x1 <- 1; x2 <- 0; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  BrCrocs <- c(BrCrocs, y)
}

#JuCrocs
JuCrocs <- c()

x1 <- 0; x2 <- 1; 
x3 <- 1; x4 <- 0
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  JuCrocs <- c(JuCrocs, y)
}

# Btur
Brtur <- c()

x1 <- 1; x2 <- 0; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  Brtur <- c(Brtur, y)
}

#Jutur
Jutur <- c()

x1 <- 0; x2 <- 1; 
x3 <- 0; x4 <- 1
for(i in 1:nrow(tax)){
  y <- B0[i] + B1[i]*x1 + B2[i]*x2 + B3[i]*x3 + B4[i]*x4 +
    B5[i]*x1*x3 + B6[i]*x2*x3 + B7[i]*x1*x4 + B8[i]*x2*x4 
  Jutur <- c(Jutur, y)
}

# Putting it all together
# put together into single data frame
tax_model <- cbind(AdGSD, BrGSD, JuGSD, AdCrocs, Adtur, BrCrocs, JuCrocs, Brtur, Jutur)
tax_model <- as.data.frame(tax_model)

# modify and add relevant columns
tax_model <- tax_model %>% 
  # get inv logit values to put them on probability space
  mutate_all(inv.logit) %>% 
  # gather all values into single column with other delineating which group they belong to
  gather(Group) %>% 
  # add life stage value
  mutate(Life.Stage = case_when( # create new column for life stage
    substr(Group, 1,1) == "A" ~ "adult.sex.ratio",
    substr(Group, 1,1) == "B" ~ "birth.sex.ratio",
    substr(Group, 1,1) == "J" ~ "juvenile.sex.ratio")) %>% 
  # add taxonomic group
  mutate(model4 = case_when( # create new column for model4 column
    substr(Group, nchar(Group) - 2, nchar(Group))  == "GSD" ~ "GSD",
    substr(Group, nchar(Group) - 4, nchar(Group))  == "Crocs" ~ "crocs",
    substr(Group, nchar(Group) - 2, nchar(Group))  == "tur" ~ "turtles"))

tax_sum <- tax_model %>% 
  group_by(Group, Life.Stage, model4) %>% 
  summarize(point_estimate(value),
            ci(value, method = "HDI"),
            mean_value = mean(value, na.rm = TRUE),
            lower_90 = quantile(value, probs = 0.05, na.rm = TRUE),
            upper_90 = quantile(value, probs = 0.95, na.rm = TRUE),
            lower_95 = quantile(value, probs = 0.025, na.rm = TRUE),
            upper_95 = quantile(value, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(grp = case_when(
    model4 == "GSD" ~ 0,
    model4 == "crocs" ~ 1,
    model4 == "turtles" ~ 3
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0*(grp == "3"))




