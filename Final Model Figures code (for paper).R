# File to create figures for paper
# difference: sources data from Model to Data Matrix file

#### Set up ####
# Load in packages
library(readr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggpubr)
library(MCMCglmm)
library(emmeans)
library(see)
library(qdapRegex)
library(MuMIn)
library(bayestestR)
library(cowplot)


# Load in models
setwd("/Users/annettew/Documents/Grad School/Mariel's paper") # SMH
source("Models to plotting data.R")

# Load in Data ####
# GitHub working directory
HH <- "https://raw.githubusercontent.com/MarielTerebiznik/TSD-Analysis/refs/heads/main/"

### Load in sex ratio data
data <- read.csv(paste0(HH, "sex%20ratio%20data%20final.csv"))

data <- df %>% 
  mutate(Sample.Size = case_when(
    N <= 100 ~ 2,
    N <= 250 ~ 3,
    N <= 500 ~ 4,
    N <= 1000 ~ 5,
    N <= 2000 ~ 6,
    N >= 2001 ~ 7))

sda <- data
conservative <- filter(data, Conservative == "Y")

#### Model 0: SSD ####
library(ggplot2)

# Function to compute predictions for a single row of newdata
str(hat_binary$model4)
hat_binary$model4 <- as.factor(hat_binary$model4)
levels(hat_binary$model4) # Verify levels exist

# Generate predict dataframe 
n_points <- 100  # Number of SDI points
#n_levels <- levels(hat_binary$model4) 
n_levels <- as.factor(c("GSD", "crocs", "turtles"))
n_lvl <- length(n_levels)

new_data <- data.frame(model4 = rep(n_levels, each = n_points))

for (model4 in n_levels){
  new_data$SDI[new_data$model4 == model4] <- rep(seq(min(hat_binary$SDI[hat_binary$model4 == model4]), 
                                                                       max(hat_binary$SDI[hat_binary$model4 == model4]),
                                                                       length.out = n_points))
}

# Add any required placeholder columns
new_data$phylo <- NA
new_data$species <- NA
new_data$population.ID <- NA

new_data$pred_mean <- NA
new_data$pred_lower <- NA
new_data$pred_upper <- NA

# Compute posterior samples
post_samples <- as.matrix(model_sdi$Sol)
colnames(post_samples)

compute_prediction <- function(row, sample) {
  if (row$model4 == "GSD") {
    pred <- sample[, "(Intercept)"] + 
      sample[, "SDI"] * row$SDI 
  } else if (row$model4 == "crocs") {
    pred <- (sample[, "(Intercept)"] + sample[, "model4crocs"]) + 
      (sample[, "SDI"] + sample[, "model4crocs:SDI"]) * row$SDI 
  } else {
    pred <- (sample[, "(Intercept)"] + sample[, "model4turtles"]) + 
      (sample[, "SDI"] + sample[, "model4turtles:SDI"]) * row$SDI
  }
  return(pred)  # Convert to probability if using logit link
}

for (i in 1:nrow(new_data)) {
  probs <- compute_prediction(new_data[i, ], post_samples)
  new_data$pred_mean[i] <- mean(probs)
  new_data$pred_lower[i] <- quantile(probs, 0.025)  # 2.5% percentile
  new_data$pred_upper[i] <- quantile(probs, 0.975)  # 97.5% percentile
}

#mat_binary$Sex_numeric <- as.numeric(mat_binary$Sex) - 1  # Male=1, Female=0
ggplot() +
  # Observed data (jittered points for visibility)
  geom_jitter(
    data = hat, 
    aes(x = SDI, y = Sex.Ratio, color = model4, shape = model4),
  ) +
  # Model predictions (mean ± 95% CI)
  geom_ribbon(
    data = new_data,
    aes(x = SDI, ymin = pred_lower, ymax = pred_upper, fill = model4),
    alpha = 0.3
  ) +
  geom_line(
    data = new_data,
    aes(x = SDI, y = pred_mean, color = model4),
    linewidth = 1.2
  ) +
  labs(
    title = "SDI vs. hatchling sex ratio for species with different SDM",
    x = "SDI", 
    y = "Hatchling Sex Ratio (proportion of male)",
    color = "Taxonomic Group",
    fill = "Taxonomic Group",
    shape = "Taxonomic Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("GSD" = "blue3", "crocs" = "coral", "turtles" = "indianred")) +
  scale_fill_manual(values = c("GSD" = "skyblue", "crocs" = "orange", "turtles" = "pink")) +
  scale_shape_manual(values = c("GSD" = 15, "crocs" = 17, "turtles" = 19))

#### Model 1: BA SMH ####
plotdata <- data %>% 
  mutate(grp = case_when(
    sex.determination == "TSD" ~ 1,
    sex.determination == "GSD" ~ 2
  )) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~ 2
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 + 0.225*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + -.225*(grp == "2"))

plotsda <- plotdata
plotcon <- filter(plotdata, Conservative == "Y")

level_order <- c("birth.sex.ratio", "adult.sex.ratio")

M1 <- ggplot(data = smh_sum, aes(colour = sex.determination, shape = sex.determination)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #raw data
  geom_point(data = plotcon, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.2)) + # actual data points 
  # summary data
  geom_point(data = smh_sum, aes(x= Life.Stage2, y = mean_value, color = sex.determination), size = 5) +
  geom_errorbar(data = smh_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  # aestheics
  #scale_size_binned("Sample Size") + # specifiying legend
  scale_colour_manual(values = c("blue3", "orange2"), labels = c("GSD", "TSD")) +
  scale_shape_manual( values = c(16, 18), labels = c("GSD", "TSD")) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Single Mechanism Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','Adult')) +
  ggtitle("Single Mechanism Model") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M1

#### Model 2: BA supertaxa #####
plotdata <- data %>% 
  mutate(grp = case_when(
    EvoH == "GSD" ~ 0,
    EvoH == "TSD crocoturtle" ~ 1,
    EvoH == "TSD squamata" ~ 2)) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "adult.sex.ratio" ~ 3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) 

plotsda <- plotdata
plotcon <- filter(plotdata, Conservative == "Y")

M2 <- ggplot(data = sut_sum, aes(colour = EvoH, shape = EvoH)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #Posterior distributions
  geom_point(data = plotcon, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  #raw data    # summary data
  geom_point(data = sut_sum, aes(x= Life.Stage2, y = mean_value, colour = EvoH), size = 5) +
  geom_errorbar(data = sut_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, size = 1) +
  # aestheics
  scale_colour_manual(values = c("blue3", "red","darkred"), labels = c("GSD", "TSD crocoturtle", "TSD squamata")) +
  scale_shape_manual(
    values = c("GSD" = 16, "TSD crocoturtle" = 18, "TSD squamata" = 15)) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c("Hatchling","", "Adult")) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Evolutionary History Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  ggtitle("TSD Type Model") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M2
#### Model 3: BA TSD #####
plotdata <- data %>% 
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

plotsda <- plotdata
plotcon <- filter(plotdata, Conservative == "Y")

level_order <- c("birth.sex.ratio","", "adult.sex.ratio")

M3 <- ggplot(data = tsd_sum, aes(colour = SDM.Type, shape = SDM.Type)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #Posterior distributions
  geom_point(data = plotcon, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  #geom_violinhalf(data = sdm_model,aes(x = factor(Life.Stage, level = level_order), y = value, color = SDM.type), size = 1.5, fill = NA)+
  #raw data    # summary data
  geom_point(data = tsd_sum, aes(x= Life.Stage2, y = mean_value, color = SDM.Type), size = 5) +
  geom_errorbar(data = tsd_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  #geom_line(data = sdm_sum, aes(x = Life.Stage2, y = avg, group = `...7`), size = 1, alpha = 0.75) +
  # aestheics
  scale_shape_manual(values = c(16, 15, 17), labels = c("GSD", "TSD Ia", "TSD II")) +
  scale_color_manual(values = c("blue3", "coral2", "yellow3"), labels = c("GSD", "TSD Ia","TSD II")) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','','Adult')) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "TSD Type Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  ggtitle("TSD Type Model") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M3
#### Model 4: BA Tax ####
plotdata <- data %>% 
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

plotsda <- plotdata
plotcon <- filter(plotdata, Conservative == "Y")

level_order <- c("birth.sex.ratio", "","adult.sex.ratio")
label2 <- c("GSD", "turtles", "crocs","lizards")

M4 <- ggplot(data = tax_sum, aes(colour = factor(model4, level=label2), shape = factor(model4, level=label2))) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", size = 1, alpha = 0.3) + # line demarking 50/50
  #Posterior distributions
  #raw data
  geom_point(data = plotcon, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  # summary data
  geom_point(data = tax_sum, aes(x= Life.Stage2, y = mean_value, color = model4), size = 5) +
  geom_errorbar(data = tax_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  # aestheics
  #facet_wrap(~sex.determination,nrow = 2) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  scale_shape_manual(values = c(16, 17, 18, 15), labels = label2) +
  scale_colour_manual(values = c("blue3", "violet", "#3CBB57FF", "purple3"), labels = label2)+
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','','Adult')) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Taxonomic Model",
    color = "Sex Determining\n Taxonomic Group",
    shape = "Sex Determining\n Taxonomic Group"
  ) +
  ggtitle("Taxonomic Model") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), 
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10))
M4



#### Model 1: full SMH ####
plotdata <- conservative %>% 
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

plotfull <- plotdata
plotnoliz <- filter(plotdata, model4 != "lizards")

level_order <- c("birth.sex.ratio", "juvenile.sex.ratio","adult.sex.ratio")

M1 <- ggplot(data = smh_sum, aes(colour = sex.determination, shape = sex.determination)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #raw data
  geom_point(data = plotfull, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.2)) + # actual data points 
  # summary data
  geom_point(data = smh_sum, aes(x= Life.Stage2, y = mean_value, color = sex.determination), size = 5) +
  geom_errorbar(data = smh_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  # aestheics
  #scale_size_binned("Sample Size") + # specifiying legend
  scale_colour_manual(values = c("blue3", "orange2"), labels = c("GSD", "TSD")) +
  scale_shape_manual( values = c(16, 18), labels = c("GSD", "TSD")) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Single Mechanism Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','Juvenile','Adult')) +
  ggtitle("Single Mechanism Model") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M1
#### Model 2: full supertaxa #####
plotdata <- conservative %>% 
  mutate(grp = case_when(
    EvoH == "GSD" ~ 0,
    EvoH == "TSD crocoturtle" ~ 1,
    EvoH == "TSD squamata" ~ 2)) %>% 
  mutate(Life.Stage2 = case_when(
    Life.Stage == "birth.sex.ratio" ~ 1,
    Life.Stage == "juvenile.sex.ratio" ~ 2,
    Life.Stage == "adult.sex.ratio" ~ 3
  )) %>% 
  mutate(Life.Stage2 = as.numeric(as.character(Life.Stage2))) %>%
  mutate(Life.Stage2 = Life.Stage2 - 0.3*(grp == "0")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "2")) 


M2 <- ggplot(data = sut_sum, aes(colour = EvoH, shape = EvoH)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #Posterior distributions
  geom_point(data = plotdata, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  #raw data    # summary data
  geom_point(data = sut_sum[1:8,], aes(x= Life.Stage2, y = mean_value, colour = EvoH), size = 5) +
  geom_errorbar(data = sut_sum[1:8,], aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, size = 1) +
  # aestheics
  scale_colour_manual(values = c("blue3", "red","darkred"), labels = c("GSD", "TSD crocoturtle", "TSD squamata")) +
  scale_shape_manual(
    values = c("GSD" = 16, "TSD crocoturtle" = 18, "TSD squamata" = 15)) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c("Hatchling","Juvenile", "Adult")) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Evolutionary History Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M2
#### Model 3: full TSD #####
plotdata <- conservative %>% 
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
  mutate(Life.Stage2 = Life.Stage2 + 0.26*(grp == "1")) %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0.26*(grp == "2"))

plotfull <- plotdata
plotnoliz <- filter(plotdata, model4 != "lizards")

level_order <- c("birth.sex.ratio", "juvenile.sex.ratio","adult.sex.ratio")

M3 <- ggplot(data = tsd_sum, aes(colour = SDM.Type, shape = SDM.Type)) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", linewidth = 1, alpha = 0.3) + # line demarking 50/50
  #Posterior distributions
  geom_point(data = plotfull, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  #geom_violinhalf(data = sdm_model,aes(x = factor(Life.Stage, level = level_order), y = value, color = SDM.type), size = 1.5, fill = NA)+
  #raw data    # summary data
  geom_point(data = tsd_sum, aes(x= Life.Stage2, y = mean_value, color = SDM.Type), size = 5) +
  geom_errorbar(data = tsd_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  #geom_line(data = sdm_sum, aes(x = Life.Stage2, y = avg, group = `...7`), size = 1, alpha = 0.75) +
  # aestheics
  scale_shape_manual(values = c(16, 15, 17), labels = c("GSD", "TSD Ia", "TSD II")) +
  scale_color_manual(values = c("blue3", "coral2", "yellow3"), labels = c("GSD", "TSD Ia","TSD II")) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','Juvenile','Adult')) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "TSD Type Model",
    color = "Sex Determining\n Mechanism",
    shape = "Sex Determining\n Mechanism"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8))
M3
#### Model 4: full Tax ####
plotdata <- conservative %>% 
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


level_order <- c("birth.sex.ratio", "juvenile.sex.ratio","adult.sex.ratio")
label2 <- c("GSD", "turtles", "crocs","lizards")

M4 <- ggplot(data = tax_sum, aes(colour = factor(model4, level=label2), shape = factor(model4, level=label2))) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", size = 1, alpha = 0.3) + # line demarking 50/50
  #raw data
  geom_point(data = plotdata, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  # summary data
  geom_point(data = filter(tax_sum, Group != "Juliz"), aes(x= Life.Stage2, y = mean_value, color = model4), size = 5) +
  geom_errorbar(data = filter(tax_sum, Group != "Juliz"), aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  # aestheics
  #facet_wrap(~sex.determination,nrow = 2) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  scale_shape_manual(values = c(16, 17, 18, 15), labels = label2) +
  scale_colour_manual(values = c("blue3", "violet", "#3CBB57FF", "purple3"), labels = label2)+
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','Juvenile','Adult')) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Taxonomic Model",
    color = "Sex Determining\n Taxonomic Group",
    shape = "Sex Determining\n Taxonomic Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), 
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10))
M4


#### Model 4: full noliz Tax ####
plotdata <- data %>% 
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
  mutate(Life.Stage2 = Life.Stage2 + 0.3*(grp == "1"))  %>% 
  mutate(Life.Stage2 = Life.Stage2 - 0*(grp == "3"))

plotnoliz <- filter(plotdata, model4 != "lizards")

level_order <- c("birth.sex.ratio", "juvenile.sex.ratio","adult.sex.ratio")
label2 <- c("GSD", "turtles", "crocs")

M4 <- ggplot(data = tax_sum, aes(colour = factor(model4, level=label2), shape = factor(model4, level=label2))) +
  # Adding the horizontal 50/50 line
  geom_hline(yintercept = 0.5, colour = "darkgrey", size = 1, alpha = 0.3) + # line demarking 50/50
  #raw data
  geom_point(data = plotnoliz, aes(x = Life.Stage2, y = Sex.Ratio, size = Sample.Size), alpha = 0.1, position = position_jitter(width = 0.08)) + # actual data points 
  # summary data
  geom_point(data = tax_sum, aes(x= Life.Stage2, y = mean_value, color = model4), size = 5) +
  geom_errorbar(data = tax_sum, aes(x = Life.Stage2, ymin = (lower_95), ymax = (upper_95)), width = 0.05, position= position_dodge(0.005), size = 1) +
  # aestheics
  #facet_wrap(~sex.determination,nrow = 2) +
  scale_size(labels = c("0-100", "101-250", "251-500", "500-1000","1000-2000","2000+")) +
  scale_shape_manual(values = c(16, 17, 18), labels = label2) +
  scale_colour_manual(values = c("blue3", "violet", "#3CBB57FF"), labels = label2)+
  theme_bw() +
  scale_x_discrete(limits = level_order, labels = c('Hatchling','Juvenile','Adult')) +
  labs(
    x = "Life Stage", 
    y = "Sex Ratio (Probability of being Male)",
    title = "Taxonomic Model",
    color = "Sex Determining\n Taxonomic Group",
    shape = "Sex Determining\n Taxonomic Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), 
    axis.text = element_text(size = 10), 
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12), 
    legend.text = element_text(size = 10))
M4



#### Combine Plots ####
library(patchwork)
M1 <- M1+labs(title = "A. Single Mechanism Model") & theme_classic(base_size = 12) + theme(
  legend.position = "none", 
  axis.title = element_text(size = 0, face = "bold"), 
  plot.title = element_text(size = 12, face = "bold"))
M2 <- M2+labs(title = "B. Evolutionary History Model") & theme_classic(base_size = 12) + theme(
  legend.position = "none", 
  axis.title = element_text(size = 0, face = "bold"), 
  plot.title = element_text(size = 12, face = "bold"))
M3 <- M3+labs(title = "C. TSD Type Model") & theme_classic(base_size = 12) + theme(
  legend.position = "none", 
  axis.title = element_text(size = 0, face = "bold"), 
  plot.title = element_text(size = 12, face = "bold"))
M4 <- M4 +
  labs(title = "D. Taxonomic Model") & theme_classic(base_size = 12) + theme(
  legend.position = "none", 
  axis.title = element_text(size = 0, face = "bold"), 
  plot.title = element_text(size = 12, face = "bold"))

# Restyle models
combined_plot1 <- M1 + M2 + M3 + M4 +  # Vertical stack
  plot_layout(ncol=2,
    guides = "collect", # Merge legends
    widths = 1,  # Adjust spacing
    heights = c(1,1)) 

combined_plot1

library(ggplot2)

# Create text grobs
sex_ratio_label <- ggplot() + 
  annotate("text", x = -5, y = 0, 
           label = "Sex Ratio (Proportion Male)", 
           angle = 90, size = 4) +
  theme_void()

life_stage_label <- ggplot() + 
  annotate("text", x = 1, y = 1, 
           label = "Life Stage", size = 4) +
  theme_void()

# Insert into layout
add_x <- combined_plot1/life_stage_label +
  plot_layout(heights = c(1,0.1))

add_y <- sex_ratio_label + add_x +
  plot_layout(widths = c(0.1,1))

add_y & theme_classic(base_size = 12) 

add_y

# Add legends
library(grid)
library(gtable)

# Create standalone legend ####
order <- c("0-100", "101-250", "251-500", "500-1000","1001-2000","2000+")
order1 <- c("GSD","TSD Crocoturtles","TSD Squamates","TSD Turtles","TSD Crocodilians","TSD Lizards")
#order2 <- c("GSD","TSD Turtles","TSD Crocodilians")
legend_SS <- cowplot::get_legend(
  ggplot(data.frame(x = factor(order, levels = order),y=1:2))+
    geom_point(aes(x, y, size = x), alpha = 0.2) +
    scale_size_manual(values = c(1,2,3,4,5,6)) +
    labs(
      size = "Sample size",
    ) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
  & theme_classic(base_size = 12)
)

legend_SDM <- cowplot::get_legend(
  ggplot(data.frame(x = c("GSD", "TSD (all)","TSD Ia","TSD II"), y = 1:2)) +
    geom_point(aes(x, y, color = x, shape = x), size=5) +
    scale_color_manual(values = c("blue3", "orange2","coral2", "yellow3")) +
    scale_shape_manual(values = c(16, 18, 15, 17)) +
    labs(
      color = "Sex Determining\n Mechanism",
      shape = "Sex Determining\n Mechanism"
    ) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
  & theme_classic(base_size = 12)
)

legend_sut <- cowplot::get_legend(
  ggplot(data.frame(x = factor(order1, levels = order1), y = 1:6)) +
    geom_point(aes(x, y, color = x, shape = x), size=5) +
    scale_color_manual(values = c("blue3", "red","darkred","violet", "#3CBB57FF", "purple3")) +
    scale_shape_manual(values = c(16, 18, 15, 17, 18, 15)) +
    labs(
      color = "Sex Determining\n Taxonomic Group",
      shape = "Sex Determining\n Taxonomic Group"
    ) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
  & theme_classic(base_size = 12)
)

legend_tax <- cowplot::get_legend(
  ggplot(data.frame(x = factor(order2, levels = order2), y = 1:3)) +
    geom_point(aes(x, y, color = x, shape = x), size=5) +
    scale_color_manual(values = c("blue3","violet", "#3CBB57FF")) +
    scale_shape_manual(values = c(16, 17, 18)) +
    labs(
      color = "Sex Determining\n Taxonomic Group",
      shape = "Sex Determining\n Taxonomic Group"
    ) +
    theme(legend.box.margin = margin(0, 0, 0, 12))
  & theme_classic(base_size = 12)
)

# Stack legends
stacked_legends <- plot_grid(
  legend_SS,
  legend_SDM,
  legend_sut,
#  legend_tax,
  ncol = 1,          # Vertical stack
  align = "v",       # Align vertically
  rel_heights = c(1,1, 1, 1)  # Adjust spacing
)

# Display
grid::grid.draw(stacked_legends)

# Insert into plot
final_plot <- add_y + stacked_legends + 
  plot_layout(nrow=1, widths = c(0.3,3,1))

final_plot & 
  theme_classic(base_size = 12)

final_plot
#### Comparing models 1-4 #####
output <- model.sel(model1.1, model1.2, model1.3, model1.4, rank = "DIC")

