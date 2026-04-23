# Terebiznik.et.al.2026
Data analysis for Terebiznik et al. (2026) Ontogenetic Variation in Sex Ratios and the Adaptive Significance of Temperature-Dependent Sex Determination.

Data and analysis of wild, natural sex ratio records of reptiles to investigate the origins of temperature-dependent sex determination. This analysis was conducted in RStudio (Version 2026.01.1+403) using the MCMCglmm package. 

**System Requirements**

The only system requirements is the program R

**Installation**

To install the MCMCglmm package, run the following code in R.

```r
install.packages(MCMCglmm)
library(MCMCglmm)
```
Instructions on how to use the MCMCglmm package can be found here: https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf

**Instructions for Use**

Included in this package is

1. a unique data set of wild, natural sex ratio records of reptiles
2. a dataset of species-level sexual size dimorphism data of reptiles
3. a dataset that includes population-level maturation age data of reptiles, obtained from Bokony et al. (2019)
4. phylogenies of species included in the dataset
5. code to reformat data, run MCMCglmm models and model selection for the sex ratio data set, run MCMCglmm models for bimaturism and sexual size dimorphism analyses, generate data matrix for model visualization, and generate plots
6. Citation Code document that links citations in the sex ratio and the Bokony et al. data set to their references

To run this code, download the data set and phylogenues and save to your working directory. Then run the code in order to first reformat the data into a binary data set, run that data set through each model, identify the best model through model selection. Running each model takes a handful of hours.

**Sex Ratio Dataset Description**

This is a unique data set of 584 sex ratio records representing 166 species of reptiles from 409 populations. The data set was assembled from previously published literature of Bokony et al. (2019), with records from ROSIE, the Reptilian Offspring Sex and Incubation Environment (Krueger and Janzen 2022; ROSIE, 2021; v1.0.0), and personal communications. Inclusion criteria for each sex ratio datum were 1) recorded from a natural, unmanipulated population using 2) non-sex biased capture methods and 3) determined by reliable sexing methods with 4) a recorded sample size and 5) a life stage of sexed individuals.

The data is further subdivided into a conservative data where each one of those criteria are clearly met and 'sda' data which includes sex ratio records where the sexing method was not explicitely stated, but the species is sexually dimorphic.

For each record, the data consists of:

- population ID: an ID value for population of each sex ratio record
- species: scientific latin name of species
- Life.Stage: the life stage of either birth (birth.sex.ratio), juvenile (juvenile.sex.ratio), or adult (adult.sex.ratio)
- Sex.Ratio: the sex ratio presented as proportion male
- N: sample size for sex ratio record
- sex.determination: type of sex determining mechanism of either genetic sex determination (GSD) or temperature-dependent sex determination (TSD)
- SDM.Type: more specific category of sex determining mechanism of either GSD, TSD Ia, or TSD II
- SuperTaxa: whether the species belongs to a sister clade of crocodilians and turtles (crocoturtle) or squamates (squamata)
- Model4: category fit for model 4 which lists GSD species as 'GSD' and any TSD species as their taxon (turtle, croc, lizard)
- Conservative: A Yes/No column indicating whether the record is part of the conservative data set (Y) or is a sexual dimorphism assumed record (N)
- Citation: citation or citation code. Citation codes can be matched to citations in the 'Citations Codes.csv' file
- Source: whether the source of the data came from Bokony et al 2018 (Bokony), from the ROSIE data base (ROSIE) or from personal communications (Other)
- sex.det.source: reference for sex determining mechanism for the species

**Phylogeny Description**

This phylogeny was created based on the phylogeny used in Bókony et al. (2019). Any species in this sex ratio data set not in the Bókony et al. (2019) were included using phytools package in R with the placement of each species was based on the same phylogenies Bókony et al. (2019) used to construct their phylogeny (Oaks 2011, Guillon et al. 2012, Pyron et al. 2013)

The phylogeny for species of the conservative data only is called 'conservative full phylogeny.phy' The phylogeny for only hatchling and adult conservative data is called 'conservative BA phylogeny.phy' The phylogeny for the sda data is called 'ultimate phylogeny.phy'. The species list was the same for the sda data with and without juveniles, therefore only one phylogeny was needed for all sda models.

**TSD Analysis Code Description**

The code for running models is attached in the file 'Terebiznik et al Code.R' This code includes the following steps:

1. Loading in of necessary packages and data sets
2. Loading in and formatting phylogenies for use in MCMCglmm models
3. Creating binary data sets for use in MCMCglmm models
4. Running models testing our 4 hypothesis using only hatchling and adult data
5. Model comparisons to identify the best model
6. Running the best model with juvenile data included

**Sexual Size Dimorphism Dataset Description**

This is a unique data set of 802 species-level maximum body length (snout-vent-length (SVL)/ carapace length (CL)) records for each sex. Exception includes 3 species for which the maximum body length data was inavailable and thus alternative body length parameters were used. Sexual dimorphism data were obtained from Oskyrko et al. (2024),  Rhodin et al. (2021),  Ceballos et al. (2013), Ceballos & Iverson (2014), Halámková et al. (2013), Deeming (2022) and further literature survey. 

For each record, the data consists of:
- species: scientific latin name of species
- F.cm.maxSVL/CL: female maximum body length in cm
- M.cm.maxSVL/CL: male maximum body length in cm
- dataset: citation
- Notes for exception: notes that indicate the nature of the body length data provided if it is not maximum SVL or CL.

**Bokony et al. (2019) Dataset Description**

This is a dataset provided by the publication below:
Bókony V, Milne G, Pipoly I, Székely T, Liker A. Sex ratios and bimaturism differ between temperature-dependent and genetic sex-determination systems in reptiles. BMC Evol Biol. 2019 Feb 18;19(1):57. doi: 10.1186/s12862-019-1386-3. 

Please consult the original paper for details.

**Citation Code Description**

In the sex ratio data set, records from the Bokony et al (2019) data have citation codes rather than direct references in the dataset. This spreadsheet matches each citation code to its corresponding reference
