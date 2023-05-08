# BIO-539-Final-Project
Repository for reproducible analysis of BIO 539 Project
The purpose of this repo is to store and centralize all material and outputs (R markdowns, graphs, raw scripts, analysis papers) for the Final Project in BIO 539

## Objectives

1. Understand association between tick-borne pathogens in terms of singular infection and co-infections in Rhode Island.
2. Understand how enviromental covariates influence infection rates.
3. Identify a possible spatial relationship for the distribution of pathogens across Rhode Island.
4. identify pathogen prevalence for tick-borne pathogens in RI.

## Materials

### Raw Data
The raw data is given in the form of various csv's which include:

Raw Data - csv detailing infection for ~ 1500 nymphs across 6 pathogens

Sampling Site Longitude and Latitude - csv of spatial informations for where we collected nymphs for this project

Specimen Date Temp Humidity - CSV of meta data which includes time of collection, data of collection, collector initals, and humidity and temperature at beggining of sampling period

Sampling Site Enviromental Parameter - CSV of enviromental covariates including # of Tick Adverse Humdity events experienced for each site in june and % forest cover for each site across mixed, decidious, or conifer stands

Please make sure when running the code that these csv's are saved in your working directory. the file paths are not hard coded, just having it in your working directory will be fine.

### Outputs

Pathogen Data Script (Model and organization) - This R script contains the entirety of all analysis done for this project including all code for modelling and graphical outputs

Graphics - This includes all outputs in the form of tables, graphs, and other visual elements. It is not coupled in the code. However, comments within the pathogen data Script that should specify which code makes which graphs

R Markdown - contains r markdown file with all outputs and thoughts/comments during analysis. The r markdown contains commentary but you will ned need to knit it yourself. The large spatial models performed near the end are computationally intensive so remember to devote plenty of resources towards it.

## Tips on use

Please remember that the raw data csv only contains information on infection of individual nymphs. If you want to join this with spatial, site data then you need to join based on the nymph ID column. That nymph ID can be used to join all of these csvs together. Instances of this are documented within the Pathogen Data Script.

In the section outlined "Spatial Data Analysis" i use a function within the package "gcKrig" called mlegc. This runs simulations to create spatial association models for covariates across multiple variables. It is resource heavy on processors. It was run on my machine using 12th Gen Intel(R) Core(TM) i7-1250U 1.10 GHz 12 core processor. Each model took ~5 minuts to run but this will depend your system's specifications.

Information on what packages and libraries can be in the pathogen data script and R markdown.
