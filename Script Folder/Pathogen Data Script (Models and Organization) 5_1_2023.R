#Install Packages
install.packages("Rtools")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages('stringr')
install.packages("installr")
install.packages("dplyr")
install.packages("gtable")
install.packages("cooccur")
install.packages("snowfall")
install.packages("parallel")
install.packages("colorspace")
install.packages("latticeExtra")
install.packages("ggfortify")
install.packages("rnaturalearth")
install.packages("tidygeocoder")
install.packages("gcKrig")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("sjlabelled")
#Libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(installr)
library(gtable)
library(purrr)
library(lubridate)
library(cooccur)
library(akima)
library(geoR) 
library(maps)
library(maptools)
library(lattice)
library(sp)
library(gstat)
library(stats)
library(corrplot)
library(snowfall)
library(gcKrig)
library(parallel)
library(colorspace)
library(latticeExtra)
library(mapdata)
library(ggfortify)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(tidygeocoder)
library(maps)
library(usmap)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

#Read Tables and import
pathogen_data_nymphs_2022 <- read_csv("././././././pathogen_data_nymphs_2022_Enhanced.csv", 
                                      locale=locale(encoding="latin1")) #this table just has the specimen and pathogen informations

Sampling_Site_Longitude_and_Lattitude <- read_csv("Sampling Site Longitude and Lattitude.csv") #longitude and lattitude

Specimen_Date_Temp_Humidity <- read_csv("Specimen_Date_Temp_Humidity.csv") #temp and humidity data for each collection event

Sampling_Site_Enviromental_Parameter <- read_csv("Sampling Site Enviromental Parameter.csv") # TAMES and % Forest Cover

#Modify Specimen environmental and collection data for ease of use
Specimen_Date_Temp_Humidity_corrected <- Specimen_Date_Temp_Humidity %>%
  separate(Date, sep='_', into =c("month","day","year")) %>%
  rename(Site_Code = `SITE CODE`) %>%
  rename(Vial_code = `VIAL CODE`)

#convert vial code column to character to allow for join later (Join wont happen unless columns are same class)
Specimen_Date_Temp_Humidity_corrected$Vial_code <- as.character(
  Specimen_Date_Temp_Humidity_corrected$Vial_code)

#Change All Instances Of Positive and Negative To 1 and 0 for analysis and to allow for summing of values
pathogen_data_nymphs_2022[pathogen_data_nymphs_2022=="Positive"] <- "1" 

pathogen_data_nymphs_2022[pathogen_data_nymphs_2022=="Negative"] <- "0"

#Remove all instances of non-ixodes and failed DNA extraction
#additionally, only include pathogens that are present

Pathogen_data_Nofails_OnlyIxodes <- pathogen_data_nymphs_2022 %>%
  filter(`Tick DNA Internal quality`=='Pass') %>%
  filter(`Dermacentor variabilis`==0) %>%
  filter(`Amblyomma americanum`==0) %>%
  select(c(`URI_Nymph_Number_on_Tube`,
           `URI_Vial_Number_on_tube`,
           'Ixodes scapularis',
           `Borrelia Genus`,
           `Borrelia burgdorferi sensu lato`,
           `Borrelia miyamotoi `,
           `Borrelia mayonii `,
           `Babesia microti`,
           `Ehrlichia muris`,
           `Anaplasma phagocytophilum`,
           `Deer Tick Virus`))

#Convert Pathogen Columns to Numeric for calcualation

Pathogen_data_numeric <- Pathogen_data_Nofails_OnlyIxodes %>%
  mutate_at(c('Borrelia Genus',
              'Borrelia burgdorferi sensu lato',
              'Borrelia miyamotoi ',
              'Borrelia mayonii ',
              'Babesia microti',
              'Ehrlichia muris',
              'Anaplasma phagocytophilum',
              'Deer Tick Virus',
              'Ixodes scapularis'), as.numeric)

#This creates a column in pathogen_data_numeric summarizing co-infections based on presence and infection amount

Pathogen_data_numeric$infection_amount <- rowSums(Pathogen_data_numeric[,5:11])

Pathogen_data_numeric <- Pathogen_data_numeric %>%
  mutate(Co_infection_present = if_else(infection_amount > 1,1,0))

#convert Vial and Specimen Number to character
Pathogen_data_numeric_pain <- Pathogen_data_numeric %>%
  mutate_at(c('URI_Nymph_Number_on_Tube', 'URI_Vial_Number_on_tube'), as.character) %>%
  mutate(Infection_present =
          if_else(`Borrelia Genus` | `Borrelia burgdorferi sensu lato` 
                  | `Borrelia miyamotoi ` | `Borrelia mayonii ` | `Babesia microti` 
                  | `Ehrlichia muris` | `Anaplasma phagocytophilum` | `Deer Tick Virus` == 1,1,0)) %>%
relocate(Infection_present, .after = `Ixodes scapularis`) 

Pathogen_data_numeric_pain <- Pathogen_data_numeric_pain[Pathogen_data_numeric_pain$`Ixodes scapularis` !=0, ]

#All of the above code will give us information on on individual nymph pathogen information

#_______________________________________________________________________________

##All of the code below this will help us summarise information at the site and time level

#Summarize total infected by vial number

Pathogen_data_total_infected_by_site <- Pathogen_data_numeric_pain %>%
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_infected_Ixodes = sum(`Infection_present`))
  

###Summarize by total positive specimens
  
#Borrelia burgdorferi
Pathogen_data_Total_Borrelia_burgdorferi <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Borrelia_burgdorferi = sum(`Borrelia burgdorferi sensu lato`))
 
#Ixodes Scapularis
Pathogen_data_Total_Ixodes_scapularis <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Ixodes_Scapularis = sum(`Ixodes scapularis`))

#Borrelia genus
Pathogen_data_Total_Borrelia_genus <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Borrelia_genus= sum(`Borrelia Genus`))

#Borrelia miyamotoi
Pathogen_data_Total_Borrelia_miyamotoi <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Borrelia_miyamotoi= sum(`Borrelia miyamotoi `))

#Borrelia Mayonii
Pathogen_data_Total_Borrelia_mayonii <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube) %>%
  summarise(total_Borrelia_Mayonii = sum(`Borrelia mayonii `))

#Babesia Microti
Pathogen_data_Total_microti <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Babesia_microti = sum(`Babesia microti`))

#Ehrlichia Muris
Pathogen_data_Total_ehrlichia <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Ehrlichia_muris = sum(`Ehrlichia muris`))

#Anaplasma Phagoctyophilum
Pathogen_data_Total_phagocytophilum <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Anaplasma_phagocytophilum = sum(`Anaplasma phagocytophilum`))

#Deer Tick Virus
Pathogen_data_Total_Powassan_DTV <- Pathogen_data_numeric_pain %>% 
  group_by(URI_Vial_Number_on_tube)  %>%
  summarise(total_Powassan_DTV = sum(`Deer Tick Virus`))

##make sure you use Co_infection_presents column, that column specifies is co-infections are present but not if it is a triple or secondary infection
# Number of Co-infections
Pathogen_data_total_coinfections <- Pathogen_data_numeric_pain %>%
  group_by(URI_Vial_Number_on_tube) %>%
  summarise(total_coinfected_individuals = sum(`Co_infection_present`))
            
#Merge all into single data
Pathogen_list <- list(Pathogen_data_Total_Ixodes_scapularis,
                      Pathogen_data_Total_Borrelia_burgdorferi,
                      Pathogen_data_Total_Borrelia_genus,
                      Pathogen_data_Total_Borrelia_miyamotoi,
                      Pathogen_data_Total_Borrelia_mayonii,
                      Pathogen_data_Total_microti,
                      Pathogen_data_Total_ehrlichia,
                      Pathogen_data_Total_phagocytophilum,
                      Pathogen_data_Total_Powassan_DTV,
                      Pathogen_data_total_infected_by_site,
                      Pathogen_data_total_coinfections)

Pathogen_combined_full <- Pathogen_list %>%
  reduce(full_join, by = 'URI_Vial_Number_on_tube') %>%
  relocate(total_infected_Ixodes, .after = total_Ixodes_Scapularis)

#Reorganize columns 
Pathogen_with_sites_codes_long_lat <- Pathogen_combined_full %>%
  rename(Vial_code = URI_Vial_Number_on_tube) %>%
  select(Vial_code:total_coinfected_individuals) 
 
Pathogen_with_sites_codes_and_all_information <- Pathogen_with_sites_codes_long_lat %>%
  left_join(Specimen_Date_Temp_Humidity_corrected, by = ('Vial_code')) %>%
  left_join(Sampling_Site_Longitude_and_Lattitude, by = ('Site_Code')) %>%
  relocate(Site_Code, .after = Vial_code) %>%
  relocate(`Longitude (x) -`, .after = Site_Code) %>%
  relocate(`Latitude (y) +`, .after = `Longitude (x) -`)


#Pathogen all by site
#Additionally, it provides a percentage of ticks at each site that are infected with pathogen

Pathogen_by_site_with_percentages <- Pathogen_with_sites_codes_and_all_information %>%
  select(c(Site_Code:Humdity)) %>%
  group_by(Site_Code) %>%
  summarise(total_Ixodes_Scapularis_nymphs = sum(total_Ixodes_Scapularis),
            total_infected_ixodes_scapularis_nymphs = sum(total_infected_Ixodes),
            total_Borrelia_burgdorferi = sum(total_Borrelia_burgdorferi),
            total_Borrelia_Mayonii = sum(total_Borrelia_Mayonii),
            total_Borrelia_miyamotoi = sum(total_Borrelia_miyamotoi),
            total_Babesia_microti = sum(total_Babesia_microti),
            total_Ehrlichia_muris = sum(total_Ehrlichia_muris),
            total_Anaplasma_phagocytophilum = sum(total_Anaplasma_phagocytophilum),
            total_Powassan_DTV = sum(total_Powassan_DTV),
            total_coinfection_individuals = sum(total_coinfected_individuals)) %>%
  mutate(percent_infected_nymphs = total_infected_ixodes_scapularis_nymphs/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Bb = total_Borrelia_burgdorferi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_miyamotoi = total_Borrelia_miyamotoi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_mayonii = total_Borrelia_Mayonii/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Babesia_microti = total_Babesia_microti/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Ehrlichia_muris = total_Ehrlichia_muris/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Anaplasma_phagocytophilum = total_Anaplasma_phagocytophilum/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Powassan_DTV = total_Powassan_DTV/total_Ixodes_Scapularis_nymphs*100,
         percent_total_coinfection_individuals = total_coinfection_individuals/total_Ixodes_Scapularis_nymphs*100)

#Pathogens all by month
#Additionally, it provides a percentage of ticks at each site that are infected with pathogen

Pathogen_by_month_with_percentages <- Pathogen_with_sites_codes_and_all_information %>%
  select(c(Vial_code:Humdity)) %>%
  group_by(month) %>%
  summarise(total_Ixodes_Scapularis_nymphs = sum(total_Ixodes_Scapularis),
            total_infected_ixodes_scapularis_nymphs = sum(total_infected_Ixodes),
            total_Borrelia_burgdorferi = sum(total_Borrelia_burgdorferi),
            total_Borrelia_Mayonii = sum(total_Borrelia_Mayonii),
            total_Borrelia_miyamotoi = sum(total_Borrelia_miyamotoi),
            total_Babesia_microti = sum(total_Babesia_microti),
            total_Ehrlichia_muris = sum(total_Ehrlichia_muris),
            total_Anaplasma_phagocytophilum = sum(total_Anaplasma_phagocytophilum),
            total_Powassan_DTV = sum(total_Powassan_DTV),
            total_coinfection_individuals = sum(total_coinfected_individuals)) %>%
  mutate(percent_infected_nymphs = total_infected_ixodes_scapularis_nymphs/total_Ixodes_Scapularis_nymphs*100,
    percent_infection_Bb = total_Borrelia_burgdorferi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_miyamotoi = total_Borrelia_miyamotoi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_mayonii = total_Borrelia_Mayonii/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Babesia_microti = total_Babesia_microti/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Ehrlichia_muris = total_Ehrlichia_muris/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Anaplasma_phagocytophilum = total_Anaplasma_phagocytophilum/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Powassan_DTV = total_Powassan_DTV/total_Ixodes_Scapularis_nymphs*100,
    percent_total_coinfection_individuals = total_coinfection_individuals/total_Ixodes_Scapularis_nymphs*100)


#Pathogens by month and by site
#Additionally, it provides a percentage of ticks at each site that are infected with pathogen

Pathogen_by_site_and_month_with_percentages <- Pathogen_with_sites_codes_and_all_information %>%
  select(c(Site_Code:Humdity)) %>%
  group_by(Site_Code, month) %>%
  summarise(total_Ixodes_Scapularis_nymphs = sum(total_Ixodes_Scapularis),
            total_infected_ixodes_scapularis_nymphs = sum(total_infected_Ixodes),
            total_Borrelia_burgdorferi = sum(total_Borrelia_burgdorferi),
            total_Borrelia_Mayonii = sum(total_Borrelia_Mayonii),
            total_Borrelia_miyamotoi = sum(total_Borrelia_miyamotoi),
            total_Babesia_microti = sum(total_Babesia_microti),
            total_Ehrlichia_muris = sum(total_Ehrlichia_muris),
            total_Anaplasma_phagocytophilum = sum(total_Anaplasma_phagocytophilum),
            total_Powassan_DTV = sum(total_Powassan_DTV),
            total_coinfection_individuals = sum(total_coinfected_individuals)) %>%
  mutate(percent_infected_nymphs = total_infected_ixodes_scapularis_nymphs/total_Ixodes_Scapularis_nymphs*100,
    percent_infection_Bb = total_Borrelia_burgdorferi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_miyamotoi = total_Borrelia_miyamotoi/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Borrelia_mayonii = total_Borrelia_Mayonii/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Babesia_microti = total_Babesia_microti/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Ehrlichia_muris = total_Ehrlichia_muris/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Anaplasma_phagocytophilum = total_Anaplasma_phagocytophilum/total_Ixodes_Scapularis_nymphs*100,
         percent_infection_Powassan_DTV = total_Powassan_DTV/total_Ixodes_Scapularis_nymphs*100,
    percent_total_coinfection_individuals = total_coinfection_individuals/total_Ixodes_Scapularis_nymphs*100)

#Join with environmental covariates
#can't do time here, you can only do site level

Pathogen_by_site_and_month_with_percentages <- Pathogen_by_site_and_month_with_percentages %>%
  left_join(Sampling_Site_Enviromental_Parameter, by = "Site_Code") 

  
Pathogen_by_site_with_percentages <- Pathogen_by_site_with_percentages %>%
  left_join(Sampling_Site_Enviromental_Parameter, by = "Site_Code")

#Relocate and select data.frames to make them more readable
Pathogen_by_site_and_month_with_percentages <- Pathogen_by_site_and_month_with_percentages %>%
  relocate(`Longitude (x) -`, .after = month) %>%
  relocate(`Latitude (y) +`, .after = `Longitude (x) -`)

Pathogen_by_site_with_percentages <- Pathogen_by_site_with_percentages %>%
  relocate(`Longitude (x) -`, .after = Site_Code) %>%
  relocate(`Latitude (y) +`, .after = `Longitude (x) -`)

#make a version of the site level counts with the buffers so it is only % forest cover for each site

Pathogen_by_site_spatial_analysis <- Pathogen_by_site_with_percentages %>%
  select(-c("percent_infected_nymphs", "percent_infection_Bb", "percent_infection_Borrelia_miyamotoi",
            "percent_infection_Borrelia_mayonii", "percent_infection_Babesia_microti","percent_infection_Ehrlichia_muris",
            "percent_infection_Anaplasma_phagocytophilum", "percent_infection_Powassan_DTV", "percent_total_coinfection_individuals",
  ))

Pathogen_by_site_spatial_analysis <- Pathogen_by_site_spatial_analysis %>%
  select(1:25) %>%
  select(-c("Site_Code"))

#only select columns we need
Pathogen_by_site_spatial_analysisz_no_buffers <- Pathogen_by_site_spatial_analysis %>%
  select(1,2,3,4,5,6,7,8,9,10,11,12,13,15,18,21,24)

#replace NA's with zeroes
Pathogen_by_site_spatial_analysisz_no_buffers[is.na(Pathogen_by_site_spatial_analysisz_no_buffers)] <- 0


#This section exports all previously collected data frames as .csv's

#CSV BY SITE regardless of month
write.csv(Pathogen_by_site_with_percentages, "../../Pathogen_by_site")

#Csv of month regardless of site
write.csv(Pathogen_by_month_with_percentages, "../../Pathogen_by_month_with_percentages")

#Csv of site and month
write.csv(Pathogen_by_site_and_month_with_percentages, "../../Pathogen_by_site_and_month_with_percentages")

##_________________________________________________________________________________________

##Let's do some exploratory analysis
#what does that entail
#how many nymphs are infected with each pathogen and overall
#how many have co-infections vs singular infections 


#We will start by counting the number of Ixodes infected with each pathogen
#for this purpose we will disregard borrelia genus as there seems to be few borrelia species beside borrelia burgdorferi
sum(Pathogen_by_month_with_percentages$total_Ixodes_Scapularis_nymphs) # 1566 amount
sum(Pathogen_by_month_with_percentages$total_infected_ixodes_scapularis_nymphs) #432 amount
sum(Pathogen_by_month_with_percentages$total_Borrelia_burgdorferi) #266 amount
sum(Pathogen_by_month_with_percentages$total_Borrelia_Mayonii) #0 amount
sum(Pathogen_by_month_with_percentages$total_Borrelia_miyamotoi) #17 amount
sum(Pathogen_by_month_with_percentages$total_Babesia_microti) #107 amount
sum(Pathogen_by_month_with_percentages$total_Ehrlichia_muris) #0 amount
sum(Pathogen_by_month_with_percentages$total_Anaplasma_phagocytophilum) #105 amount
sum(Pathogen_by_month_with_percentages$total_Powassan_DTV) #4 amount
sum(Pathogen_by_month_with_percentages$total_coinfection_individuals) #65 amount

#Bar Graph By Month with pathogen distribution

Pathogen_by_month_bar_graph_df <- Pathogen_by_month_with_percentages %>%
  rename(`Borrelia burgdorferi` = total_Borrelia_burgdorferi, 
         `Borrelia miyamotoi` = total_Borrelia_miyamotoi,
         `Babesia microti` = total_Babesia_microti,
         ` Anaplasma phagocytophilum`= total_Anaplasma_phagocytophilum,
         `Deer Tick Virus` = total_Powassan_DTV,
         `Borrelia Mayonii` = total_Borrelia_Mayonii,
         `Ehrlichia muris` = total_Ehrlichia_muris) %>%
  select(1,4,5,6,7,8,9,10)

Pathogen_by_month_bar_graph_df <- Pathogen_by_month_bar_graph_df %>%
  pivot_longer(!month, names_to = "Pathogen", values_to ="count")

Pathogen_by_month_bar_graph <- ggplot(Pathogen_by_month_bar_graph_df, aes(x = month, y = count, fill = Pathogen)) +
  geom_bar(stat="identity",position = "fill") 

Pathogen_by_month_bar_graph <- Pathogen_by_month_bar_graph + theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

Pathogen_by_month_bar_graph <- Pathogen_by_month_bar_graph + ggtitle("% of Collected Nymphs Infected With Pathogen") +
  xlab("Month") + ylab("% of I. scapularis Nymphs") +
  scale_x_discrete(breaks=c("5","6","7"),
                   labels=c("May", "June", "July"))


#Now lets see what pathogens make up those co-infections (we need to use the raw by individual data)
infection_amount_table <- table(Pathogen_data_numeric_pain$co_infections)
# 1142 No Infection
# 367 singular infection
# 63 double infection
# 2 triple infection

#bar graph for infection amount
infection_amount_df <- as.data.frame(infection_amount_table)

infection_amount_bar_graph <- ggplot(infection_amount_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity", fill = "black") 

infection_amount_bar_graph <- infection_amount_bar_graph + theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

infection_amount_bar_graph <- infection_amount_bar_graph + ggtitle("Frequency of Infection in I. scapularis nymphs") +
  xlab("Infection Type") + ylab("# of I. scapularis Nymphs") +
  scale_x_discrete(breaks=c("0","1","2","3"),
                   labels=c("No Infection",
                            "Single Infection",
                            "Double Infection",
                            "Triple Infection"))

## pathogen combinations for each co-infection

# triple infections, this was done by hand but if triple infections are high do what i did for double infections below
#but replace the infection_amount === 2 with infection_amount === 3 in the subset function

#BB/Bam/Ap = 2

#not going to make a bar graph for this since it has so few cases

#Double Infections

#make a data frame of just present pathogenes and nymphs of only double co-infections
Pathogen_data_numeric_pain_2_coinfection <- Pathogen_data_numeric_pain %>%
  select(6,7,9,11,12,13)

Pathogen_data_numeric_pain_2_coinfection <- subset(Pathogen_data_numeric_pain_2_coinfection, co_infections == 2)

#make table which combines all pathogens for double co-infection
pathogen_combinations <- with(Pathogen_data_numeric_pain_2_coinfection,
                              table(`Borrelia burgdorferi sensu lato`, 
                                    `Borrelia miyamotoi `,
                                    `Babesia microti`,
                                    `Anaplasma phagocytophilum`,
                                    `Deer Tick Virus`))

pathogen_combination_df <- as.data.frame(pathogen_combinations)

#bar graph for above data

#re-format data in pathogen_combination_df so it's easier to convert in a graph
#this isn't the most effecient way to do but it is fine for small data
Pathogen_combinations <- c("Bb/Bam","Bb/Ap","Bb/Bmy","Bam/Ap","Bmy/Ap","Ap/DTV")

Pathogen_amounts <- c(40, 15, 4, 2, 1, 1)

pathogen_combinations_and_amounts <- data.frame(Pathogen_combinations, Pathogen_amounts)

coinfection_bar_graph <- ggplot(pathogen_combinations_and_amounts, aes(x = reorder(Pathogen_combinations, -Pathogen_amounts), y = Pathogen_amounts)) +
  geom_bar(stat="identity", fill = "black")

coinfection_bar_graph <- coinfection_bar_graph + theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

coinfection_bar_graph <- coinfection_bar_graph + 
  ggtitle(label = "# of Infected Nymphs for Double Co-infections",
          subtitle = "Bb = Borrelia burgdorferi Bam = Babesia microti Ap = Anaplasmosis phagocytophilum 
          Bmy = Borrelia miyamotoi DTV = Deer Tick Virus") +
  xlab("Pathogen Combination") + 
  ylab("# of I. scapularis Nymphs") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

#the freq column in pathogen_combination_df will have the number of combinations for double infections

#Borrelia burgdorferi (Bb)/Babesia microti (Bam) = 40

#Bb/A. phagocytophilum (Ap) = 15

#Bb/B. miyamotoi = 4

#Bam/Ap = 2

#B. miyamotoi/Ap = 1

#Ap/Pow = 1

##odd how anaplasmosis and babesia have near equal infection but co-infection with burgdorferi nearly 3 times as high

#bar graph supplied below

#lets do a correlation matrix for every variable

correlation_nymphs <- cor(Pathogen_by_site_spatial_analysisz_no_buffers)

corrplot(correlation_nymphs, method = 'number', tl.cex = .5)

##there seems to be a high correlation between all pathogens in and all other pathogens
#unsuprising since we generalize based on site
#if we do a model for infection it needs to be at the individuals level rather than at the site level
#since we would be generalizing too much and losing statistical strength
#low correlation between enviromental factors and # of infected ticks by pathogen
# good candidates for covariates in model since we eliminate colinearity
#however, high correlation between all measurements for forest cover so we will only use
#forest cover percent for all forest types and not between the individual forest types

#lets do some histograms to test to see the distribution of pathogens

Burgdorferi_infected_nymphs <- Pathogen_by_site_spatial_analysis$total_Borrelia_burgdorferi
miyamotoi_infected_nymphs <- Pathogen_by_site_spatial_analysis$total_Borrelia_miyamotoi
microti_infected_nymphs <- Pathogen_by_site_spatial_analysis$total_Babesia_microti
anaplasma_infected_nymphs <- Pathogen_by_site_spatial_analysis$total_Anaplasma_phagocytophilum
coinfected_nymphs <- Pathogen_by_site_spatial_analysis$total_coinfection_individuals
infected_nymphs <- Pathogen_by_site_spatial_analysis$total_infected_ixodes_scapularis_nymphs



hist(Burgdorferi_infected_nymphs, breaks = "Sturges")
hist(miyamotoi_infected_nymphs)
hist(microti_infected_nymphs)
hist(anaplasma_infected_nymphs)
hist(coinfected_nymphs)
hist(infected_nymphs)

#they seem to be skewed to lower values
#given the binary response to infection and skew towards lower values at the site level
#models should be zero inflated or negative binomial maybe both 

###Lets do some quick descriptive map making to see if pathogen distribuion changes across rhode island

#it seems that there is indeed a consistent pattern for # of infected nymphs for each pathogen 

##Co-occurance of pathogens__________________________________________________________________________

#now that know co-infections appear with different pathogen combinations 
#lets see if there is any significance between those pathogen combinations
#make sure to add t() so you can transpose the data frame 
#the specific format for these functions need to site as columns and species as rows

#Subset data set to only pathogen co-occruance within nymphs
Pathogen_for_occur_analysis_by_individual_nymphs <- Pathogen_data_numeric %>%
  select(c(`Borrelia burgdorferi sensu lato`,`Borrelia miyamotoi `,
           `Borrelia mayonii `,`Babesia microti`,`Ehrlichia muris`,`Anaplasma phagocytophilum`,
           `Deer Tick Virus`)) %>%
  rename(`Borrelia burgdorferi` = `Borrelia burgdorferi sensu lato`) %>%
  t()

Pathogen_for_occur_analysis_by_individual_nymphs <- Pathogen_for_occur_analysis_by_individual_nymphs  %>%
  rename(`Borrelia burgdorferi` = `Borrelia burgdorferi sensu lato`)

#subset data to only measure co-occurance within site (not individuals)
#Rename column names for better presentation on graphs
Pathogen_for_occur_analysis_by_sites <- Pathogen_by_site_with_percentages %>%
  select(c(total_Borrelia_burgdorferi ,total_Borrelia_miyamotoi,total_Babesia_microti,
           total_Anaplasma_phagocytophilum,total_Powassan_DTV)) %>%
  rename(`Borrelia burgdorferi` = total_Borrelia_burgdorferi, 
         `Borrelia miyamotoi` = total_Borrelia_miyamotoi,
         `Babesia microti` = total_Babesia_microti,
         ` Anaplasma phagocytophilum`= total_Anaplasma_phagocytophilum,
         `Deer Tick Virus` = total_Powassan_DTV) %>%
  t()

##create objects for for measuring co-occurance of pathogens within seperate dimension


#presence/abscence of pathogens within individuals nymphs
pathogen_occurance_individuals <- cooccur(Pathogen_for_occur_analysis_by_individual_nymphs, spp_names = TRUE)


#prevalence of pathogens by sites

cooccur_pathogen_site_abundance <- cooccur(Pathogen_for_occur_analysis_by_sites, spp_names = TRUE)

##summary information
summary(pathogen_occurance_individuals)


summary(cooccur_pathogen_site_abundance)

#Co-infection between babesia microti and borrelia burgdorferi is significant in co-infection between individuals nymphs
#running the cooccur function on a dataframe will convert it to class coccur so when you plot
#it will autogenerate the species Co-occurrence Matrix
plot(pathogen_occurance_individuals, plotrand = TRUE)

#Presence of pathogens at the site level is signigant between all possible combinations of pathogens except deer tick virus

plot(cooccur_pathogen_site_abundance, plotrand = TRUE)

#presence/abscence of pathogens at sites 
#the code below will change all numbers above zero into 1, this will change it to abscencse/presence rather than abundance
#if you want the code line above this you need to re-run lines 507 - 510 to get the original object back
Pathogen_for_occur_analysis_by_sites[Pathogen_for_occur_analysis_by_sites > 0] <- 1 

pathogen_by_site_presence <- cooccur(Pathogen_for_occur_analysis_by_sites,spp_names = TRUE)

summary(pathogen_by_site_presence)

plot(pathogen_by_site_presence, plotrand = TRUE)


#Spatial Analysis of Pathogens______________________________________________________________________________

## Basic Plotting for each pathogen

Pathogen_by_site_spatial_analysis <- Pathogen_by_site_with_percentages %>%
  select(-c("percent_infected_nymphs", "percent_infection_Bb", "percent_infection_Borrelia_miyamotoi",
            "percent_infection_Borrelia_mayonii", "percent_infection_Babesia_microti","percent_infection_Ehrlichia_muris",
            "percent_infection_Anaplasma_phagocytophilum", "percent_infection_Powassan_DTV", "percent_total_coinfection_individuals",
  ))

Pathogen_by_site_spatial_analysis <- Pathogen_by_site_spatial_analysis %>%
  select(1:25) %>%
  select(-c("Site_Code"))

#only select columns we need which are
#spatial data, infection variables, and envriomental covariates
Pathogen_by_site_spatial_analysisz_no_buffers <- Pathogen_by_site_spatial_analysis %>%
  select(1,2,3,4,5,6,7,8,9,10,11,12,13,24)

#replace NA's with zeroes
Pathogen_by_site_spatial_analysisz_no_buffers[is.na(Pathogen_by_site_spatial_analysisz_no_buffers)] <- 0

burgdorferi_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_Borrelia_burgdorferi)) +
  ggtitle("Borrelia Burgdorferi Prevalence Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Infected Nymphs")

miyamotoi_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_Borrelia_miyamotoi)) +
  ggtitle("Borrelia Miyamotoi Prevalence Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Infected Nymphs")

microti_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_Babesia_microti)) +
  ggtitle("Babesia microti Prevalence Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Infected Nymphs")

anaplasma_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_Anaplasma_phagocytophilum)) +
  ggtitle("Anaplasma phagocytophilum Prevalence Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Infected Nymphs")

co_infected_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_coinfection_individuals)) +
  ggtitle("Co-infection Prevalence Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Co-Infected Nymphs")

infected_map <- ggplot() +
  geom_polygon(data = Rhode_Island, aes(x=long, y = lat, group = group), fill = "grey", alpha =0.8) +
  geom_point(data = Pathogen_by_site_spatial_analysisz_no_buffers,
             aes(x=Longitude, y=Latitude, size=total_infected_ixodes_scapularis_nymphs)) +
  ggtitle("Infected Nymphs Rhode Island") +
  theme(plot.title = element_text(size = 16, face = "bold", color = "black")) +
  labs(size = "# of Infected Nymphs")



#model should be either zero-inflated with binomial regression or poison regression
#im just using linear regression since i dont know how to work with those models


#Below is just simple linear models using total number of infected nymphs by pathogen for each site
#lattitude and longitude are included as covariates to account for spatial variation
#pathogen data is square root transformed to normalize data

#generalized linear model for #of infected nymphs as function of lat long using poisson distribution since data is discrete
mod0 <- glm( sqrt(total_Borrelia_burgdorferi) ~  Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod1 <- glm( sqrt(total_Borrelia_miyamotoi) ~  Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod2 <- glm( sqrt(total_Anaplasma_phagocytophilum) ~  Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod3 <- glm( sqrt(total_Babesia_microti) ~  Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod4 <- glm( sqrt(total_coinfection_individuals) ~ Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod5 <- glm( sqrt(total_coinfection_individuals) ~ Longitude + lattitude + total_Ixodes_Scapularis_nymphs, Pathogen_by_site_final, family = gaussian(link =identity))

#Pathogen prevalence as function of total ixodes scapularis nymphs collected
#here i use a gaussian distribution as the variable being predicted is continious
mod6 <- glm( sqrt(total_Borrelia_burgdorferi) ~ total_Ixodes_Scapularis_nymphs + Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod7 <- glm( sqrt(total_Borrelia_miyamotoi) ~ total_Ixodes_Scapularis_nymphs + Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod8 <- glm( sqrt(total_Anaplasma_phagocytophilum) ~ total_Ixodes_Scapularis_nymphs + Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod9 <- glm( sqrt(total_Babesia_microti) ~ total_Ixodes_Scapularis_nymphs + Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))
mod10 <-glm( sqrt(total_coinfection_individuals) ~ total_Ixodes_Scapularis_nymphs + Longitude + lattitude, Pathogen_by_site_final, family = gaussian(link =identity))


#Negative association between borrelia burgdroferi prevalence such that number of infected nymphs decreases with lattitude
summary(mod0)

#No association with borrelia miyamotoi prevalence with lattitude and longitude
summary(mod1)

#Association with anaplasma phagocytophilum and lattitude such that as lattitude increases, Ap prevalence decreases
summary(mod2)

#Significant relationship with babesia microti such that as longitude increases so does babesia microti prevalence
summary(mod3)

#this model suggests that as ixodes scapularis density increases so to does the number of co-infected nymphs 
summary(mod4)

#co-infections seems to occur more frequently in high density sites regardless of  longitude or lattitude
summary(mod5)

#However, co-nfected indivduals as a proportion of total collected specimens is not significant neither singularly 
#or with lattitude and longitude
summary(mod6)

# no signifigance with percent borrelia burgdorferi infection as function of total IS nymphs with long lat
summary(mod7)

# not signifigance with percent borrelia miyamotoi infection as function of total IS nymphs with long lat
summary(mod8)

#not signifigance with percent babesia microti infection as function of total IS nymphs with long lat
summary(mod9)

#not signifigance with percent anaplasma phagocytohilum infection as function of total IS nymohs
summary(mod10)

#start basic regression testing for each pathogen against borrelia burgdorferi to understand if borrelia burgdorferi prevalence can be associated with 
#prevalence for any other pathogen

zn.lm_1 <- glm(sqrt(total_Borrelia_burgdorferi) ~ sqrt(total_Borrelia_miyamotoi) , Pathogen_by_site_final, family = gaussian(link =identity))


zn.lm_2 <- glm(sqrt(total_Borrelia_burgdorferi) ~ sqrt(total_Anaplasma_phagocytophilum), Pathogen_by_site_final, family = gaussian(link =identity))


zn.lm_3 <- glm( sqrt(total_Borrelia_burgdorferi) ~ sqrt(total_Babesia_microti), Pathogen_by_site_final, family = gaussian(link =identity))

zn.lm_4 <- glm( sqrt(total_Borrelia_burgdorferi) ~ sqrt(total_Babesia_microti)+sqrt(total_Anaplasma_phagocytophilum)+sqrt(total_Borrelia_miyamotoi), Pathogen_by_site_final, family = gaussian(link =identity))

#no signifigance in using borrelia burgdorferi prevalence as a predictor of borrelia miyamotoi prevalence
summary(zn.lm_1)

#Borrelia burgdorferi is signifianct predictor for anaplasma phagocytophilum pathogen prevalence although R-squared is small
summary(zn.lm_2)

#no signifigance in using borrelia burgdorferi prevalence as a predictor of babesia microti prevalence by site
summary(zn.lm_3)

summary(zn.lm_4)

#_________________________________________________________________________________________________
## so far these have all been using infection when generalized across site so lets run some models
#using only the raw data by individuals

#rename columns so we can join
Specimen_Date_Temp_Humidity  <- Specimen_Date_Temp_Humidity %>%
  rename(URI_Vial_Number_on_tube = `VIAL CODE`, Site_Code = `SITE CODE`) 

#modified column so both are characters so join can work  
Specimen_Date_Temp_Humidity$URI_Vial_Number_on_tube <- as.character(Specimen_Date_Temp_Humidity$URI_Vial_Number_on_tube)
  
#left join all the individual data sheets so we can centralize into one data frame
pathogen_data_individual_enviromental_covariates <- Pathogen_data_numeric_pain %>%
  left_join(Specimen_Date_Temp_Humidity, by = join_by(URI_Vial_Number_on_tube)) %>%
  left_join(Sampling_Site_Enviromental_Parameter, by = join_by(Site_Code))

#now only select columns so their are no duplicates
pathogen_data_individual_enviromental_covariates <- pathogen_data_individual_enviromental_covariates %>%
  select(6,7,8,9,10,11,12,13,14,15,25,26,27,38)

#Now we can begin modelling

binomal_mod_1 <- glm(formula = `Borrelia burgdorferi sensu lato`~ TAME ,
            pathogen_data_individual_enviromental_covariates,
    family = binomial(link = "logit"))

binomal_mod_2 <- glm(formula = `Borrelia miyamotoi ` ~ TAME,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_3 <- glm(formula = `Babesia microti` ~ TAME,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_4 <- glm(formula = `Anaplasma phagocytophilum` ~ TAME,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_5 <- glm(formula = Co_infections_presents ~ TAME,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_6 <- glm(formula = `Borrelia burgdorferi sensu lato` ~ TAME + `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_7 <- glm(formula = `Borrelia miyamotoi ` ~ TAME + `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_8 <- glm(formula = `Babesia microti` ~ TAME + `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_9 <- glm(formula = `Anaplasma phagocytophilum` ~ TAME + `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_10 <- glm(formula = Co_infections_presents ~ TAME + `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_11 <- glm(formula = `Borrelia burgdorferi sensu lato` ~ TAME * `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_12 <- glm(formula = `Borrelia miyamotoi ` ~ TAME * `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_13 <- glm(formula = `Babesia microti` ~ TAME * `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_14 <- glm(formula = `Anaplasma phagocytophilum` ~ TAME * `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

binomal_mod_15 <- glm(formula = Co_infections_presents ~ TAME * `% Forest Cover 1 miles`,
            pathogen_data_individual_enviromental_covariates,
            family = binomial(link = "logit"))

summary(binomal_mod_1)
summary(binomal_mod_2)
summary(binomal_mod_3)
summary(binomal_mod_4)
summary(binomal_mod_5)
summary(binomal_mod_6)
summary(binomal_mod_7)
summary(binomal_mod_8)
summary(binomal_mod_9)
summary(binomal_mod_10)
summary(binomal_mod_11)
summary(binomal_mod_12)
summary(binomal_mod_13)
summary(binomal_mod_14)
summary(binomal_mod_15)


#spatial models of nymphal infection using gaussian copulas
#for this modelling we explore the spatial dependance of the data using spatial analysis principals
#we attempted to do this earlier using long and lat as covariaties but given the small size of rhode island this didn't give us any meaninful results
#this model will use a zero inflated distribution as we are focusing on infection for pathogens by site rather than in individuals annd
#based on our histograms our sites are skewed towards very low values
#additionally, given our sampling methodology the chance we miss an infected nymph is high
 
#summary tables
summary(Fitticks_miyamotoi_nobuffer)
summary(Fitticks_burgdorferi_nobuffer)
summary(Fitticks_microti_nobuffer)
summary(Fitticks_phagocytophilum_nobuffer)
summary(Fitticks_coinfection_nobuffer)
summary(Fitticks_infection_nobuffer)

summary(Fitticks_miyamotoi_nugget_nobuffer)
summary(Fitticks_burgdorferi_nugget_nobuffer)
summary(Fitticks_microti_nugget_nobuffer)
summary(Fitticks_phagocytophilum_nugget_nobuffer)
summary(Fitticks_coinfection_nugget_nobuffer)
summary(Fitticks_infection_nugget_nobuffer)

summary(Fitticks_miyamotoi_sph_nobuffer)
summary(Fitticks_burgdorferi_sph_nobuffer)
summary(Fitticks_microti_sph_nobuffer)
summary(Fitticks_phagocytophilumsph_sph_nobuffer)
summary(Fitticks_coinfection_sph_nobuffer)
summary(Fitticks_infection_sph_nobuffer)

summary(Fitticks_miyamotoi_nugget_sph_nobuffer)
summary(Fitticks_burgdorferi_nugget_sph_nobuffer)
summary(Fitticks_microti_nugget_sph_nobuffer)
summary(Fitticks_phagocytophilum_sph_nobuffer)
summary(Fitticks_coinfection_nugget_sph_nobuffer)
summary(Fitticks_infection_nugget_sph_nobuffer)

##Summary table by pathogen

#miyamotoi
summary(Fitticks_miyamotoi_nobuffer)
summary(Fitticks_miyamotoi_nugget_nobuffer)
summary(Fitticks_miyamotoi_sph_nobuffer)
summary(Fitticks_miyamotoi_nugget_sph_nobuffer)

#burgdorferi
summary(Fitticks_burgdorferi_nobuffer)
summary(Fitticks_burgdorferi_nugget_nobuffer)
summary(Fitticks_burgdorferi_sph_nobuffer)
summary(Fitticks_burgdorferi_nugget_sph_nobuffer)

#microti 
summary(Fitticks_microti_nobuffer)
summary(Fitticks_microti_nugget_nobuffer)
summary(Fitticks_microti_sph_nobuffer)
summary(Fitticks_microti_nugget_sph_nobuffer)

#anaplasmosis
summary(Fitticks_phagocytophilum_nobuffer)
summary(Fitticks_phagocytophilum_nugget_nobuffer)
summary(Fitticks_phagocytophilum_sph_nobuffer)
summary(Fitticks_phagocytophilum_nugget_sph_nobuffer)

#coinfection
summary(Fitticks_coinfection_nobuffer)
summary(Fitticks_coinfection_nugget_nobuffer)
summary(Fitticks_coinfection_sph_nobuffer)
summary(Fitticks_coinfection_nugget_sph_nobuffer)

#total infection model

summary(Fitticks_infection_nobuffer)
summary(Fitticks_infection_nugget_nobuffer)
summary(Fitticks_infection_sph_nobuffer)
summary(Fitticks_infection_nugget_sph_nobuffer)
