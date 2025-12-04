################################################################################
#####                  Overlap index and standardization                   #####
################################################################################

#Selection and standardisation of variables for risk calculation 1 ‘Risk of forest and grassland habitat degradation and associated biodiversity loss’

rm(list = ls())
#library(dplyr)
library(tidyverse)

#--- Step 1. Selection of variables for risk calculation 1 and export of the original table to CSV format ---####
all_variables <- read.csv("input/data_values.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
#head(all_variables)

# columns to select
variables_risk1_original <- all_variables [, c(
  "x",
  "y",
  "out003_land_imperviousness_density_change_2024pr",
  "out009_tree_cover_density_2024pr",
  "out010_tree_cover_density_change_2024pr",
  "out011_grassland_2024pr",
  "out012_grassland_change_2024pr",
  "out017_land_use_and_cover_nature_2k_2024pr",
  "out040_plant_phenology_index_total_productivity_2024pr",
  "out057_temperature_avg_absolute_change_2024pr_rcp85",
  "out058_potential_evapotranspiration_relative_change_2024pr_rcp85",
  "out059_precipitation_cum_relative_change_2024pr_rcp85",
  "out065_number_species_cum_all_2024pr",
  "out147_land_use_and_cover_change_2024pr"
)]   # <-- inserisci qui i codici che ti servono

# Export CSV variables_risk1_original

#cat("Exporting variables_risk1_original dataframe to CSV...\n")
#write.csv(variables_risk1_original, "output/Overlap/variables_risk1_original.csv", row.names = FALSE)
#cat("Export completed.\n")

#--- Step 2. Preprocessing: Reclassification of variable values ---####

variables_risk1_reclassified <- variables_risk1_original %>%
  mutate(
    out003_land_imperviousness_density_change_2024pr = case_when(
      out003_land_imperviousness_density_change_2024pr == 201 ~ 100, #unchanged areas with no imperviousness
      out003_land_imperviousness_density_change_2024pr == 254 ~ NA, #new imperviousness
      out003_land_imperviousness_density_change_2024pr == 255 ~ NA, #loss of imperviousness
      TRUE ~ out003_land_imperviousness_density_change_2024pr       # else remain the same
    ),
    
    out010_tree_cover_density_change_2024pr = case_when(
      out010_tree_cover_density_change_2024pr == 0  ~ 0,  #unchanged areas with no tree cover
      out010_tree_cover_density_change_2024pr == 1  ~ 1,  #new tree cover
      out010_tree_cover_density_change_2024pr == 2  ~ -1, #loss of tree cover
      out010_tree_cover_density_change_2024pr == 10 ~ 0,  #unchanged areas with tree cover
      TRUE ~ 0 #else
    ),
    
    out012_grassland_change_2024pr = case_when(
      out012_grassland_change_2024pr == 0  ~ 0,  #Unchanged non-grassland in both years
      out012_grassland_change_2024pr == 1  ~ 1,  #Grassland gain (added for completeness but not present in the data)
      out012_grassland_change_2024pr == 2  ~ -1, #Grassland loss
      out012_grassland_change_2024pr == 10 ~ 0,  #Unchanged grassland in both years
      out012_grassland_change_2024pr == 11 ~ 1,  #Unverified grassland gain
      out012_grassland_change_2024pr == 22 ~ -1, #Unverified grassland loss
      TRUE ~ 0 #else
    ),
    
    out017_land_use_and_cover_nature_2k_2024pr = case_when(
      out017_land_use_and_cover_nature_2k_2024pr %in% c(1120, 1310, 8120) ~ 1, #1120 Industrial, commercial and military units, 1310 Mineral extraction, dump and construction sites, 8120 Highly modified water courses and canals
      TRUE ~ 0 #else
    ),
    
    out147_land_use_and_cover_change_2024pr = case_when(
      out147_land_use_and_cover_change_2024pr == 1 ~ 1,   # Territori modellati artificialmente
      out147_land_use_and_cover_change_2024pr == 2 ~ 0.5, # Superfici agricole utilizzate
      out147_land_use_and_cover_change_2024pr == 5 ~ -1,  # Corpi idrici
      TRUE ~ 0 #else
    )
  )

# cat("Exporting variables_risk1_reclassified dataframe to CSV...\n")
# write.csv(variables_risk1_reclassified, "output/Overlap/variables_risk1_reclassified.csv", row.names = FALSE)
# cat("Export completed.\n")

#--- Step 3. Preprocessing: Inversion of variable values ---####

variables_risk1_reclassified_inverted <- variables_risk1_reclassified %>%
  mutate(
    out010_tree_cover_density_change_2024pr = 
      -1 * out010_tree_cover_density_change_2024pr,
    
    out012_grassland_change_2024pr = 
      -1 * out012_grassland_change_2024pr,
    
    out058_potential_evapotranspiration_relative_change_2024pr_rcp85 = 
      -1 * out058_potential_evapotranspiration_relative_change_2024pr_rcp85,
    
    out059_precipitation_cum_relative_change_2024pr_rcp85 = 
      -1 * out059_precipitation_cum_relative_change_2024pr_rcp85
  )


# print(unique(variables_risk1_reclassified$out010_tree_cover_density_change_2024pr))

# Esportazione CSV variables_risk1_reclassified

#cat("Exporting variables_risk1_reclassified_inverted dataframe to CSV...\n")
#write.csv(variables_risk1_reclassified_inverted, "output/Overlap/variables_risk1_reclassified_inverted.csv", row.names = FALSE)
#cat("Export completed.\n")

#--- Step 4. Standardisation of variables (Z-Score) ---####
# Standardizzazione delle variabili per ogni colonna: nuovo valore = (valore - media) / sd

variables_risk1_standardized <- variables_risk1_reclassified_inverted %>%
  mutate(
    across(
      c(
        "out003_land_imperviousness_density_change_2024pr",
        "out009_tree_cover_density_2024pr",
        "out040_plant_phenology_index_total_productivity_2024pr",
        "out057_temperature_avg_absolute_change_2024pr_rcp85",
        "out058_potential_evapotranspiration_relative_change_2024pr_rcp85",
        "out059_precipitation_cum_relative_change_2024pr_rcp85",
        "out065_number_species_cum_all_2024pr"
      ),
      ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)#è la formula per standardizzare ogni colonna, considerando eventuali NA (valori mancanti). Toglie la media da ogno valore e divide per la deviazione standard
    )
  )

# cat("Exporting variables_risk1_standardized dataframe to CSV...\n")
# write.csv(variables_risk1_standardized, "output/Overlap/variables_risk1_standardized.csv", row.names = FALSE)
# cat("Export completed.\n")

#--- Step 5. Calculation of the "overlap index" ---####

risk1_overlap_standardized <- variables_risk1_standardized %>%
  # Rischio totale = somma dei valori standardizzati delle variabili e poi normalizzati
  mutate(
    standard_risk_score = rowSums(across(-c("x", "y")), na.rm = TRUE) #esclude le colonne x e y e somma i valori standardizzati
    )

# Esportazione CSV
# cat("Exporting risk1_overlap_standardized dataframe to CSV...\n")
# write.csv(risk1_overlap_standardized, "output/Overlap/risk1_overlap_standardized.csv", row.names = FALSE)  # corretto qui: nome oggetto era sbagliato
# cat("Export completed.\n")

#--- Step 8. Selection of points without NA values ---####
risk1_overlap_standardized_no_na <- risk1_overlap_standardized %>%
  filter(complete.cases(.))

# Esportazione CSV risk1_no_na
cat("Exporting risk1_overlap_standardized_no_na dataframe to CSV...\n")
write.csv(risk1_overlap_standardized_no_na, "output/Overlap/risk1_overlap_standardized_no_na.csv", row.names = FALSE)
cat("Export completed.\n")

######################################################################################