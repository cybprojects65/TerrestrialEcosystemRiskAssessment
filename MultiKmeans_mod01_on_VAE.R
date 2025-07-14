rm(list = ls())
library(mclust)
library(tidyverse)
options(warn = -1)

# File to select
variables_risk1_standardized <- read.csv(file="./input/selection_high_risk_anomalies_0p_25p_10_hidden_nodes.csv", header=TRUE, sep=",")

# Remove rows with NA values
# variables_risk1_standardized <- variables_risk1_standardized %>%
#   drop_na()

# Selection of my study dataframe
# Bounding box coordinates
# bbxmin=-180
# bbxmax=180
# bbymin=-90 # latitude minimum
# bbymax=90

# Year of interest
year = "risk1_on_VAE_2024_pr"

# Selection of my study area and features
selection <- variables_risk1_standardized #[which((variables_risk1_standardized$longitude>=bbxmin)&(variables_risk1_standardized$longitude<=bbxmax)&(variables_risk1_standardized$latitude>=bbymin)&(variables_risk1_standardized$latitude<=bbymax)), ]

#select the first 300 lines to test the code=============Remember to comment
#selection <- selection %>% slice(1:100)

selected_features<-c("x",
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
)

cat(paste0("***Initialization***", "\n"))

# n_centroidi <- 7
#multi_centroidi<-  seq(from=20, to=30, by=2) #c(4)
multi_centroidi<-  seq(from=3, to=30, by=1) #c(4)
N <- 1000

bics<-c()
for (n_centroidi in multi_centroidi){
  
  selected_features_coords <- selection[,selected_features]
  # Study of centroids
  cat(paste0("Study of centroids", "\n"))
  # Prepare the dataset
  v <- as.data.frame(selected_features_coords[,3:ncol(selected_features_coords)])
  
  cat("####I'm analyzing ",n_centroidi,"centroids\n")
  
  # Create an EMPTY matrix of centroids to be filled later
  centroidi <- matrix(nrow=n_centroidi, ncol=ncol(v))
  
  for (centroide in 1:nrow(centroidi)) {
    
    centroidi[centroide,] <- as.numeric(v[centroide,])
    
  }
  
  km <- kmeans(as.matrix(v), centers = as.matrix(centroidi), iter.max = N)
  # Assign the column with the centroid value to the original dataset.
  selected_features_coords$distance_class <- km$cluster
  v$distance_class <- km$cluster
  
  centroidi<-as.matrix(km$centers)
  
  #ks.test(centroid_distribution, "punif")
  
  # Quantiles
  v_quantili <- apply(v, 2, quantile)
  
  # Prepare the centroids matrix with "M" for medium
  centroidi_labelled <- matrix("M", nrow=nrow(centroidi), ncol=ncol(centroidi))
  
  # Filling labeled centroids
  for (centroide in 1:nrow(centroidi)) {
    for (feat in 1:(ncol(v)-1)) {
      if (centroidi[centroide,feat]<v_quantili[3,feat]){
        centroidi_labelled[centroide, feat] <- "L"
      }
      else if (centroidi[centroide,feat]>v_quantili[4,feat]) {
        centroidi_labelled[centroide, feat] <- "H"
      }
      
    }
    
  }
  
  # Counting L, M, H to decide the level of risk
  c_H <- matrix(nrow = nrow(centroidi_labelled), ncol=1)
  c_M <- matrix(nrow = nrow(centroidi_labelled), ncol=1)
  c_L <- matrix(nrow = nrow(centroidi_labelled), ncol=1)
  
  # Creating empty vector for risk interpretation centroids
  centroide_interpretazione <- matrix(nrow = nrow(centroidi_labelled), ncol=1)
  
  # Counting the letters from centroidi_labelled
  for (r in 1:nrow(centroidi_labelled)) {
    c_H[r] <- sum(centroidi_labelled[r,] == "H")
    c_M[r] <- sum(centroidi_labelled[r,] == "M")
    c_L[r] <- sum(centroidi_labelled[r,] == "L")
  }  
  
  # Assignment of the risk level to the centroidi_interpretazione
  for (i in 1:nrow(centroide_interpretazione)) {
    if (c_H[i]>c_L[i] & c_H[i]>c_M[i]) {
      centroide_interpretazione[i] <- "high risk"
    }
    else if (c_L[i]>=c_H[i] & c_L[i]>c_M[i]) {
      centroide_interpretazione[i] <- "low risk"
    }
    else{ 
      centroide_interpretazione[i] <- "medium risk"
    }
  }
  
  #### Save and export the centroids with the interpretation of the risk level####
  centroidi_df <- as.data.frame(centroidi)
  centroidi_labelled_df <- as.data.frame(centroidi_labelled)
  
  # Nomina corretta delle colonne
  feature_names <- selected_features[3:length(selected_features)]
  names(centroidi_df) <- feature_names
  names(centroidi_labelled_df) <- paste0(feature_names, "_label")
  
  # ID centroidi
  centroid_id <- data.frame(centroid_id = 1:nrow(centroidi_df))
  
  # Costruiamo centroidi_annotated interponendo le colonne
  centroidi_annotated <- centroid_id  # iniziamo con ID
  
  for (i in seq_along(feature_names)) {
    num_col <- centroidi_df[[i]]
    label_col <- centroidi_labelled_df[[i]]
    
    centroidi_annotated[[feature_names[i]]] <- num_col
    centroidi_annotated[[paste0(feature_names[i], "_label")]] <- label_col
  }
  
  # Aggiungi l’interpretazione
  centroidi_annotated$risk_level <- centroide_interpretazione
  
  write.csv(centroidi_annotated, 
            file = paste0("./output/centroidi_annotated_", n_centroidi, "_", year, ".csv"), 
            row.names = FALSE)
  
  # Assign the distance class interpretation to the original dataset
  v$distance_class_interpretation <- NA
  for (i in 1:nrow(centroidi)) {
    indici <- which(v$distance_class == i)
    interpretazione <- centroide_interpretazione[i]
    v$distance_class_interpretation[indici] <- interpretazione
  }
  
  cat("Saving the dataset\n")
  nuovo_v <- cbind(selected_features_coords[,1:2],v)
  names(nuovo_v)<-c(names(selected_features_coords),"distance_class_interpretation")
  output_file<-paste0("./output/centroid_classification_assignment_",n_centroidi,"_", year, ".csv")
  write.csv(nuovo_v, output_file, row.names = F)
  
  centroid_distribution<-km$size
  #### CALCULATING ChiSqr
  if (length(which(centroid_distribution<=2))>0 || 
      ( (min(centroid_distribution)/max(centroid_distribution) ) <0.007) 
  ){
    cat("Unsuitable distribution: low uniformity:",(min(centroid_distribution)/max(centroid_distribution))," --- outliers: ",length(which(centroid_distribution<=2)),"\n")
    bic<-0
  }else{
    centroid_distribution.norm<-centroid_distribution/sum(centroid_distribution)
    reference<-rep(mean(centroid_distribution),length(centroid_distribution) )
    reference.norm<-reference/sum(reference)
    chisq<-sum((centroid_distribution.norm-reference.norm)^2/reference.norm)/length(centroid_distribution.norm)
    #high chisqr-> worse agreement with uniform distr
    #since we are selecting the maximum, let's invert the bic
    bic<-1/chisq
    
    #EXPLANATION OF THE CHI SQR CRITERION:
    #chi sqr probability calculation: for study purposes
    #the probability that the chisqr is lower than the calculcation has the inverse trend of the chisqr value
    #a small chisqr calculated means a higher probability of matching
    #a high chisqr calculated means a lower probability of matching
    #let's calculate the P(chisqr>chisqr_calculated), because the theoretical expected value of chisqr is 1
    #if this prob is high, chisqr_calculated is too far from the expected value-> non uniform distribution
    #if it is low then the chisqr_calculated is in line with the expected value-> more uniform distribution
    #uncomment for 
    #p_value <- pchisq(chisq, df = length(centroid_distribution.norm)-1, lower.tail = FALSE)
    #invert the criterion: high bic should be preferred because it corresponds to low p_value
    #bic<-1/p_value
    #cat("pvalue:",p_value,"\n")
    
    cat("Centroid distribution:",centroid_distribution.norm,"\n")
  }
  cat("Unif:",bic,"\n")
  bics<-c(bics,bic)
  cat("Done\n")
}


########################################################################################################
# Save the best clustering based on BIC
best_clusterisation<-multi_centroidi[which(bics == max(bics))]
cat("Ks: ",multi_centroidi,"\n")
cat("ChiSQRs: ",bics,"\n")
cat("Best clustering: K=",best_clusterisation,"\n")
best_clusterisation_file = paste0("./centroid_classification_assignment_",best_clusterisation,".csv")
cat("Best clustering file to take as result:",best_clusterisation_file,"\n")
plot(multi_centroidi,bics,type='l')
