################################################################################
#####                        VAE con Java                                  #####
################################################################################

#Before running this script, make sure you have Java installed on your system.
#Download the VAE model jar file from the link provided below and place it in your working directory.
#https://data.d4science.org/shub/E_RmlXSjJSbFVhZmVyT25YTFJJYlY1a3BJRWc0T0xueUVIOWNXamR3dStNV3RMZDl2WThJRE5rckY0b1cwWVU1Kw==

################################################################################
#####                        Training                                      #####
################################################################################

# Remove  all objects from R memory
rm(list=ls())

anno <- "2024_pr" #in this case, there is only one reference year 


input_file_path <-paste0("input/risk1_overlap_standardized_no_na.csv") 

variable_names<- paste0( "out003_land_imperviousness_density_change_2024pr,",
                        "out009_tree_cover_density_2024pr,",
                        "out010_tree_cover_density_change_2024pr,",
                        "out011_grassland_2024pr,",
                        "out012_grassland_change_2024pr,",
                        "out017_land_use_and_cover_nature_2k_2024pr,",
                        "out040_plant_phenology_index_total_productivity_2024pr,",
                        "out057_temperature_avg_absolute_change_2024pr_rcp85,",
                        "out058_potential_evapotranspiration_relative_change_2024pr_rcp85,",
                        "out059_precipitation_cum_relative_change_2024pr_rcp85,",
                        "out065_number_species_cum_all_2024pr,",
                        "out147_land_use_and_cover_change_2024pr"
                        )

valutazione <- FALSE

number_of_hidden_nodes <- 10
number_of_epochs <- 1000
output_folder <- paste0("./output/VAE_varational_auto_encoder/out_2024_pr_n",number_of_hidden_nodes,"_test/")

model_folder <- paste0("./output/VAE_varational_auto_encoder/out_2024_pr_n",number_of_hidden_nodes,"_test/")
number_of_reconstruction_samples <- 16
trained_model_file<-paste0(model_folder,"model_norm_49382X12_240056a46844f48936e43344446c54552153b71652737d46#12.bin")
#Before launching the testing/prediction phase, remember to rename all files in the out folder simply as ‘model’.
training_mode_active<-"false"#if true then training is performed, if false then prediction/testing is performed, training is always performed at the beginning


if(training_mode_active=="true"){
  
  
  command_training<-paste0("java -cp vae.jar it.cnr.anomaly.JavaVAE -i\"./",input_file_path,"\" -v\"",variable_names,"\" -o\"",output_folder,"\" -h",number_of_hidden_nodes," -e",number_of_epochs," -r",number_of_reconstruction_samples," -t",training_mode_active)
  
  
  VAU_execution_train<-system(command_training, intern = TRUE,
                              ignore.stdout = FALSE, ignore.stderr = FALSE,
                              wait = TRUE, input = NULL, show.output.on.console = TRUE,
                              minimized = FALSE, invisible = TRUE)
  
  
  execution_train_success<-(length(which(grepl(pattern="OK VAU Training",x=VAU_execution_train)))>0)
  log_file <- paste0(output_folder,"log_file_training_",anno,".txt")
  writeLines(VAU_execution_train, log_file)
}else{
  dir.create(output_folder)
  ################################################################################
  #####                           Test                                       #####
  ################################################################################
  
  
  
  command_test <- paste0("java -cp vae.jar it.cnr.anomaly.JavaVAE -i\"./",input_file_path,"\" -v\"",variable_names,"\" -o\"",output_folder,"\" -r",number_of_reconstruction_samples," -t",training_mode_active," -m\"",trained_model_file,"\"")
  
  
  VAU_execution_test<-system(command_test, intern = T,
                             ignore.stdout = FALSE, ignore.stderr = FALSE,
                             wait = TRUE, input = NULL, show.output.on.console = TRUE,
                             minimized = FALSE, invisible = TRUE)
  
  execution_train_success<-(length(which(grepl(pattern="OK VAU Test",x=VAU_execution_test)))>0)
  log_file <- paste0(output_folder,"log_file_test_",anno,".txt")
  writeLines(VAU_execution_test, log_file)
  
  ################################################################################
  #####                           assessment                                #####
  ################################################################################
  
  
  file_pattern <- "classification_test_"
  files <- list.files(path = output_folder, pattern = paste0("^", file_pattern))
  if (length(files) == 1) {
    # Build the complete file path
    file_path <- file.path(output_folder, files[1])
    
    # Read the CSV file
    data_projected <- read.csv(file_path,header = TRUE)
    # Display the first 6 data items
    head(data)
  } else {
    cat("Più di un file trovato o nessun file trovato.")
  }
  namelist<- unlist(strsplit(variable_names, split = ","))
  data_projected_rdx <- data_projected[,namelist]
  
  data_input<-read.csv(input_file_path,header = TRUE)
  data_input <- data_input[,namelist]
  
  vettore_differenza <- data_projected_rdx - data_input
  vettore_differenza_vector <- unlist(vettore_differenza)
  vettore_differenza_numeric <- as.numeric(vettore_differenza_vector)
  errore <- mean((as.numeric(vettore_differenza_numeric))^2)
  
  
  rec_prob_avg<-mean(data_projected$reconstruction_log_probability)
  cat(paste0("error =",errore,", average probability recostruction =",rec_prob_avg),"\n")
  
  
  file_pattern <- "classification_test_"
  files <- list.files(path = output_folder, pattern = paste0("^", file_pattern))
  
  data2 <- read.csv(paste0(output_folder,files), header=TRUE, sep=",")    # this is the file that generates the vae as output
  
  data3 <- read.csv(input_file_path, header=TRUE, sep=",")    #this is the original input
  
  data4 <- cbind(data3$longitude, data3$latitude, data2$reconstruction_log_probability)
  colnames(data4) <- c("x", "y", "reconstruction_log_probability")
  
  write.csv(data4, paste0(output_folder, "output_VAE_", anno, ".csv"), row.names = FALSE)
}