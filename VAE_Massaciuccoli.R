################################################################################
#####                        VAE con Java                                  #####
################################################################################

############################### PREREQUISITE: VAE-MODEL JAR FILE DOWNLOADABLE FROM #########################################
#https://data.d4science.org/shub/E_RmlXSjJSbFVhZmVyT25YTFJJYlY1a3BJRWc0T0xueUVIOWNXamR3dStNV3RMZDl2WThJRE5rckY0b1cwWVU1Kw==#
############################################################################################################################

################################################################################
#####                        Training                                      #####
################################################################################

# Remove  all objects from R memory
rm(list=ls())

anno <- "2024_pr" #il 2017 è il training, per passare agli anni dopo sostituisci qui, in questo caso ho solo un anno di riferimento 


#input_file_path <-"Complete_dataset_mediterranean_sea_2021_2021_2021_2021_2021_2050RCP8.5.csv"
input_file_path <-paste0("input/risk1_overlap_standardized_no_na.csv") #questo è il file che contiene le variabili di input, in questo caso sono le variabili di rischio standardizzate. Le variabili non devono contenere NA#variable_names<- "environment 2021_land_distance,environment 2021_mean_depth"

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

valutazione <- FALSE #se è true allora si fa la valutazione del modello, se è false no, in questo caso non serve

number_of_hidden_nodes <- 10 #si parte con un numero di nodi nascosti pari alle variabili di input, ma va trovato il valore che evidenzi i dati anomali considerata la grana dell'analisi (in sostanza quanti dati anomali voglio individuare?)
number_of_epochs <- 1000
output_folder <- paste0("./output/VAE_varational_auto_encoder/out_2024_pr_n",number_of_hidden_nodes,"_test/")

model_folder <- paste0("./output/VAE_varational_auto_encoder/out_2024_pr_n",number_of_hidden_nodes,"_test/")
number_of_reconstruction_samples <- 16 #numero di campioni di ricostruzione da utilizzare durante il test, che può essere regolato in base alla dimensione del dataset e alla complessità del modello (in genere si lascia 16)
trained_model_file<-paste0(model_folder,"model_norm_49382X12_240056a46844f48936e43344446c54552153b71652737d46#12.bin") #questo è il file del modello pre-addestrato da utilizzare durante il test (si genera nella cartella di output), che deve essere specificato solo se si esegue il test, in questo caso è un file di esempio, ma va sostituito con il file del modello addestrato durante la fase di training
#prima di lanciare la fase di test/predizione ricordarsi di rinominare tutti i file nella cartella out come semplicemente "model"
training_mode_active<-"false"# #se è true allora si fa il training, se è false si fa la predizione/test, il training si fa sempre all'inizio


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
  #####                           valutazione                                #####
  ################################################################################
  
  
  file_pattern <- "classification_test_"
  files <- list.files(path = output_folder, pattern = paste0("^", file_pattern))
  if (length(files) == 1) {
    # Costruire il percorso completo del file
    file_path <- file.path(output_folder, files[1])
    
    # Leggi il file CSV
    data_projected <- read.csv(file_path,header = TRUE)
    # Visualizza i primi 6 dati
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
  
  data2 <- read.csv(paste0(output_folder,files), header=TRUE, sep=",")    # questo è il file che mi genera il vae come output
  
  data3 <- read.csv(input_file_path, header=TRUE, sep=",")    #questo è l'originale input multi k means
  
  data4 <- cbind(data3$longitude, data3$latitude, data2$reconstruction_log_probability)
  colnames(data4) <- c("x", "y", "reconstruction_log_probability")
  
  write.csv(data4, paste0(output_folder, "output_VAE_", anno, ".csv"), row.names = FALSE)
}