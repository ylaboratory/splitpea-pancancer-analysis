BiocManager::install("PoisonAlien/TCGAmutations")

sp_data <- read.table("tcga_splitpea_network.largest_ccs.txt", header = TRUE, sep = "\t")

cancer_values <- sp_data$cancer


cancer_values_processed <- tolower(gsub("_.*", "", cancer_values))

unique_cancer_values <- unique(cancer_values_processed)

unique_cancer_values <- setdiff(unique_cancer_values, "background")

unique_cancer_values

tcga_studies <- unique_cancer_values


tcga_studies <- unique_cancer_values


all_data <- list()

for (study in tcga_studies) {
  print(study)
  success <- FALSE
  tryCatch({
    tcga_data <- TCGAmutations::tcga_load(study = study)
    
    data <- as.data.frame(tcga_data@data)  
    
    data$cancer <- study  
    
    all_data[[study]] <- data  
    success <- TRUE
  }, error = function(e) {
    message(paste("Error loading study:", study, " - ", e))
  })
  
  if (!success) {
    next
  }
}

combined_data <- do.call(rbind, all_data)  

combined_data  

combined_data$Participant <- substr(combined_data$Tumor_Sample_Barcode, 1, 12)

head(combined_data$Participant)

colnames(combined_data) <- make.names(colnames(combined_data), unique = TRUE)

write.csv(combined_data, "mc3_mutation_data.csv", row.names = FALSE)