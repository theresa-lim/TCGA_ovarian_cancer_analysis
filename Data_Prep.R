patients_numeric <- subset(patients_all, tcga_participant_barcode != "TCGA-59-A5PD" & !is.na(clinical_stage))
patients_numeric$PAK4dummy_var <- patients_numeric$PAK4.expression_log2
patients_numeric$lived_longer_than_3000 <- patients_numeric$days_to_death

for(row in 1:nrow(patients_numeric)){
  vital = patients_numeric[row, "vital_status"]
  if(vital == "alive"){
    patients_numeric[row, "vital_status"] = 0
  }
  else{
    patients_numeric[row, "vital_status"] = 1
  }
  
 #dummy variable for K-M and weibull diag#
  if(patients_numeric[row, "PAK4dummy_var"] < 11.34)
    patients_numeric[row,"PAK4dummy_var"] = 0
  else
    patients_numeric[row,"PAK4dummy_var"] = 1
  
  
  ###dummy variable for dependent logistic###
  if(!is.na(patients_numeric[row, "lived_longer_than_3000"]) && patients_numeric[row, "lived_longer_than_3000"] < 3000)
    patients_numeric[row, "lived_longer_than_3000"] = 0
  else
    patients_numeric[row, "lived_longer_than_3000"] = 1
  
  ###numerize stage###
  stage = patients_numeric[row, "clinical_stage"]
  if(stage == "stage iia"){
    patients_numeric[row, "clinical_stage"] = 2
  }
  else if(stage == "stage iib"){
    patients_numeric[row, "clinical_stage"] = 2
  }
  else if(stage == "stage iic"){
    patients_numeric[row, "clinical_stage"] = 2
  }
  else if(stage == "stage iiia"){
    patients_numeric[row, "clinical_stage"] = 3
  }
  else if(stage == "stage iiib"){
    patients_numeric[row, "clinical_stage"] = 3
  }
  else if(stage == "stage iiic"){
    patients_numeric[row, "clinical_stage"] = 3
  }
  else {
    patients_numeric[row, "clinical_stage"] = 4
  }

    
  ###calc days survived for alive###
  dtd = patients_numeric[row, "days_to_death"]
  if(is.na(dtd)){
    today = as.Date("01-01-2019", "%d-%m-%Y") #Case sensitive format
    day = patients_numeric[row, "day_of_form_completion"]
    month = patients_numeric[row, "month_of_form_completion"]
    year = patients_numeric[row, "year_of_form_completion"]
    date_of_diag = as.Date(paste(day, month, year, sep = "-"), "%d-%m-%Y")
    
    patients_numeric[row, "days_to_death"] = as.numeric(today - date_of_diag) #may/does skew data
  }
}
patients_numeric[, "vital_status"] = as.numeric(patients_numeric[, "vital_status"])
patients_numeric[, "clinical_stage"] = as.numeric(patients_numeric[, "clinical_stage"])

#Create ds of just dead patients
patients_dead <- patients_numeric[which(patients_numeric$vital_status == 0), ]

#to make step function work
patients_dead <- patients_dead[which(!is.na(patients_dead$clinical_stage)), ]

patients_alive <- patients_numeric[which(patients_numeric$vital_status == 1), ]
names(patients_alive)[17] <- paste("days_lived")

#patients_numeric$PAK4dummy_var = factor(patients_numeric$PAK4dummy_var, 
                                        #levels = c("0", "1"), 
                                        #labels = c("low", "high"))

#####removing outliers#####
patients_numeric <- patients_numeric[-which(patients_numeric$days_to_death > 4000), ]

write.csv(patients_numeric, file = "patient_data.csv")

