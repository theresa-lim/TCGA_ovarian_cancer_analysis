require(FirebrowseR)
cohorts = Metadata.Cohorts(format = "csv")
cancer.Type = cohorts[grep("OV", cohorts$cohort), 1]
rm(cohorts)

#gets all patients
all.Received = F
page.Counter = 1
page.size = 150
clinical_master = list()
while(all.Received == F){
  clinical_master[[page.Counter]] = Samples.Clinical(format = "csv",
                                               cohort = cancer.Type,
                                               page_size = page.size,
                                               page = page.Counter)
  if(page.Counter > 1)
    colnames(clinical_master[[page.Counter]]) = colnames(clinical_master[[page.Counter-1]])
  if(nrow(clinical_master[[page.Counter]]) < page.size){
    all.Received = T
  } else{
    page.Counter = page.Counter + 1
  }
}

clinical_master = do.call(rbind, clinical_master)
clinical_master = clinical_master[which(!clinical_master$tcga_participant_barcode == "TCGA-61-2008"), ]
clinical_master = clinical_master[which(!clinical_master$tcga_participant_barcode == "TCGA-29-1770"), ]
clinical_master = clinical_master[which(!clinical_master$tcga_participant_barcode == "TCGA-29-2414"), ]

gene_list = c("BRCA1", "BRCA2", "PTEN", "CLASP1","PAK1", "PAK2", "PAK4", "ROR1", "TP53") 
#gene_list1 = c("CDKN2A", "MAP2K4", "MAGEC1", "RIMBP2", "DIRAS3", "PEG3", "DAB2", "NF1", "ARID1A", "OPCML", "PLAGL1") from paper
vital_list = c("dead", "alive")
patients_all = list()

for(i in 1:2){
  patients_temp = clinical_master[which(clinical_master$vital_status == vital_list[[i]]), ]

  for(gene in gene_list){
    all.Received = F
    page.Counter = 1
    page.size = 600
    patients_mRNA = list()
    while(all.Received == F){
      patients_mRNA[[page.Counter]] = Samples.mRNASeq(format = "csv",
                                                    gene = gene,
                                                    cohort = cancer.Type,
                                                    tcga_participant_barcode = patients_temp$tcga_participant_barcode,
                                                    page_size = page.size,
                                                    page = page.Counter)
      if (nrow(patients_mRNA[[page.Counter]]) < page.size){
        all.Received = T
      } 
      else{
        page.Counter = page.Counter + 1
      }
    }
    patients_mRNA = do.call(rbind, patients_mRNA)
    colnames(patients_mRNA) = paste(gene, colnames(patients_mRNA), sep = ".")
    # print(paste("before", nrow(patients_temp), sep = ":"))
    patients_temp <- merge(patients_temp, patients_mRNA, by.x = "tcga_participant_barcode", 
                           by.y = paste(gene, "tcga_participant_barcode", sep = "."))
    # print(paste("after", nrow(patients_temp), sep = ":"))
  }
  patients_all[[i]] = patients_temp
}

patients_all = do.call(rbind, patients_all)

rm(patients_mRNA, patients_temp)

#IMPORT MUTATION DS???
#mutations <- Analyses.Mutation.MAF(format = "json", cohort = cancer.Type, gene = gene_list)

write.csv(patients_numeric, file = "patient_data.csv")



