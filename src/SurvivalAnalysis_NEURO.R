##################
## Got this script from MARTA for the SURVIVAL ANALYSIS

## Function to get patient barcodes

get_patient <- function(list_barcodes){
  
  IDs <- strsplit(c(as.character(list_barcodes)), "-")
  
  IDs <- ldply(IDs, rbind)
  
  IDs$patient <- apply(IDs[,c(1,2,3)],1,paste,collapse = "-" )
  
  return(IDs$patient)
  
}

# Add the label "low" and "high"

get_type_sample <-function(x,patient_up,patient_down){
  
  if(x %in% patient_down){v <- "low"}
  
  else 
    
    if(x %in% patient_up)
      
    {v <-"high"}
  
}

get_survival_table <- function(dataFilt, clinical, gene, threshDown, threshUp){
  
  
## If in your dataset, you have different types of samples, you could be interested
## in getting only TUMOR SAMPLES 
# dataSmTP <- TCGAquery_SampleTypes(colnames(dataFilt),"TP")
# length(dataSmTP)
  
## Calculate quantile to be used as threshold to define two groups of samples to compare
  
  quantile_down <- quantile(dataFilt[gene,],threshDown)
  
  quantile_up <- quantile(dataFilt[gene,],threshUp)
  
  print(paste0("the lower percentile is ",quantile_down," and the higher percentile is ",quantile_up))
  
  samples_down <- colnames(dataFilt[,])[which(dataFilt[gene,]< quantile_down)] 
  
  samples_up <- colnames(dataFilt[,])[which(dataFilt[gene,]> quantile_up)]
  
  patient_down <- get_patient(samples_down)
  patient_up <- get_patient(samples_up)

  
  print(paste0("there are ",length(intersect(patient_up,patient_down))," patients in common between two groups:",intersect(patient_up,patient_down)))
  
  
  IDs_clinic <- strsplit(c(as.character(clinical$submitter_id)), "_")
  IDs_clinic <- ldply(IDs_clinic, rbind)
  clinical$patient <- IDs_clinic$`1`

  
  clinical <- subset(clinical, clinical$patient %in% union(patient_up,patient_down),
                     
                     select = c("patient","vital_status","days_to_last_follow_up","gender","tumor_stage","age_at_diagnosis"))
  
  
  clinical$group <- as.character(lapply(clinical$patient,function(x) get_type_sample(x,patient_up,patient_down)))
  
  clinical$age_at_diagnosis <- (clinical$age_at_diagnosis)/365
  
# clinical <- subset(clinical,!is.na(clinical$age_at_diagnosis))
# median <- floor(median(clinical$age_at_diagnosis))
# clinical$age <- as.character(lapply(clinical$age_at_diagnosis,function(x) get_age_range(x,median)))
# get vital_status: dead=1, alive=0
  
  clinical$status <- as.numeric(lapply(clinical$vital_status, FUN=function(x) if(x=="dead") v<-1 else v<-0))
  
  clinical <- subset(clinical, !is.na(clinical$days_to_last_follow_up))
  
  clinical$days_to_last_follow_up <- clinical$days_to_last_follow_up/365
  
  print(paste0("there are ",length(which(clinical$group=="low"))," patients in the low group"))
  
  print(paste0("and ",length(which(clinical$group=="high"))," patients in the high group"))
  
  
  
  return(clinical)
  
}

survival_plot <- function(clinical,gene,percentile){
  
  
  
  surv_object<-survfit(Surv(days_to_last_follow_up,status)~group, data=clinical)
  
  # Log-rang test
  
  sdf <- survdiff(Surv(days_to_last_follow_up,status)~group, data = clinical)
  
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  
  
  
  #pdf(paste0("survival_plot_",gene,"_",cancer_type,"_",percentile,".pdf"), width=9, height=7)
  
  pdf(paste0("figs/survival_plot_",gene,"_",percentile,".pdf"))
  
  SurvPlot <- plot(surv_object, main = paste0 ("Survival analysis for ",gene," Expression 
                                               \n Logrank p-value= ",signif(p.val)),
                   xlab = "Time (years)", ylab = "Surviving Rate",
                   lwd = 2.5, col=c("black","purple"))
  
              par(bg = "white", col.main = "black", mar = c(5,4.5,5.5,2) + 0.4)
              
              legend(11, 0.95, legend = c("high", "low"), fill = FALSE, 
                     border = c("black", "purple"), text.col = c("black","purple"), text.font = 2)
              
  dev.off()        
                     
}
