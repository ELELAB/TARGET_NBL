## From Simon (to check if packages are loaded and/or installed, otherwise it does it)
setwd(getwd())
source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("devtools","TCGAbiolinks","SummarizedExperiment", "ggplot2", "limma",
                      "reshape2", "plyr", "survival", "graphics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

source("src/SurvivalAnalysis_NEURO.R")

# query.exp <- GDCquery(project = "TARGET-NBL", 
#                      data.category = "Transcriptome Profiling",
#                      data.type = "Gene Expression Quantification",
####                   platform = "Illumina HiSeq",
#                      experimental.strategy = "RNA-Seq",
#                      file.type = "htseq")
# 
# GDCdownload(query.exp,
#             directory = "data/GDCdata")          
# neuroblastoma.exp <- GDCprepare(query = query.exp, 
#                                 summarizedExperiment = TRUE,
#                                 save = TRUE, 
#                                 save.filename = "data/NeuroblastomaExp.rda")


neuroblastoma.exp <- get(load("data/NeuroblastomaExp.rda"))

## To see dataset's features (class, dim, metadata, assays..)
neuroblastoma.exp 

View(assay(neuroblastoma.exp))
View(colData(neuroblastoma.exp))
dim(neuroblastoma.exp)

## This will allow you to have the conversion ENSEMBLE_ID - GENE NAME 
## (that it's better NOT TO DO dealing with HARMONIZED DATA)
View(rowData(neuroblastoma.exp))
      
## To access an assay dataset (an object can also have multiple assay datasets)
assays(neuroblastoma.exp)$HTSeq 

## To retrieve clinical data
neuroClin <- GDCquery_clinic(project = "TARGET-NBL", "clinical")

## How many samples are Primary Solid Tumor?
dataSmTP <- TCGAquery_SampleTypes(colnames(neuroblastoma.exp),"TP")
length(dataSmTP)

## How many samples are Solid Tissue Normal? 
## (WE DON'T HAVE THEM HERE)
dataSmNT <- TCGAquery_SampleTypes(colnames(neuroblastoma.exp),"NT")
length(dataSmNT)

## How many samples are Recurrent Solid Tumor?
dataSmTR <- TCGAquery_SampleTypes(colnames(neuroblastoma.exp),"TR")
length(dataSmTR)

## How many samples are Blood Derived Normal)?
## (WE DON'T HAVE THEM HERE)
dataSmNB <- TCGAquery_SampleTypes(colnames(neuroblastoma.exp),"NB")
length(dataSmNB)

## To set as ROW.NAMES the GENE NAMES and not the ENSEMBLE IDs
## Be careful: for harmonized data, it is better to keep the ENSEMBLE_IDs
# mygenesnames <- data.frame(Names = c(rowData(neuroblastoma.exp[,3])))
# View(mygenesnames)
# new.variable <- as.vector(mygenesnames$Names.external_gene_name)
# row.names(neuroblastoma.exp) <- new.variable

## PreProcessing, Normalization and Filtering

dataPrep <- TCGAanalyze_Preprocessing(object = neuroblastoma.exp,
                                      cor.cut = 0.6)
dim(dataPrep)

which(rownames(dataPrep)=="ENSG00000176840")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

## To see if after the normalization the MIR7-3HG is still in the dataset or has been removed
which(rownames(dataNorm)=="ENSG00000176840")   ## This is the ENSEMBLE_ID of MIR7-3HG  

## To check how many genes you lost after the Normalization step
dim(dataNorm)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

## To check how many genes you lost after the Filtering step
dim(dataFilt)

View(dataFilt)

## The user has the possibility to save this file.
# save(dataFilt, file = "data/Neuroblastoma_PreNormFilt.rda")
# dataFilt <- get(load("data/Neuroblastoma_PreNormFilt.rda"))


design.matrix <- model.matrix(~colData(neuroblastoma.exp)$tumor_stage)
v <- voom(dataFilt, design.matrix)
dataFilt_VOOM <- v$E

## REMEMBER!! After the filtering, the dataset is no longer a SummarizedExperiment.
## This means that you cannot use functions of SumExp package on it.

tumor.stage <- data.frame(neuroblastoma.exp$barcode,neuroblastoma.exp$tumor_stage)

colData(neuroblastoma.exp)$tumor_stage

length(which(colData(neuroblastoma.exp)$tumor_stage=="stage 2b"))
length(which(colData(neuroblastoma.exp)$tumor_stage=="stage 3"))
length(which(colData(neuroblastoma.exp)$tumor_stage=="stage 4"))
length(which(colData(neuroblastoma.exp)$tumor_stage=="stage 4s"))

## In this case, we can not perform BATCH_CORRECTION
## because we do not have the information .
# batch_TSS.mirna <-TCGAbatch_Correction2(dataFilt.mirna, batch.factor="TSS", adjustment = c())

neuroblastoma.exp.me <- melt(dataFilt_VOOM)

## We now create a dataset to have also the column of the tumor stage

neuroblastoma.me.AMBRA1 <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000110497'), ]
## ENSEMBLE_ID of AMBRA1

neuro_AMBRA1.stage <- cbind(neuroblastoma.me.AMBRA1, neuroblastoma.exp$tumor_stage)

## The user has the possibility to save this file.
# save(neuroblastoma.me.AMBRA1.stage, file = "data/AMBRA1_tumorstage.rda")
# neuro_AMBRA1.stage <- get(load("data/AMBRA1_tumorstage.rda"))

## BOXPLOT: AMBRA1 expression levels depending on Tumor_Stage

png('figs/AMBRA1_tumorstage_VOOM.png')
g <- ggplot(neuro_AMBRA1.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
g + geom_boxplot() + xlab("Tumor Stage") + ylab("AMBRA1 Expression") + ggtitle("Tumor Stage and AMBRA1 expression levels")
dev.off()

####

neuroblastoma.me.MIR7.3HG <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000176840'), ]
## ENSEMBLE_ID of MIR7-3HG

neuro_MIR7.3HG.stage <- cbind(neuroblastoma.me.MIR7.3HG, neuroblastoma.exp$tumor_stage)

## The user has the possibility to save this file.
# save(neuroblastoma.me.MIR7.3HG.stage, file = "MIR7.3HG_tumorstage.rda")
# neuro_MIR7.3HG.stage <- get(load("MIR7.3HG_tumorstage.rda"))

## BOXPLOT: MIR7-3HG expression levels depending on Tumor_Stage
 
png('figs/MIR7-3HG_tumorstage_VOOM.png')
g <- ggplot(neuro_MIR7.3HG.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
g + geom_boxplot() + xlab("Tumor Stage") + ylab("MIR7-3HG Expression") + ggtitle("Tumor Stage and MIR7-3HG expression levels")
dev.off()
 
####


neuroblastoma.me.MYCN <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000134323'), ]
## ENSEMBLE_ID of MYCN

neuro_MYCN.stage <- cbind(neuroblastoma.me.MYCN, neuroblastoma.exp$tumor_stage)


# save(neuroblastoma.me.cMYC.stage, file = "cMYC_tumorstage.rda")
# neuro_cMYC.stage <- get(load("cMYC_tumorstage.rda"))
 
## BOXPLOT: c-Myc expression levels depending on Tumor_Stage

png('figs/MYCN_tumorstage_VOOM.png')
g <- ggplot(neuro_MYCN.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
g + geom_boxplot() + xlab("Tumor Stage") + ylab("MYCN Expression") + ggtitle("Tumor Stage and MYCN expression levels")
dev.off()

####


## SURVIVAL ANALYSIS: we sourced at the beginning the script 'SurvivalAnalysis_NEURO.R'
## This script is in /src directory. 

## It is possible to try different threshold values (Up,Down) for the survival analysis.
## We are going to try both 0.25-0.75 and 0.50-0.50.

## For AMBRA1
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000110497", 0.25, 0.75)
survival_plot(clinical,"AMBRA1","25")

cox.model_AMBRA1 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_AMBRA1))
print(summary(cox.model_AMBRA1))

## For MIR7-3HG
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000176840", 0.25, 0.75)
survival_plot(clinical,"MIR7-3HG","25")

cox.model_MIR7.25 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_MIR7.25))
print(summary(cox.model_MIR7.25))


## For MYCN
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000134323", 0.25, 0.75)
survival_plot(clinical,"MYCN","25")
 
cox.model_MYCN.25 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_MYCN.25))
print(summary(cox.model_MYCN.25))


#############
## For AMBRA1
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000110497", 0.5, 0.5)
survival_plot(clinical,"AMBRA1","50")

cox.model_AMBRA1 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_AMBRA1))
print(summary(cox.model_AMBRA1))

## For MIR7-3HG
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000176840", 0.5, 0.5)
survival_plot(clinical,"MIR7-3HG","50")

cox.model_MIR7.25 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_MIR7.25))
print(summary(cox.model_MIR7.25))


## For MYCN
clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000134323", 0.5, 0.5)
survival_plot(clinical,"MYCN","50")

cox.model_MYCN.50 <- coxph(Surv(days_to_last_follow_up,status)~group+tumor_stage+gender,data=clinical)
## Test the proportionality
print(cox.zph(cox.model_MYCN.50))
print(summary(cox.model_MYCN.50))



