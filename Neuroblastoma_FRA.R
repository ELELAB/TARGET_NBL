## From Simon (to check packages)
setwd(getwd())
source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("devtools","TCGAbiolinks","SummarizedExperiment", "ggplot2", "limma", "reshape2", "plyr", "survival")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
sapply(list.of.packages, require, character.only = TRUE)

# source("SurvivalAnalysis_NEURO.R")


#library(devtools)
#document()
#build()
#install()
#library(TCGAbiolinks)
#rm(list=ls(all=TRUE))
#library(SummarizedExperiment)
#library(TCGAbiolinks)

query.exp <- GDCquery(project = "TARGET-NBL", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
###                   platform = "Illumina HiSeq",
                      experimental.strategy = "RNA-Seq",
                      file.type = "htseq")
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

## Which samples are primary solid tumor
# CHAAAAAAAAAANGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE!!!!!!!
# dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP")
# length(dataSmTP)
# 
# dataSmTP <- TCGAquery_SampleTypes(colnames(neuroblastoma.exp),"TP")
# length(dataSmTP)
# 
# ## which samples are solid tissue normal (THAT WE DON'T HAVE HERE)
# dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")
# length(dataSmNT)
# 
# ## which samples are Recurrent Solid Tumor
# dataSmTR <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TR")
# length(dataSmTR)

## which samples are Blood Derived Normal (THAT WE DON'T HAVE HERE)
# dataSmNB <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NB")
# length(dataSmNB)

## To set as ROW.NAMES the GENE NAMES and not the ENSEMBLE IDs
### mygenesnames <- data.frame(Names = c(rowData(neuroblastoma.exp[,3])))
### View(mygenesnames)
### new.variable <- as.vector(mygenesnames$Names.external_gene_name)
### row.names(neuroblastoma.exp) <- new.variable

dim(neuroblastoma.exp)

## PreProcessing, Normalization and Filtering

dataPrep <- TCGAanalyze_Preprocessing(object = neuroblastoma.exp,
                                      cor.cut = 0.6)
dim(dataPrep)

which(rownames(dataPrep)=="ENSG00000176840")

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

# To see if after the normalization the MIR7-3HG is still in the dataset or has been removed
which(rownames(dataNorm)=="ENSG00000176840")   ## This is the ENSEMBLE_ID of MIR7-3HG  

## If I select geneInfoHT, I get an error:
## Error in names(y) <- 1:length(y) : 
##'names' attribute [2] must be the same length as the vector [0]
## head(TCGAbiolinks::geneInfo) ..is it possible that there are NA data?

dim(dataNorm)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataFilt)

View(dataFilt)

save(dataFilt, file = "data/Neuroblastoma_PreNormFilt.rda")

dataFilt <- get(load("data/Neuroblastoma_PreNormFilt.rda"))

library("limma")
design.matrix <- model.matrix(~colData(neuroblastoma.exp)$tumor_stage)
v <- voom(dataFilt, design.matrix)
dataFilt_VOOM <- v$E

# REMEMBER!! After the filtering, the dataset is no more a SummarizedExperiment.
## This means that you cannot use functions of SumExp package on it.

tumor.stage <- data.frame(neuroblastoma.exp$barcode,neuroblastoma.exp$tumor_stage)

colData(neuroblastoma.exp)$tumor_stage
which(colData(neuroblastoma.exp)$tumor_stage=="stage 4")

## I DON'T NEED BATCH_CORRECTION, DO I?
## batch_TSS.mirna <-TCGAbatch_Correction2(dataFilt.mirna, batch.factor="TSS", adjustment = c())

library("reshape2")

neuroblastoma.exp.me <- melt(dataFilt_VOOM)

neuroblastoma.me.AMBRA1 <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000110497'), ]
## ENSEMBLE_ID of AMBRA1

neuroblastoma.me.AMBRA1.stage <- cbind(neuroblastoma.me.AMBRA1, neuroblastoma.exp$tumor_stage)

save(neuroblastoma.me.AMBRA1.stage, file = "data/AMBRA1_tumorstage.rda")

neuro_AMBRA1.stage <- get(load("data/AMBRA1_tumorstage.rda"))

# BOXPLOT: AMBRA1 expression levels depending on Tumor_Stage

library(ggplot2)

png('figs/AMBRA1_tumorstage_VOOM.png')
g <- ggplot(neuro_AMBRA1.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
g + geom_boxplot() + xlab("Tumor Stage") + ylab("AMBRA1 Expression") + ggtitle("Tumor Stage and AMBRA1 expression levels")
dev.off()

#####

# library("reshape2")
# 
# neuroblastoma.exp.me <- melt(dataFilt_VOOM)
# 
# neuroblastoma.me.MIR7.3HG <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000176840'), ]
# ## ENSEMBLE_ID of MIR7-3HG
# 
# neuroblastoma.me.MIR7.3HG.stage <- cbind(neuroblastoma.me.MIR7.3HG, neuroblastoma.exp$tumor_stage)
# 
# save(neuroblastoma.me.MIR7.3HG.stage, file = "MIR7.3HG_tumorstage.rda")
# 
# neuro_MIR7.3HG.stage <- get(load("MIR7.3HG_tumorstage.rda"))
# 
# # BOXPLOT: AMBRA1 expression levels depending on Tumor_Stage
# 
# library(ggplot2)
# 
# png('MIR7-3HG_tumorstage_VOOM.png')
# g <- ggplot(neuro_MIR7.3HG.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
# g + geom_boxplot() + xlab("Tumor Stage") + ylab("MIR7-3HG Expression") + ggtitle("Tumor Stage and MIR7-3HG expression levels")
# dev.off()
# 
# ####
# 
# library("reshape2")
# 
# neuroblastoma.exp.me <- melt(dataFilt_VOOM)
# 
# neuroblastoma.me.cMYC <- neuroblastoma.exp.me[which(neuroblastoma.exp.me$Var1=='ENSG00000136997'), ]
# ## ENSEMBLE_ID of cMYC
# 
# neuroblastoma.me.cMYC.stage <- cbind(neuroblastoma.me.cMYC, neuroblastoma.exp$tumor_stage)
# 
# save(neuroblastoma.me.cMYC.stage, file = "cMYC_tumorstage.rda")
# 
# neuro_cMYC.stage <- get(load("cMYC_tumorstage.rda"))
# 
# # BOXPLOT: AMBRA1 expression levels depending on Tumor_Stage
# 
# library(ggplot2)
# 
# png('cMYC_tumorstage_VOOM.png')
# g <- ggplot(neuro_cMYC.stage, aes(x=neuroblastoma.exp$tumor_stage, y=value)) 
# g + geom_boxplot() + xlab("Tumor Stage") + ylab("cMYC Expression") + ggtitle("Tumor Stage and cMYC expression levels")
# dev.off()
# 
# ####
# 
# clinical <- get_survival_table(dataFilt_VOOM, neuroClin, "ENSG00000110497", 0.25, 0.75)
# survival_plot(clinical,"ENSG00000110497","25")
# 
# 
# 
# 
# 
# 
# 
