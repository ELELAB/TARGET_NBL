# Analysis of TARGET-NBL dataset to study the functional relationship among miR-7, AMBRA1 and n-Myc in neuroblastoma
We downloaded the data of neuroblastoma samples collected in the TARGET-NBL dataset, accessible through the [NCI Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/). The TARGET initiative - where TARGET stands for 'Therapeutically Applicable Research To Generate Effective Treatments' - aims at investigating the genomics of childhood cancers. The [TARGET-NBL dataset](https://ocg.cancer.gov/programs/target/projects/neuroblastoma) includes a comprehensive molecular characterization to determine the genetic changes driving the initiation and development of this aggressive cancer, accounting for ~12% of childhood cancer mortality.
In order to do this, we used an R/Bioconductor package, i.e. [TCGAbiolinks](https://www.ncbi.nlm.nih.gov/pubmed/26704973), initially developed for The Cancer Genome Atlas (TCGA) data but later adapted also to TARGET datasets.
Indeed, we:
- retrieved the TARGET-NBL dataset, integrating it with the available clinical data
- aggregated the data in the SummarizedExperiment
- pre-processed and normalized the data
- analyzed the available transcriptomics data to investigate the expression levels of AMBRA1 gene

**We developed to this aim the R script 'Neuroblastoma.R'**
## To run the analysis
### Prerequisites
# 
## Aknowledgments
Marta Lucchetta, Simo Mounir, Simon Kønig
## Author(s)
Maria Francesca Allega, Elena Papaleo
## License
This project is licensed under the MIT License - see the LICENSE.md file for details.

