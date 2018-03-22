# Analysis of TARGET-NBL dataset to study the functional relationship among miR-7, AMBRA1 and n-Myc in neuroblastoma
We downloaded the data of neuroblastoma samples collected in the TARGET-NBL dataset, accessible through the [NCI Genomic Data Commons Data Portal](https://portal.gdc.cancer.gov/). The TARGET initiative - where TARGET stands for 'Therapeutically Applicable Research To Generate Effective Treatments' - aims at investigating the genomics of childhood cancers. The [TARGET-NBL dataset](https://ocg.cancer.gov/programs/target/projects/neuroblastoma) includes a comprehensive molecular characterization to determine the genetic changes driving the initiation and development of this aggressive cancer, accounting for ~12% of childhood cancer mortality.
In order to do this, we used an R/Bioconductor package, i.e. [TCGAbiolinks](https://www.ncbi.nlm.nih.gov/pubmed/26704973), initially developed for The Cancer Genome Atlas (TCGA) data but later adapted also to TARGET datasets.
Indeed, we:
- retrieved the TARGET-NBL dataset, integrating it with the available clinical data
- aggregated the data in the SummarizedExperiment
- pre-processed and normalized the data
- analyzed the available transcriptomics data to investigate the expression levels of AMBRA1, MIR7-3HG and cMyc depending on the tumor stage
- performed the survival analysis for the same three targets based on 'high' or 'low' expression

**We developed to this aim the R script 'Neuroblastoma.R'.**
## To run the analysis
Required softwares:
```
R version 3.3.1 or higher
Rstudio version 1.1.383 or higher
Bioconductor version 3.6 or higher
```
Required packages (for the manual installation):

CRAN:
```
devtools (1.13.5)
ggplot2 (2.2.1)
graphics
plyr (1.8.4)
reshape2 (1.4.3)
survival (2.41-3)
```
Bioconductor:
```
Biobase (2.38.0)
BiocGenerics (0.24.0)
BiocInstaller (1.28.0)
DelayedArray (0.4.1)
GenomicRanges (1.30.3)
GenomInfoDb (1.14.0)
IRanges (2.12.0)
limma (3.34.9)
matrixStats (0.53.1)
SummarizedExperiment (1.8.1)
S4Vectors (0.16.0)
TCGAbiolinks (2.6.12)
```
### Effective running
Retrieve the script and the data from GitHub
```
$ git clone https://github.com/ELELAB/TARGET_NBL.git
```
Run the analysis
```
$ cd TARGET_NBL/
$ Rscript Neuroblastoma.R
```
**Running the script will perfom the analysis and source the script for the survival analysis from the /src directory.**

**In the dir /figs there will be located the output figures of both the scripts.**

## Notes 
The R script 'SurvivalAnalysis_NEURO.R' has been originally made by Marta Lucchetta.

This README file uses as a template the README file made by Simon Kønig.

## Aknowledgments
Marta Lucchetta, Simo Mounir, Simon Kønig

## Author(s)
Maria Francesca Allega, Elena Papaleo

## License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details.
