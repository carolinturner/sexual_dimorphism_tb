# Scripts for sexual_dimorphism_tb paper

The scripts provided here were written by Dr Carolin Turner, Dr Jana Jiang and Dr James Greenan-Barrett, and used for analysis and figure plotting in the manuscript *"Human in vivo immunology of tuberculosis is not affected by sex dimorphism"* (link to be added upon publication).

`Analysis scripts` 1-10 require input data available from public repositories as listed below and from the Supplementary Tables in the manuscript, and produce most of the SourceData file.
`Plotting scripts` 1-5 make all figures and tables in the manuscript and require as input the SourceData file provided alongside the manuscript.

**Data availability:**
The following data objects are available from UCL's Research Data Repository (link to be added upon publication):

 - `Processed_TPM_postSVA_Blood_ActiveTB.csv` (required for script 1)
 - `Processed_Raw_Blood_ActiveTB.csv` (required for script 2)
 - `Meta_Blood_ActiveTB.csv` (required for scripts 1-2)
 - `Rawcounts_integrated-TST_transcriptome.csv` (required for script 4)
 - `Processed_sce_TSTD2_LatentTB_h5` (required for script 5-7)
 - `GSE326212_TB.rds` (required for script 8-10)

The following data objects are available from ArrayExpress with accession number [E-MTAB-14687](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-14687?query=E-MTAB-14687%20):

 - `tpm_PC0.001_log2_genesymbol_dedup.csv` (required for script 3)
 - `E-MTAB-14687.sdrf.txt` (required for scripts 3-4)