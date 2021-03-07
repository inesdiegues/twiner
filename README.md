# twiner
Twiner is a network-based regularization parameter based on the correlation pattern between two data sets. 
In the present study, it was used to promote the selection of features that were correlated in DNA microarray and RNA-sequencing expression values data.
The higher the correlation between the genes in the two platforms, the lower the penalty term associated with it.

# Twiner code file contains:

1. Microarray data acquisition for GEO studies 
2. Microarrays pre processing steps - quantile normalization
3. Data integration for microarrays studies - ComBat for batch correction
4. RNA-sequencing data acquisition for TCGA study 
5. RNA-sequencing cross-platform integration with microarray - quantile normalization
6. Classification model with sparse logistic regression with:
    6.1. Elastic Net penalty 
    6.2. Twiner - Twin Networks Recovery - network-based regularization developed by us
7. Identification of cancer biomarkers

8. Images:
    8.1. PCA analysis
    8.2. Boxplots
    8.3. Violin plots
    8.4. Graph biological networks representation
    8.5. Survival analysis plots
    8.6. Venn diagrams
   
# Ensembl code contains the code developed to identify the genes that were protein coding
