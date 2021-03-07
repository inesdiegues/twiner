################# TWINER FOR TWO DIFFERENT GENE PROFILING TECHNIQUES ################# 
#
# Microarray studies integration 
# Cross-platform integration: microarray and RNA-seq studies
# Classification model with network-base regularization: Twiner
#
################################### LIBRARIES ################################### 

library(PRROC)
library(RColorBrewer)
library(DMwR)
library(ggfortify)
library(ggplot2)
library(gPCA)
library(preprocessCore)
library(qgraph)
library("propagate")
library("lsa")
library(glmnet)
library(survival)
library(survminer)
library("VennDiagram")
library("GEOquery")
#source("https://bioconductor.org/biocLite.R")
library(affy)
library(simpleaffy)
library("hgu133plus2.db")
library(limma)
library(sva)
library(edgeR)

################################### MICROARRAY ################################### 

## MICROARRAY DATA ACQUISITION________________________________________________________________

# Data downloaded from GEO_______________________________________________________

# STUDY 1
array_geo1 <- getGEO("GSE42568")
show(array_geo1)
array_matrix_geo1 <- exprs(array_geo1[[1]])
dim(array_matrix_geo1)

# STUDY 2
array_geo2 <- getGEO("GSE65194")
show(array_geo2)
array_matrix_geo2 <- exprs(array_geo2[[1]])
dim(array_matrix_geo2)

# 1- Create the probe to ensembl ids map (run ensembl code available for protein coding genes)

maps <- biomaRt::select(hgu133plus2.db, rownames(array_matrix_geo1), c("ENSEMBL", "GENENAME"))
maps2 <- biomaRt::select(hgu133plus2.db, rownames(array_matrix_geo2), c("ENSEMBL", "GENENAME"))

anyDuplicated(maps2$PROBEID)
anyDuplicated(maps2$ENSEMBL)

maps <- maps[!is.na(maps$ENSEMBL),]
maps2 <- maps[!is.na(maps2$ENSEMBL),]

# 2- Merge arrays with maps

######### STUDY 1
array_matrix <- as.data.frame(array_matrix_geo1)
array_matrix$PROBEID <- rownames(array_matrix)

array_matrix <- merge(array_matrix,maps, by = "PROBEID")
array_matrix$GENENAME <- NULL 
array_matrix$PROBEID <- NULL

array_ens <- aggregate(.~ENSEMBL,data=array_matrix,median) 
rownames(array_ens) <- array_ens$ENSEMBL
array_ens$ENSEMBL <- NULL
dim(array_ens) 

######### STUDY 2

array_matrix2 <- as.data.frame(array_matrix_geo2)
array_matrix2$PROBEID <- rownames(array_matrix2)

array_matrix2 <- merge(array_matrix2,maps2, by = "PROBEID")
array_matrix2$GENENAME <- NULL
array_matrix2$PROBEID <- NULL

array_ens2 <- aggregate(.~ENSEMBL,data=array_matrix2,median) 
rownames(array_ens2) <- array_ens2$ENSEMBL
array_ens2$ENSEMBL <- NULL
dim(array_ens2)

# 3 - Sample Information  

######### STUDY 1

labels <- pData(array_geo1[[1]]) 
patient_id <- labels$geo_accession
sample_identification <- labels$title
sample_identification <- substr(sample_identification,1,1)
sample_identification <- gsub("N", as.factor(0), sample_identification)
sample_identification <- gsub("B", as.factor(1), sample_identification)
sample_identification <- data.frame("Patient"= patient_id, "Label"=sample_identification)


######### STUDY 2

labels2 <- pData(array_geo2[[1]])


# 4 - Select ER+ samples

######### STUDY 1

ER_samples <- labels[labels$`er_status:ch1`== as.character(1),]
ER_samples <- ER_samples$geo_accession
ER_data <- array_ens[,colnames(array_ens) %in% ER_samples]
ER_data <- t(ER_data)

######### STUDY 2

ER_samples2 <- rbind(labels2[labels2$`sample_group:ch1` == "Luminal A",],labels2[labels2$`sample_group:ch1` == "Luminal B",]) 
ER_samples2 <- ER_samples2$geo_accession
ER_data2 <- array_ens2[,colnames(array_ens2) %in% ER_samples2]
ER_data2 <- t(ER_data2)


# 5 -  Select Normal samples

######### STUDY 1
normal_samples <- sample_identification[sample_identification$Label == "0",]
normal_samples <- normal_samples$Patient
normal_data <- array_ens[, colnames(array_ens) %in% normal_samples]
normal_data <- t(normal_data)

######### STUDY 2

normal_samples2 <- labels2[labels2$`sample_group:ch1` == "Healthy",]
normal_samples2 <- normal_samples2$geo_accession
normal_data2 <- array_ens2[,colnames(array_ens2) %in% normal_samples2]
normal_data2 <- t(normal_data2)

# 6 - Survival data [JUST FOR STUDY 1]

SA_days <- labels$`overall survival time_days:ch1`
SA_state <- labels$`overall survival event:ch1`
SA_days <- data.frame("Patient" = patient_id, "vital_status"= SA_state, "days_to_death" = as.numeric(SA_days))
SA_days <- SA_days[which(SA_days$Patient %in% as.factor(ER_samples)),]

# 7 - Join ER+ tumor and normal samples

######### STUDY 1

geo_array1 <- rbind(ER_data,normal_data)
dim(geo_array1)
y_geo_array_1<- c(rep(1,dim(ER_data)[1]),rep(0,dim(normal_data)[1]))
length(y_geo_array_1)

######### STUDY 2

geo_array2 <- rbind(ER_data2,normal_data2)
dim(geo_array2)
y_geo_array_2<- c(rep(1,dim(ER_data2)[1]),rep(0,dim(normal_data2)[1]))
length(y_geo_array_2)

# 8 - Protein coding genes and HUGO symbol 

######### STUDY 1

x_geo_array1 <- geo_array1[,which(colnames(geo_array1) %in% ensg.id)]
dim(x_geo_array1)

map_hugo <- data.frame("ensembl" = ensembl.protein.coding$gene_id,"hugo_symbol"=ensembl.protein.coding$gene_name)
map_hugo1 <- map_hugo[map_hugo$ensembl %in% colnames(x_geo_array1),]
map_hugo1 <- map_hugo1[order(map_hugo1$ensembl),]

x_geo_array1 <- x_geo_array1[,order(colnames(x_geo_array1))]
colnames(x_geo_array1) <- map_hugo1$hugo_symbol
range(x_geo_array1)
dim(x_geo_array1)
x_geo_array1 <- t(x_geo_array1)

######### STUDY 2

x_geo_array2 <- geo_array2[,which(colnames(geo_array2) %in% ensg.id)]
dim(x_geo_array2)

map_hugo2 <- map_hugo[map_hugo$ensembl %in% colnames(x_geo_array2),]
map_hugo2 <- map_hugo2[order(map_hugo2$ensembl),]

x_geo_array2 <- x_geo_array2[,order(colnames(x_geo_array2))]
colnames(x_geo_array2) <- map_hugo2$hugo_symbol
range(x_geo_array2)
dim(x_geo_array2)
x_geo_array2 <- t(x_geo_array2)


# 9 - Data visualization using Boxplots - IMAGES


n <- ncol(x_geo_array1)
col <- brewer.pal(n, "Paired")

boxplot(x_geo_array1,col=col, xlab = expression(bold("Patient")) , ylab = expression(bold("log2- gene expression value")) )

boxplot(x_geo_array2,col=col, xlab = expression(bold("Patient")) , ylab = expression(bold("log2- gene expression value")) )

# 10 - Within microarray normalization

######### STUDY 1

x_geo_array1_qn <- normalizeBetweenArrays(x_geo_array1, method="quantile")
range(log2(x_geo_array1_qn))
x_geo_array1_norm <- log2(x_geo_array1_qn)
boxplot(x_geo_array1_norm,col=col, xlab = expression(bold("Patient")) , ylab = expression(bold("log2- gene expression value")) )

######### STUDY 2

x_geo_array2_qn <- normalizeBetweenArrays(x_geo_array2, method="quantile")
range(x_geo_array2_qn)
x_geo_array2_norm <- log2(x_geo_array2_qn)
boxplot(x_geo_array2_norm,col=col, xlab = expression(bold("Patient")) , ylab = expression(bold("log2- gene expression value")) )


# CROSS-PLATFORM INTEGRATION_______________________________________________________

# Batch effect detection

# 11 - Select common genes between microarrays technologies 

comm_genes_arrays <- list()
comm_genes_arrays$geo1 <- rownames(x_geo_array1_norm)
comm_genes_arrays$geo2 <- rownames(x_geo_array2_norm)
common.symbols_arrays <- Reduce(intersect, comm_genes_arrays)

x_array_geo1 <- x_geo_array1_norm[which(rownames(x_geo_array1_norm) %in% common.symbols_arrays),]
x_array_geo1 <- x_array_geo1[order(rownames(x_array_geo1)),]
dim(x_array_geo1)
x_array_geo1<- x_array_geo1[-which(duplicated(rownames(x_array_geo1))),] 
x_array_geo1 <- t(x_array_geo1)

x_array_geo2 <- x_geo_array2_norm[which(rownames(x_geo_array2_norm) %in% common.symbols_arrays),]
x_array_geo2 <- x_array_geo2[order(rownames(x_array_geo2)),]
dim(x_array_geo2)
x_array_geo2 <- x_array_geo2[-which(duplicated(rownames(x_array_geo2))),]
x_array_geo2 <- t(x_array_geo2)

# 12 - PCAs of individualized studies - IMAGES

pca_geo1 <- as.data.frame(x_array_geo1)
group_label_geo1 <- gsub("1", "Tumor", y_geo_array_1)
group_label_geo1 <- gsub("0", "Normal", group_label_geo1)
pca_geo1$Group <- group_label_geo1

autoplot(prcomp(x_array_geo1, scale. = FALSE), data=pca_geo1, colour='Group') + theme(panel.grid.major = element_blank()) + theme_bw() + theme_classic()


pca_geo2 <- as.data.frame(x_array_geo2)
group_label_geo2 <- gsub("1", "Tumor", y_geo_array_2)
group_label_geo2 <- gsub("0", "Normal", group_label_geo2)
pca_geo2$Group <- group_label_geo2

autoplot(prcomp(x_array_geo2, scale. = FALSE),data=pca_geo2, colour='Group') +  theme(panel.grid.major = element_blank()) + theme_bw() + theme_classic()

# 13 - Cross-platform Normalization 

all_arrays <- rbind(x_array_geo1,x_array_geo2)
dim(all_arrays)

batch_label <- c(rep(1,dim(x_array_geo1)[1]), rep(2,dim(x_array_geo2)[1]))

array_batch_corrected <- ComBat(dat = t(all_arrays), batch = batch_label, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
array_batch_corrected <- t(array_batch_corrected)
dim(array_batch_corrected)


# 14 - Microarray final data

x_array <- array_batch_corrected
dim(x_array)

y_arrays <- c(rep(1, length(ER_samples)), rep(0, length(normal_samples)), rep(1, length(ER_samples2)), rep(0, length(normal_samples2)))
length(y_arrays)

# 15 - PCA for the merged data set - IMAGES

y_data_corrected <- c(rep("Tumor GEO1", length(ER_samples)), rep("Normal GEO1", length(normal_samples)), rep("Tumor GEO2", length(ER_samples2)), rep("Normal GEO2", length(normal_samples2)))

pca_plot <- as.data.frame(x_array)
pca_plot$Group <- y_data_corrected
autoplot(prcomp(x_array, scale. = FALSE),data=pca_plot, colour='Group')+ theme_bw() + theme_classic()


################################# DGE MICROARRAY ################################# 

design<- cbind(Tumor =1, Normal_vs_tumor=y_arrays)

dge_fit <- lmFit(t(x_array),design)
dge_fit <- eBayes(dge_fit)

table_arrays <- topTable(dge_fit,coef = "Normal_vs_tumor", n=dim(x_array)[2], adjust = "BH")

# DEGs 

dim(table_arrays[which(table_arrays$adj.P.Val<0.001 & abs(table_arrays$logFC)>1),])


##################################### RNASEQ #####################################
#
#
# RNA-SEQ DATA ACQUISITION____________________________________________________________

#library(devtools)
#devtools::install_url('https://github.com/averissimo/brca.data/releases/download/1.0/brca.data_1.0.tar.gz')

library(brca.data)

data('fpkm.per.tissue', 'fpkm.per.tissue.barcode', 'clinical', package = 'brca.data')
brca_tumor <- fpkm.per.tissue$primary.solid.tumor
dim(brca_tumor)

# 1 - Removing duplicated samples

brca_tumor <- t(brca_tumor[,-which(duplicated(getParticipantCode(colnames(fpkm.per.tissue$primary.solid.tumor))))]) 
dim(brca_tumor)

# 2 - ER+ samples

brca_tumor_erplus_id <- clinical$primary.solid.tumor[which(clinical$primary.solid.tumor[,51]=='Positive'),c(1,51)]

brca_tumor_erplus <- brca_tumor[which(getParticipantCode(rownames(brca_tumor)) %in% brca_tumor_erplus_id[,1]),]
dim(brca_tumor_erplus)

# 3 - Survival data 

rnaseq_survival_status <- clinical$primary.solid.tumor[,c(10,7)]
rnaseq_survival_status[which(rnaseq_survival_status[,2]=="Alive"),1] <- clinical$primary.solid.tumor[which(rnaseq_survival_status[,2]=="Alive"),11]
rnaseq_survival_status <- rnaseq_survival_status[-which(is.na(rnaseq_survival_status)),]

brca_tumor_erplus <- brca_tumor_erplus[which(getParticipantCode(rownames(brca_tumor_erplus)) %in% rownames(rnaseq_survival_status)),]
rnaseq_survival_status <- rnaseq_survival_status[which(rownames(rnaseq_survival_status) %in% getParticipantCode(rownames(brca_tumor_erplus))),]

# 4 - Normal samples

brca_normal <- fpkm.per.tissue$solid.tissue.normal
dim(brca_normal)

  # Let's consider only the normal samples taken from patients with ER+ tumors 

brca_normal_erplus_id <- clinical$solid.tissue.normal[which(clinical$solid.tissue.normal[,51]=='Positive'),c(1,51)]

brca_normal_erplus <- t(brca_normal[,which(getParticipantCode(colnames(brca_normal)) %in% brca_normal_erplus_id[,1])])
dim(brca_normal_erplus)

# 5 - Data for BRCA using RNA-SEQ

x_rnaseq <- rbind(brca_tumor_erplus,brca_normal_erplus)
dim(x_rnaseq)

# 6 - Protein coding genes and in Hugo Symbol

x_brca_rnaseq <- x_rnaseq[,which(colnames(x_rnaseq) %in% ensg.id)]
dim(x_brca_rnaseq)

mapping_ids <- data.frame("ensembl" = ensembl.protein.coding$gene_id,"hugo_symbol"=ensembl.protein.coding$gene_name)
mapping_ids <- mapping_ids[mapping_ids$ensembl %in% colnames(x_brca_rnaseq),]
mapping_ids <- mapping_ids[order(mapping_ids$ensembl),]

colnames(x_brca_rnaseq) <- mapping_ids$hugo_symbol
x_brca_rnaseq <- t(x_brca_rnaseq)
range(x_brca_rnaseq)

# 7 - Label vector 

y_brca_rnaseq <- c(rep(1,dim(brca_tumor_erplus)[1]),rep(0,dim(brca_normal_erplus)[1]))
length(y_brca_rnaseq)

################################# DGE RNA-SEQ ################################# 

design_rna <- cbind(Tumor = 1, Normal_vs_tumor = y_brca_rnaseq)
voom_transf <- voom(x_brca_rnaseq,design_rna, plot = TRUE)

dge_fit_rna <- lmFit(voom_transf,design_rna)
dge_fit_rna <- eBayes(dge_fit_rna)

table_rna <- topTable(dge_fit_rna, coef = ncol(design_rna), n=dim(x_brca_rnaseq)[2], adjust = "BH")

dim(table_arrays[which(table_rna$adj.P.Val<0.001 & abs(table_rna$logFC)>1),])

## RNASEQ DATA CROSS-PLATFORM NORMALIZATION__________________________________________

# 8 - Log2 Transformation

x_brca_rnaseq <- log2(x_brca_rnaseq+1)
dim(x_brca_rnaseq)
range(x_brca_rnaseq)


# 9 - Cross-platform Normalization - Quantile Normalization 

  # Common genes

x_array <- t(x_array)
dim(x_array)

comm_genes <- list()
comm_genes$array <- rownames(x_array)
comm_genes$seq <- rownames(x_brca_rnaseq)
common.symbols2 <- Reduce(intersect, comm_genes)

x_array <- x_array[which(rownames(x_array) %in% common.symbols2),]
x_array <- x_array[order(rownames(x_array)),]
dim(x_array)

x_brca_rnaseq <- x_brca_rnaseq[which(rownames(x_brca_rnaseq) %in% common.symbols2),]
x_brca_rnaseq <- x_brca_rnaseq[order(rownames(x_brca_rnaseq)),]
dim(x_brca_rnaseq)
x_brca_rnaseq <- x_brca_rnaseq[-which(duplicated(rownames(x_brca_rnaseq))),] #because it was given me some duplicated values...

x_array <- t(x_array)
x_brca_rnaseq <- t(x_brca_rnaseq)

# 10 - PCAs - IMAGES

  # PCA just for rna

pca_rna <- as.data.frame(x_brca_rnaseq)
group_label_rna <- gsub("1", "Tumor", y_brca_rnaseq)
group_label_rna <- gsub("0", "Normal", group_label_rna)
pca_rna$Group <- group_label_rna

autoplot(prcomp(x_brca_rnaseq, scale. = FALSE), data=pca_rna, colour="Group") +  theme(panel.grid.major = element_blank()) + theme_bw() + theme_classic()

  # Everything together 

together <- rbind(x_array,  x_brca_rnaseq)
label_color <- c(y_data_corrected, rep("Tumor RNA-seq",dim(brca_tumor_erplus)[1]),rep("Normal RNA-seq",dim(brca_normal_erplus)[1]))

pca_plot_all <- as.data.frame(together)
pca_plot_all$Group <- label_color
autoplot(prcomp(together, scale. = FALSE),data=pca_plot_all, colour='Group')+ theme(panel.grid.major = element_blank()) + theme_bw() + theme_classic()


# 11 - Quantile Normalization

rnaseq_qn <- normalize.quantiles.use.target(x = x_brca_rnaseq, target = as.vector(x_array), copy = TRUE)
rownames(rnaseq_qn) <- rownames(x_brca_rnaseq)
colnames(rnaseq_qn) <- colnames(x_brca_rnaseq)

# PCA - Normalized datasets - IMAGE

together2<- rbind(x_array, rnaseq_qn)
pca_plot_all2 <- as.data.frame(together2)
pca_plot_all2$Group <- label_color
autoplot(prcomp(together2, scale. = FALSE),data=pca_plot_all2, colour='Group')+ theme_bw() + theme_classic()


################################### 3. TWINER ################################### 
#
#
# FILTERING___________________________________________________________________________

# Remove variables with sd=0

array_sd <- x_array[,sapply(seq(ncol(x_array)), function(ix) {sd(x_array[,ix])}) != 0]
rnaseq_sd <- rnaseq_qn[,sapply(seq(ncol(rnaseq_qn)), function(ix) {sd(rnaseq_qn[,ix])}) != 0]

xarray_less <- array_sd[,which(colnames(array_sd) %in% colnames(rnaseq_sd))]
xrnaseq_less <- rnaseq_sd[,which(colnames(rnaseq_sd) %in% colnames(array_sd))]


## CLASSIFICATION PROBLEM____________________________________________________________

# Join the data from both technologies

data_join_balanced <- scale(rbind(xarray_less,xrnaseq_less))
dim(data_join_balanced)

sample_label <- c(y_arrays,y_brca_rnaseq) 

technique_label <-  c(rep("A",length(y_arrays)),rep("R",length(y_brca_rnaseq)))

# Data visualization of the scaled data

pca_plot_all2 <- as.data.frame(data_join_balanced)
pca_plot_all2$Group <- label_color
autoplot(prcomp(data_join_balanced, scale. = FALSE),data=pca_plot_all2, colour='Group')+ theme_bw() + theme_classic()

# Correlation Matrices

xarray_norm <- xarray_less
xrnaseq_norm <- xrnaseq_less

x_array_cor <- bigcor(xarray_norm, y = NULL, fun = "cor", size = 2000, verbose=FALSE)
x_array_cor <- as.data.frame(as.ffdf(x_array_cor))

x_rnaseq_cor <- bigcor(xrnaseq_norm, y = NULL, fun = "cor", size = 2000, verbose=FALSE)
x_rnaseq_cor <- as.data.frame(as.ffdf(x_rnaseq_cor))

# Angular distance vector

array_rnaseq_weights  <- vector()

for (i in 1:dim(data_join_balanced)[2]){
  array_rnaseq_weights [i] <- acos(cosine(x_array_cor[,i],x_rnaseq_cor[,i]))/pi
}

# Normalized weights (between 0 and 1)

array_rnaseq_weights_norm <- array_rnaseq_weights / max(array_rnaseq_weights)
array_rnaseq_weights <- array_rnaseq_weights_norm

# Building training and test sets (25 % assigned for testing)

set.seed(2019)

array_tumor_bl <- xarray_less[which(y_arrays == 1),] 
dim(array_tumor_bl)
array_normal_bl <- xarray_less[which(y_arrays == 0),]
dim(array_normal_bl)
rnaseq_tumor <- xrnaseq_less[which(y_brca_rnaseq == 1),]
dim(rnaseq_tumor)
rnaseq_normal <- xrnaseq_less[which(y_brca_rnaseq == 0),]
dim(rnaseq_normal)

test_id_array_tumor <- sample(1:dim(array_tumor_bl)[1],round(dim(array_tumor_bl)[1]*0.25),replace=FALSE)
test_id_array_normal <- sample(1:dim(array_normal_bl)[1],round(dim(array_normal_bl)[1]*0.25),replace=FALSE)
test_id_rnaseq_tumor <- sample(1:dim(rnaseq_tumor)[1],round(dim(rnaseq_tumor)[1]*0.25),replace=FALSE)
test_id_rnaseq_normal <- sample(1:dim(rnaseq_normal)[1],round(dim(rnaseq_normal)[1]*0.25),replace=FALSE)

y_array_tumor <- y_arrays[which(y_arrays==1)]
y_array_normal <- y_arrays[which(y_arrays==0)]
y_rnaseq_tumor <- y_brca_rnaseq[which(y_brca_rnaseq==1)]
y_rnaseq_normal <- y_brca_rnaseq[which(y_brca_rnaseq==0)]

xtrain <- rbind(array_tumor_bl[-test_id_array_tumor,],rnaseq_tumor[-test_id_rnaseq_tumor,],array_normal_bl[-test_id_array_normal,],rnaseq_normal[-test_id_rnaseq_normal,])
ytrain <- c(rep(1,length(c(y_array_tumor, y_rnaseq_tumor)) - length(c(test_id_array_tumor,test_id_rnaseq_tumor))),rep(0,length(c(y_array_normal, y_rnaseq_normal)) -  length(c(test_id_array_normal,test_id_rnaseq_normal))))

xtest <- rbind(array_tumor_bl[test_id_array_tumor,],rnaseq_tumor[test_id_rnaseq_tumor,], array_normal_bl[test_id_array_normal,],rnaseq_normal[test_id_rnaseq_normal,])
ytest <- c(rep(1,length(c(test_id_array_tumor,test_id_rnaseq_tumor))),rep(0,length(c(test_id_array_normal,test_id_rnaseq_normal))))

# Cross-validation to find the optimal alpha value for EN penalty

my_alpha <- seq(0.1,0.9,0.1)

# Optimization is going to be performed considering Logistic Regression with EN penalty

nvar_selected_EN <- matrix(0,1,length(my_alpha))
pred_EN_cv <- matrix(0,dim(xtrain)[1],length(my_alpha))
MSE_EN_cv <- matrix(0,1,length(my_alpha))
PR_EN_cv <- matrix(0,1,length(my_alpha))

# Assigning samples to folds, to be used in cross-validation when tunning alpha

set.seed(2010)
my_foldid <- sample(1:10,size=length(ytrain),replace=TRUE)

for (j in 1:length(my_alpha)){
  
  # Logistic model fitting with 10-fold cross-validation for glmnet:
  fit_EN_cv <- cv.glmnet(as.matrix(xtrain),as.factor(ytrain),family="binomial",nfolds=10,foldid=my_foldid,alpha=my_alpha[j],type.measure="mse")
  var_selected_EN <- which(fit_EN_cv$glmnet.fit$beta[,which(fit_EN_cv$cvm == min(fit_EN_cv$cvm))] != 0)
  nvar_selected_EN[j] <- length(var_selected_EN)
  
  # Predictions obtained by model i:
  pred_EN_cv[,j] <- predict(fit_EN_cv,as.matrix(xtrain),s=fit_EN_cv$lambda[which(fit_EN_cv$cvm == min(fit_EN_cv$cvm))],type="response")
  
  # Mean squared error of prediction (MSE), area under the precision recall curve:
  MSE_EN_cv[j] <- mean((ytrain-pred_EN_cv[,j])^2)
  PR_EN_cv[j] <- pr.curve(scores.class0=as.vector(round(pred_EN_cv[,j])), weights.class0=ytrain)$auc.integral
}

best_alpha_EN <- my_alpha[which(MSE_EN_cv == min(MSE_EN_cv))]


# Cross-validation to find the optimal alpha value for TWINER penalty

nvar_selected_twiner <- matrix(0,1,length(my_alpha))
pred_twiner_cv <- matrix(0,dim(xtrain)[1],length(my_alpha))
MSE_twiner_cv <- matrix(0,1,length(my_alpha))
PR_twiner_cv <- matrix(0,1,length(my_alpha))

# Assigning samples to folds, to be used in cross-validation when tunning alpha

set.seed(2010)
my_foldid <- sample(1:10,size=length(ytrain),replace=TRUE)

for (j in 1:length(my_alpha)){
  
  # Logistic model fitting with 10-fold cross-validation for glmnet:
  fit_twiner_cv <- cv.glmnet(as.matrix(xtrain),as.factor(ytrain),family="binomial",nfolds=10,foldid=my_foldid,alpha=my_alpha[j], penalty.factor=array_rnaseq_weights, type.measure="mse")
  var_selected_twiner <- which(fit_twiner_cv$glmnet.fit$beta[,which(fit_twiner_cv$cvm == min(fit_twiner_cv$cvm))] != 0)
  nvar_selected_twiner[j] <- length(var_selected_twiner)
  
  # Predictions obtained by model i:
  pred_twiner_cv[,j] <- predict(fit_twiner_cv,as.matrix(xtrain),s=fit_twiner_cv$lambda[which(fit_twiner_cv$cvm == min(fit_twiner_cv$cvm))],type="response")
  
  # Mean squared error of prediction (MSE), area under theprecision recall curve:
  MSE_twiner_cv[j] <- mean((ytrain-pred_EN_cv[,j])^2)
  PR_twiner_cv[j] <- pr.curve(scores.class0=as.vector(round(pred_EN_cv[,j])), weights.class0=ytrain)$auc.integral
}

best_alpha_twiner <- my_alpha[which(MSE_twiner_cv == min(MSE_twiner_cv))]

## Performing 100 train and test boostrap samples 

times_boot <- 100

nvar_selected_array_rnaseq_tn <- matrix(0,2,times_boot)
miscl_train_array_rnaseq_tn <- matrix(0,2,times_boot)
miscl_test_array_rnaseq_tn <- matrix(0,2,times_boot)
mse_train_array_rnaseq_tn <- matrix(0,2,times_boot)
mse_test_array_rnaseq_tn <- matrix(0,2,times_boot)
var_selected_array_rnaseq <- vector()
var_selected_array_rnaseq_netw <- vector()
var_selected_array_rnaseq_idx <- vector()
var_selected_array_rnaseq_netw_idx <- vector()
pr_auc_train_array_rnaseq_tn <- matrix(0,2,times_boot)
pr_auc_test_array_rnaseq_tn <- matrix(0,2,times_boot)


# Building random selected training and test sets (25 % assigned for testing)

set.seed(1979)

x_data_join <-data_join_balanced

new_test_id_matrix <- matrix(0,round(dim(x_data_join)[1]*0.25),times_boot) 

for (i in 1:times_boot){
  new_test_id_matrix[,i] <- sample(1:dim(x_data_join)[1],round(dim(x_data_join)[1]*0.25),replace=FALSE)
}

# Assigning samples to folds

set.seed(2010)

my_foldid_array_rnaseq_tn <- sample(1:10,size=dim(x_data_join[-new_test_id_matrix[,i],])[1]
                                    ,replace=TRUE)

#Extra analysis 0.9

best_alpha_twiner<- 0.9

for (i in 1:times_boot){
  
  new_xtrain <- x_data_join[-new_test_id_matrix[,i],]
  new_ytrain <- sample_label[-new_test_id_matrix[,i]]
  
  new_xtest <- x_data_join[new_test_id_matrix[,i],]
  new_ytest <- sample_label[new_test_id_matrix[,i]]
  
  
  # Classification by sparse logistic regression
  
  # WITH ELASTIC NET (EN) PENALTY
  
  fit_EN_cv_array_rnaseq_tn <- cv.glmnet(new_xtrain, as.factor(new_ytrain), family="binomial", nfolds=10, alpha= best_alpha_twiner, foldid=my_foldid_array_rnaseq_tn, type.measure="mse")
  
  var_selected_EN_array_rnaseq_tn_idx <- which(fit_EN_cv_array_rnaseq_tn$glmnet.fit$beta[,which(fit_EN_cv_array_rnaseq_tn$cvm == min(fit_EN_cv_array_rnaseq_tn$cvm))] != 0)
  
  nvar_selected_array_rnaseq_tn[1,i] <- length(var_selected_EN_array_rnaseq_tn_idx)
  
  var_selected_array_rnaseq_idx <- c(var_selected_array_rnaseq_idx,var_selected_EN_array_rnaseq_tn_idx)
  
  var_selected_array_rnaseq <- c(var_selected_array_rnaseq,colnames(new_xtrain[,var_selected_EN_array_rnaseq_tn_idx]))
  
  # Misclassifications, MSE, ROC and PR
  
  # train set
  pred_EN_train_array_rnaseq_tn <- predict(fit_EN_cv_array_rnaseq_tn,new_xtrain,type="response")
  table(new_ytrain,round(pred_EN_train_array_rnaseq_tn))
  miscl_train_array_rnaseq_tn[1,i] <- length(which(new_ytrain !=round(pred_EN_train_array_rnaseq_tn)))
  
  mse_train_array_rnaseq_tn[1,i] <- mean((new_ytrain - pred_EN_train_array_rnaseq_tn)^2)
  
  
  pr_auc_train_array_rnaseq_tn[1,i] <- pr.curve(scores.class0=as.vector(round(pred_EN_train_array_rnaseq_tn)), weights.class0=new_ytrain)$auc.integral
  
  # test set
  pred_EN_test_array_rnaseq_tn <- predict(fit_EN_cv_array_rnaseq_tn,new_xtest,type="response")
  table(new_ytest,round(pred_EN_test_array_rnaseq_tn))
  miscl_test_array_rnaseq_tn[1,i] <- length(which(new_ytest !=round(pred_EN_test_array_rnaseq_tn)))
  
  mse_test_array_rnaseq_tn[1,i] <- mean((new_ytest - pred_EN_test_array_rnaseq_tn)^2)
  
  pr_auc_test_array_rnaseq_tn[1,i] <- pr.curve(scores.class0=as.vector(round(pred_EN_test_array_rnaseq_tn)), weights.class0=new_ytest)$auc.integral
  
  # WITH NETWORK-BASED INFORMATION - TWINER
  
  fit_ENnetw_cv_array_rnaseq_tn <- cv.glmnet(new_xtrain, as.factor(new_ytrain), family="binomial", nfolds=10, alpha= best_alpha_twiner, foldid=my_foldid_array_rnaseq_tn, penalty.factor=array_rnaseq_weights,type.measure="mse")
  
  var_selected_ENnetw_array_rnaseq_tn_idx <- which(fit_ENnetw_cv_array_rnaseq_tn$glmnet.fit$beta[,which(fit_ENnetw_cv_array_rnaseq_tn$cvm == min(fit_ENnetw_cv_array_rnaseq_tn$cvm))] != 0)
  
  nvar_selected_array_rnaseq_tn[2,i] <- length(var_selected_ENnetw_array_rnaseq_tn_idx)
  
  var_selected_array_rnaseq_netw_idx <- c(var_selected_array_rnaseq_netw_idx,var_selected_ENnetw_array_rnaseq_tn_idx)
  
  var_selected_array_rnaseq_netw <- c(var_selected_array_rnaseq_netw,colnames(new_xtrain[,var_selected_ENnetw_array_rnaseq_tn_idx]))
  
  # Misclassifications, MSE, ROC and PR 
  
  # train set
  pred_ENnetw_train_array_rnaseq_tn <- predict(fit_ENnetw_cv_array_rnaseq_tn,new_xtrain,type="response")
  table(new_ytrain,round(pred_ENnetw_train_array_rnaseq_tn))
  miscl_train_array_rnaseq_tn[2,i] <- length(which(new_ytrain !=round(pred_ENnetw_train_array_rnaseq_tn)))
  
  mse_train_array_rnaseq_tn[2,i] <- mean((new_ytrain - pred_ENnetw_train_array_rnaseq_tn)^2)
  
  pr_auc_train_array_rnaseq_tn[2,i] <- pr.curve(scores.class0=as.vector(round(pred_ENnetw_train_array_rnaseq_tn)), weights.class0=new_ytrain)$auc.integral
  
  # test set
  pred_ENnetw_test_array_rnaseq_tn <- predict(fit_ENnetw_cv_array_rnaseq_tn,new_xtest,type="response")
  table(new_ytest,round(pred_ENnetw_test_array_rnaseq_tn))
  miscl_test_array_rnaseq_tn[2,i] <- length(which(new_ytest !=round(pred_ENnetw_test_array_rnaseq_tn)))
  
  mse_test_array_rnaseq_tn[2,i] <- mean((new_ytest - pred_ENnetw_test_array_rnaseq_tn)^2)
  
  pr_auc_test_array_rnaseq_tn[2,i] <- pr.curve(scores.class0=as.vector(round(pred_ENnetw_test_array_rnaseq_tn)), weights.class0=new_ytest)$auc.integral
  
  rm(new_xtrain,new_xtest,new_ytrain,new_ytest)
}

# STATISTIC VALUES___________________________________________________________________

# median number of variables selected

apply(nvar_selected_array_rnaseq_tn,1,median)

# median number of misclassifications in the train set

round(apply(round(miscl_train_array_rnaseq_tn),1,median))

# median number of misclassifications in the test set

round(apply(round(miscl_test_array_rnaseq_tn),1,median))

# median squared error in the train set

apply(mse_train_array_rnaseq_tn,1,median)

# median squared error in the test set

apply(mse_test_array_rnaseq_tn,1,median)

# median PR AUC in the train set

apply(pr_auc_train_array_rnaseq_tn,1,median)

# median PR AUC in the test set

apply(pr_auc_test_array_rnaseq_tn,1,median)

# improvement in model predictive performance (MSE) by incorporating network-based information

# train
median((mse_train_array_rnaseq_tn[1,]-mse_train_array_rnaseq_tn[2,])*100/mse_train_array_rnaseq_tn[1,])

#test
median((mse_test_array_rnaseq_tn[1,]-mse_test_array_rnaseq_tn[2,])*100/mse_test_array_rnaseq_tn[1,])


# GENES SELECTED___________________________________________________________________

# variables always selected by penalized logistic regression in the 100 bootstrap samples

var_selected_alw <- as.numeric(names(which(table(var_selected_array_rnaseq_idx)>99)))
var_selected_alw_code <- as.matrix(colnames(x_data_join[,var_selected_alw]))

# variables selected in more than 75% of the bootstrap runs by penalized logistic regression

var_selected_75 <- as.numeric(names(which(table(var_selected_array_rnaseq_idx)>75)))
var_selected_75_code <- as.matrix(colnames(x_data_join[,var_selected_75]))

# variables always selected by network-based penalized logistic regression in the 100 bootstrap samples

var_selected_netw_alw <- as.numeric(names(which(table(var_selected_array_rnaseq_netw_idx)>99)))
var_selected_netw_alw_code <- colnames(x_data_join[,var_selected_netw_alw])
 
# variables selected in more than 75% of the bootstrap runs by network-based penalized logistic regression

var_selected_netw_75 <- as.numeric(names(which(table(var_selected_array_rnaseq_netw_idx)>75)))
var_selected_netw_75_code <- as.matrix(colnames(x_data_join[,var_selected_netw_75]))

# variables selected in common between penalized logistic regression and network-based penalized logistic regression

common_var_selected_netw_75 <- var_selected_netw_75[which(var_selected_netw_75 %in% var_selected_75)]

common_var_selected_netw_75_code <- var_selected_netw_75_code[which(var_selected_netw_75_code %in% var_selected_75_code)]

# distinct variables selected by penalized logistic regression

distinct_var_selected_75 <- var_selected_75[-which(var_selected_75 %in% var_selected_netw_75)]

distinct_var_selected_75_code <- var_selected_75_code[-which(var_selected_75_code %in% var_selected_netw_75_code)]

# distinct variables selected by network-based penalized logistic regression

distinct_var_selected_netw_75 <- var_selected_netw_75[-which(var_selected_netw_75 %in% var_selected_75)]

distinct_var_selected_netw_75_code <- var_selected_netw_75_code[-which(var_selected_netw_75_code %in% var_selected_75_code)]

## MICROARRAY ALONE_________________________________________________________________


# Building training and test sets (25 % assigned for testing)

set.seed(2019)

test_id_array_tumor <- sample(1:dim(array_tumor_bl)[1],round(dim(array_tumor_bl)[1]*0.25),replace=FALSE)
test_id_array_normal <- sample(1:dim(array_normal_bl)[1],round(dim(array_normal_bl)[1]*0.25),replace=FALSE)

xtrain <- rbind(array_tumor_bl[-test_id_array_tumor,],array_normal_bl[-test_id_array_normal,])
ytrain <- c(rep(1,length(y_array_tumor) - length(test_id_array_tumor)),rep(0,length(y_array_normal) -  length(test_id_array_normal)))

xtest <- rbind(array_tumor_bl[test_id_array_tumor,], array_normal_bl[test_id_array_normal,])
ytest <- c(rep(1,length(test_id_array_tumor)),rep(0,length(test_id_array_normal)))

# CLASSIFICATION BY SPARSE LOGISTIC REGRESSION BASED ON EN PENALTY_______________

# Cross-validation for finding optimum alpha

my_alpha <- seq(0.5,0.9,0.1)

# Optimizing alpha and lambda for logistic regression with Elastic net regularization

nvar_selected_EN_array <- matrix(0,1,length(my_alpha))
pred_EN_cv_array <- matrix(0,dim(xtrain)[1],length(my_alpha))
MSE_EN_cv_array <- matrix(0,1,length(my_alpha))
PR_EN_cv_array <- matrix(0,1,length(my_alpha))

# Assigning samples to folds, to be used in cross-validation when tunning alpha

set.seed(2010)
my_foldid_array <- sample(1:10,size=length(ytrain),replace=TRUE)

# Classification glmnet

for (j in 1:length(my_alpha)){
  
  # Logistic model fitting with 10-fold cross-validation for glmnet:
  fit_EN_cv_array <- cv.glmnet(as.matrix(xtrain),as.factor(ytrain),family="binomial",foldid=my_foldid_array,alpha=my_alpha[j],type.measure="mse")
  var_selected_EN_array <- which(fit_EN_cv_array$glmnet.fit$beta[,which(fit_EN_cv_array$cvm == min(fit_EN_cv_array$cvm))] != 0)
  nvar_selected_EN_array[j] <- length(var_selected_EN_array)
  
  # Predictions obtained by model i:
  pred_EN_cv_array[,j] <- predict(fit_EN_cv_array,as.matrix(xtrain),s=fit_EN_cv_array$lambda[which(fit_EN_cv_array$cvm == min(fit_EN_cv_array$cvm))],type="response")
  
  # Mean squared error of prediction (MSE), area under theprecision recall curve:
  MSE_EN_cv_array[j] <- mean((ytrain-pred_EN_cv_array[,j])^2)
  PR_EN_cv_array[j] <- pr.curve(scores.class0=as.vector(round(pred_EN_cv_array[,j])), weights.class0=ytrain)$auc.integral
}

best_alpha_array <- my_alpha[which(MSE_EN_cv_array == min(MSE_EN_cv_array))]

# Sparse logistic regression

fit_EN_cv_array <- cv.glmnet(xtrain, as.factor(ytrain), family="binomial", nfolds=10, alpha=best_alpha_array,foldid=my_foldid_array, type.measure="mse")

var_selected_EN_array_idx <- which(fit_EN_cv_array$glmnet.fit$beta[,which(fit_EN_cv_array$cvm == min(fit_EN_cv_array$cvm))] != 0)

length(var_selected_EN_array_idx)

# GENES SELECTED_________________________________________________________________

var_selected_EN_array <- as.matrix(colnames(xtrain[,var_selected_EN_array_idx]))

# STATISTIC PARAMETERS___________________________________________________________

## misclassifications, MSE ROC and PR

# train set

pred_EN_train_array <- predict(fit_EN_cv_array,xtrain,s="lambda.min",type="response")
table(ytrain,round(pred_EN_train_array))

mse_EN_train_array <- mean((ytrain - pred_EN_train_array)^2)

PR_EN_train_array <- pr.curve(scores.class0=as.vector(round(pred_EN_train_array)), weights.class0=ytrain)$auc.integral

# test set 

pred_EN_test_array <- predict(fit_EN_cv_array,xtest,s="lambda.min",type="response")
table(ytest,round(pred_EN_test_array))

mse_EN_test_array <- mean((ytest - pred_EN_test_array)^2)

PR_EN_test_array <- pr.curve(scores.class0=as.vector(round(pred_EN_test_array)), weights.class0=ytest)$auc.integral

## RNA-SEQ ALONE______________________________________________________________________

xrnaseq_fit_norm <- xrnaseq_norm

# building training and test sets (25 % assigned for testing)

set.seed(2019)

test_id_rnaseq_tumor <- sample(1:dim(rnaseq_tumor)[1],round(dim(rnaseq_tumor)[1]*0.25),replace=FALSE)
test_id_rnaseq_normal <- sample(1:dim(rnaseq_normal)[1],round(dim(rnaseq_normal)[1]*0.25),replace=FALSE)

xtrain <- rbind(rnaseq_tumor[-test_id_rnaseq_tumor,],rnaseq_normal[-test_id_rnaseq_normal,])
ytrain <- c(rep(1,length(y_rnaseq_tumor) - length(test_id_rnaseq_tumor)),rep(0,length(y_rnaseq_normal) -  length(test_id_rnaseq_normal)))

xtest <- rbind(rnaseq_tumor[test_id_rnaseq_tumor,], rnaseq_normal[test_id_rnaseq_normal,])
ytest <- c(rep(1,length(test_id_rnaseq_tumor)),rep(0,length(test_id_rnaseq_normal)))

# CLASSIFICATION BY SPARSE LOGISTIC REGRESSION BASED ON EN PENALTY_______________

# Cross-validation for finding optimum alpha

# Optimizing alpha and lambda for logistic regression with Elastic net regularization

nvar_selected_EN_rnaseq <- matrix(0,1,length(my_alpha))
pred_EN_cv_rnaseq <- matrix(0,dim(xtrain)[1],length(my_alpha))
MSE_EN_cv_rnaseq <- matrix(0,1,length(my_alpha))
PR_EN_cv_rnaseq <- matrix(0,1,length(my_alpha))

# Assigning samples to folds, to be used in cross-validation when tunning alpha

set.seed(2010)
my_foldid_rnaseq <- sample(1:10,size=length(ytrain),replace=TRUE)

for (j in 1:length(my_alpha)){
  
  # Logistic model fitting with 10-fold cross-validation for glmnet:
  
  fit_EN_cv_rnaseq <- cv.glmnet(as.matrix(xtrain),as.factor(ytrain),family="binomial",foldid=my_foldid_rnaseq,alpha=my_alpha[j],type.measure="mse")
  
  var_selected_EN_rnaseq <- which(fit_EN_cv_rnaseq$glmnet.fit$beta[,which(fit_EN_cv_rnaseq$cvm == min(fit_EN_cv_rnaseq$cvm))] != 0)
  
  nvar_selected_EN_rnaseq[j] <- length(var_selected_EN_rnaseq)
  
  # Predictions obtained by model i:
  pred_EN_cv_rnaseq[,j] <- predict(fit_EN_cv_rnaseq,as.matrix(xtrain),s=fit_EN_cv_rnaseq$lambda[which(fit_EN_cv_rnaseq$cvm == min(fit_EN_cv_rnaseq$cvm))],type="response")
  
  # Mean squared error of prediction (MSE), area under the precision recall curve:
  MSE_EN_cv_rnaseq[j] <- mean((ytrain-pred_EN_cv_rnaseq[,j])^2)
  PR_EN_cv_rnaseq[j] <- pr.curve(scores.class0=as.vector(round(pred_EN_cv_rnaseq[,j])), weights.class0=ytrain)$auc.integral
}

best_alpha_rnaseq <-  my_alpha[which(MSE_EN_cv_rnaseq== min(MSE_EN_cv_rnaseq))]

# Sparse logistic regression

fit_EN_cv_rnaseq <- cv.glmnet(xtrain, as.factor(ytrain), family="binomial", nfolds=10, alpha=best_alpha_rnaseq,foldid=my_foldid_rnaseq, type.measure="mse")

var_selected_EN_rnaseq_idx <- which(fit_EN_cv_rnaseq$glmnet.fit$beta[,which(fit_EN_cv_rnaseq$cvm == min(fit_EN_cv_rnaseq$cvm))] != 0)

length(var_selected_EN_rnaseq_idx)

# GENES SELECTED___________________________________________________________________

var_selected_EN_rnaseq <- as.matrix(colnames(xtrain[,var_selected_EN_rnaseq_idx]))

# STATISTIC PARAMETERS____________________________________________________________

## misclassifications, MSE ROC and PR

# train set

pred_EN_train_rnaseq <- predict(fit_EN_cv_rnaseq,xtrain,s="lambda.min",type="response")
table(ytrain,round(pred_EN_train_rnaseq))

mse_EN_train_rnaseq <- mean((ytrain - pred_EN_train_rnaseq)^2)

PR_EN_train_rnaseq <- pr.curve(scores.class0=as.vector(round(pred_EN_train_rnaseq)), weights.class0=ytrain)$auc.integral

# test set

pred_EN_test_rnaseq <- predict(fit_EN_cv_rnaseq,xtest,s="lambda.min",type="response")
table(ytest,round(pred_EN_test_rnaseq))

mse_EN_test_rnaseq <- mean((ytest - pred_EN_test_rnaseq)^2)

PR_EN_test_rnaseq <- pr.curve(scores.class0=as.vector(round(pred_EN_test_rnaseq)), weights.class0=ytest)$auc.integral


#################################### IMAGES #################################### 


# Venn diagram_____________________________________________________________________

library("VennDiagram")

venn_diagram <- venn.diagram( x = list( "Twiner" = var_selected_netw_75_code, "EN" = var_selected_75_code, "EN-microarray" = var_selected_EN_array, "EN-RNAseq" = var_selected_EN_rnaseq), imagetype="png", filename = NULL, cat.fontface="bold", cat.cex=1.4, palette="ggplot2", label.col = c("white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "white"), fontface="bold",cex = 1.8, fill = c("#a8db72", "#89bef0",  "#e67027", "#bea4cb")) 
grid.draw(venn_diagram)


# Violin plots____________________________________________________________________ 

col1 = c("#a8db72","#e67027")

all_tumor <- rbind(array_tumor_bl,rnaseq_tumor)
all_normal <- rbind(array_normal_bl,rnaseq_normal)

twiner_gene_expression <- data.frame(c(rep("Tumor",dim(all_tumor)[1]),rep("Normal",dim(all_normal)[1])),rbind(all_tumor[,c(distinct_var_selected_netw_75_code)],all_normal[,c(distinct_var_selected_netw_75_code )]))
names(twiner_gene_expression) <- c("Samples",c(distinct_var_selected_netw_75_code))

library(reshape2)

twiner_violin <- melt(twiner_gene_expression,id.vars="Samples")
names(twiner_violin) <- c("Samples","Gene","Expression")

ggplot(data=twiner_violin, aes(x=Gene,y=Expression,fill=factor(Samples))) +
  geom_violin(aes(fill=Samples, color=Samples)) +
  theme_minimal() +
  scale_colour_manual(values=c("#a8db72","#89bef0"),guide = FALSE) +
  scale_fill_manual(values=c("#a8db72","#89bef0")) +
  guides(fill=guide_legend("Samples")) +
  theme(axis.text.x = element_text(angle = 0)) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14)) +
  theme(legend.position = "right")


# Gene networks_____________________________________________________________________

qgraph_matrix <- data_join_balanced
qgraph_vars <- c(var_selected_netw_75[-which(var_selected_netw_75 %in% var_selected_75)],var_selected_75[-which(var_selected_75 %in% var_selected_netw_75)],var_selected_netw_75[which(var_selected_netw_75 %in% var_selected_75)])

  # ARRAY

qgraph_cor_array_tumor <- cor(qgraph_matrix[technique_label=="A" & sample_label=="1",qgraph_vars], method="pearson")
qgraph_cor_array_normal <- cor(qgraph_matrix[technique_label=="A" & sample_label=="0",qgraph_vars], method="pearson")

qgraph_label <- c(distinct_var_selected_netw_75_code,distinct_var_selected_75_code,common_var_selected_netw_75_code)

qgraph_groups <- list(netEN = c(1:length(distinct_var_selected_netw_75_code)), EN = (length(distinct_var_selected_netw_75_code)+1):(length(distinct_var_selected_netw_75_code)+length(distinct_var_selected_75_code)), common = (length(distinct_var_selected_75_code)+length(distinct_var_selected_netw_75_code)+1):(length(distinct_var_selected_75_code)+length(distinct_var_selected_netw_75_code)+length(common_var_selected_netw_75_code)))

  # RNA-SEQ

qgraph_cor_rnaseq_tumor <- cor(qgraph_matrix[technique_label=="R" & sample_label=="1",qgraph_vars], method="pearson")
qgraph_cor_rnaseq_normal <- cor(qgraph_matrix[technique_label=="R" & sample_label=="0",qgraph_vars], method="pearson")

library(qgraph)

par(mfrow=c(1, 2))
qgraph(qgraph_cor_array_tumor, minimum = 0, vsize = 1.5, groups = qgraph_groups, legend=FALSE,borders = FALSE, node.width=5, color=c( "#e67027","#89bef0","#a8db72" ),labels = qgraph_label, label.cex=2, theme="TeamFortress")
title("Microarray tumor", line = 2.5)
qgraph(qgraph_cor_rnaseq_tumor, minimum = 0, vsize = 1.5, groups = qgraph_groups, legend=FALSE,borders = FALSE, node.width=4, color=c( "#e67027","#89bef0","#a8db72" ),labels = qgraph_label, label.cex=2, theme="TeamFortress")
title("RNA-seq tumor", line = 2.5)


par(mfrow=c(1, 2))
qgraph(qgraph_cor_array_normal, minimum = 0, vsize = 1.5, groups = qgraph_groups, legend=FALSE,borders = FALSE, node.width=4, color=c( "#e67027","#89bef0","#a8db72" ),labels = qgraph_label, label.cex=2, theme="TeamFortress")
title("Microarray normal", line = 2.5)
qgraph(qgraph_cor_rnaseq_normal, minimum = 0, vsize = 1.5, groups = qgraph_groups, legend=FALSE,borders = FALSE, node.width=4, color=c( "#e67027","#89bef0","#a8db72" ),labels = qgraph_label, label.cex=2, theme="TeamFortress")
title("RNA-seq normal", line = 2.5)


# Weights____________________________________________________________________________

array_rnaseq_weights_col <- matrix(1,1,length(array_rnaseq_weights))
array_rnaseq_weights_col[var_selected_75] <- "EN"
array_rnaseq_weights_col[var_selected_netw_75] <- "Twiner"
array_rnaseq_weights_col[var_selected_75[which(var_selected_75 %in% var_selected_netw_75)]] <- "Common" 
array_rnaseq_weights_col[array_rnaseq_weights_col ==1] <-"Non-selected"

array_rnaseq_weights_selected_common_plot <- data.frame(Variable = var_selected_75[which(var_selected_75 %in% var_selected_netw_75)], Weight = array_rnaseq_weights[var_selected_75[which(var_selected_75 %in% var_selected_netw_75)]], label = factor(array_rnaseq_weights_col[var_selected_75[which(var_selected_75 %in% var_selected_netw_75)]]))

array_rnaseq_weights_selected_EN_plot <- data.frame(Variable = var_selected_75, Weight = array_rnaseq_weights[var_selected_75], label = factor(array_rnaseq_weights_col[var_selected_75]))

array_rnaseq_weights_selected_netEN_plot <- data.frame(Variable = var_selected_netw_75, Weight = array_rnaseq_weights[var_selected_netw_75], label = factor(array_rnaseq_weights_col[var_selected_netw_75]))

array_rnaseq_weights_plot <- data.frame(Variable = c(1:length(array_rnaseq_weights)), Weight = array_rnaseq_weights, label = factor(array_rnaseq_weights_col))


ggplot(array_rnaseq_weights_plot, aes(Variable, Weight)) + geom_point(aes(color = label),size=1.2) + 
  scale_colour_manual(values=c("#a8db72","#89bef0","#e6e6e6","#e67027","#a8db72","#a8db72","#89bef0","#e67027","#a8db72"))+
  theme_minimal()+theme(axis.text=element_text(size=16),axis.title=element_text(size=18),legend.title = element_text(size=14, color="white"),legend.text = element_text(size=18))+  geom_point(data=array_rnaseq_weights_selected_common_plot,aes(x=Variable, y=Weight), col="#a8db72",size=1.7) +
  geom_point(data=array_rnaseq_weights_selected_netEN_plot,aes(x=Variable, y=Weight), col="#e67027",size=1.7) +  geom_point(data=array_rnaseq_weights_selected_EN_plot,aes(x=Variable, y=Weight), col="#89bef0",size=1.7) +  geom_point(data=array_rnaseq_weights_selected_common_plot,aes(x=Variable, y=Weight), col="#a8db72",size=1.7) + theme_classic()


  # variables selected by network-based sparse logistic regression only, with weights lower than 0.5, i.e., with more similar correlation pattern within array and rnaseq data sets

as.matrix(colnames(data_join_balanced)[var_selected_netw_75[which(array_rnaseq_weights[var_selected_netw_75] < 0.5)]])

  # Lowest value 

min(array_rnaseq_weights[var_selected_netw_75] )

  # Gene corresponding to the lowest weight value

genes_and_weights <- rbind(colnames(xarray_less),array_rnaseq_weights)
genes_and_weights <- genes_and_weights[,var_selected_netw_75]
genes_and_weights[,which.min(genes_and_weights[2,])] 

################################# SURVIVAL ANALYSIS #################################


# MICROARRAY

er_plus_array_2 <- x_array[which(rownames(x_array) %in% ER_samples),]
data_join <- rbind(xarray_less,xrnaseq_less)

xarray_cox <- er_plus_array_2[,colnames(data_join)]
xarray_cox <- xarray_cox[,var_selected_netw_75] #change for EN or twiner selected genes 

array_survival_status <- SA_days[,2:3]
rownames(array_survival_status) <- SA_days$Patient

xarray_cox_mat <- cbind(array_survival_status,xarray_cox)
array_fit_cox <- coxph(Surv(days_to_death,as.numeric(as.factor(vital_status))) ~ .,data=as.data.frame(xarray_cox_mat))

array_fit_cox_coefficients <- as.matrix(xarray_cox) %*% array_fit_cox$coefficients
array_fit_cox_risk_group <- array_fit_cox_coefficients
array_fit_cox_risk_group[which(array_fit_cox_coefficients > median(array_fit_cox_coefficients))] <- 2
array_fit_cox_risk_group[which(array_fit_cox_coefficients <= median(array_fit_cox_coefficients))] <- 1
array_KM_mat <- as.data.frame(cbind(array_survival_status[,1],as.numeric(as.factor(array_survival_status[,2])),array_fit_cox_risk_group))
names(array_KM_mat) <- c("survival","status","riskgroup")
fit_KM_array <- survfit(Surv(as.numeric(survival),status) ~ riskgroup, data=array_KM_mat, type="kaplan-meier", conf.type="plain")

ggsurvplot(fit_KM_array,legend="bottom", legend.labs = c("low","high"), legend.title = "Risk group", pval=TRUE,palette = c("#a8db72","#f26363"), xlab= "Time (days)")

# RNA-SEQ 

xrnaseq_cox_2<- data_join[,var_selected_netw_75] #change according to the biomarkers we want to test
xrnaseq_cox_2 <- xrnaseq_cox_2[rownames(brca_tumor_erplus_2),]
xrnaseq_cox_mat_2 <- cbind(rnaseq_survival_status,xrnaseq_cox_2)
rnaseq_fit_cox_2 <- coxph(Surv(days_to_death,as.numeric(vital_status)) ~ .,data=as.data.frame(xrnaseq_cox_mat_2))

rnaseq_fit_cox_coefficients_2 <- as.matrix(xrnaseq_cox_2) %*% rnaseq_fit_cox_2$coefficients

rnaseq_fit_cox_risk_group_2 <- rnaseq_fit_cox_coefficients_2
rnaseq_fit_cox_risk_group_2[which(rnaseq_fit_cox_coefficients_2 > median(rnaseq_fit_cox_coefficients_2))] <- 2
rnaseq_fit_cox_risk_group_2[which(rnaseq_fit_cox_coefficients_2 <= median(rnaseq_fit_cox_coefficients_2))] <- 1

rnaseq_KM_mat_2 <- as.data.frame(cbind(rnaseq_survival_status[,1],rnaseq_survival_status[,2],rnaseq_fit_cox_risk_group_2))
names(rnaseq_KM_mat_2) <- c("survival","status","riskgroup")
fit_KM_2 <- survfit(Surv(as.numeric(survival),status) ~ riskgroup, data=rnaseq_KM_mat_2, type="kaplan-meier", conf.type="plain")
ggsurvplot(fit_KM_2,legend="bottom", legend.labs = c("low","high"), legend.title = "Risk group", pval=TRUE,palette = c("#a8db72","#f26363"), xlab= "Time (days)") 

