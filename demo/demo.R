setwd("~/demo/")
library(tidyverse)
#rm(list=ls())

# import data and data wrangling
fh = read.csv(file = "malaria_feature_table.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE) %>%
  as.data.frame()
# A peak intensity table with feature on row, sample in columns.

fh <- fh[,order(fh["Label",])] ## order samples by group
labels_d1 <- as.matrix(fh[1,])

#-------------------------------------------------------------------------------
### check missing data
d1 <- fh[-1,]
sum(d1==0) ## Total missing values

d1[d1 == 0] <- NA ## Change zeros to NA
sum(is.na(d1))

# d1 is a character matrix, turn it to a numeric matrix
#d1 = as.matrix(apply(d1, 2, as.numeric)) %>%
#  `rownames<-`(., rownames(fh[-1,]))
#-------------------------------------------------------------------------------
### Assess feature presence/absence
d1$nacount <- as.matrix(apply(d1, 1, function(x) sum(is.na(x)))) #count the missing signal values by feature
d1$nacount_Native = as.matrix(apply(d1[, 1:6], 1, function(x) sum(is.na(x))))
d1$nacount_Semi = as.matrix(apply(d1[, 7:12], 1, function(x) sum(is.na(x))))

# keep a feature when it has na less than 50% of all samples.
d2 <- d1[which(d1$nacount < round((length(d1)-3)*.5)), ] %>%
        select(-c(nacount, nacount_Semi, nacount_Native))
sum(is.na(d2))
#[1] 3001
sum(is.na(d2))/(dim(d2)[1]*dim(d2)[2]) # percentage missing data
#[1] 0.0582809

# keep a feature when it has na less than 50% of CNT and PNT
d2_2 = d1[which(d1$nacount_Native < round(6*.5) & d1$nacount_Semi < round(6*.5)), ] %>%
        select(-c(nacount, nacount_Semi, nacount_Native))
sum(is.na(d2_2))
#[1] 1111
sum(is.na(d2_2))/(dim(d2_2)[1]*dim(d2_2)[2]) # percentage missing data
#[1] 0.02419847

#-------------------------------------------------------------------------------
### missing data imputation
# decide K (how many neighbor to use)
k <- d2_2 %>% length() %>% sqrt() # square root of number of features
if (k>=3){
        if (is.integer(k)) k = k else k = ceiling(k)
} else {k = 3}

imputed_meta_d1 <- impute::impute.knn(as.matrix(d2_2), k = k, rowmax = 0.5,
                              colmax = 0.8, maxp = 1500)

imputed_d1 <- as.data.frame(imputed_meta_d1[["data"]], colnames = TRUE)
MASS::truehist(unlist(imputed_d1), main = "Imputed", nbins = 50) # check data distribution
sum(imputed_d1 == 0)
sum(imputed_d1 < 0)
#-------------------------------------------------------------------------------
### log2 transformation
log2_d1 <- log2(type.convert(imputed_d1))
MASS::truehist(unlist(log2_d1), main = "log2 transformed", nbins = 50) # check data distribution
car::qqPlot(unlist(log2_d1), main = "log2 transformed",xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")

#-------------------------------------------------------------------------------
### normalization and scaling steps are optional depends on the distribution of your data

# Implement invariant normalization using  package lumi
# Implement cubic spline normalization using package affy
# Implement EigenMS using package ProteoMM

# EigenMS
library(ProteoMM)
d4 <- as.data.frame(log2_d1) #feature x sample
d4 <- cbind.data.frame(row.names(d4), d4)
rowmeta <- d4[1,]
rowmeta2 <- as.data.frame(t(colnames(d4)))
colmeta <- d4[-1,1]

intsCols <- 2:dim(d4)[2]
m_logInts = ProteoMM::make_intencities(d4, intsCols)

metaCols <- 1
m_mets.info = as.matrix(ProteoMM::make_meta(d4, metaCols))

grps = as.factor(t(labels_d1))
m_nummiss = sum(is.na(m_logInts)) #check zeros

set.seed(42)
hs_m_ints_eig1 = eig_norm1(m=m_logInts,treatment = grps,
                           prot.info=m_mets.info)
hs_m_ints_norm_1bt = eig_norm2(rv=hs_m_ints_eig1)

norm_d1 <- as.data.frame(hs_m_ints_norm_1bt[["norm_m"]])
par(mfrow=c(1,1))

MASS::truehist(unlist(t(norm_d1)), main = "EigenMS normalization", nbins = 50)
car::qqPlot(unlist(t(norm_d1)),main = "log2-EigenMS",xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")

# Implement pareto scaling using package MetabolAnalyze

#-------------------------------------------------------------------------------
### plot sample mean and boxplot
library(proBatch)

Label = as.matrix(fh[1,]) %>%
  t(.) %>%
  as.data.frame(.) %>%
  `row.names<-`(., colnames(fh)) %>%
  `colnames<-`(., "Label") %>%
  rownames_to_column(., var = "sample")

color_list = sample_annotation_to_colors(Label,sample_id_col = "sample", factor_columns = c("Label"), numeric_columns = NULL)

log2_d1_long = matrix_to_long(log2_d1, sample_id_col = "sample")
plot_boxplot(log2_d1_long, Label, sample_id_col = "sample", batch_col = "Label", color_scheme = color_list[["Label"]], ylimits = c(-10,55)) + ylab("Intensity (log2 scale)") # y-range = 45 for comparison purpose

#after normalization
norm_d1_long = matrix_to_long(norm_d1, sample_id_col = "sample")
plot_boxplot(norm_d1_long, Label, sample_id_col = "sample", batch_col = "Label", color_scheme = color_list[["Label"]], ylimits = c(-10,55)) + ylab("Intensity (after normalization)") # y-range = 45 for comparison purpose


#-------------------------------------------------------------------------------
### Univariate analysis
d5 <- rbind.data.frame(labels_d1, norm_d1)

C <- d5[-1,which(d5[1,]=="Naive")]; P <- d5[-1,which(d5[1,]=="Semi_immue")]

wt <- data.frame(Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double())

for (i in 1:dim(norm_d1)[1]){
        wt[i,1] <- t.test(type.convert(t(P[i,])),type.convert(t(C[i,])))$p.value
        wt[i,3] <- wilcox.test(type.convert(t(P[i,])),type.convert(t(C[i,])))$p.value
        az1 <- 2^(type.convert(C[i,])); az1L <- length(az1)
        az2 <- 2^(type.convert(P[i,])); az2L <- length(az2)
        wt[i,5] <- gtools::foldchange(sum(az2)/az2L, sum(az1)/az1L)
        #print(az1)
}
wt[,2] <- as.data.frame(p.adjust(wt[,1], method = 'BH'))#, n = length(p))
wt[,4] <- as.data.frame(p.adjust(wt[,3], method = 'BH'))#, n = length(p))

colnames(wt) <- c("pT", "BHT", "pW", "BHW", "FC(lin)")
row.names(wt) <- rownames(norm_d1)

# Calculate log2FoldChange
log2FoldChange =gtools::foldchange2logratio(wt$`FC(lin)`,base = 2)
wt = cbind(wt,log2FoldChange)

sum(wt[,1] < 0.05 & wt[,3] < 0.05);sum(wt[,1] < 0.05);sum(wt[,3] < 0.05)
#[1] 2895
#[1] 2971
#[1] 2927
sum(wt[,2] < 0.01 & wt[,4] < 0.01);sum(wt[,2] < 0.01);sum(wt[,4] < 0.04)
#[1] 2448
#[1] 2521
#[1] 2830
sig_wt <- wt[which(wt[,2] < 0.01 & wt[,4] < 0.01),]

#-------------------------------------------------------------------------------
### Multivariate analysis, PCA
library(pcaMethods)
pc = pcaMethods::pca(t(norm_d1), nPcs = 9, method = "svd", scale = "pareto", center = TRUE)
#pc_prcomp = prcomp(t(norm_d1), scale. = TRUE)

pca = rbind.data.frame(labels_d1, norm_d1) %>%
        t(.) %>%
        as.data.frame(.) %>%
        rownames_to_column(.) %>%
        as_tibble(.) %>%
        select(rowname, Label) %>%
        column_to_rownames(.,var = "rowname") %>%
        merge(., pcaMethods::scores(pc), by=0)
        #merge(., pc$x, by = 0)

# scree plot

var_explained_pc = data.frame(PC = paste0("PC", 1:9), var_explained = (pc@sDev)^2/sum((pc@sDev)^2), eigenvalues = pc@sDev)
var_explained_pc %>%
        ggplot(aes(x=PC,y=eigenvalues, group=1))+
        geom_point(size=4)+
        geom_line()+
        labs(title="Scree plot: PCA on scaled data")

# score plot
p_pca = ggplot(pca, aes(PC1, PC2,colour=Label)) + geom_point() + stat_ellipse() +
        xlab(paste("PC1", round(pc@R2[1] * 100, digits = 1), "% of the variance")) +
        ylab(paste("PC2", round(pc@R2[2] * 100, digits = 1), "% of the variance")) #+ ggtitle(label = "log2-Invariant-Pareto")
p_pca

pc_loadings = as.data.frame(pc@loadings)
par(mar=c(8, 3, 3, 1))
barplot( sort( pc@loadings[,1][which(abs(pc@loadings[,1]) > 0.035)], decreasing = FALSE), main = "PC1 loading plot", las = 2, cex.names=0.6)

### Multivariate analysis, PLS-DA
library(ropls)
# a == norm_d1; df == d5 (d5 is norm_d1 plus label)
classes <- as.character(labels_d1) #phenotype groups

norm_d1[,] <- as.data.frame(apply(norm_d1[,], 2,
                            function(x) as.double(as.character(x))))

my.plsda <- ropls::opls(t(norm_d1),classes)

# plots to export
#plot(my.plsda,
#     typeVc = "x-score",
#     parAsColFcVn = classes, parLabVc = rep('o', length(classes)))

my.vip <- as.data.frame(sort(ropls::getVipVn(my.plsda), decreasing = T))
colnames(my.vip)="VIP"


