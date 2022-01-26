library(tidyverse)
library(statTarget)
setwd("~/batchCorrection") # set your working directory

########### QC batch correction #################
## 1) batch correction and imputation using statTarget

## files for statTarget for drift correction
pheno <- "meta4samples.csv"
dfile <- "data4samples.csv"
labels <- read.csv("label4samples.csv", row.names = "SampleName")

## statTarget
statTarget::shiftCor(pheno, dfile,  QCspan = 0.25, Frule = 0.8,
                     degree = 2,imputeM = "KNN", ntree=500, coCV = 50) ## MLmethod = 'QCRLSC' or MLmethod = 'QCRFSC' (default)

## 2) compare before and after batch correction
meta = read.csv(file = "meta4samples.csv", header = TRUE)
QC_index = grep("QC", meta$sample)
meta_sample <- meta[-c(QC_index), ] %>% `row.names<-`(., NULL)

df_after = read.csv(file = "./statTarget/shiftCor/After_shiftCor/shift_sample_cor.csv") %>%
  column_to_rownames(., var = "sample")
#label = as.factor(df_after$class)
df_after = df_after[,-1]

df_before = read.csv(file = "./statTarget/shiftCor/Before_shiftCor/shift_sam_raw.csv")%>%
  column_to_rownames(., var = "sample")
df_before = df_before[,-c(1,2)]

log2_before <- log2(type.convert(df_before)) %>% ## log2 transformation
  t(.) %>% # feature x sample
  as.data.frame(.)
log2_after <- log2(type.convert(df_after)) %>% ## log2 transformation
  t(.) %>% # feature x sample
  as.data.frame(.)

library(proBatch)
color_list = sample_annotation_to_colors(meta, factor_columns = c("batch", "class"), numeric_columns = "order")

#plot_sample_mean(log2_after, meta, sample_id_col = "sample", order_col = "order", batch_col = "batch", color_by_batch = TRUE, ylimits = c(5,15), color_scheme = color_list[["batch"]], plot_title = "After BC") + ylab("Mean Intensity (log2 scale)") #y-range = 10 for comparison purpose
#plot_sample_mean(log2_before, meta, sample_id_col = "sample", order_col = "order", batch_col = "batch", color_by_batch = TRUE, ylimits = c(5,15), color_scheme = color_list[["batch"]], plot_title = "Before BC") + ylab("Mean Intensity (log2 scale)") #y-range = 10 for comparison purpose

log2_after_long = matrix_to_long(log2_after, sample_id_col = "sample")
log2_before_long = matrix_to_long(log2_before, sample_id_col = "sample")

plot_boxplot(log2_after_long, meta, sample_id_col = "sample", batch_col = "batch", color_scheme = color_list[["batch"]], ylimits = c(-20,25), plot_title = "After BC") + ylab("Intensity (log2 scale)") # y-range = 45 for comparison purpose
plot_boxplot(log2_before_long, meta, sample_id_col = "sample", batch_col = "batch", color_scheme = color_list[["batch"]], ylimits = c(-20,25), plot_title = "Before BC") + ylab("Intensity (log2 scale)") # y-range = 45 for comparison purpose


