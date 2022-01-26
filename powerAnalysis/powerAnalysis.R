setwd("~/power_analysis")
library(tidyverse)


# import test data
fh = read.csv(file = "human_cachexia.csv", header=TRUE) %>%
        as.data.frame()%>%
        column_to_rownames(., var = "Patient.ID")

# check missing values
d1 = fh[,-1]
sum(d1==0) ## Total missing values

# check distribution
MASS::truehist(unlist(d1), main = "raw", nbins = 50) # not normal

# log2 transformation
log2_d1 <- log2(type.convert(d1))
MASS::truehist(unlist(log2_d1), main = "log2 transformed", nbins = 50) # near normal
car::qqPlot(unlist(log2_d1), main = "log2 transformed",xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")

#-------------------------------------------------------------------------------
# Power Analysis
# this code is adapted from R package SSPA power analysis functions.
library(SSPA)
browseVignettes("SSPA") # read more about this package and how it handle power analysis
comparison_list = c("control vs. cachexic") # in the case there are more than two groups, extend this list.
class = as.factor(fh$Muscle.loss)
table(class)


## Step 1 ----
pdD = vector(mode = "list", length = length(comparison_list)) # list of pilotData
for (i in 1:length(comparison_list)) {
  grp.nms <- strsplit(comparison_list, " vs. ")[[i]]
  inx1 = which(class==grp.nms[1]); inx2 = which(class==grp.nms[2]);
  n1 <- length(inx1); n2 <- length(inx2);

  stats <- apply(log2_d1, 2, function(x) {
    tmp <- try(t.test(x[inx1], x[inx2], paired = FALSE, var.equal = T));
    if(class(tmp) == "try-error") {
      return(NA);
    }else{
      return(tmp$statistic);
    }
  })
  stats <- stats[!is.na(stats)];

  # perform power analysis using SSPA package
  pdD[[i]] <- SSPA::pilotData(statistics = stats,
                         samplesize = sqrt(n1+n2),
                         distribution="t",
                         df=n1+n2-2);
}
# create power statistics plot
SSPA::plot(pdD[[1]]);
#SSPA::plot(pdD[[2]]);
#SSPA::plot(pdD[[3]]);

## Step 2 ----
ssD = vector(mode = "list", length = length(comparison_list)) # list of sampleSize
for (i in 1:length(comparison_list)) {
  res <- round(length(pdD[[i]]@statistics)/2);
  ssD[[i]] <- SSPA::sampleSize(pdD[[i]], method="congrad", control=list(from=-6, to=6, resolution=res)); # class SampleSize
}
# create effect size plot
SSPA::plot(ssD[[1]])
#SSPA::plot(ssD[[2]])
#SSPA::plot(ssD[[3]])

## Step 3 ----
fdr.lvl = 0.1; # set fdr
smplSize = 100; # set predicted power sample size
Jpred <- GetSampleSizeLadder(smplSize);

pwrD = vector(mode = "list", length = length(comparison_list)) # list of predictedPower
for (i in 1:length(comparison_list)) {
  pi0 <- ssD[[i]]@pi0;

  if(fdr.lvl >= pi0){
    fdr.lvl <- signif(pi0-pi0/10, 3);
  }
  pwrD[[i]] <- SSPA::predictpower(ssD[[i]], samplesizes= sqrt(Jpred/2), alpha=fdr.lvl)
}

# create predicted power plot and table
pwrTable = vector(mode = "list", length = length(comparison_list)) # list of result table
pwrTable[[1]]
#pwrTable[[2]]
#pwrTable[[3]]

for (i in 1:length(comparison_list)) {
  plot(Jpred, pwrD[[i]], type="n", ylim=c(0,1), ylab="Predicted power", xlab="Sample Size (per group)", main= comparison_list[i]);
  grid(col = "lightgray", lty = "dotted", lwd = 1);
  lines(Jpred, pwrD[[i]], lwd=4, col="orange");
  points(Jpred, pwrD[[i]], pch=17);
  pwrTable[[i]] = data.frame(sampleSize = Jpred, predictedPower = pwrD[[i]])
}


#-------------------------------------------------------------------------------
GetSampleSizeLadder <- function(maxNum){
  Jpred <- c(3, 6, 10, 16, 24, 40, 60, 100, 150, seq(200, 1000, 100));
  inx <- which(Jpred == min(Jpred[Jpred>=maxNum]))
  return(Jpred[1:inx]);
}
