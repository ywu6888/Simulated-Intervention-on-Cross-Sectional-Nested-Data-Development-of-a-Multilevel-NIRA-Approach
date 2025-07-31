rm(list=ls())
library(devtools)
library(roxygen2)
library(networktools)
library(bootnet)
library(nodeIdentifyR)
library(glmmLasso)
library(qgraph)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)

devtools::install_github("ywu6888/MultilevelNIRA")

library(MultilevelNIRA)
data(MulNIRAdata)

## Fit multilevel Ising network
Mul_Isingmodel <- Multilevel_Isingfit(MulNIRAdata, group.col =  c("city"), rnd = list(city = ~1),
                                      m=100,nlambda = 100) #In our paper, since it is only for demonstration purposes, we set m = 1, nlambda = 5

#################### Multilevel Ising network result output ###########################
vars <- c("NA","UW","WTM","TR","RES","IRR","ASH")
dimnames(Mul_Isingmodel$weiadj) <- list(vars, vars)

graph2 <-qgraph::qgraph(Mul_Isingmodel$weiadj,
                        layout = 'spring',
                        groups = rep(c("GAD"), 7),
                        color  = c("#337BAC"),
                        legend = TRUE)

centralityPlot(graph2,include = c("Strength","Closeness","Betweenness"),scale = c('z-score'))

########## Bootnet stability estimation ###########
# ① Prepare purely numeric node data and grouping vector
node_data <- MulNIRAdata[, setdiff(names(MulNIRAdata), "city")]
cities    <- MulNIRAdata$city
names(cities) <- rownames(MulNIRAdata)
p <- ncol(node_data)  # number of nodes

# ② Define a 'safe' wrapper function using tryCatch to capture fitting errors
myIsingEst_safe <- function(data) {
  df <- data
  df$city <- cities[rownames(df)]
  
  out <- tryCatch({
    res <- Multilevel_Isingfit(
      df,
      group.col = "city",
      rnd       = list(city = ~1),
      gamma     = 0.25,
      m         = 100,
      nlambda   = 100,
      plot      = FALSE
    )#In our paper, since it is only for demonstration purposes, we set m = 1, nlambda = 5
    res$weiadj  
  }, error = function(e) {
    # On any error, return a p×p matrix full of NA
    mat_na <- matrix(
      NA, nrow = p, ncol = p,
      dimnames = list(colnames(data), colnames(data))
    )
    mat_na
  })
  
  out
}

# ③ Run with bootnet
bootres <- bootnet(
  data       = node_data,
  fun        = myIsingEst_safe,
  default    = "none",
  type       = "person",
  nBoots     = 1000,            
  nCores     = 1,
  statistics = c("edge", "strength","Closeness","Betweenness")
)
cs <- corStability(bootres)  
plot(bootres,statistics = c("strength","Closeness","Betweenness"))


##################### Multilevel NIRA application ###################
set.seed(2025)  ###### Set random seed

edgeWeightMatrix <- Mul_Isingmodel$weiadj  # Extract the edge weight matrix

# Extract the fixed thresholds averaged across cities
thresholdVector_fix <- Thresholds_select(Mul_Isingmodel, type = "fixed")

# Get the random thresholds for the 5th city
thresholdVector_random <- Thresholds_select(Mul_Isingmodel,
                                            type  = "random",
                                            group = "city",
                                            level = 5)    # extract logistic regression intercepts as threshold parameters; different from fixed results
thresholdVector_random <- thresholdVector_random + thresholdVector_fix


## NIRA with fixed thresholds
# Alleviating intervention
gs_IsingSamples_fixed <- simulateResponses(edgeWeightMatrix, thresholdVector_fix, "alleviating", 2)
gs_sumIsingSamples_fixed <- calculateSumScores(gs_IsingSamples_fixed)
gs_sumIsingSamplesLong_fixed <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples_fixed)
plotSumScores(sum_scores_long = gs_sumIsingSamplesLong_fixed, perturbation_type = "alleviating")

# Aggravating intervention
gs_IsingSamples_fixed <- simulateResponses(edgeWeightMatrix, thresholdVector_fix, "aggravating", 2)
gs_sumIsingSamples_fixed <- calculateSumScores(gs_IsingSamples_fixed)
gs_sumIsingSamplesLong_fixed <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples_fixed)
plotSumScores(sum_scores_long = gs_sumIsingSamplesLong_fixed, perturbation_type = "aggravating")


## NIRA with random thresholds
# Alleviating intervention
gs_IsingSamples <- simulateResponses(edgeWeightMatrix, thresholdVector_random, "alleviating", 2)
gs_sumIsingSamples <- calculateSumScores(gs_IsingSamples)
gs_sumIsingSamplesLong <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long = gs_sumIsingSamplesLong, perturbation_type = "alleviating")

# Aggravating intervention
gs_IsingSamples <- simulateResponses(edgeWeightMatrix, thresholdVector_random, "aggravating", 2)
gs_sumIsingSamples <- calculateSumScores(gs_IsingSamples)
gs_sumIsingSamplesLong <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long = gs_sumIsingSamplesLong, perturbation_type = "aggravating")

#### 5000 permutation test and plotting
devtools::install_github("kingfly51/NIRA_post")
library(NIRApost)
stat_result <- getstat(gs_sumIsingSamplesLong, method = "bonferroni")
stat_result$stat
plotmeanNIRA(stat_result, "alleviating")
plotmeanNIRA(stat_result, "aggravating")

stat_result_fixed <- getstat(gs_sumIsingSamplesLong_fixed, method = "bonferroni")
stat_result_fixed$stat
plotmeanNIRA(stat_result_fixed, "alleviating")
plotmeanNIRA(stat_result_fixed, "aggravating")
