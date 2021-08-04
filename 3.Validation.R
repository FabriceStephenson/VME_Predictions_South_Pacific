##==========================      3.Validation.R                             =====================================##

##--------------------------  Authors: Fabrice Stephenson & Owen Anderson
##--------------------------  Aim: code for assessing model fit metrics of ensemble models fitted using training 
##--------------------------  and independent evaluation data for MS "Presence-only habitat suitability models for  
##--------------------------  Vulnerable Marine Ecosystem indicator taxa in the South Pacific have reached their 
##--------------------------  predictive limit"
##--------------------------  Project : SPR201801 - Variation - SC 2019-2020
##--------------------------  Start date : 01/05/2020
##--------------------------  End date : 03/08/2021

##=============================================================================================================##
 
####=======================         PACKAGES                 =============================================####
# Load required packages
library("raster")library("pROC")

# Load functions
source("./1.Functions.R")

####=======================         LOAD FILES               =============================================####
# PROJECTION
TMERCproj<-CRS("+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs")

# Sp abbreviations
sp.names <- c("COB", "COR","DEM","ERO","GDU","HEX","MOC","PTU","SOC","SVA","Alcy")

setwd("~")

# LOAD BIOLOGICAL DATA 

# A dataframe with 15,557 rows (species records) and 34 columns (pa: 0 or 1 presence/absence;
# SpeciesName: factor with species names; X Y coordinates; FID: unique spatial identifier; 
# Update: data used in Georgian et al. 2019or evaluation data 
# and env. variables: "appox","arag","bathy","calcite","cbpm_me","cbpm_sd","dissox","sil",
# "epp_me","epp_sd","final_bathy_bpi_broad_25_125","final_bathy_bpi_fine_5_25",
# "final_bathy_slope", "final_bathy_slope_std05","final_bathy_vrm03","idwgrav",
# "idwmud","nitrate","perox","phos","poc","salinity","seamounts","sigma","temp","vgpm_sd"
load("PA.Rdata")
# Split data into training data (Georgian et al., 2019 data) and independent evaluation data
PA_T <- PA[PA$Update == "Georgian",]
PA_E <- PA[!PA$Update == "Georgian",]

####=======================    VALIDATION               =============================================####
####=======================    TRAINING DATA            =============================================####
PA <- PA_T #validation with model training data

# Loop through folders and load rasters (BRT, RF, MX, ENS - seperately)
# Test AUC, TSS, Sens, Spec for each model type

for (i  in 1:length(sp.names)){
  # ALL PRESENCES FOR TAXON I
  species <- as.character(sp.names[i])
  DF <- PA[PA$SpeciesName == species,] # presence
  
  # ALL TARGET GROUP BACKGROUND DATA FOR TAXON I
  absence <- droplevels(PA[-which(PA$SpeciesName == species), ] ) 
  absence <- absence[!absence$FID %in% DF$FID, ]
  absence$pa <- 0
  
  # IMPORT SPATIAL PREDICTIONS
  setwd(paste("~", "/TRAINING/", species, sep =""))
  
  mean.RF <- raster(paste(species,"_RF_mean.tif", sep = ""))
  mean.BRT <- raster(paste(species,"_BRT_mean.tif", sep = ""))
  mean.MX <- raster(paste(species,"_MX_mean.tif", sep = ""))
  Ensemble <- raster(paste(species,"_Ensemble.tif",sep=""))
  
  # BOOTSTRAP METRIC DATA
  rnd.boots <- 50
  Boot_RF_BRT_EN <- data.frame(RF_AUC = double(rnd.boots), RF_TSS = double(rnd.boots),RF_Sen = double(rnd.boots),RF_Spec = double(rnd.boots),
                               BRT_AUC = double(rnd.boots), BRT_TSS = double(rnd.boots),BRT_Sen = double(rnd.boots),BRT_Spec = double(rnd.boots),
                               MX_AUC = double(rnd.boots), MX_TSS = double(rnd.boots),MX_Sen = double(rnd.boots),MX_Spec = double(rnd.boots),
                               ENS_AUC = double(rnd.boots), ENS_TSS = double(rnd.boots), ENS_Sen = double(rnd.boots),ENS_Spec = double(rnd.boots))
  
  for (n in 1:rnd.boots){
    rnd.abs <- absence_unique[sample(nrow(absence_unique),nrow(DF),replace = F),]
    DF.PA <- as.data.frame(rbind(DF, rnd.abs))
    DF.PA$pa <- as.factor(DF.PA$pa)
    # RF
    preds <- raster::extract(mean.RF,DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,1] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,2] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,3] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,4] <- Sens_spec[1,2] # True negatives
    # BRT
    DF.PA <- as.data.frame(rbind(DF, rnd.abs))
    preds <- raster::extract(mean.BRT,DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,5] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,6] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,7] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,8] <- Sens_spec[1,2] # True negatives
    # Maxent
    preds <- raster::extract(mean.MX, DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,9] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,10] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,11] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,12] <- Sens_spec[1,2] # True negatives
    # Ensemble
    preds <- raster::extract(Ensemble, DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,13] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,14] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,15] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,16] <- Sens_spec[1,2] # True negatives
  }
  
  deviance_mat_EN <- apply(Boot_RF_BRT_EN, 2, function(x) c(mean = mean(x), sd = sd(x)))
  
  if (i == 1){
    Model_fit <- deviance_mat_EN
  } else {
    Model_fit <- rbind(Model_fit, deviance_mat_EN)
  }
 
  if (i == length(sp.names)){
    setwd("~")
    write.csv(Model_fit, file="Model_fit_TRAINING.csv")
    print(paste("Iteration ",i, "  ; Taxon: ", species, sep = "" ))
  } else{
    print(paste("Iteration ",i, "  ; Taxon: ", species, sep = "" ))
  }
}

####=======================    EVALUATION DATA            =============================================####
PA <- PA_E # evaluation DATA

for (i  in 1:length(sp.names)){
  # ALL PRESENCES FOR TAXON I
  species <- as.character(sp.names[i])
  DF <- PA[PA$SpeciesName == species,] # presence
  
  # ALL TARGET GROUP BACKGROUND DATA FOR TAXON I
  absence <- droplevels(PA[-which(PA$SpeciesName == species), ] ) 
  absence <- absence[!absence$FID %in% DF$FID, ]
  absence$pa <- 0
  
  # IMPORT SPATIAL PREDICTIONS
  setwd(paste("~", "/TRAINING/", species, sep =""))
  
  mean.RF <- raster(paste(species,"_RF_mean.tif", sep = ""))
  mean.BRT <- raster(paste(species,"_BRT_mean.tif", sep = ""))
  mean.MX <- raster(paste(species,"_MX_mean.tif", sep = ""))
  Ensemble <- raster(paste(species,"_Ensemble.tif",sep=""))
  
  # BOOTSTRAP METRIC DATA
  rnd.boots <- 50
  Boot_RF_BRT_EN <- data.frame(RF_AUC = double(rnd.boots), RF_TSS = double(rnd.boots),RF_Sen = double(rnd.boots),RF_Spec = double(rnd.boots),
                               BRT_AUC = double(rnd.boots), BRT_TSS = double(rnd.boots),BRT_Sen = double(rnd.boots),BRT_Spec = double(rnd.boots),
                               MX_AUC = double(rnd.boots), MX_TSS = double(rnd.boots),MX_Sen = double(rnd.boots),MX_Spec = double(rnd.boots),
                               ENS_AUC = double(rnd.boots), ENS_TSS = double(rnd.boots), ENS_Sen = double(rnd.boots),ENS_Spec = double(rnd.boots))
  
  for (n in 1:rnd.boots){
    rnd.abs <- absence_unique[sample(nrow(absence_unique),nrow(DF),replace = F),]
    DF.PA <- as.data.frame(rbind(DF, rnd.abs))
    DF.PA$pa <- as.factor(DF.PA$pa)
    # RF
    preds <- raster::extract(mean.RF,DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,1] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,2] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,3] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,4] <- Sens_spec[1,2] # True negatives
    # BRT
    DF.PA <- as.data.frame(rbind(DF, rnd.abs))
    preds <- raster::extract(mean.BRT,DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,5] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,6] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,7] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,8] <- Sens_spec[1,2] # True negatives
    # Maxent
    preds <- raster::extract(mean.MX, DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,9] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,10] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,11] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,12] <- Sens_spec[1,2] # True negatives
    # Ensemble
    preds <- raster::extract(Ensemble, DF.PA[, c("X","Y")])
    myRoc <- pROC::roc(DF.PA$pa, preds, quiet = T)
    Boot_RF_BRT_EN[n,13] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    Boot_RF_BRT_EN[n,14] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    Boot_RF_BRT_EN[n,15] <- Sens_spec[1,3] # True positives
    Boot_RF_BRT_EN[n,16] <- Sens_spec[1,2] # True negatives
  }
  
  deviance_mat_EN <- apply(Boot_RF_BRT_EN, 2, function(x) c(mean = mean(x), sd = sd(x)))
  
  if (i == 1){
    Model_fit <- deviance_mat_EN
  } else {
    Model_fit <- rbind(Model_fit, deviance_mat_EN)
  }
  
  if (i == length(sp.names)){
    setwd("~")
    write.csv(Model_fit, file="Model_fit_EVALUTAION.csv")
    print(paste("Iteration ",i, "  ; Taxon: ", species, sep = "" ))
  } else{
    print(paste("Iteration ",i, "  ; Taxon: ", species, sep = "" ))
  }
}