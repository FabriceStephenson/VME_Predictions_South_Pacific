##==========================      2.TrainingModels.R                             =====================================##

##--------------------------  Authors: Fabrice Stephenson & Owen Anderson
##--------------------------  Aim: code for ensemble modelling (Random Forest; Boosted Regression Tree and Maxent) using 
##--------------------------  training data MS "Presence-only habitat suitability models for Vulnerable Marine Ecosystem 
##--------------------------  indicator taxa in the South Pacific have reached their predictive limit"
##--------------------------  Project : SPR201801 - Variation - SC 2019-2020
##--------------------------  Start date : 01/05/2020
##--------------------------  End date : 03/08/2021

##=============================================================================================================##

####=======================         PACKAGES                 =============================================####
# Load required packages
library("raster");library("rJava")

# increase the allocation of memory for Maxent modelling - has to be done before loading DISMO
options(java.parameters = c('-Xss2560k', '-Xmx2g')) 

library("dismo");library("gbm");library("pROC")

# Packages for running bootstraps in parallel 
library(parallel);library(foreach);library(bigstatsr)

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

# LOAD ENVIRONMENTAL DATA
load("predStack1km.Rdata") 
# A raster stack of environmental variables available: "appox","arag","bathy","calcite","cbpm_me","cbpm_sd","dissox","sil",
# "epp_me","epp_sd","final_bathy_bpi_broad_25_125","final_bathy_bpi_fine_5_25",
# "final_bathy_slope", "final_bathy_slope_std05","final_bathy_vrm03","idwgrav",
# "idwmud","nitrate","perox","phos","poc","salinity","seamounts","sigma","temp","vgpm_sd"

# Environmental variables also available as a dataframe
Pred_1km <- na.omit(as.data.frame(predStack1km, xy = T))

# ENVIRONMENTAL VARIABLES FOR EACH TAXON BASED ON GEORGIAN ET AL 2019
imp.var.COB <- c("calcite","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.COR <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.DEM <- c("sil","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.ERO <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.GDU <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.HEX <- c("sil","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.MOC <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.PTU <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.SOC <- c("calcite","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")
imp.var.SVA <- c("arag","bathy","final_bathy_bpi_broad_25_125","dissox","poc","seamounts",
                 "temp","idwgrav","idwmud", "final_bathy_slope_std05","final_bathy_vrm03")

imp.var.tune <- list() # list of env preds that can be called using sp.names
for (i in 1:length(sp.names)){
  imp.var.tune[[i]] <- get(paste("imp.var.", sp.names[i], sep  = ""))
  names(imp.var.tune)[i] <- sp.names[i]
  rm(list = c(paste("imp.var.", sp.names[i], sep  = ""))) # tidy workspace
}

# learning rate for models = > 1500 trees
lr <- c(0.04, 0.025, 0.075, 0.004, 0.01, 0.1, 0.004, 0.01, 0.075,0.001,0.001)

####=======================    BOOTSTRAP: TRAINING DATA       =============================================####
n.boot <- 200
sp.length <- length(sp.names)
length.map <- nrow(Pred_1km)
PA <- PA_T

# SETUP PARRALLELISATION 
# check number of cores 
detectCores()
java <- c("rJava")

packages <- c("dismo", "gbm","pROC", "extendedForest")

# record time at start
start <- Sys.time()
# MAIN MODEL LOOP (looping through species and bootstraps)
for (i in 1:sp.length){# i = number of species
  # ALL PRESENCES FOR TAXON I
  species <- as.character(sp.names[i])
  DF <- PA[PA$SpeciesName == species,] # presence
  
  # ALL TARGET GROUP BACKGROUND DATA FOR TAXON I
  absence <- droplevels(PA[-which(PA$SpeciesName == species), ] ) 
  absence <- absence[!absence$FID %in% DF$FID, ]
  absence$pa <- 0
  
  # START PARALLEL LOOP
  cl <- parallel::makeCluster(7) # only use maximum number of cores -1 for parallelisation
  doParallel::registerDoParallel(cl)
  
  tmp3 <- foreach(k = 1:n.boot) %dopar% { 
    # LOAD PACKAGES IN FOREACH LOOP
    lapply(java, require, character.only = TRUE) # load java before DISMO
    options(java.parameters = c('-Xss2560k', '-Xmx2g')) # settings for JAVA
    lapply(packages, require, character.only = TRUE) # load other packages
    
    # RANDOM SAMPLE FOR BOOTSTRAP
    rnd.train_P <- DF[sample(nrow(DF),nrow(DF),replace = T),]
    rnd.train_A <- absence[sample(nrow(absence),nrow(DF),replace = T),]
    train_PA <- rbind(rnd.train_P,rnd.train_A)
    train_PA$pa <- as.factor(train_PA$pa)
    
    # STORAGE FOR MODEL RESULTS
    result <- multiResultClass()
    
    # RF MODEL 
    M2 <- tuneRF(x = train_PA[,imp.var.tune[[i]]],
                 y = train_PA$pa,
                 mtryStart = 3,
                 ntreeTry = 1000,
                 stepFactor = 2,
                 improve = 0.05,
                 trace = F, plot = F, doBest = T)
    
    # MODEL GOODNESS-OF-FIT METRICS
    result$SpeR2.RF <- Eval_Metrics(M2, "RF", train_PA, eval_PA)
    
    # VARIABLE CONTRIBUTION
    result$ImpVar.RF <- Var_Imp(M2, "RF")
    
    # SPATIAL PREDICTION
    result$SpPred.RF <- round(predict(M2, Pred_1km, type = "prob")[,2], digits = 2)
    
    # BRT MODEL 
    train_PA <- rbind(rnd.train_P, rnd.train_A)
    
    M1<-gbm.step(data= train_PA,
                 gbm.x = imp.var.tune[[i]],
                 gbm.y = 1,
                 family = "bernoulli", 
                 tree.complexity = 3,
                 learning.rate = lr[i],
                 bag.fraction = 0.6, n.folds = 5,
                 max.trees = 10000, step.size = 25, plot.main = F,
                 tolerance.method = "auto",
                 tolerance = 0.001, silent = T)
    
    # MODEL GOODNESS-OF-FIT METRICS
    result$SpeR2.BRT <- Eval_Metrics(M1, "BRT", train_PA, eval_PA)
    
    # VARIABLE CONTRIBUTION
    result$ImpVar.BRT <- Var_Imp(M1, "BRT")
    
    # SPATIAL PREDICTION
    result$SpPred.BRT <- round(predict.gbm(M1, Pred_1km, n.trees = M1$gbm.call$best.trees, 
                                           family = "bernoulli", type = "response"),digits = 2)
    
    # MAXENT MODEL
    # USE ALL TARGET GROUP BACKGROUND DATA FOR MAXENT
    rnd.train_A <- absence[sample(nrow(absence),nrow(absence),replace = T),]
    train_PA <- rbind(rnd.train_P, rnd.train_A)
    train_PA <- na.omit(train_PA)
    pa <- as.vector(train_PA[,1])
    Env <- train_PA[,imp.var.tune[[i]]]
    
    M3 <- dismo::maxent(x=Env, ## env conditions
                        p=pa,   ## 1:presence or 0:absence
                        maximumiterations = 5000,
                        jacknife = T)
    
    # MODEL GOODNESS-OF-FIT METRICS
    result$SpeR2.MX <- Eval_Metrics(M3, "Maxent", train_PA, eval_PA)
    
    # VARIABLE CONTRIBUTION
    result$ImpVar.MX <- Var_Imp(M3, "Maxent")
    
    # SPATIAL PREDICTION
    result$SpPred.MX <- round(predict(M3, Pred_1km),digits = 2)
    
    return(result)
  }
  parallel::stopCluster(cl)
  print(paste("Cluster stopped. Taxon ", species))
  # insert serial backend, otherwise can result in errors in repetitive tasks
  registerDoSEQ()
  
  #### BOOTSTRAP RESULTS ####
  # CREATE FOLDERS FOR SAVING OUTPUTS
  dir.create(paste("~/TRAINING/", species, sep =""))
  setwd(paste("~/TRAINING/", species, sep =""))
  
  # RANDOM FOREST MODEL SUMMARY - CALCULATE MEAN VALUES FROM BOOTSTRAPS
  # Importance of Environmental predictors
  boot_array_ImpVar <- as.data.frame(array(0, dim = c(length(imp.var), n.boot)))
  row.names(boot_array_ImpVar) <- imp.var
  for (k in 1:n.boot){boot_array_ImpVar[rownames(tmp3[[k]]$ImpVar.RF),k] <- tmp3[[k]]$ImpVar.RF}
  influence_mat_RF <- t(apply(boot_array_ImpVar, 1, function(x) c(RF.mean = mean(x), RF.sd = sd(x))))
  # write.csv(influence_mat_RF, file=paste("preds.influences_",species, ".csv", sep = ""))
  
  # Model fit metrics
  boot_array_dev <- as.data.frame(array(0, dim = c(6, n.boot)))
  for (L in 1:n.boot){boot_array_dev[,L] <- tmp3[[L]]$SpeR2.RF}
  deviance_mat_RF <-t(apply(boot_array_dev, 1, function(x) c(RF.mean = mean(x), RF.sd = sd(x))))
  rownames(deviance_mat_RF) <- row.names(tmp3[[1]]$SpeR2.RF)
  # write.csv(deviance_mat_RF, file=paste("mean.model.fit_",species, ".csv", sep = ""))
  
  # Mean spatial predictions 
  boot_mat <- array(0, dim = c(length.map, n.boot))
  for (M in 1:n.boot){boot_mat[,M] <- tmp3[[M]]$SpPred.RF}
  boot.mean.RF <-apply(boot_mat,1,mean)
  boot.sd.RF <-apply(boot_mat,1,sd)
  
  # RASTER OF MEAN SPATIAL PREDICTIONS
  mean.RF <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                      y =Pred_1km[,2], 
                                      z =boot.mean.RF),
                           crs = TMERCproj) # same proj as orginal env variables tiff files
  
  # RASTER OF SD OF THE MEAN SPATIAL PREDICTIONS
  UC.RF <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                    y =Pred_1km[,2], 
                                    z =boot.sd.RF),
                         crs = TMERCproj)
  
  # Save rasters to folder
  writeRaster(mean.RF, filename = paste(species,"_RF_mean.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  writeRaster(UC.RF, filename = paste(species,"_RF_UC.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  
  # BOOSTED REGRESSION TREE MODEL SUMMARY - CALCULATE MEAN VALUES FROM BOOTSTRAPS
  # Importance of Environmental predictors
  boot_array_ImpVar <- as.data.frame(array(0, dim = c(length(imp.var), n.boot)))
  row.names(boot_array_ImpVar) <- imp.var
  for (k in 1:n.boot){boot_array_ImpVar[rownames(tmp3[[k]]$ImpVar.BRT),k] <- tmp3[[k]]$ImpVar.BRT}
  influence_mat_BRT <- t(apply(boot_array_ImpVar, 1, function(x) c(BRT.mean = mean(x), BRT.sd = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds
  # write.csv(influence_mat_BRT, file=paste("preds.influences_",species, ".csv", sep = ""))
  
  # Model fit metrics
  boot_array_dev <- as.data.frame(array(0, dim = c(6, n.boot)))
  for (L in 1:n.boot){boot_array_dev[,L] <- tmp3[[L]]$SpeR2.BRT}
  deviance_mat_BRT <-t(apply(boot_array_dev, 1, function(x) c(BRT.mean = mean(x), BRT.sd = sd(x))))
  rownames(deviance_mat_BRT) <- row.names(tmp3[[1]]$SpeR2.BRT)
  # write.csv(deviance_mat_BRT, file=paste("mean.model.fit_",species, ".csv", sep = ""))
  
  # Mean spatial predictions 
  boot_mat <- array(0, dim = c(length.map, n.boot))
  for (M in 1:n.boot){boot_mat[,M] <- tmp3[[M]]$SpPred.BRT}
  boot.mean.BRT <-apply(boot_mat,1,mean)
  boot.sd.BRT <-apply(boot_mat,1,sd)
  
  # RASTER OF MEAN SPATIAL PREDICTIONS
  mean.BRT <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                       y =Pred_1km[,2], 
                                       z =boot.mean.BRT),
                            crs = TMERCproj) # same proj as orginal env variables tiff files
  # RASTER OF SD OF THE MEAN SPATIAL PREDICTIONS
  UC.BRT <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                     y =Pred_1km[,2], 
                                     z =boot.sd.BRT),
                          crs = TMERCproj)
  
  # save rasters to folder
  writeRaster(mean.BRT, filename = paste(species,"_BRT_mean.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  writeRaster(UC.BRT, filename = paste(species,"_BRT_UC.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  
  # MAXENT MODEL SUMMARY - CALCULATE MEAN VALUES FROM BOOTSTRAPS
  # Importance of Environmental predictors
  boot_array_ImpVar <- as.data.frame(array(0, dim = c(length(imp.var), n.boot)))
  row.names(boot_array_ImpVar) <- imp.var
  for (k in 1:n.boot){
    temp <- tmp3[[k]]$ImpVar.MX
    row.names(temp) <- sort(imp.var.tune[[i]])
    boot_array_ImpVar[rownames(temp),k] <- temp
  }
  influence_mat_MX <- t(apply(boot_array_ImpVar, 1, function(x) c(MX.mean = mean(x), MX.sd = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds
  # write.csv(influence_mat_MX, file=paste("preds.influences_",species, ".csv", sep = ""))
  
  # Model fit metrics
  boot_array_dev <- as.data.frame(array(0, dim = c(6, n.boot)))
  for (L in 1:n.boot){boot_array_dev[,L] <- tmp3[[L]]$SpeR2.MX}
  deviance_mat_MX <-t(apply(boot_array_dev, 1, function(x) c(MX.mean = mean(x), MX.sd = sd(x))))
  rownames(deviance_mat_MX) <- row.names(tmp3[[1]]$SpeR2.MX)
  # write.csv(deviance_mat_MX, file=paste("mean.model.fit_",species, ".csv", sep = ""))
  
  # mean spatial predictions 
  boot_mat <- array(0, dim = c(length.map, n.boot))
  for (M in 1:n.boot){boot_mat[,M] <- tmp3[[M]]$SpPred.MX}
  boot.mean.MX <-apply(boot_mat,1,mean)
  boot.sd.MX <-apply(boot_mat,1,sd)
  
  # RASTER OF MEAN SPATIAL PREDICTIONS
  mean.MX <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                      y =Pred_1km[,2], 
                                      z =boot.mean.MX),
                           crs = crs) # same proj as orginal env variables tiff files
  
  # RASTER OF SD OF THE MEAN SPATIAL PREDICTIONS
  UC.MX <- rasterFromXYZ(data.frame(x =Pred_1km[,1], 
                                    y =Pred_1km[,2], 
                                    z =boot.sd.MX),
                         crs = crs)
  
  # save rasters to folder
  writeRaster(mean.MX, filename = paste(species,"_MX_mean.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  writeRaster(UC.MX, filename = paste(species,"_MX_UC.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  
  # ENSEMBLE MODEL
  # Combine estimates of environmental importance 
  ENS.mean <- rowMeans(cbind(influence_mat_RF[,1], influence_mat_BRT[,1]), influence_mat_MX[,1])
  ENS.sd <- rowMeans(cbind(influence_mat_RF[,2], influence_mat_BRT[,2]), influence_mat_MX[,2])
  BRT_RF_env <- cbind(influence_mat_RF, influence_mat_BRT,influence_mat_MX, ENS.mean, ENS.sd)
  write.csv(BRT_RF_env, file=paste("Env_inf_",species, ".csv", sep = ""))
  
  # 1. model weighting 
  Model_WG <- data.frame(RF = double(1), BRT = double(1))
  if (is.na(deviance_mat_RF["AUC.eval",1])){
    Model_WG$RF<- deviance_mat_RF["AUC.train",1]/sum(deviance_mat_RF["AUC.train",1],deviance_mat_BRT["AUC.train",1], deviance_mat_MX["AUC.train",1])
    Model_WG$BRT<- deviance_mat_BRT["AUC.train",1]/sum(deviance_mat_BRT["AUC.train",1],deviance_mat_RF["AUC.train",1],deviance_mat_MX["AUC.train",1])
    Model_WG$MX<- deviance_mat_MX["AUC.train",1]/sum(deviance_mat_BRT["AUC.train",1],deviance_mat_RF["AUC.train",1],deviance_mat_MX["AUC.train",1])
  } else {
    Model_WG$RF<- deviance_mat_RF["AUC.eval",1]/sum(deviance_mat_RF["AUC.eval",1],deviance_mat_BRT["AUC.eval",1], deviance_mat_MX["AUC.eval",1])
    Model_WG$BRT<- deviance_mat_BRT["AUC.eval",1]/sum(deviance_mat_BRT["AUC.eval",1],deviance_mat_RF["AUC.eval",1], deviance_mat_MX["AUC.eval",1])
    Model_WG$MX<- deviance_mat_MX["AUC.eval",1]/sum(deviance_mat_BRT["AUC.eval",1],deviance_mat_RF["AUC.eval",1], deviance_mat_MX["AUC.eval",1])
  }
  
  # 2. Make grids of CV-based weights for BRT and RF
  BRT_WG_CV <- (1-UC.BRT) / ((1-UC.BRT) + (1-UC.RF) + (1-UC.MX))
  RF_WG_CV <- (1-UC.RF) / ((1-UC.BRT) + (1-UC.RF) + (1-UC.MX))
  MX_WG_CV <- (1-UC.MX) / ((1-UC.BRT) + (1-UC.RF) + (1-UC.MX))
  
  # 3. Combine CV-based weights with model-performance-based weights (taking the average of the two to give equal "weight" to each component of the weight)
  BRT_WG_CV <- (BRT_WG_CV + Model_WG$BRT) /2
  RF_WG_CV <- (RF_WG_CV + Model_WG$RF) /2
  MX_WG_CV <- (MX_WG_CV + Model_WG$MX) /2

  # 4. Make weighted ensemble grids
  Ensemble <- (mean.BRT*BRT_WG_CV) + (mean.RF *RF_WG_CV)+(mean.MX *MX_WG_CV)
  Ensemble.SD <- (UC.BRT*BRT_WG_CV) + (UC.RF *RF_WG_CV)+(UC.MX *MX_WG_CV)
  Ensemble_CV <- Ensemble_SD/Ensemble
  
  # save rasters to folder
  writeRaster(Ensemble, filename = paste(species,"_Ensemble.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  writeRaster(Ensemble_SD, filename = paste(species,"_Ensemble_SD.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  writeRaster(Ensemble_CV, filename = paste(species,"_Ensemble_CV.tif",sep=""),
              format = "GTiff", 
              overwrite = TRUE)
  
  print(paste("Iteration ",i, "  ; Taxon: ", species, sep = "" ))
}

end <- Sys.time()
end - start