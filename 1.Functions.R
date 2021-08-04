##==========================      1.Functions.R                             =====================================##

##--------------------------  Authors: Fabrice Stephenson & Owen Anderson
##--------------------------  Aim: functions used in subsequent R scripts for MS "Presence-only habitat 
##--------------------------  suitability models for Vulnerable Marine Ecosystem indicator taxa in the South 
##--------------------------  Pacific have reached their predictive limit"
##--------------------------  Project : SPR201801 - Variation - SC 2019-2020
##--------------------------  Start date : 01/05/2020
##--------------------------  End date : 03/08/2021
rm(list = ls())
##=============================================================================================================##

####=======================         LOAD FUNCTIONs               =============================================####

# SETTING UP FILE STRUCTURE FOR SAVING OUTPUTS FROM PARRALEL LOOP
# For each model type (RF;BRT;MX)
# ImpVar: Variables of importance (list of var, percent contribution)
# SpeR2: Model fit metrics
# SpPred: Spatial prediction of occurence

multiResultClass <- function(ImpVar.BRT=NULL, SpeR2.BRT=NULL, SpPred.BRT=NULL,
                             ImpVar.RF=NULL, SpeR2.RF=NULL, SpPred.RF=NULL,
                             ImpVar.MX =NULL, SpeR2.MX=NULL, SpPred.MX=NULL) {me <- list(
                               ImpVar.BRT = ImpVar.BRT,
                               SpeR2.BRT = SpeR2.BRT,
                               SpPred.BRT =  SpPred.BRT,
                               ImpVar.RF = ImpVar.RF,
                               SpeR2.RF = SpeR2.RF,
                               SpPred.RF =  SpPred.RF,
                               ImpVar.MX = ImpVar.MX,
                               SpeR2.MX = SpeR2.MX,
                               SpPred.MX =  SpPred.MX)
                             ## Set the name for the class
                             class(me) <- append(class(me),"multiResultClass")
                             return(me)
}


# FUNCTION TO CALCULATE EVALUTION METRICS
# model evaluation metrics calculated using training and evaluation data
# model = model object; model_type = "BRT"; "RF"; "Maxent"
# train = training dataframe - must have column "pa"
# eval = evaluation dataframe - must have column "pa"

Eval_Metrics <- function(model, model_type, train, eval) {
  if(model_type == "BRT"){
    # create blank data.frame for storing outputs
    deviance_DF <- as.data.frame(array(0, c(6,0)))
    rownames(deviance_DF) <- c("Dev.Exp.train","AUC.train","TSS.train", "Dev.Exp.eval","AUC.eval","TSS.eval")
    #internal deviance explained
    int.null.deviance <- model$self.statistics$mean.null 
    int.residual.deviance <- model$cv.statistics$deviance.mean 
    deviance_DF[1,1] <- (int.null.deviance-int.residual.deviance)/int.null.deviance
    #internal AUC
    deviance_DF[2,1] <- model$cv.statistics$discrimination.mean
    # internal TSS
    myROC <- pROC::roc(train$pa, model$fitted, quiet = T)
    Sens_spec <- pROC::coords(myROC, x="best", input="threshold", best.method="youden",transpose = FALSE)
    deviance_DF[3,] <- Sens_spec[1,3] + Sens_spec[1,2]-1 # TSS
    
    #model fit comparison with withheld evaluation data
    if (exists("eval_PA")){
      pred <- predict.gbm(model, eval, n.trees = model$gbm.call$best.trees, type = "response")
      #external deviance explained
      ext.residual.deviance <- calc.deviance(eval$pa, pred, family = "bernoulli" ,calc.mean=T) 
      ext.null.deviance <- calc.deviance(eval$pa,family = "bernoulli", rep(mean(eval$pa),nrow(eval)), calc.mean=T) 
      deviance_DF[4,]<-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance
      # external AUC
      myROC <- pROC::roc(eval$pa, pred, quiet = T)
      deviance_DF[5,] <- myROC$auc
      # external TSS
      Sens_spec <- pROC::coords(myROC, x="best", input="threshold", best.method="youden",transpose = FALSE)
      deviance_DF[6,] <- Sens_spec[1,3] + Sens_spec[1,2]-1 # TSS
    } else {
      print('no evaluation data')
    }
    colnames(deviance_DF) <- "BRT.eval"
    return(deviance_DF)
  } else if (model_type == "Maxent"){
    deviance_DF <- as.data.frame(array(0, c(6,0))) 
    rownames(deviance_DF) <- c("Dev.Exp.train","AUC.train","TSS.train", "Dev.Exp.eval","AUC.eval","TSS.eval")
    deviance_DF[2,1] <- as.numeric(model@results["Training.AUC",]) 
    pred <- predict(M3, train) 
    myROC <- pROC::roc(train$pa, pred, quiet = T)
    Sens_spec <- pROC::coords(myROC, x="best", input="threshold", best.method="youden",transpose = FALSE)
    deviance_DF[3,] <- Sens_spec[1,3] + Sens_spec[1,2]-1 # TSS
    
    if (exists("eval_PA")){
      pred <- predict(M3, eval) 
      myROC <- pROC::roc(eval$pa, pred, quiet = T)
      deviance_DF[5,] <- myROC$auc
      Sens_spec <- pROC::coords(myROC, x="best", input="threshold", best.method="youden",transpose = FALSE)
      deviance_DF[6,] <- Sens_spec[1,3] + Sens_spec[1,2]-1 # TSS
    } else {
      print('no evaluation data')
    }
    colnames(deviance_DF) <- "BRT.eval"
    return(deviance_DF)
  }  else if (model_type == "RF"){
    deviance_DF <- as.data.frame(array(0, c(6,0))) 
    rownames(deviance_DF) <- c("R2.train","AUC.train","TSS.train", "R2.eval","AUC.eval","TSS.eval")
    preds<-predict(model, type="prob")[,2]
    DATA<-data.frame(ID="ID",obs=train$pa,pred=preds)
    actual <- as.numeric(as.character(DATA$obs))
    predicted <- DATA$pred
    deviance_DF[1,1] <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
    myRoc <- pROC::roc(train$pa, preds, quiet = T)
    deviance_DF[2,1] <- myRoc$auc
    Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
    deviance_DF[3,1] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    
    if (exists("eval_PA")){
      preds<-predict(model, eval,type="prob")[,2]
      DATA<-data.frame(ID="ID",obs=eval$pa,pred=preds)
      actual <- as.numeric(as.character(DATA$obs))
      predicted <- DATA$pred
      deviance_DF[4,1] <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
      myRoc <- pROC::roc(eval$pa, preds, quiet = T)
      deviance_DF[5,1] <- myRoc$auc
      Sens_spec <- pROC::coords(myRoc, x="best", input="threshold", best.method="youden",transpose = FALSE)
      deviance_DF[6,1] <- Sens_spec[1,3]+Sens_spec[1,2]-1 # TSS
    } else {
      print('no evaluation data')
    }
    colnames(deviance_DF) <- "RF.eval"
    return(deviance_DF)
  } else {
    print("Error: enter either BRT, RF or Maxent as model_type")
  }
}

# FUNCTION TO CALCULATE ENVIRONMENTAL PREDICTOR IMPORTANCE
# Importance of environmental predictor in BRT, RF and Maxent models
# model = model object
# model_type = "BRT", "RF", "Maxent"
Var_Imp <- function(model, model_type) {
  if(model_type == "BRT"){
    M1_contrib <- as.data.frame(model$contributions) 
    env_var_ord <- M1_contrib[order(M1_contrib$var),]
    env_var_ord <- M1_contrib[,2, drop = F]
    colnames(env_var_ord) <- "BRT_rel.Inf"
    return(env_var_ord)
  } else if (model_type == "RF"){
    M1_contrib <- as.data.frame(importance(model, type=2))
    M1_contrib$norm <- 0
    for (i in 1:nrow(M1_contrib)){
      M1_contrib$norm[i] <- (M1_contrib[i,1] * 100)/sum(M1_contrib[,1])
    }
    env_var_ord <- as.data.frame(M1_contrib[order(row.names(M1_contrib[])),])
    env_var_ord <- env_var_ord[,2, drop=FALSE]
    colnames(env_var_ord) <- "RF_rel.Inf"
    return(env_var_ord)
  } else if (model_type == "Maxent"){
    env_var_ord <- as.data.frame(M3@results[c(17:26),1])
    colnames(env_var_ord) <- "RF_rel.Inf"
    return(env_var_ord)
  } else {
    print("Error: enter either BRT, Maxent or RF as model_type")
  }
}