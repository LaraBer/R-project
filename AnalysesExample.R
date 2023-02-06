#advice: in functions fitLinMod_Interactions() and fitLinMod_MainEff() you can specify
#the Mixed Models (MixedModel) for which analyses are run. 
#the function uses the defined bootstrap functions automatically
#which take as argument the completely specified ixed model



#load packages

p <- c("ggplot2", "lmerTest","BayesFactor", "standardize", "cluster", "plyr", "stats", "gridExtra", "ggplot2", "corrplot", "dplyr", "PairedData", 
       "boot","boot.pval", "lmboot","lme4", "Hmisc", "effects", "gplots", "MASS", "car", "reshape2", "rstatix", "ggpubr", "tidyverse", "data.table",
       "MuMIn", "grid", "png", "gghalves", "ggpubr", "ggtext", "colorBlindness", "latex2exp", "psych", "caTools", "corrgram", "robustbase", "robustlmm")

invisible(lapply(p, library, character.only = TRUE))


# functions ---------------------------------------------------------------


#outlier truncation: replace outliers with outlier cutoff value
chop_off_outliers <- function(var, dataframe, varname){
  out <- boxplot.stats(var)$out
  quantiles <- quantile(var)
  IQR <- as.numeric(quantiles[4]-quantiles[2])
  LB <- as.numeric(quantiles[2]-(1.5*IQR))
  UB <- as.numeric(quantiles[4]+(1.5*IQR))
  dataframe[var %in% out,varname] <- UB
  return(dataframe)
}



#plotting data: mixed effects
plot_mixedModel2Cat_woPoints <- function(data, IVwithin, DV, IVbetween, title, labelY, labelX, minValue, maxValue, colors, breaks){
  data %>% 
    ggplot() +
    aes(x = IVbetween, color = IVwithin, shape = IVwithin, y = DV, group = IVwithin, linetype = IVwithin) +
    stat_summary(fun = mean, geom = "point", size = 2.5, position = position_dodge(0.4)) + 
    stat_summary(fun = mean, geom = "line", position = position_dodge(0.4)) + theme_classic() +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.4, position = position_dodge(0.4), inherit.aes = FALSE, aes(x = IVbetween, color = IVwithin, shape = IVwithin, y = DV, group = IVwithin))+
    scale_color_manual(name = title, values=colors, breaks = breaks)+ scale_shape_manual(name = title, values = c(17, 15), breaks = breaks)+
    scale_linetype_manual(name = title, values=c("dotted", "solid"), breaks = breaks)+
    ylab(labelY) + scale_y_continuous(limits = c(minValue, maxValue))+ xlab(labelX)+
    theme(legend.position = "top", legend.title = element_text(size = 25), legend.text = element_text(size = 20), axis.title.y = element_text(size = 30, color = "gray12"), axis.text.y = element_text(size = 25, color = "gray12"), axis.title.x = element_text(size = 30, color = "gray12"), axis.text.x = element_text(size = 25, color = "gray12")) 
}


#fit and plot linear mixed models, write them into txt and pgn files in designated folders
#Interaction effects (returns only significant interaction effects)
fitLinMod_Interactions <- function(dat, dependent,FileName, to_test, to_plot, wd){
  data <- dat
  for (var in dependent){ #do this for each dependent variable you want to test
    setwd(wd)
    if(dir.exists(var)){ 
      unlink(var)
    }
    dir.create(var) #create a folder for results for each DV
    setwd(var) #go into this directory
    DV = data[,var] #select column corresponding to DV in data set
    for(varset in to_test){ #for each set of predictors
      varsetCat_name <- to_plot[which(to_test == varset)] #get the corresponding categorized predictor set
      varsetCat <- eval(as.name(varsetCat_name))
      varsetCont <- eval(as.name(varset))
      for (variable in names(varsetCont)){ #loop through all variables in set of predictors 
        curVarCont <- scale(data[,variable])
        MixedModel <- lmer(DV ~ curVarCont*condition + (1|id) , data = data) #specify mixed model with selected variables
        cat(paste(var, "\n", sep = ""), capture.output(MixedModel), file=FileName, sep="\n", append=TRUE) #write un-bootstrapped output to file
        p_val <- bootstrapAnova(variable, FileName, MixedModel, data)
        p_val <- as.numeric(unlist(p_val))
        if (p_val[3]<=0.05){ 
          curVarCat <- names(varsetCat[which(names(varsetCont) == variable)])
          curVarCat <- data[,curVarCat]
          p <- (plot_mixedModel2Cat_woPoints(data, as.factor(curVarCat), DV, data$condition, variable, 
                                             var, "condition", min(DV), max(DV), c("#009E73","#D55E00"), c("low", "high")))
          png(paste(variable, ".png", sep=""))
          plot(p)
          dev.off()
          
        }
      }
      
    }
    setwd(wd)
    unlink(var)
    
  }
}


#main effects (returns only significant main effects)
fitLinMod_MainEff <- function(dat, dependent,FileName, to_test, to_plot, wd){
  data <- dat
  for (var in dependent){ #do this for each dependent variable you want to test
    setwd(wd)
    if(dir.exists(var)){ 
      unlink(var)
    }
    dir.create(var) #create a folder for results for each DV
    setwd(var) #go into this directory
    DV = data[,var] #select column corresponding to DV in data set
    
    
    for(varset in to_test){ #for each set of predictors
      varsetCat_name <- to_plot[which(to_test == varset)] #get the corresponding categorized predictor set
      varsetCat <- eval(as.name(varsetCat_name))
      varsetCont <- eval(as.name(varset))
      for (variable in names(varsetCont)){ #loop through all variables in set of predictors 
        curVarCont <- scale(data[,variable])
        MixedModel <- lmer(DV ~ curVarCont*condition + (1|id), data = data) #fit linear model with this predictor
        cat(paste(var, "\n", sep = ""), capture.output(MixedModel), file=FileName, sep="\n", append=TRUE)
        p_val <- bootstrapAnova(variable, FileName, MixedModel, data) #feed mixed model into bootsrap anova function
        p_val <- as.numeric(unlist(p_val))
        if (p_val[1]<=0.05){ #only for significant main effects
          
          xlabel <- variable
          ylabel <- var 
          maxi <- max(DV)
          
          p <- ggplot(data, aes(x = curVarCont, y = DV, color = curVarCont)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
            geom_smooth(method=lm, color = "blue") + labs(x=xlabel, y = ylabel)+
            theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
            geom_density_2d(color = "blue") + scale_y_continuous(limits = c(0, maxi))
          
          
          png(paste(variable, ".png", sep=""))
          plot(p)
          dev.off()
          
        }
      }
      
    }
    setwd(wd)
    unlink(var)
    
  }
}

#bootstrap functions

#resampling for bootstrapping
performResampling <- function(data){
  #these are the different condition orders that could occur
  data_012 <- data[data$conditionOrder == "012",]
  data_021 <- data[data$conditionOrder == "021",]
  data_102 <- data[data$conditionOrder == "102",]
  data_120 <- data[data$conditionOrder == "120",]
  data_201 <- data[data$conditionOrder == "201",]
  data_210 <- data[data$conditionOrder == "210",]
  samplesize = nrow(data)/6 #we have six different orders and want to have them equally represented
  #from each data subset (six = 1 for each condition order) draw sample with standardized size
  dat_012 <- data[sample(x = nrow(data_012), size = samplesize, replace = TRUE),]
  dat_021 <- data[sample(x = nrow(data_021), size = samplesize, replace = TRUE),]
  dat_102 <- data[sample(x = nrow(data_102), size = samplesize, replace = TRUE),]
  dat_120 <- data[sample(x = nrow(data_120), size = samplesize, replace = TRUE),]
  dat_201 <- data[sample(x = nrow(data_201), size = samplesize, replace = TRUE),]
  dat_210 <- data[sample(x = nrow(data_210), size = samplesize, replace = TRUE),]
  #bind the resulting resampled data subframes together
  d <- cbind.data.frame(dat_012, dat_021, dat_102, dat_120, dat_201, dat_210)
  return(d)
}

#prepare empty data frame to store bootstrap results over which will be averaged to obtain results
prepareBootstrap <- function(listOfOutputPars,listOfCoeffs, numberPars, numberCoeffs){
  Pr <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  F.val <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  SumSq <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  MeanSumSq <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  NumDF <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  DenDF <- as.data.frame(rep(NA, numberPars), row.names = listOfOutputPars)
  coefficients_model <- as.data.frame(rep(NA,numberCoeffs), listOfCoeffs)
  estimate_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  stErr_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  tVal_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  pVal_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  lowerCI_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  upperCI_model <- as.data.frame(rep(NA, numberCoeffs), listOfCoeffs)
  return(list(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model, estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model))
}

#function to update data frame (add result from new bootstrap)
updateBootstrap <- function(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model,estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model, a, model, c, ci){
  Pr <- cbind.data.frame(Pr, a$`Pr(>F)`)
  F.val <- cbind.data.frame(F.val, a$`F value`)
  SumSq <- cbind.data.frame(SumSq, a$`Sum Sq`)
  MeanSumSq <- cbind.data.frame(MeanSumSq, a$`Mean Sq`)
  NumDF <- cbind.data.frame(NumDF, a$`NumDF`)
  DenDF <- cbind.data.frame(DenDF, a$`DenDF`)
  coefficients_model <- cbind.data.frame(coefficients_model, coef(summary(model))[,"Estimate"])
  estimate_model <- cbind.data.frame(estimate_model, c$Estimate)
  stErr_model <- cbind.data.frame(stErr_model, c$`Std. Error`)
  tVal_model <- cbind.data.frame(tVal_model, c$`t value`)
  pVal_model <- cbind.data.frame(pVal_model, c$`Pr(>|t|)`)
  lowerCI_model <- cbind.data.frame(lowerCI_model, ci$`2.5 %`)
  upperCI_model <- cbind.data.frame(upperCI_model, ci$`97.5 %`)
  return(list(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model, estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model))
}

#compute the mean over bootstrap results
meanBootstrap <- function(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model, estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model){
  Pr <- as.data.frame(Pr)
  Pr_F <- rowMeans(Pr[,-1])
  F.val <- as.data.frame(F.val)
  F_val <- rowMeans(F.val[,-1])
  SumSq <- as.data.frame(SumSq)
  Sum_Sq <- rowMeans(SumSq[,-1])
  MeanSumSq <- as.data.frame(MeanSumSq)
  MeanSum_Sq <- rowMeans(MeanSumSq[,-1])
  NumDF <- as.data.frame(NumDF)
  Num_DF <- rowMeans(NumDF[,-1])
  DenDF <- as.data.frame(DenDF)
  Den_DF <- rowMeans(DenDF[,-1])
  coefficients_model <- as.data.frame(coefficients_model)
  coeffi <- rowMeans(coefficients_model[,-1])
  estimate_model <- as.data.frame(estimate_model)
  esti <- rowMeans(estimate_model[,-1])
  stErr_model <- as.data.frame(stErr_model)
  stErr <- rowMeans(stErr_model[,-1])
  tVal_model <- as.data.frame(tVal_model)
  tVal <- rowMeans(tVal_model[,-1])
  pVal_model <- as.data.frame(pVal_model)
  pVal <- rowMeans(pVal_model[,-1])
  lowerCI_model <- as.data.frame(lowerCI_model)
  lowerCI <- rowMeans(lowerCI_model[,-1])
  upperCI_model <- as.data.frame(upperCI_model)
  upperCI <- rowMeans(upperCI_model[,-1])
  return(list(Pr_F, F_val, Sum_Sq, MeanSum_Sq, Num_DF, Den_DF, coeffi, esti, stErr, tVal, pVal, lowerCI, upperCI))
}

#function writing bootstrap results (average) to file
writeBootResults <- function(F_val, Pr_F, Sum_Sq, MeanSum_Sq, Num_DF, Den_DF, coeffi, esti, stErr, tVal, pVal, lowerCI, upperCI, pathToFile, var){
  cat(paste(var, "\n F-val", sep = ""), as.character(F_val), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n", "Statistical significance", "\n", sep = ""), as.character(Pr_F), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","SumSq", "\n", sep = ""), as.character(Sum_Sq), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","MeanSumSq", "\n", sep = ""), as.character(MeanSum_Sq), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","NumDF", "\n", sep = ""), as.character(Num_DF), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","DenDF", "\n", sep = ""), as.character(Den_DF), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","Coefficients", "\n", sep = ""), as.character(coeffi), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","Estimate", "\n", sep = ""), as.character(esti), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","stdErr", "\n", sep = ""), as.character(stErr), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","tVal", "\n", sep = ""), as.character(tVal), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","pVal", "\n", sep = ""), as.character(pVal), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","lowerCi", "\n", sep = ""), as.character(lowerCI), file=pathToFile, sep="\n", append=TRUE)
  cat(paste("\n","upperCI", "\n", sep = ""), as.character(upperCI), file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
}

#bootstrap anova makes use of the predefined functions: bootraps anova results and writes significant results to file
bootstrapAnova <- function(var, FileName, MixedModel, data){
  pathToFile <- FileName
  dat <- data
  DV <-  dat[,var] #select column corresponding to DV in data set
  model <- MixedModel #fit model
  a <- anova(model) #anova of teh model, returning things Liek SS, p, F...
  ci <- as.data.frame(confint(model)) #confidence interval of model coefficients
  ci <- ci[-c(1,2),] #delete sigmal rows
  c <- summary(model) #summarize the model for model coefficients
  c <- as.data.frame(c$coefficients) #need it in data frame format
  listOfCoeffs<-rownames(c) #these are the coefficients, include intercept
  listOfOutputPars <- rownames(a) #excludes intercept, as not included in anova
  numberPars <- length(listOfOutputPars) #number of anova parameters -1 compared to coefficients as excludes intercept
  numberCoefs <- length(listOfCoeffs) 
  #listOfCoeffNames <- rownames(coef(summary(model)))
  numberCoeffs <- length(rownames(coef(summary(model))))
  bootParameters <- prepareBootstrap(listOfOutputPars,listOfCoeffs, numberPars, numberCoeffs)
  Pr <- bootParameters[1] 
  F.val <- bootParameters[2]
  SumSq <- bootParameters[3]
  MeanSumSq <- bootParameters[4]
  NumDF <- bootParameters[5]
  DenDF <- bootParameters[6]
  coefficients_model <- bootParameters[7]
  estimate_model <- bootParameters[8]
  stErr_model <- bootParameters[9]
  tVal_model <- bootParameters[10]
  pVal_model <- bootParameters[11]
  lowerCI_model <- bootParameters[12]
  upperCI_model <- bootParameters[13]
  for(i in 1:100){
    dat <- as.data.frame(performResampling(data))
    DV <-  dat[,var] #select column corresponding to DV in data set
    model <- MixedModel
    a <- anova(model)
    bootParameters <- updateBootstrap(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model, estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model, a, model, c, ci)
    Pr <- bootParameters[1] 
    F.val <- bootParameters[2]
    SumSq <- bootParameters[3]
    MeanSumSq <- bootParameters[4]
    NumDF <- bootParameters[5]
    DenDF <- bootParameters[6]
    coefficients_model <- bootParameters[7]
    estimate_model <- bootParameters[8]
    stErr_model <- bootParameters[9]
    tVal_model <- bootParameters[10]
    pVal_model <- bootParameters[11]
    lowerCI_model <- bootParameters[12]
    upperCI_model <- bootParameters[13]
    
  }
  means <- meanBootstrap(Pr, F.val, SumSq, MeanSumSq, NumDF, DenDF, coefficients_model, estimate_model, stErr_model, tVal_model, pVal_model, lowerCI_model, upperCI_model)
  Pr_F <- means[1]
  F_val <- means[2]
  Sum_Sq <- means[3]
  MeanSum_Sq <- means[4]
  Num_DF <- means[5]
  Den_DF <- means[6]
  coeffi <- means[7]
  esti <- means[8]
  stErr <- means[9]
  tVal <- means[10]
  pVal <- means[11]
  lowerCI <- means[12]
  upperCI <- means[13]
  p <- as.array(as.numeric(unlist(Pr_F)))
  if (any(p<=0.05)){ 
    writeBootResults(F_val, Pr_F, Sum_Sq, MeanSum_Sq, Num_DF, Den_DF, coeffi, esti, stErr, tVal, pVal, lowerCI, upperCI, pathToFile, var)
  }
  return(Pr_F)
}

#applies bootstrap anova function to array of dependent variables
multiplePredictorModel <- function(dependentVar, FileName, MixedModel, data){
    bootstrapAnova(dependentVar, FileName, MixedModel, data)
}

#does the same as multiplePredictorModel but only takes single dependent variable as argument
multiplePredictorModel_woEff <- function(dependentVar, FileName, MixedModel, data){
  
  bootstrapAnova(dependentVar, FileName, MixedModel, data)
  
}


#function
ML_Regression <- function(data_train, data_test, RegressionModel, FileName, predictedVar){
  #split data into train and test set
  pathToFile <- FileName
  trainSet <- data_train
  testSet <- data_test
  #train model on training set
  model <- RegressionModel
  #summary(model)
  #vif(model)
  
  #Plot model residuals
  modelResiduals <- as.data.frame(residuals(model)) 
  #Here we expect to see something approximately normally distributed.
  ggplot(modelResiduals, aes(residuals(model))) +
    geom_histogram(fill='deepskyblue', color='black')
  
  #predict points based on trained model and compare with actual points in test set
  preds <- predict(model, testSet)
  modelEval <- cbind(testSet[,predictedVar], preds)
  colnames(modelEval) <- c('Actual', 'Predicted')
  modelEval <- as.data.frame(modelEval)
  R2 <- r.squaredGLMM(model)
  
  #evaluate model
  mae <- mean(modelEval$Actual - modelEval$Predicted)
  mse <- mean((modelEval$Actual - modelEval$Predicted)^2)
  rmse <- sqrt(mse)
  
  significance <- capture.output(anova(model))
  results <- capture.output(summary(model))
  
  cat(paste(predictedVar, "\n", sep = ""), results, file=pathToFile, sep="\n", append=TRUE)
  cat(paste("Statistical significance", "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
  cat(paste("MSE", "", sep = ""), mse, file=pathToFile, sep="\n", append=TRUE)
  cat(paste("RMSE", "", sep = ""), rmse, file=pathToFile, sep="\n", append=TRUE)
  cat(paste("R2", "\n", sep = ""), R2, file=pathToFile, sep="\n", append=TRUE)
  write("\n\n",file=pathToFile,append=TRUE) 
}


compute_regressions_Generic <- function(dat, dependentVars, FileName, Predictors){
  for(var in dependentVars){
    sampleSplit <- sample.split(Y=dat[,var], SplitRatio=0.7)
    trainSet <- as.data.frame(subset(x=dat, sampleSplit==TRUE))
    testSet <- as.data.frame(subset(x=dat, sampleSplit==FALSE))
    dependentVar <- trainSet[,var]
    flm <- reformulate(Predictors, var)
    RegressionModel <- lm(flm, data <- trainSet)
    
    ML_Regression(trainSet, testSet, RegressionModel, FileName, var)
  }
  
}


#summary function for bootstrap
mySummary1 <- function(.){s <- sigma(.)
c(beta = getME(., "beta"), sigma=s, sig01 = unname(s*getME(.,"theta")))}


  
 # model <- MixedModel
 # a <- anova(model)
  #pathToFile = FileName
  #if (any(a$`Pr(>F)`<=0.05)){ 
   # significance <- capture.output(anova(model))
    #results <- capture.output(summary(model))
    #R2 <- r.squaredGLMM(model)
    #cat(paste("secondsForBonusThisRound", "\n", sep = ""), results, file=pathToFile, sep="\n", append=TRUE)
    #cat(paste("Statistical significance", "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    #cat(paste("R2", "\n", sep = ""), R2, file=pathToFile, sep="\n", append=TRUE)
    #write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  #}
#}




#load data
data <- as.data.frame(read.csv("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/DataAnalysis/Analyses/UeberarbeiteteScoringFiles/Spef/2022-feb-experiment-extended-data-withEIG.csv", na.strings = "", stringsAsFactors = FALSE, header=T, sep = ",")) 

#id must be factor
data$id <- as.factor(data$id)

#create categorized variables for visualization
continuous <- c("DominanceTrait", "ValenceTrait", "ArousalTrait", "ValenceState", "DominanceState", "ArousalState", "NeedForCognition", "FaithInIntuition",
                "Numeracy", "WorkingMemory", "Agreeableness", "Conscientiousness", "Extraversion", "OpenMind", "NegEmotionality")
categorical <- c("DominanceTraitCat", "ValenceTraitCat", "ArousalTraitCat", "ValenceStateCat", "DominanceStateCat", "ArousalStateCat", "NfCCat", "FiICat", 
                 "NumCat", "WMCat", "AgrCat", "ConsCat", "ExtraCat", "OpenMindCat", "NegEmoCat")

for (i in 1:length(continuous)){
  cont = data[,continuous[i]]
  cat = categorical[i]
  data[,cat] = rep(0, nrow(data))
  data[, cat][cont >= median(cont)] <- "high"
  data[, cat][cont < median(cont)] <- "low"
}


#demographics

#get demographic information from long data set
sample_size <- length(unique(data$id))

uniqueIDs <- unique(data$id)
ages <- rep(0, sample_size)
genders <- rep(0, sample_size)

demographics <- as.data.frame(cbind(uniqueIDs, ages, genders))

for (id in uniqueIDs){
  cur = data[data$id == id,]
  cur = cur[1,]
  age = 2022 - cur$YoB
  gender = cur$gender
  demographics[demographics$uniqueIDs==id,"ages"]= age
  demographics[demographics$uniqueIDs==id,"genders"]= gender
  
}

#number of males, females, queers, NA in sample
females = nrow(demographics[demographics$genders == "f",])
males = nrow(demographics[demographics$genders == "m",])
queer = nrow(demographics[demographics$genders == "q",])

#mean and sd age
M_age = mean(demographics$ages)
SD_age = sd(demographics$ages)


#Attention check
#load the raw data
attentionCheck <- as.data.frame(read.csv("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/DataAnalysis/Analyses/CompleteDataFiles/Spef/data_speedefficiency_2022-02-23_12-55.csv", na.strings = "", stringsAsFactors = FALSE, header=T, sep = ";")) 

AC1 <- attentionCheck$CASE[attentionCheck$RE01_11 != 1]

#check if participants passed attention check
for(i in 1:nrow(attentionCheck)){
  cur <- attentionCheck$AC01_01[i]
  id <- attentionCheck$CASE[i]
  check1 <- grepl("astermind", cur, fixed=TRUE)
  check2 <- grepl("ASTERMIND", cur, fixed=TRUE)
  if((! check1) & (! check2)){
    append(AC1, id)
  }
}

data$AttentionCheck <- rep(1, nrow(data))
data$AttentionCheck[data$id %in% AC1] <- 0


data <- data[data$AttentionCheck==1,]


#check validity of game play data: outliers in number of guesses, time spent on games
#detect outliers
out_numGuess <- boxplot.stats(data$numguesses)$out
out_rt <- boxplot.stats(data$rttotal)$out

#all outliers are replaced by the value that defines outlier cutoff value (kind of truncation)
data <- chop_off_outliers(data$numguesses,data, "numguesses")
data <- chop_off_outliers(data$rttotal, data, "rttotal")

#compute step by step parameters

#first create new column codejarML with the most likely item in codejar
data$codejarML <- gsub("]","",data$codejar)
data$codejarML <- sub(".","",data$codejarML)
data$codejarML <- gsub(";", ",", data$codejarML)

len <- length(data$codejarML)
for (i in 1:len){
  one <- data$codejarML[i] 
  one <- as.numeric(unlist(strsplit(one, ",")))
  data$codejarML [i] <- which.max(one)
  
}



#new colum for guesses, rts and feedback (happy, neutral, sad) new in which we delete all non-integer symbols

data$guessesNew <- gsub("]","",data$guesses, fixed = TRUE)
data$guessesNew <- gsub(",", "", data$guessesNew)
data$guessesNew <- gsub("[","",data$guessesNew, fixed = TRUE)


data$rts <- gsub("]","",data$rt_by_step)
data$rts <- gsub(",","",data$rts)
data$rts <- gsub("[","",data$rts, fixed = TRUE)


data$f_p <- gsub("]","",data$feedback_smiley)
data$f_p <- gsub(",","",data$f_p)
data$f_p <- gsub("[","",data$f_p, fixed = TRUE)

data$f_m <- gsub("]","",data$feedback_neutral)
data$f_m <- gsub(",","",data$f_m)
data$f_m <- gsub("[","",data$f_m, fixed = TRUE)


data$f_n <- gsub("]","",data$feedback_frownie)
data$f_n <- gsub(",","",data$f_n)
data$f_n <- gsub("[","",data$f_n, fixed = TRUE)

#extract guess, feedback and the kind of changes participants made to guess in response to feedback (stepwise)
for (j in 1:nrow(data)){
  g = str_split(data$guessesNew[j], pattern = " ")
  rt = str_split(data$rts[j], pattern = " ")
  fp = str_split(data$f_p[j], pattern = " ")
  fm = str_split(data$f_m[j], pattern = " ")
  fn = str_split(data$f_n[j], pattern = " ")
  ml = data$codejarML[j]
  g_l = g[[1]]
  rt_l = rt[[1]]
  fp_l = fp[[1]]
  fm_l = fm[[1]]
  fn_l = fn[[1]]
  counter = 0
  for(i in 1:length(g_l)){
    if(i %% 3 == 1){
      new_g = g_l[i]
    }
    if(i %% 3 == 2){
      new_g = append(new_g, g_l[i])
    }
    if(i %% 3 == 0){
      counter = counter + 1
      new_rt = rt_l[counter]
      new_feedback = c(fp_l[counter], fm_l[counter], fn_l[counter])
      new_g = append(new_g, g_l[i])
      nam <- paste("guess", counter, sep = "")
      nam2 <- paste("nr_ML", counter, sep = "")
      nam3 <- paste("rt", counter, sep = "")
      nam4 <- paste("fb", counter, sep = "")
      nam5 <- paste("nr_diff", counter, sep = " ")
      data[j, nam] = str_flatten(new_g)
      data[j, nam2] = str_count(str_flatten(new_g), ml)
      data[j, nam3] = str_flatten(new_rt)
      data[j, nam4] = str_flatten(new_feedback)
      data[j, nam5] = length((unique(new_g)))
      if (counter > 1 && !is.na(new_g)){
        g2 = str_split(data[j, nam], "")[[1]]
        id_prev_guess = paste("guess", counter-1, sep = "")
        g1 = str_split(data[j, id_prev_guess], "")[[1]]
        same_position = g1[g1 == g2]
        same_kind <- 0
        for (k in 1:length(g1)){
          item = g1[k]
          for (l in 1:length(g2)){
            if (item == g2[l]){
              same_kind = same_kind + 1
              g2 = g2[-l]
              break 
            }
          }
        }
        color_changes = paste("C", 3-same_kind, sep = "")
        position_changes = paste("L",same_kind-length(same_position), sep = "")
        q_type = paste(color_changes, position_changes, sep = "")
        change_index = paste(counter-1, counter, sep = "")
        change_index = paste("question_fb_type", change_index, sep = "")
        data[j, change_index] = q_type
      }
    }
  }
  
}

data$fb1 <- as.factor(data$fb1)
data$question_fb_type12 <- as.factor(data$question_fb_type12)
data$fb1 <- as.factor(data$fb2)
data$question_fb_type12 <- as.factor(data$question_fb_type23)

dataGuesses <- data[,c("nr_diff 1", "nr_diff 2", "nr_diff 3")]

data$meanNrDiffItems <- rowSums(dataGuesses, na.rm = TRUE)


#define the set of predictor variables you want to run analyses for
#continuous (for analyses)
cognitive <- data[,c("Numeracy", "WorkingMemory")]

thinkingstyles <- data[,c("FaithInIntuition", "NeedForCognition")]

personality <- data[,c("Agreeableness", "Conscientiousness", "Extraversion", 
                     "NegEmotionality", "OpenMind")]

emotions <- data[,c("DominanceTrait", "DominanceState", "ValenceTrait", "ValenceState", "ArousalTrait", 
                       "ArousalState")]

#second-order array with array of parameters
to_test <- c("personality", "emotions", "cognitive", "thinkingstyles")


#categorical (for visualization)
cognitiveCat <- data[,c("WMCat","NumCat")]

thinkingstylesCat <- data[,c("FiICat", "NfCCat")]

personalityCat <- data[,c("AgrCat", "ConsCat", "ExtraCat",
                        "NegEmoCat", "OpenMindCat")]

emotionsCat <- data[,c("DominanceTraitCat", "DominanceStateCat", "ValenceTraitCat", "ValenceStateCat", "ArousalTraitCat", 
                          "ArousalStateCat")]

#second order array
to_plot <- c("personalityCat", "emotionsCat", "cognitiveCat", "thinkingstylesCat")


#define dependent variables
dependent_allConditions <- c("pointsThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", "rt_first", "distance_to_previous_query_mean", "EIG0", "nr_diff 1", "meanNrDiffItems")



#get parameters required for resampling (bootstrap) to account for potential order effects
#counts for each condition
ns <- nrow(data[data$condition == "speed",])
ne <- nrow(data[data$condition == "efficiency",])
nm <- nrow(data[data$condition == "mixed",])

#proportion of games where speed/eff/mixed came 1erst/2nd/3rd
p1s <- nrow(data[data$indexInRound == 0 & data$condition == "speed",])/ns
p1e <- nrow(data[data$indexInRound == 0 & data$condition == "efficiency",])/ne
p1m <- nrow(data[data$indexInRound == 0 & data$condition == "mixed",])/nm

p2s <- nrow(data[data$indexInRound == 1 & data$condition == "speed",])/ns
p2e <- nrow(data[data$indexInRound == 1 & data$condition == "efficiency",])/ne
p2m <- nrow(data[data$indexInRound == 1 & data$condition == "mixed",])/nm

p3s <- nrow(data[data$indexInRound == 2 & data$condition == "speed",])/ns
p3e <- nrow(data[data$indexInRound == 2 & data$condition == "efficiency",])/ne
p3m <- nrow(data[data$indexInRound == 2 & data$condition == "mixed",])/nm


#idea: as (1/x)*x = 1, we assign the inverse of the occurrence probability for resampling
data$weighting <- rep(NA, nrow(data))
data$weighting[data$indexInRound == 0 & data$condition == "speed"] <- 1/p1s
data$weighting[data$indexInRound == 0 & data$condition == "efficiency"] <- 1/p1e
data$weighting[data$indexInRound == 0 & data$condition == "mixed"] <- 1/p1m

data$weighting[data$indexInRound == 1 & data$condition == "speed"] <- 1/p2s
data$weighting[data$indexInRound == 1 & data$condition == "efficiency"] <- 1/p2e
data$weighting[data$indexInRound == 1 & data$condition == "mixed"] <- 1/p2m

data$weighting[data$indexInRound == 2 & data$condition == "speed"] <- 1/p3s
data$weighting[data$indexInRound == 2 & data$condition == "efficiency"] <- 1/p3e
data$weighting[data$indexInRound == 2 & data$condition == "mixed"] <- 1/p3m


#get new clumn with condition order as we need this information for resampling function
data$conditionOrder <- as.factor(paste(data$indexInRound_speed, data$indexInRound_efficiency, data$indexInRound_mixed, sep = ""))


#rt mean wo first and distance to previous query mean contain NAs, replace them with 0 as these cases did not need more than 1 query
data$rt_mean_wo_first[is.na(data$rt_mean_wo_first)] <- 0
data$distance_to_previous_query_mean[is.na(data$distance_to_previous_query_mean)] <- 0


data_backup <- data

set.seed(22)

# Interaction and Main Effects for all predictors including visualization -----------

#create directory to store interaction effects
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Interactions")
#dir.create("withEIG") #if you get a warning, you may have already created the directory

#interaction effects for all dependent variables and all predictors
fitLinMod_Interactions(data, dependent_allConditions,"resultsSingleLMs", to_test, to_plot, "C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Interactions")

#run the same analysis for seconds, but here the efficiency condition must be excluded (time bonus did not vary here)
dependent_woEfficiency <- c("secondsForBonusThisRound")

fitLinMod_Interactions(data[data$condition!="efficiency",], dependent_woEfficiency,"resultsSingleLMs", to_test, to_plot, "C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Interactions")


#create directory for main effects
dependent <- c("pointsThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", "rt_first", "distance_to_previous_query_mean", "EIG0")
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/MainEffects")

#dir.create("withEIG")

#
fitLinMod_MainEff(data, dependent,"resultsSingleLMs", to_test, to_plot, "C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Main Effects")

#seconds for Bonus
fitLinMod_MainEff(data[data$condition!="efficiency",], dependent_woEfficiency,"resultsSingleLMs", to_test, to_plot, "C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Main Effects")





# Incremental variance ----------------------------------------------------
set.seed(22)
#check for which psychological variables emotions explain additional variance beyond numeracy and working memory
dependent <- c("pointsThisRound","secondsForBonusThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first", "EIG0", "nr_ML1") #nr ml 1 muss erst später im skript berechnet werden




#create new directory
#setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models")
#dir.create("Incremental Variance Explained by Psychological Variables")

#State emotions

setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models/Incremental Variance Explained by Psychological Variables")

for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarEmoState.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + DominanceState*condition + ValenceState*condition +ArousalState*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woEmo, model_withEmo)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + DominanceState*condition + ValenceState*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF)  # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

#explore state emo model further
set.seed(22)
#base vs. state emotion only
for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarEmoState_vs_base.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woEmo <- lmer(DV ~ condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withEmo <- lmer(DV ~ DominanceState*condition + ValenceState*condition +ArousalState*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woEmo, model_withEmo)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ DominanceState*condition + ValenceState*condition + ArousalState*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF)  # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

set.seed(22)

#state emotion only vs cognitive only
for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "EmoState_vs_Cognitive.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withEmo <- lmer(DV ~ DominanceState*condition + ValenceState*condition +ArousalState*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woEmo, model_withEmo)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ DominanceState*condition + ValenceState*condition + ArousalState*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF)  # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

set.seed(22)

#state emotion only vs. cognitive + state emotion model
for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "EmoState_vs_EmoStateCognitive.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woEmo <- lmer(DV ~ DominanceState*condition + ValenceState*condition +ArousalState*condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + DominanceState*condition + ValenceState*condition +ArousalState*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woEmo, model_withEmo)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + DominanceState*condition + ValenceState*condition +ArousalState*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ DominanceState*condition + ValenceState*condition +ArousalState*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF)  # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

#Trait emotions
set.seed(22)
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models/Incremental Variance Explained by Psychological Variables")

for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarEmoTrait.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withEmo <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + DominanceTrait*condition + ValenceTrait*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woEmo, model_withEmo)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + DominanceTrait*condition + ValenceTrait*condition +ArousalTrait*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF)  # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}


#Personality
set.seed(22)
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models/Incremental Variance Explained by Psychological Variables")

for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarPersonality.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woPers <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id) , data = dat) #fit linear model with this predictor
  model_withPers <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + Agreeableness*condition + Conscientiousness*condition + OpenMind*condition + NegEmotionality*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woPers, model_withPers)
  
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + Agreeableness*condition + Conscientiousness*condition + OpenMind*condition + NegEmotionality*condition + id + conditionOrder, data = dat, whichRandom = c('id', "conditionOrder"))
  null_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF) # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}




#Thinking Styles
set.seed(22)
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models/Incremental Variance Explained by Psychological Variables")


for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarThinkingStyles.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_woTS <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  model_withTS <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + NeedForCognition*condition + FaithInIntuition*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woTS, model_withTS)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + NeedForCognition*condition + FaithInIntuition*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF) # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models")


#Base model = cognitive variables only
set.seed(22)
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models/Incremental Variance Explained by Psychological Variables")

for (var in dependent){ #do this for each dependent variable you want to test
  pathToFile <- "incVarCognVars.txt" #create text file where significant results will be stored
  if(var == "secondsForBonusThisRound"){
    dat = data[data$condition!="efficiency",]
  }
  else{
    dat = data 
  }
  DV = dat[,var] #select column corresponding to DV in data set
  model_withC <- lmer(DV ~ Numeracy*condition + WorkingMemory*condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  model_woC <- lmer(DV ~ condition + conditionOrder + (1|id), data = dat) #fit linear model with this predictor
  a <- anova(model_woC, model_withC)
  dat$curVar <- dat[,var]
  full_BF = lmBF(curVar ~ Numeracy*condition + WorkingMemory*condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  null_BF = lmBF(curVar ~ condition + id + conditionOrder, data = dat, whichRandom = c('id', 'conditionOrder'))
  BF <- as.data.frame(full_BF / null_BF) # The Bayes factor in favor of the full model
  
  if (a$`Pr(>Chisq)`[2]<=0.05){ 
    significance <- capture.output(a)
    cat(paste(var, "\n", sep = ""), significance, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
    cat(paste("Bayes Factor:", sep = ""), BF$bf, file=pathToFile, sep="\n", append=TRUE)
    cat(paste("Error:", sep = ""), BF$error, file=pathToFile, sep="\n", append=TRUE)
    write("\n\n",file=pathToFile,append=TRUE) #write results (anova, summary) into results file
  }
}

setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience")




# Linear Mixed Models with multiple predictors: only extended model with state emotions was better than base model with WM, Num, cond ----------------------------
#all models are run with bootstrapping
#oddly, we have to specify a very precise model first for the syntax to run,
#bin then in the multiplePredictorModel function we use this generally specified model
#as a template and loop over dependent variables


#setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results")
#dir.create("Multiple Predictor Models")



setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Multiple Predictor Models")


dependent_all <- c("pointsThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first", "EIG0")

set.seed(22)
dependent <- c("pointsThisRound", "guess1jointp_perc_of_maxp")
FileName = "resultsComplexLM_Cognitive+StateEmotion"
  
for (var in dependent){
  DV <- data[,var]
  MixedModel = lmer(DV ~ scale(Numeracy)*condition + scale(WorkingMemory)*condition + scale(DominanceState)*condition + scale(ValenceState)*condition + scale(ArousalState)*condition + (1|id), data = data)
  multiplePredictorModel(var, FileName, MixedModel, data)
  significance <- capture.output(anova(MixedModel))
  results <- capture.output(summary(MixedModel))
  cat(paste(var, "\n", sep = ""), results, file=FileName, sep="\n", append=TRUE)
  cat(paste("Statistical significance", "\n", sep = ""), significance, file=FileName, sep="\n", append=TRUE)
}

var = "secondsForBonusThisRound"
DV <- data[data$condition!="efficiency",var]
MixedModel = lmer(DV ~ scale(Numeracy)*condition + scale(WorkingMemory)*condition + scale(DominanceState)*condition + scale(ValenceState)*condition + scale(ArousalState)*condition + (1|id), data = data[data$condition!="efficiency",])
multiplePredictorModel(var, FileName, MixedModel, data)
significance <- capture.output(anova(MixedModel))
results <- capture.output(summary(MixedModel))
cat(paste(var, "\n", sep = ""), results, file=FileName, sep="\n", append=TRUE)
cat(paste("Statistical significance", "\n", sep = ""), significance, file=FileName, sep="\n", append=TRUE)


#we can do teh same for all other complex models but maybe no need as only 
#extended model with state emotions is better than base model with only num and WM

#for numuguesses & EIG0 (EIG of first guess), the cognitive model seemed best choice

FileName = "resultsComplexLM_Cognitive"
var = "numguesses"
DV <- data[,var]
MixedModel = lmer(DV ~ scale(Numeracy)*condition + scale(WorkingMemory)*condition  + (1|id), data = data)
multiplePredictorModel(var, FileName, MixedModel, data)
significance <- capture.output(anova(MixedModel))
results <- capture.output(summary(MixedModel))
cat(paste(var, "\n", sep = ""), results, file=FileName, sep="\n", append=TRUE)
cat(paste("Statistical significance", "\n", sep = ""), significance, file=FileName, sep="\n", append=TRUE)

var = "EIG0"
DV <- data[,var]
MixedModel = lmer(DV ~ scale(Numeracy)*condition + scale(WorkingMemory)*condition  + (1|id), data = data)
multiplePredictorModel(var, FileName, MixedModel, data)
significance <- capture.output(anova(MixedModel))
results <- capture.output(summary(MixedModel))
cat(paste(var, "\n", sep = ""), results, file=FileName, sep="\n", append=TRUE)
cat(paste("Statistical significance", "\n", sep = ""), significance, file=FileName, sep="\n", append=TRUE)




#for complexity, distance, mean time, and time first guess the base model is sufficient
data$rt_first <- data$rt_first/1000
data$rt_mean_wo_first <- data$rt_mean_wo_first/1000
dependent <- c("complexity_of_query_mean", "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first")
FileName = "resultsComplexLM_Base"

for (var in dependent){
  DV <- data[,var]
  MixedModel = lmer(DV ~ condition + (1|id), data = data)
  multiplePredictorModel(var, FileName, MixedModel, data)
  significance <- capture.output(anova(MixedModel))
  results <- capture.output(summary(MixedModel))
  cat(paste(var, "\n", sep = ""), results, file=FileName, sep="\n", append=TRUE)
  cat(paste("Statistical significance", "\n", sep = ""), significance, file=FileName, sep="\n", append=TRUE)
}

# Supervised Learning: Regression --------------------------------------------------------


#Emotions
set.seed(22)
variables <- data[,c("pointsThisRound","secondsForBonusThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", 
                     "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first", "EIG0",
                     "ValenceTrait", "ValenceState", "DominanceTrait", "DominanceState", 
                     "WorkingMemory", "Numeracy")]
#variables <- data[,c("pointsThisRound", "ValenceState", "DominanceState")]

any(is.na(variables))
nor <- as.data.frame(scale(variables))
nor <- cbind(nor, data$id, data$condition)
colnames(nor)[which(names(nor) =="data$id")] <- "id"
colnames(nor)[which(names(nor) =="data$condition")] <- "condition"



dat_summarized <- nor %>%
        group_by(id) %>%
        summarise(
          count = n(),
          Numeracy = mean(Numeracy),
          WorkingMemory = mean(WorkingMemory),
          pointsThisRound = mean(pointsThisRound),
          secondsForBonusThisRound =mean(secondsForBonusThisRound),
          numguesses = mean(numguesses),
          guess1jointp_perc_of_maxp = mean(guess1jointp_perc_of_maxp),
          complexity_of_query_mean = mean(complexity_of_query_mean),
          rt_first = mean(rt_first),
          distance_to_previous_query_mean = mean(distance_to_previous_query_mean),
          rt_mean_wo_first = mean(rt_mean_wo_first),
          EIG0_mean = mean(EIG0),
          ValenceTrait = mean(ValenceTrait),
          ValenceState = mean(ValenceState),
          DominanceTrait = mean(DominanceTrait),
          DominanceState = mean(DominanceState)
        )


dat_mixed <- nor[nor$condition == "mixed",]


#setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/MLRegression")
#dir.create("ML_Regression")
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/ML_Regression")

#cognitive + state emotion model
dependentVars <- c("pointsThisRound","secondsForBonusThisRound","guess1jointp_perc_of_maxp")
Predictors <-  c("condition*Numeracy", "condition*WorkingMemory","condition*DominanceState","condition*ValenceState")
compute_regressions_Generic(nor, dependentVars, "ML_cognitive+stateEmotio",Predictors)

#cognitive model
dependentVars <- c("numguesses", "EIG0")
Predictors <-  c("condition*Numeracy", "condition*WorkingMemory")
compute_regressions_Generic(nor, dependentVars, "ML_cognitive",Predictors)

#base model
dependentVars <- c("complexity_of_query_mean", "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first")
Predictors <-  c("condition")
compute_regressions_Generic(nor, dependentVars, "ML_base",Predictors)


#all dependent variables with cognitive + emo model
dependentVars <- c("pointsThisRound","secondsForBonusThisRound","numguesses","complexity_of_query_mean", "guess1jointp_perc_of_maxp", 
                     "rt_first", "distance_to_previous_query_mean", "rt_mean_wo_first", "EIG0")

Predictors <-  c("condition*Numeracy", "condition*WorkingMemory","condition*DominanceState","condition*ValenceState")
compute_regressions_Generic(nor, dependentVars, "ML_cognitive+stateEmotio_allDepVars",Predictors)


# Unsupervised Learning: Cluster Analysis ---------------------------------



distance = dist(nor)

wss <- (nrow(nor)-1)*sum(apply(nor,2,var))
for (i in 2:20) {
  wss[i] <- sum(kmeans(nor, centers=i)$withinss)
  plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
}




# Correlation Matrix ------------------------------------------------------
#pretest variables
dat_cor <- data[,c("ValenceTrait", "ValenceState","DominanceTrait", "DominanceState", "WorkingMemory", "Numeracy", "Conscientiousness",
                  "Agreeableness", "NegEmotionality", "OpenMind", "Extraversion", "FaithInIntuition", "NeedForCognition")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("ValTr", "ValSt", "DomTr", "DomSt", "WM", "Num", "Cons", "Agree", "NegEmo", "OpMi", "Extra", "FiI", "NfC")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.4, alpha = 0.01, stars = TRUE)

#game play variables
dat_cor <- data[,c("pointsThisRound", "secondsForBonusThisRound", "numguesses",
                  "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp", "EIG0")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst", "EIG" )

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.6, alpha = 0.01, stars = TRUE)


#overall emotion and game play
dat_cor <- data[,c("ValenceTrait", "ValenceState","DominanceTrait", "DominanceState", "pointsThisRound", "secondsForBonusThisRound", "numguesses",
                   "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp", "EIG0")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("ValTr", "ValSt", "DomTr", "DomSt","points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst", "EIG" )

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2,cex.labels = 1.4, alpha = 0.01, stars = TRUE)


#overall Num, WM and game play
dat_cor <- data[,c("Numeracy", "WorkingMemory","pointsThisRound", "secondsForBonusThisRound", "numguesses",
                   "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp", "EIG0")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("Num", "WM","points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst", "EIG" )

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2,cex.labels = 1.4, alpha = 0.01, stars = TRUE)


#BIS HIER MIT EIG0 GERECHNET


#per condition pretest and game play

#speed, emo + cognitive
dat_cor <- data[dat$condition=="speed",c("ValenceTrait", "ValenceState","DominanceTrait", "DominanceState", "WorkingMemory", "Numeracy", "pointsThisRound", "secondsForBonusThisRound", "numguesses",
                  "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("ValTr", "ValSt", "DomTr", "DomSt", "WM", "Num", "points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)


#efficiency
dat_cor <- dat[dat$condition=="efficiency",c("ValenceTrait", "ValenceState","DominanceTrait", "DominanceState", "WorkingMemory", "Numeracy", "pointsThisRound", "numguesses",
                                        "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("ValTr", "ValSt", "DomTr", "DomSt", "WM", "Num", "points", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)


#mixed

dat_cor <- dat[dat$condition=="mixed",c("ValenceTrait", "ValenceState","DominanceTrait", "DominanceState", "WorkingMemory", "Numeracy", "pointsThisRound","secondsForBonusThisRound", "numguesses",
                                        "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("ValTr", "ValSt", "DomTr", "DomSt", "WM", "Num", "points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)



#thinking styles, personality

#speed, emo + cognitive
dat_cor <- data[dat$condition=="speed",c("Conscientiousness",
                                         "Agreeableness", "NegEmotionality", "OpenMind", "Extraversion", "FaithInIntuition", "NeedForCognition", "pointsThisRound", "secondsForBonusThisRound", "numguesses",
                                         "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("Cons", "Agree", "NegEmo", "OpMi", "Extra", "FiI", "NfC", "points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)


#efficiency
dat_cor <- dat[dat$condition=="efficiency",c("Conscientiousness",
                                             "Agreeableness", "NegEmotionality", "OpenMind", "Extraversion", "FaithInIntuition", "NeedForCognition", "pointsThisRound", "numguesses",
                                             "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("Cons", "Agree", "NegEmo", "OpMi", "Extra", "FiI", "NfC", "points", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)


#mixed

dat_cor <- dat[dat$condition=="mixed",c("Conscientiousness",
                                        "Agreeableness", "NegEmotionality", "OpenMind", "Extraversion", "FaithInIntuition", "NeedForCognition", "pointsThisRound","secondsForBonusThisRound", "numguesses",
                                        "complexity_of_query_mean", "rt_first","rt_mean_wo_first", "distance_to_previous_query_mean", "guess1jointp_perc_of_maxp")]

dat_cor <- sapply(dat_cor, scale)

colnames(dat_cor) <- c("Cons", "Agree", "NegEmo", "OpMi", "Extra", "FiI", "NfC", "points", "time", "nrGuesses", "complexity","rtFirst","rtMean", "queryDist", "jointPr1rst")

CorrMatr <- pairs.panels(dat_cor, density = TRUE, smoother = TRUE, lm = TRUE, scale = FALSE, method = "spearman", ci = TRUE, cex.cor = 1.2, cex.labels = 1.1, alpha = 0.01, stars = TRUE)



# robust and bootstrapped linear models for the significant emotion gameplay relationships --------
#which we actually don't really need as we already bootstrapped all other analyses
#summary statistics
mySumm <- function(.) { s <- sigma(.)
c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }
(t0 <- mySumm(fm01ML)) # just three parameters
## alternatively:
mySumm2 <- function(.) {
  c(beta=fixef(.),sigma=sigma(.), sig01=sqrt(unlist(VarCorr(.))))
}
#points
#robust
linmod <- rlmer(pointsThisRound~condition*DominanceState+condition*ValenceState + (1|id), data = data)
summary(linmod)

#linearity
plot(linmod)

#bootstrapped
lm <- lmer(pointsThisRound~condition*DominanceState+condition*ValenceState + (1|id), data = data)


lm_boot <- bootMer(lm, mySumm2, nsim = 100)

lm_boot_results <- boot_summary(
  lm,
  type = "perc",
  method = "parametric",
  conf.level = 0.95,
  R = 1000,
  coef = "raw",
  adjust.method = "bonferroni",

)

#seconds for bonus calculation
#robust
linmod <- rlmer(secondsForBonusThisRound~condition*DominanceState+condition*ValenceState + (1|id), data = data)
summary(linmod)

#linearity
plot(linmod)

#bootstrapped
lm <- lmer(secondsForBonusThisRound~condition*DominanceState+condition*ValenceState + (1|id), data = data)


lm_boot <- bootMer(lm, mySumm2, nsim = 100)

lm_boot_results <- boot_summary(
  lm,
  type = "perc",
  method = "parametric",
  conf.level = 0.95,
  R = 1000,
  coef = "raw",
  adjust.method = "bonferroni",
  
)

print(lm_boot_results)

#complexity
#robust
linmod <- rlmer(guess1jointp_perc_of_maxp~condition*DominanceState+condition*ValenceState + (1|id), data = data)
summary(linmod)

#linearity
plot(linmod)

#bootstrapped
lm <- lmer(guess1jointp_perc_of_maxp~condition*DominanceState+condition*ValenceState + (1|id), data = data)


lm_boot <- bootMer(lm, mySumm2, nsim = 100)

lm_boot_results <- boot_summary(
  lm,
  type = "perc",
  method = "parametric",
  conf.level = 0.95,
  R = 1000,
  coef = "raw",
  adjust.method = "bonferroni",
  
)

print(lm_boot_results)

# psychologically rational model of subjective information gain -----------

#nr of guesses as outcome variable: fewer queries means more information gained in single query
#predictors: complexity, joint prob, nr smileys, nr frownies, nr neutral -> are these plausible predictors of subjective information gain?
linmod <- lm(numguesses~complexity_of_query_mean*WorkingMemory + complexity_of_query_mean*Numeracy + guess1jointp_perc_of_maxp*Numeracy + guess1jointp_perc_of_maxp*WorkingMemory , data = data)
summary(linmod)

#linearity
plot(linmod, 1) #linear fit seems appropriate
#normally distributed residuals
plot(linmod, 2) #non-normality
shapiro.test(studres(linmod)) #confirmation that residuals are not normally distributed
#high leverage points
plot(linmod, 5)
plot(linmod, 4) #cook's distance

#homoskedasticity
#scale-location plot
plot(linmod, 3)
ncvTest(linmod)#we have homoskedasticity in the data

#as assumptions are violated, we fit robust linear model
#at least
lm_robust <- lmrob(numguesses~complexity_of_query_mean*WorkingMemory + complexity_of_query_mean*Numeracy + guess1jointp_perc_of_maxp*Numeracy + guess1jointp_perc_of_maxp*WorkingMemory , data = dat)
summary(lm_robust)

lm_boot <- ANOVA.boot(data$numguesses~data$complexity_of_query_mean*data$WorkingMemory + data$complexity_of_query_mean*data$Numeracy + data$guess1jointp_perc_of_maxp*data$Numeracy + data$guess1jointp_perc_of_maxp*data$WorkingMemory, B = 10000, type = "residual", wild.dist = "normal", 
                      seed = NULL, data = NULL, keep.boot.resp = FALSE )

lm_boot$`p-values`
lm_boot$terms




# Step by step analysis of game play --------------------------------------
setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Modal Strategies")


#create data frame subsets for first and second guesses + modal strategies (counts of strategies people used)

subset_firstGuess_feedback <- data[, c("guess1", "guess2", "fb1","fb2", "rt1", "rt2", "question_fb_type12")]
subset_secondGuess_feedback <- data[, c("guess2", "guess3", "fb2","fb3", "rt2", "rt3", "question_fb_type23")]
modal_hypotheses12 <- table(subset_firstGuess_feedback$fb1, subset_firstGuess_feedback$question_fb_type12)
modal_hypotheses23 <- table(subset_secondGuess_feedback$fb2, subset_secondGuess_feedback$question_fb_type23)

#significant difference in counts between strategies per feedback type
for (i in 1:(nrow(modal_hypotheses12))){
  one <- modal_hypotheses12[i,]
  chi2test <- chisq.test(one)
  cat(paste(rownames(modal_hypotheses12)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess12", sep="\n", append=TRUE) 
}


for (i in 1:(nrow(modal_hypotheses23))){
  one <- modal_hypotheses23[i,]
  chi2test <- chisq.test(one)
  cat(paste(rownames(modal_hypotheses23)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess23", sep="\n", append=TRUE) 
}


#differences in strategies between 12 and 23

for (i in 1:(nrow(modal_hypotheses12))){
  one <- modal_hypotheses12[i,]
  two <- modal_hypotheses23[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c, simulate.p.value = TRUE)
  cat(paste(rownames(modal_hypotheses12)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess12vs23", sep="\n", append=TRUE) 
}

one <- modal_hypotheses12["021",]
two <- modal_hypotheses23["021",]
c <- rbind(one, two)
colnames(c) <- colnames(modal_hypotheses12)
chisq.test(c[,-1])
chi2test <- fisher.test(c, simulate.p.value = TRUE)


subset_firstGuess_feedback_speed <- data[data$condition=="speed", c("guess1", "guess2", "fb1","fb2", "rt1", "rt2", "question_fb_type12")]
subset_secondGuess_feedback_speed <- data[data$condition=="speed", c("guess2", "guess3", "fb2","fb3", "rt2", "rt3", "question_fb_type23")]
modal_hypotheses12_speed <- table(subset_firstGuess_feedback_speed$fb1, subset_firstGuess_feedback_speed$question_fb_type12)
modal_hypotheses23_speed <- table(subset_secondGuess_feedback_speed$fb2, subset_secondGuess_feedback_speed$question_fb_type23)

subset_firstGuess_feedback_efficiency <- data[data$condition=="efficiency", c("guess1", "guess2", "fb1","fb2", "rt1", "rt2", "question_fb_type12")]
subset_secondGuess_feedback_efficiency <- data[data$condition=="efficiency", c("guess2", "guess3", "fb2","fb3", "rt2", "rt3", "question_fb_type23")]
modal_hypotheses12_efficiency <- table(subset_firstGuess_feedback_efficiency$fb1, subset_firstGuess_feedback_efficiency$question_fb_type12)
modal_hypotheses23_efficiency <- table(subset_secondGuess_feedback_efficiency$fb2, subset_secondGuess_feedback_efficiency$question_fb_type23)

subset_firstGuess_feedback_mixed <- data[data$condition=="mixed", c("guess1", "guess2", "fb1","fb2", "rt1", "rt2", "question_fb_type12")]
subset_secondGuess_feedback_mixed <- data[data$condition=="mixed", c("guess2", "guess3", "fb2","fb3", "rt2", "rt3", "question_fb_type23")]
modal_hypotheses12_mixed <- table(subset_firstGuess_feedback_mixed$fb1, subset_firstGuess_feedback_mixed$question_fb_type12)
modal_hypotheses23_mixed <- table(subset_secondGuess_feedback_mixed$fb2, subset_secondGuess_feedback_mixed$question_fb_type23)

#significant differences in strategies per condition
#speed vs efficiency
for (i in 1:(nrow(modal_hypotheses12))){
  one <- modal_hypotheses12_speed[i,]
  two <- modal_hypotheses12_efficiency[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses12)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess12_speedvsefficiency", sep="\n", append=TRUE) 
}


for (i in 1:(nrow(modal_hypotheses23))){
  one <- modal_hypotheses23_speed[i,]
  two <- modal_hypotheses23_efficiency[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses23)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess23_speedvsefficiency", sep="\n", append=TRUE) 
}


#speed vs mixed
for (i in 1:(nrow(modal_hypotheses12))){
  one <- modal_hypotheses12_speed[i,]
  two <- modal_hypotheses12_mixed[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses12)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess12_speedvsmixed", sep="\n", append=TRUE) 
}


for (i in 1:(nrow(modal_hypotheses23))){
  one <- modal_hypotheses23_speed[i,]
  two <- modal_hypotheses23_mixed[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses23)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess23_speedvsmixed", sep="\n", append=TRUE) 
}


#speed vs mixed
for (i in 1:(nrow(modal_hypotheses12))){
  one <- modal_hypotheses12_mixed[i,]
  two <- modal_hypotheses12_efficiency[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses12)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess12_mixedvsefficiency", sep="\n", append=TRUE) 
}


for (i in 1:(nrow(modal_hypotheses23))){
  one <- modal_hypotheses23_mixed[i,]
  two <- modal_hypotheses23_efficiency[i,]
  c <- rbind(one, two)
  colnames(c) <- colnames(modal_hypotheses12)
  chi2test <- fisher.test(c)
  cat(paste(rownames(modal_hypotheses23)[i], "\n", sep = ""), capture.output(chi2test), file="chi2_feedback_guess23_mixedvsefficiency", sep="\n", append=TRUE) 
}

#setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results")

#dir.create("Modal Strategies")

setwd("C:/Users/larab/OneDrive - University of Surrey/LaraiCloud/Promotion/Publications/PsychScience/Results/Modal Strategies")

write.csv(modal_hypotheses12, "modal_strategies_12.csv")
write.csv(modal_hypotheses23, "modal_strategies_23.csv")

write.csv(modal_hypotheses12_speed, "modal_strategies_12_speed.csv")
write.csv(modal_hypotheses23_speed, "modal_strategies_23_speed.csv")

write.csv(modal_hypotheses12_efficiency, "modal_strategies_12_efficiency.csv")
write.csv(modal_hypotheses23_efficiency, "modal_strategies_23_efficiency.csv")


write.csv(modal_hypotheses12_mixed, "modal_strategies_12_mixed.csv")
write.csv(modal_hypotheses23_mixed, "modal_strategies_23_mixed.csv")
#normalized counts for modal strategies

#load modal strategies for first + second guess

m12 <- read.csv2("modal_strategies_12.csv", na.strings = "", stringsAsFactors = FALSE, header=T, sep = ",")
rownames(m12) <- m12[,1]
m12 <- m12[,-1]
m12$total_fbType <- rowSums(m12) #add rowsums for normalization

#normalize counts
for (i in 1:nrow(m12)){
  for (j in 1:(length(m12)-1)){
    m12[i, j] <- m12[i, j]/m12[i, "total_fbType"]*100
  }
}


#load modal strategies for second + third guess
m23 <- read.csv2("modal_strategies_23.csv", na.strings = "", stringsAsFactors = FALSE, header=T, sep = ",")
rownames(m23) <- m23[,1]
m23 <- m23[,-1]
m23$total_fbType <- rowSums(m23)

#normalize counts
for (i in 1:nrow(m23)){
  for (j in 1:(length(m23)-1)){
    m23[i, j] <- m23[i, j]/m23[i, "total_fbType"]*100
  }
}


write.csv(m12, "modal_strategies_12_normalized.csv")
write.csv(m23, "modal_strategies_23_normalized.csv")


# guess + feedback combinations by pretest variables ----------------

#guess 12
#select relevant pretest variables and create data frame with means by guess-feedback combinations
num_by_guessFbCombi <- data %>%
  group_by(fb1, question_fb_type12) %>%
  summarise(mean(Numeracy), mean(WorkingMemory), mean(DominanceState), mean(ValenceState))

num_by_guessFbCombi <-  data.frame(num_by_guessFbCombi)

#create data frame for guess1-feedback-guess2 combinations
n <- data.frame(matrix(nrow = length(unique(num_by_guessFbCombi$fb1)), ncol = length(unique(num_by_guessFbCombi$question_fb_type12))))
rownames(n) <- unique(num_by_guessFbCombi$fb1)
colnames(n) <- unique(num_by_guessFbCombi$question_fb_type12)
n <- n[,-10]
wm <- n
ds <- n
vs <- n

for (i in unique(num_by_guessFbCombi$fb1)){
  for (j in unique(num_by_guessFbCombi$question_fb_type12)[-length(unique(num_by_guessFbCombi$question_fb_type12))]){
    nm <- num_by_guessFbCombi$mean.Numeracy.[which(num_by_guessFbCombi$fb1 == i & num_by_guessFbCombi$question_fb_type12==j)]
    if(length(nm)>0){
      n[i, j] <- nm
    }
    wmm <- num_by_guessFbCombi$mean.WorkingMemory.[which(num_by_guessFbCombi$fb1 == i & num_by_guessFbCombi$question_fb_type12==j)]
    if(length(wmm)>0){
      wm[i, j] <- wmm
    }
    dsm <- num_by_guessFbCombi$mean.DominanceState.[which(num_by_guessFbCombi$fb1 == i & num_by_guessFbCombi$question_fb_type12==j)]
    if(length(dsm)>0){
      ds[i, j] <- dsm
    }
    vsm <- num_by_guessFbCombi$mean.ValenceState.[which(num_by_guessFbCombi$fb1 == i & num_by_guessFbCombi$question_fb_type12==j)]
    if(length(vsm)>0){
      vs[i, j] <- vsm
    }
  
  }
}


write.csv(n, "numeracy_modalStrategies12.csv")
write.csv(wm, "workingMemory_modalStrategies12.csv")
write.csv(ds, "dominanceState_modalStrategies12.csv")
write.csv(vs, "valenceState_modalStrategies12.csv")




#select relevant pretest variables and create data frame with means by guess-feedback combinations


num_by_guessFbCombi <- data %>%
  group_by(fb2, question_fb_type23) %>%
  summarise(mean(Numeracy), mean(WorkingMemory), mean(DominanceState), mean(ValenceState))

num_by_guessFbCombi <-  data.frame(num_by_guessFbCombi)


n <- data.frame(matrix(nrow = (length(unique(num_by_guessFbCombi$fb2))-1), ncol = (length(unique(num_by_guessFbCombi$question_fb_type23))-1)))
rownames(n) <- unique(num_by_guessFbCombi$fb2)[-length(unique(num_by_guessFbCombi$fb2))]
colnames(n) <- unique(num_by_guessFbCombi$question_fb_type23)[-9]

wm <- n
ds <- n
vs <- n

for (i in unique(num_by_guessFbCombi$fb2)[-length(unique(num_by_guessFbCombi$fb2))]){
  for (j in unique(num_by_guessFbCombi$question_fb_type23)[-9]){
    nm <- num_by_guessFbCombi$mean.Numeracy.[which(num_by_guessFbCombi$fb2 == i & num_by_guessFbCombi$question_fb_type23==j)]
    if(length(nm)>0){
      n[i, j] <- nm
    }
    wmm <- num_by_guessFbCombi$mean.WorkingMemory.[which(num_by_guessFbCombi$fb2 == i & num_by_guessFbCombi$question_fb_type23==j)]
    if(length(wmm)>0){
      wm[i, j] <- wmm
    }
    dsm <- num_by_guessFbCombi$mean.DominanceState.[which(num_by_guessFbCombi$fb2 == i & num_by_guessFbCombi$question_fb_type23==j)]
    if(length(dsm)>0){
      ds[i, j] <- dsm
    }
    vsm <- num_by_guessFbCombi$mean.ValenceState.[which(num_by_guessFbCombi$fb2 == i & num_by_guessFbCombi$question_fb_type23==j)]
    if(length(vsm)>0){
      vs[i, j] <- vsm
    }
  }
}


write.csv(n, "numeracy_modalStrategies23.csv")
write.csv(wm, "workingMemory_modalStrategies23.csv")
write.csv(ds, "dominanceState_modalStrategies23.csv")
write.csv(vs, "valenceState_modalStrategies23.csv")


#chi2 test for numeracy and other pretest variables

#no chi2 test possible

#does numeracy differ between strategies for each feedback type: guess 12



#numeracy 
fbs <- unique(data$fb1[!is.na(data$fb1) & data$fb1 != "300"])

filename <- "Numeracy_by_Strategy"

for (fb in fbs){
  d <- data[data$fb1 == fb,]
  l_m <- lm(d$Numeracy~as.factor(d$question_fb_type12), data = d)
  a <- anova(l_m)
  s <- summary(l_m)
  cat(paste(fb, "\n", sep = ""), capture.output(a), file=filename, sep="\n", append=TRUE) 
  cat(paste("summary", "\n", sep = ""), capture.output(s), file=filename, sep="\n", append=TRUE) 
}

d <- data[!is.na(data$fb1[!is.na(data$fb1) & data$fb1 != "300"]),]
o_lm <- lm(d$Numeracy~as.factor(d$fb1)*as.factor(d$question_fb_type12), data = d)
summary(o_lm)
anova(o_lm)
aov(o_lm)

filename <- "WorkingMemory_by_Strategy"

for (fb in fbs){
  d <- data[data$fb1 == fb,]
  l_m <- lm(d$WorkingMemory~as.factor(d$question_fb_type12), data = d)
  a <- anova(l_m)
  s <- summary(l_m)
  cat(paste(fb, "\n", sep = ""), capture.output(a), file=filename, sep="\n", append=TRUE) 
  cat(paste("summary", "\n", sep = ""), capture.output(s), file=filename, sep="\n", append=TRUE) 
}

d <- data[!is.na(data$fb1[!is.na(data$question_fb_type12) & data$fb1 != "300"]),]
o_lm <- lm(d$WorkingMemory~d$fb1*d$question_fb_type12, data = d)
summary(o_lm)
anova(o_lm)


filename <- "Dominance_by_Strategy"

for (fb in fbs){
  d <- data[data$fb1 == fb,]
  l_m <- lm(d$DominanceState~as.factor(d$question_fb_type12), data = d)
  a <- anova(l_m)
  s <- summary(l_m)
  cat(paste(fb, "\n", sep = ""), capture.output(a), file=filename, sep="\n", append=TRUE) 
  cat(paste("summary", "\n", sep = ""), capture.output(s), file=filename, sep="\n", append=TRUE) 
}

d <- data[!is.na(data$fb1[!is.na(data$fb1) & data$fb1 != "300"]),]
o_lm <- lm(d$DominanceState~as.factor(d$fb1)*as.factor(d$question_fb_type12), data = d)
summary(o_lm)
anova(o_lm)

filename <- "Valence_by_Strategy"

for (fb in fbs){
  d <- data[data$fb1 == fb,]
  l_m <- lm(d$ValenceState~as.factor(d$question_fb_type12), data = d)
  a <- anova(l_m)
  s <- summary(l_m)
  cat(paste(fb, "\n", sep = ""), capture.output(a), file=filename, sep="\n", append=TRUE) 
  cat(paste("summary", "\n", sep = ""), capture.output(s), file=filename, sep="\n", append=TRUE) 
}

d <- data[!is.na(data$fb1[!is.na(data$fb1) & data$fb1 != "300"]),]
o_lm <- lm(d$ValenceState~as.factor(d$fb1)+as.factor(d$question_fb_type12), data = d)
summary(o_lm)
anova(o_lm)


#M and SD for significantly different groups
#003 feedback
mean(data$Numeracy[data$fb1 == "003" & data$question_fb_type12 == "C3L0" ])
mean(data$Numeracy[data$fb1 == "003" & data$question_fb_type12 != "C3L0" & !is.na(data$Numeracy)])

sd(data$Numeracy[data$fb1 == "003" & data$question_fb_type12 == "C3L0" ])
sd(data$Numeracy[data$fb1 == "003" & data$question_fb_type12 != "C3L0" & !is.na(data$Numeracy) ])


#021 feedback
mean(data$Numeracy[data$fb1 == "021" & data$question_fb_type12 == "C3L0" ])
mean(data$Numeracy[data$fb1 == "021" & data$question_fb_type12 != "C3L0" & !is.na(data$Numeracy)])

sd(data$Numeracy[data$fb1 == "021" & data$question_fb_type12 == "C3L0" ])
sd(data$Numeracy[data$fb1 == "021" & data$question_fb_type12 != "C3L0"& !is.na(data$Numeracy)] )


#M and SD for significantly different groups: Working Memory
#030 feedback
mean(data$WorkingMemory[!is.na(data$WorkingMemory[data$fb1 == "030" & data$question_fb_type12 == "C0L3"])])
mean(data$WorkingMemory[!is.na(data$WorkingMemory[data$fb1 == "030" & data$question_fb_type12 != "C0L3"])])

sd(!is.na(data$WorkingMemory[data$fb1 == "030" & data$question_fb_type12 == "C0L3" ]))
sd(!is.na(data$WorkingMemory[data$fb1 == "030" & data$question_fb_type12 != "C0L3" ]))


#021 feedback
mean(data$WorkingMemory[data$fb1 == "021" & data$question_fb_type12 == "C3L0" ])
mean(data$WorkingMemory[data$fb1 == "021" & data$question_fb_type12 != "C3L0" & !is.na(data$Numeracy)])

sd(data$WorkingMemory[data$fb1 == "021" & data$question_fb_type12 == "C3L0" ])
sd(data$WorkingMemory[data$fb1 == "021" & data$question_fb_type12 != "C3L0"& !is.na(data$Numeracy)] )




#multinomial logistic regression: predicting strategy per feedback type from pretest variables

install.packages("nnet")
library(nnet)

d <- data[data$fb1[!is.na(data$fb1) & data$fb1 == "003" ],]
m <- multinom(d$question_fb_type12 ~ d$Numeracy, data = d)
summary(m)

o_lm <- lm(d$WorkingMemory~as.factor(d$fb1)*as.factor(d$question_fb_type12), data = d)
summary(o_lm)
anova(o_lm)








#more positiv efeedback


dataMeans <- ddply(data, .(id), summarise, points = mean(pointsThisRound), secPoints = mean(secondsForBonusThisRound), numGuesses = mean(numguesses),
                   rtFirst = mean(rt_first), rtMean = mean(rt_mean), jpFirst = mean(guess1jointp_perc_of_maxp), complexityMean = mean(complexity_of_query_mean),
                   distanceMean = mean(distance_to_previous_query_mean), EIGFQMean = mean(EIG0),ValT = mean(ValenceTrait), ValS = mean(ValenceState), 
                   DomT = mean(DominanceTrait), DomS = mean(DominanceState), Numeracy = mean(Numeracy), WorkingMemory = mean(WorkingMemory))

plot_mixedModel2Cat_withPoints <- function(data, IVwithin, DV, IVbetween, title, labelY, labelX, minValue, maxValue, colors, breaks){
  data %>% 
    ggplot() +
    aes(x = IVbetween, color = IVwithin, shape = IVwithin, y = DV, group = IVwithin, linetype = IVwithin) +
    geom_point(aes (color = as.factor(IVwithin)), position = position_jitterdodge(), alpha = 0.2)+ 
    stat_summary(fun = mean, geom = "point", size = 2.5, position = position_dodge(0.4)) + 
    stat_summary(fun = mean, geom = "line", position = position_dodge(0.4)) + theme_classic() +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.4, position = position_dodge(0.4), inherit.aes = FALSE, aes(x = IVbetween, color = IVwithin, shape = IVwithin, y = DV, group = IVwithin))+
    scale_color_manual(name = title, values=colors, breaks = breaks)+ scale_shape_manual(name = title, values = c(17, 15), breaks = breaks)+
    scale_linetype_manual(name = title, values=c("dotted", "solid"), breaks = breaks)+
    ylab(labelY) + scale_y_continuous(limits = c(minValue, maxValue))+ xlab(labelX)+
    theme(legend.position = "top", legend.title = element_text(size = 25), legend.text = element_text(size = 20), axis.title.y = element_text(size = 30, color = "gray12"), axis.text.y = element_text(size = 25, color = "gray12"), axis.title.x = element_text(size = 30, color = "gray12"), axis.text.x = element_text(size = 25, color = "gray12")) 
}

data_woEff <- data[data$condition!="efficiency",]

#First Query strategies

#Histograms

dat$condition[dat$condition == "efficiency"] <- "accuracy"


plotEIG = ggplot(dat, aes(x = as.factor(condition), y = EIG0, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+
  labs(x="", y = "EIG FQ") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotEIG

#dat$rt_first_sec <- dat$rt_first/1000

plotrtFQ = ggplot(dat, aes(x = as.factor(condition), y = rt_first_sec, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+ ylim(c(0, 60))+
  labs(x="", y = "time FQ (sec)") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotrtFQ


plotjpFQ = ggplot(dat, aes(x = as.factor(condition), y = guesses_joint_prob_mean_perc_of_maxp, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+ ylim(c(0, 60))+
  labs(x="", y = "jp FQ") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotjpFQ

hist_fq <- grid.arrange(plotEIG, plotjpFQ, plotrtFQ,  nrow = 1)


#EIG

data$condition[data$condition == "efficiency"] <- "accuracy"

Numeracy_EIG <- plot_mixedModel2Cat_withPoints(data, data$NumCat, data$EIG0, data$condition, "Numeracy", "EIG FQ","", min(data$EIG0), max(data$EIG0), c("#0072B2", "#D55E00"), c("low", "high"))
Numeracy_EIG


WM_EIG <- plot_mixedModel2Cat_withPoints(data, data$WMCat, data$EIG0, data$condition, "Working Memory", "EIG FQ","", min(data$EIG0), max(data$EIG0), c("#0072B2", "#D55E00"), c("low", "high"))
WM_EIG

int_EIG <- grid.arrange(Numeracy_EIG, WM_EIG, nrow = 1)

baseCogn <- grid.arrange(hist_fq, int_EIG, nrow = 2)

#JP FQ

DomS_JP <- plot_mixedModel2Cat_withPoints(data, data$DominanceStateCat, data$guess1jointp_perc_of_maxp, data$condition, "Dominance State", "joint prob FG","", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
DomS_JP


ValS_JP <- plot_mixedModel2Cat_withPoints(data, data$ValenceStateCat, data$guess1jointp_perc_of_maxp, data$condition, "Valence State", "joint prob FG","", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
ValS_JP


Numeracy_JP <- plot_mixedModel2Cat_withPoints(data, data$NumCat, data$guess1jointp_perc_of_maxp, data$condition, "Numeracy", "joint prob FG","", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
Numeracy_JP


WM_JP <- plot_mixedModel2Cat_withPoints(data, data$WMCat, data$guess1jointp_perc_of_maxp, data$condition, "Working Memory", "joint prob FG","", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
WM_JP


emoCog <- grid.arrange(Numeracy_JP, WM_JP, DomS_JP, ValS_JP, nrow = 2)

all <- grid.arrange(baseCogn, emoCog, nrow = 1)

#cognitive + state emotion model for points, seconds for bonus, joint probability

#seconds: main effects of numeracy, working memory and interactions between condition and numeracy, working memory, state dominance, state valence

#Main effects
#Base model

plotqcMean = ggplot(dat, aes(x = as.factor(condition), y = complexity_of_query_mean, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+ 
  labs(x="", y = "mean query complexity") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotqcMean



plotdMean = ggplot(dat, aes(x = as.factor(condition), y = dat$distance_to_previous_query_mean, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+ 
  labs(x="", y = "mean query distance") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotdMean

dat$rt_mean_wo_first_sec <- dat$rt_mean_wo_first/1000
plotrtMean = ggplot(dat, aes(x = as.factor(condition), y = rt_mean_wo_first_sec, fill = as.factor(condition))) + geom_boxplot() + theme_classic()+
  #geom_point(alpha = 0.2)+
  stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2)+ ylim(c(0,60))+
  labs(x="", y = "mean time on query (sec)") + scale_fill_brewer(palette="PuBu") + theme(legend.position='none', axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 20))

plotrtMean


base <- grid.arrange(plotqcMean, plotdMean, plotrtMean, nrow = 1)

Base_all <- grid.arrange(hist_fq, base, nrow = 2)


#cognitive

Numeracy_secs_main <- ggplot(dataMeans, aes(x = Numeracy, y = secPoints, color = Numeracy)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Numeracy", y = "seconds")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$secPoints), max(dataMeans$secPoints)))

Numeracy_secs_main

WorkingMemory_secs_main <- ggplot(dataMeans, aes(x = WorkingMemory, y = secPoints, color = WorkingMemory)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Working Memory", y = "seconds")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$secPoints), max(dataMeans$secPoints)))

WorkingMemory_secs_main

#Interactions
DomS_secs <- plot_mixedModel2Cat_withPoints(data_woEff, data_woEff$DominanceStateCat, data_woEff$secondsForBonusThisRound, data_woEff$condition, "Dominance (state)", "sec","", min(data_woEff$secondsForBonusThisRound), max(data_woEff$secondsForBonusThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
DomS_secs

ValS_secs <- plot_mixedModel2Cat_withPoints(data_woEff, data_woEff$ValenceStateCat, data_woEff$secondsForBonusThisRound, data_woEff$condition, "Valence (state)", "sec","", min(data_woEff$secondsForBonusThisRound), max(data_woEff$secondsForBonusThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
ValS_secs

Numeracy_secs <- plot_mixedModel2Cat_withPoints(data_woEff, data_woEff$NumCat, data_woEff$secondsForBonusThisRound, data_woEff$condition, "Numeracy", "sec","", min(data_woEff$secondsForBonusThisRound), max(data_woEff$secondsForBonusThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
Numeracy_secs

WM_secs <- plot_mixedModel2Cat_withPoints(data_woEff, data_woEff$WMCat, data_woEff$secondsForBonusThisRound, data_woEff$condition, "Working Memory", "sec","", min(data_woEff$secondsForBonusThisRound), max(data_woEff$secondsForBonusThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
WM_secs

#points: main effects of numeracy, working memory and interactions between condition and numeracy, working memory, state dominance, state valence
#main effects:

Numeracy_points_main <- ggplot(dataMeans, aes(x = Numeracy, y = points, color = Numeracy)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Numeracy", y = "points")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$points), max(dataMeans$points)))

Numeracy_points_main

WorkingMemory_points_main <- ggplot(dataMeans, aes(x = WorkingMemory, y = points, color = WorkingMemory)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Working Memory", y = "points")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$points), max(dataMeans$points)))

WorkingMemory_points_main

#interactions

DomS_points <- plot_mixedModel2Cat_withPoints(data, data$DominanceStateCat, data$pointsThisRound, data$condition, "Dominance State", "points","", min(data$pointsThisRound), max(data$pointsThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
DomS_points


ValS_points <- plot_mixedModel2Cat_withPoints(data, data$ValenceStateCat, data$pointsThisRound, data$condition, "Valence State", "points","", min(data$pointsThisRound), max(data$pointsThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
ValS_points


Numeracy_points <- plot_mixedModel2Cat_withPoints(data, data$NumCat, data$pointsThisRound, data$condition, "Numeracy", "points","", min(data$pointsThisRound), max(data$pointsThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
Numeracy_points


WM_points <- plot_mixedModel2Cat_withPoints(data, data$WMCat, data$pointsThisRound, data$condition, "Working Memory", "points","", min(data$pointsThisRound), max(data$pointsThisRound), c("#0072B2", "#D55E00"), c("low", "high"))
WM_points


Numeracy_numGuess <- plot_mixedModel2Cat_withPoints(data, data$NumCat, data$numguesses, data$condition, "Numeracy", "nr guesses","", min(data$numguesses), max(data$numguesses), c("#0072B2", "#D55E00"), c("low", "high"))
Numeracy_numGuess


WM_numGuess <- plot_mixedModel2Cat_withPoints(data, data$WMCat, data$numguesses, data$condition, "Working Memory", "nr guesses","", min(data$numguesses), max(data$numguesses), c("#0072B2", "#D55E00"), c("low", "high"))
WM_numGuess



#Joint Probability FG: main effects of numeracy and state valence; interaction between condition and state dominance and working memory
#main effects:

Numeracy_JP_main <- ggplot(dataMeans, aes(x = Numeracy, y = jpFirst, color = Numeracy)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Numeracy", y = "joint prob FG")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$jpFirst), max(dataMeans$jpFirst)))

Numeracy_JP_main


Valence_JP_main <- ggplot(dataMeans, aes(x = ValS, y = jpFirst, color = ValS)) + geom_point(position = position_jitterdodge(), alpha = 0.3) + theme_classic()+
  geom_smooth(method=lm, color = "blue") + labs(x="Valence (state)", y = "joint prob FG")+
  theme(legend.position='none', axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
  geom_density_2d(color = "blue") + scale_y_continuous(limits = c(min(dataMeans$jpFirst), max(dataMeans$jpFirst)))

Valence_JP_main


#Interactions

DomS_JP <- plot_mixedModel2Cat_withPoints(data, data$DominanceStateCat, data$guess1jointp_perc_of_maxp, data$condition, "Dominance State", "joint prob FG","condition", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
DomS_JP


#ValS_JP <- plot_mixedModel2Cat_withPoints(data, data$ValenceStateCat, data$guess1jointp_perc_of_maxp, data$condition, "Valence State", "joint prob FG","condition", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
#ValS_JP


#Numeracy_JP <- plot_mixedModel2Cat_withPoints(data, data$NumCat, data$guess1jointp_perc_of_maxp, data$condition, "Numeracy", "joint prob FG","condition", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
#Numeracy_JP


WM_JP <- plot_mixedModel2Cat_withPoints(data, data$WMCat, data$guess1jointp_perc_of_maxp, data$condition, "Working Memory", "joint prob FG","condition", min(data$guess1jointp_perc_of_maxp), max(data$guess1jointp_perc_of_maxp), c("#0072B2", "#D55E00"), c("low", "high"))
WM_JP

#numguesses: cognitive model


#main effects:

grid.arrange(Numeracy_secs_main, WorkingMemory_secs_main, Numeracy_points_main, WorkingMemory_points_main, 
             Numeracy_JP_main, Valence_JP_main, nrow = 3)

grid.arrange(DomS_points, ValS_points, Numeracy_points, WM_points, nrow = 2)
grid.arrange(DomS_secs, ValS_secs, Numeracy_secs, WM_secs, nrow = 2)
grid.arrange(DomS_JP, WM_JP, nrow = 1)

grid.arrange(Numeracy_numGuess, WM_numGuess, nrow = 1)
#game play variables


#check mean time spent on games, figure out how to weight conditions

#game play (learning) over time by pretest variables