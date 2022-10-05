#######################
### Small2 scenario ###
#######################

library(DynForest) # Install with DynForest_devlopment_version/DynForest_0.1.0.tar.gz
library(JMbayes)
library(doParallel)
library(survival)
library(rstpm2)
library(splines)
library(lcmm)
library(pec)
library(timeROC)

source("simu_function.R")

load("marker2_splines.RData")

evt.link <- list(model = ~ I(marker1_class == 1) + I(marker1_class == 2) + I(marker1_class == 3) + I(marker1_class == 4),
                 gamma = c(-2.5, -1, 1.5, 3))

set.seed(1000)

testData <- simuData(n = 500, marker = marker, class.lat = 4,
                     evt.link = evt.link, Y.type = "surv",
                     censoring.value = exp(-2.5),
                     var.timevis = TRUE,
                     truncation = TRUE, trunctime = 10,
                     baseline.surv.param = c(0.1, 2), baseline.surv.fct = "Weibull",
                     scale_marker = FALSE,
                     true.prob = TRUE, tLMs = c(2,4), tHors = c(1,2))

testData_long <- testData$data.long
testData_surv <- subset(testData$data.surv, select = c(id, cont_covar1, bin_covar1, time_event, evt))

# starting values mixed models

marker1.best <- hlme(fixed = marker1 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                     random = ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                     subject = "id", data = testData_long)$best

marker2.best <- hlme(fixed = marker2 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                     random = ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                     subject = "id", data = testData_long)$best

#######################
### Parallelisation ###
#######################

replicates.min <- 1
replicates.max <- 250
ncores_fct <- 12

res <- list()

for (replicate in replicates.min:replicates.max){

  cat(paste0(replicate,"/",replicates.max),"\n")

  cat("Simu data train","\n")

  set.seed(replicate)

  trainData <- simuData(n = 500, marker = marker, class.lat = 4,
                        evt.link = evt.link, Y.type = "surv",
                        censoring.value = exp(-2.5),
                        var.timevis = TRUE,
                        truncation = TRUE, trunctime = 10,
                        baseline.surv.param = c(0.1, 2), baseline.surv.fct = "Weibull",
                        scale_marker = FALSE)

  trainData_long <- trainData$data.long
  trainData_surv <- subset(trainData$data.surv, select = c(id, cont_covar1, bin_covar1, time_event, evt))

  ####################
  #### Estimation ####
  ####################

  # DynForest

  Y <- list(type = "surv", Y = Surv(trainData_surv$time_event, trainData_surv$evt), id = trainData_surv$id)
  scalar <- list(X = trainData_surv[,c("cont_covar1"), drop = FALSE], id = trainData_surv$id)
  facto <- list(X = trainData_surv[,c("bin_covar1"), drop = FALSE], id = trainData_surv$id)
  curve <- list(X = trainData_long[,c("marker1","marker2"), drop = FALSE], id = trainData_long$id,
                time = trainData_long$time,
                model = list(marker1 = list(fixed = marker1 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                                            random = ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                                            init.param = marker1.best),
                             marker2 = list(fixed = marker2 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                                            random = ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)),
                                            init.param = marker2.best)))

  # mtry 1

  set.seed(replicate)

  cat("Estimation DynForest","\n")

  res_dyn1 <- tryCatch(DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                 Y=Y, ntree = 250, imp = TRUE,
                                 mtry = 1, nodesize = 3, minsplit = 5,
                                 ncores = ncores_fct),
                       error = function(e) return(NULL))

  if (!is.null(res_dyn1)){
    dyn_depth1 <- var_depth(res_dyn1)
    vimp_dyn1 <- list(res_dyn1$Importance, res_dyn1$Inputs)
  }else{
    dyn_depth1 <- NULL
    vimp_dyn1 <- NULL
  }



  # mtry 2

  set.seed(replicate)

  cat("Estimation DynForest","\n")

  res_dyn2 <- tryCatch(DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                 Y=Y, ntree = 250, imp = TRUE,
                                 mtry = 2, nodesize = 3, minsplit = 5,
                                 ncores = ncores_fct),
                       error = function(e) return(NULL))

  if (!is.null(res_dyn2)){
    dyn_depth2 <- var_depth(res_dyn2)
    vimp_dyn2 <- list(res_dyn2$Importance, res_dyn2$Inputs)
  }else{
    dyn_depth2 <- NULL
    vimp_dyn2 <- NULL
  }

  # mtry 3

  set.seed(replicate)

  cat("Estimation DynForest","\n")

  res_dyn3 <- tryCatch(DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                 Y=Y, ntree = 250, imp = TRUE,
                                 mtry = 3, nodesize = 3, minsplit = 5,
                                 ncores = ncores_fct),
                       error = function(e) return(NULL))

  if (!is.null(res_dyn1)){
    dyn_depth3 <- var_depth(res_dyn3)
    vimp_dyn3 <- list(res_dyn3$Importance, res_dyn3$Inputs)
  }else{
    dyn_depth3 <- NULL
    vimp_dyn3 <- NULL
  }

  # mtry 4

  set.seed(replicate)

  cat("Estimation DynForest","\n")

  res_dyn4 <- tryCatch(DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                 Y=Y, ntree = 250, imp = TRUE,
                                 mtry = 4, nodesize = 3, minsplit = 5,
                                 ncores = ncores_fct),
                       error = function(e) return(NULL))

  if (!is.null(res_dyn4)){
    dyn_depth4 <- var_depth(res_dyn4)
    vimp_dyn4 <- list(res_dyn4$Importance, res_dyn4$Inputs)
  }else{
    dyn_depth4 <- NULL
    vimp_dyn4 <- NULL
  }

  # JMbayes intercept + slope

  cat("Estimation JMbayes","\n")

  mm_fit <- mvglmer(list(marker1 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)) +
                           (splines::ns(time, knots = 5, Boundary.knots = c(0,10))|id),
                         marker2 ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10)) +
                           (splines::ns(time, knots = 5, Boundary.knots = c(0,10))|id)),
                    data = trainData_long,
                    families = list(gaussian, gaussian),
                    control = list(n.processors = ncores_fct))

  cat("Cox...","\n")

  cox_fit <- coxph(Surv(time_event, evt) ~ cont_covar1 + bin_covar1,
                   data = trainData_surv, model = TRUE)

  Forms_intslope <- list("marker1" = "value",
                         "marker1" = list(fixed = ~ -1 + rstpm2::nsxD(time, knots = 5, Boundary.knots = c(0,10)),
                                          random = ~ -1 + rstpm2::nsxD(time, knots = 5, Boundary.knots = c(0,10)),
                                          indFixed = c(2,3),
                                          indRandom = c(2,3), name = "slope"),
                         "marker2" = "value",
                         "marker2" = list(fixed = ~ -1 + rstpm2::nsxD(time, knots = 5, Boundary.knots = c(0,10)),
                                          random = ~ -1 + rstpm2::nsxD(time, knots = 5, Boundary.knots = c(0,10)),
                                          indFixed = c(2,3),
                                          indRandom = c(2,3), name = "slope"))

  set.seed(replicate)

  JM_fit_intslope <- mvJointModelBayes(mm_fit, cox_fit, timeVar = "time",
                                       Formulas = Forms_intslope,
                                       control = list(n_cores = ncores_fct))

  ####################
  #### Prediction ####
  ####################

  tLMs <- c(2,4)

  pred_dyn1 <- pred_dyn2 <- pred_dyn3 <- pred_dyn4 <- list()
  surv_JM_intslope <- list()
  pred_KM <- list()

  BS.fit <- AUC.fit <- msep <- biais <- list()

  for (i in seq(length(tLMs))){

    tLM <- tLMs[i]

    predTimes <- tLM + c(1,2)

    testData_long_tLM <- testData_long[which(testData_long$time<=tLM&testData_long$time_event>tLM),]
    testData_long_tLM <- testData_long_tLM[order(testData_long_tLM$id, testData_long_tLM$time),]

    testData_surv_tLM <- testData_surv[which(testData_surv$time_event>tLM),]

    # DynForest

    scalar_test <- list(X = testData_surv_tLM[,c("cont_covar1"), drop = FALSE], id = testData_surv_tLM$id)
    facto_test <- list(X = testData_surv_tLM[,c("bin_covar1"), drop = FALSE], id = testData_surv_tLM$id)
    curve_test <- list(X = testData_long_tLM[,c("marker1","marker2"), drop = FALSE], id = testData_long_tLM$id,
                       time = testData_long_tLM$time)

    cat("Prediction DynForest1","\n")

    if (!is.null(res_dyn1)){
      pred_dyn1[[i]] <- 1 - predict(res_dyn1, Curve = curve_test, Factor = facto_test, Scalar = scalar_test,
                                    predTimes = predTimes, t0 = tLM,
                                    ncores = ncores_fct)
    }else{
      pred_dyn1[[i]] <- NA
    }


    cat("Prediction DynForest2","\n")

    if (!is.null(res_dyn2)){
      pred_dyn2[[i]] <- 1 - predict(res_dyn2, Curve = curve_test, Factor = facto_test, Scalar = scalar_test,
                                    predTimes = predTimes, t0 = tLM,
                                    ncores = ncores_fct)
    }else{
      pred_dyn2[[i]] <- NA
    }

    cat("Prediction DynForest3","\n")

    if (!is.null(res_dyn3)){
      pred_dyn3[[i]] <- 1 - predict(res_dyn3, Curve = curve_test, Factor = facto_test, Scalar = scalar_test,
                                    predTimes = predTimes, t0 = tLM,
                                    ncores = ncores_fct)
    }else{
      pred_dyn3[[i]] <- NA
    }

    cat("Prediction DynForest4","\n")

    if (!is.null(res_dyn4)){
      pred_dyn4[[i]] <- 1 - predict(res_dyn4, Curve = curve_test, Factor = facto_test, Scalar = scalar_test,
                                    predTimes = predTimes, t0 = tLM,
                                    ncores = ncores_fct)
    }else{
      pred_dyn4[[i]] <- NA
    }

    # JMbayes intercept + slope

    cat("Prediction JMbayes 3","\n")

    pred_JM_intslope <- survfitJM(JM_fit_intslope, newdata = testData_long_tLM, idVar="id", last.time=tLM,
                                  survTimes = predTimes)

    if (length(predTimes)>1){
      surv_list <- lapply(pred_JM_intslope$full.results, function(x) do.call(rbind, x))
      surv_JM_intslope[[i]] <- Reduce("+", surv_list) / length(surv_list)
    }else{
      surv_JM_intslope[[i]] <- Reduce("+", pred_JM_intslope$full.results) / length(pred_JM_intslope$full.results)
    }

    # KM

    cat("Kaplan-Meier","\n")

    KM_fit <- survfit(Surv(time_event, evt) ~ 1, data = testData_surv_tLM)

    pred_KM[[i]] <- sapply(predTimes, function(x){
      id.time <- sum(KM_fit$time <= x)
      rep(KM_fit$surv[id.time], nrow(testData_surv_tLM))
    })

    ####################
    #### Assessment ####
    ####################

    pred_method <- c("pred_dyn1","pred_dyn2","pred_dyn3","pred_dyn4",
                     "surv_JM_intslope",
                     "pred_KM")

    # BS

    BS.fit[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      if(is.null(nrow(pred_x))){
        length_x <- length(pred_x)
      }else{
        length_x <- nrow(pred_x)
      }

      out <- sapply(predTimes, function(y){

        testData_surv_tLM_tHor <- testData_surv_tLM
        testData_surv_tLM_tHor$evt[which(testData_surv_tLM_tHor$time_event>y)] <- 0
        testData_surv_tLM_tHor$time_event[which(testData_surv_tLM_tHor$time_event>y)] <- y

        tryCatch(pec(object = cbind(rep(1, length_x), pred_x[,which(y==predTimes)]), formula = Surv(time_event, evt) ~ 1,
                     data = testData_surv_tLM_tHor, cens.model = "marginal",
                     exact = FALSE, times = c(0, y))$AppErr$matrix[-1],
                 error = function(e) return(NA))

      })

      return(t(out))

    })

    # AUC

    AUC.fit[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      if(is.null(nrow(pred_x))){
        length_x <- length(pred_x)
      }else{
        length_x <- nrow(pred_x)
      }

      out <- sapply(predTimes, function(y){

        testData_surv_tLM_tHor <- testData_surv_tLM
        testData_surv_tLM_tHor$evt[which(testData_surv_tLM_tHor$time_event>y)] <- 0
        testData_surv_tLM_tHor$time_event[which(testData_surv_tLM_tHor$time_event>y)] <- y

        tryCatch(timeROC(T = testData_surv_tLM_tHor$time_event, delta = testData_surv_tLM_tHor$evt,
                         marker = 1 - pred_x[,which(y==predTimes)],
                         cause = 1, iid = TRUE,
                         times = y - 0.001)$AUC[-1],
                 error = function(e) return(NA))

      })

      return(t(out))

    })

    # msep

    id_final <- unique(testData_long_tLM$id)

    msep[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      out <- sapply(seq(length(predTimes)), FUN = function(y){

        return(sum((pred_x[,y] - testData$true.prob[[i]][id_final,y])^2)/nrow(pred_x))

      })

      return(t(out))

    })

    # biais

    id_final <- unique(testData_long_tLM$id)

    biais[[i]] <- lapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      out <- sapply(seq(length(predTimes)), FUN = function(y){

        return(pred_x[,y] - testData$true.prob[[i]][id_final,y])

      })

      return(out)

    })

    names(biais[[i]]) <- pred_method

  }

  cat("AUC BS DONE !","\n")

  res[[replicate]] <- list(oob.err = list(res_dyn1$oob.err, res_dyn2$oob.err,
                                          res_dyn3$oob.err, res_dyn4$oob.err),
                           Var_depth = list(dyn_depth1, dyn_depth2, dyn_depth3, dyn_depth4),
                           VIMP = list(vimp_dyn1, vimp_dyn2, vimp_dyn3, vimp_dyn4),
                           pred = list(pred_dyn1, pred_dyn2, pred_dyn3, pred_dyn4,
                                       surv_JM_intslope,
                                       pred_KM),
                           BS = BS.fit, AUC = AUC.fit, msep = msep, biais = biais)


  cat("Saving...","\n")

  save(res, testData, trainData, file = "Scenarios/Small2/Small2_resu.RData")

}