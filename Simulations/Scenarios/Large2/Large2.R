#######################
### Large2 scenario ###
#######################

library(DynForest) # Install with DynForest_devlopment_version/DynForest_0.1.0.tar.gz
library(doParallel)
library(survival)
library(lcmm)
library(pec)
library(timeROC)

source("simu_function.R")

load("marker20.RData")

evt.link <- list(model = ~ I(marker1_class == 1) + I(marker1_class == 2) + I(marker3_class == 1) + I(marker3_class == 3)
                 ,
                 gamma = c(-2, 2, -1, 3))

res <- list()

set.seed(1000)

testData <- simuData(n = 500, marker = marker, class.lat = 4,
                     evt.link = evt.link, Y.type = "surv",
                     censoring.value = exp(-2.5),
                     var.timevis = TRUE,
                     truncation = TRUE, trunctime = 10,
                     baseline.surv.param = c(0.1, 2), baseline.surv.fct = "Weibull",
                     scale_marker = FALSE,
                     true.prob = T, tLMs = c(2,4), tHors = c(1,2))

testData_long <- testData$data.long
testData_surv <- subset(testData$data.surv, select = c(id, cont_covar1, bin_covar1, time_event, evt))

#######################
### Parallelisation ###
#######################

replicates.min <- 1
replicates.max <- 250
ncores_fct <- 12


for (replicate in replicates.min:replicates.max){

  cat(paste0(replicate,"/",replicates.max),"\n")

  cat("Simu data train","\n")

  set.seed(replicate)

  trainData <- simuData(n = 500, marker = marker, class.lat = 4,
                        evt.link = evt.link, Y.type = "surv",
                        censoring.value = exp(-2.5),
                        var.timevis = TRUE,
                        truncation = TRUE, trunctime = 10,
                        scale_marker = FALSE,
                        baseline.surv.param = c(0.1, 2), baseline.surv.fct = "Weibull")

  trainData_long <- trainData$data.long
  trainData_surv <- subset(trainData$data.surv, select = c(id, cont_covar1, bin_covar1, time_event, evt))

  ####################
  #### Estimation ####
  ####################


  # starting values mixed models

  model <- list()
  RE.df <- data.frame(id = seq(500))
  marker.fit <- list()

  for (i in seq(20)){

    mm.fixed.formula <- as.formula(paste0("marker",i,"~","splines::ns(time, knots = c(5),
                             Boundary.knots = c(0, 10))"))
    mm.random.formula <- as.formula(paste0("~","splines::ns(time, knots = c(5),
                             Boundary.knots = c(0, 10))"))

    marker.fit[[i]] <- hlme(fixed = mm.fixed.formula,
                            random = mm.random.formula,
                            subject = "id",
                            data = trainData_long)

    init.param <- marker.fit[[i]]$best

    model[[i]] <- list(fixed = mm.fixed.formula,
                       random = mm.random.formula,
                       init.param = init.param)

    RE <- marker.fit[[i]]$predRE
    colnames(RE)[-1] <- paste0("marker",i,"_RE", seq(ncol(RE)-1)-1)

    RE.df <- merge(RE.df, RE, by = "id")

  }

  names(model) <- paste0("marker",seq(20))

  # DynForest

  Y <- list(type = "surv", Y = Surv(trainData_surv$time_event, trainData_surv$evt), id = trainData_surv$id)
  scalar <- list(X = trainData_surv[,c("cont_covar1"), drop = FALSE], id = trainData_surv$id)
  facto <- list(X = trainData_surv[,c("bin_covar1"), drop = FALSE], id = trainData_surv$id)
  curve <- list(X = trainData_long[,paste0("marker", seq(20)), drop = FALSE], id = trainData_long$id,
                time = trainData_long$time,
                model = model)

  # tuning mtry

  if (length(res)==0){

    mtry.max <- ncol(curve$X)+ncol(scalar$X)+ncol(facto$X)
    mtrys <- seq(1, mtry.max, 4)
    #mtrys <- mtry.max
    OOB.err_Dyn <- vector("numeric", length(mtrys))
    best.OOB.err_Dyn <- Inf

    for (mtry in mtrys){

      cat(mtry, "\n")

      ind.mtry <- which(mtrys == mtry)

      res_dyn_mtry <- DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                Y=Y, ntree = 250, imp = FALSE,
                                mtry = mtry, nodesize = 3, minsplit = 5,
                                ncores = ncores_fct)

      OOB.err_Dyn[ind.mtry] <- mean(res_dyn_mtry$oob.err, na.rm = TRUE)

      if (OOB.err_Dyn[ind.mtry] < best.OOB.err_Dyn){

        best.OOB.err_Dyn <- OOB.err_Dyn[ind.mtry]
        best.mtry_Dyn <- mtry

      }

    }

    save(OOB.err_Dyn, best.mtry_Dyn, file = "Scenarios/Large2/Large2_mtry_DynForest.RData")
  }

  # best mtry

  set.seed(replicate)

  cat("Estimation DynForest","\n")

  # mtry opt

  res_dyn <- tryCatch(DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                Y=Y, ntree = 250, OOB_error = FALSE, imp = FALSE,
                                mtry = best.mtry_Dyn, nodesize = 3, minsplit = 5,
                                ncores = ncores_fct),
                      error = function(e) return(NULL))

  if (!is.null(res_dyn)){
    vimp.Dyn <- list(vimp = res_dyn$Importance,
                     var = res_dyn$Inputs)
  }else{
    vimp.Dyn <- NULL
  }

  ################################################
  # 2-steps mixed effect + DynForest

  cat("Estimation 2step MM","\n")

  scalar_2step_X <- merge(trainData_surv, RE.df, by = "id")

  scalar_2steps <- list(X = subset(scalar_2step_X, select = -c(id, bin_covar1, time_event, evt)),
                        id = scalar_2step_X$id)

  # tuning mtry

  if (length(res)==0){

    mtry.max <- ncol(scalar_2steps$X)+ncol(facto$X)
    mtrys <- seq(1, mtry.max, 5)
    #mtrys <- mtry.max
    OOB.err_2steps <- vector("numeric", length(mtrys))
    best.OOB.err_2steps <- Inf

    for (mtry in mtrys){

      cat(mtry, "\n")

      ind.mtry <- which(mtrys == mtry)

      res_dyn_mtry <- DynForest(Curve = curve, Factor = facto, Scalar = scalar,
                                Y=Y, ntree = 250, imp = FALSE,
                                mtry = mtry, nodesize = 3, minsplit = 5,
                                ncores = ncores_fct)

      OOB.err_2steps[ind.mtry] <- mean(res_dyn_mtry$oob.err, na.rm = TRUE)

      if (OOB.err_2steps[ind.mtry] < best.OOB.err_2steps){

        best.OOB.err_2steps <- OOB.err_2steps[ind.mtry]
        best.mtry_2steps <- mtry

      }

    }

    save(OOB.err_2steps, best.mtry_2steps, file = "Scenarios/Large2/Large2_mtry_2steps.RData")

  }

  # best mtry

  set.seed(replicate)

  cat("Estimation 2step DynForest","\n")

  res_2step_dyn <- tryCatch(DynForest(Curve = NULL, Factor = facto, Scalar = scalar_2steps,
                                      Y=Y, ntree = 250, OOB_error = FALSE, imp = FALSE,
                                      mtry = best.mtry_2steps,
                                      nodesize = 3, minsplit = 5,
                                      ncores = ncores_fct),
                            error = function(e) return(NULL))

  if (!is.null(res_2step_dyn)){
    dyn_2step_depth <- var_depth(res_2step_dyn)

    vimp.2steps <- list(vimp = res_2step_dyn$Importance,
                        var = res_2step_dyn$Inputs)
  }else{
    dyn_2step_depth <- NULL
    vimp.2steps <- NULL
  }

  ####################
  #### Prediction ####
  ####################

  tLMs <- c(2,4)

  pred_dyn <- list()
  pred_2step_dyn <- list()
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
    curve_test <- list(X = testData_long_tLM[, paste0("marker", seq(20)), drop = FALSE], id = testData_long_tLM$id,
                       time = testData_long_tLM$time)

    cat("Prediction DynForest","\n")

    if (!is.null(res_dyn)){
      pred_dyn[[i]] <- 1 - predict(res_dyn, Curve = curve_test, Factor = facto_test, Scalar = scalar_test,
                                   predTimes = predTimes, t0 = tLM,
                                   ncores = ncores_fct)
    }else{
      pred_dyn[[i]] <- NA
    }

    # 2-steps MM + DynForest

    cat("Prediction 2step DynForest","\n")

    pred.RE.df <- data.frame(id = seq(500))

    for (j in seq(20)){

      mm.fixed.formula <- as.formula(paste0("marker",j,"~","splines::ns(time, knots = c(5),
                             Boundary.knots = c(0, 10))"))
      mm.random.formula <- as.formula(paste0("~","splines::ns(time, knots = c(5),
                             Boundary.knots = c(0, 10))"))

      pred.RE <- hlme(fixed = mm.fixed.formula,
                      random = mm.random.formula,
                      subject = "id", data = testData_long_tLM,
                      #maxiter = 0,
                      B = marker.fit[[j]]$best)$predRE

      colnames(pred.RE)[-1] <- paste0("marker",j,"_RE", seq(ncol(pred.RE)-1)-1)

      pred.RE.df <- merge(pred.RE.df, pred.RE, by = "id")

    }

    scalar_2step_X_test <- merge(testData_surv_tLM, pred.RE.df, by = "id")

    scalar_2steps_test <- list(X = subset(scalar_2step_X_test, select = -c(id, bin_covar1, time_event, evt)),
                               id = scalar_2step_X_test$id)

    if (!is.null(res_2step_dyn)){
      pred_2step_dyn[[i]] <- 1 - predict(res_2step_dyn, Curve = NULL, Factor = facto_test,
                                         Scalar = scalar_2steps_test,
                                         predTimes = predTimes, t0 = tLM,
                                         ncores = ncores_fct)
    }else{
      pred_2step_dyn[[i]] <- NA
    }

    # KM

    cat("Kaplan-Meier","\n")

    KM_fit <- survfit(Surv(time_event, evt) ~ 1, data = testData_surv_tLM)

    pred_KM[[i]] <- sapply(predTimes, function(x){
      id.time <- sum(KM_fit$time <= x)
      rep(KM_fit$surv[id.time], nrow(testData_surv_tLM))
    })

    rownames(pred_KM[[i]]) <- testData_surv_tLM$id

    ####################
    #### Assessment ####
    ####################

    pred_method <- c("pred_dyn", "pred_2step_dyn",
                     "pred_KM")

    id_final <- Reduce(intersect, list(seq(500),
                                       rownames(get("pred_dyn")[[i]]),
                                       rownames(get("pred_2step_dyn")[[i]]),
                                       rownames(get("pred_KM")[[i]])
    )
    )

    # BS

    BS.fit[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      pred_x <- pred_x[which(rownames(pred_x)%in%id_final),]

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
                     data = testData_surv_tLM_tHor[which(testData_surv_tLM_tHor$id%in%id_final),], cens.model = "marginal",
                     exact = FALSE, times = c(0, y))$AppErr$matrix[-1],
                 error = function(e) return(NA))

      })

      return(t(out))

    })

    # AUC

    AUC.fit[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      pred_x <- pred_x[which(rownames(pred_x)%in%id_final),]

      if(is.null(nrow(pred_x))){
        length_x <- length(pred_x)
      }else{
        length_x <- nrow(pred_x)
      }

      out <- sapply(predTimes, function(y){

        testData_surv_tLM_tHor <- testData_surv_tLM
        testData_surv_tLM_tHor$evt[which(testData_surv_tLM_tHor$time_event>y)] <- 0
        testData_surv_tLM_tHor$time_event[which(testData_surv_tLM_tHor$time_event>y)] <- y

        tryCatch(timeROC(T = testData_surv_tLM_tHor$time_event[which(testData_surv_tLM_tHor$id%in%id_final)],
                         delta = testData_surv_tLM_tHor$evt[which(testData_surv_tLM_tHor$id%in%id_final)],
                         marker = 1 - pred_x[,which(y==predTimes)],
                         cause = 1, iid = TRUE,
                         times = y - 0.001)$AUC[-1],
                 error = function(e) return(NA))

      })

      return(t(out))

    })

    # msep

    msep[[i]] <- sapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      pred_x <- pred_x[which(rownames(pred_x)%in%id_final),]

      out <- sapply(seq(length(predTimes)), FUN = function(y){

        return(sum((pred_x[,y] - testData$true.prob[[i]][as.integer(id_final),y])^2)/nrow(pred_x))

      })

      return(t(out))

    })

    # biais

    biais[[i]] <- lapply(pred_method, function(x){

      pred_x <- get(x)
      pred_x <- pred_x[[i]]

      pred_x <- pred_x[which(rownames(pred_x)%in%id_final),]

      out <- sapply(seq(length(predTimes)), FUN = function(y){

        return(pred_x[,y] - testData$true.prob[[i]][as.integer(id_final),y])

      })

      return(t(out))

    })

    names(biais[[i]]) <- pred_method

  }

  cat("AUC BS DONE !","\n")

  res[[replicate]] <- list(oob.err = list(res_dyn$oob.err, res_2step_dyn$oob.err),
                           pred = list(pred_dyn,
                                       pred_2step_dyn,
                                       pred_KM),
                           BS = BS.fit, AUC = AUC.fit, msep = msep, biais = biais)

  best.mtry <- list(best.mtry_Dyn, best.mtry_2steps)

  cat("Saving...","\n")

  save(res, best.mtry, testData, trainData, file = "Scenarios/Large2/Large2_resu.RData")

}