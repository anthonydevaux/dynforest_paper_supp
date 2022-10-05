simuData <- function(n = 100, marker = NULL, nb.marker = 1,
                     evt.link = NULL, form.predlin = "linear", nb.predlin.var = 1,
                     Y.type = "scalar", Y.class.lat = FALSE,
                     pNA = rep(0,length(marker)),
                     vis = seq(0,9,1), var.timevis = TRUE,
                     nb.cont.covar = 1, param.cont.covar = list(c(0,1)),
                     nb.bin.covar = 1, param.bin.covar = 0.5, class.lat = 4,
                     baseline.surv.fct = "Weibull", baseline.surv.param = NULL,
                     censoring.value = exp(-3),
                     truncation = TRUE, trunctime = 10,
                     scaling = FALSE, scale_marker = TRUE,
                     true.prob = FALSE, tLMs = NULL, tHors = NULL){

  library(mvtnorm)

  if (length(param.cont.covar)!=nb.cont.covar){
    stop("param.cont.covar's length should be equal to nb.cont.covar", "\n")
  }

  if (length(param.bin.covar)!=nb.bin.covar){
    stop("nb.bin.covar's length should be equal to param.bin.covar", "\n")
  }

  nbvis <- length(vis)

  if (is.null(marker)){

    marker <- list()

    for (i in seq(nb.marker)){

      marker[[i]] <- list()

      #complexity <- floor(runif(1, 2, 5))
      complexity <- 2
      marker_courant <- paste0("marker",i)

      for (num_class in seq(class.lat)){

        if (complexity==2){
          fixed <- reformulate(termlabels = "time",
                               response = marker_courant)
          random <- ~ time
          params <- list(beta = round(c(runif(1, 2, 5), runif(1, -4, 4)), 2),
                         sd.re = rep(1, 2),
                         cor.re = 0.3,
                         sigmae = 1)
        }

        if (complexity==3){
          fixed <- reformulate(termlabels = "splines::ns(time, knots = 5, Boundary.knots = c(0,10))",
                               response = marker_courant)
          random <- ~ splines::ns(time, knots = 5, Boundary.knots = c(0,10))
          params <- list(beta = round(c(runif(1, 2, 5), runif(2, -4, 4)), 2),
                         sd.re = rep(1, 3),
                         cor.re = rep(0.3, 3),
                         sigmae = 1)
        }

        if (complexity==4){
          fixed <- reformulate(termlabels = "splines::ns(time, knots = c(3,6), Boundary.knots = c(0,10))",
                               response = marker_courant)
          random <- ~ splines::ns(time, knots = c(3,6), Boundary.knots = c(0,10))
          params <- list(beta = round(c(runif(1, 2, 5), runif(3, -4, 4)), 2),
                         sd.re = rep(1, 4),
                         cor.re = rep(0.3, 6),
                         sigmae = 1)
        }

        marker[[i]][[num_class]] <- list(fixed = fixed, random = random,
                                         params = params)

      }

      marker[[i]][[num_class+1]] <- "cont"
      names(marker[[i]]) <- c(paste0("class",seq(class.lat)), "type.var")

    }

    names(marker) <- paste0("marker", seq(nb.marker))

  }else{

    nb.marker <- length(marker)

  }

  ###############################

  donnees <- matrix(NA, nrow = n * nbvis,
                    ncol = 2 + nb.cont.covar + nb.bin.covar + nb.marker)

  colnames(donnees) <- c("id", "time", paste0("cont_covar", seq(nb.cont.covar)),
                         paste0("bin_covar", seq(nb.bin.covar)),
                         names(marker))

  donnees <- as.data.frame(donnees)

  bi.mat.marker <- lapply(marker, function(x){

    tirage <- sample(seq(class.lat), n, replace = TRUE) # structure latente specifique a chaque marqueur

    marker_class <- list()

    for (num_class in seq(class.lat)){

      if (length(x[[num_class]]$params$sd.re)==1){

        B <- x[[num_class]]$params$sd.re[1]
        bi.mat <- rnorm(length(which(tirage==num_class)), 0, B)

      }

      if (length(x[[num_class]]$params$sd.re)==2){

        covar <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$cor.re

        B <- matrix(c(x[[num_class]]$params$sd.re[1]^2, covar,
                      covar, x[[num_class]]$params$sd.re[2]^2), ncol = 2)

        if (det(B)==0){
          B <- B + diag(ncol(B))*0.01
        }

        bi.mat <- rmvnorm(length(which(tirage==num_class)), mean = rep(0, ncol(B)), sigma = B)

      }

      if (length(x[[num_class]]$params$sd.re)==3){

        covar01 <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$cor.re[1]
        covar02 <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[3] * x[[num_class]]$params$cor.re[2]
        covar12 <- x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$sd.re[3] * x[[num_class]]$params$cor.re[3]

        B <- matrix(c(x[[num_class]]$params$sd.re[1]^2, covar01, covar02,
                      covar01, x[[num_class]]$params$sd.re[2]^2, covar12,
                      covar02, covar12, x[[num_class]]$params$sd.re[3]^2), ncol = 3)

        if (det(B)==0){
          B <- B + diag(ncol(B))*0.01
        }

        bi.mat <- rmvnorm(length(which(tirage==num_class)), mean = rep(0, ncol(B)), sigma = B)

      }

      if (length(x[[num_class]]$params$sd.re)==4){

        covar01 <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$cor.re[1]
        covar02 <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[3] * x[[num_class]]$params$cor.re[2]
        covar03 <- x[[num_class]]$params$sd.re[1] * x[[num_class]]$params$sd.re[4] * x[[num_class]]$params$cor.re[3]
        covar12 <- x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$sd.re[3] * x[[num_class]]$params$cor.re[4]
        covar13 <- x[[num_class]]$params$sd.re[2] * x[[num_class]]$params$sd.re[4] * x[[num_class]]$params$cor.re[5]
        covar23 <- x[[num_class]]$params$sd.re[3] * x[[num_class]]$params$sd.re[4] * x[[num_class]]$params$cor.re[6]

        B <- matrix(c(x[[num_class]]$params$sd.re[1]^2, covar01, covar02, covar03,
                      covar01, x[[num_class]]$params$sd.re[2]^2, covar12, covar13,
                      covar02, covar12, x[[num_class]]$params$sd.re[3]^2, covar23,
                      covar03, covar13, covar23, x[[num_class]]$params$sd.re[4]^2), ncol = 4)

        if (det(B)==0){
          B <- B + diag(ncol(B))*0.01
        }

        bi.mat <- rmvnorm(length(which(tirage==num_class)), mean = rep(0, ncol(B)), sigma = B)

      }

      rownames(bi.mat) <- which(tirage==num_class)
      colnames(bi.mat) <- paste0(as.character(x[[num_class]]$fixed)[2], "_bi", c(0:(ncol(bi.mat)-1)))
      #marker_class[[num_class]] <- bi.mat[, 1, drop = FALSE] # on garde que intercept
      marker_class[[num_class]] <- bi.mat

    }

    bi.mat.class <- do.call(rbind, marker_class)
    bi.mat.class <- bi.mat.class[order(as.integer(rownames(bi.mat.class))), , drop = FALSE]

    return(list(bi.mat.class = bi.mat.class, tirage = tirage))

  })

  tirage_list <- lapply(bi.mat.marker, FUN = function(x) return(x$tirage))
  tirage <- do.call(cbind, tirage_list)
  colnames(tirage) <- paste(colnames(tirage), "class", sep = "_")
  bi.mat.marker <- lapply(bi.mat.marker, FUN = function(x) return(x$bi.mat.class))

  for (i in seq(n)){

    donnees[((nbvis*(i-1))+1):(nbvis*i),1] <- i # id
    donnees[((nbvis*(i-1))+1),2] <- vis[1] # t0

    if (var.timevis){
      donnees[((nbvis*(i-1))+2):(nbvis*i),2] <- vis[-1] + rexp(length(vis[-1]), rate = 5) # variation temps de vis
    }else{
      donnees[((nbvis*(i-1))+2):(nbvis*i),2] <- vis[-1] + rep(0, length(vis[-1]))
    }

    # generation covariables

    for (ind.nb.cont.covar in 1:nb.cont.covar){

      mean <- param.cont.covar[[ind.nb.cont.covar]][1]
      sd <- param.cont.covar[[ind.nb.cont.covar]][2]

      donnees[((nbvis*(i-1))+1):(nbvis*i),2+ind.nb.cont.covar] <- rnorm(1, mean, sd)

    }

    for (ind.nb.bin.covar in 1:nb.bin.covar){

      p <- param.bin.covar[ind.nb.bin.covar]

      donnees[((nbvis*(i-1))+1):(nbvis*i),2+nb.cont.covar+ind.nb.bin.covar] <- rbinom(1, 1, p)

    }

    for (j in seq(nbvis)){

      Y <- mapply(function(x, idx){

        num_class <- tirage_list[[idx]][i]

        X <- model.matrix(x[[num_class]]$fixed[-2], data = donnees[(nbvis*(i-1))+j,])
        Z <- model.matrix(x[[num_class]]$random, data = donnees[(nbvis*(i-1))+j,])

        beta <- x[[num_class]]$params$beta

        bi <- bi.mat.marker[[idx]][i,]

        eta <- X%*%beta + Z%*%bi

        if (x$type.var=="bin"){eta <- rbinom(1, 1 , prob = exp(eta)/(1+exp(eta)))}
        if (x$type.var=="cont"){eta <- eta + rnorm(1,0,x[[num_class]]$params$sigmae)}

        return(eta)

      }, marker, seq_along(marker))

      donnees[(nbvis*(i-1))+j,names(Y)] <- Y

    }

  }

  if (any(pNA<0 | pNA>1)){
    stop("pNA takes values in [0;1]")
  }

  if (any(pNA > 0 & pNA < 1)){
    censures <- lapply(pNA,
                       function(x){
                         ifelse(c(t(cbind(rep(0,n),
                                          matrix(rbinom((nbvis-1)*n,1,x),nrow = n, ncol = nbvis-1)
                         )))==1,NA,0)})

    donnees[,(3+nb.cont.covar+nb.bin.covar):(2+nb.cont.covar+nb.bin.covar+length(censures))] <-
      mapply(function(x, idx){ifelse(is.na(x),NA,donnees[,2+nb.cont.covar+nb.bin.covar+idx])},
             censures, seq_along(censures))
  }

  marker.bin <- unlist(lapply(marker, function(x){x$type.var=="binaire"}))

  var.bin <- c(names(marker.bin)[which(marker.bin==TRUE)],paste0("bin_covar", seq(nb.bin.covar)))

  for (varbin in var.bin){

    donnees[,varbin] <- as.factor(donnees[,varbin])

  }

  data.surv <- cbind(unique(donnees[,which(!colnames(donnees)%in%c("time", names(marker)))]),
                     do.call(cbind, bi.mat.marker)) # ajout des resumes

  rownames(data.surv) <- 1:nrow(data.surv)

  # Scaling ?

  if (scaling){

    data.surv[,(nb.cont.covar+nb.bin.covar+3):ncol(data.surv)] <-
      scale(data.surv[,(nb.cont.covar+nb.bin.covar+3):ncol(data.surv)])

  }

  # generation outcome

  #####################
  # matrice des parametres marqueur/class latente specifique pour chaque individu

  mat.beta <- mapply(function(x, idx){

    test <- sapply(tirage_list[[idx]], FUN = function(y){

      return(x[[y]]$params$beta)

    })

    return(t(test))

  }, marker, seq_along(marker), USE.NAMES = FALSE, SIMPLIFY = FALSE)

  mat.beta <- do.call(cbind, mat.beta)

  seq.name.beta <- unlist(lapply(marker, FUN = function(x) return(seq(x[[1]]$params$beta)-1)))
  length.beta <- lapply(marker, FUN = function(x) return(length(x[[1]]$params$beta)))
  nb.beta.param <- sum(unlist(length.beta))

  colnames(mat.beta) <- paste0("beta", seq.name.beta, "_", rep(1:length(marker), unlist(length.beta)))
  data.surv <- cbind(data.surv, mat.beta, tirage)

  #####################
  # calcul vrai bi

  if (true.prob){

    bi.mat0 <- mapply(function(x, idx){

      bi.mat0.marker.list <- vector("list", length = class.lat)

      for (class in 1:class.lat){

        id_class <- which(tirage_list[[idx]]==class)
        donnees_class <- donnees[which(donnees$id%in%id_class),]
        bi.mat0 <- predRE_simu(x[[class]], donnees_class)$b_i
        colnames(bi.mat0) <- paste0(as.character(x[[class]]$fixed)[2], "_bi", c(0:(ncol(bi.mat0)-1)))
        bi.mat0.marker.list[[class]] <- bi.mat0

      }

      bi.mat0.marker <- do.call(rbind, bi.mat0.marker.list)
      bi.mat0.marker <- bi.mat0.marker[order(as.integer(rownames(bi.mat0.marker))), , drop = FALSE]

      return(bi.mat0.marker)

    }, marker, seq_along(marker), SIMPLIFY = FALSE)

    bi.mat0 <- do.call(cbind, bi.mat0)
    data.surv0 <- cbind(data.surv[,c(1:3)], bi.mat0, mat.beta, tirage)
  }

  #####################

  if (is.null(evt.link)){

    evt.link <- list()

    if (Y.class.lat){

      listoffactors <- paste0("I(class == ", seq(class.lat), ")")

      evt.link$model <- reformulate(termlabels = listoffactors,
                                    response = NULL)

      evt.link$gamma <- round(runif(length(listoffactors), -5, 5), 2)

      evt.link$sigmae <- round(runif(1), 2)

    }else{

      num.predlin.var <- sort(sample(seq(nb.beta.param), nb.predlin.var))
      listoffactors <- colnames(data.surv[,(nb.cont.covar+nb.bin.covar+2):
                                            (nb.cont.covar+nb.bin.covar+2+nb.beta.param-1)])[num.predlin.var]

      listofbetas <- colnames(data.surv[,(nb.cont.covar+nb.bin.covar+2+nb.beta.param+1):
                                          ncol(data.surv)])[num.predlin.var]

      if (form.predlin=="linear"){

        evt.link$model <- reformulate(termlabels = listoffactors,
                                      response = NULL)



      }

      if (form.predlin=="non-linear"){

        listoffactors2 <- paste0("I(", listoffactors[1:floor(length(listoffactors)/2)], "^2)")
        listoffactorsbin <- paste0("I(", listoffactors[(floor(length(listoffactors)/2)+1):length(listoffactors)],
                                   " > median(", listoffactors[(floor(length(listoffactors)/2)+1):length(listoffactors)], "))")

        listoffactors <- c(listoffactors2, listoffactorsbin)

        evt.link$model <- reformulate(termlabels = listoffactors,
                                      response = NULL)

      }

      evt.link$gamma <- round(runif(length(listoffactors), -5, 5), 2)

    }

    if (Y.type == "scalar"){

      evt.link$gamma <- c(round(runif(1, 5, 10), 2), evt.link$gamma) # Y intercept
      evt.link$sigmae <- round(runif(1), 2)

    }

  }

  if (Y.type == "scalar"){ # Y scalar

    X <- model.matrix(evt.link$model, data = data.surv)
    Y_res <- X%*%evt.link$gamma + rnorm(nrow(X), 0, evt.link$sigmae)

    data.surv <- cbind(data.surv, Y_res)

    donnees <- merge(donnees, data.surv[,c("id","Y_res")])

  }

  if (Y.type == "surv"){

    X <- as.matrix(model.matrix(evt.link$model, data = data.surv)[,-1]) # no intercept

    if (scale_marker){
      X <- scale(X) # scaling
    }

    predlin <- X%*%evt.link$gamma

    if (baseline.surv.fct=="Weibull"){

      b <- baseline.surv.param[1]
      c <- baseline.surv.param[2]

      u <- runif(nrow(X))

      time_event <- as.vector((1/b)*(-log(u)/exp(predlin))^(1/c))
      #time_event <- as.vector((-log(u)/(b*exp(predlin)))^(1/c)) # Bender et al. (2005)

      if (true.prob){

        X0 <- as.matrix(model.matrix(evt.link$model, data = data.surv0)[,-1]) # no intercept

        if (scale_marker){
          X0 <- scale(X0) # scaling
        }

        predlin0 <- X0%*%evt.link$gamma

        prob0 <- lapply(tLMs, function(x){

          Ss <- exp(-exp(as.vector(predlin0)) * (b*x)^c) # S(S)
          Fs <- 1 - Ss # F(S)

          prob0_tLM <- sapply(tHors, function(y){

            Fst <- 1 - exp(-exp(as.vector(predlin0)) * (b*(x+y))^c) # F(S+t)

            return((Fst-Fs)/Ss)

          })

          return(1 - prob0_tLM)

        })

      }

    }

    if (baseline.surv.fct=="Exponential"){

      theta <- baseline.surv.param[1]
      u <- runif(nrow(X))

      time_event <- as.vector(-log(u)/(theta*exp(predlin))) # Bender et al. (2005)

    }

    if (baseline.surv.fct=="Gompertz"){

      alpha <- baseline.surv.param[1]
      gamma <- baseline.surv.param[2]
      u <- runif(nrow(X))

      time_event <- as.vector((1/alpha)*log(1- (alpha*log(u))/(gamma*exp(predlin)))) # Bender et al. (2005)

    }

    if (!is.null(evt.link$censoring.value)){

      predlin.censo <- X%*%evt.link$censoring.value

      b <- baseline.surv.param[1]
      c <- baseline.surv.param[2]

      u <- runif(nrow(X))

      censure <- as.vector((1/b)*(-log(u)/exp(predlin.censo))^(1/c)) # weibull

    }else{
      censure <- -log(runif(n)) / censoring.value # loi exponentiel
      #censure <- runif(n, 1, 30) # generation censure selon loi uniforme
    }


    if (truncation) {
      censure <- pmin(censure, trunctime)
    }

    evt <- ifelse(time_event<censure, 1, 0)
    time_event <- pmin(time_event, censure)

    data.surv <- cbind(data.surv, time_event, evt)
    donnees <- merge(donnees, data.surv[,c("id","time_event","evt")])
    donnees <- donnees[which(donnees$time<donnees$time_event),]

    cat(paste0(round(100 * sum(evt) / n, 1), "% experienced event\n"))

  }

  if (true.prob){
    res <- list(data.long = donnees,
                data.surv = data.surv,
                class = tirage,
                true.prob = prob0,
                marker = marker,
                evt.link = evt.link)
  }else{
    res <- list(data.long = donnees,
                data.surv = data.surv,
                class = tirage,
                marker = marker,
                evt.link = evt.link)
  }

  return(res)

}
