# IRT probabilities
pr <- function(abilities, diff, discr, guess, eps=1.4e-07) 
{
  lp <- discr * (abilities - diff)
  prob <- plogis(lp)
  prob[prob >= 1-eps] <- 1 - 1-eps
  prob[prob <= eps] <- eps
  return(guess + (1 - guess) * prob)
}



# simulation of a dataset from a 3PL model
# IRT parameterization
simtpl <- function(n, diff, discr, guess)
{
  abilities <- rnorm(n, 0, 1) # generation of abilities
  # generation of probabilities
  nitems <- length(diff)
  prob <- matrix(NA, n, nitems)
  for(j in 1:nitems) prob[,j] <- pr(abilities, diff[j], discr[j], guess[j])
  # generation of responses
  data <- (prob > matrix(runif(n * nitems), n, nitems)) * 1 
  colnames(data) <- paste0("I", formatC(1:nitems, flag = "0", width = 2))
  return(data)
}


# shrinkage estimation of a 3PL model
tpl <- function(data, penalized = FALSE, BR = FALSE, lambda = 0, InfCrit = FALSE, cv = FALSE,
                K = 5, nq = 21, empirical = TRUE, M = 5000, ncores = NULL, trace = TRUE, init_BR = FALSE)
{
  if (penalized & all(lambda==0)) stop("If argument penalized is TRUE lambda should take a positive value.")
  if (length(lambda) > 1 & lambda[1]==0) lambda <- lambda[-1]
  # maximum likelihood estimation of a 3PL model 
  #if (trace) message("=============================================\n")
  if (trace) message("* running maximum likelihood estimation")
  #if (trace) message("=============================================\n")
  res_MLE <- fit3pl_MLE(data = data, nq = nq, InfCrit = InfCrit)
  est_MLE <- res_MLE$est_MLE
  I <- ncol(data)
  init <- if(init_BR) est_MLE[c((I + 1):(3 * I), 1:I)] else NULL 
  lik <- res_MLE$lik
  conv_MLE <- res_MLE$conv_MLE
  aic_MLE <- res_MLE$aic
  bic_MLE <- res_MLE$bic
  cv_res <- NULL
  est_pen <- NULL
  conv_pen <- NULL
  est_pen_sel_cv <- est_pen_sel_aic <- est_pen_sel_bic <- NULL
  cv_res <- list(est = NULL, lik = NULL, sel = NULL)
  lambda_sel_aic <- lambda_sel_bic <- lambda_sel_cv <- NULL
  aic_pen <- bic_pen <- NULL
  if (penalized)
    # ridge-type penalty
  {
    #if (trace) message("=============================================\n")
    if (trace) message("* running penalized likelihood estimation")
    #if (trace) message("=============================================\n")
    res_ridgepen <- fit3pl_ridgepen(est_MLE = est_MLE, lambda = lambda, data = data, nq = nq, K = K, trace = trace, cv = cv, InfCrit = InfCrit)
    est_pen <- res_ridgepen$est_pen
    cv_res <- res_ridgepen$cv_res
    conv_pen <- res_ridgepen$conv_pen
    if (cv) 
    {
      est_pen_sel_cv <- est_pen[, cv_res$sel]
      lambda_sel_cv <- c(0, lambda)[cv_res$sel]
    }
    if (InfCrit) 
    {
      aic_pen <- res_ridgepen$aic
      bic_pen <- res_ridgepen$bic
      sel_aic <- which.min(c(aic_MLE, aic_pen))
      lambda_sel_aic <- c(0,lambda)[sel_aic]
      sel_bic <- which.min(c(bic_MLE, bic_pen))
      lambda_sel_bic <- c(0,lambda)[sel_bic]
      est_pen_sel_aic <- est_pen[, sel_aic]
      est_pen_sel_bic <- est_pen[, sel_bic]
    }
  }
  
  est_BR <- NULL
  conv_BR <- NULL
  # bias reduction
  if (BR)
  {
    #if (trace) message("=============================================\n")
    if (trace) message("* running estimation with bias reduction methodology")
    #if (trace) message("=============================================\n")
    res_BR <- fit3pl_BR(data = data, nq = nq, ncores = ncores, empirical = empirical, M = M, trace = trace, init)
    est_BR <- res_BR$ris.br
    conv_BR <- res_BR$conv.br / nrow(data)
  }
  out <- list(est_MLE = est_MLE, logLik = lik, conv_MLE = conv_MLE, AIC_MLE = aic_MLE, BIC_MLE = bic_MLE,
              est_pen = est_pen,conv_pen = conv_pen, AIC_pen = aic_pen, BIC_pen = bic_pen,
              lambda_sel_AIC = lambda_sel_aic, lambda_sel_BIC = lambda_sel_bic,
              est_pen_sel_AIC = est_pen_sel_aic, est_pen_sel_BIC = est_pen_sel_bic,
              cv_est = cv_res$est, cv_lik = cv_res$lik, lambda_sel_cv = lambda_sel_cv,
              est_pen_sel_cv = est_pen_sel_cv,
              lambda = lambda, est_BR = est_BR,conv_BR = conv_BR, itemnames = colnames(data))
  class(out) <- "tpl"
  return(out)
}


fit3pl_MLE <- function(data, nq, InfCrit) {
  suppressMessages(mod <- try(mirt(data, 1, itemtype = '3PL',verbose = FALSE), silent = TRUE))
  mirtest <- import.mirt(mod, display = FALSE) # estimates, parameterization used for estimation
  gq <- gauss.quad.prob(nq, dist = "normal")
  opt <- optim(par = as.vector(mirtest$coef), fn = IRTlikRcppA, gr = gradIRTlikRcppA, data = data,
               pen = FALSE, nodes = gq$nodes, weights = gq$weights, hessian = FALSE,
               control = list(fnscale = -1), method = "BFGS")
  est_MLE <- opt$par
  lik <- opt$value
  conv_MLE <- opt$convergence
  aic <- bic <- NULL
  if (InfCrit)
  {
    deri <- grad_i_IRTlikRcppA(par = opt$par,data = data, pen = FALSE,
                               nodes = gq$nodes, weights = gq$weights)
    hess <- jacobian(func = gradIRTlikRcppA, x = opt$par,data = data, pen = FALSE,
                     nodes = gq$nodes, weights = gq$weights,
                     method="simple")
    df <- try(sum(diag(solve(- hess) %*% t(deri) %*% deri)))
    if (inherits(df, "try-error")) df <- NA
    #print(df)
    aic <- 2 * df - 2 * lik
    bic <- df * log (nrow(data)) - 2 * lik
  }
  return(list(est_MLE = est_MLE, lik = lik, conv_MLE = conv_MLE, aic = aic, bic = bic))
}


fit3pl_ridgepen <- function(est_MLE, lambda, data, nq, K, trace, cv, InfCrit)
{
  gq <- statmod::gauss.quad.prob(nq, dist = "normal")
  est_pen <- matrix(est_MLE, ncol = 1)
  conv<-c()
  aic <- bic <- NULL
  dfout<-c()
  likout<-c()
  for (i in seq_along(lambda) + 1)
  {
    ini <- est_pen[,i-1]
    lambda_i <- lambda[i-1]
    opt1 <- optim(par = ini, fn = IRTlikRcppA, gr = gradIRTlikRcppA, method = "BFGS", data = data, pen = TRUE,   
                  lambda = lambda_i, nodes = gq$nodes, weights = gq$weights, 
                  control = list(fnscale = -1, maxit = 10000))
    est_pen <- cbind(est_pen, opt1$par)
    conv <- c(conv, opt1$convergence)
    if (InfCrit)
    {
      deri <- grad_i_IRTlikRcppA(par = opt1$par,data = data, pen = TRUE,
                                 lambda = lambda_i, nodes = gq$nodes, weights = gq$weights)
      hess <- jacobian(func = gradIRTlikRcppA, x = opt1$par,data = data, pen = TRUE,
                       lambda = lambda_i, nodes = gq$nodes, weights = gq$weights,
                       method="simple")
      df <- try(sum(diag(solve(- hess) %*% t(deri) %*% deri)))
      if (inherits(df, "try-error")) df <- NA
      dfout<-c(dfout,df)
      #print(df)
      lik <- IRTlikRcppA(par = opt1$par, data = data, pen = FALSE, nodes = gq$nodes, weights = gq$weights)
      #print(lik)
      likout<-c(likout,lik)
      aic <- c(aic, 2 * df - 2 * lik)
      bic <- c(bic, df * log (nrow(data)) - 2 * lik)
    }
  }
  cv_res <- NULL
  if (cv) cv_res <- cv_tpl_ridge(data = data, lambda = lambda, K = K, est_MLE = est_MLE, gq = gq, trace = trace)
  return(list(est_pen = est_pen, conv_pen = conv, cv_res = cv_res, 
              aic = aic, bic = bic, df=dfout,lik=likout))
}


fit3pl_BR <- function(data, nq, ncores, empirical, M, trace, init) {
  if (!is.matrix(data)) data<-as.matrix(data)
  trace <- trace * 1
  I <-  ncol(data) #number of items
  parameters <- list(alpha = rep(0, I), beta=rep(1, I), gamma = rep(0, I))
  gh <- statmod::gauss.quad(nq, "hermite")
  if (is.null(ncores) & .Platform$OS.type != "windows") ncores <- parallel::detectCores()-1
  if (is.null(ncores) & .Platform$OS.type == "windows") ncores <- 1
  if (.Platform$OS.type == "windows" & ncores > 1)
  {
    message("'mc.cores' > 1 is not supported on Windows; ncores set to 1")
    ncores <- 1
  }
  data.listM <- list(y = as.matrix(data), y1 = as.matrix(1 - data),
                     nodes=as.vector(gh$nodes), weights=as.vector(gh$weights))
  obj3 <- TMB::MakeADFun(data = data.listM, parameters = parameters, DLL = "BR3PL", silent = TRUE)
  start.val <- if(is.null(init))  obj3$par else init
  
  
  if (empirical)
  {
    # empirical approximation
    brJ <- nleqslv::nleqslv(start.val, gsolv_pack, ADobj = obj3, gh = gh, ncores = ncores, M = NULL,
                            seed = NULL, y = data,  control=list(trace = trace), method = "Newton", jac = jsolv)
    ris.br <- brJ$x[c((2 * I + 1):(3 * I), 1:(2 * I))] # put guessing first
    conv.br <- max(abs(brJ$fvec))
  }
  # Monte Carlo approximation
  else
  {
    brJ <- nleqslv::nleqslv(start.val, gsolv_pack, ADobj = obj3, gh = gh, ncores = ncores, M = M, 
                            seed = NULL, y = data,  control = list(trace = trace), method = "Newton", 
                            jac = jsolv) 
    ris.br <- brJ$x[c((2 * I + 1):(3 * I), 1:(2 * I))] # put guessing first
    conv.br <- max(abs(brJ$fvec))
  }
  return(list(ris.br = ris.br, conv.br = conv.br))
}


print.tpl <- function(x, IRTparam = TRUE,...)
{
  est_MLE <- matrix(x$est_MLE, ncol = 3)
  rownames(est_MLE) <- x$itemnames
  if (IRTparam) 
  {
    est_MLE[,1] <-plogis(est_MLE[, 1])
    est_MLE[,2] <- -est_MLE[, 2] / est_MLE[, 3]
    colnames(est_MLE) <- c("c", "b", "a")
  }
  else colnames(est_MLE) <- c("beta3", "beta1", "beta2")
  message("Maximum likelihood estimates")
  print.default(est_MLE, quote = FALSE, right = TRUE)
  message("\n")
  
  if (!is.null(x$est_pen))
  {
    message("Penalized estimates \n")
    if (length(x$lambda) == 1)
    {
      est_pen_sel <- matrix(x$est_pen[,2], ncol = 3)
      message("(only 1 value for lambda)")
    }
    else
    {
      if (!is.null(x$cv_est))
      {
        message("selected by cross-validation")
        est_pen_sel<-matrix(x$est_pen_sel_cv,ncol=3)
        rownames(est_pen_sel)<-x$itemnames
        if (IRTparam) 
        {
          est_pen_sel[,1] <- plogis(est_pen_sel[, 1])
          est_pen_sel[,2] <- - est_pen_sel[, 2] / est_pen_sel[, 3]
          colnames(est_pen_sel) <- c("c", "b", "a")
        }
        else colnames(est_pen_sel) <- c("beta3", "beta1", "beta2")
        print.default(est_pen_sel, quote = FALSE, right = TRUE)
        message("\n")
      }
      if (!is.null(x$est_pen_sel_AIC))
      {
        message("selected by AIC")
        est_pen_sel<-matrix(x$est_pen_sel_AIC,ncol=3)
        rownames(est_pen_sel)<-x$itemnames
        if (IRTparam) 
        {
          est_pen_sel[,1] <- plogis(est_pen_sel[, 1])
          est_pen_sel[,2] <- - est_pen_sel[, 2] / est_pen_sel[, 3]
          colnames(est_pen_sel) <- c("c", "b", "a")
        }
        else colnames(est_pen_sel) <- c("beta3", "beta1", "beta2")
        print.default(est_pen_sel, quote = FALSE, right = TRUE)
        message("\n")
        message("selected by BIC")
        est_pen_sel<-matrix(x$est_pen_sel_BIC,ncol=3)
        rownames(est_pen_sel)<-x$itemnames
        if (IRTparam) 
        {
          est_pen_sel[,1] <- plogis(est_pen_sel[, 1])
          est_pen_sel[,2] <- - est_pen_sel[, 2] / est_pen_sel[, 3]
          colnames(est_pen_sel) <- c("c", "b", "a")
        }
        else colnames(est_pen_sel) <- c("beta3", "beta1", "beta2")
        print.default(est_pen_sel, quote = FALSE, right = TRUE)
        message("\n")
      }
      if (is.null(x$cv_est) & is.null(x$est_pen_sel_AIC)) 
        message("Cross-validation or information criteria are not available to select the penalized estimates.")
    }
  }
  
  if (!is.null(x$est_BR))
  {
    est_BR <- matrix(x$est_BR, ncol = 3)
    rownames(est_BR) <- x$itemnames
    if (IRTparam) 
    {
      est_BR[,1] <- plogis(est_BR[,1])
      est_BR[,2] <- -est_BR[,2] / est_BR[,3]
      colnames(est_BR) <- c("c", "b", "a")
    }
    else colnames(est_BR) <- c("beta3", "beta1", "beta2")
    message("Bias reduction estimates")
    print.default(est_BR, quote = FALSE, right = TRUE)
    message("\n")
  }
}


cv_tpl_ridge<-function(data, lambda, K, est_MLE, gq, trace)
{
  #if (trace) message("=============================================\n")
  if (trace) message("* running cross validation")
  #if (trace) message("=============================================\n")
  n <- nrow(data)
  # the following for excluding subsets with columns all 0 or all 1
  all01 <- TRUE
  while (all01) {
    all01_tmp <- FALSE
    gr <- split(sample(n, n, replace = FALSE), as.factor(1:K)) # generation of subsets for CROSS VALIDATION
    for (k in 1:K)
    {
      data_k <- data[-gr[[k]],]
      if (any(colMeans(data_k, na.rm = TRUE) == 0 | colMeans(data_k, na.rm = TRUE) == 1)) all01_tmp <- TRUE
    }
    all01 <- all01_tmp
  }
  nitems <- ncol(data)
  nlambda <- length(lambda) + 1
  est <- vector("list", nlambda) # estimates
  lik <- vector("list", nlambda) # likelihood index
  
  # ===============================
  # not penalized
  # ===============================
  
  if (trace) message("lambda = 0      ", appendLF = FALSE)
  est[[1]] <- matrix(NA, K, nitems * 3) # estimates on the training set
  lik[[1]] <- rep(NA, K) # not penalized likelihood on the validation
  if (trace) message("folds: ", appendLF = FALSE)
  for (k in 1:K) # k = subset
  {
    if (trace) message(k ,appendLF = FALSE)
    data_k<-data[-gr[[k]],] # training set
    suppressMessages(mod <- try(mirt(data_k, 1, itemtype = '3PL', verbose = FALSE), silent = TRUE))
    mirtest <- import.mirt(mod, display = FALSE)
    opt <- optim(par = est_MLE, fn = IRTlikRcppA, gr= gradIRTlikRcppA, data = data_k, 
                 pen = FALSE, nodes = gq$nodes, weights = gq$weights, hessian = FALSE,
                 control = list(fnscale = -1), method = "BFGS")
    estik <- opt$par
    est[[1]][k,] <- estik
    lik[[1]][k] <- IRTlikRcppA(par = estik, data = data[gr[[k]],], pen = FALSE, nodes = gq$nodes, weights = gq$weights)
  }
  if (trace) message("")
  
  # ===============================
  # penalized
  # ===============================
  
  for (i in 2:nlambda) {
    est[[i]] <- matrix(NA, K, nitems * 3)
    lik[[i]] <- rep(NA, K)
    lambda_i <- lambda[i - 1]
    if (trace) message("lambda = ", formatC(lambda_i, digits = 4, flag = "-") ,"  ", appendLF = FALSE)
    if (trace) message("folds: ", appendLF = FALSE)
    for (k in 1:K) # k = subset
    {
      if (trace) message(k, appendLF = FALSE)
      data_k <- data[-gr[[k]],] # training set
      # estimate for increasing amounts of shrikage
      ini <- est[[i-1]][k,]
      opt1 <- optim(par = ini, fn = IRTlikRcppA, gr = gradIRTlikRcppA, method = "BFGS", data = data_k, pen = TRUE, 
                    lambda = lambda_i, nodes = gq$nodes, weights = gq$weights,
                    control=list(fnscale = -1, maxit = 10000))
      estik <- opt1$par
      est[[i]][k,] <- estik
      lik[[i]][k] <- IRTlikRcppA(par = estik, data = data[gr[[k]],], pen = FALSE, nodes = gq$nodes, weights = gq$weights)
    }
    if (trace) message("")
  }
  # select the maximum likelihood
  sel <- which.max(sapply(lik, mean))
  names(lik) <- names(est) <- paste("lambda =", c(0, lambda))
  return(list(est = est, lik = lik, sel = sel))
}


plot.tpl <- function(x, onlyGuess = TRUE, IRTparam = TRUE, ...)
{
  est_pen <- x$est_pen
  if (is.null(est_pen)) stop("Penalized estimation was not performed.")
  cv_lik <- x$cv_lik
  lambda <- c(0, x$lambda)
  lambda_sel_cv <- x$lambda_sel_cv
  lambda_sel_AIC <- x$lambda_sel_AIC
  lambda_sel_BIC <- x$lambda_sel_BIC
  nitems <- nrow(est_pen) / 3
  
  #if (onlyGuess) par(mfrow = c(1, 2)) else par(mfrow = c(2, 2))
  par(ask = TRUE)
  
  # likelihood
  if (!is.null(lambda_sel_cv))
  {
    lik_mean <- sapply(cv_lik, mean)
    plot(lambda, -lik_mean, xlab = expression(lambda), ylab = "CV error", col = 1, type = "l")
    abline(v = lambda_sel_cv, lty = 2)
  }
  if (!is.null(lambda_sel_AIC))
  {
    aic <- c(x$AIC_MLE, x$AIC_pen)
    plot(lambda, aic, xlab = expression(lambda), ylab = "AIC", col = 1, type = "l")
    abline(v = lambda_sel_AIC, lty = 2)
  }
  if (!is.null(lambda_sel_AIC))
  {
    bic <- c(x$BIC_MLE, x$BIC_pen)
    plot(lambda, bic, xlab = expression(lambda), ylab = "BIC", col = 1, type = "l")
    abline(v = lambda_sel_BIC, lty = 2)
  }
  if (IRTparam)
  {
    est_pen[1:nitems,]<-plogis(est_pen[1:nitems, ]) # guessing
    est_pen[(nitems + 1):(2 * nitems), ]<- -est_pen[(nitems + 1):(2 * nitems),] / est_pen[(nitems * 2 + 1):(3 * nitems), ] #difficulty
  }
  
  # guessing
  est1 <- est_pen[1:nitems, ]
  plot(lambda,est1[1, ], xlim = c(0, max(lambda)), ylim = c(min(est1), max(est1)), 
       xlab = expression(lambda), ylab = "", col = 1, type = "l")
  if (IRTparam) title(ylab = expression(hat(c)[j]), line = 2.4)
  else title(ylab = expression(hat(beta)[3][j]),line = 2.4)
  for (i in 2:nitems) lines(lambda, est1[i, ])
  abline(v = lambda_sel_cv, lty = 2)
  abline(v = lambda_sel_AIC, lty = 2)
  abline(v = lambda_sel_BIC, lty = 2)
  
  if (!onlyGuess) 
    
  {
    # difficulty
    est2 <- est_pen[(nitems + 1):(2 * nitems), ]
    plot(lambda, est2[1, ], xlim = c(0, max(lambda)), ylim = c(min(est2), max(est2)), xlab = expression(lambda), 
         ylab = "", col = 1, type = "l")
    if (IRTparam) title(ylab = expression(hat(b)[j]), line = 2.4)
    else title(ylab = expression(hat(beta)[1][j]), line = 2.4)
    for (i in 2:nitems) lines(lambda, est2[i,])
    abline(v = lambda_sel_cv, lty = 2)
    abline(v = lambda_sel_AIC, lty = 3)
    abline(v = lambda_sel_BIC, lty = 4)
    
    # discrimination
    est3 <- est_pen[(nitems * 2 + 1):(3 * nitems),]
    plot(lambda,est3[1, ], xlim = c(0, max(lambda)), ylim = c(min(est3), max(est3)), xlab = expression(lambda),
         ylab = "", col = 1, type = "l")
    if (IRTparam) title(ylab = expression(hat(a)[j]), line = 2.4)
    else title(ylab = expression(hat(beta)[2][j]),line = 2.4)
    for (i in 2:nitems) lines(lambda, est3[i, ])
    abline(v = lambda_sel_cv, lty = 2)
    abline(v = lambda_sel_AIC, lty = 3)
    abline(v = lambda_sel_BIC, lty = 4)
  }
  par(ask = FALSE)
}




