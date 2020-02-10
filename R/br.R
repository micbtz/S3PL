gendata3PL <- function(beta, S)
{
  I <- length(beta) / 3
  parameters <-  list(alpha = beta[1:I], beta=beta[1:I + I], gamma = beta[1:I + 2 * I])
  data.y <- matrix(0, nrow=S, ncol=I)
  theta <- rnorm(S)
  for(i in 1:I)
         {
          ci <- plogis(parameters$gamma[i])
	      probi <- ci + (1 - ci) * plogis(parameters$alpha[i] + theta * parameters$beta[i])
	      data.y[,i] <- rbinom(S, size = 1, prob = probi)
	      }
   return(data.y)      
}

tr <- function(x) sum(diag(x))

A_star_pack <- function(para, y, gh, ncores)
{
  S <- nrow(y)
  I <- ncol(y)
  p <- 3 * I
  parameters <-  list(alpha = para[1:I], beta = para[1:I+I], gamma = para[1:I+2*I]) 
  a <- rep(0, p)
  listy <- splitdata(y, ncores) 
  myf <- function(v) {
     mat <- matrix(v, ncol = I, byrow = FALSE)
     E <- array(0, dim=c(p, p, p))
     B.size <- nrow(mat)
     out1 <- matrix(0, B.size, p) 
     for(i in 1:B.size){
     data.list1 <- list(y = matrix(mat[i,], nrow = 1), y1 = matrix(1-mat[i,], nrow = 1), 
                      nodes = as.vector(gh$nodes), weights = as.vector(gh$weights / sqrt(pi)))
     uni1 <- TMB::MakeADFun(data = data.list1, parameters = parameters, DLL = "BR3PL", silent = TRUE) 
     out1[i,] <- -uni1$gr(para)
     out2 <- uni1$he(para)
     t1 <- tcrossprod(out1[i,])
     JS <- out2 - t1
     E <- E + outer(out1[i,], JS) 
     }
     list(E = E, out1 = out1)
    }
  ogg <-  parallel::mclapply(listy, myf, mc.cores = ncores)
  startind <- 0
  E <- array(0, dim=c(p, p, p))
  out1 <- matrix(0, S, p)
  for(i in 1:length(ogg)) 
    {
        B.size <- nrow(ogg[[i]]$out1)
        ind <- startind + 1:B.size
        startind <- startind + B.size 
        E <- E + ogg[[i]]$E / S 
        out1[ind,] <- ogg[[i]]$out1
   } 
  J <- cov(out1)
  F1 <- try(solve(J)) 
  if(is.numeric(F1)) for(j in 1:p) a[j] <- -0.5 * tr(F1 %*% E[,,j]) 
  return(list(a = a, J = J))
 } 
  
  
 
##split the data in ncores groups  
splitdata <- function(y, ncores)
{
  M <- nrow(y)
  rest <- M%%ncores   
  B.size <- floor(M / ncores)
  ind <-  if(ncores>1) sort(c(rep(1:(ncores-1), B.size), rep(ncores, B.size+rest)))
           else rep(1,B.size)
  listy <- split(y, ind)
  return(listy)
}  
  


exp.info1 <- function(para, yi, gh, ncores)
{
  listy <- splitdata(yi, ncores)
  M <- nrow(yi)
  I <- ncol(yi)
  p <- 3 * I
  parameters <-  list(alpha = para[1:I], beta=para[1:I + I], gamma = para[1:I + 2 * I]) 
  myf <- function(v) {
     mat <- matrix(v, ncol=I, byrow=FALSE)
     out1 <- matrix(0, nrow(mat), p) 
     for(i in 1:nrow(mat)){
     data.list1 <- list(y=matrix(mat[i,], nrow = 1), y1=matrix(1 - mat[i,], nrow = 1), 
                      nodes=as.vector(gh$nodes), weights=as.vector(gh$weights))
     uni1 <- TMB::MakeADFun(data = data.list1, parameters = parameters, DLL = "BR3PL", silent = TRUE) 
     out1[i,] <- -uni1$gr(para)
     }
     list(out1=out1)
  }
  ogg <-  parallel::mclapply(listy, myf, mc.cores = ncores)
  out1 <- matrix(0, M, p)
  startind <- 0
  for(i in 1:length(ogg)) 
    {
    	B.size <- nrow(ogg[[i]]$out1)
        ind <- startind + 1:B.size
        startind <- startind + B.size 
        out1[ind,] <- ogg[[i]]$out1
   } 
  J <- crossprod(out1) / M 
  return(J)
 } 


gsolv_pack <- function(para, y, gh, ncores, ADobj, M, seed=NULL)
{
  if(!is.null(seed)){ 
         set.seed(seed)
         yi <- gendata3PL(para, M)
      } else yi <- y  
  a <- A_star_pack(para, yi, gh, ncores)$a 
  score <- as.vector(-ADobj$gr(para))
  out <- score + a 
  return(out) 
}

 
### uses exp info 
jsolv <- function(para, y, gh, ncores, ADobj, M, seed=NULL)
{
  if(!is.null(seed)){ 
         set.seed(seed)
         yi <- gendata3PL(para, M)
      } else yi <- y  
  out <- exp.info1(para, yi, gh, ncores) * nrow(y) #####
  return(-out)
 }
