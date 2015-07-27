




####################################################################
############### Truncated estimators (1st paper) ##################
####################################################################



### estimation of gamma

trHill <- function(data, r = 1, tol = 1e-8, maxiter = 100, plot = FALSE, add = FALSE, 
                         main="Estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  gamma <- numeric(n)
  H <- numeric(n)
  H1 <- numeric(n)
  R <- numeric(n)
  K <- r:(n-1)
  
  ### Truncated Hill estimator 
  
  # Trimmed Hill
  H[K] <- cumsum(log(X[(n-r+1):(n-tail(K,1)+1)]))/(K-r+1) - log(X[n-K]) 
  K1 <- 1:(n-1)
  # Hill
  H1[K] <- (cumsum(log(X[n-K1+1]))/K1 - log(X[n-K1]))[K] 
  
  R[K] <- X[n-K]/X[n-r+1]
  
  ### Newton-Raphson to obtain estimates of gamma
  for (k in (n-1):r) {  
    # Note that y = gamma = 1/alpha
    #y <- H1[k]
    #Use Hill's trimmed estimator as initial approximation
    y <- H[k]
    i <- 1
    while (i <= maxiter){
      #temp[i] <- y
      ym <- 1/y
      z <- y + (H[k]-y-R[k]^ym*log(R[k])/(1-R[k]^ym))/(1-ym^2*R[k]^ym*log(R[k])^2/((1-R[k]^ym)^2))
      
      if (abs(z-y) < tol | identical(abs(z-y),NaN)) { break }
      
      i <- i+1
      y <- z
    }
    gamma[k] <- z
  }
  
  gamma[gamma<=0] <- NA
  
  ### plots if TRUE  
  plotfun(K, gamma[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  
  ### output list with values of k and
  ### corresponding estimates for gamma
  
  output(list(k=K, gamma=gamma[K], H=H[K]),plot=plot,add=add)
  
}


### estimation of DT

trDT <- function(data, r = 1, gamma, plot=FALSE, add=FALSE, 
                      main="Estimates of DT", ...) {
  
  # Check input arguments
  checkInput(data,gamma=gamma,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  DT <- numeric(n)
  K <- r:(n-1)
  
  if(length(gamma)!=length(K)) {
    stop(paste("gamma should have length", length(K)))
  }
  
  
  A <- 1/gamma
  R <- X[n-K]/X[n-r+1]
  
  DT[K] <-  pmax( (K+1)/(n+1) * (R^A-r/(K+1)) / (1-R^A), 0)
  
  ### plots if TRUE  
  plotfun(K, DT[K], type="l", xlab="k", ylab="DT", main=main, plot=plot, add=add, ...)
  
  
  ### output list with values of k and
  ### corresponding estimates for DT
  
  output(list(k=K, DT=DT[K]),plot=plot,add=add)
  
}

trEndpoint <- function(data, r = 1, gamma, DT, plot = FALSE, add = FALSE, 
                      main = "Estimates of Endpoint", ...) {
  
  # Check input arguments
  checkInput(data,gamma=gamma,DT=DT,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  Tk <- numeric(n)
  K <- r:(n-1)
  

  
  if(length(gamma)!=length(K)) {
    stop(paste("gamma should have length", length(K)))
  }
  
  if(length(DT)!=length(K)) {
    stop(paste("DT should have length", length(K)))
  }
  
  
  Tk[K] <-  exp( pmax( log(X[n-K]) + gamma * log(1+(K+1)/((n+1)*DT)), log(X[n])) )
  
  ### plots if TRUE  
  plotfun(K, Tk[K], type="l", xlab="k", ylab="Tk", main=main, plot=plot, add=add, ...)
  
  ### output list with values of k and
  ### corresponding estimates for DT
  
  output(list(k=K, Tk=Tk[K]),plot=plot,add=add)
  
}

### quantile estimation
# rough indicates if we are in the rough truncation case
trQuant <- function(data, r = 1, rough = TRUE, gamma, DT, p, plot = FALSE, add = FALSE,
                         main="Estimates of extreme quantile", ...) {
  
  # Check input arguments
  checkInput(data,gamma=gamma,DT=DT,r=r)
  checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  quant <- numeric(n)
  K <- r:(n-1)
  
  ### Estimator for extreme quantiles

  if (rough) {
    quant[K] <- X[n-K] * ((K+1)/((n+1)*p))^gamma * ((1+DT*(n+1)/(K+1))/(1+DT/p))^gamma
    
  } else {
    quant[K] <- X[n-K] * ((K+1)/(p*(n+1)))^gamma
  }


  
  ### plots if TRUE
  
  plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  ### output list with values of k, corresponding quantile estimates 
  ### and the considered small tail probability p
  
  output(list(k=K, Q=quant[K], p=p),plot=plot,add=add)
  
}



### exceedance probability estimation
trProb <- function(data, r = 1, gamma, q, warnings = TRUE, plot = FALSE, add = FALSE,
                         main="Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  checkInput(data,gamma=gamma,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  prob <- numeric(n)
  K <- r:(n-1)
  
  R <- X[n-K]/X[n-r+1]
  
  ### Estimator for exceedance probability
  prob[K] <- (K+1)/(n+1) * ( (q/X[n-K])^(-1/gamma) - R^(1/gamma) ) / (1-R^(1/gamma))
  
  if (q > max(data) & warnings) {
    warning("q should be smaller than the largest observation.")
  } 

  
  prob[prob<0 | prob>1] <- NA
  
  ### plots if TRUE
  
  if( !all(is.na(prob[K])) ) {
    plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  }
  
  
  ### output list with values of k, exceedance probabilitye estimates 
  ### and the considered high quantile q
  
  output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
  
}


trParetoQQ <- function(data, r = 1, DT, kstar = NULL, main = "TPa QQ-plot") {
  
  # Check input arguments
  checkInput(data,DT=DT)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  j <- 1:n
  K <- r:(n-1) 
  
  if(length(DT)==1) {
    DT = rep(0,length(K))
  }
  

  
  # calculate theoretical and empirical quantiles
    

  pqq.emp <- log(X[n-j+1])
  
  # Select kstar based on maximising correlation
  if(is.null(kstar)) {
    k = ifelse(length(K)>10,10,1):length(K)
    cors <- numeric(length(k)) 
    for(i in 1:length(k)) {
      cors[i] <- abs(cor(log(X[n-(1:k[i])+1]),log(DT[K==k[i]]+(1:k[i])/n)))
    }
    kstar <- k[which.max(cors)]
  } 
  pqq.the <- log(DT[K==kstar]+j/n)

  #Old margins
  oldmar <-  par("mar")
  
  #Increase margins to the left
  par(mar=c(5.1, 4.6, 4.1, 2.1))
  
  # plots if TRUE
  plot(pqq.emp, pqq.the, type="p", xlab=bquote(log(X["n-j+1,n"])), 
       ylab=bquote(log(hat(D)[.(paste0("T,",r,",",kstar))]+j/n)), 
          main=main)
  
  #Reset margins
  par(mar=oldmar)
  
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  # invisible makes sure that the returned list is not printed, but can still be assigned
  invisible(list(pqq.the=pqq.the, pqq.emp=pqq.emp, DT_opt = DT[K==kstar],kstar=kstar))
}


#Test for non-truncated vs. truncated Pareto-type tails
trTest <- function(data, alpha = 0.05, plot = TRUE, main = "Test for truncation", ...) {
  
  # Check input arguments
  checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  K <- 1:(n-1)
  
  H <- Hill(X)$gamma
  
  #Compute L(H_{k,n})
  E <- numeric(n)

  for(k in K) {
    E[k] <- sum( (X[n-k]/X[n-(1:k)+1])^(1/H[k]) ) / k
  }  
    
  L <- (E[K]-0.5)/(1-E[K])
  
  #Test value
  tv <- sqrt(12*K)*L
  
  #Critical value
  cv <- qnorm(alpha)
  
  #Reject
  reject <- (tv<=cv)
  
  #P-values
  Pval <- pnorm(tv)
  
  #Plot P-values 
  
  plotfun(K, Pval[K], ylim=c(0,1), type="l", xlab="k", ylab="P-value", main=main, plot=plot, add=FALSE, ...)
  if (plot) abline(h=alpha, col="blue", ...)
  
  output(list(k=K, testVal=tv, critVal=cv, Pval=Pval, reject=reject), 
         plot=plot, add=FALSE)
  
}


# ####################################################################
# ################ Truncated estimators (proposal) ###################
# ####################################################################
# 
# ### estimation of DT
# 
# tr.DT <- function(data, eps=1e-8, plot=FALSE, add=FALSE, 
#                   main="Estimates of DT", ...) {
#   
#   # Check input arguments
#   checkInput(data)
#   
#   X <- as.numeric(sort(data))
#   n <- length(X)
#   
#   DT <- numeric(n)
#   H <- numeric(n)
#   EH <- numeric(n)
#   K <- 1:(n-1)
#   
# 
#   H[K] <- cumsum(log(X[n-K+1])) / K - log(X[n-K])
#   K2 <- 1:(n-2) 
#   EH[K2] <- (1/K2)*cumsum(H[K2])/H[K2+1]
#   
# #   EH2 <- numeric(n)
# #   for (k in 1:(n-1)) {
# #       Z1.values <- log(X[n:(n-k+1)]/X[n-k])
# #       H[k] <- (1/k)*sum(Z1.values)
# #       EH2[k-1] <- (1/(k-1))*sum(H[1:(k-1)])/H[k]
# #   }
# #   print(max(abs(EH2-EH)/EH2,na.rm=TRUE))
#   
#   ### estimates of DT
#   f <- function(dt,eh,n_f,k_f) {
#     logtemp <- log(1+(k_f+1)/(n_f*dt))
#     sumtemp <- sum(log(1+(1:k_f)/(n_f*dt))/(1:k_f))
#     return((1-n_f*dt/(k_f+1)*sumtemp)/(1-n_f*dt/(k_f+1)*logtemp) - eh)
#   }  
#   
#   for (k in (n-2):1) {  
#     
#     end = 2
#     if(sign(f(eps,EH[k],n,k))!=sign(f(end,EH[k],n,k))){
#       DT[k] <- uniroot(f, interval=c(eps,end), eh=EH[k], n_f=n, k_f=k)$root
#     }
#   }
#   
#   plotfun(K, DT[K], type="l", xlab="k", ylab="DT", main=main, plot=plot, add=FALSE,...)
#   
#   ### output list with values of k and
#   ### corresponding estimates for DT
#   
#   output(list(k=K, DT=DT[K]),plot=plot,add=add)
#   
# }
# 
# 
# ### estimation of beta
# 
# tr.beta <- function(data, DT, plot=FALSE, add=FALSE, 
#                     main="Estimates of beta", ...) {
#   
#   # Check input arguments
#   checkInput(data)
#   
#   X <- as.numeric(sort(data))
#   n <- length(X)
#   
#   Dt <- numeric(n)
#   beta <- numeric(n)
#   H <- numeric(n)
#   K <- 1:(n-1)
#   
#   
#   if(length(DT)!=length(K)) {
#     stop(paste("DT should have length", length(K)))
#   }
#   
#   
#   Dt[K] <- DT
#   
#   ### estimates of beta
#   H[K] <- cumsum(log(X[n-K+1])) / K - log(X[n-K])
#   K2 <- 1:(n-2)
#   beta[K2] <- H[K2]/(1-n*Dt[K2]/(K2+1)*log(1+(K2+1)/(n*Dt[K2])))
#   
#   
#   plotfun(K, DT[K], type="l", xlab="k", ylab=bquote(beta), main=main, plot=plot, add=FALSE,...)
#  
#   ### output list with values of k and
#   ### corresponding estimates for DT
#   
#   output(list(k=K, beta=beta[K]),plot=plot,add=add)
#   
# }
# 
# 
# ### estimation of gamma
# 
# tr.gamma <- function(data, beta, DT, plot=FALSE, add=FALSE, ...) {
#   
#   X <- sort(data)
#   # Check input arguments
#   checkInput(data)
#   
#   X <- as.numeric(sort(data))
#   n <- length(X)
#   
#   Dt <- numeric(n)
#   Beta <- numeric(n)
#   gamma <- numeric(n)
#   R <- numeric(n)
#   E <- numeric(n)
#   K <- 1:(n-1)
#   
#   Dt[K] <- DT
#   Beta[K] <- beta
#   R[K] <- X[n-K+1]/X[n-K]
#   
#   ### obtain estimates of gamma
#   h <- function(x,gamma) {
#     ifelse(gamma==0,log(x),(x^gamma-1)/gamma)
#   }
#   
#   
#   f <- function(gamma,dt,beta,e,n_f,k_f) {
# 
#     sumtemp <- sum(h(((n_f)*dt)/((n_f)*dt+(1:k_f)+1),gamma))/k_f
#     temp <- h(((n_f)*dt)/((n_f)*dt+k_f+1),gamma)
#     f <- (1+beta*sumtemp)/(1+beta*temp) - e
#   }
#   
#   for (k in (n-2):1) {
#     if (!identical(Beta[k],NaN)) {
#       E[k] <- sum(R[1:k])/k
#       
#       if(sign(f(0,Dt[k],Beta[k],E[k],n,k))!=sign(f(25,Dt[k],Beta[k],E[k],n,k))){
#         gamma[k] <- uniroot(f, interval=c(0,25), dt=Dt[k], beta=Beta[k], e=E[k], n_f=n, k_f=k)$root
#       }
#     }
#     else {
#       gamma[k] = NaN
#     }
#   }
#   
#   ### plots if TRUE  
#   if (plot || add){
#     if ( !add ) {   ### plot estimates
#       plot(K, gamma[K], type="l", ylab="alpha", xlab="k", main="Estimates of gamma", ...)
#     }
#     else { ### adds estimates to egammasting plot
#       lines(K, gamma[K], ...)
#     }
#   }
#   
#   ### output list with values of k and
#   ### corresponding estimates for DT
#   
#   list(k=K, gamma=gamma[K])
#   
# }
# 
# 
# tr.quant <- function(data, gamma, beta, DT, p, plot=FALSE, add=FALSE, ...) {
#   
#   # Check input arguments
#   checkInput(data)
#   checkInput(p)
#   
#   X <- as.numeric(sort(data))
#   n <- length(X)
#   
#   Dt <- numeric(n)
#   gamma <- numeric(n)
#   Beta <- numeric(n)
#   quant <- numeric(n)
#   K <- 1:(n-1)
#   
#   gamma[K] <- gamma
#   Dt[K] <- DT
#   Beta[K] <- beta
#   
#   ### Estimator for extreme quantiles
#   
#   h <- function(x,gamma) {
#     ifelse(gamma==0,log(x),(x^gamma-1)/gamma)
#   }
#   
#   K2 <- K[!identical(gamma[K],NaN)]
#   quant[K2] <- X[n-K2]*(1+Beta[K2]*h(Dt[K2]/(Dt[K2]+p),gamma[K2]))/(1+Beta[K2]*h(n*Dt[K2]/(n*Dt[K2]+K2),gamma[K2]))
#   
#   
# #   h <- function(x,gamma) {
# #     h <- (x^gamma-1)/gamma
# #   }
# #   quant2 <- numeric(n)
# #   for (k in (n-1):1) {
# #     if (1-identical(gamma[k],NaN)) {
# #       if (gamma[k]!=0) {
# #         #quant[k] <- X[n-k]*(gamma[k]+Beta[k]*((Dt[k]/(Dt[k]+p))^gamma[k]-1))/(gamma[k]+Beta[k]*(((n+1)*Dt[k]/((n+1)*Dt[k]+k+1))^gamma[k]-1))
# #         quant2[k] <- X[n-k]*(1+Beta[k]*h(Dt[k]/(Dt[k]+p),gamma[k]))/(1+Beta[k]*h(n*Dt[k]/(n*Dt[k]+k),gamma[k]))
# #       }
# #       else {
# #         quant2[k] <- X[n-k]*(1+Beta[k]*log(Dt[k]/(Dt[k]+p)))/(1+Beta[k]*log(n*Dt[k]/(n*Dt[k]+k)))
# #       }
# #     }  
# #   }
# #   print(max(abs(quant-quant2),na.rm=TRUE))
# 
#   
#   ### plots if TRUE
#   
#   if (plot || add){
#     if ( !add ) {     ### plot estimates
#       plot(K, quant[K], type="l", ylab="Q(1-p)", xlab="k", main="Estimates of extreme quantile", ...)
#     }
#     else {     	### adds estimates to egammasting plot
#       lines(K, quant[K], ...)
#     }
#   }
#   
#   ### output list with values of k, corresponding quantile estimates 
#   ### and the considered small tail probability p
#   
#   list(k=K, Q=quant[K], p=p)
#   
# }



