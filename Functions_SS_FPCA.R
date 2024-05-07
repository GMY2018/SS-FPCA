## -----------------------------------------------------------
## Functions for the AECM algorithm for implementing SS-FPCA
## -----------------------------------------------------------


library(sm)
library(splines)
library(nlme)
library(matrixcalc)
library(Matrix)


## -----------------------------------------------------
## PART A : the cycle for state space component
## -----------------------------------------------------

## (I) A function for computing the score function of sigma
score.sigma <- function(sigma, xi, Lambda, B, N){
  # sigma is essentially sigma^2
  # B can be obtained after the smoothing is done, see Xu & Wikle (2007)
  # N is the number of time points (images)
  # n is the number of observations per time point
  n <- nrow(B)   
  Score <- 0
  Heissian <- 0
  for (i in 1:ncol(xi)) {
    Ai <- xi[,i]%*%t(xi[,i])
    Score <- Score + N*(1/(Lambda[i,i] + sigma) - 1/sigma)*sum(diag(Ai)) -
      (1/(Lambda[i,i] + sigma)^2 - 1/(sigma^2))*sum(diag(Ai%*%B))
    Heissian <- Heissian - N*(1/(Lambda[i,i] + sigma)^2 - 1/(sigma^2))*sum(diag(Ai)) + 
      (1/(Lambda[i,i] + sigma)^3 - 1/(sigma^3))*sum(diag(Ai%*%B))
  }
  # the Score and Heissian
  Score <- Score + N*n/sigma - sum(diag(B))/(sigma^2)
  Heissian <- Heissian - N*n/(sigma^2) + sum(diag(B))/(sigma^3)
  list(Score=Score, Heissian=Heissian)
}


## (II) Functions for computing the likelihood for the Kalman filter
## 1. For system states
AECM.KFa.lik <- function(Phi, KS, Plag.s, mu0, Sigma0, HHt) {
  ## preparations
  Phi <- as.matrix(Phi)
  a.s <- KS$xs
  P.s <- KS$Ps
  a0 <- KS$x0n
  P0 <- KS$P0n
  
  N <- ncol(a.s)
  df <- nrow(a.s)
  
  ## The likelihood w.r.t. the system states
  S00 <- S01 <- S11 <- matrix(0, nrow=df, ncol=df)
  cHt <- chol(HHt)
  HHt.inv <- solve(cHt)%*%t(solve(cHt))
  
  # The smoothed state from time point 2 to N
  for (i in 2:N){
    at0 <- as.matrix(a.s[,i-1])
    at1 <- as.matrix(a.s[,i])
    Pt0 <- as.matrix(P.s[,,i-1])
    Pt1 <- as.matrix(P.s[,,i])
    Pt01 <-  as.matrix(Plag.s[,,i])
    
    S11 <- S11 + Pt1 + at1%*%t(at1)
    S01 <- S01 + Pt01 + at1%*%t(at0)
    S00 <- S00 + Pt0 + at0%*%t(at0)
  }
  
  # Add on the first smoothed state (T = 1)
  S11 <- S11 + as.matrix(P.s[,,1]) + as.matrix(a.s[,1])%*%t(as.matrix(a.s[,1]))
  S01 <- S01 + as.matrix(Plag.s[,,1]) + as.matrix(a.s[,1])%*%t(a0)
  S00 <- S00 + P0 + a0%*%t(a0)
  
  # The likelihood component
  at.lik <- S11 - S01%*%t(Phi) - Phi%*%t(S01) + Phi%*%S00%*%t(Phi)
  a0.lik <- P0 + (a0 - mu0)%*%t(a0 - mu0)
  at.sum <- sum(diag(at.lik%*%HHt.inv)) + 2*N*sum(log(diag(cHt)))
  a0.sum <- sum(diag(a0.lik%*%solve(Sigma0))) + sum(log(diag(Sigma0)))
  
  list(a0.sum=a0.sum, at.sum=at.sum, S00=S00, S01=S01, S11=S11)
}


## 2. For observation
# will need to compute 'GGt' and 'GGt.inv' first, using
# temp.inv <- diag(sigma/diag(Lambda) + 1)
# GGt.inv <- 1/sigma*diag(n) - xi%*%temp.inv%*%t(xi)
AECM.KFy.lik <- function(y, A, Phi, KS, GGt, GGt.inv, Pmiss=1) {
  ## preparations
  y.mat <- as.matrix(y)
  Phi <- as.matrix(Phi)
  cGt <- chol(GGt)  # this is not diagonal
  a.s <- KS$xs
  P.s <- KS$Ps
  
  n <- nrow(y)
  N <- ncol(y)
  df <- ncol(A)
  
  Yt.lik <- matrix(0, nrow=n, ncol=n)
  
  for (i in 1:N){
    whereNA <- is.na(y.mat[,i])
    PNA <- sum(whereNA)/length(whereNA) 
    
    if (!any(whereNA)) {
      rt <- y.mat[,i] - A %*% as.matrix(a.s[,i])
      Yt.lik <- Yt.lik + A %*% as.matrix(P.s[,,i]) %*%t(A) + rt%*%t(rt)
    }
    else {
      if((all(whereNA))|(PNA >= Pmiss)){
        Yt.lik <- Yt.lik + 0
      }
      else {
        yes <- as.vector(!whereNA)
        no <- as.vector(whereNA) 
        nyes <- sum(yes)
        
        y.obs <- rep(0, times=n)
        y.obs[1:nyes] <- y.mat[yes,i]
        A.obs <- matrix(0, nrow=n, ncol=df)
        A.obs[1:nyes,] <- A[yes,]
        G.mis <- matrix(0, nrow=n, ncol=n)
        G.mis[(nyes+1):n,(nyes+1):n] <- GGt[no,no]
        
        rt <- y.obs - A.obs%*%as.matrix(a.s[,i])
        Yt.lik <- Yt.lik + A.obs%*%as.matrix(P.s[,,i])%*%t(A.obs) + rt%*%t(rt) + G.mis
      }
    }
  }
  
  Yt.sum <- sum(diag(Yt.lik %*% GGt.inv)) + 2*N*sum(log(diag(cGt)))
  
  list(Yt.lik=Yt.lik, Yt.sum=Yt.sum)
}


## (III) the filtering/smoothing functions
## 1. The Kalman filter with threshold
sparse.Pfilter.G <- function (y, A, mu0, Sigma0, Phi, cHt, Theta, Lambda, sigma, Pmiss, NA.mat) {
  ## (1) build the model
  # the observation equation
  y.mat <- as.matrix(y)
  # NA.mat is the matrix indicating the locations of missing observations
  n <- nrow(y.mat)
  N <- ncol(y.mat)  
  
  EOF <- Theta %*% Lambda %*% t(Theta)
  G <- A %*% EOF %*% t(A) + sigma*diag(n)
  # cGt <- chol(G)
  # cGt.inv <- solve(cGt)
  # G.inv <- cGt.inv %*% t(cGt.inv)
  
  # the state transition equation
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)
  H <- t(cHt) %*% cHt
  cHt.inv <- solve(cHt)
  H.inv <- cHt.inv %*% t(cHt.inv)
  
  # the initial state
  x00 <- as.matrix(mu0, nrow=df, ncol=1)
  P00 <- as.matrix(Sigma0, nrow=df, ncol=df)
  
  ## (2) For storing the results
  xp <- array(NA, dim=c(df, N))
  Pp <- array(NA, dim=c(df, df, N))
  xf <- array(NA, dim=c(df, N+1))
  Pf <- array(NA, dim=c(df, df, N+1))
  xf[,1] <- x00
  Pf[,,1] <- P00
  
  ## (3) Begin the filter!!
  miss <- NULL
  st <- 0
  
  for (i in 1:N) {
    whereNA <- NA.mat[,i]
    PNA <- sum(whereNA)/n 
    miss <- c(miss, PNA >= Pmiss)
    
    if (!any(whereNA)) {
      # NO missing
      xp[,i] <- Phi %*% xf[,i]
      Pptemp <- Phi %*% Pf[,,i] %*% t(Phi) + H
      Pp[,,i] <- (t(Pptemp) + Pptemp)/2
      # The 'sigtemp <- A %*% Pp[,,1] %*% t(A) + G' needs to be 
      # changed by using Sherman-Morrison-Woodbury
      # Will also need to adapt to the FPCA component
      # May also considering Pseudo Inverse
      
      P.temp <- Pp[,,i] + EOF
      P.inv <- svd.inverse(P.temp) * sigma
      # P.inv <- svd.inverse(Pp[,,i])  for the standard filtering function
      gain.inv <- svd.inverse(P.inv + diag(df))
      gain <- (diag(n) - A %*% gain.inv %*% t(A)) / sigma  
      # Using the Woodbury identity
      K <- Pp[,,i] %*% t(A) %*% gain
      
      innov <- y.mat[,i] - A %*% xp[,i]
      xf[,i+1] <- xp[,i] + K %*% innov
      Pftemp <- Pp[,,i] - K %*% A %*% Pp[,,i]
      Pf[,,i+1] <- (t(Pftemp) + Pftemp)/2
    }
    else {
      if ((all(whereNA))|(PNA >= Pmiss)) {
        # ALL missing or smaller than a threshold
        xp[,i] <- Phi %*% xf[,i]
        Pp[,,i] <- Phi %*% Pf[,,i] %*% t(Phi) + H
        xf[,i+1] <- xp[,i]
        Pf[,,i+1] <- Pp[,,i]
      }
      else {
        # Partly missing and greater than a threshold
        yes <- as.vector(!whereNA)
        no <- as.vector(whereNA) 
        nyes <- sum(yes)
        
        y.obs <- rep(0, times=n)
        y.obs[1:nyes] <- y.mat[yes,i]
        A.obs <- matrix(0, nrow=n, ncol=df)
        A.obs[1:nyes,] <- A[yes,]
        
        xp[,i] <- Phi %*% xf[,i]
        Pptemp <- Phi %*% Pf[,,i] %*% t(Phi) + H
        Pp[,,i] <- (t(Pptemp) + Pptemp)/2
        
        # Will this be faster?? 
        P.temp <- Pp[,,i] + EOF
        P.inv <- svd.inverse(P.temp) * sigma
        
        # the adjustment used here is
        I.mis.inv <- diag(n)
        adjust <- G[no,no] / sigma
        adj.inv <- solve(chol(adjust))
        I.adj.inv <- adj.inv %*% t(adj.inv)
        I.mis.inv[(nyes+1):n, (nyes+1):n] <- I.adj.inv
        
        AI.mis <- t(A.obs) %*% I.mis.inv   # for simplification
        # sigtemp <- (P.inv + t(A.obs) %*% I.mis.inv %*% A.obs)
        sigtemp <- (P.inv + AI.mis %*% A.obs)
        sig <- (t(sigtemp) + sigtemp)/2
        siginv <- solve(sig)
        
        # gain <- (I.mis.inv - I.mis.inv %*% A.obs %*% siginv %*% t(A.obs) %*% I.mis.inv) / sigma
        gain <- (I.mis.inv - t(AI.mis) %*% siginv %*% AI.mis) / sigma
        # Using the Woodbury identity
        K <- Pp[,,i] %*% t(A.obs) %*% gain
        
        innov <- y.obs - A.obs %*% xp[,i]
        xf[,i+1] <- xp[,i] + K %*% innov
        Pftemp <- Pp[,,i] - K %*% A.obs %*% Pp[,,i]
        Pf[,,i+1] <- (t(Pftemp) + Pftemp)/2
      }
    }
    
    # The n-step BACKWARD SMOOTHING!
    if ((i >= 2) & (miss[i]==0)) {
      if (miss[i-1]==1){
        # apply the smoother
        smallKS <- sparse.Back(y.mat[,(i-st):i], A=A, Phi=Phi, H=H, xf=xf[,(i-st+1):(i+1)], 
                               Pf=Pf[,,(i-st+1):(i+1)], step=st)
        xf[,(i-st+1):i] <- smallKS$xs
        Pf[,,(i-st+1):i] <- smallKS$Ps
        xp[,(i-st+1):i] <- smallKS$xp
        Pp[,,(i-st+1):i] <- smallKS$Pp
        # start the new counting
        st <- 0 
      }
    }
    else st <- st + 1
  }
  
  ## The final results
  list(xp = xp, Pp = Pp, xf = xf[,2:(N+1)], Pf = Pf[,,2:(N+1)], 
       innov = innov, K = K)
}


## 2. The smoother corresponding to the Kalman filter
sparse.Psmooth.G <- function(y, A, mu0, Sigma0, Phi, cHt, Theta, Lambda, sigma, Pmiss, NA.mat) {
  y.mat <- as.matrix(y)
  n <- nrow(y.mat)
  N <- ncol(y.mat)
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)
  
  ## run the filter
  kf <- sparse.Pfilter.G(y=y.mat, A=A, mu0=mu0, Sigma0=Sigma0, Phi=Phi, cHt=cHt, Theta=Theta, 
                         Lambda=Lambda, sigma=sigma, Pmiss=Pmiss, NA.mat=NA.mat)
  
  ## The smoother
  xs <- array(NA, dim=c(df, N))
  Ps <- array(NA, dim=c(df, df, N))
  J <- array(NA, dim=c(df, df, N))
  
  xs[,N] <- kf$xf[,N]
  Ps[,,N] <- kf$Pf[,,N]
  for (i in N:2) {
    J[,,i-1] <- kf$Pf[,,i-1] %*% t(Phi) %*% solve(kf$Pp[,,i])
    xs[,i-1] <- kf$xf[,i-1] + J[,,i-1] %*% (xs[,i] - kf$xp[,i])
    Pstemp <- kf$Pf[,,i-1] + J[,,i-1] %*% (Ps[,,i] - kf$Pp[,,i]) %*% t(J[,,i-1])
    Ps[,,i-1] <- (t(Pstemp) + Pstemp)/2
  }
  
  x00 <- mu0
  P00 <- Sigma0
  J0 <- as.matrix((P00 %*% t(Phi)) %*% solve(kf$Pp[,,1]), nrow=df, ncol=df)
  x0n <- as.matrix(x00 + J0 %*% (xs[,1] - kf$xp[,1]), nrow=df, ncol=1)
  P0n <- P00 + J0 %*% (Ps[,,1] - kf$Pp[,,1]) %*% t(J0)
  
  ## The results
  list(xs = xs, Ps = Ps, x0n = x0n, P0n = P0n, J0 = J0, J = J, 
       xp = kf$xp, Pp = kf$Pp, xf = kf$xf, Pf = kf$Pf, K = kf$K)
}



## 3. The n-step backward smoother
sparse.Back <- function(y, A, Phi, H, xf, Pf, step=1) {
  # throw in the data from missing time points, t, t-1, ...
  # also the data from t+1
  y.mat <- as.matrix(y)
  n <- nrow(y.mat)
  N <- ncol(y.mat)
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)
  
  if (N != (step+1)) print('Incorrect number of steps')
  else if (N == (step+1)) {
    ## The n-step back smoothing (n >= 1)
    xs <- array(NA, dim=c(df, N))
    Ps <- array(NA, dim=c(df, df, N))
    J <- array(NA, dim=c(df, df, step))
    xp <- array(NA, dim=c(df, step))
    Pp <- array(NA, dim=c(df, df, step))
    
    # note that xf = xp and Pf = Pp for missing data
    xs[,N] <- xf[,N]
    Ps[,,N] <- Pf[,,N]
    for (i in N:2) {
      Pp.old <- Phi %*% Pf[,,(i-1)] %*% t(Phi) + H
      J[,,i-1] <- Pf[,,i-1] %*% t(Phi) %*% solve(Pp.old)
      xs[,i-1] <- xf[,i-1] + J[,,i-1] %*% (xs[,i] - Phi%*%xf[,i-1])
      Pstemp <- Pf[,,i-1] + J[,,i-1] %*% (Ps[,,i] - Pp.old) %*% t(J[,,i-1])
      Ps[,,i-1] <- (t(Pstemp) + Pstemp)/2
      xp[,i-1] <- Phi %*% xs[,i-1]
      Pptemp <- Phi %*% Ps[,,i-1] %*% t(Phi) + H
      Pp[,,i-1] <- (t(Pptemp) + Pptemp)/2
    }
  }
  
  ## The results, excluding t+1
  list(xs = xs[,-N], Ps = Ps[,,-N], xp=xp, Pp=Pp)
}


## 4. The lag-1 smoothed covariance
sparse.Plag1 <- function(y, A, Phi, Pf, J, J0, K){
  y.mat <- as.matrix(y)
  n <- nrow(y.mat)
  N <- ncol(y.mat)
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)
  
  # backward recursion
  Plag <- array(NA, dim=c(df, df, N))
  whereNA <- is.na(y.mat[,N])
  yes <- as.vector(!whereNA)
  no <- as.vector(whereNA) 
  D.mat <- rbind(diag(n)[yes,], diag(n)[no,])
  A.temp <- A
  A.temp[no,] <- 0
  A.obs <- D.mat%*%A.temp
  
  # the following doesn't work if nyes=0
  # nyes <- sum(yes)
  # A.obs <- matrix(0, nrow=n, ncol=df)
  # A.obs[1:nyes,] <- A[yes,]
  
  Plag[,,N] <- (diag(df) - K %*% A.obs) %*% Phi %*% Pf[,,N-1]
  for (i in N:3) {
    # from N-1 to 2
    Plag[,,i-1] <- Pf[,,i-1] %*% t(J[,,i-2]) + J[,,i-1] %*% (Plag[,,i] -                                                               Phi %*% Pf[,,i-1]) %*% t(J[,,i-2])
  }
  Plag[,,1] <- Pf[,,1] %*% t(J0) + J[,,1] %*% (Plag[,,2] - Phi %*% Pf[,,1]) %*% t(J0)
  
  return(Plag)
}



# 5. The function for compute the MLEs of the state space component
sparse.MLE.H <- function(y, A, Phi, mu0, Sigma0, cHt, Theta, Lambda, sigma, Pmiss, NA.mat) {
  y.mat <- as.matrix(y)
  n <- nrow(y.mat)    # this is qdim
  N <- ncol(y.mat)    # this is num
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)     # this is pdim
  
  ## (1) The filter and smoother
  KS <- sparse.Psmooth.G(y=y, A=A, Phi=Phi, mu0=mu0, Sigma0=Sigma0, cHt=cHt, Theta=Theta, 
                         Lambda=Lambda, sigma=sigma, Pmiss=Pmiss, NA.mat=NA.mat)
  a.f <- KS$xf
  P.f <- KS$Pf
  a.s <- KS$xs
  P.s <- KS$Ps
  a.p <- KS$xp
  P.p <- KS$Pp
  Jmat <- KS$J
  J0 <- KS$J0
  Kgain <- KS$K   # only the last Kalman gain
  
  # the lag-1 covariance smoother
  Plag.s <- sparse.Plag1(y=y, A=A, Phi=Phi, Pf=P.f, J=Jmat, J0=J0, K=Kgain)
  
  # the initial state
  a00 <- as.matrix(mu0, nrow=df, ncol=1)
  P00 <- as.matrix(Sigma0, nrow=df, ncol=df)
  
  ## (2) The elements of the likelihood
  # the likelihood of the smoothed state
  S00 <- S01 <- S11 <- matrix(0, nrow=df, ncol=df)
  cHt.inv <- solve(cHt)
  HHt.inv <- cHt.inv %*% t(cHt.inv)
  
  for (i in 2:N){
    # compute the smoothed states and their variance
    at0 <- as.matrix(a.s[,i-1])
    at1 <- as.matrix(a.s[,i])
    Pt0 <- as.matrix(P.s[,,i-1])
    Pt1 <- as.matrix(P.s[,,i])
    Pt01 <-  as.matrix(Plag.s[,,i])
    
    S11 <- S11 + Pt1 + at1%*%t(at1)
    S01 <- S01 + Pt01 + at1%*%t(at0)
    S00 <- S00 + Pt0 + at0%*%t(at0)
  }
  
  # Add on the first smoothed state (T = 1)
  # for this, we need the initial states (T = 0)
  a0 <- KS$x0n
  P0 <- KS$P0n
  S11 <- S11 + as.matrix(P.s[,,1]) + as.matrix(a.s[,1])%*%t(as.matrix(a.s[,1]))
  S01 <- S01 + as.matrix(Plag.s[,,1]) + as.matrix(a.s[,1])%*%t(a0)
  S00 <- S00 + P0 + a0%*%t(a0)
  
  at.lik <- S11 - S01%*%t(Phi) - Phi%*%t(S01) + Phi%*%S00%*%t(Phi)
  a0.lik <- P0 + (a0-a00)%*%t(a0-a00)
  at.sum <- sum(diag(at.lik%*%HHt.inv)) + 2*N*sum(log(diag(cHt)))
  a0.sum <- sum(diag(a0.lik%*%solve(Sigma0))) + sum(log(diag(Sigma0)))
  
  xi <- A %*% Theta
  lambda <- diag(Lambda)
  eta <- lambda/(sigma*(sigma + lambda))
  GGt <- xi %*% Lambda %*% t(xi) + sigma*diag(n)
  cGt <- chol(GGt)
  GGt.inv <- 1/sigma*diag(n) - xi %*% diag(eta) %*%t(xi)
  Yt.lik <- matrix(0, nrow=n, ncol=n)
  
  for (i in 1:N){
    whereNA <- NA.matrix[,i]
    PNA <- sum(whereNA)/length(whereNA) 
    
    if (!any(whereNA)) {
      rt <- y.mat[,i] - A %*% as.matrix(a.s[,i])
      Yt.lik <- Yt.lik + A %*% as.matrix(P.s[,,i]) %*%t(A) + rt%*%t(rt)
    }
    else {
      if((all(whereNA))|(PNA >= Pmiss)){
        Yt.lik <- Yt.lik + 0
      }
      else {
        yes <- as.vector(!whereNA)
        no <- as.vector(whereNA) 
        nyes <- sum(yes)
        
        y.obs <- rep(0, times=n)
        y.obs[1:nyes] <- y.mat[yes,i]
        A.obs <- matrix(0, nrow=n, ncol=df)
        A.obs[1:nyes,] <- A[yes,]
        G.mis <- matrix(0, nrow=n, ncol=n)
        G.mis[(nyes+1):n,(nyes+1):n] <- GGt[no,no]
        
        rt <- y.obs - A.obs%*%as.matrix(a.s[,i])
        Yt.lik <- Yt.lik + A.obs%*%as.matrix(P.s[,,i])%*%t(A.obs) + rt%*%t(rt) + G.mis
      }
    }
  }
  Yt.sum <- sum(diag(Yt.lik%*%GGt.inv)) + 2*N*sum(log(diag(cGt)))
  
  # So the log-likelihood is
  like <- 0.5*(a0.sum + at.sum + Yt.sum)
  
  # (3) The MLEs (analytical solution)
  # Phi.mle <- S01%*%solve(S00)
  # Htemp <- 1/N * (S11 - S01%*%solve(S00)%*%t(S01))
  Htemp <- 1/N * (S11 - 2*S01 + S00)
  H.mle <- (t(Htemp) + Htemp)/2
  
  list(H.mle=H.mle, loglike=like)
}


## If the covariance matrix is G = Phi*Theta*Lambda*Theta*Phi + sigma*I
## and sigma = c is the parameter to be estiamted
## For the one-step Newton-Raphson
st.findc <- function(c, xi, Lambda, B, N){
  # B can be obtained after the smoothing is done
  # N is the number of time points (images)
  n <- nrow(B)  # this is qdim
  Score <- 0
  Heissian <- 0
  for (i in 1:ncol(xi)) {
    Ai <- xi[,i] %*% t(xi[,i])
    Score <- Score + N*(1/(Lambda[i,i]+c) - 1/c)*sum(diag(Ai)) -
      (1/(Lambda[i,i] + c)^2 - 1/(c^2))*sum(diag(Ai %*% B))
    Heissian <- Heissian - N*(1/(Lambda[i,i] + c)^2 - 1/(c^2))*sum(diag(Ai)) + 
      (1/(Lambda[i,i] + c)^3 - 1/(c^3))*sum(diag(Ai %*% B))
  }
  # the Score and Heissian
  Score <- Score + N*n/c - sum(diag(B))/(c^2)
  Heissian <- Heissian - N*n/(c^2) + sum(diag(B))/(c^3)
  list(Score=Score, Heissian=Heissian)
}


## For the MLE of c
sparse.MLE.c <- function(y, A, Phi, mu0, Sigma0, cHt, c, Theta, Lambda, Pmiss) {
  y.mat <- as.matrix(y)
  n <- nrow(y.mat)    # this is num
  N <- ncol(y.mat)    # this is qdim
  Phi <- as.matrix(Phi)
  df <- nrow(Phi)     # this is pdim
  
  NA.matrix <- is.na(y.mat)
  
  ## (1) The filter and smoother
  KS <- sparse.Psmooth.G(y=y, A=A, Phi=Phi, mu0=mu0, Sigma0=Sigma0, cHt=cHt, Theta=Theta, 
                         Lambda=Lambda, sigma=c, Pmiss=Pmiss, NA.mat=NA.matrix)
  a.f <- KS$xf
  P.f <- KS$Pf
  a.s <- KS$xs
  P.s <- KS$Ps
  a.p <- KS$xp
  P.p <- KS$Pp
  Jmat <- KS$J
  J0 <- KS$J0
  Kgain <- KS$K   # only the last Kalman gain
  
  # the lag-1 covariance smoother
  Plag.s <- sparse.Plag1(y=y, A=A, Phi=Phi, Pf=P.f, J=Jmat, J0=J0, K=Kgain)
  
  # the initial state
  a00 <- as.matrix(mu0, nrow=df, ncol=1)
  P00 <- as.matrix(Sigma0, nrow=df, ncol=df)
  
  ## (2) The elements of the likelihood
  # the likelihood of the smoothed state
  S00 <- S01 <- S11 <- matrix(0, nrow=df, ncol=df)
  HHt.inv <- solve(cHt)%*%t(solve(cHt))
  
  for (i in 2:N){
    # compute the smoothed states and their variance
    at0 <- as.matrix(a.s[,i-1])
    at1 <- as.matrix(a.s[,i])
    Pt0 <- as.matrix(P.s[,,i-1])
    Pt1 <- as.matrix(P.s[,,i])
    Pt01 <-  as.matrix(Plag.s[,,i])
    
    S11 <- S11 + Pt1 + at1%*%t(at1)
    S01 <- S01 + Pt01 + at1%*%t(at0)
    S00 <- S00 + Pt0 + at0%*%t(at0)
  }
  
  # Add on the first smoothed state (T = 1)
  # for this, we need the initial states (T = 0)
  a0 <- KS$x0n
  P0 <- KS$P0n
  S11 <- S11 + as.matrix(P.s[,,1]) + as.matrix(a.s[,1])%*%t(as.matrix(a.s[,1]))
  S01 <- S01 + as.matrix(Plag.s[,,1]) + as.matrix(a.s[,1])%*%t(a0)
  S00 <- S00 + P0 + a0%*%t(a0)
  
  at.lik <- S11 - S01%*%t(Phi) - Phi%*%t(S01) + Phi%*%S00%*%t(Phi)
  a0.lik <- P0 + (a0-a00)%*%t(a0-a00)
  at.sum <- sum(diag(at.lik%*%HHt.inv)) + 2*N*sum(log(diag(cHt)))
  a0.sum <- sum(diag(a0.lik%*%solve(Sigma0))) + sum(log(diag(Sigma0)))
  
  # the likelihood of the observations
  xi <- A %*% Theta
  lambda <- diag(Lambda)
  eta <- lambda/(c*(c + lambda))
  GGt <- xi %*% Lambda %*% t(xi) + c*diag(n)
  cGt <- chol(GGt)
  GGt.inv <- 1/c*diag(n) - xi %*% diag(eta) %*%t(xi)
  Yt.lik <- matrix(0, nrow=n, ncol=n)
  
  for (i in 1:N){
    whereNA <- NA.matrix[,i]
    PNA <- sum(whereNA)/length(whereNA) 
    
    if (!any(whereNA)) {
      rt <- y.mat[,i] - A %*% as.matrix(a.s[,i])
      Yt.lik <- Yt.lik + A %*% as.matrix(P.s[,,i]) %*%t(A) + rt%*%t(rt)
    }
    else {
      if((all(whereNA))|(PNA >= Pmiss)){
        Yt.lik <- Yt.lik + 0
      }
      else {
        yes <- as.vector(!whereNA)
        no <- as.vector(whereNA) 
        nyes <- sum(yes)
        
        y.obs <- rep(0, times=n)
        y.obs[1:nyes] <- y.mat[yes,i]
        A.obs <- matrix(0, nrow=n, ncol=df)
        A.obs[1:nyes,] <- A[yes,]
        G.mis <- matrix(0, nrow=n, ncol=n)
        G.mis[(nyes+1):n,(nyes+1):n] <- GGt[no,no]
        
        rt <- y.obs - A.obs%*%as.matrix(a.s[,i])
        Yt.lik <- Yt.lik + A.obs%*%as.matrix(P.s[,,i])%*%t(A.obs) + rt%*%t(rt) + G.mis
      }
    }
  }
  Yt.sum <- sum(diag(Yt.lik%*%GGt.inv)) + 2*N*sum(log(diag(cGt)))
  
  # So the log-likelihood is
  like <- 0.5*(a0.sum + at.sum + Yt.sum)
  
  # (3) The MLEs (analytical solution)
  # Phi.mle <- S01%*%solve(S00)
  # Htemp <- 1/N * (S11 - S01%*%solve(S00)%*%t(S01))
  Htemp <- 1/N * (S11 - 2*S01 + S00)
  H.mle <- (t(Htemp) + Htemp)/2
  fc <- function(x){
    st.findc(c=x, xi=xi, Lambda=Lambda, B=Yt.lik, N=N)$Score
  }
  c.mle <- uniroot(fc, interval=c(1e-4, 1e-1), extendInt='yes', trace=1)$root
  
  list(H.mle=H.mle, c.mle=c.mle, loglike=like)
}




## -----------------------------------------------------
## PART B : the cycle for FPCA component
## -----------------------------------------------------
## (I) Initialization
fpca.ini <- function(data, basis, K, take.mean=T, pert=0.01){
  ## data are stored in a matrix
  ## base is the orthogonalized basis matrix
  N <- dim(data)[2]
  n <- dim(data)[1]
  dimb <- dim(basis)
  if (take.mean==T) data.m <- apply(data, MARGIN=1, FUN=mean, na.rm=T)
  else data.m <- rep(0, times=n)
  
  ## The basis array
  B.arr <- array(NA, dim=c(dimb[1],dimb[2],N))
  for (i in 1:N){
    index <- which(!is.na(data[,i]))
    B.arr[index,,i] <- basis[index,]
  }
  
  ## Mean function using pooled data
  data.p <- as.vector(data)
  base.p <- apply(basis, MARGIN=2, FUN=rep, times=N)
  data.mv <- rep(data.m, times=N)
  # remove 'NA's to get the complete data
  # then center the observations
  index2 <- which(!is.na(data.p))
  data.c <- data.p[index2] - data.mv[index2]
  base.c <- base.p[index2,]
  base.inv <- solve(crossprod(base.c))
  theta.ini <- base.inv %*% t(base.c) %*% data.c
  meanfun <- basis %*% theta.ini + data.m
  
  ## Residuals and sigma
  resi <- data
  for (i in 1:N){
    Z <- data[,i]
    resi[,i] <- Z - meanfun
  }
  sigma.ini <- sum(resi^2, na.rm=T)/length(resi)
  
  ## Random coefficient Theta and covariance D
  ## Follow the methods in 'fpca'
  Gamma <- matrix(0, nrow=dimb[2], ncol=N)
  for(i in 1:N){
    index <- which(!is.na(data[,i]))
    B <- B.arr[index,,i]
    Z <- resi[index,i]
    Gamma[,i] <- solve(crossprod(B) + pert*diag(dimb[2]))%*%t(B)%*%Z
    # crude estimation of random effects
  }
  Gamma.cov <- cov(t(Gamma))
  d.Gamma <- eigen(Gamma.cov, symmetric=T)
  Theta.ini <- d.Gamma$vectors[,1:K]
  D.ini <- diag(d.Gamma$values[1:K])
  
  ## return the results
  list(basis=B.arr, theta=theta.ini, Theta=Theta.ini, D=D.ini,
       sigma=sigma.ini)
}



## (II) E-step function for the new FPCA algorithm
AECM.fpca.Estep <- function(data, basis, P, meanfun,  a.s, P.s, Theta, sigma, D){
  # data here is n*N
  # D is Lambda matrix
  N <- ncol(data)
  n <- nrow(data)
  df <- ncol(basis)
  
  alpha <- matrix(0, nrow=P, ncol=N)
  alpha2 <- array(0, dim=c(P, P, N))
  alphabeta <- array(0, dim=c(P, df, N))
  
  # the E-step estimations
  for(i in 1:N){
    index <- which(!is.na(data[,i]))
    B <- basis[index,]
    Z <- data[index,i]
    
    # the following are the general results (without orthogonal assumption on the basis)
    mix.varinv <- svd.inverse(P.s[,,i] + Theta%*%D%*%t(Theta))
    mix.inv <- svd.inverse(sigma*mix.varinv + crossprod(B))
    Estep.inv <- 1/sigma*diag(length(index)) - 1/sigma*B%*%mix.inv%*%t(B)
    
    alpha[,i] <- D%*%t(Theta)%*%t(B) %*% Estep.inv %*% (Z - meanfun[index,i])
    alpha2[,,i] <- D - D%*%t(Theta)%*%t(B) %*% Estep.inv %*% B%*%Theta%*%D + 
      alpha[,i]%*%t(alpha[,i])
    alphabeta[,,i] <- alpha[,i]%*%t(a.s[,i])
    
  }
  
  list(alpha=alpha, alpha2=alpha2, alphabeta=alphabeta)
}


## (III) M-step function for the new FPCA algorithm
AECM.fpca.Mstep <- function(data, basis, P, meanfun, a.s, P.s, Theta, 
                            alpha, alpha2, alphabeta, tol=0.0001){
  N <- ncol(data)
  n <- nrow(data)
  df <- ncol(basis)
  
  ## Estimate variance matrix D
  D <- matrix(0, nrow=P, ncol=P)
  for (p in 1:P){
    D[p,p] <- sum(alpha2[p,p,])/N
  }
  
  ## Estimate Theta 
  # A loop is used here because this is essentially an ECM step
  rss.old <- 0
  rss.new <- 1
  loop <- 1
  while ((abs(rss.new - rss.old)/rss.new > tol)&(loop < 20)){
    rss.old <- rss.new
    
    # BB for t(B)*B, BZ for t(B)*Z
    # Bc for t(B)*B*alpha*beta
    # Br for t(B)*Z*alpha
    # Ba for t(B)*B*alpha2 
    BB <- array(0, dim=c(df, df, N))
    BZ <- matrix(0, nrow=df, ncol=N)
    Ba <- array(0, dim=c(df, df, P, N))
    Bc <- Br <- BZ
    
    for(p in 1:P){
      for(i in 1:N){
        index <- which(!is.na(data[,i]))
        B <- basis[index,]
        Z <- data[index,i]
        BB[,,i] <- crossprod(B)
        BZ[,i] <- t(B)%*%Z
        Br[,i] <- BZ[,i]*alpha[p,i]
        Bc[,i] <- BB[,,i]%*%alphabeta[p,,i]
        Ba[,,,i] <- outer(BB[,,i], alpha2[,p,i])
      }
      
      # To compute the estimation equations, sum over time first (i), 
      # multiply the corresponding column in Theta, then sum over PCs (p)
      Ba.sum <- apply(Ba, MARGIN=c(1,2,3), FUN=sum)
      Ba.p.sum <- Ba.sum[,,p]
      Ba.notp <- matrix(0, nrow=df, ncol=P)
      j.notp <- (1:P)[-p]
      for(j in j.notp){
        Ba.notp[,j] <- Ba.sum[,,j]%*%Theta[,j]
        # The p-th column of Ba.notp.Theta is zero (not updated)
      }
      Ba.notp.sum <- apply(Ba.notp, MARGIN=1, FUN=sum)
      Br.p.sum <- apply(Br, MARGIN=1, FUN=sum) 
      Bc.p.sum <- apply(Bc, MARGIN=1, FUN=sum)
      Theta[,p] <- solve(Ba.p.sum)%*%(Br.p.sum - Bc.p.sum - Ba.notp.sum)
    }
    
    ## (3) If I also estimate sigma in this cycle
    sig <- 1:N
    err <- 1:N
    n.obs <- 1:N
    for (i in 1:N){
      index <- which(!is.na(data[,i]))
      B <- basis[index,]
      Z <- data[index,i]
      mTime <- meanfun[index,i]  # the state space component
      n.obs[i] <- length(index)
      var.temp <- P.s[,,i] + a.s[,i]%*%t(a.s[,i]) + 2*Theta%*%alphabeta[,,i] + 
        Theta%*%alpha2[,,i]%*%t(Theta)
      sig[i] <- t(Z) %*% (Z - 2 * (mTime + B%*%Theta%*%alpha[,i])) + 
        sum(diag(B%*%var.temp%*%t(B)))
      err[i] <- crossprod(Z - mTime - B%*%Theta%*%alpha[,i])
    }
    sigma <- sum(sig)/sum(n.obs)
    
    # update the criterion and count of loops
    rss.new <- sum(err)
    loop <- loop + 1
  }
  
  list(D=D, Theta=Theta, loop=loop, rss=rss.new, sigma=sigma)
}


## (IV) The function for computing the final joint likelihood
## Will not use this...
AECM.fpca.lik <- function(data, basis, Phi, KS, Plag.s, mu0, Sigma0, 
                          cHt, Theta, D, sigma, alpha, rss){
  N <- ncol(data)
  n <- nrow(data)
  df <- ncol(basis)
  
  # computing the likelihood for alpha and part of the observation
  lik.each <- 1:N
  nobs <- 1:N
  for (i in 1:N){
    index <- which(!is.na(data[,i]))
    nobs[i] <- length(index)
    Z <- data[index,i]
    B <- basis[index,]
    
    lik.obs <- - 0.5*nobs[i]*log(sigma)
    lik.alpha <- - 0.5*sum(log(diag(D))) - 0.5*t(alpha[,i])%*%solve(D)%*%alpha[,i]
    lik.each[i] <-  lik.obs + lik.alpha
  }
  
  beta.lik <- AECM.KFa.lik(Phi=Phi, KS=KS, Plag.s=Plag.s, mu0=mu0, Sigma0=Sigma0, cHt=cHt, est=F)
  lik.beta <- - 0.5*(beta.lik$at.sum + beta.lik$a0.sum)
  lik.joint <- - 0.5/sigma*rss + sum(lik.each) + lik.beta
  # the constant, - 0.5*(sum(nobs) + df + P)*log(2*pi), has been omitted here
  # I might need it later as it involves df and P
  lik.joint
}



## (V) The final orthogonalization
AECM.fpca.orth <- function(data, basis, P, meanfun, a.s, P.s, Theta, sigma, D){
  
  # (1) Get eigenfunctions, using decomposition
  n <- dim(data)[1]  # 'grid length = n' in the mixFPCA code
  N <- dim(data)[2]
  
  Gamma <- Theta %*% D %*% t(Theta)
  decomp <- eigen(Gamma, symmetric=T)
  eigenval <- decomp$values[1:P] # / grid.length
  eigenvec <- decomp$vectors[,1:P]
  eigenfun <- basis %*% eigenvec
  Lambda <- diag(eigenval)  # * grid.length
  
  # (2) Update the alpha (scores)
  # this is slightly different as now there is the State space component
  # affecting the variance of the model, hence the conditional distribution
  alpha.est <- matrix(0, nrow=P, ncol=N)
  
  for(i in 1:N){
    index <- which(!is.na(data[,i]))
    B <- basis[index,]
    Z <- data[index,i]
    E <- eigenfun[index,]
    
    gamma <- eigenvec %*% Lambda %*% t(eigenvec)
    mix.varinv <- svd.inverse(P.s[,,i] + gamma)
    mix.inv <- svd.inverse(sigma*mix.varinv + crossprod(B))
    Estep.inv <- 1/sigma*diag(length(index)) - 1/sigma*B%*%mix.inv%*%t(B)
    alpha.est[,i] <- Lambda%*%t(E) %*% Estep.inv %*% (Z - meanfun[index,i])
    
  }
  
  # (3) Return the results
  list(newTheta=eigenvec, Lambda=Lambda, eigenfun=eigenfun, score=alpha.est)
  
}




