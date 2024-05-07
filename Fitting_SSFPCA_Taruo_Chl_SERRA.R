## -------------------------------------------
## Fitting the model using AECM algorithm
## -------------------------------------------

## NOT RUN! FOR INFORMATION ONLY!
# load('Taro_Lake_Chl.RData')
# 
# # Mask file
# maskfile <- '00000383_reducedmask.tif'
# mask  <- raster(maskfile)
# 
# # Log transfomation
# Taro.arr[Taro.arr < 0.005] <- 0.005
# logTaro.arr <- log(Taro.arr)
# 
# # Trimming
# lon.trim <- Lon[8:110]
# lat.trim <- Lat[10:54]
# Taro.trim <- array(NA, dim=c(length(lon.trim), length(lat.trim), 101))
# for (i in 1:101) {
#   Taro.trim[,,i] <- t(Taro.arr[10:54, 8:110, i])
# }
# 
# Taro.mean <- apply(Taro.trim[,,-1], MARGIN=3, FUN=mean, na.rm=T)
# Taro.data <- apply(Taro.trim[,,-1], MARGIN=3, FUN=as.vector)
# Taro.data <- Taro.data - matrix(rep(Taro.mean, each=nrow(Taro.data)), nrow=nrow(Taro.data), ncol=100)
# 
# mask.mat <- matrix(flip(mask, direction='y'), nrow=length(Lon), ncol=length(Lat))
# mask.mat <- mask.mat[8:110, 10:54]
# sum(mask.mat == 1, na.rm=T) / length(mask.mat)
# # 94.43% data for the lake
# # 84.97% data for the grid



## IMPLEMENTING THE SSFPCA

source('Functions_AECM_G.R')

load('Taro_Chl_trimmed.RData')



## Step 0: Preparations
# Center the data
Taro.mean <- apply(Taro.trim[,,-1], MARGIN=3, FUN=mean, na.rm=T)
Taro.data <- apply(Taro.trim[,,-1], MARGIN=3, FUN=as.vector)
Taro.data <- Taro.data - matrix(rep(Taro.mean, each=nrow(Taro.data)), nrow=nrow(Taro.data), ncol=100)

# create the basis matrix
k1 <- 4   # number of interior knots along longitude
k2 <- 2   # number of interior knots along latitude

knot1 <- quantile(lon.trim, prob=(1:k1)/(k1+1))
knot2 <- quantile(lat.trim, prob=(1:k2)/(k2+1))
bs.lon <- bs(lon.trim, knots=knot1, degree=3, intercept=T)  # longitude
bs.lat <- bs(lat.trim, knots=knot2, degree=3, intercept=T)  # latitude
bs.grid <- kronecker(X=as.matrix(bs.lon), Y=as.matrix(bs.lat)) 

temp <- t(bs.lon)%*%bs.lon
R.lon <- t(chol(temp)) 
B.lon <- t(solve(R.lon, t(bs.lon)))
temp <- t(bs.lat)%*%bs.lat
R.lat <- t(chol(temp)) 
B.lat <- t(solve(R.lat, t(bs.lat)))
Base <- kronecker(X=B.lat, Y=B.lon)    # The tensor spline basis
df <- ncol(Base)   # the basis dimension

# some other quantities
P <- 6   # the number of functional PCs
Pmiss <- 0.9   # the Kalman filtering threshold

n <- nrow(Taro.data)
N <- ncol(Taro.data)
n.obs <- sum(!is.na(Taro.data))
NA.matrix <- is.na(Taro.data)
index.lk <- mask.mat == 1   # indeces of lake pixels (mask.mat gives the lake mask)
index.lk[is.na(index.lk)] <- FALSE



## Step 1: Initialization
# the state space componnt
beta0 <- mu0 <- rep(0, times=df)
P0 <- Sigma0 <-  100*diag(df)
Mt <- 1*diag(df)
sigma.h <- sum(apply(Taro.data, MARGIN=1, FUN=var, na.rm=T), na.rm=T)/n
HHt <- sigma.h*diag(df)
cHt <- chol(HHt)

# the FPCA component
par.zero <- fpca.ini(data=Taro.data, basis=Base, P=P, take.mean=F, pert=0.001)
Theta <- par.zero$Theta
xi <- Base %*% Theta
Lambda <- par.zero$D
sigma <- par.zero$sigma



## Step 2: The AECM iterations
loglike0 <- 0
loglike1 <- 1
it <- 0

while ((abs((loglike1 - loglike0) / loglike0) > 0.0005)&&(loglike1 - loglike0 >= 0)&&(it <= 10)) {

  begin <- proc.time()

  it <- it + 1
  loglike0 <- loglike1

  ## CYCLE 1 : the state space model
  ## Using a set of notation consistent with the KFilter functions

  ## Apply the filter and smoother (E-step)
  KS <- sparse.Psmooth.G(y=Taro.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
                         Theta=Theta, Lambda=Lambda, sigma=sigma,
                         Pmiss=Pmiss, NA.mat=NA.matrix)
  a.f <- KS$xf
  P.f <- KS$Pf
  a.s <- KS$xs
  P.s <- KS$Ps
  a0 <- KS$x0n
  P0 <- KS$P0n
  Jmat <- KS$J
  J0 <- KS$J0
  Kgain <- KS$K
  Plag.s <- sparse.Plag1(y=Taro.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)

  ## Compute the likelihood and get the MLEs (M-step)
  state.lik <- AECM.KFa.lik(Phi=Mt, KS=KS, Plag.s=Plag.s, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  a0.sum <- state.lik$a0.sum
  at.sum <- state.lik$at.sum
  S00 <- state.lik$S00
  S01 <- state.lik$S01
  S11 <- state.lik$S11

  ## The MLEs (including the scoring method)
  Mt.temp <- diag(S01) / diag(S00)
  Mt.mle <- diag(Mt.temp)

  Htemp <- 1/N * (S11 - 2*S01%*%t(Mt.mle) + Mt.mle%*%S00%*%t(Mt.mle))
  H.mle <- (t(Htemp) + Htemp)/2

  ## Update the smoothed state with the new MLEs
  Mt <- Mt.mle
  HHt <- H.mle
  cHt <- chol(HHt)

  KS.cycle1 <- sparse.Psmooth.G(y=Taro.data, A=Base, Phi=Mt, mu0=mu0, Sigma0=Sigma0, cHt=cHt,
                                Theta=Theta, Lambda=Lambda, sigma=sigma,
                                Pmiss=Pmiss, NA.mat=NA.matrix)
  a.f <- KS.cycle1$xf
  P.f <- KS.cycle1$Pf
  a.s <- KS.cycle1$xs
  P.s <- KS.cycle1$Ps
  a0 <- KS.cycle1$x0n
  P0 <- KS.cycle1$P0n
  a.p <- KS.cycle1$xp
  P.p <- KS.cycle1$Pp
  Jmat <- KS.cycle1$J
  J0 <- KS.cycle1$J0
  Kgain <- KS.cycle1$K
  Plag.cycle1 <- sparse.Plag1(y=Taro.data, A=Base, Phi=Mt, Pf=P.f, J=Jmat, J0=J0, K=Kgain)


  ## CYCLE2 : the mixed model FPCA
  meanfun <- Base%*%a.s

  ## The E-step
  par.E <- AECM.fpca.Estep(data=Taro.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
                           Theta=Theta, sigma=sigma, D=Lambda)
  alpha <- par.E$alpha
  alpha2 <- par.E$alpha2
  alphabeta <- par.E$alphabeta

  ## M-step
  par.M <- AECM.fpca.Mstep(data=Taro.data, basis=Base, P=P, meanfun=meanfun, a.s=a.s, P.s=P.s,
                           Theta=Theta, alpha=alpha, alpha2=alpha2, alphabeta=alphabeta)
  Lambda.mle <- par.M$D
  Theta.mle <- par.M$Theta
  sigma.mle <- par.M$sigma
  rss <- par.M$rss

  # update the parameters and related quantities
  Theta <- Theta.mle
  Lambda <- Lambda.mle
  sigma <- sigma.mle


  ## Update the joint log-likelihood
  lik.each <- 1:N
  nobs <- 1:N
  for (i in 1:N){
    index <- which(!is.na(Taro.data[,i]))
    nobs[i] <- length(index)
    Z <- Taro.data[index,i]
    B <- Base[index,]

    lik.obs <- - 0.5*nobs[i]*log(sigma)
    lik.alpha <- - 0.5*sum(log(diag(Lambda))) - 0.5*t(alpha[,i])%*%solve(Lambda)%*%alpha[,i]
    lik.each[i] <-  lik.obs + lik.alpha
  }

  beta.lik <- AECM.KF.lik(Phi=Mt, KS=KS.cycle1, mu0=mu0, Sigma0=Sigma0, HHt=HHt)
  lik.beta <- - 0.5*(beta.lik$at.sum + beta.lik$a0.sum)
  lik.joint <- - 0.5/sigma*rss + sum(lik.each) + lik.beta   # - 0.5*(sum(nobs) + df + P)*log(2*pi)
  loglike1 <- lik.joint

  end <- proc.time() - begin
  print(paste('Iteration', it, '; Timer', end[1], '; loglike = ', loglike1))

}



## Step 3: Keep the results
AECM.sigma <- sigma
AECM.Theta <- Theta

AECM.Mt <- Mt
AECM.HHt <- HHt
AECM.beta <- a.s
AECM.vbeta <- P.s
AECM.D <- Lambda

AECM.tmean <- Base %*% AECM.beta
final <- AECM.fpca.orth(data=Taro.data, basis=Base, P=P, meanfun=AECM.tmean, a.s=AECM.beta,
                        P.s=AECM.vbeta, Theta=AECM.Theta, sigma=AECM.sigma, D=AECM.D)
AECM.newTheta <- final$newTheta
AECM.Lambda <- final$Lambda
AECM.alpha <- final$score

