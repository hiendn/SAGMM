#' @import stats MixSim mvtnorm mclust lowmemtkmeans



#' Data for simulations.
#' 
#' @description  Wrapper for MixSim
#' @param Groups Number of mixture components. Default 5.
#' @param Dimensions number of Dimensions. Default 5.
#' @param Number number of samples. Default 10^4.
#' @return List of results: X, Y, simobject.
#' @examples
#' sims<-generateSimData(Groups=10, Dimensions=10, Number=10^4)
#' sims<-generateSimData()
#'
#'@export
generateSimData<-function(Groups=5, Dimensions=5, Number=10^4,BAR=0.01){
    MS <- MixSim::MixSim(BarOmega=BAR,K=Groups,p=Dimensions,PiLow=(0.1/Groups)) # Simulation code
    Z <- MixSim::simdataset(n=Number,Pi=MS$Pi,Mu=MS$Mu,S=MS$S) # More simulation code, look at package vignette
    X <- Z[[1]] # Extract features
    Y <- Z[[2]] # Extract Labels
    SAMPLE <- sample(1:Number) # Randomize
    X <- X[SAMPLE,] # Randomize
    Y <- Y[SAMPLE] # Randomize
    return(list(X=X, Y=Y, MS=MS))
}

main_loop_R<-function(Number,Groups, PISTAR_O, MU_O, LAMBDA_O, GAMMA, X, Dimensions, SIGMA){
    
    LogLike <- 0
    MU <- MU_O
    LAMBDA <- LAMBDA_O
    PISTAR <- PISTAR_O
    
    ### Main loop
    for (ii in 1:Number) {
        ### PI
        PI <- exp(PISTAR_O)/sum(exp(PISTAR_O)) # Compute Pi from Pistar olds
        
        ### Tau # this is the component/overall_distribution that goes at the start of each gradient component
        Tau <- c()
        for (gg in 1:Groups) {
            Tau[gg] <- PI[gg]*mvtnorm::dmvnorm(X[ii,],MU_O[gg,],diag(Dimensions)*LAMBDA_O[gg]^2/2)
        }
        Tau <- Tau/sum(Tau) 
        
        ### Mean
        for (gg in 1:Groups) {
            MU[gg,] <- MU_O[gg,] + GAMMA[ii]*Tau[gg]*2/LAMBDA_O[gg]^2*(X[ii,]-MU_O[gg,])
        }
        
        LAMBDASTAR <- log(LAMBDA)
        # Sigma
        for (gg in 1:Groups) {
            LAMBDASTAR[gg] <- log(LAMBDA_O[gg]) - GAMMA[ii]*Tau[gg]*(Dimensions-2/LAMBDA_O[gg]^2*t((X[ii,]-MU_O[gg,]))%*%(X[ii,]-MU_O[gg,]))
        }
        LAMBDA <- exp(LAMBDASTAR)
        for (gg in 1:Groups) {
            SIGMA[[gg]] <- diag(Dimensions)*LAMBDA[gg]^2/2
        }
        
        # ### Pi
        for (gg in 1:(Groups-1)) {
            PISTAR[gg] <- PISTAR_O[gg] + GAMMA[ii]*Tau[gg]*(2 - exp(PISTAR_O[gg])/sum(exp(PISTAR_O)))
        }
        
        Comps <- c()
        for (gg in 1:Groups) {
            Comps[gg] <- log(PISTAR[gg]/sum(PISTAR)) + dmvnorm(X[ii,],MU[gg,],SIGMA[[gg]],log=T)
        }
        LogLike <- LogLike + max(Comps) + log(sum(exp(Comps-max(Comps))))
        
        
        ### Old stuff
        LAMBDA_O <- LAMBDA # Put old stuff back
        MU_O <- MU
        PISTAR_O <- PISTAR
        
        #message(c(ii))
    }
    
    
    return(list(PI=exp(PISTAR)/sum(exp(PISTAR)), MU=MU, LAMBDA=LAMBDA))
}


#' Clustering via Stochastic Approximation and Gaussian Mixture Models
#' 
#' @description  Fit a GMM with SA. 
#' @param X numeric maxtrix of the data.
#' @param Y True classes.
#' @param MS MixSim Simulation object.
#' @param BURNIN Ratio of observations to use as a burn in before algorithm begins.
#' @param Groups Number of mixture components.
#' @param kstart number of kmeans starts to initialise.
#' @param mode for testing "C" uses C++, "R" uses R code.
#' @param plot If TRUE generates a plot of the clustering.
#' @return List of results:
#' @examples
#' sims<-generateSimData(Groups=10, Dimensions=10, Number=10^4)
#' res1<-SAGMMFit(sims$X, sims$Y, sims$MS, mode="C")
#'
#'@export
SAGMMFit <- function (X, Y, BURNIN = 10, Groups = 5, kstart = 10, mode = "C") 
{
  Number <- nrow(X)
  Dimensions <- ncol(X)
  KM <- suppressWarnings(stats::kmeans(X[1:round(Number/BURNIN), 
                                         ], Groups, nstart = kstart))
  MU <- KM[[2]]
  LAMBDA <- rep(max(sqrt(KM$withinss/KM$size)), Groups)
  SIGMA <- list()
  SIGMA_C <- array(0, c(Dimensions, Dimensions, Groups))
  for (gg in 1:Groups) {
    SIGMA[[gg]] <- diag(Dimensions) * LAMBDA[gg]^2/2
    SIGMA_C[, , gg] <- SIGMA[[gg]]
  }
  PI <- rep(1/Groups, Groups)
  PISTAR <- rep(0, Groups)
  for (it in 1:100) {
    for (gg in 1:(Groups - 1)) {
      PISTAR[gg] <- log((1/PI[gg] - 1)^-1 * sum(exp(PISTAR[-gg])))
    }
  }
  GAMMA <- log(Number)/Number*c(rep(1,
                                    round(Number/BURNIN)),
                                1/(1+(1:(Number-round(Number/BURNIN))/
                                        (Number-round(Number/BURNIN)))))
  LAMBDA_O <- LAMBDA
  MU_O <- MU
  PISTAR_O <- PISTAR
  if (mode == "C") {
    results <- main_loop_C(Number, Groups, PISTAR_O, MU_O, 
                           LAMBDA_O, GAMMA, X, Dimensions, SIGMA_C)
    TauMAT <- results$Tau
    PI <- results$PI
    MU <- results$MU
    SIGMA <- list()
    for (gg in 1:Groups) {
      SIGMA[[gg]] <- diag(Dimensions) * results$LAMBDA[gg]^2/2
    }
  }
  if (mode == "R") {
    results <- main_loop_R(Number, as.integer(Groups), PISTAR_O, 
                           MU_O, LAMBDA_O, GAMMA, X, Dimensions, SIGMA)
    PI <- results$PI
    MU <- results$MU
    SIGMA <- list()
    for (gg in 1:Groups) {
      SIGMA[[gg]] <- diag(Dimensions) * results$LAMBDA[gg]^2/2
    }
    TauMAT <- matrix(NA, Number, Groups)
    for (ii in 1:Number) {
      Tau <- c()
      for (gg in 1:Groups) {
        Tau[gg] <- log(PI[gg]) + dmvnorm(X[ii, ], MU[gg, 
                                                     ], SIGMA[[gg]]/2, log = T)
      }
      TauMAT[ii, ] <- Tau
    }
  }
  retList <- list(fit = results)
  return(retList)
}















