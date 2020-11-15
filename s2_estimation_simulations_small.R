# Simulation setting 1
# Joint ICA, separate LNGCA and SING with different rho values

# Source all the functions
source("jngcaFunctions.R")
source("mCCAjointICA.R")

# Use R libraries for parallel environment
library(doParallel)
library(doRNG)


# Source for C++ functions - makes them then usable within R
####################################################################################
# Libraries needed to run C functions
library(Rcpp)
library(RcppArmadillo)

sourceCpp("CfunctionsCurvilinear.cpp")


# Number of replications for each variable type
rep = 100
# Levels for signal to noise ratio
snr_levels = c(0.2, 5)
n_snr = length(snr_levels)
snr1levels = rep(snr_levels, n_snr)
snr2levels = rep(snr_levels, each = n_snr)

# Level for variance of inactive points in components
var_level = 0.005

# Snr for each replication
snr1 = rep(snr1levels, rep)
snr2 = rep(snr2levels, rep)

nrep = length(snr1)
nworkers <- detectCores() 
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

set.seed(234623)

#results <- foreach(i=1:nrep, .errorhandling='pass', .noexport = "updateUboth_c", .packages = "Rcpp") %dorng% {

results <- foreach(i=1:nrep, .errorhandling='pass') %dorng% {
  #sourceCpp("CfunctionsCurvilinear.cpp")
  
  # Generate the data with the specified signal and noise ratio
  data <- generateData_v2(nsubject = 48, snr = c(snr1[i], snr2[i]), vars = c(var_level, var_level))
  trueJx <- data$mj %*% data$sjX
  trueJy <- t(t(data$mj)*c(-5,2)) %*% data$sjY
  trueJxF2 <- sum(trueJx^2) # Frobenius norm squared of joint part of X
  trueJyF2 <- sum(trueJy^2) # Frobenius norm squared of joint part of Y
  
  # Prepare the output list
  output <- list(snr1 = snr1[i], snr2 = snr2[i], JBsep = NULL, JBjoint = NULL, jointICA = NULL, R2x = data$R2x, R2y = data$R2y)
  
  ##################################################################################################
  # Create whitened data
  n = nrow(data$dX) # sample size
  pX = ncol(data$dX)
  pY = ncol(data$dY)
  
  # For X
  # Center
  dXsA <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
  # Scale rowwise
  est.sigmaXA = tcrossprod(dXsA)/(pX-1)  ## since already centered, can just take tcrossprod
  whitenerXA = est.sigmaXA%^%(-0.5)
  invLx = est.sigmaXA%^%(0.5)
  xDataA = whitenerXA %*% dXsA
  
  # For Y
  # Center
  dYsA <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)
  # Scale rowwise
  est.sigmaYA = tcrossprod(dYsA)/(pY-1)  ## since already centered, can just take tcrossprod
  whitenerYA = est.sigmaYA%^%(-0.5)
  yDataA = whitenerYA %*% dYsA
  invLy = est.sigmaYA%^%(0.5)
  
  ##################################################################################################
  # Apply separate JB to each (LNGCA model), calculate error
  ##################################################################################################
  # JB on X
  estX_JB = mlcaFP(xData = t(data$dX), n.comp = 3, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

  # Mixing matrix for X
  Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX))
  Uxfull = tcrossprod(whitenerXA, Mx_JB)

  # JB on Y
  estY_JB = mlcaFP(xData = t(data$dY), n.comp = 4, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

  # Mixing matrix for Y
  My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))
  Uyfull = tcrossprod(whitenerYA, My_JB)
  
  # Get 2 joint components out
  # Match the mixing matrices, the first two should be joint
  matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))
  
  # Reorder components
  Sx = estX_JB$S[, matchMxMy$mapX]
  Sy = estY_JB$S[, matchMxMy$mapY]

  # Component errors
  errorSx = frobICA(S2 = t(data$sjX), S1 = Sx[ , 1:2], standardize  = T)
  errorSy = frobICA(S2 = t(data$sjY), S1 = Sy[ , 1:2], standardize  = T)
  
  # Mixing matrices errors
  errorMx = frobICA(M1 = t(matchMxMy$Mx[,1:2]), M2 = t(data$mj), standardize  = T)*sqrt(nrow(data$mj))
  errorMy = frobICA(M1 = t(matchMxMy$My[,1:2]), M2 = t(data$mj), standardize  = T)*sqrt(nrow(data$mj))
  
  Mdist = frobICA(M1 = t(matchMxMy$Mx[,1:2]), M2 = t(matchMxMy$My[,1:2]), standardize  = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((tcrossprod(matchMxMy$Mx[,1:2], Sx[ , 1:2]) - trueJx)^2)/trueJxF2
  errorJy = sum((tcrossprod(matchMxMy$My[,1:2], Sy[ , 1:2]) - trueJy)^2)/trueJyF2
  
  output$JBsep <- list(errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, Mdist = Mdist, errorJx = errorJx, errorJy = errorJy)

  # Apply separate JB to each (LNGCA model), but then average and recompute
  ##################################################################################################
  
  # Error from averaged mixing matrix
  newM = aveM(matchMxMy$Mx[,1:2], matchMxMy$My[,1:2])
  errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize  = T)*sqrt(nrow(data$mj))
  
  # Back projection on the mixing matrix
  outx <- est.S.backproject(t(Sx[ , 1:2]), matchMxMy$Mx[,1:2], newM)
  outy <- est.S.backproject(t(Sy[ , 1:2]), matchMxMy$My[,1:2], newM)
  
  errorSxAve = frobICA(S2 = t(data$sjX), S1 = t(outx$S), standardize  = T)
  errorSyAve = frobICA(S2 = t(data$sjY), S1 = t(outy$S), standardize  = T)
  
  # Joint signal Frobenius norm reconstruction error with average
  errorJxave = sum((newM %*% diag(outx$D) %*% outx$S/sqrt(pX-1) - trueJx)^2)/trueJxF2
  errorJyave = sum((newM %*% diag(outy$D) %*% outy$S/sqrt(pY-1) - trueJy)^2)/trueJyF2
  
  output$JBsepAve <- list(errorSx = errorSxAve, errorSy = errorSyAve, errorMx = errorM, errorMy = errorM, errorJx = errorJxave, errorJy = errorJyave)
  

  ##################################################################################################
  # Apply joint approach with 3 different values of rho, calculate error
  ##################################################################################################
  # Calculate JB values
  JBall = calculateJB(matchMxMy$Ux[1:2, ], X = xDataA) + calculateJB(matchMxMy$Uy[1:2, ], X = yDataA)
  
  # Decide on rho values based on JB values: small, medium and large

  # Small rho (JB/10)
  ##################################################################################################
  rho = JBall/10
  
  out1 = updateUboth(Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-10, r0 = 2)
  #out1 = updateUboth_c(Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-10, r0 = 2)
  
  # Calculate mixing matrices and components
  Mxjoint = tcrossprod(invLx, out1$Ux[1:2, ])
  Sxjoint = out1$Ux[1:2, ] %*% xDataA
  Myjoint = tcrossprod(invLy, out1$Uy[1:2, ])
  Syjoint = out1$Uy[1:2, ] %*% yDataA
  
  
  # What are the errors?
  # Components
  errorSx = frobICA(S2 = t(data$sjX), S1 = t(Sxjoint), standardize  = T)
  errorSy = frobICA(S2 = t(data$sjY), S1 = t(Syjoint), standardize  = T)
  # Mixing matrices
  errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  # Error from averaged mixing matrix
  newM = aveM(Mxjoint, Myjoint)
  errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  
  Mdist = frobICA(M1 = t(Mxjoint), M2 = t(Myjoint), standardize  = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((Mxjoint %*% Sxjoint - trueJx)^2)/trueJxF2
  errorJy = sum((Myjoint %*% Syjoint - trueJy)^2)/trueJyF2
  
  

  small = list(rho = rho, errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, errorM = errorM, Mdist = Mdist, errorJx = errorJx, errorJy = errorJy)

  # Medium rho (JB)
  ##################################################################################################
  rho = JBall

  out2 = updateUboth(Ux = out1$Ux, Uy = out1$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-10, r0 = 2)
  #out2 = updateUboth_c(Ux = out1$Ux, Uy = out1$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-10, r0 = 2)
  
  # Calculate mixing matrices and components
  Mxjoint = tcrossprod(invLx, out2$Ux[1:2, ])
  Sxjoint = out2$Ux[1:2, ] %*% xDataA
  Myjoint = tcrossprod(invLy, out2$Uy[1:2, ])
  Syjoint = out2$Uy[1:2, ] %*% yDataA
  
  
  # What are the errors?
  # Components
  errorSx = frobICA(S2 = t(data$sjX), S1 = t(Sxjoint), standardize  = T)
  errorSy = frobICA(S2 = t(data$sjY), S1 = t(Syjoint), standardize  = T)
  # Mixing matrices
  errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  
  # Error from averaged mixing matrix
  newM = aveM(Mxjoint, Myjoint)
  errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))

  Mdist = frobICA(M1 = t(Mxjoint), M2 = t(Myjoint), standardize  = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((Mxjoint %*% Sxjoint - trueJx)^2)/trueJxF2
  errorJy = sum((Myjoint %*% Syjoint - trueJy)^2)/trueJyF2

  medium = list(rho = rho, errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, errorM = errorM, Mdist = Mdist, errorJx = errorJx, errorJy = errorJy)

  # Large rho (JB X 2)
  ##################################################################################################
  rho = JBall*20 

  out3 = updateUboth(Ux = out2$Ux, Uy = out2$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0 = 2)
  #out3 = updateUboth_c(Ux = out2$Ux, Uy = out2$Uy, xData = xDataA, yData = yDataA, invLx = invLx, invLy = invLy, rho = rho, tau = 0.01, alpha = 0.8, maxiter = 1000, tol = 1e-6, r0 = 2)
  
  # Calculate mixing matrices and components
  Mxjoint = tcrossprod(invLx, out3$Ux[1:2, ])
  Sxjoint = out3$Ux[1:2, ] %*% xDataA
  Myjoint = tcrossprod(invLy, out3$Uy[1:2, ])
  Syjoint = out3$Uy[1:2, ] %*% yDataA
  
  
  # What are the errors?
  # Components
  errorSx = frobICA(S2 = t(data$sjX), S1 = t(Sxjoint), standardize  = T)
  errorSy = frobICA(S2 = t(data$sjY), S1 = t(Syjoint), standardize  = T)
  # Mixing matrices
  errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  
  # Error from averaged mixing matrix
  newM = aveM(Mxjoint, Myjoint)
  errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))

  Mdist = frobICA(M1 = t(Mxjoint), M2 = t(Myjoint), standardize  = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((Mxjoint %*% Sxjoint - trueJx)^2)/trueJxF2
  errorJy = sum((Myjoint %*% Syjoint - trueJy)^2)/trueJyF2

  large = list(rho = rho, errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, errorM = errorM, Mdist = Mdist, errorJx = errorJx, errorJy = errorJy)

  output$JBjoint <- list(small = small, medium = medium, large = large)
  
  ##################################################################################################
  # Apply Joint ICA, and calculate the errors
  ##################################################################################################
  out_jointICA <- jointICA(dXsA, dYsA, r0 = 2)
  
  # What are the errors?
  errorSx = frobICA(S1 = t(data$sjX), S2 = out_jointICA$S[1:pX, ], standardize  = T)
  errorSy = frobICA(S1 = t(data$sjY), S2 = out_jointICA$S[(pX+1):(pX+pY), ], standardize  = T)
  errorM = frobICA(M1 = out_jointICA$Mjoint, M2 = t(data$mj), standardize = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((tcrossprod(t(out_jointICA$Mjoint), out_jointICA$S[1:pX, ]) * sqrt(mean(dXsA^2)) - trueJx)^2)/trueJxF2
  errorJy = sum((tcrossprod(t(out_jointICA$Mjoint), out_jointICA$S[(pX+1):(pX+pY), ])* sqrt(mean(dYsA^2)) - trueJy)^2)/trueJyF2
  
  output$jointICA <- list(errorSx = errorSx, errorSy = errorSy, errorM = errorM, errorJx = errorJx, errorJy = errorJy)
  
  ###############################################################################
  # Apply mCCA + Joint ICA, and calculate the errors
  ###############################################################################
  out_mcca <- mCCA_jointICA(dXsA, dYsA, Mx = 3, My = 4, M = 2)
  
  # What are the errors?
  errorSx = frobICA(S1 = t(data$sjX), S2 = out_mcca$S[1:pX, ], standardize  = T)
  errorSy = frobICA(S1 = t(data$sjY), S2 = out_mcca$S[(pX+1):(pX+pY), ], standardize  = T)
  errorMx = frobICA(M1 = out_mcca$Mx, M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  errorMy = frobICA(M1 = out_mcca$My, M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
  
  Mdist = frobICA(M1 = out_mcca$Mx, M2 = out_mcca$My, standardize  = T) * sqrt(nrow(data$mj))
  
  # Joint signal Frobenius norm reconstruction error
  errorJx = sum((tcrossprod(t(out_mcca$Mx), out_mcca$S[1:pX, ]) - trueJx)^2)/trueJxF2
  errorJy = sum((tcrossprod(t(out_mcca$My), out_mcca$S[(pX+1):(pX+pY), ]) - trueJy)^2)/trueJyF2
  
  output$mCCA <- list(errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, Mdist = Mdist, errorJx = errorJx, errorJy = errorJy)
  
  output
}
# stop cluster
stopCluster(cl)

save(results, file = "SimResults_estimation_small.Rdata")
