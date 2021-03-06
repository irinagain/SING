# Simulation Setting 3 (True components are sparse)
# Identification of the number of components + methods application 


# Load the simulated data as well as true components
load("SimulatedLargeScaleindiv10_Sparse.Rdata")

# Source functions for estimation
source("jngcaFunctions.R")
source("mCCAjointICA.R")

# Center X and Y
n = nrow(dX)
pX = ncol(dX)
pY = ncol(dY)
dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
dYcentered <- dY - matrix(rowMeans(dY), n, pY, byrow = F)


## Apply Joint ICA
#############################################################
out_jointICA <- jointICA(dXcentered, dYcentered, r0 = 2)

# Prepare the output list for joint ICA
save(out_jointICA, file = "jointICA_LargeScaleindiv10_Sparse.Rda")

# Joint rank estimation based on oversaturated models
#####################################################

estX_JB = mlcaFP(xData = t(dX), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(dX))

estY_JB = mlcaFP(xData = t(dY), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

# Get joint components out
My_JB = est.M.ols(sData = estY_JB$S, xData = t(dY))

alpha = 0.05
nperms = 1000
matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
joint_rank = min(which(permuteJoint$pvalues > alpha)) - 1
pval_joint = permuteJoint$pvalues
joint_rank # selects rank 2


## Apply separate JB 
#############################################################
# JB on X
estX_JB = mlcaFP(xData = t(dX), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uxfull <- estX_JB$Ws

# Match mixing matrices to get correct ordering, then can get starting points
Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(dX))

errorMx = frobICA(M1 = Mx_JB, M2 = t(mj), standardize = T)*sqrt(ncol(Mx_JB))
errorMx # 0.225

errorSx = frobICA(S1 = t(Sx[1:2, ]), S2 = estX_JB$S, standardize  = T)
errorSx # 0.029

# JB on Y
estY_JB = mlcaFP(xData = t(dY), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uyfull <- estY_JB$Ws

# Get joint components out
My_JB = est.M.ols(sData = estY_JB$S, xData = t(dY))

errorMy = frobICA(M1 = My_JB, M2 = t(mj), standardize = T)*sqrt(ncol(My_JB))
errorMy # 0.168

errorSy = frobICA(S1 = t(Sy[1:2, ]), S2 = estY_JB$S, standardize  = T)
errorSy # 0.015


# Prepare the output list for separate JB
output_sepJB <- list(errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, estX_JB = estX_JB, estY_JB = estY_JB, mj = mj)
save(output_sepJB, file = "sepJB_LargeScale_Sparse.Rda")


# Joint rank estimation based on separate models fitted with true number of components
#####################################################
# Rank estimation via permutation
alpha = 0.05
nperms = 1000
matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
joint_rank = min(which(permuteJoint$pvalues > alpha)) - 1
pval_joint = permuteJoint$pvalues
joint_rank # selects rank 2


# Whiten dX and dY
###################################

# Scale rowwise
est.sigmaXA = tcrossprod(dXcentered)/(pX-1)  ## since already centered, can just take tcrossprod
whitenerXA = est.sigmaXA%^%(-0.5)
invLx = est.sigmaXA%^%(0.5)
xDataA = whitenerXA %*% dXcentered

# For Y
# Scale rowwise
est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
whitenerYA = est.sigmaYA%^%(-0.5)
yDataA = whitenerYA %*% dYcentered
invLy = est.sigmaYA%^%(0.5)

# Starting points for the algorithm are Us
###########################################

matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))

cor(matchMxMy$Mx[, 1:2], matchMxMy$My[, 1:2]) # 0.97 and 0.95

# Calculate JB values
JBall = calculateJB(matchMxMy$Ux[1:2, ], X = xDataA) + calculateJB(matchMxMy$Uy[1:2, ], X = yDataA)


# SING on shared + individual
###################################
# Libraries needed to run C functions
library(Rcpp)
library(RcppArmadillo)

sourceCpp("CfunctionsCurvilinear.cpp")

# JB and tolerance parameters
alpha = 0.8
tol = 1e-10

# small rho
rho = JBall/10
out_indiv_small <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_small, file = "out_indiv_small_Sparse.Rda")

# medium rho
rho = JBall
out_indiv_medium <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = out_indiv_small$Ux, Uy = out_indiv_small$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_medium, file = "out_indiv_medium_Sparse.Rda")

# large rho
rho = JBall * 10
out_indiv_large <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = out_indiv_medium$Ux, Uy = out_indiv_medium$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_large, file = "out_indiv_large_Sparse.Rda")

###############################################################################
# Apply mCCA + Joint ICA, and calculate the errors
###############################################################################
out_mcca <- mCCA_jointICA(dXcentered, dYcentered, Mx = 12, My = 12, M = 2)

save(out_mcca, file = "mCCA_LargeScale_Sparse.Rda")
