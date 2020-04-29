# Simulation setting 1
# Rank identification


# Source all the functions
source("jngcaFunctions.R")

# Use R libraries for parallel environment
library(doParallel)
library(doRNG)


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

######################################################################################
# Both individual and joint rank estimation
######################################################################################

# Threshold for p-value significance
alphavec = c(0.01, 0.05, 0.1)
# Number of permutations
nperms = 1000

nworkers <- detectCores() 
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

set.seed(234623)

results <- foreach(i=1:nrep, .errorhandling='pass') %dorng% {
  # Generate the data with the specified signal and noise ratio
  data <- generateData_v2(nsubject = 48, snr = c(snr1[i], snr2[i]), vars = c(var_level, var_level))
  
  # Prepare the output list
  output <- list(snr1 = snr1[i], snr2 = snr2[i], alpha = list())
  
  ##################################################################################################   
  # Individual rank estimation
  ##################################################################################################   
  permuteX = permmatRank_sequential_JB(t(data$dX), maxcompdata = 6, ncomplist = c(1:6), nperms = nperms, ninitialization_data = 5, ninitialization_perm=5)
    
  permuteY = permmatRank_sequential_JB(t(data$dY), maxcompdata = 6, ncomplist = c(1:6), nperms = nperms, ninitialization_data = 5, ninitialization_perm=5)
  
  for (l in 1:length(alphavec)){
    # Get alpha value
    alpha = alphavec[l]
    # Figure out the rank of X
    compX = permuteX$pvalues > alpha
    if(sum(compX) > 0){
      rankX = min(which(compX)) - 1
    }else{
      rankX = length(compX)
    }
    if (rankX < 0){rankX  = 0}
    # Figure out the rank of Y
    compY = permuteY$pvalues > alpha
    if (sum(compY) > 0){
      rankY = min(which(compY)) - 1
    }else{
      rankY = length(compY)
    }
    if (rankY < 0){rankY  = 0}
    # Figure out the joint rank
    if ((rankX > 0)&(rankY > 0)){
      # JB on X
      estX_JB = mlcaFP(xData = t(data$dX), n.comp = rankX, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
      # Mixing matrix for X
      Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX))
      # JB on Y
      estY_JB = mlcaFP(xData = t(data$dY), n.comp = rankY, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
      # Mixing matrix for Y
      My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))
      # Make sure both Mx and My are in the matrix form
      if (rankX == 1){
        Mx_JB = matrix(Mx_JB, 1, length(Mx_JB))
      }
      if (rankY == 1){
        My_JB = matrix(My_JB, 1, length(My_JB))
      }
      # Match the mixing matrices, the first two should be joint
      matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
      permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
      compJ = permuteJoint$pvalues > alpha
      if (sum(compJ) > 0){
        joint_rank = min(which(compJ)) - 1
      }else{
        joint_rank = length(compJ)
      }
      pval_joint = permuteJoint$pvalues
    }else{
      # Joint rank is zero when one of the ranks is zero
      joint_rank = 0
      pval_joint = NA
    }
    output$alpha[[l]] = list(rankX = rankX, rankY = rankY, joint_rank = joint_rank, pvalX = permuteX$pvalues, pvalY = permuteY$pvalues, pval_joint = pval_joint, alpha = alpha)
  }
  output
}
# stop cluster
stopCluster(cl)

save(results, file = "SimResults_ForRanks.Rdata")


######################################################################################
# Only joint rank estimation based on the full saturated model
######################################################################################

# Threshold for p-value significance
alpha = 0.01
# Number of permutations
nperms = 1000

nworkers <- detectCores() 
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)

set.seed(234623)

results <- foreach(i=1:nrep, .errorhandling='pass') %dorng% {
  # Generate the data with the specified signal and noise ratio
  data <- generateData_v2(nsubject = 48, snr = c(snr1[i], snr2[i]), vars = c(var_level, var_level))
  
  # Prepare the output list
  output <- list(snr1 = snr1[i], snr2 = snr2[i], pval_joint = NULL, rank_joint = NULL)
    
    #########################################################################################
    # Apply separate JB to each (NGCA model) with full number of components
    #########################################################################################
    # JB on X
    estX_JB = mlcaFP(xData = t(data$dX), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
    
    # Mixing matrix for X
    Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX))
    
    # JB on Y
    estY_JB = mlcaFP(xData = t(data$dY), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
    
    # Mixing matrix for Y
    My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))
    
    
    #########################################################################################
    # Determine joint rank
    ###################################
    
    # Match the mixing matrices, the first two should be joint
    matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
    
    permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
    
    joint_rank = min(which(permuteJoint$pvalues > alpha)) - 1
    if (joint_rank < 0){joint_rank = 0}
    output$joint_rank = joint_rank
    output$pval_joint = permuteJoint$pvalues
  
  output
}
# stop cluster
stopCluster(cl)

save(results, file = "SimResults_ForRanks_JointFull.Rdata")
