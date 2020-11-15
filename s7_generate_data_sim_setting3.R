# Simulation Setting 3
# Generate data from SING model with components based on HCP data
# Aadd additional sparsity to the components

# Load components to use for data generation
load("ComponentsForSim_setting2.Rda")
nindiv = 10 # number of inidividual components in each dataset

# Check the current component distribution
# threshold
prop_zero <- function(mat, threshold){
  rowMeans(abs(mat) <= threshold)
}


prop_zero(rbind(sj1x, sj2x, sIxs), 5) # gives around 99.9% zeroes

prop_zero(rbind(sj1y, sj2y, sIys), 5) # gives around 99.8% zeroes

# CAREFUL: if just set to 0, may not longer be uncorrelated, so will have to adjust for that manually. Plus the mean will no longer be zero.
threshold = 5
Sx = rbind(sj1x, sj2x, sIxs)
Sy = rbind(sj1y, sj2y, sIys)
Sx[abs(Sx) < threshold] = 0
Sy[abs(Sy) < threshold] = 0
sum(colSums(Sx) == 0) #58710 are completely zero
sum(colSums(Sy) == 0) #58710 are completely zero

# Adjust the mean of S to be zero without changing the sparsity pattern
mean_adjust <- function(S){
  sumX = rowSums(S)
  nonzeroX = rowSums(S != 0)
  whichzerox = which(S == 0)
  S = S - matrix(sumX/nonzeroX, nrow(S), ncol(S))
  S[whichzerox] = 0
  return(S)
}

# Apply the above function
Sx = mean_adjust(Sx)
Sy = mean_adjust(Sy)

# Adjust the scale 1st time
scale_adjust <- function(S){
  p = ncol(S)
  currents = rowSums(S^2)
  S = S * sqrt(p) / sqrt(currents)
  return(S)
}

Sx = scale_adjust(Sx)
Sy = scale_adjust(Sy)

# Adjust so that they are uncorrelated again without changing the sparsity pattern


# For Sx, the only remaining correlation is between s1 and s10
# Adjust S1 to be orthogonal to S10?
nonzeroj = which(Sx[1, ] != 0)
Sx[1, nonzeroj] = Sx[1, nonzeroj, drop = F] - Sx[1, nonzeroj, drop = F] %*% crossprod(Sx[12, nonzeroj, drop = F])/ sum(Sx[12, nonzeroj]^2)

Sx = scale_adjust(Sx)


cor_adjust <- function(S){
  r = nrow(S)
  for (j in 1:(r - 1)){
    nonzeroj = which(S[j, ] != 0)
    indexcor = which(S[j, nonzeroj] %*% t(S[, nonzeroj]) != 0)
    indexcor = indexcor[indexcor != j]
    if (length(indexcor) > 0){
      for (l in indexcor){
        #S[j, nonzeroj] = S[j, nonzeroj] - S[j, nonzeroj] %*% t(S[indexcor, nonzeroj, drop = F])%*% solve(tcrossprod(S[indexcor, nonzeroj, drop = F]), S[indexcor, nonzeroj, drop = F])
        S[j, nonzeroj] = S[j, nonzeroj] - S[j, nonzeroj, drop = F] %*% crossprod(S[l, nonzeroj, drop = F])/ sum(S[l, nonzeroj]^2)
      }
    }
  }
  return(S)
}
Sy = cor_adjust(Sy)

Sy = scale_adjust(Sy)

# What is the difference between original and new
# Source all the functions
source("jngcaFunctions.R")

errorSxJ = frobICA(S1 = t(Sx[1:2, ]), S2 = t(rbind(sj1x, sj2x)), standardize  = T)
errorSxJ # 1.25

errorSyJ = frobICA(S1 = t(Sy[1:2, ]), S2 = t(rbind(sj1y, sj2y)), standardize  = T)
errorSyJ # 1.14

set.seed(2358234)
# Generate mixing matrices
nsubject = 48
n1 = round(nsubject/2)
# Mixing for joint components
mj1 = c(rep(1, n1), rep(-1, nsubject - n1)) + rnorm(nsubject)
mj2 = c(rep(-1, n1), rep(1, nsubject - n1)) + rnorm(nsubject)
mj = cbind(mj1, mj2)


# Generate data with 2 joint components, and 10 individual component for each dataset
#################################################################################
# Set signal to noise ratio
snr = c(0.5, 0.5) # Ratio of Non-gaussian to Gaussian
# Create joint structure for X and Y
djX = mj%*%rbind(Sx[1, ], Sx[2, ])
scalemj = t(t(mj)*c(-5,2))
djY = scalemj%*%rbind(Sy[1, ], Sy[2, ])

# Create individual structure for X and Y
# Here the individual mixing matrices are approximately uncorrelated with joint -> adjust?
set.seed(928364)
miX = matrix(rnorm(nsubject * nindiv, sd = 10), nsubject, nindiv)
diX = miX%*%Sx[-c(1:2), ]
miY = matrix(rnorm(nsubject * nindiv, sd = 30), nsubject, nindiv)
diY = miY%*%Sy[-c(1:2), ]

# Calculate Frobenius norm of the signal
signalXF2 = sum((djX + diX)^2)
signalYF2 = sum((djY + diY)^2)

# Generate noise for X and Y 
set.seed(120836)
px = ncol(djX)
nX = t(scale(matrix(rnorm((nsubject-12)*px),px)))
mnX = matrix(rnorm((nsubject-12)*nsubject),nsubject) # change to gaussian
dnX = mnX%*%nX

py = ncol(djY)
nY = t(scale(matrix(rnorm((nsubject-12)*py),py)))
mnY = matrix(rnorm((nsubject-12)*nsubject),nsubject)
dnY = mnY%*%nY

# Adjust the noise with snr ratio
# Wrt to Frobenius norm
dnX = dnX * sqrt(signalXF2/(sum(dnX^2)*snr[1]))
dnY = dnY * sqrt(signalYF2/(sum(dnY^2)*snr[2]))

# Create data matrix X and Y
dX = djX + diX + dnX
dY = djY + diY + dnY

# Signal all to noise
signalXF2/sum(dnX^2)
signalYF2/sum(dnY^2)

# Signal joint to all (aim to have 0.0017 for X and 0.0021 for Y which mimics real data)
sum(djX^2)/sum(dX^2) # 0.0014 (close enough)
sum(djY^2)/sum(dY^2) # 0.0021 (close enough)

# Save generated data
save(dX, dY, mj, Sx, Sy, file = "SimulatedLargeScaleindiv10_Sparse.Rdata")
