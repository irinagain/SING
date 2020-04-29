# Simulation Setting 2
# Generate data from SING model with components based on HCP data

# Load components to use for data generation
load("ComponentsForSim_setting2.Rda")
nindiv = 10 # number of inidividual components in each dataset


set.seed(2358234)
# Generate mixing matrices
nsubject = 48
n1 = round(nsubject/2)
# Mixing for joint components
mj1 = c(rep(1,n1),rep(-1,nsubject-n1))+rnorm(nsubject)
mj2 = c(rep(-1,n1),rep(1,nsubject - n1))+rnorm(nsubject)
mj = cbind(mj1,mj2)


# Generate data with 2 joint components, and 10 individual component for each dataset
#################################################################################
# Set signal to noise ratio
snr = c(0.5, 0.5) # Ratio of Non-gaussian to Gaussian
# Create joint structure for X and Y
djX = mj%*%rbind(sj1x, sj2x)
scalemj = t(t(mj)*c(-5,2))
djY = scalemj%*%rbind(sj1y, sj2y)

# Create individual structure for X and Y
# Here the individual mixing matrices are approximately uncorrelated with joint -> adjust?
set.seed(928364)
miX = matrix(rnorm(nsubject * nindiv, sd = 10), nsubject, nindiv)
diX = miX%*%sIxs
miY = matrix(rnorm(nsubject * nindiv, sd = 30), nsubject, nindiv)
diY = miY%*%sIys

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
save(dX, dY, mj, sj1x, sj1y, sj2x, sj2y, sIxs, sIys, file = "SimulatedLargeScaleindiv10.Rdata")
