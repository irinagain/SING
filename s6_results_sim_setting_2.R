# Simulation setting 2
# Error calculation for all methods, and creation of figures

# Load the results of joint ICA
load("jointICA_LargeScaleindiv10.Rda")
# Load the results from separate approach (rho = 0)
load("sepJB_LargeScale.Rda")
# Load the results from SING (from small rho to large rho)
load("out_indiv_small.Rda")
load("out_indiv_medium.Rda")
load("out_indiv_large.Rda")

# Load true underlying components
load("SimulatedLargeScaleindiv10.Rdata")

# Source all the functions
source("jngcaFunctions.R")

n = nrow(dX)
pX = ncol(dX)
pY = ncol(dY)
# For X
# Center
dXsA <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
# Scale rowwise
est.sigmaXA = tcrossprod(dXsA)/(pX-1)  ## since already centered, can just take tcrossprod
whitenerXA = est.sigmaXA%^%(-0.5)
invLx = est.sigmaXA%^%(0.5)
xDataA = whitenerXA %*% dXsA

# For Y
# Center
dYsA <- dY - matrix(rowMeans(dY), n, pY, byrow = F)
# Scale rowwise
est.sigmaYA = tcrossprod(dYsA)/(pY-1)  ## since already centered, can just take tcrossprod
whitenerYA = est.sigmaYA%^%(-0.5)
yDataA = whitenerYA %*% dYsA
invLy = est.sigmaYA%^%(0.5)

# Starting points for the algorithm are Us
###########################################
Uxfull <- output_sepJB$estX_JB$Ws
Mx_JB = est.M.ols(sData = output_sepJB$estX_JB$S, xData = t(dX))
Uyfull <- output_sepJB$estY_JB$Ws
My_JB = est.M.ols(sData = output_sepJB$estY_JB$S, xData = t(dY))

matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))

cor(matchMxMy$Mx[, 1:2], matchMxMy$My[, 1:2]) # 0.95 and 0.93

# Calculate errors for Sx
############################
#joint ICA
errorSx = frobICA(S1 = t(rbind(sj1x, sj2x)), S2 = output_jointICA$estXY_JB$S[1:pX, ], standardize  = T)
errorSx # 1.408

# rho = 0 with all 2 matched components
Sxmatched = matchMxMy$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(rbind(sj1x, sj2x)), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.051

# small rho
Sxmatched = out_indiv_small$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(rbind(sj1x, sj2x)), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.035

# medium rho
Sxmatched = out_indiv_medium$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(rbind(sj1x, sj2x)), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.035

# large rho
Sxmatched = out_indiv_large$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(rbind(sj1x, sj2x)), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.035

# Calculate errors for Sy
############################
#joint ICA
errorSy = frobICA(S1 = t(rbind(sj1y, sj2y)), S2 = output_jointICA$estXY_JB$S[(pX+1):(pX+pY), ], standardize  = T)
errorSy # 1.409


# rho = 0 with all 2 matched components
Symatched = matchMxMy$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(rbind(sj1y, sj2y)), S2 = t(Symatched), standardize  = T)
errorSy # 0.031, same

# small rho 
Symatched = out_indiv_small$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(rbind(sj1y, sj2y)), S2 = t(Symatched), standardize  = T)
errorSy # 0.023


# medium rho
Symatched = out_indiv_medium$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(rbind(sj1y, sj2y)), S2 = t(Symatched), standardize  = T)
errorSy # 0.023

# large rho
Symatched = out_indiv_large$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(rbind(sj1y, sj2y)), S2 = t(Symatched), standardize  = T)
errorSy # 0.023

# Calculate errors for M
############################
# joint ICA
dXcentered <- dX - matrix(rowMeans(dX), n, pX, byrow = F)
dYcentered <- dY - matrix(rowMeans(dY), n, pY, byrow = F)
# Normalization step here. Divide by \|X\|_F^2 each part
dXsA <- dXcentered/sqrt(mean(dXcentered^2))
dYsA <- dYcentered/sqrt(mean(dYcentered^2))
# Concatenate them together [X, Y] and perform SVD (PCA)
dXY <- cbind(dXsA, dYsA) # [X, Y] ~ UDV'
Mjoint = est.M.ols(sData = output_jointICA$estXY_JB$S, xData = t(dXY))
errorM = frobICA(M1 = Mjoint, M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 1.362


# rho = 0 with on matched
Mxjoint = tcrossprod(invLx, matchMxMy$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.295

Myjoint = tcrossprod(invLy, matchMxMy$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.208

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag # 0.95 and 0.93
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.171


# small rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_small$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.118

Myjoint = tcrossprod(invLy, out_indiv_small$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.121

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.9994; 0.9999
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.12

# medium rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_medium$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.115

Myjoint = tcrossprod(invLy, out_indiv_medium$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.115

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.99999; 0.999998
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.115

# large rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_large$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.114

Myjoint = tcrossprod(invLy, out_indiv_large$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.114

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.999999999; 1
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.114



# Save true components, estimated from jointICA, estimated with rho=0, estimated with small rho
Sxtrue = t(rbind(sj1x, sj2x))
Sytrue = t(rbind(sj1y, sj2y))
SxjointICA = output_jointICA$estXY_JB$S[1:pX, ]
SyjointICA = output_jointICA$estXY_JB$S[(pX+1):(pX + pY), ]
Sx_rho0 = t(matchMxMy$Ux[1:2, ] %*% xDataA)
Sy_rho0 = t(matchMxMy$Uy[1:2, ] %*% yDataA)
Sx_rhoSmall = t(out_indiv_small$Ux[1:2, ] %*% xDataA)
Sy_rhoSmall = t(out_indiv_small$Uy[1:2, ] %*% yDataA)

Sxtrue = signchange(Sxtrue)
Sytrue = signchange(Sytrue)
SxjointICA = signchange(SxjointICA)
SyjointICA = signchange(SyjointICA)
Sx_rhoSmall = signchange(Sx_rhoSmall)
Sy_rhoSmall = signchange(Sy_rhoSmall)

# Save the components
save(Sxtrue, Sytrue, SxjointICA, SyjointICA, Sx_rhoSmall, Sy_rhoSmall, file = "EstimatedComponentsLarge.Rda")

# Create all the plots for joint components
out_true1 = plotNetwork(Sytrue[,1], title='Truth',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv') 

out_true2 = plotNetwork(Sytrue[,2], title='Truth',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv') 

out_joint1 = plotNetwork(SyjointICA[,1], title='joint ICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv') 

out_joint2 = plotNetwork(SyjointICA[,2], title='joint ICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv') 

out_rhoSmall1 = plotNetwork(Sy_rhoSmall[,1], title='small~rho',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv') 

out_rhoSmall2 = plotNetwork(Sy_rhoSmall[,2], title='small~rho',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus.csv')



# Compare loadings from the components of Y shared structure
############################################################################
mmp_modules = read.csv('community_affiliation_mmpplus.csv')
mmp_order = order(mmp_modules$Community_Vector)
Community = factor(mmp_modules$Community_Label)[mmp_order]

# 1st component loadings
dataLoad <- data.frame(x = rep(c(1:379), 3), loadingsum2 = c(out_true1$loadingsummary[mmp_order], out_rhoSmall1$loadingsummary[mmp_order], out_joint1$loadingsummary[mmp_order]), Community = rep(Community, 3), Method = rep(c("Truth", "Small~rho", "Joint~ICA"), each = 379))

dataLoad$Method <- relevel(dataLoad$Method, "Truth")

# 2nd component loadings
dataLoad2 <- data.frame(x = rep(c(1:379), 3), loadingsum2 = c(out_true2$loadingsummary[mmp_order], out_rhoSmall2$loadingsummary[mmp_order], out_joint2$loadingsummary[mmp_order]), Community = rep(Community, 3), Method = rep(c("Truth", "Small~rho", "Joint~ICA"), each = 379))

dataLoad2$Method <- relevel(dataLoad$Method, "Truth")

dataLoad$component <- "Component~1"
dataLoad2$component <- "Component~2"

dataAll <- rbind(dataLoad, dataLoad2)

pdf(file = "LargeScaleBothComponentsY.pdf", width = 14, height = 4)
p = ggplot(dataAll, aes(x = x, y = loadingsum2, col = Community))+geom_point(size = 4) + facet_grid(component~Method, labeller = label_parsed)+xlab('MMP Index')+ylab('L1 Norm of the Rows') + theme(text = element_text(size=20))
print(p)
dev.off()


# Compare visual networks from the components of Y shared structure
############################################################################
qmin=0.005
qmax=0.995
labels = c('VI','SM','DS','VS','DM','CE','SC')
coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)
zmin = -2.5
zmax = 2.5

meltsub = create.graph.long(out_true1$netmat,mmp_order)
meltsub$Method = "Truth"
meltsub$component = "Component~1"
meltsubT = create.graph.long(out_true2$netmat,mmp_order)
meltsubT$Method = "Truth"
meltsubT$component = "Component~2"
meltsub = rbind(meltsub, meltsubT)


meltsub2 = create.graph.long(out_joint1$netmat, mmp_order)
meltsub2$Method = "Joint~ICA"
meltsub2$component = "Component~1"
meltsub2T = create.graph.long(out_joint2$netmat,mmp_order)
meltsub2T$Method = "Joint~ICA"
meltsub2T$component = "Component~2"
meltsub2 = rbind(meltsub2, meltsub2T)


meltsub3 = create.graph.long(out_rhoSmall1$netmat, mmp_order)
meltsub3$Method = "Small~rho"
meltsub3$component = "Component~1"
meltsub3T = create.graph.long(out_rhoSmall2$netmat,mmp_order)
meltsub3T$Method = "Small~rho"
meltsub3T$component = "Component~2"
meltsub3 = rbind(meltsub3, meltsub3T)


meltsubAll = rbind(meltsub, meltsub2, meltsub3)
meltsubAll$Method <- as.factor(meltsubAll$Method)
meltsubAll$Method <- relevel(meltsubAll$Method, "Truth")

g2 = ggplot(meltsubAll, aes(X1, X2,fill=value)) + geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red", limits=c(zmin,zmax), oob=squish)+labs(x = "Node 1", y = "Node 2") + coord_cartesian(clip='off',xlim=c(-0,390)) + facet_grid(component~Method, labeller = label_parsed)

for (i in 1:7) {
  if (i!=3) {
    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
  } else{
    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
  }
}

g2  = g2 + theme(text = element_text(size=24))

# pdf(file = "LargeScaleComponentsYnetwork.pdf", width = 14, height = 6)
# print(g2)
# dev.off()

png(filename = "LargeScaleComponentsYnetwork.png",
    width = 1200, height = 1000)
print(g2)
dev.off()


# Create cifti files:
#############################################################################################
source('makecifti.R')

# match joint ICA:
SxjointICA_matched = matchICA(S = SxjointICA,template = Sxtrue)
cor(SxjointICA_matched,Sxtrue)
# joint components do not contain the joint signal
cor(Sx_rhoSmall,Sxtrue)

apply(Sx_rhoSmall,2,function(x) mean(x^3))
apply(SxjointICA_matched,2,function(x) mean(x^3))

# If necessary run the following:
#Sxtrue = signchange(Sxtrue)
#Sx_rhoLarge = signchange(Sx_rhoLarge)
#SxjointICA = signchange(SxjointICA)
#apply(t(sIxs),2,function(x) mean(x^3))
#SIxs = signchange(t(sIxs))
# This function has dependencies that will need to be modified. It requires matlab and wb_command:
makecifti(Sxtrue,"Sxtrue.dtseries.nii")
makecifti(Sx_rhoSmall,"Sx_rhoSmall.dtseries.nii")
makecifti(SxjointICA_matched,"SxjointICA.dtseries.nii")

# Alternative showing how to specify pathnames on other systems:
#makecifti(Sxtrue,"Sxtrue.dtseries.nii", wbloc='/home/benjamin/Applications2/workbench/bin_linux64/wb_command', matlabloc='/usr/local/bin/matlab')
