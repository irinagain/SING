# Simulation setting 1
# Generates all the figures and results


######################################################
# Component pictures for Simulation Setting 1 for SING paper
######################################################

# Data generation functions
source("jngcaFunctions.R")

# Generate data, will return both components and mixing matrix
data <- generateData_v2(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))

# Components for X
lgrid = 33

pdf(file = "ComponentsX.pdf", width = 12, height = 5)
par(mfrow = c(1,3))
image(matrix(data$siX[1,], lgrid, lgrid), xaxt = "n", yaxt = "n")
image(matrix(data$sjX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(matrix(data$sjX[2,], lgrid, lgrid), xaxt = "n", yaxt = "n")
dev.off()

# Components for Y
pdf(file = "ComponentsY.pdf", width = 12, height = 5)
par(mfrow = c(1,4))
image(vec2net(data$sjY[1,]), xaxt = "n", yaxt = "n")
image(vec2net(data$sjY[2,]), xaxt = "n", yaxt = "n")
image(vec2net(data$siY[1,]), xaxt = "n", yaxt = "n")
image(vec2net(data$siY[2,]), xaxt = "n", yaxt = "n")
dev.off()

######################################################
# Component selection results
#####################################################
load("SimResults_ForRanks.Rdata")
# The design repeated
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

# Threshold for p-value significance
alphavec = c(0.01, 0.05, 0.1)

# Create dataframe with all results

data_ranks = c()

for (i in 1:length(alphavec)){
  data_ranks = rbind(data_ranks, data.frame(snr1 = sapply(results, function(x) x$snr1), snr2 = sapply(results, function(x) x$snr2), alpha = rep(alphavec[i], length(results)), rankX = sapply(results, function(x) x$alpha[[i]]$rankX), rankY = sapply(results, function(x) x$alpha[[i]]$rankY), joint_rank = sapply(results, function(x) x$alpha[[i]]$joint_rank)))
}

# Adjust for Inf joint rank, it should be 1 (mistake in calculation, fixed now)
data_ranks$joint_rank[data_ranks$joint_rank == Inf] = 1

# Now need to add full rank results from oversaturated model
load("SimResults_ForRanks_JointFull.Rdata")


data_ranks = rbind(data_ranks, data.frame(snr1 = sapply(results, function(x) x$snr1), snr2 = sapply(results, function(x) x$snr2), alpha = rep("Full", length(results)), rankX = rep(48, length(results)), rankY = rep(48, length(results)), joint_rank = sapply(results, function(x) x$joint_rank)))


data_ranks$snr1 <- paste("SNRx = ", data_ranks$snr1, sep="")
data_ranks$snr2 <- paste("SNRy = ", data_ranks$snr2, sep="")
data_ranks$var1 <- paste("var1 = ", data_ranks$var1, sep="")
data_ranks$var2 <- paste("var2 = ", data_ranks$var2, sep="")

# Count how many times the rank identification based on saturated model made mistakes
library(tidyverse)
test = data_ranks %>% filter(alpha == "Full")
table(test$snr1, test$snr2, test$joint_rank)

# Create rank figure

data_melt = reshape2::melt(data_ranks, measure.vars = c("rankX", "rankY", "joint_rank"))

library(ggplot2)



levels(data_melt$variable) = c(expression("r"["x"]), expression("r"["y"]), expression("r"["J"]))
data_melt$snr1 <- as.factor(data_melt$snr1)
data_melt$snr2 <- as.factor(data_melt$snr2)
levels(data_melt$snr1) = c(expression(Low~SNR~X), expression(High~SNR~X))
levels(data_melt$snr2) = c(expression(Low~SNR~Y), expression(High~SNR~Y))



pdf(file = "Boxplots_Ranks_ForPaper.pdf", width = 14, height = 12)
p = data_melt%>%
  filter(!is.na(value))%>%
  filter(value < 48)%>%
  ggplot(aes(x = value, group = alpha, fill = alpha)) +geom_histogram(position = "dodge", binwidth = 0.7, col = "black") +facet_grid(variable ~ snr1 + snr2, labeller = label_parsed) + xlab("") + theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(), legend.text.align = 0)+ scale_fill_brewer(palette="GnBu", labels = expression(paste(alpha, '=', 0.01), paste(alpha, '=', 0.05), paste(alpha, '=', 0.1), Full))+scale_x_continuous(breaks = c(0,1,2,3,4,5,6))+xlab("Selected rank value") + ylab("Frequency")
print(p)
dev.off()


######################################################
# Estimation results
######################################################
load("SimResults_estimation_small.Rdata")

# Create dataframe with all results

dataSep <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBsep$errorSx), errorSy = sapply(results, function(x) x$JBsep$errorSy), errorMx = sapply(results, function(x) x$JBsep$errorMx), errorMy = sapply(results, function(x) x$JBsep$errorMy), method = rep("Sep", nrep), errorM = sapply(results, function(x) x$JBsep$errorM), Mdist = sapply(results, function(x) x$JBsep$Mdist))

dataSmall <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$small$errorSx), errorSy = sapply(results, function(x) x$JBjoint$small$errorSy), errorMx = sapply(results, function(x) x$JBjoint$small$errorMx), errorMy = sapply(results, function(x) x$JBjoint$small$errorMy), method = rep("Small rho", nrep), errorM = sapply(results, function(x) x$JBjoint$small$errorM), Mdist = sapply(results, function(x) x$JBjoint$small$Mdist))

dataMed <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$medium$errorSx), errorSy = sapply(results, function(x) x$JBjoint$medium$errorSy), errorMx = sapply(results, function(x) x$JBjoint$medium$errorMx), errorMy = sapply(results, function(x) x$JBjoint$medium$errorMy), method = rep("Medium rho", nrep), errorM = sapply(results, function(x) x$JBjoint$medium$errorM), Mdist = sapply(results, function(x) x$JBjoint$medium$Mdist))

dataLarge <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$large$errorSx), errorSy = sapply(results, function(x) x$JBjoint$large$errorSy), errorMx = sapply(results, function(x) x$JBjoint$large$errorMx), errorMy = sapply(results, function(x) x$JBjoint$large$errorMy), method = rep("Large rho", nrep), errorM = sapply(results, function(x) x$JBjoint$large$errorM), Mdist = sapply(results, function(x) x$JBjoint$large$Mdist))

# Dataframe for joint ICA
dataJointICA <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$jointICA$errorSx), errorSy = sapply(results, function(x) x$jointICA$errorSy), errorMx = sapply(results, function(x) x$jointICA$errorM), errorMy = sapply(results, function(x) x$jointICA$errorM), method = rep("Joint ICA", nrep), errorM = sapply(results, function(x) x$jointICA$errorM), Mdist = rep(0, length(snr1)))

dataAll <- rbind(dataSep, dataSmall, dataMed, dataLarge, dataJointICA)

dataAll$snr1 <- paste("SNRx = ", dataAll$snr1, sep="")
dataAll$snr2 <- paste("SNRy = ", dataAll$snr2, sep="")

dataAll$snr1 <- as.factor(dataAll$snr1)
dataAll$snr2 <- as.factor(dataAll$snr2)
levels(dataAll$snr1) = c(expression(Low~SNR~X), expression(High~SNR~X))
levels(dataAll$snr2) = c(expression(Low~SNR~Y), expression(High~SNR~Y))


colnames(dataAll)[3] = "Sx"
colnames(dataAll)[4] = "Sy"
colnames(dataAll)[5] = "Mx"
colnames(dataAll)[6] = "My"

levels(dataAll$method) <- c(expression(rho~"="~0), expression(small~rho), expression(medium~rho), expression(large~rho), expression(joint~ICA))

# All methods together
library(ggplot2)

dataAllmelt = dataAll%>%
  reshape2::melt(variable.name = "type")%>%
  filter((type != "errorM")&(type != "Mdist"))

levels(dataAllmelt$type) <- c(expression("S"["Jx"]), expression("S"["Jy"]), expression("M"["Jx"]), expression("M"["Jy"]), expression(M), expression(distM))

p = dataAllmelt%>%
  filter((type != "M")&(type != "distM"))%>%
  ggplot(aes(x = method, y = value, fill = method)) +geom_boxplot() +facet_grid(type ~ snr1 + snr2, labeller = label_parsed) + xlab("") + ylab("Error") + theme(legend.position="none", text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1, size = 24))+scale_fill_brewer(palette="GnBu") + scale_x_discrete(labels = parse(text = levels(dataAllmelt$method)))
print(p)

pdf(file = "Boxplots_WithJoint_ForPaper.pdf", width = 14, height = 12)
print(p)
dev.off()


# All methods without Joint ICA
p2 = dataAllmelt%>%
  filter(method != "joint ~ ICA")%>%
  ggplot(aes(x = method, y = value, fill = method)) +geom_boxplot() +facet_grid(type ~ snr1 + snr2, labeller = label_parsed) + xlab("") + ylab("Error") + theme(legend.position="none", text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1, size = 24))+scale_fill_brewer(palette="GnBu")+ scale_x_discrete(labels = parse(text = levels(dataAllmelt$method)))
print(p2)

pdf(file = "Boxplots_WithoutJoint_ForPaper.pdf", width = 14, height = 12)
print(p2)
dev.off()

# Get numbers for R^2_{jx} and R^2_{jy}

dataR2 <- data.frame(snr1 = snr1, snr2 = snr2, R2x = sapply(results, function(x) x$R2x), R2y = sapply(results, function(x) x$R2y))

summary(dataR2$R2x[dataR2$snr1 == 0.2]) #[0.09, 0.13]
summary(dataR2$R2x[dataR2$snr1 == 5]) #[0.4, 0.65]
summary(dataR2$R2y[dataR2$snr2 == 0.2]) #[0.15, 0.16] 
summary(dataR2$R2y[dataR2$snr2 == 5]) #[0.72, 0.80]
