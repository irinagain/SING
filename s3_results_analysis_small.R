# Simulation setting 1
# Generates all the figures and results


######################################################
# Component pictures for Simulation Setting 1 for SING paper
######################################################

# Data generation functions
source("jngcaFunctions.R")

# Generate data, will return both components and mixing matrix
set.seed(0573452)
data <- generateData_v2(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))

# Components for X
lgrid = 33

pdf(file = "ComponentsX.pdf", width = 12, height = 5)
par(mfrow = c(1,3))
image(matrix(data$siX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(matrix(data$sjX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(matrix(data$sjX[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n")
dev.off()

# Components for Y
pdf(file = "ComponentsY.pdf", width = 12, height = 5)
par(mfrow = c(1,4))
image(vec2net(data$sjY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(vec2net(data$sjY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(vec2net(data$siY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n")
image(vec2net(data$siY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n")
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
load("SimResults_estimation_small_rho20.Rdata")

# Create dataframe with all results

dataSep <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBsep$errorSx), errorSy = sapply(results, function(x) x$JBsep$errorSy), errorMx = sapply(results, function(x) x$JBsep$errorMx), errorMy = sapply(results, function(x) x$JBsep$errorMy), method = rep("Sep", nrep), errorJx = sapply(results, function(x) x$JBsep$errorJx), errorJy = sapply(results, function(x) x$JBsep$errorJy))

dataSepAve <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBsepAve$errorSx), errorSy = sapply(results, function(x) x$JBsepAve$errorSy), errorMx = sapply(results, function(x) x$JBsepAve$errorMx), errorMy = sapply(results, function(x) x$JBsepAve$errorMy), method = rep("SepAve", nrep), errorJx = sapply(results, function(x) x$JBsepAve$errorJx), errorJy = sapply(results, function(x) x$JBsepAve$errorJy))

dataSmall <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$small$errorSx), errorSy = sapply(results, function(x) x$JBjoint$small$errorSy), errorMx = sapply(results, function(x) x$JBjoint$small$errorMx), errorMy = sapply(results, function(x) x$JBjoint$small$errorMy), method = rep("Small rho", nrep), errorJx = sapply(results, function(x) x$JBjoint$small$errorJx), errorJy = sapply(results, function(x) x$JBjoint$small$errorJy))

dataMed <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$medium$errorSx), errorSy = sapply(results, function(x) x$JBjoint$medium$errorSy), errorMx = sapply(results, function(x) x$JBjoint$medium$errorMx), errorMy = sapply(results, function(x) x$JBjoint$medium$errorMy), method = rep("Medium rho", nrep), errorJx = sapply(results, function(x) x$JBjoint$medium$errorJx), errorJy = sapply(results, function(x) x$JBjoint$medium$errorJy))

dataLarge <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$JBjoint$large$errorSx), errorSy = sapply(results, function(x) x$JBjoint$large$errorSy), errorMx = sapply(results, function(x) x$JBjoint$large$errorMx), errorMy = sapply(results, function(x) x$JBjoint$large$errorMy), method = rep("Large rho", nrep), errorJx = sapply(results, function(x) x$JBjoint$large$errorJx), errorJy = sapply(results, function(x) x$JBjoint$large$errorJy))

# Dataframe for joint ICA
dataJointICA <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$jointICA$errorSx), errorSy = sapply(results, function(x) x$jointICA$errorSy), errorMx = sapply(results, function(x) x$jointICA$errorM), errorMy = sapply(results, function(x) x$jointICA$errorM), method = rep("Joint ICA", nrep), errorJx = sapply(results, function(x) x$jointICA$errorJx), errorJy = sapply(results, function(x) x$jointICA$errorJy))

# Dataframe for mCCA + joint ICA
datamCCA <- data.frame(snr1 = snr1, snr2 = snr2, errorSx = sapply(results, function(x) x$mCCA$errorSx), errorSy = sapply(results, function(x) x$mCCA$errorSy), errorMx = sapply(results, function(x) x$mCCA$errorMx), errorMy = sapply(results, function(x) x$mCCA$errorMy), method = rep("mCCA+jICA", nrep), errorJx = sapply(results, function(x) x$mCCA$errorJx), errorJy = sapply(results, function(x) x$mCCA$errorJy))

dataAll <- rbind(dataSep, dataSepAve, dataSmall, dataMed, dataLarge, dataJointICA, datamCCA)

# Do square root error for Jx and Jy for consistency
dataAll$errorJx = sqrt(dataAll$errorJx)
dataAll$errorJy = sqrt(dataAll$errorJy)

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
colnames(dataAll)[8] = "Jx"
colnames(dataAll)[9] = "Jy"

dataAll$method <- factor(dataAll$method, levels = c("Sep", "Small rho", "Medium rho", "Large rho", "SepAve", "Joint ICA", "mCCA+jICA"))

levels(dataAll$method) <- c(expression(rho~"="~0), expression(Small~rho), expression(Medium~rho), expression(Large~rho), expression(SING-averaged), expression(Joint~ICA), expression(mCCA~"+"~jICA))



# All methods together
library(ggplot2)
library(tidyverse)

dataAllmelt = dataAll%>%
  reshape2::melt(variable.name = "type")



###################################################################################
# Plots of errors
###################################################################################
# no hats
#levels(dataAllmelt$type) <- c(expression("S"["Jx"]), expression("S"["Jy"]), expression("M"["Jx"]), expression("M"["Jy"]), expression("J"["x"]), expression("J"["y"]))

# do hats
levels(dataAllmelt$type) <- c(expression(hat(S)["Jx"]), expression(hat(S)["Jy"]), expression(hat(M)["Jx"]), expression(hat(M)["Jy"]), expression(hat(J)["x"]), expression(hat(J)["y"]))

# All methods without rho = 0, and without medium/small rho
p = dataAllmelt%>%
  filter(method != "rho ~ \"=\" ~ 0")%>%
  filter(method != "Medium ~ rho")%>%
  filter(method != "Small ~ rho")%>%
  ggplot(aes(x = method, y = value, fill = method)) +geom_boxplot() +facet_grid(type ~ snr1 + snr2, labeller = label_parsed) + xlab("") + ylab("Error") + theme(legend.position="none", text = element_text(size=28), axis.text.x = element_text(angle = 45, hjust = 1, size = 30))+scale_fill_brewer(palette="GnBu") 
p = p + scale_x_discrete(labels = parse(text = levels(dataAllmelt$method)[-c(1:3)]))
print(p)

pdf(file = "Boxplots_WithJoint_ForPaper_new.pdf", width = 14, height = 12)
print(p)
dev.off()

# All methods without Joint ICA
p2 = dataAllmelt%>%
  filter(method != "Joint ~ ICA")%>%
  filter(method != "mCCA ~ \"+\" ~ jICA")%>%
  ggplot(aes(x = method, y = value, fill = method)) +geom_boxplot() +facet_grid(type ~ snr1 + snr2, labeller = label_parsed) + xlab("") + ylab("Error") + theme(legend.position="none", text = element_text(size=24), axis.text.x = element_text(angle = 45, hjust = 1, size = 24))+scale_fill_brewer(palette="GnBu")
p2 = p2 + scale_x_discrete(labels = parse(text = levels(dataAllmelt$method)[-c(6,7)]))
print(p2)

pdf(file = "Boxplots_WithoutJoint_ForPaper_new.pdf", width = 14, height = 16)
print(p2)
dev.off()


# Zoomed in on just Ms


# All methods without Joint ICA
p3 = dataAllmelt%>%
  filter(method != "Joint ~ ICA")%>%
  filter(method != "mCCA ~ \"+\" ~ jICA")%>%
  filter(type != "hat(S)[\"Jx\"]")%>%
  filter(type != "hat(S)[\"Jy\"]")%>%
  filter(type != "hat(J)[\"x\"]")%>%
  filter(type != "hat(J)[\"y\"]")%>%
  ggplot(aes(x = method, y = value, fill = method)) +geom_boxplot() +facet_grid(type ~ snr1 + snr2, labeller = label_parsed) + xlab("") + ylab("Error") + theme(legend.position="none", text = element_text(size=28), axis.text.x = element_text(angle = 45, hjust = 1, size = 30))+scale_fill_brewer(palette="GnBu")
p3 = p3 + scale_x_discrete(labels = parse(text = levels(dataAllmelt$method)[-c(6,7)]))
print(p3)

pdf(file = "Boxplots_ZoomedInMs.pdf", width = 14, height = 12)
print(p3)
dev.off()



# Get numbers for R^2_{jx} and R^2_{jy}

dataR2 <- data.frame(snr1 = snr1, snr2 = snr2, R2x = sapply(results, function(x) x$R2x), R2y = sapply(results, function(x) x$R2y))

summary(dataR2$R2x[dataR2$snr1 == 0.2]) #[0.09, 0.13]
summary(dataR2$R2x[dataR2$snr1 == 5]) #[0.4, 0.65]
summary(dataR2$R2y[dataR2$snr2 == 0.2]) #[0.15, 0.16] 
summary(dataR2$R2y[dataR2$snr2 == 5]) #[0.72, 0.80]

###################################################################################
# Alternative - Table of errors
###################################################################################

table_sum = dataAllmelt %>%
  dplyr::group_by(snr1, snr2, method, type)%>%
  summarize(errormean = round(mean(value), 2), errorse = round(sd(value)/10, 3), errorsd = round(sd(value), 2))

# table_sum$error = paste(sprintf('%.3f', table_sum$errormean), "(", sprintf('%.3f', table_sum$errorse), ")", sep = "")
table_sum$error = paste(sprintf('%.2f', table_sum$errormean), "(", sprintf('%.2f', table_sum$errorsd), ")", sep = "")

output_table = table_sum %>%
  select(snr1, snr2, method, type, error)%>%
  pivot_wider(id_cols = c(snr1, snr2, method), names_from = type, values_from = error)


levels(output_table$snr1) = c("Low SNR $\\bX$", "High SNR $\\bX$")
levels(output_table$snr2) = c("Low SNR $\\bY$", "High SNR $\\bY$")

levels(output_table$method) = c("$\\rho = 0$", "small $\\rho$", "medium $\\rho$", "large $\\rho$", "SING-averaged$", "Joint ICA", "mCCA $+$ jICA")

output_table$snr1 = as.character(output_table$snr1)
output_table$snr2 = as.character(output_table$snr2)

# Leave only the first line in each snr regime
output_table[rep(2:7, 4) + rep(c(0:3)*(nrow(output_table)/4), each = 6), 1:2] = ""

colnames(output_table) = c("", "", "Method", "$\\bS_x$", "$\\bS_y$", "$\\bM_x$", "$\\bM_y$", "$\\bJ_x$", "$\\bJ_y$")

#Combine SNRS

output_table[rep(2, 4) + c(0:3)*(nrow(output_table)/4), 1] = output_table[rep(1, 4) + c(0:3)*(nrow(output_table)/4), 2]

library(xtable)
xoutput = xtable(output_table[, -2])

align(xoutput) <-"l|l|l|cc|cc|cc|"

hlines = c(-1, 0, c(1:3)*(nrow(xoutput)/4), nrow(xoutput))

print(xoutput, include.rownames = FALSE, hline.after = hlines, latex.environments = "center", sanitize.text.function = function(x) {x})

