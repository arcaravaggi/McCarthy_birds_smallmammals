#### Small Mammal Models ####

# Clear R's memory
rm(list=ls())

# Load required packages
library(dplyr)
library(caret)
library(e1071)
library(ggplot2)

# Set working directory
setwd("C:\\Users\\mccarthya\\OneDrive - University College Cork\\Hen Harriers\\Hen Harrier MSc\\Prey Abundance Surveys\\Data\\Data Analysis\\Data & code for publication\\Data")

# Load vegetation data
veg.b <- read.csv("vegetation_data.csv")

## Principle Component Analysis ##
veg.pca <- preProcess(veg.b[3:19],
                         method=c("BoxCox", "center", "scale", "pca"),# Correct for skew and standardise 
                         thresh = .9) # retain components which account for 90% of variance
c500 <- predict(veg.pca, veg.b[3:19]) # components

# Calculate % variation explained by each PC
pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}

pcaCharts(prcomp(c500, center=FALSE))

# Loadings
l500 <- data.frame(veg.pca$rotation)

sm.veg <- data.frame(veg.b[1:2], c500)

# Plot the PCs against treatment and vegetation variables
l500$var <- rownames(l500)

# Load dplyr (function overwritten by previously loaded package)
library(dplyr)

# Plot treatment against PC1 and PC2
sm.veg %>% group_by(treatment) %>% select (-c(location)) %>%
  summarize_all(funs(mean, sd)) %>%
  ggplot(aes(x=PC1_mean, y=PC2_mean)) +
  geom_point() +
  xlim(-6,6) + ylim(-5,5) +
  theme_classic() +
  geom_errorbar(aes(ymin=PC2_mean-PC2_sd, ymax=PC2_mean+PC2_sd), width=.1) +
  geom_errorbarh(aes(xmin=PC1_mean-PC1_sd, xmax=PC1_mean+PC1_sd)) +
  geom_point(data=sm.veg, aes(x=PC1, y=PC2, shape = treatment, colour = treatment, size = 0.5))


# Load small mammal vegetation surveyor
veg.surv <- read.csv("Small_mammal_veg_surveyors.csv")

# Merge back in vegetation surveyor
names(veg.surv) <- c("location", "treatment", "veg.surveyor")
veg.c <- merge(sm.veg, veg.surv, by = c("location", "treatment"))
# Delete replicates of each row
veg.c <- veg.c[!duplicated(veg.c), ]

# Read in small mammal trapping data
cap.hist.sum <- read.csv("Capture_history.csv")
names(cap.hist.sum)
names(veg.c)
# Merge with vegetation data
veg_cap.hist <- merge(veg.c, cap.hist.sum, by = c("location", "treatment"))

# Make vegetation surveyor and mammal surveyor column
names(veg_cap.hist)
veg_cap.hist$Surveyor <- paste(veg_cap.hist[,11], veg_cap.hist[,14], sep="_")

# Tidy duplicated surveyor names
veg_cap.hist$Surveyor <- as.character(veg_cap.hist$Surveyor)
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "IS_IS"] <- "IS"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "AMcC_AMcC"] <- "AMc"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "MOC_MOC"] <- "MOC"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "JK_JK"] <- "JK"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "TK_TK"] <- "TK"


###################
#### Modelling ####

# Zero-inflated negative binomial general linear mixed model
# Captures = Response variable
# PCs = Fixed effects
# Treatment, location, surveyor = Random effects

# Load required packages 
library(glmmTMB)
library(TMB)
library(DHARMa)


# Make new dataframes with each species for all sites
bvcaptures <- veg_cap.hist[veg_cap.hist$species == "BV",]
wmcaptures <- veg_cap.hist[veg_cap.hist$species == "FM",]
gscaptures <- veg_cap.hist[veg_cap.hist$species == "GS",]


###############################
## Treatment as fixed effect ##

## Basic model of captures vs. treatment ##
ziglmm_basic <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                    (1|location/treatment), data = veg_cap.hist, family = nbinom2)
summary(ziglmm_basic)
# Simulate the residuals
simres_basic <- simulateResiduals(ziglmm_basic)
# Plot the simulated residuals
plot(simres_basic, rank = T)

# Run model checks
testUniformity(simres_basic)
testDispersion(simres_basic)
testZeroInflation(simres_basic)

# Check for overdispersion
overdisp_fun <- function(ziglmm_basic) {
  rdf <- df.residual(ziglmm_basic)
  rp <- residuals(ziglmm_basic,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmm_basic)


# Null model comparison
ziglmm_basic_null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                         (1|location/treatment), data = veg_cap.hist, family = nbinom2)
summary(ziglmm_basic_null)
# Simulate the residuals
simres_basic_null <- simulateResiduals(ziglmm_basic_null)
# Plot the simulated residuals
plot(simres_basic_null, rank = T)
# Compare the model with the null model
anova(ziglmm_basic, ziglmm_basic_null)
# Model model is a better fit

# Re-order to compare early pre-thicket and late pre-thicket
# Do this by renaming Y as A, therefore A will be the intercept in the output
veg_cap.hist.reorder <- veg_cap.hist
veg_cap.hist.reorder$treatment <- gsub("Y", "A", veg_cap.hist.reorder$treatment)

# Compare late pre-thicket and early pre-thicket (now A)
ziglmm_basic2 <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                          (1|location/treatment), data = veg_cap.hist.reorder, family = nbinom2)
summary(ziglmm_basic2)




#### Model by each species ####
## Bank vole only ##
ziglmmbv_tr <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                      (1|location/treatment), data = bvcaptures, family = nbinom2)
summary(ziglmmbv_tr)
# Simulate the residuals
simres.bv_tr <- simulateResiduals(ziglmmbv_tr)
# Plot the simulated residuals
plot(simres.bv_tr, rank = T)

# Run model checks
testUniformity(simres.bv_tr)
testDispersion(simres.bv_tr)
testZeroInflation(simres.bv_tr)

# Check for overdispersion
overdisp_fun <- function(ziglmmbv_tr) {
  rdf <- df.residual(ziglmmbv_tr)
  rp <- residuals(ziglmmbv_tr,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmmbv_tr)


# Null model comparison
ziglmmbv_tr_null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                           (1|location/treatment), data = bvcaptures, family = nbinom2)
summary(ziglmmbv_tr_null)
# Simulate the residuals
simres.bv_tr_null <- simulateResiduals(ziglmmbv_tr_null)
# Plot the simulated residuals
plot(simres.bv_tr_null, rank = T)
# Compare the model with the null model
anova(ziglmmbv_tr,ziglmmbv_tr_null)
# Main model is a better fit

# Re-order to compare early pre-thicket and late pre-thicket
bvcaptures.reorder <- bvcaptures
bvcaptures.reorder$treatment <- gsub("Y", "A", bvcaptures.reorder$treatment)
# Compare late pre-thicket and early pre-thicket (now A)
ziglmmbv_tr2 <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                         (1|location/treatment), data = bvcaptures.reorder, family = nbinom2)
summary(ziglmmbv_tr2)




## Wood mouse only ##
ziglmmwm_tr <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                         (1|location/treatment), data = wmcaptures, family = nbinom2)
summary(ziglmmwm_tr)
# Simulate the residuals
simres.wm_tr <- simulateResiduals(ziglmmwm_tr)
# Plot the simulated residuals
plot(simres.wm_tr, rank = T)

# Run model checks
testUniformity(simres.wm_tr)
testDispersion(simres.wm_tr)
testZeroInflation(simres.wm_tr)

# Check for overdispersion
overdisp_fun <- function(ziglmmwm_tr) {
  rdf <- df.residual(ziglmmwm_tr)
  rp <- residuals(ziglmmwm_tr,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmmwm_tr)


# Null model comparison
ziglmmwm_tr_null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                              (1|location/treatment), data = wmcaptures, family = nbinom2)
summary(ziglmmwm_tr_null)
# Simulate the residuals
simres.wm_tr_null <- simulateResiduals(ziglmmwm_tr_null)
# Plot the simulated residuals
plot(simres.wm_tr_null, rank = T)
# Compare model with null model
anova(ziglmmwm_tr,ziglmmwm_tr_null)
# Main model is a better fit

# Re-order to compare early pre-thicket and late pre-thicket
wmcaptures.reorder <- wmcaptures
wmcaptures.reorder$treatment <- gsub("Y", "A", wmcaptures.reorder$treatment)
# Compare late pre-thicket and early pre-thicket (now A)
ziglmmwm_tr2 <- glmmTMB(captures ~ treatment + (1|Surveyor) +
                         (1|location/treatment), data = wmcaptures.reorder, family = nbinom2)
summary(ziglmmwm_tr2)




##########################################
## Principle Components as fixed effect ##

# Run stepwise regression to determine which model is the best fit
step(glmmTMB(captures ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + (1|Surveyor) +
                      (1|location/treatment), data = veg_cap.hist, family = nbinom2), direction = "both")
# PC1 + PC2 are best fit
# Run model with just these two PCs

ziglmm <- glmmTMB(captures ~ PC1 + PC2 + (1|Surveyor) +
                    (1|location/treatment), data = veg_cap.hist, family = nbinom2)
summary(ziglmm)
# Simulate the residuals
simres <- simulateResiduals(ziglmm)
# Plot the simulated residuals
plot(simres, rank = T)

# Run model checks
testUniformity(simres)
testDispersion(simres)
testZeroInflation(simres)

# Check for overdispersion
overdisp_fun <- function(ziglmm) {
  rdf <- df.residual(ziglmm)
  rp <- residuals(ziglmm,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmm)


# Null model comparison
ziglmm.null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                    (1|location/treatment), data = veg_cap.hist, family = nbinom2)
summary(ziglmm.null)
# Simulate the residuals
simres.null <- simulateResiduals(ziglmm.null)
# Plot the simulated residuals
plot(simres.null, rank = T)
# Compare model with null model
anova(ziglmm, ziglmm.null)
# Main model is a better fit




#### Model by each species ####
## Bank vole only ##
ziglmmbv <- glmmTMB(captures ~ PC1 + PC2 + (1|Surveyor) +
                    (1|location/treatment), data = bvcaptures, family = nbinom2)
summary(ziglmmbv)
# Simulate the residuals
simres.bv <- simulateResiduals(ziglmmbv)
# Plot the simulated residuals
plot(simres.bv, rank = T)

# Run model checks
testUniformity(simres.bv)
testDispersion(simres.bv)
testZeroInflation(simres.bv)

# Check for overdispersion
overdisp_fun <- function(ziglmmbv) {
  rdf <- df.residual(ziglmmbv)
  rp <- residuals(ziglmmbv,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmmbv)


# Null model comparison
ziglmmbv_null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                      (1|location/treatment), data = bvcaptures, family = nbinom2)
summary(ziglmmbv_null)
# Simulate the residuals
simres.bv_null <- simulateResiduals(ziglmmbv_null)
# Plot the simulated residuals
plot(simres.bv_null, rank = T)
# Compare model with the null model
anova(ziglmmbv,ziglmmbv_null)
# Near significant and lower AIC value in the actual model
# Main model is better fit


## Wood mouse only ##
ziglmmwm <- glmmTMB(captures ~ PC1 + PC2 + (1|Surveyor) +
                      (1|location/treatment), data = wmcaptures, family = nbinom2)
summary(ziglmmwm)
# Simulate the residuals
simres.wm <- simulateResiduals(ziglmmwm)
# Plot the simulated residuals
plot(simres.wm, rank = T)

# Run model checks
testUniformity(simres.wm)
testDispersion(simres.wm)
testZeroInflation(simres.wm)

# Check for overdispersion
overdisp_fun <- function(ziglmmwm) {
  rdf <- df.residual(ziglmmwm)
  rp <- residuals(ziglmmwm,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(ziglmmwm)



# Null model comparison
ziglmmwm_null <- glmmTMB(captures ~ 1 + (1|Surveyor) +
                           (1|location/treatment), data = wmcaptures, family = nbinom2)
summary(ziglmmwm_null)
# Simulate the residuals
simres.wm_null <- simulateResiduals(ziglmmwm_null)
# Plot the simulated residuals
plot(simres.wm_null, rank = T)
# Compare the model with the null model
anova(ziglmmwm,ziglmmwm_null)
# Main model is a better fit
