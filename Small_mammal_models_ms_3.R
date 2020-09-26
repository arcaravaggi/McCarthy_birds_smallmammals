#### Small Mammal Vegetation Model ####

# Clear R's memory
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Small Mammal Models ####

# Load required packages
library(dplyr)
library(caret)
library(e1071)
library(ggplot2)
library(factoextra)
library(glmmTMB)
library(MuMIn)
library(car)
library(DHARMa)
library(performance)

# Load vegetation data
veg.b <- read.csv("vegetation_data_4.csv")
names(veg.b)

# Assess correlation between brash measures
cor.test(veg.b$brash.vol,veg.b$brash.dec)
cor.test(veg.b$brash.vol,veg.b$brash.veg)
cor.test(veg.b$brash.dec,veg.b$brash.veg)
# Decay and vegetation are highly correlated

# Remove decay
veg.b$brash.dec <- NULL

# Remove moisture as it's a constant
veg.b$moisture <- NULL
names(veg.b)

# Vegetation data PCA
v.pca <- prcomp(veg.b[c(7,10:16)], scale = TRUE)
fviz_eig(v.pca)
# Eigenvalues
eig.val_v <- get_eigenvalue(v.pca)
# Results for Variables
v.var <- get_pca_var(v.pca)
v.var$coord          # Coordinates
v.var$contrib        # Contributions to the PCs
v.var$cos2           # Quality of representation 
# Results for individuals
v.ind <- get_pca_ind(v.pca)
v.ind$coord          # Coordinates
v.ind$contrib        # Contributions to the PCs
v.ind$cos2           # Quality of representation

veg <- v.ind$coord[,c(1:4)]


# Tree data PCA
names(veg.b)
t.pca <- prcomp(veg.b[c(3:6,8,17)], scale = TRUE)
fviz_eig(t.pca)
# Eigenvalues
eig.val_t <- get_eigenvalue(t.pca)
# Results for Variables
t.var <- get_pca_var(t.pca)
t.var$coord
t.var$contrib
t.var$cos2
# Results for individuals
t.ind <- get_pca_ind(t.pca)

trees <- t.ind$coord[,c(1:4)]


# Build data frame
names(veg.b)
veg.b2 <- data.frame(veg.b[c(1:2,7:9,16)], veg, trees)
names(veg.b2)
names(veg.b2)[7:14] <- c("vPC.1", "vPC.2", "vPC.3", "vPC.4", "tPC.1", "tPC.2", "tPC.3", "tPC.4")
names(veg.b2)



# Load small mammal vegetation surveyor
veg.surv <- read.csv("Small_mammal_veg_surveyors.csv")

# Merge back in vegetation surveyor
names(veg.surv) <- c("location", "treatment", "veg.surveyor")
veg.c <- merge(veg.b2, veg.surv, by = c("location", "treatment"))
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
veg_cap.hist$Surveyor <- paste(veg_cap.hist$veg.surveyor, veg_cap.hist$mam.surveyor, sep="_")
veg_cap.hist$veg.surveyor <- NULL
veg_cap.hist$mam.surveyor <- NULL

# Tidy duplicated surveyor names
veg_cap.hist$Surveyor <- as.character(veg_cap.hist$Surveyor)
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "IS_IS"] <- "IS"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "AMcC_AMcC"] <- "AMc"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "MOC_MOC"] <- "MOC"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "JK_JK"] <- "JK"
veg_cap.hist$Surveyor[veg_cap.hist$Surveyor == "TK_TK"] <- "TK"

#############################################################################
## Plot the Principle Components ##

# Vegetation PCA
# PCs 1 and 2
veg_cap.hist %>%
  group_by(treatment) %>%
  summarise(v1.m = mean(vPC.1), v1.s = sd(vPC.1), v2.m = mean(vPC.2), v2.s = sd(vPC.2)) %>%
  ggplot(aes(x = v1.m, y = v2.m, color = treatment)) + geom_point()  +
  geom_errorbar(aes(ymin = v2.m - v2.s, ymax = v2.m + v2.s), width = 0.2) +
  geom_errorbar(aes(xmin = v1.m - v1.s, xmax = v1.m + v1.s), width = 0.2) +
  xlim(-3,3) + ylim(-3,3) +
  xlab("Vegetation PC1") +  ylab("Vegetation PC2") +
  theme_bw(base_family = "serif") %+replace%
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),  
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16,face="bold"),
    legend.key=element_rect(colour=NA, fill =NA),
    legend.text = element_text(colour="black", size = 14),
    legend.title = element_text(colour="black", size = 14),
    legend.title.align = 0,
    legend.position = c(.85, .15),
    panel.border = element_rect(fill = NA, colour = "black", size=0),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = NA)) +
    scale_colour_manual(name="Treatment",
                      breaks=c("C", "O", "Y"),
                      labels=c("Moorland", "Late pre-thicket\nforest", "Early pre-thicket\nforest"),
                      values = c("black", "blue","orange"))


# Vegetation PCA
# PCs 3 and 4
veg_cap.hist %>%
  group_by(treatment) %>%
  summarise(v3.m = mean(vPC.3), v3.s = sd(vPC.3), v4.m = mean(vPC.4), v4.s = sd(vPC.4)) %>%
  ggplot(aes(x = v3.m, y = v4.m, color = treatment)) + geom_point()  +
  geom_errorbar(aes(ymin = v4.m - v4.s, ymax = v4.m + v4.s), width = 0.2) +
  geom_errorbar(aes(xmin = v3.m - v3.s, xmax = v3.m + v3.s), width = 0.2) +
  xlim(-2.5, 2.5) + ylim(-2,2) +
  xlab("Vegetation PC3") +  ylab("Vegetation PC4") +
  theme_bw(base_family = "serif") %+replace%
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),  
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16,face="bold"),
    legend.key=element_rect(colour=NA, fill =NA),
    legend.text = element_text(colour="black", size = 14),
    legend.title = element_text(colour="black", size = 14),
    legend.title.align = 0,
    legend.position = c(.85, .15),
    panel.border = element_rect(fill = NA, colour = "black", size=0),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = NA)) +
  scale_colour_manual(name="Treatment",
                      breaks=c("C", "O", "Y"),
                      labels=c("Moorland", "Late pre-thicket\nforest", "Early pre-thicket\nforest"),
                      values = c("black", "blue","orange"))






# Tree PCA
# PCs 1 and 2
veg_cap.hist %>%
  group_by(treatment) %>%
  summarise(t1.m = mean(tPC.1), t1.s = sd(tPC.1), t2.m = mean(tPC.2), t2.s = sd(tPC.2)) %>%
  ggplot(aes(x = t1.m, y = t2.m, color = treatment)) + geom_point()  +
  geom_errorbar(aes(ymin = t2.m - t2.s, ymax = t2.m + t2.s), width = 0.2) +
  geom_errorbar(aes(xmin = t1.m - t1.s, xmax = t1.m + t1.s), width = 0.2) +
  xlim(-4, 4) + ylim(-2.5,2.5) +
  xlab("Tree PC1") +  ylab("Tree PC2") +
  theme_bw(base_family = "serif") %+replace%
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),  
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16,face="bold"),
    legend.key=element_rect(colour=NA, fill =NA),
    legend.text = element_text(colour="black", size = 14),
    legend.title = element_text(colour="black", size = 14),
    legend.title.align = 0,
    legend.position = c(.85, .15),
    panel.border = element_rect(fill = NA, colour = "black", size=0),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = NA)) +
  scale_colour_manual(name="Treatment",
                      breaks=c("C", "O", "Y"),
                      labels=c("Moorland", "Late pre-thicket\nforest", "Early pre-thicket\nforest"),
                      values = c("black", "blue","orange"))


# Tree PCA
# PCs 3 and 4
veg_cap.hist %>%
  group_by(treatment) %>%
  summarise(t3.m = mean(tPC.3), t3.s = sd(tPC.3), t4.m = mean(tPC.4), t4.s = sd(tPC.4)) %>%
  ggplot(aes(x = t3.m, y = t4.m, color = treatment)) + geom_point()  +
  geom_errorbar(aes(ymin = t4.m - t4.s, ymax = t4.m + t4.s), width = 0.2) +
  geom_errorbar(aes(xmin = t3.m - t3.s, xmax = t3.m + t3.s), width = 0.2) +
  xlim(-1, 1) + ylim(-1,1) +
  xlab("Tree PC3") +  ylab("Tree PC4") +
  theme_bw(base_family = "serif") %+replace%
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),  
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16,face="bold"),
    legend.key=element_rect(colour=NA, fill =NA),
    legend.text = element_text(colour="black", size = 14),
    legend.title = element_text(colour="black", size = 14),
    legend.title.align = 0,
    legend.position = c(.85, .15),
    panel.border = element_rect(fill = NA, colour = "black", size=0),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = NA)) +
  scale_colour_manual(name="Treatment",
                      breaks=c("C", "O", "Y"),
                      labels=c("Moorland", "Late pre-thicket\nforest", "Early pre-thicket\nforest"),
                      values = c("black", "blue","orange"))


#######################################


# Make new dataframes with each species for all sites
bvcaptures <- veg_cap.hist[veg_cap.hist$species == "BV",]
wmcaptures <- veg_cap.hist[veg_cap.hist$species == "FM",]
gscaptures <- veg_cap.hist[veg_cap.hist$species == "GS",]


###################################
#### Treatment as fixed effect ####

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
# Main model is a better fit

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


