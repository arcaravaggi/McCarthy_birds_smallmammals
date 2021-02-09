#### McCarthy et al. 2021 ####
#### Bird Models ####

# Clear R's memory
rm(list = ls())

# Load required packages
library(ggplot2)
library(lme4)
library(dplyr)
library(FSA)
library(emmeans)
library(tidyr)
library(PMCMRplus)
library(car)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load in the 3 datafiles
all_spp <- read.csv("all_spp_complete.csv")
open_spp <- read.csv("open_spp_complete.csv")
scrub_spp <- read.csv("scrub_spp_complete.csv")

# Read in and merge time
time_data <- read.csv("Point_count_times.csv")
names(time_data)

# Create a STPC (site,treatment,point,count) column
names(time_data)
time_data$STPC <- paste(time_data[,1], time_data[,2], time_data[,6], time_data[,3], sep="_")
names(time_data)
time_data <- time_data[,-c(1:3,6)]
# Delete duplicated STPCs
time_data <- time_data[!duplicated(time_data$STPC),]

# Merge relevant times with the three dataframes based on STPC
all_spp <- merge(all_spp, time_data, by = "STPC", all = FALSE)
open_spp <- merge(open_spp, time_data, by = "STPC", all = FALSE)
scrub_spp <- merge(scrub_spp, time_data, by = "STPC", all = FALSE)

# Graph start times of counts for each treatment
# This will ensure treatments were surveyed at similar times
all_enc_graph <- ggplot(data = all_spp, aes(x= Start.time, y = Treatment)) +
  stat_summary(geom = "point")
all_enc_graph

open_enc_graph <- ggplot(data = open_spp, aes(x= Start.time, y = Treatment)) +
  stat_summary(geom = "point")
open_enc_graph

scrub_enc_graph <- ggplot(data = scrub_spp, aes(x= Start.time, y = Treatment)) +
  stat_summary(geom = "point")
scrub_enc_graph
# No difference in time of surveys between treatments




#### Graphing Bird Densities ####
#### All Species Group ####
# Exploring the data by graphing it

# Graph density for each treatment
ggplot(data = all_spp, aes(x= Treatment, y = Estimate, fill = Treatment)) +
  geom_boxplot() + stat_summary(fun.y = mean)

# All treatments together
hist(all_spp$Estimate)

# Moorland
all_spp_con <- all_spp[all_spp$Treatment == "C",]
hist(all_spp_con$Estimate)

# Early pre-thicket
all_spp_young <- all_spp[all_spp$Treatment == "Y",]
hist(all_spp_young$Estimate)

# Late pre-thicket
all_spp_old <- all_spp[all_spp$Treatment == "O",]
hist(all_spp_old$Estimate)




#### Open Species Group ####

# Do histograms of density and biomass for each treatment
ggplot(data = open_spp, aes(x= Treatment, y = Estimate, fill = Treatment)) +
  geom_boxplot() + stat_summary(fun.y = mean)

# Moorland
open_spp_con <- open_spp[open_spp$Treatment == "C",]
hist(open_spp_con$Estimate)

# Early pre-thicket
open_spp_young <- open_spp[open_spp$Treatment == "Y",]
hist(open_spp_young$Estimate)




#### Scrub Species Group ####

# Do histograms of density and biomass for each treatment
ggplot(data = scrub_spp, aes(x= Treatment, y = Estimate, fill = Treatment)) +
  geom_boxplot() + stat_summary(fun.y = mean)

# Late pre-thicket
scrub_spp_old <- scrub_spp[scrub_spp$Treatment == "O",]
hist(scrub_spp_old$Estimate)

# Early pre-thicket
scrub_spp_young <- scrub_spp[scrub_spp$Treatment == "Y",]
hist(scrub_spp_young$Estimate)




###################################
#### Statistical Tests- Models ####
###################################

#### All Species Group Stats ####
## Genereal Linear Mixed Effects Models ##

# Test for normality
hist(all_spp$Estimate)
shapiro.test(all_spp$Estimate)
# Density not normally distributed

# Outlier present in the data
# Outlier produced by two birds being v. close to observer
# and the remaining birds far away
# Atypical disribution of birds, therefore remove outlier
all_spp1<-subset(all_spp,Estimate<15)
# Test for normality again
hist(all_spp1$Estimate)
shapiro.test(all_spp1$Estimate)
# Now normally distributed


# All Species Group mixed effects model
all_spp_lmer <- lmer(Estimate ~ Treatment + (1|Site), data = all_spp1)
summary(all_spp_lmer)
Anova(all_spp_lmer)
plot(all_spp_lmer)
qqnorm(resid(all_spp_lmer)); qqline(resid(all_spp_lmer))
# QQplot slightly tailed though acceptable

# Conduct a post-hoc test
emmeans(all_spp_lmer, list(pairwise ~ Treatment), adjust = "tukey")
# Significant difference between all treatments


################################################

#### Open Species Group Stats ####
## Genereal Linear Mixed Effects Models ##

# Test for normality
hist(open_spp$Estimate)
shapiro.test(open_spp$Estimate)
# Density is not normally distributed

# Delete outlier as above
open_spp1<-subset(open_spp,Estimate<15)

# Test for normality
hist(open_spp1$Estimate)
shapiro.test(open_spp1$Estimate)
# Now normally distributed (density)

## Genereal Linear Mixed Effects Model ##
# Open Species Group mixed effects model
open_spp_lmer <- lmer(Estimate ~ Treatment + (1|Site), data = open_spp1)
summary(open_spp_lmer)
anova(open_spp_lmer)
plot(open_spp_lmer)
qqnorm(resid(open_spp_lmer)); qqline(resid(open_spp_lmer))
# Residuals and QQ plot look okay

# Conduct a post-hoc test
emmeans(open_spp_lmer, list(pairwise ~ Treatment), adjust = "tukey")
# Significant difference between the two treatments


###############################################

#### Scrub Species Group Stats ####
## Genereal Linear Mixed Effects Models ##

# General Linear Mixed Effects Model
scrub_spp_lmer <- lmer(Estimate ~ Treatment + (1|Site), data = scrub_spp)
summary(scrub_spp_lmer)
plot(scrub_spp_lmer)
qqnorm(resid(scrub_spp_lmer)); qqline(resid(scrub_spp_lmer))
# qqplot and residuals are okay

# Conduct a post-hoc test
emmeans(scrub_spp_lmer, list(pairwise ~ Treatment), adjust = "tukey")
# Significant difference between the two treatments


############################################## ##############################################
############################################## ##############################################

#### DESCRIPTIVES ####

# Obtain mean and standard error (SE) of density within each treatment
# Write a function that calculates SE
ste <- function(x) sd(x)/sqrt(length(x))

# ALL SPECIES #
all_spp1_C <- all_spp1[all_spp1$Treatment == "C",]
mean(all_spp1_C$Estimate)
# 8.228
ste(all_spp1_C$Estimate)
#0.427

all_spp1_Y <- all_spp1[all_spp1$Treatment == "Y",]
mean(all_spp1_Y$Estimate)
#3.412
ste(all_spp1_Y$Estimate)
#0.416

all_spp1_O <- all_spp1[all_spp1$Treatment == "O",]
mean(all_spp1_O$Estimate)
# 6.492
ste(all_spp1_O$Estimate)
# 0.413


# OPEN SPECIES #
open_spp1_C <- open_spp1[open_spp1$Treatment == "C",]
mean(open_spp1_C$Estimate)
# 8.151
ste(open_spp1_C$Estimate)
# 0.450

open_spp1_Y <- open_spp1[open_spp1$Treatment == "Y",]
mean(open_spp1_Y$Estimate)
# 2.740
ste(open_spp1_Y$Estimate)
# 0.387


# SCRUB SPECIES #
scrub_spp_Y <- scrub_spp[scrub_spp$Treatment == "Y",]
mean(scrub_spp_Y$Estimate)
# 0.884
ste(scrub_spp_Y$Estimate)
# 0.204

scrub_spp_O <- scrub_spp[scrub_spp$Treatment == "O",]
mean(scrub_spp_O$Estimate)
# 5.226
ste(scrub_spp_O$Estimate)
# 0.453

##########################
#### GRAPHING RESULTS ####
##########################

# Set common theme for graphs
theme_ac1 <- function(base_family = "serif", base_size_a = 18, base_size_t = 18){
  theme_bw(base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      plot.title = element_text(size=10),
      panel.grid = element_blank(), 
      axis.text = element_text(size = base_size_a),
      axis.title = element_text(size=base_size_t,face="bold"),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.border = element_rect(fill = NA, colour = "black", size=0),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}

# Make a palette with 3 colours of grey for graphs
greypalette_3 <- c("gray99", "gray65","grey35")

# Get rid of unecessary columns
names(all_spp1)
all_spp_1 <- all_spp1[,-c(1:3,7,8,10,11)]
names(all_spp_1)

# Get rid of unecessary columns
names(open_spp1)
open_spp_1 <- open_spp1[,-c(1:3,7,8,10,11)]
names(open_spp_1)

# Get rid of unecessary columns
names(scrub_spp)
scrub_spp_1 <- scrub_spp[,-c(1:3,7,8,10,11:13)]
names(scrub_spp_1)

# Combine the three dataframes above into one dataframe
graph_data <- rbind(all_spp_1, open_spp_1, scrub_spp_1)
names(scrub_spp_1)
# In order to make bars the same width,
# add in rows for open species in late pre-thicket
# and scrub species in moorland
# These are in a separate csv
# Read in and merge
graph_data_blank <- read.csv("graph_data_blanks.csv")
graph_data <- rbind(graph_data, graph_data_blank)

# Now change the bird species group labels to full titles
graph_data$Analysis <- gsub("Open", "Open country bird species", graph_data$Analysis)
graph_data$Analysis <- gsub("Scrub", "Scrub bird species", graph_data$Analysis)
graph_data$Analysis <- gsub("All", "All bird species", graph_data$Analysis)


## OVERALL BIRD DENSITY GRAPH ##
#overview_density_graph <- 
graph_data %>%
ggplot(aes(x= Treatment, y = Estimate, fill = Analysis)) +
 stat_summary(geom = "bar", fun = mean, position = "dodge", width = 0.8, colour = "black") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.8), width = .4)  +
  scale_fill_manual(values = greypalette_3) +
  scale_x_discrete(breaks=unique(graph_data$Treatment),
                   limits = c("C", "Y", "O"),
                   labels = c("Moorland", "Late\npre-thicket\nforest", "Early\npre-thicket\nforest")) + #the order of labels is different but it comes out correct
  ylab("Bird density (birds/ha)") +
  xlab("") +
  theme_ac1(base_size_a = 14, base_size_t = 14) +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14)) +
  guides(fill=guide_legend(title="Group"))

# Save graph at 300dpi resolution
#ggsave(
#  width = 8, height = 5,
#  "Bird density_graph.tiff",
#  dpi = 300
#)

