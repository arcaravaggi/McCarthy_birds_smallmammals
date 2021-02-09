#### McCarthy et al. 2021 ####
#### Diversity Analysis ####

# Clear R's memory
rm(list = ls())

# Load packages
library(vegan)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in data
diversity_data_raw <- read.csv("Diversity_data.csv")
names(diversity_data_raw)

# Make new column for location + treatment + count, e.g. BHE_C_1
diversity_data_raw$st_tr_co <- paste(diversity_data_raw[,2], diversity_data_raw[,4], diversity_data_raw[,6],sep="_")

# Make new dataframe of totals for each site (location + treatment + count period)
templs1 <- split(diversity_data_raw, diversity_data_raw$st_tr_co)
templs2 <- lapply(templs1, function(x){x <- x[,7:23]})
tempdf1 <- do.call(rbind, lapply(templs1, function(x)(x <- x[1,])))
tempdf2 <- do.call(rbind, lapply(templs2, colSums))
tempdf1[,7:23] <- tempdf2

# Remove unwanted dataframes and lists
rm(tempdf2, templs1, templs2)

# Rename dataframe
diversity_totals <- tempdf1

# Remove old name
rm(tempdf1)

# Delete sites where nothing was recorded
diversity_totals <- diversity_totals[!(diversity_totals$CHAFF == "0" &
                                         diversity_totals$MEAPI == "0" &
                                         diversity_totals$REDGR == "0" &
                                         diversity_totals$SKYLA == "0" &
                                         diversity_totals$WILWA == "0" &
                                         diversity_totals$WREN == "0" &
                                         diversity_totals$BLABI == "0" &
                                         diversity_totals$BLUTI == "0" &
                                         diversity_totals$COATI == "0" &
                                         diversity_totals$DUNNO == "0" &
                                         diversity_totals$LESRE == "0" &
                                         diversity_totals$PHEAS == "0" &
                                         diversity_totals$ROBIN == "0" &
                                         diversity_totals$WHITE == "0" &
                                         diversity_totals$BLACA == "0" &
                                         diversity_totals$GOLDC == "0" &
                                         diversity_totals$REEBU == "0"),]

# Seperate out the count periods
div_count_1 <- diversity_totals[diversity_totals$Count == "1",]
div_count_2 <- diversity_totals[diversity_totals$Count == "2",]


## Count period 1
# Delete non-species columns
names(div_count_1)
diversity_1 <- div_count_1[,-c(1:6, 24)]
names(diversity_1)

# Calculate diversity index
diversity_index_1 <- diversity(diversity_1,index = "simpson")
div_count_1[,(ncol(div_count_1) + 1)] <- diversity_index_1


## Count period 2
# Delete non-species columns
names(div_count_2)
diversity_2 <- div_count_2[,-c(1:6, 24)]
names(diversity_2)

# Calculate diversity index
diversity_index_2 <- diversity(diversity_2,index = "simpson")
div_count_2[,(ncol(div_count_2) + 1)] <- diversity_index_2




#### Descriptives ####

# Obtain mean +/- SE of diversity indices for each treatment during each count period
# Load package for calculating standard error
library(plotrix)

## Count period 1 means/SE per treatment ##
# Moorland
div_1_c_mean <- div_count_1[div_count_1$Treatment == "C",]; mean(div_1_c_mean$V25)
# 0.2286798
std.error(div_1_c_mean$V25)
# 0.06636877

# Early pre-thicket
div_1_y_mean <- div_count_1[div_count_1$Treatment == "Y",]; mean(div_1_y_mean$V25)
# 0.4286771
std.error(div_1_y_mean$V25)
# 0.1209766

# Late pre-thicket
div_1_o_mean <- div_count_1[div_count_1$Treatment == "O",]; mean(div_1_o_mean$V25)
# 0.6922271
std.error(div_1_o_mean$V25)
# 0.03478615




## Count period 2 means/SE per treatment ##
# Moorland
div_2_c_mean <- div_count_2[div_count_2$Treatment == "C",]; mean(div_2_c_mean$V25)
# 0.07844592
std.error(div_2_c_mean$V25)
# 0.0298865

# Early pre-thicket
div_2_y_mean <- div_count_2[div_count_2$Treatment == "Y",]; mean(div_2_y_mean$V25)
# 0.3730296
std.error(div_2_y_mean$V25)
# 0.09368495

# Late pre-thicket
div_2_o_mean <- div_count_2[div_count_2$Treatment == "O",]; mean(div_2_o_mean$V25)
# 0.6158634
std.error(div_2_o_mean$V25)
# 0.05103559



####################################################################################

# Read in small mammal data
sm_data <- read.csv("Small_mammal_diversity.csv")

# Delete sites where nothing was recorded
sm_data_1 <- sm_data[!(sm_data$BV == "0" &
                           sm_data$FM == "0" &
                           sm_data$GS == "0"),]

# Delete non-species columns
names(sm_data_1)
sm_data_2 <- sm_data_1[,-c(1:2)]
names(sm_data_2)


# Calculate diversity index
sm_div <- diversity(sm_data_2,index = "simpson")
sm_data_1[,(ncol(sm_data_1) + 1)] <- sm_div



#### Descriptives ####

# Moorland
sm_data_c_mean <- sm_data_1[sm_data_1$Treatment == "C",]; mean(sm_data_c_mean$V6)
# 0.4184362
std.error(sm_data_c_mean$V6)
# 0.1361153

# Early pre-thicket
sm_data_y_mean <- sm_data_1[sm_data_1$Treatment == "Y",]; mean(sm_data_y_mean$V6)
# 0.406773
std.error(sm_data_y_mean$V6)
# 0.08798551

# Late pre-thicket
sm_data_o_mean <- sm_data_1[sm_data_1$Treatment == "O",]; mean(sm_data_o_mean$V6)
# 0.203393
std.error(sm_data_o_mean$V6)
# 0.08443713



