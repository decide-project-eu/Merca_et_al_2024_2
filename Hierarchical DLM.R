#### Hierarchical DLM - with 2 harmonic waves for country level

#libraries
library(dplyr)

#get data
setwd("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Datasets")
load("df_salmon_20.RData")
df <- df_salmon_20


### Create Learning and Test sets ----

source("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Functions for Salmon study - Hierarchical.R")

#get Learning and Test sets
#use the first ~3/4 (that will be 75%) of time (dates) for Learning.set
N <- round(3*(length(unique(df[order(df[,"date"]), "date"]))/4))
sets <- get.learning.test.sets(df, relevant.names="log0.mortality.rel.20", hierarchical=TRUE, N=N)
Learning.set.i <- sets[["Learning.set"]]
Test.set.i <- sets[["Test.set"]]

### Re-arrange datasets to work with this approach - the datasets have to be each column one farm, and each row one date
#Learning.set
col.needed <- Learning.set.i[,c(1,2,3,29)] 
col.needed <- col.needed %>% # - order by region
  arrange(local.authority)
col.needed[,"site_region"] <- paste(col.needed$site, col.needed$local.authority, sep="_")
col.needed <- col.needed[,c(1,4,5)]

# get re-arranged Learning.set
Learning.set <- as.data.frame(matrix(NA, length(unique(col.needed$date)), length(unique(col.needed$site_region))+1))
colnames(Learning.set) <- c("date", unique(col.needed$site_region))
Learning.set$date <- unique(col.needed$date)
Learning.set <- Learning.set %>% # - order by date
  arrange(date)
for (i in 1:dim(col.needed)[1]){ 
  Learning.set[which(Learning.set$date == col.needed[i,"date"]), col.needed[i,"site_region"]] <- col.needed[i,"log0.mortality.rel.20"]
}

#Test.set
col.needed <- Test.set.i[,c(1,2,3,29)] 
col.needed <- col.needed %>% # - order by region
  arrange(local.authority)
col.needed[,"site_region"] <- paste(col.needed$site, col.needed$local.authority, sep="_")
col.needed <- col.needed[,c(1,4,5)]

# get re-arranged Test.set
Test.set <- as.data.frame(matrix(NA, length(unique(col.needed$date)), length(unique(col.needed$site_region))+1))
colnames(Test.set) <- c("date", unique(col.needed$site_region))
Test.set$date <- unique(col.needed$date)
Test.set <- Test.set %>% # - order by date
  arrange(date)
for (i in 1:dim(col.needed)[1]){ 
  Test.set[which(Test.set$date == col.needed[i,"date"]), col.needed[i,"site_region"]] <- col.needed[i,"log0.mortality.rel.20"]
}

### Create a dataset with information about months since start for each production cycle - the dataset has to be each column one farm, and each row one date
#Learning.set
months.learning.set <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
months.learning.set[,1] <- Learning.set$date
colnames(months.learning.set) <- colnames(Learning.set)

for (farm in unique(Learning.set.i$site)){ 
  df.farm <- subset(Learning.set.i, Learning.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      months.learning.set[which(months.learning.set$date==date), 
                          paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      months.learning.set[which(months.learning.set$date==date), 
                          paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "months.since.start"]
    }
  }
}

#Test.set
months.test.set <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
months.test.set[,1] <- Test.set$date
colnames(months.test.set) <- colnames(Test.set)

for (farm in unique(Test.set.i$site)){ 
  df.farm <- subset(Test.set.i, Test.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      months.test.set[which(months.test.set$date==date), 
                      paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      months.test.set[which(months.test.set$date==date), 
                      paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "months.since.start"]
    }
  }
}

### Create a dataset with information about the number of the production cycle (ID) - the dataset has to be each column one farm, and each row one date
#Learning.set
nseq.learning.set <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
nseq.learning.set[,1] <- Learning.set$date
colnames(nseq.learning.set) <- colnames(Learning.set)

for (farm in unique(Learning.set.i$site)){
  df.farm <- subset(Learning.set.i, Learning.set.i$site == farm)
  
  for(date in unique(df.farm$date)){
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      nseq.learning.set[which(nseq.learning.set$date==date), 
                        paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      nseq.learning.set[which(nseq.learning.set$date==date), 
                        paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "nseq"]
    }
  }
}

#Test.set
nseq.test.set <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
nseq.test.set[,1] <- Test.set$date
colnames(nseq.test.set) <- colnames(Test.set)

for (farm in unique(Test.set.i$site)){ 
  df.farm <- subset(Test.set.i, Test.set.i$site == farm)
  
  for(date in unique(df.farm$date)){ 
    if(is.na(df.farm[which(df.farm$date==date), "months.since.start"])){ #if fallow period (NaN)
      nseq.test.set[which(nseq.test.set$date==date), 
                    paste(farm, unique(df.farm$local.authority), sep="_")] <- NaN
    }else{
      nseq.test.set[which(nseq.test.set$date==date), 
                    paste(farm, unique(df.farm$local.authority), sep="_")] <- df.farm[which(df.farm$date==date), "nseq"]
    }
  }
}


### See if all dates have information at least for one farm - if not remove them
dataset <- Learning.set

dataset$sum <- rowMeans(dataset[,2:ncol(dataset)], na.rm=T)
row.remove <- which(is.na(dataset$sum)) #remove this row from the Learning sets

Learning.set <- Learning.set[-row.remove,]
months.learning.set <- months.learning.set[-row.remove,]
nseq.learning.set <- nseq.learning.set[-row.remove,]


### Get the optimum values of Country, Region and Farms variances and discount factors  ----

# Make sure you avoid scientific notation in your output - it will make things easier!
options(scipen=999)

# Standardize the learning set data (force the data to be standard normal distributed)
all.values <- unlist(Learning.set[,c(2:ncol(Learning.set))])
Mean <- mean(na.omit(all.values))
SD <- sd(na.omit(all.values))
Standarized.factors <- c(Mean=Mean, SD=SD)

Learning.set.stand <- data.frame(matrix(NA, nrow = dim(Learning.set)[1], ncol = dim(Learning.set)[2]))
Learning.set.stand[,1] <- Learning.set$date
colnames(Learning.set.stand) <- colnames(Learning.set)

Learning.set.stand[1:dim(Learning.set.stand)[1], 2:dim(Learning.set.stand)[2]] <- (unlist(Learning.set[,c(2:ncol(Learning.set))]) - Mean)/SD

# Get mu0 and C0
regions.names <- c("Argyll and Bute", "Eilean Siar", "Highland", "Orkney Islands", "Shetland Islands")
mu <- get.m0(D=Learning.set.stand, D.months=months.learning.set, expected.start.time=0, H.w.country=2, regions=regions.names)
mu0 <- mu$mu0
Spline.list  <- mu$Spline.list

C <- get.C0(D=Learning.set.stand, D.months=months.learning.set, expected.start.time=0, expected.finish.time=6, H.w.country=2, regions=regions.names)
C0 <- C$C0

# Run the optimum function with DLM and mu0 and C0 to get the best estimates for the 3 V's and 3 deltas (Country, Region and Farm), priorVs and priorDeltas are values from 0 to 1 and are of minor importance (arbitrarily chosen)
## This can take some time!
Start.time <- Sys.time()

est = estimateDiscountModel(priorVs=c(0.01, 0.01, 0.5),
                            priorDeltas=c(0.99, 0.99, 0.99),
                            D=Learning.set.stand, D.months=months.learning.set, D.nseq=nseq.learning.set,
                            mu0, C0, H.w.country=2, Spline.list, regions=regions.names)

print(Sys.time()-Start.time)

# - Save the 3 V's and 3 deltas
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM"
saveRDS(object = est, file = paste(model.dir, '/estimates.Vs.deltas.RDS', sep=''))
Vs.deltas <- readRDS(file = paste(model.dir, '/estimates.Vs.deltas.RDS', sep=''))

# - Get the transformed back values
Vs.deltas.t <- transformResults(Vs.deltas, vars=3)


### Apply Univariate hierarchical DLM ----

# Get the functions needed
source("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Functions for Salmon study - Hierarchical.R")

# Standardize the test set data (force the data to be standard normal distributed)
Test.set.stand <- data.frame(matrix(NA, nrow = dim(Test.set)[1], ncol = dim(Test.set)[2]))
Test.set.stand[,1] <- Test.set$date
colnames(Test.set.stand) <- colnames(Test.set)

Test.set.stand[1:dim(Test.set.stand)[1], 2:dim(Test.set.stand)[2]] <- (unlist(Test.set[,c(2:ncol(Test.set))]) - Standarized.factors["Mean"])/Standarized.factors["SD"]

# Get the optimum V's and deltas
countryVar = Vs.deltas.t[1]
regionVar = Vs.deltas.t[2]
farmVar = Vs.deltas.t[3]
deltas = Vs.deltas.t[4:6]

# Run the DLM with the optimum V's and deltas
out <- runDiscountDLM(D=Test.set.stand, D.months=months.test.set, mu0, C0, countryVar, regionVar, farmVar, deltas, H.w.country=2, Spline.list, regions=regions.names)

# Run the Smoother on the DLM outcomes
smot <- runSmoother(out)

# - Save the outcomes and the Smoothed values
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM"
saveRDS(object = out, file = paste(model.dir, '/outDLM.RDS', sep=''))
saveRDS(object = smot, file = paste(model.dir, '/smotDLM.RDS', sep=''))

# Get the results for out and smot
results <- extract.res(out, smot, D=Test.set.stand, D.nseq=nseq.test.set, H.w.country=2, regions=regions.names)

# - Save the results
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM"
saveRDS(object = results, file = paste(model.dir, '/resultsDLM.RDS', sep=''))

# - Get RMSE (applied to the Test set!)
#### Get the results
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM"
results <- readRDS(file = paste(model.dir, '/resultsDLM.RDS', sep=''))

#we can get the RMSE from the DLM on test set
et.results <- results[["et"]]
et.all <- na.omit(unlist(et.results[,c(2:ncol(et.results))]))
RMSE <- round(sqrt(mean(et.all^2)),4) 


### Plot Univariate hierarchical DLM results ----
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM"
results <- readRDS(file = paste(model.dir, '/resultsDLM.RDS', sep=''))
out <- readRDS(file = paste(model.dir, '/outDLM.RDS', sep=''))

# D is now the standardize test set
D <- Test.set.stand
D.nseq <- nseq.test.set

# - plot to pdf - using ggplot
setwd("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Hierarchical DLM/plots")
pdf("Hiearchical_mortality.pdf")

## Country ----
# - Plot with log transformed and standardized values
min <- min(c(results$mt$Country-1.96*sqrt(results$Ct$Country),
             results$mts$Country-1.96*sqrt(results$Cts$Country)), na.rm=T)
max <- max(c(results$mt$Country+1.96*sqrt(results$Ct$Country),
             results$mts$Country+1.96*sqrt(results$Cts$Country)), na.rm=T)

#mt - darkgreen
low.limit.mt <- (results$mt$Country)-1.96*sqrt(results$Ct$Country)
high.limit.mt <- (results$mt$Country)+1.96*sqrt(results$Ct$Country)

mt.df <- data.frame(date=D[,"date"], mt=results$mt$Country,
                    low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)

plot(mt.df$mt ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
     main="Country", ylim=c(min,max), lwd=2, col='darkgreen',
     cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean

lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance

legend(x = "bottomright",                                # Position
       legend = c("mt", "95% CI"),                       # Legend texts
       col = c("darkgreen", "darkgreen"),                # Line colors
       lwd = 1,                                          # Line thickness
       lty = c(1,2),
       cex=0.7)

#mts (smoother) - blue
low.limit.mts <- results$mts$Country-1.96*sqrt(results$Cts$Country)
high.limit.mts <- results$mts$Country+1.96*sqrt(results$Cts$Country)

mts.df <- data.frame(date=D[,"date"], mts=results$mts$Country,
                     low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)

plot(mts.df$mts ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
     main="Country", ylim=c(min,max), lwd=2, col='blue',
     cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean

lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance

legend(x = "bottomright",                                # Position
       legend = c("mts", "95% CI"),                      # Legend texts
       col = c("blue", "blue"),                          # Line colors
       lwd = 1,                                          # Line thickness
       lty = c(1,2),
       cex=0.7)

# - remove standardization from mortality
mt_no_s <- (results$mt$Country * Standarized.factors["SD"]) + Standarized.factors["Mean"]
mts_no_s <- (results$mts$Country * Standarized.factors["SD"]) + Standarized.factors["Mean"]

## - For the CI's we have to remove the standardization on the all CI calculation, 
## - not individually on Ct, Qt and Cts
low.limit.mt_no_s <- (low.limit.mt * Standarized.factors["SD"]) + Standarized.factors["Mean"]
high.limit.mt_no_s <- (high.limit.mt * Standarized.factors["SD"]) + Standarized.factors["Mean"]
low.limit.mts_no_s <- (low.limit.mts * Standarized.factors["SD"]) + Standarized.factors["Mean"]
high.limit.mts_no_s <- (high.limit.mts * Standarized.factors["SD"]) + Standarized.factors["Mean"]

# - remove log transformation from mortality (solve in https://www.wolframalpha.com/)
mt <- exp(1)^mt_no_s - 1/20000
mts <- exp(1)^mts_no_s - 1/20000

## - For the CI's we have again to remove log transformation on the all CI calculation, 
## - not individually on Ct, Qt and Cts
low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000

# - plot with the real values
min <- min(c(low.limit.mt,low.limit.mt), na.rm=T)
max <- max(c(high.limit.mt,high.limit.mts), na.rm=T)

#mt - darkgreen
mt.df <- data.frame(date=D[,"date"], mt=mt,
                    low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)

plot(mt.df$mt ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
     main="Country", ylim=c(min,max), lwd=2, col='darkgreen',
     cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean

lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance

legend(x = "topright",                                   # Position
       legend = c("mt", "95% CI"),                       # Legend texts
       col = c("darkgreen", "darkgreen"),                # Line colors
       lwd = 1,                                          # Line thickness
       lty = c(1,2),
       cex=0.7)

#mts - blue
mts.df <- data.frame(date=D[,"date"], mts=mts,
                     low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)

plot(mts.df$mts ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
     main="Country", ylim=c(min,max), lwd=2, col='blue',
     cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean

lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance

legend(x = "topright",                                   # Position
       legend = c("mts", "95% CI"),                      # Legend texts
       col = c("blue", "blue"),                          # Line colors
       lwd = 1,                                          # Line thickness
       lty = c(1,2),
       cex=0.7)


## Region ----

regions <- c("Argyll and Bute", "Eilean Siar", "Highland", "Orkney Islands", "Shetland Islands")

for (region in regions){ 
  
  min <- -1
  max <- 2
  
  #mt - darkgreen
  mt=results$mt[region]
  colnames(mt) <- NULL
  
  low.limit.mt <- unlist((results$mt[region])-1.96*sqrt(results$Ct[region]))
  high.limit.mt <- unlist((results$mt[region])+1.96*sqrt(results$Ct[region]))
  
  mt.df <- data.frame(date=D[,"date"], mt=mt,
                      low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)
  
  plot(mt.df$mt ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
       main=region, ylim=c(min,max), lwd=2, col='darkgreen',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean
  
  lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
  lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance
  
  legend(x = "bottomright",                                # Position
         legend = c("mt", "95% CI"),                       # Legend texts
         col = c("darkgreen", "darkgreen"),                # Line colors
         lwd = 1,                                          # Line thickness
         lty = c(1,2),
         cex=0.7)
  
  #mts (smoother) - blue
  mts=results$mts[region]
  colnames(mts) <- NULL
  
  low.limit.mts <- unlist((results$mts[region])-1.96*sqrt(results$Cts[region]))
  high.limit.mts <- unlist((results$mts[region])+1.96*sqrt(results$Cts[region]))
  
  mts.df <- data.frame(date=D[,"date"], mts=mts,
                       low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)
  
  plot(mts.df$mts ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
       main=region, ylim=c(min,max), lwd=2, col='blue',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean
  
  lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
  lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance
  
  legend(x = "bottomright",                                # Position
         legend = c("mts", "95% CI"),                      # Legend texts
         col = c("blue", "blue"),                          # Line colors
         lwd = 1,                                          # Line thickness
         lty = c(1,2),
         cex=0.7)
  
  # - remove standardization from mortality
  mt_no_s <- (results$mt[region] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  mts_no_s <- (results$mts[region] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (unlist(low.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mt_no_s <- (unlist(high.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.mts_no_s <- (unlist(low.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mts_no_s <- (unlist(high.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  # - remove log transformation from mortality (solve in https://www.wolframalpha.com/)
  mt <- exp(1)^mt_no_s - 1/20000
  colnames(mt) <- NULL
  mts <- exp(1)^mts_no_s - 1/20000
  colnames(mts) <- NULL
  
  ## - For the CI's we have again to remove log transformation on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
  high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
  low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
  high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000
  
  # - plot with the real values
  min <- min(c(low.limit.mt, low.limit.mts), na.rm=T)
  max <- max(c(high.limit.mt, high.limit.mts), na.rm=T)
  
  #mt - darkgreen
  mt.df <- data.frame(date=D[,"date"], mt=mt,
                      low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)
  
  plot(mt.df$mt ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
       main=region, ylim=c(min,max), lwd=2, col='darkgreen',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean
  
  lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
  lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance
  
  legend(x = "topright",                                   # Position
         legend = c("mt", "95% CI"),                       # Legend texts
         col = c("darkgreen", "darkgreen"),                # Line colors
         lwd = 1,                                          # Line thickness
         lty = c(1,2),
         cex=0.7)
  
  #mts - blue
  mts.df <- data.frame(date=D[,"date"], mts=mts,
                       low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)
  
  plot(mts.df$mts ~ mt.df$date, type="l", xlab="date", ylab="log.mortality", 
       main=region, ylim=c(min,max), lwd=2, col='blue',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #filtered mean
  
  lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
  lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance
  
  legend(x = "topright",                                   # Position
         legend = c("mts", "95% CI"),                      # Legend texts
         col = c("blue", "blue"),                          # Line colors
         lwd = 1,                                          # Line thickness
         lty = c(1,2),
         cex=0.7)
}


## Farms ----

farms <- colnames(D[2:ncol(D)])

for (farm in farms){
  
  # - Plot with log transformed and standardized values
  min <- min(c(unlist(results$ft[farm]-1.96*sqrt(results$Qt[farm])), D[,farm]), na.rm=T)
  max <- max(c(unlist(results$ft[farm]+1.96*sqrt(results$Qt[farm])), D[,farm]), na.rm=T)
  
  # - information about when each production cycle starts and ends
  cycles <- unique(D.nseq[,farm])[!is.na(unique(D.nseq[,farm]))]
  rects <- as.data.frame(matrix(NA, ncol = 2, nrow = length(cycles)))
  colnames(rects) <- c("xstart", "xend")
  loop=0
  for(cycle in cycles){
    loop=loop+1
    dates <- D.nseq[which(D.nseq[,farm]==cycle), "date"]
    rects[loop, "xstart"] <- dates[1]
    rects[loop, "xend"] <- dates[length(dates)]
  }
  rects$xstart <- as.Date(rects$xstart, format = "%Y-%m-%d", origin = "1970-01-01")
  rects$xend <- as.Date(rects$xend, format = "%Y-%m-%d", origin = "1970-01-01")
  
  #mt - darkgreen
  mt=results$mt[farm]
  colnames(mt) <- NULL
  
  low.limit.mt <- unlist((results$mt[farm])-1.96*sqrt(results$Ct[farm]))
  high.limit.mt <- unlist((results$mt[farm])+1.96*sqrt(results$Ct[farm]))
  
  mt.df <- data.frame(date=D[,"date"], obs=D[,farm], mt=mt,
                      low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)
  
  plot(mt.df$mt ~ mt.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='darkgreen',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -10, xright = rects$xend[i], ytop = 10, 
         col = "lightgrey", border = NA)
  }
  
  lines(mt.df$mt ~ mt.df$date, col='darkgreen', lwd=2)               #filtered mean
  lines(mt.df$obs ~ mt.df$date, col='black', lwd=2)                  #observations
  lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
  lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance
  
  legend(x = "bottomright",                                # Position
         legend = c("obs", "mt", "95% CI", "cycle"),       # Legend texts
         col = c("black", "darkgreen", "darkgreen", NA),   # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)

  #ft - red
  ft=results$ft[farm]
  colnames(ft) <- NULL
  
  low.limit.ft <- unlist((results$ft[farm])-1.96*sqrt(results$Qt[farm]))
  high.limit.ft <- unlist((results$ft[farm])+1.96*sqrt(results$Qt[farm]))
  
  ft.df <- data.frame(date=D[,"date"], obs=D[,farm], ft=ft,
                      low.limit.ft=low.limit.ft, high.limit.ft=high.limit.ft)
  
  plot(ft.df$ft ~ ft.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='red',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)           #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -10, xright = rects$xend[i], ytop = 10, 
         col = "lightgrey", border = NA)
  }
  
  lines(ft.df$ft ~ ft.df$date, col='red', lwd=2)               #filtered mean
  lines(ft.df$obs ~ ft.df$date, col='black', lwd=2)            #observations
  lines(low.limit.ft ~ ft.df$date, col='red', lty = "dashed")  #filtered variance
  lines(high.limit.ft ~ ft.df$date, col='red', lty = "dashed") #filtered variance
  
  legend(x = "bottomright",                                # Position
         legend = c("obs", "ft", "95% CI", "cycle"),       # Legend texts
         col = c("black", "red", "red", NA),               # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)
  
  #mts (smoother) - blue
  mts=results$mts[farm]
  colnames(mts) <- NULL
  
  low.limit.mts <- unlist((results$mts[farm])-1.96*sqrt(results$Cts[farm]))
  high.limit.mts <- unlist((results$mts[farm])+1.96*sqrt(results$Cts[farm]))
  colnames(low.limit.mts) <- NULL
  colnames(high.limit.mts) <- NULL
  
  mts.df <- data.frame(date=D[,"date"], obs=D[,farm], mts=mts,
                       low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)
  
  plot(mts.df$mts ~ mts.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='blue',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)            #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -10, xright = rects$xend[i], ytop = 10, 
         col = "lightgrey", border = NA)
  }
  
  lines(mts.df$mts ~ mts.df$date, col='blue', lwd=2)              #filtered mean
  lines(mts.df$obs ~ mts.df$date, col='black', lwd=2)             #observations
  lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
  lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance
  
  legend(x = "bottomright",                                # Position
         legend = c("obs", "mts", "95% CI", "cycle"),      # Legend texts
         col = c("black", "blue", "blue", NA),             # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)
    
  # - remove standardization from mortality
  mt_no_s <- (results$mt[farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  ft_no_s <- (results$ft[farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  mts_no_s <- (results$mts[farm] * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (unlist(low.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mt_no_s <- (unlist(high.limit.mt) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.ft_no_s <- (unlist(low.limit.ft) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.ft_no_s <- (unlist(high.limit.ft) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  low.limit.mts_no_s <- (unlist(low.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  high.limit.mts_no_s <- (unlist(high.limit.mts) * Standarized.factors["SD"]) + Standarized.factors["Mean"]
  
  # - remove log transformation from mortality (solve in https://www.wolframalpha.com/)
  Obs <- exp(1)^(Test.set[,farm]) - 1/20000
  mt <- exp(1)^mt_no_s - 1/20000
  colnames(mt) <- NULL
  ft <- exp(1)^ft_no_s - 1/20000
  colnames(ft) <- NULL
  mts <- exp(1)^mts_no_s - 1/20000
  colnames(mts) <- NULL
  
  ## - For the CI's we have again to remove log transformation on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
  high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
  low.limit.ft <- exp(1)^low.limit.ft_no_s - 1/20000
  high.limit.ft <- exp(1)^high.limit.ft_no_s - 1/20000
  low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
  high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000
  
  # - plot with the real values
  min <- min(c(low.limit.ft, Obs), na.rm=T)
  max <- max(c(high.limit.ft, Obs), na.rm=T)
  
  #mt - darkgreen
  mt.df <- data.frame(date=D[,"date"], obs=Obs, mt=mt,
                      low.limit.mt=low.limit.mt, high.limit.mt=high.limit.mt)
  
  plot(mt.df$mt ~ mt.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='darkgreen',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)       #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -20, xright = rects$xend[i], ytop = 20, 
         col = "lightgrey", border = NA)
  }
  
  lines(mt.df$mt ~ mt.df$date, col='darkgreen', lwd=2)               #filtered mean
  lines(mt.df$obs ~ mt.df$date, col='black', lwd=2)                  #observations
  lines(low.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed")  #filtered variance
  lines(high.limit.mt ~ mt.df$date, col='darkgreen', lty = "dashed") #filtered variance
  
  legend(x = "topright",                                   # Position
         legend = c("obs", "mt", "95% CI", "cycle"),       # Legend texts
         col = c("black", "darkgreen", "darkgreen", NA),   # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)
  
  #ft - red
  ft.df <- data.frame(date=D[,"date"], obs=Obs, ft=ft,
                      low.limit.ft=low.limit.ft, high.limit.ft=high.limit.ft)
  
  plot(ft.df$ft ~ ft.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='red',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)           #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -20, xright = rects$xend[i], ytop = 20, 
         col = "lightgrey", border = NA)
  }
  
  lines(ft.df$ft ~ ft.df$date, col='red', lwd=2)               #filtered mean
  lines(ft.df$obs ~ ft.df$date, col='black', lwd=2)            #observations
  lines(low.limit.ft ~ ft.df$date, col='red', lty = "dashed")  #filtered variance
  lines(high.limit.ft ~ ft.df$date, col='red', lty = "dashed") #filtered variance
  
  legend(x = "topright",                                   # Position
         legend = c("obs", "ft", "95% CI", "cycle"),       # Legend texts
         col = c("black", "red", "red", NA),               # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)
  
  #mts - blue
  mts.df <- data.frame(date=D[,"date"], obs=Obs, mts=mts,
                       low.limit.mts=low.limit.mts, high.limit.mts=high.limit.mts)
  
  plot(mts.df$mts ~ mts.df$date, type="n", xlab="date", ylab="log.mortality", 
       main=farm, ylim=c(min,max), lwd=2, col='blue',
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)             #just the plot
  
  for(i in 1:length(rects$xstart)){
    rect(xleft = rects$xstart[i], ybottom = -20, xright = rects$xend[i], ytop = 20, 
         col = "lightgrey", border = NA)
  }
  
  lines(mts.df$mts ~ mts.df$date, col='blue', lwd=2)              #filtered mean
  lines(mts.df$obs ~ mts.df$date, col='black', lwd=2)             #observations
  lines(low.limit.mts ~ mts.df$date, col='blue', lty = "dashed")  #filtered variance
  lines(high.limit.mts ~ mts.df$date, col='blue', lty = "dashed") #filtered variance
  
  legend(x = "topright",                                   # Position
         legend = c("obs", "mts", "95% CI", "cycle"),      # Legend texts
         col = c("black", "blue", "blue", NA),             # Line colors
         lwd = c(1, 1, 1, NA),                             # Line thickness
         lty = c(1, 1, 2, NA),                             # Line type
         fill = c(NA, NA, NA, "lightgrey"),                # Fill for rectangles 
         border = NA,                                      # No border for rectangles
         cex=0.7)
}
#turn off PDF plotting
dev.off() 

