### Codes for the Salmon study - Hierarchical DLM
#### With 2 harmonic waves for country level

#libraries
library(lubridate)
library(dplyr)


################################################################################################################
### CREATE LEARNING AND TEST SETS ###

# Function to create Learning and Test sets for DLMs ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step where the datasets are divided into learning and test sets
learning.test.sets <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){ #what to do with the nseqs that are cut
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){ #if it's in the middle goes to Learning.set
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{#put the others in learning or test set
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100  
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100       
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)             
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){ 
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){ 
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){ #if has less than 6 observations != NA
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ #if has only 1 nseq - remove farm
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 always equal to NA and 
  # - farms with only 1 nseq in Learning.set  
  all.sites.exclude <- c(sites.exclude.all.na, sites.exclude.1.nseq) #52
  Learning.set <- subset(Learning.set, Learning.set$site %in% all.sites.exclude == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% all.sites.exclude == FALSE)            
  
  
  #remove the farms that all environmental data is NA - none
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)            
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) #sites that are in Learning.set but not in Test.set
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) #sites in Test.set that are not in Learning.set
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)            
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}



# Function to create Learning and Test sets for EM algorithm ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step where the datasets are divided into learning and test sets
## Difference is that doesn't remove farms on the Learning set with only 1 nseq (!=0)
learning.test.sets.EM <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){ 
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){ 
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){ #if it's in the middle goes to Learning.set
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{#put the others in learning or test set
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100 
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100      
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)             
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){ 
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){ 
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){ #if has less than 6 observations != NA
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ #if has only 1 nseq - don't remove farm
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with mortality.rel.20/log0.mortality.rel.10 always equal to NA
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.exclude.all.na == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.exclude.all.na == FALSE)            
  
  
  #remove the farms that all environmental data is NA - none
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)            
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) #sites that are in Learning.set but not in Test.set
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) #sites in Test.set that are not in Learning.set
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE)
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)            
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}



# Function to get the final Learning and Test sets ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not, here hierarchical=TRUE
# - N is the time step where the datasets are divided into learning and test sets
## Use only the sites that will be there after dividing the Learning set twice
## (so that the sites in the train.set used for the EM algorithm will be the same as in the Learning.set)
get.learning.test.sets <- function(df, relevant.names, hierarchical, N){
  
  df_0 <- df[-c(which(df$nseq==0)),] 
  
  sets <- learning.test.sets(df_0, relevant.names, N)
  Learning.set_0 <- sets[["Learning.set"]]
  Test.set_0 <- sets[["Test.set"]]
  
  N.EM <- round(3*(length(unique(Learning.set_0[order(Learning.set_0[,"date"]), "date"]))/4))
  
  sets.EM <- learning.test.sets.EM(df=Learning.set_0, relevant.names, N=N.EM)
  train.set.EM_0 <- sets.EM[["Learning.set"]]
  test.set.EM_0 <- sets.EM[["Test.set"]]
  
  sites.to.remove <- setdiff(Learning.set_0$site, train.set.EM_0$site)
  sites.to.keep <- unique(train.set.EM_0$site)
  
  if(hierarchical==FALSE){
    Learning.set.final <- subset(Learning.set_0, Learning.set_0$site %in% sites.to.remove == FALSE)
    Test.set.final <- subset(Test.set_0, Test.set_0$site %in% sites.to.remove == FALSE)
  }
  if(hierarchical==TRUE){
    sets <- learning.test.sets(df, relevant.names, N)
    Learning.set <- sets[["Learning.set"]]
    Test.set <- sets[["Test.set"]]
    
    Learning.set.final <- subset(Learning.set, Learning.set$site %in% sites.to.keep == TRUE)
    Test.set.final <- subset(Test.set, Test.set$site %in% sites.to.keep == TRUE)
  }
  return(list(Learning.set=Learning.set.final, Test.set=Test.set.final))
}


################################################################################################################
### MODELLING ###

# Function to estimate the mu0 (initial parameter vector) ----
# - D is the learning set with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle for the learning set - each column should be a farm and each row a date
# - expected.start.time is the observation time when we expect the model to start
# - H.w.country is the number of harmonic waves we want to model at country level
# - regions is a vector containing the name of the regions
get.m0 <- function(D, D.months, expected.start.time, H.w.country, regions){
  
  ### Elements for the initial parameter vector (mu0)
  mu0 <- c()
  Spline.list <- list()
  
  # Remove outliers, based on overall mean and moving SD
  all.values <- unlist(D[,c(2:ncol(D))])
  ylim <- range(na.omit(all.values))
  Mean <- mean(na.omit(all.values))
  SD <- sd(na.omit(all.values))
  upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
  lower <- Mean - 3*SD
  for(farm in colnames(D[2:ncol(D)])){
    remove.i <- which(D[,farm] > upper | D[,farm] < lower)
    D[remove.i, farm] <- NA
  }
  
  #Country - with harmonic waves
  ## required packages
  library('MuMIn')
  lmer <- lme4::lmer
  
  ## define the base frequency
  w <- (2*pi)/12
  
  ## we make an empty vector, into which we will add the relevant text for the linear function
  lm.vector <- c()
  
  ## we also make an empty vector, into which we will add the names for the parameter vector
  names.vector <- c()
  
  ## we add the target variable (the variable we model with this collection of harmonics)
  lm.vector <- c(lm.vector, "var", '~')
  names.vector <- c(names.vector, 'Country')
  
  ## now we add the relevant number of harmonics, one at a time
  for(i in 1:H.w.country){
    if(i == 1){
      a <- paste('cos(', i, '*w*', "month", ') + sin(', i, '*w*', "month", ')')
      lm.vector <- c(lm.vector, a)
      names <- c(paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
      
      names.vector <- c(names.vector, paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    }else{
      a <- paste('+ cos(', i, '*w*', "month", ') + sin(', i, '*w*', "month", ')')
      lm.vector <- c(lm.vector, a)
      names.vector <- c(names.vector, paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    }
  }
  
  ## add monthly country mean and month of the year in D
  library(lubridate)
  D2 <- D
  D2$var <- rowMeans(D[,2:ncol(D)], na.rm=T) #national average per date
  D2$month <- month(as.POSIXlt(D$date, format="%Y/%m/%d"))
  
  ## now we make it into a linear function - without random effect
  b1 <- as.formula(paste(lm.vector, collapse=' '))
  lm.1 <- lm(b1, data = D2)
  S1 <- summary(lm.1)
  
  ## now we directly have the initial parameter vector
  m0.H.w <- matrix(S1$coefficients[,'Estimate'])
  rownames(m0.H.w) <- names.vector
  mu0 <- rbind(mu0, m0.H.w)
  
  
  #Regions
  regions.means <- c()
  
  for (region in regions){
    ## get the farms that belong to the current region
    farms.per.region <- grep(region, names(D), value=TRUE)
    
    ## mean of mortality of the farms that belong to that region
    m0.R <- mean(unlist(D[,farms.per.region]), na.rm=T)
    
    ## save the mean to a vector for latter
    regions.means <- c(regions.means, m0.R)
    names(regions.means)[length(regions.means)] <- region
    
    ## get the difference between the mean of the region and the harmonic waves of the country
    for(i in 1:length(m0.H.w)){
      if(i == 1){
        x <-  m0.R - m0.H.w[i,]
      } else{
        m0.R.C <- x - m0.H.w[i,]
        x <- m0.R.C
      }
    }
    rm(x)
    
    ## add that value to mu0
    mu0 <- rbind(mu0, m0.R.C)
    row.names(mu0)[nrow(mu0)] <- region
  }
  
  #Farms
  for(farm in colnames(D[2:ncol(D)])){
    
    ## Level is 0 and trend is 1 for all farms
    mu <- c(0, 1)
    
    ## Give names to columns and rows
    mu.names <- c(farm, paste('d.', farm, sep=''))
    mu <- matrix(mu)
    row.names(mu) <- mu.names
    
    ## Initial parameter vector for all variables (mu0)
    mu0 <- rbind(mu0, mu)
    colnames(mu0) <- "mu0"
    
    ## Get the spline function to be used in Gt
    y.i <- D[,farm]
    y <- y.i[-which(is.na(y.i))]
    
    x.i <- D.months[,farm]
    x <- x.i[-which(is.na(y.i))] #because what is NA in mortality (not known or fallow period) those months have to be removed
    
    Spline = smooth.spline(x = x, y = y)
    
    plot(y~x, xlab="month", ylab=farm, main = paste("Spline", '=', farm))
    lines(Spline, col='red', lwd=3)
    
    ## Save the spline function - we will need it for the Gt matrix
    spline.name <- paste(farm, '_Spline', sep='')
    Spline.list[[length(Spline.list)+1]] <- Spline
    names(Spline.list)[length(Spline.list)] <- spline.name
  }
  
  return(list('mu0'=mu0, 'Spline.list'=Spline.list))
}



# Function to estimate C0 (prior variance) ----
# - D is the learning set with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle for the learning set - each column should be a farm and each row a date
# - expected.start.time is the observation time when we expect the model to start
# - expected.finish.time is the number of time steps we want to include to estimate C0 (we used expected.finish.time=6)
# - H.w.country is the number of harmonic waves we want to model at country level
# - regions is a vector containing the name of the regions
get.C0 <- function(D, D.months, expected.start.time, expected.finish.time, H.w.country, regions){
  
  # Remove outliers, based on overall mean and moving SD
  all.values <- unlist(D[,c(2:ncol(D))])
  ylim <- range(na.omit(all.values))
  Mean <- mean(na.omit(all.values))
  SD <- sd(na.omit(all.values))
  upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
  lower <- Mean - 3*SD
  for(farm in colnames(D[2:ncol(D)])){
    remove.i <- which(D[,farm] > upper | D[,farm] < lower)
    D[remove.i, farm] <- NA
  }
  
  # Make the empty C0 matrix
  C.all <- list()
  n.rows.cols <- 1+(2*H.w.country)+length(regions)+length(colnames(D[2:ncol(D)]))*2
  C0 <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
  
  # Give row and column names to C0
  country.names <- c()
  for(i in 1:H.w.country){
    names <- c(paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    country.names <- c(country.names, names)
  }
  country.names <- c("Country", country.names)
  regions.names <- regions
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){ 
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.col.names <- c(country.names, regions.names, level.trend.names) 
  row.names(C0) <- row.col.names
  colnames(C0) <- row.col.names
  
  # Make an empty table for saving regions and farms means
  means <- as.data.frame(matrix(NA, ncol=length(regions)+length(colnames(D[2:ncol(D)])), nrow=length((expected.start.time):(expected.start.time+expected.finish.time))))
  rownames(means) <- (expected.start.time):(expected.start.time+expected.finish.time)
  
  colnames(means) <- c(regions, colnames(D[2:ncol(D)]))
  
  
  # Elements for the prior variance (C0)
  #Country - harmonic waves
  
  ## required packages
  library('MuMIn')
  lmer <- lme4::lmer
  
  ## define the base frequency
  w <- (2*pi)/12
  
  ## we make an empty vector, into which we will add the relevant text for the linear function
  lm.vector <- c()
  
  ## we also make an empty vector, into which we will add the names for the parameter vector
  names.vector <- c()
  
  ## we add the target variable (the variable we model with this collection of harmonics)
  lm.vector <- c(lm.vector, "var", '~')
  names.vector <- c(names.vector, 'Country')
  
  ## now we add the relevant number of harmonics, one at a time
  
  for(i in 1:H.w.country){
    if(i == 1){
      a <- paste('cos(', i, '*w*', "month", ') + sin(', i, '*w*', "month", ')')
      lm.vector <- c(lm.vector, a)
      names <- c(paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
      
      names.vector <- c(names.vector, paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    }else{
      a <- paste('+ cos(', i, '*w*', "month", ') + sin(', i, '*w*', "month", ')')
      lm.vector <- c(lm.vector, a)
      names.vector <- c(names.vector, paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    }
  }
  
  ## add monthly country mean and month of the year in D
  library(lubridate)
  D2 <- D
  D2$var <- rowMeans(D[,2:ncol(D)], na.rm=T)
  D2$month <- month(as.POSIXlt(D$date, format="%Y/%m/%d"))
  
  ## now we make it into a linear function - without random effect
  b1 <- as.formula(paste(lm.vector, collapse=' '))
  lm.1 <- lm(b1, data = D2)
  S1 <- summary(lm.1)
  
  ## now we directly have the initial parameter vector
  m0.H.w <- matrix(S1$coefficients[,'Estimate'])
  
  ## and we directly have the initial prior variance matrix
  C <- as.matrix(vcov(lm.1))
  
  ## finalize and save in C0
  C0[names.vector, names.vector] <- C
  
  ## save all C matrices to see where there is NA values and replace them afterwards
  C.all[["Country"]] <- C
  
  
  #Regions
  for (region in regions){
    
    ## Get the farms that belong to the current region
    indx.region <- grepl(region, colnames(D[,c(2:ncol(D))]))
    farms.region <- colnames(D[,c(2:ncol(D))])[indx.region]
    
    ## Create table with the relevant information
    D.A <- data.frame(region=as.vector(unlist(D[,farms.region])), months=as.vector(unlist(D.months[,farms.region])))
    #only consider the first n months
    D.A <- subset(D.A, D.A[,"months"] %in% (expected.start.time):(expected.start.time+expected.finish.time))
    colnames(D.A)[1] <- region
    
    ## Get the mean of each region for the first 7 months (in our case) 
    for(i in (expected.start.time):(expected.start.time+expected.finish.time)){ 
      D.A.i <- subset(D.A, D.A$months == i)
      means[row.names(means)[row.names(means)==i], region] <- mean(D.A.i[, region], na.rm=T)
    }
    
    ## Get the difference between the means of the region and the harmonic waves of the country
    R.C.all <- c()
    for(r in (expected.start.time):(expected.start.time+expected.finish.time)){
      for(i in 1:length(m0.H.w)){
        if(i == 1){
          x <-  means[row.names(means)[row.names(means)==r], region] - m0.H.w[i,]
        } else{
          R.C <- x - m0.H.w[i,]
          x <- R.C
        }
      }
      R.C.all <- c(R.C.all, R.C)
      rm(x)
    }
    
    ## Finalize and save C0
    C <- var(R.C.all)
    C0[region, region] <- C
    
    ## Save all C matrices to see where there is NA values and replace them afterwards
    C.all[[region]] <- C
  }
  
  
  #Farms
  for(farm in colnames(D[2:ncol(D)])){
    
    ## Create table with the relevant information
    D.A <- data.frame(Farm=D[,farm], months=D.months[,farm])
    #only consider the first n months
    D.A <- subset(D.A, D.A[,"months"] %in% (expected.start.time):(expected.start.time+expected.finish.time))
    colnames(D.A)[1] <- farm
    
    ## Get the mean of each farm for the first 7 months (in pour case)
    for(i in (expected.start.time):(expected.start.time+expected.finish.time)){
      D.A.i <- subset(D.A, D.A$months == i)
      means[row.names(means)[row.names(means)==i], farm] <- mean(D.A.i[, farm], na.rm=T)
    }
    
    ## Get the region which this farm belongs to
    region <- sapply(strsplit(farm, '_'), `[`, 2) #get what is after _
    
    ## Get the difference between the means of the farm, region and the harmonic waves of the country
    F.R.C.all <- c()
    for(r in (expected.start.time):(expected.start.time+expected.finish.time)){
      for(i in 1:length(m0.H.w)){
        if(i == 1){
          x <-  means[row.names(means)[row.names(means)==r], farm] - means[row.names(means)[row.names(means)==r], region] - m0.H.w[i,]
        } else{
          F.R.C <- x - m0.H.w[i,]
          x <- F.R.C
        }
      }
      F.R.C.all <- c(F.R.C.all, F.R.C)
      rm(x)
    }
    
    ## Difference between the present observation and the previous observation (deviations)
    Diff <- diff(F.R.C.all)
    F.R.C.Diff <- c(NA, Diff)
    
    ## Combine the differences and the deviations
    df <- data.frame(F.R.C=F.R.C.all, F.R.C.Diff=F.R.C.Diff)
    df <- na.omit(df)
    
    ## Finalize and save C0
    C <- cov(df[,c("F.R.C", "F.R.C.Diff")])
    colnames(C) <- c(farm, paste("d.", farm, sep=""))
    row.names(C) <- c(farm, paste("d.", farm, sep=""))
    C0[c(farm, paste("d.", farm, sep="")), c(farm, paste("d.", farm, sep=""))] <- C
    C0 <- as.matrix(C0)
    
    ## Save all C matrices to see where there is NA values and replace them afterwards
    C.all[[farm]] <- C
    
  }
  # - Test if C0 is legal: if length(negs) is 0, the matrix is positive semi definite.
  eig = eigen(C0)
  eigVal = eig$values
  negs = eigVal[eigVal < 0]
  if(length(negs)!=0){
    warning("C0 is not positive semi definite.") 
  }
  return(list('C0'=C0, 'C.all'=C.all))
}



# Function to estimate the Ft (design matrix) ----
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - H.w.country is the number of harmonic waves we want to model at country level
# - regions is a vector containing the name of the regions
get.Ft <- function(D, H.w.country, regions){
  
  # Make the empty Ft matrix with the correct dimensions
  n.rows <- 1+(2*H.w.country)+length(regions)+length(colnames(D[2:ncol(D)]))*2
  n.cols <- length(colnames(D[2:ncol(D)]))
  Ft <- as.data.frame(matrix(0, ncol = n.cols, nrow = n.rows))
  
  # Give row and column names to Ft
  # - Rows
  country.names <- c()
  for(i in 1:H.w.country){
    names <- c(paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    country.names <- c(country.names, names)
  }
  country.names <- c("Country", country.names)
  regions.names <- regions
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){ 
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.names <- c(country.names, regions.names, level.trend.names) 
  row.names(Ft) <- row.names
  # - Columns
  colnames(Ft) <- colnames(D[2:ncol(D)])
  
  # Ft matrix
  for(farm in colnames(Ft)){ 
    # - Put the value 1 in the country part for all farms
    Ft["Country", farm] <- 1
    # - Put the value 1 in the respective region of the current farm
    region <- sapply(strsplit(farm, '_'), `[`, 2) #get what is after _
    Ft[region, farm] <- 1
    # - Put the value 1 in the current farm
    Ft[farm, farm] <- 1
  }
  # - Put 1 and 0 for the harmonic waves elements for country
  for(i in 1:H.w.country){
    Ft[paste("Country.", i, ".H.w", sep=""),] = 1
    Ft[paste("d.Country.", i, ".H.w", sep=""),] = 0
  }
  return(as.matrix(Ft))
}



# Function to estimate the Gt (system matrix) ----
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - i is the time step you are modelling in the DLM, you can also give it equal to NA if you want to use dates instead
# - Date is the date you are currently modelling in the DLM
# - Date_1 is the previously date to what you are currently modelling in the DLM
# - H.w.country is the number of harmonic waves we want to model at country level
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector containing the name of the regions
get.Gt <- function(D, D.months, i, Date, Date_1, H.w.country, Spline.list, regions){ 
  
  # Make the empty Gt matrix
  n.rows.cols <- 1+(2*H.w.country)+length(regions)+length(colnames(D[2:ncol(D)]))*2
  Gt <- as.data.frame(diag(1, ncol = n.rows.cols, nrow = n.rows.cols)) #put 1 in all diagonal
  
  # Give row and column names to Gt
  country.names <- c()
  for(i in 1:H.w.country){
    names <- c(paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep=""))
    country.names <- c(country.names, names)
  }
  country.names <- c("Country", country.names)
  regions.names <- regions
  level.trend.names <- c()
  for(farm in colnames(D[2:ncol(D)])){
    level.trend.names <- c(level.trend.names, farm, paste('d.', farm, sep=''))      
  }
  row.col.names <- c(country.names, regions.names, level.trend.names) 
  row.names(Gt) <- row.col.names
  colnames(Gt) <- row.col.names
  
  
  # Elements for the system matrix (Gt)
  
  ## Country
  # - define the base frequency
  w <- (2*pi)/12
  
  # - add the harmonic waves for the country
  for(i in 1:H.w.country){
    Gt[paste("Country.", i, ".H.w", sep=""), paste("Country.", i, ".H.w", sep="")] = cos(i*w)
    Gt[paste("Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep="")] = sin(i*w)
    Gt[paste("d.Country.", i, ".H.w", sep=""), paste("Country.", i, ".H.w", sep="")] = -sin(i*w)
    Gt[paste("d.Country.", i, ".H.w", sep=""), paste("d.Country.", i, ".H.w", sep="")] = cos(i*w)
  }
  
  ## Region
  # - nothing needs to be done since it's only the value 1 (already added in line 428)
  
  ## Farm
  
  for(farm in colnames(D[2:ncol(D)])){ 
    
    # - get the relevant spline function
    spline.name <- paste(farm, 'Spline', sep='_') #1 spline per farm
    Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
    
    # - get the date index
    i <- which(unique(D$date) == Date)
    i_1 <- which(unique(D$date) == Date_1)
    
    if(isTRUE(length(i)!=0)){
      
      time <- D.months[i,farm] #The observation time of the current observations
      time_1 <- D.months[i_1,farm] #The observation time of the previous observations
      if(length(time_1) == 0){
        time_1 <- 0
      }
      
      # - get trend for Gt matrix
      if(!is.na(time)){
        pred <- predict(object = Spline, x = time)$y
        if(isTRUE(time != time_1)){
          pred_1 <- predict(object = Spline, x = time_1)$y
          Trend <- pred - pred_1
        }else{
          Trend <- pred 
        }
      }else{
        Trend <- 0
      }
      
      Gt[farm, paste("d.", farm, sep="")] <- Trend
      
    }else{
      Gt[farm, paste("d.", farm, sep="")] <- 0
    }
  }
  return(as.matrix(Gt))
}



# Function to estimate the V (observation variance) ----
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - countryVar is the estimated observation variance for the country level
# - regionVar is the estimated observation variance for the region level
# - farmVar is the estimated observation variance for the farm level
get.Vt <- function(D, countryVar, regionVar, farmVar){
  
  # Make the empty Vt matrix
  n.rows.cols <- length(colnames(D[2:ncol(D)]))
  Vt <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
  
  # Give row and column names to Ft
  rownames(Vt) <- colnames(D[2:ncol(D)])
  colnames(Vt) <- colnames(D[2:ncol(D)])
  
  for(farm in rownames(Vt)){ 
    
    #Country
    Vt[farm, ] <- countryVar #all farms belong to the same country so all have countryVar 
    
    #Region
    region <- sapply(strsplit(farm, '_'), `[`, 2) #get what is after _
    indx.region <- grepl(region, rownames(Vt))
    farms.region <- rownames(Vt)[indx.region]
    Vt[farm, farms.region] <- countryVar + regionVar #covariances of the farms that belong to the same region is countryVar + regionVar
    
    #Farm
    Vt[farm, farm] <- countryVar + regionVar + farmVar #diagonal is countryVar + regionVar + farmVar
  }
  Vt <- as.matrix(Vt)
  return(Vt)
}



# Function to update mu0 ----
## when a dataset does not start in January we have to update the country levels (harmonic waves) in mu0 for the correspondent month
# - mu0 is the initial parameter vector
# - Gt is the system matrix
# - start.time is the month number of the dataset's first date
update.mu0 <- function(mu0, Gt, start.time){
  for(i in 1:(start.time-1)){
    mu0 <- Gt %*% mu0
  }
  return(mu0)
}



# Function to initialize at ----
## when a new production cycle starts the at for that farm should be initialized (= Gt %*% mu0) 
# - at is the prior mean
# - mu0 is the initial parameter vector
# - Gt is the system matrix
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - Date is the date you are currently modelling in the DLM
get.at.hierarchical <- function(at, mu0, Gt, D.months, Date){
  
  #Observation vector only for that Date
  Date.set.months <- D.months[which(unique(D.months$date) == Date),2:ncol(D.months)]
  
  #Farms that started a new production cycle on that Date
  Date.first.set <- which(Date.set.months==0)
  farms.first <- colnames(Date.set.months[Date.first.set])
  
  at.h <- at
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){
      
      #get level and trends names for those farms
      farm.name.mu0.l <- farm
      farm.name.mu0.t <- paste("d", farm, sep=".")
      
      #replace at for those farms by Gt %*% mu0 (initialize)
      at.h[c(farm.name.mu0.l, farm.name.mu0.t),] <- as.matrix(Gt[c(farm.name.mu0.l, farm.name.mu0.t),c(farm.name.mu0.l, farm.name.mu0.t)]) %*% 
        as.matrix(mu0[c(farm.name.mu0.l, farm.name.mu0.t),])
    }
  }
  return(at.h)
}



# Function to initialize Gt ----
## when a new production cycle starts the Gt for that farm should be initialized (=0) 
# - Gt is the system matrix
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - Date is the date you are currently modelling in the DLM
get.Gt.hierarchical <- function(Gt, D.months, Date){
  
  #Observation vector only for that Date
  Date.set.months <- D.months[which(unique(D.months$date) == Date),2:ncol(D.months)]
  
  #Farms that started a new production cycle on that Date
  Date.first.set <- which(Date.set.months==0)
  farms.first <- colnames(Date.set.months[Date.first.set])
  
  Gt.h <- Gt
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){
      
      #get level and trends names for those farms
      farm.name.Gt.l <- farm
      farm.name.Gt.t <- paste("d", farm, sep=".")
      
      #replace Gt for those farms by 0
      Gt.h[farm.name.Gt.l, farm.name.Gt.l] <- 0
      Gt.h[farm.name.Gt.l, farm.name.Gt.t] <- 0
      Gt.h[farm.name.Gt.t, farm.name.Gt.l] <- 0
      Gt.h[farm.name.Gt.t, farm.name.Gt.t] <- 0
    }
  }
  return(Gt.h)
}



# Function to create an extraVariance matrix ----
## when a new production cycle starts the we add extraVariance to Wt (equal to C0 values for that farm) 
# - C0 is the prior variance
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - Date is the date you are currently modelling in the DLM
add.extraVariance <- function(C0, D.months, Date){
  
  #Observation vector only for that Date
  Date.set.months <- D.months[which(unique(D.months$date) == Date),2:ncol(D.months)]
  
  #Farms that started a new production cycle on that Date
  Date.first.set <- which(Date.set.months==0)
  farms.first <- colnames(Date.set.months[Date.first.set])
  
  #Create a zero matrix with the same dimensions as Wt/C0
  zeroMatrix <- matrix(0, ncol = length(C0[,1]), nrow = length(C0[,1]))
  rownames(zeroMatrix) <- rownames(C0)
  colnames(zeroMatrix) <- colnames(C0)
  
  extraVariance <- zeroMatrix
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){
      
      #get level and trends names for those farms
      farm.name.C0.l <- farm
      farm.name.C0.t <- paste("d", farm, sep=".")
      
      #Add C0 for that farm to extraVariance
      extraVariance[farm.name.C0.l, farm.name.C0.l] <- C0[farm.name.C0.l, farm.name.C0.l]
      extraVariance[farm.name.C0.l, farm.name.C0.t] <- C0[farm.name.C0.l, farm.name.C0.t]
      extraVariance[farm.name.C0.t, farm.name.C0.l] <- C0[farm.name.C0.t, farm.name.C0.l]
      extraVariance[farm.name.C0.t, farm.name.C0.t] <- C0[farm.name.C0.t, farm.name.C0.t]
    }
  }
  return(extraVariance)
}



# Function to get a (hard coded) factor to multiply to the extraVariance (if needed) ----
getFactor = function() {
  return(1.0)
}



# Function to check if matrices are positive semi definite ----
# - text is the name of the matrix we are checking
# - step is the time step you are modelling in the DLM
# - varMat is the matrix being checked
checkValidity = function(text, step, varMat) {
  eig = eigen(varMat)
  eigVal = eig$values
  negs = eigVal[eigVal < 0]
  if(length(negs)!=0){
    print(paste(text, " is not positive semi definite at step ", step, sep = ""))
    return(FALSE)
  }
  return(TRUE)
}



# Function for the DLM for hierarchical models with deltas ----
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - countryVar is the estimated observation variance for the country level
# - regionVar is the estimated observation variance for the region level
# - farmVar is the estimated observation variance for the farm level
# - deltas is a vector with the estimated discount factors (Country, Region and Farm)
# - H.w.country is the number of harmonic waves we want to model at country level
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector containing the name of the regions
runDiscountDLM <- function(D, D.months, mu0, C0, countryVar, regionVar, farmVar, deltas, H.w.country, Spline.list, regions){
  
  # Get the total number of time steps in the dataset D and the first date available
  n <- 1:length(unique(D$date))
  first.date <- as.Date(D$date)[1]
  
  # Create empty lists
  Yt.list <- list()
  at.list <- list()		 
  Rt.list <- list()
  ft.list <- list()
  Qt.list <- list()
  At.list <- list()
  et.list <- list()
  ut.list <- list()
  mt.list <- list()
  Ct.list <- list()
  Ft.list <- list()
  VSE.list <- list()
  Gt.list <- list()
  fullFt.list <- list()
  
  # Create a list for discount groups to be used in creating Wt
  discount.groups <- list()
  discount.groups[[1]] <- 1:(1+(H.w.country*2))                                                        #country
  discount.groups[[2]] <- (last(discount.groups[[1]])+1):(last(discount.groups[[1]])+length(regions))  #regions
  discount.groups[[3]] <- (last(discount.groups[[2]])+1):dim(C0)[1]                                    #all farms (levels and trends)
  
  # Update country level of mu0 if the dataset does not start in January
  if(month(first.date)!=1){
    m <- month(first.date)
    Gt <- get.Gt(D, D.months, i=NA, Date=first.date, Date_1=0, H.w.country, Spline.list, regions)
    mu0.all <- update.mu0(mu = mu0, Gt = Gt, start.time = m)
    mu0[1:(1+(H.w.country*2)),] <- mu0.all[1:(1+(H.w.country*2)),]
  }
  mt <- mu0				
  Ct <- C0
  
  for(i in n){
    
    # Get dates (current and previous)
    Date <- sort(unique(D$date))[i]
    Date_1 <- sort(unique(D$date))[i-1]
    if(length(Date_1)==0){
      Date_1 <- 0
    }
    
    # Get the observation vector (Yt)
    Yt <- D[which(unique(D$date) == Date),2:ncol(D)]
    
    # Get the observational variance (Vt)
    Vt <- get.Vt(D, countryVar, regionVar, farmVar)
    
    # Get design matrix (Ft) 
    Ft <- get.Ft(D, H.w.country, regions)
    
    # Get the system matrix (Gt) 
    Gt <- get.Gt(D, D.months, i=NA, Date, Date_1, H.w.country, Spline.list, regions)
    
    
    # Remove missing values from Yt, Vt and Ft
    missing = which(is.na(Yt))
    fullYt <- Yt
    fullFt <- Ft
    fullVt <- Vt
    
    # - If length of missing is > 0 there is at least one missing
    if (length(missing) > 0) {
      
      # remove from Yt
      a <- colnames(fullYt)[-missing]
      Yt <- t(as.matrix(fullYt[-missing]))
      
      # Remove from Ft
      Ft <- as.matrix(fullFt[, -missing])
      
      # Remove from Vt
      Vt <- fullVt[-missing, -missing]
    }
    
    
    ## Initialize at and Gt when a new production cycle starts or when starts to have information for one farm
    if(Date > first.date){
      at <- Gt %*% mt   
      at <- get.at.hierarchical(at, mu0, Gt, D.months, Date)
      Gt <- get.Gt.hierarchical(Gt, D.months, Date)
    }else{
      at <- Gt %*% mt                                
    }
    
    # Start with the Kalman filter
    Pt <- Gt %*% Ct %*% t(Gt)      
    
    # Define system variance
    Wt <- matrix(0, ncol = length(C0[,1]), nrow = length(C0[,1]))
    for (b in 1:3) {
      Wt[discount.groups[[b]], discount.groups[[b]]] = ((1-deltas[b])/deltas[b]) * Pt[discount.groups[[b]], discount.groups[[b]]]
    }
    
    # Add extravariance to Wt
    extraVariance <- add.extraVariance(C0, D.months, Date)
    Wt = Wt + getFactor()*extraVariance
    
    # Return to Kalman Filter
    Rt = Pt + Wt                                   # Prior Variance
    Rt = (Rt + t(Rt))/2                            # Make sure Rt is symmetrical
    ft = t(Ft) %*% at                              # One-step Forecast mean
    ft2 = t(fullFt) %*% at                         # Save this one to have estimates of ft when we don't have observations
    Qt = t(Ft) %*% Rt %*% Ft + Vt                  # One-step Forecast variance
    Qt2 = t(fullFt) %*% Rt %*% fullFt + fullVt     # Save this one to have estimates of Qt when we don't have observations
    At = Rt %*% Ft %*% solve(Qt)                   # Adaptative Coef. matrix, also called Kalman gain
    et = Yt - ft	                                 # One-step forecast error
    ut = et / sqrt(diag(Qt))                       # Standardized forecast error
    
    # Update the parameter vector and variance matrix
    mt = at + At %*% et                            # Filtered mean
    Ct = Rt - At  %*% Qt %*% t(At)	               # Filtered variance
    
    # Make sure Ct is symmetrical
    Ct = (Ct + t(Ct))/2
    
    # Check if Rt and Qt are positive semidefinite
    checkValidity("Rt", i, Rt)
    checkValidity("Qt", i, Qt)
    
    # Save the values in lists
    Yt.list[[i]] <- Yt
    at.list[[i]] <- at
    Rt.list[[i]] <- Rt
    ft.list[[i]] <- ft2
    Qt.list[[i]] <- Qt2
    At.list[[i]] <- At
    et.list[[i]] <- et
    ut.list[[i]] <- ut
    mt.list[[i]] <- mt
    Ct.list[[i]] <- Ct
    Ft.list[[i]] <- t(Ft)
    Gt.list[[i]] <- Gt
    fullFt.list[[i]] <- t(fullFt)
  }
  
  return(list(
    Yt=Yt.list,
    at=at.list,
    Rt=Rt.list,
    ft=ft.list,
    Qt=Qt.list,
    At=At.list,
    et=et.list,
    ut=ut.list,
    mt=mt.list,
    Ct=Ct.list,
    F=Ft.list,
    Gt.list=Gt.list,
    fullFt.list=fullFt.list))
}



# Function of the Kalman Smoother ----
# - res is the result returned from the DLM filter (runDiscountDLM)
runSmoother <- function(res) {
  
  # Create empty mts (smoothed mean) and Cts (smoothed variance) with the correct dimensions
  n = length(res$mt)
  p = length(res$mt[[1]])
  mts <- array(NA,dim=c(p,1,n));
  Cts <- array(NA,dim=c(p,p,n));
  
  # Put last value equal to filtered
  mts[,,n] <- res$mt[[n]]
  Cts[,,n] <- res$Ct[[n]]
  
  # These are useful
  Bt <- array(NA,dim=c(p,p,n))
  Lt <- array(NA,dim=c(p,p,n))
  
  # Iterate backwards over months
  for(i in ((n-1):1))   {
    
    # Get Gt and Rt
    Gt <- res$Gt.list[[i+1]]
    res$Rt[[i+1]] <- as.matrix(res$Rt[[i+1]])
    
    # Kalman smoother
    Bt[,,i] <- as.matrix( res$Ct[[i]] %*% t(Gt) %*% solve(res$Rt[[i+1]]) )
    mts[,,i] <- res$mt[[i]] + Bt[,,i] %*% (mts[,,i+1] - res$at[[i+1]])
    Cts[,,i] <- as.matrix( res$Ct[[i]] + Bt[,,i] %*% (Cts[,,i+1] - res$Rt[[i+1]]) %*% t(Bt[,,i]) )
    
    # Make sure Ct is symmetrical
    Cts[,,i] <- (Cts[,,i] + t(Cts[,,i]))/2
    
    # Check for being positive semidefinite
    checkValidity("Cts", i, Cts[,,i])
    
    mts[,,i] <- as.matrix(mts[,,i])
  }
  # give names
  rownames(mts) <- rownames(res$mt[[i]])
  colnames(mts) <- "mts"
  rownames(Cts) <- rownames(res$Ct[[i]])
  colnames(Cts) <- colnames(res$Ct[[i]])
  
  # Now when we are at it: Find L and store it you want to use it for the EM algorithm (we didn't use it)
  for(i in ((n):2))  {
    Lt[,,i] <- Cts[,,i] + Gt%*%Cts[,,i-1]%*%t(Gt) - Cts[,,i]%*%t(Bt[,,i-1]) - Bt[,,i-1]%*%Cts[,,i]
  }
  rownames(Lt) <- rownames(res$Ct[[i]])
  colnames(Lt) <- colnames(res$Ct[[i]])
  
  return(list(mts=mts,
              Cts=Cts,
              Lt=Lt,
              D=res$D));
}



# Function to extract the relevant information ----
# - res is the result returned from the DLM filter (runDiscountDLM)
# - smot is the result returned from the Kalman Smoother (runSmoother)
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - D.nseq is a dataset with information about the number of the production cycle (ID) - each column should be a farm and each row a date
# - H.w.country is the number of harmonic waves we want to model at country level
# - regions is a vector containing the name of the regions
extract.res <- function(res, smot, D, D.nseq, H.w.country, regions){
  
  # Get farms names and create an empty list to save the results
  farms <- colnames(D[2:ncol(D)])
  results.list <- list()
  
  # Get the filtered mean (mt)
  ## create empty data frame for mt results
  df.mt <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1)) #regions and country
  df.mt[,1] <- D$date
  colnames(df.mt) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## for farms 
  ##(sum of country, region and farm (use Ftfull))
  for(i in 1:length(res$mt)){ 
    f <- res$fullFt.list[[i]] %*% res$mt[[i]] #multiply: Ft * mt
    df.mt[i, rownames(f)] <- t(f)
  }
  
  ## for region 
  ##(sum of country and 1 region levels and the average of the levels of the farms for each region)
  for(region in regions){
    
    # - relevant vector to extract the country level
    vector.country <- rep(0, dim(res$mt[[1]])[1])
    names(vector.country) <- rownames(res$mt[[1]])
    
    if(is.null(H.w.country)){
      vector.country[1] <- 1
    }else{
      Fti <- 1
      for(n in 1:H.w.country){
        Fti <- c(Fti, c(1,0))
      }
      vector.country[c(1:(1+2*H.w.country))] <- Fti
    }
    
    # - relevant vector to extract the region and country levels
    vector.country.region <- vector.country
    vector.country.region[region] <- 1
    
    # - farms that belong to the current region
    farms.region <- farms[grepl(region, farms)]
    
    # - relevant vector to extract the region and country levels and the average of the levels of the farms for the current region
    for(i in 1:length(res$mt)){ 
      vector.country.region.farm <- vector.country.region
      farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
      vector.country.region.farm[farms.region.i] <- 1/length(farms.region.i)
      
      # - extract complete region level
      r <- vector.country.region.farm %*% res$mt[[i]]
      df.mt[i, region] <- r
    }
  }
  
  ## for country 
  ##(sum of the country and average of the region levels and the average of the levels of the farms divided by the number of regions)
  
  # - relevant vector to extract the country and average of the region levels
  vector.sum <- vector.country
  vector.sum[regions] <- 1/length(regions)
  
  # - relevant vector to extract the country and average of the region levels, and the average of the levels of the farms divided by the number of regions
  for(i in 1:length(res$mt)){ 
    vector <- vector.sum
    for(region in regions){
      # - farms that belong to that region
      farms.region <- farms[grepl(region, farms)]
      farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
      vector[farms.region.i] <- 1/length(farms.region.i)/length(regions)
    }
    
    # - extract complete country level
    c <- vector %*% res$mt[[i]]
    df.mt[i, "Country"] <- c
  }
  results.list[["mt"]] <- df.mt
  
  
  # Get the forecasts (ft)
  ## create empty data frame for ft results
  df.ft <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1)) 
  df.ft[,1] <- D$date
  colnames(df.ft) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## get ft
  for(i in 1:length(res$ft)){ 
    for(a in rownames(res$ft[[i]])){ #put the observations on the right date
      df.ft[i, a] <- res$ft[[i]][a,1]
    }
  }
  results.list[["ft"]] <- df.ft
  
  
  # Get the raw forecasts errors (et)
  ## create empty data frame for et results
  df.et <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1)) #regions and country
  df.et[,1] <- D$date
  colnames(df.et) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## get et
  for(i in 1:length(res$et)){ 
    for(a in rownames(res$et[[i]])){ #put the et that exist on the right date, everything else is NA
      df.et[i, a] <- res$et[[i]][a,1]
    }
  }
  results.list[["et"]] <- df.et
  
  
  # Get the standardized forecasts errors (ut)
  ## create empty data frame for ut results
  df.ut <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1)) #regions and country
  df.ut[,1] <- D$date
  colnames(df.ut) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## get ut
  for(i in 1:length(res$ut)){ 
    for(a in rownames(res$ut[[i]])){ #put the ut that exist on the right date, everything else is NA
      df.ut[i, a] <- res$ut[[i]][a,1]
    }
  }
  results.list[["ut"]] <- df.ut
  
  
  # Get the filtered variance (Ct) 
  ## create empty data frame for Ct results
  df.Ct <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1))
  df.Ct[,1] <- D$date
  colnames(df.Ct) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ##for farms 
  ##(considering country, region and farm (use Ftfull))
  for(i in 1:length(res$Ct)){ 
    f <- res$fullFt.list[[i]] %*% res$Ct[[i]] %*% t(res$fullFt.list[[i]]) #multiply: Ft * Ct * t(Ft)
    for(farm in colnames(f)){ 
      df.Ct[i, farm] <- f[farm, farm] #now get just the diagonals
    }
  }
  
  ## for region
  ##(sum of country and 1 region levels and the average of the levels of the farms for each region)
  for(region in regions){
    
    # - relevant vector to extract the region and country levels
    vector.country.region <- vector.country
    vector.country.region[region] <- 1
    
    # - farms that belong to the current region
    farms.region <- farms[grepl(region, farms)]
    
    # - relevant vector to extract the region and country levels and the average of the levels of the farms for the current region
    for(i in 1:length(res$Ct)){
      vector.country.region.farm <- vector.country.region
      farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
      vector.country.region.farm[farms.region.i] <- 1/length(farms.region.i)
      
      # - extract complete region level
      r <- vector.country.region.farm %*% res$Ct[[i]] %*% vector.country.region.farm 
      df.Ct[i, region] <- r
    }
  }
  
  ## for country 
  ##(sum of the country and average of the region levels and the average of the levels of the farms divided by the number of regions)

  # - relevant vector to extract the country and average of the region levels
  vector.sum <- vector.country
  vector.sum[regions] <- 1/length(regions)
  
  # - relevant vector to extract the country and average of the region levels, and the average of the levels of the farms divided by the number of regions
  for(i in 1:length(res$Ct)){ 
    vector <- vector.sum
    for(region in regions){
      # - farms that belong to the current region
      farms.region <- farms[grepl(region, farms)]
      farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
      vector[farms.region.i] <- 1/length(farms.region.i)/length(regions)
    }
    
    # - extract complete country level
    c <- vector %*% res$Ct[[i]] %*% vector
    df.Ct[i, "Country"] <- c
  }
  results.list[["Ct"]] <- df.Ct
  
  
  # Get the forecast variance (Qt)
  ## create empty data frame for Qt results
  df.Qt <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1)) #regions and country
  df.Qt[,1] <- D$date
  colnames(df.Qt) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
  
  ## get Qt
  for(i in 1:length(res$Qt)){ 
    for(a in rownames(res$Qt[[i]])){  #put the Qt that exist on the right date, everything else is NA
      df.Qt[i, a] <- res$Qt[[i]][a,a] #diagonal values
    }
  }
  results.list[["Qt"]] <- df.Qt
  
  
  # Get the smoothed values
  if (!is.null(smot)) {
    
    # Get the smoothed mean (mts)
    ## create empty data frame for mts results
    df.mts <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1))
    df.mts[,1] <- D$date
    colnames(df.mts) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
    
    ## for farms 
    ##(sum of country, region and farm (use Ftfull))
    for(i in 1:dim(smot$mts)[3]){ 
      f <- res$fullFt.list[[i]] %*% smot$mts[,,i] #multiply: Ft * mts
      df.mts[i, rownames(f)] <- t(f)
    }
    
    ## for region 
    ##(sum of country and 1 region levels and the average of the levels of the farms for each region)
    for(region in regions){
      
      # - relevant vector to extract the region and country levels
      vector.country.region <- vector.country
      vector.country.region[region] <- 1
      
      # - farms that belong to the current region
      farms.region <- farms[grepl(region, farms)]
      
      # - relevant vector to extract the region and country levels and the average of the levels of the farms for the current region
      for(i in 1:dim(smot$mts)[3]){ 
        vector.country.region.farm <- vector.country.region
        farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
        vector.country.region.farm[farms.region.i] <- 1/length(farms.region.i)
        
        # - extract complete region level
        r <- vector.country.region.farm %*% smot$mts[,,i]
        df.mts[i, region] <- r
      }
    }
    
    ## for country 
    ##(sum of the country and average of the region levels and the average of the levels of the farms divided by the number of regions)

    # - relevant vector to extract the country and average of the region levels
    vector.sum <- vector.country
    vector.sum[regions] <- 1/length(regions)
    
    # - relevant vector to extract the country and average of the region levels, and the average of the levels of the farms divided by the number of regions
    for(i in 1:dim(smot$mts)[3]){ 
      vector <- vector.sum
      for(region in regions){
        # - farms that belong to the current region
        farms.region <- farms[grepl(region, farms)]
        farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
        vector[farms.region.i] <- 1/length(farms.region.i)/length(regions)
      }
      
      # - extract complete country level
      c <- vector %*% smot$mts[,,i]
      df.mts[i, "Country"] <- c
    }
    results.list[["mts"]] <- df.mts
    
    
    # Get the smoothed variance (Cts) 
    ## create empty data frame for Cts results
    df.Cts <- data.frame(matrix(NA, nrow = dim(D)[1], ncol = dim(D)[2] + length(regions) + 1))
    df.Cts[,1] <- D$date
    colnames(df.Cts) <- c("date", "Country", regions, colnames(D[2:ncol(D)]))
    
    ## for farms
    ##(considering country, region and farm (use Ftfull))
    for(i in 1:dim(smot$Cts)[3]){ 
      f <- res$fullFt.list[[i]] %*% smot$Cts[,,i] %*% t(res$fullFt.list[[i]]) #multiply: Ft * Cts * t(Ft)
      for(farm in colnames(f)){ 
        df.Cts[i, farm] <- f[farm, farm] #now get just the diagonals
      }
    }
    
    ## for region
    ##(sum of country and 1 region levels and the average of the levels of the farms for each region)
    for(region in regions){
      
      # - relevant vector to extract the region and country levels
      vector.country.region <- vector.country
      vector.country.region[region] <- 1
      
      # - farms that belong to the current region
      farms.region <- farms[grepl(region, farms)]
      
      # - relevant vector to extract the region and country levels and the average of the levels of the farms for the current region
      for(i in 1:dim(smot$Cts)[3]){ 
        vector.country.region.farm <- vector.country.region
        farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
        vector.country.region.farm[farms.region.i] <- 1/length(farms.region.i)
        
        # - extract complete region level
        r <- vector.country.region.farm %*% smot$Cts[,,i] %*% vector.country.region.farm 
        df.Cts[i, region] <- r
      }
    }
    
    ## for country 
    ##(sum of the country and average of the region levels and the average of the levels of the farms divided by the number of regions)

    # - relevant vector to extract the country and average of the region levels
    vector.sum <- vector.country
    vector.sum[regions] <- 1/length(regions)
    
    # - relevant vector to extract the country and average of the region levels, and the average of the levels of the farms divided by the number of regions
    for(i in 1:dim(smot$Cts)[3]){ 
      vector <- vector.sum
      for(region in regions){
        # - farms that belong to the current region
        farms.region <- farms[grepl(region, farms)]
        farms.region.i <- farms.region[!is.na(D.nseq[i,farms.region])] #farms that have observations on that time step
        vector[farms.region.i] <- 1/length(farms.region.i)/length(regions)
      }
      
      # - extract complete country level
      c <- vector %*% smot$Cts[,,i] %*% vector
      df.Cts[i, "Country"] <- c
    }
    results.list[["Cts"]] <- df.Cts
  }
  return(results.list)
}


################################################################################################################
### Estimate discount factors by standard optimization  ###

# Function to be minimized ----
# - parms is a vector with logarithmic transforms of the observation variances and the logit transforms of the discount factors
# - vars is the number of observation variances and discount factors we want, in this case 3 because we want one per country, one per region and one per farm
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - D.nseq is a dataset with information about the number of the production cycle (ID) - each column should be a farm and each row a date
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - H.w.country is the number of harmonic waves we want to model at country level
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector containing the name of the regions
## We wish to find the observation variances and the discount factors that minimize the RMSE
getFit = function(parms, vars=3, D, D.months, D.nseq, mu0, C0, H.w.country, Spline.list, regions) {
  
  # Transform parameters back
  # - Variances 
  countryVar = exp(parms[1])
  regionVar = exp(parms[2])
  farmVar = exp(parms[3])
  
  # - Discount factors  
  deltas = c()
  for (i in (vars+1):length(parms)) {
    delta = 1.0/(1 + exp(-parms[i]))
    deltas = c(deltas, delta)
  }
  
  # Call the filtering function
  resDis = runDiscountDLM(D, D.months, mu0, C0, countryVar, regionVar, farmVar, deltas, H.w.country, Spline.list, regions)
  
  # Calculate the RMSE
  results <- extract.res(resDis, smot=NULL, D, D.nseq, H.w.country, regions)
  et.results <- results[["et"]]
  et.all <- na.omit(unlist(et.results[,c(2:ncol(et.results))]))
  RMSE <- round(sqrt(mean(na.omit(et.all)^2)),4)
  
  # The next line can be enabled if you wish to follow the values 
  # during optimization (but it delays the optimization)
  print(paste("Parms:", countryVar, regionVar, farmVar, toString(deltas), RMSE))
  
  return(RMSE)
}



# Function to transform the V's and deltas to log and logit, respectively ----
# - priorVs is a vector of observation variances (country, region and farm)
# - priorDeltas is a vector of discount factors (country, region and farm)
createPriorTransform = function(priorVs, priorDeltas) {
  transPrior = c()
  vars = length(priorVs)
  N = vars + length(priorDeltas)
  transPrior[1:vars] = log(priorVs[1:vars])
  transPrior[(vars+1):N] = log(priorDeltas/(1 - priorDeltas))
  return(transPrior)
}



# Function to transform the estimated values back to normal ----
# - optimRes is a vector with the optimized observation variances and discount factors with log and logit transformations, respectively
# - vars is the number of observation variances and discount factors we want, in this case 3 because we want one per country, one per region and one per farm
transformResults = function(optimRes, vars) {
  pars = c()
  N = length(optimRes$par)
  pars[1:vars] = exp(optimRes$par[1:vars])
  pars[(vars+1):N] = 1.0/(1 + exp(-optimRes$par[(vars+1):N]))
  return(pars)
}



# Function to estimate the discount factors and observation variances ----
# - priorVs is a vector of observation variances (country, region and farm)
# - priorDeltas is a vector of discount factors (country, region and farm)
# - D is the dataset with the mortality data - each column should be a farm and each row a date
# - D.months is a dataset with information about the months since start for each production cycle - each column should be a farm and each row a date
# - D.nseq is a dataset with information about the number of the production cycle (ID) - each column should be a farm and each row a date
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - H.w.country is the number of harmonic waves we want to model at country level
# - Spline.list is a list with the splines calculated for each farm
# - regions is a vector containing the name of the regions
estimateDiscountModel = function(priorVs, priorDeltas, D, D.months, D.nseq, mu0, C0, H.w.country, Spline.list, regions) {
  
  # Transform the observation variances and the discount factors with log and logit transformations, respectively
  initialParms = createPriorTransform(priorVs, priorDeltas)
  
  # Estimated values
  estim = optim(initialParms, getFit, gr = NULL, vars=3, D, D.months, D.nseq, mu0, C0, H.w.country, Spline.list, regions)
  
  # The "optim object" is returned. The parameters can be extracted by transformResults(optimRes, vars)
  return(estim)
}
