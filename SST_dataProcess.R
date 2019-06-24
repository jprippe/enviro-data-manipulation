# Adapted from code by Fernando Lima (CIBIO/InBIO, Universidado do Porto)
# Last updated 19 Jun 2019

#############################################
####### Generate summary SST datasets #######
#############################################
setwd('')

# Import SST datasets
hadi <- read.csv('', header = T)
oi <- read.csv('', header = T)
mur <- read.csv('', header = T)

# Add date column (these commands may or may not be applicable depending on input data)
hadi$date <- as.Date(do.call(sprintf, c(hadi[,1:3], fmt = c("%02d-%02d-%4d"))), format = "%m-%d-%Y")
hadi <- hadi[,-c(1:3)]
oi$date <- as.POSIXct(oi$date, format = '%m/%d/%y')
oi <- oi[format(oi$date,'%Y') != "1981", ]
mur$date <- as.POSIXct(mur$date, format = '%m/%d/%y')

# Convert data to long form for downstream analyses (if not in long form already)
mur.long <- melt(mur, id = 'date', value.name = 'sst', variable.name = 'siteID')
oi.long <- melt(oi, id = 'date', value.name = 'sst', variable.name = 'siteID')
hadi.long <- melt(hadi, id = 'date', value.name = 'sst', variable.name = 'siteID')

weekly.mean.all <- weekly.max.all <- weekly.max.all <- weekly.min.all <- weekly.range.all <- 
  annual.weekly.range.all <- monthly.max.all <- monthly.min.all <- monthly.range.all <- 
  summer.mean.all <- days.27to29.all <- days.above28.all <- DHW.all <- annual.DHW.all <- bleaching.hotspots.all <- bleaching.hotspots.higher.than.one.all <- 
  annual.cumulative.stress.all <- n.stress.events.all <- duration.stress.events.all <- 
  between.stress.events.all <- n.mod.stress.events.all <- duration.mod.stress.events.all <- 
  between.mod.stress.events.all <- n.sev.stress.events.all <- duration.sev.stress.events.all <-
  between.sev.stress.events.all <- monthly.mean.all <- annual.mean.all <- annual.max.all <-
  annual.min.all <- annual.range.all <- month.clim.all <- site.statistics.all <- list()

for(x in c('hadi','oi','mur')){data <- get(x);
for(reg in c('y','g','w')){
  
  sub.data <- cbind.data.frame(data$date,data[,which(substr(names(data),1,1)==reg)])
  colnames(sub.data) <- c('date',colnames(data[which(substr(names(data),1,1)==reg)]))
  
  # Convert data to extensible time series for temporal subsetting and calculations
  xts.data <- xts(sub.data, sub.data$date)[,-1]
  
  # Temperature averages, extremes, and ranges over weekly, monthly and annual time periods
  if(x!='hadi'){
    weekly.mean <- apply.weekly(xts.data,mean)
    weekly.max <- apply.weekly(xts.data,apply,2,max)
    weekly.min <- apply.weekly(xts.data,apply,2,min)
    weekly.range <- data.frame(matrix(as.numeric(weekly.max)-as.numeric(weekly.min), ncol = ncol(weekly.max)), row.names = index(weekly.max))
    annual.weekly.range <- data.frame(apply.yearly(weekly.range,mean))
    monthly.max <- apply.monthly(xts.data,apply,2,max)
    monthly.min <- apply.monthly(xts.data,apply,2,min)
    monthly.range <- data.frame(matrix(as.numeric(monthly.max)-as.numeric(monthly.min), ncol = ncol(monthly.max)), row.names = index(monthly.max))
    days.27to29 <- data.frame(apply.yearly(xts.data,apply,2,function(x) sum(27 < x & x < 29)))
    days.above28 <- data.frame(apply.yearly(xts.data,apply,2,function(x) sum(x > 28)))
    colnames(weekly.range) <- colnames(annual.weekly.range) <- colnames(monthly.range) <- colnames(weekly.max)
    rownames(annual.weekly.range) <- rownames(days.27to29) <- rownames(days.above28) <- substr(rownames(annual.weekly.range),1,4)
  }
  monthly.mean <- apply.monthly(xts.data,mean)
  annual.mean <- data.frame(apply.yearly(xts.data,mean))
  annual.max <- data.frame(apply.yearly(xts.data,apply,2,max))
  annual.min <- data.frame(apply.yearly(xts.data,apply,2,min))
  annual.range <- data.frame(matrix(as.numeric(apply.yearly(xts.data,apply,2,max))-as.numeric(apply.yearly(xts.data,apply,2,min)), ncol = ncol(annual.max)), row.names = substr(rownames(annual.max),1,4))
  summer.mean <- data.frame(apply.yearly(xts.data[.indexmon(xts.data) %in% c(6,7,8)], mean))
  colnames(annual.range) <- colnames(annual.max)
  rownames(annual.mean) <- rownames(annual.max) <- rownames(annual.min) <- rownames(summer.mean) <- substr(rownames(annual.mean),1,4)
  
  # Monthly climatologies
  month.clim <- data.frame(aggregate(monthly.mean,factor(substr(index(monthly.mean),6,7)),mean,na.rm=T))
  
  # Hottest and coldest month averages
  mean.hottest.month <- apply(month.clim,2,max)
  mean.coldest.month <- apply(month.clim,2,min)
  which.hottest.month <- apply(month.clim,2,which.max)
  which.coldest.month <- apply(month.clim,2,which.min)
  regional.hottest.month <- max(mean.hottest.month,na.rm=T)
  regional.coldest.month <- min(mean.coldest.month,na.rm=T)
  
  if(x!='hadi'){
    # Bleaching hotspots from http://coralreefwatch.noaa.gov/satellite/methodology/methodology.php#hotspot
    bleaching.hotspots <- data.frame(matrix(as.numeric(xts.data)-as.numeric(mean.hottest.month), ncol = ncol(xts.data)), row.names = index(xts.data))
    bleaching.hotspots.higher.than.one <- bleaching.hotspots
    bleaching.hotspots.higher.than.one[bleaching.hotspots.higher.than.one<1] <- NA
    colnames(bleaching.hotspots) <- colnames(bleaching.hotspots.higher.than.one) <- names(xts.data)
    annual.cumulative.stress <- aggregate(bleaching.hotspots.higher.than.one, list(substr(rownames(bleaching.hotspots.higher.than.one),1,4)), sum, na.rm=T)
    colnames(annual.cumulative.stress)[1] <- 'Year'
    
    # Degree heating weeks from http://coralreefwatch.noaa.gov/satellite/methodology/methodology.php#dhw
    DHW <- data.frame(rollapply(bleaching.hotspots.higher.than.one, 84, align="right", sum, na.rm=T)/7) # accumulation of hotspots for 84 day (12 week) intervals, values divided by 7 because we have 7 daily values in a week
    rownames(DHW) <- rownames(bleaching.hotspots)[-c(1:83)]
    annual.DHW <- apply.yearly(DHW, mean)[-1,]
    
    # Stress frequency from http://coralreefwatch.noaa.gov/satellite/thermal_history/th_stress_frequency.php
    # Number of thermal stress events (Hotspot > 1 degC)
    count.sep.events <- function(y) {sum(rle(!is.na(y))$values)}
    stress.events <- bleaching.hotspots.higher.than.one
    n.stress.events <- aggregate(stress.events, list(substr(rownames(stress.events),1,4)), function(y) count.sep.events(y))
    # Duration of stress events (Hotspot > 1 degC)
    stress.duration <- function(y) {
      temp <- data.frame(values=rle(!is.na(y))$values,length=rle(!is.na(y))$length)
      avg.duration <- mean(temp[temp$values == TRUE,]$length)
      return(avg.duration)
    }
    duration.stress.events <- aggregate(stress.events, list(substr(rownames(stress.events),1,4)), function(y) stress.duration(y))
    # Average period between stress events (Hotspot > 1 degC) starting with first event and ending with last
    # Years with only one stress event return NA
    no.stress.duration <- function(y) {
      if(sum(rle(!is.na(y))$values) <= 1){
        return(NA)
      } else {
        first.event <- min(which(!is.na(y))); last.event <- max(which(!is.na(y)))
        first.to.last <- y[c(first.event:last.event)]
        temp <- data.frame(values=rle(!is.na(first.to.last))$values,length=rle(!is.na(first.to.last))$length)
        period.between <- mean(temp[temp$values == FALSE,]$length)
        return(period.between)
      }
    }
    between.stress.events <- aggregate(stress.events, list(substr(rownames(stress.events),1,4)), function(y) no.stress.duration(y))
    colnames(n.stress.events)[1] <- colnames(duration.stress.events)[1] <- colnames(between.stress.events)[1] <- 'Year'
    
    # Number of moderate thermal stress events (DHW > 4)
    count.mod.events <- function(y) {sum(rle(y>=4)$values)}
    n.mod.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) count.mod.events(y))
    # Duration of stress events (DHW > 4)
    mod.stress.duration <- function(y) {
      temp <- data.frame(values=rle(y>=4)$values,length=rle(y>=4)$length)
      avg.duration <- mean(temp[temp$values == TRUE,]$length)
      return(avg.duration)
    }
    duration.mod.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) mod.stress.duration(y))
    # Average period between stress events (DHW > 4) starting with first event and ending with last
    # Years with only one stress event return NA
    no.mod.stress.duration <- function(y) {
      if(sum(rle(y>=4)$values) <= 1){
        return(NA)
      } else {
        first.event <- min(which(y>=4)); last.event <- max(which(y>=4))
        first.to.last <- y[c(first.event:last.event)]
        temp <- data.frame(values=rle(first.to.last>=4)$values,length=rle(first.to.last>=4 )$length)
        period.between <- mean(temp[temp$values == FALSE,]$length)
        return(period.between)
      }
    }
    between.mod.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) no.mod.stress.duration(y))
    colnames(n.mod.stress.events)[1] <- colnames(duration.mod.stress.events)[1] <- colnames(between.mod.stress.events)[1] <- 'Year'
    
    # Number of severe thermal stress events (DHW > 8)
    count.sev.events <- function(y) {sum(rle(y>=8)$values)}
    n.sev.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) count.sev.events(y))
    # Duration of stress events (DHW > 8)
    sev.stress.duration <- function(y) {
      temp <- data.frame(values=rle(y>=8)$values,length=rle(y>=8)$length)
      avg.duration <- mean(temp[temp$values == TRUE,]$length)
      return(avg.duration)
    }
    duration.sev.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) sev.stress.duration(y))
    # Average period between stress events (DHW > 8) starting with first event and ending with last
    # Years with only one stress event return NA
    no.sev.stress.duration <- function(y) {
      if(sum(rle(y>=8)$values) <= 1){
        return(NA)
      } else {
        first.event <- min(which(y>=8)); last.event <- max(which(y>=8))
        first.to.last <- y[c(first.event:last.event)]
        temp <- data.frame(values=rle(first.to.last>=8)$values,length=rle(first.to.last>=8 )$length)
        period.between <- mean(temp[temp$values == FALSE,]$length)
        return(period.between)
      }
    }
    between.sev.stress.events <- aggregate(DHW, list(substr(rownames(DHW),1,4)), function(y) no.sev.stress.duration(y))
    colnames(n.sev.stress.events)[1] <- colnames(duration.sev.stress.events)[1] <- colnames(between.sev.stress.events)[1] <- 'Year'
  }
  
  ###############################
  #Bleaching Alert Level
  
  #Stress Level	   	          Definition	   	                Effect
  #No Stress                  HotSpot <= 0
  #Bleaching Watch            0 < HotSpot < 1
  #Bleaching Warning          1 <= HotSpot and 0 < DHW < 4    Possible Bleaching
  #Bleaching Alert Level 1    1 <= HotSpot and 4 <= DHW < 8   Bleaching Likely
  #Bleaching Alert Level 2    1 <= HotSpot and 8 <= DHW       Mortality Likely
  
  ###############################
  
  # Warming rate
  # function to find the significance given the correlation coefficient (R) and the number of observations (N)
  r.to.sign<-function(R,N){T<-R/(sqrt((1-(R^2))/(N-2))); return((1-pt(T,df=N-2))*2)}
  # function to calculate the effective size of an autocorrelated data series (quenouille correction)
  eff.sample<-function(corr.coef,N){return(round(N*((1-corr.coef)/(1+corr.coef)),0))}
  
  slope.list <- vector('list')
  for(i in names(xts.data)){
    if(sum(!is.na(sub.data[,i])) < 3){
      slope.temp <- prob <- AC.coef <- AC.sign <- NA
    } else {
      lm.temp <- lm(get(i) ~ date, sub.data)
      if(x!='hadi'){
        slope.temp <- lm.temp$coefficients[2]*3600*24*365.25*10 #ºC/decade (POSIXct dates stored numerically in units of seconds)
      } else {slope.temp <- lm.temp$coefficients[2]*365.25*10} #ºC/decade
      prob <- sqrt(summary(lm.temp)[[8]])
      # Autocorrelation
      AC.coef <- cor(lm.temp$residuals[1:(length(lm.temp$residuals)-1)],lm.temp$residuals[2:length(lm.temp$residuals)])
      # Test if autocorrelation is significant
      AC.sign <- r.to.sign(AC.coef,length(lm.temp$residuals)-1)
      # and adjusted correlation probability
      if(AC.sign < 0.05){
        eff.n <- eff.sample(AC.coef,nrow(sub.data))
        prob <- r.to.sign(prob,eff.n)
      }
    }
    vec <- as.vector(c(slope.temp,prob,AC.coef,AC.sign))
    slope.list[[i]] <- vec
  }
  slope.df <- do.call('rbind', slope.list)
  site.statistics <- cbind(mean.hottest.month, which.hottest.month, mean.coldest.month, which.coldest.month, slope.df)
  colnames(site.statistics) <- c('mean.hot.month', 'hot.month', 'mean.cold.month', 'cold.month', 'slope','prob','AC.coef','AC.sign')
  
  if(x!='hadi'){
    assign(paste0('weekly.mean.',reg), weekly.mean)
    assign(paste0('weekly.max.',reg), data.frame(weekly.max))
    assign(paste0('weekly.min.',reg), data.frame(weekly.min))
    assign(paste0('weekly.range.',reg), weekly.range)
    assign(paste0('annual.weekly.range.',reg), annual.weekly.range)
    assign(paste0('days.27to29.',reg), days.27to29)
    assign(paste0('days.above28.',reg), days.above28)
    assign(paste0('monthly.max.',reg), data.frame(monthly.max))
    assign(paste0('monthly.min.',reg), data.frame(monthly.min))
    assign(paste0('monthly.range.',reg), monthly.range)
    assign(paste0('DHW.',reg), DHW)
    assign(paste0('annual.DHW.',reg), annual.DHW)
    assign(paste0('bleaching.hotspots.',reg), bleaching.hotspots)
    assign(paste0('bleaching.hotspots.higher.than.one.',reg), bleaching.hotspots.higher.than.one)
    assign(paste0('annual.cumulative.stress.',reg), annual.cumulative.stress)
    assign(paste0('n.stress.events.',reg), n.stress.events)
    assign(paste0('duration.stress.events.',reg), duration.stress.events)
    assign(paste0('between.stress.events.',reg), between.stress.events)
    assign(paste0('n.mod.stress.events.',reg), n.mod.stress.events)
    assign(paste0('duration.mod.stress.events.',reg), duration.mod.stress.events)
    assign(paste0('between.mod.stress.events.',reg), between.mod.stress.events)
    assign(paste0('n.sev.stress.events.',reg), n.sev.stress.events)
    assign(paste0('duration.sev.stress.events.',reg), duration.sev.stress.events)
    assign(paste0('between.sev.stress.events.',reg), between.sev.stress.events)
  }
  
  assign(paste0('monthly.mean.',reg), monthly.mean)
  assign(paste0('annual.mean.',reg), annual.mean)
  assign(paste0('annual.max.',reg), data.frame(annual.max))
  assign(paste0('annual.min.',reg), data.frame(annual.min))
  assign(paste0('annual.range.',reg), annual.range)
  assign(paste0('summer.mean.',reg), summer.mean)
  assign(paste0('month.clim.',reg), month.clim)
  assign(paste0('site.statistics.',reg), site.statistics)
}

if(x!='hadi'){
  weekly.mean.all[[x]] <- cbind(weekly.mean.y, weekly.mean.g, weekly.mean.w)
  weekly.max.all[[x]] <- cbind(weekly.max.y, weekly.max.g, weekly.max.w)
  weekly.min.all[[x]] <- cbind(weekly.min.y, weekly.min.g, weekly.min.w)
  weekly.range.all[[x]] <- cbind(weekly.range.y, weekly.range.g, weekly.range.w)
  annual.weekly.range.all[[x]] <- cbind(annual.weekly.range.y, annual.weekly.range.g, annual.weekly.range.w)
  days.27to29.all[[x]] <- cbind(days.27to29.y, days.27to29.g, days.27to29.w)
  days.above28.all[[x]] <- cbind(days.above28.y, days.above28.g, days.above28.w)
  monthly.max.all[[x]] <- cbind(monthly.max.y, monthly.max.g, monthly.max.w)
  monthly.min.all[[x]] <- cbind(monthly.min.y, monthly.min.g, monthly.min.w)
  monthly.range.all[[x]] <- cbind(monthly.range.y, monthly.range.g, monthly.range.w)
  DHW.all[[x]] <- cbind(DHW.y, DHW.g, DHW.w)
  annual.DHW.all[[x]] <- cbind(annual.DHW.y, annual.DHW.g, annual.DHW.w)
  bleaching.hotspots.all[[x]] <- cbind(bleaching.hotspots.y, bleaching.hotspots.g, bleaching.hotspots.w)
  bleaching.hotspots.higher.than.one.all[[x]] <- cbind(bleaching.hotspots.higher.than.one.y, bleaching.hotspots.higher.than.one.g, bleaching.hotspots.higher.than.one.w)
  annual.cumulative.stress.all[[x]] <- cbind(annual.cumulative.stress.y, annual.cumulative.stress.g[,-1], annual.cumulative.stress.w[,-1])
  n.stress.events.all[[x]] <- cbind(n.stress.events.y, n.stress.events.g[,-1], n.stress.events.w[,-1])
  duration.stress.events.all[[x]] <- cbind(duration.stress.events.y, duration.stress.events.g[,-1], duration.stress.events.w[,-1])
  between.stress.events.all[[x]] <- cbind(between.stress.events.y, between.stress.events.g[,-1], between.stress.events.w[,-1])
  n.mod.stress.events.all[[x]] <- cbind(n.mod.stress.events.y, n.mod.stress.events.g[,-1], n.mod.stress.events.w[,-1])
  duration.mod.stress.events.all[[x]] <- cbind(duration.mod.stress.events.y, duration.mod.stress.events.g[,-1], duration.mod.stress.events.w[,-1])
  between.mod.stress.events.all[[x]] <- cbind(between.mod.stress.events.y, between.mod.stress.events.g[,-1], between.mod.stress.events.w[,-1])
  n.sev.stress.events.all[[x]] <- cbind(n.sev.stress.events.y, n.sev.stress.events.g[,-1], n.sev.stress.events.w[,-1])
  duration.sev.stress.events.all[[x]] <- cbind(duration.sev.stress.events.y, duration.sev.stress.events.g[,-1], duration.sev.stress.events.w[,-1])
  between.sev.stress.events.all[[x]] <- cbind(between.sev.stress.events.y, between.sev.stress.events.g[,-1], between.sev.stress.events.w[,-1])
  
  write.zoo(weekly.mean.all[[x]], paste0('weekly.mean_',x,'.csv'), sep = ',')
  write.csv(weekly.max.all[[x]], paste0('weekly.max_',x,'.csv'))
  write.csv(weekly.min.all[[x]], paste0('weekly.min_',x,'.csv'))
  write.csv(weekly.range.all[[x]], paste0('weekly.range_',x,'.csv'))
  write.csv(annual.weekly.range.all[[x]], paste0('annual.weekly.range_',x,'.csv'))
  write.csv(days.27to29.all[[x]], paste0('days.27to29_',x,'.csv'))
  write.csv(days.above28.all[[x]], paste0('days.above28_',x,'.csv'))
  write.csv(monthly.max.all[[x]], paste0('monthly.max_',x,'.csv'))
  write.csv(monthly.min.all[[x]], paste0('monthly.min_',x,'.csv'))
  write.csv(monthly.range.all[[x]], paste0('monthly.range_',x,'.csv'))
  write.csv(DHW.all[[x]], paste0('DHW_',x,'.csv'))
  write.csv(annual.DHW.all[[x]], paste0('annual.DHW_',x,'.csv'))
  write.csv(bleaching.hotspots.all[[x]], paste0('bleaching.hotspots_',x,'.csv'))
  write.csv(bleaching.hotspots.higher.than.one.all[[x]], paste0('bleaching.hotspots.higher.than.one_',x,'.csv'))
  write.csv(annual.cumulative.stress.all[[x]], paste0('annual.cumulative.stress_',x,'.csv'), row.names = F)
  write.csv(n.stress.events.all[[x]], paste0('n.stress.events_',x,'.csv'), row.names = F)
  write.csv(duration.stress.events.all[[x]], paste0('duration.stress.events_',x,'.csv'), row.names = F)
  write.csv(between.stress.events.all[[x]], paste0('between.stress.events_',x,'.csv'), row.names = F)
  write.csv(n.mod.stress.events.all[[x]], paste0('n.mod.stress.events_',x,'.csv'), row.names = F)
  write.csv(duration.mod.stress.events.all[[x]], paste0('duration.mod.stress.events_',x,'.csv'), row.names = F)
  write.csv(between.mod.stress.events.all[[x]], paste0('between.mod.stress.events_',x,'.csv'), row.names = F)
  write.csv(n.sev.stress.events.all[[x]], paste0('n.sev.stress.events_',x,'.csv'), row.names = F)
  write.csv(duration.sev.stress.events.all[[x]], paste0('duration.sev.stress.events_',x,'.csv'), row.names = F)
  write.csv(between.sev.stress.events.all[[x]], paste0('between.sev.stress.events_',x,'.csv'), row.names = F)
}

monthly.mean.all[[x]] <- cbind(monthly.mean.y, monthly.mean.g, monthly.mean.w)
annual.mean.all[[x]] <- cbind(annual.mean.y, annual.mean.g, annual.mean.w)
annual.max.all[[x]] <- cbind(annual.max.y, annual.max.g, annual.max.w)
annual.min.all[[x]] <- cbind(annual.min.y, annual.min.g, annual.min.w)
annual.range.all[[x]] <- cbind(annual.range.y, annual.range.g, annual.range.w)
summer.mean.all[[x]] <- cbind(summer.mean.y, summer.mean.g, summer.mean.w)
month.clim.all[[x]] <- cbind(month.clim.y, month.clim.g, month.clim.w)
site.statistics.all[[x]] <- rbind(site.statistics.y, site.statistics.g, site.statistics.w)

write.zoo(monthly.mean.all[[x]], paste0('monthly.mean_',x,'.csv'), sep = ',')
write.csv(annual.mean.all[[x]], paste0('annual.mean_',x,'.csv'))
write.csv(annual.max.all[[x]], paste0('annual.max_',x,'.csv'))
write.csv(annual.min.all[[x]], paste0('annual.min_',x,'.csv'))
write.csv(annual.range.all[[x]], paste0('annual.range_',x,'.csv'))
write.csv(summer.mean.all[[x]], paste0('summer.mean_',x,'.csv'))
write.csv(month.clim.all[[x]], paste0('month.clim_',x,'.csv'))
write.csv(site.statistics.all[[x]], paste0('site.statistics_',x,'.csv'))
}

# Create list of lists containing annual summaries for all three SST datasets
move.Year <- function(y){rownames(y) <- y[['Year']]; y <- y[,-1]}
annual.cumulative.stress.all <- lapply(annual.cumulative.stress.all, function(y) move.Year(y))
n.stress.events.all <- lapply(n.stress.events.all, function(y) move.Year(y))
duration.stress.events.all <- lapply(duration.stress.events.all, function(y) move.Year(y))
between.stress.events.all <- lapply(between.stress.events.all, function(y) move.Year(y))

annual.data.lol <- list(annual.mean.all, annual.max.all, annual.min.all, annual.range.all, annual.weekly.range.all, summer.mean.all, days.27to29.all, days.above28.all, annual.cumulative.stress.all, n.stress.events.all, duration.stress.events.all, between.stress.events.all)
names(annual.data.lol) <- c('annual.mean.all', 'annual.max.all', 'annual.min.all', 'annual.range.all', 'annual.weekly.range.all', 'summer.mean.all', 'days.27to29.all', 'days.above28.all', 'annual.cumulative.stress.all', 'n.stress.events.all', 'duration.stress.events.all', 'between.stress.events.all')
for(z in 1:length(annual.data.lol)){
  melt.hack <- function(y) {
    df <- data.frame(rep(rownames(y),ncol(y)), unlist(y), rep(colnames(y), rep(nrow(y),ncol(y))))
    colnames(df) <- c('Year', names(annual.data.lol)[z], 'site'); return(df)}
  annual.data.lol[[z]] <- lapply(annual.data.lol[[z]], function(y) melt.hack(y))
}