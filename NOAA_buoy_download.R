require(data.table)
require(dplyr)
require(ggplot2)
require(scales)

# Specify which buoys to compare based on their identifiers on the NDBC website
buoy.names <- c('chav3','fwyf1')

# Create empty list to fill with buoy data records
site.list <- list()
for (buoy in buoy.names){
  temp.list <- list()
  #Specify the range of years to include (make sure the chosen years are actually available on the NDBC website)
  for (i in 2007:2018) {
    dat <- fread(paste('https://www.ndbc.noaa.gov/view_text_file.php?filename=', buoy, 'h', as.character(i),'.txt.gz&dir=data/historical/stdmet/', sep=''))
    temp.list[[i]] <- dat
  }
  raw.df <- do.call(rbind, temp.list)
  filt.df <- data.frame(date = as.POSIXct(paste(with(raw.df[-1,], paste(`#YY`, MM, DD, sep = '-')), with(raw.df[-1,], paste(hh, mm, sep = ':')), sep = ' '), format = '%Y-%m-%d %H:%M'),
                        ATMP = raw.df$ATMP[-1],
                        WTMP = raw.df$WTMP[-1],
                        buoy = as.factor(buoy))
  site.list[[buoy]] <- filt.df
}
all.site.df <- do.call(rbind, site.list)
all.site.df$WTMP <- as.numeric(as.character(all.site.df$WTMP)); all.site.df$ATMP <- as.numeric(as.character(all.site.df$ATMP))
# Remove any anomalous data points (change bounds based on expected temperatures)
all.site.filt <- subset(all.site.df, WTMP > 10 & WTMP < 40)
# Calculate daily climatology (average temp for each day of the year) for each location
daily.clim <- aggregate(all.site.filt$WTMP, list(substr(all.site.filt$date,6,10), all.site.filt$buoy), mean, na.rm=T)
colnames(daily.clim) <- c('date','buoy','temp')
daily.clim$date <- as.Date(paste0('2016-', daily.clim$date)) #for plotting purposes, 2016 is a leap year so Feb 29 is included
daily.clim <- daily.clim %>% group_by(buoy) %>% mutate(avgTemp = mean(temp))

ggplot(data = daily.clim, aes(col = buoy))+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = 'none', axis.title = element_blank())+
  geom_line(aes(x=date, y=temp), size=0.3)+
  geom_hline(yintercept = unique(daily.clim$avgTemp), linetype = 'dashed')+
  scale_x_date(labels = date_format('%b'))
