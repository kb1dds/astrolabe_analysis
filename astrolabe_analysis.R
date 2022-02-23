# Astrolabe data exploration
#
# Copyright (c) 2022 Michael Robinson
# Available under the MIT license

library(tidyverse)
library(lubridate)
library(modelr)

#### Convenience functions

rad2ce <- function(lat,ra,dec,st){
  # Convert latitude (deg),right ascension (hr), declination (deg), 
  # and local sidereal time (rad) to elevation (radians)
  latrd <- lat*pi/180.0
  rard <- ra*2*pi/24
  decrd <- dec*pi/180.0

  lha <- st-rard
  
  num <- sin(latrd)*sin(decrd)+cos(latrd)*cos(decrd)*cos(lha)
  den <- sqrt(1-num^2)
  return(atan2(num,den))
}

radc2a <- function(lat,ra,dec,st){
  # Convert latitude (deg),right ascension (hr), declination (deg), 
  # local sidereal time (rad) to azimuth (radians)
  
  latrd <- lat*pi/180.0
  rard <- ra*2*pi/24
  decrd <- dec*pi/180.0
  
  lha <- st-rard

  el <- rad2ce(lat,ra,dec,st)
	
  num <- -sin(lha)*cos(decrd)/cos(el)
  den <- (sin(decrd)-sin(el)*sin(latrd))/(cos(el)*cos(latrd))
  if( atan2(num,den) < 0){
    return(atan2(num,den)+2*pi)
  }
  return(atan2(num,den))
}

utc2s <- function(lon,date){
  # Convert Longitude (degrees), Julian date past 1 Jan 2000,
  # Source https://en.wikipedia.org/wiki/Sidereal_time

  n <- as.duration(ymd('2000-01-01') %--% date)/ddays(1)
  
  UTC2S <- (0.77905727+1.0027378*(n-0.5))+lon/360.0
  return((UTC2S-as.integer(UTC2S))*2*pi)
}

sun_ra <- function(date){
  # Compute right ascension (in radians) of the sun from its day n past 1 January 2000.
  #  Source https://en.wikipedia.org/wiki/Position_of_the_Sun
  n <- as.duration(ymd('2000-01-01') %--% date)/ddays(1)
  L <- 280.460+0.9856003*n
  g <- 357.528+0.9856003*n
  lam <- L + 1.915*sin(g*pi/180)+0.02*sin(2*pi/180*g)
  eps <- 23.439+0.00000004*n
  ra <- atan2(cos(eps*pi/180)*sin(lam*pi/180),cos(lam*pi/180))
  return(if_else(ra<0,(ra+2*pi),ra)*12/pi)
}

sundc <- function(n){
  # Compute declination of the sun (in degrees) of the sun from its day
  # past 1 Jan 2000
  # Source https://en.wikipedia.org/wiki/Position_of_the_Sun

  l <- 280.460+0.9856003*n
  g <- 357.528+0.9856003*n
  lam <- l+1.915*sin(g*pi/180)+0.02*sin(2*pi/180*g)
  eps <- 23.439+0.00000004*n

  num <- sin(eps*pi/180)*sin(lam*pi/180)
  den <- sqrt(1-num^2)
  return(atan2(num,den)*180/pi)
}

#### Actual code begins!

data <- read_csv('~/astrolabe_analysis/observations.csv') %>%
  mutate(`UTC date`=mdy(`UTC date`,tz='UTC'),
         date=`UTC date`+`UTC time`,
         `Local true apparent solar time`=hms::as_hms(with_tz(date,tzone='EST')),
         `Sun right ascension`=sun_ra(date),
         `Jupiter right ascension`=hour(`Jupiter right ascension`)+
           minute(`Jupiter right ascension`)/60,
         `Saturn right ascension`=hour(`Saturn right ascension`)+
           minute(`Saturn right ascension`)/60,
         `Venus right ascension`=hour(`Venus right ascension`)+
           minute(`Venus right ascension`)/60,
         `Mars right ascension`=hour(`Mars right ascension`)+
           minute(`Mars right ascension`)/60,
         `Moon right ascension`=hour(`Moon right ascension`)+
           minute(`Moon right ascension`)/60
         )

#### Time measurement

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local apparent solar time`,tzone='EST')) %>%
  ggplot(aes(x=date,y=time_diff/60,color=Instrument)) + 
  geom_point()+
  ylim(-120,120) + 
  ylab('Error (minutes)') +
  theme(legend.position = 'top')

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')) %>%
  ggplot(aes(x=date,y=time_diff/60,color=Instrument)) + 
  geom_point() +
  ylim(-120,120) + 
  ylab('Error (minutes)') +
  theme(legend.position = 'top')

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST'),
         daytime=is.na(`Sun elevation`)) %>%
  ggplot(aes(x=date,y=time_diff/60,color=daytime)) + 
  geom_point() +
  ylim(-120,120) + 
  ylab('Error (minutes)') +
  theme(legend.position = 'top')

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')) %>%
  ggplot(aes(y=time_diff/60)) + 
  geom_histogram(bins=100) +
  ylim(-120,120) + 
  ylab('Error (minutes)')

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')/60) %>%
  summarise(mean(time_diff,na.rm=TRUE))

# Q-normal plot of errors; looks fairly normal to me!
data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')) %>%
  with(qqnorm(time_diff/60,ylim=c(-120,120)))

# Time estimation error, split by solar/stellar
solar_time_error_anova <- data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST'),
         daytime=!is.na(`Sun elevation`)) %>%
  filter(!is.na(time_diff)) %>%
  with(aov(daytime~time_diff))

summary(solar_time_error_anova) # Not significant!

data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')/60,
         daytime=!is.na(`Sun elevation`)) %>%
  ggplot(aes(x=time_diff,y=daytime)) + 
  geom_boxplot() +
  xlab('Error (minutes)') +
  ylab('Observation of the sun?') +
  scale_x_continuous(limits=c(-120,120))

# Solar time estimation error, versus sun elevation
data %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST'),
         daytime=!is.na(`Sun elevation`),
         morning=hour(`Local mean time`)<12) %>%
  filter(daytime==TRUE) %>%
  ggplot(aes(x=cut(`Sun elevation`,breaks=10),y=time_diff/60,color=morning)) +
  geom_boxplot()+
  ylim(-120,120) + 
  ylab('Error (minutes)') +
  xlab('Sun elevation (degrees)') +
  theme(axis.text.x = element_text(angle=90))

#### Pivot to handling planets and individual stars

# This takes several steps because stars don't have varying right ascensions,
# and therefore are not listed in the observation spreadsheet
# Objects in the solar system have other properties that can be observed,
# and therefore these need to be pivoted in separately
data_individual <- data %>% 
  mutate(daytime=!is.na(`Sun elevation`)) %>%
  pivot_longer(cols=ends_with('elevation'),
               names_to='object',
               values_to='elevation') %>%
  mutate(object=str_remove(object,' elevation')) %>%
  filter(!is.na(elevation)) %>%
  mutate(right_ascension=NA, # These are the properties that solar system objects
         azimuth=NA,  # have and need pivoting in separately
         phase=NA)

# All objects
objects <- data_individual$object %>% unique()

# Pivot in objects that have their own right ascension columns (solar system objects)
for( obj in objects ){
  if( paste(obj,'right ascension') %in% colnames(data_individual) ){
    data_individual$right_ascension[data_individual$object==obj] <- 
      data_individual[[paste(obj,'right ascension')]][data_individual$object==obj]
  }
  if( paste(obj,'azimuth') %in% colnames(data_individual) ){
    data_individual$azimuth[data_individual$object==obj] <- 
      data_individual[[paste(obj,'azimuth')]][data_individual$object==obj]
  }
  if( paste(obj,'phase') %in% colnames(data_individual) ){
    data_individual$phase[data_individual$object==obj] <- 
      data_individual[[paste(obj,'phase')]][data_individual$object==obj]
  }
}

# Remove solar system right ascension columns
data_individual <- data_individual %>%
  select(-ends_with(' right ascension')) %>%
  select(-ends_with(' azimuth')) %>%
  select(-ends_with(' phase')) %>%
  mutate(solar_system_object=!is.na(right_ascension)) %>%
  mutate(sun_right_ascension=sun_ra(date)) # Yes, it looks like we recomputed this.. oh well

# Set the phases in proper order
# This is a small hack because I record Moon phase as waxing or waning, 
# but I don't record wax/wane status on other objects... perhaps I should!
phases<-c(sapply(1:99,function(x){sprintf('Wax%02d',x)}),
          'Full',
          sapply(99:1,function(x){sprintf('Wane%02d',x)}),
          1:99) # Bit of a hack here...

data_individual <- data_individual %>% 
  mutate(phase=factor(phase,phases))

# Join in the fixed geocentric locations of the stars
starchart <- read_csv('~/astrolabe_analysis/stars.csv')

data_individual <- left_join(data_individual,starchart,by='object') %>% 
  mutate(right_ascension=coalesce(right_ascension.x,right_ascension.y)) %>% 
  select(-right_ascension.x,-right_ascension.y)

# What's in the data?
data_individual %>% 
  ggplot(aes(y=object,fill=daytime)) + 
  geom_bar() +
  theme(legend.position = 'top')

data_individual %>% 
  group_by(date,Instrument) %>% 
  count() %>% 
  ggplot(aes(as.factor(n),fill=Instrument)) + 
  stat_count() + 
  xlab('Sightings per observing session')

data_individual %>% 
  group_by(date,daytime) %>% 
  count() %>% 
  ggplot(aes(as.factor(n),fill=daytime)) + 
  stat_count() + 
  xlab('Sightings per observing session')

#### Elevation estimation error

# Sun elevation error
sun_elevation_data <- data_individual %>% 
  filter(object=='Sun') %>%
  mutate(jd=as.integer(as.duration(ymd('2000-01-01') %--% date)/ddays(1)),
         sidereal_time=utc2s(Longitude,date),
         declination=sundc(jd),
         expected_elevation=rad2ce(Latitude,
                                   right_ascension,
                                   declination,
                                   sidereal_time)*180/pi,
         elevation_difference=elevation-expected_elevation) 

mean(sun_elevation_data$elevation_difference,na.rm=TRUE)
sd(sun_elevation_data$elevation_difference,na.rm=TRUE)

sun_elevation_data %>% ggplot(aes(elevation_difference)) + geom_histogram(bins=100)

sun_elevation_data %>%
  with(qqnorm(elevation_difference,ylim=c(-5,5)))

sun_elevation_data %>% 
  ggplot(aes(elevation,elevation_difference)) +
  geom_point()

# Stellar elevation error
elevation_data <- 
  data_individual %>% 
  filter(!solar_system_object,
         object!="Moon", object!='Mars',object !='Saturn') %>% # Silliness due to some missing data
  mutate(sidereal_time=utc2s(Longitude+360,date),
         expected_elevation=rad2ce(Latitude,
                                   right_ascension,
                                   declination,
                                   sidereal_time)*180/pi,
         elevation_difference=elevation-expected_elevation)

sd(elevation_data$elevation_difference,na.rm=TRUE)
mean(elevation_data$elevation_difference,na.rm=TRUE)

elevation_data %>% 
  ggplot(aes(elevation_difference)) + 
  geom_histogram(binwidth = 0.2) +
  scale_x_continuous(limits = c(-10,10))

elevation_data %>%
  with(qqnorm(elevation_difference,ylim=c(-5,5)))

elevation_data %>% 
  group_by(object) %>%
  ggplot(aes(elevation_difference,object)) + 
  geom_boxplot()

elevation_data %>% 
  mutate(year=as.factor(year(date))) %>%
  group_by(year) %>%
  ggplot(aes(elevation_difference,year)) + 
  geom_boxplot()

elevation_data %>% 
  group_by(Instrument) %>%
  ggplot(aes(elevation_difference,Instrument)) + 
  geom_boxplot()

elevation_data %>% 
  group_by(Instrument) %>%
  ggplot(aes(elevation,elevation_difference)) +
  geom_point() + 
  geom_smooth(method='lm') +
  scale_y_continuous(limits=c(-5,5))

#### Planets

data_individual %>% 
  filter(solar_system_object,object!='Moon') %>%
  ggplot(aes(x=date,y=right_ascension,color=object)) +
  geom_point() +
  ylab('Right ascension (hours)') +
  theme(legend.position = 'top')

data_individual %>% 
  filter(solar_system_object,
         object!='Moon',
         object!='Sun') %>%
  mutate(right_ascension_diff=right_ascension-sun_right_ascension) %>%
  ggplot(aes(x=date,y=right_ascension_diff,color=object)) +
  geom_point() +
  geom_ref_line(h=0) +
  geom_ref_line(h=12) +
  geom_ref_line(h=-12) +
  geom_ref_line(h=-24) +
  ylab('Right ascension offset from sun (hours)') +
  theme(legend.position = 'top')

#### Determining orbital radii

### Determining Venus's orbital radius by fitting a sinusoid
venus_sinusoid <- nls(diff~a*sin(2*pi*n/b+ph)+12,
    start=list(a=3,b=600,ph=0), # Eyeballed off graph
    data=data_individual %>% 
      filter(object=='Venus') %>%
      mutate(diff=(right_ascension-sun_right_ascension-12)%%24,
             n=as.duration(ymd('2000-01-01') %--% date)/ddays(1)))

summary(venus_sinusoid)

data_individual %>% filter(object=='Venus') %>%
  mutate(diff=(right_ascension-sun_right_ascension-12)%%24,
         n=as.duration(ymd('2000-01-01') %--% date)/ddays(1)) %>%
  ggplot() +
  geom_point(aes(x=n,y=diff,color=daytime)) +
  geom_line(data=tibble(n=seq_range(c(7400,8100),n=100)) %>%
              add_predictions(venus_sinusoid),
            aes(x=n,y=pred)) +
  geom_ref_line(h=0) + 
  geom_ref_line(h=12) + 
  geom_ref_line(h=24) + 
  xlab('Julian date') +
  ylab('Venus-to-sun right ascension difference (hours)') +
  theme(legend.position = 'top')

venus_radius<-(coef(venus_sinusoid)[[2]]/(coef(venus_sinusoid)[[2]]+365))^(2/3)

# Mars orbital radius (unlikely to work because Mars orbit is eccentric)
data_individual %>% filter(object=='Mars') %>%
  mutate(diff=(right_ascension-sun_right_ascension)%%24,
         n=as.duration(ymd('2000-01-01') %--% date)/ddays(1)) %>%
  ggplot() +
  geom_point(aes(x=n,y=diff)) +
  geom_ref_line(h=12)+
  xlab('Julian date') +
  ylab('Mars-to-sun right ascension difference (hours)') +
  theme(legend.position = 'top')

mars_radius<-(800/(800-365))^(2/3)

# Determining Jupiter's orbital radius using time between oppositions. 
# I have observed three oppositions
odata<-data_individual %>% filter(object=='Jupiter') %>%
  mutate(diff = right_ascension-sun_right_ascension-12,
         year=as.factor(year(date)))

# Determining opposition by lines of best fit
# Exclude observations that are too far from opposition
lfit_2019 <- lm(diff~as_date(date),odata%>%filter(year==2019))
opp_2019<- -lfit_2019$coefficients[[1]]/lfit_2019$coefficients[[2]]

lfit_2020 <- lm(diff~as_date(date),odata%>%filter(year==2020,diff>-5))
opp_2020<- -lfit_2020$coefficients[[1]]/lfit_2020$coefficients[[2]]

lfit_2021 <- lm(diff~as_date(date),odata%>%filter(year==2021))
opp_2021<- -lfit_2021$coefficients[[1]]/lfit_2021$coefficients[[2]]

odata %>%
  ggplot(aes(x=as_date(date),y=diff,color=year)) +
  geom_point() +
  geom_ref_line(h=0) + 
  geom_abline(slope=coef(lfit_2019)[2],intercept=coef(lfit_2019)[1]) +
  geom_abline(slope=coef(lfit_2020)[2],intercept=coef(lfit_2020)[1]) +
  geom_abline(slope=coef(lfit_2021)[2],intercept=coef(lfit_2021)[1]) +
  xlab('Date') +
  ylab('Jupiter-to-sun right ascension difference (hours)') +
  theme(legend.position = 'top')

# Determine radius of Jupiter's orbit in astronomical units, using
# Newton's circular orbit period model
jupiter_radius <- ((opp_2020-opp_2019)/((opp_2020-opp_2019)-365))^(2/3)
jupiter_radius <- ((opp_2021-opp_2020)/((opp_2021-opp_2020)-365))^(2/3)
jupiter_radius <- ((opp_2021-opp_2019)/((opp_2021-opp_2019)-365*2))^(2/3)

# Determining Saturn's orbital radius using time between oppositions. 
# I have observed three oppositions
odata<-data_individual %>% filter(object=='Saturn') %>%
  mutate(diff = right_ascension-sun_right_ascension-12,
         year=as.factor(year(date)))

# Determining opposition by lines of best fit
# Exclude observations that are too far from opposition
lfit_2019 <- lm(diff~as_date(date),odata%>%filter(year==2019))
opp_2019<- -lfit_2019$coefficients[[1]]/lfit_2019$coefficients[[2]]

lfit_2020 <- lm(diff~as_date(date),odata%>%filter(year==2020))
opp_2020<- -lfit_2020$coefficients[[1]]/lfit_2020$coefficients[[2]]

lfit_2021 <- lm(diff~as_date(date),odata%>%filter(year==2021))
opp_2021<- -lfit_2021$coefficients[[1]]/lfit_2021$coefficients[[2]]

# Determine radius of Saturn's orbit in astronomical units, using
# Newton's circular orbit period model
saturn_radius <- ((opp_2020-opp_2019)/((opp_2020-opp_2019)-365))^(2/3)
saturn_radius <- ((opp_2021-opp_2020)/((opp_2021-opp_2020)-365))^(2/3)
saturn_radius <- ((opp_2021-opp_2019)/((opp_2021-opp_2019)-365*2))^(2/3)

odata %>%
  ggplot(aes(x=as_date(date),y=diff,color=year)) +
  geom_point() +
  geom_ref_line(h=0) + 
  geom_abline(slope=coef(lfit_2019)[2],intercept=coef(lfit_2019)[1]) +
  geom_abline(slope=coef(lfit_2020)[2],intercept=coef(lfit_2020)[1]) +
  geom_abline(slope=coef(lfit_2021)[2],intercept=coef(lfit_2021)[1]) +
  xlab('Date') +
  ylab('Saturn-to-sun right ascension difference (hours)') +
  theme(legend.position = 'top')

## Trying to get multiple right ascension differences

data_individual %>% 
  filter(object=='Jupiter') %>% 
  mutate(diff=right_ascension-sun_right_ascension) %>% 
  group_by(object,round(diff/2)*2) %>% 
  summarise(mindate=min(date),maxdate=max(date)) %>%
  mutate(synodic=as.duration(mindate%--%maxdate)/ddays(1)) %>%
  filter(synodic>10)%>%
  mutate(synodic=ifelse(synodic>500,synodic/2,synodic),
         radius=(synodic/(synodic-365))^(2/3))

data_individual %>% 
  filter(object=='Saturn') %>% 
  mutate(diff=right_ascension-sun_right_ascension) %>% 
  group_by(object,round(diff/2)*2) %>% 
  summarise(mindate=min(date),maxdate=max(date)) %>%
  mutate(synodic=as.duration(mindate%--%maxdate)/ddays(1)) %>%
  filter(synodic>10)%>%
  mutate(synodic=ifelse(synodic>500,synodic/2,synodic),
         radius=(synodic/(synodic-365))^(2/3))

data_individual %>% 
  filter(object=='Mars') %>% 
  mutate(diff=right_ascension-sun_right_ascension) %>% 
  group_by(object,round(diff/2)*2) %>% 
  summarise(mindate=min(date),maxdate=max(date)) %>%
  mutate(synodic=as.duration(mindate%--%maxdate)/ddays(1)) %>%
  filter(synodic>10)%>%
  mutate(#synodic=ifelse(synodic>500,synodic/2,synodic),
         radius=(synodic/(synodic-365))^(2/3))

data_individual %>% 
  filter(object=='Venus') %>% 
  mutate(diff=right_ascension-sun_right_ascension) %>% 
  group_by(object,round(diff/2)*2) %>% 
  summarise(mindate=min(date),maxdate=max(date)) %>%
  mutate(synodic=as.duration(mindate%--%maxdate)/ddays(1)) %>%
  filter(synodic>10)%>%
  mutate(#synodic=ifelse(synodic>500,synodic/2,synodic),
         radius=(synodic/(synodic+365))^(2/3))

#### Phase analysis

data_individual %>% 
  filter(object=='Moon') %>%
  mutate(diff=(right_ascension-sun_right_ascension)%%24) %>%
  group_by(phase) %>%
  ggplot(aes(y=phase,x=diff)) +
  geom_ref_line(v=0)+
  geom_ref_line(v=12)+
  geom_ref_line(v=24)+
  geom_boxplot()

data_individual %>% 
  filter(object=='Venus') %>%
  filter(!is.na(phase)) %>%
  mutate(diff=(right_ascension-sun_right_ascension)%%24) %>%
  group_by(phase) %>%
  ggplot(aes(y=phase,x=diff)) +
  geom_ref_line(v=0)+
  geom_ref_line(v=12)+
  geom_ref_line(v=24)+
  geom_point()
