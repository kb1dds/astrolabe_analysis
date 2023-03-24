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
  return(atan2(num,den))
}

utc2s <- function(lon,date){
  # Convert Longitude (degrees), Julian date past 1 Jan 2000, to sidereal time (radians)
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

heliocentric <- function(NN,II,WW,AA,EE,MM){
  # Compute heliocentric position of a body orbiting the sun
  #
  # Inputs:
  # Orbital elements of planet (about the sun)
  # NN = longitude of ascending node (radians)
  # II = inclination to the ecliptic (radians)
  # WW = Argument of perihelion (radians)
  # AA = mean distance to the sun (AU)
  # EE = eccentricity (dimensionless)
  # MM = Mean anomaly (radians)
  #
  # Returns a tibble with columns: xh, yh, zh (AU)
  #
  # Reference: [https://stjarnhimlen.se/comp/ppcomp.html]

  # Computing current radius to the sun (AU)
  E <- MM+EE*sin(MM)*(1.0+EE*cos(MM))
	XV <- AA*(cos(E)-EE)
	YV <- AA*sqrt(1.0-EE*EE)*sin(E)
	V <- atan2(YV,XV)
	R <- sqrt(XV*XV+YV*YV)
	
  #  Planet's position in space (centered on the sun in AU)
	return(tibble(xh=R*(cos(NN)*cos(V+WW)-sin(NN)*sin(V+WW)*cos(II)),
	              yh=R*(sin(NN)*cos(V+WW)+cos(NN)*sin(V+WW)*cos(II)),
	              zh=R*(sin(V+WW)*sin(II))))
}

planet_ra <- function(n,xh,yh,zh){
  # Compute right ascension of a planet based on its heliocentric 
  # coordinates (xh,yh,zh) and day past 1 January 2000
  #
  # Reference: [https://stjarnhimlen.se/comp/ppcomp.html]
  
  # Orbital elements of the sun, as a geocentric object
  EES <- (0.016709-1.151E-9*n)*pi/180
  MMS <- (356.0470+0.9856002585*n)*pi/180
  WWS <- (282.9404+4.70935E-5*n)*pi/180
  
  # Geocentric position of the sun
  ES <- MMS+EES*sin(MMS)*(1.0+EES*cos(MMS))
  XVS <- cos(ES)-EES
  YVS <- sqrt(1.0-EES*EES)*sin(ES)
  VS <- atan2(YVS,XVS)	
  RS <- sqrt(XVS*XVS+YVS*YVS)
  LONSUN <- VS+WWS
  XS <- RS*cos(LONSUN)
  YS <- RS*sin(LONSUN)
  
  # Geocentric coordinates of planet
  XG <- xh+XS
  YG <- yh+YS
  ZG <- zh
  
  # Ecliptic coordinates
  ECL <- (23.4393-3.563E-7*n)*pi/180
  XE <- XG
  YE <- YG*cos(ECL)-ZG*sin(ECL)
  ZE <- YG*sin(ECL)+ZG*cos(ECL)
  
  # Final output in equatorial coordinates
  return(atan2(YE,XE)*24/6.28319)
}

planet_dc <- function(n,xh,yh,zh){
  # Compute right ascension of a planet based on its heliocentric 
  # coordinates (xh,yh,zh) and day past 1 January 2000
  #
  # Reference: [https://stjarnhimlen.se/comp/ppcomp.html]
  
  # Orbital elements of the sun, as a geocentric object
  EES <- (0.016709-1.151E-9*n)*pi/180
  MMS <- (356.0470+0.9856002585*n)*pi/180
  WWS <- (282.9404+4.70935E-5*n)*pi/180
  
  # Geocentric position of the sun
  ES <- MMS+EES*sin(MMS)*(1.0+EES*cos(MMS))
  XVS <- cos(ES)-EES
  YVS <- sqrt(1.0-EES*EES)*sin(ES)
  VS <- atan2(YVS,XVS)	
  RS <- sqrt(XVS*XVS+YVS*YVS)
  LONSUN <- VS+WWS
  XS <- RS*cos(LONSUN)
  YS <- RS*sin(LONSUN)
  
  # Geocentric coordinates of planet
  XG <- xh+XS
  YG <- yh+YS
  ZG <- zh
  
  # Ecliptic coordinates
  ECL <- (23.4393-3.563E-7*n)*pi/180
  XE <- XG
  YE <- YG*cos(ECL)-ZG*sin(ECL)
  ZE <- YG*sin(ECL)+ZG*cos(ECL)
  
  # Final output in equatorial coordinates
  return(atan2(ZE,sqrt(XE^2+YE^2))*180/pi)
}

planet_xy_from_data <- function(n,ra,r) {
  # Determine planet xh,yh heliocentric coordinates from a radius estimate (r) in AU 
  # and date past 1 Jan 2000
  
  # Orbital elements of the earth around the sun
  EES <- (0.016709-1.151E-9*n)*pi/180
  MMS <- (356.0470+0.9856002585*n)*pi/180
  WWS <- (282.9404+4.70935E-5*n)*pi/180
  
  # Heliocentric position of the earth
  ES <- MMS+EES*sin(MMS)*(1.0+EES*cos(MMS))
  XVS <- cos(ES)-EES
  YVS <- sqrt(1.0-EES*EES)*sin(ES)
  VS <- atan2(YVS,XVS)	
  RS <- sqrt(XVS*XVS+YVS*YVS)
  LONSUN <- VS+WWS
  XE <- -RS*cos(LONSUN)
  YE <- -RS*sin(LONSUN)
  
  tmp <- XE*cos(ra*pi/12)+YE*sin(ra*pi/12)
  
  t <- -tmp+sqrt(tmp^2-(1-r^2))

  return(tibble(xhp=XE+t*cos(ra*pi/12),
                yhp=YE+t*sin(ra*pi/12)))
}

#### Actual code begins!

# Set the planets in orbit order (I can't see Pluto anyway, so no need to quibble here!)
planets <- c('Mercury','Venus','Mars','Jupiter','Saturn','Uranus','Neptune')
planet_colors <- c('#555555','#ff0000','#ff9900','#00ff00') # Note: only the planets in the data are colored...
planet_shapes <- c(1,3,4,8) # Note: only the planets in the data get symbols...

# Set the phases in proper order
# This is a small hack because I record Moon phase as waxing or waning, 
# but I don't record wax/wane status on other objects... perhaps I should!
phases<-c(sapply(1:99,function(x){sprintf('Wax%02d',x)}),
          'Full',
          sapply(99:1,function(x){sprintf('Wane%02d',x)}),
          1:99) # Bit of a hack here...

data <- read_csv('observations.csv') %>%
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
  mutate(daytime=factor(if_else(daytime,'Day','Night'),c('Day','Night'))) %>%
  ggplot(aes(x=date,y=time_diff/60,color=daytime,shape=daytime)) + 
  geom_point() +
  scale_color_discrete() +
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

data_individual <- data_individual %>% 
  mutate(phase=factor(phase,phases))

# Join in the fixed geocentric locations of the stars
starchart <- read_csv('~/astrolabe_analysis/stars.csv')

data_individual <- left_join(data_individual,starchart,by='object') %>% 
  mutate(right_ascension=coalesce(right_ascension.x,right_ascension.y)) %>% 
  select(-right_ascension.x,-right_ascension.y)

# What's in the data?
data_individual %>% 
  mutate(daytime=factor(if_else(daytime,'Day','Night'),c('Night','Day'))) %>%
  ggplot(aes(y=object,fill=daytime)) + 
  geom_bar() +
  scale_fill_grey() +
  theme(legend.position = 'top')

data_individual %>% 
  group_by(date,Instrument) %>% 
  count() %>% 
  ggplot(aes(as.factor(n),fill=Instrument)) + 
  stat_count() + 
  xlab('Sightings per observing session')

data_individual %>% 
  group_by(date,daytime) %>% 
  mutate(daytime=factor(if_else(daytime,'Day','Night'),c('Night','Day'))) %>%
  count() %>% 
  ggplot(aes(as.factor(n),fill=daytime)) + 
  stat_count() + 
  scale_fill_grey() +
  xlab('Sightings per observing session')

data_individual %>% 
  filter(!is.na(right_ascension),object!='Moon') %>%
  mutate(daytime=factor(if_else(daytime,'Day','Night'),c('Night','Day'))) %>%
  group_by(date,daytime) %>% 
  count() %>% 
  ggplot(aes(as.factor(n),fill=daytime)) + 
  stat_count() + 
  scale_fill_grey() +
  xlab('Sightings per observing session')

# Time estimation accuracy versus number of sightings
data_individual %>% 
  mutate(time_diff = interval(`Local true apparent solar time`,
                              `Local mean time`,tzone='EST')/60) %>%
  filter(!is.na(right_ascension),object!='Moon') %>%
  group_by(date,daytime) %>% 
  add_count() %>%
  ggplot(aes(as.factor(n),time_diff)) + 
  geom_boxplot() +
  xlab('Sightings per observing session') +
  ylab('Time estimation error (minutes)') +
  scale_y_continuous(limits=c(-120,120))

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

sun_elevation_data %>% ggplot(aes(elevation_difference)) + 
  geom_histogram(binwidth = 0.2) +
  scale_x_continuous(limits = c(-10,10))

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
  ggplot(aes(elevation,elevation_difference)) +
  geom_point() + 
  geom_smooth(method='lm') +
  scale_y_continuous(limits=c(-5,5))

elevation_data %>%
  group_by(Instrument) %>% 
  summarize(mean=mean(elevation_difference,na.rm=TRUE),
            sd=sd(elevation_difference,na.rm=TRUE),
            n=n())

# Difference between sun/star elevation error?
elevation_anova <- with(bind_rows(mutate(elevation_data,
             solar=FALSE),
        mutate(sun_elevation_data,
             solar=TRUE)),
     aov(elevation_difference~solar))

summary(elevation_anova)

#### Occlusion analysis, basically list all observed az/el pairs

data_individual_azel <- data_individual %>% 
  mutate(declination=if_else(solar_system_object,0,declination),
         sidereal_time=utc2s(Longitude+360,date),
         expected_elevation=rad2ce(Latitude,
                                   right_ascension,
                                   declination,
                                   sidereal_time)*180/pi,
         expected_azimuth=radc2a(Latitude,
                               right_ascension,
                               declination,
                               sidereal_time)*180/pi)

data_individual_azel %>%
  ggplot(aes(x=expected_azimuth,
             y=expected_elevation,
             color=object)) + 
  geom_point()

#### Planets

data_individual %>% 
  filter(solar_system_object,object!='Moon',object != 'Sun') %>%
  mutate(object=ordered(object,planets)) %>%
  ggplot(aes(x=date,y=right_ascension,color=object)) +
  geom_point() +
  scale_color_manual(values=planet_colors)+
  geom_line(aes(x=date,y=right_ascension),
            data = data_individual %>% 
              filter(object=='Sun') %>% 
              filter(date<mdy('3/22/2019')),
            color='black') +
  geom_line(aes(x=date,y=right_ascension),
            data = data_individual %>% 
              filter(object=='Sun') %>% 
              filter(date>mdy('3/22/2019'),date<mdy('3/22/2020')),
            color='black') +
  geom_line(aes(x=date,y=right_ascension),
            data = data_individual %>% 
              filter(object=='Sun') %>% 
              filter(date>mdy('3/22/2020'),date<mdy('3/20/2021')),
            color='black') +
  geom_line(aes(x=date,y=right_ascension),
            data = data_individual %>% 
              filter(object=='Sun') %>% 
              filter(date>mdy('3/20/2021'),date<mdy('3/22/2022')),
            color='black') +
  ylab('Right ascension (hours)') +
  theme(legend.position = 'top')

data_individual %>% 
  filter(solar_system_object,
         object!='Moon',
         object!='Sun') %>%
  mutate(object=ordered(object,planets)) %>%
  mutate(right_ascension_diff=(right_ascension-sun_right_ascension)%%24) %>%
  mutate(right_ascension_diff=if_else(right_ascension_diff>12,
                                      right_ascension_diff-24,right_ascension_diff)) %>%
  ggplot(aes(x=date,y=right_ascension_diff,color=object,shape=object)) +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  geom_ref_line(h=0) +
  geom_ref_line(h=12) +
  geom_ref_line(h=-12) +
  geom_point() +
  ylab('Sun-relative right ascension (hours)') +
  theme(legend.position = 'top')

#### Determining orbital radii

# Load orbital elements and first order corrections
# Note: orbital elements from [https://stjarnhimlen.se/comp/ppcomp.html]
orbital_elements <- read_csv('~/astrolabe_analysis/orbital_elements.csv') %>% # Only planets orbit the sun
  mutate(object=ordered(object,planets)) 

### Planetary orbital radii all at once, by looking at elapsed time between 
### right ascension recurrences

data_planetary <- data_individual %>%
  filter(solar_system_object&object!='Sun'&object!='Moon') %>% # Only planets orbit the sun
  mutate(object=ordered(object,planets)) # Order planets by orbit position

periods_radii <- data_planetary %>% 
  mutate(diff=right_ascension-sun_right_ascension) %>%  # Referencing against the sun
  group_by(object,round(diff/2)*2) %>% # Does the recurrence grouping: 2 hour bins 
  nest(data=date) %>%
  filter(length(data)>2) %>% # Get rid of bins with too few observations
  unnest(cols=data) %>%
  summarize(synodic=combn(date,2,function(x){as.duration(x[[1]]%--%x[[2]])/ddays(1)})) %>%
  mutate(synodic=map(synodic,function(x){x[[1]]})%>%unlist) %>% # Compute synodic periods (days)
  mutate(synodic=if_else((object=='Jupiter'|object=='Saturn')&synodic>600, # I observed two full periods for Jupiter and Saturn only, so...
                         synodic/2,  # ...if the period is unusually large we're looking at doubled periods, so split the period in half
                         synodic)) %>%   # ...otherwise it's fine as is
  filter(synodic>200) %>% # Short synodic periods aren't actually synodic periods; they're the planet being observed sequentially
  mutate(radius=if_else(object=='Venus',  # Compute orbital radii (AU) using Kepler
                        (synodic/(synodic+365))^(2/3),  # Venus is an inner planet
                        (synodic/(synodic-365))^(2/3))) %>% # all others are outer
  ungroup() %>%
  select(object,synodic,radius) %>%
  group_by(object)

periods_radii %>%
  ggplot(aes(synodic,object)) + geom_boxplot() +
  xlab('Synodic period (days)')

periods_radii %>% 
  left_join(orbital_elements,by='object') %>%
  filter(radius < 20) %>%
  ggplot() + geom_boxplot(aes(radius,object)) +
  geom_point(aes(a0,object),color='red',shape='square') + # Note: a0 is mean orbital radius (AU)
  xlab('Orbital semi-major axis (AU)')

# Summary table
periods_radii %>%
  left_join(orbital_elements,by='object') %>%
  filter(radius < 20) %>%
  group_by(object) %>%
  summarize(n=n(),
            true=mean(a0),
            mean=mean(radius),
            median=median(radius),
            sd=sd(radius))

### Planet actual locations

data_planetary_extended <- data_planetary %>% 
  left_join(orbital_elements,by='object')  %>%
  mutate(n=as.duration(ymd('2000-01-01') %--% date)/ddays(1), # Compute corrections to orbital elements
         NN=(n0+n*n1)*pi/180, 
         II=(i0+n*i1)*pi/180,
         WW=(w0+n*w1)*pi/180,
         AA=a0+n*a1,
         EE=e0+n*e1,
         MM=(m0+n*m1)*pi/180,
         pos=heliocentric(NN,II,WW,AA,EE,MM)) %>% # Predict planet heliocentric positions
  unnest(pos) %>%
  mutate(sidereal_time=utc2s(Longitude+360,date),
         true_right_ascension=planet_ra(n,xh,yh,zh)%%24,
         true_declination=planet_dc(n,xh,yh,zh),
         true_elevation=rad2ce(Latitude,
                               true_right_ascension,
                               true_declination,
                               sidereal_time)*180/pi,
         true_azimuth=radc2a(Latitude,
                             true_right_ascension,
                             true_declination,
                             sidereal_time)*180/pi)

# Right ascension error
data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  summarize(n=n(),
            mean=mean(right_ascension_error),
            sd=sd(right_ascension_error))

data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  group_by(object) %>%
  summarize(n=n(),
            mean=mean(right_ascension_error),
            sd=sd(right_ascension_error))

data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  ggplot(aes(right_ascension_error,object)) + 
  geom_boxplot() +
  xlab('Right ascension error (hours)')

data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  ggplot(aes(right_ascension_error)) + 
  geom_histogram(bins=50)

data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  ggplot(aes(date,right_ascension_error,color=object,shape=object)) +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  geom_point()

# Checking if Chaucer is right: when a planet is near due south, 
# right ascension error increases
data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  ggplot(aes(true_azimuth%%360,right_ascension_error,color=object,shape=object)) +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  geom_ref_line(v=180) +
  geom_point() +
  xlab('Azimuth (degrees)') +
  ylab('Right ascension error (hours)') +
  theme(legend.position = 'top')

# Seeing if declination impacts right ascension error.  
# It seems that Jupiter is generally worse, but otherwise not much of an impact
data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6) %>%
  ggplot(aes(true_declination,right_ascension_error,color=object,shape=object)) +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  geom_point() +
  xlab('Declination (degrees)') +
  ylab('Right ascension error (hours)') +
  theme(legend.position = 'top')

# Testing East/West azimuth dependence of positive/negative right ascension errors
data_planetary_extended %>%
  mutate(right_ascension_error=(right_ascension-true_right_ascension+6)%%12-6,
         easterly=true_azimuth%%360<180,
         pos_ra_error=right_ascension_error>0) %>%
  group_by(easterly,pos_ra_error) %>%
  summarize(n=n(),
            mean=mean(right_ascension_error),
            sd=sd(right_ascension_error)) 

# True right ascensions
data_planetary_extended %>%
  ggplot(aes(x=date,y=true_right_ascension,color=object,shape=object)) +
  geom_point() +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  ylab('Right ascension (hours)') +
  theme(legend.position = 'top')

# True planet locations
data_planetary_extended %>% 
  ggplot(aes(x=xh,y=yh,color=object,shape=object)) + geom_point() +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)

# True planet locations
data_planetary_extended %>% 
  mutate(pos=planet_xy_from_data(n,right_ascension,AA)) %>%
  unnest(pos) %>% # view()
  ggplot(aes(x=xhp,y=yhp,color=object,shape=object)) + geom_point() +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)

# Orrery plot of true positions (circles) connected to observed positions (triangles)
data_planetary_extended %>% 
  left_join(periods_radii %>% group_by(object) %>% summarize(radius=median(radius,na.rm = TRUE)),
            by='object') %>%
  mutate(pos=planet_xy_from_data(n,right_ascension,radius)) %>%
  unnest(pos) %>%
  ggplot() + 
  geom_segment(aes(x=xh,y=yh,
                   xend=xhp,yend=yhp,color=object,shape=object)) +
  geom_point(aes(x=xh,y=yh),shape=1) +
  geom_point(aes(x=xhp,y=yhp),shape=2) +
  scale_color_manual(values=planet_colors)+
  scale_shape_manual(values=planet_shapes)+
  theme(legend.position = 'top') +
  coord_fixed()

### Determining Venus's orbital radius by parallax

data_planetary %>%
  filter(object=='Venus') %>%
  mutate(diff=(right_ascension-sun_right_ascension+6)%%12-6) %>%  # Referencing against the sun
  ggplot(aes(date,diff)) +
  geom_point() +
  geom_ref_line(h=0) +
  geom_ref_line(h=12) +
  geom_ref_line(h=-12) +
  ylab('Sun-relative right ascension of Venus (hours)') 
  
data_planetary %>%
  filter(object=='Venus') %>%
  mutate(diff=(right_ascension-sun_right_ascension+6)%%12-6) %>%  # Referencing against the sun
  ggplot(aes(date,diff,color=phase)) +
  geom_point() +
  geom_ref_line(h=0) +
  geom_ref_line(h=12) +
  geom_ref_line(h=-12) +
  ylab('Sun-relative right ascension of Venus (hours)') 

data_planetary %>%
  filter(object=='Venus') %>%
  mutate(diff=right_ascension-sun_right_ascension) %>%  # Referencing against the sun
  mutate(radial_estimate=abs(sin(diff*pi/12))) %>%
  ggplot(aes(radial_estimate)) + 
  geom_histogram(bins=20) +
  geom_ref_line(v=orbital_elements[orbital_elements$object=='Venus',]$a0) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10)) +
  xlab('Venus orbital semi-major axis (AU)')

data_planetary %>%
  filter(object=='Venus') %>%
  mutate(diff=right_ascension-sun_right_ascension) %>%  # Referencing against the sun
  mutate(radial_estimate=abs(sin(diff*pi/12))) %>%
  summarize(n=n(),
            semimajor_mean=mean(radial_estimate),
            semimajor_median=median(radial_estimate),
            stddev=sd(radial_estimate))

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
  ylab('Sun-relative right ascension of Venus (hours)') +
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
  ylab('Sun-relative right ascension of Mars (hours)') +
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
  ylab('Sun-relative right ascension of Jupiter (hours)') +
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
  ylab('Sun-relative right ascension of Saturn (hours)') +
  theme(legend.position = 'top')

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

