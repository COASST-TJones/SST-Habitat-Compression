################### LOAD LIBRARIES #####################################
# SSTa maps
library(maps)
library(mapdata)
library(maptools)
library(RColorBrewer)
library(lattice)
library(RNetCDF)
library(PBSmapping)
library(ncdf4)
library(PBSmapping)
library(dplyr) 
library(lubridate) 
library(tidync) 
library(rerddap) 
library(doParallel) 

mapfun <- function(xlim, ylim){
  
  data(worldLLhigh) # load the nepacLL data set
  plotMap(worldLLhigh, # plot the nepacLL data set
          xlim=xlim, # limit the region horizontally
          ylim=ylim, # limit the region vertically
          plt=c(0.05, 0.95, 0.1, 0.95), # specify the plot region size
          col=rgb(0.9,0.9,0.9,1),
          bg=rgb(1,1,1,1),
          tck=c(-0.03), # set the tick mark length
          cex = 1.2, # adjust the font size
          mgp=c(1.9, 0.7, 0)) 
}

########################################################################

################### ACCESS DATA ########################################
# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", 
              url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# Function for downloading data for the Channel Islands by time-period
OISST_sub_dl <- function(time_df){
  OISST_dat <- griddap(datasetx = "ncdcOisst21Agg_LonPM180", 
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                       time = c(time_df$start, time_df$end), 
                       zlev = c(0, 0),
                       latitude = c(31, 37),
                       longitude = c(-125, -116),
                       fields = "sst")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(t = time, temp = sst) %>% 
    select(longitude, latitude, t, temp) %>% 
    na.omit()
}

# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:5,
                       start = as.Date(c("1982-01-01", "1990-01-01", 
                                         "1998-01-01", "2006-01-01", "2014-01-01")),
                       end = as.Date(c("1989-12-31", "1997-12-31", 
                                       "2005-12-31", "2013-12-31", "2019-12-31")))

# Access data
OISST_data <- dl_years %>% 
    group_by(date_index) %>% 
    group_modify(~OISST_sub_dl(.x)) %>% 
    ungroup() %>% 
    select(longitude, latitude, t, temp)

# add doy and year labels
OISST_data$doy <- yday(OISST_data$t)
OISST_data$year <- year(OISST_data$t)

####################################################################

############## PLOTTING DAILY TEMPERATURE MAPS #####################
library(colorRamps)

# Function for a day of year and year
plotSST <- function(doy, year){
  # Now extract SST data
  sstselect <- OISST_data[OISST_data$doy == doy &
                            OISST_data$year == year,]
  
  # define color palette
  colpal <- matlab.like2(35)
  
  # Identify colors
  colmat <- floor(sstselect$temp*2)-16
  
  # Plot on map
  mapfun(xlim = c(360-126, 360-116), ylim = c(30,38))
  
  # For the given selection add polygons
  for(i in 1:nrow(sstselect)){
    latD <- sstselect$latitude[i]
    lonD <- sstselect$longitude[i]
    
    # Add polygon
    polygon(x=360+lonD+c(-0.125,0.125,0.125,-0.125),
              y=latD+c(-0.125,-0.125,0.125,0.125), 
              col=colpal[colmat[i]], 
              border=NA)      
  }
  
  l.vals <- seq(from=8.5, to=25, by=0.5)
  u.vals <- l.vals + 0.5
  legV <- paste(l.vals, u.vals, sep="-")
  legend(x=244, y=37.5, col=colpal, pch=15, xpd=NA,
         legend=legV, bty="n", pt.cex=2, ncol=2, title="SST (C)")
}

plotSST(doy=210, year=2013)

######################################################################

############# CALCULATE CLIMATOLOGY ##################################
# Create a climatology based on 1982 to 2011
OISST_data$latlon <- paste(OISST_data$latitude,
                           OISST_data$longitude,
                           sep="_")

sst.clim <- list()

for(i in 1:365){
  # Find days of year within +/- 5 days of doy
  daycent <- as.Date(paste(1982:2011, "-01-01", sep=""))+i-1
  
  # construct date lists
  datemin <- daycent-5
  datemax <- daycent+5
  
  datewant <- seq(from=datemin[1], to=datemax[1], by=1)
  for(j in 2:length(datemin)){
    datewant <- c(datewant,seq(from=datemin[j], to=datemax[j], by=1))
  }

  # extract all IDs corresponding to those dates
  Dwant <- OISST_data[OISST_data$t %in% datewant,]
  
  # Calculate mean by lat long
  sst.clim[[i]] <- Dwant %>% 
    group_by(latlon) %>% 
    summarise(mean=mean(temp), doy=i)
}

# Combine sst clim lists into a single data frame
sst_climatology <- bind_rows(sst.clim, .id = "doy")

#####################################################################

############ CALCULATE ANOMALIES #####################################
# Now merge climatology with the original dataset
OISST_data$mID <- paste(OISST_data$latlon, OISST_data$doy,sep="_")
sst_climatology$mID <- paste(sst_climatology$latlon, sst_climatology$doy,sep="_")

sst_df <- merge(x=OISST_data, 
                y=sst_climatology,
                by="mID", all.x=T)
sst_df <- sst_df[,c("t","latitude","longitude","doy.x","year",
                    "temp","mean")]
names(sst_df)[4] <- "doy"
sst_df <- sst_df[order(sst_df$t),]

# calculate anomaly
sst_df$anom <- sst_df$temp - sst_df$mean

# Create anomaly plot
colmatch <- function(vv){
  if(is.na(vv)){
    return(rgb(0,0,0,0))
  } else {
    bb <- ceiling(abs(vv)*2)
    if(bb >=9){
      bb <- 9 
    }
    if(vv < 0){
      return(brewer.pal(n=9, name="Blues")[bb])
    } else {
      return(brewer.pal(n=9, name="Reds")[bb])
    }
  }
}

# Function for a day of year and year
plotANOM <- function(doy, year){
  # Now extract SST data
  sstselect <- sst_df[sst_df$doy==doy & sst_df$year == year,]
  
  # Plot on map
  mapfun(xlim = c(360-126, 360-116), ylim = c(30,38))
  
  # For the given selection add polygons
  
  for(i in 1:nrow(sstselect)){
    latD <- sstselect$latitude[i]
    lonD <- sstselect$longitude[i]
    
    # Add polygon
    polygon(x=360+lonD+c(-0.125,0.125,0.125,-0.125),
            y=latD+c(-0.125,-0.125,0.125,0.125), 
            col=colmatch(sstselect$anom[i]), 
            border=NA)      
  }
  
  pal <- c(brewer.pal(n=9, name="Blues")[9:1],brewer.pal(n=9, name="Reds")[1:9])
  
  l.vals <- seq(from=-4.5, to=4, by=0.5)
  u.vals <- l.vals + 0.5
  legV <- paste(l.vals, u.vals, sep=" to ")
  legend(x=244, y=37.5, col=pal, pch=15, xpd=NA,
         legend=legV, bty="n", pt.cex=2, ncol=2, title="SSTa (C)")
  
}

plotANOM(doy=210, year=2019)

#######################################################################

######### CALCULATE TIME-SERIES OF COMPRESSION ########################
# Go through and calculate average SSTa and prop above 1C
sst_df$sst_gt1 <- 0
sst_df$sst_gt1[sst_df$anom > 1] <- 1

# Summarise
ssta.ts <- sst_df %>% 
  group_by(t) %>%
  summarise(ave.SSTa = mean(anom), prop.SSTa.gt1 = sum(sst_gt1)/595)

# subsample every 5th day 

ssta.ts.sub <- ssta.ts[seq(from=15, to=13879, by=15),]

pp1 <- brewer.pal(n=9, name="YlOrRd")[3:9]
plot(prop.SSTa.gt1 ~ t, data=ssta.ts.sub, type="h", ylab="Habitat compression (Area SSTa > 1C)", xlab="", 
     xaxt="n", lwd=2,
     col=pp1[floor(ssta.ts.sub$prop.SSTa.gt1*5)+1])

axis(side=1, at=as.Date(paste(1982:2020, "-01-01", sep="")), labels=1982:2020)

################################################################

