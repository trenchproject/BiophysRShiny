#source biophysical model
setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\CAREER\\Rcode\\MESA\\biophys\\")
source("BiophysModel_mat.R")
#source zenith angle calculating function
source("ZenithAngleFunction.R")
#source dtr function
source("DTRfunction_Feb2013.R")

library("RAtmosphere")
library(sp)
library(raster)
#----------------------------
#PARAM

#svl= 60 #mm
#mass= 3.55*10^-5*(svl)^3.00 #Tinkle and Ballinger 1972 (g)


#BIOPHYSICAL MODELING ***************************************

# read climate data
#10' DATA
setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\CAREER\\Rcode\\MESA\\climate\\")
climdata=read.csv("ClimateData_2160x684latlon.csv")

#read surface temperature data
tsurfdata=read.csv("Tsurface_NARRall.csv");
#match to climate data
match1= match(climdata$GridID, tsurfdata$GridID)
tsurfdata= tsurfdata[match1,]

lat= climdata[,"lat"] ; #latitude in degrees
lon= climdata[,"lon"] ; #longitude in degrees

#use surface T data 
Tm= tsurfdata[,2:13]; #mean of an average day of each month
Tr= tsurfdata[,14:25]; #diurnal temperature range of an average day of each month
Wind= climdata[,4]; #mean wind speed m/s
Wind[]=0.1 #assume 0.1 m/s wind
Albedo= climdata[,5:8]; #Albedo percentage for Jan / Apr / Jul / Oct
Elev= climdata[,33]; #mean elevation (m)
TSmax=climdata[,34:45];  #max soil T, mean 14:00hr temperature for five days in the middle of each month (K), from LDAS
#TSr=climdata[,46:57];  #range between 02:00 and 14:00hr temperature for five days in the middle of each month (K), from LDAS
Tairm= climdata[,9:20]; #mean of an average day of each month
Tairr= climdata[,21:32]; #diurnal temperature range of an average day of each month

#---------------------------------------
#Calculate hourly temperature trend

Thourly= function(monthk, hourk){

J.all=c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
J= J.all[monthk]  # Julian calendar day

Trise.set= suncalc(J, Lat = lat, Long = lon, UTC = FALSE)
set= Trise.set$sunset
rise= Trise.set$sunrise

psi_rad=zenith.rad(J, lat, lon, hourk) #Calculate zenith in radians
psi_rad[psi_rad>=pi/2]=pi/2 #set zenith position below horizon

#dayhrs= hourk[which(hourk<sunset & hourk>sunrise)]
#nighthrs= hourk[-dayhrs]	

# Calculate the diurnal temperature trend, Parton and Logan 1981
#surface temp
Tx= Tm[, monthk] + Tr[, monthk] / 2  # maximum daily temperature
Tn= Tm[, monthk] - Tr[, monthk] / 2 # minimum daily temperature

Ta=apply(cbind(Tx,Tn,rise, set),FUN=Thour.mat, MARGIN=1, Hr=hourk,  alpha=1.86, beta= -0.17, gamma=2.20)

#air temp
Tairx= Tairm[, monthk] + Tairr[, monthk] / 2  # maximum daily temperature
Tairn= Tairm[, monthk] - Tairr[, monthk] / 2 # minimum daily temperature

Taira  = apply(cbind(Tairx,Tairn,rise, set),FUN=Thour.mat, MARGIN=1, Hr=hourk,  alpha=1.86, beta= -0.17, 
               gamma=2.20)

return(cbind(Taira,Ta))
}

#---------------------------------------
#Write Te raster

TeRaster= function(monthk, hourk, svl, mass1, absorb){

  #Calculate Te 
  albk= floor (monthk/3.01)+1 # compute season for albedo measurement
  rho_S= Albedo[, albk]/100 # (Table 12.2)
  
  J.all=c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  J= J.all[monthk]  # Julian calendar day
  
  psi_rad=zenith.rad(J, lat, lon, hourk) #Calculate zenith in radians
  psi_rad[psi_rad>=pi/2]=pi/2 #set zenith position below horizon
  
  Tas= Thourly(monthk, hourk)
  
  Tmat= cbind(Tas,Wind,psi_rad,rho_S,Elev)
  
  Te  = apply(Tmat,FUN=biophys.mat, MARGIN=1, SVL=svl, MASS=mass1, J=J, absorb=absorb)
  
  Te.dat= as.data.frame(cbind(lat,lon,Te, Tas[,1]))
  
  #MAP
  ys= c(-40,50)
  xs= c(-180,100) 
  Te.dat.sub= subset(Te.dat, lat>ys[1] & lat<ys[2] & lon>xs[1] & lon<xs[2])
  
  xy.sp= cbind( Te.dat.sub$lon, Te.dat.sub$lat)
  xy.cc= coordinates(xy.sp)
  bbox(xy.sp)
  
  grd1 <- SpatialPixelsDataFrame(points=xy.sp, data = as.data.frame(as.numeric(Te.dat.sub[,c(3)])), tolerance=0.2)
  Te.raster= raster(grd1)
  # replacing NA's by zero
#  Te.raster[is.na(Te.raster[])] <- 0 
  
  projection(Te.raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  writeRaster(Te.raster, "Tegrid.grd", format='raster', overwrite=TRUE)
  
  return(Te.dat.sub[,3])
}  

