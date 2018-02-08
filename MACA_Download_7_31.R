#Download MACA Data through OPeNDAP
#Alec Nelson
# 7/31/16

############################
#setwd("~/Documents/Clemson_Files/Spring_2016/MastersThesis/DataAnalysis")
setwd("/Volumes/Logos/DataAnalysis")

list.of.packages <- c("raster","pracma", "reshape","car","compute.es","effects","rgdal","fields","chron", "ff","downloader",
                      "multcomp","pastecs","lme4","MuMIn", "ncdf4","sp","dismo","stringr", "data.table","RCurl","rio")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)} 

#Load all packages
lapply(list.of.packages, require, character.only = TRUE)
############################
FileList <- read.table("macav2livneh_urls-2.txt")
nrFiles<-nrow(FileList)
FailedList<-list()
for(i in 1:length(FileList[,1])){
  URL.i<-str_replace(as.character(FileList[i,]),"fileServer/NWCSC_INTEGRATED_SCENARIOS_ALL_CLIMATE","dodsC" )
  if(i==1){URL.list<-URL.i}else{URL.list<-c(URL.list,URL.i)}}
FileList<-cbind(FileList,URL.list)
str(FileList)

#Load NetCDF Information
for(j in 1:length(FileList[,1])){
  # Start the clock!
  ptm <- proc.time()
  #####################
  url.j<-as.character(FileList[j,1])
  filename.j<-substr(FileList[j,1],(gregexpr('macav2livneh', FileList[j,1])[[1]][2])+13,nchar(as.character(FileList[j,1])))
  destfile.j<-paste0(getwd(),"/",filename.j)
  maxretryattempts <- 5 #If there is an error downloading data how many times to retry
  ######## Try the download:
  for(t in 1:maxretryattempts){
    error.t<-NULL
    tryCatch(nc_open(filename.j),error = function(e) error.t<<-(e))
      if(is.null(error.t)){
        nc<-tryCatch(nc_open(filename.j, write=FALSE, readunlim=FALSE),error = function(e) print(e))
        try(print(j))
        try(print(nc[[1]]))
        cat("File successfully downloaded")
          break}
        #The data wasn't previously downloaded so let's attempt to download it
        cat("(",j,"/",nrFiles,") ","Downloading ", filename.j , "\t\t Attempt: ", t , "/", maxretryattempts,"\n")
        tryCatch(download.file(url.j,destfile=destfile.j,mode="wb"),error = function(e) print(e))
    ######## Load and Check on Max Retry attempt:
      if(t==maxretryattempts){
        error.t<-NULL
        tryCatch(nc_open(filename.j),error = function(e) error.t<<-(e))
        if(is.null(error.t)){
          #The variable exists so confirm download
          nc<-tryCatch(nc_open(filename.j, write=FALSE, readunlim=FALSE),error = function(e) print(e))
          try(print(j))
          try(print(nc[[1]]))
          cat("File successfully downloaded on last attempt")}
        else{
          FailedList<-c(FailedList,filename.j)
          cat("Error in downloading",filename.j)}
      }  
  }
  # Stop the clock
  proc.time() - ptm
  Sys.sleep(2)
}

#CHECK FAILED LIST!!:
print(FailedList)

## SHOW SOME METADATA 
nc[[1]]
nc[[11]]
class(nc)
names(nc)
(nc$var)$vals
print(nc)

## DISPLAY INFORMATION ABOUT AN ATTRIBUTE
ncatt_get(nc,names(nc$var))

######################################
# Testing methods of downloading data from WebURL
######################################
as.character(FileList[1,1])
# f.url <- fread(as.character(FileList[1,1]))
# #fname <- file.choose()
# ff.data<-read.delim.ffdf(as.character(FileList[1,1]))
# myfile <- getURL(as.character(FileList[1,1]), ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
# data <- import(as.character(FileList[1,1]))
# ncin <- open.connection(as.character(FileList[1,1]))

error.t<-NULL
tryCatch(nc_open(filename.j),error = function(e) error.t<<-(e))
tryCatch(nc_open(as.character(FileList[1,2])),error = function(e) error.t<<-(e))
nc <- nc_open(as.character(FileList[1,2]))

!is.null(error.t)

if(nc){
  t<-TRUE
}else{
  t<-FALSE
}
t


# This illustrates how to read data one timestep at a time.  In this
# example we will elaborately show how to deal with a variable whose
# shape is completely unknown (i.e., how many dimensions, and what their
# sizes are).  We will also, for illustration of a common case, show how
# to read in the values of the time dimension at each timestep.
v3      <- nc$var[[1]]
varsize <- v3$varsize
ndims   <- v3$ndims
nt      <- varsize[ndims]  # Remember timelike dim is always the LAST dimension!
## DEFINE OUR VARIABLE NAME 
var=names(nc$var)

for( i in 1:nt ) {
  # Initialize start and count to read one timestep of the variable.
  start <- rep(1,ndims)    # begin with start=(1,1,1,...,1)
  start[ndims] <- i    # change to start=(1,1,1,...,i) to read timestep i
  count <- varsize    # begin w/count=(nx,ny,nz,...,nt), reads entire var
  count[ndims] <- 1    # change to count=(nx,ny,nz,...,1) to read 1 tstep
  
  data3 <- ncvar_get( nc, var, start=start, count=count )
  
  # Now read in the value of the timelike dimension
  timeval <- ncvar_get( nc, v3$dim[[ndims]]$name, start=i, count=1 )
  
  print(paste("Data for variable",v3$name,"at timestep",i,
              " (time value=",timeval,v3$dim[[ndims]]$units,"):"))
  print(data3)
}

nc_close(nc)


## GET DATA SIZES: http://www.inside-r.org/packages/cran/ncdf4/docs/ncvar_get
## NOTE: FILE DIMENSIONS ARE lon,lat,time
v3 <- nc$var[[1]]
lonsize <- v3$varsize[1]
latsize <- v3$varsize[2] 
endcount <- v3$varsize[3] 


### DEFINE OUR POINT OF INTEREST 
## NOTE: MAKE SURE TO CHECK WHETHER YOUR SOURCE STARTS COUNTING AT 0 OR 1
## e.g. ncdf4 PACKAGE STARTS COUNTING AT 1 BUT OPeNDAP DATASET ACCESS FORM STARTS AT 0:
## http://convection.meas.ncsu.edu:8080/thredds/dodsC/pinemap/maca/past/macav2livneh_pr_bcc-csm1-1-m_historical_1970_1989_CONUS.nc.html
lons = ncvar_get(nc,varid="lon")
lats = ncvar_get(nc,varid="lat")

lon=478
lat=176

## DEFINE OUR VARIABLE NAME 
var=names(nc$var)

## READ THE DATA VARIABLE (e.g. precipitation IN THIS CASE): http://www.inside-r.org/packages/cran/ncdf4/docs/ncvar_get 
## AND http://stackoverflow.com/questions/19936432/faster-reading-of-time-series-from-netcdf
## ORDER OF start= AND count= IS BASED ON ORDER IN BRACKETS AFTER VARIABLE NAME (SHOWN WHEN DISPLAYING YOUR METADATA)
## FROM THE DOCUMENTATION... "If [start] not specified, reading starts at the beginning of the file (1,1,1,...)."
## AND "If [count] not specified and the variable does NOT have an unlimited dimension, the entire variable is read. 
## As a special case, the value "-1" indicates that all entries along that dimension should be read."
data <- ncvar_get(nc, var, start=c(lon,lat,1),count=c(1,1,endcount))
class(data)
data <- ncvar_get(nc, var,count=c(1,1,endcount))
## READ THE TIME VARIABLE
time <- ncvar_get(nc, "time", start=c(1),count=c(endcount))
## CONVERT TIME FROM "days since 1900-01-01" TO YYYY-MM-DD
time=as.Date(time, origin="1900-01-01") ##note: assumes leap years! http://stat.ethz.ch/R-manual/R-patched/library/base/html/as.Date.html
# PUT EVERYTHING INTO A DATA FRAME
c <- data.frame(time,data)

## CLOSE THE FILE
nc_close(nc)

## PLOT THE DATA 
plot(c$time,c$data,main=paste("Daily ",var," for ",c$time[1]," through ",c$time[nrow(c)], " at ",lat,",",lon,sep=""),xlab="Date",ylab="Precipitation (mm)")

## EXTRCT RASTER AND ATTEMPT TO PLOT IT
data.i <- ncvar_put( nc, var, data)

variable <- raster(data.i, varname=names(nc$var))
plot(variable)


# Create dimensions lon, lat, level and time
dim_lon  <- ncdim_def('longitude', 'degrees_east', seq(-179.75,179.75,by=0.5))
dim_lat  <- ncdim_def('latitude', 'degrees_north', seq(89.75,-89.75,by=-0.5))
dim_lev  <- ncdim_def('level', 'level/index', 1)
dim_time <- ncdim_def('time', "years since 1990-01-01", 1:12, unlim=T)

# Create a new variable "precipitation", create netcdf file, put updated contents on it and close file
# Note that variable "data" is the actual contents of the original netcdf file
var_out <- ncvar_def('precipitation', 'mm/day', list(dim_lon,dim_lat,dim_lev,dim_time), 9.e20)
ncid_out <- nc_create('updated_netcdf.nc', var_out)
ncvar_put(ncid_out, var_out, data, start=c(1, 1, 1, 1), count=c(720, 360, 1, 12))
nc_close(ncid_out)





## TEST RASTER AND PLOTTING MECHANICS
writeRaster(control, "Data2.nc", "CDF", overwrite=TRUE,
            varname="specific_humidity",varunit="kg kg-1",longname="Daily Mean Near-Surface Specific Humidity",
            xname="lon",yname="lat",bylayer=FALSE,NAflag=-9999)

control = nc_open("Test_File1.nc", write=FALSE, readunlim=FALSE)
cat(paste(control$filename,"has",control$nvars,"variables"), fill=TRUE)

lonmat  = ncvar_get(nc=control,varid="lon")   # reads entire matrix
latmat  = ncvar_get(nc=control,varid="lat")    # ditto
timearr = ncvar_get(nc=control,varid="time")        # reads entire time array


# Get date index from the file
idx <- getZ(control)

# Put coordinates and extract values for all time steps
coords <- matrix(c(-55, -31), ncol = 2) # here -55 is longitude and -31 is latitude
vals <- extract(b, coords, df=T)

# Merge dates and values and fix data frame names
df <- data.frame(idx, t(vals)[-1,])
rownames(df) <- NULL
names(df) <- c('date','value')

# Take a peek on extracted values:
head(df)




plot( lonmat, latmat, main='The (MM5) grid locations')
US( add=T )

summary(control)

