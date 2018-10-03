#Data Processing for MACA Data
#Alec Nelson
# 10/17/16

############################
basedirectory <- "C:\\Users\\ahn11803\\Documents\\GitHub\\Spatial-R_Code"
inputdata_path <- "D:/DataAnalysis/TestData"
countryshape_path <- "D:/DataAnalysis/App_Boundary_SHP"
#outputdata_path <- "D:/DataAnalysis/TestResults5"
#monthlydata_path<-"D:/DataAnalysis/MonthlyData3"
sampledata_path<-"D:/DataAnalysis/TestSample"
SampledCSV_path<-"D:/DataAnalysis/SampledClimateData"
EnvirData_path<-"D:/DataAnalysis/EnvironmentData"
CombinIter_NoSeas_path<-"D:\\DataAnalysis\\MultiVarLimit_Iter_NO_SAdj"
CombinIter_SeasAdj_path<-"D:\\DataAnalysis\\MultiVarLimit_Iter_SeasonAdj"
DataAnalysis_path<- "D:/DataAnalysis"

setwd(inputdata_path)

list.of.packages <- c("raster","pracma", "reshape","car","compute.es","effects","rgdal","fields",
                      "chron", "ff","downloader","magrittr","maptools","GSIF","rgeos","ggplot2",
                      "multcomp","pastecs","data.table","MuMIn", "ncdf4","sp","dismo","stringr", 
                      "data.table","RCurl","rio","RNetCDF","parallel","qicharts","qcc","zoo","dplyr","purrr","plyr",
                      "geoR","geoRglm","MSQC","coda","MASS","relaimpo","arcgisbinding","lme4","glmm","nlme","arm","rms")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)} 

#Load all packages
lapply(list.of.packages, require, character.only = TRUE)
########################################################
########################################################
# Process data into a Raster format and store for a smaller (study area) extent
setwd(inputdata_path)
raster_file_list <- list.files(pattern = ".nc" , all.files = FALSE , full.names = FALSE )
countryshape<-rgdal::readOGR( countryshape_path, "AppalachianLvl2_GCS84" )
sample_points<-rgdal::readOGR("D:/DataAnalysis", "SamplePoints1")

# print(nc.i)
# summary(nc.i)
# str(nc.i)
# nc.i[14]
# lon <- ncvar_get(nc.i,"lon")
# nlon <- dim(lon)
# head(lon)
# lat <- ncvar_get(nc.i,"lat",verbose=F)
# nlat <- dim(lat)
# head(lat)
# print(c(nlon,nlat))
# t <- ncvar_get(nc.i,"time")
# tunits <- ncatt_get(nc.i,"time","units")
# nt <- dim(t)
#tmp_array <- ncvar_get(nc.i,var.i)
# dlname <- ncatt_get(nc.i,var.i,"long_name")
# dunits <- ncatt_get(nc.i,var.i,"units")
# fillvalue <- ncatt_get(nc.i,var.i,"_FillValue")
# dim(tmp_array)

##################
# Start the clock!
ptm <- proc.time()
##################

#Iterate process through raster file list

for(f in 1:length(raster_file_list)){
  
  raster_file_name<-raster_file_list[f]
  pos = regexpr('_', raster_file_name)[1] 
  var_name.f<-substr(raster_file_name, 1, (pos-1))
  raster_years<-substr(raster_file_name, nchar(raster_file_name)-23, nchar(raster_file_name)-15)
  outputpath_name.f<-paste0("D:/DataAnalysis/DailyData_MACA/Daily_",var_name.f,"_",raster_years)
  
  dir.create(outputpath_name.f)
  outputdata_path <- outputpath_name.f
  
  ########################################################
  #for(j in 1:length(raster_file_list)){
    setwd(inputdata_path)
    nc.i<-nc_open(as.character(raster_file_list[f]))
    var.i<-ncatt_get(nc.i,names(nc.i$var))$standard_name
    v3  <- ncatt_get(nc.i,names(nc.i$var))$standard_name
    v2 <-nc.i$var[[1]]
    varsize <- v2$varsize
    ndims   <- v2$ndims
    nt      <- varsize[ndims]
    nc_close(nc.i)
    
    print(paste0("Loaded parameters for: ",as.character(raster_file_list[f])))
    
    for( i in 1:nt ) {
      setwd(inputdata_path)
      data.r<-raster( as.character(raster_file_list[f]) ,ymn = 25.125, ymx = 52.875, xmn = 235.375, xmx = 293,
                      ncdf=TRUE, varname=v3,lvar=3,level=1,band=i)
      timeval<-getZ(data.r)
      #proj4string( data.r )
      shift.raster <- raster::shift(data.r,-360)
      data.projected<-projectRaster( shift.raster , crs = proj4string( countryshape ),method = "bilinear" )
      data.proj.r<-raster::crop( data.projected , extent( countryshape ) )
      #   plot(data.proj.r)
      #   plot(countryshape,add=T)    
      filename.i<-paste0(substr(raster_file_list[f],1,(nchar(raster_file_list[f])-24)),timeval,".tif")
      print(paste0("Saving file: ",filename.i))
      setwd(outputdata_path)
      rf <- writeRaster(data.proj.r, filename=filename.i, format="GTiff", overwrite=TRUE)
      setwd(inputdata_path)
      #data.mean.r<- raster::extract( data.proj.r , countryshape ,f=TRUE, fun = mean, na.rm = TRUE  )
    }
  #}
  
  ####################################################################################
  ####################################################################################
  
  setwd(outputdata_path)
  #Z-score calculation in R
  
  #Stack all of the raster layers
  clipped_file_list <- list.files(pattern = ".tif" , all.files = FALSE , full.names = FALSE )
  #clipped_file_list <- list.files(pattern = "1956" , all.files = FALSE , full.names = FALSE )
  stack.g <- stack() 
  for(g in 1:length(clipped_file_list)){
    stack.g<-stack(stack.g,raster(clipped_file_list[g]))
    print(paste("Added file",clipped_file_list[g],"to stack"))
  }
  
  
  ########################################################
  ########################################################
  
  # Summarize data by monthly raster layers
  
  stats.m<-c("Mean","Max","Min")
  output.stat.months<-c("D:/DataAnalysis/MonthlyData_MACA/MonthlyData_Mean","D:/DataAnalysis/MonthlyData_MACA/MonthlyData_Max","D:/DataAnalysis/MonthlyData_MACA/MonthlyData_Min")
  
  output.stat.months.f1<-paste0(output.stat.months[1],"_",var_name.f,"_",raster_years)
  dir.create(output.stat.months.f1, showWarnings = FALSE)
  output.stat.months.f2<-paste0(output.stat.months[2],"_",var_name.f,"_",raster_years)
  dir.create(output.stat.months.f2, showWarnings = FALSE)
  output.stat.months.f3<-paste0(output.stat.months[3],"_",var_name.f,"_",raster_years)
  dir.create(output.stat.months.f3, showWarnings = FALSE)
  
  output.stat.months<-c(output.stat.months.f1,output.stat.months.f2,output.stat.months.f3)
  
  stack.j<-stack.g
  
  for(m in 1:length(stats.m)){
    month.prev<-0
    year.prev<-1950
    stack.month.i<-stack()
    for(k in 1:nlayers(stack.j)){
      
      Date.k<-(substr(as.character(names(stack.j)[k]),nchar(names(stack.j)[k])-9,nchar(names(stack.j)[k])))
      Date.k<-gsub("[.]", "-", Date.k)
      #Dates.all<-c(Dates.all,Date.k)
      Month.k<-as.numeric(format(as.Date(Date.k, origin = "1900-01-01"), "%m"))
      Month.char.k<-format(as.Date(Date.k, origin = "1900-01-01"), "%m")
      Year.k<-as.numeric(format(as.Date(Date.k, origin = "1900-01-01"), "%Y"))
      Year.char.k<-format(as.Date(Date.k, origin = "1900-01-01"), "%Y")
      
      ## Change this to save layer to a particular month stack, 
      # then perform summary operation, then save to monthly stack
      
      if(Month.k > month.prev | Year.k > year.prev | k == (nlayers(stack.j))){
        if(k>1){
          #stack.month.stat<- cellStats(stack.month.i,stats.m[m],na.rm=TRUE)
          if(stats.m[m]=="Mean"){
            stack.month.stat<- mean(stack.month.i,na.rm=TRUE)
            print(paste0("Calculated ", stats.m[m], " of Previous Month"))
          }else if(stats.m[m]=="Max"){
            stack.month.stat<- max(stack.month.i,na.rm=TRUE)
            print(paste0("Calculated ", stats.m[m], " of Previous Month"))
          }else if(stats.m[m]=="Min"){
            stack.month.stat<- min(stack.month.i,na.rm=TRUE)
            print(paste0("Calculated ", stats.m[m], " of Previous Month"))
          }else{
            print("ERROR: NO STATISTICAL FUNCTION PERFORMED!")
          }
          
          filename.k<-paste0(stats.m[m],"_",(substr(as.character(names(stack.j)[k]),1,nchar(names(stack.j)[k])-10)),year.prev.char.k,"_",month.prev.char.k)
          setwd(output.stat.months[m])
          r.output <- writeRaster(stack.month.stat, filename=filename.k, format="GTiff", overwrite=TRUE)
          setwd(outputdata_path)
        }
        stack.month.i<-stack()
        stack.month.i<-stack(stack.month.i,stack.j[[k]])
        print(paste("Began new stack w/",names(stack.j[[k]]),"as first layer"))
      }else{
        stack.month.i<-stack(stack.month.i,stack.j[[k]])
        print(paste("Added layer",names(stack.j[[k]]),"to Month Stack"))
      }
      month.prev<-Month.k
      month.prev.char.k<-Month.char.k
      year.prev<-Year.k
      year.prev.char.k<-Year.char.k
    }
  }
  
  ################################################################
  
  stack.list<-list()
  
  for(m in 1:length(stats.m)){
    
    setwd(output.stat.months[m])
    
    #Stack all of the monthly raster layers
    monthly_file_list <- list.files(pattern = ".tif" , all.files = FALSE , full.names = FALSE )
    #clipped_file_list <- list.files(pattern = "1956" , all.files = FALSE , full.names = FALSE )
    
    stack.m <- stack() 
    for(g in 1:length(monthly_file_list)){
      stack.m<-stack(stack.m,raster(monthly_file_list[g]))
      print(paste("Added file",monthly_file_list[g],"to stack"))
    }
    
    assign(paste("stack.m.", stats.m[m], sep=""),stack.m)
    
    stack.add<-paste("stack.m.", stats.m[m], sep="")
    stack.list<-c(stack.list,eval(parse(text = stack.add)))
    
    setwd(basedirectory)
  }
  
  # stack.list[[1]]
  # stack.list[[2]]
  
  ################################################################
  # Test iterations of sampling over multiple layers
  # setwd(sampledata_path)
  # 
  # layer.t<-stack.m[[1]]
  # xres.t<-xres(layer.t)
  # yres.t<-yres(layer.t)
  # 
  # layer.SPdf<-as(layer.t,"SpatialPointsDataFrame")
  # layer.SPdf.sample<-GSIF::sample.grid(layer.SPdf, cell.size = c((xres.t*10), (yres.t*10)), n = 1)
  # layer.grid<-layer.SPdf.sample$grid
  # layer.grid.raster<-raster(layer.grid)
  # layer.resamp<-resample(layer.grid.raster,layer.t,method='ngb')
  # samp_strat.rand<- sampleStratified(layer.resamp, 10, xy = TRUE, sp=TRUE, na.rm = TRUE)
  # samp_strat.points<-as(samp_strat.rand,"SpatialPoints")
  #writeSpatialShape(samp_strat.rand,"D:/DataAnalysis/SamplePoints1.shp")
  # samplepoints.df<-as.data.frame(coordinates(samp_strat.points))
  # samplepointnum<-c(1:nrow(samplepoints.df))
  # samplepoints.df<-cbind(samplepointnum,samplepoints.df)
  ################################################################
  
  samp_strat.points<-sample_points
  samplepoints.df<-as.data.frame(coordinates(samp_strat.points))
  samplepointnum<-c(1:nrow(samplepoints.df))
  samplepoints.df<-cbind(samplepointnum,samplepoints.df)
  
  
  # samples.list<-list()
  
  for(x in 1:length(stack.list)){
    
    Samples.stack.m<-data.frame()
    stack.m<-stack.list[[x]]
    
    for(t in 1:nlayers(stack.m)){
      layer.t<-stack.m[[t]]
      layer.t.name<-names(layer.t)
      month.t<-substr(layer.t.name,nchar(layer.t.name)-1,nchar(layer.t.name))
      year.t<-substr(layer.t.name,nchar(layer.t.name)-6,nchar(layer.t.name)-3)
      var.t<-substr(layer.t.name,1,4)
      
      extract.sample<-raster::extract(layer.t,samp_strat.points)
      
      extract.df<-as.data.frame(extract.sample)
      sample.names<-rep(layer.t.name,nrow(samplepoints.df))
      sample.month<-rep(month.t,nrow(samplepoints.df))
      sample.year<-rep(year.t,nrow(samplepoints.df))
      sample.coordinates.t<-cbind(sample.names,sample.year,sample.month,samplepoints.df,extract.df)
      
      Samples.stack.m<-rbind(Samples.stack.m,sample.coordinates.t)
      print(paste("Sampled layer",layer.t.name,", layer",as.character(t),"out of", as.character(nlayers(stack.m))))
    }
    
    assign(paste("Samples.stack.m.", stats.m[x], sep=""),Samples.stack.m)
    
    # sample.add<-paste("Samples.stack.m.", stats.m[x], sep="")
    # samples.list<-list(samples.list,eval(parse(text = sample.add)))
    
  }

  #str(samples.list)
  
  str(Samples.stack.m.Mean)
  summary(Samples.stack.m.Mean)
  str(Samples.stack.m.Max)
  summary(Samples.stack.m.Max)
  str(Samples.stack.m.Min)
  summary(Samples.stack.m.Min)
  
  Samples.stack.m.Mean1<-Samples.stack.m.Mean[,-1]
  Samples.stack.m.Max1<-Samples.stack.m.Max[,-1]
  Samples.stack.m.Min1<-Samples.stack.m.Min[,-1]
  
  colnames(Samples.stack.m.Mean1)[colnames(Samples.stack.m.Mean1)=="extract.sample"] <- paste0(var_name.f,"_Mean")
  colnames(Samples.stack.m.Max1)[colnames(Samples.stack.m.Max1)=="extract.sample"] <- paste0(var_name.f,"_Max")
  colnames(Samples.stack.m.Min1)[colnames(Samples.stack.m.Min1)=="extract.sample"] <- paste0(var_name.f,"_Min")
  
  Samples.stack.total <- merge(Samples.stack.m.Mean1,Samples.stack.m.Max1,by=c("sample.year","sample.month","samplepointnum","coords.x1","coords.x2"),all=TRUE)
  Samples.stack.total <- merge(Samples.stack.total,Samples.stack.m.Min1,by=c("sample.year","sample.month","samplepointnum","coords.x1","coords.x2"),all=TRUE)
  
  filename.f<-paste0("D:/DataAnalysis/SampleStack_",var_name.f,"_",raster_years,".csv")
  
  write.csv(Samples.stack.total,file=filename.f,row.names=FALSE)

}

##################
#Check the time elapsed
proc.time() - ptm
##################

########################################################


str(Samples.stack.total)
summary(Samples.stack.total)
#unlink(output.stat.months[1])

setwd(SampledCSV_path)
csv_file_list <- list.files(pattern = ".csv" , all.files = FALSE , full.names = FALSE )

lastyear<-as.numeric(substr(csv_file_list[1],nchar(csv_file_list[1])-7,nchar(csv_file_list[1])-4))

filename.c<-paste0()

Samples.1<-read.csv(csv_file_list[12])
Samples.2<-read.csv(csv_file_list[15])
Samples.3<-read.csv(csv_file_list[9])

Samples.TasTotal <- merge(Samples.1,Samples.2,by=c("sample.year","sample.month","samplepointnum","coords.x1","coords.x2"),all=TRUE)
Samples.TasTotal <- merge(Samples.TasTotal,Samples.3,by=c("sample.year","sample.month","samplepointnum","coords.x1","coords.x2"),all=TRUE)

Samples.1950_1969<-Samples.TasTotal
Samples.1970_1989<-Samples.TasTotal
Samples.1990_2005<-Samples.TasTotal

Samples.TasTotal<-rbind(Samples.1950_1969,Samples.1970_1989,Samples.1990_2005)

newdata <- Samples.TasTotal[order(Samples.TasTotal$sample.year,Samples.TasTotal$sample.month,Samples.TasTotal$samplepointnum),]
write.csv(newdata,file="D:/DataAnalysis/SampleStack_ClimateTotal.csv",row.names=FALSE)


################################################################################################################
################################################################################################################
setwd(DataAnalysis_path)
SampleTotal<-read.csv("SampleStack_ClimateTotal.csv")
str(SampleTotal)
summary(SampleTotal)

YearMonth<-as.yearmon(paste(SampleTotal$sample.year[1],SampleTotal$sample.month[1],sep="-"))

for(y in 2:nrow(SampleTotal)){
  YearMonth.y<-as.yearmon(paste(SampleTotal$sample.year[y],SampleTotal$sample.month[y],sep="-"))
  YearMonth<-c(YearMonth,YearMonth.y)
  percent<-(y/nrow(SampleTotal))*100
  print(paste("Added YearMonth for",YearMonth.y,", row",as.character(y),"out of", as.character(nrow(SampleTotal)), "(",percent,"%)" ))
  }

SampleTotal.y<-cbind(SampleTotal,YearMonth)

row.ha.na<-apply(SampleTotal.y,1,function(x){any(is.na(x))})
sum(row.ha.na)
SampleTotal.final<-SampleTotal.y[!row.ha.na,]

str(SampleTotal.final)

write.csv(SampleTotal.final,file="D:/DataAnalysis/SampleStack_ClimateTotal_Final.csv",row.names=FALSE)

##############################################
##############################################
# Subtract Mean by Month for each point
setwd(basedirectory)
SampleTotal.final<-read.csv("SampleStack_ClimateTotal_Final.csv")
# 
# y<-rpois(24,16)
# qic(y)
# test.sample<-SampleTotal.final[SampleTotal.final$samplepointnum == 1,]
# str(test.sample)
# summary(test.sample)

SampleTotal.final.pqrs<-SampleTotal.final


samplepointnum.unique<-unique(SampleTotal.final$samplepointnum)
samplepoint.month.unique<-unique(SampleTotal.final$sample.month)
vars.final<-names(SampleTotal.final)[6:20]
col.nums<-c(6:20)

SampleTotal.final.pqrs %>%
  filter(sample.month==1 & samplepointnum==1) %>%
  select(tasmax_Mean) %>%
  
summarise()

purrr::map(col.nums,function(i){
  SampleTotal.final.pqrs %>%
    filter(sample.month==1 & samplepointnum==1) %>%
    select(i) %>%
    head
    #mean(.,na.rm=TRUE)
})
##############################

Season.Adj<-purrr::map(samplepointnum.unique,function(j){
    purrr::map(samplepoint.month.unique,function(i){
      SampleTotal.final.pqrs %>%
        filter(sample.month==i & samplepointnum==j) %>%
        select(col.nums) %>%
        mutate_each(funs(.-mean(.)))
    })
    })

class(Season.Adj)

Season.Adj[1][[1]][[1]]
Season.Adj[1][[1]][1]

unique(SampleTotal.final.pqrs$YearMonth)


# Unlist.test1<- do.call("rbind",lapply(Season.Adj,as.data.frame))
# Unlist.test2<- do.call("rbind",lapply(Season.Adj,as.data.frame,stringsAsFactors=FALSE))
# Unlist.test3<- data.frame(Reduce(rbind,Season.Adj))
# Unlist.test4<- data.frame(t(sapply(Season.Adj,c)))
# Unlist.test5<- do.call(rbind.data.frame,Season.Adj)
# 
# Unlist.test6<- data.frame(matrix(Season.Adj[1][[1]][[1]],nrow=nrow(SampleTotal.final.pqrs %>% 
#                           filter(sample.month== 1 & samplepointnum==1)),byrow=F),stringsAsFactors = FALSE)
# 
# Unlist.test6<- data.frame(matrix(Season.Adj[1][[1]][[1]],nrow=nrow(Season.Adj[1][[1]][[1]]),byrow=F),stringsAsFactors = FALSE)

Unlist.test6<- data.frame((Season.Adj[1][[1]][[1]]),stringsAsFactors = FALSE)


# Unlist.test7<- data.frame(lNames=rep(names(Season.Adj),lapply(Season.Adj,length)),lVal=unlist(Season.Adj))
# Unlist.test8<- ldply(Season.Adj,data.frame)

str(Unlist.test6)
#names(Unlist.test6)[1:15]<-vars.final

SampleTotal.final.Adj<-data.frame(matrix(nrow=nrow(SampleTotal.final.pqrs),ncol=length(vars.final)))
names(SampleTotal.final.Adj)[1:15]<-vars.final

Unlist.test6<- data.frame((Season.Adj[1][[1]][[1]]),stringsAsFactors = FALSE)
Unlist.test6.1<- data.frame((Season.Adj[2][[1]][[1]]),stringsAsFactors = FALSE)

# for(i in 1:)



##############################

for(p in 1:length(samplepointnum.unique)){
  sample.point.p<-SampleTotal.final[SampleTotal.final$samplepointnum == samplepointnum.unique[p],]

for(q in 1:length(samplepoint.month.unique)){
  sample.month.p<-sample.point.p[sample.point.p$sample.month == samplepoint.month.unique[q],]

for(r in 1:length(vars.final)){
  sample.mean.var<-mean(sample.month.p[,vars.final[[r]]])
  
  test.sample.rows<-which(SampleTotal.final$samplepointnum==samplepointnum.unique[p] & SampleTotal.final$sample.month ==  samplepoint.month.unique[q])
  
for(s in 1:length(test.sample.rows)){
  SampleTotal.final.pqrs[test.sample.rows[s],vars.final[[r]]]<-(SampleTotal.final[s,vars.final[[r]]]-sample.mean.var)

}
  #print(paste("Adjusting variable",vars.final[r],"(Sample Month #",samplepoint.month.unique[q],")" ))
  }
  print(paste("Adjusted values for month",q,"(Sample Point #",samplepointnum.unique[p],")" ))
}
  percent<-(p/length(samplepointnum.unique))*100
  print(paste("Completed Adjusting for Sample Point",samplepointnum.unique[p],"","(",percent,"% of Total)" ))
}

write.csv(SampleTotal.final.pqrs,file="D:/DataAnalysis/SampleTotal_SeasonAdj.csv",row.names=FALSE)

##############################################



qic(tasmax_Mean,x=YearMonth,ylab="tasmax_Mean",xlab="YearMonth",data=test.sample,chart="i")
qic(tasmax_Mean,ylab="tasmax_Mean",xlab="Date",data=test.sample,chart="i")
qic(tasmax_Mean,ylab="tasmax_Mean",xlab="Date",data=test.sample,chart="c")
qic(tasmax_Mean,x=YearMonth,ylab="tasmax_Mean",xlab="Date",data=test.sample,chart="xbar")


xbar_chart1<-qcc(data=test.sample$tasmax_Max,type="xbar.one",plot=TRUE,digits=5)

xbar_chart1<-qcc(data=test.sample$tasmax_Max,type="g",plot=TRUE,digits=5,confidence.level = 0.95)


summary(xbar_chart1)

##############################################
setwd(basedirectory)
SampleTotal.SeasonAdj<-read.csv("SampleTotal_SeasonAdj.csv")
str(SampleTotal.SeasonAdj)

SampleTotal.final<-SampleTotal.SeasonAdj
############
xbar.samples<-data.frame(matrix(ncol=31,nrow=length(unique(SampleTotal.final$samplepointnum))))
colnames(xbar.samples)<-c("samplepointnum","tasmax_Mean.v.limits","tasmax_Mean.v.runs","tasmax_Max.v.limits","tasmax_Max.v.runs","tasmax_Min.v.limits","tasmax_Min.v.runs"
                          ,"tasmin_Mean.v.limits","tasmin_Mean.v.runs","tasmin_Max.v.limits","tasmin_Max.v.runs","tasmin_Min.v.limits","tasmin_Min.v.runs"
                          ,"huss_Mean.v.limits","huss_Mean.v.runs","huss_Max.v.limits","huss_Max.v.runs","huss_Min.v.limits","huss_Min.v.runs"
                          ,"pr_Mean.v.limits","pr_Mean.v.runs","pr_Max.v.limits","pr_Max.v.runs","pr_Min.v.limits","pr_Min.v.runs"
                          ,"rsds_Mean.v.limits","rsds_Mean.v.runs","rsds_Max.v.limits","rsds_Max.v.runs","rsds_Min.v.limits","rsds_Min.v.runs")
samplepointnum.unique<-unique(SampleTotal.final$samplepointnum)
vars.final<-names(SampleTotal.final)[6:20]

for(x in 1:length(samplepointnum.unique)){
  
  test.sample<-SampleTotal.final[SampleTotal.final$samplepointnum == samplepointnum.unique[x],]
  xbar.samples[x,1]<-samplepointnum.unique[x]
  
  for(z in 1:length(vars.final)){
    z.var<-z+5
    xbar_chart1<-qcc(data=test.sample[z.var],type="xbar.one",plot=FALSE,digits=5)
    sample.violation.limits<-length(xbar_chart1$violations$beyond.limits)
    sample.violation.runs<-length(xbar_chart1$violations$violating.runs)
    
    xbar.samples[x,(z*2)]<-sample.violation.limits
    xbar.samples[x,(z*2+1)]<-sample.violation.runs
  }
  percent<-(x/length(unique(SampleTotal.final$samplepointnum)))*100
  print(paste("Tested Limits for",samplepointnum.unique[x],"samplenum ","(",percent,"%)" ))
  }


str(xbar.samples)
summary(xbar.samples)


plot(xbar.samples$samplepointnum,xbar.samples[,3],ylim=c(0,80))
points(xbar.samples$samplepointnum,xbar.samples[,5],col="red")
points(xbar.samples$samplepointnum,xbar.samples[,7],col="blue")

plot(xbar.samples$samplepointnum,xbar.samples$tasmax_Mean.v.limits)

str(SampleTotal.final)
SamplePointNums<-SampleTotal.final[3:5]
table.d<-SamplePointNums[1:length(unique(SamplePointNums$samplepointnum)),]

Xbar_SampleTotal <- merge(xbar.samples,table.d,by="samplepointnum")
write.csv(Xbar_SampleTotal,file="D:/DataAnalysis/Xbar_SampleTotal_SeasonAdj.csv",row.names=FALSE)


################################################################################################################
################################################################################################################
# Convert both original and seasonally-adjusted datasets into a Multivariate Control Chart

setwd(basedirectory)
SampleTotal.final<-read.csv("SampleStack_ClimateTotal_Final.csv")
SampleTotal.SeasonAdj<-read.csv("SampleTotal_SeasonAdj.csv")
summary(SampleTotal.final)

#SampleTotal.final<-SampleTotal.SeasonAdj
SampleTotal.final$pr_Min<-NULL

#samplepointnum.unique<-unique(SampleTotal.final$samplepointnum)
#xbar.samples[x,1]<-samplepointnum.unique[x]

#xbar_chart1<-qcc(data=test.sample[z.var],type="xbar.one",plot=FALSE,digits=5)
# sample.violation.limits<-length(xbar_chart1$violations$beyond.limits)
# sample.violation.runs<-length(xbar_chart1$violations$violating.runs)
# cor(test.sample[,6:19])
# test.sample.1<-SampleTotal.final[SampleTotal.final$samplepointnum == samplepointnum.unique[x],]
# test.sample<-test.sample.1[,c(6,7,8,9,10,11,12,13,14,15,16,17,18,19)]
# test.sample<-test.sample.1[,6:19]
# T2_single_chart1<-mqcc(test.sample,type="T2.single",limits=TRUE, pred.limits = FALSE)
# sample.violation.limits<-length(T2_single_chart1$violations$beyond.limits)


############
mqcc.samples<-data.frame(matrix(ncol=2,nrow=length(unique(SampleTotal.final$samplepointnum))))

samplepointnum.unique<-unique(SampleTotal.final$samplepointnum)
vars.final<-names(SampleTotal.final)[6:19]

for(x in 1:length(samplepointnum.unique)){
  
  test.sample<-SampleTotal.final[SampleTotal.final$samplepointnum == samplepointnum.unique[x],]
  mqcc.samples[x,1]<-samplepointnum.unique[x]
  
  test.sample.x<-test.sample[,6:19]
  
    T2_single_chart1<-mqcc(test.sample.x,type="T2.single",limits=TRUE, pred.limits = FALSE)
    sample.violation.limits<-length(T2_single_chart1$violations$beyond.limits)
    
    mqcc.samples[x,2]<-sample.violation.limits

  percent<-(x/length(unique(SampleTotal.final$samplepointnum)))*100
  print(paste("Tested Limits for",samplepointnum.unique[x],"samplenum ","(",percent,"%)" ))
}


str(mqcc.samples)
summary(mqcc.samples)
names(mqcc.samples)[1]<-"samplepointnum"
names(mqcc.samples)[2]<-"multivar.v.limits"

SamplePointNums<-SampleTotal.final[3:5]
table.d<-SamplePointNums[1:length(unique(SamplePointNums$samplepointnum)),]
MQCC_SampleTotal <- merge(table.d,mqcc.samples,by="samplepointnum")

#write.csv(MQCC_SampleTotal,file="D:/DataAnalysis/MultiQCC_SampleTotal.csv",row.names=FALSE)
write.csv(MQCC_SampleTotal,file="D:/DataAnalysis/MultiQCC_SampleTotal_SeasonAdj.csv",row.names=FALSE)


########################################################

setwd(basedirectory)
SampleTotal.final<-read.csv("SampleStack_ClimateTotal_Final.csv")

summary(SampleTotal.final)
cor(SampleTotal.final[6:20])
pairs(SampleTotal.final[6:20])
str(SampleTotal.final)

SampleTotal.final.2<-subset(SampleTotal.final,select=c(tasmax_Max,tasmin_Min,rsds_Mean))
str(SampleTotal.final.2)

pca=prcomp(SampleTotal.final.2)

par(mar=rep(2,4))
plot(pca)

pca$rotation=-pca$rotation
pca$x=pca$x
biplot(pca,scale=0)

################################################################################################################
################################################################################################################
#Import and process Raster Data
setwd(EnvirData_path)

NED.raster<-raster("AppLCC_ned_30_proj.tif")
#proj4string( data.r )
#shift.raster <- raster::shift(data.r,-360)
#NED.projected<-projectRaster( NED.raster , crs = proj4string( countryshape ),method = "bilinear" )
#NED.proj.r<-raster::crop( NED.raster , extent( countryshape ) )

#Apply the terrain functions to the DEM dataset

Slope.NED<-terrain(NED.raster,opt='slope',unit='degrees',neighbors=8,filename = "Slope_ned_proj.tif")
removeTmpFiles()
gc()
Aspect.NED<-terrain(NED.raster,opt='aspect',unit='degrees',neighbors=8,filename = "Aspect_ned_proj.tif")
removeTmpFiles()
gc()
TPI.NED<-terrain(NED.raster,opt='TPI',filename = "TPI_ned_proj.tif")
removeTmpFiles()
gc()
TRI.NED<-terrain(NED.raster,opt='TRI',filename = "TRI_ned_proj.tif")
removeTmpFiles()
gc()
Roughness.NED<-terrain(NED.raster,opt='roughness',filename = "Roughness_ned_proj.tif")
removeTmpFiles()
gc()


Slope.raster<-raster("Slope_ned_proj.tif")
Aspect.raster<-raster("Aspect_ned_proj.tif")
TPI.raster<-raster("TPI_ned_proj.tif")
TRI.raster<-raster("TRI_ned_proj.tif")
Roughness.raster<-raster("Roughness_ned_proj.tif")


samp_strat.points<-sample_points
samplepoints.df<-as.data.frame(coordinates(samp_strat.points))
samplepointnum<-c(1:nrow(samplepoints.df))
samplepoints.df<-cbind(samplepointnum,samplepoints.df)


NED.extract<-raster::extract(NED.raster,samp_strat.points)
Slope.extract<-raster::extract(Slope.raster,samp_strat.points)
Aspect.extract<-raster::extract(Aspect.raster,samp_strat.points)
TPI.extract<-raster::extract(TPI.raster,samp_strat.points)
TRI.extract<-raster::extract(TRI.raster,samp_strat.points)
Roughness.extract<-raster::extract(Roughness.raster,samp_strat.points)


##################
# Start the clock!
ptm <- proc.time()
##################
NEDfocal <- focal(NED.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max, filename='NED_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
##################
#Check the time elapsed
proc.time() - ptm
##################
##################
# Start the clock!
ptm <- proc.time()
##################
NEDfocal <- focal(NED.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max, filename='NED_focal_3x3max_na_rmF.tif',na.rm=FALSE) 
removeTmpFiles()
gc()
##################
#Check the time elapsed
proc.time() - ptm
##################
##################
# Start the clock!
ptm <- proc.time()
##################
NEDfocal <- focal(NED.raster, w=matrix(1/21025,nrow=145,ncol=145),fun=mean, filename='NED_focal_145x145mean_na_rmF.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
##################
#Check the time elapsed
proc.time() - ptm
##################
Slopefocal <- focal(Slope.raster, w=matrix(1/9,nrow=3,ncol=3),fun=mean,  filename='Slope_focal_3x3mean.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
Aspectfocal <- focal(Aspect.raster, w=matrix(1/9,nrow=3,ncol=3),fun=mean,  filename='Aspect_focal_3x3mean.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
TPIfocal <- focal(TPI.raster, w=matrix(1/9,nrow=3,ncol=3),fun=mean,  filename='TPI_focal_3x3mean.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
TRIfocal <- focal(TRI.raster, w=matrix(1/9,nrow=3,ncol=3),fun=mean,  filename='TRI_focal_3x3mean.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
Roughnessfocal <- focal(Roughness.raster, w=matrix(1/9,nrow=3,ncol=3),fun=mean,  filename='Roughness_focal_3x3mean.tif',na.rm=TRUE) 
removeTmpFiles()
gc()

Slopefocal <- focal(Slope.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max,  filename='Slope_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
Aspectfocal <- focal(Aspect.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max,  filename='Aspect_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
TPIfocal <- focal(TPI.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max,  filename='TPI_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
TRIfocal <- focal(TRI.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max,  filename='TRI_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()
Roughnessfocal <- focal(Roughness.raster, w=matrix(1/9,nrow=3,ncol=3),fun=max,  filename='Roughness_focal_3x3max.tif',na.rm=TRUE) 
removeTmpFiles()
gc()

# NEDfocal<-raster("NED_focal_3x3mean.tif")
# NEDfocal<-raster("NED_focal_5x5mean.tif")
# Slopefocal<-raster("Slope_focal_5x5mean.tif")
# Aspectfocal<-raster("Aspect_focal_5x5mean.tif")
# TPIfocal<-raster("TPI_focal_5x5mean.tif")
# TRIfocal<-raster("TRI_focal_5x5mean.tif")
# Roughnessfocal<-raster("Roughness_focal_5x5mean.tif")

setwd(EnvirData_path)

NEDfocal<-raster("NEDFocal_150x150_mean.tif")
Slopefocal<-raster("SlopeFocal2_150x150_mean.tif")
Aspectfocal<-raster("AspectFocal_150x150_mean.tif")
TPIfocal<-raster("focaltpi_150")
TRIfocal<-raster("focaltri_150")
Roughnessfocal<-raster("focalrou_150")



samp_strat.points<-sample_points
samplepoints.df<-as.data.frame(coordinates(samp_strat.points))
samplepointnum<-c(1:nrow(samplepoints.df))
samplepoints.df<-cbind(samplepointnum,samplepoints.df)


NED.extract<-raster::extract(NEDfocal,samp_strat.points)
Slope.extract<-raster::extract(Slopefocal,samp_strat.points)
Aspect.extract<-raster::extract(Aspectfocal,samp_strat.points)
TPI.extract<-raster::extract(TPIfocal,samp_strat.points)
TRI.extract<-raster::extract(TRIfocal,samp_strat.points)
Roughness.extract<-raster::extract(Roughnessfocal,samp_strat.points)


NED.extract.df<-as.data.frame(NED.extract)
Slope.extract.df<-as.data.frame(Slope.extract)
Aspect.extract.df<-as.data.frame(Aspect.extract)
TPI.extract.df<-as.data.frame(TPI.extract)
TRI.extract.df<-as.data.frame(TRI.extract)
Roughness.extract.df<-as.data.frame(Roughness.extract)


#sample.names<-rep(layer.t.name,nrow(samplepoints.df))
#sample.month<-rep(month.t,nrow(samplepoints.df))
#sample.year<-rep(year.t,nrow(samplepoints.df))
NED.samples<-cbind(samplepoints.df,NED.extract.df,Slope.extract.df,Aspect.extract.df,TPI.extract.df,TRI.extract.df,Roughness.extract.df)
str(NED.samples)

# setwd(basedirectory)
# 
# Xbar_SampleTotal<-read.csv("Xbar_SampleTotal.csv")
# Xbar_SampleTotal<-read.csv("Xbar_SampleTotal_SeasonAdj.csv")
# str(Xbar_SampleTotal)
# 
# Xbar_SampleTotal1 <- merge(NED.samples,Xbar_SampleTotal,by="samplepointnum")

#write.csv(Xbar_SampleTotal1,file="D:/DataAnalysis/Topo_Xbar_SampleTotal.csv",row.names=FALSE)
#write.csv(Xbar_SampleTotal1,file="D:/DataAnalysis/Topo_Xbar_SampleTotal_SeasonAdj.csv",row.names=FALSE)

setwd(basedirectory)
MQCC_SampleTotal<-read.csv("MultiQCC_SampleTotal.csv")
MQCC_SampleTotal<-read.csv("MultiQCC_SampleTotal_SeasonAdj.csv")
str(MQCC_SampleTotal)

MQCC_SampleTotal1 <- merge(NED.samples,MQCC_SampleTotal,by=c("samplepointnum","coords.x1","coords.x2"))
MQCC_SampleTotal1<- MQCC_SampleTotal1[order(MQCC_SampleTotal1[,1]), ]

names(MQCC_SampleTotal1)[2]<-"x1"
names(MQCC_SampleTotal1)[3]<-"x2"
names(MQCC_SampleTotal1)[4]<-"NEDextract"
names(MQCC_SampleTotal1)[5]<-"Slopeextract"
names(MQCC_SampleTotal1)[6]<-"Aspectextract"
names(MQCC_SampleTotal1)[7]<-"TPIextract"
names(MQCC_SampleTotal1)[8]<-"TRIextract"
names(MQCC_SampleTotal1)[9]<-"Roughnessextract"
names(MQCC_SampleTotal1)[10]<-"multivarvlimits"

PointNums<-read.csv("AppLCC_Points.csv")
App.points<-PointNums[,1]

selectedRows <- (MQCC_SampleTotal1$samplepointnum %in% App.points)

MQCC_SampleTotal.App <- MQCC_SampleTotal1[selectedRows,]

MQCC_SampleTotal.App<-MQCC_SampleTotal.App[complete.cases(MQCC_SampleTotal.App),]

#write.csv(MQCC_SampleTotal.App,file="D:/DataAnalysis/Topo_MQCC_Focal_SampleTotal.csv",row.names=FALSE)
#write.csv(MQCC_SampleTotal.App,file="D:/DataAnalysis/Topo_MQCC__Focal_SampleTotal_SeasonAdj.csv",row.names=FALSE)

#write.csv(MQCC_SampleTotal,file="D:/DataAnalysis/MQCC_SampleTotal_complete.csv",row.names=FALSE)


########################################################

setwd(basedirectory)
MQCC_SampleTotal<-read.csv("Topo_MQCC_Focal_SampleTotal_SeasonAdj.csv")
str(MQCC_SampleTotal)

var.col<-ncol(MQCC_SampleTotal)
names(MQCC_SampleTotal)

cor(MQCC_SampleTotal[,4:9])

#Simple Linear Regression
lm.out<-lm(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract, data=MQCC_SampleTotal)
lm.out<-lm(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract+x1, data=MQCC_SampleTotal)
summary(lm.out)

plot(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract, data=MQCC_SampleTotal,main="Topo Climate Plot")
abline(lm.out,col="red")

par(mfrow=c(2,2))
plot(lm.out)

#Poisson Regression
glm.out<-glm(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract, data=MQCC_SampleTotal,family=poisson())
summary(glm.out)

coefficients(glm.out)
#fitted(glm.out)
vcov(glm.out)
#influence(glm.out)

lm.out<-lm(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract, data=MQCC_SampleTotal)

glm.out<-glm(multivarvlimits ~ NEDextract+Slopeextract+Aspectextract+TPIextract, data=MQCC_SampleTotal,family=poisson())

step.out<-stepAIC(glm.out,direction="backward")
step.out$anova

calc.relimp(lm.out,type=c("lmg","last","first","pratt"),rela=TRUE)


str(MQCC_SampleTotal)
names(MQCC_SampleTotal)

#############################
#Convert data to geodata form
MQCC_Sample.data<-as.geodata(MQCC_SampleTotal,coords.col = 2:3,data.col = var.col,covar.col = 4:7)

#MCMC SIMULATION
model2 <- list(cov.pars = c(1, 1), beta = 1, family = "poisson")

mcmc2.tune <- mcmc.control(S.scale = 0.21, thin = 1)
test2.tune <- glsm.mcmc(MQCC_Sample.data, model = model2, mcmc.input = mcmc2.tune)


test2.tune.c <- create.mcmc.coda(test2.tune, mcmc.input = mcmc2.tune)
test2.tune.c <- create.mcmc.coda(test2.tune$simulations[45,], mcmc.input = list(S.scale = 0.5, thin = 1))
par(mfrow = c(1, 2))
plot(test2.tune.c, density = FALSE, ask = FALSE, auto.layout = FALSE)
autocorr.plot(test2.tune.c, ask = FALSE, auto.layout = FALSE)

mcmc2 <- mcmc.control(S.scale = 0.21)
test2 <- glsm.mcmc(MQCC_Sample.data, model = model2, mcmc.input = mcmc2)
summary(test2)
test2$acc.rate
test2$model

#SPATIAL PREDICTION
out2 <- output.glm.control(sim.predict = TRUE)
pred.test2 <- glsm.krige(test2, locations = cbind(as.numeric(MQCC_SampleTotal[10,2:3]), as.numeric(MQCC_SampleTotal[2,2:3])), output = out2)
cbind(pred.test2$predict, pred.test2$mcmc.error)


# BAYESIAN ANALYSIS
prior5 <- prior.glm.control(phi.prior = "fixed", phi = 0.1)
mcmc5.tune <- mcmc.control(S.scale = 0.001, thin = 10)
test5.tune <- pois.krige.bayes(MQCC_Sample.data, prior = prior5, mcmc.input = mcmc5.tune)

as.numeric(test5.tune$posterior$acc.rate[,2])

# test.scale.i<-seq(0.001,0.999,0.001)
# Acc.test.list<-list(length=0)
# for(i in 1:length(test.scale.i)){
#   prior5 <- prior.glm.control(phi.prior = "fixed", phi = 0.1)
#   mcmc5.tune <- mcmc.control(S.scale = test.scale.i[i], thin = 10)
#   test5.tune <- pois.krige.bayes(MQCC_Sample.data, prior = prior5, mcmc.input = mcmc5.tune)
#   Acc.rate.i<-as.numeric(test5.tune$posterior$acc.rate[,2])
#   if(any(Acc.rate.i>0)){
#     Acc.test.list<-c(Acc.test.list,test.scale.i[i])
#   }
# }
# as.numeric(test5.tune$posterior$acc.rate[,2])

# Prediction of model?
mcmc5 <- mcmc.control(S.scale = 0.001, thin = 100)
out5 <- output.glm.control(threshold = 10, quantile = c(0.05,0.99))
test5 <- pois.krige.bayes(MQCC_Sample.data, locations = t(cbind(as.numeric(MQCC_SampleTotal[10,2:3]), as.numeric(MQCC_SampleTotal[2,2:3]))), prior = prior5, mcmc.input = mcmc5,output = out5)

cbind(test5$predict, test5$mcmc.error)

mcmc6.tune <- mcmc.control(S.scale = 0.001, n.iter = 2000,thin = 100, phi.scale = 0.01)
prior6 <- prior.glm.control(phi.prior = "uniform", phi.discrete = seq(0.02, 1, 0.02), tausq.rel = 0.05)
test6.tune <- pois.krige.bayes(MQCC_Sample.data, prior = prior6, mcmc.input = mcmc6.tune)

mcmc6 <- mcmc.control(S.scale = 0.001, n.iter = 4e+05, thin = 200, burn.in = 5000, phi.scale = 0.12, phi.start = 0.5)
test6 <- pois.krige.bayes(MQCC_Sample.data, locations = t(cbind(as.numeric(MQCC_SampleTotal[10,2:3]), as.numeric(MQCC_SampleTotal[2,2:3]))), prior = prior6, mcmc.input = mcmc6)

test6$predict
test6$krige.var
test6$mcmc.error

par(mfrow = c(1,3))
hist(test6$posterior$beta$sample, main ="beta")
hist(test6$posterior$sigmasq$sample, main = "sigmasq")
hist(test6$posterior$phi$sample, main = "phi")

sqrt(asympvar(test6$posterior$beta$sample)/2000) 
sqrt(asympvar(test6$posterior$sigmasq$sample)/2000) 
sqrt(asympvar(test6$posterior$phi$sample)/2000) 
sqrt(asympvar(log(test6$posterior$simulations))/2000) 

#test.coda <- mcmc(t(log(test6$intensity)),thin = 1)

#############################

setwd(basedirectory)
MQCC_SampleTotal<-read.csv("Topo_MQCC_Focal_SampleTotal_SeasonAdj.csv")
str(MQCC_SampleTotal)

var.col<-ncol(MQCC_SampleTotal)
names(MQCC_SampleTotal)

summary(MQCC_SampleTotal)

MQCC_SampleTotal$test.vars<-MQCC_SampleTotal$multivarvlimits
Third.Qu.x<-as.numeric(quantile(MQCC_SampleTotal$test.vars)[4])

for(x in 1:nrow(MQCC_SampleTotal)){
  test.var.x<-MQCC_SampleTotal$test.vars[x]
  if(test.var.x>=Third.Qu.x){
    MQCC_SampleTotal$test.vars[x]<-1
  }else{
    MQCC_SampleTotal$test.vars[x]<-0
  }
}


Model.test.1<-lmer(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract + (1|test.vars), data = MQCC_SampleTotal, REML = FALSE)
Model.test.2<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, data = MQCC_SampleTotal,method="REML")
Model.test.3<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corAR1(form=~1),data = MQCC_SampleTotal,method="REML")
Model.test.4<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corLin(form=~1),data = MQCC_SampleTotal,method="REML")


summary(Model.test.3)
coef(Model.test.2)
coef(Model.test.3)
coef(Model.test.4)
r.squaredGLMM(Model.test.1)
AICc(Model.test.1)
AICc(Model.test.2)
AICc(Model.test.3)
AICc(Model.test.4)
anova(Model.test.1)
anova(Model.test.2)
anova(Model.test.3)
anova(Model.test.4)

anova(Model.test.2,Model.test.3,Model.test.4)

#############################
# Using GLMM package

Model.test.2<-glmm(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, 
                   random = test.vars,varcomps.names = c(1,0),data=MQCC_SampleTotal,family.glmm = bernoulli,debug=TRUE)

#############################
#############################
# Iterating through Control Chart parameters


setwd(basedirectory)
SampleTotal.final<-read.csv("SampleStack_ClimateTotal_Final.csv")
SampleTotal.SeasonAdj<-read.csv("SampleTotal_SeasonAdj.csv")
summary(SampleTotal.final)

SampleTotal.final<-SampleTotal.SeasonAdj
SampleTotal.final$pr_Min<-NULL

vars.total.all<-names(SampleTotal.final)[6:19]
############

iter.comb<-c(11:14)

samplepointnum.unique<-unique(SampleTotal.final$samplepointnum)
SamplePointNums<-SampleTotal.final[3:5]
table.d<-SamplePointNums[1:length(unique(SamplePointNums$samplepointnum)),]
PointNums<-read.csv("AppLCC_Points.csv")
App.points<-PointNums[,1]


for(i in 1:length(iter.comb)){
  
  vars.list<-combn(names(SampleTotal.final)[6:19],iter.comb[i])
  
  for(c in 1:ncol(vars.list)){

    mqcc.samples<-data.frame(matrix(ncol=2,nrow=length(unique(SampleTotal.final$samplepointnum))))
    vars.final<-vars.list[,c]
    
    for(x in 1:length(samplepointnum.unique)){
      
      test.sample<-SampleTotal.final[SampleTotal.final$samplepointnum == samplepointnum.unique[x],]
      mqcc.samples[x,1]<-samplepointnum.unique[x]
      
      test.sample.x<-test.sample[,vars.final]
      
      T2_single_chart1<-mqcc(test.sample.x,type="T2.single",limits=TRUE, pred.limits = FALSE)
      sample.violation.limits<-length(T2_single_chart1$violations$beyond.limits)
      
      mqcc.samples[x,2]<-sample.violation.limits
      
      names(mqcc.samples)[1]<-"samplepointnum"
      names(mqcc.samples)[2]<-"multivarvlimits"
      MQCC_SampleTotal <- merge(table.d,mqcc.samples,by="samplepointnum")
      selectedRows <- (MQCC_SampleTotal$samplepointnum %in% App.points)
      MQCC_SampleTotal.App <- MQCC_SampleTotal[selectedRows,]
      MQCC_SampleTotal.App<-MQCC_SampleTotal.App[complete.cases(MQCC_SampleTotal.App),]
      # percent<-(x/length(unique(SampleTotal.final$samplepointnum)))*100
      # print(paste("Tested Limits for",samplepointnum.unique[x],"samplenum ","(",percent,"%)" ))
    
      index.nums<-which(vars.total.all %in% vars.final)
      index.coll<-paste(index.nums,collapse="~")
      
      filename.f<-paste0("D:/DataAnalysis/MultiVarLimit_Iter_SeasonAdj/MQCC_Limits_",index.coll,".csv")
      write.csv(MQCC_SampleTotal.App,file=filename.f,row.names=FALSE)
      }
    ############
    print(paste("Measured Limits for",c," iteration out of ",ncol(vars.list),"(Combination Length",iter.comb[i],")" ))
    print(paste("Saved to:",filename.f))
  }
}





str(mqcc.samples)
summary(mqcc.samples)


#setwd(basedirectory)
MQCC_SampleTotal<-read.csv("Topo_MQCC_Focal_SampleTotal_SeasonAdj.csv")
NED.samples<-MQCC_SampleTotal[,-10]

setwd(CombinIter_SeasAdj_path)
csv_file_list <- list.files(pattern = ".csv" , all.files = FALSE , full.names = FALSE )

Model.compare.AIC<-data.frame(matrix(ncol=5,nrow=length(csv_file_list)))
names(Model.compare.AIC)<-c("Model.num","Model.vars","Model1.AIC","Model2.AIC","Model3.AIC")

# Model.compare.AIC<-data.frame(matrix(ncol=4,nrow=length(csv_file_list)))
# names(Model.compare.AIC)<-c("Model.num","Model.vars","Model1.AIC","Model2.AIC")

for(i in 1:length(csv_file_list)){
  
  mqcc.samples<-read.csv(csv_file_list[i])
  names(mqcc.samples)[2]<-"x1"
  names(mqcc.samples)[3]<-"x2"
  
  MQCC_SampleTotal.merge <- merge(NED.samples,mqcc.samples,by=c("samplepointnum","x1","x2"))
  MQCC_SampleTotal.merge<- MQCC_SampleTotal.merge[order(MQCC_SampleTotal.merge[,1]), ]
  
  Model.AIC.1<-9999
  Model.AIC.2<-9999
  Model.AIC.3<-9999
  
  Model.base.1<-NA
  Model.corAR1.2<-NA
  Model.corLIN.3<-NA
  
  if(min(MQCC_SampleTotal.merge$multivarvlimits)>0){
    
    tryCatch(   Model.base.1<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, 
                                  data = MQCC_SampleTotal.merge,method="REML"),error=function(e) {NA})
    tryCatch(   Model.corAR1.2<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corAR1(form=~1),
                                    data = MQCC_SampleTotal.merge,method="REML"),error=function(e) {NA})
    tryCatch(   Model.corLIN.3<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corLin(form=~ x1 + x2),
                                    data = MQCC_SampleTotal.merge,method="REML"),error=function(e) {NA})

    Model.AIC.1<-tryCatch(AICc(Model.base.1),error=function(e) {print(paste0("Error in Model 1 of iteration ",i));NA})
    Model.AIC.2<-tryCatch(AICc(Model.corAR1.2),error=function(e) {print(paste0("Error in Model 2 of iteration ",i));NA})
    Model.AIC.3<-tryCatch(AICc(Model.corLIN.3),error=function(e) {print(paste0("Error in Model 3 of iteration ",i));NA})
    
    if(!is.na(Model.AIC.1)){Model.compare.AIC[i,3]<-AICc(Model.base.1)}
    if(!is.na(Model.AIC.2)){Model.compare.AIC[i,4]<-AICc(Model.corAR1.2)}
    if(!is.na(Model.AIC.3)){Model.compare.AIC[i,5]<-AICc(Model.corLIN.3)}
    }
  
  filename.f<-csv_file_list[i]
  col.nums.char<-substr(strsplit(filename.f,"_")[[1]][3],1,nchar(strsplit(filename.f,"_")[[1]][3])-4)
  # col.nums.unique<-unique(na.omit(as.numeric(unlist(strsplit(unlist(col.nums.char),"[^0-9+]")))))
  # vars.result<-vars.total.all[col.nums.unique]
  
  Model.compare.AIC[i,1]<-i
  Model.compare.AIC[i,2]<-col.nums.char
  print(paste("Models tested for",i,"iteration out of",length(csv_file_list) ))
  print(Model.compare.AIC[i,])
}


summary(Model.compare.AIC)

#write.csv(Model.compare.AIC,file="D:/DataAnalysis/Model.compare.AIC_SeasonAdj.csv",row.names=FALSE)
Model.compare.AIC<-read.csv("D:/DataAnalysis/Model.compare.AIC_SeasonAdj.csv")

barplot(Model.compare.AIC[order(Model.compare.AIC$Model3.AIC), ]$Model3.AIC)

Model.compare.AIC.sort<- Model.compare.AIC[order(Model.compare.AIC$Model3.AIC), ]

top.100.AIC<-Model.compare.AIC.sort[1:100,2]

#col.nums.unique<-unique(na.omit(as.numeric(unlist(strsplit(unlist(as.character(top.100.AIC[1])),"[^0-9+]")))))

top.100.AIC.t <- vector(mode="numeric", length=0)

for(t in 1:length(top.100.AIC)){
  top.100.AIC.t<-c(top.100.AIC.t,unique(na.omit(as.numeric(unlist(strsplit(unlist(as.character(top.100.AIC[t])),"[^0-9+]"))))))
}

order(-count(top.100.AIC.t)$freq)

min(Model.compare.AIC[,3:5], na.rm=T)
min(Model.compare.AIC[,3], na.rm=T)
min(Model.compare.AIC[,4], na.rm=T)
min(Model.compare.AIC[,5], na.rm=T)

which.min(Model.compare.AIC[,3:5])
minOfRows=apply(Model.compare.AIC, 1, function(x) min(x[x!=0])) 

Model.compare.AIC[which.min(Model.compare.AIC$Model1.AIC),]
Model.compare.AIC[which.min(Model.compare.AIC$Model2.AIC),]
Model.compare.AIC[which.min(Model.compare.AIC$Model3.AIC),]

Model.compare.AIC.sort<- Model.compare.AIC[order(Model.compare.AIC$Model3.AIC), ]


col.nums.unique<-unique(na.omit(as.numeric(unlist(strsplit(unlist(Model.compare.AIC[420,2]),"[^0-9+]")))))
vars.result<-vars.total.all[col.nums.unique]

Vars.result_NoSeas_170<-c("tasmax_Mean","tasmax_Max","tasmax_Min","tasmin_Mean","huss_Max","huss_Min","pr_Mean","pr_Max","rsds_Mean","rsds_Max","rsds_Min")
#Vars.result_SeasAdj_405<-c("tasmax_Max","tasmax_Min","tasmin_Min","tasmin_Mean","tasmin_Max","huss_Mean","huss_Min","pr_Mean","pr_Max","rsds_Max","rsds_Min")
Vars.result_NoSeas_218<-c("tasmax_Mean","tasmax_Max","tasmax_Min","tasmin_Max","huss_Mean","huss_Max","huss_Min","pr_Mean","pr_Max","rsds_Max","rsds_Min")

#Final Model
Vars.result_NoSeas_420<-c("tasmax_Max","tasmax_Min","tasmin_Mean","tasmin_Max","huss_Mean","huss_Max","huss_Min","pr_Mean","pr_Max","rsds_Max","rsds_Min")

Vars.result_NoSeas_Final<-c("tasmax_Max","tasmax_Min","tasmin_Mean","tasmin_Max","huss_Mean","huss_Max","huss_Min","pr_Max","rsds_Mean","rsds_Max","rsds_Min")


index.nums<-which(vars.total.all %in% Vars.result_NoSeas_420)
index.nums<-which(vars.total.all %in% Vars.result_NoSeas_218)
index.nums<-which(vars.total.all %in% Vars.result_NoSeas_Final)

index.coll<-paste(index.nums,collapse="~")
filename.f<-paste0("D:/DataAnalysis/MultiVarLimit_Iter_SeasonAdj/MQCC_Limits_",index.coll,".csv")

filename.f<-paste0("D:/DataAnalysis/MultiVarLimit_Iter_SeasonAdj/MQCC_Limits_",index.coll,".csv")


mqcc.samples<-read.csv(filename.f)
names(mqcc.samples)[2]<-"x1"
names(mqcc.samples)[3]<-"x2"
MQCC_SampleTotal.merge <- merge(NED.samples,mqcc.samples,by=c("samplepointnum","x1","x2"))
MQCC_SampleTotal.merge<- MQCC_SampleTotal.merge[order(MQCC_SampleTotal.merge[,1]), ]

#write.csv(MQCC_SampleTotal.merge,file="/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/Topo.Model.AIC_SeasonAdj420.csv",row.names=FALSE)

Model.corAR1.2<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corAR1(form=~1),
                    data = MQCC_SampleTotal.merge,method="REML")
Model.corLIN.3<-gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corLin(form=~ x1 + x2),
                    data = MQCC_SampleTotal.merge,method="REML",verbose = TRUE)
Model.corLIN.GLs<-Gls(multivarvlimits ~ NEDextract + Slopeextract + Aspectextract + TPIextract, correlation=corLin(form=~ x1 + x2),data = MQCC_SampleTotal.merge,method="REML")
summary(Model.corAR1.2)
summary(Model.corLIN.3)
summary(Model.corLIN.GLs)

1-(0.4994/0.7251)

################################################################################################################
################################################################################################################

#Test of taking means

r.mean<- cellStats(stack.j[[1]],"mean",na.rm=TRUE)

if(Month.k %in% Winter){
  if(max.min.k == "tasmax"){
    Winter.k_max<-stack(Winter.k_max,stack.j[[k]])
    print(paste("Added layer",names(stack.j[[k]]),"to Winter_max Stack"))
  }else if(max.min.k=="tasmin"){
    Winter.k_min<-stack(Winter.k_min,stack.j[[k]])
    print(paste("Added layer",names(stack.j[[k]]),"to Winter_min Stack"))
  }
}



################################################################################################################
################################################################################################################


#1. Generate a Vector Grid of an appropriate Size over the study area
layer.t<-stack.g[[1]]
plot(layer.t)

layer.SPdf<-as(layer.t,"SpatialPointsDataFrame")
layer.SPdf.sample<-GSIF::sample.grid(layer.SPdf, cell.size = c((0.0625*10), (0.0625*10)), n = 1)

layer.grid<-layer.SPdf.sample$grid
plot(layer.grid)
layer.SGdf<-as(layer.t,"SpatialGridDataFrame")
layer.grid.raster<-raster(layer.grid)
pp <- rasterToPolygons(layer.grid.raster, dissolve=TRUE)
#writePolyShape(pp,"D:/DataAnalysis/Grid_Sample.shp")
plot(layer.t)
plot(pp,add=T)

#2. Generate random SpatialPoints within that sampling vector grid
layer.resamp<-resample(layer.grid.raster,layer.t,method='ngb')

samp_strat.rand<- sampleStratified(layer.resamp, 10, xy = TRUE, sp=TRUE, na.rm = TRUE)
samp_strat.points<-as(samp_strat.rand,"SpatialPoints")

plot(layer.t)
plot(pp,add=T)
points(samp_strat.points,col="black",pch=16,cex=0.5)
#points(samp_strat.rand)
#writePointsShape(samp_strat.rand,"D:/DataAnalysis/Point_Sample.shp")

#3. Extract values from raster layer based on SpatialPoints generated in the previous step
extract.sample<-raster::extract(layer.t,samp_strat.points)
extract.df<-as.data.frame(extract.sample)
samplepoints.df<-as.data.frame(coordinates(samp_strat.points))
sample.coordinates<-cbind(samplepoints.df,extract.df)


# r<-raster(ncol=30,nrow=20)
# r[]<-1:(30*20)
# #r[runif(30*20)>0.30]<-NA
# plot(r)
# SpatialPixels(SpatialPoints(coordinates(r)[!is.na(values(r)),]))
# 
# SP.r<-SpatialPoints(coordinates(r)[!is.na(values(r)),])
# proj4string(SP.r) <- CRS("+proj=longlat")
# class(SP.r)
# g<-as(r,"SpatialPointsDataFrame")
# p<-as(r,"SpatialPixels")
# plot(g)
# proj4string(g) <- CRS("+proj=longlat +datum=WGS84")
# class(g)



################################################################################################################
################################################################################################################
################################################################################################################

#Z-score calculation in R
##################
# Start the clock!
ptm <- proc.time()
##################
r.mean<- cellStats(stack.j[[1]],"mean",na.rm=TRUE)
r.sd<- cellStats(stack.j[[1]],"sd",na.rm=TRUE)
r.j<-stack.j[[1]]
r.j[is.na(r.j)] <- -9999
r.Z<-(r.j-r.mean)/r.sd
r.Z[r.Z < -500] <- NA
##################
#Check the time elapsed
proc.time() - ptm
##################

layer.t<-stack.g[[1]]
plot(layer.t)


#Convert values to Z-scores for entire stack:
stack.j<-stack.g
for(i in 1:nlayers(stack.j)){
  r.mean<- cellStats(stack.j[[i]],"mean",na.rm=TRUE)
  r.sd<- cellStats(stack.j[[i]],"sd",na.rm=TRUE)
  r.j<-stack.j[[i]]
  r.j[is.na(r.j)] <- -9999
  r.Z<-(r.j-r.mean)/r.sd
  r.Z[r.Z < -500] <- NA
  stack.j[[i]]<-r.Z
  print(paste("Converted layer",names(stack.j[[i]]),"in stack"))
}








#Generate six different grids of climate variability:

# Subset raster stack by season based on the month of the date since 1900

Winter<-c(12,1,2)
Summer<-c(6,7,8)
Winter.k_max<-stack()
Summer.k_max<-stack()
Winter.k_min<-stack()
Summer.k_min<-stack()
Dates.all<-vector(mode="character", length=0)

for(k in 1:nlayers(stack.j)){
  Date.k<-(substr(as.character(names(stack.j)[k]),nchar(names(stack.j)[k])-9,nchar(names(stack.j)[k])))
  Date.k<-gsub("[.]", "-", Date.k)
  Dates.all<-c(Dates.all,Date.k)
  Month.k<-as.numeric(format(as.Date(Date.k, origin = "1900-01-01"), "%m"))
  max.min.k<-(substr(as.character(names(stack.j)[k]),1,6))
  if(Month.k %in% Winter){
    if(max.min.k == "tasmax"){
      Winter.k_max<-stack(Winter.k_max,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Winter_max Stack"))
    }else if(max.min.k=="tasmin"){
      Winter.k_min<-stack(Winter.k_min,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Winter_min Stack"))
    }
  }
  if(Month.k %in% Summer){
    if(max.min.k == "tasmax"){
      Summer.k_max<-stack(Summer.k_max,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Summer_max Stack"))
    }else if(max.min.k=="tasmin"){
      Summer.k_min<-stack(Summer.k_min,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Summer_min Stack"))
    }
  }
  }

#SUBSET stack.j by YEAR to input into year calculation (MAX and MIN)

Dates.all

Years.all<-c(1950,1951)

Year1.k_max<-stack()
Year2.k_max<-stack()
Year1.k_min<-stack()
Year2.k_min<-stack()

for(k in 1:nlayers(stack.j)){
  Date.k<-(substr(as.character(names(stack.j)[k]),nchar(names(stack.j)[k])-9,nchar(names(stack.j)[k])))
  Date.k<-gsub("[.]", "-", Date.k)
  Year.k<-as.numeric(format(as.Date(Date.k, origin = "1900-01-01"), "%Y"))
  max.min.k<-(substr(as.character(names(stack.j)[k]),1,6))
  if(Year.k == Years.all[1]){
    if(max.min.k == "tasmax"){
      Year1.k_max<-stack(Year1.k_max,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Year1_max Stack"))
    }else if(max.min.k=="tasmin"){
      Year1.k_min<-stack(Year1.k_min,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Year1_min Stack"))
    }
  }
  if(Year.k == Years.all[2]){
    if(max.min.k == "tasmax"){
      Year2.k_max<-stack(Year2.k_max,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Year2_max Stack"))
    }else if(max.min.k=="tasmin"){
      Year2.k_min<-stack(Year2.k_min,stack.j[[k]])
      print(paste("Added layer",names(stack.j[[k]]),"to Year2_min Stack"))
    }
  }
}


#1) intra-seasonal variation in maximum temperatures, calculated as the 95th percentile of 
# summer (December-February) maximums minus the 5th percentile of summer maximums

Q95 <- calc(Summer.k_max, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q05 <- calc(Summer.k_max, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

Var_Summer_max<-(Q95-Q05)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Summer_max, filename="Var_Summer_max", format="GTiff", overwrite=TRUE)
Var_Summer_max<-raster("Var_Summer_max.tif")
setwd(outputdata_path)

plot(Q95)
plot(Q05)
plot(Var_Data_1)

#2) intra-annual variation in maximum temperatures, calculated as the 95th percentile of 
# winter maximum temperatures minus the 95th percentile of winter (June–August) maximum temperatures

Q95 <- calc(Winter.k_max, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q05 <- calc(Winter.k_max, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

Var_Winter_max<-(Q95-Q05)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Winter_max, filename="Var_Winter_max", format="GTiff", overwrite=TRUE)
Var_Winter_max<-raster("Var_Winter_max.tif")
setwd(outputdata_path)

#3) inter-annual variation in maximum temperatures, calculated as the difference 
# in the 95th percentile of maximum temperatures between the 2 years


Q95_y1 <- calc(Year1.k_max, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q95_y2 <- calc(Year2.k_max, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )

Var_Annual_max<-(Q95_y1-Q95_y2)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Annual_max, filename="Var_Annual_max", format="GTiff", overwrite=TRUE)
Var_Annual_max<-raster("Var_Annual_max.tif")
setwd(outputdata_path)

#4) intra-seasonal variation in minimum temperatures, calculated as the 95th percentile of 
# winter minimum temperatures minus the 5th percentile of winter minimum temperatures

Q95 <- calc(Summer.k_min, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q05 <- calc(Summer.k_min, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

Var_Summer_min<-(Q95-Q05)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Summer_min, filename="Var_Summer_min", format="GTiff", overwrite=TRUE)
Var_Summer_min<-raster("Var_Summer_min.tif")
setwd(outputdata_path)

#5) intra-annual variation in minimum tempera- tures, calculated as the 5th percentile of summer minimum
# temperatures minus the 5th percentile of winter minimum temperatures

Q95 <- calc(Winter.k_min, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q05 <- calc(Winter.k_min, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

Var_Winter_min<-(Q95-Q05)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Winter_min, filename="Var_Winter_min", format="GTiff", overwrite=TRUE)
Var_Winter_min<-raster("Var_Winter_min.tif")
setwd(outputdata_path)

#6) inter-annual variation in minimum temperatures, calculated as the difference in the 
# 5th percentile of minimum temperatures between the 2 years

Q95_y1 <- calc(Year1.k_min, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )
Q95_y2 <- calc(Year2.k_min, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

Var_Annual_min<-(Q95_y1-Q95_y2)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Annual_min, filename="Var_Annual_min", format="GTiff", overwrite=TRUE)
Var_Annual_min<-raster("Var_Annual_min.tif")
setwd(outputdata_path)

############################
#Generate Z-score analysis of data over study area

#1) #Z scores for the coldest/hottest locations:
# Calculated for both 5th percentile min temp and 95th percentile max temp

Q95_max <- calc(stack.j, fun = function(x) {quantile(x,probs = c(.95),na.rm=TRUE)} )
Q05_min <- calc(stack.j, fun = function(x) {quantile(x,probs = c(.05),na.rm=TRUE)} )

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Q95_max, filename="Total_Q95_max", format="GTiff", overwrite=TRUE)
Q95_max<-raster("Q95_max.tif")
setwd(outputdata_path)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Q05_min, filename="Total_Q05_min", format="GTiff", overwrite=TRUE)
Q05_min<-raster("Q05_min.tif")
setwd(outputdata_path)

#2) #Z scores for the most climatically stable locations:
# "Z-scores for three variability grids associated with 95th percentile 
# of maximum temperatures were averaged to estimate overall climate stability, 
# and same was done for three variability grids associated with 5th percentile of minimum temperatures"

Var_Q95_max <- ((Var_Summer_max+Var_Winter_max+Var_Annual_max)/3)
Var_Q05_min <- ((Var_Summer_min+Var_Winter_min+Var_Annual_min)/3)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Q95_max, filename="Var_Q95_max", format="GTiff", overwrite=TRUE)
Var_Q95_max<-raster("Var_Q95_max.tif")
setwd(outputdata_path)

setwd("/Users/alecnelson/Documents/Clemson_Files/MastersThesis/DataAnalysis/TestProcessingOutput")
#r.output <- writeRaster(Var_Q05_min, filename="Var_Q05_min", format="GTiff", overwrite=TRUE)
Var_Q05_min<-raster("Var_Q05_min.tif")
setwd(outputdata_path)

#3) #Z scores for the difference from the matrix:
# Calculated as difference between site temperature and average temperature 
  # within 5km radius moving window

test.j1<-stack.g[[1]]

cfw <- focalWeight(test.j1, 0.15, "circle")
r.mean <- focal(test.j1,  w=cfw, mean, na.rm=TRUE)
r.mean <- focal(test.j1, w=matrix(1/25, nc=5, nr=5))
r.diff<-(test.j1-r.mean)
plot(test.j1)
plot(r.mean)
plot(r.diff)

sign.matrix<-r.diff
sign.matrix[ sign.matrix > 0 ] <- 1
sign.matrix[ sign.matrix < 0 ] <- -1
plot(sign.matrix)



#Calculate final Refugia Index for Max and Min temperatures:

RI_max <- (Q95_max + (sign.matrix.max*Var_Q95_max) + Matrix_max)/3
RI_min <- (Q95_min + (sign.matrix.min*Var_Q95_min) + Matrix_min)/3



