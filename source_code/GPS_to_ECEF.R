#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2024/5/16.
if( dev.cur() > 1 ) dev.off()

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################


TIMEFORMAT=1
ORDER_longi_lati=1

options(digits=22)
basedir=dirname(here::here())
Day=c("Day1_Pre_drone","Day3_Post_drone")
FL=c("FL01","FL02","FL03")

dir_list=c()
for(i in Day){
  for(j in FL){
    d=paste0(basedir,"/csv/",i,"/",j)
    dir_list=append(dir_list,d)
  }
}


for(dir in dir_list){
  
  
  
  setwd(paste0(dir,"/original GPS data"))
  
  files  <- list.files(pattern=".csv", recursive = T)
  # files=subset(files,regexpr("position", files)<0)
  # files=subset(files,regexpr("xyz", files)<0)
  # files=subset(files,regexpr("#", files)<0)
  
  # Individual index is fixed by the "_#N_" token in each raw file name
  # (same bird = same number on every day; #1 is the drone-experienced bird,
  # whose name also keeps the "_drone" tag). Ordering by N makes list element
  # i always bird #i, independent of GPS-logger serial numbers.
  bird_no=as.integer(sub(".*_#([0-9]+)_.*","\\1",files))
  if(any(is.na(bird_no))||any(duplicated(bird_no))){
    stop("raw GPS file names must carry unique _#N_ tags: ",paste(files,collapse=", "))
  }
  files=files[order(bird_no)]
  if(!grepl("_drone",files[1])){
    stop("bird #1 must be the drone-experienced individual (_drone): ",files[1])
  }
  data=list()
  for (file.name in files) {
    # raw files carry a header: timestamp,latitude,longitude,altitude
    data<-c(data,list(read.csv(file.name,header=T)[,1:4]))
  }
  
  if(ORDER_longi_lati){
    for(i in 1:length(data)){
      data[[i]]=data[[i]][,c(1,3,2,4)]
    }
  }
  
  for(i in 1:length(data)){
    colnames(data[[i]])=c("HMS","Longitude","Latitude","Altitude")
  }
  
  result<-list()
  basepoint=as.numeric(data[[1]][1,2:4])
  
  for(i in 1:length(data)){
    xyz=transform(data[[i]],basepoint)
    set=data.frame(sec=addSubSec(data[[i]])$sec.double, xyz)
    result=c(result,list(set))
  }
  
  #smoothing&sync
  output_fps=5
  timerange=OverlappedTimeRange(result)
  xyz_sync=getKinematicInfo(result,timerange)
  savedir=paste0(dir,"/earth_frame data")
  dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
  save(xyz_sync, file=paste0(savedir,"/xyz_sync.Rdata"))
  
  result<-list()
  for(i in 1:length(data)){
    set<-data.frame(sec=addSubSec(data[[i]])$sec.double,data[[i]][,2:4])
    result=c(result,list(set))
  }
  timerange=OverlappedTimeRange(result)
  lonlat_sync=getKinematicInfo(result,timerange,1)
  for(i in 1:length(lonlat_sync)){
    colnames(lonlat_sync[[i]]$x)=c("sec","Longitude","Latitude","Altitude")
  }
  save(lonlat_sync, file=paste0(savedir,"/lonlat_sync.Rdata"))
}
