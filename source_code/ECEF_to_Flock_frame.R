#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2022/12/01.


if( dev.cur() > 1 ) dev.off()
#############################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################


output_dpi=300
im_size=5
arrow_size=0.5
arrowhead_size=0.01

ElapsedTime=T
vlim=15

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
  setwd(dir)
  
  load("./earth_frame data/xyz_sync.Rdata")
  load("./earth_frame data/lonlat_sync.Rdata")
  
  
  
  if(ElapsedTime){
    basetime=xyz_sync[[1]]$x$sec[1]
    xyz_sync=resetTime(xyz_sync,basetime)
    lonlat_sync=resetTime(lonlat_sync,basetime)
  }
  
  
  xyz_sync2=addMean(xyz_sync)
  lonlat_sync2=addMean(lonlat_sync)
  
  objects=resetNorm(xyz_sync2)
  objects=getFormationSD(objects)
  rot_objects=getLocalcoordinate(objects)[[1]]
  
  object=c(list(objects),list(rot_objects), list(lonlat_sync2))
  names(object)=c("ori","rot","ori_geo")
  
  # add distance to the center of #2-#4, "dist_center_of_three"
  center_flockmates=getMean(object$rot[2:4])
  object$rot=addDistance(object$rot,center_flockmates$x,"dist_center_of_three")
  
  # #5 data update 
  # $ave.dist_center_of_all is the average distance of individuals #2–#4 to the flock center  
  # $ave.dist_center_of_three is the average distance of individuals #2–#4 to the center of #2-#4

  object$rot[[5]]$x=data.frame(object$rot[[5]]$x,
                               ave.dist_center_of_all=getMean(object$rot[2:4])$x$norm,
                               ave.dist_center_of_three=getMean(object$rot[2:4])$x$dist_center_of_three)
    
    
  # save(object, file="./transformed data/final product/parameters_no_crop.Rdata")
  
  
  
  #select data in cruise, faster than a flight velocity threshold
  object_crop=crop_longestDuration(object,vlim)
  
  # check_croprange(list(object,object_crop),"./visualized data/croprange.png")
  parts <- strsplit(dir, "/")[[1]]
  day <- strsplit(parts[length(parts)-1], "_")[[1]][1]
  prop <- paste0(day,"_",tail(parts, 1))
  savedir=paste0(dir,"/flock_frame data/")
  dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
  save(object, object_crop, file=paste0(savedir,"flockframe_",prop,".Rdata"))
  for(individual in 1:length(object_crop$rot)){
    suffix=NULL
    if(individual==length(object_crop$rot)){
      suffix="(centerpoint)"
    }
    write.csv(object_crop$rot[[individual]]$x,paste0(savedir,"flockframe_",prop,"_#",individual,suffix,".csv"))
  }
  
  filepath=paste0(basedir,"/derived_Rdata/")
  dir.create(filepath, recursive = TRUE, showWarnings = FALSE)
  save(object, object_crop, file=paste0(filepath,"flockframe_",prop,".Rdata"))
  
  
}

