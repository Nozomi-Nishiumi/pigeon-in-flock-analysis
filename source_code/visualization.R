#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2024/5/16.
if( dev.cur() > 1 ) dev.off()



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

# ---- Shared density scale for the 2D distribution plots --------------------
# geom_density2d() by default picks both its KDE bandwidth and its contour
# levels from each series alone, so a contour line means a different density
# in every panel. To make the distribution plots visually comparable across
# days AND across flocks, every series is smoothed with one fixed common
# bandwidth, and the contours are probability-mass isopleths (same stay-time
# fractions for every bird), so region size directly encodes dispersion.
# The KDE-mode markers in draw2d() use the same shared bandwidth, keeping
# them inside the innermost isopleth.
h_shared=c(6.5,6.5)

# Probability-mass isopleths: each contour encloses the given fraction of a
# bird's positions (5-80%: the 5% line hugs the density peak, the 80% line
# shows the outer extent).
mass_probs=c(0.05,0.25,0.45,0.65,0.8)

list_F_B=list()
list_H_L=list()

for(dir in dir_list){
  setwd(dir)
  parts <- strsplit(dir, "/")[[1]]
  day <- strsplit(parts[length(parts)-1], "_")[[1]][1]
  prop <- paste0(day,"_",tail(parts, 1))
  load(paste0("./flock_frame data/","flockframe_",prop,".Rdata"))
    obj_num=length(object$ori_geo)
  if(length(object$rot)!=obj_num){
    print("num of objects is not consistent")
    exit()
  }
    
    
    
    dir.create("visualized data", showWarnings = FALSE)
    localScene(300,7000,object_crop$ori,"./visualized data/entirescene")
    localScene(300,7000,object_crop$rot,"./visualized data/formation")
    draw2d(object_crop$rot,c("x","y"),15,F,T,F,F,NULL,"./visualized data/distribution_xy.png","backON",h=h_shared,mass_probs=mass_probs)
    draw2d(object_crop$rot,c("x","z"),15,F,T,F,F,NULL,"./visualized data/distribution_xz.png","backON",h=h_shared,mass_probs=mass_probs)

    savedir=paste0(basedir,"/figures/distribution")
    dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
    draw2d(object_crop$rot,c("x","y"),15,F,T,F,F,NULL,paste0(savedir,"/distribution_xy_",prop,".png"),"backOFF",h=h_shared,mass_probs=mass_probs)
    draw2d(object_crop$rot,c("x","z"),15,F,T,F,F,NULL,paste0(savedir,"/distribution_xz_",prop,".png"),"backOFF",h=h_shared,mass_probs=mass_probs)
    
    
    
  
  cp=object$ori[[length(object$ori)]]$x
  path_list=list()
  for (i in 1:obj_num) {
    geom_path(
      aes(x = Longitude, y = Latitude), color =hsv(110/180,0.5,1,alpha=0.1), linewidth = 1, data = object$ori_geo[[i]]$x
    )|> blend("screen")->path
    
    path_list=c(path_list,list(path))
  }
  
  axis_F_B=object$rot[[1]]$x$y
  axis_H_L=object$rot[[1]]$x$z
  gps=object$ori_geo[[1]]$x
  
  point_F_B=coloring2(axis_F_B,gps,6,c(31/360,180/360))
  point_H_L=coloring2(axis_H_L,gps,25,c(0/360,240/360))
  list_F_B=c(list_F_B,list(point_F_B))
  list_H_L=c(list_H_L,list(point_H_L))
}

all<-list(list_F_B,list_H_L)
names(all)=c("FB","HL")

centerplace=c(-83.384415,33.792843)
drone_engage=c(-83.3928,33.8147)
drone_disengage=c(-83.39208221,33.82188416)
engage_dist=calc_xy(drone_disengage[1],drone_disengage[2],
                      drone_engage[1],drone_engage[2])
engage_radius=sqrt(engage_dist$x^2+engage_dist$y^2)
engage_midpoint=(drone_engage+drone_disengage)/2

theta=seq(0,2*pi, length.out=100)
r=engage_radius
x <- r*sin(theta)
y <- r*cos(theta)
drone_area <- calc_lambda(x, y, drone_engage)

engage_p<-geom_point( aes(x = drone_engage[1], y = drone_engage[2]),
                      shape = 21, size = 1, stroke = 1, 
                      colour = hsv(1, 1, 1), 
                      fill = hsv(0,0, 1, alpha=0.08))

interception_area<-geom_polygon( aes(x = lambda_deg, y = phi_deg), data=drone_area, 
                     colour = hsv(0, 0, 1,alpha=1), linewidth=0.8,
                     fill = hsv(0,0, 1, alpha=0.12))


scalebar=scalebar_gmap(map,0.8,0.03,1000,"1000 m")

savedir=paste0(basedir,"/figures/trajectories")
dir.create(savedir, recursive = TRUE, showWarnings = FALSE)

for(j in 1:length(all)){
for(i in 1:length(FL)){
ggsave(paste0(savedir,"/map_",names(all)[j],"_FL0",i,".png"), map+interception_area+scalebar+all[[j]][[i]]+all[[j]][[i+3]],
       bg = "transparent",dpi=400,width=5,height=5)
}
}






