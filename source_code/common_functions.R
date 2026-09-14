if( dev.cur() > 1 ) dev.off()

targetPackages <- c(
  "tidyverse",
  "ggplot2", 
  "ggmap",
  "RgoogleMaps",
  "ggblend",
  "Rmpfr",
  
  "rgl",
  "lubridate",
  "stringr",
  "ggthemes",
  "ggExtra",
  "ggblend",
  "here",
  
  "rmarkdown",
  "knitr",
  "MASS"
) 

newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")

# rgl in headless (null) mode: drawScene()/localScene() only need to export
# HTML via rglwidget(); we do not want an OpenGL window (avoids X11/GLX
# errors on macOS without XQuartz, and lets the pipeline run in RStudio /
# Rscript without a display). Must be set BEFORE library(rgl) below.
options(rgl.useNULL = TRUE)

for(package in targetPackages) library(package, character.only = T)



coloring<-function(forward_back,gps,range){
  colorlist=c(31/360,180/360)
  point_list=list()
  high_low=object$rot[[1]]$x$z
  hue_i=(forward_back<1)+1
  sat=abs(forward_back)/range
  for(k in seq(1, nrow(gps), length = 400)){
    s=min(max(sat[k],0),0.9)
    v=min(max(sat[k],0.9),1)
    geom_point(
      aes(x = Longitude, y = Latitude), color =hsv(colorlist[hue_i[k]],s,1,alpha=0.8), size = 1, data = gps[k,]
    )->point
    
    point_list=c(point_list,list(point))
  }
  return(point_list)
}

coloring2<-function(forward_back,gps,range,colorlist){
  point_list=list()
  high_low=object$rot[[1]]$x$z
  hue_i=(forward_back<1)+1
  hue=colorlist[hue_i]
  sat=abs(forward_back)/range
  for(s in 1:length(sat)){
    sat[s]=min(max(sat[s],0),1)
  }
  visualset=data.frame(gps,hue,sat)
  head(visualset)
  
  geom_point(
    aes(x = Longitude, y = Latitude), color =hsv(visualset$hue,visualset$sat,1,alpha=0.8), size = 1, data = visualset
  )->testpoint
  
  return(testpoint)
}

CurrentTime<-function(srcHMS,srcMS=0){
  if(TIMEFORMAT==0){
    timedata=hms(srcHMS)
  }else{
    timedata=ymd_hms(srcHMS)
  }
  srcMS=round(srcMS)
  curr_time=hour(timedata)*60*60+minute(timedata)*60+second(timedata)+srcMS/1000
  base_time=hour(timedata)[1]*60*60+minute(timedata)[1]*60+second(timedata)[1]+srcMS[1]/1000
  elapsed_time=curr_time-base_time
  return(round(curr_time,digits=5))
}

addSubSec<-function(d){
  sec=CurrentTime(d$HMS)
  d=cbind(d,sec=CurrentTime(d$HMS))
  n_distinct(sec)
  timegroup=group_by(d,sec)
  summary=summarise(timegroup, n = length(HMS))
  
  sec.double=numeric()
  for(i in 1:nrow(summary)){
    for(j in 1:summary$n[i]){
      subsec=summary$sec[i]+1/summary$n[i]*(j-1)
      sec.double=c(sec.double,subsec)
    }
  }
  
  d=cbind(d,sec.double)
  return(d)
}

dataOmitting<-function(d){
  data_omit=d[1,]
  for(i in 1:nrow(d)){
    prev=data_omit[nrow(data_omit),]
    curr=d[i,]
    diff=(prev-curr)[1,c("Longitude","Latitude","Altitude")]
    
    if(diff$Longitude==0&&diff$Latitude==0&&diff$Altitude==0){
      next
    }
    data_omit=rbind(data_omit,curr)
  }
  rownames(data_omit)=NULL
  return(data_omit)
}


# 補助関数
A_array<-function(n){
  A0 = 1 + (n^2)/4.0+ (n^4)/64.0
  A1 = -(3.0/2)*( n - (n^3)/8.0- (n^5)/64.0) 
  A2 =  (15.0/16)*( n^2 - (n^4)/4.0)
  A3 = -(35.0/48)*( n^3 - (5.0/16)*(n^5) )
  A4 =  (315.0/512)*( n^4 )
  A5 = -(693.0/1280)*( n^5 )
  return(c(A0, A1, A2, A3, A4, A5))
}

alpha_array<-function(n){
  a0 = 0 # dummy
  a1 = (1.0/2)*n - (2.0/3)*(n^2) + (5.0/16)*(n^3) + (41.0/180)*(n^4) - (127.0/288)*(n^5)
  a2 = (13.0/48)*(n^2) - (3.0/5)*(n^3) + (557.0/1440)*(n^4) + (281.0/630)*(n^5)
  a3 = (61.0/240)*(n^3) - (103.0/140)*(n^4) + (15061.0/26880)*(n^5)
  a4 = (49561.0/161280)*(n^4) - (179.0/168)*(n^5)
  a5 = (34729.0/80640)*(n^5)
  return(c(a0, a1, a2, a3, a4, a5))
}

calc_xy<-function(lambda_deg,phi_deg, lambda0_deg,phi0_deg){
  
  phi_rad = phi_deg*pi/180
  lambda_rad = lambda_deg*pi/180
  phi0_rad = phi0_deg*pi/180
  lambda0_rad = lambda0_deg*pi/180
  
  
  # 定数 (a, F: 世界測地系-測地基準系1980（GRS80）楕円体)
  m0 = 0.9999 
  a = 6378137.0
  f = 298.257222101
  
  # (1) n, A_i, alpha_iの計算
  n = 1.0/ (2*f - 1)
  A_array = A_array(n)
  alpha_array = alpha_array(n)
  
  # (2), S, Aの計算
  A_ = ( (m0*a)/(1.0+n) )*A_array[1] # [m]
  S_ = ( (m0*a)/(1.0+n) )*( A_array[1]*phi0_rad + as.numeric(A_array[2:6]%*%sin(2*phi0_rad*c(1:5))) ) # [m]
  
  # (3) lambda_c, lambda_sの計算
  lambda_c = cos(lambda_rad - lambda0_rad)
  lambda_s = sin(lambda_rad - lambda0_rad)
  
  # (4) t, t_の計算
  t = sinh( atanh(sin(phi_rad)) - ((2*sqrt(n)) / (1+n))*atanh(((2*sqrt(n)) / (1+n)) * sin(phi_rad)) )
  t_ = sqrt(1 + t*t)
  
  # (5) xi', eta'の計算
  xi2  = atan(t / lambda_c) # [rad]
  eta2 = atanh(lambda_s / t_)
  
  # (6) x, yの計算
  y = mapply(function(xi2,eta2){A_ * (xi2 + sum(alpha_array[2:6]*(sin(2*xi2*c(1:5)))*(cosh(2*eta2*c(1:5))))) - S_ },
             xi2,eta2)
  x = mapply(function(xi2,eta2){A_ * (eta2 + sum(alpha_array[2:6]*(cos(2*xi2*c(1:5)))*(sinh(2*eta2*c(1:5)))))},
             xi2,eta2) 
  # return
  return(data.frame(x, y)) # [m]
}

# 逆展開係数 β（標準的な5次までの係数）
beta_array <- function(n) {
  b0 = 0
  b1 = (1.0/2)*n - (2.0/3)*(n^2) + (37.0/96)*(n^3) - (1.0/360)*(n^4) - (81.0/512)*(n^5)
  b2 = (1.0/48)*(n^2) + (1.0/15)*(n^3) - (437.0/1440)*(n^4) + (46.0/105)*(n^5)
  b3 = (17.0/480)*(n^3) - (37.0/840)*(n^4) - (209.0/4480)*(n^5)
  b4 = (4397.0/161280)*(n^4) - (11.0/504)*(n^5)
  b5 = (4583.0/161280)*(n^5)
  return(c(b0, b1, b2, b3, b4, b5))
}


# --- 逆変換（修正版） ---
inverse_transform <- function(xyz, base) {
  latlon = inverse_calc_xy(xyz$x, xyz$y, base[1], base[2])
  data = data.frame(
    Latitude = latlon$Latitude,
    Longitude = latlon$Longitude,
    Altitude = xyz$z + base[3]
  )
  return(data)
}

inverse_calc_xy <- function(x, y, lambda0_deg, phi0_deg) {
  # 定数
  m0 = 0.9999 
  a = 6378137.0
  f = 298.257222101
  n = 1.0 / (2*f - 1)
  
  A_array_val = A_array(n)
  alpha_array_val = alpha_array(n)
  beta_array_val = beta_array(n)
  
  phi0_rad = phi0_deg * pi / 180
  lambda0_rad = lambda0_deg * pi / 180
  
  A_ = (m0*a)/(1.0+n) * A_array_val[1]
  S_ = (m0*a)/(1.0+n) * (A_array_val[1]*phi0_rad + sum(A_array_val[2:6]*sin(2*phi0_rad*(1:5))))
  
  N = length(x)
  Lat = numeric(N)
  Lon = numeric(N)
  term = (2*sqrt(n)) / (1+n)
  
  for (i in seq_len(N)) {
    xi_val = (y[i] + S_) / A_   # = ξ (forwardのξ)
    eta_val = x[i] / A_        # = η (forwardのη)
    js = 1:5
    # ξ' (xi_prime) と η' (eta_prime) を求める（逆補正）
    xi_prime = xi_val - sum( beta_array_val[js+1] * sin(2*js*xi_val) * cosh(2*js*eta_val) )
    eta_prime = eta_val - sum( beta_array_val[js+1] * cos(2*js*xi_val) * sinh(2*js*eta_val) )
    
    # chi（conformal latitude）
    chi = asin( sin(xi_prime) / cosh(eta_prime) )
    
    # chi -> phi をニュートン法で解く:
    #    atanh(sin(phi)) - term*atanh(term*sin(phi)) = atanh(sin(chi))
    phi = chi  # 初期値
    rhs = atanh( sin(chi) )
    for (k in 1:50) {
      fval = atanh( sin(phi) ) - term * atanh( term * sin(phi) ) - rhs
      # f' = 1/cos(phi) - term^2 * cos(phi) / (1 - term^2 * sin(phi)^2)
      fp = 1.0 / cos(phi) - (term^2) * cos(phi) / (1.0 - (term^2) * sin(phi)^2)
      delta = fval / fp
      phi = phi - delta
      if (abs(delta) < 1e-12) break
    }
    
    lambda_rad = lambda0_rad + atan2( sinh(eta_prime), cos(xi_prime) )
    
    Lat[i] = phi * 180.0 / pi
    Lon[i] = lambda_rad * 180.0 / pi
  }
  
  return(data.frame(Latitude = Lat, Longitude = Lon))
}

GPStoXY<-function(d){
  base=d[1,]
  ang=90-base$latitude
  
  
  longi=d$Longitude*pi/180
  lati=d$Latitude*pi/180
  base_lati=lati[1]
  
  alpha=6378137
  e=0.081819191042815791
  x=alpha*cos(base_lati)/sqrt(1-e^2*sin(base_lati)^2)*longi
  y=alpha*(1-e^2)/((1-e^2*sin(base_lati)^2)^(3/2))*lati
  return(data.frame(x,y))
}

smoothed_position<-function(positiondata,func, output_fps, time_range=NULL){
  if(is.null(time_range)){
    time_range=c(min(positiondata$sec), max(positiondata$sec))
  }
  TIME=seq(time_range[1],time_range[2], by = 1/output_fps)
  result=data.frame(TIME)
  for(funcindex in 1:3){
    predicted <- predict(func[[funcindex]],TIME, deriv = 0)$y
    result=cbind(result,predicted)
  }
  
  colnames(result)=c("sec","x","y","z")
  return(result)
}

diff.calibration<-function(src,base){
  result=src[,2:4]-base[,2:4]
  return(result)
}

transform<-function(data,base){
  xy=calc_xy(data$Longitude,data$Latitude,base[1],base[2])
  xyz=data.frame(x=xy$x,y=xy$y,z=data$Altitude-base[3])
  return(xyz)
}

horizontal_distance<-function(data,base){
  xy=calc_xy(data$Longitude,data$Latitude,base[1],base[2])
  xyz=data.frame(x=xy$x,y=xy$y,0)
  dist=getNorm(xyz)
  return(dist)
}

getMotionFunc<-function(xyz){
  motionFunc=list()
  for(funcindex in 1:(length(xyz)-1)){
    func_each <- smooth.spline(xyz$sec, xyz[,funcindex+1])
    motionFunc=c(motionFunc,list(func_each))
  }
  return(motionFunc)
}

OverlappedTimeRange<-function(src){
  timeset=data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
  for(i in 1:length(src)){
    trange=c(min(src[[i]]$sec),max(src[[i]]$sec))
    timeset=rbind(timeset, trange)
    colnames(timeset)=c("init","end")
  }
  
  timerange=c(max(timeset$init),min(timeset$end))
  return(timerange)
}


drawScene<-function(data,fov=60,res=0,name="scene"){
  
  
  plotdata=list()
  for(i in 1:length(data)){
    #switch axis
    df=data.frame(x=data[[i]]$y,
                  y=data[[i]]$z,
                  z=data[[i]]$x,
                  sec=data[[i]]$sec)
    
    #set resolution
    if(res!=0){
      df=df[seq(1,nrow(data[[i]]),length.out=res),]
    }
    plotdata=c(plotdata,list(df))
  }
  
  
  data=plotdata
  num=length(data)
  
  range<-data.frame(matrix(rep(NA, length(2)), ncol=3))[,numeric(0)]
  
  for(i in 1:3){
    maxvalue=numeric()
    minvalue=numeric()
    for(j in 1:num){
      maxvalue=c(maxvalue,max(na.omit(data[[j]][,i])))
      minvalue=c(minvalue,min(na.omit(data[[j]][,i])))
    }
    r=c(min(minvalue),max(maxvalue))
    range=cbind(range,r)
  }
  d_omit<-as.data.frame(range[is.finite(rowSums(range)),])
  if(nrow(d_omit)==0){
    print("ERROR")
    return()
  }
  colnames(range)=c("x","y","z")
  
  diff=numeric()
  for(i in 1:3){
    diff=c(diff,range[2,i]-range[1,i])  
  }
  
  ratio=c(diff[1]/max(diff),diff[2]/max(diff),diff[3]/max(diff))
  
  
  
  
  if(rgl.cur()!=0){
    close3d()
  }
  plot3d(x=c(0),y=c(0),z=c(0),type="n",
         xlim=range$x,ylim=range$y,zlim=range$z,
         xlab="y",ylab="z",zlab="x")
  par3d(windowRect = c(0, 0, 700, 700))
  aspect3d(ratio)
  
  t_range=max(data[[1]]$sec)-min(data[[1]]$sec)
  candidatestep=c(0.01,0.1,1,10,60)
  timedif=abs(candidatestep-t_range/5)
  timestep=candidatestep[which.min(timedif)]
  
  
  
  for(i in 1:num){
    col_list=huelisting(num)
    sat_list=c(rep(0.65,num-1),0)
    val_list=c(rep(0.8,num-1),0.6)
    colour = hsv(col_list[i], sat_list[i], val_list[i])
    
    lines3d(data[[i]]$x,data[[i]]$y,data[[i]]$z,col=colour)
    start=(min(data[[i]]$sec)%/%timestep)*timestep
    ext=seq(start, max(data[[i]]$sec), by = timestep)
    ext_index <- sapply(ext, function(z)which.min(abs(z - data[[i]]$sec)))
    text3d(data[[i]]$x[ext_index],data[[i]]$y[ext_index],data[[i]]$z[ext_index]+diff[1]/30, texts=ext, col=colour,alpha=0.2)
  }
  
  bgplot3d({plot.new(); title(name)})
  view3d(theta=-90, phi=90, fov=fov, zoom=0.8, scale=par3d("scale"), interactive=TRUE)
  
  # setwd("./graphics/scene/")
  HTML <- rglwidget( width=1000,height=600)
  fullfilepath=paste0(name,".html")
  htmlwidgets::saveWidget(HTML, fullfilepath,selfcontained = T)
  # setwd("../../")
  close3d()
} 

getKinematicInfo<-function(objects,t_range=NULL, depth=3){
  
  res=list()
  for(i in 1:length(objects)){
    obj=objects[[i]]
    func=getMotionFunc(obj)
    if(!is.null(t_range)){
      TIME=seq(t_range[1],t_range[2], by = 1/output_fps)
    }else{
      TIME=obj$sec
    }
    
    
    
    x=data.frame(TIME)
    for(funcindex in 1:length(func)){
      predicted <- predict(func[[funcindex]],TIME, deriv = 0)$y
      x=cbind(x,predicted)
    }
    # norm=getNorm(x[,2:4])
    # x=cbind(x,norm)
    colnames(x)=c("sec","x","y","z")  
    object=list(x=x)
    
    if(depth>1){
      v=data.frame(TIME)
      for(funcindex in 1:3){
        predicted <- predict(func[[funcindex]],TIME, deriv = 1)$y
        v=cbind(v,predicted)
      }
      # norm=getNorm(v[,2:4])
      # v=cbind(v,norm)
      colnames(v)=c("sec","x","y","z")  
      object=c(object,list(v=v))
    }
    
    if(depth>2){
      a=data.frame(TIME)
      for(funcindex in 1:3){
        predicted <- predict(func[[funcindex]],TIME, deriv = 2)$y
        a=cbind(a,predicted)
      }
      # norm=getNorm(a[,2:4])
      # a=cbind(a,norm)
      colnames(a)=c("sec","x","y","z")  
      object=c(object,list(a=a))
    }
    
    res=c(res,list(object))
  }
  return(res)
}




####


data2obj<-function(dataset){
  obj_num=(ncol(dataset)-1)/3
  objects=c()
  for(i in 1:obj_num){
    obj=data.frame(sec=dataset[,1],x=dataset[,(i-1)*3+2],y=dataset[,(i-1)*3+3],z=dataset[,(i-1)*3+4])
    objects=c(objects,list(obj))
  }
  return(objects)
}

localScene<-function(time,range,objects,name="scene",res=0){
  obj_num=length(objects)
  
  data=list()
  for(i in 1:obj_num){
    df=data.frame(sec=objects[[i]]$x$sec,
                  x=objects[[i]]$x$x,
                  y=objects[[i]]$x$y,
                  z=objects[[i]]$x$z)
    
    if(res!=0){
      df=df[seq(1,nrow(df),length.out=res),]
    }
    timedif=abs(df$sec-time)
    timestep=which.min(timedif)
    
    loc=objects[[obj_num]]$x[timestep,]
    area=data.frame(x=c(loc$x-range,loc$x+range),
                    y=c(loc$y-range,loc$y+range),
                    z=c(loc$z-range,loc$z+range))
    area_frame=which(df$x>area$x[1]&df$x<area$x[2]&
                       df$y>area$y[1]&df$y<area$y[2]&
                       df$z>area$z[1]&df$z<area$z[2])
    
    if(length(area_frame)>0){
      df=df[area_frame,]
      data=c(data,list(df))
    }
  }
  
  drawScene(data,60,0,name)
}



Steel.Dwass<- function(data,group){
  OK <- complete.cases(data, group)
  data <- data[OK]
  group <- group[OK]
  n.i <- table(group)
  ng <- length(n.i)
  t <- combn(ng, 2, function(ij) {
    i <- ij[1]
    j <- ij[2]
    r <- rank(c(data[group == i], data[group == j]))
    R <- sum(r[1:n.i[i]])
    N <- n.i[i]+n.i[j]
    E <- n.i[i]*(N+1)/2
    V <- n.i[i]*n.i[j]/(N*(N-1))*(sum(r^2)-N*(N+1)^2/4)
    return(abs(R-E)/sqrt(V))
  })
  p <- ptukey(t*sqrt(2), ng, Inf, lower.tail=FALSE)
  result <- cbind(t, p)
  rownames(result) <- combn(ng, 2, paste, collapse=":")
  return(result)
}

getCenterpoint<-function(objects){
  CP=data.frame(matrix(rep(0, ncol(objects[[1]])*nrow(objects[[1]])), nrow=nrow(objects[[1]])))
  for(i in 1:length(objects)){
    CP=CP+objects[[i]]
  }
  CP=CP/length(objects)
  colnames(CP)=c("sec","x","y","z") 
  objects=c(objects,list(CP))
  return(objects)
}

getMean<-function(objects){
  
  obj_num=length(objects)
  param_num=length(objects[[1]])
  
  names(objects[[1]])
  
  centerpoint=list()
  for(i in 1:param_num){
    CP=data.frame(matrix(rep(0, ncol(objects[[1]][[i]])*nrow(objects[[1]][[i]])), nrow=nrow(objects[[1]][[1]])))
    for(j in 1:obj_num){
      CP=CP+objects[[j]][[i]]
    }
    CP=CP/obj_num
    colnames(CP)=colnames(objects[[1]][[i]])
    centerpoint=c(centerpoint,list(CP))
    
  }
  names(centerpoint)=names(objects[[1]])
  return(centerpoint)
}

addMean<-function(src, add=src){
  
  obj_num=length(add)
  param_num=length(add[[1]])
  
  names(add[[1]])
  
  centerpoint=list()
  for(i in 1:param_num){
    CP=data.frame(matrix(rep(0, ncol(add[[1]][[i]])*nrow(add[[1]][[i]])), nrow=nrow(add[[1]][[i]])))
    for(j in 1:obj_num){
      CP=CP+add[[j]][[i]]
    }
    CP=CP/obj_num
    colnames(CP)=colnames(add[[1]][[i]])
    centerpoint=c(centerpoint,list(CP))
    
  }
  names(centerpoint)=names(add[[1]])
  return(c(src,list(centerpoint)))
}

resetNorm<-function(objects){
  for(i in 1:length(objects)){
    for(j in 1:length(objects[[1]])){
      df=objects[[i]][[j]][,2:4]
      norm=getNorm(df)
      objects[[i]][[j]]=cbind(objects[[i]][[j]][,-5],norm)
    }
  }
  return(objects)
}

resetTime<-function(objects,subset_t){
  for(i in 1:length(objects)){
    for(j in 1:length(objects[[i]])){
      objects[[i]][[j]]$sec=objects[[i]][[j]]$sec-subset_t
    }
  }
  return(objects)
}

huelisting<-function(num){
  colis=numeric()
  for(i in 0:(num-1)){
    colis=c(colis,i/num)
  }
  return(colis)
}

# Common density scale for draw2d(): one KDE bandwidth and one set of contour
# levels shared by every series, so that the same contour line means the same
# density in every panel (fair comparison across days and flocks).
#   series_list - list of data.frames holding the columns dim1, dim2
#   Returns list(h=c(h1,h2), breaks=contour levels)
sharedDensityScale<-function(series_list,dim1="x",dim2="y",n_breaks=10){
  ok=sapply(series_list,function(d)sd(d[[dim1]])>0.0001&sd(d[[dim2]])>0.0001)
  series_list=series_list[ok]
  # median of the per-series normal-reference bandwidths: a representative
  # smoothing level that is not inflated by between-flock spread (pooling
  # all series into one bandwidth.nrd() call would oversmooth)
  h=c(median(sapply(series_list,function(d)bandwidth.nrd(d[[dim1]]))),
      median(sapply(series_list,function(d)bandwidth.nrd(d[[dim2]]))))
  # the highest density any series reaches under the shared bandwidth sets
  # the top of the common contour ladder
  maxd=0
  for(d in series_list){
    k=kde2d(d[[dim1]],d[[dim2]],h=h,n=100)
    maxd=max(maxd,max(k$z))
  }
  breaks=pretty(c(0,maxd),n_breaks)
  breaks=breaks[breaks>0&breaks<maxd]
  return(list(h=h,breaks=breaks))
}

# Density levels enclosing given probability masses for one kde2d() result.
# For each p in probs, returns the density value whose superlevel set
# (highest-density region) contains fraction p of the total mass.
massToDensityLevels<-function(k,probs){
  cell=diff(k$x[1:2])*diff(k$y[1:2])
  z=sort(as.vector(k$z),decreasing=T)
  cm=cumsum(z)*cell
  sapply(probs,function(p){
    hit=which(cm>=p)
    if(length(hit)==0)return(min(z))
    z[hit[1]]
  })
}

# axes: which flock-frame coordinates go on the horizontal and vertical axis,
# e.g. c("x","y") (top view) or c("x","z") (side view). Positions and
# velocities are projected explicitly onto these two columns as a1/a2 (no
# column relabelling), so a1/a2 always mean "plotted axis 1 / axis 2".
draw2d<-function(objects,axes=c("x","y"),range,path=F,dnsty=F,mn=F,drawcp=F,time=NULL,filename=NULL,backOFF="backON",h=NULL,breaks=NULL,mass_probs=NULL){
  if(length(axes)!=2||!all(axes%in%c("x","y","z"))||axes[1]==axes[2]){
    stop("axes must be two distinct names among x, y, z")
  }
  project<-function(d)data.frame(sec=d$sec,a1=d[[axes[1]]],a2=d[[axes[2]]])
  
  obj_num=length(objects)
  timedif=abs(objects[[obj_num]]$x$sec-time)
  timestep=which.min(timedif)
  
  if(!drawcp){
    objects=objects[-length(objects)]
  }
  plist=list()
  slist=list()
  llist=list()
  dlist=list()
  mlist=list()
  mslist=list()
  
  loop=1
  col_list=huelisting(obj_num)
  sat_list=c(rep(0.65,obj_num-1),0)
  sat_list=c(0.65,rep(0,obj_num-1))
  val_list=c(rep(0.8,obj_num-1),0.6)
  val_list=c(0.8,rep(0.6,obj_num-1))
  
  draw<-ggplot()
  
  for(i in 1:length(objects)){
    pos=project(objects[[i]]$x)   # positions on the plotted axes
    vel=project(objects[[i]]$v)   # velocities on the plotted axes
    
    line<-geom_path(aes(x = a1, y = a2), 
                    colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=0.5),data = pos)
    llist=c(llist,list(line))
    
    point<-geom_point(aes(x = a1, y = a2), 
                      colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=1),size=5, 
                      data=pos[timestep,])
    plist<-c(plist,list(point))
    
    m=data.frame(a1=mean(pos$a1),a2=mean(pos$a2))
    mp<-geom_point(aes(x = a1, y = a2), 
                   colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=1),size=2, 
                   data=m)
    mlist<-c(mlist,list(mp))
    
    segment=data.frame(p1=pos$a1,p2=pos$a2,v1=vel$a1,v2=vel$a2)
    seg<-geom_segment(aes(x = p1, y = p2, xend = p1+v1/3, yend = p2+v2/3), 
                      arrow = arrow(length = unit(0.05, "inches"),type = "closed"), linejoin = "mitre",
                      colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=1),data = segment[timestep,])
    slist<-c(slist,list(seg))
    
    if(!is.null(h)){
      # position mode from the SAME shared-bandwidth KDE as the isopleths,
      # so the marker always falls inside the innermost contour
      x.top=max.location_KDE(pos$a1,pos$a2,h=h,
              lims=c(min(pos$a1)-3*h[1],max(pos$a1)+3*h[1],
                     min(pos$a2)-3*h[2],max(pos$a2)+3*h[2]))
    }else{
      x.top=max.location_KDE(pos$a1,pos$a2)
    }
    v.top=max.location_KDE(vel$a1,vel$a2)
    top=data.frame(p1=x.top$x,p2=x.top$y,v1=v.top$x,v2=v.top$y)
    mseg<-geom_segment(aes(x = p1, y = p2, xend = p1+v1/100, yend = p2+v2/100), 
                       arrow = arrow(length = unit(arrowhead_size, "inches")), linejoin = "mitre", linewidth=arrow_size,
                       colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=1), data=top)
    mslist<-c(mslist,list(mseg))
    
    if(sd(pos$a1)>0.0001){
      if(!is.null(mass_probs)&&!is.null(h)){
        # probability-mass isopleths: each line encloses the fraction of this
        # bird's positions given in mass_probs (e.g. 0.9 = smallest region
        # holding 90% of its time). The bandwidth h is shared across panels,
        # so region sizes are directly comparable between days and flocks.
        k=kde2d(pos$a1,pos$a2,h=h,n=100,
                lims=c(min(pos$a1)-3*h[1],max(pos$a1)+3*h[1],
                       min(pos$a2)-3*h[2],max(pos$a2)+3*h[2]))
        grid_df=expand.grid(a1=k$x,a2=k$y)
        grid_df$dens=as.vector(k$z)
        density<-geom_contour(aes(x=a1,y=a2,z=dens),
                              breaks=unique(massToDensityLevels(k,mass_probs)),
                              colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=0.5),linewidth=0.5,data=grid_df)
      }else if(is.null(h)||is.null(breaks)){
        # per-series automatic bandwidth and contour levels (not comparable
        # between panels)
        density<-geom_density2d(aes(x=a1, y=a2),
                                colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=0.5),size=0.5,data=pos)
      }else{
        # shared scale from sharedDensityScale(): identical contour lines
        # represent identical densities in every panel
        density<-geom_density2d(aes(x=a1, y=a2),
                                colour = hsv(col_list[i], sat_list[i], val_list[i], alpha=0.5),size=0.5,data=pos,
                                h=h,breaks=breaks)
      }
      dlist=c(dlist,list(density))
    }
    
    loop=loop+1
  }
  
  
  draw<-draw+geom_hline(yintercept=0,alpha=0.3)+geom_vline(xintercept=0,alpha=0.2)
  
  if(dnsty){draw<-draw+dlist}
  if(path){draw<-draw+llist}
  if(mn){draw<-draw+mlist}
  
  if(!is.null(time)){
    draw<-draw+plist+slist
    loc=project(objects[[length(objects)]]$x)[timestep,]
  }else{
    draw<-draw+mslist
    loc=data.frame(a1=c(0),a2=c(0))
  }
  
  area=data.frame(a1=c(loc$a1-range,loc$a1+range),
                  a2=c(loc$a2-range,loc$a2+range))
  draw<-draw+coord_cartesian(xlim = area$a1, ylim = area$a2)
  draw<-draw+xlab(paste0(axes[1]," (m)"))+ylab(paste0(axes[2]," (m)"))
  
  
  draw<-draw+theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
  
  if(!is.null(filename)){
    if(backOFF=="backOFF"){
      draw<-draw+theme(
        panel.background = element_rect(fill = "transparent",color = hsv(0,0,0, alpha=0.5)),
        panel.grid.minor = element_line(color = NA),
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA)
        # axis.text= element_blank(),
        # plot.title = element_blank(),
        # axis.title = element_blank()
      )
    }
    
    ggsave(filename, draw, bg = "transparent",dpi=output_dpi,width=im_size,height=im_size)
  }
  
  print(draw)
  
}

getRvector<-function(objects, ref){
  
  obj_num=length(objects)
  res=list()
  for(i in 1:(obj_num)){
    obj=objects[[i]]$x
    R=obj-ref
    R$sec=ref$sec
    add=c(objects[[i]],r=list(R))
    res=c(res,list(add))
  }
  return(res)
}

addDistance<-function(objects, ref, name="dist"){
  
  obj_num=length(objects)
  res=list()
  for(i in 1:(obj_num)){
    obj=objects[[i]]$x
    R=obj-ref
    dist=getNorm(data.frame(R$x,R$y,R$z))
    names=c(colnames(obj),name)
    objects[[i]]$x=data.frame(obj, dist)
    colnames(objects[[i]]$x)=names
  }
  return(objects)
}



direction2<-function(target, base=c(0,0,0)){
  R=target-base
  h=atan2(R$y,R$z)
  v=atan2(sqrt(R$y^2+R$z^2),R$x)
  angle=atan2(R$x,R$y)
  dist=sqrt(R$x^2+R$y^2+R$z^2)  
  distXY=sqrt(R$x^2+R$y^2)  
  
  x=sin(angle)
  y=cos(angle)
  E=data.frame(x,y)
  
  H=h*180/pi
  V=v*180/pi
  H=normalizeRotation(H)
  
  directionData<-list(data.frame(sec=target$sec, angle,norm=dist),data.frame(sec=target$sec,x=R$x,y=R$y,z=R$z,norm=dist))  
  
  names(directionData)=c("bearing","xyz")
  return(directionData)
}

angleForAlign<-function(vector){
  
  
  x=vector$x
  y=vector$y
  z=vector$z
  
  
  aZ=atan2(y,x)
  rotatedZ=data.frame(cos(aZ-aZ)*sqrt(x^2+y^2),sin(aZ-aZ)*sqrt(x^2+y^2),z)
  colnames(rotatedZ)=c("x","y","z")
  
  
  aY=atan2(rotatedZ$z,rotatedZ$x)
  rotatedY=data.frame(cos(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2),rotatedZ$y,sin(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2))
  colnames(rotatedY)=c("x","y","z")
  
  
  a<-data.frame(aZ,aY)
  colnames(a)=c("Z","Y")
  return(a)
}

rotationPolar<-function(vector,angle,xyz){
  
  x=vector$x
  y=vector$y
  z=vector$z
  
  aX=atan2(y,z)
  aY=atan2(z,x)
  aZ=atan2(y,x)
  
  if(xyz==1){
    rotated=data.frame(x,sin(aX-angle)*sqrt(z^2+y^2),cos(aX-angle)*sqrt(z^2+y^2))
  }
  if(xyz==2){
    rotated=data.frame(cos(aY-angle)*sqrt(z^2+x^2),y,sin(aY-angle)*sqrt(z^2+x^2))
  }
  if(xyz==3){
    rotated=data.frame(cos(aZ-angle)*sqrt(x^2+y^2),sin(aZ-angle)*sqrt(x^2+y^2),z)
  }
  
  colnames(rotated)=c("x","y","z")
  
  return(rotated)
}

normalizeRotation<-function(x){
  while(length(x[which(x>=360)])>0){
    x[which(x>=360)]=x[which(x>=360)]-360
  }
  
  while(length(x[which(x<0)])>0){
    
    x[which(x<0)]=x[which(x<0)]+360
  }
  
  return(x)
}

getLocalcoordinate<-function(objects){
  obj_num=length(objects)
  cp=objects[[obj_num]]
  objects=getRvector(objects, cp$x)
  ang=angleForAlign(cp$v)
  res0=list()
  for(i in 1:obj_num){
    obj=objects[[i]][2:length(objects[[i]])]
    
    res2=list()
    for(j in 1:length(obj)){
      rotated=cbind(sec=objects[[i]]$x$sec,rotationPolar(obj[[j]],ang$Z-pi/2,3))
      polar=direction2(rotated)
      name=c(names(objects[[i]]),paste0("rot_",names(obj)[j]))
      objects[[i]]=c(objects[[i]],add=list(polar))
      names(objects[[i]])=name
    }
    res1=c(x=list(objects[[i]]$rot_r$xyz),
           v=list(objects[[i]]$rot_v$xyz),
           a=list(objects[[i]]$rot_a$xyz))
    res0=c(res0,list(res1))
    
  }
  return(list(res0,objects))
}

getNorm<-function(df){
  res=apply(df,1,function(x){return(sqrt(x[1]^2+x[2]^2+x[3]^2))})
  return(res)
}

getFormationSD<-function(objects){
  
  sd_df=data.frame(sec=objects[[1]]$x$sec)
  
  for(j in 1:3){
    x=data.frame(sec=objects[[1]]$x$sec)
    for(i in 1:length(objects)){
      x=cbind(x,objects[[i]]$x[,1+j])
    }
    x=x[,-1]
    sd_x=apply(x,1,sd)
    sd_df=cbind(sd_df,sd_x)
  }
  norm=getNorm(sd_df[,2:4])
  sd_df=cbind(sd_df,norm)
  
  colnames(sd_df)=c("sec","x","y","z","norm")  
  
  
  obj_num=length(objects)
  objects[[obj_num]]=c(objects[[obj_num]],sd=list(sd_df))
  return(objects)
}


crop<-function(obj,t){
  if(is.list(obj)&!is.data.frame(obj)){
    for(i in 1:length(obj)){
      obj[[i]]=crop(obj[[i]],t)
    }
  }else{
    tgtindex=which(obj$sec>=t[1]&obj$sec<=t[2])
    obj=obj[tgtindex,]
  }
  return(obj)
}

swap<-function(obj){
  if(is.list(obj)&!is.data.frame(obj)){
    for(i in 1:length(obj)){
      obj[[i]]=swap(obj[[i]])
    }
  }else{
    if(ncol(obj)>=4){
      name=colnames(obj)
      index=1:length(name)
      index[2:4]=c(3,4,2)
      obj=obj[,index]
      colnames(obj)=name
    }
  }
  return(obj)
}


labeling<-function(obj,lab){
  if(is.list(obj)&!is.data.frame(obj)){
    for(i in 1:length(obj)){
      obj[[i]]=labeling(obj[[i]],lab)
    }
  }else{
    for(j in 1:length(lab)){
      label=rep(lab[j],nrow(obj))
      obj=cbind(obj,label)
    }
  }
  return(obj)
}

find_longestDuration<-function(index){
  if(length(index)<=1){return (c(1,1))}
  j=1
  after<-numeric()
  blanklength<-numeric()
  for(i in 1:(length(index)-1)){
    if(index[i+1]-index[i]>1){
      after[j]=i
      j=j+1
    }
  }
  
  blank=c(0,after,length(index))
  for(i in 1:(length(blank)-1)){
    blanklength[i]=blank[i+1]-blank[i]
  }
  start=(blank+1)[-length(blank)]
  end=start+blanklength-1
  
  tgtindex=which(blanklength==max(blanklength))
  return(c(index[start[tgtindex]],index[end[tgtindex]]))
}

check_croprange<-function(candidates,filepath){
  p=candidates[[1]]$ori_geo[[5]]
  
  v=candidates[[1]]$ori[[5]]$v$norm
  time=candidates[[1]]$ori[[5]]$v$sec
  
  path_list=list()
  range_list=list()
  
  num=length(candidates)
  for(c in 2:num){
    focused=candidates[[c]]
    cp=focused$ori_geo[[5]]
    vf=focused$ori[[5]]$v$norm
    
    range<-annotate(
      "rect",xmin=min(cp$x$sec),xmax=max(cp$x$sec), 
      ymin=min(vf), ymax=max(vf),
      color=hsv(1/num*(c-1),1,1),alpha=0.2, fill=hsv(1/num*(c-1),1,1))
    
    geom_path(
      aes(x = Longitude, y = Latitude), color =hsv(1/num*(c-1),1,1,alpha=1), linewidth = 1, data = cp$x
    )->path
    
    range_list=c(range_list,list(range))
    path_list=c(path_list,list(path))
  }
  
  v_line<-geom_line(aes(x=time,y = v), linewidth = 0.5)
  entire_path<-geom_path(
    aes(x = Longitude, y = Latitude), color =hsv(1,0,1,alpha=1), linewidth = 1, data = p$x
  )
  
  
  rangeeffect<-ggplot()+theme_bw()+
    theme(panel.grid = element_blank())+
    range_list+v_line+
    labs(
      x = "sec",
      y = "m/s"
    )
  
  print(rangeeffect)
  # print(map+entire_path+path_list)
  
  
  ggsave(
    bg="white",
    paste0(filepath,"/time_crop_area.png"),
    plot = rangeeffect,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  ggsave(paste0(filepath,"/map_crop_area.png"), map+entire_path+path_list,
         bg = "transparent",dpi=400,width=5,height=5)
  
  
}

crop_longestDuration<-function(object, vlim){
  cp=object$ori[[length(objects)]]
  time=cp$x$sec
  v=cp$v$norm
  dataset=data.frame(time,v)
  
  indexrange=which(v>vlim)
  range=find_longestDuration(indexrange)
  start=time[range[1]]
  end=time[range[2]]
  object_crop=crop(object,c(start,end))  
  
  dir.create("visualized data", showWarnings = FALSE)
  check_croprange(list(object,object_crop), "./visualized data")
  
  return(object_crop)
}


max.location_KDE<-function(x,y,n=100,h=NULL,lims=NULL){
  # h/lims default to the legacy per-series behaviour; pass the shared
  # bandwidth (and matching lims) to make the mode consistent with the
  # shared-scale isopleths drawn in draw2d()
  if(is.null(h))h=c(bandwidth.nrd(x),bandwidth.nrd(y))
  if(is.null(lims))lims=c(range(x),range(y))
  rd=kde2d(x=x,y=y,h=h,n=n,lims=lims)
  m=max(rd$z)
  rowcol_ <- which(rd$z==m)
  row_=(rowcol_ - 1) %/% nrow(rd$z) + 1
  col_=(rowcol_ - 1) %% nrow(rd$z) +1
  
  res=data.frame(x=rd$x[col_],y=rd$y[row_])
  return(res)
}



##################################################################################################################
####


calc_lambda_phi <- function(x, y, base_deg, precBits = 256) {
  lambda0_deg=base_deg[1]
  phi0_deg=base_deg[2]
  # 高精度数値型を設定
  x <- mpfr(x, precBits)
  y <- mpfr(y, precBits)
  lambda0_deg <- mpfr(lambda0_deg, precBits)
  phi0_deg <- mpfr(phi0_deg, precBits)
  
  # 定数
  m0 <- mpfr(0.9999, precBits)
  a <- mpfr(6378137.0, precBits)
  f <- mpfr(298.257222101, precBits)
  
  phi0_rad <- phi0_deg * mpfr(pi, precBits) / 180
  lambda0_rad <- lambda0_deg * mpfr(pi, precBits) / 180
  
  n <- mpfr(1.0, precBits) / (2 * f - 1)
  A_array <- A_array_mp(n, precBits)
  alpha_array <- alpha_array_mp(n, precBits)
  
  A_ <- ((m0 * a) / (1.0 + n)) * A_array[1]
  S_ <- ((m0 * a) / (1.0 + n)) * (A_array[1] * phi0_rad + 
                                    sum(A_array[2:6] * sin(2 * phi0_rad * mpfr(1:5, precBits))))
  
  # (1) xi', eta' の計算
  xi2 <- (y + S_) / A_
  eta2 <- x / A_
  
  # (2) xi, eta の計算（補正項の適用 - ループ方式で明示的に計算）
  xi <- xi2
  eta <- eta2
  for (i in 1:10) {  # 精度向上のために反復処理
    xi_correction <- mpfr(0, precBits)
    eta_correction <- mpfr(0, precBits)
    
    for (k in 1:5) {
      k_mpfr <- mpfr(k, precBits)
      xi_correction <- xi_correction + alpha_array[k + 1] * sin(2 * k_mpfr * xi) * cosh(2 * k_mpfr * eta)
      eta_correction <- eta_correction + alpha_array[k + 1] * cos(2 * k_mpfr * xi) * sinh(2 * k_mpfr * eta)
    }
    
    xi <- xi2 - xi_correction
    eta <- eta2 - eta_correction
  }
  
  # (3) t の計算
  t <- sinh(atanh(sin(xi)) + ((2 * sqrt(n)) / (1 + n)) * atanh(((2 * sqrt(n)) / (1 + n)) * sin(xi)))
  
  # (4) phi_rad の計算
  phi_rad <- asin(t / sqrt(1 + t^2))
  
  # (5) lambda_rad の計算
  lambda_rad <- lambda0_rad + atan2(sinh(eta), cos(xi))
  
  # 結果を度に変換
  phi_deg <- phi_rad * 180 / mpfr(pi, precBits)
  lambda_deg <- lambda_rad * 180 / mpfr(pi, precBits)
  
  # 結果を数値に戻して返す
  return(data.frame(lambda_deg = as.numeric(lambda_deg), phi_deg = as.numeric(phi_deg)))
}

calc_lambda<-function(x, y, base_deg, precBits = 256){
  result <- calc_lambda_phi(x, y, base_deg, precBits)
  offset=base_deg-calc_lambda_phi(0,0, base_deg, precBits) #計算過程でどうしても誤差が生じてしまうので、誤差値を使って補正
  return(as.data.frame(mapply(`+`, result,offset)))
}
# 高精度版 A_array
A_array_mp <- function(n, precBits) {
  n <- mpfr(n, precBits)
  c(
    1 + (n^2) / 4 + (n^4) / 64,
    -(3 / 2) * (n - (n^3) / 8 - (n^5) / 64),
    (15 / 16) * (n^2 - (n^4) / 4),
    -(35 / 48) * (n^3 - (5 / 16) * (n^5)),
    (315 / 512) * (n^4),
    -(693 / 1280) * (n^5)
  )
}

# 修正版 alpha_array
alpha_array_mp <- function(n, precBits) {
  n <- mpfr(n, precBits)
  c(
    mpfr(0, precBits),  # dummy
    (1.0 / 2) * n - (2.0 / 3) * (n^2) + (5.0 / 16) * (n^3) + (41.0 / 180) * (n^4) - (127.0 / 288) * (n^5),
    (13.0 / 48) * (n^2) - (3.0 / 5) * (n^3) + (557.0 / 1440) * (n^4) + (281.0 / 630) * (n^5),
    (61.0 / 240) * (n^3) - (103.0 / 140) * (n^4) + (15061.0 / 26880) * (n^5),
    (49561.0 / 161280) * (n^4) - (179.0 / 168) * (n^5),
    (34729.0 / 80640) * (n^5)
  )
}

scalebar_gmap<-function(map,locx,locy,norm,text=""){
  gg_build <- ggplot_build(map)
  panel_limits <- gg_build$layout$panel_params[[1]]$x.range # x軸範囲
  x_range <- panel_limits
  y_range <- gg_build$layout$panel_params[[1]]$y.range # y軸範囲
  
  origin=c(x_range[1],y_range[1])
  x1y0=c(x_range[2],y_range[1])
  x0y1=c(x_range[1],y_range[2])
  aspect=as.data.frame(rbind(x1y0,x0y1))
  diff_lambda=c((x1y0-origin)[1],(x0y1-origin)[2])
  diff_xy=calc_xy(aspect[,1],aspect[,2],origin[1],origin[2])
  angle_h=atan2(diff_xy$y[1],diff_xy$x[1])
  angle_v=atan2(diff_xy$x[2],diff_xy$y[2])
  r=seq(0,norm,length.out=20)
  x <- r*cos(angle_h)
  y <- r*sin(angle_h)
  
  loc<-function(v1,v2){
    return(c(origin[1]+diff_lambda[1]*v1, origin[2]+diff_lambda[2]*v2))
  }
  base=loc(locx,locy)
  bar <- calc_lambda(x, y, base)
  
  midpoint=c(
    (bar[nrow(bar),1]+bar[1,1])/2,
    (bar[nrow(bar),2]+bar[1,2])/2
    )
  
  textloc=midpoint+c(0,diff_lambda[2]*0.02)
  scale_bar=list(
    geom_path(aes(x = lambda_deg, y =phi_deg ),data=bar,col="white"),
    annotate("text",x = textloc[1],y = textloc[2],
      label = text, family = "Helvetica", color = "white"
    )
  )
  
  
  
  return(scale_bar)
}


