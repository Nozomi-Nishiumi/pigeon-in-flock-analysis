#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2025/10/8.
if( dev.cur() > 1 ) dev.off()


####### For Flight velocity analysis
drone_engage=c(-83.3928,33.8147)

options(digits=22)
basedir="/Users/nozomi/Dropbox/Nozomi data share/Homing data Analysis"
scandir=paste0(basedir,"/csv")
Day=c("#1_Pre_drone","#3_Post_drone")
FL=c("FL01","FL02","FL03")



list_perDay=list()
for(i in Day){
  list_perFL=list()
  for(j in FL){
    dir=paste0(scandir,"/",i,"/",j)
    setwd(dir)
    load("./transformed data/final product/parameters.Rdata")
    
    obj_num=length(object$ori_geo)
    if(length(object$rot)!=obj_num){
      print("num of objects is not consistent")
      exit()
    }
    
    cp_x=object$ori[[length(object$ori)]]$x
    cp_v=object$ori[[length(object$ori)]]$v
    cp_geo=object$ori_geo[[length(object$ori_geo)]]$x
    
    dist_drone=horizontal_distance(cp_geo,drone_engage)
    cp_v=cbind(cp_v,dist_drone)
    closest_index=which(cp_v$dist_drone==min(cp_v$dist_drone))
    app_cp_v=cp_v[1:closest_index,]#ドローン最近接地点までのデータを対象とする
  
    plot(app_cp_v$sec,app_cp_v$dist_drone)
        
    list_perFL=c(list_perFL,setNames(list(app_cp_v),j))
  }
  list_perDay=c(list_perDay,setNames(list(list_perFL),i))
}

df_plot <- bind_rows(
  list_perDay$`#1_Pre_drone`$FL01%>% mutate(group = "Day1"),
  list_perDay$`#3_Post_drone`$FL01%>% mutate(group = "Day2")
)

# 折れ線グラフを描画
ggplot(df_plot, aes(x = -1*dist_drone, y = mag*3600/1000, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "FL01",
    x = "distance to drone area",
    y = "flight velocity (km/h)",
    color = "Dataset"
  )



setwd(basedir)
if(!dir.exists("fig")){
  dir.create("fig")
}

# ループ処理
for (fl in FL) {
  
  # データをまとめる
  df_plot <- bind_rows(
    list_perDay$`#1_Pre_drone`[[fl]] %>% mutate(group = "Day1"),
    list_perDay$`#3_Post_drone`[[fl]] %>% mutate(group = "Day3")
  )
  
  # 折れ線グラフ_横軸はsec
  p <- ggplot(df_plot, aes(x = sec, y = mag*3600/1000, color = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    theme_minimal(base_size = 14) +
    labs(
      title = fl,
      x = "sec",
      y = "flight velocity (km/h)",
      color = ""
    )
  
  print(p)  # RStudioなどでループ中にも確認できるように出力
  
  
  # PNGで保存
  ggsave(
    # bg="white",
    filename = paste0("fig/velocity_",fl, "_sec.png"),
    plot = p,
    width = 7,
    height = 4,
    dpi = 300
  )


# 折れ線グラフ_横軸はsec
p <- ggplot(df_plot, aes(x = -dist_drone, y = mag*3600/1000, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = fl,
    x = "dist",
    y = "flight velocity (km/h)",
    color = ""
  )

print(p)  # RStudioなどでループ中にも確認できるように出力


# PNGで保存
ggsave(
  # bg="white",
  filename = paste0("fig/velocity_",fl, "_dist.png"),
  plot = p,
  width = 7,
  height = 4,
  dpi = 300
)
}




for (fl in FL) {
  
  # 1. データ結合
  df_box <- bind_rows(
    list_perDay$`#1_Pre_drone`[[fl]] %>% mutate(group = "Day1"),
    list_perDay$`#3_Post_drone`[[fl]] %>% mutate(group = "Day3")
  ) 
  # 2. 箱ひげ図
  p <- ggplot(df_box, aes(x = group, y = mag*3600/1000, fill = group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "gray60") +
    # geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +  # 個別データも薄く表示
    theme_minimal(base_size = 14) +
    theme(
      # panel.background = element_rect(fill = "white", color = NA),
      # plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    labs(
      title = fl,
      x = "",
      y = "Flight velocity (km/h)"
    )
  
  print(p)  # RStudioで確認
  
  # 3. 画像保存
  ggsave(
    filename = paste0("fig/velocity_",fl, "_boxplot.png"),
    plot = p,
    width = 2,
    height = 4,
    dpi = 300
  )
}
