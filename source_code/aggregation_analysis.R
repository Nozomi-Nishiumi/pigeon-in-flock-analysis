#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2025/10/8.
if( dev.cur() > 1 ) dev.off()


####### For Flight velocity analysis
drone_engage=c(-83.3928,33.8147)

options(digits=22)
basedir=dirname(here::here())
scandir=paste0(basedir,"/csv")
Day=c("Day1_Pre_drone","Day3_Post_drone")
FL=c("FL01","FL02","FL03")



list_perDay=list()
for(i in Day){
  list_perFL=list()
  for(j in FL){
    dir=paste0(scandir,"/",i,"/",j)
    setwd(dir)
    parts <- strsplit(dir, "/")[[1]]
    day <- strsplit(parts[length(parts)-1], "_")[[1]][1]
    prop <- paste0(day,"_",tail(parts, 1))
    load(paste0("./flock_frame data/","flockframe_",prop,".Rdata"))
    obj_num=length(object_crop$ori_geo)
    if(length(object_crop$rot)!=obj_num){
      print("num of objects is not consistent")
      exit()
    }
  
      
    
    dev<-data.frame(matrix(rep(NA, length(2)), ncol=3))[,numeric(0)]
    for(individual in 1:(obj_num-1)){
      dev=cbind(dev, object_crop$rot[[individual]]$x$dist_center_of_three)
      colnames(dev)[individual]=paste0("P",individual)
    }
    ave <- apply(dev, 1,mean)
    s.dev <- apply(dev, 1,sd)
    
    cp_geo=object_crop$ori_geo[[length(object_crop$ori_geo)]]$x
    dist_drone=horizontal_distance(cp_geo,drone_engage)
    
    dev=cbind(sec=object_crop$rot[[1]]$x$sec,dev,mean=ave, sd=s.dev, dist_drone)
    plot(dev$sec,dev$dist_drone)
    plot(object_crop$ori[[5]]$v$sec,object_crop$ori[[5]]$v$norm)
    
    
    list_perFL=c(list_perFL,setNames(list(dev),j))
  }
  list_perDay=c(list_perDay,setNames(list(list_perFL),i))
}

savedir=paste0(basedir,"/derived_Rdata")
dir.create(savedir, showWarnings = FALSE)
save(list_perDay, file=paste0(savedir,"/list_perDay.Rdata"))

df_plot <- bind_rows(
  list_perDay$`Day1_Pre_drone`$FL01%>% mutate(group = "Day1"),
  list_perDay$`Day3_Post_drone`$FL01%>% mutate(group = "Day3")
)




savedir=paste0(basedir,"/figures/closeness")
dir.create(savedir, showWarnings = FALSE)

for (fl in FL) {
  
  df_plot <- bind_rows(
    list_perDay$`Day1_Pre_drone`[[fl]] %>% mutate(group = "Day1"),
    list_perDay$`Day3_Post_drone`[[fl]] %>% mutate(group = "Day3")
  )
  
  p <- ggplot(df_plot, aes(x = sec, y = mean, color = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    theme_minimal(base_size = 14) +
    labs(
      title = fl,
      x = "sec",
      y = "Mean distance to center (m)",
      color = ""
    )
  
  print(p)  
  
  
  ggsave(
    bg="white",
    filename = paste0(savedir,"/deviation_",fl, "_sec.png"),
    plot = p,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  p <- ggplot(df_plot, aes(x = -dist_drone, y = mean, color = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    theme_minimal(base_size = 14) +
    labs(
      title = fl,
      x = "Distance to drone area (m)",
      y = "Mean distance to center (m)",
      color = ""
    )
  
  print(p)  
  
  
  ggsave(
    bg="white",
    filename = paste0(savedir,"/deviation_",fl, "_dist.png"),
    plot = p,
    width = 7,
    height = 4,
    dpi = 300
  )
}


for (fl in FL) {
  
  # 1. data binding
  cols_to_compare <- c("P2","P3","P4") #set target columns for comparison. Here we choose flockmates other than drone-experienced one, P1. 
  
  df_box <- bind_rows(
    list_perDay$`Day1_Pre_drone`[[fl]] %>% mutate(group = "Day1")%>%
      pivot_longer(cols = all_of(cols_to_compare), names_to = "variable", values_to = "value"),
    list_perDay$`Day3_Post_drone`[[fl]] %>% mutate(group = "Day3")%>%
      pivot_longer(cols = all_of(cols_to_compare), names_to = "variable", values_to = "value")
  ) 
  
  
  # 2. box plot
  # plot_each <- ggplot(df_box, aes(x = variable, y = value, fill = group)) +
  #   geom_boxplot(width = 0.4, position = position_dodge(width = 0.4)) +
  #  
  #   theme_minimal(base_size = 14) +
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background = element_rect(fill = "white", color = NA),
  #     legend.position = "none"
  #   ) +
  #   labs(
  #     title = fl,
  #     x = "",
  #     y = "Distance to center (m)"
  #   )
  
  plot_each <- ggplot(
    df_box,
    aes(
      x = variable,
      y = value,
      group = interaction(variable, group)
    )
  ) +
    geom_boxplot(
      width = 0.4,
      position = position_dodge(width = 0.4),
      fill = "gray"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    labs(
      title = fl,
      x = "",
      y = "Distance to center (m)"
    )
  
  print(plot_each)  
  
  
  # plot_total=ggplot(df_box, aes(x = group, y = value, fill = group)) +
  #   geom_boxplot() +
  #   theme_minimal(base_size = 14) +
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background = element_rect(fill = "white", color = NA),
  #     legend.position = "none"
  #   ) +
  #   labs(
  #     title = fl,
  #     x = "",
  #     y = "Distance to center (m)"
  #   )
  
  plot_total <- ggplot(
    df_box,
    aes(x = group, y = value, group = group)
  ) +
    geom_boxplot(
      fill = "gray"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    labs(
      title = fl,
      x = "",
      y = "Distance to center (m)"
    )
  
  print(plot_total)  
  
  
  # plot_data <- bind_rows(df_box, df_box_total)
  # 
  # ggplot(plot_data, aes(x = group, y = value, fill = group)) +
  #   geom_boxplot(width = 0.45, position = position_dodge(width = 0.45)) +
  #   facet_wrap(~variable, scales = "free_x") +  # 列ごと + total を横並びに
  #   theme_minimal(base_size = 14) +
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background = element_rect(fill = "white", color = NA),
  #     legend.position = "none"
  #   ) +
  #   labs(
  #     title = fl,
  #     x = "",
  #     y = "Distance to center (m)"
  #   )
  # 
  
  
  # 3. save image
  ggsave(
    filename = paste0(savedir, "/deviation_",fl, "_boxplot_indivitual.png"),
    plot = plot_each,
    width = 8,
    height = 4,
    dpi = 300
  )


ggsave(
  filename = paste0(savedir, "/deviation_",fl, "_boxplot_total.png"),
  plot = plot_total,
  width = 3,
  height = 4,
  dpi = 300
)
}
# 
# 
# for (fl in FL) {
#   
#   # 1. data binding
#   df_box <- bind_rows(
#     list_perDay$`Day1_Pre_drone`[[fl]] %>% mutate(group = "Day1"),
#     list_perDay$`Day3_Post_drone`[[fl]] %>% mutate(group = "Day3")
#   ) 
#   # 2. box plot
#   p <- ggplot(df_box, aes(x = group, y = mean, fill = group)) +
#     geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "gray60") +
#     # geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +  # 個別データも薄く表示
#     theme_minimal(base_size = 14) +
#     theme(
#       panel.background = element_rect(fill = "white", color = NA),
#       plot.background = element_rect(fill = "white", color = NA),
#       legend.position = "none"
#     ) +
#     labs(
#       title = fl,
#       x = "",
#       y = "Mean distance to center (m)"
#     )
#   
#   print(p)  
#   
# 3. save image
#   ggsave(
#     filename = paste0("fig/deviation_",fl, "_boxplot.png"),
#     plot = p,
#     width = 2,
#     height = 4,
#     dpi = 300
#   )
# }
# 



