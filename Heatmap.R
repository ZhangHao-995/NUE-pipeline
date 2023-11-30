library(ComplexHeatmap)
library(BuenColors)

h1 <- Heatmap(plot_ac[,6:9],name="H3K27ac",cluster_columns = F,
  row_split = plot_ac$cluster,col=jdb_palette("solar_extra", type = "continuous"),
  border=T,column_title = "Promoter H3K27ac",
  heatmap_legend_param = list(direction = "horizontal"),column_names_rot = 45)
ha = rowAnnotation(Marker = anno_mark(at = 1:11, labels = plot_exp$name))
h2 <- Heatmap(plot_exp[,6:9],name="Exp",cluster_columns = F,
  row_split = plot_exp$cluster,row_title = NULL,right_annotation = ha,
  show_row_dend = F,border=T,column_title = "DEG Expression",
  col=jdb_palette("wolfgang_extra", type = "continuous"),
  heatmap_legend_param = list(direction = "horizontal"),column_names_rot = 45)
draw(h1+h2,main_heatmap="Exp",heatmap_legend_side = "bottom")
