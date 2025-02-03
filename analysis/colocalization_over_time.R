devtools::load_all(path = "/Users/jacobchang/Lab/metadisco")
library(dplyr)
library(ggplot2)
library(ggpubr)

# Starting analysis with assembloids_02 since it has the most complete time course
base_path_osi <- "data/coculture/osimertinib/assembloids_02"
base_path_control <- "data/coculture/control/assembloids_02"



time_points <- c("0h", "4h", "12h", "24h", "72h", "4dw", "7dw")
time_points_hours <- c(0, 4, 12, 24, 72, 168, 240)


osi_clq_list <- list()
control_clq_list <- list()
for(i in seq_along(time_points)){
  spomic <- createSpomic(p = paste0(base_path_osi, "/Assembloid2_osi_", time_points[i], ".csv"))
  spomic <- setSpomicHypers(
    spomic_obj = spomic, 
    tile_size = NULL, 
    window_size = NULL, 
    k_neigh = 15, 
    bandwidth = 100, 
    n_bootstrap = NULL, 
    weight_scheme = "linear", 
    precompute_neighbors = TRUE
  )
  osi_clq_list[[i]] <- getCLQs(spomic) |> dplyr::mutate(t = time_points_hours[i])
  # print(plotSpomic(spomic, point_size = 0.1))
  # print(plotCellProportions(spomic, show_pct = TRUE))
  
  spomic <- createSpomic(p = paste0(base_path_control, "/Assembloid2_ctl_", time_points[i], ".csv"))
  spomic <- setSpomicHypers(
    spomic_obj = spomic, 
    tile_size = 250, 
    window_size = 10, 
    k_neigh = 15, 
    bandwidth = 100, 
    n_bootstrap = 100, 
    weight_scheme = "linear", 
    precompute_neighbors = TRUE
  )
  control_clq_list[[i]] <- getCLQs(spomic) |> dplyr::mutate(t = time_points_hours[i])
  # print(plotSpomic(spomic, point_size = 0.1))
  # print(plotCellProportions(spomic, show_pct = TRUE))
}
osi_clqs <- dplyr::bind_rows(osi_clq_list)
control_clqs <- dplyr::bind_rows(control_clq_list)

# Generate and save plots
a_b <- unique(c(osi_clqs$A_B, control_clqs$A_B))
osi_df <- osi_clqs |> dplyr::mutate(treatment = "osi")
control_df <- control_clqs |> dplyr::mutate(treatment = "control")
total_df <- rbind(osi_df, control_df)

pdf(file = file.path("output", 
                     "colocalization_vs_time", 
                     "assembloids_02", 
                     "colocalization_plots.pdf"), 
    width = 6.5, 
    height = 3)

for(i in seq_along(a_b)) {
  df <- total_df |> filter(A_B == a_b[i])
  p <- ggplot(df, aes(x = t, y = colocalization_stat, color = treatment, group = treatment)) + 
    geom_point() + 
    geom_line() + 
    geom_vline(xintercept =  72, linetype = "dotted", color = "darkgray") + 
    labs(title = a_b[i], x = "Time (hours)", y = "CLQ", color = "Group") + 
    theme_pubr() + 
    theme_publication_text()
  print(p)
}
dev.off()
