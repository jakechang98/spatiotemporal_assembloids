devtools::load_all(path = "/Users/jacobchang/Lab/metadisco")
source("utils/plotting_utils.R")
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

base_path <- "data/assembloids/assembloids_02"
time_points <- c("0h", "4h", "12h", "24h", "72h", "4dw", "7dw")
time_points_numeric <- c(0, 4, 12, 24, 72, 96, 168)

process_sample <- function(fp){
  sample = tools::file_path_sans_ext(basename(fp))
  df <- readxl::read_xlsx(path = fp) |> 
    dplyr::mutate(sample = sample) |> 
    dplyr::select(sample, X, Y, "Final cell type") |> 
    dplyr::rename(x = X, y = Y, cell_type = "Final cell type")
  
  return(df)
}

clq_list <- list()
for(i in seq_along(time_points)){
  df <- process_sample(fp = paste0(base_path, "/Assembloid2_osi_", time_points[i], ".xlsx"))
  spomic <- createSpomic(df)
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
  clq_list[[i]] <- getCLQs(spomic) |> dplyr::mutate(t = time_points_numeric[i], )
  print(plotSpomic(spomic, point_size = 0.1))
  print(plotCellProportions(spomic, show_pct = TRUE))
}
clqs <- dplyr::bind_rows(clq_list)

# Generate and save plots
pdf(file = file.path("output", 
                     "colocalization_vs_time", 
                     "assembloids_02", 
                     "colocalization_plots.pdf"), 
    width = 6.5, 
    height = 3)

a_b <- unique(clqs$A_B)
for(i in seq_along(a_b)) {
  df <- clqs |> filter(A_B == a_b[i])
  p <- ggplot(df, aes(x = t, y = colocalization_stat, group = 1)) + 
    geom_point() + 
    geom_line() + 
    geom_vline(xintercept =  72, linetype = "dotted", color = "darkgray") + 
    labs(title = a_b[i], x = "Time (hours)", y = "CLQ") + 
    theme_pubr() + 
    theme_publication_text()
  print(p)
}
dev.off()













