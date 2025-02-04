# NOTE: I want to expand this script to combine all samples.
# Should also test this over multiple sets of colocalization hyperparameters...

devtools::load_all(path = "/Users/jacobchang/Lab/metadisco")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(beepr)

# Starting analysis with assembloids_02 since it has the most complete time course
base_path_osi <- "data/coculture/osimertinib/assembloids_02"
base_path_control <- "data/coculture/control/assembloids_02"

time_points <- c("0h", "4h", "12h", "24h", "72h", "4dw", "7dw")
time_points_hours <- c(0, 4, 12, 24, 72, 168, 240)

osi_clq_list <- list()
control_clq_list <- list()
for(i in seq_along(time_points)){
  print(time_points[i])
  df <- read.csv(paste0(base_path_osi, "/Assembloid2_osi_", time_points[i], ".csv"))
  spomic <- createSpomic(p = df)
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
  
  permuted_clqs <- list()
  pb <- txtProgressBar(min = 1, max = 100, style = 3)
  for(j in 1:100) {
      df2 <- df
      df2$cell_type <- sample(df$cell_type)
      permuted_spomic <- createSpomic(p = df2)
      permuted_spomic <- setSpomicHypers(
        spomic_obj = permuted_spomic, 
        tile_size = NULL, 
        window_size = NULL, 
        k_neigh = 15, 
        bandwidth = 100, 
        n_bootstrap = NULL, 
        weight_scheme = "linear", 
        precompute_neighbors = TRUE
      )
      permuted_clqs[[j]] <- getCLQs(permuted_spomic)
      setTxtProgressBar(pb, j)
    }
    close(pb)
    permutatation_df <- bind_rows(permuted_clqs) |> 
      group_by(A_B) |> 
      summarize(lower_bound = quantile(colocalization_stat, probs = c(0.025), na.rm = TRUE),
                upper_bound = quantile(colocalization_stat, probs = c(0.975), na.rm = TRUE)) |> 
      ungroup() |> 
      select(A_B, lower_bound, upper_bound)    
    osi_clq_list[[i]] <- getCLQs(spomic) |> 
      dplyr::mutate(t = time_points_hours[i]) |> 
      left_join(permutatation_df)
    
}
beep()

for(i in seq_along(time_points)){
  print(time_points[i])
  df <- read.csv(paste0(base_path_control, "/Assembloid2_ctl_", time_points[i], ".csv"))
  spomic <- createSpomic(p = df)
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
  
  permuted_clqs <- list()
  pb <- txtProgressBar(min = 1, max = 100, style = 3)
  for(j in 1:100) {
    df2 <- df
    df2$cell_type <- sample(df$cell_type)
    permuted_spomic <- createSpomic(p = df2)
    permuted_spomic <- setSpomicHypers(
      spomic_obj = permuted_spomic, 
      tile_size = NULL, 
      window_size = NULL, 
      k_neigh = 15, 
      bandwidth = 100, 
      n_bootstrap = NULL, 
      weight_scheme = "linear", 
      precompute_neighbors = TRUE
    )
    permuted_clqs[[j]] <- getCLQs(permuted_spomic)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  permutatation_df <- bind_rows(permuted_clqs) |> 
    group_by(A_B) |> 
    summarize(lower_bound = quantile(colocalization_stat, probs = c(0.025), na.rm = TRUE),
              upper_bound = quantile(colocalization_stat, probs = c(0.975), na.rm = TRUE)) |> 
    ungroup() |> 
    select(A_B, lower_bound, upper_bound)    
  control_clq_list[[i]] <- getCLQs(spomic) |> 
    dplyr::mutate(t = time_points_hours[i]) |> 
    left_join(permutatation_df)
  
}
beep()

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
  p <- ggplot(total_df|> filter(A_B == a_b[i]), 
              aes(x = t, y = colocalization_stat, color = treatment, group = treatment)) +
    geom_point(aes(color = treatment), size = 1) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = treatment), color = NA, alpha = 0.1) + # Shaded region for bounds
    geom_vline(xintercept =  72, linetype = "dotted", color = "darkgray") + 
    labs(title = a_b[i],
         x = "Time",
         y = "Colocalization Statistic",
         color = "Treatment Group", 
         fill = "Random Colocalization Band") +
    theme_pubr() + 
    theme_publication_text()
print(p)
}
dev.off()

write.csv(total_df, "output/colocalization_vs_time/assembloids_02/colocalization_vs_time.csv", row.names = FALSE)

# ==============================================================================
