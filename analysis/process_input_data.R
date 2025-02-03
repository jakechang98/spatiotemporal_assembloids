devtools::load_all(path = "/Users/jacobchang/Lab/metadisco")
library(dplyr)

# Restructure data frames for use with Spomic
process_sample <- function(fp){
  sample <- tools::file_path_sans_ext(basename(fp))
  df <- readxl::read_xlsx(path = fp) |> 
    dplyr::mutate(sample = sample) |> 
    dplyr::select(sample, X, Y, "Final cell type", "Cell type number") |> 
    dplyr::rename(x = X, y = Y, cell_type = "Final cell type", cell_type_cluster_id = "Cell type number")
  
  return(df)
}


# Process all samples and rewrite to new directory
main <- function(base_directory, output_directory) {
  if(!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  files <-  list.files(base_directory, full.names = TRUE)
  files_short_name <- list.files(base_directory, full.names = FALSE)
  pb <- txtProgressBar(min = 1, max = length(files), style = 3)
  for(i in seq_along(files)) {
    df <- process_sample(files[i])
    write.csv(df, 
              file = file.path(
                output_directory, 
                paste0(tools::file_path_sans_ext(files_short_name[i]), ".csv")),
              row.names = FALSE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
}


conditions <- c("coculture")
assembloids <- c("assembloids_01", "assembloids_02", "assembloids_03")

for(c in conditions) {
  for(a in assembloids) {
    main(base_directory = paste0("data/original_celesta_data/", c,"/osimertinib/", a), 
         output_directory = paste0("data/", c, "/osimertinib/", a))
    main(base_directory = paste0("data/original_celesta_data/", c,"/control/", a), 
         output_directory = paste0("data/", c, "/control/", a))
  }
}

# ==============================================================================
