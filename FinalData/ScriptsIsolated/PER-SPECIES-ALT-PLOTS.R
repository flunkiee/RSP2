library(tidyverse)
library(patchwork)
library(tcltk)    # for tk_choose.files()
library(tools)    # for file_path_sans_ext()

# Step 1: Choose multiple files using tcltk
file_paths <- tk_choose.files(caption = "Select gene data files", multi = TRUE)

# Step 2: Function to read and label each file
read_gene_file <- function(file_path) {
  gene_name <- file_path_sans_ext(basename(file_path))  # extract gene name from filename
  df <- read_csv(file_path)
  df <- df %>%
    mutate(
      Gene = gene_name,
      Pair = paste(`Species 1`, `Species 2`, sep = " vs ")
    )
  return(df)
}

# Step 3: Read and combine all files
all_data <- map_dfr(file_paths, read_gene_file)

# Step 4: Plot function
plot_gene <- function(data, gene_name) {
  ggplot(data, aes(x = Pair, y = Dist)) +
    geom_point(shape = 21, color = "black", fill = "white") +
    geom_errorbar(aes(ymin = Dist - `Std. Err`, ymax = Dist + `Std. Err`), width = 0.2) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = gene_name,
      y = "n-s"
    )
}

# Step 5: Generate plots by gene
plots <- all_data %>%
  split(.$Gene) %>%
  map2(names(.), ~ plot_gene(.x, .y))

# Step 6: Combine and display plots
wrap_plots(plots) +
  plot_annotation(title = "Metschnikowia pulcherrima pN-pS")

