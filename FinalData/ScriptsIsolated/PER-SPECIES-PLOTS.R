library(tidyverse)
library(patchwork)

# Step 1: File chooser (opens a file selection dialog)
file_path <- file.choose()  # this opens your system file explorer

# Step 2: Load the selected CSV file
df <- read_csv(file_path)

# Step 3: Extract gene name and prepare plotting labels
df <- df %>%
  mutate(
    Gene = str_extract(`Species 1`, "(?<=:)\\w+$"),
    Pair = paste(`Species 1`, `Species 2`, sep = " vs ")
  )

plot_gene <- function(data, gene_name) {
  ggplot(data, aes(x = Pair, y = Dist)) +
    geom_point(shape = 21, color = "black", fill = "white") +  # White fill, black border
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




# Step 5: Generate and combine all plots by gene
plots <- df %>%
  split(.$Gene) %>%
  map2(names(.), ~ plot_gene(.x, .y))

# Step 6: Display all plots together
wrap_plots(plots) +
  plot_annotation(title = "Bacillus subtilis group dN-dS")
