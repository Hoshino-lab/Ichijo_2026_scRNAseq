

library(tidyr)

# Output directory (relative path; created if it does not exist)
out_dir <- "results"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# 1. Extract metadata and assign cell types and mouse IDs

metadata <- ASlin0.5_sc5_0.1[[]]
cluster_map <- c("0" = "Ast", "1" = "BG", "2" = "BGLP")

library_key <- expand_grid(
  stage = c("P3", "P10"),
  genotype = c("ctr", "cKO"),
  suffix = c("n1b1", "n1l1", "n1l2")
) %>%
  mutate(
    library = paste(stage, genotype, suffix, sep = "_"),
    bio_rep = if_else(suffix == "n1b1", "rep1", "rep2")
  )

library_index <- match(as.character(metadata$library), library_key$library)
metadata_rep <- sub("^.*_(rep[12])$", "\\1", as.character(metadata$bio_rep))

# Verify the cluster labels and the correspondence between libraries and mice.
stopifnot(
  !anyNA(metadata$seurat_clusters),
  setequal(as.character(metadata$seurat_clusters), names(cluster_map)),
  !anyNA(library_index),
  all(as.character(metadata$stage) == library_key$stage[library_index]),
  all(as.character(metadata$genotype) == library_key$genotype[library_index]),
  all(metadata_rep == library_key$bio_rep[library_index])
)

cell_metadata <- tibble(
  stage = as.character(metadata$stage),
  genotype = as.character(metadata$genotype),
  bio_rep = library_key$bio_rep[library_index],
  cell_type = unname(cluster_map[as.character(metadata$seurat_clusters)])
)

# 2. Pool split libraries and calculate proportions for each mouse

cell_composition <- cell_metadata %>%
  count(stage, genotype, bio_rep, cell_type, name = "n_cells") %>%
  complete(
    stage = c("P3", "P10"),
    genotype = c("ctr", "cKO"),
    bio_rep = c("rep1", "rep2"),
    cell_type = c("Ast", "BG", "BGLP"),
    fill = list(n_cells = 0)
  ) %>%
  group_by(stage, genotype, bio_rep) %>%
  mutate(
    total_lineage_cells = sum(n_cells),
    proportion = n_cells / total_lineage_cells,
    percentage = 100 * proportion
  ) %>%
  ungroup() %>%
  mutate(
    stage = factor(stage, levels = c("P3", "P10")),
    genotype = factor(genotype, levels = c("ctr", "cKO")),
    cell_type = factor(cell_type, levels = c("Ast", "BG", "BGLP")),
    mouse_id = paste(stage, genotype, bio_rep, sep = "__")
  )

# An absent cell type is represented by zero; an absent mouse is not allowed.
stopifnot(all(cell_composition$total_lineage_cells > 0))

# 3. Calculate mean proportions and SEM across mice
# Each mouse receives equal weight, irrespective of its cell count.

composition_summary <- cell_composition %>%
  group_by(stage, genotype, cell_type) %>%
  summarise(
    n_mice = n_distinct(bio_rep),
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion),
    .groups = "drop"
  ) %>%
  mutate(
    sem_proportion = sd_proportion / sqrt(n_mice),
    mean_percentage = 100 * mean_proportion,
    sd_percentage = 100 * sd_proportion,
    sem_percentage = 100 * sem_proportion
  )

stopifnot(all(composition_summary$n_mice == 2L))

# 4. Calculate SEM for cumulative boundaries of the stacked bars
# BGLP boundary = BGLP proportion.
# BGLP + BG boundary = cumulative BGLP + BG proportion, not BG alone.
# Calculate cumulative values within each mouse before calculating SEM.

boundary_by_mouse <- cell_composition %>%
  select(stage, genotype, bio_rep, cell_type, proportion) %>%
  pivot_wider(names_from = cell_type, values_from = proportion) %>%
  transmute(
    stage, genotype, bio_rep,
    BGLP_boundary = BGLP,
    BGLP_BG_boundary = BGLP + BG
  ) %>%
  pivot_longer(
    cols = c(BGLP_boundary, BGLP_BG_boundary),
    names_to = "boundary",
    values_to = "boundary_proportion"
  )

boundary_summary <- boundary_by_mouse %>%
  group_by(stage, genotype, boundary) %>%
  summarise(
    n_mice = n_distinct(bio_rep),
    mean_boundary = mean(boundary_proportion),
    sd_boundary = sd(boundary_proportion),
    .groups = "drop"
  ) %>%
  mutate(
    sem_boundary = sd_boundary / sqrt(n_mice),
    upper_sem = mean_boundary + sem_boundary,
    x_position = as.integer(genotype)
  )

# 5. Shared plot settings

cell_colors <- c(Ast = "#A6CEE3", BG = "#1F78B4", BGLP = "#B2DF8A")

fill_scale <- function() {
  scale_fill_manual(
    values = cell_colors,
    breaks = c("Ast", "BG", "BGLP"),
    labels = c("Astrocyte", "BG", "BGLP"),
    drop = FALSE
  )
}

y_scale <- function() {
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  )
}

shared_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "top",
    plot.margin = margin(10, 15, 10, 10)
  )

# 6. Individual-mouse stacked bars (eight mice; no error bars)
# Default stacking places BGLP at the bottom, followed by BG and Ast.

mouse_levels <- c(
  "P3__ctr__rep1", "P3__ctr__rep2",
  "P3__cKO__rep1", "P3__cKO__rep2",
  "P10__ctr__rep1", "P10__ctr__rep2",
  "P10__cKO__rep1", "P10__cKO__rep2"
)

mouse_labels <- setNames(
  c(
    "P3 Ctrl\nrep1", "P3 Ctrl\nrep2",
    "P3 cKO\nrep1", "P3 cKO\nrep2",
    "P10 Ctrl\nrep1", "P10 Ctrl\nrep2",
    "P10 cKO\nrep1", "P10 cKO\nrep2"
  ),
  mouse_levels
)

cell_composition <- cell_composition %>%
  mutate(mouse_id = factor(mouse_id, levels = mouse_levels)) %>%
  arrange(mouse_id, cell_type)

p_by_mouse <- ggplot(
  cell_composition,
  aes(x = mouse_id, y = proportion, fill = cell_type)
) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.45) +
  facet_grid(. ~ stage, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = mouse_labels, drop = TRUE) +
  fill_scale() +
  y_scale() +
  coord_cartesian(ylim = c(0, 1.03)) +
  labs(x = NULL, y = "Astroglial lineage cell proportion", fill = NULL) +
  shared_theme

# 7. Mean stacked bars with upper SEM at the two cumulative boundaries
# Explicit rectangles set the order from bottom to top: BGLP, BG, Ast.

mean_rects <- composition_summary %>%
  mutate(
    stack_order = match(as.character(cell_type), c("BGLP", "BG", "Ast")),
    x_position = as.integer(genotype)
  ) %>%
  arrange(stage, genotype, stack_order) %>%
  group_by(stage, genotype) %>%
  mutate(
    y_max = cumsum(mean_proportion),
    y_min = lag(y_max, default = 0)
  ) %>%
  ungroup()

p_mean_upper_sem <- ggplot() +
  geom_rect(
    data = mean_rects,
    aes(
      xmin = x_position - 0.36, xmax = x_position + 0.36,
      ymin = y_min, ymax = y_max, fill = cell_type
    ),
    colour = "black", linewidth = 0.45
  ) +
  geom_segment(
    data = boundary_summary,
    aes(x = x_position, xend = x_position, y = mean_boundary, yend = upper_sem),
    colour = "black", linewidth = 0.7
  ) +
  geom_segment(
    data = boundary_summary,
    aes(
      x = x_position - 0.06, xend = x_position + 0.06,
      y = upper_sem, yend = upper_sem
    ),
    colour = "black", linewidth = 0.7
  ) +
  facet_grid(. ~ stage) +
  scale_x_continuous(
    breaks = c(1, 2), labels = c("Control", "En1::Meis1 cKO"),
    limits = c(0.4, 2.6), expand = expansion(mult = c(0, 0))
  ) +
  fill_scale() +
  y_scale() +
  coord_cartesian(ylim = c(0, 1.03)) +
  labs(x = NULL, y = "Astroglial lineage cell proportion", fill = NULL) +
  shared_theme

#sup7.j
print(p_by_mouse)
ggsave(file.path(fig_dir, "01_composition_by_mouse.png"), p_by_mouse,
       width = 12, height = 6.5, units = "in", dpi = 300, bg = "white")

#fig2.d
print(p_mean_upper_sem)
ggsave(file.path(fig_dir, "02_composition_mean_upper_SEM.png"), p_mean_upper_sem,
       width = 10, height = 6.5, units = "in", dpi = 300, bg = "white")



