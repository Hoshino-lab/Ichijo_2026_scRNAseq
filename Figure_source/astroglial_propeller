
library(speckle)
library(limma)
library(dplyr)
library(tidyr)

# 1. Prepare mouse-level cell counts from the composition analysis

composition_input <- cell_composition


da_counts <- composition_input %>%
  transmute(
    stage = as.character(stage),
    genotype = as.character(genotype),
    bio_rep = as.character(bio_rep),
    pb_celltype = as.character(cell_type),
    n_cells = n_cells
  )

# Accept short replicate labels or matching full mouse IDs.
rep_short <- sub("^.*_(rep[12])$", "\\1", da_counts$bio_rep)
stopifnot(
  !anyNA(da_counts),
  is.numeric(da_counts$n_cells),
  all(is.finite(da_counts$n_cells)),
  all(da_counts$n_cells >= 0),
  all(da_counts$n_cells == round(da_counts$n_cells)),
  all(da_counts$stage %in% c("P3", "P10")),
  all(da_counts$genotype %in% c("ctr", "cKO")),
  all(da_counts$pb_celltype %in% c("Ast", "BG", "BGLP")),
  all(rep_short %in% c("rep1", "rep2")),
  all(da_counts$bio_rep == rep_short |
        da_counts$bio_rep == paste(
          da_counts$stage, da_counts$genotype, rep_short, sep = "_"
        ))
)
da_counts$bio_rep <- rep_short
stopifnot(!anyDuplicated(
  da_counts[c("stage", "genotype", "bio_rep", "pb_celltype")]
))

# 2. Complete zero-count cell types and calculate proportions

da_counts <- da_counts %>%
  complete(
    stage = c("P3", "P10"),
    genotype = c("ctr", "cKO"),
    bio_rep = c("rep1", "rep2"),
    pb_celltype = c("Ast", "BG", "BGLP"),
    fill = list(n_cells = 0)
  ) %>%
  mutate(mouse_id = paste(stage, genotype, bio_rep, sep = "_")) %>%
  group_by(stage, genotype, bio_rep, mouse_id) %>%
  mutate(
    total_lineage_cells = sum(n_cells),
    proportion = n_cells / total_lineage_cells,
    percentage = 100 * proportion
  ) %>%
  ungroup()

# Missing cell types are zero; a missing mouse must not enter the analysis.
stopifnot(all(da_counts$total_lineage_cells > 0))

# 3. Summarize observed proportions and effect sizes
# Each mouse receives equal weight, irrespective of its cell count.

da_summary <- da_counts %>%
  group_by(stage, genotype, pb_celltype) %>%
  summarise(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion),
    n_mice = n_distinct(mouse_id),
    .groups = "drop"
  ) %>%
  mutate(sem_proportion = sd_proportion / sqrt(n_mice))

da_effects <- da_summary %>%
  select(stage, pb_celltype, genotype, mean_proportion) %>%
  pivot_wider(
    names_from = genotype,
    values_from = mean_proportion,
    names_prefix = "mean_"
  ) %>%
  mutate(
    difference_percentage_points = 100 * (mean_cKO - mean_ctr)
  )

# 4. Test genotype differences within each age
# Use the same settings for both transformations; do not select a
# transformation based on statistical significance.

da_results <- list()

for (transform_i in c("logit", "asin")) {
  for (stage_i in c("P3", "P10")) {
    a <- as.data.frame(filter(da_counts, stage == stage_i))
    stopifnot(all(tapply(a$n_cells, a$pb_celltype, sum) > 0))

    # Reconstruct cell-type/sample labels exactly from observed counts.
    # This does not simulate cells or create biological replicates.
    idx <- rep(seq_len(nrow(a)), times = a$n_cells)
    props <- speckle::getTransformedProps(
      clusters = factor(a$pb_celltype[idx], levels = c("Ast", "BG", "BGLP")),
      sample = a$mouse_id[idx],
      transform = transform_i
    )

    # Align sample metadata to the columns of the proportion matrix.
    info <- unique(a[, c("mouse_id", "genotype", "bio_rep")])
    info <- info[match(colnames(props$Proportions), info$mouse_id), , drop = FALSE]
    info$genotype <- factor(info$genotype, levels = c("ctr", "cKO"))
    stopifnot(
      !anyNA(info$mouse_id),
      identical(info$mouse_id, colnames(props$Proportions)),
      identical(as.integer(table(info$genotype)), c(2L, 2L))
    )

    design <- model.matrix(~ 0 + genotype, data = info)
    rownames(design) <- info$mouse_id
    contrast <- limma::makeContrasts(
      genotypecKO - genotypectr, levels = design
    )

    # Empirical Bayes moderated t-tests, without robust or trend estimation.
    tab <- speckle::propeller.ttest(
      prop.list = props,
      design = design,
      contrasts = contrast,
      robust = FALSE,
      trend = FALSE,
      sort = FALSE
    )
    tab <- data.frame(
      pb_celltype = rownames(tab), tab,
      row.names = NULL, check.names = FALSE
    )
    stopifnot(
      "P.Value" %in% names(tab),
      nrow(tab) == 3L,
      all(is.finite(tab$P.Value))
    )

    # Supplementary correction across the three cell types within each age.
    tab$FDR_within_stage <- p.adjust(tab$P.Value, method = "BH")
    tab$FDR <- NULL
    tab$stage <- stage_i
    tab$transform <- transform_i
    tab <- left_join(tab, da_effects, by = c("stage", "pb_celltype"))
    da_results[[paste(stage_i, transform_i, sep = "_")]] <- tab
  }
}

# 5. Correct across three cell types x two ages (six tests)


da_results <- bind_rows(da_results) %>%
  group_by(transform) %>%
  mutate(FDR_across_6_tests = p.adjust(P.Value, method = "BH")) %>%
  ungroup()

primary_results <- filter(da_results, transform == "logit")
sensitivity_results <- filter(da_results, transform == "asin")

# 6. Export observed proportions, summary statistics, and test results

write.csv(da_counts, "mouse_level_cell_counts_and_proportions.csv", row.names = FALSE)
write.csv(da_summary, "group_mean_SD_SEM.csv", row.names = FALSE)
write.csv(primary_results, "PRIMARY_logit_P3_P10_propeller.csv", row.names = FALSE)
write.csv(sensitivity_results, "SENSITIVITY_asin_P3_P10_propeller.csv", row.names = FALSE)

print(primary_results)

# Display the current software environment without writing another file.
sessionInfo()
#R version 4.4.2 (2024-10-31)
#Platform: x86_64-apple-darwin20
#Running under: macOS Ventura 13.0
#[1] speckle_1.6.0            scales_1.4.0             tidyr_1.3.2              stringr_1.6.0           
#[5] RColorBrewer_1.1-3       dplyr_1.2.0              future_1.69.0            hdf5r_1.3.12            
#[9] Seurat_5.4.0             variancePartition_1.42.0 BiocParallel_1.40.2      limma_3.62.2            
#[13] ggplot2_4.0.2            SeuratObject_5.3.0       sp_2.2-1                
