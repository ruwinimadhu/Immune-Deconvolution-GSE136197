setwd("C:/Users/user/OneDrive/Desktop/e")


# Install devtools 
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install EPIC from GitHub
devtools::install_github("GfellerLab/EPIC")

library(EPIC)

# Load TPM data 
tpm <- read.csv("GSE136197_Canine_RNASeq.csv", row.names = 1)
tpm <- t(tpm)

# Run EPIC
res <- EPIC(bulk = tpm)

# Save output
write.csv(res$cellFractions, "EPIC_CellFractions.csv")


# After running EPIC
res <- EPIC(bulk = tpm)


head(res$fit.gof)


# Load libraries
library(tidyverse)
library(ggpubr)
library(effectsize)

# Load EPIC output
epic_raw <- read_csv("EPIC_CellFractions.csv", show_col_types = FALSE)

# Preprocess: rename first column and reshape
epic_clean <- epic_raw %>%
  rename(SampleID = 1) %>%
  mutate(group = if_else(str_detect(SampleID, "C$|E$"), "Normal", "Tumor")) %>%
  pivot_longer(
    cols = -c(SampleID, group),
    names_to = "CellType",
    values_to = "CellFraction"
  )



# group assignment: mark all samples starting with a digit (1A, 2B...) as Tumor, and N01, N02... as Normal
epic <- epic %>%
  rename(Sample = 1) %>%
  mutate(group = case_when(
    str_detect(Sample, "N\\d+") ~ "Normal",
    str_detect(Sample, "^\\d") ~ "Tumor",
    TRUE ~ NA_character_
  ))


epic <- epic %>% filter(!is.na(group))


epic_long <- epic %>%
  pivot_longer(cols = -c(Sample, group), names_to = "CellType", values_to = "CellFraction") %>%
  filter(!is.na(CellFraction))


library(readr)
library(dplyr)
library(tidyr)

# Load expression fractions 
epic <- read_csv("EPIC_CellFractions.csv")  # column 1 = sample ID

# Load metadata 
meta <- read_tsv("GSE136197_clean_metadata.txt")

# Rename first column in EPIC to match metadata
colnames(epic)[1] <- "sample"

# Merge with metadata
epic_merged <- left_join(epic, meta, by = "sample")

# Convert to long format
epic_long <- epic_merged %>%
  pivot_longer(
    cols = Bcells:otherCells,
    names_to = "CellType",
    values_to = "CellFraction"
  )

library(ggplot2)
library(ggpubr)

# Remove NAs
epic_filtered <- epic_long %>%
  filter(!is.na(group))

# Plot with Wilcoxon test
ggplot(epic_filtered, aes(x = group, y = CellFraction, fill = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Normal", "Tumor")),
    label = "p.format",
    label.y = 1.05
  ) +
  facet_wrap(~ CellType, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Immune Cell Fractions (EPIC): Normal vs Tumor",
    x = "Group", y = "Cell Fraction"
  ) +
  scale_fill_manual(values = c("Normal" = "#F8766D", "Tumor" = "#00BFC4"))


# Save as PNG (300 dpi)
ggsave("Immune_Cell_Fractions_EPIC.png", dpi = 300, width = 10, height = 7, units = "in")

# Save as PDF (high quality, dpi ignored)
ggsave("Immune_Cell_Fractions_EPIC.pdf", width = 10, height = 7, units = "in")


library(dplyr)
library(ggplot2)
library(ggpubr)

# Calculate p-values for each cell type
pvals <- epic_clean %>%
  group_by(CellType) %>%
  summarise(p = wilcox.test(CellFraction ~ group)$p.value) %>%
  mutate(p_label = paste0("p = ", signif(p, 3)))

# Merge back with main data
epic_annot <- epic_clean %>%
  left_join(pvals, by = "CellType") %>%
  mutate(CellType_facet = paste0(CellType, "\n", p_label))

# Plot with custom facet labels
ggplot(epic_annot, aes(x = group, y = CellFraction, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
  facet_wrap(~ CellType_facet, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#F8766D", "Tumor" = "#00BFC4")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Immune Cell Fractions (EPIC): Normal vs Tumor",
    y = "Cell Fraction",
    x = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )


library(dplyr)
library(tidyr)
library(ggpubr)

# Summary table with group-wise means and Wilcoxon p-values
summary_table <- epic_long %>%
  filter(group %in% c("Normal", "Tumor")) %>%
  group_by(CellType, group) %>%
  summarise(
    mean_fraction = mean(CellFraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = group, values_from = mean_fraction) %>%
  left_join(
    epic_long %>%
      filter(group %in% c("Normal", "Tumor")) %>%
      group_by(CellType) %>%
      summarise(
        p_value = wilcox.test(CellFraction ~ group)$p.value,
        .groups = "drop"
      ),
    by = "CellType"
  )

# View it
print(summary_table)
# save it
write.csv(summary_table, "EPIC_celltype_means_pvalues.csv", row.names = FALSE)





# Load libraries
library(tidyverse)
library(readr)
library(pheatmap)

# Load EPIC fractions CSV
epic <- read_csv("EPIC_CellFractions.csv")

# Check column names
print("EPIC column names BEFORE renaming:")
print(colnames(epic))

# Rename first column if it's unnamed (...1)
if ("...1" %in% colnames(epic)) {
  colnames(epic)[which(colnames(epic) == "...1")] <- "sample"
}

# Confirm fix
print("EPIC column names AFTER renaming:")
print(colnames(epic))

# Load metadata
metadata <- read_delim("GSE136197_clean_metadata.txt", delim = "\t")

# Show number of matching samples
common_samples <- intersect(epic$sample, metadata$sample)
print(paste("Number of common samples:", length(common_samples)))

# Filter to common samples only
epic_filtered <- epic %>% filter(sample %in% common_samples)
metadata_filtered <- metadata %>% filter(sample %in% common_samples)

# Match order
epic_filtered <- epic_filtered %>% arrange(match(sample, metadata_filtered$sample))
metadata_filtered <- metadata_filtered %>% arrange(match(sample, epic_filtered$sample))

# Final check
print("Are sample orders matching?")
print(all(epic_filtered$sample == metadata_filtered$sample))

# Transpose matrix
epic_mat <- epic_filtered %>% column_to_rownames("sample") %>% as.matrix()
fractions <- t(epic_mat)

# Create annotation for heatmap
annotation_col <- data.frame(Group = metadata_filtered$group)
rownames(annotation_col) <- metadata_filtered$sample


pheatmap(
  fractions,
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 6,
  main = "EPIC Immune Cell Fractions (Z-score)",
  filename = "EPIC_Immune_Heatmap_Scaled.pdf",  # Save as file
  width = 8,
  height = 6
)

pheatmap(
  fractions,
  scale = "row",
  annotation_col = annotation_col,
  main = "EPIC Immune Cell Fractions (Z-score)"
)


# Save PDF
pheatmap(
  fractions,
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 6,
  main = "EPIC Immune Cell Fractions (Z-score)",
  filename = "EPIC_Immune_Heatmap_Scaled.pdf",
  width = 8,
  height = 6
)

# Save JPG
jpeg("EPIC_Immune_Heatmap_Scaled.jpg", width = 1200, height = 900, res = 150)
pheatmap(
  fractions,
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 6,
  main = "EPIC Immune Cell Fractions (Z-score)"
)
dev.off()


library(effsize)

# EPIC output 
epic <- read.csv("EPIC_CellFractions.csv")
colnames(epic)[1] <- "sample"  # Rename first column

# Metadata file 
metadata <- read.delim("GSE136197_clean_metadata.txt")

# Filter EPIC to only include samples in metadata
epic <- epic %>% filter(sample %in% metadata$sample)

# Match sample order
metadata <- metadata %>% filter(sample %in% epic$sample)
epic <- epic[match(metadata$sample, epic$sample), ]
stopifnot(all(metadata$sample == epic$sample))  # Ensure perfect match

# Set sample names as rownames and remove column
rownames(epic) <- epic$sample
epic_mat <- epic[, -1]  # Only cell types

# CALCULATE COHEN’S D (Effect Size)
cohen_result <- cohen.d(epic$Bcells ~ metadata$group)
print(cohen_result)

# loop over all cell types
cell_types <- colnames(epic_mat)
effect_sizes <- lapply(cell_types, function(cell) {
  d <- cohen.d(epic[[cell]] ~ metadata$group)
  data.frame(CellType = cell, d = d$estimate, magnitude = d$magnitude)
})
effect_sizes_df <- do.call(rbind, effect_sizes)
print(effect_sizes_df)

# Save to CSV
write.csv(effect_sizes_df, "EPIC_EffectSizes_CohenD.csv", row.names = FALSE)


library(tidyverse)

epic <- read.csv("EPIC_CellFractions.csv")  # 89 rows, immune cell fractions
metadata <- read.delim("GSE136197_clean_metadata.txt")  # 48 rows, histology info



colnames(epic)[1] <- "sample"
merged <- inner_join(metadata, epic, by = "sample")

long_data <- merged %>%
  pivot_longer(cols = Bcells:otherCells, names_to = "cell_type", values_to = "fraction")

ggplot(long_data, aes(x = histology, y = fraction, fill = histology)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), size = 0.8, alpha = 0.4) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1.2, "lines")
  ) +
  labs(
    title = "Immune Cell Infiltration by Tumor Histology",
    x = "Tumor Histology",
    y = "Estimated Cell Fraction"
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave("Immune_Histology_Boxplot.pdf", width = 10, height = 7)
ggsave("Immune_Histology_Boxplot.jpg", width = 10, height = 7, dpi = 300)






library(dplyr)
library(FSA)

# Transform to long format 
epic_long <- epic_mat %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "CellType", values_to = "Fraction") %>%
  left_join(metadata[, c("sample", "histology")], by = "sample")

# Kruskal-Wallis test only for cell types with ≥2 histology groups
kruskal_results <- epic_long %>%
  group_by(CellType) %>%
  filter(n_distinct(histology) > 1) %>%
  summarise(
    p_value = tryCatch(
      kruskal.test(Fraction ~ histology)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH"))

# View results
print(kruskal_results)



