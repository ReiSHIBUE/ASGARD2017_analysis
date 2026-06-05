### S25_transect_compare.R
### ASGARD 2017 Survey — Class composition compared across transects
### トランセクト別 Alpha/Gamma/Bacteroidia 比較（DBO3 vs 他）
###
### REQUIRES (from 00_setup.R, S01, S02):
###   asgard_filtered  - 181×258 ASV proportion matrix
###   meta_asgard      - metadata
###   clusnum11        - 11-cluster assignments
###   hier_levels_11   - ordered cluster names
###   shorternames     - taxonomy short labels
###   fullnameboot     - taxonomy + bootstrap df
###
### OUTPUT:
###   output/survey/transect_compare/ASGARD_transect_class_compare.pdf
###   output/survey/transect_compare/ASGARD_transect_class_compare.csv
###   output/survey/transect_compare/ASGARD_dbo3_vs_others_test.csv

library(tidyverse)

# ==============================================================================
# Section 1: Assign each sample to a transect
# ==============================================================================

m <- meta_asgard %>% rownames_to_column("Sample")
m$cluster11 <- factor(as.character(clusnum11[m$Sample]), levels = hier_levels_11)

m$transect <- ifelse(grepl("^DBO3", m$station), "DBO3",
              ifelse(grepl("^DBO2", m$station), "DBO2",
              ifelse(grepl("^DBO1", m$station), "DBO1",
              ifelse(grepl("^CBE",  m$station), "CBE",
              ifelse(grepl("^CBW",  m$station), "CBW",
              ifelse(grepl("^CB[0-9]", m$station), "CB",
              ifelse(grepl("^CNL",  m$station), "CNL",
              ifelse(grepl("^CPL",  m$station), "CPL",
              ifelse(grepl("^IL",   m$station), "IL",
              ifelse(grepl("^KL",   m$station), "KL",
              ifelse(grepl("^CL",   m$station), "CL",
              ifelse(grepl("^SLY",  m$station), "SLY", "Other"))))))))))))
m$transect <- factor(m$transect, levels = c("SLY","CB","CBE","CBW","CPL",
                                            "DBO1","DBO2","DBO3","CNL","IL",
                                            "KL","CL","Other"))

cat("=== Samples per transect ===\n")
print(table(m$transect))

# Mark DBO3 vs others
m$is_dbo3 <- factor(ifelse(m$transect == "DBO3", "DBO3", "Other transects"),
                    levels = c("Other transects", "DBO3"))

# ==============================================================================
# Section 2: Compute per-sample Class composition
# ==============================================================================

asv_idx   <- match(colnames(asgard_filtered), shorternames)
strip_boot <- function(x) sub("\\s*\\(\\d+\\)\\s*$", "", x)
asv_class <- strip_boot(fullnameboot$Class[asv_idx])
names(asv_class) <- colnames(asgard_filtered)

long_c <- as.data.frame(asgard_filtered) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "Abundance") %>%
  mutate(class = asv_class[ASV]) %>%
  inner_join(m %>% select(Sample, transect, is_dbo3, cluster11), by = "Sample") %>%
  group_by(Sample) %>% mutate(RA = Abundance / sum(Abundance)) %>% ungroup()

# Per-sample, per-class total
sample_class <- long_c %>%
  group_by(Sample, transect, is_dbo3, cluster11, class) %>%
  summarise(sample_pct = sum(RA) * 100, .groups = "drop")

# Keep Big-3 only for the focused comparison
big3 <- c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia")
big3_df <- sample_class %>% filter(class %in% big3)
big3_df$class <- factor(big3_df$class, levels = big3)

# ==============================================================================
# Section 3: Plots — boxplots per transect
# ==============================================================================

dir.create(here::here("output", "survey", "transect_compare"),
           showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "transect_compare",
                     "ASGARD_transect_class_compare.pdf"),
    width = 14, height = 10)

# Page 1: Big-3 by transect
print(
  ggplot(big3_df, aes(x = transect, y = sample_pct, fill = transect)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.4) +
    facet_wrap(~ class, nrow = 1, scales = "free_y") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    labs(title = "Big-3 class composition by transect",
         subtitle = "DBO3 highlighted; each point = one sample",
         x = "Transect", y = "Relative abundance (%)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
)

# Page 2: Big-3 — DBO3 vs Other transects (binary)
print(
  ggplot(big3_df, aes(x = is_dbo3, y = sample_pct, fill = is_dbo3)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.6, alpha = 0.4) +
    facet_wrap(~ class, nrow = 1) +
    scale_fill_manual(values = c("Other transects" = "gray70", "DBO3" = "#E41A1C"),
                      guide = "none") +
    labs(title = "Big-3 class composition: DBO3 vs other transects",
         x = NULL, y = "Relative abundance (%)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12))
)

# Page 3: Alpha share (Alpha / (Alpha + Gamma)) per transect
alpha_share <- sample_class %>%
  filter(class %in% c("Alphaproteobacteria", "Gammaproteobacteria")) %>%
  pivot_wider(names_from = class, values_from = sample_pct, values_fill = 0) %>%
  rename(Alpha = Alphaproteobacteria, Gamma = Gammaproteobacteria) %>%
  mutate(alpha_share = ifelse((Alpha + Gamma) > 0,
                              Alpha / (Alpha + Gamma) * 100, NA))

print(
  ggplot(alpha_share, aes(x = transect, y = alpha_share, fill = transect)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.4, alpha = 0.4) +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    labs(title = "Alpha share (= Alpha / (Alpha + Gamma) %) by transect",
         subtitle = "High = Thalassiosira-like signal; Low = Chaetoceros-like signal",
         x = "Transect", y = "Alpha share (%)") +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11),
          axis.text.x   = element_text(angle = 45, hjust = 1))
)

# Page 4: Alpha share — DBO3 vs Other (binary)
print(
  ggplot(alpha_share, aes(x = is_dbo3, y = alpha_share, fill = is_dbo3)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.6, alpha = 0.4) +
    scale_fill_manual(values = c("Other transects" = "gray70", "DBO3" = "#E41A1C"),
                      guide = "none") +
    labs(title = "Alpha share: DBO3 vs other transects",
         x = NULL, y = "Alpha share (%)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))
)

dev.off()

# ==============================================================================
# Section 4: Statistical tests — DBO3 vs other transects
# ==============================================================================

test_results <- list()
for (cl in big3) {
  d <- sample_class %>% filter(class == cl)
  wt <- wilcox.test(sample_pct ~ is_dbo3, data = d)
  test_results[[cl]] <- data.frame(
    variable = cl,
    n_DBO3   = sum(d$is_dbo3 == "DBO3"),
    n_other  = sum(d$is_dbo3 == "Other transects"),
    median_DBO3   = round(median(d$sample_pct[d$is_dbo3 == "DBO3"]), 2),
    median_other  = round(median(d$sample_pct[d$is_dbo3 == "Other transects"]), 2),
    W_statistic   = round(wt$statistic, 1),
    p_value       = round(wt$p.value, 4)
  )
}
# Alpha share test
as_d <- alpha_share %>% filter(!is.na(alpha_share))
wt <- wilcox.test(alpha_share ~ is_dbo3, data = as_d)
test_results[["Alpha_share"]] <- data.frame(
  variable = "Alpha_share",
  n_DBO3   = sum(as_d$is_dbo3 == "DBO3"),
  n_other  = sum(as_d$is_dbo3 == "Other transects"),
  median_DBO3   = round(median(as_d$alpha_share[as_d$is_dbo3 == "DBO3"]), 2),
  median_other  = round(median(as_d$alpha_share[as_d$is_dbo3 == "Other transects"]), 2),
  W_statistic   = round(wt$statistic, 1),
  p_value       = round(wt$p.value, 4)
)

test_df <- bind_rows(test_results) %>%
  mutate(sig = ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01,  "**",
              ifelse(p_value < 0.05,  "*", "ns"))))

cat("\n=== Wilcoxon rank-sum: DBO3 vs other transects ===\n")
print(test_df, row.names = FALSE)

write.csv(test_df,
  here::here("output", "survey", "transect_compare",
             "ASGARD_dbo3_vs_others_test.csv"),
  row.names = FALSE)

# Save per-sample data
write.csv(sample_class,
  here::here("output", "survey", "transect_compare",
             "ASGARD_transect_class_compare.csv"),
  row.names = FALSE)

message("\nS25_transect_compare.R: done.")
message("  PDF: output/survey/transect_compare/ASGARD_transect_class_compare.pdf")
message("  CSV: output/survey/transect_compare/ASGARD_transect_class_compare.csv")
message("  CSV: output/survey/transect_compare/ASGARD_dbo3_vs_others_test.csv")
