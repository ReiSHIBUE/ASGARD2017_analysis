### S20_stratification.R
### ASGARD 2017 Survey — Water column stratification analysis
### 水柱成層解析（温度・塩分・密度・DO・NO3・Chl の鉛直分布）
###
### REQUIRES (from 00_setup.R, S01, S02):
###   meta_asgard      - metadata 181 samples
###   clusnum11        - 11-cluster assignments
###   hier_levels_11   - ordered cluster names
###   cc11             - 11-colour palette
###
### PRODUCES:
###   stratification_summary - df of Delta sigma-t per station
###
### OUTPUT:
###   output/survey/stratification/ASGARD_vertical_profiles.pdf
###   output/survey/stratification/ASGARD_stratification_summary.csv

library(tidyverse)

# ==============================================================================
# Section 1: Data preparation and sigma-t calculation
# データ準備 + ポテンシャル密度 (sigma-t, EOS-80) の計算
# ==============================================================================

m <- meta_asgard %>% rownames_to_column("Sample")
m$cluster11  <- factor(as.character(clusnum11[m$Sample]), levels = hier_levels_11)
m$division   <- factor(sub("^(.).*", "\\1", as.character(m$cluster11)),
                       levels = c("A", "B", "C"))
m$depth_type <- factor(m$depth_type, levels = c("surf", "mid", "bottom"))

# EOS-80 potential density (sigma-t, kg/m^3) — surface pressure
sigma_theta <- function(T, S) {
  999.842594 + 6.793952e-2*T - 9.095290e-3*T^2 + 1.001685e-4*T^3 -
    1.120083e-6*T^4 + 6.536332e-9*T^5 +
    (0.824493 - 4.0899e-3*T + 7.6438e-5*T^2 - 8.2467e-7*T^3 + 5.3875e-9*T^4)*S +
    (-5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T^2)*S^1.5 +
    4.8314e-4*S^2 - 1000
}
m$sigma_t <- with(m, sigma_theta(temp, salinity))

# Region from station prefix
m$region <- ifelse(grepl("^CBE", m$station),    "CBE",
            ifelse(grepl("^CBW", m$station),    "CBW",
            ifelse(grepl("^DBO", m$station),    "DBO",
            ifelse(grepl("^CB[0-9]", m$station),"CB",
            ifelse(grepl("^SLY", m$station),    "SLY",
            ifelse(grepl("^CPL", m$station),    "CPL", "Other"))))))

# ==============================================================================
# Section 2: Stratification metric per station (Delta sigma-t)
# ステーション毎の成層度 (sigma-t 差)
# ==============================================================================

stratification_summary <- m %>%
  group_by(station) %>%
  summarise(
    n_depths     = n(),
    depth_min    = round(min(depth_m, na.rm = TRUE), 1),
    depth_max    = round(max(depth_m, na.rm = TRUE), 1),
    depth_range  = round(diff(range(depth_m, na.rm = TRUE)), 1),
    temp_surf    = round(temp[which.min(depth_m)], 2),
    temp_bottom  = round(temp[which.max(depth_m)], 2),
    sal_surf     = round(salinity[which.min(depth_m)], 2),
    sal_bottom   = round(salinity[which.max(depth_m)], 2),
    sigma_surf   = round(sigma_t[which.min(depth_m)], 3),
    sigma_bottom = round(sigma_t[which.max(depth_m)], 3),
    delta_sigma  = round(sigma_bottom - sigma_surf, 3),
    .groups = "drop"
  ) %>%
  filter(n_depths >= 2) %>%
  mutate(strat_class = case_when(
    delta_sigma >  1   ~ "strongly stratified",
    delta_sigma >  0.5 ~ "weakly stratified",
    TRUE               ~ "well-mixed"
  ))

cat("\n=== Stratification class per station ===\n")
print(table(stratification_summary$strat_class))

# ==============================================================================
# Section 3: Vertical profile plots
# ==============================================================================

plot_df <- m %>%
  select(Sample, station, lat, lon, depth_m, depth_type, cluster11, division, region,
         temp, salinity, DO, sigma_t, `NO3(uM)`, `chl (ug/l)`) %>%
  rename(NO3 = `NO3(uM)`, chl = `chl (ug/l)`) %>%
  pivot_longer(cols = c(temp, salinity, sigma_t, DO, NO3, chl),
               names_to = "variable", values_to = "value") %>%
  filter(!is.na(value), !is.na(depth_m))

plot_df$variable <- factor(plot_df$variable,
  levels = c("temp", "salinity", "sigma_t", "DO", "NO3", "chl"),
  labels = c("Temperature (°C)", "Salinity (PSU)", "Sigma-t (kg/m³)",
             "DO (µmol/kg)", "NO3 (µM)", "Chl a (µg/L)"))

dir.create(here::here("output", "survey", "stratification"),
           showWarnings = FALSE, recursive = TRUE)

pdf(file = here::here("output", "survey", "stratification",
                     "ASGARD_vertical_profiles.pdf"),
    width = 22, height = 12)

# Page 1: vertical profiles colored by Division (A/B/C)
print(
  ggplot(plot_df, aes(x = value, y = depth_m, color = division)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_y_reverse() +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +
    scale_color_manual(values = c("A" = "#E31A1C", "B" = "#33A02C", "C" = "#1F78B4")) +
    labs(title = "Vertical profiles by Division (A/B/C)",
         x = NULL, y = "Depth (m)", color = "Division") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
)

# Page 2: vertical profiles colored by 11-cluster
print(
  ggplot(plot_df, aes(x = value, y = depth_m, color = cluster11)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_y_reverse() +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +
    scale_color_manual(values = cc11) +
    labs(title = "Vertical profiles by 11-cluster",
         x = NULL, y = "Depth (m)", color = "Cluster") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
)

# Page 3: per-station profile lines, colored by transect region
print(
  ggplot(plot_df,
         aes(x = value, y = depth_m, group = station, color = region)) +
    geom_path(alpha = 0.5) +
    geom_point(alpha = 0.7, size = 1) +
    scale_y_reverse() +
    facet_wrap(~ variable, scales = "free_x", nrow = 2) +
    labs(title = "Per-station vertical profiles by transect region",
         x = NULL, y = "Depth (m)", color = "Region") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
)

# Page 4: stratification scatter (Delta sigma-t per station)
print(
  ggplot(stratification_summary, aes(x = depth_range, y = delta_sigma,
                                     color = strat_class)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = station), size = 2.5, vjust = -1, color = "black") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("well-mixed"            = "#377EB8",
                                  "weakly stratified"    = "#FFB000",
                                  "strongly stratified"  = "#E41A1C")) +
    labs(title = "Stratification index per station (Δ sigma-t = bottom − surface)",
         subtitle = "Dashed lines at 0.5 and 1.0 kg/m³",
         x = "Water-column depth range (m)",
         y = "Δ Sigma-t (kg/m³)",
         color = "Stratification class") +
    theme_bw(base_size = 12)
)

# Page 5: vertical profiles faceted by 11-cluster × variable
# 11 クラスター × 6 変数のグリッド (各セルが独立した鉛直分布)
print(
  ggplot(plot_df, aes(x = value, y = depth_m)) +
    geom_point(aes(color = cluster11), alpha = 0.7, size = 1.3) +
    scale_y_reverse() +
    facet_grid(variable ~ cluster11, scales = "free_x") +
    scale_color_manual(values = cc11, guide = "none") +
    labs(title = "Vertical profiles by 11-cluster (rows: variable, cols: cluster)",
         x = NULL, y = "Depth (m)") +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(face = "bold", size = 14))
)

# Page 6: points per variable × cluster, x-axis = depth_type (surf/mid/bottom)
# 6 variables × 11 clusters のグリッド、各セル内で depth_type 別の点プロット
print(
  ggplot(plot_df %>% filter(!is.na(depth_type)),
         aes(x = depth_type, y = value, color = cluster11)) +
    geom_jitter(width = 0.2, size = 1.2, alpha = 0.7) +
    facet_grid(variable ~ cluster11, scales = "free_y") +
    scale_color_manual(values = cc11, guide = "none") +
    labs(title = "Variable distribution by depth-type, within each cluster",
         subtitle = "Rows: variable, cols: cluster; x-axis = surface/mid/bottom",
         x = "Depth type", y = NULL) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(face = "bold", size = 14))
)

dev.off()

# ==============================================================================
# Section 4: CSV summary
# ==============================================================================

write.csv(stratification_summary,
          here::here("output", "survey", "stratification",
                     "ASGARD_stratification_summary.csv"),
          row.names = FALSE)

message("\nS20_stratification.R: done.")
message("  Stations with n>=2 depths: ", nrow(stratification_summary))
message("  Strongly stratified (Δ > 1.0): ",
        sum(stratification_summary$strat_class == "strongly stratified"))
message("  Weakly stratified (0.5 < Δ <= 1.0): ",
        sum(stratification_summary$strat_class == "weakly stratified"))
message("  Well-mixed (Δ <= 0.5): ",
        sum(stratification_summary$strat_class == "well-mixed"))
message("  PDF: output/survey/stratification/ASGARD_vertical_profiles.pdf")
message("  CSV: output/survey/stratification/ASGARD_stratification_summary.csv")
