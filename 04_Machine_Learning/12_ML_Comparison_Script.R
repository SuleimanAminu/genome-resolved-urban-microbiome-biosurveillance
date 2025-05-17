# plot_model_performance.R
set.seed(43)
# === Load Libraries ===
library(tidyverse)
library(ggplot2)

#########################################################################
# Load your long-format results
df <- read.csv("Env_model_foldwise_scores.csv")
names(df) <- make.names(names(df))

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = c(Accuracy, F1.Macro, Kappa),
               names_to = "Metric", 
               values_to = "Score")

# Optional: rename metrics to match desired visual
df_long$Metric <- recode(df_long$Metric,
                         "F1.Macro" = "F1 Macro",
                         "Accuracy" = "Mean Balanced Accuracy")

# Reorder models by Kappa for cleaner plot
model_order <- df_long %>%
  filter(Metric == "Kappa") %>%
  group_by(Model) %>%
  summarise(avg = mean(Score)) %>%
  arrange(desc(avg)) %>%
  pull(Model)

df_long$Model <- factor(df_long$Model, levels = model_order)

top_models <- df_long %>%
  filter(Metric == "Kappa") %>%
  group_by(Model) %>%
  summarise(mean_kappa = mean(Score)) %>%
  arrange(desc(mean_kappa)) %>%
  slice_head(n = 3) %>%
  pull(Model)

df_long$Highlight <- ifelse(df_long$Model %in% top_models, "Top 3", "Other")


# Identify top 3 models by average Kappa
top_models <- df_long %>%
  filter(Metric == "Kappa") %>%
  group_by(Model) %>%
  summarise(mean_kappa = mean(Score)) %>%
  arrange(desc(mean_kappa)) %>%
  slice_head(n = 3) %>%
  pull(Model)

# Add a Highlight column
df_long$Highlight <- ifelse(df_long$Model %in% top_models, "Top 3", "Other")

# === Save to TIFF ===
tiff("model_comparison_plot.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")  ####Suppelementary Figure 5

ggplot(df_long, aes(x = Score, y = Model)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", fill = "transparent") +
  stat_summary(aes(fill = Highlight), fun = mean, geom = "point", shape = 22, size = 4) +
  scale_fill_manual(values = c("Top 3" = "deeppink", "Other" = "blue")) +
  facet_wrap(~ Metric, scales = "free_x") +
  theme_minimal(base_size = 14) +
  labs(x = "Score", y = "") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8, color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 10)),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 1),
    legend.position = "none"
  )

dev.off()


