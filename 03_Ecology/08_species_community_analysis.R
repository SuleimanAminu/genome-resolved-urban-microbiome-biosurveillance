################    This script was developed by Suleiman and AbdulAziz; 4th December, 2024   ################################

#######################################  Setting Seed for the Analysis
set.seed(45)
### Loading Libraries
library(vegan)
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(dplyr)
library(FSA) 
library(tidyverse)
library(reshape2)
library(parallel)
library(ggrepel)
library(ComplexUpset)
library(tidyr)
library(UpSetR)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(indicspecies)
library(Matrix)

##############################################################################################

# Load data files
species_contaminants_corrected <- read.csv("Env_Pathogenic_species.csv", row.names = 1)
metadata <- read.csv("Merged_metadata_all.csv")

# Transpose the species data so that samples are rows
species_contaminants_corrected <- t(species_contaminants_corrected)

rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[rownames(species_contaminants_corrected), ]


#  Save the Aligned Metadata
metadata_export <- metadata
metadata_export$Sample_ID <- rownames(metadata_export)

# Save
write.csv(metadata_export, "Metadata_Aligned_to_species_corrected.csv", row.names = FALSE)

#################################################################
# Calculate Shannon diversity index
alpha_diversity_shannon <- diversity(species_contaminants_corrected, index = "shannon")

# Calculate Simpson diversity index
alpha_diversity_simpson <- diversity(species_contaminants_corrected, index = "simpson")

# Calculate Richness (number of observed species per sample)
alpha_diversity_richness <- specnumber(species_contaminants_corrected)

# Combine diversity indices into a data frame for plotting
alpha_diversity_df <- data.frame(
  Sample_ID = rownames(species_contaminants_corrected),
  Shannon = alpha_diversity_shannon,
  Simpson = alpha_diversity_simpson,
  Richness = alpha_diversity_richness
)

# Merge with metadata for grouping information
alpha_diversity_df <- merge(alpha_diversity_df, metadata, by.x = "Sample_ID", by.y = "Sample_ID")


# Print the first few rows of the data frame to verify the content
print(head(alpha_diversity_df))
summary(alpha_diversity_df)
# Print the unique group names to verify available groups
print(unique(alpha_diversity_df$Group))

# Get the unique group names
unique_groups <- unique(alpha_diversity_df$Group)


# Print the first few rows of the data frame to verify the content
print(head(alpha_diversity_df))

# Print the unique group names to verify available groups
print(unique(alpha_diversity_df$Group))

# Get the unique group names
unique_groups <- unique(alpha_diversity_df$Group)

###Save plots ################
group_colors <- c(
  "Ambulance" = "#1f77b4",             # Deep blue
  "Hospital_environment" = "#aec7e8",  # Light blue
  "Hospital_sewage" = "#ff7f0e",       # Orange
  "Public_transport" = "#fdbf6f"       # Light orange/yellow
)

library(ggplot2)

# Start TIFF device
tiff("Species_Shannon.tiff", width = 4, height = 4, units = "in", res = 600) ##### Figure 2H

# Plot
p1 <- ggplot(alpha_diversity_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +  # Apply consistent colors
  theme_minimal() +
  labs(title = "", x = "", y = "Shannon Index") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.ticks.length = unit(0.25, "cm")
  )

# Print & save
print(p1)
dev.off()





tiff("Species_Richness.tiff", width = 4, height = 4, units = "in", res = 600)   ###### Figure 2I
# Updated Richness plot
p2 <- ggplot(alpha_diversity_df, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  labs(title = "", x = "", y = "Richness") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.2),  # Increase tick thickness
    axis.ticks.length = unit(0.25, "cm")  # Increase tick length
  )

p2

dev.off()


tiff("Species_Simpson.tiff", width = 4, height = 4, units = "in", res = 600)   ##### Figure 2J

# Updated Simpson Diversity plot
p3 <- ggplot(alpha_diversity_df, aes(x = Group, y = Simpson, fill = Group)) +
  scale_fill_manual(values = group_colors) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "", x = "", y = "Simpson Index") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.2),  # Increase tick thickness
    axis.ticks.length = unit(0.25, "cm")  # Increase tick length
  )

p3
dev.off()

##########STATISTICS FOR THE ECOLOGY DIVERSITY##########

### 1. Kruskal-Wallis Results
kw_shannon <- kruskal.test(Shannon ~ Group, data = alpha_diversity_df)
kw_simpson <- kruskal.test(Simpson ~ Group, data = alpha_diversity_df)
kw_richness <- kruskal.test(Richness ~ Group, data = alpha_diversity_df)

# Format Kruskal-Wallis results
kw_df <- data.frame(
  Test = "Kruskal-Wallis",
  Metric = c("Shannon", "Simpson", "Richness"),
  Comparison = NA,
  Z = NA,
  P.unadj = NA,
  P.adj = NA,
  Chi_squared = c(kw_shannon$statistic, kw_simpson$statistic, kw_richness$statistic),
  df = c(kw_shannon$parameter, kw_simpson$parameter, kw_richness$parameter),
  p_value = c(kw_shannon$p.value, kw_simpson$p.value, kw_richness$p.value)
)

### 2. Dunn’s Post-Hoc Results
dunn_shannon <- dunnTest(Shannon ~ Group, data = alpha_diversity_df, method = "bh")$res
dunn_shannon$Metric <- "Shannon"

dunn_simpson <- dunnTest(Simpson ~ Group, data = alpha_diversity_df, method = "bh")$res
dunn_simpson$Metric <- "Simpson"

dunn_richness <- dunnTest(Richness ~ Group, data = alpha_diversity_df, method = "bh")$res
dunn_richness$Metric <- "Richness"

# Combine all Dunn results
dunn_all <- bind_rows(dunn_shannon, dunn_simpson, dunn_richness)
dunn_all <- dunn_all[, c("Metric", "Comparison", "Z", "P.unadj", "P.adj")]
dunn_all$Test <- "Dunn's Test"

# Add missing columns to Dunn to match Kruskal-Wallis structure
dunn_all$Chi_squared <- NA
dunn_all$df <- NA
dunn_all$p_value <- NA

# Reorder to match
dunn_all <- dunn_all[, colnames(kw_df)]

# Combine into one dataframe
combined_output <- bind_rows(kw_df, dunn_all)

# Save
write.csv(combined_output, "AlphaDiversity_KW_Dunn_Combined.csv", row.names = FALSE)  ##### Supplementary Table 4 (Statistical_Analysis_Species_Diversity)


############################################################## ANOSIM ANALYSIS ##################
# Replace 'grouping' with the actual column name in your 'md' metadata object
grouping <- metadata$Group

# Compute dissimilarity matrix using Bray-Curtis distance
dist_matrix <- vegdist(species_contaminants_corrected, method = "bray", binary = TRUE, parallel = detectCores() - 1) # The default distance method is bray(preferable), if it doesnt give good values then use gower or any other.


# Perform ANOSIM test
anosim_result <- anosim(dist_matrix, grouping, permutations = 999, parallel = detectCores() - 1)

# Print ANOSIM result
print(anosim_result)

# Set color vector based on ordering of the groups + "Between"
ordered_groups <- levels(as.factor(grouping))  # Ensure consistent order
anosim_colors <- c("black",                     # Between
                   group_colors["Ambulance"],
                   group_colors["Hospital_environment"],
                   group_colors["Hospital_sewage"],
                   group_colors["Public_transport"])
tiff("ANOSIM_species.tiff", width = 9, height = 5, units = "in", res = 600) 

plot(anosim_result,
     col = anosim_colors,
     ylab = "Dissimilarity Rank Value", xlab = "",
     cex.lab = 0.8, cex.axis = 0.8,
     lwd = 1)

box(lwd = 5)  # Border thickness
dev.off()





###################################################################################################
# Perform a Non-metric Multidimensional Scaling (NMDS) which allows you to collapse all these species axes into 2 to plot in cartesian space in order to visualize the differences between samples and sites
NMDS <- metaMDS(species_contaminants_corrected, k = 2, try = 10, trymax = 30, distance = "bray", no.share = TRUE, parallel = detectCores() - 1) # change the k and the try's to optimize the dimensions for better representation when necessary.

NMDS_result <- NMDS

NMDS_result$stress  # take the value of the stress result here. if below 0.05 = excellent. if >0.2 then bad if  between 0.05 to 0.1 = good fit of the NMDS model.
# Make a stress plot for the NMDS.The stressplot evaluates how well the ordination represented the complexity in the data.They help to show you how closely the ordination (y-axis) represents the dissimilarities calculated (x-axis). The points around the red stair steps are the communities, and the distance from the line represents the “stress”, or how they are pulled from their original position to be represented in their ordination.
stressplot(NMDS) # visualize how well the model fits here.

# Customize the plot
par(
  font.axis = 2,   # Bold axis text
  font.lab = 2,    # Bold axis labels
  lwd = 2          # Thickness of the axis lines and border
)

# Re-plot with customizations
stressplot(NMDS)

# Add custom axis labels and lines if necessary
#title(xlab = "Discrepancy in Ordination", ylab = "Stress Value", font.lab = 2)
box(lwd = 2)  # Add thick border around the plot
#Close TIFF device
#dev.off()

plot(NMDS)




####################Plot NMDS######

# Extract NMDS coordinates
nmds_points <- as.data.frame(NMDS$points)
nmds_points$Sample_ID <- rownames(nmds_points)

# Merge with metadata
nmds_plot_df <- merge(nmds_points, metadata, by = "Sample_ID")

# Define colors
group_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hospital_environment" = "#aec7e8",
  "Hospital_sewage" = "#ff7f0e",
  "Public_transport" = "#fdbf6f"
)


# Compute convex hulls per group (optional)
hulls <- nmds_plot_df %>%
  group_by(Group) %>%
  slice(chull(MDS1, MDS2))

# Compute group centroids
centroids <- nmds_plot_df %>%
  group_by(Group) %>%
  summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2))




tiff("Species_NMDS_Ordination.tiff", width = 4.5, height = 4, units = "in", res = 600)   ### Figure 2K

ggplot(nmds_plot_df, aes(x = MDS1, y = MDS2, color = Group)) +
  # Hulls for group clusters
  geom_polygon(data = hulls, aes(fill = Group), alpha = 0.1, show.legend = FALSE, color = NA) +
  
  # Sample points
  geom_point(size = 1, alpha = 0.7) +
  
  # Colors
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  
  # Group labels at centroids
  geom_text(data = centroids, aes(label = Group), color = "black", fontface = "bold", size = 2.5) +
  
  #  Add NMDS stress value annotation
  annotate("text", x = -1, y = -3, label = "Stress = 0.123", color = "RED", size = 3, fontface = "bold", fontface = "italic") +
  
  # Theme and styling
  theme_minimal() +
  labs(title = "", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    panel.grid = element_blank()
  )


dev.off()



##################################################################
### BEFORE RUNNING PERMANOVA CHECK THE HOMOGENEITY OF MULTIVARIATE DISPERSION (the assumption) ##############


# Calculate Bray-Curtis dissimilarity matrix (you’ve done this)
Perm_dist <- vegdist(species_contaminants_corrected, method = "bray", parallel = detectCores() - 1)

# Run betadisper
dispersion <- betadisper(Perm_dist, grouping, type = "centroid")

# Run ANOVA on dispersion
anova_dispersion <- anova(dispersion)

# Convert ANOVA table to a data frame
anova_df <- as.data.frame(anova_dispersion)
anova_df$Term <- rownames(anova_df)
anova_df <- anova_df %>%
  select(Term, everything()) 

anova_df <- anova_df[, c("Term", setdiff(names(anova_df), "Term"))]

# Save ANOVA result
write.csv(anova_df, "MultivariateDispersion_ANOVA.csv", row.names = FALSE)   ##### Supplementary Table 4 (Multivariate_Dispersion_Analysis)

# Extract and save individual sample scores (distances to centroid)
dispersion_scores <- data.frame(Sample = rownames(scores(dispersion, display = "sites")),
                                scores(dispersion, display = "sites"))
dispersion_scores$Group <- grouping

write.csv(dispersion_scores, "MultivariateDispersion_SampleScores.csv", row.names = FALSE) 

# Extract and save centroids
centroids_df <- data.frame(Group = rownames(scores(dispersion, display = "centroids")),
                           scores(dispersion, display = "centroids"))
write.csv(centroids_df, "MultivariateDispersion_Centroids.csv", row.names = FALSE)    ### Supplementary Table 4 (Multivariate_Dispersion_Analysis)




################################################Run permanova if the above condition is satisfied ##################
###################################PERMANOVA ANALYSIS#######################################################################
# PERMANOVA : Assessing differences in community composition is done with permutational Multivariate Analysis of Variance, or perMANOVA.
Env_perm <- adonis2(species_contaminants_corrected ~ Group, data = metadata, distance = "bray", binary = TRUE, parallel = detectCores() - 1) #change this distance to "gower" if the default(bray) doesnt work
Env_perm
write.csv(Env_perm, "Permanova(ADONIS (TABLE 5)).csv", row.names = TRUE)  ####### Supplementary Table 4 (overall_Community_Composition_P)
#################################
# Load vegan if not already


# Custom function for pairwise adonis2
pairwise_adonis2 <- function(data, grouping, method = "bray", binary = TRUE, permutations = 999) {
  groups <- unique(grouping)
  comparisons <- combn(groups, 2, simplify = FALSE)
  results <- data.frame()
  
  for (comp in comparisons) {
    subset_rows <- grouping %in% comp
    dist_sub <- vegdist(data[subset_rows, ], method = method, binary = binary)
    meta_sub <- data.frame(Group = factor(grouping[subset_rows]))
    
    ad <- adonis2(dist_sub ~ Group, data = meta_sub, permutations = permutations)
    
    results <- rbind(results, data.frame(
      Comparison = paste(comp[1], "vs", comp[2]),
      R2 = round(ad$R2[1], 4),
      F_model = round(ad$F[1], 3),
      p_value = ad$`Pr(>F)`[1]
    ))
  }
  
  results$Significance <- cut(results$p_value,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                              labels = c("***", "**", "*", "ns"))
  return(results)
}


# Example run
pairwise_results <- pairwise_adonis2(
  data = species_contaminants_corrected,
  grouping = metadata$Group,
  method = "bray",
  binary = TRUE,
  permutations = 999
)

# View result
print(pairwise_results)

# Save to CSV
write.csv(pairwise_results, "PERMANOVA_Pairwise_Comparisons.csv", row.names = FALSE)  #### Supplementary Table 4 (Pairwise_Comparison_Across_Envi)
#######################################################################################################################################################

# Assuming species_contaminants_corrected is samples x species
# and metadata contains Group info

# Assuming species_contaminants_corrected is samples x species
# and metadata contains Group info



# Add Group info to species matrix
species_df <- as.data.frame(species_contaminants_corrected)
species_df$Group <- metadata$Group

# Convert to long format
species_long <- pivot_longer(species_df, 
                             cols = -Group, 
                             names_to = "Species", 
                             values_to = "Abundance")

# Convert abundance to binary presence/absence
species_long <- species_long %>%
  mutate(Present = ifelse(Abundance > 0, 1, 0))

# Calculate prevalence within each environment
prevalence_by_group <- species_long %>%
  group_by(Group, Species) %>%
  summarise(Prevalence = sum(Present) / n() * 100, .groups = "drop")

# Classify core/secondary/peripheral
prevalence_by_group <- prevalence_by_group %>%
  mutate(Category = case_when(
    Prevalence >= 80 ~ "Core",
    Prevalence >= 50 & Prevalence < 80 ~ "Secondary",
    TRUE ~ "Peripheral"
  ))

# Save the output
write.csv(prevalence_by_group, "Core_Pathogens_By_Environment.csv", row.names = FALSE)




# Count number of species in each category per group
core_counts <- prevalence_by_group %>%
  group_by(Group, Category) %>%
  summarise(n = n(), .groups = "drop")

# Define color palette
category_colors <- c(
  "Core" = "#e377c2",
  "Secondary" = "#969696",
  "Peripheral" = "#6baed6"
)

tiff("Prevalence_corePlot.tiff", width = 8, height = 6, units = "in", res = 600)    ### Supplementary 1G
# Updated plot with flipped coordinates
core_plot <- ggplot(core_counts, aes(x = Group, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = category_colors) +
  labs(title = "",
       x = "", y = "Number of Species") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  # Horizontal axes now
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.x = element_text(size = 12),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2)
  ) +
  coord_flip()  # 💥 Flip coordinates


# Print the plot
print(core_plot)

dev.off()
###############################################
######Calculating the species prevalence

### We need to make the Species in a row to really get their numbers for our Categorization
species_contaminants_corrected1 = t(species_contaminants_corrected)

# Assuming presence_matrix has already been defined
presence_matrix <- species_contaminants_corrected1  # Adjust this to your actual data
presence_matrix[presence_matrix > 0] <- 1  # Convert abundance data to presence (1) or absence (0)

# Get the number of samples and species
n <- nrow(presence_matrix)  # Number of samples
m <- ncol(presence_matrix)  # Number of species

# Calculate the prevalence of each species across all samples
x <- rep(0, n)
for (i in 1:n) {
  x[i] <- sum(presence_matrix[i, ] == 1)
}
y <- (x / m) * 100  # Convert to percentage

# Categorize species based on prevalence
Core_species <- row.names(presence_matrix[which(y > 80),])  # Prevalence greater than 80%.
Secondary_species <- row.names(presence_matrix[which(y > 50 & y <= 80),])  # Prevalence between 50% and 80%.
Peripheral_species <- row.names(presence_matrix[which(y <= 50),])  # Prevalence less than 50


# Create a dataframe for visualization
species_categories <- data.frame(
  Species = row.names(presence_matrix),
  Prevalence = y,
  Category = ifelse(y > 80, "Core", 
                    ifelse(y > 50 & y <= 80, "Secondary", "Peripheral"))
)



# Define refined color palette
prevalence_colors <- c(
  "Core" = "#e377c2",         # Muted violet-pink
  "Secondary" = "#969696",    # Elegant mid-gray
  "Peripheral" = "#6baed6"    # Desaturated steel blue
)

# Save the updated plot
tiff("Species_prevalence_mature.tiff", width = 5, height = 4, units = "in", res = 600)   ###### Figure 3B

P1 <- ggplot(species_categories, aes(x = Prevalence, fill = Category)) +
  geom_density(alpha = 0.25, color = "black", adjust = 1.5) +
  geom_rug(aes(color = Category), sides = "b", alpha = 0.4, length = unit(0.02, "npc")) +
  
  # Apply refined colors
  scale_fill_manual(values = prevalence_colors) +
  scale_color_manual(values = prevalence_colors) +
  
  # Thresholds
  geom_vline(xintercept = c(50, 80), linetype = "dashed", color = "red", size = 0.6) +
  
  # Annotations
  annotate("text", x = 94, y = 0.02, label = paste0("Core (>80%)\n(", length(Core_species), " species)"),
           size =1.8, fontface = "bold", hjust = 1) +
  annotate("text", x = 68, y = 0.02, label = paste0("Secondary (50–80%)\n(", length(Secondary_species), " species)"),
           size = 1.8, fontface = "bold", hjust = 0.5) +
  annotate("text", x = 16, y = 0.017, label = paste0("Peripheral (≤50%)\n(", length(Peripheral_species), " species)"),
           size = 1.8, fontface = "bold", hjust = 0) +
  
  # Labels and theme
  labs(
    x = "Species Prevalence (%)",
    y = "Density"
  ) +
  coord_cartesian(ylim = c(0, 0.18)) +  # Optional: remove excess white space
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    plot.margin = margin(10, 10, 10, 10)
  )

# Print and save
print(P1)
dev.off()


######################## Upset Plot
# Convert species data to long format
species_long <- species_contaminants_corrected %>%
  as.data.frame() %>%
  mutate(Sample_ID = rownames(.)) %>%
  left_join(metadata[, c("Sample_ID", "Group")], by = "Sample_ID") %>%
  pivot_longer(-c(Sample_ID, Group), names_to = "Species", values_to = "Abundance") %>%
  mutate(Present = ifelse(Abundance > 0, 1, 0)) %>%
  filter(Present == 1)


# Transform long-format back to species x group matrix
species_env_matrix <- species_long %>%
  filter(Present == 1) %>%
  distinct(Group, Species) %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = Group, values_from = Present, values_fill = 0)

# Preview
head(species_env_matrix)


# If needed, convert tibble to data.frame
species_env_matrix <- as.data.frame(species_env_matrix)


# Make sure species names are the rownames
rownames(species_env_matrix) <- species_env_matrix$Species
species_env_matrix$Species <- NULL


group_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hospital_environment" = "#aec7e8",
  "Hospital_sewage" = "#ff7f0e",
  "Public_transport" = "#fdbf6f"
)


# Adjust margins and line width for high-res plotting
par(lwd = 2, mar = c(5, 5, 2, 2))

# Your set of columns representing presence (1) / absence (0) of species in each environment
set_cols <- c("Ambulance", "Hospital_environment", "Hospital_sewage", "Public_transport")
bar_colors <- group_colors[set_cols]

# Create UpSet plot
tiff("UpSet_Plot(Figure_11).tiff", width = 10, height = 6, units = "in", res = 600)  #### Figure 3A
UpSetR::upset(
  species_env_matrix,
  sets = set_cols,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = bar_colors,      # Color per environment
  main.bar.color = "black",         # Black bars for intersections
  text.scale = 2,
  point.size = 5,
  line.size = 1
)
dev.off()




#######################################################Top 30 Pathogenic Species per Environment

# Limit to Core species only
top_core_species <- prevalence_by_group %>%
  filter(Category == "Core") %>%
  group_by(Group) %>%
  arrange(desc(Prevalence)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Define consistent fill colors per environment
group_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hospital_environment" = "#aec7e8",
  "Hospital_sewage" = "#ff7f0e",
  "Public_transport" = "#fdbf6f"
)

tiff("Core_species_per_Env.tiff", width = 8, height = 6, units = "in", res = 600)   ##### Supplementary Figure 3
# Plot: Top 20 Core Pathogens per Environment
ggplot(top_core_species, aes(x = reorder(Species, Prevalence), y = Prevalence, fill = Group)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(title = "",
       x = "Species", y = "Prevalence (%)") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()
################################################### Top 20 Gneral Core species Heatmap

# Get total prevalence of each core species across environments
top_core_20 <- prevalence_by_group %>%
  filter(Category == "Core") %>%
  group_by(Species) %>%
  summarise(TotalPrevalence = sum(Prevalence, na.rm = TRUE)) %>%
  arrange(desc(TotalPrevalence)) %>%
  slice_head(n = 20)

# Filter the main prevalence table to just top 20
top20_matrix <- prevalence_by_group %>%
  filter(Species %in% top_core_20$Species) %>%
  select(Group, Species, Prevalence) %>%
  pivot_wider(names_from = Group, values_from = Prevalence, values_fill = 0) %>%
  column_to_rownames("Species")


# Create improved color palette (e.g., blue gradient)
heat_colors <- colorRampPalette(brewer.pal(n = 8, name = "Blues"))(100)

# Plot heatmap with updated theme
tiff("Top20_Core_Heatmap(Figure13).tiff", width = 6, height = 8, units = "in", res = 600)   ##### Figure 3C

pheatmap(top20_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,  # Don't cluster environments to keep logical order
         color = heat_colors,
         fontsize_row = 9,
         fontsize_col = 11,
         main = "",
         angle_col = 45,
         border_color = NA,
         na_col = "white",
         display_numbers = FALSE,
         legend = TRUE,
         cellwidth = 25,  # adjust if too wide/narrow
         cellheight = 12  # adjust for label legibility
)

dev.off()

###################################################
####Indicator Species Analysis

# Convert abundance data to binary (presence/absence)
species_binary <- species_contaminants_corrected
species_binary[species_binary > 0] <- 1


# Check alignment
all(rownames(species_binary) == metadata$Sample_ID)  # Should be TRUE

# If not, re-order metadata
metadata <- metadata[match(rownames(species_binary), metadata$Sample_ID), ]

# Check for NAs in metadata group
sum(is.na(metadata$Group))  # should be 0

# Check for NAs in species matrix
sum(is.na(species_binary))  # should also be 0



# Filter out rows where Group is NA
valid_rows <- !is.na(metadata$Group)

# Apply filter to both metadata and species matrix
metadata_clean <- metadata[valid_rows, ]
species_binary_clean <- species_binary[valid_rows, ]


set.seed(123)

indval_result <- multipatt(species_binary_clean, metadata_clean$Group,
                           func = "IndVal.g", control = how(nperm = 999))

# View summary of significant results
summary(indval_result, alpha = 0.05)

# Extract the indicator statistics
indval_df <- as.data.frame(indval_result$sign)
indval_df$Species <- rownames(indval_df)

# Add p-values directly from the result object
indval_df$p.value <- indval_result$p.value  # <-- this is where the p-values live

names(indval_result)
# View column names inside the sign object
colnames(indval_result$sign)



# Extract indicator species results
indval_df <- as.data.frame(indval_result$sign)
indval_df$Species <- rownames(indval_df)

significant_indicators <- indval_df %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)



significant_indicators <- significant_indicators %>%
  rowwise() %>%
  mutate(Environment = names(.)[which.max(c_across(starts_with("s.")))]) %>%
  ungroup() %>%
  mutate(Environment = gsub("s\\.", "", Environment))  # Clean up label

write.csv(significant_indicators, "Indicator_Species_Significant.csv", row.names = FALSE)    #### Supplementary Table 6
head(significant_indicators[, c("Species", "Environment", "stat", "p.value")])


# Load and prep
indicators <- read.csv("Indicator_Species_Significant.csv")

top_indicators <- indicators %>%
  arrange(desc(stat)) %>%
  slice_head(n = 60)

# Use your environment color palette
env_colors <- c(
  "Ambulance" = "#1f77b4",
  "Hospital_environment" = "#aec7e8",
  "Hospital_sewage" = "#ff7f0e",
  "Public_transport" = "#fdbf6f"
)


# Plot Indicator Species with updated theme
tiff("Top_indicators.tiff", width = 6, height = 8, units = "in", res = 600)   #### Figure 3F 

# Plot
ggplot(top_indicators, aes(x = stat, y = reorder(Species, stat), color = Environment)) +
  geom_point(size = 3.5) +
  labs(
    title = "",
    x = "Indicator Value",
    y = ""
  ) +
  scale_color_manual(values = env_colors) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    
    # Legend styling — inside figure
    legend.title = element_blank(),
    legend.position = c(0.6, 0.02),  #  Adjust this for best position
    legend.justification = c("left", "bottom"),
    #legend.box.background = element_rect(color = "black", size = 0.3),
    #legend.background = element_rect(fill = "white", color = NA),
    #  NEW legend sizing tweaks
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()


##########################################nETWORK aNALYSIS For Priority Pathogens ######

# Define exact species
exact_matches <- c(
  "Escherichia coli", "Acinetobacter baumannii", "Mycobacterium tuberculosis", "Salmonella Typhii",
  "Enterococcus faecium", "Pseudomonas aeruginosa", "Neisseria gonorrhoeae",
  "Staphylococcus aureus", "Streptococcus pneumoniae",
  "Haemophilus influenza", "Helicobacter pylori", "Clostridioides difficile",
  "Klebsiella oxytoca", "Enterococcus faecalis", "Enterococcus avium",
  "Candida albicans", "Candida parapsilosis", "Candida glabrata", "Candida dubliniensis",
  "Proteus mirabilis", "Klebsiella pneumoniae", "Stenotrophomonas maltophilia"
)

# Define genus-level group matches
group_matches <- list(
  "Shigella" = "Shigella",
  "Enterobacter" = "Enterobacter",
  "Citrobacter" = "Citrobacter",
  "Proteus" = "Proteus",
  "Serratia" = "Serratia",
  "Streptococcus" = "Streptococcus",
  "Morganella" = "Morganella",
  "Providencia" = "Providencia",
  "Campylobacter" = "Campylobacter"
)

# Combine both: species or genus match
all_priority_species <- colnames(species_contaminants_corrected)[
  colnames(species_contaminants_corrected) %in% exact_matches |
    sapply(colnames(species_contaminants_corrected), function(x) {
      any(sapply(group_matches, function(pattern) grepl(pattern, x, ignore.case = TRUE)))
    })
]

# Subset matrix for WHO + co-occurring species
library(Matrix)
library(vegan)


# Optional: remove rare species
filtered_matrix <- species_contaminants_corrected[, colSums(species_contaminants_corrected) > 4]

# Use Jaccard or Spearman for binary data
cor_matrix <- cor(filtered_matrix, method = "spearman")  # for abundance
# OR
# jaccard_dist <- vegdist(t(filtered_matrix), method = "jaccard")


library(igraph)
set.seed(43)

# Threshold to include edges
threshold <- 0.8

# Threshold to include edges and convert to absolute
adj_matrix <- ifelse(abs(cor_matrix) > threshold, abs(cor_matrix), 0)
diag(adj_matrix) <- 0

# Build igraph object
network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Annotate priority species
V(network)$type <- ifelse(V(network)$name %in% all_priority_species, "WHO Priority", "Other")

# Set visual attributes
V(network)$color <- ifelse(V(network)$type == "WHO Priority", "#d62728", "#1f77b4")
V(network)$size <- ifelse(V(network)$type == "WHO Priority", 6, 3)

plot(
  network,
  vertex.label = NA,
  edge.color = "gray80",
  edge.width = E(network)$weight * 2,
  layout = layout_with_fr,
  main = ""
)

# Add labels for WHO priority pathogens only
V(network)$label <- ifelse(V(network)$name %in% all_priority_species, V(network)$name, NA)

# Optional: Customize label size and font
V(network)$label.cex <- 0.8
V(network)$label.color <- "black"
V(network)$label.font <- 2  # Bold

# Save plot to TIFF
tiff("WHO_Priority_Network.tiff", width = 8, height = 8, units = "in", res = 600)    #### Figure 3G

# Re-plot with labels
plot(
  network,
  vertex.label = V(network)$label,
  edge.color = "gray80",
  edge.width = E(network)$weight * 4,
  vertex.size = V(network)$size,
  vertex.color = V(network)$color,
  layout = layout_with_fr,
  main = ""
)

dev.off()
######################################################## End ##########################
