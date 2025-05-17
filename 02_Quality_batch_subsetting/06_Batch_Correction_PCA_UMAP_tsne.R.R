# -----------------------------------------------
# Description: Batch correction, PCA, UMAP, t-SNE, PERMANOVA on Hellinger-transformed species matrix
# -----------------------------------------------
#############################################
# Load libraries
library(vegan)
library(ggplot2)
library(car)
library(gridExtra)
library(lme4)
library(MuMIn)
library(lattice)
library(scales)
library(limma)
library(umap)
library(Rtsne)
library(shadowtext)
library(dplyr)
##############################################################

# Set seed for reproducibility
set.seed(45)



###############################################################
# Load data files
exp_matrix <- read.csv("/home/rachid_lab/Downloads/Species_expression_matrix.csv", row.names = 1)
metadata <- read.csv("/home/rachid_lab/Downloads/Species_metadata.csv")

# Data Normalisation using vegan
counts <- decostand(exp_matrix, method =  "hellinger")


# (Optional) Write to a new CSV file
write.csv(counts, "Species_Helinger_transformed.csv", row.names = TRUE)

# Ensure the columns of counts match the Sample_IDs in metadata
counts <- counts[, metadata$Sample_ID]

###########################################################################################

#####Inspecting the original Data

# Perform PCA (Before)

pca <- prcomp(t(counts), scale. = TRUE)

# Combine PCA results with metadata for plotting

pca_data <- as.data.frame(pca$x)
pca_data$Project_ID <- metadata$Project_ID
pca_data$Country <- metadata$Country
pca_data$Region <- metadata$Region
pca_data$Continent <- metadata$Continent  
pca_data$Instrument <- metadata$Instrument
pca_data$Sequencing_Center <- metadata$Sequencing_Center
pca_data$Group <- metadata$Group
pca_data$environmental_material <- metadata$environmental_material


# Create PCA plots

pca_plot_project_ID <- xyplot(PC2 ~ PC1, data = pca_data, groups = Project_ID,
                              auto.key = FALSE, main = "Project", xlab = "PC1", ylab = "PC2",
                              par.settings = list(axis.line = list(col = "black")))

pca_plot_country <- xyplot(PC2 ~ PC1, data = pca_data, groups = Country,
                           auto.key = FALSE, main = "Country", xlab = "PC1", ylab = "PC2",
                           par.settings = list(axis.line = list(col = "black")))

pca_plot_region <- xyplot(PC2 ~ PC1, data = pca_data, groups = Region,
                          auto.key = FALSE, main = "Region", xlab = "PC1", ylab = "PC2",
                          par.settings = list(axis.line = list(col = "black")))

pca_plot_continent <- xyplot(PC2 ~ PC1, data = pca_data, groups = Continent,
                             auto.key = FALSE, main = "Continent", xlab = "PC1", ylab = "PC2",
                             par.settings = list(axis.line = list(col = "black")))

pca_plot_instrument <- xyplot(PC2 ~ PC1, data = pca_data, groups = Instrument,
                              auto.key = FALSE, main = "Instrument", xlab = "PC1", ylab = "PC2",
                              par.settings = list(axis.line = list(col = "black")))

pca_plot_sequencing_center <- xyplot(PC2 ~ PC1, data = pca_data, groups = Sequencing_Center,
                                     auto.key = FALSE, main = "Sequencing Center", xlab = "PC1", ylab = "PC2",
                                     par.settings = list(axis.line = list(col = "black")))

pca_plot_group <- xyplot(PC2 ~ PC1, data = pca_data, groups = Group,
                         auto.key = FALSE, main = "Group", xlab = "PC1", ylab = "PC2",
                         par.settings = list(axis.line = list(col = "black")))

pca_plot_environmental_material <- xyplot(PC2 ~ PC1, data = pca_data, groups = environmental_material,
                                          auto.key = FALSE, main = "Environmental Material", xlab = "PC1", ylab = "PC2",
                                          par.settings = list(axis.line = list(col = "black")))


# Helper function to create plots

create_pca_plot <- function(data, group_var, title) {
  xyplot(PC2 ~ PC1, data = data, groups = group_var,
         auto.key = FALSE,
         xlab = list("PC1", cex = 0.6, font = 2),  # Reduce x-axis title font size
         ylab = list("PC2", cex = 0.6, font = 2),  # Reduce y-axis title font size
         scales = list(
           x = list(rot = 45, cex = 0.6, font = 2),  # Rotate x-axis labels and reduce font size
           y = list(cex = 0.6, font = 2)  # Reduce y-axis labels font size
         ),
         par.settings = list(axis.line = list(col = "black")),
         axis.components = list(
           top = list(tck = 0.5),  # No tick marks on top axis
           right = list(tck = 0.5),  # No tick marks on right axis
           left = list(tck = 0.5),  # Maintain tick marks on left axis
           bottom = list(tck = 1)  # Maintain tick marks on bottom axis
         ),
         xlab.top = list(title, cex = 0.6, col = "black", border = TRUE, font = 2),  # Reduce font size and make not bold
         panel = function(...) {
           panel.grid(h = -1, v = -1, lty = 2, col = alpha("black", 0.3))  # Dashed grid lines with increased alpha
           panel.superpose(...)
         })
}

# Save the plot as a tiff file

tiff("Technical_Variation_PCA.tiff", width = 8, height = 6, units = "in", res = 600)   ###### Supplementary Figure 2


# Create PCA plots

pca_plot_project_ID <- create_pca_plot(pca_data, pca_data$Project_ID, "Project")
pca_plot_country <- create_pca_plot(pca_data, pca_data$Country, "Country")
pca_plot_region <- create_pca_plot(pca_data, pca_data$Region, "Region")
pca_plot_continent <- create_pca_plot(pca_data, pca_data$Continent, "Continent")
pca_plot_instrument <- create_pca_plot(pca_data, pca_data$Instrument, "Instrument")
pca_plot_sequencing_center <- create_pca_plot(pca_data, pca_data$Sequencing_Center, "Sequencing Center")
pca_plot_group <- create_pca_plot(pca_data, pca_data$Group, "Group")
pca_plot_environmental_material <- create_pca_plot(pca_data, pca_data$environmental_material, "Environmental Material")

# Arrange plots in a grid

print(pca_plot_group, split = c(1, 1, 4, 2), more = TRUE)
print(pca_plot_environmental_material, split = c(2, 1, 4, 2), more = TRUE)
print(pca_plot_project_ID, split = c(3, 1, 4, 2), more = TRUE)
print(pca_plot_continent, split = c(4, 1, 4, 2), more = TRUE)
print(pca_plot_country, split = c(1, 2, 4, 2), more = TRUE)
print(pca_plot_region, split = c(2, 2, 4, 2), more = TRUE)
print(pca_plot_sequencing_center, split = c(3, 2, 4, 2), more = TRUE)
print(pca_plot_instrument, split = c(4, 2, 4, 2), more = FALSE)


#Close
dev.off()

##############################################################################################
#Perform PERMANOVA (adonis) for to understand batch contribution before correction

# Prepare counts data for PERMANOVA
counts_t <- t(counts)  # Transpose counts matrix for proper alignment

# Load necessary library
library(vegan)

# Run PERMANOVA for different factors
adonis_result_group <- adonis2(counts_t ~ metadata$Group, data = metadata, method = "euclidean")
adonis_result_environmental_material <- adonis2(counts_t ~ metadata$environmental_material, data = metadata, method = "euclidean")
adonis_result_Project_ID <- adonis2(counts_t ~ metadata$Project_ID, data = metadata, method = "euclidean")
adonis_result_Country <- adonis2(counts_t ~ metadata$Country, data = metadata, method = "euclidean")
adonis_result_Region <- adonis2(counts_t ~ metadata$Region, data = metadata, method = "euclidean")
adonis_result_Continent <- adonis2(counts_t ~ metadata$Continent, data = metadata, method = "euclidean")
adonis_result_Instrument <- adonis2(counts_t ~ metadata$Instrument, data = metadata, method = "euclidean")
adonis_result_Sequencing_Center <- adonis2(counts_t ~ metadata$Sequencing_Center, data = metadata, method = "euclidean")

# Convert the PERMANOVA results to data frames
adonis_df_group <- as.data.frame(adonis_result_group)
adonis_df_environmental_material <- as.data.frame(adonis_result_environmental_material)
adonis_df_Project_ID <- as.data.frame(adonis_result_Project_ID)
adonis_df_Country <- as.data.frame(adonis_result_Country)
adonis_df_Region <- as.data.frame(adonis_result_Region)
adonis_df_Continent <- as.data.frame(adonis_result_Continent)
adonis_df_Instrument <- as.data.frame(adonis_result_Instrument)
adonis_df_Sequencing_Center <- as.data.frame(adonis_result_Sequencing_Center)

# Combine all PERMANOVA results into a single data frame
combined_adonis_results <- rbind(
  cbind(Factor = "Group", adonis_df_group),
  cbind(Factor = "Environmental_Material", adonis_df_environmental_material),
  cbind(Factor = "Project_ID", adonis_df_Project_ID),
  cbind(Factor = "Country", adonis_df_Country),
  cbind(Factor = "Region", adonis_df_Region),
  cbind(Factor = "Continent", adonis_df_Continent),
  cbind(Factor = "Instrument", adonis_df_Instrument),
  cbind(Factor = "Sequencing_Center", adonis_df_Sequencing_Center)
)

# Save the combined PERMANOVA results to a CSV file
write.csv(combined_adonis_results, file = "Technical_variation_Adonis.csv", row.names = FALSE)   ##### Supplementary Table 3 (Contribution of metadata variables)




#################################################
# Remove Batch Effects using Limma
counts_matrix <- as.matrix(counts)

# Create the design matrix for biological covariates (e.g., Group)
design <- model.matrix(~ metadata$Group)

# Create a covariate matrix excluding Country and Region
covariates <- model.matrix(~ metadata$Instrument + metadata$Sequencing_Center)[, -1]


# Apply removeBatchEffect to remove the impact of multiple batch effects
counts_corrected <- removeBatchEffect(
  counts_matrix,
  batch = metadata$Project_ID,
  covariates = covariates,
  design = design
)

# Save the corrected counts
write.csv(counts_corrected, "Species_batch_corrected.csv", row.names = TRUE)

# PCA Before and After Batch Correction
pca_before <- prcomp(t(counts_matrix), scale. = TRUE)
pca_after <- prcomp(t(counts_corrected), scale. = TRUE)

# PCA before correction
pca_data_before <- as.data.frame(pca_before$x)
pca_data_before$Group <- metadata$Group
p_before <- ggplot(pca_data_before, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA Before Batch Correction")

# PCA after correction
pca_data_after <- as.data.frame(pca_after$x)
pca_data_after$Group <- metadata$Group
p_after <- ggplot(pca_data_after, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA After Batch Correction")

# Arrange the PCA plots side by side
grid.arrange(p_before, p_after, ncol = 2)

######################################################
#Perform UMAP (Before and After Batch Correction)

# UMAP for data before batch correction
umap_before <- umap(t(counts_matrix), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data_before <- as.data.frame(umap_before$layout)
umap_data_before$Group <- metadata$Group
colnames(umap_data_before) <- c("UMAP1", "UMAP2", "Group")

# UMAP for data after batch correction
umap_after <- umap(t(counts_corrected), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data_after <- as.data.frame(umap_after$layout)
umap_data_after$Group <- metadata$Group
colnames(umap_data_after) <- c("UMAP1", "UMAP2", "Group")

# Plot UMAP results before and after batch correction

# UMAP plot before correction
umap_plot_before <- ggplot(umap_data_before, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "UMAP Before Batch Correction")

# UMAP plot after correction
umap_plot_after <- ggplot(umap_data_after, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "UMAP After Batch Correction")

# Arrange the UMAP plots side by side
grid.arrange(umap_plot_before, umap_plot_after, ncol = 2)
############################################################### 
##### Perform permanova to asses the biological variability before and after 

# Load necessary library for PERMANOVA
library(vegan)

# Run PERMANOVA on UMAP coordinates before batch correction
permanova_result_before <- adonis2(umap_data_before[, 1:2] ~ Group, data = umap_data_before, method = "euclidean")
print("PERMANOVA result for UMAP coordinates before batch correction:")
print(permanova_result_before)

# Convert the PERMANOVA result to a data frame for exporting
permanova_result_before_df <- as.data.frame(permanova_result_before)

# Add an identifier column to indicate 'Before Batch Correction'
permanova_result_before_df$Correction_Status <- "Before Batch Correction"

# Run PERMANOVA on UMAP coordinates after batch correction
permanova_result_after <- adonis2(umap_data_after[, 1:2] ~ Group, data = umap_data_after, method = "euclidean")
print("PERMANOVA result for UMAP coordinates after batch correction:")
print(permanova_result_after)

# Convert the PERMANOVA result to a data frame for exporting
permanova_result_after_df <- as.data.frame(permanova_result_after)

# Add an identifier column to indicate 'After Batch Correction'
permanova_result_after_df$Correction_Status <- "After Batch Correction"

# Combine both data frames
combined_permanova_results <- rbind(permanova_result_before_df, permanova_result_after_df)

# Save the combined PERMANOVA results to a single CSV file
write.csv(combined_permanova_results, file = "UMAP_Adonis_Biological_variation_corrected.csv", row.names = TRUE)  ### Supplementary Table 3 (Post_Correction PERMANOVA Analysis)

# Print a message indicating successful save
print("Combined PERMANOVA results have been saved to 'PERMANOVA_combined_results.csv'")


################## Preferably use this to visualize the before and after correction ################
set.seed(43)

# Remove duplicates and Run t-SNE on the transposed data matrix (Before Correction)
counts_matrix_unique <- t(counts_matrix)  # Transpose the matrix for t-SNE
counts_matrix_unique <- counts_matrix_unique[!duplicated(counts_matrix_unique), ]  # Remove duplicate rows

tsne_result_before <- Rtsne(counts_matrix_unique, dims = 2, perplexity = 50)

# Convert t-SNE result to a data frame for plotting
tsne_data_before <- as.data.frame(tsne_result_before$Y)
colnames(tsne_data_before) <- c("tSNE1", "tSNE2")  # Set the correct column names
group_labels_before <- metadata$Group[!duplicated(t(counts_matrix))]  # Match Group labels to unique rows
tsne_data_before$Group <- group_labels_before

# Remove duplicates and Run t-SNE on the transposed data matrix (After Correction)
counts_corrected_unique <- t(counts_corrected)  # Transpose the matrix for t-SNE
counts_corrected_unique <- counts_corrected_unique[!duplicated(counts_corrected_unique), ]  # Remove duplicate rows

tsne_result_after <- Rtsne(counts_corrected_unique, dims = 2, perplexity = 50)

# Convert t-SNE result to a data frame for plotting
tsne_data_after <- as.data.frame(tsne_result_after$Y)
colnames(tsne_data_after) <- c("tSNE1", "tSNE2")  # Set the correct column names
group_labels_after <- metadata$Group[!duplicated(t(counts_corrected))]  # Match Group labels to unique rows
tsne_data_after$Group <- group_labels_after

# Plot t-SNE result before and after batch correction

# Save the plot as a tiff file
tiff("Biological_variability_TSNE.tiff", width = 8, height = 4, units = "in", res = 600)

# t-SNE plot before correction
tsne_plot_before <- ggplot(tsne_data_before, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 0.8, stroke = 1, shape = 20) +
  theme_minimal() +
  labs(title = "t-SNE Before Batch Correction") + 
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_blank(),  # Remove individual axis lines
    axis.ticks = element_line(color = "black", size = 0.4),  # Add tick marks
    axis.ticks.length = unit(0.2, "cm"),  # Length of the tick marks
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the entire plot
    plot.title = element_text(hjust = 0.5, size = 12),  # Center the plot title
    #legend.position = "none"  # Remove the legend from this plot
  )

# t-SNE plot after correction
tsne_plot_after <- ggplot(tsne_data_after, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 0.8, stroke = 1, shape = 20) +
  theme_minimal() +
  labs(title = "t-SNE After Batch Correction") +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_blank(),  # Remove individual axis lines
    axis.ticks = element_line(color = "black", size = 0.4),  # Add tick marks
    axis.ticks.length = unit(0.2, "cm"),  # Length of the tick marks
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border around the entire plot
    plot.title = element_text(hjust = 0.5, size = 12),  # Center the plot title
    #legend.position = "none"  # Add the legend to the right side of the plot
  )

# Arrange the t-SNE plots side by side
grid.arrange(tsne_plot_before, tsne_plot_after, ncol = 2)

dev.off()


#####################################################################################################
#### Figures used in Manuscript




# Define group colors
# Color palette 
group_colors <- c(
  "Ambulance" = "#1f77b4",            
  "Hospital_environment" = "#aec7e8", 
  "Hospital_sewage" = "#ff7f0e",      
  "Public_transport" = "#fdbf6f"      
)


###################
#Before Correction
#################
# t-SNE plot before correction
tsne_plot_before <- ggplot(tsne_data_before, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 2, shape = 16, alpha = 1.6, stroke = 0) +  # ← Softer, rounder points
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Hides the legend completely
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.ticks.length = unit(0.1, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    plot.title = element_text(hjust = 0.5, size = 6)
  )



# Save the plot
ggsave("tSNE_Before_Batch_Correction.tiff",
       plot = tsne_plot_before,
       width = 4, height = 4, dpi = 600,
       bg = "white")  # ← Ensures white background

# Show the plot

print(tsne_plot_before)  #### Figure 2C (Before correction)


###################
#After Correction
#################

# Optional: define group colors
# New softened color palette from the violin plot
group_colors <- c(
  "Ambulance" = "#1f77b4",            # Soft deep blue
  "Hospital_environment" = "#aec7e8", # Light blue
  "Hospital_sewage" = "#ff7f0e",      # Orange
  "Public_transport" = "#fdbf6f"      # Light orange/yellow
)

# --- Compute Centroids ---
centroids <- tsne_data_after %>%
  group_by(Group) %>%
  summarize(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2))

# --- Compute Convex Hulls for each group ---
hulls <- tsne_data_after %>%
  group_by(Group) %>%
  slice(chull(tSNE1, tSNE2))  # get boundary points

# --- Create Plot with Points + Hulls + Centroids ---
tsne_plot_with_hulls <- ggplot(tsne_data_after, aes(x = tSNE1, y = tSNE2, color = Group)) +
  # Convex hull polygons
  geom_polygon(data = hulls, aes(fill = Group), alpha = 0.1, color = NA, show.legend = FALSE) +
  # Main points
  geom_point(size = 2, shape = 16, alpha = 1.6, stroke = 0) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  # Centroid text labels
  geom_shadowtext(data = centroids, aes(label = Group), 
                  size = 2.4, fontface = "bold", color = "black",
                  bg.color = "white", bg.r = 0.4) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

ggsave("tSNE_After_Bach_Correction.tiff",
       plot = tsne_plot_with_hulls,
       width = 4.3, height = 4, dpi = 600, bg = "white")

# Show the plot
print(tsne_plot_with_hulls)      ##### Figure 2c (After correction)


#########################################################
# -----------------------------------------------
# Final Cleanup
# -----------------------------------------------

species_corrected <- read.csv("Species_batch_corrected.csv", row.names = 1)
species_corrected[species_corrected < 0] <- 0
write.csv(species_corrected, "species_final.csv", row.names = TRUE)

cat("Batch correction and dimensionality analysis complete.\n")
#######################################################################################################################################################