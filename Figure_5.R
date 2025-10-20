# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)

Bacteria_ps_filtered <- readRDS("Bacteria_ps_filtered_clean.rds")
# 1. **Normalize Data**
#Cumulative Sum Scaling (CSS):
# CSS: Adjusts for varying sequencing depths and compositional bias.
# Another option that normalizes count data while accounting for the overabundance of a few dominant taxa, mitigating compositional bias.
# Adjusts for differences in library sizes and overrepresentation of high-abundance taxa.

library(metagenomeSeq)
# Convert phyloseq to metagenomeSeq format
physeq_metagenomeSeq <- phyloseq_to_metagenomeSeq(Bacteria_ps_filtered)
# Apply CSS normalization
pCSS <- cumNorm(physeq_metagenomeSeq)
# Convert back to phyloseq object
otu_table_css <- MRcounts(pCSS, norm=TRUE)
physeq_css <- phyloseq(otu_table(otu_table_css, taxa_are_rows=TRUE), sample_data(Bacteria_ps_filtered), tax_table(Bacteria_ps_filtered), phy_tree(Bacteria_ps_filtered))



# Subset data for pH 6.2
subset_physeq_6.2 <- subset_samples(physeq_css, pH == 6.2)

# Compute Bray-Curtis distance for the subset
bray_dist_subset_6.2 <- distance(subset_physeq_6.2, method = "bray")

# Perform PCoA
pcoa_bray_subset_6.2 <- ordinate(subset_physeq_6.2, method = "PCoA", distance = "bray")

# Variance explained by axes
variance_explained_subset_6.2 <- pcoa_bray_subset_6.2$values$Relative_eig

# Convert sample data to data frame for the subset
metadata_subset_6.2 <- as.data.frame(as.matrix(sample_data(subset_physeq_6.2)))
metadata_subset_6.2$pH <- as.factor(metadata_subset_6.2$pH)
metadata_subset_6.2$Species <- as.factor(metadata_subset_6.2$Species)

# PERMANOVA analysis for the subset
adonis_results_subset_6.2 <- adonis2(bray_dist_subset_6.2 ~ Species, data = metadata_subset_6.2, permutations = 999)

# Perform pairwise PERMANOVA
pairwise_results_subset_6.2 <- pairwise.adonis(
  bray_dist_subset_6.2,  # Your Bray-Curtis distance matrix
  factors = metadata_subset_6.2$Species,  # Species as factors
  perm = 999  # Number of permutations
)

# Save pairwise PERMANOVA results as CSV
write.csv(pairwise_results_subset_6.2, "pairwise_results_pH_6.2.csv")

# Create legend text for PERMANOVA results
#legend_text_subset <- paste0("PERMANOVA: F = ", round(adonis_results_subset_6.2$F[1], 2),
    #                         ", R² = ", round(adonis_results_subset_6.2$R2[1], 3),
 #                            ", p = ", format.pval(adonis_results_subset_6.2$Pr[1], digits = 3))
legend_text_subset_6_2 <- paste0("F = ", round(adonis_results_subset_6.2$F[1], 2),
                             ", R² = ", round(adonis_results_subset_6.2$R2[1], 3),
                             ", p = ", format.pval(adonis_results_subset_6.2$Pr[1], digits = 3))


# Prepare ordination data for plotting
ordination_data_subset_6.2 <- plot_ordination(subset_physeq_6.2, pcoa_bray_subset_6.2, color = "pH", shape = "Species", justDF = TRUE)

# Define color and shape mappings
color_map <- c("Empty" = "black", "Tomato" = "brown", "Petunia" = "darkgreen", "Geranium" = "goldenrod", "Marigold" = "cyan3")
shape_mapping <- c("Empty" = 21, "Geranium" = 16, "Marigold" = 7, "Petunia" = 3, "Tomato" = 17)

# Plot for pH 6.2 with PERMANOVA results, custom color mapping, and shape customization
plot_6_2 <- ggplot(ordination_data_subset_6.2) +
  # Plot points with customized stroke based on species
  geom_jitter(
    aes(x = Axis.1, y = Axis.2, color = Species, shape = Species),
    size = 7, width = 0.05, height = 0.05,
    stroke = ifelse(ordination_data_subset_6.2$Species == "Empty" | ordination_data_subset_6.2$Species == "Petunia", 3, 1)  # Thicker border for "Empty" and "Petunia"
  ) +
  # Ellipses for each species group
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = Species), type = "norm", linetype = 2, size = 1) +
  # Axis labels with explained variance
  labs(
    x = paste0("PCoA1 (", round(variance_explained_subset_6.2[1] * 100, 1), "%)"), 
    y = paste0("PCoA2 (", round(variance_explained_subset_6.2[2] * 100, 1), "%)")
  ) +
  # Apply color and shape mappings
  scale_shape_manual(values = shape_mapping, name = "Species") +
  scale_color_manual(values = color_map, name = "Species") +
  # Customize the theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "none"
    
  ) +
  # Override legend to apply thicker borders for "Empty" and "Petunia" in the legend
  guides(
    shape = guide_legend(
      override.aes = list(
        size = 12,
        stroke = c(3, 1, 1, 3, 1)  # Specify thicker borders only for "Empty" and "Petunia"
      )
    ),
    color = guide_legend(override.aes = list(size = 10))
  ) +
  # Add the PERMANOVA annotation
  annotate("text", x = Inf, y = Inf, label = legend_text_subset_6_2, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE)+
  labs(title = "", tag = "C") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))

# Display the plot
plot_6_2


# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)

# Subset data for pH 7
subset_physeq_7 <- subset_samples(physeq_css, pH == 7)

# Compute Bray-Curtis distance for the subset
bray_dist_subset_7 <- distance(subset_physeq_7, method = "bray")

# Perform PCoA
pcoa_bray_subset_7 <- ordinate(subset_physeq_7, method = "PCoA", distance = "bray")

# Variance explained by axes
variance_explained_subset_7 <- pcoa_bray_subset_7$values$Relative_eig

# Convert sample data to data frame for the subset
metadata_subset_7 <- as.data.frame(as.matrix(sample_data(subset_physeq_7)))
metadata_subset_7$pH <- as.factor(metadata_subset_7$pH)
metadata_subset_7$Species <- as.factor(metadata_subset_7$Species)

# PERMANOVA analysis for the subset
adonis_results_subset_7 <- adonis2(bray_dist_subset_7 ~ Species, data = metadata_subset_7, permutations = 999)

# Perform pairwise PERMANOVA
pairwise_results_subset_7 <- pairwise.adonis(
  bray_dist_subset_7,  # Your Bray-Curtis distance matrix
  factors = metadata_subset_7$Species,  # Species as factors
  perm = 999  # Number of permutations
)

# Save pairwise PERMANOVA results as CSV
write.csv(pairwise_results_subset_7, "pairwise_results_pH_7.csv")

# Create legend text for PERMANOVA results
legend_text_subset_7 <- paste0("F = ", round(adonis_results_subset_7$F[1], 2),
                               ", R² = ", round(adonis_results_subset_7$R2[1], 3),
                               ", p = ", format.pval(adonis_results_subset_7$Pr[1], digits = 3))

# Prepare ordination data for plotting
ordination_data_subset_7 <- plot_ordination(subset_physeq_7, pcoa_bray_subset_7, color = "pH", shape = "Species", justDF = TRUE)

# Plot for pH 6.2 with PERMANOVA results, custom color mapping, and shape customization
plot_7_0 <- ggplot(ordination_data_subset_7) +
  # Plot points with customized stroke based on species
  geom_jitter(
    aes(x = Axis.1, y = Axis.2, color = Species, shape = Species),
size = 7, width = 0.05, height = 0.05,
    stroke = ifelse(ordination_data_subset_7$Species == "Empty" | ordination_data_subset_7$Species == "Petunia", 3, 1)  # Thicker border for "Empty" and "Petunia"
  ) +
  # Ellipses for each species group
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = Species), type = "norm", linetype = 2, size = 1) +
  # Axis labels with explained variance
  labs(
    x = paste0("PCoA1 (", round(variance_explained_subset_7[1] * 100, 1), "%)"), 
    y = paste0("PCoA2 (", round(variance_explained_subset_7[2] * 100, 1), "%)")
  ) +
  # Apply color and shape mappings
  scale_shape_manual(values = shape_mapping, name = "Species") +
  scale_color_manual(values = color_map, name = "Species") +
  # Customize the theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 20, face = "bold"),
legend.text = element_text(size = 22),  # Increase legend text size    
legend.title = element_blank(),     
legend.key.size = unit(1.3, "cm"),     
legend.spacing.x = unit(0.5, "cm"),     
legend.spacing.y = unit(0.3, "cm"),     
legend.position = "none",
  ) +
  # Override legend to apply thicker borders for "Empty" and "Petunia" in the legend
  guides(
    shape = guide_legend(
      override.aes = list(
        size = 12,
        stroke = c(3, 1, 1, 3, 1)  # Specify thicker borders only for "Empty" and "Petunia"
      )
    ),
    color = guide_legend(override.aes = list(size = 10))
  ) +
  # Add the PERMANOVA annotation
  annotate("text", x = Inf, y = Inf, label = legend_text_subset_7, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE)+
  labs(title = "", tag = "D") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))

# Display the plot
plot_7_0
# Save the plot
ggsave("Figures/PCoA_of_Bray_Curtis_Dissimilarity_pH_7_with_PERMANOVA.jpeg", width = 9, height = 8, dpi = 600, units = "in")
ggsave("Figures/PCoA_of_Bray_Curtis_Dissimilarity_pH_7_with_PERMANOVA.pdf", width = 9, height = 8, dpi = 600, units = "in")
ggsave("Figures/PCoA_of_Bray_Curtis_Dissimilarity_pH_7_with_PERMANOVA.tiff", width = 9, height = 8, dpi = 600, units = "in")



# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)

# Subset data for pH 4.5
subset_physeq_4_5 <- subset_samples(physeq_css, pH == 4.5)

# Compute Bray-Curtis distance for the subset
bray_dist_subset_4_5 <- distance(subset_physeq_4_5, method = "bray")

# Perform PCoA
pcoa_bray_subset_4_5 <- ordinate(subset_physeq_4_5, method = "PCoA", distance = "bray")

# Variance explained by axes
variance_explained_subset_4_5 <- pcoa_bray_subset_4_5$values$Relative_eig

# Convert sample data to data frame for the subset
metadata_subset_4_5 <- as.data.frame(as.matrix(sample_data(subset_physeq_4_5)))
metadata_subset_4_5$pH <- as.factor(metadata_subset_4_5$pH)
metadata_subset_4_5$Species <- as.factor(metadata_subset_4_5$Species)

# PERMANOVA analysis for the subset
adonis_results_subset_4_5 <- adonis2(bray_dist_subset_4_5 ~ Species, data = metadata_subset_4_5, permutations = 999)

# Perform pairwise PERMANOVA
pairwise_results_subset_4_5 <- pairwise.adonis(
  bray_dist_subset_4_5,  # Your Bray-Curtis distance matrix
  factors = metadata_subset_4_5$Species,  # Species as factors
  perm = 999  # Number of permutations
)

# Save pairwise PERMANOVA results as CSV
write.csv(pairwise_results_subset_4_5, "pairwise_results_pH_4.5.csv")

# Create legend text for PERMANOVA results
legend_text_subset_4_5 <- paste0("F = ", round(adonis_results_subset_4_5$F[1], 2),
                                 ", R² = ", round(adonis_results_subset_4_5$R2[1], 3),
                                 ", p = ", format.pval(adonis_results_subset_4_5$Pr[1], digits = 3))

# Prepare ordination data for plotting
ordination_data_subset_4_5 <- plot_ordination(subset_physeq_4_5, pcoa_bray_subset_4_5, color = "pH", shape = "Species", justDF = TRUE)

# Define shapes for species
shape_mapping <- c("Empty" = 21, "Geranium" = 16, "Marigold" = 7, "Petunia" = 3, "Tomato" = 17)

# Plot for pH 6.2 with PERMANOVA results, custom color mapping, and shape customization
plot_4_5 <- ggplot(ordination_data_subset_4_5) +
  # Plot points with customized stroke based on species
  geom_jitter(
    aes(x = Axis.1, y = Axis.2, color = Species, shape = Species),
size = 7, width = 0.05, height = 0.05,
    stroke = ifelse(ordination_data_subset_4_5$Species == "Empty" | ordination_data_subset_4_5$Species == "Petunia", 3, 1)  # Thicker border for "Empty" and "Petunia"
  ) +
  # Ellipses for each species group
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = Species), type = "norm", linetype = 2, size = 1) +
  # Axis labels with explained variance
  labs(
    x = paste0("PCoA1 (", round(variance_explained_subset_4_5[1] * 100, 1), "%)"), 
    y = paste0("PCoA2 (", round(variance_explained_subset_4_5[2] * 100, 1), "%)")
  ) +
  # Apply color and shape mappings
  scale_shape_manual(values = shape_mapping, name = "Species") +
  scale_color_manual(values = color_map, name = "Species") +
  # Customize the theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 20, face = "bold"),
legend.text = element_text(size = 22),  # Increase legend text size     
legend.title = element_blank(),    
legend.key.size = unit(1.3, "cm"),     
legend.spacing.x = unit(0.5, "cm"),     
legend.spacing.y = unit(0.3, "cm"),     
legend.position = "none"
  ) +
  # Override legend to apply thicker borders for "Empty" and "Petunia" in the legend
  guides(
    shape = guide_legend(
      override.aes = list(
        size = 12,
        stroke = c(3, 1, 1, 3, 1)  # Specify thicker borders only for "Empty" and "Petunia"
      )
    ),
    color = guide_legend(override.aes = list(size = 10))
  ) +
  # Add the PERMANOVA annotation
  annotate("text", x = Inf, y = Inf, label = legend_text_subset_4_5, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE)+
  labs(title = "", tag = "A") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))

# Display the plot
plot_4_5


# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)

# Subset data for pH 5.5
subset_physeq_5_5 <- subset_samples(physeq_css, pH == 5.5)

# Compute Bray-Curtis distance for the subset
bray_dist_subset_5_5 <- distance(subset_physeq_5_5, method = "bray")

# Perform PCoA
pcoa_bray_subset_5_5 <- ordinate(subset_physeq_5_5, method = "PCoA", distance = "bray")

# Variance explained by axes
variance_explained_subset_5_5 <- pcoa_bray_subset_5_5$values$Relative_eig

# Convert sample data to data frame for the subset
metadata_subset_5_5 <- as.data.frame(as.matrix(sample_data(subset_physeq_5_5)))
metadata_subset_5_5$pH <- as.factor(metadata_subset_5_5$pH)
metadata_subset_5_5$Species <- as.factor(metadata_subset_5_5$Species)

# PERMANOVA analysis for the subset
adonis_results_subset_5_5 <- adonis2(bray_dist_subset_5_5 ~ Species, data = metadata_subset_5_5, permutations = 999)

# Perform pairwise PERMANOVA
pairwise_results_subset_5_5 <- pairwise.adonis(
  bray_dist_subset_5_5,  # Your Bray-Curtis distance matrix
  factors = metadata_subset_5_5$Species,  # Species as factors
  perm = 999  # Number of permutations
)

# Save pairwise PERMANOVA results as CSV
write.csv(pairwise_results_subset_5_5, "pairwise_results_pH_5.5.csv")

# Create legend text for PERMANOVA results
legend_text_subset_5_5 <- paste0("F = ", round(adonis_results_subset_5_5$F[1], 2),
                                 ", R² = ", round(adonis_results_subset_5_5$R2[1], 3),
                                 ", p = ", format.pval(adonis_results_subset_5_5$Pr[1], digits = 3))

# Prepare ordination data for plotting
ordination_data_subset_5_5 <- plot_ordination(subset_physeq_5_5, pcoa_bray_subset_5_5, color = "pH", shape = "Species", justDF = TRUE)

# Define shapes for species
shape_mapping <- c("Empty" = 21, "Geranium" = 16, "Marigold" = 7, "Petunia" = 3, "Tomato" = 17)

# Plot for pH 6.2 with PERMANOVA results, custom color mapping, and shape customization
plot_5_5 <- ggplot(ordination_data_subset_5_5) +
  # Plot points with customized stroke based on species
  geom_jitter(
    aes(x = Axis.1, y = Axis.2, color = Species, shape = Species),
size = 7, width = 0.05, height = 0.05,
    stroke = ifelse(ordination_data_subset_5_5$Species == "Empty" | ordination_data_subset_5_5$Species == "Petunia", 3, 1)  # Thicker border for "Empty" and "Petunia"
  ) +
  # Ellipses for each species group
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = Species), type = "norm", linetype = 2, size = 1) +
  # Axis labels with explained variance
  labs(
    x = paste0("PCoA1 (", round(variance_explained_subset_5_5[1] * 100, 1), "%)"), 
    y = paste0("PCoA2 (", round(variance_explained_subset_5_5[2] * 100, 1), "%)")
  ) +
  # Apply color and shape mappings
  scale_shape_manual(values = shape_mapping, name = "Species") +
  scale_color_manual(values = color_map, name = "Species") +
  # Customize the theme
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 28, color = "black"),
    plot.title = element_text(size = 20, face = "bold"),
legend.text = element_text(size = 22),  # Increase legend text size     
legend.title = element_blank(),     
legend.key.size = unit(1.3, "cm"),     
legend.spacing.x = unit(0.5, "cm"),     
legend.spacing.y = unit(0.3, "cm"),     
legend.position = "none"
  ) +
  # Override legend to apply thicker borders for "Empty" and "Petunia" in the legend
  guides(
    shape = guide_legend(
      override.aes = list(
        size = 12,
        stroke = c(3, 1, 1, 3, 1)  # Specify thicker borders only for "Empty" and "Petunia"
      )
    ),
    color = guide_legend(override.aes = list(size = 10))
  ) +
  # Add the PERMANOVA annotation
  annotate("text", x = Inf, y = Inf, label = legend_text_subset_5_5, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE)+
  labs(title = "", tag = "B") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))

# Display the plot
plot_5_5

# Extract the legend from endpH
#legend <- get_legend(plot_7_0 + 
#                       theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0)))

##
library(cowplot)
# Overlay "pH 4.5" label exactly near tag "A"
plot_4_5 <- ggdraw(plot_4_5) +
  draw_label("pH 4.5", x = 0.15, y = 0.975, hjust = 0, fontface = "bold", size = 28)
plot_4_5
plot_5_5 <- ggdraw(plot_5_5) +
  draw_label("pH 5.5", x = 0.15, y = 0.975, hjust = 0, fontface = "bold", size = 28)
plot_5_5
plot_6_2 <- ggdraw(plot_6_2) +
  draw_label("pH 6.2", x = 0.15, y = 0.975, hjust = 0, fontface = "bold", size = 28)
plot_6_2
plot_7_0 <- ggdraw(plot_7_0) +
  draw_label("pH 7.0", x = 0.15, y = 0.975, hjust = 0, fontface = "bold", size = 28) + theme(legend.position = "none")
plot_7_0

# Extract the legend from one of the plots
plot_4_5 <- plot_4_5 + theme(legend.position = "none")
plot_5_5 <- plot_5_5 + theme(legend.position = "none")
plot_6_2 <- plot_6_2 + theme(legend.position = "none")
plot_7_0 <- plot_7_0 + theme(legend.position = "none")

# Create a row of plots first
plots_row <- plot_grid(
  plot_4_5, 
  plot_5_5,
  plot_6_2,
  plot_7_0,
  ncol = 2,
  align = 'h',
  axis = 'b'
)

# Combine plots row with centered legend
combined_plot <- plot_grid(
  plots_row,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1),  # Adjust the 0.1 to control legend space
  align = 'v',
  axis = 'l'
)
combined_plot
# Save as high-resolution JPEG
ggsave("combined_plot.jpg", combined_plot, 
       width = 16, height = 14,  # 16:9 aspect ratio
       dpi = 300, quality = 100)



# Save the combined plot
ggsave("Figures/Fig3_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.jpeg", plot = combined_plot, width = 16, height = 14, dpi = 600)
ggsave("Figures/Fig3N_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.jpeg", plot = combined_plot, width = 16, height = 14,dpi = 600, quality = 100)
ggsave("Figures/Fig3_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.pdf", plot = combined_plot, width = 16, height = 14, dpi = 600)
ggsave("Figures/Fig3_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeightedl.tiff", plot = combined_plot, width = 16, height = 14, dpi = 600)


