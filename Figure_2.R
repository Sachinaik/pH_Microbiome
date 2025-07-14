# Load necessary libraries
library(vegan)          # For beta diversity analysis
library(phyloseq)       # Microbiome data handling
library(ggplot2)        # For visualization
library(dplyr)          # Data manipulation
library(ape)            # For phylogenetic trees, if needed

# Example microbiome data (you'll replace this with your actual data)
# Assume 'physeq' is your phyloseq object containing your microbiome data
# 'meta_data' is your sample metadata

Bacteria_ps_filtered <- readRDS("Bacteria_ps_filtered_clean.rds")

Bacteria_ps_filtered
# Create a filtered phyloseq object without mitochondria and chloroplast
Bacteria_ps_filtered_clean <- Bacteria_ps_filtered %>%
  subset_taxa(
    !grepl("mitochondria", Family, ignore.case = TRUE) &
      !grepl("mitochondria", Genus, ignore.case = TRUE) &
      !grepl("mitochondria", Order, ignore.case = TRUE) &
      !grepl("chloroplast", Family, ignore.case = TRUE) &
      !grepl("chloroplast", Order, ignore.case = TRUE) &
      !grepl("chloroplast", Class, ignore.case = TRUE) &
      !grepl("chloroplast", Genus, ignore.case = TRUE)
  )
Bacteria_ps_filtered <- Bacteria_ps_filtered_clean
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


# 2. **Beta Diversity Calculation**
# Compute beta diversity using Bray-Curtis dissimilarity (abundance-based)
bray_dist <- distance(physeq_css, method = "bray")

# 3. **Ordination - PCoA**
# Perform Principal Coordinate Analysis (PCoA) using Bray-Curtis distance
pcoa_bray <- ordinate(physeq_css, method = "PCoA", distance = "bray")
pcoa_var1 <- round(100 * pcoa_bray$values$Relative_eig[1], 1)  # PCoA1 variance
pcoa_var2 <- round(100 * pcoa_bray$values$Relative_eig[2], 1)  # PCoA2 variance

# Convert the sample_data into a proper data frame
metadata <- as.data.frame(as.matrix(sample_data(physeq_css)))

# Check if 'pH' and 'Species' columns exist and are correct
str(metadata)  # This should show 'pH' and 'Species'

# Ensure 'pH' is numeric and 'Species' is a factor
metadata$pH <- as.numeric(metadata$pH)
metadata$Species <- as.factor(metadata$Species)

# Convert pH to a factor
metadata$pH <- as.factor(metadata$pH)

# Rerun the adonis analysis
adonis_results <- adonis2(bray_dist ~ pH + Species, data = metadata, permutations = 999)


# Run PERMANOVA with pH and Species as predictors
adonis_results <- adonis2(bray_dist ~ pH + Species, data = metadata, 
                          permutations = 999)
# PERMANOVA with only pH as the predictor
adonis_pH <- adonis2(bray_dist ~ pH, data = metadata, permutations = 999)

# Display results
print(adonis_pH)
# PERMANOVA with only Species as the predictor
adonis_Species <- adonis2(bray_dist ~ Species, data = metadata, permutations = 999)

# Display results
print(adonis_Species)

# PERMANOVA with interaction between pH and Species
adonis_interaction <- adonis2(bray_dist ~ pH * Species, data = metadata, permutations = 999)

# Display results
print(adonis_interaction)

# Calculate PCoA variances
variance_explained <- pcoa_bray$values$Relative_eig

# Create legend text for PERMANOVA results
#legend_text <- paste0("PERMANOVA: F = ", round(adonis_results$F[1], 2),
#                      ", R² = ", round(adonis_results$R2[1], 3),
#                      ", p = ", format.pval(adonis_results$Pr[1], digits = 3))
# Create legend text for PERMANOVA results
legend_text <- paste0(" F = ", round(adonis_pH$F[1], 2),
                      ", R² = ", round(adonis_pH$R2[1], 3),
                      ", p = ", format.pval(adonis_pH$Pr[1], digits = 3))

# Suppress default point layer in plot_ordination
ordination_data <- plot_ordination(physeq_css, pcoa_bray, color = "pH", shape = "Species", justDF = TRUE)

# Convert pH in ordination_data to numeric (if not done already)
ordination_data$pH <- as.numeric(as.character(ordination_data$pH))

# Define shapes for species using the specified shapes
shape_mapping <- c("Empty" = 21,        # Empty pot shape (open square)
                   "Geranium" = 16,    # Filled circle
                   "Marigold" = 7,    # Filled triangle
                   "Petunia" = 3,      # Plus
                   "Tomato" = 17)       # Cross mark (filled X)

# Plot PCoA with jitter and no duplication of points
# Define custom colors for each pH level
custom_colors2 <- c("4.5" = "#E69F00", "5.5" = "#56B4E9", "6.2" = "#009E73", "7" = "#d62728")

bray <- ggplot(ordination_data) +
  geom_jitter(aes(x = Axis.1, y = Axis.2, color = pH, shape = Species), size = 7, width = 0.05, height = 0.05) +  # Use continuous pH
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = factor(pH)), type = "norm", linetype = 2, size = 1) +  # Group by factor(pH)
  labs(
    title = "", 
    x = paste0("PCoA1 (", round(variance_explained[1] * 100, 1), "%)"), 
    y = paste0("PCoA2 (", round(variance_explained[2] * 100, 1), "%)")
  ) +
  scale_color_gradientn(
    colors = c("#E69F00", "#56B4E9", "#009E73", "#d62728"),  # Gradient colors
    values = scales::rescale(c(4.5, 5.5, 6.2, 7)),  # Rescale pH values
    breaks = c(4.5, 5.5, 6.2, 7),  # Specify which values to show in the legend
    name = "pH"  # Legend title
  ) +
  scale_shape_manual(
    values = shape_mapping,  # Shape mapping for Species
    name = "Species"  # Legend title for shapes
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Box around plot
    axis.text = element_text(size = 18),  # Larger text for axis labels
    axis.title = element_text(size = 28),  # Larger text for axis titles
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 18),
    legend.position = "bottom"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8))  # Adjust size of shapes in legend
  ) +
  annotate("text", x = Inf, y = Inf, label = legend_text, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE) +
  labs(title = "", tag = "A") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))

bray

# 2. Perform Betadisper (to check homogeneity of dispersion across groups)
betadisper_result <- betadisper(bray_dist, sample_data(physeq_css)$pH)
# Print betadisper results (this shows the dispersion for each group)
print(betadisper_result)
# Perform ANOVA to test for homogeneity of dispersion
anova_betadisperpH <- anova(betadisper_result)
print(anova_betadisperpH)

betadisper_result <- betadisper(bray_dist, sample_data(physeq_css)$Species)
# Print betadisper results (this shows the dispersion for each group)
print(betadisper_result)
# Perform ANOVA to test for homogeneity of dispersion
anova_betadisperSpecies <- anova(betadisper_result)
print(anova_betadisperSpecies)

betadisper_result <- betadisper(bray_dist, sample_data(physeq_css)$Treatments)
# Print betadisper results (this shows the dispersion for each group)
print(betadisper_result)
# Perform ANOVA to test for homogeneity of dispersion
anova_betadisperTreatments <- anova(betadisper_result)
print(anova_betadisperTreatments)

#########

# 1. **Beta Diversity Calculation using Weighted UniFrac**
# Calculate the weighted UniFrac distance
uni_frac_dist <- distance(physeq_css, method = "wunifrac")

# 2. **Ordination - PCoA**
# Perform Principal Coordinate Analysis (PCoA) using Weighted UniFrac distance
pcoa_uni_frac <- ordinate(physeq_css, method = "PCoA", distance = "wunifrac")

# Calculate variances explained by the axes
pcoa_var1 <- round(100 * pcoa_uni_frac$values$Relative_eig[1], 1)  # PCoA1 variance
pcoa_var2 <- round(100 * pcoa_uni_frac$values$Relative_eig[2], 1)  # PCoA2 variance

# Convert the sample_data into a proper data frame
metadata <- as.data.frame(as.matrix(sample_data(physeq_css)))

# Check if 'pH' and 'Species' columns exist and are correct
str(metadata)  # This should show 'pH' and 'Species'

# Ensure 'pH' is numeric and 'Species' is a factor
metadata$pH <- as.numeric(metadata$pH)
metadata$Species <- as.factor(metadata$Species)

# Convert pH to a factor
metadata$pH <- as.factor(metadata$pH)

# 3. **Plot PCoA with Weighted UniFrac**
# Create a data frame for plotting
ordination_data <- plot_ordination(physeq_css, pcoa_uni_frac, color = "pH", shape = "Species", justDF = TRUE)

# Assuming you have performed the PERMANOVA and stored the results
# Run PERMANOVA with UniFrac distance and pH, Species, and their interaction as predictors
adonis_uni_pH <- adonis2(uni_frac_dist ~ pH, data = metadata, permutations = 999)
adonis_uni_Species <- adonis2(uni_frac_dist ~ Species, data = metadata, permutations = 999)
adonis_uni_interaction <- adonis2(uni_frac_dist ~ pH * Species, data = metadata, permutations = 999)
# Run betadisper for pH and Species (homogeneity of dispersion)
anova_betadisper_uni_pH <- anova(betadisper(uni_frac_dist, metadata$pH))
anova_betadisper_uni_Species <- anova(betadisper(uni_frac_dist, metadata$Species))

# Create a summary of the results for plotting
permanova_text <- paste0(" F = ", round(adonis_uni_pH$F[1], 2),
                         ", R² = ", round(adonis_uni_pH$R2[1], 3),
                         ", p = ", format.pval(adonis_uni_pH$Pr[1], digits = 3))
# Plot PCoA with jitter and no duplication of points
uniW <- ggplot(ordination_data) +
  geom_jitter(aes(x = Axis.1, y = Axis.2, color = pH, shape = Species), size = 7, width = 0.05, height = 0.05) +
  stat_ellipse(aes(x = Axis.1, y = Axis.2, group = factor(pH)), type = "norm", linetype = 2, size = 1) +
  labs(
    title = "", 
    x = paste0("PCoA1 (", pcoa_var1, "%)"), 
    y = paste0("PCoA2 (", pcoa_var2, "%)")
  ) +
  scale_color_gradientn(
    colors = c("#E69F00", "#56B4E9", "#009E73", "#d62728"),  # Gradient colors
    values = scales::rescale(c(4.5, 5.5, 6.2, 7)),  # Rescale pH values
    breaks = c(4.5, 5.5, 6.2, 7),  # Specify which values to show in the legend
    name = "pH"  # Legend title
  ) +
  scale_shape_manual(
    values = shape_mapping,  # Use the shape mapping defined above
    name = "Species"  # Set a name for the shape legend
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Box around plot
    axis.text = element_text(size = 18),  # Larger text for axis labels
    axis.title = element_text(size = 28),  # Larger text for axis titles
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 22),  # Increase legend text size
    legend.title = element_blank(),
    legend.key.size = unit(1.3, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.position = "bottom"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 10))  # Keep shapes for the Species legend
  ) +
  annotate("text", x = Inf, y = Inf, label = permanova_text, 
           hjust = 1.1, vjust = 1.2, size = 7, color = "black", 
           fontface = "italic", parse = FALSE)+
  labs(title = "", tag = "B") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))
uniW


# Extract the legend from endpH
legend <- get_legend(uniW + 
                       theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0)))

# Extract the legend from one of the plots
bray <- bray + theme(legend.position = "none")
uniW <- uniW + theme(legend.position = "none")
# Create a row of plots first
plots_row <- plot_grid(
  bray, 
  uniW,
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
       width = 16, height = 9,  # 16:9 aspect ratio
       dpi = 300, quality = 100)

# Save the combined plot
ggsave("Figures/Fig2_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.jpeg", plot = combined_plot, width = 16, height = 12, dpi = 600)
ggsave("Figures/Fig2N_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.jpeg", plot = combined_plot, width = 16, height = 9,dpi = 300, quality = 100)
ggsave("Figures/Fig2_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeighted.pdf", plot = combined_plot, width = 16, height = 12, dpi = 600)
ggsave("Figures/Fig2_Combined_PCoA_Bray_Curtis_Dissimilarity_UnifracWeightedl.tiff", plot = combined_plot, width = 16, height = 12, dpi = 600)



