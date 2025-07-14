# Load required libraries
library(ggplot2)
library(dplyr)

# Step 1: Transform counts to relative abundance
Bacteria_ps_filtered <- readRDS("Bacteria_ps_filtered_clean.rds")
Bacteria_ps_relabund <- transform_sample_counts(Bacteria_ps_filtered, function(x) x / sum(x))
composite <- psmelt(Bacteria_ps_relabund)

# Step 2: Filter for the top 100% phyla
get_ALL_Phylum <- function(data) {
  data %>%
    group_by(Phylum) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    top_n(n = floor(1 * n()), wt = mean_abundance) %>%
    pull(Phylum)
}

top_50_percent_Phylum <- get_ALL_Phylum(composite)

composite_top50 <- composite %>%
  filter(Phylum %in% top_50_percent_Phylum)

# Step 3: Calculate mean relative abundance for jitter points
composite_mean <- composite_top50 %>%
  group_by(Phylum, pH, sample_Species) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE) * 100) %>%  # Convert to percentage
  ungroup()
# Save the composite_mean data frame as a CSV file
write.csv(composite_mean, file = "Phylum_composite_mean_relative_abundance.csv", row.names = FALSE)

# Define custom colors for pH and shapes for species
custom_colors2 <- c("4.5" = "#E69F00", "5.5" = "#56B4E9", "6.2" = "#009E73", "7" = "#d62728")
custom_shapes <- c("Empty" = 21, "Geranium" = 16, "Marigold" = 7, "Petunia" = 3, "Tomato" = 17)

# Load required libraries
library(ggplot2)
library(dplyr)

# ... (previous code remains unchanged)

# Create the boxplot with jitter points
p1 <- ggplot(composite_mean, aes(x = Phylum, y = Mean_Abundance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(data = composite_mean, aes(x = Phylum, y = Mean_Abundance, color = as.numeric(pH), shape = sample_Species), 
              width = 0.2, height = 0, size = 3, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#E69F00", "#56B4E9", "#009E73", "#d62728"),
    values = scales::rescale(c(4.5, 5.5, 6.2, 7)),
    breaks = c(4.5, 5.5, 6.2, 7),
    labels = c("4.5", "5.5", "6.2", "7"),
    name = "pH"
  ) +
  scale_shape_manual(values = custom_shapes) +
  labs(
    title = "",
    x = "Phylum",
    y = "Relative abundance (%)",
    color = "pH",
    shape = "Species"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.title = element_text(size = 24),
    axis.ticks = element_line(color = "black", size = 0.5),  # Add axis ticks
    axis.ticks.length = unit(0.1, "cm"),  # Set tick length
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.key.size = unit(1.3, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.position = "bottom"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 10)),
    color = guide_colorbar(title = "pH")
  )
p1
####Genus

# Step 2: Filter for the top 5% genera
get_top_5_percent_Genus <- function(data) {
  data %>%
    group_by(Genus) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    top_n(n = floor(0.05 * n()), wt = mean_abundance) %>%
    pull(Genus)
}

top_5_percent_Genus <- get_top_5_percent_Genus(composite)

composite_top5 <- composite %>%
  filter(Genus %in% top_5_percent_Genus)

# Step 2.5: Shorten genus names
composite_top5 <- composite_top5 %>%
  mutate(
    Genus = case_when(
      Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" ~ "Rhizobium",
      Genus == "Burkholderia-Caballeronia-Paraburkholderia" ~ "Burkholderia",
      TRUE ~ Genus
    )
  )

# Step 3: Calculate mean relative abundance for jitter points
composite_mean <- composite_top5 %>%
  group_by(Genus, pH, sample_Species) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE) * 100) %>%  # Convert to percentage
  ungroup()

# Define custom colors for pH and shapes for species
custom_colors2 <- c("4.5" = "#E69F00", "5.5" = "#56B4E9", "6.2" = "#009E73", "7" = "#d62728")
custom_shapes <- c("Empty" = 21, "Geranium" = 16, "Marigold" = 7, "Petunia" = 3, "Tomato" = 17)



# Create the boxplot with jitter points
p2 <- ggplot(composite_mean, aes(x = Genus, y = Mean_Abundance)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(data = composite_mean, aes(x = Genus, y = Mean_Abundance, color = as.numeric(pH), shape = sample_Species), 
              width = 0.2, height = 0, size = 3, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#E69F00", "#56B4E9", "#009E73", "#d62728"),
    values = scales::rescale(c(4.5, 5.5, 6.2, 7)),
    breaks = c(4.5, 5.5, 6.2, 7),
    labels = c("4.5", "5.5", "6.2", "7"),
    name = "pH"
  ) +
  scale_shape_manual(values = custom_shapes) +  # Custom shapes for species
  labs(
    title = "",
    x = "Genus",
    y = "Relative abundance(%)",
    color = "pH",
    shape = "Species"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.title = element_text(size = 24),
    axis.ticks = element_line(color = "black", size = 0.5),  # Add axis ticks
    axis.ticks.length = unit(0.1, "cm"),  # Set tick length
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.size = unit(1.3, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    legend.position = "bottom"
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 10)),
    color = guide_colorbar(title = "pH")
  )
p2


legend_components <- get_plot_component(p1, "guide-box", return_all = TRUE)
selected_legend <- legend_components[[3]]

# Create a grid of plots without legends
plots_row <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  labels = c("A", "B"),
  label_size = 20,
  nrow = 2
)

# Add the extracted legend at the bottom
final_plot <- plot_grid(
  plots_row, 
  selected_legend, 
  ncol = 1,
  rel_heights = c(1, 0.1)
)

# Display the final plot
final_plot

# Save the plot
ggsave("Figures/Fig4_Phylum_Genus_Boxplot_With_Jitter.jpeg", plot = final_plot, width = 14, height = 16, dpi = 600)
ggsave("Figures/Fig4_Phylum_Genus_Boxplot_With_Jitter.pdf", plot = final_plot, width = 14, height = 16, dpi = 600)
ggsave("Figures/Fig4_Phylum_Genus_Boxplot_With_Jitter.tiff", plot = final_plot, width = 14, height = 16, dpi = 600, device = "tiff")


