# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(ggrepel)
library(multcompView)
library(car)
library(ggpubr)
# Read the rarefied phyloseq object
physeq_rarefied <- readRDS("physeq_rarefied.rds")

# Calculate alpha diversity measures
alpha_diversity <- estimate_richness(physeq_rarefied, measures = c("Observed", "Shannon", "Simpson", "Chao1"))
rownames(alpha_diversity) <- gsub("^X", "", rownames(alpha_diversity))

# Extract sample data and merge with alpha diversity metrics
meta <- data.frame(sample_data(physeq_rarefied))
alpha_diversity$SampleID <- rownames(alpha_diversity)
alpha_diversity_merged <- merge(alpha_diversity, meta, by = "row.names")

# Define treatment mapping and custom colors
treatment_map <- c(
  "Empty_4.5" = "E4.5", "Petunia_4.5" = "P4.5", "Tomato_4.5" = "T4.5", "Marigold_4.5" = "M4.5", "Geranium_4.5" = "G4.5",
  "Empty_5.5" = "E5.5", "Petunia_5.5" = "P5.5", "Tomato_5.5" = "T5.5", "Marigold_5.5" = "M5.5", "Geranium_5.5" = "G5.5",
  "Empty_6.2" = "E6.2", "Petunia_6.2" = "P6.2", "Tomato_6.2" = "T6.2", "Marigold_6.2" = "M6.2", "Geranium_6.2" = "G6.2",
  "Empty_7" = "E7.0", "Petunia_7" = "P7.0", "Tomato_7" = "T7.0", "Marigold_7" = "M7.0", "Geranium_7" = "G7.0"
)

custom_colors2 <- c("4.5" = "#E69F00", "5.5" = "#56B4E9", "6.2" = "#009E73", "7" = "#d62728")

# Apply treatment mapping and convert pH to factor
alpha_diversity_merged1 <- alpha_diversity_merged %>%
  mutate(Treatments = str_replace_all(Treatments, treatment_map),
         pH = as.factor(pH))

# Calculate means for each treatment
means <- alpha_diversity_merged1 %>%
  group_by(Treatments) %>%
  summarize(mean_pH = mean(endpH, na.rm = TRUE))

# Create endpH plot
endpH <- ggplot(alpha_diversity_merged1, aes(x = Treatments, y = endpH, fill = pH)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6, color = "black") +
  theme_classic() +
  labs(x = "", y = "pH") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.6, vjust = 0.6, size = 28, color = "black"),
    axis.text.y = element_text(size = 28, color = "black"),
    axis.title = element_text(size = 32),
    axis.ticks = element_line(size = 1, color = "black", lineend = "square"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.8),
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 34),
    legend.key.size = unit(1.5, "cm")
  ) +
  scale_y_continuous(limits = c(3, 8), breaks = seq(3, 8, 0.5)) +
  scale_fill_manual(values = custom_colors2) +
  geom_text(data = means, aes(x = Treatments, y = mean_pH + 0.7, label = round(mean_pH, 2)),
            inherit.aes = FALSE, size = 9, color = "darkblue")+
labs(title = "", tag = "A") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))
endpH
# Perform statistical tests
anova_model <- aov(Shannon ~ Treatments, data = alpha_diversity_merged)
shapiro.test(residuals(anova_model))
leveneTest(Shannon ~ Treatments, data = alpha_diversity_merged)
kruskal_test <- kruskal.test(Shannon ~ Treatments, data = alpha_diversity_merged)

# Pairwise Wilcoxon test
pairwise_test <- compare_means(Shannon ~ Treatments, data = alpha_diversity_merged, method = "wilcox.test", p.adjust.method = "BH")

# Generate significance letters
pairwise_pvals <- setNames(pairwise_test$p.adj, paste(pairwise_test$group1, pairwise_test$group2, sep = "-"))
signif_letters <- multcompLetters(pairwise_pvals, threshold = 0.05)$Letters

# Create a dataframe mapping treatments to letters
treatment_letters <- data.frame(Treatments = unique(c(pairwise_test$group1, pairwise_test$group2)),
                                Letters = signif_letters)

# Merge significant letters with alpha_diversity_merged dataset
alpha_diversity_merged <- merge(alpha_diversity_merged, treatment_letters, by = "Treatments", all.x = TRUE)

# Filter and prepare data for Shannon plot
alpha_diversity_filtered <- alpha_diversity_merged %>%
  filter(!is.na(Shannon)) %>%
  mutate(Treatments = str_replace_all(Treatments, treatment_map),
         pH = as.factor(pH))

treatment_letters_unique <- alpha_diversity_filtered %>%
  group_by(Treatments) %>%
  summarise(
    Letters = dplyr::first(Letters),
    Shannon = max(Shannon, na.rm = TRUE),
    pH = dplyr::first(pH)
  )

# Create Shannon plot
Shano <- ggplot(alpha_diversity_filtered, aes(x = Treatments, y = Shannon, fill = pH)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6, color = "black") +
  theme_classic() +
  labs(x = "Treatments", y = "Shannon diversity") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6, size = 28, color = "black"),
    axis.text.y = element_text(size = 28, color = "black"),
    axis.title.x = element_text(size = 32, margin = margin(t = 15)),  # <-- This line adds space
    axis.title.y = element_text(size = 32),
    panel.border = element_rect(color = "black", fill = NA, size = 1.8),
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 34),
    legend.key.size = unit(1.5, "cm")
  ) +
  scale_y_continuous(limits = c(2.5, 7.5), breaks = seq(2.5, 7.5, 0.5)) +
  scale_fill_manual(values = custom_colors2) +
  geom_text_repel(
    data = treatment_letters_unique,
    aes(x = Treatments, y = Shannon + 0.3, label = Letters, color = pH),
    size = 9.5, box.padding = 0.6, point.padding = 0.5,
    nudge_y = 0.3, direction = "y", max.overlaps = 20, segment.color = NA
  ) +
  scale_color_manual(values = custom_colors2) +
  theme(legend.position = "bottom") +
  labs(title = "", tag = "B") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text(size = 28, face = "bold"),
        plot.tag.position = c(0.020, 0.985))
Shano

# Remove legends from both plots
endpH_no_legend <- endpH + theme(legend.position = "none")
Shano_no_legend <- Shano + theme(legend.position = "none")


# Extract the legend from endpH
legend <- get_legend(endpH + 
                       theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0)))

# Combine the plots and add the legend
combined_plot <- plot_grid(
  endpH_no_legend,
  Shano_no_legend,
  legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.2),
  align = 'v',
  axis = 'l'
)


# Save the plot in different formats
ggsave("Figures/Fig1_Shannon_Diversity_Plot_Wilcoxon_pH.jpeg", plot = combined_plot, width = 16, height = 14, dpi = 600, units = "in")
ggsave("Figures/Fig1_Shannon_Diversity_Plot_Wilcoxon_pH.pdf", plot = combined_plot, width = 16, height = 14, dpi = 600, units = "in")
ggsave("Figures/Fig1_Shannon_Diversity_Plot_Wilcoxon_pH.tiff", plot = combined_plot, width = 16, height = 14, dpi = 600, units = "in")
