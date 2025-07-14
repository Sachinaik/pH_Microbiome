########################################
# Load necessary libraries
########################################
library(phyloseq)
library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)

########################################
# 1) Read Data and Prepare Phyloseq
########################################
physeq_filterednochloromito <- readRDS("Bacteria_ps_filtered_nochlormito.rds")

# Replace with your actual metadata CSV file
sample_info_tab <- read.csv("Leachate_MD.csv", 
                            header = TRUE, 
                            row.names = 1,
                            check.names = FALSE)

# Convert sample_info_tab to a sample_data object
sample_data_object <- sample_data(sample_info_tab)

# Assuming you already have the Phyloseq object 'physeq_filterednochloromito'
# Merge the sample data with the existing phyloseq object
merged_physeq <- merge_phyloseq(physeq_filterednochloromito, sample_data_object)

# Normalize Data
merged_physeq <- transform_sample_counts(merged_physeq, function(x) x / sum(x))

########################################
# 2) Extract Genus-level data
########################################
Genus_data <- tax_glom(merged_physeq, taxrank = "Genus")
Genus_rel_abund <- transform_sample_counts(Genus_data, function(x) x / sum(x))

# Extract the Genus names from the tax_table
tax_table_df <- as.data.frame(tax_table(Genus_rel_abund))
Genus_names <- tax_table_df$Genus

# Shorten specific genus names
Genus_names <- ifelse(
  Genus_names %in% c("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"), 
  "Rhizobium", 
  ifelse(
    Genus_names %in% c("Burkholderia-Caballeronia-Paraburkholderia"), 
    "Burkholderia", 
    Genus_names
  )
)

# Create a new otu_table with Genus names as column names
Genus_table <- otu_table(Genus_rel_abund)
Genus_table <- t(Genus_table) # Transpose to have samples as rows
colnames(Genus_table) <- Genus_names

########################################
# 3) Select all Genera
########################################
Genus_sums <- colSums(Genus_table)
all_percent <- Genus_sums >= quantile(Genus_sums, 0)
all_phyla <- Genus_table[, all_percent, drop = FALSE]

########################################
# 4) Combine with Nutrient Data
########################################
# List of nutrients in your metadata
nutrients <- c(
  "NPOC", "TN", "Al", "B", "Cu", "Fe", 
  "K", "Mn", "Mo","P", "Zn", "Ca"
)

# Extract nutrient data from sample_data
nutrient_data <- data.frame(sample_data(merged_physeq)[, nutrients])

# Combine Genus and nutrient data with Treatment info
combined_data <- cbind(
  as.data.frame(all_phyla),
  nutrient_data,
  Treatment = sample_data(merged_physeq)$Treatment
)

########################################
# 5) Define correlation functions
########################################

# A function to compute correlation and p-values
correlation_with_pvalues <- function(Genus_data, nutrient_data, method = "spearman") {
  # Correlation matrix
  cor_matrix <- cor(Genus_data, nutrient_data, method = method, use = "complete.obs")
  
  # Initialize p-value matrix
  p_matrix <- matrix(
    NA, 
    nrow = ncol(Genus_data), 
    ncol = ncol(nutrient_data),
    dimnames = list(colnames(Genus_data), colnames(nutrient_data))
  )
  
  # Loop through each pair, run cor.test, save p-values
  for (i in 1:ncol(Genus_data)) {
    for (j in 1:ncol(nutrient_data)) {
      test <- cor.test(Genus_data[, i], nutrient_data[, j], method = method)
      p_matrix[i, j] <- test$p.value
    }
  }
  
  list(correlation = cor_matrix, pvalues = p_matrix)
}

# A function to create correlation heatmaps
# Only plots correlations that are both significant (p < 0.05) AND strong (|r| >= 0.5)
create_correlation_plot <- function(cor_matrix, p_matrix, title) {
  plot_data <- reshape2::melt(cor_matrix) %>%
    dplyr::rename(Genus = Var1, Nutrient = Var2, Correlation = value) %>%
    dplyr::left_join(
      reshape2::melt(p_matrix) %>%
        dplyr::rename(Genus = Var1, Nutrient = Var2, Pvalue = value),
      by = c("Genus", "Nutrient")
    )
  
  # Define significance levels and labels
  plot_data <- plot_data %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        Pvalue < 0.001 ~ "***",
        Pvalue < 0.01  ~ "**",
        Pvalue < 0.05  ~ "*",
        TRUE           ~ ""
      ),
      Label = ifelse(
        Significance != "",
        paste0(format(round(Correlation, 2), nsmall = 2), Significance),
        ""
      )
    )
  
  # Filter for both significant correlations (p < 0.05) AND |correlation| >= 0.5
  plot_data_significant <- dplyr::filter(plot_data, Pvalue < 0.05 & abs(Correlation) >= 0.5)
  
  # First convert to character to ensure proper string sorting
  sorted_genera <- sort(as.character(unique(plot_data_significant$Genus)))
  
  # Create a factor with explicit levels in alphabetical order
  plot_data_significant$Genus <- factor(plot_data_significant$Genus, 
                                        levels = sorted_genera)
  
  # Sort data frame by Genus to ensure correct order
  plot_data_significant <- plot_data_significant %>%
    dplyr::arrange(Genus)
  
  # Now plot only the strong significant correlations
  ggplot(plot_data_significant, aes(x = Nutrient, y = Genus, fill = Correlation)) +
    # Explicitly set y-axis order
    scale_y_discrete(limits = sorted_genera) +
    
    # Visible borders around tiles
    geom_tile(color = "black", size = 0.3) +
    
    # Two-color gradient (blue-white-red) for correlation
    scale_fill_gradient2(
      low = "#241571",
      mid = "white",
      high = "#BF0A30",
      midpoint = 0,
      limit = c(-1, 1),
      name = "Spearman\nCorrelation"
    ) +
    
    # Add correlation values and significance
    geom_text(aes(label = Label), color = "white", size = 3) +
    
    labs(title = title) +
    
    # Use theme_bw to get visible grid lines/borders
    theme_bw(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      axis.title         = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 16)
    )
}

########################################
# 6) Generate Correlation Plot for Top 25% Genera
########################################

# Extract just the data for the top 25% genera
genus_data <- as.data.frame(all_phyla)

# Calculate correlations for top 25% genera
top_correlations <- correlation_with_pvalues(genus_data, nutrient_data, method = "spearman")

# Create the correlation plot with Strong correlations only
top_genera_correlation_plot <- create_correlation_plot(
  top_correlations$correlation, 
  top_correlations$pvalues, 
  ""
)

# Display the plot
print(top_genera_correlation_plot)

# Save plots
ggsave("RDA_plots/Figure5s_Spearman_Corelation_05.jpeg", top_genera_correlation_plot, device = "jpeg", dpi = 600, 
       width = 14, height = 26, units = "in")
ggsave("RDA_plots/Figure5s_Spearman_Corelation_05.pdf", top_genera_correlation_plot, device = "pdf", dpi = 600,
       width = 14, height = 26, units = "in")
ggsave("RDA_plots/Figure5s_Spearman_Corelation_05.tiff", top_genera_correlation_plot, device = "tiff", dpi = 600,
       width = 14, height = 26, units = "in")


# Save plots
# Save plots
ggsave("RDA_plots/Supplementary_Figure_S8.jpeg", top_genera_correlation_plot, device = "jpeg", dpi = 600, 
       width = 14, height = 26, units = "in")
ggsave("RDA_plots/Supplementary_Figure_S8.pdf", top_genera_correlation_plot, device = "pdf", dpi = 600,
       width = 14, height = 26, units = "in")
ggsave("RDA_plots/Supplementary_Figure_S8.tiff", top_genera_correlation_plot, device = "tiff", dpi = 600,
       width = 14, height = 26, units = "in")
