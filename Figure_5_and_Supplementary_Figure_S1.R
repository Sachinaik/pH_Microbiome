library(cowplot)
library(stringr)
library(ggtext)  # For rich text styling in ggplot2
library(ggplot2)
library(FSA)  # For Dunn's test
library(purrr)
library(dplyr)      # For data manipulation
library(purrr)      # For map() function
library(tidyr)      # For nest() and unnest()
library(broom) 

###. PHYLUM
# Step 1: Transform counts to relative abundance
Bacteria_ps_filtered <- readRDS("Bacteria_ps_filtered_nochlormito.rds")
Bacteria_ps_relabund <- transform_sample_counts(Bacteria_ps_filtered, function(x) x / sum(x))
composite <- psmelt(Bacteria_ps_relabund)

# Step 3: Filter data for only the top 100% phyla across all pH levels
#i am using all phylum
composite_top5 <- composite

# Step 4: Perform Kruskal-Wallis test for each Phylum
sig_phyla <- composite_top5 %>%
  nest(data = -Phylum) %>%
  mutate(test = map(data, ~ kruskal.test(Abundance ~ Treatments, data = .x) %>% tidy())) %>%
  unnest(test) %>%
  filter(p.value < 0.001) %>%
  select(Phylum, p.value)

# Step 5: Merge significant phyla results with composite data
significant_phyla_data <- composite_top5 %>%
  inner_join(sig_phyla, by = "Phylum") %>%
  mutate(Relative_Abundance = Abundance * 100)

write.csv(significant_phyla_data, "Signifcant_taxa/significant_phyla_dataall_Phylum(p0.01).csv", row.names = FALSE)

# Step 6: Perform Dunn's Test for pairwise comparisons for each significant Phylum
pairwise_results <- significant_phyla_data %>%
  nest(data = -Phylum) %>%
  mutate(pairwise_test = map(data, ~ {
    test_result <- tryCatch({
      dunnTest(Abundance ~ Treatments, data = .x, method = "bonferroni")$res
    }, error = function(e) NULL)  # Return NULL if error occurs
    test_result
  })) %>%
  unnest(pairwise_test) %>%
  filter(!is.na(P.adj) & P.adj < 0.01) %>%
  select(Phylum, Comparison, Z, P.adj)
write.csv(pairwise_results, "Signifcant_taxa/pairwise_resultsphylum(p0.01).csv", row.names = FALSE)

pairwise_results <- pairwise_results[abs(pairwise_results$Z) >= 5, ]


# Change format of Comparison labels (E4.5 - E5.5 to E4.5_E5.5)
# Change format of Comparison labels using base R instead of stringr
pairwise_results <- pairwise_results %>%
  mutate(Comparison = gsub(" - ", "/", Comparison))

# Adding a column to mark significant comparisons with asterisks
pairwise_results <- pairwise_results %>%
  mutate(Significance = case_when(
    P.adj < 0.001 ~ "***",
    P.adj < 0.01  ~ "**",
    P.adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

# Define treatment mapping
treatment_map <- c(
  "Empty_4.5" = "E4.5", "Petunia_4.5" = "P4.5", "Tomato_4.5" = "T4.5", "Marigold_4.5" = "M4.5", "Geranium_4.5" = "G4.5",
  "Empty_5.5" = "E5.5", "Petunia_5.5" = "P5.5", "Tomato_5.5" = "T5.5", "Marigold_5.5" = "M5.5", "Geranium_5.5" = "G5.5",
  "Empty_6.2" = "E6.2", "Petunia_6.2" = "P6.2", "Tomato_6.2" = "T6.2", "Marigold_6.2" = "M6.2", "Geranium_6.2" = "G6.2",
  "Empty_7" = "E7.0", "Petunia_7" = "P7.0", "Tomato_7" = "T7.0", "Marigold_7" = "M7.0", "Geranium_7" = "G7.0"
)

# Apply the mapping to the Comparison column using base R
pairwise_results$Comparison <- as.character(pairwise_results$Comparison)

# Loop through each treatment pattern and replace it
for (pattern in names(treatment_map)) {
  pairwise_results$Comparison <- gsub(pattern, treatment_map[pattern], pairwise_results$Comparison)
}

# View the updated dataframe
pairwise_results

# Define colors for each prefix
color_map <- c("E" = "black", "T" = "brown", "P" = "darkgreen", "G" = "goldenrod", "M" = "cyan3")

# Function to style each label part based on the prefix using ggtext-compatible HTML tags
style_comparison_label <- function(label) {
  parts <- strsplit(label, "/")[[1]]
  styled_parts <- sapply(parts, function(part) {
    prefix <- substr(trimws(part), 1, 1)
    color <- color_map[prefix]
    sprintf("<span style='color:%s;'>%s</span>", color, trimws(part))
  })
  paste(styled_parts, collapse = "/")
}

# Apply this styling function to the Comparison column
pairwise_results <- pairwise_results %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label))


# 1. Remove all comparisons involving Empty pots
pairwise_results_filtered <- pairwise_results 

# 2. Extract actual comparisons that exist in your data
actual_comparisons <- unique(as.character(pairwise_results_filtered$Comparison))

# 3. Create a custom ordering for comparisons grouped by plant species
plant_types <- c("E", "P", "M", "T","G")
pH_values <- c("4.5", "5.5", "6.2", "7.0")

# Create all possible pairwise comparisons within each plant type
ordered_comparisons <- c()
for(plant in plant_types) {
  for(i in 1:(length(pH_values)-1)) {
    for(j in (i+1):length(pH_values)) {
      # Create both possible orders of the comparison
      forward <- paste0(plant, pH_values[i], "/", plant, pH_values[j])
      reverse <- paste0(plant, pH_values[j], "/", plant, pH_values[i])
      
      # Add the comparison that exists in the data (if any)
      if(forward %in% actual_comparisons) {
        ordered_comparisons <- c(ordered_comparisons, forward)
      } else if(reverse %in% actual_comparisons) {
        ordered_comparisons <- c(ordered_comparisons, reverse)
      }
    }
  }
}

# 4. Filter to include only the within-species comparisons we want
pairwise_results_species_grouped <- pairwise_results_filtered %>%
  filter(Comparison %in% ordered_comparisons)
write.csv(pairwise_results_species_grouped, "Signifcant_taxa/pairwise_results_species_grouped_Phylum_All_z5(p0.01).csv", row.names = FALSE)


# 5. Apply styling to the Comparison column
pairwise_results_species_grouped <- pairwise_results_species_grouped %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label)) %>%
  mutate(Comparison = factor(Comparison, levels = ordered_comparisons))

# 6. Split for plotting if needed
split_index <- ceiling(length(ordered_comparisons) / 2)
comparisons_part1 <- ordered_comparisons[1:split_index]
comparisons_part2 <- ordered_comparisons[(split_index+1):length(ordered_comparisons)]

# 7. Filter data for each part
pairwise_results_part1 <- pairwise_results_species_grouped %>%
  filter(Comparison %in% comparisons_part1)
pairwise_results_part2 <- pairwise_results_species_grouped %>%
  filter(Comparison %in% comparisons_part2)

# 8. Create the heatmap plots
heatmap_plotP1 <- ggplot(pairwise_results_species_grouped, aes(x = StyledComparison, y = Phylum, fill = Z)) +
  geom_tile(color = "white", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.6, hjust = 1, size = 17.5),
    axis.text.y = element_text(size = 20, face = "italic", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_discrete(labels = setNames(pairwise_results_part1$StyledComparison, pairwise_results_part1$Comparison))

#####pH grouped

pairwise_results_filtered <- pairwise_results 

# 2. Extract actual comparisons that exist in your data
actual_comparisons <- unique(as.character(pairwise_results_filtered$Comparison))

# 3. Create a custom ordering for comparisons grouped by pH
pH_values <- c("4.5", "5.5", "6.2", "7.0")
plant_types <- c("E", "P", "M", "T","G")

# Create all possible pairwise comparisons within each pH value
ordered_comparisons <- c()
for(pH in pH_values) {
  # Store comparisons that exist for this pH
  pH_comparisons <- c()
  
  for(i in 1:(length(plant_types)-1)) {
    for(j in (i+1):length(plant_types)) {
      # Check both possible orders of the comparison
      forward <- paste0(plant_types[i], pH, "/", plant_types[j], pH)
      reverse <- paste0(plant_types[j], pH, "/", plant_types[i], pH)
      
      if(forward %in% actual_comparisons) {
        pH_comparisons <- c(pH_comparisons, forward)
      } else if(reverse %in% actual_comparisons) {
        pH_comparisons <- c(pH_comparisons, reverse)
      }
    }
  }
  
  # Add this pH's comparisons to the ordered list
  ordered_comparisons <- c(ordered_comparisons, pH_comparisons)
}

# 4. Filter to include only the pH-grouped comparisons
pairwise_results_pH_grouped <- pairwise_results_filtered %>%
  filter(Comparison %in% ordered_comparisons)

# 5. Apply styling to the Comparison column
pairwise_results_pH_grouped <- pairwise_results_pH_grouped %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label)) %>%
  mutate(Comparison = factor(Comparison, levels = ordered_comparisons))

# 6. Split for plotting if needed
split_index <- ceiling(length(ordered_comparisons) / 2)
comparisons_part1 <- ordered_comparisons[1:split_index]
comparisons_part2 <- ordered_comparisons[(split_index+1):length(ordered_comparisons)]

# 7. Filter data for each part
pairwise_results_part1 <- pairwise_results_pH_grouped %>%
  filter(Comparison %in% comparisons_part1)
pairwise_results_part2 <- pairwise_results_pH_grouped %>%
  filter(Comparison %in% comparisons_part2)

# 8. Create the heatmap plots using StyledComparison for colored x-axis text
heatmap_plotP2 <- ggplot(pairwise_results_pH_grouped, aes(x = StyledComparison, y = Phylum, fill = Z)) +
  geom_tile(color = "white", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.6, hjust = 1, size = 17.5),
    axis.text.y = element_text(size = 20, face = "italic", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_discrete(labels = setNames(pairwise_results_part1$StyledComparison, pairwise_results_part1$Comparison))


########genus top 25%
library(cowplot)
library(stringr)
library(ggtext)  # For rich text styling in ggplot2
library(ggplot2)
library(FSA)  # For Dunn's test
library(purrr)
library(dplyr)      # For data manipulation
library(purrr)      # For map() function
library(tidyr)      # For nest() and unnest()
library(broom) 

###. Genus
# Step 1: Transform counts to relative abundance
Bacteria_ps_filtered <- readRDS("Bacteria_ps_filtered_nochlormito.rds")
Bacteria_ps_relabund <- transform_sample_counts(Bacteria_ps_filtered, function(x) x / sum(x))
composite <- psmelt(Bacteria_ps_relabund)


# Step 2: Filter for the top 25% genera
get_top_25_percent_Genus <- function(data) {
  data %>%
    group_by(Genus) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    top_n(n = floor(0.25 * n()), wt = mean_abundance) %>%
    pull(Genus)
}

top_25_percent_Genus <- get_top_25_percent_Genus(composite)

composite_top5 <- composite %>%
  filter(Genus %in% top_25_percent_Genus)

# Step 2.5: Shorten genus names
composite_top5 <- composite_top5 %>%
  mutate(
    Genus = case_when(
      Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" ~ "Rhizobium",
      Genus == "Burkholderia-Caballeronia-Paraburkholderia" ~ "Burkholderia",
      TRUE ~ Genus
    )
  )




# Step 3: Filter data for only the top 100% phyla across all pH levels

# Step 4: Perform Kruskal-Wallis test for each Genus
sig_phyla <- composite_top5 %>%
  nest(data = -Genus) %>%
  mutate(test = map(data, ~ kruskal.test(Abundance ~ Treatments, data = .x) %>% tidy())) %>%
  unnest(test) %>%
  filter(p.value < 0.001) %>%
  select(Genus, p.value)

# Step 5: Merge significant phyla results with composite data
significant_phyla_data <- composite_top5 %>%
  inner_join(sig_phyla, by = "Genus") %>%
  mutate(Relative_Abundance = Abundance * 100)

write.csv(significant_phyla_data, "Signifcant_taxa/significant_phyla_dataall_Genus(p0.01).csv", row.names = FALSE)

# Step 6: Perform Dunn's Test for pairwise comparisons for each significant Genus
pairwise_results <- significant_phyla_data %>%
  nest(data = -Genus) %>%
  mutate(pairwise_test = map(data, ~ {
    test_result <- tryCatch({
      dunnTest(Abundance ~ Treatments, data = .x, method = "bonferroni")$res
    }, error = function(e) NULL)  # Return NULL if error occurs
    test_result
  })) %>%
  unnest(pairwise_test) %>%
  filter(!is.na(P.adj) & P.adj < 0.01) %>%
  select(Genus, Comparison, Z, P.adj)
write.csv(pairwise_results, "Signifcant_taxa/pairwise_resultsGenus(p0.01).csv", row.names = FALSE)

pairwise_results <- pairwise_results[abs(pairwise_results$Z) >= 5, ]
write.csv(pairwise_results, "Signifcant_taxa/pairwise_resultsGenus(p0.01)_TOP25_P_Z_5.csv", row.names = FALSE)


# Change format of Comparison labels (E4.5 - E5.5 to E4.5_E5.5)
# Change format of Comparison labels using base R instead of stringr
pairwise_results <- pairwise_results %>%
  mutate(Comparison = gsub(" - ", "/", Comparison))

# Adding a column to mark significant comparisons with asterisks
pairwise_results <- pairwise_results %>%
  mutate(Significance = case_when(
    P.adj < 0.001 ~ "***",
    P.adj < 0.01  ~ "**",
    P.adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

# Define treatment mapping
treatment_map <- c(
  "Empty_4.5" = "E4.5", "Petunia_4.5" = "P4.5", "Tomato_4.5" = "T4.5", "Marigold_4.5" = "M4.5", "Geranium_4.5" = "G4.5",
  "Empty_5.5" = "E5.5", "Petunia_5.5" = "P5.5", "Tomato_5.5" = "T5.5", "Marigold_5.5" = "M5.5", "Geranium_5.5" = "G5.5",
  "Empty_6.2" = "E6.2", "Petunia_6.2" = "P6.2", "Tomato_6.2" = "T6.2", "Marigold_6.2" = "M6.2", "Geranium_6.2" = "G6.2",
  "Empty_7" = "E7.0", "Petunia_7" = "P7.0", "Tomato_7" = "T7.0", "Marigold_7" = "M7.0", "Geranium_7" = "G7.0"
)

# Apply the mapping to the Comparison column using base R
pairwise_results$Comparison <- as.character(pairwise_results$Comparison)

# Loop through each treatment pattern and replace it
for (pattern in names(treatment_map)) {
  pairwise_results$Comparison <- gsub(pattern, treatment_map[pattern], pairwise_results$Comparison)
}

# View the updated dataframe
pairwise_results

exclude_phyla <- c(NA)
#"Firmicutes"
pairwise_results <- pairwise_results %>%
  filter(!Genus %in% exclude_phyla)


# Define colors for each prefix
color_map <- c("E" = "black", "T" = "brown", "P" = "darkgreen", "G" = "goldenrod", "M" = "cyan3")

# Function to style each label part based on the prefix using ggtext-compatible HTML tags
style_comparison_label <- function(label) {
  parts <- strsplit(label, "/")[[1]]
  styled_parts <- sapply(parts, function(part) {
    prefix <- substr(trimws(part), 1, 1)
    color <- color_map[prefix]
    sprintf("<span style='color:%s;'>%s</span>", color, trimws(part))
  })
  paste(styled_parts, collapse = "/")
}

# Apply this styling function to the Comparison column
pairwise_results <- pairwise_results %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label))


# 1. Remove all comparisons involving Empty pots
pairwise_results_filtered <- pairwise_results 

# 2. Extract actual comparisons that exist in your data
actual_comparisons <- unique(as.character(pairwise_results_filtered$Comparison))

# 3. Create a custom ordering for comparisons grouped by plant species
plant_types <- c("E", "P", "M", "T","G")
pH_values <- c("4.5", "5.5", "6.2", "7.0")

# Create all possible pairwise comparisons within each plant type
ordered_comparisons <- c()
for(plant in plant_types) {
  for(i in 1:(length(pH_values)-1)) {
    for(j in (i+1):length(pH_values)) {
      # Create both possible orders of the comparison
      forward <- paste0(plant, pH_values[i], "/", plant, pH_values[j])
      reverse <- paste0(plant, pH_values[j], "/", plant, pH_values[i])
      
      # Add the comparison that exists in the data (if any)
      if(forward %in% actual_comparisons) {
        ordered_comparisons <- c(ordered_comparisons, forward)
      } else if(reverse %in% actual_comparisons) {
        ordered_comparisons <- c(ordered_comparisons, reverse)
      }
    }
  }
}

# 4. Filter to include only the within-species comparisons we want
pairwise_results_species_grouped <- pairwise_results_filtered %>%
  filter(Comparison %in% ordered_comparisons)
write.csv(pairwise_results_species_grouped, "Signifcant_taxa/pairwise_results_species_grouped_Genus_All_z5(p0.01).csv", row.names = FALSE)


# 5. Apply styling to the Comparison column
pairwise_results_species_grouped <- pairwise_results_species_grouped %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label)) %>%
  mutate(Comparison = factor(Comparison, levels = ordered_comparisons))

# 6. Split for plotting if needed
split_index <- ceiling(length(ordered_comparisons) / 2)
comparisons_part1 <- ordered_comparisons[1:split_index]
comparisons_part2 <- ordered_comparisons[(split_index+1):length(ordered_comparisons)]

# 7. Filter data for each part
pairwise_results_part1 <- pairwise_results_species_grouped %>%
  filter(Comparison %in% comparisons_part1)
pairwise_results_part2 <- pairwise_results_species_grouped %>%
  filter(Comparison %in% comparisons_part2)

# 8. Create the heatmap plots
heatmap_plotG1 <- ggplot(pairwise_results_species_grouped, aes(x = StyledComparison, y = Genus, fill = Z)) +
  geom_tile(color = "white", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.6, hjust = 1, size = 17.5),
    axis.text.y = element_text(size = 20, face = "italic", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_discrete(labels = setNames(pairwise_results_part1$StyledComparison, pairwise_results_part1$Comparison))

heatmap_plotG1

###genus pH grouped

# 1. Remove all comparisons involving Empty pots
pairwise_results_filtered <- pairwise_results 

# 2. Extract actual comparisons that exist in your data
actual_comparisons <- unique(as.character(pairwise_results_filtered$Comparison))

# 3. Create a custom ordering for comparisons grouped by pH
pH_values <- c("4.5", "5.5", "6.2", "7.0")
plant_types <- c("E", "P", "M", "T","G")

# Create all possible pairwise comparisons within each pH value
ordered_comparisons <- c()
for(pH in pH_values) {
  # Store comparisons that exist for this pH
  pH_comparisons <- c()
  
  for(i in 1:(length(plant_types)-1)) {
    for(j in (i+1):length(plant_types)) {
      # Check both possible orders of the comparison
      forward <- paste0(plant_types[i], pH, "/", plant_types[j], pH)
      reverse <- paste0(plant_types[j], pH, "/", plant_types[i], pH)
      
      if(forward %in% actual_comparisons) {
        pH_comparisons <- c(pH_comparisons, forward)
      } else if(reverse %in% actual_comparisons) {
        pH_comparisons <- c(pH_comparisons, reverse)
      }
    }
  }
  
  # Add this pH's comparisons to the ordered list
  ordered_comparisons <- c(ordered_comparisons, pH_comparisons)
}

# 4. Filter to include only the pH-grouped comparisons
pairwise_results_pH_grouped <- pairwise_results_filtered %>%
  filter(Comparison %in% ordered_comparisons)

# 5. Apply styling to the Comparison column
pairwise_results_pH_grouped <- pairwise_results_pH_grouped %>%
  mutate(StyledComparison = sapply(Comparison, style_comparison_label)) %>%
  mutate(Comparison = factor(Comparison, levels = ordered_comparisons))

# 6. Split for plotting if needed
split_index <- ceiling(length(ordered_comparisons) / 2)
comparisons_part1 <- ordered_comparisons[1:split_index]
comparisons_part2 <- ordered_comparisons[(split_index+1):length(ordered_comparisons)]

# 7. Filter data for each part
pairwise_results_part1 <- pairwise_results_pH_grouped %>%
  filter(Comparison %in% comparisons_part1)
pairwise_results_part2 <- pairwise_results_pH_grouped %>%
  filter(Comparison %in% comparisons_part2)

# 8. Create the heatmap plots using StyledComparison for colored x-axis text
heatmap_plotG2 <- ggplot(pairwise_results_pH_grouped, aes(x = StyledComparison, y = Genus, fill = Z)) +
  geom_tile(color = "white", size = 0.05) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-score") +
  theme_minimal() +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.6, hjust = 1, size = 17.5),
    axis.text.y = element_text(size = 20, face = "italic", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_discrete(labels = setNames(pairwise_results_part1$StyledComparison, pairwise_results_part1$Comparison))

heatmap_plotG2


# Combine plots for Figure 5
combined_heatmapsPhylum <- plot_grid(heatmap_plotP1, heatmap_plotG1, ncol = 1, rel_heights = c(1,1.65), align = "v", labels = c("A", "B"), label_size = 20)
combined_heatmapsPhylum

# Combine plots for Figure 5s
combined_heatmapsPhylum2 <- plot_grid(heatmap_plotP2, heatmap_plotG2, ncol = 1, rel_heights = c(1.25,1), align = "v", labels = c("A", "B"), label_size = 20)
combined_heatmapsPhylum2

heatmap_plotP2 <- heatmap_plotP2 + theme(plot.margin = margin(20, 10, 10, 30))  # top, right, bottom, left
heatmap_plotG2 <- heatmap_plotG2 + theme(plot.margin = margin(20, 10, 10, 30))

# Combine plots with adjusted label positions
combined_heatmapsPhylum2 <- plot_grid(
  heatmap_plotP2, 
  heatmap_plotG2, 
  ncol = 1, 
  rel_heights = c(1.35, 1), 
  align = "v", 
  labels = c("A", "B"), 
  label_size = 20,
  label_x = 0.02,  # Move labels slightly to the right
  hjust = 0.45,     # Adjust horizontal justification to move labels away from y-axis
  label_y = 1.0 # Move labels slightly to the right
)

# Display the combined plot
combined_heatmapsPhylum2


# Save plot
ggsave("Figures/Fig5_pH_grouped_heatmap_Phylumallpz5_Genustop25.jpeg", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)
ggsave("Figures/Fig5_pH_grouped_heatmap_Phylumallpz5_Genustop25.pdf", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)
ggsave("Figures/Fig5_pH_grouped_heatmap_Phylumallpz5_Genustop25.tiff", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)

# Save plot
ggsave("Figures/Fig5s_species_grouped_heatmap_Phylumallpz5_Genustop25.jpeg", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)
ggsave("Figures/Fig5s_species_grouped_heatmap_Phylumallpz5_Genustop25.pdf", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)
ggsave("Figures/Fig5s_species_grouped_heatmap_Phylumallpz5_Genustop25.tiff", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)


# Save plot
ggsave("Figures/Figure_5.jpeg", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)
ggsave("Figures/Figure_5.pdf", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)
ggsave("Figures/Figure_5.tiff", plot = combined_heatmapsPhylum, width = 16, height = 18, dpi = 600)

# Save plot
ggsave("Figures/Figure_S1.jpeg", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)
ggsave("Figures/Figure_S1.pdf", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)
ggsave("Figures/Figure_S1.tiff", plot = combined_heatmapsPhylum2, width = 8, height = 10, dpi = 600)


