# Load necessary packages
library(ggplot2)
library(FSA)  # for Dunn's test
library(dunn.test)
library(cowplot)  # For combining the plots

# Convert 'Treatment' to a factor explicitly (assuming pH values are in 'Treatment')
data$Treatment <- as.factor(data$Treatment)
data$pH <- as.factor(data$pH)

# Define variables for each category
macronutrients <- c("K", "Ca", "Mg", "P",)
macronutrients1 <- c("S", "TN", "Mg", "NPOC")
micronutrients <- c("Fe", "Mn", "Cu", "Zn")
micronutrients1 <- c("B", "Mo", "Al","Na")
npoc_tn <- c("NPOC", "TN", "K", "Zn")
anions <- c("Flouride_CD1", "Chloride_CD1", "Bromide_CD1", "Nitrate_CD1" )
anions1 <- c("Sulfate_CD1", "Phosphate_CD1", "Bromide_CD1", "Nitrate_CD1")
cations <- c("Calcium_CD2", "Magnesium_CD2", "Potassium_CD2", "Sodium_CD2")

# Custom color scale for pH
custom_colors <- c("4.5" = "#E69F00", "5.5" = "#56B4E9", "6.2" = "#009E73", "7" = "#d62728")

# Function to generate plots for a given group of variables
generate_plot <- function(variables) {
  plots <- list()  # Store plots here
  
  for (var in variables) {
    # Kruskal-Wallis test
    kruskal_test <- kruskal.test(as.formula(paste(var, "~ Treatment")), data = data)
    
    # If Kruskal-Wallis is significant, proceed with Dunn's test
    if (kruskal_test$p.value < 0.05) {
      # Perform Dunn's test
      dunn_results <- dunnTest(as.formula(paste(var, "~ Treatment")), data = data, method = "bonferroni")
      
      # Access the correct results from the 'res' data frame
      dunn_results_df <- dunn_results$res
      
      # Filter to keep only significant results (p.adj < 0.05)
      significant_results <- dunn_results_df[dunn_results_df$P.adj < 0.05, ]
      
      # Add a new column with the variable name
      significant_results$variable <- var
      
      # Plot for each variable with treatment groups
      if (nrow(significant_results) > 0) {
        # Ensure the data for each plot is correctly mapped using the 'get(var)' function
        plot <- ggplot(data, aes(x = Treatment, y = .data[[var]])) +
          geom_boxplot(aes(fill = pH), alpha = 0.5, color = "black") +  # Use Treatment for fill
          geom_jitter(width = 0.1, size = 2, aes(color = pH)) +  # Use Treatment for color
          scale_fill_manual(values = custom_colors) +  # Apply custom color scale for pH
          scale_color_manual(values = custom_colors) +  # Apply custom color scale for jitter points
          labs(title = paste("", var), x = "Treatment", y = paste(var, "(ppm)")) +  # Add PPM to the y-axis title
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20, color = "black"),  # Set size and color of x-axis text
            axis.text.y = element_text(size = 20, color = "black"),  # Set size and color of y-axis text
            axis.title.x = element_text(size = 22, color = "black"),  # Set size and color of x-axis title
            axis.title.y = element_text(size = 22, color = "black"),  # Set size and color of y-axis title
            axis.ticks = element_line(color = "black", size = 1),  # Add ticks to the axes
            panel.grid = element_blank(),  # Remove grid
            panel.background = element_blank(),  # Remove panel background
            plot.background = element_blank(),
            plot.title = element_text(size = 26),  # Remove plot background
            axis.line = element_line(color = "black"),  # Add axis lines (box around the plot)
            legend.position = "right",  # Position the legend on the right
            legend.title = element_text(size = 16),  # Increase the legend title size
            legend.text = element_text(size = 14),  # Increase the legend text size
            legend.key.size = unit(0.5, "cm")  # Increase the size of the legend keys
          )
        
        plots[[var]] <- plot  # Save the plot
      }
    }
  }
  
  # Combine all the individual plots using cowplot
  combined_plot <- plot_grid(plotlist = plots, ncol = 2, align = "v", label_size = 20)
  
  # Save the combined plot
  ggsave(paste0("combined_plot_", paste(variables, collapse = "_"), ".jpeg"), plot = combined_plot, height = 12, width = 16, dpi = 300)
}

# Generate the cowplots for each category
generate_plot(macronutrients)
generate_plot(macronutrients1)
generate_plot(micronutrients)
generate_plot(micronutrients1)
generate_plot(npoc_tn)
generate_plot(anions)
generate_plot(anions1)
generate_plot(cations)

# Save the combined results to a CSV file
write.csv(final_results, "All_Dunn_test_results.csv", row.names = FALSE)
