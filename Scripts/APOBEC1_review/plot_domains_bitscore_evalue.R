# Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")

# 1. Load the data
# Assuming the file is named 'for_plotting_domain_full_hits.tsv'
data <- read.table("for_plotting_domain_full_hits.tsv",
                   header = TRUE,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   quote = "")

# 2. Preprocess the data
plot_data <- data %>%
  # Filter out 'superfamily' hits to focus on specific/non-specific hits
  filter(Hit_type != "superfamily") %>%
  mutate(
    # Remove underscores for clean species names
    Species_Clean = str_replace_all(Query, "_", " "),
    # Calculate -log10(E-value) for better visual scaling
    neg_log10_evalue = -log10(E_Value),
    # Handle E_Value of 0 to avoid Infinity
    neg_log10_evalue = ifelse(is.infinite(neg_log10_evalue), 
                              max(neg_log10_evalue[is.finite(neg_log10_evalue)], na.rm = TRUE) + 10, 
                              neg_log10_evalue)
  )

# 3. Calculate Average Statistics and create formatted X-axis labels
# We use symbols: '*' for Bitscore and 'ε' (or 'e') for E-value
domain_labels <- plot_data %>%
  group_by(Short_name) %>%
  summarise(
    avg_bit = mean(Bitscore, na.rm = TRUE),
    avg_e = mean(E_Value, na.rm = TRUE)
  ) %>%
  arrange(desc(avg_bit)) %>%
  mutate(
    # Custom label format: DomainName \n *Bitscore  ~E-value
    custom_label = paste0(
      Short_name, 
      "\n*", round(avg_bit, 1), 
      "  ~", formatC(avg_e, format = "e", digits = 1)
    )
  )

# Merge the labels back and set order based on Bitscore significance
plot_data <- plot_data %>%
  left_join(domain_labels, by = "Short_name") %>%
  mutate(custom_label = factor(custom_label, levels = domain_labels$custom_label))

# 4. Create the plot
p <- ggplot(plot_data, aes(x = custom_label, y = Species_Clean)) +
  # Bubble points: Size = Bitscore, Color = -log10(E-value)
  geom_point(aes(size = Bitscore, color = neg_log10_evalue), alpha = 0.8) +
  
  # Color gradient with high variation (Blue -> Green -> Yellow -> Red)
  scale_color_gradientn(
    colors = c("blue", "cyan", "green", "yellow", "orange", "red"), 
    name = "-log10(E-value)"
  ) +
  
  # Adjust bubble sizes
  scale_size_continuous(range = c(2, 10), name = "Bitscore") +
  
  theme_bw() +
  labs(
    title = "Conserved Domain Hits for Avian APOBEC1",
    subtitle = "X-axis: Domain (*Avg Bitscore, ~Avg E-value)",
    x = "Conserved Domain (Summary Stats Across Species)",
    y = "Species",
    caption = "Data: NCBI Batch CD-Search | Symbol Key: * = Average Bitscore, ~ = Average E-Value"
  ) +
  theme(
    # Bold and angled X-axis labels for domain summaries
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 9),
    # Italics for species names
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey95")
  )

# 5. Save and display the plot
print(p)
ggsave("APOBEC1_domain_hits_refined.png", width = 12, height = 15, dpi = 300)
