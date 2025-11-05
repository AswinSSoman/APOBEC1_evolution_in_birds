# Load required package
library(ggplot2)

# Read input file
data <- read.table("ERV_DNA_edit_sites_comparison", header = TRUE)

# Create a new simplified color category
data$Group <- ifelse(data$Order == "Passeriformes", "Passeriformes", "Others")

# Convert Loss_info to factor for shape mapping
data$Loss_info <- factor(data$Loss_info, levels = c(0, 1), labels = c("Present", "Lost"))

# Create scatter plot
p <- ggplot(data, aes(x = GA_edit_sites, y = ERV_content,
                      color = Group, shape = Loss_info)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  scale_color_manual(values = c("Passeriformes" = "#d95f02", "Others" = "black")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Number of APOBEC1-associated G-A Edit Sites",
    y = "ERV Content (bp)",
    color = "Order Group",
    shape = "APOBEC1 Status",
    title = "Higher GA Edit Sites are Associated with Higher ERV Content"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

# Display in R
print(p)


# Save as high-quality PDF
ggsave(
  filename = "GA_edit_vs_ERV_content.pdf",
  plot = p,
  width = 18,    # wider canvas
  height = 15,    # taller canvas
  dpi = 600,     # high-resolution for export
  device = cairo_pdf  # ensures vector quality text and shapes
)
