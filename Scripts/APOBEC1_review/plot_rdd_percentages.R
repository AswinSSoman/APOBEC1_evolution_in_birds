#!/usr/bin/env Rscript

# =========================================================
# Creative RDD summary plot (Refined Labels)
# =========================================================

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1){
    stop("Usage: Rscript plot_rdd_summary.R <input_file>")
}

input_file <- args[1]

# Read table
df <- read.table(input_file,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE,
                 check.names = FALSE)

# RDD columns
rdd_cols <- c("ADAR", "APOBEC1",
              "A>C", "A>T", "C>A", "C>G",
              "G>C", "G>T", "T>A", "T>G")

# Colors
cols <- c(
    "ADAR"    = "#1565C0",   # blue
    "APOBEC1" = "#EF6C00",   # orange
    "Others"  = "#7B1FA2"    # purple
)

# Statistics table
stats <- data.frame(
    RDD_Type = rdd_cols,
    Min = NA, Q1 = NA, Median = NA, Mean = NA, Q3 = NA, Max = NA, Color = NA
)

# Calculate stats
for(i in seq_along(rdd_cols)) {
    vals <- df[[ rdd_cols[i] ]]
    stats$Min[i]    <- min(vals, na.rm = TRUE)
    stats$Q1[i]     <- quantile(vals, 0.25, na.rm = TRUE)
    stats$Median[i] <- median(vals, na.rm = TRUE)
    stats$Mean[i]   <- mean(vals, na.rm = TRUE)
    stats$Q3[i]     <- quantile(vals, 0.75, na.rm = TRUE)
    stats$Max[i]    <- max(vals, na.rm = TRUE)

    if(stats$RDD_Type[i] == "ADAR"){
        stats$Color[i] <- cols["ADAR"]
    } else if(stats$RDD_Type[i] == "APOBEC1"){
        stats$Color[i] <- cols["APOBEC1"]
    } else {
        stats$Color[i] <- cols["Others"]
    }
}

# PDF
pdf("rdd_percentage_summary_plot.pdf", width = 8.5, height = 6)

par(mar = c(8,5,4,2), bg = "white")

# Increase ylim slightly more to accommodate labels at the top
plot(1:nrow(stats),
     stats$Mean,
     type = "n",
     ylim = c(0, max(stats$Max) + 20), 
     xaxt = "n",
     xlab = "",
     ylab = "Percentage (%)",
     main = "RDD Percentage Distribution Across Samples",
     cex.main = 1.6, cex.lab = 1.4, font.lab = 2, font.main = 2)

# Grid
abline(h = pretty(c(0, max(stats$Max))), col = "grey90", lty = "dashed")

# X-axis labels
axis(1, at = 1:nrow(stats), labels = FALSE)
text(1:nrow(stats), par("usr")[3] - 3, labels = stats$RDD_Type,
     srt = 45, adj = 1, xpd = TRUE, cex = 1.2, font = 2, col = stats$Color)

# Draw each category
for(i in 1:nrow(stats)) {
    x <- i
    colx <- stats$Color[i]

    # Min-Max whisker
    arrows(x, stats$Min[i], x, stats$Max[i],
           angle = 90, code = 3, length = 0.05, lwd = 2, col = colx)

    # Q1-Q3 thick bar
    segments(x, stats$Q1[i], x, stats$Q3[i],
             lwd = 14, col = adjustcolor(colx, alpha.f = 0.75))

    # Mean point
    points(x, stats$Mean[i], pch = 21, bg = "white", col = colx, lwd = 2, cex = 2)

    # Median point
    points(x, stats$Median[i], pch = 18, col = "black", cex = 1.3)

    # --- LABEL LOGIC ---
    # Set label height relative to the Max whisker to avoid overlap
    label_y <- stats$Max[i] + 8 

    # Dotted connector line: from Max whisker to just below the text
    segments(x, stats$Max[i] + 1, 
             x, label_y - 2, 
             lty = "dotted", lwd = 1.2, col = "grey40")

    # Mean label: Smaller (cex=1.0), not bold (font=1)
    text(x, label_y,
         labels = sprintf("%.2f", stats$Mean[i]),
         cex = 1.0, 
         font = 1, 
         col = "black")
}

# Legend
legend("topright",
       legend = c("ADAR", "APOBEC1", "Other RDD types"),
       fill = c(cols["ADAR"], cols["APOBEC1"], cols["Others"]),
       border = NA, bty = "n", cex = 1.1)

dev.off()
cat("\nOutput written to: rdd_percentage_summary_plot.pdf\n")
