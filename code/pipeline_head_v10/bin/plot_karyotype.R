#!/usr/bin/env Rscript
# =============================================================================
# plot_karyotype.R — TE density karyoplot (chromosome painting), Nature-styled.
#
# Published next to the figure so you can adjust the number of scaffolds shown,
# the colour ramp, window size, etc. and re-run:
#
#   Rscript plot_karyotype.R <density.tsv> <lengths.tsv> <prefix> <formats> <top_n>
#
#   <density.tsv>  columns: chrom  start  end  frac   (frac = masked fraction of window)
#   <lengths.tsv>  columns: chrom  length
#   <prefix>       output basename (e.g. sampleX.karyoplot)
#   <formats>      comma list png,svg,pdf  (default png,svg)
#   <top_n>        number of longest scaffolds to draw (default 30)
#
# Only base R + ggplot2 are required (ragg/svglite used if present).
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
den_f   <- if (length(args) >= 1) args[1] else "density.tsv"
len_f   <- if (length(args) >= 2) args[2] else "lengths.tsv"
prefix  <- if (length(args) >= 3) args[3] else "karyoplot"
formats <- if (length(args) >= 4) trimws(strsplit(args[4], ",")[[1]]) else c("png", "svg")
top_n   <- if (length(args) >= 5) as.integer(args[5]) else 30L

theme_nature <- function(base_size = 7, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid       = element_blank(),
      axis.text        = element_text(colour = "black", size = base_size - 1),
      axis.text.y      = element_text(size = base_size - 2),
      axis.ticks.x     = element_line(linewidth = 0.3, colour = "black"),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(6, "mm"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1),
      legend.position  = "bottom",
      plot.title       = element_text(size = base_size, face = "bold")
    )
}

save_plot <- function(p, base, formats, width = 120, height = 90, units = "mm", dpi = 600) {
  for (fmt in formats) {
    f <- paste0(base, ".", fmt)
    tryCatch({
      if (fmt == "png") {
        if (requireNamespace("ragg", quietly = TRUE))
          ragg::agg_png(f, width = width, height = height, units = units, res = dpi)
        else png(f, width = width, height = height, units = units, res = dpi, type = "cairo")
      } else if (fmt == "svg") {
        if (requireNamespace("svglite", quietly = TRUE))
          svglite::svglite(f, width = width / 25.4, height = height / 25.4)
        else svg(f, width = width / 25.4, height = height / 25.4)
      } else if (fmt == "pdf") {
        grDevices::cairo_pdf(f, width = width / 25.4, height = height / 25.4)
      } else stop("unknown format")
      print(p); dev.off()
    }, error = function(e) message("WARN: could not write ", f, ": ", conditionMessage(e)))
  }
}

empty_plot <- function(msg) ggplot() + annotate("text", 0, 0, label = msg, size = 3) + theme_void()

ok <- file.exists(den_f) && file.info(den_f)$size > 0 &&
      file.exists(len_f) && file.info(len_f)$size > 0
if (!ok) { save_plot(empty_plot("No annotation to plot"), prefix, formats); quit(save = "no", status = 0) }

den <- read.table(den_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
len <- read.table(len_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
if (nrow(den) == 0 || nrow(len) == 0) {
  save_plot(empty_plot("No annotation to plot"), prefix, formats); quit(save = "no", status = 0)
}

# Keep the top_n longest scaffolds, ordered longest → shortest (top of the plot).
len   <- len[order(-len$length), ]
keep  <- head(len$chrom, top_n)
den   <- den[den$chrom %in% keep, ]
den$chrom <- factor(den$chrom, levels = rev(keep))
den$mid_mb <- (den$start + den$end) / 2 / 1e6

p <- ggplot(den, aes(x = mid_mb, y = chrom, fill = 100 * frac)) +
  geom_tile(height = 0.75) +
  scale_fill_gradient(low = "#F0F0F0", high = "#31548F",
                      name = "TE density (% of window)",
                      limits = c(0, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Position (Mb)", y = NULL) +
  theme_nature() +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

# Height grows with the number of scaffolds so rows stay legible.
h <- max(45, 8 + 3.2 * length(keep))
save_plot(p, prefix, formats, width = 130, height = h)
