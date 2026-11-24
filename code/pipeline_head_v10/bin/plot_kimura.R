#!/usr/bin/env Rscript
# =============================================================================
# plot_kimura.R — TE divergence (Kimura) landscape, Nature-styled.
#
# This is the exact script the pipeline used to draw the plot. It is published
# next to the figure so you can tweak colours, sizes, labels, etc. and re-run:
#
#   Rscript plot_kimura.R <kimura_matrix.tsv> <prefix> <formats> <faceted>
#
#   <kimura_matrix.tsv>  RepeatMasker/calcDivergenceFromAlign coverage matrix
#                        (first column "Div" = Kimura %, then one bp column per class;
#                        "PercDiv" = the pipeline's fallback from the .out's plain
#                        percent divergence, labelled as such)
#   <prefix>             output basename (e.g. sampleX.kimura_landscape)
#   <formats>            comma list: png,svg,pdf   (default png,svg)
#   <faceted>            1/true = one panel per class; 0 = single stacked area
#
# Only base R + ggplot2 are required (ragg/svglite used if present).
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

args     <- commandArgs(trailingOnly = TRUE)
mat_file <- if (length(args) >= 1) args[1] else "kimura_matrix.tsv"
prefix   <- if (length(args) >= 2) args[2] else "kimura_landscape"
formats  <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("png", "svg")
faceted  <- if (length(args) >= 4) (args[4] == "1" || tolower(args[4]) == "true") else FALSE
formats  <- trimws(formats)

# ── Nature-style theme ───────────────────────────────────────────────────────
# Set base_family = "Arial" (or "Helvetica") if the font is installed on your system.
theme_nature <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size - 1),
      legend.key.size  = unit(3, "mm"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 1),
      legend.position  = "right",
      plot.title       = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 1, face = "bold")
    )
}

# Colour-blind-safe palette keyed by top-level RepeatMasker class.
te_colours <- c(
  LTR = "#4E79A7", LINE = "#F28E2B", SINE = "#E15759", DNA = "#76B7B2",
  RC = "#59A14F", Retroposon = "#9C755F", Satellite = "#EDC948",
  Simple_repeat = "#B07AA1", Low_complexity = "#FF9DA7", Unknown = "#BAB0AC"
)
fallback <- c("#8CD17D", "#B6992D", "#499894", "#D37295", "#79706E", "#D7B5A6")

save_plot <- function(p, base, formats, width = 120, height = 75, units = "mm", dpi = 600) {
  for (fmt in formats) {
    f  <- paste0(base, ".", fmt)
    ok <- tryCatch({
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
      print(p); dev.off(); TRUE
    }, error = function(e) { message("WARN: could not write ", f, ": ", conditionMessage(e)); FALSE })
  }
}

# ── Load + reshape ───────────────────────────────────────────────────────────
empty_plot <- function(msg) ggplot() + annotate("text", 0, 0, label = msg, size = 3) +
  theme_void()

ok <- file.exists(mat_file) && file.info(mat_file)$size > 0
m  <- if (ok) tryCatch(read.table(mat_file, header = TRUE, check.names = FALSE,
                                  comment.char = ""), error = function(e) NULL) else NULL

if (is.null(m) || nrow(m) == 0 || ncol(m) < 2) {
  save_plot(empty_plot("No divergence data"), prefix, formats)
  quit(save = "no", status = 0)
}

div     <- m[[1]]
classes <- colnames(m)[-1]
long <- data.frame(
  Div   = rep(div, times = length(classes)),
  Class = factor(rep(classes, each = length(div)), levels = classes),
  bp    = as.numeric(as.vector(as.matrix(m[, -1, drop = FALSE])))
)
long$Mb <- long$bp / 1e6

# Drop classes with no coverage; order by total abundance (largest at the bottom).
tot         <- tapply(long$bp, long$Class, sum)
keep        <- names(tot)[tot > 0]
long        <- long[long$Class %in% keep, ]
ord         <- names(sort(tot[keep], decreasing = TRUE))
long$Class  <- factor(long$Class, levels = ord)

if (nrow(long) == 0) {
  save_plot(empty_plot("No divergence data"), prefix, formats); quit(save = "no", status = 0)
}

# Assign colours (known classes from palette, the rest from the fallback ramp).
cols <- setNames(rep(NA_character_, length(ord)), ord)
for (cl in ord) if (!is.na(te_colours[cl])) cols[cl] <- te_colours[cl]
missing <- ord[is.na(cols)]
if (length(missing)) cols[missing] <- rep(fallback, length.out = length(missing))

p <- ggplot(long, aes(Div, Mb, fill = Class)) +
  geom_area(colour = NA, linewidth = 0) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = if (colnames(m)[1] == "PercDiv")
             "Divergence from consensus (%, RepeatMasker, not Kimura-corrected)"
           else "Kimura substitution level (CpG-adjusted, %)",
       y = "Sequence (Mb)") +
  theme_nature()

if (faceted) {
  p <- p + facet_wrap(~Class, scales = "free_y") +
    guides(fill = "none")
  save_plot(p, prefix, formats, width = 170, height = 120)
} else {
  save_plot(p, prefix, formats)
}
