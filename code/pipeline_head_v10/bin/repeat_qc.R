#!/usr/bin/env Rscript
# =============================================================================
# repeat_qc.R — repeat detection & masking quality-control figure, Nature-styled.
#
# Published next to the figure/table so you can re-style or extend it:
#
#   Rscript repeat_qc.R <composition.tsv> <summary.tsv> <prefix> <formats>
#
#   <composition.tsv>  columns: class  bp  pct    (pct = % of total genome length)
#   <summary.tsv>      columns: metric  value  status   (QC metrics + PASS/WARN)
#   <prefix>           output basename (e.g. sampleX.repeat_qc)
#   <formats>          comma list png,svg,pdf  (default png,svg)
#
# Only base R + ggplot2 are required (ragg/svglite used if present).
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
comp_f  <- if (length(args) >= 1) args[1] else "composition.tsv"
summ_f  <- if (length(args) >= 2) args[2] else "summary.tsv"
prefix  <- if (length(args) >= 3) args[3] else "repeat_qc"
formats <- if (length(args) >= 4) trimws(strsplit(args[4], ",")[[1]]) else c("png", "svg")

theme_nature <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line   = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.3, colour = "black"),
      axis.text   = element_text(colour = "black", size = base_size - 1),
      plot.title  = element_text(size = base_size, face = "bold"),
      legend.position = "none"
    )
}

save_plot <- function(p, base, formats, width = 110, height = 75, units = "mm", dpi = 600) {
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

ok <- file.exists(comp_f) && file.info(comp_f)$size > 0
comp <- if (ok) tryCatch(read.table(comp_f, header = TRUE, sep = "\t",
                                    stringsAsFactors = FALSE, comment.char = ""),
                         error = function(e) NULL) else NULL

# Colour-blind-safe palette keyed by top-level class.
te_colours <- c(
  LTR = "#4E79A7", LINE = "#F28E2B", SINE = "#E15759", DNA = "#76B7B2",
  RC = "#59A14F", Retroposon = "#9C755F", Satellite = "#EDC948",
  Simple_repeat = "#B07AA1", Low_complexity = "#FF9DA7", Unknown = "#BAB0AC")

if (is.null(comp) || nrow(comp) == 0) {
  save_plot(empty_plot("No repeat composition"), paste0(prefix, "_composition"), formats)
} else {
  # Without a genome length (no 'total length' in the .tbl) pct is NA: plot Mb instead.
  comp$pct <- suppressWarnings(as.numeric(comp$pct))
  by_pct <- any(!is.na(comp$pct))
  comp$value <- if (by_pct) comp$pct else comp$bp / 1e6
  comp <- comp[order(comp$value), ]
  comp$class <- factor(comp$class, levels = comp$class)
  cols <- sapply(as.character(comp$class), function(x)
    if (!is.na(te_colours[x])) te_colours[x] else "#9AA0A6")
  # Per class: the union of its masked bases (classes may overlap, so bars can sum
  # to more than the masked fraction of the genome).
  p <- ggplot(comp, aes(value, class, fill = class)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = if (by_pct) "% of genome" else "Masked sequence (Mb)", y = NULL,
         title = "Repeat composition") +
    theme_nature()
  save_plot(p, paste0(prefix, "_composition"), formats)
}

# Echo the QC summary to the log so it is visible in .command.out too.
if (file.exists(summ_f) && file.info(summ_f)$size > 0) {
  s <- read.table(summ_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  comment.char = "", quote = "")
  message("── Repeat detection & masking QC ─────────────────────────")
  for (i in seq_len(nrow(s)))
    message(sprintf("  %-28s %-14s %s", s$metric[i], s$value[i], s$status[i]))
}
