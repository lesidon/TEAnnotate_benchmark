#!/usr/bin/env Rscript
# =============================================================================
# plot_kimura.R — TE divergence (Kimura) landscape.
#
# This is the exact script the pipeline used to draw the plot. It is published
# next to the figure so you can tweak colours, sizes, labels, etc. and re-run:
#
#   Rscript plot_kimura.R <kimura_matrix.tsv> <prefix> <formats> <layout> [genome_bp] [max_per_class] [min_share]
#
#   <kimura_matrix.tsv>  RepeatMasker/calcDivergenceFromAlign coverage matrix
#                        (first column "Div" = Kimura %, then one bp column per
#                        repeat label, e.g. "DNA/hAT-Charlie"; "PercDiv" = the
#                        pipeline's fallback from the .out's plain percent
#                        divergence, labelled as such)
#   <prefix>             output basename (e.g. sampleX.kimura_landscape)
#   <formats>            comma list: png,svg,pdf   (default png,svg)
#   <layout>             both (default) | stacked | faceted. The stacked plot is
#                        written to <prefix>.<fmt>, the one-panel-per-class plot
#                        to <prefix>_faceted.<fmt>  (1/0 = faceted/stacked)
#   [genome_bp]          genome length in bp; the y axis is then % of the genome
#                        (0 or empty = Mb of sequence)
#   [max_per_class]      superfamilies shown separately per class (default 4)
#   [min_share]          ...and only if they cover at least this fraction of all
#                        TE bases (default 0.01); the rest of a class is "Other"
#
# Labels are collapsed to superfamilies (DNA/hAT-Charlie -> DNA/hAT), each class
# gets one hue and its superfamilies shades of it (darker = more abundant), and
# simple, low-complexity, satellite and RNA repeats are left out, so the legend
# stays short and every colour on the plot is a different element.
# Only base R + ggplot2 are required (ragg/svglite used if present).
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))

args          <- commandArgs(trailingOnly = TRUE)
mat_file      <- if (length(args) >= 1) args[1] else "kimura_matrix.tsv"
prefix        <- if (length(args) >= 2) args[2] else "kimura_landscape"
formats       <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("png", "svg")
layout        <- if (length(args) >= 4) tolower(args[4]) else "both"
genome_bp     <- if (length(args) >= 5) suppressWarnings(as.numeric(args[5])) else NA
max_per_class <- if (length(args) >= 6) as.integer(args[6]) else 4L
min_share     <- if (length(args) >= 7) as.numeric(args[7]) else 0.01
formats       <- trimws(formats)
layout        <- switch(layout, "1" = , "true" = "faceted", "0" = , "false" = "stacked", layout)
if (!layout %in% c("both", "stacked", "faceted")) layout <- "both"
draw_stacked  <- layout %in% c("both", "stacked")
draw_faceted  <- layout %in% c("both", "faceted")
use_pct       <- !is.na(genome_bp) && genome_bp > 0
facet_prefix  <- paste0(prefix, "_faceted")

# ── Theme ────────────────────────────────────────────────────────────────────
# Set base_family = "Arial" (or "Helvetica") if the font is installed on your system.
theme_landscape <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size - 1),
      legend.key.size  = unit(2.6, "mm"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.position  = "right",
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 1, face = "bold")
    )
}

# One hue per TE class: ramp from light (rare / "Other") to dark (most abundant).
class_ramps <- list(
  DNA     = c("#d9f99d", "#2f8f2f", "#0b5d1e"),
  LINE    = c("#e4d9f2", "#8a5cc2", "#4b2a7b"),
  LTR     = c("#fde3b5", "#f68b2c", "#b8400b"),
  SINE    = c("#f7c6d3", "#c2185b", "#6e0e33"),   # rose, kept apart from the LTR oranges
  RC      = c("#c6e0f5", "#2a7fb7", "#164a6e"),
  Unknown = c("#d9d9d9", "#a6a6a6", "#7f7f7f")
)
class_order <- names(class_ramps)            # bottom-to-top stacking order

save_plot <- function(p, base, formats, width = 120, height = 75, units = "mm", dpi = 600) {
  for (fmt in formats) {
    f  <- paste0(base, ".", fmt)
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
save_empty <- function(msg) {
  if (draw_stacked) save_plot(empty_plot(msg), prefix, formats)
  if (draw_faceted) save_plot(empty_plot(msg), facet_prefix, formats)
  quit(save = "no", status = 0)
}

# ── Load ─────────────────────────────────────────────────────────────────────
ok <- file.exists(mat_file) && file.info(mat_file)$size > 0
m  <- if (ok) tryCatch(read.table(mat_file, header = TRUE, check.names = FALSE,
                                  comment.char = ""), error = function(e) NULL) else NULL
if (is.null(m) || nrow(m) == 0 || ncol(m) < 2) save_empty("No divergence data")
div_kind <- colnames(m)[1]
div      <- m[[1]]
labels   <- colnames(m)[-1]

# ── Label -> class and superfamily ───────────────────────────────────────────
top   <- sub("/.*$", "", labels)
cls   <- ifelse(top %in% c("Retroposon"), "SINE",
         ifelse(top %in% c("PLE"), "LINE",
         ifelse(top %in% c("DIRS"), "LTR",
         ifelse(top %in% c("Unknown", "Unspecified", "Unclassified"), "Unknown", top))))
sub_  <- ifelse(grepl("/", labels), sub("^[^/]*/", "", labels), "")
sfam  <- sub("-.*$", "", sub_)                      # hAT-Charlie -> hAT, ERVL-MaLR -> ERVL
sfam[tolower(sfam) %in% c("", "other", "unknown", "na")] <- "Other"
keep_lab <- cls %in% class_order                    # drops Simple_repeat, Satellite, rRNA, ...

long <- do.call(rbind, lapply(which(keep_lab), function(j)
  data.frame(Div = div, Class = cls[j], Sfam = sfam[j],
             bp = as.numeric(m[[j + 1]]), stringsAsFactors = FALSE)))
if (is.null(long) || sum(long$bp, na.rm = TRUE) == 0) save_empty("No TE divergence data")
long$bp[is.na(long$bp)] <- 0
long$Class[long$Class == "Unknown"] <- "Unknown"
long$Sfam[long$Class == "Unknown"]  <- "Other"      # Unknown is one element

# ── Keep the major superfamilies of each class; the rest becomes <class>/Other
tot_te <- sum(long$bp)
sf_tot <- aggregate(bp ~ Class + Sfam, data = long, FUN = sum)
sf_tot <- sf_tot[order(match(sf_tot$Class, class_order), -sf_tot$bp), ]
sf_tot$keep <- FALSE
for (cl in unique(sf_tot$Class)) {
  i <- which(sf_tot$Class == cl & sf_tot$Sfam != "Other" & sf_tot$bp / tot_te >= min_share)
  sf_tot$keep[head(i, max_per_class)] <- TRUE
}
long <- merge(long, sf_tot[, c("Class", "Sfam", "keep")], by = c("Class", "Sfam"))
long$Sfam[!long$keep] <- "Other"
long <- aggregate(bp ~ Div + Class + Sfam, data = long, FUN = sum)
long$Label <- ifelse(long$Class == "Unknown", "Unknown",
              ifelse(long$Sfam == "Other", paste0(long$Class, "/Other"),
                     paste0(long$Class, "/", long$Sfam)))
long$y  <- if (use_pct) 100 * long$bp / genome_bp else long$bp / 1e6

# ── Order and colours ────────────────────────────────────────────────────────
lab_tot <- aggregate(bp ~ Class + Sfam + Label, data = long, FUN = sum)
lab_tot <- lab_tot[lab_tot$bp > 0, ]
lab_tot <- lab_tot[order(match(lab_tot$Class, class_order), lab_tot$Sfam == "Other", -lab_tot$bp), ]
cols <- character(0)
for (cl in unique(lab_tot$Class)) {
  labs_cl <- lab_tot$Label[lab_tot$Class == cl]
  named   <- labs_cl[!grepl("/Other$", labs_cl)]
  ramp    <- colorRampPalette(rev(class_ramps[[cl]]))(max(length(named), 1) + 1)
  cl_cols <- setNames(ramp[seq_along(named)], named)
  if (any(grepl("/Other$", labs_cl)))          # "Other": the lightest tint of the class
    cl_cols <- c(cl_cols, setNames(class_ramps[[cl]][1], labs_cl[grepl("/Other$", labs_cl)]))
  if (cl == "Unknown") cl_cols <- setNames(class_ramps$Unknown[2], "Unknown")
  cols <- c(cols, cl_cols)
}
lev <- lab_tot$Label
long <- long[long$Label %in% lev, ]
long$Label <- factor(long$Label, levels = lev)
long$Class <- factor(long$Class, levels = class_order[class_order %in% long$Class])

# x range: up to the last 1% bin that still holds sequence, rounded up to 10, at most 60
has  <- tapply(long$bp, long$Div, sum)
xmax <- min(60, max(10, ceiling(max(as.numeric(names(has))[has > 0]) / 10) * 10))
long <- long[long$Div <= xmax, ]

p <- ggplot(long, aes(Div, y, fill = Label)) +
  geom_col(width = 1, colour = NA, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols, name = "Repeat class / superfamily",
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_x_continuous(limits = c(-0.5, xmax + 0.5), breaks = seq(0, xmax, 10), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = if (div_kind == "PercDiv")
             "Divergence from consensus (%, RepeatMasker, not Kimura-corrected)"
           else "Kimura substitution level (CpG-adjusted, %)",
       y = if (use_pct) "Genome (%)" else "Sequence (Mb)") +
  theme_landscape()

if (draw_stacked) save_plot(p, prefix, formats)
if (draw_faceted) {
  # one panel per class, each with its own y range so small classes stay readable
  save_plot(p + facet_wrap(~Class, scales = "free_y", ncol = 2),
            facet_prefix, formats, width = 170, height = 120)
}
