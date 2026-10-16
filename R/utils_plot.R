# =============================================================================
# utils_plot.R —— shared plotting style for all main and supplementary figures.
# Principle: every number in a figure comes from a result table; no statistical
# decisions are taken inside figure scripts.
# =============================================================================

PAL <- list(
  AD = "#C2553A",          # adipose
  BM = "#3A6EA5",          # bone marrow
  UC = "#4E9A6B",          # umbilical cord
  claimed = "#3A6EA5",     # reproduced in the claimed direction
  reversed = "#C2553A",    # reproduced in the reverse direction
  not_rep = "#B8B8B8",     # not reproduced
  null = "#9E9E9E",
  accent = "#222222"
)

theme_msc <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key.size = grid::unit(3, "mm"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1)
    )
}

# save both PDF (vector graphics for submission) and PNG (preview); sizes in mm
save_fig <- function(p, name, width_mm = 180, height_mm = 120) {
  dir.create(proj_path("figures"), showWarnings = FALSE, recursive = TRUE)
  f_pdf <- proj_path("figures", paste0(name, ".pdf"))
  f_png <- proj_path("figures", paste0(name, ".png"))
  ggplot2::ggsave(f_pdf, p, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf)
  ggplot2::ggsave(f_png, p, width = width_mm, height = height_mm, units = "mm", dpi = 300)
  log_msg("saved figure ", name)
  invisible(c(f_pdf, f_png))
}
