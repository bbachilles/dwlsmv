# phi2_tail_grid.R -- the grids and the cell design of the Phi2 tail study
# (docs/notes/phi2-tail/prereg.md §3), shared by tools/make_phi2_tail_ref.R and
# docs/notes/phi2-tail/design_phi2.R so that the two cannot drift apart.
# Expects `smoke` (logical) to be set by the caller.

# ---- grids (docs/notes/phi2-tail/prereg.md §3) -----------------------------
L <- 4; B <- 0.999
if (smoke) {
  G_cert <- c(-3.7, 0.3, 3.3); G_env <- 5.3
  rho_grid <- c(-0.95, 0.15, 0.97)
} else {
  G_cert <- seq(-L, L, by = 0.25); G_env <- c(-6, -5, -4.5, 4.5, 5, 6)
  rho_tail <- c(0.9999, 0.9995, 0.999, 0.998, 0.9975, 0.995, 0.99, 0.985,
                0.98, 0.975, 0.97, 0.95)
  rho_grid <- sort(unique(c(-rho_tail, rho_tail, -0.9, -0.8, -0.6, -0.3, 0,
                            0.3, 0.6, 0.8, 0.9)))
}
G <- sort(c(G_cert, G_env))
in_cert <- function(v, rho) all(!is.finite(v) | abs(v) <= L) && abs(rho) <= B

# ---- Arm B: the cell design ------------------------------------------------
# Item templates (thresholds) on the 0.25 lattice, slopes in {-1.5, 0, 1.5}:
# every shifted bound lies in G and within +-L.
if (smoke) {
  tmpl <- list(c(-1.1, 1.3), 0.3); slopes <- c(0, 1.2)
} else {
  tmpl <- list(c(-1, 0, 1), c(-2.5, 0, 2.5), c(2, 2.5), c(-2.5, -2), 0.25, 2.5)
  slopes <- c(-1.5, 0, 1.5)
}
# the slope pairs of a template pair: x = 0 makes the slopes irrelevant to the
# bounds and gradients, so it gets one combination
slope_pairs <- function(x) {
  if (x == 0) return(list(c(0, 0)))
  lapply(seq_len(length(slopes)^2), function(q)
    c(slopes[(q - 1) %/% length(slopes) + 1], slopes[(q - 1) %% length(slopes) + 1]))
}
