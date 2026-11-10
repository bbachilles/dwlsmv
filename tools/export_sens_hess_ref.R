# export_sens_hess_ref.R -- references for tests/unit/test_sens_hess.do
#
# Writes tests/unit/sens_hess_ref.mata: at 12 FIXED parameter points of design E
# (G = 2, K = 4 4 4 2 2), theta, sigma(theta) and the nonzero entries of Result
# MS6's Hessian H (@sec-mg-hessian), in dwlsmv_sens_ddelta's layout (column
# (l - 1) q + i is d2 sigma / d theta_i d theta_l). The functions are not
# transcribed: every definition is evaluated from
# spec/derivations/dwlsmv-derivations.qmd's own chunks (mg2_*, mg2h_*, inv_tau,
# inv_lam, inv_designs), picked by name. The points are fixed values, not fits,
# so the references do not depend on Stata's version or processor count.
#
# Points, per level (configural, thresholds, loadings = Stata's metric,
# intercepts = Stata's scalar): g = mg2_point (MS6's preregistered generic
# point), s = mg2_start (Module 3's start), c = g + 0.03 cos(1:q) (a third,
# with no two theta elements equal).
#
# Run from the repository root:  Rscript tools/export_sens_hess_ref.R
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

qmd <- "spec/derivations/dwlsmv-derivations.qmd"
want <- c("inv_tau", "inv_lam", "inv_designs", "inv_levels",
          "mg2_block", "mg2_row", "mg2_fidx", "mg2_table", "mg2_tpos", "mg2_Dg", "mg2_off", "mg2_pars",
          "mg2_sigma_g", "mg2_sigma", "mg2_dcol_g", "mg2_dcols", "mg2_delta",
          "mg2_level_table", "mg2_point", "mg2_start", "mg2h_domain", "mg2h_d2", "mg2h_hess")
x <- readLines(qmd, encoding = "UTF-8")
open <- grep("^```\\{r", x); close <- grep("^```$", x)
env <- new.env()
for (o in open) {
  e <- close[close > o][1]
  ex <- tryCatch(parse(text = x[(o + 1):(e - 1)], keep.source = FALSE), error = function(err) NULL)
  for (z in ex) {
    if (is.call(z) && identical(z[[1]], as.name("<-")) && is.name(z[[2]]) && as.character(z[[2]]) %in% want)
      eval(z, env)
  }
}
miss <- setdiff(want, ls(env))
if (length(miss)) stop("not found in the qmd: ", paste(miss, collapse = ", "))

out <- with(env, {
  KE <- inv_designs$E
  stata <- c(configural = "configural", thresholds = "thresholds", loadings = "metric", intercepts = "scalar")
  pts <- list()
  for (L in inv_levels) {
    tab <- mg2_level_table(KE, L); g <- mg2_point(tab)
    pts[[length(pts) + 1L]] <- list(lev = stata[[L]], nm = "g", tab = tab, th = g)
    pts[[length(pts) + 1L]] <- list(lev = stata[[L]], nm = "s", tab = tab, th = mg2_start(tab, inv_tau(KE)))
    pts[[length(pts) + 1L]] <- list(lev = stata[[L]], nm = "c", tab = tab, th = g + 0.03 * cos(seq_along(g)))
  }
  lapply(pts, function(p) {
    H <- mg2h_hess(p$tab, p$th); D <- dim(H)[1]; q <- dim(H)[2]
    M <- matrix(H, D, q * q); nz <- which(M != 0, arr.ind = TRUE)
    c(p[c("lev", "nm", "th")], list(sig = mg2_sigma(p$tab, p$th), nz = cbind(nz, M[nz]), D = D, q = q))
  })
})

f17 <- function(v) paste(sprintf("%.17g", v), collapse = ", ")
L <- c("*! sens_hess_ref.mata -- written by tools/export_sens_hess_ref.R; do not edit by hand",
       sprintf("*! Written %s by %s from %s's own chunks (mg2_*, mg2h_*).", Sys.Date(), R.version.string, qmd),
       "*! Fixed points, not fits: independent of Stata's version and processor count.",
       "*!",
       "*! SPDX-License-Identifier: GPL-3.0-or-later",
       "*! SPDX-FileCopyrightText: 2026 Erik Ruzek",
       "", "mata:", "",
       "// sh_ref(k): k-th point. Returns (level, name) and fills theta, sigma, and H's",
       "// nonzero entries as rows (row, column, value) of nz; D and q in dq.",
       "string rowvector sh_ref(real scalar k, real colvector th, real colvector sig, real matrix nz, real rowvector dq)",
       "{")
for (k in seq_along(out)) {
  p <- out[[k]]
  L <- c(L, sprintf("    if (k == %d) {", k),
         sprintf("        th  = (%s)'", f17(p$th)),
         sprintf("        sig = (%s)'", f17(p$sig)),
         sprintf("        dq  = (%d, %d)", p$D, p$q),
         sprintf("        nz  = J(%d, 3, .)", nrow(p$nz)))
  rws <- split(seq_len(nrow(p$nz)), ceiling(seq_len(nrow(p$nz)) / 50))
  for (r in rws)
    L <- c(L, sprintf("        nz[(%d::%d), .] = colshape((%s), 3)", min(r), max(r), f17(t(p$nz[r, , drop = FALSE]))))
  L <- c(L, sprintf("        return((\"%s\", \"%s\"))", p$lev, p$nm), "    }")
}
L <- c(L, "    _error(\"sh_ref: no such point\")", "}", "",
       "real scalar sh_nref()", "{", sprintf("    return(%d)", length(out)), "}", "", "end")
con <- file("tests/unit/sens_hess_ref.mata", "wb"); writeLines(L, con, sep = "\n"); close(con)
cat(sprintf("wrote %d points; nonzeros per point %s\n", length(out),
            paste(sapply(out, function(p) nrow(p$nz)), collapse = " ")))
