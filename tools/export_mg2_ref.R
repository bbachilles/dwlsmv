# ---------------------------------------------------------------------------
# dwlsmv -- export the multiple-group sigma/Delta reference values of
# derivation @sec-mg-sigma-delta into tests/unit/mg_model_ref.mata
#
# Runs the derivation's chunks in document order through mg2-delta (base R,
# no lavaan needed past the chunks themselves), then writes, for each case,
# a Mata function returning the parameter table (with its group key), theta,
# the per-group category counts, and mg2_sigma / mg2_delta at that theta.
# Cases: tabS at thS (two factors, two groups, every row type, ties, beta,
# K differing between groups), and design E at all four levels at mg2_point.
# The values are the derivation's closed forms, which mg2-delta asserts
# against central differences and mg2-lavaan / mg2-lavdelta against lavaan.
#
# Run from the repository root (writes one file, LF line endings):
#     Rscript tools/export_mg2_ref.R
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek
# ---------------------------------------------------------------------------

root <- normalizePath(".")
qmd  <- file.path(root, "spec", "derivations", "dwlsmv-derivations.qmd")
.mg2_out <- file.path(root, "tests", "unit", "mg_model_ref.mata")

# The chunks assign at top level, so the loop's own names live in a function.
.run_chunks <- function(qmd, last) {
  pur <- tempfile(fileext = ".R")
  knitr::purl(qmd, output = pur, quiet = TRUE, documentation = 1)
  x <- readLines(pur, encoding = "UTF-8")
  hdr <- grep("^## ----[^ -]", x)
  lbl <- sub("^## ----([^-,]+(-[^-,]+)*).*$", "\\1", x[hdr])
  stop_at <- which(lbl == last)
  stopifnot(length(stop_at) == 1)
  old <- setwd(dirname(qmd))
  on.exit(setwd(old))
  for (k in seq_len(stop_at)) {
    j <- c(hdr[hdr > hdr[k]], length(x) + 1)[1]
    eval(parse(text = x[hdr[k]:(j - 1)], encoding = "UTF-8"), envir = globalenv())
  }
}
.run_chunks(qmd, "mg2-delta")


ptype <- c(lambda = 1, psi = 2, beta = 3, tau = 4, nu = 5, omega = 6, alpha = 7)
num <- function(v) paste(sprintf("%.17g", v), collapse = ", ")
mvec <- function(v) if (length(v)) paste0("(", gsub(", ", " \\\\ ", num(v)), ")") else "J(0, 1, .)"
mmat <- function(M) paste0("(", paste(apply(M, 1, num), collapse = " \\ "), ")")
# One statement per row: a whole Delta as one literal exceeds Mata's token limit.
mrows <- function(nm, M) c(sprintf("    %s = J(%d, %d, .)", nm, nrow(M), ncol(M)),
                          sprintf("    %s[%d, .] = (%s)", nm, seq_len(nrow(M)), apply(M, 1, num)))

case_fn <- function(name, tab, th) {
  R <- tab$rows
  nm <- sprintf("%s g%d %d %d", R$type, R$g, R$i, R$j)
  fx <- ifelse(R$free, R$fidx, 0L)
  vl <- ifelse(R$free, 0, R$val)
  sig <- mg2_sigma(tab, th); Dl <- mg2_delta(tab, th)
  K <- sapply(tab$K, identity)
  c(sprintf("void mgref_%s(struct dwlsmv_parmap scalar par, real colvector th,", name),
    "                real matrix K, real colvector sig, real matrix Dl)",
    "{",
    sprintf("    par.npar  = %d", tab$npar),
    sprintf("    par.nall  = %d", tab$nall),
    sprintf("    par.name  = (%s)", paste0('"', nm, '"', collapse = " \\ ")),
    sprintf("    par.ptype = %s", mvec(ptype[R$type])),
    sprintf("    par.row   = %s", mvec(R$i)),
    sprintf("    par.col   = %s", mvec(R$j)),
    sprintf("    par.free  = %s", mvec(as.integer(R$free))),
    sprintf("    par.value = %s", mvec(vl)),
    sprintf("    par.fidx  = %s", mvec(fx)),
    sprintf("    par.grp   = %s", mvec(R$g)),
    sprintf("    th  = %s", mvec(th)),
    sprintf("    K   = %s", mmat(K)),
    sprintf("    sig = %s", mvec(sig)),
    mrows("Dl", Dl),
    "}", "")
}

lines <- c(
  "*! Reference values for tests/unit/test_mg_model.do, WRITTEN BY tools/export_mg2_ref.R",
  "*! from derivation @sec-mg-sigma-delta (mg2_sigma, mg2_delta). Do not edit by hand.",
  sprintf("*! Generated %s, %s.", format(Sys.Date()), R.version.string),
  "*! K is p x G: column g is group g's category counts.",
  "*!", "*! SPDX-License-Identifier: GPL-3.0-or-later",
  "*! SPDX-FileCopyrightText: 2026 Erik Ruzek", "", "version 16", "", "mata:", "",
  case_fn("S", tabS, thS))
KE <- inv_designs$E
for (l in inv_levels) {
  tab <- mg2_level_table(KE, l)
  lines <- c(lines, case_fn(paste0("E_", l), tab, mg2_point(tab)))
}
lines <- c(lines, "end")
con <- file(.mg2_out, open = "wb"); writeLines(lines, con, sep = "\n"); close(con)
cat("wrote", .mg2_out, length(lines), "lines\n")
