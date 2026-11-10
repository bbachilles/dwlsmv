# ---------------------------------------------------------------------------
# dwlsmv -- the R reference of tests/unit/test_sens_t3.do: B_T3's terms by the
# subsection's own functions (@sec-mg-sensitivity-t3), on dwlsmv_sens_t3()'s
# own inputs, written to tests/unit/sens_t3_ref.mata
#
# Executed 2026-10-16 (R 4.6.0) on export_sens_t3_in.do's output; wrote tests/unit/sens_t3_ref.mata.
#
# Run from the repository root, after tools/export_sens_t3_in.do:
#     Rscript tools/export_sens_t3_ref.R <dir with in_<case>.csv> [<out .mata> <provenance>]
# With three arguments it writes a per-version, per-count reference (policy 2026-10-19:
# tests/unit/sens_t3_ref_s<version>_p<processors>.mata; the Stata 18 six-processor set stays
# tests/unit/sens_t3_ref.mata).
#
# The t3_* functions are taken from spec/derivations/dwlsmv-derivations.qmd
# itself: every R chunk is parsed, and only the top-level assignments of a
# function to a name starting "t3_" are evaluated (no chunk is run, so no
# lavaan and no earlier chunk is needed). t3_u is set here to 2^-53, the
# chunk's value. The functions then see Stata's inputs:
#   - sigma is the constant function returning Stata's sigma(theta-hat), since
#     every t3_ function evaluates it at theta-hat only;
#   - dD, Ebar, J are Stata's (dwlsmv_sens_ddelta, _ebar, _jmat), which R
#     cannot rebuild without the model; the test pins each in Stata;
#   - Theta is recomputed in R by t3_Theta() from Stata's P, A, c, r, step and
#     K, and written beside the terms, so the test can compare Stata's Theta.
# What is compared is the assembly: t3_grad, t3_gth_abs, t3_beta, t3_beval
# (with t3_eta_q4 (SE6; production since 2026-10-20), c_ch = 6(q + 1)), t3_belem, and TS7's formula.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1, 3))   # <indir> [<out> <provenance>]: a per-version, per-count reference
indir <- args[1]
qmd <- "spec/derivations/dwlsmv-derivations.qmd"
out <- if (length(args) == 3) args[2] else "tests/unit/sens_t3_ref.mata"
prov <- if (length(args) == 3) args[3] else "StataMP 18, six processors"

## the subsection's functions
L <- readLines(qmd, encoding = "UTF-8")
st <- grep("^```\\{r", L); en <- grep("^```\\s*$", L)
code <- unlist(lapply(st, function(a) { b <- en[en > a][1]; if (b > a + 1) L[(a + 1):(b - 1)] else character(0) }))
ex <- parse(text = code, keep.source = FALSE)
env <- new.env(parent = globalenv())
assign("t3_u", 2^-53, envir = env)
nfun <- 0L
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      startsWith(as.character(e[[2]]), "t3_") && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = env); nfun <- nfun + 1L
  }
}
need <- c("t3_grad", "t3_at", "t3_tr", "t3_of", "t3_inv", "t3_gth_abs", "t3_beta", "t3_nT_of",
          "t3_beval", "t3_evgrad", "t3_eta_q4", "t3_belem", "t3_Theta", "t3_rT")
stopifnot(all(vapply(need, exists, NA, envir = env, inherits = FALSE)))
cat(sprintf("ST3REF|%d t3_ function definitions taken from the qmd\n", nfun))

## %21x (docs/notes/sens-grid/analyze.R's decoder; no missing values expected here)
HXRE <- "^([+-])([01])\\.([0-9a-fA-F]{13})X([+-][0-9a-fA-F]{3})$"
hx <- function(s) {
  s <- trimws(s); stopifnot(all(grepl(HXRE, s)))
  e <- sub(HXRE, "\\4", s); stopifnot(all(e != "+3ff"))
  sg <- ifelse(sub(HXRE, "\\1", s) == "-", -1, 1); lead <- as.numeric(sub(HXRE, "\\2", s))
  fr <- sub(HXRE, "\\3", s)
  ex <- strtoi(substring(e, 2), 16L) * ifelse(substr(e, 1, 1) == "-", -1, 1)
  ex <- ifelse(lead == 0, -1022, ex)
  sg * vapply(seq_along(s), function(i) {
    d <- strtoi(strsplit(fr[i], "")[[1]], 16L); (lead[i] + sum(d * 16^-(seq_along(d)))) * 2^ex[i] }, 0)
}
stopifnot(identical(hx(c("+1.0000000000000X+000", "-1.8000000000000X-001", "+0.0000000000000X-3ff")), c(1, -0.75, 0)))

rd <- function(f) {
  d <- read.csv(f, header = FALSE, colClasses = "character", col.names = c("nm", "i", "j", "hex"))
  d$i <- as.integer(d$i); d$j <- as.integer(d$j); d$v <- hx(d$hex)
  lapply(split(d, d$nm), function(z) { M <- matrix(NA_real_, max(z$i), max(z$j)); M[cbind(z$i, z$j)] <- z$v
    stopifnot(!anyNA(M)); M })
}

cases <- c("metric_full", "scalar_full", "metric_none", "metric_bigeps")
res <- list()
with(env, {
  for (cs in cases) {
    f <- file.path(indir, paste0("in_", cs, ".csv"))
    z <- rd(f)
    th <- drop(z$th); Dl <- z$Dl; Om <- z$Om; s <- drop(z$s); sgv <- drop(z$sg); q <- length(th)
    nT <- drop(z$N) - drop(z$G); T3s <- drop(z$T3); epsV <- drop(z$epsV); rhoT <- drop(z$rhoT)
    sig <- function(x) sgv
    g <- t3_grad(sig, th, Dl, Om, s = s, nT = nT)
    dD <- lapply(seq_len(q), function(l) z$dD[, ((l - 1) * q + 1):(l * q), drop = FALSE])
    emax <- max(abs(z$eps))
    gtha <- t3_gth_abs(g, dD, g$gh, z$Eb, Dl, eps_max = emax)
    o <- list(P = crossprod(Dl, drop(z$wf) * Dl), A = z$A, r = s - sgv, D = Dl, c = drop(z$c))
    ThR <- t3_Theta(o, drop(z$wf), z$Eb, eps = drop(z$eps), ds = drop(z$ds), K = drop(z$K), stK = drop(z$step))
    ThS <- (1 + abs(th)) * ((if (drop(z$full) == 1) z$sens[, 5] else 0) + z$sens[, 3] + z$sens[, 6])
    bb <- t3_beta(g, z$Eb, gtha, ThS, z$Jm, dOm = z$dOm, ds = drop(z$ds))
    be <- t3_beval(Dl, Om, g$r, nT, t3_eta_q4)              # SE6: production uses Q4 (2026-10-20)
    bel <- t3_belem(g, Dl, sig, th, eps_s5 = epsV)
    dn <- 1 + abs(T3s)
    v <- c(T3s, 2 * bb[["b1"]] / dn, bb[["b2"]] / dn, bb[["b3"]] / dn, be / dn, bel / dn, NA, NA)
    v[7] <- (1 + rhoT) * (v[2] + v[3] + v[4]) + v[5] + v[6] + t3_rT(T3s)
    v[8] <- v[5] + t3_rT(T3s)
    stata <- drop(z$out)
    cat(sprintf("ST3REF|%s|R T3 %.17g, Stata T3 %.17g, rel %.3g\n", cs, g$v[["T3"]], T3s, abs(g$v[["T3"]] / T3s - 1)))
    cat(sprintf("ST3REF|%s|col %d R %.6e Stata %.6e rel %.3g\n", cs, 1:8, v, stata, abs(stata - v) / abs(v)), sep = "")
    cat(sprintf("ST3REF|%s|Theta max rel (R t3_Theta vs Stata) %.3g\n", cs, max(abs(ThR - ThS) / abs(ThR))))
    res[[cs]] <- list(v = v, Th = ThR, md5 = unname(tools::md5sum(f)), T3R = g$v[["T3"]])
  }
})
res <- get("res", envir = env)

## write the reference
fmt <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
o <- c(sprintf("*! %s -- written by tools/export_sens_t3_ref.R; do not edit by hand", basename(out)),
       sprintf("*! Inputs from %s.", prov),
       "*!",
       sprintf("*! Written %s by R %s from %s's t3_ functions, on the inputs", format(Sys.Date()), getRversion(), qmd),
       "*! tools/export_sens_t3_in.do wrote (md5 of each input file below). Values in %.17g.",
       "*!",
       "*! SPDX-License-Identifier: GPL-3.0-or-later",
       "*! SPDX-FileCopyrightText: 2026 Erik Ruzek",
       "", "mata:", "")
for (cs in cases) {
  o <- c(o, sprintf("// %s: input md5 %s; R's own T3 %.17g", cs, res[[cs]]$md5, res[[cs]]$T3R))
}
o <- c(o, "", "real rowvector st3_ref(string scalar cs)", "{")
for (cs in cases) o <- c(o, sprintf("    if (cs == \"%s\") return((%s))", cs, fmt(res[[cs]]$v)))
o <- c(o, "    _error(\"st3_ref: unknown case\")", "}", "", "real colvector st3_ref_theta(string scalar cs)", "{")
for (cs in cases) o <- c(o, sprintf("    if (cs == \"%s\") return((%s)')", cs, fmt(res[[cs]]$Th)))
o <- c(o, "    _error(\"st3_ref_theta: unknown case\")", "}", "", "end")
con <- file(out, "wb"); writeLines(o, con, sep = "\n"); close(con)
cat("ST3REF|wrote", out, "\n")
