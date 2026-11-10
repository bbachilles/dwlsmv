# make_rootfind_embed.R -- reference values for Module 2 option-B unit tests.
#
# Sources the R code of spec/derivations/dwlsmv-derivations.qmd (whole file,
# extracted with knitr::purl, run with cwd spec/derivations so that every chunk,
# every stopifnot of @sec-rootfind included, runs as written), then writes the
# cases of @sec-rootfind to CSV. No mathematics is copied here: every number
# comes from the .qmd's own rf_* functions and objects.
#
# Drafted by a read-only agent (2026-09-30); adapted by the main session after
# the .qmd's aborts and new controls (4b11c87), and run there.
#
# RERUNNABLE. Two consecutive runs give byte-identical CSVs (set.seed is inside
# the .qmd; nothing else is random). About 3.3 minutes. Needs the gitignored
# replicates in fixtures/replication/ (fixtures/generate/50_replicates.R).
# Then tools/emit_rootfind_embed.R turns the CSVs into the committed
# tests/unit/rootfind_embed.mata. From the repo root, in Git Bash:
#   R="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
#   "$R" tools/make_rootfind_embed.R . "$TEMP/optB/embed"
#   "$R" tools/emit_rootfind_embed.R "$TEMP/optB/embed" .
# Outputs (LF, comma, unquoted, %.17g, NA = undefined): cases.csv, tables.csv,
# thresholds.csv, aborts.csv, constants.csv, provenance.txt, and the purled code.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[1] else ".", "/")
outdir <- if (length(args) >= 2) args[2] else file.path(tempdir(), "rootfind_embed")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, "/")
qmd <- file.path(repo, "spec", "derivations", "dwlsmv-derivations.qmd")
stopifnot(file.exists(qmd))

# ---- 1. extract and run the whole .qmd --------------------------------------
purled <- file.path(outdir, "derivations_purled.R")
knitr::purl(qmd, output = purled, documentation = 0, quiet = TRUE)
E <- globalenv()  # adaptation: some chunk resolves names in globalenv (a private env failed on dP_drho)
old <- setwd(dirname(qmd))
t0 <- Sys.time()
sys.source(purled, envir = E, keep.source = FALSE)
# cwd stays at spec/derivations: rf_Y() reads replicates by a relative path
cat("qmd code ran to completion in", format(Sys.time() - t0), "\n")

# ---- 2. writers: %.17g, comma, LF, no quoting --------------------------------
num <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
wcsv <- function(df, name) {
  cols <- lapply(df, function(v) if (is.numeric(v) && !is.integer(v)) num(v) else
    ifelse(is.na(v), "NA", as.character(v)))
  txt <- c(paste(names(df), collapse = ","), do.call(paste, c(cols, sep = ",")))
  stopifnot(!any(grepl(",", unlist(lapply(df, function(v) if (is.character(v)) v))))) # text never holds a comma
  con <- file(file.path(outdir, name), "wb"); on.exit(close(con))
  writeLines(txt, con, sep = "\n")
}

with(E, {
  cases <- list(); tabs <- list(); thrs <- list()
  add <- function(id, set, source, status_of, n, a, b, r, B = 0.999, cap = 40, capc = 54,
                  mod = "none", bracket = c(NA, NA), rk = NA, budget = NA,
                  group = NA, cond = NA, rep = NA, pair = NA, note = "none") {
    e <- r$e
    cases[[length(cases) + 1]] <<- data.frame(
      case_id = id, set = set, source = source, status_tested = status_of,
      group = group, cond = cond, rep = rep, pair = pair,
      K1 = nrow(n), K2 = ncol(n), B = B, cap = cap, capc = capc, modification = mod,
      bracket_lo = bracket[1], bracket_hi = bracket[2],
      expected_status = r$case, expected_converged = as.integer(r$converged),
      expected_rho = r$rho, S_over_N = if (is.null(e$S)) NA else e$S, c_scaler = if (is.null(e$c)) NA else e$c,
      tau_S = if (is.null(e$tauS)) NA else e$tauS,
      tau_rho = if (is.null(e$tauS)) NA else rf_taur(e),
      width = if (!is.null(r$L) && !is.null(r$R)) r$R$x - r$L$x else NA,
      L = if (is.null(r$L)) NA else r$L$x, R = if (is.null(r$R)) NA else r$R$x,
      evals = r$evals, evals_check = if (is.null(r$evals_check)) NA else r$evals_check,
      end_lo = if (is.null(r$ends)) NA else r$ends[[1]]$x,
      end_hi = if (is.null(r$ends)) NA else r$ends[[2]]$x,
      rk_rho = rk, rk_budget = budget, note = note, stringsAsFactors = FALSE)
    tabs[[length(tabs) + 1]] <<- data.frame(case_id = id, r = as.integer(row(n)), s = as.integer(col(n)),
                                            n = as.vector(n))
    thrs[[length(thrs) + 1]] <<- rbind(
      data.frame(case_id = id, var = "a", k = seq_along(a), value = a),
      data.frame(case_id = id, var = "b", k = seq_along(b), value = b))
  }

  # ---- Level R1: FAIL, CTRL, BOUND (section) ----
  P <- rf_L1[rf_L1$group %in% c("FAIL", "CTRL", "BOUND"), ]
  P <- P[order(match(P$group, c("FAIL", "CTRL", "BOUND")), P$cond, P$rep, P$pair), ]
  for (k in seq_len(nrow(P))) {
    p <- P[k, ]; ij <- as.integer(strsplit(p$pair, ",")[[1]]); q <- rf_pair(p$cond, p$rep, ij[1], ij[2])
    r <- rf_root(q$n, q$a, q$b)
    stopifnot(identical(r$rho, p$rf), r$case == p$case)          # same fit as Level R1
    add(sprintf("R1_%s_%s_r%d_%d%d", p$group, p$cond, p$rep, ij[1], ij[2]), "pairs", "section L1",
        r$case, q$n, q$a, q$b, r, rk = p$rk, budget = p$budget, group = p$group,
        cond = p$cond, rep = p$rep, pair = paste0(ij[1], "-", ij[2]))
  }

  # ---- statuses and negative controls (Level R5, section) ----
  bd1 <- rf_root(rf_bd, a2, b2); bd2 <- rf_root(rf_bd, a2, b2, B = 1 - 2^-30)
  add("R5_emptycell_2x2_B0999", "status", "section L5", "no sign change up to the bound", rf_bd, a2, b2, bd1,
      note = "empty-cell boundary; score exactly 0 from 1-2^-16")
  add("R5_emptycell_2x2_Bnear1", "status", "section L5", "no sign change up to the bound", rf_bd, a2, b2, bd2,
      B = 1 - 2^-30, note = "same table at B = 1-2^-30")
  p1 <- F1[1, ]
  isl <- function(x) if (x >= 0.97 && x <= 0.985) list(x = x, res = FALSE) else rf_eval(q1$n, q1$a, q1$b, x)
  add("R5_island", "status", "section L5", "gap in the usable region", q1$n, q1$a, q1$b, rf_L5_isl,
      mod = "evaluator declares [0.97;0.985] unresolvable (res=0)", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_second_sign_change", "status", "section L5", "more than one sign change", q1$n, q1$a, q1$b, rf_L5_two,
      mod = "S replaced by -|S| where x < -0.5", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_zero_at_end", "status", "section L5", "ambiguous endpoint sign", q1$n, q1$a, q1$b, rf_L5_amb,
      mod = "S replaced by 0 where x < -0.5", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_cap5", "status", "section L5", "evaluation cap", q1$n, q1$a, q1$b, rf_L5_cap, cap = 5,
      cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_endpoint_cap10", "status", "section L5", "evaluation cap", q1$n, q1$a, q1$b, rf_L5_cc, capc = 10,
      mod = "endpoint check's cap set to 10 per side", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_s0_touch", "status", "section L5", "root", q1$n, q1$a, q1$b, rf_L5_s0t,
      mod = "S set to 0 at x == 0 exactly", cond = p1$cond, rep = p1$rep, pair = "1-2",
      rk = p1$rk, budget = p1$budget, note = "touching zero at 0; walks on from the probe")
  add("R5_s0_cross", "status", "section L5", "root", q1$n, q1$a, q1$b, rf_L5_s0c,
      mod = "S replaced by -sign(x)*max(|x|*c; 2*tau_S)", cond = p1$cond, rep = p1$rep, pair = "1-2",
      note = "certified crossing at 0")
  add("R5_s0_flat", "status", "section L5", "ambiguous endpoint sign", q1$n, q1$a, q1$b, rf_L5_s0f,
      mod = "S replaced by 0 where |x| < 1e-6", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_s0_side", "status", "section L5", "root", q1$n, q1$a, q1$b, rf_L5_s0s,
      mod = "S set to 0 at x == 0; (0.4999;0.5) unusable", cond = p1$cond, rep = p1$rep, pair = "1-2",
      rk = p1$rk, budget = p1$budget, note = "only a walk from the wrong probe meets the hole")
  add("R5_s0_neg", "status", "section L5", "no sign change up to the usable end", q1$n, q1$a, q1$b, rf_L5_s0n,
      mod = "dP_drho negated; S set to 0 at x == 0; (-0.5;-0.4999) unusable", cond = p1$cond, rep = p1$rep, pair = "1-2")
  add("R5_flat_stretch", "status", "section L5", "evaluation cap", q1$n, q1$a, q1$b, rf_L5_flat,
      mod = "S set to tau_S/2 on [0.996;0.9975]", cond = p1$cond, rep = p1$rep, pair = "1-2",
      note = "the width condition: without it the loop accepts 0.99609375")

  # ---- aborts (production rc 3498), not statuses: the expected message only ----
  stopifnot(rf_L5_inf == "abort: unresolvable bracket end",
            rf_L5_ori == "abort: bracket not maximum-oriented",
            rf_L5_z0 == "abort: rho = 0 is not resolvable")
  aborts <- data.frame(
    case_id = c("A_unresolvable_bracket_end", "A_bracket_not_max_oriented", "A_rho0_unresolvable"),
    modification = c("bracket supplied; walk skipped", "bracket supplied; walk skipped",
                     "evaluator declares x == 0 unresolvable (res=0)"),
    bracket_lo = c(-0.5, 0, NA), bracket_hi = c(x_inf, 0.5, NA),
    expected = c(rf_L5_inf, rf_L5_ori, rf_L5_z0),
    note = c("x_inf = first -1+2^-k with an observed cell P <= 0", "S > 0 at both ends", "none"),
    stringsAsFactors = FALSE)
  for (id in aborts$case_id) {
    tabs[[length(tabs) + 1]] <- data.frame(case_id = id, r = as.integer(row(q1$n)), s = as.integer(col(q1$n)),
                                           n = as.vector(q1$n))
    thrs[[length(thrs) + 1]] <- rbind(data.frame(case_id = id, var = "a", k = seq_along(q1$a), value = q1$a),
                                      data.frame(case_id = id, var = "b", k = seq_along(q1$b), value = q1$b))
  }

  # ---- CONSTRUCTED here (not in the section): unresolvable point inside the bracket ----
  # The section's island evaluator, applied to the single point of the plain fit's
  # first evaluation after the bracket closed (a Fisher/bisection trial inside (l, h)).
  r0 <- rf_root(q1$n, q1$a, q1$b)
  # index of first evaluation strictly inside the final walk bracket: walk evals are e0 + rf_side's k
  wk <- rf_side(function(x) rf_eval(q1$n, q1$a, q1$b, x), rf_eval(q1$n, q1$a, q1$b, 0),
                sign(rf_eval(q1$n, q1$a, q1$b, 0)$S), 0.999)
  xin <- r0$seen[[1 + wk$k + 1]]$x
  stopifnot(xin > wk$inner$x, xin < wk$outer$x)
  hole <- function(x) if (x == xin) list(x = x, res = FALSE) else rf_eval(q1$n, q1$a, q1$b, x)
  rh <- rf_root(q1$n, q1$a, q1$b, f = hole)
  stopifnot(!rh$converged, rh$case == "unresolvable point inside the bracket")
  add("C_unresolvable_inside_bracket", "status", "constructed by agent",
      "unresolvable point inside the bracket", q1$n, q1$a, q1$b, rh,
      mod = sprintf("evaluator declares the single point x = %s unresolvable", num(xin)),
      cond = p1$cond, rep = p1$rep, pair = "1-2",
      note = "hole = first trial after the walk bracket; constructed not taken from the section")

  # ---- negative controls on every FAIL pair (section) ----
  for (k in seq_len(nrow(F1))) {
    p <- F1[k, ]; ij <- as.integer(strsplit(p$pair, ",")[[1]]); q <- rf_pair(p$cond, p$rep, ij[1], ij[2])
    tag <- sprintf("%s_r%d_%d%d", p$cond, p$rep, ij[1], ij[2])
    rv <- rf_root(q$n, q$a, q$b, dS = function(a, b, x) -dP_drho(a, b, x))
    dr <- rf_root(q$n, q$a, q$b, dS = rf_dP_drop)
    add(paste0("R5_reversed_", tag), "negcontrol", "section L5", "no sign change up to the usable end",
        q$n, q$a, q$b, rv, mod = "dP_drho negated", rk = p$rk, budget = p$budget,
        group = "FAIL", cond = p$cond, rep = p$rep, pair = paste0(ij[1], "-", ij[2]))
    add(paste0("R5_dropped_", tag), "negcontrol", "section L5", "wrong root: |rho - rk_rho| > rk_budget",
        q$n, q$a, q$b, dr, mod = "dP_drho without its +phi2(a_{r-1};b_{s-1}) corner", rk = p$rk,
        budget = p$budget, group = "FAIL", cond = p$cond, rep = p$rep, pair = paste0(ij[1], "-", ij[2]),
        note = "caught by comparison with rk_rho not by the algorithm")
  }

  C <- do.call(rbind, cases)
  wcsv(C, "cases.csv"); wcsv(do.call(rbind, tabs), "tables.csv"); wcsv(do.call(rbind, thrs), "thresholds.csv")
  wcsv(aborts, "aborts.csv")
  wcsv(data.frame(name = c("eps", "delta_Phi", "delta_P", "B", "cap", "cap_endpoint_per_side", "frontier_ulps",
                           "tau_rho_factor", "tau_rho_ulps", "tauD_const", "rk_stop", "H_step"),
                  value = c(rf_eps, rf_dPhi, rf_dP, 0.999, 40, 54, 4, 16, 8, 7, 1e-15, 1e-5)), "constants.csv")
  print(C[, c("case_id", "expected_status", "expected_converged", "expected_rho", "evals")], digits = 17)
})

# ---- 3. provenance ----------------------------------------------------------
head <- system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)
dirty <- system2("git", c("-C", repo, "status", "--porcelain", "--", "spec/derivations/dwlsmv-derivations.qmd"), stdout = TRUE)
pk <- c("knitr", "haven", "mvtnorm", "lavaan")
prov <- c(paste0("repo_head=", head), paste0("qmd_dirty=", length(dirty) > 0),
          paste0("qmd_md5=", unname(tools::md5sum(qmd))), paste0("R=", R.version.string),
          paste0(pk, "=", vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "absent"), "")),
          paste0("loaded=", paste(sort(loadedNamespaces()), collapse = " ")))
con <- file(file.path(outdir, "provenance.txt"), "wb"); writeLines(prov, con, sep = "\n"); close(con)
cat("done\n")
