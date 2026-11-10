# ---------------------------------------------------------------------------
# Verify committed fixtures against their manifests.
#
#   Rscript fixtures/verify.R        (from the repository root)
#
# Exits 0 if every recorded file is present and unchanged, 1 otherwise.
#
# Both tiers are checked:
#
#   fixtures/MANIFEST.json         tier 1, 10_fixtures.R          (required)
#   fixtures/MANIFEST_lavaan.json  tier 2, 20_lavaan_reference.R  (optional)
#   fixtures/MANIFEST_mplus.json   tier 3, 30_mplus_reference.R   (optional)
#   fixtures/MANIFEST_mc.json      lavaan Monte Carlo SE reference,
#                                  70_lavaan_mc_reference.R        (optional)
#   fixtures/MANIFEST_lavaan_std.json  lavaan standardized solution (D7),
#                                  80_lavaan_standardized.R        (optional)
#   fixtures/MANIFEST_lavaan_emptycell.json  lavaan, one-empty-cell pair (O6),
#                                  90_lavaan_emptycell.R           (optional)
#   fixtures/MANIFEST_marker_trap.json  tier 1, marker-trap data and population,
#                                  95_marker_trap.R, fixture mode  (optional)
#   fixtures/MANIFEST_lavaan_marker_trap.json  lavaan, marker-trap reference,
#                                  96_lavaan_marker_trap.R         (optional)
#   fixtures/MANIFEST_lavaan_difftest.json  lavaan, DIFFTEST nested pairs,
#                                  97_lavaan_difftest.R            (optional)
#   fixtures/MANIFEST_lavaan_multigroup.json  lavaan, two-group design E,
#                                  98_lavaan_multigroup.R          (optional)
#   fixtures/MANIFEST_lavaan_mg_beta.json  lavaan, two-group design M with beta,
#                                  99_lavaan_mg_beta.R             (optional)
#   fixtures/MANIFEST_lavaan_mg_betafix.json  lavaan, the same data, beta fixed .4,
#                                  99b_lavaan_mg_betafix.R         (optional)
#   fixtures/MANIFEST_lavaan_mg_betag2.json  lavaan, the same data, beta g2only,
#                                  99c_lavaan_mg_betag2.R          (optional)
#   fixtures/MANIFEST_lavaan_mg_betag1.json  lavaan, the same data, beta g1only,
#                                  99d_lavaan_mg_betag1.R          (optional)
#   fixtures/MANIFEST_lavaan_mg_difftest.json  lavaan, DIFFTEST across generated
#                                  levels, designs E, B, D,
#                                  99e_lavaan_mg_difftest.R        (optional)
#   fixtures/MANIFEST_lavaan_fitindices.json  lavaan scaled fit indices, c1
#                                  baseline, 99f_lavaan_fitindices.R (optional)
#   fixtures/MANIFEST_lavaan_condx.json  lavaan conditional.x Stage 1, one binary x,
#                                  99h_lavaan_condx.R (optional)
#   fixtures/MANIFEST_lavaan_condx_fit.json  lavaan conditional.x full fit,
#                                  99i_lavaan_condx_fit.R (optional)
#   fixtures/MANIFEST_lavaan_mg_fitidx.json  lavaan G = 2 fit indices across levels,
#                                  99j_lavaan_mg_fitidx.R (optional)
#   fixtures/MANIFEST_phi2_corners.json  Phi2 corners, MPFR M1 reference
#                                  certified by M4, 99j_phi2_corners.R (optional)
#   fixtures/MANIFEST_efa_data.json  tier 1, geomin EFA data, Cases F and G,
#                                  91_efa_data.R                   (optional)
#   fixtures/MANIFEST_lavaan_efa_geomin.json  lavaan oblique geomin, E/F/G,
#                                  99h_lavaan_efa_geomin.R         (optional)
#
# Tiers 2 and 3 are optional because they need lavaan and Mplus respectively,
# and the whole point of the tier split is that tier 1 regenerates on any
# machine forever. An absent optional manifest is reported and skipped. A
# present one is checked exactly as strictly, because a silently rewritten
# reference is no less damaging than a rewritten tier-1 one -- it just fails
# later, inside a module comparison, where the cause is harder to see.
#
# What this actually guards against. The fixtures are committed files, so any
# process that rewrites them silently -- git's line-ending conversion on
# Windows checkout being the usual culprit -- would invalidate every
# comparison the Stata test suite makes, without any test failing. Run this
# after a fresh clone, after changing .gitattributes, and in CI. A mismatch
# means either a fixture was regenerated (expected: regenerate the manifest
# too, in the same commit) or something rewrote it behind your back (not
# expected: find out what).
#
# Base R only, like the rest of fixtures/generate/.
# ---------------------------------------------------------------------------
#
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: 2026 Erik Ruzek

check_manifest <- function(manifest, label) {
  lines <- readLines(manifest, warn = FALSE)

  # the files block holds lines of the form:   "fixtures/...": "<32 hex>",
  pat  <- paste0('"(fixtures/[^"]+)"[[:space:]]*:',
                 '[[:space:]]*"([0-9a-f]{32})"')
  m    <- regmatches(lines, regexec(pat, lines))
  hits <- Filter(function(x) length(x) == 3L, m)
  if (!length(hits))
    stop("no file hashes found in ", manifest, " -- is it the expected format?")

  paths <- vapply(hits, `[`, character(1), 2L)
  want  <- vapply(hits, `[`, character(1), 3L)

  missing <- !file.exists(paths)
  got     <- rep(NA_character_, length(paths))
  got[!missing] <- unname(tools::md5sum(paths[!missing]))

  status <- ifelse(missing, "MISSING",
                   ifelse(got == want, "ok", "CHANGED"))

  cat(sprintf("\n%s (%s)\n", label, manifest))
  for (i in seq_along(paths))
    cat(sprintf("  %-7s %s\n", status[i], paths[i]))

  c(checked = length(paths), bad = sum(status != "ok"))
}

tier1 <- file.path("fixtures", "MANIFEST.json")
if (!file.exists(tier1))
  stop("no manifest at ", tier1,
       " -- run fixtures/generate/10_fixtures.R first")

totals <- check_manifest(tier1, "tier 1: base R, regenerable anywhere")

optional <- list(
  list(file = "MANIFEST_lavaan.json", label = "tier 2: lavaan reference",
       gen = "20_lavaan_reference.R"),
  list(file = "MANIFEST_mplus.json", label = "tier 3: Mplus reference",
       gen = "30_mplus_reference.R"),
  # The .dta files are BINARY and committed, so they need this more than the
  # CSVs do, not less: a line-ending rewrite would corrupt them outright rather
  # than merely changing them. .gitattributes marks fixtures/data/** as -text,
  # and this is what proves it held.
  list(file = "MANIFEST_dta.json", label = "transfer channel: .dta fixtures",
       gen = "40_dta_roundtrip.R"),
  # Its own manifest, not a merge into MANIFEST_lavaan.json: 20_ rewrites that
  # one from its own list. Only the condition-level summaries are committed;
  # the replicates themselves are transient (50_replicates.R's header).
  list(file = "MANIFEST_mc.json", label = "lavaan Monte Carlo SE reference",
       gen = "70_lavaan_mc_reference.R (after 50_replicates.R)"),
  # Module 11 D7: lavaan's standardizedSolution() on Case E. Own manifest for
  # the same reason as MANIFEST_mc.json.
  list(file = "MANIFEST_lavaan_std.json", label = "lavaan standardized solution (D7)",
       gen = "80_lavaan_standardized.R"),
  # O6: lavaan on data with a one-empty-cell pair, its own data file included.
  # Own manifest for the same reason as MANIFEST_mc.json.
  list(file = "MANIFEST_lavaan_emptycell.json", label = "lavaan one-empty-cell reference (O6)",
       gen = "90_lavaan_emptycell.R"),
  # Marker-metric start trap (E-010, 2026-10-06): tier-1 data and
  # population values, then lavaan's fits of those data. Own manifests.
  list(file = "MANIFEST_marker_trap.json", label = "tier 1: marker-trap data and population",
       gen = "95_marker_trap.R with DWLSMV_TRAP_FIXTURE=1"),
  list(file = "MANIFEST_lavaan_marker_trap.json", label = "lavaan marker-trap reference",
       gen = "96_lavaan_marker_trap.R (after 95_ in fixture mode)"),
  # DIFFTEST (0.3.0): lavaan's scaled-shifted difference test, Case E pairs.
  # Own manifest for the same reason as MANIFEST_mc.json.
  list(file = "MANIFEST_lavaan_difftest.json", label = "lavaan difference-test reference (DIFFTEST)",
       gen = "97_lavaan_difftest.R"),
  # Multiple groups (0.3.0): design E, two unequal groups, four levels.
  list(file = "MANIFEST_lavaan_multigroup.json", label = "lavaan two-group reference (pooled fit)",
       gen = "98_lavaan_multigroup.R"),
  # Standardized solution at G = 2 (0.3.0, O17): designs E, B, D, B2F, four levels.
  list(file = "MANIFEST_lavaan_mg_std.json", label = "lavaan two-group standardized solution",
       gen = "99k_lavaan_mg_standardized.R"),
  # Multiple groups with a structural path (0.3.0, O2 lift): design M, beta free.
  list(file = "MANIFEST_lavaan_mg_beta.json", label = "lavaan two-group reference with beta",
       gen = "99_lavaan_mg_beta.R"),
  # The same data, beta fixed at 0.4 in both groups (O2 widened 2026-10-11).
  list(file = "MANIFEST_lavaan_mg_betafix.json", label = "lavaan two-group reference, beta fixed",
       gen = "99b_lavaan_mg_betafix.R"),
  # The same data, beta free in group 2 and absent in group 1 (O2 widened 2026-10-11, g2only).
  list(file = "MANIFEST_lavaan_mg_betag2.json", label = "lavaan two-group reference, beta g2only",
       gen = "99c_lavaan_mg_betag2.R"),
  # The same data, beta free in group 1 and absent in group 2 (O2 widened 2026-10-11, g1only).
  list(file = "MANIFEST_lavaan_mg_betag1.json", label = "lavaan two-group reference, beta g1only",
       gen = "99d_lavaan_mg_betag1.R"),
  # DIFFTEST between generated invariance levels (exemption gate, 2026-10-13): designs E, B, D.
  list(file = "MANIFEST_lavaan_mg_difftest.json", label = "lavaan two-group difference test across levels",
       gen = "99e_lavaan_mg_difftest.R"),
  # Scaled fit indices and the c1 baseline (0.3.0, option A, 2026-10-21): case E, design E.
  list(file = "MANIFEST_lavaan_fitindices.json", label = "lavaan scaled fit indices and c1 baseline",
       gen = "99f_lavaan_fitindices.R"),
  # Conditional-covariate mode (v0.4.0, item 4a, 2026-10-24): lavaan conditional.x Stage 1.
  list(file = "MANIFEST_lavaan_condx.json", label = "lavaan conditional.x Stage 1 (s and gamma)",
       gen = "99h_lavaan_condx.R"),
  # Conditional-covariate mode (v0.4.0, item 4b, 2026-10-25): lavaan conditional.x full fit.
  list(file = "MANIFEST_lavaan_condx_fit.json", label = "lavaan conditional.x full fit (estimates, SEs, T3)",
       gen = "99i_lavaan_condx_fit.R"),
  # G > 1 fit indices and their differences (mg-fitidx-t2, 2026-10-31).
  list(file = "MANIFEST_lavaan_mg_fitidx.json", label = "lavaan G = 2 fit indices across invariance levels",
       gen = "99j_lavaan_mg_fitidx.R"),
  # Numerical reference layer (2026-10-28): Phi2 corners, M1 certified by M4.
  list(file = "MANIFEST_phi2_corners.json", label = "Phi2 corners, MPFR M1/M4 reference",
       gen = "99j_phi2_corners.R"),
  # ESEM groundwork: unrotated two-factor EFA on Case E. Own manifest, as above.
  list(file = "MANIFEST_lavaan_efa.json", label = "lavaan unrotated EFA (ESEM step 0)",
       gen = "99g_lavaan_efa.R"),
  # 0.3.0 geomin EFA (Unit 4, 2026-10-29): tier-1 data for Cases F and G, then
  # lavaan's oblique geomin on E, F, G. Own manifests, as above.
  list(file = "MANIFEST_efa_data.json", label = "tier 1: geomin EFA data (Cases F, G)",
       gen = "91_efa_data.R"),
  list(file = "MANIFEST_lavaan_efa_geomin.json", label = "lavaan oblique geomin EFA reference",
       gen = "99h_lavaan_efa_geomin.R")
)

for (t in optional) {
  path <- file.path("fixtures", t$file)
  if (file.exists(path)) {
    totals <- totals + check_manifest(path, t$label)
  } else {
    cat(sprintf("\n%s: %s not present -- skipping.\n", t$label, t$file))
    cat(sprintf("  Generate: Rscript fixtures/generate/%s\n", t$gen))
  }
}

cat(sprintf("\n%d file(s) checked, %d problem(s)\n",
            totals[["checked"]], totals[["bad"]]))

if (totals[["bad"]]) {
  cat("\nIf you regenerated fixtures on purpose, regenerate the manifest\n",
      "in the same commit. If you did not, check .gitattributes: on\n",
      "Windows, git rewrites line endings on checkout unless the fixture\n",
      "paths are marked -text.\n",
      "\n",
      "A lavaan-generated mismatch has one innocent explanation the\n",
      "others do not: lavaan was upgraded. Check the package versions\n",
      "recorded in that file's own manifest -- MANIFEST_lavaan.json for\n",
      "tier 2, MANIFEST_mc.json for the Monte Carlo SE reference --\n",
      "against the installed ones before hunting a bug. A reference output\n",
      "that moved because the package moved is not a test failure, but it\n",
      "does need the manifest regenerated and the move noted.\n",
      sep = "")
  quit(status = 1L)
}

cat("fixtures verified against the manifest\n")
