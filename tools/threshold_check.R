# dwlsmv -- threshold_check.R: are the DWLS/ULS optimum's thresholds Stage 1's
# sample thresholds exactly? lavaan, delta and theta parameterizations, plus a
# shared-threshold negative control. Executed 2026-10-07, R 4.6.0, lavaan 0.7.2;
# figures in E-075. Run from the repository root: Rscript tools/threshold_check.R
suppressPackageStartupMessages(library(lavaan))
R <- "fixtures/data/"   # run from the repository root
mt <- read.csv(paste0(R, "marker_trap.csv")); mt[] <- lapply(mt, ordered)
ce <- read.csv(paste0(R, "case_E.csv"));      ce[] <- lapply(ce, ordered)
q <- sprintf("q%d", 1:8); a <- sprintf("att%d", 1:5)
m_uv <- paste0("MathAb =~ ", paste0("NA*", q, collapse=" + "), "\n",
  "MathAtt =~ ", paste0("NA*", a, collapse=" + "), "\nMathAb ~ MathAtt\n",
  "MathAtt ~~ 1*MathAtt\nMathAb ~~ 1*MathAb\n")
m_mk <- paste0("MathAb =~ ", paste(q, collapse=" + "), "\nMathAtt =~ ",
  paste(a, collapse=" + "), "\nMathAb ~ MathAtt\n")
m_E <- "f1 =~ y1 + y2 + y3\n f2 =~ y4 + y5 + y6\n f2 ~ f1"
run <- function(lbl, model, dat, est, param = "delta") {
  f <- sem(model, data = dat, ordered = names(dat), estimator = est,
           parameterization = param, bad.marker.crit = 0)
  sth <- lavInspect(f, "sampstat")$th; fth <- lavInspect(f, "fitted")$th
  pe <- parameterEstimates(f); pth <- pe$est[pe$op == "|"]
  cat(sprintf("%-10s %-4s %-6s conv=%s nth=%d max|fitted-s|=%.3e max|est-s|=%.3e iter=%d\n",
      lbl, est, param, lavInspect(f, "converged"), length(sth),
      max(abs(fth - sth)), max(abs(pth - as.numeric(sth))),
      lavInspect(f, "iterations")))
}
for (e in c("DWLS", "ULS")) {
  run("trap-uv", m_uv, mt, e); run("trap-mk", m_mk, mt, e); run("caseE", m_E, ce, e)
  run("trap-uv", m_uv, mt, e, "theta"); run("caseE", m_E, ce, e, "theta")
}
# Shared-threshold breaker: equal thresholds across y1,y2 in case E (DWLS)
f <- sem(paste0(m_E, "\n y1 | a*t1\n y2 | a*t1"), data = ce, ordered = names(ce),
         estimator = "DWLS")
cat(sprintf("caseE shared y1|t1=y2|t1: max|fitted-s|=%.3e\n",
    max(abs(lavInspect(f, "fitted")$th - lavInspect(f, "sampstat")$th))))
