# Step 3: measurement invariance of every questionnaire, congenital aphantasia against each reference group.
#
# Per questionnaire: multi-group ordinal confirmatory factor models (WLSMV, delta parameterisation) with the published
# factor structure, specified with semTools::measEq.syntax and the identification of Wu & Estabrook (2016):
#   configural -> thresholds (equal thresholds) -> metric (+ equal loadings) -> scalar (+ equal intercepts, which
#   frees the latent means of the aphantasia group)
# Each step is compared with the one before (scaled CFI may fall by at most .010 and scaled RMSEA rise by at most
# .015; Chen, 2007). If the scalar step fails, the item with the largest score test on its threshold and intercept
# constraints is freed, one item at a time, up to MAXFREE items, and freeing stops as soon as it no longer brings the
# model closer to the criteria -> partial scalar model. The latent mean difference (aphantasia minus reference, in SD
# of the reference group) comes from the final model.
# Comparisons: congenital vs Prolific controls; congenital vs Aphantasia Network members with imagery.
# The Psi-Q is not tested: it defines the groups and is at the floor in aphantasia.
#
# Input : output/scored.csv
# Output: output/03_invariance_fit_<ref>.csv, 03_invariance_freed_<ref>.csv, 03_latent_means_<ref>.csv,
#         03_invariance_summary_<ref>.md
# Run   : Rscript scripts/03_invariance.R [control | an_with_imagery]
suppressPackageStartupMessages({library(lavaan); library(semTools)})
args <- commandArgs(FALSE)
me <- sub("^--file=", "", args[grep("^--file=", args)])
ROOT <- if (length(me)) normalizePath(file.path(dirname(me), ".."), winslash = "/") else getwd()   # repository root
MAXFREE <- 6
d <- read.csv(file.path(ROOT, "output/scored.csv"), check.names = FALSE)
d$grp <- ifelse(d$grp_ext %in% c("cong_floor", "cong_above"), "congenital", d$grp_ext)

f <- function(name, prefix, idx) paste0(name, " =~ ", paste0(prefix, idx, collapse = " + "))
MODELS <- list(
  DES  = f("dis", "DES_Q", 1:28),
  MAIA = paste(f("noticing", "MAIA_", 1:4), f("notdistract", "MAIA_", 5:10), f("notworry", "MAIA_", 11:15),
               f("attention", "MAIA_", 16:22), f("emotional", "MAIA_", 23:27), f("selfreg", "MAIA_", 28:31),
               f("listening", "MAIA_", 32:34), f("trusting", "MAIA_", 35:37), sep = "\n"),
  # IAMI: 10 future and 10 past items, the same split script 01 uses for the subscale scores
  IAMI = paste(f("future", "IAMI_Q", c(1, 3, 5, 6, 7, 10, 12, 15, 16, 18)),
               f("past", "IAMI_Q", c(2, 4, 8, 9, 11, 13, 14, 17, 19, 20)), sep = "\n"),
  CTQ  = paste(f("emoabuse", "CTQ_A_Q", 1:5), f("phyabuse", "CTQ_B_Q", 1:5), f("sexabuse", "CTQ_C_Q", 1:5),
               f("emoneglect", "CTQ_D_Q", 1:5), f("phyneglect", "CTQ_E_Q", 1:5), sep = "\n"),
  PHQ  = f("dep", "PHQ_Q", 1:9),
  GAD  = f("anx", "GAD_Q", 1:7),
  PCL  = paste(f("intrusion", "PCL_Q", 1:5), f("avoidance", "PCL_Q", 6:7), f("negmood", "PCL_Q", 8:14),
               f("arousal", "PCL_Q", 15:20), sep = "\n"))
STEPS <- list(configural = "configural", thresholds = "thresholds", metric = c("thresholds", "loadings"),
              scalar = c("thresholds", "loadings", "intercepts"))

collapse <- function(x, g, minn = 5) {
  # integer-code an ordinal item; DES 0-100 is first binned to 0, 10, 20, 30-40 and 50+ (% of the time); then merge
  # adjacent categories until every category has at least minn responses in both groups
  if (length(unique(x)) > 7) x <- cut(x, c(-Inf, 0, 10, 20, 40, Inf), labels = FALSE)
  x <- as.integer(factor(x))
  repeat {
    tab <- table(factor(x, levels = seq_len(max(x))), g)
    low <- apply(tab, 1, min)
    if (all(low >= minn) || max(x) <= 2) break
    k <- which.min(low)
    if (k == 1) x[x >= 2] <- x[x >= 2] - 1 else x[x >= k] <- x[x >= k] - 1
  }
  x
}

fm <- function(fit) fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.scaled", "rmsea.scaled", "srmr"))
fits <- list(); freed <- list(); means <- list()

args <- commandArgs(TRUE)
REFS <- if (length(args)) args else c("control", "an_with_imagery")   # one comparison per process runs both in parallel
TAG <- paste(REFS, collapse = "_")
for (ref in REFS) for (ins in names(MODELS)) {
  items <- unique(unlist(regmatches(MODELS[[ins]], gregexpr("[A-Z]+_[A-Z_]*Q?[0-9]+", MODELS[[ins]]))))
  s <- d[d$grp %in% c(ref, "congenital"), c("grp", items)]
  s <- s[complete.cases(s), ]
  s$grp <- factor(s$grp, levels = c(ref, "congenital"))
  s <- s[order(s$grp), ]                                  # lavaan numbers groups by first appearance: reference = group 1
  for (i in items) s[[i]] <- collapse(s[[i]], s$grp)
  n <- table(s$grp)
  run <- function(eq, partial = character(0)) tryCatch({
    syn <- as.character(measEq.syntax(configural.model = MODELS[[ins]], data = s, group = "grp", ordered = items,
                                      parameterization = "delta", ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
                                      group.equal = eq, group.partial = if (length(partial)) partial else NULL))
    cfa(syn, data = s, group = "grp", ordered = items, estimator = "WLSMV", parameterization = "delta")
  }, error = function(e) { message(ref, " ", ins, ": ", conditionMessage(e)); NULL })
  steps <- lapply(STEPS, run)
  if (any(vapply(steps, is.null, TRUE))) next
  m <- do.call(rbind, lapply(steps, fm))
  # distance from the criteria for the scalar step (0 = invariant)
  gap <- function(fit) { a <- fm(fit); max(0, -.01 - (a["cfi.scaled"] - m["metric", "cfi.scaled"])) +
                                       max(0, (a["rmsea.scaled"] - m["metric", "rmsea.scaled"]) - .015) }
  part <- steps$scalar; out <- character(0); partial <- character(0)
  while (gap(part) > 0 && length(out) < MAXFREE) {
    pt <- parTable(part)
    uni <- suppressWarnings(lavTestScore(part)$uni)                       # ordinary score tests of the equality constraints
    row <- match(uni$lhs, pt$plabel)
    uni$item <- pt$lhs[row]; uni$op <- pt$op[row]
    keep <- uni$op %in% c("|", "~1") & uni$item %in% setdiff(items, out)
    x2 <- sort(tapply(uni$X2[keep], uni$item[keep], sum), decreasing = TRUE)
    if (!length(x2)) break
    it <- names(x2)[1]
    nth <- sum(pt$op == "|" & pt$lhs == it & pt$group == 1)
    try_partial <- c(partial, paste0(it, "|t", seq_len(nth)), paste0(it, "~1"))
    nxt <- run(STEPS$scalar, try_partial)
    if (is.null(nxt) || gap(nxt) >= gap(part)) break                     # freeing no longer helps
    part <- nxt; out <- c(out, it); partial <- try_partial
  }
  fails <- gap(part) > 0
  if (length(out)) m <- rbind(m, partial_scalar = fm(part))
  m <- data.frame(contrast = ref, instrument = ins, n_ref = n[1], n_cong = n[2], model = rownames(m), m, row.names = NULL)
  m$d_cfi <- c(NA, diff(m$cfi.scaled)); m$d_rmsea <- c(NA, diff(m$rmsea.scaled))
  if (length(out)) {                                                      # partial scalar is compared with metric
    m$d_cfi[5] <- m$cfi.scaled[5] - m$cfi.scaled[3]; m$d_rmsea[5] <- m$rmsea.scaled[5] - m$rmsea.scaled[3]
  }
  fits[[paste(ref, ins)]] <- m
  freed[[paste(ref, ins)]] <- data.frame(contrast = ref, instrument = ins, n_items = length(items),
                                         n_freed = length(out), freed = paste(out, collapse = " "),
                                         still_fails = fails)
  # latent mean of congenital aphantasia; the reference group has mean 0 and variance 1 (ID.fac = "std.lv")
  lm <- function(fit) { e <- lavInspect(fit, "est"); as.numeric(e[[2]]$alpha) }
  means[[paste(ref, ins)]] <- data.frame(contrast = ref, instrument = ins, factor = lavNames(steps$scalar, "lv"),
                                         d_scalar = lm(steps$scalar), d_partial = lm(part))
  cat(ref, ins, "done; freed:", paste(out, collapse = " "), "\n")
  # save after every model, so a stopped run keeps what it finished
  write.csv(do.call(rbind, fits), file.path(ROOT, paste0("output/03_invariance_fit_", TAG, ".csv")), row.names = FALSE)
  write.csv(do.call(rbind, freed), file.path(ROOT, paste0("output/03_invariance_freed_", TAG, ".csv")), row.names = FALSE)
  write.csv(do.call(rbind, means), file.path(ROOT, paste0("output/03_latent_means_", TAG, ".csv")), row.names = FALSE)
}

fits <- do.call(rbind, fits); freed <- do.call(rbind, freed); means <- do.call(rbind, means)
num <- sapply(fits, is.numeric); fits[num] <- round(fits[num], 3); means[4:5] <- round(means[4:5], 2)
sink(file.path(ROOT, paste0("output/03_invariance_summary_", TAG, ".md")))
cat("# 03 - measurement invariance (ordinal multi-group CFA, WLSMV, Wu & Estabrook 2016 identification)\n\n",
    "Criterion per step: dCFI >= -.010 and dRMSEA <= .015.\n\n## Fit\n\n", sep = "")
print(knitr::kable(fits, row.names = FALSE))
cat("\n\n## Items whose thresholds and intercept were freed\n\n"); print(knitr::kable(freed, row.names = FALSE))
cat("\n\n## Latent mean of congenital aphantasia (SD units of the reference group)\n\n"); print(knitr::kable(means, row.names = FALSE))
sink()
