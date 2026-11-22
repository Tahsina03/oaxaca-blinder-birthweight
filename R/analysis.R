# Blinder-Oaxaca decomposition of the birth weight gap between babies
# of mothers who smoked during pregnancy and those who did not.
#
# Data: MASS::birthwt (Baystate Medical Center, Springfield MA, 1986; n = 189)
# Methods: threefold decomposition (endowments, coefficients, interaction),
#          twofold decomposition with pooled coefficients (Jann 2008),
#          detailed contribution of each covariate, bootstrap 95% CIs.
# The decomposition is coded by hand from its formulas, then checked
# against the identity  gap = E + C + I.
#
# Run from the repository root:  Rscript R/analysis.R

suppressPackageStartupMessages({
  library(MASS)
  library(ggplot2)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
d <- MASS::birthwt
d$smoker   <- d$smoke == 1
d$race_blk <- as.integer(d$race == 2)
d$race_oth <- as.integer(d$race == 3)
d$prev_ptl <- as.integer(d$ptl > 0)      # any previous premature labour
d$any_ftv  <- as.integer(d$ftv > 0)      # any physician visit in 1st trimester
d$lwt_10   <- d$lwt / 10                 # mother's weight, per 10 lb
d$age_5    <- d$age / 5                  # mother's age, per 5 years

covars <- c("age_5", "lwt_10", "race_blk", "race_oth", "ht", "ui", "prev_ptl", "any_ftv")
labels <- c(age_5 = "Mother's age", lwt_10 = "Mother's weight", race_blk = "Race: Black",
            race_oth = "Race: other", ht = "Hypertension", ui = "Uterine irritability",
            prev_ptl = "Previous premature labour", any_ftv = "1st-trimester doctor visit")
fml <- as.formula(paste("bwt ~", paste(covars, collapse = " + ")))

desc_vars <- c("bwt", "age", "lwt", covars[-(1:2)])
desc <- data.frame(
  variable   = c("Birth weight (g)", "Mother's age (years)", "Mother's weight (lb)",
                 unname(labels[covars[-(1:2)]])),
  non_smoker = colMeans(d[!d$smoker, desc_vars]),
  smoker     = colMeans(d[d$smoker, desc_vars]),
  row.names = NULL
)
desc[, 2:3] <- round(desc[, 2:3], 2)
write.csv(desc, "output/tables/group_means.csv", row.names = FALSE)
cat("n non-smokers =", sum(!d$smoker), " n smokers =", sum(d$smoker), "\n")
print(desc)

# ---------------------------------------------------------------------------
# 2. Decomposition function
#    A = non-smokers (advantaged group), B = smokers
#    gap = mean(Y_A) - mean(Y_B)
#    Threefold (reference = group B coefficients):
#      E = (Xbar_A - Xbar_B)' beta_B
#      C = Xbar_B' (beta_A - beta_B)
#      I = (Xbar_A - Xbar_B)' (beta_A - beta_B)
#    Twofold (pooled reference beta*, from a pooled model that includes a group dummy):
#      explained   = (Xbar_A - Xbar_B)' beta*
#      unexplained = Xbar_A'(beta_A - beta*) + Xbar_B'(beta* - beta_B)
# ---------------------------------------------------------------------------
decompose <- function(data) {
  A <- data[!data$smoker, ]; B <- data[data$smoker, ]
  fA <- lm(fml, data = A); fB <- lm(fml, data = B)
  pooled <- lm(update(fml, . ~ . + smoker), data = data)
  bA <- coef(fA); bB <- coef(fB)
  bP <- coef(pooled)[names(bA)]
  xA <- c(1, colMeans(A[, covars])); xB <- c(1, colMeans(B[, covars]))
  dx <- xA - xB
  gap <- mean(A$bwt) - mean(B$bwt)

  E <- sum(dx * bB); C <- sum(xB * (bA - bB)); I <- sum(dx * (bA - bB))
  expl <- sum(dx * bP)
  unexpl <- sum(xA * (bA - bP)) + sum(xB * (bP - bB))
  detail_expl <- (dx * bP)[-1]                      # per-covariate explained part

  list(summary = c(gap = gap, endowments = E, coefficients = C, interaction = I,
                   explained = expl, unexplained = unexpl),
       detail = detail_expl)
}

est <- decompose(d)
stopifnot(abs(sum(est$summary[c("endowments", "coefficients", "interaction")]) -
              est$summary["gap"]) < 1e-8,
          abs(sum(est$summary[c("explained", "unexplained")]) - est$summary["gap"]) < 1e-8)

# ---------------------------------------------------------------------------
# 3. Bootstrap confidence intervals (stratified by group)
# ---------------------------------------------------------------------------
B_boot <- 2000
idxA <- which(!d$smoker); idxB <- which(d$smoker)
boot <- replicate(B_boot, {
  bd <- d[c(sample(idxA, replace = TRUE), sample(idxB, replace = TRUE)), ]
  r <- decompose(bd); c(r$summary, r$detail)
})
# A few resamples have a binary covariate with no variation in one group
# (e.g. no hypertensive smokers), so lm() returns NA for it. Those are dropped.
ok <- colSums(is.na(boot)) == 0
cat("Usable bootstrap resamples:", sum(ok), "of", B_boot, "\n")
boot <- boot[, ok]
ci <- t(apply(boot, 1, quantile, c(0.025, 0.975)))

sum_tab <- data.frame(
  component = names(est$summary),
  estimate  = round(est$summary, 1),
  lower95   = round(ci[names(est$summary), 1], 1),
  upper95   = round(ci[names(est$summary), 2], 1),
  pct_of_gap = round(100 * est$summary / est$summary["gap"], 1),
  row.names = NULL
)
write.csv(sum_tab, "output/tables/decomposition_summary.csv", row.names = FALSE)
print(sum_tab)

det_tab <- data.frame(
  covariate = unname(labels[names(est$detail)]),
  explained = round(est$detail, 1),
  lower95   = round(ci[names(est$detail), 1], 1),
  upper95   = round(ci[names(est$detail), 2], 1),
  row.names = NULL
)
det_tab <- det_tab[order(-abs(det_tab$explained)), ]
write.csv(det_tab, "output/tables/detailed_explained.csv", row.names = FALSE)
print(det_tab)

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
three <- sum_tab[sum_tab$component %in% c("endowments", "coefficients", "interaction"), ]
three$component <- factor(three$component, levels = c("endowments", "coefficients", "interaction"),
                          labels = c("Endowments\n(differences in characteristics)",
                                     "Coefficients\n(differences in returns)",
                                     "Interaction"))
p1 <- ggplot(three, aes(component, estimate)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(fill = c("#1b6ca8", "#e67e22", "grey60"), width = 0.6) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.15) +
  labs(x = NULL, y = "Grams",
       title = sprintf("Threefold decomposition of the %.0f g birth weight gap",
                       est$summary["gap"]),
       subtitle = sprintf("Non-smokers minus smokers, with bootstrap 95%% CIs (%s resamples)",
                         format(ncol(boot), big.mark = ","))) +
  theme_minimal(base_size = 12)
ggsave("output/figures/threefold.png", p1, width = 7.5, height = 5, dpi = 150, bg = "white")

det_tab$covariate <- factor(det_tab$covariate, levels = rev(det_tab$covariate))
p2 <- ggplot(det_tab, aes(explained, covariate)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = lower95, xmax = upper95), colour = "#1b6ca8") +
  labs(x = "Contribution to the explained gap (grams)", y = NULL,
       title = "Which characteristics explain part of the gap?",
       subtitle = "Detailed explained component (pooled coefficients), bootstrap 95% CIs") +
  theme_minimal(base_size = 12)
ggsave("output/figures/detailed_explained.png", p2, width = 7.5, height = 5, dpi = 150, bg = "white")

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("\nDone. Results written to output/\n")
