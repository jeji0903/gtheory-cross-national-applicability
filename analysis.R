pkgs <- c("dplyr", "tidyr", "purrr", "lme4", "lavaan", "semTools", "readr")
invisible(lapply(pkgs, library, character.only = TRUE))
source("simulation.R")

items_d1 <- c("st094q01na", "st094q02na", "st094q03na",
              "st094q04na", "st094q05na")
items_d2 <- c("st095q04na", "st095q07na", "st095q08na",
              "st095q13na", "st095q15na")


enjoy_long <- dat2 |>
  select(cntryid, person_id, all_of(items_d1)) |>
  pivot_longer(all_of(items_d1), names_to = "item", values_to = "y") |>
  drop_na() |>
  rename(country = cntryid) |>
  mutate(across(c(country, person_id, item), factor))

interest_long <- dat2 |>
  select(cntryid, person_id, all_of(items_d2)) |>
  pivot_longer(all_of(items_d2), names_to = "item", values_to = "y") |>
  drop_na() |>
  rename(country = cntryid) |>
  mutate(across(c(country, person_id, item), factor))

fit_enjoy    <- fit_gtheory(enjoy_long)
fit_interest <- fit_gtheory(interest_long)

vc_pisa_enjoy    <- extract_vc(fit_enjoy)
vc_pisa_interest <- extract_vc(fit_interest)

cat("\n===== PISA Variance Components: Enjoyment =====\n")
print(vc_pisa_enjoy)
cat("\n===== PISA Variance Components: Interest =====\n")
print(vc_pisa_interest)


table3_pisa <- expand_grid(
  Dimension          = c("Enjoyment", "Interest"),
  StudentsPerCountry = c(100, 250, 500)
) |>
  mutate(
    ItemsPerDimension = 5,
    vc = ifelse(Dimension == "Enjoyment",
                list(vc_pisa_enjoy), list(vc_pisa_interest))
  ) |>
  rowwise() |>
  mutate(rel = list(dstudy_reliability(vc,
                                       n_p = StudentsPerCountry,
                                       n_i = ItemsPerDimension))) |>
  ungroup() |>
  unnest(rel) |>
  transmute(Dimension, `Students per Country` = StudentsPerCountry,
            `Person G` = round(Person_G, 4),
            `Person Φ` = round(Person_Phi, 4),
            `Country G` = round(Country_G, 4),
            `Country Φ` = round(Country_Phi, 4))

cat("\n===== Table 3: PISA D-study Reliabilities =====\n")
print(table3_pisa, n = Inf, width = Inf)

cfa_dat <- dat2 |>
  select(all_of(c("cntryid", items_d1, items_d2))) |>
  drop_na() |>
  mutate(cntryid = as.character(cntryid))

cfa_model <- paste0(
  "Enjoy =~ ",    paste(items_d1, collapse = " + "), "\n",
  "Interest =~ ", paste(items_d2, collapse = " + ")
)

fit_pooled <- cfa(model = cfa_model, data = cfa_dat,
                  ordered = c(items_d1, items_d2), estimator = "WLSMV")

cat("\n===== Pooled CFA Fit =====\n")
print(round(fitMeasures(fit_pooled,
      c("chisq","df","pvalue","cfi","tli","rmsea",
        "rmsea.ci.lower","rmsea.ci.upper","srmr")), 4))

cat("\n===== Standardized Factor Loadings =====\n")
std_sol <- standardizedSolution(fit_pooled)
print(std_sol |> filter(op == "=~") |>
        select(lhs, rhs, est.std, se, pvalue) |>
        mutate(across(c(est.std, se), ~round(.x, 3))))

fit_1f <- cfa(
  model = paste0("F =~ ", paste(c(items_d1, items_d2), collapse = " + ")),
  data = cfa_dat, ordered = c(items_d1, items_d2), estimator = "WLSMV")

cat("\n===== One-Factor CFA Fit =====\n")
print(round(fitMeasures(fit_1f, c("cfi","tli","rmsea","srmr")), 4))

MIN_N <- 200
country_n <- cfa_dat |> count(cntryid, name = "n")
countries_keep <- country_n |> filter(n >= MIN_N) |> pull(cntryid)
cfa_dat_mg <- cfa_dat |> filter(cntryid %in% countries_keep)

cat(sprintf("\nMulti-group invariance: %d countries (n >= %d), N = %s\n",
            length(countries_keep), MIN_N,
            format(nrow(cfa_dat_mg), big.mark = ",")))

fit_config <- cfa(model = cfa_model, data = cfa_dat_mg, group = "cntryid",
                  ordered = c(items_d1, items_d2), estimator = "WLSMV")

fit_thresh <- cfa(model = cfa_model, data = cfa_dat_mg, group = "cntryid",
                  ordered = c(items_d1, items_d2), estimator = "WLSMV",
                  group.equal = "thresholds")

fit_scalar <- cfa(model = cfa_model, data = cfa_dat_mg, group = "cntryid",
                  ordered = c(items_d1, items_d2), estimator = "WLSMV",
                  group.equal = c("thresholds", "loadings"))

fit_strict <- tryCatch(
  cfa(model = cfa_model, data = cfa_dat_mg, group = "cntryid",
      ordered = c(items_d1, items_d2), estimator = "WLSMV",
      group.equal = c("thresholds", "loadings", "residuals")),
  error = function(e) { cat("Strict model did not converge.\n"); NULL }
)

get_fit_row <- function(fit, label) {

  fm <- fitMeasures(fit, c("chisq.scaled","df.scaled","cfi.scaled",
                            "tli.scaled","rmsea.scaled","srmr"))
  tibble(Model = label,
         `χ²` = round(fm["chisq.scaled"], 1),
         df   = round(fm["df.scaled"], 0),
         CFI  = round(fm["cfi.scaled"], 4),
         TLI  = round(fm["tli.scaled"], 4),
         RMSEA = round(fm["rmsea.scaled"], 4),
         SRMR  = round(fm["srmr"], 4))
}

comparison <- bind_rows(
  get_fit_row(fit_config, "1. Configural"),
  get_fit_row(fit_thresh, "2. Threshold"),
  get_fit_row(fit_scalar, "3. Scalar")
)
if (!is.null(fit_strict)) {
  comparison <- bind_rows(comparison, get_fit_row(fit_strict, "4. Strict"))
}

cat("\n===== Measurement Invariance Comparison =====\n")
print(comparison, width = Inf)

# Delta-fit (Chen, 2007 criteria)
cat("\n===== Change in Fit =====\n")
for (i in 2:nrow(comparison)) {
  d_cfi   <- comparison$CFI[i] - comparison$CFI[i-1]
  d_rmsea <- comparison$RMSEA[i] - comparison$RMSEA[i-1]
  cat(sprintf("  %s vs %s: ΔCFI = %+.4f (%s), ΔRMSEA = %+.4f (%s)\n",
              comparison$Model[i], comparison$Model[i-1],
              d_cfi, ifelse(abs(d_cfi) <= 0.010, "pass", "fail"),
              d_rmsea, ifelse(abs(d_rmsea) <= 0.015, "pass", "fail")))
}

cat("\n===== Scaled χ² Difference Tests =====\n")
print(lavTestLRT(fit_config, fit_thresh))
print(lavTestLRT(fit_thresh, fit_scalar))
if (!is.null(fit_strict)) print(lavTestLRT(fit_scalar, fit_strict))

write_csv(std_sol |> filter(op == "=~") |>
            select(Factor = lhs, Item = rhs, Loading = est.std, SE = se, p = pvalue),
          "CFA_factor_loadings.csv")
write_csv(comparison, "CFA_invariance_comparison.csv")

cat("\nDone. Results saved to CSV.\n")
