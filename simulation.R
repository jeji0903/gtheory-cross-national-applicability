pkgs <- c("dplyr", "tidyr", "purrr", "lme4", "ggplot2")
invisible(lapply(pkgs, library, character.only = TRUE))


extract_vc <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  get_v <- function(grp) {
    x <- vc$vcov[vc$grp == grp]
    if (length(x) == 0) 0 else as.numeric(x[1])
  }
  list(
    country              = get_v("country"),
    person_within_country = get_v("country:person_id"),
    item                 = get_v("item"),
    country_item         = get_v("country:item"),
    residual             = sigma(fit)^2
  )
}

dstudy_reliability <- function(vc, n_p, n_i) {
  s_c  <- vc$country
  s_p  <- vc$person_within_country
  s_i  <- vc$item
  s_ci <- vc$country_item
  s_e  <- vc$residual

  person_G   <- s_p / (s_p + s_e / n_i)
  person_Phi <- s_p / (s_p + s_i / n_i + s_e / n_i)

  country_G   <- s_c / (s_c + s_ci / n_i + s_p / n_p + s_e / (n_p * n_i))
  country_Phi <- s_c / (s_c + s_i / n_i + s_ci / n_i + s_p / n_p + s_e / (n_p * n_i))

  tibble(Person_G = person_G, Person_Phi = person_Phi,
         Country_G = country_G, Country_Phi = country_Phi)
}

fit_gtheory <- function(dat_long) {
  lmer(y ~ 1 + (1|country) + (1|country:person_id) +
         (1|item) + (1|country:item),
       data = dat_long, REML = TRUE)
}


simulate_dimension <- function(seed,
                               n_country = 58, n_person = 500, n_item = 5,
                               mu = 0,
                               var_country, var_person,
                               var_item, var_cxitem, var_resid,
                               dim_name = "Enjoyment") {
  set.seed(seed)

  countries <- sprintf("C%02d", 1:n_country)
  items     <- sprintf("%s_Item%02d", dim_name, 1:n_item)

  u_c  <- rnorm(n_country, 0, sqrt(var_country)); names(u_c) <- countries
  u_i  <- rnorm(n_item,    0, sqrt(var_item));     names(u_i) <- items

  persons <- expand_grid(country = countries, person = 1:n_person) |>
    mutate(person_id = paste(country, sprintf("P%04d", person), sep = ":"))
  u_p <- rnorm(nrow(persons), 0, sqrt(var_person)); names(u_p) <- persons$person_id

  cx <- expand_grid(country = countries, item = items) |>
    mutate(cx_id = paste(country, item, sep = ":"))
  u_ci <- rnorm(nrow(cx), 0, sqrt(var_cxitem)); names(u_ci) <- cx$cx_id

  expand_grid(country = countries, person = 1:n_person, item = items) |>
    mutate(person_id = paste(country, sprintf("P%04d", person), sep = ":"),
           cx_id     = paste(country, item, sep = ":"),
           y = mu + u_c[country] + u_p[person_id] + u_i[item] +
             u_ci[cx_id] + rnorm(n(), 0, sqrt(var_resid))) |>
    select(country, person_id, item, y) |>
    mutate(across(c(country, person_id, item), factor))
}

SHARED <- list(var_item = 0.09, var_cxitem = 0.15, var_resid = 0.42)

conditions <- tibble(
  condition   = c("A: Strongly applicable",
                  "B: Applicable (baseline)",
                  "C: Borderline",
                  "D: Weakly inapplicable",
                  "E: Strongly inapplicable"),
  var_country = c(0.05, 0.20, 0.40, 0.60, 0.75),
  var_person  = c(0.80, 0.60, 0.40, 0.30, 0.25),
  seed        = c(2024, 2025, 2026, 2027, 2028)
)

cat("Running simulations...\n")

sim_results <- conditions |>
  rowwise() |>
  mutate(
    data = list(simulate_dimension(
      seed = seed, var_country = var_country, var_person = var_person,
      var_item = SHARED$var_item, var_cxitem = SHARED$var_cxitem,
      var_resid = SHARED$var_resid, dim_name = condition)),
    fit = list(fit_gtheory(data)),
    vc  = list(extract_vc(fit))
  ) |>
  ungroup()

cat("All simulations complete.\n")

table1 <- sim_results |>
  rowwise() |>
  mutate(Country = vc$country, `Person (w/in C)` = vc$person_within_country,
         Item = vc$item, `Country × Item` = vc$country_item,
         Residual = vc$residual) |>
  ungroup() |>
  transmute(Condition = condition,
            `Set c:p` = sprintf("%.2f : %.2f", var_country, var_person),
            across(c(Country, `Person (w/in C)`, Item, `Country × Item`, Residual),
                   ~round(.x, 4)))

cat("\n===== Table 1: Variance Components =====\n")
print(table1, n = Inf, width = Inf)


table2 <- sim_results |>
  select(condition, vc) |>
  crossing(StudentsPerCountry = c(100, 250, 500, 1000),
           ItemsPerDimension  = c(3, 5, 7)) |>
  rowwise() |>
  mutate(rel = list(dstudy_reliability(vc, n_p = StudentsPerCountry,
                                       n_i = ItemsPerDimension))) |>
  ungroup() |>
  unnest(rel) |>
  transmute(Condition = condition,
            `Students per Country` = StudentsPerCountry,
            `Items per Dimension`  = ItemsPerDimension,
            across(c(Person_G, Person_Phi, Country_G, Country_Phi),
                   ~round(.x, 4)))

cat("\n===== Table 2: D-study Reliabilities =====\n")
print(table2, n = Inf, width = Inf)


theme_pub <- theme_minimal(base_size = 12, base_family = "serif") +
  theme(plot.title       = element_text(face = "bold", size = 13, hjust = 0),
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text       = element_text(face = "bold", size = 10),
        axis.title       = element_text(size = 11),
        plot.margin      = margin(10, 15, 10, 10))

cond_labels <- c("A: Strongly\napplicable", "B: Applicable\n(baseline)",
                 "C: Borderline", "D: Weakly\ninapplicable",
                 "E: Strongly\ninapplicable")

comp_colors <- c("Country" = "#1565C0", "Person (w/in)" = "#42A5F5",
                 "Item" = "#A5D6A7", "Country×Item" = "#FFCC80",
                 "Residual" = "#BDBDBD")

fig1_data <- sim_results |>
  rowwise() |>
  mutate(Country = vc$country, `Person (w/in)` = vc$person_within_country,
         Item = vc$item, `Country×Item` = vc$country_item,
         Residual = vc$residual) |>
  ungroup() |>
  select(condition, Country, `Person (w/in)`, Item, `Country×Item`, Residual) |>
  pivot_longer(-condition, names_to = "Component", values_to = "Variance") |>
  mutate(condition = factor(condition, levels = conditions$condition),
         Component = factor(Component,
                            levels = c("Residual","Country×Item","Item",
                                       "Person (w/in)","Country")))

fig1 <- ggplot(fig1_data, aes(x = condition, y = Variance, fill = Component)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = comp_colors, guide = guide_legend(reverse = TRUE)) +
  scale_x_discrete(labels = cond_labels) +
  labs(x = NULL, y = "Variance", fill = "Component") +
  theme_pub + theme(axis.text.x = element_text(size = 9, lineheight = 0.9))

ggsave("Figure1_variance_components.pdf", fig1, width = 8, height = 5)

fig2_data <- sim_results |>
  rowwise() |>
  mutate(`Country (σ²c)` = vc$country,
         `Person (σ²p)`  = vc$person_within_country) |>
  ungroup() |>
  select(condition, `Country (σ²c)`, `Person (σ²p)`) |>
  pivot_longer(-condition, names_to = "Source", values_to = "Variance") |>
  mutate(condition = factor(condition, levels = rev(conditions$condition)))

fig2 <- ggplot(fig2_data, aes(x = condition, y = Variance, fill = Source)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Country (σ²c)" = "#1565C0",
                                "Person (σ²p)"  = "#EF5350")) +
  coord_flip() +
  labs(x = NULL, y = "Estimated Variance", fill = NULL) +
  theme_pub + theme(legend.position = "top")

ggsave("Figure2_country_vs_person.pdf", fig2, width = 7, height = 4.5)

fig3_data <- sim_results |>
  select(condition, vc) |>
  crossing(n_p = seq(50, 1500, by = 50)) |>
  rowwise() |>
  mutate(rel = list(dstudy_reliability(vc, n_p = n_p, n_i = 5))) |>
  ungroup() |>
  unnest(rel) |>
  select(condition, n_p, Person_G, Country_G) |>
  pivot_longer(c(Person_G, Country_G),
               names_to = "Level", values_to = "Value") |>
  mutate(condition = factor(condition, levels = conditions$condition),
         Level = gsub("_G$", "", Level))

fig3 <- ggplot(fig3_data, aes(x = n_p, y = Value, color = Level, linetype = Level)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ condition, nrow = 1) +
  scale_color_manual(values = c("Person" = "#EF5350", "Country" = "#1565C0")) +
  scale_linetype_manual(values = c("Person" = "solid", "Country" = "dashed")) +
  scale_y_continuous(limits = c(0.3, 1.0), breaks = seq(0.3, 1.0, 0.1)) +
  scale_x_continuous(breaks = c(250, 750, 1250)) +
  labs(x = "Number of Students per Country", y = "G Coefficient") +
  theme_pub + theme(strip.text = element_text(size = 8),
                    axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

ggsave("Figure3_reliability_by_np.pdf", fig3, width = 12, height = 4.5)

fig4_data <- sim_results |>
  select(condition, vc) |>
  crossing(n_p = seq(50, 1500, by = 50)) |>
  rowwise() |>
  mutate(rel = list(dstudy_reliability(vc, n_p = n_p, n_i = 5))) |>
  ungroup() |>
  unnest(rel) |>
  pivot_longer(c(Person_G, Person_Phi, Country_G, Country_Phi),
               names_to = "Coefficient", values_to = "Value") |>
  mutate(condition = factor(condition, levels = conditions$condition),
         Level = ifelse(grepl("Person", Coefficient), "Person", "Country"),
         Type  = ifelse(grepl("_G$", Coefficient), "G (relative)", "Φ (absolute)"))

fig4 <- ggplot(fig4_data, aes(x = n_p, y = Value, color = Level, linetype = Level)) +
  geom_line(linewidth = 0.8) +
  facet_grid(Type ~ condition) +
  scale_color_manual(values = c("Person" = "#EF5350", "Country" = "#1565C0")) +
  scale_linetype_manual(values = c("Person" = "solid", "Country" = "dashed")) +
  scale_y_continuous(limits = c(0.3, 1.0)) +
  labs(x = "Number of Students per Country", y = "Reliability Coefficient") +
  theme_pub + theme(strip.text = element_text(size = 8),
                    axis.text.x = element_text(size = 7, angle = 45, hjust = 1))

ggsave("Figure4_reliability_G_Phi.pdf", fig4, width = 13, height = 5.5)

fig5_data <- sim_results |>
  rowwise() |>
  mutate(rel = list(dstudy_reliability(vc, n_p = 500, n_i = 5)),
         ratio = var_country / var_person,
         ratio_label = sprintf("%.2f:%.2f", var_country, var_person)) |>
  ungroup() |>
  unnest(rel) |>
  select(condition, ratio, ratio_label, Person_G, Country_G) |>
  pivot_longer(c(Person_G, Country_G),
               names_to = "Coefficient", values_to = "Value") |>
  mutate(Coefficient = recode(Coefficient,
                              Person_G = "Person G", Country_G = "Country G"),
         condition = factor(condition, levels = conditions$condition))

fig5 <- ggplot(fig5_data, aes(x = ratio, y = Value,
                               color = Coefficient, shape = Coefficient)) +
  geom_point(size = 4) + geom_line(linewidth = 1) +
  geom_text(aes(label = condition), size = 2.5, vjust = -1.2,
            show.legend = FALSE, color = "grey30") +
  scale_color_manual(values = c("Person G" = "#EF5350", "Country G" = "#1565C0")) +
  scale_x_continuous(
    breaks = unique(fig5_data$ratio),
    labels = (fig5_data |> distinct(ratio, ratio_label))$ratio_label) +
  scale_y_continuous(limits = c(0.4, 1.0)) +
  labs(x = "Country-to-Person Variance Ratio (σ²c : σ²p)",
       y = "G Coefficient (np = 500, ni = 5)", color = NULL, shape = NULL) +
  theme_pub + theme(legend.position = "top")

ggsave("Figure5_inverse_relationship.pdf", fig5, width = 7, height = 5)

cat("\nDone. All figures saved as PDF.\n")
