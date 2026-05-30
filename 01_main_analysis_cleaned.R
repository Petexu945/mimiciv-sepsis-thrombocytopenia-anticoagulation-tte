CFG <- list(
  db = list(
    dbname = Sys.getenv("MIMIC_DBNAME", "mimiciv"),
    host = Sys.getenv("MIMIC_HOST", "localhost"),
    port = as.integer(Sys.getenv("MIMIC_PORT", "5432")),
    user = Sys.getenv("MIMIC_USER"),
    password = Sys.getenv("MIMIC_PASSWORD"),
    analysis_schema = Sys.getenv("MIMIC_ANALYSIS_SCHEMA", "anticoagulant")
  ),
  out_dir = Sys.getenv("TTE_OUTPUT_DIR", "tte_outputs_csv"),
  export_backup_csv = TRUE,
  export_balance_plots = TRUE,
  save_analysis_rds = FALSE,
  seed = 2026,
  n_boot = 1000
)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  DBI, RPostgres,
  survival, survminer,
  cobalt, WeightIt, MatchIt,
  tableone, survey,
  broom
)

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
arrange   <- dplyr::arrange

options(dplyr.summarise.inform = FALSE)
set.seed(CFG$seed)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

to_time <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  as.POSIXct(x, tz = "UTC")
}

write_backup_csv <- function(x, file_name, ...) {
  if (isTRUE(CFG$export_backup_csv)) {
    readr::write_csv(x, file.path(CFG$out_dir, file_name), ...)
  }
  invisible(x)
}

ess <- function(w) {
  w <- w[is.finite(w) & !is.na(w)]
  if (length(w) == 0) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

weighted_risk <- function(y, w) {
  ok <- !is.na(y) & !is.na(w) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * y[ok]) / sum(w[ok])
}

cox_robust_row <- function(fit, term = "treat_fProphylaxis") {
  s <- summary(fit)$coefficients
  b <- s[term, "coef"]
  se <- if ("robust se" %in% colnames(s)) s[term, "robust se"] else s[term, "se(coef)"]
  tibble(
    HR  = exp(b),
    lcl = exp(b - 1.96 * se),
    ucl = exp(b + 1.96 * se),
    p   = 2 * pnorm(abs(b / se), lower.tail = FALSE)
  )
}

svy_fit_row <- function(fit, term = "treat_fProphylaxis") {
  cf <- summary(fit)$coefficients
  b  <- cf[term, "Estimate"]
  se <- cf[term, "Std. Error"]
  tibble(
    effect = b,
    lcl = b - 1.96 * se,
    ucl = b + 1.96 * se,
    p = 2 * pnorm(abs(b / se), lower.tail = FALSE)
  )
}

risk_rr_rd_point <- function(dat, y, w_col) {
  w  <- dat[[w_col]]
  tr <- dat$treat

  p1 <- weighted_risk(dat[[y]][tr == 1], w[tr == 1])
  p0 <- weighted_risk(dat[[y]][tr == 0], w[tr == 0])
  rr <- if (is.na(p0) || p0 <= 0) NA_real_ else p1 / p0

  note <- dplyr::case_when(
    is.na(p0) | is.na(p1) ~ "missing risk",
    p0 == 0 & p1 == 0     ~ "no events both groups (RR undefined)",
    p0 == 0 & p1 > 0      ~ "no events in control (RR undefined)",
    TRUE                  ~ NA_character_
  )

  tibble(
    outcome = y, weight = w_col,
    risk_control = p0, risk_prophy = p1,
    RD = p1 - p0, RR = rr, note = note
  )
}

risk_rr_rd_boot <- function(dat, y, w_col, R = CFG$n_boot, seed = CFG$seed) {
  set.seed(seed)
  n <- nrow(dat)

  stat <- replicate(R, {
    idx <- sample.int(n, n, replace = TRUE)
    tmp <- risk_rr_rd_point(dat[idx, , drop = FALSE], y = y, w_col = w_col)
    c(RD = tmp$RD, RR = tmp$RR)
  })

  pt <- risk_rr_rd_point(dat, y = y, w_col = w_col)

  tibble(
    outcome = y, weight = w_col,
    risk_control = pt$risk_control,
    risk_prophy  = pt$risk_prophy,
    RD = pt$RD,
    RD_l = quantile(stat["RD", ], 0.025, na.rm = TRUE),
    RD_u = quantile(stat["RD", ], 0.975, na.rm = TRUE),
    RR = pt$RR,
    RR_l = quantile(stat["RR", ], 0.025, na.rm = TRUE),
    RR_u = quantile(stat["RR", ], 0.975, na.rm = TRUE),
    R = R,
    note = pt$note
  )
}

ncb_point <- function(dat, w_col, w_bleed = 1) {
  w  <- dat[[w_col]]
  tr <- dat$treat

  p_vte_T <- weighted_risk(dat$vte_28_flag[tr == 1], w[tr == 1])
  p_vte_C <- weighted_risk(dat$vte_28_flag[tr == 0], w[tr == 0])
  p_bld_T <- weighted_risk(dat$bleed_28_flag[tr == 1], w[tr == 1])
  p_bld_C <- weighted_risk(dat$bleed_28_flag[tr == 0], w[tr == 0])

  vte_benefit <- p_vte_C - p_vte_T
  bleed_harm  <- p_bld_T - p_bld_C
  ncb         <- vte_benefit - w_bleed * bleed_harm

  tibble(
    weight = w_col,
    w_bleed = w_bleed,
    risk_vte_control   = p_vte_C,
    risk_vte_prophy    = p_vte_T,
    risk_bleed_control = p_bld_C,
    risk_bleed_prophy  = p_bld_T,
    vte_benefit = vte_benefit,
    bleed_harm  = bleed_harm,
    NCB = ncb
  )
}

ncb_boot <- function(dat, w_col, w_bleed = 1, R = CFG$n_boot, seed = CFG$seed) {
  set.seed(seed)
  n <- nrow(dat)

  stat <- replicate(R, {
    idx <- sample.int(n, n, replace = TRUE)
    pt  <- ncb_point(dat[idx, , drop = FALSE], w_col = w_col, w_bleed = w_bleed)
    c(vte_benefit = pt$vte_benefit, bleed_harm = pt$bleed_harm, NCB = pt$NCB)
  })

  pt <- ncb_point(dat, w_col = w_col, w_bleed = w_bleed)

  dplyr::bind_cols(
    pt,
    tibble(
      vte_benefit_l = quantile(stat["vte_benefit", ], 0.025, na.rm = TRUE),
      vte_benefit_u = quantile(stat["vte_benefit", ], 0.975, na.rm = TRUE),

      bleed_harm_l  = quantile(stat["bleed_harm", ],  0.025, na.rm = TRUE),
      bleed_harm_u  = quantile(stat["bleed_harm", ],  0.975, na.rm = TRUE),
      NCB_l         = quantile(stat["NCB", ],         0.025, na.rm = TRUE),
      NCB_u         = quantile(stat["NCB", ],         0.975, na.rm = TRUE),
      R = R
    )
  )
}

wald_joint_p <- function(fit, pattern) {
  b <- coef(fit)
  V <- vcov(fit)
  idx <- grep(pattern, names(b))
  if (length(idx) == 0) return(NA_real_)

  b_int <- b[idx]
  V_int <- V[idx, idx, drop = FALSE]
  if (any(!is.finite(b_int)) || any(!is.finite(V_int))) return(NA_real_)
  if (qr(V_int)$rank < ncol(V_int)) return(NA_real_)

  W <- as.numeric(t(b_int) %*% solve(V_int) %*% b_int)
  pchisq(W, df = length(idx), lower.tail = FALSE)
}

con <- dbConnect(
  RPostgres::Postgres(),
  dbname   = CFG$db$dbname,
  host     = CFG$db$host,
  port     = CFG$db$port,
  user     = CFG$db$user,
  password = CFG$db$password
)

fq_table <- function(table_name) {
  as.character(DBI::dbQuoteIdentifier(con, DBI::Id(schema = CFG$db$analysis_schema, table = table_name)))
}

if (dbExistsTable(con, Id(schema = CFG$db$analysis_schema, table = "analysis_wide_v2"))) {
  df_raw <- dbGetQuery(con, sprintf("
    SELECT *
    FROM %s
    WHERE primary_analysis_flag = 1
      AND exposure_group IN ('A_prophylactic','B_none');
  ", fq_table("analysis_wide_v2")))
  data_source <- paste0(CFG$db$analysis_schema, ".analysis_wide_v2")
} else if (dbExistsTable(con, Id(schema = CFG$db$analysis_schema, table = "analysis_primary"))) {
  df_raw <- dbGetQuery(con, sprintf("
    SELECT *
    FROM %s
    WHERE exposure_group IN ('A_prophylactic','B_none');
  ", fq_table("analysis_primary")))
  data_source <- paste0(CFG$db$analysis_schema, ".analysis_primary")
} else {
  df_raw <- dbGetQuery(con, sprintf("
    SELECT *
    FROM %s
    WHERE primary_analysis_flag = 1
      AND exposure_group IN ('A_prophylactic','B_none');
  ", fq_table("analysis_wide")))
  data_source <- paste0(CFG$db$analysis_schema, ".analysis_wide")
}

df_careunit <- dbGetQuery(con, sprintf("
  SELECT stay_id, first_careunit
  FROM %s;
", fq_table("base_first_stay")))

dbDisconnect(con)

cohort_raw <- df_raw %>%
  left_join(df_careunit, by = "stay_id")

cat(sprintf("读取数据源: %s\n", data_source))
cat(sprintf("主分析集（A/B）样本量: %d\n", nrow(cohort_raw)))
stopifnot(all(cohort_raw$exposure_group %in% c("A_prophylactic", "B_none")))

if ("sofa_24h_t0" %in% names(cohort_raw) && !"sofa_t0" %in% names(cohort_raw)) {
  cohort_raw <- cohort_raw %>% rename(sofa_t0 = sofa_24h_t0)
}
if (!"tL" %in% names(cohort_raw) && "tl" %in% names(cohort_raw)) {
  cohort_raw <- cohort_raw %>% rename(tL = tl)
}

has_t0 <- "t0" %in% names(cohort_raw)
has_tL <- "tL" %in% names(cohort_raw)
has_time0 <- "time0" %in% names(cohort_raw)
has_censor_time <- "censor_time" %in% names(cohort_raw)
has_obs_end_28d <- "obs_end_28d" %in% names(cohort_raw)
has_vte_flag_28d <- "vte_flag_28d" %in% names(cohort_raw)
has_bleed_flag_28d <- "major_bleed_flag_28d" %in% names(cohort_raw)
has_death_28_flag <- "death_28_flag_28d" %in% names(cohort_raw)
has_rbc_any_28d <- "rbc_any_28d" %in% names(cohort_raw)
has_plt_any_28d <- "plt_any_28d" %in% names(cohort_raw)
has_hfd_28 <- "hfd_28" %in% names(cohort_raw)
has_vfd_28 <- "vfd_28" %in% names(cohort_raw)
has_los_28 <- "hosp_days_28d_from_tl" %in% names(cohort_raw)
has_ventdays_28 <- "vent_days_28d" %in% names(cohort_raw)

cohort <- cohort_raw %>%
  mutate(
    t0          = if (has_t0) to_time(t0) else as.POSIXct(NA),
    tL          = if (has_tL) to_time(tL) else as.POSIXct(NA),
    time0       = if (has_time0) to_time(time0) else tL,
    censor_time = if (has_censor_time) to_time(censor_time) else as.POSIXct(NA),
    obs_end_28d = if (has_obs_end_28d) to_time(obs_end_28d) else as.POSIXct(NA),

    treat   = as.integer(treat),
    treat_f = factor(treat, levels = c(0, 1), labels = c("Control", "Prophylaxis")),

    gender          = as.factor(gender),
    cancer_flag     = as.factor(cancer_flag),
    vent_on_t0_flag = as.factor(vent_on_t0_flag),
    rrt_on_t0_flag  = as.factor(rrt_on_t0_flag),

    race_group = case_when(
      is.na(race) ~ "Unknown",
      str_detect(race, regex("WHITE", ignore_case = TRUE)) ~ "White",
      str_detect(race, regex("BLACK|AFRICAN", ignore_case = TRUE)) ~ "Black",
      str_detect(race, regex("ASIAN", ignore_case = TRUE)) ~ "Asian",
      str_detect(race, regex("HISPANIC|LATIN", ignore_case = TRUE)) ~ "Hispanic",
      TRUE ~ "Other"
    ) %>% factor(),

    vte_28_flag   = if (has_vte_flag_28d) as.integer(vte_flag_28d) else as.integer(vte_flag_eff),
    bleed_28_flag = if (has_bleed_flag_28d) as.integer(major_bleed_flag_28d) else as.integer(major_bleed_flag_eff),

    tte_vte_days_28   = clamp(as.numeric(tte_vte_days),   0.001, 28),
    tte_bleed_days_28 = clamp(as.numeric(tte_bleed_days), 0.001, 28),

    death_28_flag = if (has_death_28_flag) as.integer(death_28_flag_28d) else as.integer(death_flag_eff),
    rbc_any_28    = if (has_rbc_any_28d) as.integer(rbc_any_28d) else NA_integer_,
    plt_any_28    = if (has_plt_any_28d) as.integer(plt_any_28d) else NA_integer_,

    hfd_28      = as.numeric(if (has_hfd_28) hfd_28 else NA_real_),
    vfd_28      = as.numeric(if (has_vfd_28) vfd_28 else NA_real_),
    los_28      = as.numeric(if (has_los_28) hosp_days_28d_from_tl else NA_real_),
    ventdays_28 = as.numeric(if (has_ventdays_28) vent_days_28d else NA_real_)
  )

cat("\n=== Landmark / 28d validation ===\n")
if (all(c("time0", "tL") %in% names(cohort))) {
  dt_sec <- abs(as.numeric(difftime(cohort$time0, cohort$tL, units = "secs")))
  if (any(dt_sec > 1, na.rm = TRUE)) {
    bad <- which(dt_sec > 1)[1]
    stop(sprintf(
      "time0 != tL（示例 stay_id=%s, time0=%s, tL=%s）",
      cohort$stay_id[bad], cohort$time0[bad], cohort$tL[bad]
    ))
  }
}
if (all(c("censor_time", "time0") %in% names(cohort))) {
  if (any(cohort$censor_time < cohort$time0, na.rm = TRUE)) {
    bad <- which(cohort$censor_time < cohort$time0)[1]
    stop(sprintf(
      "censor_time < time0（示例 stay_id=%s, censor_time=%s, time0=%s）",
      cohort$stay_id[bad], cohort$censor_time[bad], cohort$time0[bad]
    ))
  }
}
if ("fup_days" %in% names(cohort) && any(cohort$fup_days > 28.001, na.rm = TRUE)) {
  warning("发现 fup_days > 28 天，请复核 SQL 端 28d 行政删失逻辑。")
}
if (any(cohort$tte_vte_days_28 <= 0, na.rm = TRUE))   stop("tte_vte_days_28 <= 0")
if (any(cohort$tte_bleed_days_28 <= 0, na.rm = TRUE)) stop("tte_bleed_days_28 <= 0")
if (any(cohort$vte_28_flag == 1 & cohort$tte_vte_days_28 > 28, na.rm = TRUE)) {
  stop("发现 vte_28_flag=1 但 tte_vte_days_28>28，请复核 SQL。")
}
if (any(cohort$bleed_28_flag == 1 & cohort$tte_bleed_days_28 > 28, na.rm = TRUE)) {
  stop("发现 bleed_28_flag=1 但 tte_bleed_days_28>28，请复核 SQL。")
}

cohort_desc <- cohort %>%
  summarise(
    n = n(),
    n_control = sum(treat == 0, na.rm = TRUE),
    n_prophy  = sum(treat == 1, na.rm = TRUE),
    n_vte     = sum(vte_28_flag == 1, na.rm = TRUE),
    n_bleed   = sum(bleed_28_flag == 1, na.rm = TRUE),
    n_death   = sum(death_28_flag == 1, na.rm = TRUE),
    vte_rate    = mean(vte_28_flag == 1, na.rm = TRUE),
    bleed_rate  = mean(bleed_28_flag == 1, na.rm = TRUE),
    death_rate  = mean(death_28_flag == 1, na.rm = TRUE),
    rbc_any_rate = mean(rbc_any_28 == 1, na.rm = TRUE),
    plt_any_rate = mean(plt_any_28 == 1, na.rm = TRUE)
  )
print(cohort_desc)
write_backup_csv(cohort_desc, "00_desc_primary_set.csv")

vars_to_impute <- c(
  "age","weight_kg","sofa_t0","cci","plt_k_t0","inr","aptt","hb",
  "creat_mgdl","lactate","egfr_ckdepi2021"
)

missing_summary <- cohort %>%
  select(all_of(vars_to_impute)) %>%
  summarise(across(everything(), ~ sum(is.na(.)) / n() * 100)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Missing_Rate_Percent"
  ) %>%
  arrange(desc(Missing_Rate_Percent))

cat("\n=== Variable-level missingness (%) ===\n")
print(missing_summary, n = Inf)

incomplete_rows <- cohort %>%
  select(all_of(vars_to_impute)) %>%
  filter(!complete.cases(.)) %>%
  nrow()

total_rows <- nrow(cohort)
row_missing_rate <- (incomplete_rows / total_rows) * 100

row_missing_summary <- tibble(
  total_rows = total_rows,
  incomplete_rows = incomplete_rows,
  row_missing_rate_percent = row_missing_rate
)

cat("\n=== Row-wise missingness ===\n")
cat(sprintf("Total patients: %d\n", total_rows))
cat(sprintf("Patients with at least one missing value: %d\n", incomplete_rows))
cat(sprintf("Overall row-wise missingness: %.2f%%\n", row_missing_rate))

write_backup_csv(missing_summary, "05_missing_summary_variable_level.csv")
write_backup_csv(row_missing_summary, "06_missing_summary_row_level.csv")

cohort <- cohort %>%
  mutate(across(all_of(vars_to_impute),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .),
                .names = "{.col}_imp"))

cohort <- cohort %>%
  mutate(across(all_of(vars_to_impute),
                ~ ifelse(is.na(.), median(., na.rm = TRUE), .),
                .names = "{.col}_imp"))

ps_vars <- c(
  "age_imp","gender","race_group","weight_kg_imp",
  "sofa_t0_imp","cci_imp","plt_k_t0_imp",
  "inr_imp","aptt_imp","hb_imp","creat_mgdl_imp","lactate_imp",
  "cancer_flag","egfr_ckdepi2021_imp",
  "vent_on_t0_flag","rrt_on_t0_flag"
)

adj_vars_formula <- "age_imp + sofa_t0_imp + plt_k_t0_imp + cancer_flag + inr_imp + egfr_ckdepi2021_imp"
missing_ps <- setdiff(ps_vars, names(cohort))
if (length(missing_ps) > 0) stop("PS variables missing: ", paste(missing_ps, collapse = ", "))

ps_formula <- as.formula(paste("treat_f ~", paste(ps_vars, collapse = " + ")))
ps_model <- glm(ps_formula, data = cohort, family = binomial(link = "logit"))
cohort$ps <- predict(ps_model, type = "response")
cohort$ps <- pmin(pmax(cohort$ps, 1e-6), 1 - 1e-6)

p_treat <- mean(cohort$treat == 1, na.rm = TRUE)
cohort$sw <- ifelse(cohort$treat == 1, p_treat / cohort$ps, (1 - p_treat) / (1 - cohort$ps))

q1  <- quantile(cohort$sw, 0.01, na.rm = TRUE)
q99 <- quantile(cohort$sw, 0.99, na.rm = TRUE)
cohort$sw_trunc <- pmin(pmax(cohort$sw, q1), q99)

weight_diag_ipw <- tibble(
  weight = "IPW_sw_trunc",
  mean = mean(cohort$sw_trunc, na.rm = TRUE),
  sd   = sd(cohort$sw_trunc, na.rm = TRUE),
  min  = min(cohort$sw_trunc, na.rm = TRUE),
  p1   = quantile(cohort$sw_trunc, 0.01, na.rm = TRUE),
  p50  = quantile(cohort$sw_trunc, 0.50, na.rm = TRUE),
  p99  = quantile(cohort$sw_trunc, 0.99, na.rm = TRUE),
  max  = max(cohort$sw_trunc, na.rm = TRUE),
  ESS  = ess(cohort$sw_trunc)
)
print(weight_diag_ipw)
write_backup_csv(weight_diag_ipw, "10_weight_summary_ipw.csv")

bal_ipw <- cobalt::bal.tab(ps_formula, data = cohort, weights = cohort$sw_trunc, un = TRUE, m.threshold = 0.1)
bal_ipw_df <- bal_ipw$Balance %>% as.data.frame() %>% rownames_to_column("variable")
write_backup_csv(bal_ipw_df, "11_balance_ipw_table.csv")

if (isTRUE(CFG$export_balance_plots)) {
  p_love_ipw <- cobalt::love.plot(
    ps_formula, data = cohort, weights = cohort$sw_trunc,
    binary = "std", thresholds = c(m = 0.1), abs = TRUE,
    title = "Balance (IPW stabilized, truncated)"
  )
  png(file.path(CFG$out_dir, "11_loveplot_ipw.png"), width = 8, height = 6, units = "in", res = 600)
  print(p_love_ipw)
  graphics.off()
  svg(file.path(CFG$out_dir, "11_loveplot_ipw.svg"), width = 8, height = 6)
  print(p_love_ipw)
  graphics.off()
}

cox_dat_ipw <- cohort %>%
  filter(is.finite(sw_trunc), sw_trunc > 0) %>%
  drop_na(
    stay_id, treat_f, sw_trunc,
    tte_vte_days_28, vte_28_flag,
    tte_bleed_days_28, bleed_28_flag,
    age_imp, sofa_t0_imp, plt_k_t0_imp, cancer_flag, inr_imp, egfr_ckdepi2021_imp
  )

fit_vte_ipw <- coxph(
  as.formula(paste0("Surv(tte_vte_days_28, vte_28_flag) ~ treat_f + ", adj_vars_formula)),
  data = cox_dat_ipw, weights = sw_trunc, ties = "efron", cluster = stay_id
)
fit_bleed_ipw <- coxph(
  as.formula(paste0("Surv(tte_bleed_days_28, bleed_28_flag) ~ treat_f + ", adj_vars_formula)),
  data = cox_dat_ipw, weights = sw_trunc, ties = "efron", cluster = stay_id
)

res_primary <- bind_rows(
  cox_robust_row(fit_vte_ipw)   %>% mutate(outcome = "VTE (28d)", method = "IPW (robust, cluster)"),
  cox_robust_row(fit_bleed_ipw) %>% mutate(outcome = "ISTH major bleed (28d)", method = "IPW (robust, cluster)")
) %>% select(outcome, method, HR, lcl, ucl, p)
print(res_primary)
write_backup_csv(res_primary, "01_primary_cox_ipw.csv")

cat("\n=== E-value sensitivity analysis (based on primary outcomes) ===\n")

if (!requireNamespace("EValue", quietly = TRUE)) install.packages("EValue")
library(EValue)

get_evalue_row <- function(HR, lcl, ucl) {
  ev <- as.data.frame(EValue::evalues.HR(est = HR, lo = lcl, hi = ucl, rare = FALSE))

  tibble(
    E_value_Estimate = unname(ev["E-values", "point"]),
    E_value_CI_limit = if_else(
      HR < 1,
      unname(ev["E-values", "upper"]),
      unname(ev["E-values", "lower"])
    )
  )
}

evalue_extra <- purrr::pmap_dfr(
  list(res_primary$HR, res_primary$lcl, res_primary$ucl),
  get_evalue_row
)

evalue_res <- bind_cols(
  res_primary %>% select(outcome, method, HR, lcl, ucl),
  evalue_extra
)

print(evalue_res)
write_backup_csv(evalue_res, "50_evalue_sensitivity.csv")

eb_out <- weightit(ps_formula, data = cohort, method = "ebal", estimand = "ATE", moments = 2)
cohort$w_eb <- eb_out$weights

weight_diag_eb <- tibble(
  weight = "EB_w_eb",
  mean = mean(cohort$w_eb, na.rm = TRUE),
  sd   = sd(cohort$w_eb, na.rm = TRUE),
  min  = min(cohort$w_eb, na.rm = TRUE),
  p1   = quantile(cohort$w_eb, 0.01, na.rm = TRUE),
  p50  = quantile(cohort$w_eb, 0.50, na.rm = TRUE),
  p99  = quantile(cohort$w_eb, 0.99, na.rm = TRUE),
  max  = max(cohort$w_eb, na.rm = TRUE),
  ESS  = ess(cohort$w_eb)
)
print(weight_diag_eb)
write_backup_csv(weight_diag_eb, "12_weight_summary_eb.csv")

cox_dat_eb <- cohort %>%
  filter(is.finite(w_eb), w_eb > 0) %>%
  drop_na(
    stay_id, treat_f, w_eb,
    tte_vte_days_28, vte_28_flag,
    tte_bleed_days_28, bleed_28_flag,
    age_imp, sofa_t0_imp, plt_k_t0_imp, cancer_flag, inr_imp, egfr_ckdepi2021_imp
  )

fit_vte_eb <- coxph(
  as.formula(paste0("Surv(tte_vte_days_28, vte_28_flag) ~ treat_f + ", adj_vars_formula)),
  data = cox_dat_eb, weights = w_eb, ties = "efron", cluster = stay_id
)
fit_bleed_eb <- coxph(
  as.formula(paste0("Surv(tte_bleed_days_28, bleed_28_flag) ~ treat_f + ", adj_vars_formula)),
  data = cox_dat_eb, weights = w_eb, ties = "efron", cluster = stay_id
)

res_eb <- bind_rows(
  cox_robust_row(fit_vte_eb)   %>% mutate(outcome = "VTE (28d)", method = "EB (robust, cluster)"),
  cox_robust_row(fit_bleed_eb) %>% mutate(outcome = "ISTH major bleed (28d)", method = "EB (robust, cluster)")
) %>% select(outcome, method, HR, lcl, ucl, p)
print(res_eb)
write_backup_csv(res_eb, "02_primary_cox_eb.csv")

cat("\n=== Sensitivity: PSM 1:1 (manual logit-PS distance) ===\n")

stopifnot("ps" %in% names(cohort))
stopifnot(all(is.finite(cohort$ps)))

eps <- 1e-6
cohort <- cohort %>%
  mutate(
    ps_clamped = pmin(pmax(ps, eps), 1 - eps),
    logit_ps   = qlogis(ps_clamped)
  )
stopifnot(all(is.finite(cohort$logit_ps)))

cal_logit <- 0.2 * sd(cohort$logit_ps, na.rm = TRUE)
cat(sprintf("Calculated logit caliper (0.2*SD): %.4f\n", cal_logit))

m_out <- NULL
cohort_matched <- NULL
res_psm <- tibble(outcome = character(), method = character(), HR = numeric(), lcl = numeric(), ucl = numeric(), p = numeric())

tryCatch({
  m_out <- MatchIt::matchit(
    as.formula(paste("treat_f ~", paste(ps_vars, collapse = " + "))),
    data = cohort,
    method      = "nearest",
    distance    = cohort$logit_ps,
    caliper     = cal_logit,
    std.caliper = FALSE,
    ratio       = 1,
    replace     = FALSE,
    m.order     = "closest"
  )

  bal_psm <- cobalt::bal.tab(m_out, un = TRUE, m.threshold = 0.1)
  bal_psm_df <- bal_psm$Balance %>% as.data.frame() %>% rownames_to_column("variable")
  write_backup_csv(bal_psm_df, "13_balance_psm_table.csv")

  if (isTRUE(CFG$export_balance_plots)) {
    p_love_psm <- cobalt::love.plot(
      m_out, binary = "std", thresholds = c(m = 0.1), abs = TRUE,
      title = "Balance (PSM 1:1, Austin caliper on logit-PS)"
    )
    png(file.path(CFG$out_dir, "13_loveplot_psm.png"), width = 8, height = 6, units = "in", res = 600)
    print(p_love_psm)
    graphics.off()
    svg(file.path(CFG$out_dir, "13_loveplot_psm.svg"), width = 8, height = 6)
    print(p_love_psm)
    graphics.off()
  }

  cohort_matched <- MatchIt::match.data(m_out)
  cat(sprintf("Matched sample: N=%d (Pairs≈%d)\n", nrow(cohort_matched), nrow(cohort_matched)/2))

  fit_vte_psm <- coxph(
    as.formula(paste0("Surv(tte_vte_days_28, vte_28_flag) ~ treat_f + ", adj_vars_formula)),
    data = cohort_matched, ties = "efron", cluster = subclass
  )
  fit_bleed_psm <- coxph(
    as.formula(paste0("Surv(tte_bleed_days_28, bleed_28_flag) ~ treat_f + ", adj_vars_formula)),
    data = cohort_matched, ties = "efron", cluster = subclass
  )

  res_psm <<- bind_rows(
    cox_robust_row(fit_vte_psm)   %>% mutate(outcome = "VTE (28d)", method = "PSM 1:1 (robust, cluster)"),
    cox_robust_row(fit_bleed_psm) %>% mutate(outcome = "ISTH major bleed (28d)", method = "PSM 1:1 (robust, cluster)")
  ) %>% select(outcome, method, HR, lcl, ucl, p)

  print(res_psm)
  write_backup_csv(res_psm, "03_primary_cox_psm.csv")

}, error = function(e) {
  message("PSM step failed: ", e$message)
})

sec_bin_ipw <- bind_rows(
  risk_rr_rd_boot(cohort, "death_28_flag", "sw_trunc", R = CFG$n_boot),
  risk_rr_rd_boot(cohort, "rbc_any_28",    "sw_trunc", R = CFG$n_boot),
  risk_rr_rd_boot(cohort, "plt_any_28",    "sw_trunc", R = CFG$n_boot)
) %>% mutate(method = "IPW standardized risk (bootstrap)")
print(sec_bin_ipw)
write_backup_csv(sec_bin_ipw, "20_secondary_binary_ipw_boot.csv")

sec_bin_eb <- bind_rows(
  risk_rr_rd_boot(cohort, "death_28_flag", "w_eb", R = CFG$n_boot),
  risk_rr_rd_boot(cohort, "rbc_any_28",    "w_eb", R = CFG$n_boot),
  risk_rr_rd_boot(cohort, "plt_any_28",    "w_eb", R = CFG$n_boot)
) %>% mutate(method = "EB standardized risk (bootstrap)")
write_backup_csv(sec_bin_eb, "21_secondary_binary_eb_boot.csv")

ncb_ipw_w1 <- ncb_boot(cohort, w_col = "sw_trunc", w_bleed = 1, R = CFG$n_boot)
ncb_eb_w1  <- ncb_boot(cohort, w_col = "w_eb",     w_bleed = 1, R = CFG$n_boot)

ncb_all <- bind_rows(
  ncb_ipw_w1 %>% mutate(method = "IPW standardized risk (bootstrap)"),
  ncb_eb_w1  %>% mutate(method = "EB standardized risk (bootstrap)")
) %>%
  mutate(outcome = "NCB (28d)") %>%
  select(
    outcome, method, weight, w_bleed,
    risk_vte_control, risk_vte_prophy, risk_bleed_control, risk_bleed_prophy,
    vte_benefit, vte_benefit_l, vte_benefit_u,
    bleed_harm,  bleed_harm_l,  bleed_harm_u,
    NCB, NCB_l, NCB_u, R
  )
print(ncb_all)
write_backup_csv(ncb_all, "23_secondary_ncb_boot.csv")

des_ipw <- svydesign(ids = ~1, weights = ~sw_trunc, data = cohort)
f_cont <- as.formula(paste0(" ~ treat_f + ", adj_vars_formula))

fit_hfd  <- svyglm(update(f_cont, hfd_28      ~ .), design = des_ipw)
fit_vfd  <- svyglm(update(f_cont, vfd_28      ~ .), design = des_ipw)
fit_los  <- svyglm(update(f_cont, los_28      ~ .), design = des_ipw)
fit_vent <- svyglm(update(f_cont, ventdays_28 ~ .), design = des_ipw)

sec_cont_ipw <- bind_rows(
  svy_fit_row(fit_hfd)  %>% mutate(outcome = "HFD-28", method = "IPW svyglm"),
  svy_fit_row(fit_vfd)  %>% mutate(outcome = "VFD-28", method = "IPW svyglm"),
  svy_fit_row(fit_los)  %>% mutate(outcome = "LOS-28 (from tl)", method = "IPW svyglm"),
  svy_fit_row(fit_vent) %>% mutate(outcome = "Vent-days-28", method = "IPW svyglm")
) %>% select(outcome, method, effect, lcl, ucl, p)
print(sec_cont_ipw)
write_backup_csv(sec_cont_ipw, "22_secondary_continuous_ipw_svyglm.csv")

cohort <- cohort %>%
  mutate(
    sub_age  = factor(ifelse(age >= 65, "≥65", "<65"), levels = c("<65", "≥65")),
    sub_plt  = case_when(
      is.na(plt_k_t0) ~ NA_character_,
      plt_k_t0 >= 40 ~ "40,000–50,000/µL",
      TRUE ~ "30,000–39,999/µL"
    ) %>% factor(levels = c("30,000–39,999/µL", "40,000–50,000/µL")),
    sub_inr  = case_when(
      is.na(inr) ~ NA_character_,
      inr > 1.5 ~ ">1.5",
      TRUE ~ "≤1.5"
    ) %>% factor(levels = c("≤1.5", ">1.5")),
    sub_sofa = factor(ifelse(sofa_t0 >= 10, "≥10", "<10"), levels = c("<10", "≥10")),
    sub_egfr = case_when(
      is.na(egfr_ckdepi2021) ~ NA_character_,
      egfr_ckdepi2021 >= 30 ~ "≥30",
      TRUE ~ "<30"
    ) %>% factor(levels = c("<30", "≥30")),
    sub_icu = case_when(
      str_detect(first_careunit, "MICU") ~ "MICU",
      str_detect(first_careunit, "SICU|TSICU") ~ "SICU/TSICU",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("MICU", "SICU/TSICU"))
  )

sub_defs <- tibble(
  var = c("sub_age","sub_plt","sub_inr","sub_sofa","sub_egfr","sub_icu"),
  label = c("年龄","血小板","INR","SOFA","eGFR","ICU类型（MICU vs SICU/TSICU）")
)

run_subgroup <- function(data, time_var, status_var, subgroup_var, subgroup_label, outcome_label) {
  d <- data %>% filter(!is.na(.data[[subgroup_var]]))
  if (subgroup_var == "sub_icu") d <- d %>% filter(sub_icu %in% c("MICU", "SICU/TSICU"))

  if (nrow(d) < 50) {
    p_row <- tibble(outcome = outcome_label, subgroup = subgroup_label, p_interaction = NA_real_, n_used = nrow(d), err = "n < 50", detail = NA_character_)
    hr_rows <- tibble(outcome = outcome_label, subgroup = subgroup_label, level = character(), n = integer(), events = integer(), HR = numeric(), lcl = numeric(), ucl = numeric(), p = numeric(), p_interaction = numeric())
    return(list(p_row = p_row, hr_rows = hr_rows))
  }

  tab <- d %>%
    group_by(.data[[subgroup_var]]) %>%
    summarise(n = n(), events = sum(.data[[status_var]] == 1, na.rm = TRUE), .groups = "drop")
  names(tab)[1] <- "level"

  f_full <- as.formula(paste0(
    "Surv(", time_var, ", ", status_var, ") ~ treat_f * ", subgroup_var, " + ", adj_vars_formula
  ))

  fit_full <- tryCatch(
    coxph(f_full, data = d, weights = sw_trunc, ties = "efron", cluster = stay_id),
    error = function(e) NULL
  )

  p_int <- if (is.null(fit_full)) NA_real_ else wald_joint_p(fit_full, paste0("^treat_fProphylaxis:", subgroup_var))
  detail_txt <- paste(paste(tab$level, tab$n, tab$events, sep = ":n/events="), collapse = " | ")

  hr_rows <- map_dfr(levels(d[[subgroup_var]]), function(lv) {
    dd <- d %>% filter(.data[[subgroup_var]] == lv)
    ne <- sum(dd[[status_var]] == 1, na.rm = TRUE)
    if (nrow(dd) < 30 || ne < 5) {
      return(tibble(level = lv, n = nrow(dd), events = ne, HR = NA_real_, lcl = NA_real_, ucl = NA_real_, p = NA_real_))
    }
    fit_lv <- suppressWarnings(
      coxph(
        as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ treat_f + ", adj_vars_formula)),
        data = dd, weights = sw_trunc, ties = "efron", cluster = stay_id
      )
    )
    rr <- cox_robust_row(fit_lv)
    tibble(level = lv, n = nrow(dd), events = ne, HR = rr$HR, lcl = rr$lcl, ucl = rr$ucl, p = rr$p)
  })

  p_row <- tibble(
    outcome = outcome_label,
    subgroup = subgroup_label,
    p_interaction = p_int,
    n_used = nrow(d),
    err = ifelse(is.na(p_int), "interaction NA / sparse or singular", NA_character_),
    detail = detail_txt
  )

  list(
    p_row = p_row,
    hr_rows = hr_rows %>% mutate(outcome = outcome_label, subgroup = subgroup_label, p_interaction = p_int)
  )
}

vte_sub <- map2(sub_defs$var, sub_defs$label,
                ~run_subgroup(cohort, "tte_vte_days_28", "vte_28_flag", .x, .y, "VTE (28d)"))
bleed_sub <- map2(sub_defs$var, sub_defs$label,
                  ~run_subgroup(cohort, "tte_bleed_days_28", "bleed_28_flag", .x, .y, "ISTH major bleed (28d)"))

sub_p_all <- bind_rows(bind_rows(map(vte_sub, "p_row")), bind_rows(map(bleed_sub, "p_row")))
sub_hr_all <- bind_rows(bind_rows(map(vte_sub, "hr_rows")), bind_rows(map(bleed_sub, "hr_rows")))

print(sub_p_all)
print(sub_hr_all)
write_backup_csv(sub_p_all,  "30_subgroup_interaction_p_ipw.csv")
write_backup_csv(sub_hr_all, "31_subgroup_level_hr_ipw.csv")

baseline_vars <- c(
  "age","gender","race_group","weight_kg","sofa_t0","cci","plt_k_t0",
  "inr","aptt","hb","creat_mgdl","lactate","cancer_flag","egfr_ckdepi2021",
  "vent_on_t0_flag","rrt_on_t0_flag"
)
factor_vars <- c("gender","race_group","cancer_flag","vent_on_t0_flag","rrt_on_t0_flag")

tab_unw <- CreateTableOne(vars = baseline_vars, strata = "treat_f", data = cohort, factorVars = factor_vars)
tab_unw_df <- print(tab_unw, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>% as.data.frame()
write_backup_csv(rownames_to_column(tab_unw_df, "variable"), "40_baseline_table_unweighted.csv")

des_ipw_base <- svydesign(ids = ~1, data = cohort, weights = ~sw_trunc)
tab_w <- svyCreateTableOne(vars = baseline_vars, strata = "treat_f", data = des_ipw_base, factorVars = factor_vars)
tab_w_df <- print(tab_w, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>% as.data.frame()
write_backup_csv(rownames_to_column(tab_w_df, "variable"), "41_baseline_table_ipw.csv")

plot_df <- cohort %>%
  select(
    subject_id, hadm_id, stay_id, t0, tL, time0, censor_time, obs_end_28d,
    exposure_group, treat, treat_f,
    ps, sw_trunc, w_eb,
    tte_vte_days_28, vte_28_flag,
    tte_bleed_days_28, bleed_28_flag,
    death_28_flag, rbc_any_28, plt_any_28,
    hfd_28, vfd_28, los_28, ventdays_28,
    first_careunit,
    sub_age, sub_plt, sub_inr, sub_sofa, sub_egfr, sub_icu
  )
write_backup_csv(plot_df, "00_plot_ready_cohort.csv")

if (!is.null(cohort_matched) && nrow(cohort_matched) > 0) {
  psm_plot_df <- cohort_matched %>%
    select(
      stay_id, subject_id, hadm_id,
      exposure_group, treat, treat_f,
      weights, subclass,
      tte_vte_days_28, vte_28_flag,
      tte_bleed_days_28, bleed_28_flag,
      age, gender, race_group, weight_kg,
      sofa_t0, cci, plt_k_t0, inr, aptt, hb, creat_mgdl, lactate, egfr_ckdepi2021,
      cancer_flag, vent_on_t0_flag, rrt_on_t0_flag,
      first_careunit
    )
  write_backup_csv(psm_plot_df, "04_cohort_matched_plot_data_psm.csv")
}

event_summary <- cohort %>%
  summarise(
    n = n(),
    n_control = sum(treat == 0, na.rm = TRUE),
    n_prophy  = sum(treat == 1, na.rm = TRUE),
    n_vte     = sum(vte_28_flag == 1, na.rm = TRUE),
    n_bleed   = sum(bleed_28_flag == 1, na.rm = TRUE),
    n_death   = sum(death_28_flag == 1, na.rm = TRUE),
    vte_rate  = mean(vte_28_flag == 1, na.rm = TRUE),
    bleed_rate= mean(bleed_28_flag == 1, na.rm = TRUE),
    death_rate= mean(death_28_flag == 1, na.rm = TRUE)
  )

weight_diag_all <- bind_rows(weight_diag_ipw, weight_diag_eb)
write_backup_csv(event_summary, "03_event_summary.csv")
write_backup_csv(weight_diag_all, "04_weight_diagnostics.csv")

final_summary <- bind_rows(
  res_primary %>% mutate(block = "Primary Cox"),
  res_eb      %>% mutate(block = "Primary Cox (EB)"),
  res_psm     %>% mutate(block = "Primary Cox (PSM)"),
  sec_bin_ipw %>% mutate(block = "Secondary binary (IPW)") %>%
    transmute(block, outcome, method, estimate = RD, lcl = RD_l, ucl = RD_u, p = NA_real_, extra = paste0("RR=", RR, " [", RR_l, ", ", RR_u, "]"), note),
  sec_cont_ipw %>% mutate(block = "Secondary continuous (IPW)") %>%
    transmute(block, outcome, method, estimate = effect, lcl, ucl, p, extra = NA_character_, note = NA_character_)
)
write_backup_csv(final_summary, "99_final_summary_all_blocks.csv")

analysis_objects <- list(
  config = CFG,
  data_source = data_source,
  cohort = cohort,
  cohort_desc = cohort_desc,
  ps_model = ps_model,
  bal_ipw = bal_ipw,
  eb_out = eb_out,
  m_out = m_out,
  cohort_matched = cohort_matched,
  res_primary = res_primary,
  evalue_res = evalue_res,
  res_eb = res_eb,
  res_psm = res_psm,
  sec_bin_ipw = sec_bin_ipw,
  sec_bin_eb = sec_bin_eb,
  sec_cont_ipw = sec_cont_ipw,
  ncb_all = ncb_all,
  sub_p_all = sub_p_all,
  sub_hr_all = sub_hr_all,
  tab_unw_df = tab_unw_df,
  tab_w_df = tab_w_df,
  plot_df = plot_df,
  event_summary = event_summary,
  weight_diag_all = weight_diag_all,
  final_summary = final_summary
)

if (isTRUE(CFG$save_analysis_rds)) {
  saveRDS(analysis_objects, file.path(CFG$out_dir, "98_analysis_objects.rds"))
}

cat("\nMain analysis completed.\n")
cat("输出目录：", CFG$out_dir, "\n")

pkgs <- c("tidyverse", "survival", "survminer", "survey", "cobalt",
          "tableone", "patchwork", "scales", "svglite", "DiagrammeR",
          "DiagrammeRsvg", "rsvg", "broom")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

INPUT_DIR <- if (exists("CFG") && !is.null(CFG$out_dir)) CFG$out_dir else "tte_outputs_csv"
OUT_ROOT <- Sys.getenv("TTE_RESULTS_DIR", file.path("results", "publication_outputs"))
FIG_DIR   <- file.path(OUT_ROOT, "figures")
TAB_DIR   <- file.path(OUT_ROOT, "tables")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

read_if_exists <- function(path) {
  if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else NULL
}

obj_if_exists <- function(name, envir = .GlobalEnv) {
  if (exists(name, envir = envir, inherits = TRUE)) get(name, envir = envir, inherits = TRUE) else NULL
}

analysis_rds_path <- file.path(INPUT_DIR, "98_analysis_objects.rds")
analysis_objects_local <- NULL

if (exists("analysis_objects", inherits = TRUE)) {
  analysis_objects_local <- get("analysis_objects", inherits = TRUE)
  message(">>> 检测到当前会话中已有 analysis_objects，优先直接使用对象。")
}

if (is.null(analysis_objects_local) && file.exists(analysis_rds_path)) {
  analysis_objects_local <- readRDS(analysis_rds_path)
  message("Analysis objects loaded.")
}

get_from_pipeline <- function(obj_name = NULL, csv_name = NULL, analysis_list = analysis_objects_local) {

  if (!is.null(obj_name)) {
    x <- obj_if_exists(obj_name)
    if (!is.null(x)) return(x)
  }

  if (!is.null(obj_name) && !is.null(analysis_list) && obj_name %in% names(analysis_list)) {
    return(analysis_list[[obj_name]])
  }

  if (!is.null(csv_name)) {
    x <- read_if_exists(file.path(INPUT_DIR, csv_name))
    if (!is.null(x)) return(x)
  }

  return(NULL)
}

df <- get_from_pipeline(obj_name = "plot_df", csv_name = "00_plot_ready_cohort.csv")

if (is.null(df)) {
  stop("既没有找到 plot_df 对象，也没有找到 00_plot_ready_cohort.csv，请先运行主分析脚本。")
}

df <- df %>%
  mutate(
    treat_f = case_when(
      as.character(treat_f) %in% c("B_none", "Control", "0") ~ "Control",
      as.character(treat_f) %in% c("A_prophylactic", "Prophylaxis", "1") ~ "Prophylaxis",
      TRUE ~ as.character(treat_f)
    ),
    treat_f = factor(treat_f, levels = c("Control", "Prophylaxis")),
    treat   = suppressWarnings(as.numeric(treat))
  )

t_pri <- get_from_pipeline(obj_name = "res_primary",  csv_name = "01_primary_cox_ipw.csv")
t_bin <- get_from_pipeline(obj_name = "sec_bin_ipw",  csv_name = "20_secondary_binary_ipw_boot.csv")
t_con <- get_from_pipeline(obj_name = "sec_cont_ipw", csv_name = "22_secondary_continuous_ipw_svyglm.csv")
sub_hr <- get_from_pipeline(obj_name = "sub_hr_all",  csv_name = "31_subgroup_level_hr_ipw.csv")

BASE_FAMILY <- "Times New Roman"
PAL_GROUP <- c("Control" = "#1F77B4", "Prophylaxis" = "#D62728")

theme_pub <- function(base_size = 18) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      axis.title = element_text(size = base_size),
      axis.text  = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size),
      plot.title   = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle= element_text(size = base_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_gg <- function(p, name, w = 7, h = 4.5,
                    png_dpi = 300, tiff_dpi = 300,
                    keep_svg = FALSE) {

  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".pdf")),
    plot = p,
    width = w, height = h, units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".png")),
    plot = p,
    width = w, height = h, units = "in",
    dpi = png_dpi,
    device = "png",
    bg = "white"
  )

  ggsave(
    filename = file.path(FIG_DIR, paste0(name, ".tiff")),
    plot = p,
    width = w, height = h, units = "in",
    dpi = tiff_dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  if (isTRUE(keep_svg)) {
    ggsave(
      filename = file.path(FIG_DIR, paste0(name, ".svg")),
      plot = p,
      width = w, height = h, units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  }
}

save_dot <- function(dot, name, w_in = 8, h_in = 9, dpi = 600) {
  g <- DiagrammeR::grViz(dot)
  svg_txt <- DiagrammeRsvg::export_svg(g)
  writeLines(svg_txt, con = file.path(FIG_DIR, paste0(name, ".svg")), useBytes = TRUE)
  rsvg::rsvg_png(charToRaw(svg_txt), file = file.path(FIG_DIR, paste0(name, ".png")), width = w_in * dpi, height = h_in * dpi)
}

message(">>> 模块 0 运行完毕。已优先对接主分析对象；若对象不存在，则自动 fallback 到 CSV。")

save_km_pair <- function(g1, g2, name, w = 12, h = 6.8) {
  p_combined <- (
    (g1$plot | g2$plot) /
      (g1$table | g2$table)
  ) + patchwork::plot_layout(heights = c(3.6, 1.1))

  save_gg(p_combined, name, w = w, h = h)
}

make_km <- function(dat, time_col, event_col, weight_col = NULL,
                    title = "", subtitle = "") {

  d <- dat %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(treat_f)
    )

  form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ treat_f"))

  sf <- if (is.null(weight_col)) {
    survival::survfit(form, data = d)
  } else {
    survival::survfit(form, data = d, weights = d[[weight_col]])
  }

  sf$call$formula <- form

  g <- suppressWarnings(
    survminer::ggsurvplot(
      fit = sf,
      data = d,
      conf.int = TRUE,
      censor = TRUE,
      risk.table = TRUE,
      risk.table.title = "Number at risk",
      risk.table.height = 0.22,
      risk.table.y.text = TRUE,
      risk.table.y.text.col = FALSE,
      legend.title = "",
      legend.labs = c("Control", "Prophylaxis"),
      xlim = c(0, 28),
      break.time.by = 7
    )
  )

  g$plot <- g$plot +
    scale_color_manual(values = PAL_GROUP) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Days",
      y = "Survival probability"
    ) +
    theme_pub()

  g$table <- g$table +
    survminer::theme_cleantable(base_size = 12, base_family = BASE_FAMILY) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    )

  return(g)
}

p3a <- make_km(df, "tte_vte_days_28", "vte_28_flag",
               weight_col = NULL,
               title = "Venous thromboembolism",
               subtitle = "Unweighted")

p3b <- make_km(df, "tte_vte_days_28", "vte_28_flag",
               weight_col = "sw_trunc",
               title = "Venous thromboembolism",
               subtitle = "IPW-weighted")

save_km_pair(p3a, p3b, "Figure3_VTE_KM", w = 12, h = 6.8)

p3sa <- make_km(df, "tte_bleed_days_28", "bleed_28_flag",
                weight_col = NULL,
                title = "ISTH major bleeding",
                subtitle = "Unweighted")

p3sb <- make_km(df, "tte_bleed_days_28", "bleed_28_flag",
                weight_col = "sw_trunc",
                title = "ISTH major bleeding",
                subtitle = "IPW-weighted")

save_km_pair(p3sa, p3sb, "Figure3S_Bleed_KM", w = 12, h = 6.8)

message(">>> 模块 2 运行完毕。KM曲线已保存。")

if (!requireNamespace("forestploter", quietly = TRUE)) install.packages("forestploter")
library(forestploter)
library(grid)

if (is.null(sub_hr)) {
  sub_hr <- get_from_pipeline(obj_name = "sub_hr_all", csv_name = "31_subgroup_level_hr_ipw.csv")
}

if (!is.null(sub_hr)) {

  dat_clean <- sub_hr
  if("estimate" %in% names(dat_clean) & !"HR" %in% names(dat_clean)) dat_clean <- rename(dat_clean, HR = estimate)
  if("conf.low" %in% names(dat_clean) & !"lcl" %in% names(dat_clean)) dat_clean <- rename(dat_clean, lcl = conf.low)
  if("conf.high" %in% names(dat_clean) & !"ucl" %in% names(dat_clean)) dat_clean <- rename(dat_clean, ucl = conf.high)
  if("p_interaction" %in% names(dat_clean) & !"p_int" %in% names(dat_clean)) dat_clean <- rename(dat_clean, p_int = p_interaction)

  dat_clean <- dat_clean %>%
    mutate(
      subgroup = case_when(
        subgroup == "年龄" ~ "Age",
        subgroup == "血小板" ~ "Platelet count",
        subgroup == "INR" ~ "INR",
        subgroup == "SOFA" ~ "SOFA score",
        subgroup == "eGFR" ~ "eGFR",
        stringr::str_detect(subgroup, "ICU类型") ~ "ICU type",
        TRUE ~ subgroup
      )
    )

  make_publication_forest <- function(data_subset, file_name, x_lims = c(0.1, 10)) {

    df_fmt <- data_subset %>%
      mutate(
        HR_CI_Str = sprintf("%.2f (%.2f-%.2f)", HR, lcl, ucl),
        P_Int_Str = ifelse(is.na(p_int), "", sprintf("%.3f", p_int))
      )

    plot_data <- tibble()
    subgroups <- unique(df_fmt$subgroup)

    for (sg in subgroups) {
      group_data <- df_fmt %>% filter(subgroup == sg)
      p_val <- group_data$P_Int_Str[1]

      header_row <- tibble(
        Subgroup = sg,
        ` ` = paste(rep(" ", 30), collapse = ""),
        `HR (95% CI)` = "",
        `P-interaction` = p_val,
        HR = NA_real_, lcl = NA_real_, ucl = NA_real_,
        is_summary = TRUE
      )
      plot_data <- bind_rows(plot_data, header_row)

      data_rows <- group_data %>%
        transmute(
          Subgroup = paste0("    ", level),
          ` ` = "",
          `HR (95% CI)` = HR_CI_Str,
          `P-interaction` = "",
          HR = HR, lcl = lcl, ucl = ucl,
          is_summary = FALSE
        )
      plot_data <- bind_rows(plot_data, data_rows)
    }

    tm <- forest_theme(
      base_size = 10,

      core = list(
        bg_params = list(fill = c("#F5F5F5", "#FFFFFF")),
        align = c("l", "c", "c", "c"),

        padding = unit(c(16, 4), "mm")
      ),

      colhead = list(
        align = c("l", "c", "c", "c"),

        padding = unit(c(2, 6), "mm")
      ),

      summary_col = "black",
      summary_fill = "transparent",
      refline_col = "#A00000",
      refline_lty = "solid",
      refline_lwd = 1.2,
      ci_col = "black",
      ci_fill = "black",
      ci_alpha = 1,
      ci_lty = 1,
      ci_lwd = 1.2,
      ci_Theight = 0.2
    )

    p <- forest(
      data = plot_data[, 1:4],
      est = plot_data$HR,
      lower = plot_data$lcl,
      upper = plot_data$ucl,
      ci_column = 2,
      ref_line = 1,
      xlim = x_lims,
      ticks_at = c(0.1, 0.5, 1, 2, 5, 10),
      x_trans = "log10",
      arrow_lab = c("Favors Control", "Favors Prophylaxis"),
      is_summary = plot_data$is_summary,
      theme = tm
    )

    plot_height <- 1.5 + nrow(plot_data) * 0.25

    save_gg(p, file_name, w = 7.5, h = plot_height)}

  vte_dat <- dat_clean %>% filter(stringr::str_detect(outcome, "(?i)vte"))
  if(nrow(vte_dat) > 0) make_publication_forest(vte_dat, "Figure5_Forest_VTE_Publication", c(0.1, 10))

  bld_dat <- dat_clean %>% filter(stringr::str_detect(outcome, "(?i)bleed"))
  if(nrow(bld_dat) > 0) make_publication_forest(bld_dat, "Figure5S_Forest_Bleed_Publication", c(0.1, 10))

  message(">>> 模块 5 (极简顶刊版) 运行完毕。请去 figures 文件夹查看最新结果。")
}

des_ipw <- survey::svydesign(ids = ~1, weights = ~sw_trunc, data = df)

fmt_mean_sd <- function(m, s) sprintf("%.2f ± %.2f", m, s)
fmt_est_ci  <- function(est, lcl, ucl) sprintf("%.2f (%.2f, %.2f)", est, lcl, ucl)
fmt_p       <- function(p) if(is.na(p)) "NA" else if(p < 0.001) "<0.001" else sprintf("%.3f", p)

calc_pri <- function(o_var, t_var, o_name, m_name) {
  c_n <- sum(df[[o_var]][df$treat==0]==1, na.rm=T); c_N <- sum(df$treat==0, na.rm=T)
  t_n <- sum(df[[o_var]][df$treat==1]==1, na.rm=T); t_N <- sum(df$treat==1, na.rm=T)
  s_unw <- summary(survival::coxph(as.formula(paste0("Surv(",t_var,",",o_var,")~treat_f")), data=df))

  m_c <- survey::svymean(as.formula(paste0("~",o_var)), subset(des_ipw, treat==0), na.rm=T)[1]
  m_t <- survey::svymean(as.formula(paste0("~",o_var)), subset(des_ipw, treat==1), na.rm=T)[1]
  r_ipw <- t_pri %>% filter(outcome == m_name)

  tibble(Outcome = o_name, Type = "Primary Time-to-Event Outcomes (HR)",
         Unw_Ctrl = sprintf("%d/%d (%.3f)", c_n, c_N, c_n/c_N), Unw_Prophy = sprintf("%d/%d (%.3f)", t_n, t_N, t_n/t_N),
         Unw_Eff = fmt_est_ci(s_unw$conf.int[1,1], s_unw$conf.int[1,3], s_unw$conf.int[1,4]), Unw_P = fmt_p(s_unw$coefficients[1,5]),
         IPW_Ctrl = sprintf("%.3f", m_c), IPW_Prophy = sprintf("%.3f", m_t),
         IPW_Eff = if(nrow(r_ipw)>0) fmt_est_ci(r_ipw$HR, r_ipw$lcl, r_ipw$ucl) else "NA", IPW_P = if(nrow(r_ipw)>0) fmt_p(r_ipw$p) else "NA")
}

calc_bin <- function(o_var, o_name, m_name) {
  c_n <- sum(df[[o_var]][df$treat==0]==1, na.rm=T); c_N <- sum(df$treat==0, na.rm=T)
  t_n <- sum(df[[o_var]][df$treat==1]==1, na.rm=T); t_N <- sum(df$treat==1, na.rm=T)
  pt <- prop.test(c(t_n, c_n), c(t_N, c_N), correct=F)

  r_ipw <- t_bin %>% filter(outcome == m_name)
  ipw_p <- "NA"
  if(nrow(r_ipw)>0) {
    se <- (r_ipw$RD_u - r_ipw$RD_l) / (2 * 1.96)
    if(!is.na(se) && se > 0) ipw_p <- fmt_p(2 * pnorm(abs(r_ipw$RD / se), lower.tail = F))
  }

  tibble(Outcome = o_name, Type = "Secondary Dichotomous Outcomes (RD)",
         Unw_Ctrl = sprintf("%d/%d (%.3f)", c_n, c_N, c_n/c_N), Unw_Prophy = sprintf("%d/%d (%.3f)", t_n, t_N, t_n/t_N),
         Unw_Eff = sprintf("%.3f (%.3f, %.3f)", (t_n/t_N)-(c_n/c_N), pt$conf.int[1], pt$conf.int[2]), Unw_P = fmt_p(pt$p.value),
         IPW_Ctrl = if(nrow(r_ipw)>0) sprintf("%.3f", r_ipw$risk_control) else "NA",
         IPW_Prophy = if(nrow(r_ipw)>0) sprintf("%.3f", r_ipw$risk_prophy) else "NA",
         IPW_Eff = if(nrow(r_ipw)>0) sprintf("%.3f (%.3f, %.3f)", r_ipw$RD, r_ipw$RD_l, r_ipw$RD_u) else "NA", IPW_P = ipw_p)
}

calc_con <- function(o_var, o_name, m_name) {
  d_c <- df[[o_var]][df$treat==0]; d_t <- df[[o_var]][df$treat==1]
  tt <- t.test(d_t, d_c)

  m_c <- survey::svymean(as.formula(paste0("~",o_var)), subset(des_ipw, treat==0), na.rm=T)[1]
  m_t <- survey::svymean(as.formula(paste0("~",o_var)), subset(des_ipw, treat==1), na.rm=T)[1]
  r_ipw <- t_con %>% filter(outcome == m_name)

  tibble(Outcome = o_name, Type = "Secondary Continuous Outcomes (MD)",
         Unw_Ctrl = fmt_mean_sd(mean(d_c,na.rm=T), sd(d_c,na.rm=T)), Unw_Prophy = fmt_mean_sd(mean(d_t,na.rm=T), sd(d_t,na.rm=T)),
         Unw_Eff = sprintf("%.2f (%.2f, %.2f)", diff(tt$estimate), tt$conf.int[1], tt$conf.int[2]), Unw_P = fmt_p(tt$p.value),
         IPW_Ctrl = sprintf("%.2f", m_c), IPW_Prophy = sprintf("%.2f", m_t),
         IPW_Eff = if(nrow(r_ipw)>0) fmt_est_ci(r_ipw$effect, r_ipw$lcl, r_ipw$ucl) else "NA", IPW_P = if(nrow(r_ipw)>0) fmt_p(r_ipw$p) else "NA")
}

table2 <- bind_rows(
  tryCatch(calc_pri("vte_28_flag", "tte_vte_days_28", "Venous thromboembolism", "VTE (28d)"), error=function(e) NULL),
  tryCatch(calc_pri("bleed_28_flag", "tte_bleed_days_28", "ISTH major bleeding", "ISTH major bleed (28d)"), error=function(e) NULL),
  tryCatch(calc_bin("death_28_flag", "All-cause mortality", "death_28_flag"), error=function(e) NULL),
  tryCatch(calc_bin("rbc_any_28", "Red blood cell transfusion", "rbc_any_28"), error=function(e) NULL),
  tryCatch(calc_bin("plt_any_28", "Platelet transfusion", "plt_any_28"), error=function(e) NULL),
  tryCatch(calc_con("hfd_28", "Hospital-free days", "HFD-28"), error=function(e) NULL),
  tryCatch(calc_con("vfd_28", "Ventilator-free days", "VFD-28"), error=function(e) NULL),
  tryCatch(calc_con("los_28", "Hospital length of stay, days", "LOS-28 (from tl)"), error=function(e) NULL),
  tryCatch(calc_con("ventdays_28", "Duration of mechanical ventilation, days", "Vent-days-28"), error=function(e) NULL)
)

readr::write_excel_csv(table2, file.path(TAB_DIR, "Table2_Formatted_Final.csv"))

message(">>> 模块 6 运行完毕！完美的 Table 2 已经导出到 tables 文件夹。")

cat("\n>>> 开始生成 Table S1 (EBW) 与 Table S2 (PSM)...\n")

df_full <- get_from_pipeline(obj_name = "cohort")

if (is.null(df_full)) {
  stop("未能获取完整的 cohort 对象！请确认主分析输出的 98_analysis_objects.rds 文件在对应目录下。因为 CSV 版本的绘图数据剔除了基线变量。")
}

df_full <- df_full %>%
  mutate(
    treat_f = case_when(
      as.character(treat_f) %in% c("B_none", "Control", "0") ~ "Control",
      as.character(treat_f) %in% c("A_prophylactic", "Prophylaxis", "1") ~ "Prophylaxis",
      TRUE ~ as.character(treat_f)
    ),
    treat_f = factor(treat_f, levels = c("Control", "Prophylaxis"))
  )

df_psm <- get_from_pipeline(obj_name = "cohort_matched", csv_name = "04_cohort_matched_plot_data_psm.csv")

if (is.null(df_psm)) {
  stop("未能获取到 PSM 匹配后的数据集，请确保主分析输出了 matched_cohort。")
}

df_psm <- df_psm %>%
  mutate(
    treat_f = case_when(
      as.character(treat_f) %in% c("B_none", "Control", "0") ~ "Control",
      as.character(treat_f) %in% c("A_prophylactic", "Prophylaxis", "1") ~ "Prophylaxis",
      TRUE ~ as.character(treat_f)
    ),
    treat_f = factor(treat_f, levels = c("Control", "Prophylaxis"))
  )

baseline_vars <- c(
  "age","gender","race_group","weight_kg","sofa_t0","cci","plt_k_t0",
  "inr","aptt","hb","creat_mgdl","lactate","cancer_flag","egfr_ckdepi2021",
  "vent_on_t0_flag","rrt_on_t0_flag"
)
factor_vars <- c("gender","race_group","cancer_flag","vent_on_t0_flag","rrt_on_t0_flag")

des_eb <- survey::svydesign(ids = ~1, data = df_full, weights = ~w_eb)
tab_eb <- tableone::svyCreateTableOne(vars = baseline_vars, strata = "treat_f", data = des_eb, factorVars = factor_vars)
tab_eb_df <- print(tab_eb, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>% as.data.frame() %>% rownames_to_column("variable")
readr::write_excel_csv(tab_eb_df, file.path(TAB_DIR, "TableS1_Baseline_EBW.csv"))

tab_psm <- tableone::CreateTableOne(vars = baseline_vars, strata = "treat_f", data = df_psm, factorVars = factor_vars)
tab_psm_df <- print(tab_psm, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>% as.data.frame() %>% rownames_to_column("variable")
readr::write_excel_csv(tab_psm_df, file.path(TAB_DIR, "TableS2_Baseline_PSM.csv"))

cat(">>> Table S1 & Table S2 已输出至 tables 文件夹。\n")

cat("\n>>> 开始生成 Figure S1 (合并 Love Plot, 使用插补后变量)...\n")

plot_vars_imp <- c(
  "age_imp", "gender", "race_group", "weight_kg_imp",
  "sofa_t0_imp", "cci_imp", "plt_k_t0_imp",
  "inr_imp", "aptt_imp", "hb_imp", "creat_mgdl_imp", "lactate_imp",
  "cancer_flag", "egfr_ckdepi2021_imp",
  "vent_on_t0_flag", "rrt_on_t0_flag"
)

clean_var_names <- data.frame(
  old = c(
    "age_imp", "gender", "race_group", "weight_kg_imp",
    "sofa_t0_imp", "cci_imp", "plt_k_t0_imp",
    "inr_imp", "aptt_imp", "hb_imp", "creat_mgdl_imp", "lactate_imp",
    "cancer_flag", "egfr_ckdepi2021_imp",
    "vent_on_t0_flag", "rrt_on_t0_flag"
  ),
  new = c(
    "Age", "Gender", "Race", "Weight",
    "SOFA Score", "CCI", "Platelets",
    "INR", "APTT", "Hemoglobin", "Creatinine", "Lactate",
    "Cancer", "eGFR",
    "Ventilation at t0", "RRT at t0"
  )
)

df_full <- df_full %>%
  mutate(
    w_psm = ifelse(stay_id %in% df_psm$stay_id, 1, 0)
  )

covs_formula_imp <- as.formula(paste("treat_f ~", paste(plot_vars_imp, collapse = " + ")))

p_love_merged <- cobalt::love.plot(
  covs_formula_imp,
  data = df_full,
  weights = list(
    IPW = df_full$sw_trunc,
    PSM = df_full$w_psm,
    EBW = df_full$w_eb
  ),
  drop.distance = TRUE,
  binary = "std",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  var.names = clean_var_names,
  colors = c("#8C8C8C", "#1F77B4", "#2CA02C", "#D62728"),
  shapes = c(16, 17, 15, 18),
  title = "Covariate Balance Across Causal Inference Frameworks",
  sample.names = c("Unadjusted", "IPW", "PSM", "EBW")
) +
  theme_pub(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(
  filename = file.path(FIG_DIR, "FigureS1_Merged_LovePlot.pdf"),
  plot = p_love_merged,
  width = 9, height = 7, units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(FIG_DIR, "FigureS1_Merged_LovePlot.tiff"),
  plot = p_love_merged,
  width = 9, height = 7, units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

cat(">>> Figure S1 (插补后变量合并 Love Plot) 已输出至 figures 文件夹。\n")

cat("\n>>> 开始生成 Figure S2: 4-Panel KM Curves (PSM & EBW)...\n")

if (!"weights" %in% names(df_psm)) df_psm$weights <- 1

pA <- make_km(df_psm, "tte_vte_days_28", "vte_28_flag", weight_col = "weights",
              title = "Panel A: VTE (PSM)", subtitle = "Propensity Score Matching (1:1)")

pC <- make_km(df_psm, "tte_bleed_days_28", "bleed_28_flag", weight_col = "weights",
              title = "Panel C: Major Bleeding (PSM)", subtitle = "Propensity Score Matching (1:1)")

pB <- make_km(df_full, "tte_vte_days_28", "vte_28_flag", weight_col = "w_eb",
              title = "Panel B: VTE (EBW)", subtitle = "Entropy Balancing Weighting")

pD <- make_km(df_full, "tte_bleed_days_28", "bleed_28_flag", weight_col = "w_eb",
              title = "Panel D: Major Bleeding (EBW)", subtitle = "Entropy Balancing Weighting")

library(patchwork)

p_4panel <- (
  (pA$plot | pB$plot) /
    (pA$table | pB$table) /
    (pC$plot | pD$plot) /
    (pC$table | pD$table)
) + plot_layout(heights = c(3, 1, 3, 1))

ggsave(
  filename = file.path(FIG_DIR, "FigureS2_4Panel_KM_Sensitivity.pdf"),
  plot = p_4panel,
  width = 16, height = 14, units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(FIG_DIR, "FigureS2_4Panel_KM_Sensitivity.tiff"),
  plot = p_4panel,
  width = 16, height = 14, units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white"
)

cat(">>> 补充 KM 曲线 4-Panel 拼图已输出至 figures 文件夹。\n")
