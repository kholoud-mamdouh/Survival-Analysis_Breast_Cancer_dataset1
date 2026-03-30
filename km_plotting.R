rm(list = ls())
library(survival)
library(survminer)
library(ggplot2)

surv_data = read.csv('survival_data1.csv')
surv_data$Status = ifelse(surv_data$Status == "Dead", 1, 0)

fit <- survfit(Surv(time = OS, event = Status)~group, data = surv_data)
p.value <- surv_pvalue(fit = fit, data = surv_data)

# ---- 12. Plot ----
ggsurvplot(
  fit,
  data = surv_data,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  ggtheme = theme_minimal(),
  palette = c("#E41A1C", "#377EB8"),
  risk.table.height = 0.25
)
