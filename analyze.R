analyze_oncology_data <- function(data, timepoint) {
  require(dplyr)
  require(survival)
  require(pander)
  
  # 获取时间点列名
  time_col <- paste0("Month_", timepoint)
  
  # 1. 分析各组各状态人数
  status_counts <- data %>%
    group_by(Group, .data[[time_col]]) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Group, values_from = Count) %>%
    mutate(Total = rowSums(select(., where(is.numeric)), na.rm = TRUE))
  
  # 2. 排除未入组患者
  enrolled_data <- data %>%
    filter(.data[[time_col]] != "Not Enrolled")
  
  # 计算ORR
  orr_data <- enrolled_data %>%
    group_by(Group) %>%
    summarise(
      Total = n(),
      Responders = sum(.data[[time_col]] %in% c("CR", "PR")),
      ORR = paste0(round(Responders/Total * 100, 1), "%")  # 在 ORR 后面加上百分号
    )
  
  
  # Fisher精确检验
  resp_table <- matrix(
    c(orr_data$Responders[1], orr_data$Total[1] - orr_data$Responders[1],
      orr_data$Responders[2], orr_data$Total[2] - orr_data$Responders[2]),
    nrow = 2
  )
  fisher_test <- fisher.test(resp_table)
  
  # 3. 生存分析
  # 确定生存状态
  enrolled_data <- enrolled_data %>%
    mutate(
      # DoR状态
      dor_event = case_when(
        .data[[time_col]] %in% c("CR", "PR") ~ 1,
        TRUE ~ 0
      ),
      # TTR状态
      ttr_event = case_when(
        .data[[time_col]] %in% c("CR", "PR") ~ 1,
        TRUE ~ 0
      ),
      # PFS状态
      pfs_event = case_when(
        .data[[time_col]] %in% c("PD","Exited/Died") ~ 1,
        TRUE ~ 0
      )
    )
  
  # DoR分析 - 仅包含有DoR值的患者
  dor_data <- enrolled_data %>%
    filter(!is.na(DoR))
  
  dor_surv <- Surv(dor_data$DoR, dor_data$dor_event)
  dor_fit <- survfit(dor_surv ~ dor_data$Group)
  dor_cox <- coxph(dor_surv ~ dor_data$Group)
  dor_logrank <- survdiff(dor_surv ~ dor_data$Group)
  
  dor_median <- summary(dor_fit)$table[, "median"]
  dor_ci <- summary(dor_fit)$table[, c("0.95LCL", "0.95UCL")]
  dor_hr <- exp(dor_cox$coefficients)
  dor_hr_ci <- exp(confint(dor_cox))
  dor_p <- 1 - pchisq(dor_logrank$chisq, df = 1)
  
  # TTR分析 - 仅包含有TTR值的患者
  ttr_data <- enrolled_data %>%
    filter(!is.na(TTR))
  
  ttr_surv <- Surv(ttr_data$TTR, ttr_data$ttr_event)
  ttr_fit <- survfit(ttr_surv ~ ttr_data$Group)
  ttr_cox <- coxph(ttr_surv ~ ttr_data$Group)
  ttr_logrank <- survdiff(ttr_surv ~ ttr_data$Group)
  
  ttr_median <- summary(ttr_fit)$table[, "median"]
  ttr_ci <- summary(ttr_fit)$table[, c("0.95LCL", "0.95UCL")]
  ttr_hr <- exp(ttr_cox$coefficients)
  ttr_hr_ci <- exp(confint(ttr_cox))
  ttr_p <- 1 - pchisq(ttr_logrank$chisq, df = 1)
  
  # PFS分析
  pfs_surv <- Surv(enrolled_data$PFS, enrolled_data$pfs_event)
  pfs_fit <- survfit(pfs_surv ~ enrolled_data$Group)
  pfs_cox <- coxph(pfs_surv ~ enrolled_data$Group)
  pfs_logrank <- survdiff(pfs_surv ~ enrolled_data$Group)
  
  pfs_median <- summary(pfs_fit)$table[, "median"]
  pfs_ci <- summary(pfs_fit)$table[, c("0.95LCL", "0.95UCL")]
  pfs_hr <- exp(pfs_cox$coefficients)
  pfs_hr_ci <- exp(confint(pfs_cox))
  pfs_p <- 1 - pchisq(pfs_logrank$chisq, df = 1)
  
  # 使用pander输出结果
  cat("\nStatus：\n")
  pander(status_counts)
  cat("\nObjective Response Rate Analysis：\n")
  pander(orr_data)
  cat("Fisher Exact Test p value:", round(fisher_test$p.value, 4), "\n")
  cat("\nSurvival Analysis Results：\n")
  survival_table <- data.frame(
    Endpoint = c("DoR", "TTR", "PFS"),
    `Treatment Median` = sprintf("%.1f (%.1f-%.1f)", 
                                          c(dor_median[1], ttr_median[1], pfs_median[1]),
                                          c(dor_ci[1,1], ttr_ci[1,1], pfs_ci[1,1]),
                                          c(dor_ci[1,2], ttr_ci[1,2], pfs_ci[1,2])),
    `Control Median` = sprintf("%.1f (%.1f-%.1f)", 
                                        c(dor_median[2], ttr_median[2], pfs_median[2]),
                                        c(dor_ci[2,1], ttr_ci[2,1], pfs_ci[2,1]),
                                        c(dor_ci[2,2], ttr_ci[2,2], pfs_ci[2,2])),
    HR = sprintf("%.2f (%.2f-%.2f)",
                 c(dor_hr, ttr_hr, pfs_hr),
                 c(dor_hr_ci[1], ttr_hr_ci[1], pfs_hr_ci[1]),
                 c(dor_hr_ci[2], ttr_hr_ci[2], pfs_hr_ci[2])),
    `P value` = round(c(dor_p, ttr_p, pfs_p), 4)
  )
  pander(survival_table)
  cat("The P value was calculated using the Log-Rank Test.")
  # 返回结果
  results <- list(
    status_counts = status_counts,
    orr = orr_data,
    fisher_p = fisher_test$p.value,
    dor = list(median = dor_median, ci = dor_ci, hr = dor_hr, hr_ci = dor_hr_ci, p_value = dor_p),
    ttr = list(median = ttr_median, ci = ttr_ci, hr = ttr_hr, hr_ci = ttr_hr_ci, p_value = ttr_p),
    pfs = list(median = pfs_median, ci = pfs_ci, hr = pfs_hr, hr_ci = pfs_hr_ci, p_value = pfs_p)
  )
  
  return(invisible(results))
}
# Perform the analysis
analyze_oncology_data(results, 6)
analyze_oncology_data(results, 12) 
analyze_oncology_data(results, 24) 

