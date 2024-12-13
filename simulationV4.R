#' 肿瘤临床试验数据模拟函数
#' 
#' @description 
#' 模拟双臂肿瘤临床试验数据，包括患者响应状态、生存时间等关键指标
#' 
#' @param p1 试验组客观缓解率(ORR)
#' @param p2 对照组客观缓解率(ORR)
#' @param alpha 显著性水平
#' @param power 检验效能
#' @param dropout 预计脱落率
#' @param enrollment_period 入组期（月）
#' @param followup_period 随访期（月）
#' @param assessment_interval 疾病评估间隔（月）
#' @param cr_pr_ratio CR/PR的比例向量
#' @param cr_pfs CR患者的中位PFS时间
#' @param pr_pfs PR患者的中位PFS时间
#' @param non_responder_pfs 非响应者的中位PFS时间
#' @param ttr 达到响应的中位时间
#' @param seed 随机数种子
#'
#' @return 包含模拟临床试验数据的数据框
simulate_oncology_trial <- function(
    p1 = 0.5,                  # 试验组ORR
    p2 = 0.4,                  # 对照组ORR
    alpha = 0.05,              # 显著性水平
    power = 0.8,               # 检验效能
    dropout = 0.2,             # 脱落率
    enrollment_period = 12,    # 入组期（月）
    followup_period = 24,      # 随访期（月）
    assessment_interval = 6,   # 评估间隔（月）
    cr_pr_ratio = c(0.3, 0.7), # CR/PR比例
    cr_pfs = 15,               # CR的PFS时间均值
    pr_pfs = 10,               # PR的PFS时间均值
    non_responder_pfs = 5,     # 非响应者的PFS时间均值
    ttr = 2,                   # ttr的均值
    seed = 123                 # 随机种子
) {
  set.seed(seed)
  
  # ============= 1. 样本量计算 =============
  # 计算合并的响应率
  p_pooled <- mean(c(p1, p2))
  
  # 使用标准公式计算每组所需样本量
  n_per_group <- ceiling(
    (qnorm(1-alpha/2) * sqrt(2*p_pooled*(1-p_pooled)) + 
       qnorm(power) * sqrt(p1*(1-p1) + p2*(1-p2)))^2 / (p1-p2)^2
  )
  
  # 考虑脱落率调整样本量
  n_adjusted <- ceiling(n_per_group / (1 - dropout))
  total_n <- 2 * n_adjusted
  
  # ============= 2. 创建基础数据框 =============
  df <- data.frame(
    SubjID = 1:total_n,
    Group = rep(c("Treatment", "Control"), each = n_adjusted),
    stringsAsFactors = FALSE
  )
  
  # ============= 3. 生成初始响应状态 =============
  #' 生成患者响应状态的内部函数
  #' @param n 患者数量
  #' @param p_resp 响应概率
  #' @param cr_pr_ratio CR和PR的比例
  generate_response <- function(n, p_resp, cr_pr_ratio) {
    # 先生成响应/非响应状态
    response <- rbinom(n, 1, p_resp)
    status <- rep("PD", n)
    
    # 对响应者随机分配CR/PR状态
    responders <- which(response == 1)
    if(length(responders) > 0) {
      status[responders] <- sample(
        c("CR", "PR"), 
        length(responders),
        prob = cr_pr_ratio,
        replace = TRUE
      )
    }
    return(status)
  }
  
  # 为试验组和对照组分别生成响应状态
  df$InitialResponse <- ifelse(
    df$Group == "Treatment",
    generate_response(n_adjusted, p1, cr_pr_ratio),
    generate_response(n_adjusted, p2, cr_pr_ratio)
  )
  
  # ============= 4. 定义状态生成函数 =============
  #' 生成患者疾病进展和事件状态的内部函数
  #' @param response_type 响应类型（CR/PR/PD）
  #' @param pfs_mean 预期PFS时间
  #' @param enrollment_time 入组时间
  #' @param ttr_time 达到响应的时间（可选）
  #' @param death_prob 死亡概率
  #' @param ltfu_prob 失访概率
  generate_patient_status <- function(response_type, pfs_mean, enrollment_time, 
                                      ttr_time = NULL,
                                      death_prob = 0.6,    # 调整死亡概率
                                      ltfu_prob = 0.4) {   # 调整失访概率
    # 生成PFS时间
    if(is.null(ttr_time)) {
      pfs_time <- rexp(1, 1/pfs_mean)
    } else {
      pfs_time <- ttr_time + rexp(1, 1/pfs_mean)
      while(pfs_time <= ttr_time) {
        pfs_time <- ttr_time + rexp(1, 1/pfs_mean)
      }
    }
    
    # 调整事件类型概率分布
    event_prob <- runif(1)
    event_type <- if(event_prob < death_prob) {
      "Exited/Died"
    } else if(event_prob < (death_prob + ltfu_prob)) {
      "Lost to Follow-up"
    } else {
      "PD"
    }
    
    return(list(
      pfs_time = pfs_time,
      event_type = event_type
    ))
  }
  
  # ============= 5. 初始化评估时间点 =============
  # 生成所有随访时间点
  assessment_times <- seq(0, followup_period, by = assessment_interval)
  status_cols <- paste0("Month_", assessment_times)
  df[status_cols] <- ""
  
  # 生成入组时间和初始化时间指标
  df$EnrollmentTime <- runif(total_n, 0, enrollment_period)
  df$TTR <- df$DoR <- df$PFS <- NA
  
  # ============= 6. 生成患者状态 =============
  for(i in 1:nrow(df)) {
    response <- df$InitialResponse[i]
    enrollment_time <- df$EnrollmentTime[i]
    
    # 根据初始响应状态处理
    if(response == "PD") {
      # 处理PD患者
      status <- generate_patient_status("PD", non_responder_pfs, enrollment_time)
      df$PFS[i] <- status$pfs_time
      
    } else {
      # 处理CR/PR患者
      ttr_time <- rexp(1, 1/ttr)
      df$TTR[i] <- ttr_time
      
      # 根据响应类型选择PFS时间
      pfs_mean <- if(response == "CR") cr_pfs else pr_pfs
      status <- generate_patient_status(response, pfs_mean, enrollment_time, ttr_time)
      
      df$PFS[i] <- status$pfs_time
      df$DoR[i] <- status$pfs_time - ttr_time
    }
    
    # 填充每个评估时间点的疾病状态
    for(t in assessment_times) {
      col <- paste0("Month_", t)
      time_since_enrollment <- t - enrollment_time
      
      # 根据时间点确定状态
      df[i, col] <- if(time_since_enrollment < 0) {
        "Not Enrolled"
      } else if(response != "PD" && time_since_enrollment < df$TTR[i]) {
        "SD"
      } else if(time_since_enrollment >= df$PFS[i]) {
        status$event_type
      } else {
        if(response == "PD") "SD" else response
      }
    }
  }
  
  # ============= 7. 整理最终数据框 =============
  # 按照指定顺序排列列名
  desired_order <- c("SubjID", "Group", "InitialResponse", "EnrollmentTime", 
                     "TTR", "DoR", "PFS", status_cols)
  df <- df[, desired_order]
  
  return(df)
}
results <- simulate_oncology_trial()

