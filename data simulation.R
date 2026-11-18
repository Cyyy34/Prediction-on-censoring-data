library(readxl)
library(dplyr)
library(tidyr)
library(copula)
library(VineCopula)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(reshape2)
library(pROC)  # 用于ROC曲线

####################create generated dataset###########################################################################################################
# 读取数据
df <- read_excel("C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/cleaned_dataset_with_indicators.xlsx")

# 仅选择 numeric columns
df_numeric <- df %>% select(where(is.numeric))

# 计算 numeric 变量的缺失率（从少到多排序）
missing_rate_numeric <- df_numeric %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "missing_rate") %>%
  arrange(missing_rate)

# 选出缺失最少的前 20 个 numeric variables
top20_numeric <- missing_rate_numeric %>% slice(1:20)

# 查看结果
print(top20_numeric)

vars <- c(
  "U3CC07","U3CD07","Height","BMI_V3","BMI","SFLT1_V1","fbHCG_V1","V2BA01_KG",
  "Age_at_V1","ADAM12_V1","AFP_V1","PLGF_V1","INHIBINA_V1","BMI_V2","U2BB02",
  "PAPPA_V1","V3BA01_KG","SFLT1_V2","PAPPA_V2","fbHCG_V2","AFP_V2","PLGF_V2",
  "ENDOGLIN_V1","INHIBINA_V2","U3BB02","ADAM12_V2","ENDOGLIN_V2","U2CD07",
  "U2CC07","VEGF_V2","VITD_MCG"
)

# 依次绘制直方图
for (v in vars) {
  if (v %in% names(df)) {
    if (is.numeric(df[[v]])) {
      print(
        ggplot(df, aes(x = .data[[v]])) +
          geom_histogram(bins = 30, fill = "skyblue", color = "black") +
          ggtitle(paste("Histogram of", v)) +
          theme_minimal()
      )
    } else {
      message(v, " is not numeric — skipped.")
    }
  } else {
    message(v, " is not found in dataset — skipped.")
  }
}

vars_top15 <- c("U3CC07","U3CD07","Height","BMI","SFLT1_V1","fbHCG_V1","V2BA01_KG",
  "Age_at_V1","ADAM12_V1","AFP_V1","PLGF_V1","INHIBINA_V1","U2BB02",
  "PAPPA_V1","SFLT1_V2")

vars_top15

#提取子数据集，并去掉有缺失的行（fit copula 时常用 complete cases）
dat15 <- df %>%
  select(all_of(vars_top15))

dat15_cc <- dat15[complete.cases(dat15), ]

# 如果 complete cases 太少，也可以只对部分变量做，或者先做简单插补

# 转换到 (0,1) 上的 pseudo-observations （经验分布 -> U(0,1)）
U <- pobs(as.matrix(dat15_cc))   # copula 包里的函数

# 拟合 t-copula
d <- ncol(U)
t.cop.start <- tCopula(dim = d, dispstr = "un")  # "un" = 一般相关结构
# ... (前面的代码保持不变，直到 fit.tcop 那一行) ...

# 确保数据转换正确
# 注意：如果 complete cases 很少，后续会报错。建议先检查 nrow(dat15_cc)
if(nrow(dat15_cc) < 50) warning("样本量可能太少，建议检查缺失值处理")

# ==========================================
# 修改部分开始
# ==========================================

# --- 方案：快速拟合 t-Copula (推荐用于保留原模型结构) ---
message("开始拟合 t-Copula (使用 itau.mpl 方法)...")

# 模拟数据 
# 1. 重新定义模拟样本量（通常设为和原始数据 dat15_cc 行数一致）
n_sim <- nrow(dat15_cc)  

# 2. 再次运行模拟
set.seed(123)
U_sim <- rCopula(n_sim, tcop)

# 3. 检查是否成功
dim(U_sim) 
# 应该显示: [1] 1756   15   (假设你的 dat15_cc 是 1756 行)

# method = "itau.mpl" 策略：
# 1. 用 Kendall's Tau 估计相关矩阵 (itau)
# 2. 用 最大似然 估计自由度 (mpl)
# 这将把 106 个参数的优化问题简化为 1 个参数的优化问题
fit.tcop <- fitCopula(t.cop.start, data = U, method = "itau.mpl")

summary(fit.tcop) # 查看拟合结果
tcop <- fit.tcop@copula


# ==========================================
# 修改部分结束
# ==========================================

# 反变换回原始尺度：用经验分布做 inverse CDF
sim_mat <- sapply(seq_len(d), function(j) {
  x <- dat15_cc[[j]]
  # 用经验分布的分位数作为 inverse CDF
  quantile(x, probs = U_sim[, j], na.rm = TRUE, type = 8)
})

sim_dat15 <- as.data.frame(sim_mat)
colnames(sim_dat15) <- vars_top15

# 查看结果
head(sim_dat15)

# 可视化对比：对比原始数据和模拟数据的某个变量（例如 BMI）
par(mfrow=c(1,2))
hist(dat15_cc$BMI, main="Original BMI", col="skyblue", breaks=20)
hist(sim_dat15$BMI, main="Simulated BMI", col="lightgreen", breaks=20)
par(mfrow=c(1,1))

# ==========================================
# 1. 小修正：清理行名
# ==========================================
# 现在的行名是百分数（因为是 quantile 函数生成的），这里重置一下
rownames(sim_dat15) <- NULL

# ==========================================
# 2. 核心验证：检查相关性结构是否保留
# ==========================================

# 计算原始数据和模拟数据的相关系数矩阵 (使用 Spearman 相关系数，更能捕捉非线性关系)
cor_orig <- cor(dat15_cc, method = "spearman")
cor_sim  <- cor(sim_dat15, method = "spearman")

# 计算两者之差
cor_diff <- cor_sim - cor_orig

# --- 可视化对比 A: 相关性矩阵热图 ---
# 这一步能让你一眼看出模拟数据是否“跑偏”了
# 如果颜色很深（红或蓝），说明模拟的数据丢失了相关性
# 如果颜色很浅（接近白色），说明模拟非常成功

library(corrplot)

par(mfrow = c(1, 2))
corrplot(cor_orig, method = "color", type = "upper", 
         tl.col = "black", tl.cex = 0.6, title = "Original Correlation", mar=c(0,0,1,0))
corrplot(cor_sim, method = "color", type = "upper", 
         tl.col = "black", tl.cex = 0.6, title = "Simulated Correlation", mar=c(0,0,1,0))

# --- 可视化对比 B: 差值矩阵 ---
# 专门画出“差异”，看看误差都在哪里
par(mfrow = c(1, 1))
corrplot(cor_diff, method = "color", type = "upper", 
         tl.col = "black", tl.cex = 0.6, 
         title = "Difference (Simulated - Original)", mar=c(0,0,2,0),
         cl.lim = c(-0.2, 0.2)) # 限制颜色范围，突出显示差异大的地方

# 打印平均误差
mean_abs_diff <- mean(abs(cor_diff[upper.tri(cor_diff)]))
message("平均相关系数绝对误差: ", round(mean_abs_diff, 4))

# 设置输出文件名为 comparison_histograms.pdf
pdf("comparison_histograms.pdf", width = 10, height = 6)

# 循环画图
par(mfrow = c(1, 2)) # 一页放左右两张图

for (v in colnames(sim_dat15)) {
  
  # 获取当前变量的数据范围，保证两张图坐标轴一致，方便对比
  rng <- range(c(dat15_cc[[v]], sim_dat15[[v]]), na.rm = TRUE)
  
  # 画原始数据
  hist(dat15_cc[[v]], 
       main = paste("Original:", v), 
       xlab = v, col = "skyblue", 
       xlim = rng, breaks = 30)
  
  # 画模拟数据
  hist(sim_dat15[[v]], 
       main = paste("Simulated:", v), 
       xlab = v, col = "lightgreen", 
       xlim = rng, breaks = 30)
}

# 关闭设备，完成保存
dev.off()

message("图片已保存到你的工作目录下的 comparison_histograms.pdf")

# 定义两组变量名
group_top10 <- c(
  "Height",       # 正态分布 (Normal)
  "U2BB02",       # 正态分布 (Normal, 尖峰不同)
  "Age_at_V1",    # 均匀分布 (Uniform-like)
  "U3CC07",       # 极度右偏 (Extreme Right Skew, 零膨胀)
  "PAPPA_V1",     # 长尾分布 (Long Tail)
  "BMI",          # 对数正态 (Log-normal)
  "SFLT1_V1",     # 偏移偏态 (Shifted Skew)
  "V2BA01_KG",    # 宽跨度右偏 (Wide Span Skew)
  "fbHCG_V1",     # 伽马分布 (Gamma-like)
  "INHIBINA_V1"   # 多峰/不规则 (Irregular)
)

group_remaining5 <- c(
  "U3CD07",       # 重复的极度右偏
  "ADAM12_V1",    # 离散/稀疏右偏
  "AFP_V1",       # 标准右偏
  "PLGF_V1",      # 重复的偏移偏态
  "SFLT1_V2"      # 重复的生物标志物形态
)

# 将它们放入一个 List 中，方便管理
selected_vars_list <- list(
  top10 = group_top10,
  remaining5 = group_remaining5
)

# 查看内容
print(selected_vars_list)

# --- 使用示例 ---
# 如果你想直接从模拟数据中提取这 Top 10 的数据：
# sim_top10_data <- sim_dat15 %>% select(all_of(selected_vars_list$top10))

# ============================================================
# 1. 准备数据
# ============================================================
# 提取之前选出的 Top 10 变量
X_inf <- sim_dat15 %>% select(all_of(selected_vars_list$top10))

# **重要**：标准化数据 (Z-score)，确保 Beta 的大小能真实反映“信号强度”
# 如果不标准化，PAPPA_V1 (数值上千) 的影响会彻底淹没 Height (数值100多)
X_scaled <- scale(X_inf) 

# ============================================================
# 2. 定义 "Signal Gradient" (预设的真实系数 Beta)
# ============================================================
# 假设我们按顺序分配以下信号强度：
# 1-2: 强正相关 (Strong)
# 3-5: 中正相关 (Moderate)
# 6-7: 弱正相关 (Weak)
# 8-10: 负相关 (Negative - 混合强度)

# 你可以根据实际生物学意义调整这个向量
betas <- c(
  # Strong (e.g., log-odds ratio ~ 1.0 - 1.5)
  1.5,  1.2,   
  # Moderate (e.g., ~ 0.5 - 0.8)
  0.8,  0.6,  0.5,
  # Weak (e.g., ~ 0.1 - 0.3)
  0.2,  0.1,
  # Negative (e.g., ~ -0.5 to -1.0)
  -0.5, -0.8, -1.2
)

# 给系数命名，方便查看
names(betas) <- colnames(X_scaled)
print("预设的真实系数 (True Betas):")
print(betas)

# 计算线性预测部分 (X * Beta)
# 矩阵乘法: [n x 10] * [10 x 1] = [n x 1]
linear_pred_no_intercept <- X_scaled %*% betas

# ============================================================
# 3. 寻找截距 beta0 以匹配目标事件率 (32.94%)
# ============================================================
target_prevalence <- 0.3294

# 定义一个函数：给定 beta0，计算当前的平均预测概率与目标的差值
calc_prevalence_diff <- function(beta0) {
  # Logit 公式: p = 1 / (1 + exp(-(beta0 + XB)))
  probs <- plogis(beta0 + linear_pred_no_intercept)
  return(mean(probs) - target_prevalence)
}

# 使用 uniroot 在区间内寻找让差值为 0 的 beta0
# 这里的区间 c(-10, 10) 通常足够覆盖 logit 截距的范围
opt_res <- uniroot(calc_prevalence_diff, interval = c(-10, 10))
beta0_final <- opt_res$root

message("调整后的截距 (beta0): ", round(beta0_final, 4))

# ============================================================
# 4. 生成最终概率和二分类结果
# ============================================================
# 计算最终概率向量 p_i
true_probs <- plogis(beta0_final + linear_pred_no_intercept)

# 伯努利采样生成 Y (0 或 1)
set.seed(2025) # 设定种子保证可复现
Y_sim <- rbinom(n = length(true_probs), size = 1, prob = true_probs)

# ============================================================
# 5. 合并结果并验证
# ============================================================
# 将生成的 Y 和真实概率 p 拼回到模拟数据中
sim_dataset_final <- sim_dat15 %>%
  mutate(
    true_prob = as.vector(true_probs), # 记录上帝视角的真实概率
    outcome = Y_sim                    # 最终生成的二分类 Target
  )

# 验证 1: 检查事件率是否接近 32.94%
actual_rate <- mean(sim_dataset_final$outcome)
message("模拟数据的最终事件率: ", round(actual_rate * 100, 2), "% (目标: 32.94%)")

# 验证 2: 简单看一下几个变量跟 Outcome 的关系 (Boxplot)
# 选取一个“强正相关”变量和一个“负相关”变量看看趋势是否正确
par(mfrow=c(1,2))
boxplot(sim_dataset_final[[names(betas)[1]]] ~ sim_dataset_final$outcome, 
        main=paste("Strong Pos:", names(betas)[1]), col="lightblue")
boxplot(sim_dataset_final[[names(betas)[10]]] ~ sim_dataset_final$outcome, 
        main=paste("Strong Neg:", names(betas)[10]), col="salmon")
par(mfrow=c(1,1))

# 将数据保存到当前工作目录
write.csv(sim_dataset_final, "simulated_dataset_final.csv", row.names = FALSE)

message("文件已保存为 simulated_dataset_final.csv")


###############数据集基本信息#######################
df <- read_csv("C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/simulated_dataset_final.csv")
summary(df)

################################make the truncated datasets###############
# 读取数据
df <- read.csv("C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/simulated_dataset_final.csv")

# 定义需要进行人为截断的 8 个变量
vars_to_truncate <- c(
  # --- Informative (5个) ---
  "BMI",          # Log-normal, Weak Signal
  "PAPPA_V1",     # Long-tail, Moderate Signal
  "fbHCG_V1",     # Gamma-like, Negative Signal
  "Age_at_V1",    # Uniform-like, Moderate Signal
  "Height",       # Normal, Strong Signal
  
  # --- Remaining/Noise (3个) ---
  "U3CD07",       # Extreme Skew
  "PLGF_V1",      # Shifted Skew (Likely correlated with SFLT1)
  "ADAM12_V1"     # Discrete/Sparse
)

# 验证一下这些变量是否都在数据集里
all(vars_to_truncate %in% names(sim_dataset_final)) 
# 应该返回 TRUE

print(vars_to_truncate)

library(dplyr)
library(tidyr)

# 0. 准备工作
# -------------------------------------------------------
# 原始数据
df_orig <- df 

# 之前选定的 8 个变量
vars_to_truncate <- c(
  "BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", 
  "U3CD07", "PLGF_V1", "ADAM12_V1"
)

# 定义截断比例
truncation_levels <- c(0.10, 0.20, 0.30, 0.40)

# 用于存储结果的 List，以及存储 Cutoff 信息的 Data Frame
datasets_list <- list()
cutoff_info_list <- list()

# 1. 循环生成 4 种截断数据集
# -------------------------------------------------------
for (p in truncation_levels) {
  
  # 复制一份原始数据
  df_current <- df_orig
  
  # 用于收集当前比例下的 cutoff 信息
  current_cutoffs <- data.frame(
    Variable = character(),
    Truncation_Level = numeric(),
    Lower_Cutoff = numeric(),
    Upper_Cutoff = numeric(),
    stringsAsFactors = FALSE
  )
  
  message(paste0("正在处理截断比例: Left ", p*100, "% / Right ", p*100, "% ..."))
  
  for (v in vars_to_truncate) {
    
    # 获取原始向量
    x <- df_orig[[v]]
    
    # 计算分位数 (Cutoff Points)
    # 左侧 p% 分位点，右侧 (1-p)% 分位点
    q_lower <- quantile(x, probs = p, na.rm = TRUE)
    q_upper <- quantile(x, probs = 1 - p, na.rm = TRUE)
    
    # 记录 Cutoff 值
    current_cutoffs <- rbind(current_cutoffs, data.frame(
      Variable = v,
      Truncation_Level = p,
      Lower_Cutoff = q_lower,
      Upper_Cutoff = q_upper
    ))
    
    # --- 创建指示列 (Indicator Column) ---
    # 命名规则：变量名 + "_ind" (例如 BMI_ind)
    # 0 = Observed (未缺失)
    # 1 = Below Lower Cutoff (左侧截断)
    # 2 = Above Upper Cutoff (右侧截断)
    # (或者用字符标记，这里用字符更直观)
    
    ind_col_name <- paste0(v, "_ind")
    
    # 初始化为 "Observed"
    indicator_vec <- rep("Observed", nrow(df_current))
    
    # 标记低于下界的
    indicator_vec[x < q_lower] <- "Below_Lower"
    
    # 标记高于上界的
    indicator_vec[x > q_upper] <- "Above_Upper"
    
    # 将指示列加入数据框
    df_current[[ind_col_name]] <- indicator_vec
    
    # --- 实施截断 (设为 NA) ---
    # 将超出范围的值设为 NA
    df_current[[v]][x < q_lower | x > q_upper] <- NA
  }
  
  # 保存当前截断比例的数据集到 List
  dataset_name <- paste0("df_trunc_", p * 100)
  datasets_list[[dataset_name]] <- df_current
  
  # 保存 Cutoff 信息
  cutoff_info_list[[dataset_name]] <- current_cutoffs
}

# 2. 提取结果到独立对象 (可选)
# -------------------------------------------------------
df_trunc_10 <- datasets_list[["df_trunc_10"]]
df_trunc_20 <- datasets_list[["df_trunc_20"]]
df_trunc_30 <- datasets_list[["df_trunc_30"]]
df_trunc_40 <- datasets_list[["df_trunc_40"]]

# 合并所有的 Cutoff 信息为一个大表，方便查看
all_cutoffs <- bind_rows(cutoff_info_list)

# 3. 检查与验证
# -------------------------------------------------------

# 查看 Cutoff 表
print("截断阈值记录表 (部分展示):")
print(head(all_cutoffs, 10))

# 验证 df_trunc_10 中的 BMI 缺失率
# 理论上应该是 10% (左) + 10% (右) = 20% 总缺失
missing_rate_bmi_10 <- mean(is.na(df_trunc_10$BMI))
message("df_trunc_10 中 BMI 的实际缺失率: ", missing_rate_bmi_10 * 100, "% (预期约 20%)")

# 查看指示列 (Indicator) 是否正确生成
print("查看 df_trunc_10 中 BMI 及其指示列的前 10 行:")
print(df_trunc_10 %>% select(BMI, BMI_ind) %>% head(10))

# 4. 保存这些数据集 (可选)
# 1. 创建一个新文件夹来存放这些文件 (如果文件夹不存在的话)
output_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 2. 保存 4 个截断数据集
write.csv(df_trunc_10, file.path(output_dir, "dataset_trunc_10.csv"), row.names = FALSE)
write.csv(df_trunc_20, file.path(output_dir, "dataset_trunc_20.csv"), row.names = FALSE)
write.csv(df_trunc_30, file.path(output_dir, "dataset_trunc_30.csv"), row.names = FALSE)
write.csv(df_trunc_40, file.path(output_dir, "dataset_trunc_40.csv"), row.names = FALSE)

# 3. 保存截断阈值记录表 (Cutoff Info) - 这个非常重要，后续分析要用
write.csv(all_cutoffs, file.path(output_dir, "cutoff_info.csv"), row.names = FALSE)

message("所有文件已成功保存到文件夹: ", file.path(getwd(), output_dir))



##############################开始填充#########################
# ============================================================
# 手动创建截断阈值列表 (Ground Truth Cutoffs)
# ============================================================

# 1. 创建 10% 截断的阈值表
cutoff_10 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.1,
  Lower_Cutoff = c(20.061273, 200.055249, 8.495351, 19.000000, 155.000000, 0.630000, 20.662977, 2.811597),
  Upper_Cutoff = c(36.500505, 3093.325368, 43.406748, 34.000000, 172.700000, 1.052760, 76.117738, 7.467111),
  stringsAsFactors = FALSE
)

# 2. 创建 20% 截断的阈值表
cutoff_20 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.2,
  Lower_Cutoff = c(21.5316497, 341.2159561, 11.8236356, 21.0000000, 158.0000000, 0.7862489, 26.5367602, 3.3568065),
  Upper_Cutoff = c(31.4288844, 2044.8760798, 35.2850673, 32.0000000, 170.0000000, 0.9451719, 60.6895559, 6.4242105),
  stringsAsFactors = FALSE
)

# 3. 创建 30% 截断的阈值表
cutoff_30 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.3,
  Lower_Cutoff = c(22.6611810, 508.1810942, 14.1212004, 23.0000000, 160.2000000, 0.8218406, 30.6946994, 3.8344199),
  Upper_Cutoff = c(28.5752868, 1518.1500077, 28.0546296, 30.0000000, 167.5000000, 0.9198789, 53.0351045, 5.7293256),
  stringsAsFactors = FALSE
)

# 4. 创建 40% 截断的阈值表
cutoff_40 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.4,
  Lower_Cutoff = c(23.8219026, 711.0014683, 16.8910435, 25.0000000, 162.5000000, 0.8434391, 35.6720577, 4.2460000),
  Upper_Cutoff = c(26.7537356, 1229.4455148, 23.5820478, 29.0000000, 165.1000000, 0.8948467, 45.8456378, 5.1540000),
  stringsAsFactors = FALSE
)

# 5. 将这四个表打包进一个 List
cutoff_list <- list(
  df_trunc_10 = cutoff_10,
  df_trunc_20 = cutoff_20,
  df_trunc_30 = cutoff_30,
  df_trunc_40 = cutoff_40
)

# ============================================================
# 调用示例
# ============================================================

# 查看 10% 截断时的 BMI 阈值
# bmi_cutoffs_10 <- cutoff_list$df_trunc_10[cutoff_list$df_trunc_10$Variable == "BMI", ]
# print(bmi_cutoffs_10)

# 打印整个列表结构
print(str(cutoff_list))

message("cutoff_list 已创建完成。")


##########################random imputation###########################
library(dplyr)
library(readr)

# ============================================================
# 0. 准备工作
# ============================================================

# 1. 定义物理/理论极限 (作为唯一的填充范围依据)
limits <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Lower_limit = c(15.35315791, 50.38720449, 1.900280425, 13.4904681089308, 140, 0.26, 3.12588676962321, 1.20355732717561),
  Upper_limit = c(58.8941687473194, 14175.77793, 106.20096951712, 43, 196.5, 4.47, 301.970355243691, 15.9587598514815),
  stringsAsFactors = FALSE
)

# 2. 定义文件夹路径
dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

# 3. 创建输出文件夹 (命名为 Random_Baseline 以示区分)
output_base_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/imputed_Random_Baseline"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir)

# 设置随机种子
set.seed(2025)

# ============================================================
# 1. 循环处理四个数据集
# ============================================================

for (i in seq_along(dirs)) {
  
  level_name <- names(dirs)[i]     
  dir_path   <- dirs[[i]]          
  
  # 读取 CSV
  csv_file <- list.files(dir_path, pattern = "*.csv", full.names = TRUE)[1]
  
  if (is.na(csv_file)) {
    warning(paste("在文件夹", dir_path, "中未找到 CSV 文件，跳过。"))
    next
  }
  
  message(paste("正在处理 (Baseline Random):", level_name, "| 文件:", basename(csv_file)))
  
  df <- read_csv(csv_file, show_col_types = FALSE)
  
  # ============================================================
  # 2. 执行完全随机填充 (Full Range Uniform)
  # ============================================================
  
  # 遍历 limits 表中的每一个变量
  for (j in 1:nrow(limits)) {
    
    var_name <- limits$Variable[j]
    phys_min <- limits$Lower_limit[j]
    phys_max <- limits$Upper_limit[j]
    
    ind_col <- paste0(var_name, "_ind")
    
    # 检查变量是否存在
    if (var_name %in% names(df) && ind_col %in% names(df)) {
      
      # 找出所有因为截断而缺失的行 (不管它是 Below 还是 Above)
      # 逻辑：只要是 NA，且指示列标记了原因，就视作需要填充的目标
      mask_missing <- is.na(df[[var_name]]) & 
        (df[[ind_col]] %in% c("Below_Lower", "Above_Upper"))
      
      n_missing <- sum(mask_missing, na.rm = TRUE)
      
      if (n_missing > 0) {
        # --- 核心修改 ---
        # 不考虑截断方向，直接在 [Phys_Min, Phys_Max] 全域内随机生成
        random_vals <- runif(n_missing, min = phys_min, max = phys_max)
        
        df[[var_name]][mask_missing] <- random_vals
      }
    }
  }
  
  # ============================================================
  # 3. 保存文件
  # ============================================================
  
  save_dir <- file.path(output_base_dir, paste0("imputed_", level_name))
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  save_name <- file.path(save_dir, paste0("imputed_RandomBase_", level_name, ".csv"))
  write_csv(df, save_name)
  
  message(paste("  -> 已保存至:", save_name))
}

message("所有 Baseline Random 填充任务完成！")

#########################fixed value imputation 1#####################

library(dplyr)
library(readr)

# ============================================================
# 0. 准备工作：定义路径和 Cutoff 信息
# ============================================================

# 定义四个数据集的文件夹路径
dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

# 对应 cutoff_list 中的名字 (确保你之前运行过 cutoff_list 的定义代码)
cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

# 创建输出文件夹 (保存填充后的数据)
output_base_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/imputed_LOD"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir)

# ============================================================
# 1. 循环处理四个数据集
# ============================================================

# 使用索引来同时遍历文件夹和 Cutoff 列表
for (i in seq_along(dirs)) {
  
  level_name <- names(dirs)[i]     # e.g., "10%"
  dir_path   <- dirs[[i]]          # 文件夹路径
  cutoff_key <- cutoff_keys[i]     # e.g., "df_trunc_10"
  
  # 获取当前层级的 Cutoff 标准答案
  current_cutoffs <- cutoff_list[[cutoff_key]]
  
  # 读取数据集 (假设文件夹里有一个 dataset_trunc_XX.csv)
  # 这里我们需要找到该文件夹下的 csv 文件
  csv_file <- list.files(dir_path, pattern = "*.csv", full.names = TRUE)[1]
  
  if (is.na(csv_file)) {
    warning(paste("在文件夹", dir_path, "中未找到 CSV 文件，跳过。"))
    next
  }
  
  message(paste("正在处理:", level_name, "| 文件:", basename(csv_file)))
  
  # 读取数据
  df <- read_csv(csv_file, show_col_types = FALSE)
  
  # ============================================================
  # 2. 执行 LOD 填充 (Cutoff Imputation)
  # ============================================================
  
  # 遍历 Cutoff 表中的每一个变量进行填充
  for (j in 1:nrow(current_cutoffs)) {
    
    var_name <- current_cutoffs$Variable[j]
    lower_val <- current_cutoffs$Lower_Cutoff[j]
    upper_val <- current_cutoffs$Upper_Cutoff[j]
    
    # 构造指示列的名字 (例如 BMI_ind)
    ind_col <- paste0(var_name, "_ind")
    
    # 检查变量和指示列是否存在于数据中
    if (var_name %in% names(df) && ind_col %in% names(df)) {
      
      # --- 核心填充逻辑 ---
      
      # 1. 填充下界 (Below Lower -> 填入 Lower Cutoff)
      # 只有当它是 NA 且 指示列显示 "Below_Lower" 时才填
      mask_lower <- is.na(df[[var_name]]) & (df[[ind_col]] == "Below_Lower")
      df[[var_name]][mask_lower] <- lower_val
      
      # 2. 填充上界 (Above Upper -> 填入 Upper Cutoff)
      # 只有当它是 NA 且 指示列显示 "Above_Upper" 时才填
      mask_upper <- is.na(df[[var_name]]) & (df[[ind_col]] == "Above_Upper")
      df[[var_name]][mask_upper] <- upper_val
      
    }
  }
  
  # ============================================================
  # 3. 保存填充后的数据
  # ============================================================
  
  # 为这个比例创建一个子文件夹
  save_dir <- file.path(output_base_dir, paste0("imputed_", level_name))
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  # 保存文件
  save_name <- file.path(save_dir, paste0("imputed_LOD_", level_name, ".csv"))
  write_csv(df, save_name)
  
  message(paste("  -> 已保存至:", save_name))
}

message("所有 LOD 填充任务完成！")

###############################fixed value imputation 2#############################################

library(dplyr)
library(readr)
library(survival) # 必须加载，用于 survreg 函数

# ============================================================
# 0. 准备工作
# ============================================================

# 定义输入文件夹 (截断数据集路径)
input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

# 对应 cutoff_list 的 key
cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

# 定义输出文件夹
output_base_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/imputed_LOD_Half_Survival"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir)

# ============================================================
# 1. 循环处理四个截断比例
# ============================================================

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  dir_path   <- input_dirs[[i]]
  cutoff_key <- cutoff_keys[i]
  
  # 获取当前的 Cutoff 标准
  current_cutoffs <- cutoff_list[[cutoff_key]]
  
  # 读取 CSV 文件
  csv_file <- list.files(dir_path, pattern = "*.csv", full.names = TRUE)[1]
  if (is.na(csv_file)) next
  
  message(paste("正在处理:", level_name, "| 文件:", basename(csv_file)))
  df <- read_csv(csv_file, show_col_types = FALSE)
  
  # ============================================================
  # 2. 逐变量进行混合填充
  # ============================================================
  
  for (j in 1:nrow(current_cutoffs)) {
    
    var_name  <- current_cutoffs$Variable[j]
    lower_val <- current_cutoffs$Lower_Cutoff[j]
    upper_val <- current_cutoffs$Upper_Cutoff[j]
    ind_col   <- paste0(var_name, "_ind")
    
    # 确保变量存在
    if (var_name %in% names(df) && ind_col %in% names(df)) {
      
      # --- 步骤 A: 左侧截断 (Left-Censored) -> LOD/2 替代 ---
      mask_lower <- is.na(df[[var_name]]) & (df[[ind_col]] == "Below_Lower")
      
      # 替换为 LOD_lower / 2
      df[[var_name]][mask_lower] <- lower_val / 2
      
      # --- 步骤 B: 右侧截断 (Right-Censored) -> Survival Conditional Mean ---
      mask_upper <- is.na(df[[var_name]]) & (df[[ind_col]] == "Above_Upper")
      
      # 如果存在右侧缺失，才进行模型拟合
      if (sum(mask_upper) > 0) {
        
        # 1. 准备生存分析数据
        # 构造临时向量 y_temp:
        # - 已观测值 (Observed): 保持原值
        # - 左侧填充值 (Left Imputed): 视为已观测 (Status=1)，值为 LOD/2
        # - 右侧截断值 (Right Censored): 暂时设为 Upper Cutoff，且 Status=0
        
        y_vec <- df[[var_name]]
        status_vec <- rep(1, length(y_vec)) # 初始化状态为 1 (Event/Observed)
        
        # 对于右侧截断的点：
        # Time = LOD_upper
        # Status = 0 (Censored)
        y_vec[mask_upper] <- upper_val
        status_vec[mask_upper] <- 0 
        
        # 2. 拟合 Log-normal AFT 模型
        # survreg(Surv(time, status) ~ 1, dist = "lognormal")
        # 只用截距模型 (~ 1)，因为我们是拟合边缘分布 (Marginal Distribution)
        tryCatch({
          fit <- survreg(Surv(y_vec, status_vec) ~ 1, dist = "lognormal")
          
          # 3. 提取参数
          # survreg 的 Intercept 对应 log-normal 的 mu
          # survreg 的 Scale 对应 log-normal 的 sigma
          mu_hat <- coef(fit)["(Intercept)"]
          sigma_hat <- fit$scale
          
          # 4. 计算条件均值 E[X | X > LOD_upper]
          # 公式：
          # num = 1 - Phi( (ln(LOD) - mu - sigma^2) / sigma )
          # den = 1 - Phi( (ln(LOD) - mu) / sigma )
          # factor = exp(mu + 0.5 * sigma^2)
          # Result = factor * (num / den)
          
          log_LOD <- log(upper_val)
          
          term1 <- (log_LOD - mu_hat - sigma_hat^2) / sigma_hat
          term2 <- (log_LOD - mu_hat) / sigma_hat
          
          numerator   <- 1 - pnorm(term1)
          denominator <- 1 - pnorm(term2)
          pre_factor  <- exp(mu_hat + 0.5 * sigma_hat^2)
          
          cond_mean <- pre_factor * (numerator / denominator)
          
          # 5. 填充右侧缺失值
          df[[var_name]][mask_upper] <- cond_mean
          
        }, error = function(e) {
          message(paste("  [Warning] 模型拟合失败:", var_name, "-", e$message))
          # 如果拟合失败（极少情况），可以用备选方案（如直接填 Cutoff）
          df[[var_name]][mask_upper] <- upper_val 
        })
      }
    }
  }
  
  # ============================================================
  # 3. 保存结果
  # ============================================================
  
  save_dir <- file.path(output_base_dir, paste0("imputed_", level_name))
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  save_name <- file.path(save_dir, paste0("imputed_LOD_Half_Survival_", level_name, ".csv"))
  write_csv(df, save_name)
  
  message(paste("  -> 已保存至:", save_name))
}

message("所有混合填充 (LOD/2 + Survival Mean) 任务完成！")

##################################distribution based imputation##################
library(dplyr)
library(readr)
library(survival) # 用于 survreg

library(dplyr)
library(readr)
library(AER) # 核心包：包含 tobit 函数

# ============================================================
# 0. 准备工作
# ============================================================

input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

# 定义输出路径
out_gauss_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_Gaussian_Draw"
out_lnorm_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_Lognormal_Draw"

if (!dir.exists(out_gauss_dir)) dir.create(out_gauss_dir)
if (!dir.exists(out_lnorm_dir)) dir.create(out_lnorm_dir)

# ============================================================
# 基于 Tobit 的随机填充函数
# ============================================================
impute_with_tobit <- function(df, cutoffs, target_dist) {
  
  df_out <- df 
  
  for (j in 1:nrow(cutoffs)) {
    var_name  <- cutoffs$Variable[j]
    lower_val <- cutoffs$Lower_Cutoff[j]
    upper_val <- cutoffs$Upper_Cutoff[j]
    ind_col   <- paste0(var_name, "_ind")
    
    # 确保变量存在
    if (var_name %in% names(df) && ind_col %in% names(df)) {
      
      vals <- df[[var_name]]
      inds <- df[[ind_col]]
      
      # 1. 确定掩码 (Masks)
      mask_obs   <- !is.na(vals)
      mask_lower <- is.na(vals) & (inds == "Below_Lower")
      mask_upper <- is.na(vals) & (inds == "Above_Upper")
      
      mask_lower[is.na(mask_lower)] <- FALSE
      mask_upper[is.na(mask_upper)] <- FALSE
      
      if (sum(mask_lower) == 0 && sum(mask_upper) == 0) next
      
      # 2. 准备 Tobit 需要的数据 (Temp Value)
      # Tobit 要求截断点必须有值（等于边界值），不能是 NA
      temp_vec <- vals
      temp_vec[mask_lower] <- lower_val
      temp_vec[mask_upper] <- upper_val
      
      # 3. 处理分布 (Gaussian vs Lognormal)
      # 如果是 Lognormal，我们在 log 尺度上做 Tobit (即 Log-Tobit)
      use_log <- (target_dist == "lognormal")
      
      work_vec   <- temp_vec
      work_left  <- lower_val
      work_right <- upper_val
      
      if (use_log) {
        # 数据清洗：确保 > 0，防止 log 报错
        if (any(work_vec <= 0)) work_vec <- pmax(work_vec, 1e-6)
        work_left  <- max(lower_val, 1e-6)
        work_right <- max(upper_val, 1e-6)
        
        # 取对数
        work_vec   <- log(work_vec)
        work_left  <- log(work_left)
        work_right <- log(work_right)
      } else {
        # Gaussian 下界可以用 -Inf，但在 Tobit 中如果数据等于 left 参数就被视为截断
        # 所以这里不需要改动，只要 work_vec[mask_lower] == work_left 即可
      }
      
      # 4. 拟合 Tobit 模型
      tryCatch({
        # 这里的关键：temp_value ~ 1 表示只估计均值和方差
        # left 和 right 参数告诉函数哪些值是截断的
        tobit_model <- tobit(work_vec ~ 1, left = work_left, right = work_right)
        
        # 提取参数 (Tobit 默认假设高斯分布)
        mu_hat    <- coef(tobit_model)["(Intercept)"]
        sigma_hat <- tobit_model$scale
        
        # 5. DB-Draw: 随机抽样填充
        # 使用 Inverse-CDF 法：runif -> qnorm
        
        set.seed(123 + j + as.numeric(charToRaw(target_dist))[1])
        
        # >>> Left Side <<<
        if (sum(mask_lower) > 0) {
          n_miss <- sum(mask_lower)
          # F(L) 在当前尺度 (Log or Original)
          p_limit <- pnorm(work_left, mean = mu_hat, sd = sigma_hat)
          p_limit <- max(min(p_limit, 1 - 1e-9), 1e-9)
          
          u_draws <- runif(n_miss, min = 0, max = p_limit)
          imputed <- qnorm(u_draws, mean = mu_hat, sd = sigma_hat)
          
          # 如果是 Log 尺度，还原
          if (use_log) imputed <- exp(imputed)
          
          df_out[[var_name]][mask_lower] <- imputed
        }
        
        # >>> Right Side <<<
        if (sum(mask_upper) > 0) {
          n_miss <- sum(mask_upper)
          # F(U) 在当前尺度
          p_limit <- pnorm(work_right, mean = mu_hat, sd = sigma_hat)
          p_limit <- max(min(p_limit, 1 - 1e-9), 1e-9)
          
          u_draws <- runif(n_miss, min = p_limit, max = 1)
          imputed <- qnorm(u_draws, mean = mu_hat, sd = sigma_hat)
          
          if (use_log) imputed <- exp(imputed)
          
          df_out[[var_name]][mask_upper] <- imputed
        }
        
      }, error = function(e) {
        message(paste("  [Error]", target_dist, "Tobit 拟合失败:", var_name, "-", e$message))
        # 失败兜底：填 Cutoff
        df_out[[var_name]][mask_lower] <- lower_val
        df_out[[var_name]][mask_upper] <- upper_val
      })
    }
  }
  return(df_out)
}

# ============================================================
# 主循环
# ============================================================

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  dir_path   <- input_dirs[[i]]
  cutoff_key <- cutoff_keys[i]
  
  current_cutoffs <- cutoff_list[[cutoff_key]]
  
  csv_file <- list.files(dir_path, pattern = "*.csv", full.names = TRUE)[1]
  if (is.na(csv_file)) next
  
  message(paste("正在处理:", level_name, "| 文件:", basename(csv_file)))
  df_orig <- read_csv(csv_file, show_col_types = FALSE)
  
  # 1. Gaussian 填充 (使用 Tobit)
  message("  -> 开始 Gaussian (Tobit) 拟合...")
  df_gauss <- impute_with_tobit(df_orig, current_cutoffs, target_dist = "gaussian")
  
  save_dir_g <- file.path(out_gauss_dir, paste0("imputed_", level_name))
  if (!dir.exists(save_dir_g)) dir.create(save_dir_g)
  write_csv(df_gauss, file.path(save_dir_g, paste0("imputed_Gaussian_Draw_", level_name, ".csv")))
  
  # 2. Lognormal 填充 (使用 Log-Tobit)
  message("  -> 开始 Lognormal (Log-Tobit) 拟合...")
  df_lnorm <- impute_with_tobit(df_orig, current_cutoffs, target_dist = "lognormal")
  
  save_dir_l <- file.path(out_lnorm_dir, paste0("imputed_", level_name))
  if (!dir.exists(save_dir_l)) dir.create(save_dir_l)
  write_csv(df_lnorm, file.path(save_dir_l, paste0("imputed_Lognormal_Draw_", level_name, ".csv")))
  
  message("  -> 保存完成。")
}

message("所有 Tobit 任务圆满完成！")

#################################MICE imputation##############################
library(mice)
library(truncnorm) # 用于截断正态分布采样
library(dplyr)

# 确保之前的 cutoff_list 和 truncated 数据集路径已经准备好
# 如果没有加载，请先运行你之前的 cutoff_list 定义代码

# ============================================================
# 辅助函数：获取每个缺失值的动态边界 [L_i, U_i]
# ============================================================
# 这个函数会根据指示列 (e.g., BMI_ind) 判断每一行缺失值的允许范围
get_dynamic_bounds <- function(y_name, n_rows, data_df, cutoff_info) {
  
  # 获取当前变量的阈值标准
  var_cutoffs <- cutoff_info %>% filter(Variable == y_name)
  
  # 如果找不到阈值信息（比如非截断变量），返回无限制 (-Inf, Inf)
  if (nrow(var_cutoffs) == 0) {
    return(list(lower = rep(-Inf, n_rows), upper = rep(Inf, n_rows)))
  }
  
  l_limit <- var_cutoffs$Lower_Cutoff
  u_limit <- var_cutoffs$Upper_Cutoff
  
  # 获取指示列 (例如 BMI -> BMI_ind)
  ind_col <- paste0(y_name, "_ind")
  
  # 初始化边界为无限制
  lower_vec <- rep(-Inf, n_rows)
  upper_vec <- rep(Inf, n_rows)
  
  # 如果数据中有指示列，根据指示列设置边界
  if (ind_col %in% names(data_df)) {
    indicators <- data_df[[ind_col]]
    
    # 情况 1: Below Lower -> 允许范围 (-Inf, Lower_Limit]
    # 注意：我们是对"潜在值"建模，所以必须小于下限
    idx_below <- which(indicators == "Below_Lower")
    upper_vec[idx_below] <- l_limit
    
    # 情况 2: Above Upper -> 允许范围 [Upper_Limit, Inf)
    idx_above <- which(indicators == "Above_Upper")
    lower_vec[idx_above] <- u_limit
  }
  
  return(list(lower = lower_vec, upper = upper_vec))
}
# ============================================================
# 方法 1: Post-hoc Truncation (事后截断)
# ============================================================
mice.impute.posthoc <- function(y, ry, x, wy = NULL, ...) {
  # 1. 使用标准正态回归 (Norm) 进行预测
  # mice.impute.norm 会返回贝叶斯线性回归的预测值
  imputed_values <- mice.impute.norm(y, ry, x, wy = wy, ...)
  
  # 2. 获取当前上下文信息
  # mice 也就是在这个函数里不知道变量名，我们需要从 parent.frame 获取数据
  # 这是一个高级技巧，用于在 mice 内部获取完整数据框
  # 我们约定：调用 mice 时，必须传递 custom_cutoff_info 和 data_ref
  
  # 从 ... 中提取参数，或者从全局环境获取 (为了简单，这里假设全局有 global_df_ref)
  # 在实际调用循环中，我们会设置这个全局变量
  
  # 识别当前是哪个变量 (通过对比 y 和 global_df_ref 中的列)
  # 这种匹配稍微有点风险，但对于固定流程是可行的
  # 更稳健的方法是利用 ... 传递 varname，但 mice 接口不支持
  
  # 简单方案：我们假设外部循环会设置一个全局变量 CURRENT_VAR_NAME
  # (虽然不优雅，但在 mice 自定义函数中非常实用)
  
  bounds <- get_dynamic_bounds(CURRENT_VAR_NAME, length(y), global_df_ref, global_cutoff_info)
  
  # 提取缺失值对应的边界 (wy 是缺失值的逻辑向量)
  # mice.impute.norm 返回的是 wy=TRUE 部分的插补值
  # 所以我们需要 wy=TRUE 部分的边界
  l_bounds_miss <- bounds$lower[wy]
  u_bounds_miss <- bounds$upper[wy]
  
  # 3. 实施截断逻辑 (Post-hoc Logic)
  # 如果 pred < L (且 L > -Inf)，设为 L? 
  # 不，根据你的文本: 
  #   - 如果 pred < L (对于左删失点，L是可行上限)，无需操作
  #   - 如果 pred > L (对于左删失点，这是违规的)，设为 L
  
  # 逻辑修正：
  # 左删失点 (Below Lower): 范围 (-Inf, LOD_L]. 如果 pred > LOD_L，强设为 LOD_L.
  # 右删失点 (Above Upper): 范围 [LOD_U, Inf). 如果 pred < LOD_U，强设为 LOD_U.
  
  # 对比 Lower Bound (即 u_bounds_miss 对于左删失点是 LOD_L? 不，get_dynamic_bounds 返回的是区间 [L, U])
  # 左删失: lower=-Inf, upper=LOD_L
  # 右删失: lower=LOD_U, upper=Inf
  
  # 修正小于下界的 (针对右删失点预测偏低的情况)
  imputed_values <- pmax(imputed_values, l_bounds_miss)
  
  # 修正大于上界的 (针对左删失点预测偏高的情况)
  imputed_values <- pmin(imputed_values, u_bounds_miss)
  
  return(imputed_values)
}

# ============================================================
# 方法 2: Truncated Normal Sampling (截断正态采样)
# ============================================================
mice.impute.truncnorm <- function(y, ry, x, wy = NULL, ...) {
  # 1. 手动执行贝叶斯线性回归步骤 (仿照 mice.impute.norm)
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw(y, ry, x) # mice 内部函数，抽取 beta 和 sigma
  
  # 2. 计算预测均值 (mu)
  x_miss <- x[wy, , drop = FALSE]
  mu_miss <- x_miss %*% parm$beta
  sigma <- parm$sigma
  
  # 3. 获取边界
  bounds <- get_dynamic_bounds(CURRENT_VAR_NAME, length(y), global_df_ref, global_cutoff_info)
  l_bounds_miss <- bounds$lower[wy]
  u_bounds_miss <- bounds$upper[wy]
  
  # 4. 从截断正态分布采样
  # rtruncnorm(n, a, b, mean, sd)
  imputed_values <- rtruncnorm(sum(wy), 
                               a = l_bounds_miss, 
                               b = u_bounds_miss, 
                               mean = mu_miss, 
                               sd = sigma)
  
  return(imputed_values)
}

# 必须导出 .norm.draw 供我们使用 (它通常是隐藏的)
.norm.draw <- mice:::.norm.draw

# ============================================================
# 主执行循环
# ============================================================

# 输入输出设置
input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)
cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

# 输出目录
out_mice_posthoc <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_Posthoc"
out_mice_trunc   <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_TruncNormal"
if(!dir.exists(out_mice_posthoc)) dir.create(out_mice_posthoc)
if(!dir.exists(out_mice_trunc)) dir.create(out_mice_trunc)

# MICE 参数
M_IMP <- 5   # 多重插补次数
MAXIT <- 5   # 迭代次数

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  cutoff_key <- cutoff_keys[i]
  
  # 1. 读取数据
  csv_file <- list.files(input_dirs[[i]], pattern = "*.csv", full.names = TRUE)[1]
  if(is.na(csv_file)) next
  
  message(paste("处理:", level_name, "| Reading:", basename(csv_file)))
  df_orig <- read.csv(csv_file) # 使用 read.csv 保持 data.frame 格式
  
  # 设置全局变量供自定义函数使用
  global_cutoff_info <<- cutoff_list[[cutoff_key]]
  global_df_ref <<- df_orig
  
  # 2. 配置 MICE
  # 初始化 mice
  ini <- mice(df_orig, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  # 确定需要插补的变量 (在 Cutoff 列表里的)
  target_vars <- global_cutoff_info$Variable
  
  # 移除指示列 (_ind) 作为预测变量，避免多重共线性或泄露
  # (指示列完全决定了截断状态，不能直接放进回归方程)
  ind_cols <- grep("_ind$", names(df_orig), value = TRUE)
  pred[, ind_cols] <- 0 
  # 确保指示列不被插补
  meth[ind_cols] <- ""
  
  # --- 执行方法 1: MICE Post-hoc ---
  meth_posthoc <- meth
  # 动态修改方法
  # mice 允许我们直接把函数名赋值给 method
  for(v in target_vars) {
    if(v %in% names(df_orig)) {
      meth_posthoc[v] <- "posthoc" # 对应 mice.impute.posthoc
    }
  }
  
  message("  -> Running MICE Post-hoc...")
  # 我们需要 hook 每次插补，因为 mice 每次只处理一个变量
  # 我们可以利用 visitSequence 或者是 mice 的内部机制
  # 实际上，只要 meth[v] == "posthoc"，mice 就会去找 mice.impute.posthoc
  # 关键点：我们需要让函数知道当前处理的是哪个变量。
  # 这里的 trick 是：我们无法直接传参。
  # 所以我们必须写死函数吗？或者，我们可以为每个变量生成一个 wrapper？
  # 最简单的 Hack：修改 mice.impute.posthoc，让它去匹配 y 和 global_df_ref 的列
  # (上面的 get_dynamic_bounds 已经预留了 CURRENT_VAR_NAME，我们需要在这里设置)
  
  # 修正：mice 不支持在迭代中传递变量名。
  # 我们可以使用 `mice` 的 `post` 处理参数？
  # 不，最稳妥的方法是为每个变量定义一个特定的函数，或者更高级的：
  # 在 mice.impute.posthoc 内部，通过 `match.call()` 或者比较 `y` 值来反推变量名。
  # 下面是一个自动反推变量名的增强版函数逻辑 (放在循环内定义，利用闭包特性)
  
  make_imputer <- function(method_type) {
    function(y, ry, x, wy = NULL, ...) {
      # 反推变量名：找到 global_df_ref 中与 y (非缺失部分) 一致的那一列
      var_name <- NA
      for(col in names(global_df_ref)) {
        # 比较非缺失值是否完全相等
        if(is.numeric(global_df_ref[[col]]) && 
           length(y) == length(global_df_ref[[col]]) &&
           all(y[ry] == global_df_ref[[col]][ry])) {
          var_name <- col
          break
        }
      }
      
      # 设置全局变量给 helper 使用
      assign("CURRENT_VAR_NAME", var_name, envir = .GlobalEnv)
      
      if(method_type == "posthoc") {
        return(mice.impute.posthoc(y, ry, x, wy, ...))
      } else {
        return(mice.impute.truncnorm(y, ry, x, wy, ...))
      }
    }
  }
  
  # 注册临时函数到 mice 命名空间 (或者全局)
  mice.impute.posthoc_wrapper <<- make_imputer("posthoc")
  meth_posthoc[target_vars] <- "posthoc_wrapper"
  
  imp_posthoc <- mice(df_orig, method = meth_posthoc, predictorMatrix = pred, 
                      m = M_IMP, maxit = MAXIT, printFlag = FALSE)
  
  # 保存 Posthoc 结果 (保存第1个数据集作为示例，或者保存 pooled mean，通常保存 complete stacked)
  # 这里我们保存所有 m 个数据集的 stack
  df_posthoc_out <- complete(imp_posthoc, action = "long", include = TRUE)
  write.csv(df_posthoc_out, 
            file.path(out_mice_posthoc, paste0("mice_posthoc_long_", level_name, ".csv")), 
            row.names = FALSE)
  
  # --- 执行方法 2: MICE Truncated Normal ---
  meth_trunc <- meth
  mice.impute.trunc_wrapper <<- make_imputer("trunc")
  meth_trunc[target_vars] <- "trunc_wrapper"
  
  message("  -> Running MICE Truncated Normal...")
  imp_trunc <- mice(df_orig, method = meth_trunc, predictorMatrix = pred, 
                    m = M_IMP, maxit = MAXIT, printFlag = FALSE)
  
  df_trunc_out <- complete(imp_trunc, action = "long", include = TRUE)
  write.csv(df_trunc_out, 
            file.path(out_mice_trunc, paste0("mice_truncnorm_long_", level_name, ".csv")), 
            row.names = FALSE)
  
  message("  -> Done.")
}

library(mice)
library(truncnorm)
library(dplyr)

# ============================================================
# 辅助函数：将 MICE 的 m 个数据集“平均”为一个
# ============================================================
pool_mice_to_single <- function(mice_obj) {
  # 1. 提取所有 m 次插补的数据 (Long Format)
  # include = FALSE 表示不包含含有 NA 的原始数据，只包含 m 个插补好的
  long_dat <- complete(mice_obj, action = "long", include = FALSE)
  
  # 2. 按原始行号 (.id) 分组
  # 对于数值型变量：计算 m 次插补的均值 (Mean)
  # 对于非数值型变量 (如 _ind 指示列)：直接取第一个值 (因为它们在插补中不会变)
  pooled_dat <- long_dat %>%
    group_by(.id) %>%
    summarise(across(everything(), function(x) {
      if(is.numeric(x)) mean(x) else dplyr::first(x)
    })) %>%
    select(-.id, -.imp) # 移除 MICE 自动生成的辅助列
  
  return(as.data.frame(pooled_dat))
}

# ... (这里省略之前的 make_imputer 和辅助函数定义，请保持之前的一样) ...
# 请确保此时环境中已经定义了: input_dirs, cutoff_keys, global_cutoff_info 等

# ============================================================
# 主循环 (包含 averaging 步骤)
# ============================================================

# 定义输出路径
out_mice_posthoc <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_Posthoc"
out_mice_trunc   <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_TruncNormal"

if(!dir.exists(out_mice_posthoc)) dir.create(out_mice_posthoc)
if(!dir.exists(out_mice_trunc)) dir.create(out_mice_trunc)

M_IMP <- 5
MAXIT <- 5

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  cutoff_key <- cutoff_keys[i]
  
  # 读取数据
  csv_file <- list.files(input_dirs[[i]], pattern = "*.csv", full.names = TRUE)[1]
  if(is.na(csv_file)) next
  
  message(paste("处理:", level_name, "| Reading:", basename(csv_file)))
  df_orig <- read.csv(csv_file)
  
  # 设置全局变量
  global_cutoff_info <<- cutoff_list[[cutoff_key]]
  global_df_ref <<- df_orig
  
  # 配置 MICE 基础
  ini <- mice(df_orig, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  target_vars <- global_cutoff_info$Variable
  ind_cols <- grep("_ind$", names(df_orig), value = TRUE)
  pred[, ind_cols] <- 0
  meth[ind_cols] <- ""
  
  # --------------------------------------------------------
  # 1. Post-hoc Truncation -> Averaged Single Dataset
  # --------------------------------------------------------
  meth_posthoc <- meth
  # 重新注册 wrapper (防止环境丢失)
  mice.impute.posthoc_wrapper <<- make_imputer("posthoc")
  meth_posthoc[target_vars] <- "posthoc_wrapper"
  
  message("  -> Running MICE Post-hoc (m=5)...")
  imp_posthoc <- mice(df_orig, method = meth_posthoc, predictorMatrix = pred, 
                      m = M_IMP, maxit = MAXIT, printFlag = FALSE)
  
  # >>> 关键步骤：取平均生成单一数据集 <<<
  df_posthoc_single <- pool_mice_to_single(imp_posthoc)
  
  # 保存单一数据集
  save_path <- file.path(out_mice_posthoc, paste0("mice_posthoc_averaged_", level_name, ".csv"))
  write.csv(df_posthoc_single, save_path, row.names = FALSE)
  message(paste("     [Saved] Single Averaged Dataset:", basename(save_path)))
  
  # --------------------------------------------------------
  # 2. Truncated Normal -> Averaged Single Dataset
  # --------------------------------------------------------
  meth_trunc <- meth
  mice.impute.trunc_wrapper <<- make_imputer("trunc")
  meth_trunc[target_vars] <- "trunc_wrapper"
  
  message("  -> Running MICE Truncated Normal (m=5)...")
  imp_trunc <- mice(df_orig, method = meth_trunc, predictorMatrix = pred, 
                    m = M_IMP, maxit = MAXIT, printFlag = FALSE)
  
  # >>> 关键步骤：取平均生成单一数据集 <<<
  df_trunc_single <- pool_mice_to_single(imp_trunc)
  
  # 保存单一数据集
  save_path <- file.path(out_mice_trunc, paste0("mice_truncnorm_averaged_", level_name, ".csv"))
  write.csv(df_trunc_single, save_path, row.names = FALSE)
  message(paste("     [Saved] Single Averaged Dataset:", basename(save_path)))
}

message("所有 MICE (Averaged) 填充任务完成！")

#####################################TOBIT imputation############################
library(dplyr)
library(readr)
library(AER)       # 核心包：包含 tobit 函数
library(truncnorm) # 用于从截断正态分布中采样

# ============================================================
# 0. 准备工作
# ============================================================

input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

# 定义输出路径
output_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_Tobit"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ============================================================
# 核心函数：单次 Tobit 插补 (包含迭代过程)
# ============================================================
impute_tobit_single_run <- function(df, cutoff_info, max_iter = 5, seed_val = 123) {
  
  set.seed(seed_val)
  
  # 1. 初始化：先用均值填充所有 NA，作为起始点 (predictor 不能有 NA)
  df_filled <- df
  target_vars <- cutoff_info$Variable
  
  # 简单的初始化循环
  for (v in names(df_filled)) {
    if (is.numeric(df_filled[[v]]) && any(is.na(df_filled[[v]]))) {
      df_filled[[v]][is.na(df_filled[[v]])] <- mean(df_filled[[v]], na.rm = TRUE)
    }
  }
  
  # 2. 迭代优化 (Gibbs Sampling style)
  for (iter in 1:max_iter) {
    
    # 对每个需要插补的变量进行循环
    for (j in 1:nrow(cutoff_info)) {
      
      var_name  <- cutoff_info$Variable[j]
      lower_val <- cutoff_info$Lower_Cutoff[j]
      upper_val <- cutoff_info$Upper_Cutoff[j]
      ind_col   <- paste0(var_name, "_ind")
      
      if (!var_name %in% names(df)) next
      
      # 识别缺失位置
      # 注意：我们使用原始 df 的 _ind 来确定哪些是真的缺失
      orig_inds <- df[[ind_col]]
      mask_lower <- (orig_inds == "Below_Lower")
      mask_upper <- (orig_inds == "Above_Upper")
      
      mask_lower[is.na(mask_lower)] <- FALSE
      mask_upper[is.na(mask_upper)] <- FALSE
      
      if (sum(mask_lower) == 0 && sum(mask_upper) == 0) next
      
      # 3. 准备 Tobit 回归数据
      # Tobit 要求 Y 变量在截断点必须有值，不能是 NA
      # 我们使用当前 df_filled 中的值，但在截断点处，必须强制设为边界值
      # 以便 tobit 函数能识别出它们是 censored observation
      
      y_temp <- df_filled[[var_name]]
      
      # 关键步骤：AER::tobit 需要将 censored values 设为边界值
      y_temp[mask_lower] <- lower_val
      y_temp[mask_upper] <- upper_val
      
      # 构建回归公式：Y ~ 所有其他数值型变量
      # 排除指示列和其他非数值列
      predictors <- setdiff(names(df_filled), c(var_name, grep("_ind$", names(df_filled), value = TRUE)))
      # 只保留数值型预测变量
      predictors <- predictors[sapply(df_filled[predictors], is.numeric)]
      
      formula_str <- paste(var_name, "~", paste(predictors, collapse = " + "))
      
      # 构造临时数据框用于拟合 (Y 列是处理过的 y_temp)
      data_for_fit <- df_filled
      data_for_fit[[var_name]] <- y_temp
      
      # 4. 拟合 Tobit 模型
      tryCatch({
        # left 和 right 参数定义了截断阈值
        fit <- tobit(as.formula(formula_str), left = lower_val, right = upper_val, data = data_for_fit)
        
        # 提取参数
        sigma_hat <- fit$scale
        # 获取线性预测值 (mu = X * beta)
        # predict() 默认返回 linear predictor
        mu_hat <- predict(fit, newdata = data_for_fit)
        
        # 5. 从截断正态分布中采样 (Stochastic Imputation)
        
        # >>> Left Censored <<<
        if (sum(mask_lower) > 0) {
          # 在 (-Inf, Lower] 采样
          draws <- rtruncnorm(sum(mask_lower), 
                              a = -Inf, 
                              b = lower_val, 
                              mean = mu_hat[mask_lower], 
                              sd = sigma_hat)
          df_filled[[var_name]][mask_lower] <- draws
        }
        
        # >>> Right Censored <<<
        if (sum(mask_upper) > 0) {
          # 在 [Upper, Inf) 采样
          draws <- rtruncnorm(sum(mask_upper), 
                              a = upper_val, 
                              b = Inf, 
                              mean = mu_hat[mask_upper], 
                              sd = sigma_hat)
          df_filled[[var_name]][mask_upper] <- draws
        }
        
      }, error = function(e) {
        # 如果回归失败 (例如共线性)，保持上一次迭代的值或均值
        # 这里不做操作，直接跳过
      })
    }
  }
  
  return(df_filled)
}

# ============================================================
# 主执行循环
# ============================================================

M_IMP <- 5  # 生成 5 个数据集
MAXIT <- 5  # 每个数据集内部迭代 5 次以稳定结果

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  cutoff_key <- cutoff_keys[i]
  
  # 读取数据
  csv_file <- list.files(input_dirs[[i]], pattern = "*.csv", full.names = TRUE)[1]
  if(is.na(csv_file)) next
  
  message(paste("处理:", level_name, "| Reading:", basename(csv_file)))
  df_orig <- read.csv(csv_file)
  
  current_cutoffs <- cutoff_list[[cutoff_key]]
  
  # 存储 m 次插补的结果
  imputed_datasets <- list()
  
  message(paste("  -> Running Tobit Imputation (m =", M_IMP, ")..."))
  
  for (m in 1:M_IMP) {
    # 每次使用不同的随机种子
    df_imp <- impute_tobit_single_run(df_orig, current_cutoffs, max_iter = MAXIT, seed_val = 1000 + m)
    imputed_datasets[[m]] <- df_imp
  }
  
  # --------------------------------------------------------
  # Pooling: 取平均生成单一数据集
  # --------------------------------------------------------
  message("  -> Pooling (Averaging) datasets...")
  
  # 初始化为第一个数据集的结构
  df_pooled <- imputed_datasets[[1]]
  cols_numeric <- names(df_pooled)[sapply(df_pooled, is.numeric)]
  
  # 对数值列取平均
  for (col in cols_numeric) {
    # 提取所有 m 个数据集中该列的值，组成矩阵 (n x m)
    vals_matrix <- sapply(imputed_datasets, function(d) d[[col]])
    # 计算行均值
    df_pooled[[col]] <- rowMeans(vals_matrix)
  }
  
  # 保存最终结果
  save_path <- file.path(output_dir, paste0("tobit_averaged_", level_name, ".csv"))
  write.csv(df_pooled, save_path, row.names = FALSE)
  
  message(paste("     [Saved]:", basename(save_path)))
}

message("所有 Tobit 填充任务完成！")


##########################单次填充，对MICE和TOBIT，更新时间09072026############################
# ============================================================

# 1. 创建 10% 截断的阈值表
cutoff_10 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.1,
  Lower_Cutoff = c(20.061273, 200.055249, 8.495351, 19.000000, 155.000000, 0.630000, 20.662977, 2.811597),
  Upper_Cutoff = c(36.500505, 3093.325368, 43.406748, 34.000000, 172.700000, 1.052760, 76.117738, 7.467111),
  stringsAsFactors = FALSE
)

# 2. 创建 20% 截断的阈值表
cutoff_20 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.2,
  Lower_Cutoff = c(21.5316497, 341.2159561, 11.8236356, 21.0000000, 158.0000000, 0.7862489, 26.5367602, 3.3568065),
  Upper_Cutoff = c(31.4288844, 2044.8760798, 35.2850673, 32.0000000, 170.0000000, 0.9451719, 60.6895559, 6.4242105),
  stringsAsFactors = FALSE
)

# 3. 创建 30% 截断的阈值表
cutoff_30 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.3,
  Lower_Cutoff = c(22.6611810, 508.1810942, 14.1212004, 23.0000000, 160.2000000, 0.8218406, 30.6946994, 3.8344199),
  Upper_Cutoff = c(28.5752868, 1518.1500077, 28.0546296, 30.0000000, 167.5000000, 0.9198789, 53.0351045, 5.7293256),
  stringsAsFactors = FALSE
)

# 4. 创建 40% 截断的阈值表
cutoff_40 <- data.frame(
  Variable = c("BMI", "PAPPA_V1", "fbHCG_V1", "Age_at_V1", "Height", "U3CD07", "PLGF_V1", "ADAM12_V1"),
  Truncation_Level = 0.4,
  Lower_Cutoff = c(23.8219026, 711.0014683, 16.8910435, 25.0000000, 162.5000000, 0.8434391, 35.6720577, 4.2460000),
  Upper_Cutoff = c(26.7537356, 1229.4455148, 23.5820478, 29.0000000, 165.1000000, 0.8948467, 45.8456378, 5.1540000),
  stringsAsFactors = FALSE
)

# 5. 将这四个表打包进一个 List
cutoff_list <- list(
  df_trunc_10 = cutoff_10,
  df_trunc_20 = cutoff_20,
  df_trunc_30 = cutoff_30,
  df_trunc_40 = cutoff_40
)

#################### MICE imputation (SINGLE imputation) ####################
library(mice)
library(truncnorm)
library(dplyr)

# 依赖：环境中需已定义 cutoff_list（各 truncation 级别的 cutoff 信息）

# ============================================================
# 辅助函数：获取每个缺失值的动态边界 [L_i, U_i]
# ============================================================
get_dynamic_bounds <- function(y_name, n_rows, data_df, cutoff_info) {
  var_cutoffs <- cutoff_info %>% filter(Variable == y_name)
  if (nrow(var_cutoffs) == 0) {
    return(list(lower = rep(-Inf, n_rows), upper = rep(Inf, n_rows)))
  }
  l_limit <- var_cutoffs$Lower_Cutoff
  u_limit <- var_cutoffs$Upper_Cutoff
  ind_col <- paste0(y_name, "_ind")
  
  lower_vec <- rep(-Inf, n_rows)
  upper_vec <- rep(Inf,  n_rows)
  
  if (ind_col %in% names(data_df)) {
    indicators <- data_df[[ind_col]]
    # 左删失 (Below_Lower)：潜在真值 <= 下限  -> 区间 (-Inf, L]
    idx_below <- which(indicators == "Below_Lower")
    upper_vec[idx_below] <- l_limit
    # 右删失 (Above_Upper)：潜在真值 >= 上限  -> 区间 [U, Inf)
    idx_above <- which(indicators == "Above_Upper")
    lower_vec[idx_above] <- u_limit
  }
  return(list(lower = lower_vec, upper = upper_vec))
}

# ============================================================
# 方法 1: Post-hoc Truncation
# ============================================================
mice.impute.posthoc <- function(y, ry, x, wy = NULL, ...) {
  imputed_values <- mice.impute.norm(y, ry, x, wy = wy, ...)
  
  bounds <- get_dynamic_bounds(CURRENT_VAR_NAME, length(y), global_df_ref, global_cutoff_info)
  l_bounds_miss <- bounds$lower[wy]
  u_bounds_miss <- bounds$upper[wy]
  
  # 左删失: (-Inf, L] -> 预测超过 L 的强设为 L
  # 右删失: [U, Inf)  -> 预测低于 U 的强设为 U
  imputed_values <- pmax(imputed_values, l_bounds_miss)
  imputed_values <- pmin(imputed_values, u_bounds_miss)
  return(imputed_values)
}

# ============================================================
# 方法 2: Truncated Normal Sampling
# ============================================================
mice.impute.truncnorm <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  x <- cbind(1, as.matrix(x))
  parm <- .norm.draw(y, ry, x)
  
  x_miss  <- x[wy, , drop = FALSE]
  mu_miss <- x_miss %*% parm$beta
  sigma   <- parm$sigma
  
  bounds <- get_dynamic_bounds(CURRENT_VAR_NAME, length(y), global_df_ref, global_cutoff_info)
  l_bounds_miss <- bounds$lower[wy]
  u_bounds_miss <- bounds$upper[wy]
  
  imputed_values <- rtruncnorm(sum(wy),
                               a = l_bounds_miss,
                               b = u_bounds_miss,
                               mean = mu_miss,
                               sd = sigma)
  return(imputed_values)
}

.norm.draw <- mice:::.norm.draw

# ============================================================
# 闭包：自动反推当前正在插补的变量名
# ============================================================
make_imputer <- function(method_type) {
  function(y, ry, x, wy = NULL, ...) {
    var_name <- NA
    for (col in names(global_df_ref)) {
      if (is.numeric(global_df_ref[[col]]) &&
          length(y) == length(global_df_ref[[col]]) &&
          all(y[ry] == global_df_ref[[col]][ry])) {
        var_name <- col
        break
      }
    }
    assign("CURRENT_VAR_NAME", var_name, envir = .GlobalEnv)
    if (method_type == "posthoc") {
      return(mice.impute.posthoc(y, ry, x, wy, ...))
    } else {
      return(mice.impute.truncnorm(y, ry, x, wy, ...))
    }
  }
}

# ============================================================
# 主执行循环（单次填充）
# ============================================================
input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)
cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

out_mice_posthoc <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_Posthoc"
out_mice_trunc   <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_MICE_TruncNormal"
if (!dir.exists(out_mice_posthoc)) dir.create(out_mice_posthoc)
if (!dir.exists(out_mice_trunc))   dir.create(out_mice_trunc)

M_IMP <- 1   # <<< 单次填充（原来是 5）
MAXIT <- 5   # 链式方程迭代次数（收敛用，保持不变）
SEED  <- 123 # 保证可复现

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  cutoff_key <- cutoff_keys[i]
  
  csv_file <- list.files(input_dirs[[i]], pattern = "*.csv", full.names = TRUE)[1]
  if (is.na(csv_file)) next
  
  message(paste("处理:", level_name, "| Reading:", basename(csv_file)))
  df_orig <- read.csv(csv_file)
  
  global_cutoff_info <<- cutoff_list[[cutoff_key]]
  global_df_ref      <<- df_orig
  
  ini  <- mice(df_orig, maxit = 0, print = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  target_vars <- global_cutoff_info$Variable
  
  #ind_cols <- grep("_ind$", names(df_orig), value = TRUE)
  #pred[, ind_cols] <- 0
  #meth[ind_cols]   <- ""
  df_orig$true_prob <- NULL                     # 彻底删掉上帝变量
  ini  <- mice(df_orig, maxit = 0, print = FALSE)  # 删完之后再 init（pred 矩阵才不含 true_prob）
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  ind_cols <- grep("_ind$", names(df_orig), value = TRUE)
  pred[, c(ind_cols, "outcome")] <- 0           # _ind 和 outcome 都不许当 predictor
  meth[ind_cols] <- ""                          # _ind 不插补；outcome 没缺失，mice 不会动它
  
  # ---------- 方法 1: MICE Post-hoc（单次）----------
  meth_posthoc <- meth
  mice.impute.posthoc_wrapper <<- make_imputer("posthoc")
  meth_posthoc[target_vars] <- "posthoc_wrapper"
  
  message("  -> Running MICE Post-hoc (m=1)...")
  imp_posthoc <- mice(df_orig, method = meth_posthoc, predictorMatrix = pred,
                      m = M_IMP, maxit = MAXIT, printFlag = FALSE, seed = SEED)
  
  df_posthoc_single <- complete(imp_posthoc, 1)   # <<< 直接取唯一的填充数据集
  save_path <- file.path(out_mice_posthoc, paste0("mice_posthoc_single_", level_name, ".csv"))
  write.csv(df_posthoc_single, save_path, row.names = FALSE)
  message(paste("     [Saved]:", basename(save_path)))
  
  # ---------- 方法 2: MICE Truncated Normal（单次）----------
  meth_trunc <- meth
  mice.impute.trunc_wrapper <<- make_imputer("trunc")
  meth_trunc[target_vars] <- "trunc_wrapper"
  
  message("  -> Running MICE Truncated Normal (m=1)...")
  imp_trunc <- mice(df_orig, method = meth_trunc, predictorMatrix = pred,
                    m = M_IMP, maxit = MAXIT, printFlag = FALSE, seed = SEED)
  
  df_trunc_single <- complete(imp_trunc, 1)        # <<< 直接取唯一的填充数据集
  save_path <- file.path(out_mice_trunc, paste0("mice_truncnorm_single_", level_name, ".csv"))
  write.csv(df_trunc_single, save_path, row.names = FALSE)
  message(paste("     [Saved]:", basename(save_path)))
}

message("所有 MICE (Single) 填充任务完成！")
#################### TOBIT imputation (SINGLE imputation) ####################
library(dplyr)
library(readr)
library(AER)
library(truncnorm)

# ============================================================
# 0. 路径设置
# ============================================================
input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)
cutoff_keys <- c("df_trunc_10", "df_trunc_20", "df_trunc_30", "df_trunc_40")

output_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/Imputed_Tobit"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ============================================================
# 核心函数：单次 Tobit 插补（内部含 Gibbs 式迭代稳定）
# ============================================================
impute_tobit_single_run <- function(df, cutoff_info, max_iter = 5, seed_val = 123) {
  
  set.seed(seed_val)
  
  df_filled   <- df
  target_vars <- cutoff_info$Variable
  
  # 初始化：均值填充作为回归起点（predictor 不能有 NA）
  for (v in names(df_filled)) {
    if (is.numeric(df_filled[[v]]) && any(is.na(df_filled[[v]]))) {
      df_filled[[v]][is.na(df_filled[[v]])] <- mean(df_filled[[v]], na.rm = TRUE)
    }
  }
  
  # 迭代优化（Gibbs sampling style，收敛用，保持 max_iter）
  for (iter in 1:max_iter) {
    for (j in 1:nrow(cutoff_info)) {
      
      var_name  <- cutoff_info$Variable[j]
      lower_val <- cutoff_info$Lower_Cutoff[j]
      upper_val <- cutoff_info$Upper_Cutoff[j]
      ind_col   <- paste0(var_name, "_ind")
      
      if (!var_name %in% names(df)) next
      
      orig_inds  <- df[[ind_col]]
      mask_lower <- (orig_inds == "Below_Lower")
      mask_upper <- (orig_inds == "Above_Upper")
      mask_lower[is.na(mask_lower)] <- FALSE
      mask_upper[is.na(mask_upper)] <- FALSE
      
      if (sum(mask_lower) == 0 && sum(mask_upper) == 0) next
      
      # Tobit 要求 censored 点设为边界值
      y_temp <- df_filled[[var_name]]
      y_temp[mask_lower] <- lower_val
      y_temp[mask_upper] <- upper_val
      
      #predictors <- setdiff(names(df_filled),
                            #c(var_name, grep("_ind$", names(df_filled), value = TRUE)))
      #predictors <- predictors[sapply(df_filled[predictors], is.numeric)]
      
      predictors <- setdiff(names(df_filled),
                            c(var_name, "outcome", "true_prob",
                              grep("_ind$", names(df_filled), value = TRUE)))
      predictors <- predictors[sapply(df_filled[predictors], is.numeric)]
      
      formula_str <- paste(var_name, "~", paste(predictors, collapse = " + "))
      
      data_for_fit <- df_filled
      data_for_fit[[var_name]] <- y_temp
      
      tryCatch({
        fit <- tobit(as.formula(formula_str), left = lower_val, right = upper_val,
                     data = data_for_fit)
        
        sigma_hat <- fit$scale
        mu_hat    <- predict(fit, newdata = data_for_fit)
        
        # 从截断正态抽样（stochastic imputation）
        if (sum(mask_lower) > 0) {
          draws <- rtruncnorm(sum(mask_lower), a = -Inf, b = lower_val,
                              mean = mu_hat[mask_lower], sd = sigma_hat)
          df_filled[[var_name]][mask_lower] <- draws
        }
        if (sum(mask_upper) > 0) {
          draws <- rtruncnorm(sum(mask_upper), a = upper_val, b = Inf,
                              mean = mu_hat[mask_upper], sd = sigma_hat)
          df_filled[[var_name]][mask_upper] <- draws
        }
      }, error = function(e) {
        # 拟合失败：保持上一轮的值，跳过
      })
    }
  }
  
  return(df_filled)
}

# ============================================================
# 主执行循环（单次填充，无 pooling）
# ============================================================
MAXIT <- 5   # 单个数据集内部迭代次数（收敛用，保持不变）
SEED  <- 1000

for (i in seq_along(input_dirs)) {
  
  level_name <- names(input_dirs)[i]
  cutoff_key <- cutoff_keys[i]
  
  csv_file <- list.files(input_dirs[[i]], pattern = "*.csv", full.names = TRUE)[1]
  if (is.na(csv_file)) next
  
  message(paste("处理:", level_name, "| Reading:", basename(csv_file)))
  df_orig <- read.csv(csv_file)
  
  current_cutoffs <- cutoff_list[[cutoff_key]]
  
  message("  -> Running Tobit Imputation (single)...")
  df_imp <- impute_tobit_single_run(df_orig, current_cutoffs,
                                    max_iter = MAXIT, seed_val = SEED)
  
  save_path <- file.path(output_dir, paste0("tobit_single_", level_name, ".csv"))
  write.csv(df_imp, save_path, row.names = FALSE)
  message(paste("     [Saved]:", basename(save_path)))
}

message("所有 Tobit (Single) 填充任务完成！")
############################建模和评估###########################

# ============================================================
# 1. 配置 & 文件映射
# ============================================================

# !!! 请确认你的路径 !!!
work_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"

# 15个变量
predictors <- c(
  "U3CC07", "U3CD07", "Height", "BMI", "SFLT1_V1", 
  "fbHCG_V1", "V2BA01_KG", "Age_at_V1", "ADAM12_V1", "AFP_V1", 
  "PLGF_V1", "INHIBINA_V1", "U2BB02", "PAPPA_V1", "SFLT1_V2"
)

outcome_col <- "outcome"

# 文件名映射表
file_pattern_map <- list(
  "simulated_dataset_final" = "Complete_Data",
  "imputed_RandomBase" = "Random",
  "dataset_trunc" = "Raw_Trunc",
  "imputed_LOD_Half" = "LOD_Half_Sub",
  "imputed_LOD_" = "LOD_Fill",
  "imputed_Gaussian" = "Gauss_Draw",
  "imputed_Lognormal" = "Lognorm_Draw",
  "mice_posthoc" = "MICE_Posthoc",
  "mice_truncnorm" = "MICE_Trunc",
  "tobit_averaged" = "Tobit_Avg",
  "imputed_DB_Mean" = "DB_Mean",
  "imputed_DB_Draw" = "DB_Draw"
)

# ============================================================
# 2. 核心函数定义
# ============================================================

# --- 2.1 辅助函数：计算指标 (AUC + Accuracy + CI) ---
calculate_metrics_with_ci <- function(preds_df, threshold = 0.5) {
  # AUC & DeLong CI
  roc_obj <- roc(preds_df$truth, preds_df$prob, quiet = TRUE)
  ci_res <- ci.auc(roc_obj, method = "delong")
  
  # Accuracy & Binomial CI
  pred_class <- ifelse(preds_df$prob >= threshold, 1, 0)
  pred_factor  <- factor(pred_class, levels = c(0, 1))
  truth_factor <- factor(preds_df$truth, levels = c(0, 1))
  cm <- confusionMatrix(pred_factor, truth_factor, positive = "1")
  
  return(c(
    AUC = as.numeric(ci_res[2]), 
    AUC_Lower = as.numeric(ci_res[1]), 
    AUC_Upper = as.numeric(ci_res[3]),
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Acc_Lower = as.numeric(cm$overall["AccuracyLower"]),
    Acc_Upper = as.numeric(cm$overall["AccuracyUpper"])
  ))
}

# --- 2.2 XGBoost 训练函数 (Top 10 Frequency + Score) ---
train_xgboost_cv <- function(df, predictors, outcome_col, n_folds = 5) {
  X <- as.matrix(df[, predictors])
  y <- as.numeric(df[[outcome_col]])
  
  set.seed(123)
  folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
  pooled_preds <- data.frame(row_id = integer(), prob = numeric(), truth = numeric())
  
  # 5-Fold CV
  for (f in 1:n_folds) {
    idx_test <- folds[[f]]
    idx_train <- setdiff(1:nrow(df), idx_test)
    
    dtrain <- xgb.DMatrix(data = X[idx_train, ], label = y[idx_train])
    dtest  <- xgb.DMatrix(data = X[idx_test, ], label = y[idx_test])
    
    params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05, max_depth = 4)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
    preds <- predict(model, dtest)
    
    pooled_preds <- rbind(pooled_preds, data.frame(row_id = idx_test, prob = preds, truth = y[idx_test]))
  }
  pooled_preds <- pooled_preds[order(pooled_preds$row_id), ]
  
  # --- 修改 Bootstrap 部分 ---
  message("   -> XGBoost Feature Importance with Shadow Feature Threshold (B=1000)...")
  B_BOOT <- 1000
  imp_list <- list()
  
  # 1. 定义包含 Shadow 的预测变量列表
  # 我们不改变原始 predictors，只在训练时临时加一列
  
  X_raw <- as.matrix(df[, predictors])
  y_raw <- as.numeric(df[[outcome_col]])
  
  for (b in 1:B_BOOT) {
    # 2. 重采样索引
    boot_idx <- sample(1:nrow(X_raw), nrow(X_raw), replace = TRUE)
    
    # 3. 构造数据：添加一列 "Shadow" (纯随机噪音)
    # 这一列完全是随机数，没有任何预测能力
    X_boot <- X_raw[boot_idx, ]
    y_boot <- y_raw[boot_idx]
    
    # 添加 Shadow 列
    shadow_col <- runif(nrow(X_boot))
    X_boot_with_shadow <- cbind(X_boot, Shadow = shadow_col)
    
    # 4. 训练 XGBoost
    dboot <- xgb.DMatrix(data = X_boot_with_shadow, label = y_boot)
    # 保持参数一致
    params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05, max_depth = 4)
    bst_boot <- xgb.train(params = params, data = dboot, nrounds = 100, verbose = 0)
    
    # 5. 提取重要性
    imp_df <- xgb.importance(model = bst_boot)
    
    # 6. 【关键步骤】找到 Shadow 变量的重要性分数
    # 如果 Shadow 没被用到（Importance里没有它），那它的分数就是 0
    shadow_row <- imp_df %>% filter(Feature == "Shadow")
    shadow_threshold <- if(nrow(shadow_row) > 0) shadow_row$Gain else 0
    
    # 7. 定义 "Selected"：Gain 必须严格大于 Shadow 的 Gain
    full_imp <- data.frame(Feature = predictors) # 原始变量列表
    
    if (!is.null(imp_df)) {
      full_imp <- left_join(full_imp, imp_df, by = "Feature")
    }
    full_imp[is.na(full_imp)] <- 0 # 没出现的变量 Gain 补 0
    
    # 这里的 Total 通常用 Gain 就够了，或者 Gain + Cover
    # 建议只用 Gain，因为这代表它对 Loss 的贡献
    full_imp$Metric_Used <- full_imp$Gain 
    
    # 动态筛选：只有打败了随机噪音的，才算被选中
    full_imp$Selected <- ifelse(full_imp$Metric_Used > shadow_threshold, 1, 0)
    
    imp_list[[b]] <- full_imp
  }
  
  imp_combined <- bind_rows(imp_list)
  
  final_imp <- imp_combined %>% 
    group_by(Feature) %>%
    summarise(
      Score_Mean = mean(Total), 
      Score_Lower = quantile(Total, 0.025), 
      Score_Upper = quantile(Total, 0.975),
      Selection_Frequency = mean(Selected) 
    ) %>%
    arrange(desc(Score_Mean))
  
  return(list(preds = pooled_preds, importance = final_imp))
}

# --- 2.3 Elastic Net 训练函数 (True Coef + Frequency) ---
train_glmnet_cv <- function(df, predictors, outcome_col, n_folds = 5) {
  X_raw <- df[, predictors]
  X_scaled <- scale(X_raw)
  X_scaled[is.na(X_scaled)] <- 0 
  y <- as.numeric(df[[outcome_col]])
  
  set.seed(123)
  folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
  pooled_preds <- data.frame(row_id = integer(), prob = numeric(), truth = numeric())
  
  # Part A: CV for Refit Performance
  for (f in 1:n_folds) {
    idx_test <- folds[[f]]
    idx_train <- setdiff(1:nrow(df), idx_test)
    X_train <- X_scaled[idx_train, ]; y_train <- y[idx_train]
    X_test  <- X_scaled[idx_test, ]
    
    # 1. Select
    cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 0.5, type.measure = "auc")
    coefs <- coef(cv_fit, s = cv_fit$lambda.1se)
    selected_vars <- rownames(coefs)[which(coefs != 0)]
    selected_vars <- selected_vars[selected_vars != "(Intercept)"]
    
    # 2. Refit
    if (length(selected_vars) == 0) {
      preds <- rep(mean(y_train), length(idx_test))
    } else {
      df_train_sel <- data.frame(y = y_train, X_train[, selected_vars, drop=FALSE])
      df_test_sel  <- data.frame(X_test[, selected_vars, drop=FALSE])
      fit_refit <- glm(y ~ ., data = df_train_sel, family = binomial)
      preds <- predict(fit_refit, newdata = df_test_sel, type = "response")
    }
    pooled_preds <- rbind(pooled_preds, data.frame(row_id = idx_test, prob = preds, truth = y[idx_test]))
  }
  pooled_preds <- pooled_preds[order(pooled_preds$row_id), ]
  
  # Feature Importance Bootstrap
  # 修复：设置 B=1000
  message("   -> GLMNet Feature Importance & Frequency Bootstrap (B=1000)...")
  B_BOOT <- 1000
  imp_list <- list()
  
  for (b in 1:B_BOOT) {
    boot_idx <- sample(1:nrow(X_scaled), nrow(X_scaled), replace = TRUE)
    cv_boot <- cv.glmnet(X_scaled[boot_idx, ], y[boot_idx], family = "binomial", alpha = 0.5)
    coef_boot <- as.matrix(coef(cv_boot, s = cv_boot$lambda.1se))
    
    df_coef <- data.frame(Feature = rownames(coef_boot), Coef = as.numeric(coef_boot)) %>%
      filter(Feature != "(Intercept)") 
    
    df_coef$Selected <- ifelse(df_coef$Coef != 0, 1, 0)
    
    imp_list[[b]] <- df_coef
  }
  
  imp_combined <- bind_rows(imp_list)
  
  final_imp <- imp_combined %>%
    group_by(Feature) %>%
    summarise(
      Score_Mean = mean(Coef),
      Score_Lower = quantile(Coef, 0.025),
      Score_Upper = quantile(Coef, 0.975),
      Selection_Frequency = mean(Selected)
    ) %>%
    right_join(data.frame(Feature = predictors), by = "Feature") %>%
    replace_na(list(Score_Mean = 0, Score_Lower = 0, Score_Upper = 0, Selection_Frequency = 0)) %>%
    arrange(desc(abs(Score_Mean)))
  
  return(list(preds = pooled_preds, importance = final_imp))
}

# --- 2.4 Agreement & DeLong 计算函数 (修复命名匹配问题) ---
calculate_agreement_pair <- function(probs1, probs2, threshold = 0.5, B = 1000) {
  y1 <- ifelse(probs1 >= threshold, 1, 0); y2 <- ifelse(probs2 >= threshold, 1, 0)
  n <- length(y1)
  boot_vals <- numeric(B)
  set.seed(123)
  for (b in 1:B) {
    idx <- sample(1:n, n, replace = TRUE)
    boot_vals[b] <- mean(y1[idx] == y2[idx])
  }
  # 修复：确保返回名字是规范的，不依赖 quantile 的默认命名
  q_low <- quantile(boot_vals, 0.025)
  q_up  <- quantile(boot_vals, 0.975)
  return(c(Agreement = mean(y1 == y2), CI_Lower = as.numeric(q_low), CI_Upper = as.numeric(q_up)))
}

compute_agreement_matrices_with_ci <- function(preds_list) {
  methods <- names(preds_list); n <- length(methods)
  mat_mean  <- matrix(NA, n, n, dimnames = list(methods, methods))
  mat_lower <- matrix(NA, n, n, dimnames = list(methods, methods))
  mat_upper <- matrix(NA, n, n, dimnames = list(methods, methods))
  
  # 修复：使用 B=1000
  message("   Running pairwise agreement bootstrap (B=1000)...")
  
  for(i in 1:n) { for(j in 1:n) { if(i <= j) {
    p1 <- preds_list[[methods[i]]]; p2 <- preds_list[[methods[j]]]
    common_ids <- intersect(p1$row_id, p2$row_id)
    if(length(common_ids) > 0) {
      sub1 <- p1[p1$row_id %in% common_ids, ]; sub2 <- p2[p2$row_id %in% common_ids, ]
      
      res <- calculate_agreement_pair(sub1$prob, sub2$prob, B = 1000)
      
      mat_mean[i, j]  <- mat_mean[j, i]  <- res["Agreement"]
      # 修复：现在名字匹配了
      mat_lower[i, j] <- mat_lower[j, i] <- res["CI_Lower"] 
      mat_upper[i, j] <- mat_upper[j, i] <- res["CI_Upper"] 
    }
  }}}
  return(list(Mean = mat_mean, Lower = mat_lower, Upper = mat_upper))
}

compute_delong_matrices <- function(preds_list) {
  methods <- names(preds_list); n <- length(methods)
  mat_p <- matrix(NA, n, n, dimnames = list(methods, methods))
  
  for(i in 1:n) { for(j in 1:n) {
    if(i < j) {
      p1 <- preds_list[[methods[i]]]; p2 <- preds_list[[methods[j]]]
      common_ids <- intersect(p1$row_id, p2$row_id)
      if(length(common_ids) > 0) {
        sub1 <- p1[p1$row_id %in% common_ids, ]; sub2 <- p2[p2$row_id %in% common_ids, ]
        t_res <- roc.test(roc(sub1$truth, sub1$prob, quiet=T), roc(sub2$truth, sub2$prob, quiet=T), method="delong")
        mat_p[i, j] <- mat_p[j, i] <- t_res$p.value
      }
    } else if (i == j) { mat_p[i, j] <- 1.0 }
  }}
  return(mat_p)
}

# ============================================================
# 3. 主程序循环
# ============================================================

all_files <- list.files(work_dir, pattern = "*.csv", full.names = TRUE)
results_xgboost_perf <- list(); results_glmnet_perf  <- list()
results_xgboost_imp <- list(); results_glmnet_imp   <- list()
preds_storage_xgb    <- list(); preds_storage_glm    <- list()

message(paste("Found", length(all_files), "files. Starting..."))

for (file_path in all_files) {
  file_name <- basename(file_path)
  method_label <- NA
  for (pattern in names(file_pattern_map)) {
    if (grepl(pattern, file_name, fixed = TRUE)) {
      method_label <- file_pattern_map[[pattern]]
      break 
    }
  }
  if (is.na(method_label)) next
  
  message(paste("Processing:", method_label))
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # --- XGBoost ---
  res_xgb <- train_xgboost_cv(df, predictors, outcome_col)
  metrics_xgb <- calculate_metrics_with_ci(res_xgb$preds)
  
  results_xgboost_perf[[method_label]] <- metrics_xgb
  res_xgb$importance$Method <- method_label
  results_xgboost_imp[[method_label]] <- res_xgb$importance
  preds_storage_xgb[[method_label]] <- res_xgb$preds
  
  # --- Elastic Net (Skip Raw_Trunc) ---
  if (method_label != "Raw_Trunc") {
    res_glm <- train_glmnet_cv(df, predictors, outcome_col)
    metrics_glm <- calculate_metrics_with_ci(res_glm$preds)
    
    results_glmnet_perf[[method_label]] <- metrics_glm
    res_glm$importance$Method <- method_label
    results_glmnet_imp[[method_label]] <- res_glm$importance
    preds_storage_glm[[method_label]] <- res_glm$preds 
  }
}

# ============================================================
# 4. 统计计算 (Agreement CI & DeLong)
# ============================================================
message("Calculating Agreement Matrices (with CI) and DeLong Tests...")
agree_list_xgb <- compute_agreement_matrices_with_ci(preds_storage_xgb)
delong_mat_xgb <- compute_delong_matrices(preds_storage_xgb)

agree_list_glm <- compute_agreement_matrices_with_ci(preds_storage_glm)
delong_mat_glm <- compute_delong_matrices(preds_storage_glm)

# ============================================================
# 5. 保存结果
# ============================================================
res_dir <- file.path(work_dir, "Evaluation_Results_Full")
if(!dir.exists(res_dir)) dir.create(res_dir)

message(paste("Saving all results to:", res_dir))

# Performance
write.csv(bind_rows(results_xgboost_perf, .id = "Method"), file.path(res_dir, "Performance_XGBoost.csv"), row.names = FALSE)
write.csv(bind_rows(results_glmnet_perf, .id = "Method"), file.path(res_dir, "Performance_ElasticNet.csv"), row.names = FALSE)

# Importance
write.csv(bind_rows(results_xgboost_imp), file.path(res_dir, "Importance_XGBoost.csv"), row.names = FALSE)
write.csv(bind_rows(results_glmnet_imp), file.path(res_dir, "Importance_ElasticNet.csv"), row.names = FALSE)

# Agreement
write.csv(agree_list_xgb$Mean,  file.path(res_dir, "Agreement_Mean_XGBoost.csv"))
write.csv(agree_list_xgb$Lower, file.path(res_dir, "Agreement_Lower_XGBoost.csv"))
write.csv(agree_list_xgb$Upper, file.path(res_dir, "Agreement_Upper_XGBoost.csv"))

write.csv(agree_list_glm$Mean,  file.path(res_dir, "Agreement_Mean_ElasticNet.csv"))
write.csv(agree_list_glm$Lower, file.path(res_dir, "Agreement_Lower_ElasticNet.csv"))
write.csv(agree_list_glm$Upper, file.path(res_dir, "Agreement_Upper_ElasticNet.csv"))

# DeLong P-values
write.csv(delong_mat_xgb, file.path(res_dir, "DeLong_Pvalues_XGBoost.csv"))
write.csv(delong_mat_glm, file.path(res_dir, "DeLong_Pvalues_ElasticNet.csv"))

message("✅ All results saved successfully!")





































# ============================================================
# 极速补救版：只算 Agreement CI 和 DeLong Test
# ============================================================

library(dplyr)
library(readr)
library(xgboost)
library(glmnet)
library(pROC)
library(caret)
library(tidyr)

# !!! 修改你的路径 !!!
work_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"

# 15个变量
predictors <- c(
  "U3CC07", "U3CD07", "Height", "BMI", "SFLT1_V1", 
  "fbHCG_V1", "V2BA01_KG", "Age_at_V1", "ADAM12_V1", "AFP_V1", 
  "PLGF_V1", "INHIBINA_V1", "U2BB02", "PAPPA_V1", "SFLT1_V2"
)
outcome_col <- "outcome"

# 映射表
file_pattern_map <- list(
  "simulated_dataset_final" = "Complete_Data",
  "imputed_RandomBase" = "Random",
  "dataset_trunc" = "Raw_Trunc",
  "imputed_LOD_Half" = "LOD_Half_Sub",
  "imputed_LOD_" = "LOD_Fill",
  "imputed_Gaussian" = "Gauss_Draw",
  "imputed_Lognormal" = "Lognorm_Draw",
  "mice_posthoc" = "MICE_Posthoc",
  "mice_truncnorm" = "MICE_Trunc",
  "tobit_averaged" = "Tobit_Avg",
  "imputed_DB_Mean" = "DB_Mean",
  "imputed_DB_Draw" = "DB_Draw"
)

# ============================================================
# 1. 极速训练函数 (只为了拿预测值，不须算Importance)
# ============================================================

# XGBoost: 只做 5-Fold CV，不做 Bootstrap
train_xgb_fast <- function(df, predictors, outcome_col, n_folds = 5) {
  X <- as.matrix(df[, predictors])
  y <- as.numeric(df[[outcome_col]])
  
  set.seed(123)
  folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
  pooled_preds <- data.frame()
  
  for (f in 1:n_folds) {
    idx_test <- folds[[f]]
    idx_train <- setdiff(1:nrow(df), idx_test)
    
    dtrain <- xgb.DMatrix(data = X[idx_train, ], label = y[idx_train])
    dtest  <- xgb.DMatrix(data = X[idx_test, ], label = y[idx_test])
    
    params <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05, max_depth = 4)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
    
    preds <- predict(model, dtest)
    pooled_preds <- rbind(pooled_preds, data.frame(row_id = idx_test, prob = preds, truth = y[idx_test]))
  }
  return(pooled_preds[order(pooled_preds$row_id), ])
}

# Elastic Net: 只做 5-Fold CV
train_glm_fast <- function(df, predictors, outcome_col, n_folds = 5) {
  X_raw <- df[, predictors]
  X_scaled <- scale(X_raw); X_scaled[is.na(X_scaled)] <- 0 
  y <- as.numeric(df[[outcome_col]])
  
  set.seed(123)
  folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
  pooled_preds <- data.frame()
  
  for (f in 1:n_folds) {
    idx_test <- folds[[f]]
    idx_train <- setdiff(1:nrow(df), idx_test)
    
    cv_fit <- cv.glmnet(X_scaled[idx_train,], y[idx_train], family="binomial", alpha=0.5)
    coefs <- coef(cv_fit, s = cv_fit$lambda.1se)
    sel_vars <- rownames(coefs)[which(coefs != 0)]; sel_vars <- sel_vars[sel_vars != "(Intercept)"]
    
    if (length(sel_vars) == 0) {
      preds <- rep(mean(y[idx_train]), length(idx_test))
    } else {
      df_train <- data.frame(y=y[idx_train], X_scaled[idx_train, sel_vars, drop=FALSE])
      df_test  <- data.frame(X_scaled[idx_test, sel_vars, drop=FALSE])
      fit_refit <- glm(y ~ ., data = df_train, family = binomial)
      preds <- predict(fit_refit, newdata = df_test, type = "response")
    }
    pooled_preds <- rbind(pooled_preds, data.frame(row_id = idx_test, prob = preds, truth = y[idx_test]))
  }
  return(pooled_preds[order(pooled_preds$row_id), ])
}

# ============================================================
# 2. 主循环：快速生成预测值
# ============================================================

all_files <- list.files(work_dir, pattern = "*.csv", full.names = TRUE)
preds_storage_xgb <- list()
preds_storage_glm <- list()

message(">>> 1. Starting Fast Prediction Generation (Skipping Importance)...")

for (file_path in all_files) {
  file_name <- basename(file_path)
  method_label <- NA
  for (p in names(file_pattern_map)) { if (grepl(p, file_name, fixed=TRUE)) { method_label <- file_pattern_map[[p]]; break } }
  
  if (is.na(method_label)) next
  message(paste("Processing:", method_label))
  
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # XGBoost
  preds_storage_xgb[[method_label]] <- train_xgb_fast(df, predictors, outcome_col)
  
  # Elastic Net (跳过 Raw_Trunc)
  if (method_label != "Raw_Trunc") {
    preds_storage_glm[[method_label]] <- train_glm_fast(df, predictors, outcome_col)
  }
}

# ============================================================
# 3. 统计计算 (Agreement B=100 & DeLong)
# ============================================================

# --- 辅助函数：Agreement (B=100) ---
calc_agree_b100 <- function(p1, p2) {
  y1 <- ifelse(p1 >= 0.5, 1, 0); y2 <- ifelse(p2 >= 0.5, 1, 0)
  obs <- mean(y1 == y2)
  # 向量化 Bootstrap，极速
  n <- length(y1)
  idx_mat <- matrix(sample(1:n, n * 100, replace = TRUE), nrow = n, ncol = 100)
  boot_vals <- colMeans(matrix((y1 == y2)[idx_mat], nrow = n))
  return(c(obs, quantile(boot_vals, 0.025), quantile(boot_vals, 0.975)))
}

# --- 批量 Agreement ---
compute_matrix_b100 <- function(p_list) {
  methods <- names(p_list); n <- length(methods)
  m_mean <- matrix(NA, n, n, dimnames=list(methods, methods))
  m_low  <- matrix(NA, n, n, dimnames=list(methods, methods))
  m_up   <- matrix(NA, n, n, dimnames=list(methods, methods))
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i <= j) {
        sub1 <- p_list[[methods[i]]]; sub2 <- p_list[[methods[j]]]
        c_ids <- intersect(sub1$row_id, sub2$row_id)
        if(length(c_ids)>0) {
          s1 <- sub1[sub1$row_id %in% c_ids, ]; s2 <- sub2[sub2$row_id %in% c_ids, ]
          res <- calc_agree_b100(s1$prob, s2$prob)
          
          m_mean[i,j] <- m_mean[j,i] <- res[1]
          m_low[i,j]  <- m_low[j,i]  <- res[2]
          m_up[i,j]   <- m_up[j,i]   <- res[3]
        }
      }
    }
  }
  return(list(Mean=m_mean, Lower=m_low, Upper=m_up))
}

# --- 批量 DeLong ---
compute_delong <- function(p_list) {
  methods <- names(p_list); n <- length(methods)
  m_p <- matrix(NA, n, n, dimnames=list(methods, methods))
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i < j) {
        sub1 <- p_list[[methods[i]]]; sub2 <- p_list[[methods[j]]]
        c_ids <- intersect(sub1$row_id, sub2$row_id)
        if(length(c_ids)>0) {
          s1 <- sub1[sub1$row_id %in% c_ids, ]; s2 <- sub2[sub2$row_id %in% c_ids, ]
          # pROC 包自带 delong 很快
          t_res <- roc.test(roc(s1$truth, s1$prob, quiet=T), roc(s2$truth, s2$prob, quiet=T), method="delong")
          m_p[i,j] <- m_p[j,i] <- t_res$p.value
        }
      } else if (i==j) { m_p[i,j] <- 1 }
    }
  }
  return(m_p)
}

message(">>> 2. Calculating Statistics (Agreement B=100 & DeLong)...")

# 计算 XGBoost
res_agree_xgb <- compute_matrix_b100(preds_storage_xgb)
res_delong_xgb <- compute_delong(preds_storage_xgb)

# 计算 Elastic Net
res_agree_glm <- compute_matrix_b100(preds_storage_glm)
res_delong_glm <- compute_delong(preds_storage_glm)

# ============================================================
# 4. 保存结果
# ============================================================
res_dir <- file.path(work_dir, "Evaluation_Fast_Results")
if(!dir.exists(res_dir)) dir.create(res_dir)

# 保存 XGBoost
write.csv(res_agree_xgb$Mean, file.path(res_dir, "Agreement_Mean_XGBoost.csv"))
write.csv(res_agree_xgb$Lower, file.path(res_dir, "Agreement_Lower_XGBoost.csv"))
write.csv(res_agree_xgb$Upper, file.path(res_dir, "Agreement_Upper_XGBoost.csv"))
write.csv(res_delong_xgb, file.path(res_dir, "DeLong_Pvalues_XGBoost.csv"))

# 保存 Elastic Net
write.csv(res_agree_glm$Mean, file.path(res_dir, "Agreement_Mean_ElasticNet.csv"))
write.csv(res_agree_glm$Lower, file.path(res_dir, "Agreement_Lower_ElasticNet.csv"))
write.csv(res_agree_glm$Upper, file.path(res_dir, "Agreement_Upper_ElasticNet.csv"))
write.csv(res_delong_glm, file.path(res_dir, "DeLong_Pvalues_ElasticNet.csv"))

message("✅ 搞定！结果已保存到 Evaluation_Fast_Results 文件夹。")
# ============================================================
# XGBoost Feature Selection - CPI 方法 v5
# B=100 bootstrap, n_perm=100, 70/30 split, 并行4核
# 入选标准：CI 下界 > 0
# ============================================================

library(dplyr)
library(readr)
library(xgboost)
library(parallel)
library(doParallel)
library(foreach)

# ============================================================
# 配置
# ============================================================

input_dirs <- list(
  "10%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_10",
  "20%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_20",
  "30%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_30",
  "40%" = "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/trunc_40"
)

output_dir <- "C:/Users/20210607/OneDrive - TU Eindhoven/Desktop/Task 2 2025-2026/CPI_Results"

predictors <- c(
  "U3CC07", "U3CD07", "Height", "BMI", "SFLT1_V1",
  "fbHCG_V1", "V2BA01_KG", "Age_at_V1", "ADAM12_V1", "AFP_V1",
  "PLGF_V1", "INHIBINA_V1", "U2BB02", "PAPPA_V1", "SFLT1_V2"
)

outcome_col <- "outcome"

file_pattern_map <- list(
  "simulated_dataset_final" = "Complete_Data",
  "imputed_RandomBase"      = "Random",
  "dataset_trunc"           = "Raw_Trunc",
  "imputed_LOD_Half"        = "LOD_Half_Sub",
  "imputed_LOD_"            = "LOD_Fill",
  "imputed_Gaussian"        = "Gauss_Draw",
  "imputed_Lognormal"       = "Lognorm_Draw",
  "mice_posthoc"            = "MICE_Posthoc",
  "mice_truncnorm"          = "MICE_Trunc",
  "tobit_averaged"          = "Tobit_Avg",
  "imputed_DB_Mean"         = "DB_Mean",
  "imputed_DB_Draw"         = "DB_Draw"
)

N_CORES <- 4
B_BOOT  <- 100
N_PERM  <- 100
TRAIN_RATIO <- 0.7

# ============================================================
# 辅助函数
# ============================================================

logloss <- function(y, p) {
  p <- pmax(pmin(p, 1 - 1e-10), 1e-10)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

conditional_permute <- function(X, var_idx) {
  xi      <- X[, var_idx]
  X_other <- X[, -var_idx, drop = FALSE]
  complete_idx <- complete.cases(cbind(xi, X_other))
  
  if (sum(complete_idx) < 10) {
    xi_perm <- X[, var_idx]
    xi_perm[!is.na(xi_perm)] <- sample(xi_perm[!is.na(xi_perm)])
    X_perm            <- X
    X_perm[, var_idx] <- xi_perm
    return(X_perm)
  }
  
  fit <- lm(xi ~ X_other, subset = complete_idx, na.action = na.exclude)
  fitted_vals <- rep(NA, length(xi))
  resid_vals  <- rep(NA, length(xi))
  fitted_vals[complete_idx] <- fitted(fit)
  resid_vals[complete_idx]  <- residuals(fit)
  non_na_idx             <- which(!is.na(resid_vals))
  resid_vals[non_na_idx] <- sample(resid_vals[non_na_idx])
  xi_perm           <- fitted_vals + resid_vals
  X_perm            <- X
  X_perm[, var_idx] <- xi_perm
  return(X_perm)
}

# ============================================================
# 单次 bootstrap 函数（并行的最小单元）
# ============================================================

one_bootstrap <- function(b, X, y, predictors, train_ratio, n_perm) {
  
  set.seed(b)
  n         <- nrow(X)
  train_idx <- sample(1:n, floor(n * train_ratio), replace = FALSE)
  test_idx  <- setdiff(1:n, train_idx)
  
  X_train <- X[train_idx, ]; y_train <- y[train_idx]
  X_test  <- X[test_idx,  ]; y_test  <- y[test_idx]
  
  # 训练
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  params <- list(objective = "binary:logistic", eval_metric = "logloss",
                 eta = 0.05, max_depth = 4)
  model  <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
  
  # ── Feature Importance (Gain) ──
  imp_df   <- xgb.importance(model = model)
  gain_vec <- setNames(rep(0, length(predictors)), predictors)
  if (!is.null(imp_df) && nrow(imp_df) > 0) {
    for (feat in imp_df$Feature) {
      if (feat %in% predictors)
        gain_vec[feat] <- imp_df$Gain[imp_df$Feature == feat]
    }
  }
  
  # ── CPI: baseline loss 在测试集上 ──
  dtest     <- xgb.DMatrix(data = X_test, label = y_test)
  loss_base <- logloss(y_test, predict(model, dtest))
  
  # 对每个变量做条件置换
  cpi_deltas <- list()
  for (i in seq_along(predictors)) {
    perm_losses <- numeric(n_perm)
    for (k in 1:n_perm) {
      X_perm         <- conditional_permute(X_test, i)
      dperm          <- xgb.DMatrix(data = X_perm, label = y_test)
      perm_losses[k] <- logloss(y_test, predict(model, dperm))
    }
    cpi_deltas[[i]] <- perm_losses - loss_base
  }
  
  return(list(gain = gain_vec, cpi = cpi_deltas))
}

# ============================================================
# 核心函数：对单个数据集跑 CPI
# ============================================================

compute_xgb_cpi_boot <- function(df, predictors, outcome_col,
                                 B = 100, n_perm = 100,
                                 train_ratio = 0.7, n_cores = 4) {
  
  X <- as.matrix(df[, predictors])
  y <- as.numeric(df[[outcome_col]])
  
  # 并行跑 B 次 bootstrap
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  boot_results <- foreach(
    b = 1:B,
    .packages  = c("xgboost"),
    .export    = c("one_bootstrap", "conditional_permute", "logloss"),
    .errorhandling = "pass"
  ) %dopar% {
    one_bootstrap(b, X, y, predictors, train_ratio, n_perm)
  }
  
  stopCluster(cl)
  
  # 汇总：Importance + Selection 一并输出
  results <- list()
  
  for (i in seq_along(predictors)) {
    
    # ── Feature Importance (Gain) ──
    # 100次bootstrap每次的Gain值
    gain_vals  <- sapply(boot_results, function(br) br$gain[predictors[i]])
    gain_mean  <- mean(gain_vals)
    gain_lower <- quantile(gain_vals, 0.025)
    gain_upper <- quantile(gain_vals, 0.975)
    
    # ── CPI Selection ──
    # 每次bootstrap的n_perm个delta
    per_boot_deltas <- lapply(boot_results, function(br) br$cpi[[i]])
    
    # 每次bootstrap单独判断入选（该次CI下界>0）
    per_boot_selected <- sapply(per_boot_deltas, function(deltas) {
      as.integer(quantile(deltas, 0.025) > 0)
    })
    selection_freq <- mean(per_boot_selected)
    
    # 合并所有delta算总体CPI估计
    all_deltas <- unlist(per_boot_deltas)
    cpi_score  <- mean(all_deltas)
    cpi_lower  <- quantile(all_deltas, 0.025)
    cpi_upper  <- quantile(all_deltas, 0.975)
    selected   <- as.integer(cpi_lower > 0)
    
    results[[i]] <- data.frame(
      Feature        = predictors[i],
      # Importance
      Gain_Mean      = gain_mean,
      Gain_Lower     = as.numeric(gain_lower),
      Gain_Upper     = as.numeric(gain_upper),
      # Selection
      CPI_Score      = cpi_score,
      CPI_Lower      = as.numeric(cpi_lower),
      CPI_Upper      = as.numeric(cpi_upper),
      Selected       = selected,        # 合并所有delta的最终判断
      Selection_Freq = selection_freq   # 100次bootstrap里入选的频率
    )
  }
  
  bind_rows(results) %>% arrange(desc(CPI_Score))
}

# ============================================================
# 主循环
# ============================================================

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_results <- list()

for (trunc_label in names(input_dirs)) {
  
  work_dir  <- input_dirs[[trunc_label]]
  all_files <- list.files(work_dir, pattern = "*.csv", full.names = TRUE)
  
  message(paste0("\n========== Truncation: ", trunc_label,
                 " | Found ", length(all_files), " files =========="))
  
  for (file_path in all_files) {
    
    file_name    <- basename(file_path)
    method_label <- NA
    for (pattern in names(file_pattern_map)) {
      if (grepl(pattern, file_name, fixed = TRUE)) {
        method_label <- file_pattern_map[[pattern]]
        break
      }
    }
    if (is.na(method_label)) next
    
    message(paste("\n--- Method:", method_label, "| Truncation:", trunc_label, "---"))
    t_start <- proc.time()
    
    df  <- read_csv(file_path, show_col_types = FALSE)
    cpi <- compute_xgb_cpi_boot(df, predictors, outcome_col,
                                B = B_BOOT, n_perm = N_PERM,
                                train_ratio = TRAIN_RATIO, n_cores = N_CORES)
    
    cpi$Method     <- method_label
    cpi$Truncation <- trunc_label
    
    t_elapsed <- round((proc.time() - t_start)["elapsed"] / 60, 1)
    message(paste("   完成，耗时", t_elapsed, "分钟"))
    message(paste("   入选变量:", paste(cpi$Feature[cpi$Selected == 1], collapse = ", ")))
    
    key                <- paste0(trunc_label, "_", method_label)
    all_results[[key]] <- cpi
  }
}

# ============================================================
# 保存结果
# ============================================================

final_table <- bind_rows(all_results)

write.csv(
  final_table,
  file.path(output_dir, "XGBoost_CPI_Selection_AllTrunc.csv"),
  row.names = FALSE
)

message("\n✅ 完成！结果保存至: ", file.path(output_dir, "XGBoost_CPI_Selection_AllTrunc.csv"))
message("   输出列说明：")
message("   Gain_Mean      : 100次bootstrap的Gain均值（feature importance）")
message("   Gain_Lower/Upper: Gain的95% CI")
message("   CPI_Score      : 置换后误差升高量均值（越大越重要）")
message("   CPI_Lower/Upper: CPI的95% CI")
message("   Selected       : 1=入选（合并所有delta后CI下界>0）")
message("   Selection_Freq : 100次bootstrap里入选的频率（稳定性指标）")
message("   Truncation     : 截断比例")

