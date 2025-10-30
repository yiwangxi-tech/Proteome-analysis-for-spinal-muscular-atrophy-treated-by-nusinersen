library(openxlsx)
library(tidyverse)
library(gtsummary)
# 连续变量及分类变量的Wlicoxon秩和检验及Fisher精确检验
Pat_info <- read.xlsx("Patient_info.xlsx", colNames = TRUE)
Pat_info <- Pat_info %>%
  mutate(
    Group = factor(Type),
    Gender = factor(Gender),
    Classification = factor(Classification),
    SMN2.Copy.Number = factor(SMN2.Copy.Number) 
  ) %>%
  select(
    Group, Age, Gender, Classification, SMN2.Copy.Number
  )

continuous_vars <- c("Age", "SMN2.Copy.Number")
categorical_vars <- c("Gender", "Classification")
test_results <- list()

for (var in continuous_vars) {
  Pat_info[[var]] <- as.numeric(Pat_info[[var]])
  test <- wilcox.test(Pat_info[[var]] ~ Pat_info$Group, data = Pat_info)
  test_results[[var]] <- test
  cat(paste("变量:", var, "\n"))
  cat(paste("  P 值:", sprintf("%.4f", test$p.value), "\n"))
  print(test)
  cat("\n")
}

for (var in categorical_vars) {
  contingency_table <- table(Pat_info[[var]], Pat_info$Group)
  chi_sq_test <- tryCatch({
    chisq.test(contingency_table)
  }, warning = function(w) {
    if (grepl("Chi-squared approximation may be incorrect", w$message) || any(chisq.test(contingency_table)$expected < 5)) {
      cat(paste("  注意:", var, "的卡方检验期望频数过低，自动使用 Fisher 精确检验。\n"))
      return(fisher.test(contingency_table))
    } else {
      return(chisq.test(contingency_table))
    }
  })
  test_results[[var]] <- chi_sq_test
  cat(paste("变量:", var, "\n"))
  if ("Fisher's Exact Test" %in% chi_sq_test$method) {
      cat(paste("  检验类型: Fisher's Exact Test\n"))
      cat(paste("  P 值:", sprintf("%.4f", chi_sq_test$p.value), "\n"))
  } else {
      cat(paste("  检验类型: Pearson's Chi-squared Test\n"))
      cat(paste("  P 值:", sprintf("%.4f", chi_sq_test$p.value), "\n"))
  }
  print(chi_sq_test)
  cat("\n")
}
# 评分变化的Wlicoxon秩和检验
score_vars <- c("HFMSE_Pre", "RULM_Pre", "HFMSE_1y", "RULM_1y", 
                "HFMSE_2y", "RULM_2y", "HFMSE_2.5y", "RULM_2.5y")
Pat_info <- Pat_info %>%
  mutate(
    Group = factor(Type), 
    across(all_of(score_vars), as.numeric) 
  ) %>%
  mutate(
    # HFMSE 变化量
    HFMSE_Delta_1y   = HFMSE_1y - HFMSE_Pre,
    HFMSE_Delta_2y   = HFMSE_2y - HFMSE_Pre,
    HFMSE_Delta_2.5y = HFMSE_2.5y - HFMSE_Pre,
    
    # RULM 变化量
    RULM_Delta_1y   = RULM_1y - RULM_Pre,
    RULM_Delta_2y   = RULM_2y - RULM_Pre,
    RULM_Delta_2.5y = RULM_2.5y - RULM_Pre
  )
baseline_vars <- c("HFMSE_Pre", "RULM_Pre")
delta_vars <- c(
  "HFMSE_Delta_1y", "HFMSE_Delta_2y", "HFMSE_Delta_2.5y",
  "RULM_Delta_1y", "RULM_Delta_2y", "RULM_Delta_2.5y"
)
groups <- c("Respond", "Non-respond")
all_continuous_vars <- c(baseline_vars, delta_vars)
for (var in all_continuous_vars) {
  if (var %in% baseline_vars) {
      cat(paste("【指标 (基线):", var, "】\n"))
  } else {
      cat(paste("【指标 (变化量):", var, "】\n"))
  }
  for (group in groups) {
    data_subset <- Pat_info[[var]][Pat_info$Type == group]
    quantiles <- quantile(data_subset, na.rm = TRUE)
    cat(paste("  分组:", group, "\n"))
    print(quantiles)
    median_val <- quantiles["50%"]
    iqr_val <- paste0("(", quantiles["25%"], ", ", quantiles["75%"], ")")
    cat(paste("  **中位数 (IQR):**", median_val, iqr_val, "\n\n"))
  }
  test_result <- wilcox.test(Pat_info[[var]] ~ Pat_info$Group, data = Pat_info, exact = FALSE)
  cat("--- 组间差异检验 (Wilcoxon 秩和检验) ---\n")
  cat(paste("  P 值:", sprintf("%.4f", test_result$p.value), "\n"))
  print(test_result)
  cat("========================================================================\n\n")
}
# 不同时间点箱线图绘制
time_points <- c("Pre", "1y", "2y", "2.5y")
df_long <- Pat_info %>%
  select(
    Patient.ID, Group, 
    HFMSE_Pre, RULM_Pre, HFMSE_1y, RULM_1y, HFMSE_2y, RULM_2y, HFMSE_2.5y, RULM_2.5y
  ) %>%
  pivot_longer(
    cols = -c(Patient.ID, Group),
    names_to = "Time_Score_Combined",
    values_to = "Score_Value"
  ) %>%
  separate(
    col = Time_Score_Combined,
    into = c("Score_Type", "Time_Point"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    Time_Point = factor(Time_Point, levels = time_points)
  )
plot_score_trajectories_final <- function(data) {
  score_labels <- c(
    "HFMSE" = "HFMSE Score", 
    "RULM" = "RULM Score"
  )
  ggplot(data, aes(x = Time_Point, y = Score_Value, fill = Group)) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.6,
      width = 0.7,
      position = position_dodge(0.8)
    ) +
    geom_point(
      aes(color = Group),
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
      size = 2,
      alpha = 0.7
    ) +
    facet_wrap(~ Score_Type, 
               scales = "free_y", 
               ncol = 2,
               labeller = as_labeller(score_labels), 
               strip.position = "left") + 
    labs(
      title = "评分随时间的变化和组间比较"
    ) +
    scale_fill_manual(values = c("Respond" = "#F48892", "Non-respond" = "#91CAE8")) +
    scale_color_manual(values = c("Respond" = "#F48892", "Non-respond" = "#91CAE8")) +
    theme_bw(base_size = 14, base_family = "Times New Roman") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"), 
      axis.title = element_text(face = "bold", family = "Times New Roman"),
      axis.text = element_text(family = "Times New Roman"),
      strip.background = element_blank(), 
      strip.placement = "outside", 
      strip.text = element_text(face = "bold", hjust = 0.5, family = "Times New Roman"),
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Times New Roman"),
      legend.position = "bottom"
    )
}
# 蛋白质组PCA降维分析
pro <- read.xlsx("蛋白定量列表.xlsx", colNames = TRUE)
row.names(pro) <- pro[, 1]
pro_data <- pro %>% select(-1)
all_samples <- colnames(pro_data)
# log2转化+数据处理
data_matrix_log2 <- pro_data %>% 
  mutate(across(everything(), as.numeric)) %>% 
  mutate(across(everything(), log2)) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  as.matrix()
sample_id <- as.numeric(gsub("[AB]", "", all_samples))
sample_metadata <- data.frame(
  Sample = all_samples,
  ID = sample_id,
  Treatment_Group = ifelse(sample_id %in% c(5, 6, 9), "Respond", "Non-respond"),
  Time_Point = ifelse(grepl("A$", all_samples), "Pre-Treatment", "Post-Treatment")
) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group),
    Time_Point = factor(Time_Point, levels = c("Pre-Treatment", "Post-Treatment"))
  )
pca_result <- prcomp(t(data_matrix_log2), center = TRUE, scale. = TRUE)
pca_df <- data.frame(
     PC1 = pca_result$x[, 1],
     PC2 = pca_result$x[, 2],
     Sample = rownames(pca_result$x) 
) %>%
  left_join(sample_metadata, by = "Sample")
var_explained <- summary(pca_result)$importance[2, 1:2] * 100
print(var_explained)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment_Group, shape = Time_Point)) +
  geom_text_repel(
    aes(label = ID, color = Treatment_Group), 
    size = 3,                                
    segment.alpha = 0.5,                     
    family = "Times New Roman",              
    show.legend = FALSE                      
  ) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Treatment_Group, color = Treatment_Group), geom = "polygon", alpha = 0.1, level = 0.95, show.legend = FALSE) +
  geom_line(aes(group = ID, color = Treatment_Group), alpha = 0.4) +
  scale_color_manual(
        values = c("Respond" = "#F48892", "Non-respond" = "#91CAE8")
  ) +
  scale_shape_manual(
        values = c("Pre-Treatment" = 16, "Post-Treatment" = 17) 
  ) +
  labs(
    title = "Proteome PCA",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "Treatment response",
    shape = "Pre- and Post-treatment"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"),
    legend.title = element_text(size = 14, family = "Times New Roman"),
    legend.text = element_text(size = 12, family = "Times New Roman"),
    panel.grid.major = element_line(linetype = "dotted", color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12, family = "Times New Roman"), 
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman") 
  )

#蛋白质组差异分析(10A/B样本)
library(openxlsx)
library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)

Pro <- read.xlsx("蛋白定量列表.xlsx", colNames = TRUE)
protein_matrix <- as.matrix(Pro[, -1])
rownames(protein_matrix) <- Pro$Protein
log2_data <- log2(protein_matrix)

respond_before <- c('5A', '6A', '9A')
respond_after <- c('5B', '6B', '9B')
non_respond_before <- c('1A', '2A', '3A', '4A', '7A', '8A', '10A')
non_respond_after <- c('1B', '2B', '3B', '4B', '7B', '8B', '10B')

output_dir <- "D:/博士/参考文献/蛋白质组学/许婷婷-协和-20个脑脊液ddia-241029/R_output1013"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

perform_de_analysis <- function() {
  results_list <- list()
  
  # Analysis 1: Respond vs Non-respond (Before treatment)
  before_samples <- c(respond_before, non_respond_before)
  before_data <- log2_data[, before_samples]
  before_group <- factor(c(rep("Respond", 3), rep("Non_respond", 7)))
  before_design <- model.matrix(~0 + before_group)
  colnames(before_design) <- c("Non_respond", "Respond")
  before_fit <- lmFit(before_data, before_design)
  before_contrast <- makeContrasts(Respond - Non_respond, levels = before_design)
  before_fit2 <- contrasts.fit(before_fit, before_contrast)
  before_fit2 <- eBayes(before_fit2)
  results_list$before_results <- topTable(before_fit2, n = Inf)
  results_list$before_results$Protein <- rownames(results_list$before_results)
  
  # Analysis 2: Respond vs Non-respond (After treatment)
  after_samples <- c(respond_after, non_respond_after)
  after_data <- log2_data[, after_samples]
  after_group <- factor(c(rep("Respond", 3), rep("Non_respond", 7)))
  after_design <- model.matrix(~0 + after_group)
  colnames(after_design) <- c("Non_respond", "Respond")
  after_fit <- lmFit(after_data, after_design)
  after_contrast <- makeContrasts(Respond - Non_respond, levels = after_design)
  after_fit2 <- contrasts.fit(after_fit, after_contrast)
  after_fit2 <- eBayes(after_fit2)
  results_list$after_results <- topTable(after_fit2, n = Inf)
  results_list$after_results$Protein <- rownames(results_list$after_results)
  
  # Analysis 3: Respond group (Before vs After treatment)
  respond_samples <- c(respond_before, respond_after)
  respond_data <- log2_data[, respond_samples]
  respond_time_group <- factor(c(rep("Before", 3), rep("After", 3)))
  respond_design <- model.matrix(~0 + respond_time_group)
  colnames(respond_design) <- c("After", "Before")
  respond_fit <- lmFit(respond_data, respond_design)
  respond_contrast <- makeContrasts(After - Before, levels = respond_design)
  respond_fit2 <- contrasts.fit(respond_fit, respond_contrast)
  respond_fit2 <- eBayes(respond_fit2)
  results_list$respond_ba_results <- topTable(respond_fit2, n = Inf)
  results_list$respond_ba_results$Protein <- rownames(results_list$respond_ba_results)
  
  # Analysis 4: Non-respond group (Before vs After treatment)
  nonrespond_samples <- c(non_respond_before, non_respond_after)
  nonrespond_data <- log2_data[, nonrespond_samples]
  nonrespond_time_group <- factor(c(rep("Before", 7), rep("After", 7)))
  nonrespond_design <- model.matrix(~0 + nonrespond_time_group)
  colnames(nonrespond_design) <- c("After", "Before")
  nonrespond_fit <- lmFit(nonrespond_data, nonrespond_design)
  nonrespond_contrast <- makeContrasts(After - Before, levels = nonrespond_design)
  nonrespond_fit2 <- contrasts.fit(nonrespond_fit, nonrespond_contrast)
  nonrespond_fit2 <- eBayes(nonrespond_fit2)
  results_list$nonrespond_ba_results <- topTable(nonrespond_fit2, n = Inf)
  results_list$nonrespond_ba_results$Protein <- rownames(results_list$nonrespond_ba_results)
  
  return(results_list)
}

all_results <- perform_de_analysis()

create_volcano_plot <- function(results, title, filename) {
    results$significance <- "Not Significant"
    results$significance[results$P.Value < 0.05 & abs(results$logFC) > 0.5] <- "Significant"
    results$significance[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5] <- "Highly Significant"
    
    results$logP <- -log10(results$P.Value)
    
    colors <- c("#2472A3", "#7EA4D1", "#82AD7F", "#F8F2A4", "#DCA96A", "#C1565E", "#B02425")
    
    p <- ggplot(results, aes(x = logFC, y = logP)) +
        geom_point(
            data = subset(results, significance == "Not Significant"),
            aes(color = logFC, size = logP),
            alpha = 0.6,
            shape = 16
        ) +
        geom_point(
            data = subset(results, significance != "Not Significant"),
            aes(color = logFC, size = logP),
            shape = 16
        ) +
        scale_color_gradientn(colors = colors, name = expression(log[2]*FoldChange)) +
        scale_size_continuous(range = c(1, 4), name = expression(-log[10](P[value]))) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        labs(
            title = title,
            x = expression(log[2]*FoldChange),
            y = expression(-log[10](P[value]))
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.title = element_text(size = 14),   
            axis.text = element_text(size = 12)     
        )

    top_proteins <- head(results[results$P.Value < 0.05 & abs(results$logFC) > 0.5, ], 5)
    if (nrow(top_proteins) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_proteins,
            aes(label = Protein),
            size = 3,
            box.padding = 0.3,
            max.overlaps = 10
        )
    }
    
    ggsave(filename, p, width = 10, height = 8, dpi = 300)
    print(paste("Volcano plot saved:", filename))
    return(p)
}

create_volcano_plot(all_results$before_results, "Respond vs Non-respond (Before Treatment)", 
                   file.path(output_dir, "volcano_before_respond_vs_nonrespond.png"))
create_volcano_plot(all_results$after_results, "Respond vs Non-respond (After Treatment)", 
                   file.path(output_dir, "volcano_after_respond_vs_nonrespond.png"))
create_volcano_plot(all_results$respond_ba_results, "Respond Group: After vs Before Treatment", 
                   file.path(output_dir, "volcano_respond_after_vs_before.png"))
create_volcano_plot(all_results$nonrespond_ba_results, "Non-respond Group: After vs Before Treatment", 
                   file.path(output_dir, "volcano_nonrespond_after_vs_before.png"))

SMA1_before <- c('8A')
SMA1_after <- c('8B')
SMA2_before <- c('4A', '5A', '10A')
SMA2_after <- c('4B', '5B', '10B')
SMA3_before <- c('1A', '2A', '3A', '6A', '7A', '9A')
SMA3_after <- c('1B', '2B', '3B', '6B', '7B', '9B')

perform_de_analysis <- function() {
     results_list <- list()
     
     # Analysis 1: SMA-1 vs SMA-2 (Before treatment)
     before_samples_1 <- c(SMA1_before, SMA2_before)
     before_data_1 <- log2_data[, before_samples_1]
     before_group_1 <- factor(c(rep("SMA1", 1), rep("SMA2", 3)))
     before_design_1 <- model.matrix(~0 + before_group_1)
     colnames(before_design_1) <- c("SMA1", "SMA2")
     before_fit_1 <- lmFit(before_data_1, before_design_1)
     before_contrast_1 <- makeContrasts(SMA1 - SMA2, levels = before_design_1)
     before_fit2_1 <- contrasts.fit(before_fit_1, before_contrast_1)
     before_fit2_1 <- eBayes(before_fit2_1)
     results_list$before_results_1 <- topTable(before_fit2_1, n = Inf)
     results_list$before_results$Protein_1 <- rownames(results_list$before_results_1)
     
     # Analysis 2: SMA-1 vs SMA-3 (Before treatment)
     before_samples_2 <- c(SMA1_before, SMA3_before)
     before_data_2 <- log2_data[, before_samples_2]
     before_group_2 <- factor(c(rep("SMA1", 1), rep("SMA3", 6)))
     before_design_2 <- model.matrix(~0 + before_group_2)
     colnames(before_design_2) <- c("SMA1", "SMA3")
     before_fit_2 <- lmFit(before_data_2, before_design_2)
     before_contrast_2 <- makeContrasts(SMA1 - SMA3, levels = before_design_2)
     before_fit2_2 <- contrasts.fit(before_fit_2, before_contrast_2)
     before_fit2_2 <- eBayes(before_fit2_2)
     results_list$before_results_2 <- topTable(before_fit2_2, n = Inf)
     results_list$before_results$Protein_2 <- rownames(results_list$before_results_2)
     
     # Analysis 3: SMA-2 vs SMA-3 (Before treatment)
     before_samples_3 <- c(SMA2_before, SMA3_before)
     before_data_3 <- log2_data[, before_samples_3]
     before_group_3 <- factor(c(rep("SMA2", 3), rep("SMA3", 6)))
     before_design_3 <- model.matrix(~0 + before_group_3)
     colnames(before_design_3) <- c("SMA2", "SMA3")
     before_fit_3 <- lmFit(before_data_3, before_design_3)
     before_contrast_3 <- makeContrasts(SMA2 - SMA3, levels = before_design_3)
     before_fit2_3 <- contrasts.fit(before_fit_3, before_contrast_3)
     before_fit2_3 <- eBayes(before_fit2_3)
     results_list$before_results_3 <- topTable(before_fit2_3, n = Inf)
     results_list$before_results$Protein_3 <- rownames(results_list$before_results_3)
     
     # Analysis 4: SMA-1 vs SMA-2 (After treatment)
     after_samples_1 <- c(SMA1_after, SMA2_after)
     after_data_1 <- log2_data[, after_samples_1]
     after_group_1 <- factor(c(rep("SMA1", 1), rep("SMA2", 3)))
     after_design_1 <- model.matrix(~0 + after_group_1)
     colnames(after_design_1) <- c("SMA1", "SMA2")
     after_fit_1 <- lmFit(after_data_1, after_design_1)
     after_contrast_1 <- makeContrasts(SMA1 - SMA2, levels = after_design_1)
     after_fit2_1 <- contrasts.fit(after_fit_1, after_contrast_1)
     after_fit2_1 <- eBayes(after_fit2_1)
     results_list$after_results_1 <- topTable(after_fit2_1, n = Inf)
     results_list$after_results$Protein_1 <- rownames(results_list$after_results_1)
     
     # Analysis 5: SMA-1 vs SMA-3 (After treatment)
     after_samples_2 <- c(SMA1_after, SMA3_after)
     after_data_2 <- log2_data[, after_samples_2]
     after_group_2 <- factor(c(rep("SMA1", 1), rep("SMA3", 6)))
     after_design_2 <- model.matrix(~0 + after_group_2)
     colnames(after_design_2) <- c("SMA1", "SMA3")
     after_fit_2 <- lmFit(after_data_5, after_design_2)
     after_contrast_2 <- makeContrasts(SMA1 - SMA3, levels = after_design_2)
     after_fit2_2 <- contrasts.fit(before_fit_2, after_contrast_2)
     after_fit2_2 <- eBayes(after_fit2_2)
     results_list$after_results_2 <- topTable(after_fit2_2, n = Inf)
     results_list$after_results$Protein_2 <- rownames(results_list$after_results_2)
     
     # Analysis 6: SMA-2 vs SMA-3 (After treatment)
     after_samples_3 <- c(SMA2_after, SMA3_after)
     after_data_3 <- log2_data[, after_samples_3]
     after_group_3 <- factor(c(rep("SMA2", 3), rep("SMA3", 6)))
     after_design_3 <- model.matrix(~0 + after_group_3)
     colnames(after_design_3) <- c("SMA2", "SMA3")
     after_fit_3 <- lmFit(after_data_3, after_design_3)
     after_contrast_3 <- makeContrasts(SMA2 - SMA3, levels = after_design_3)
     after_fit2_3 <- contrasts.fit(after_fit_3, after_contrast_3)
     after_fit2_3 <- eBayes(after_fit2_3)
     results_list$after_results_3 <- topTable(after_fit2_3, n = Inf)
     results_list$after_results$Protein_3 <- rownames(results_list$after_results_3)
     
     # Analysis 7: Before vs After (SMA-1)
     sma1_samples <- c(SMA1_before, SMA1_after)
     sma1_data <- log2_data[, sma1_samples]
     sma1_group <- factor(c(rep("Before", 1), rep("After", 1)))
     sma1_design <- model.matrix(~0 + sma1_group)
     colnames(sma1_design) <- c("Before", " After")
     sma1_fit <- lmFit(sma1_data_7, sma1_design)
     sma1_contrast <- makeContrasts(After - Before, levels = sma1_design)
     sma1_fit2 <- contrasts.fit(sma1_fit, sma1_contrast)
     sma1_fit2 <- eBayes(sma1_fit2)
     results_list$sma1_results <- topTable(sma1_fit2, n = Inf)
     results_list$sma1_results$Protein <- rownames(results_list$sma1_results)
     
     # Analysis 8: Before vs After (SMA-2)
     sma2_samples <- c(SMA2_before, SMA2_after)
     sma2_data <- log2_data[, sma2_samples]
     sma2_group <- factor(c(rep("Before", 3), rep("After", 3)))
     sma2_design <- model.matrix(~0 + sma2_group)
     colnames(before_design) <- c("Before", " After")
     sma2_fit <- lmFit(sma2_data, sma2_design)
     sma2_contrast <- makeContrasts(After - Before, levels = sma2_design)
     sma2_fit2 <- contrasts.fit(sma2_fit, sma2_contrast)
     sma2_fit2 <- eBayes(sma2_fit2)
     results_list$sma2_results <- topTable(sma2_fit2, n = Inf)
     results_list$sma2_results$Protein <- rownames(results_list$sma2_results_)
     
     # Analysis 9: Before vs After (SMA-3)
     sma3_samples <- c(SMA3_before, SMA3_after)
     sma3_data <- log2_data[, sma3_samples]
     sma3_group <- factor(c(rep("Before", 6), rep("After", 6)))
     sma3_design <- model.matrix(~0 + sma3_group)
     colnames(sma3_design) <- c("Before", " After")
     sma3_fit <- lmFit(sma3_data, sma3_design)
     sma3_contrast <- makeContrasts(After - Before, levels = sma3_design)
     sma3_fit2 <- contrasts.fit(sma3_fit, sma3_contrast)
     sma3_fit2 <- eBayes(sma3_fit2)
     results_list$sma3_results <- topTable(sma3_fit2, n = Inf)
     results_list$sma3_results$Protein <- rownames(results_list$sma3_results)
     
     return(results_list)
}

