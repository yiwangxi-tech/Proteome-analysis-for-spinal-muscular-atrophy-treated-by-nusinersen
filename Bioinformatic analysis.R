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

# 相关性热图绘制
library(pheatmap) 
library(tibble)
library(openxlsx)
library(dplyr)
Pro <- read.xlsx("蛋白定量列表.xlsx", colNames = TRUE)
protein_matrix <- as.matrix(Pro[, -1])
rownames(protein_matrix) <- Pro$Protein
log2_data <- log2(protein_matrix)
zscore_data <- t(scale(t(log2_data)))
effective_patients <- c(5, 6, 9)
sample_names <- colnames(zscore_data)
patients <- as.numeric(gsub("[AB]", "", sample_names))
conditions <- gsub("[0-9]", "", sample_names)
effectiveness <- ifelse(patients %in% effective_patients, "Respond", "Non-respond")

sample_info <- data.frame(
  Sample = sample_names,
  Patient = patients,
  Condition = conditions,
  Effectiveness = effectiveness
)
effective_idx <- which(effectiveness == "Respond")
ineffective_idx <- which(effectiveness == "Non-respond")

corr_all <- cor(zscore_data, method = "pearson")
corr_eff <- cor(zscore_data[, effective_idx], method = "pearson")  
corr_ineff <- cor(zscore_data[, ineffective_idx], method = "pearson")

mean_all <- mean(corr_all[upper.tri(corr_all)])
mean_eff <- mean(corr_eff[upper.tri(corr_eff)])
mean_ineff <- mean(corr_ineff[upper.tri(corr_ineff)])

col_annotation <- data.frame(
  row.names = sample_names,
  Effectiveness = factor(effectiveness),
  Condition = factor(conditions)
)

ann_colors <- list(
  Effectiveness = c("Respond" = "#F48892", "Non-respond" = "#91CAE8"),
  Condition = c("A" = "#84BA42", "B" = "#D4562E")
)
# Heatmap 1: All samples
pheatmap(corr_all,
         color = colorRampPalette(c("#E63946", "white", "#2E86AB"))(100),
         breaks = seq(-1, 1, length.out = 101),
         annotation_col = col_annotation,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         main = "Sample-to-Sample Correlations (All Samples)",
         fontsize = 14, 
         fontsize_row = 12, 
         fontsize_col = 12,
         cellwidth = 35,
         cellheight = 35,
         border_color = "grey90")

# Heatmap 2: Grouped by treatment effectiveness
sample_order <- c(sample_names[effective_idx], sample_names[ineffective_idx])
corr_grouped <- corr_all[sample_order, sample_order]
ann_grouped <- col_annotation[sample_order, , drop = FALSE]

pheatmap(corr_grouped,
         color = colorRampPalette(c("#E63946", "white", "#2E86AB"))(100),
         breaks = seq(-1, 1, length.out = 101),
         annotation_col = ann_grouped,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Sample Correlations: Respond vs Non-respond Groups",
         fontsize = 14,
         cellwidth = 30,
         cellheight = 30,
         border_color = "grey90")

# Protein analysis and biomarker identification
protein_vars <- apply(zscore_data, 1, var, na.rm = TRUE)
protein_mean_eff <- rowMeans(zscore_data[, effective_idx], na.rm = TRUE)
protein_mean_ineff <- rowMeans(zscore_data[, ineffective_idx], na.rm = TRUE)
protein_diff <- protein_mean_eff - protein_mean_ineff

# Create biomarker data frame
biomarkers <- data.frame(
  Protein = names(protein_diff),
  Mean_Effective = protein_mean_eff,
  Mean_Ineffective = protein_mean_ineff,
  Response_Difference = protein_diff,
  Abs_Response_Diff = abs(protein_diff),
  Variance = protein_vars
) %>% arrange(desc(Abs_Response_Diff))

top10 <- head(biomarkers, 10)
for(i in 1:nrow(top10)) {
  cat("  ", sprintf("%2d", i), ". ", top10$Protein[i], " (Response Diff: ", 
      sprintf("%6.2f", top10$Response_Difference[i]), ")\n", sep="")
}

# Protein-protein correlation analysis
n_top_proteins <- 150
top_variable_proteins <- names(sort(protein_vars, decreasing = TRUE)[1:n_top_proteins])
zscore_subset <- zscore_data[top_variable_proteins, ]

prot_corr_eff <- cor(t(zscore_subset[, effective_idx]), method = "pearson")
prot_corr_ineff <- cor(t(zscore_subset[, ineffective_idx]), method = "pearson")
corr_difference <- prot_corr_eff - prot_corr_ineff

# Create protein correlation heatmaps
pheatmap(prot_corr_eff,
         color = colorRampPalette(c("#E63946", "white", "#2E86AB"))(100),
         breaks = seq(-1, 1, length.out = 101),
         clustering_method = "average",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = paste("Protein Network Correlations: Respond Group (", n_top_proteins, " proteins)"),
         fontsize = 14,
         border_color = NA)

pheatmap(prot_corr_ineff,
         color = colorRampPalette(c("#E63946", "white", "#2E86AB"))(100),
         breaks = seq(-1, 1, length.out = 101),
         clustering_method = "average",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = paste("Protein Network Correlations: Non-respond Group (", n_top_proteins, " proteins)"),
         fontsize = 14,
         border_color = NA)

pheatmap(corr_difference,
         color = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
         breaks = seq(-2, 2, length.out = 101),
         clustering_method = "average",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Network Correlation Differences (Respond - Non-respond)",
         fontsize = 14,
         border_color = NA)

# WGCNA分析
setwd("D:/博士/参考文献/蛋白质组学/许婷婷-协和-20个脑脊液ddia-241029/R_output1013")
library(WGCNA)
library(pheatmap)
library(openxlsx)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)
## 全局WGCNA
Pro <- read.xlsx("蛋白定量列表.xlsx", colNames = TRUE)
expr_matrix <- as.matrix(Pro[,-1])
rownames(expr_matrix) <- Pro$Protein
datExpr <- as.data.frame(t(expr_matrix))
datExpr <- log2(datExpr)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
gene_vars <- apply(datExpr, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])
datExpr <- datExpr[, top_genes]
sample_names <- rownames(datExpr)
sample_info <- data.frame(
     sample_name = sample_names,
     sample_number = as.numeric(gsub("([0-9]+)[AB]", "\\1", sample_names)),
     treatment = gsub("[0-9]+([AB])", "\\1", sample_names)
 )
datTraits <- data.frame(
     Treatment_Response = ifelse(sample_info$sample_number %in% c(5,6,9), 1, 0),
     Treatment_Status = ifelse(sample_info$treatment == "A", 0, 1)
 )
rownames(datTraits) <- sample_names
sft <- pickSoftThreshold(datExpr, powerVector = c(1:15), verbose = 5)
softPower <- 6
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                        minModuleSize = 10, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, saveTOMs = FALSE,
                        verbose = 3)
moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
final_results <- data.frame(
     Module = rownames(moduleTraitCor),
     Module_Size = as.numeric(table(moduleColors)[gsub("ME", "", rownames(moduleTraitCor))]),
     Treatment_Response_Correlation = round(moduleTraitCor[, "Treatment_Response"], 4),
     Treatment_Response_Pvalue = round(moduleTraitPvalue[, "Treatment_Response"], 4),
     Treatment_Status_Correlation = round(moduleTraitCor[, "Treatment_Status"], 4),
     Treatment_Status_Pvalue = round(moduleTraitPvalue[, "Treatment_Status"], 4)
 )
module_assignment <- data.frame(
     Protein_ID = colnames(datExpr),
     Module_Color = moduleColors,
     stringsAsFactors = FALSE
 )
output_dir <- "D:/博士/参考文献/蛋白质组学/许婷婷-协和-20个脑脊液ddia-241029/R_output1013/WGCNA"
write.xlsx(final_results, paste0(output_dir, "/WGCNA_Final_Results.xlsx"), rowNames = FALSE)
write.xlsx(module_assignment, paste0(output_dir, "/Protein_Module_Assignment.xlsx"), rowNames = FALSE)

module <- "green"
moduleGenes <- (moduleColors == module)
trait_name <- "Treatment_status"  
trait <- datTraits[, trait_name]
MM <- cor(datExpr[, moduleGenes], MEs[, paste0("ME", module)], use="p")
MM_pvalue <- corPvalueStudent(MM, nSamples = nrow(datExpr))
GS <- cor(datExpr[, moduleGenes], trait, use="p")
GS_pvalue <- corPvalueStudent(GS, nSamples = nrow(datExpr))

## Respond及Non-respond WGNCA
library(dplyr)
library(stringr)
library(tibble)
Pro_matrix <- Pro
Pro_matrix <- Pro %>% column_to_rownames(var = colnames(Pro)[1])

respond_ids <- c("5", "6", "9")
non_respond_ids <- c("1", "2", "3", "4", "7", "8", "10")
filter_columns <- function(matrix, ids, time_suffix) {
     target_cols <- paste0(ids, time_suffix)
     selected_cols <- intersect(target_cols, colnames(matrix))
     if (length(selected_cols) == 0) {
         warning(paste("警告：未找到匹配", paste(ids, time_suffix, sep=""), "格式的列名。"))
         return(NULL)
     }
     return(matrix[, selected_cols, drop = FALSE])
 }
Resp_Before <- filter_columns(Pro_matrix, respond_ids, "A")
Resp_After <- filter_columns(Pro_matrix, respond_ids, "B")
NonResp_Before <- filter_columns(Pro_matrix, non_respond_ids, "A")
NonResp_After <- filter_columns(Pro_matrix, non_respond_ids, "B")
if (!is.null(Resp_Before) && !is.null(Resp_After)) {
     Respond_WGCNA_Matrix <- cbind(Resp_Before, Resp_After)
     cat("Respond 组 WGCNA 矩阵样本数:", ncol(Respond_WGCNA_Matrix), "\n")
 } else {
     Respond_WGCNA_Matrix <- NULL
 }
if (!is.null(NonResp_Before) && !is.null(NonResp_After)) {
     NonRespond_WGCNA_Matrix <- cbind(NonResp_Before, NonResp_After)
     cat("\nNon-respond 组 WGCNA 矩阵样本数:", ncol(NonRespond_WGCNA_Matrix), "\n")
 } else {
     NonRespond_WGCNA_Matrix <- NULL
 }

Respond_WGCNA_Matrix <- log2(Respond_WGCNA_Matrix)
NonRespond_WGCNA_Matrix <- log2(NonRespond_WGCNA_Matrix)
Respond_WGCNA_Matrix <- t(Respond_WGCNA_Matrix)
NonRespond_WGCNA_Matrix <- t(NonRespond_WGCNA_Matrix)

sample_names_R <- rownames(Respond_WGCNA_Matrix)
Time_Status_R <- ifelse(str_detect(sample_names_R, "A$"), 0, 1)
Respond_datTraits <- data.frame(
     Treatment_Status = Time_Status_R,
     Treatment_Response = 1, 
     row.names = sample_names_R
 )

sample_names_N <- rownames(NonRespond_WGCNA_Matrix)
Time_Status_N <- ifelse(str_detect(sample_names_N, "A$"), 0, 1)
NonRespond_datTraits <- data.frame(
     Treatment_Status = Time_Status_N,
     Treatment_Response = 0, 
     row.names = sample_names_N
 )

gene_vars_R <- apply(Respond_WGCNA_Matrix, 2, var, na.rm = TRUE)
top_genes_R <- names(sort(gene_vars_R, decreasing = TRUE)[1:1000])
Respond_WGCNA_Matrix <- Respond_WGCNA_Matrix[, top_genes_R]
sft_R <- pickSoftThreshold(Respond_WGCNA_Matrix, powerVector = c(1:30), verbose = 5)
softPower_R <- 25
net_R <- blockwiseModules(Respond_WGCNA_Matrix, power = softPower_R, TOMType = "unsigned",
                          minModuleSize = 10, reassignThreshold = 0,
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, saveTOMs = FALSE,
                          verbose = 3)
moduleColors_R <- labels2colors(net_R$colors)
MEs_R <- moduleEigengenes(Respond_WGCNA_Matrix, moduleColors_R)$eigengenes
MEs_R <- orderMEs(MEs_R)
moduleTraitCor_R <- cor(MEs_R, Respond_datTraits, use = "p")
moduleTraitPvalue_R <- corPvalueStudent(moduleTraitCor_R, nrow(Respond_WGCNA_Matrix))
final_results_R <- data.frame(
     Module = rownames(moduleTraitCor_R),
     Module_Size = as.numeric(table(moduleColors_R)[gsub("ME", "", rownames(moduleTraitCor_R))]),
     Treatment_Status_Correlation = round(moduleTraitCor_R[, "Treatment_Status"], 4),
     Treatment_Status_Pvalue = round(moduleTraitPvalue_R[, "Treatment_Status"], 4)
 )
module_assignment_R <- data.frame(
          Protein_ID = colnames(Respond_WGCNA_Matrix),
          Module_Color = moduleColors_R,
          stringsAsFactors = FALSE
          )
write.xlsx(final_results_R, paste0(output_dir, "/WGCNA_Final_Results_R_log2.xlsx"), rowNames = FALSE)
write.xlsx(module_assignment_R, paste0(output_dir, "/Protein_Module_Assignment_R_log2.xlsx"), rowNames = FALSE)

gene_vars_N <- apply(NonRespond_WGCNA_Matrix, 2, var, na.rm = TRUE)
top_genes_N <- names(sort(gene_vars_N, decreasing = TRUE)[1:1000])
NonRespond_WGCNA_Matrix <- NonRespond_WGCNA_Matrix[, top_genes_N]
sft_NR <- pickSoftThreshold(NonRespond_WGCNA_Matrix, powerVector = c(1:30), verbose = 5)
softPower_NR <- 9
net_NR <- blockwiseModules(NonRespond_WGCNA_Matrix, power = softPower_NR, TOMType = "unsigned",
                          minModuleSize = 10, reassignThreshold = 0,
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, saveTOMs = FALSE,
                          verbose = 3)
moduleColors_NR <- labels2colors(net_NR$colors)
MEs_NR <- moduleEigengenes(NonRespond_WGCNA_Matrix, moduleColors_NR)$eigengenes
MEs_NR <- orderMEs(MEs_NR)
moduleTraitCor_NR <- cor(MEs_NR, NonRespond_datTraits, use = "p")
moduleTraitPvalue_NR <- corPvalueStudent(moduleTraitCor_NR, nrow(NonRespond_WGCNA_Matrix))
final_results_NR <- data.frame(
               Module = rownames(moduleTraitCor_NR),
               Module_Size = as.numeric(table(moduleColors_NR)[gsub("ME", "", rownames(moduleTraitCor_NR))]),
               Treatment_Status_Correlation = round(moduleTraitCor_NR[, "Treatment_Status"], 4),
               Treatment_Status_Pvalue = round(moduleTraitPvalue_NR[, "Treatment_Status"], 4)
 )
module_assignment_NR <- data.frame(
          Protein_ID = colnames(NonRespond_WGCNA_Matrix),
          Module_Color = moduleColors_NR,
          stringsAsFactors = FALSE
          )
write.xlsx(final_results_NR, paste0(output_dir, "/WGCNA_Final_Results_NR_log2.xlsx"), rowNames = FALSE)
write.xlsx(module_assignment_NR, paste0(output_dir, "/Protein_Module_Assignment_NR_log2.xlsx"), rowNames = FALSE)

