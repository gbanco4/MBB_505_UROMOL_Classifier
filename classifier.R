#INSTALL PACKAGES IF NOT ALREADY
packages <- c(
  "tidyverse",
  "survival",
  "glmnet",
  "survminer",
  "survcomp",
  "rms",
  "timeROC"
)

bioc_packages <- c(
  "survcomp"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

missing_bioc <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)
}

# LOAD PACKAGES
library(tidyverse)
library(survival)
library(glmnet)
library(survminer)
library(survcomp)
library(rms)
library(timeROC)

# SET SEED TO ENSURE REPRODUCIBILITY AND PREVENT ANYTHING RANDOM FROM HAPPENING
set.seed(111)

# READ IN DATASETS
knowles_full <- readRDS("knowles_matched_TaLG_final.rds") # THIS COHORT IS FOR EXTERNAL VALIDATION

uromol_full <- readRDS("UROMOL_TaLG.teachingcohort.rds") # THIS COHORT IS FOR TRAINING, TESTING, AND INTERNAL VALIDATION

# REMOVE PATIENTS WITH NA RECURRENCE DATA AND FORMAT COLUMNS
uromol <- uromol_full %>%
  filter(!is.na(Recurrence)) %>%
  filter(!is.na(RFS_time)) %>%
  filter(RFS_time > 0) %>%
  mutate(
    Sex = as.factor(Sex),
    Smoking = as.factor(Smoking),
    Tumor.stage = as.factor(Tumor.stage),
    Tumor.grade = as.factor(Tumor.grade),
    Concomitant.CIS = as.factor(Concomitant.CIS),
    Tumor.size = as.factor(Tumor.size),
    Incident.tumor = as.factor(Incident.tumor),
    EAU.risk = as.factor(EAU.risk),
    BCG = as.factor(BCG),
    UROMOL2021.classification = as.factor(UROMOL2021.classification),
    Age = as.numeric(Age),
    FUtime_days. = as.numeric(FUtime_days.),
    Recurrence = as.numeric(Recurrence)
  )
rownames(uromol) <- uromol$UROMOL.ID

knowles <- knowles_full %>%
  filter(!is.na(Recurrence)) %>%
  filter(!is.na(RFS_time)) %>%
  filter(RFS_time > 0)

# CREATE OBJECT WITH CLINICAL VARIABLES TO SEE WHICH ARE INFORMATIVE FOR PROGNOSIS
uromol_clinical_variables_full <- uromol %>%
  select(
    Age,
    Sex,
    Smoking,
    Tumor.stage,
    Tumor.grade,
    Tumor.size,
    Concomitant.CIS,
    Incident.tumor,
    EAU.risk,
    BCG,
    UROMOL2021.classification
  ) %>% select(where(~ !is.factor(.x) || nlevels(.x) > 1)) # REMOVE COLUMNS THAT ONLY HAVE ONE VALUE

clinical_vars <- colnames(uromol_clinical_variables_full)

uromol_clinical_variables_no_na <- uromol_clinical_variables_full %>%
  filter(complete.cases(.)) # REMOVE ROWS WITH NA VALUES

# SUBSET UROMOL TO ROWS WITHOUT ANY NA VALUES FOR THESE VARIABLES FOR CREATING SURVIVAL OBJECT
uromol_subset_by_na_clin <- uromol[rownames(uromol_clinical_variables_no_na), ]

# CREATE SURVIVAL OBJECT FOR SUBSET DATA
uromol_survival_object_subset <- with(uromol_subset_by_na_clin, Surv(RFS_time, Recurrence))

# CREATE SURVIVAL OBJECT FOR FULL DATA
uromol_survival_object_full <- with(uromol, Surv(RFS_time, Recurrence))

# COX UNIVARIABLE TESTING ON SUBSET UROMOL DATA
univariable_cox_subset <- lapply(clinical_vars, function(v) {
  f <- as.formula(paste("Surv(RFS_time, Recurrence) ~", v))
  summary(coxph(f, data = uromol_subset_by_na_clin))
})

names(univariable_cox_subset) <- clinical_vars

# EXTRACT CLINICAL VARIABLES THAT HAVE RELATIONSHIP TO RISK
clin_var_pvals_subset <- sapply(univariable_cox_subset, function(x) x$logtest["pvalue"])
clin_var_pvals_subset

clin_var_kept_subset <- names(clin_var_pvals_subset)[clin_var_pvals_subset < 0.20]
clin_var_kept_subset

clin_var_hr_subset <- lapply(univariable_cox_subset, function(x) {
  exp(x$coefficients[, "coef"])
})
clin_var_hr_subset

# CRTEATE TABLE WITH P-VALUE FOR RETAINED CLINICAL VARIABLES
clinical_results_subset <- data.frame(
  Variable = names(univariable_cox_subset),
  P_value = clin_var_pvals_subset,
  stringsAsFactors = FALSE
)
clinical_results_subset

# COX UNIVARIABLE TESTING ON FULL UROMOL DATA
univariable_cox_full <- lapply(clinical_vars, function(v) {
  f <- as.formula(paste("Surv(RFS_time, Recurrence) ~", v))
  summary(coxph(f, data = uromol))
})

names(univariable_cox_full) <- clinical_vars

clin_var_pvals_full <- sapply(univariable_cox_full, function(x) x$logtest["pvalue"])

clin_var_kept_full <- names(clin_var_pvals_full)[clin_var_pvals_full < 0.2]

clin_var_p_val_table <- data.frame(
  Variable = clinical_vars,
  P_subset = clin_var_pvals_subset,
  P_full = clin_var_pvals_full
)

clin_var_p_val_table_sig <- clin_var_p_val_table %>% filter(P_subset < 0.2 | P_full < 0.2)

# COX MODEL OF JUST EAU.risk IN SUBSET UROMOL DATASET (BASELINE FOR COMPARISON)
cox_eau <- coxph(uromol_survival_object_subset ~ EAU.risk, data = uromol_subset_by_na_clin)
cox_eau

# EVALUATE ABILITY OF EAU.risk TO PREDICT RFS_time and Recurrence IN FULL UROMOL DATASET USING COX (BASELINE FOR COMPARISON)
cox_eau_full <- coxph(uromol_survival_object_full ~ EAU.risk, data = uromol)
cox_eau_full

# EVALUATE ABILITY OF EAU.RISK TO PREDICT RFS_time and Recurrence IN FULL UROMOL DATASET USING KAPLAN-MEIER CURVE
eau_km <- survfit(Surv(RFS_time, Recurrence) ~ EAU.risk, data = uromol)

ggsurvplot(
  eau_km,
  data = uromol,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "EAU Risk Group",
  legend.labs = c("High", "Intermediate", "Low"),
  palette = c("Red", "Blue", "Green"),
  xlab = "Time (months)",
  ylab = "Recurrence-free survival probability by \n EAU risk group",
  ggtheme = theme_minimal()
)

# EVALUATE ABILITY OF EAU.RISK TO PREDICT RFS_time and Recurrence IN SUBSET UROMOL DATASET USING KAPLAN-MEIER CURVE
eau_km_subset <- survfit(Surv(RFS_time, Recurrence) ~ EAU.risk, data = uromol_subset_by_na_clin)

ggsurvplot(
  eau_km_subset,
  data = uromol_subset_by_na_clin,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "EAU Risk Group",
  legend.labs = c("High", "Intermediate", "Low"),
  palette = c("Red", "Blue", "Green"),
  xlab = "Time (months)",
  ylab = "Recurrence-free survival probability /n by EAU risk group (subset \n to no NA clinical variables)",
  ggtheme = theme_minimal()
)

# EVALUATE ABILITY OF EAU.RISK TO PREDICT RFS_time and Recurrence IN FULL UROMOL DATASET USING CONCORDANCE-INDEX
# IDENTIFY ROWS USED IN COX MODEL (NO MISSING VALUES FOR EAU.risk)
rows_in_cox <- rownames(model.frame(cox_eau_full))

linear_predictor_eau <- predict(cox_eau_full, type = "lp")

# C-INDEX
cindex_eau <- concordance.index(
  linear_predictor_eau,
  surv.time = uromol[rows_in_cox, ]$RFS_time,
  surv.event = uromol[rows_in_cox, ]$Recurrence
)$c.index
cindex_eau

# CREATE OBJECT WITH EXPRESSION DATA
uromol_expression <- uromol$exprs
if (is.matrix(uromol_expression) && is.null(rownames(uromol_expression))) {
  rownames(uromol_expression) <- rownames(uromol)
}

# TRANSPOSE
uromol_expr_mat <- t(uromol$exprs)

# CHECK IF UROMOL RNASEQ CONTAINS BATCH EFFECTS BY GENERATING PCA PLOTS OF exprs WITH OTHER VARIABLES AS col = VARIABLE
pca <- prcomp(t(uromol$exprs))
plot(pca$x[,1],
     pca$x[,2],
     col = as.factor(uromol$Smoking),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2")
# THE PCA DOES NOT APPEAR TO SHOW BATCH EFFECTS. I GENERATED PCA PLOTS USING MULTIPLE CLINICAL VARIABLES (NOT SHOWN HERE)

# FILTER OUT LOW EXPRESSION GENES (BELOW ARBITRARY THRESHOLD OF MEAN EXPRESSION = 3, WHICH I DETERMINED BY LOOKING AT A HISTOGRAM OF uromol$exprs)
high_expression_genes <- rowMeans(uromol_expr_mat) > 3

expr_high <- uromol_expr_mat[high_expression_genes, ]

# FILTER OUT LOW VARIANCE GENES (RETAIN TOP 75% OF GENES BY VARIANCE)
gene_var <- apply(expr_high, 1, var)

expr_high_var <- expr_high[gene_var > quantile(gene_var, 0.25), ]

# Z-SCORE FILTERED EXPRESSION DATA (SINCE IT ALREADY APPEARS TO HAVE BEEN LOG-SCALED BEFORE PUBLICATION IN THIS DATASET, MEANING QUANTILE NORMALIZATION DOESN'T MAKE SENSE)
genes_z_scored <- t(scale(t(expr_high_var)))

# ENSURE UROMOL IDs ARE ALL THE SAME AS IN THE ORIGINAL UROMOL OBJECT AND CREATE SURVIVAL OBJECT FROM THIS
common_ids <- intersect(
  rownames(uromol),
  colnames(genes_z_scored)
)

expression_z_subset <- genes_z_scored[, common_ids]

uromol_matched <- uromol[common_ids, ]

uromol_survival_object_full <- with(uromol_matched, Surv(RFS_time, Recurrence))

# UNIVARIABLE COX FOR EACH GENE TO REDUCE DIMENSIONALITY AND IDENTIFY GOOD CANDIDATE GENES
gene_cox <- apply(expression_z_subset, 1, function(g) {
  summary(coxph(uromol_survival_object_full ~ g))
})

# SELECT GENES BASED ON P VALUE (P-VALUE < 0.02)
gene_pvals <- sapply(gene_cox, function(x) x$coefficients[, "Pr(>|z|)"])
genes_kept <- names(gene_pvals)[gene_pvals < 0.02]
expression_screened <- expression_z_subset[genes_kept, ]

# BUILD RNASEQ ONLY MODEL
# TRANSPOSE
x_rna <- t(expression_screened)

y_surv <- uromol_survival_object_full

fit_rna <- cv.glmnet(
  x = x_rna,
  y = y_surv,
  family = "cox",
  alpha = 1
)

coef_rna <- coef(fit_rna, s = "lambda.min")

genes_rna <- rownames(coef_rna)[coef_rna[,1] != 0]

# PRINT GENES

genes_rna

# DETERMINE RNASEQ MOLECULAR RISK SCORE
rna_score <- x_rna[, genes_rna] %*% coef_rna[genes_rna, ]
uromol$rna_score <- as.numeric(rna_score)

# EVALUATE RNASEQ ONLY MODEL
# C-INDEX
cindex_rna <- concordance.index(
  uromol$rna_score,
  surv.time = uromol$RFS_time,
  surv.event = uromol$Recurrence
)$c.index
cindex_rna

# KM CURVE BY RNA-SCORE RISK GROUP
uromol$rna_group <- cut(
  uromol$rna_score,
  breaks = quantile(uromol$rna_score, c(0, 1/2, 1)),
  include.lowest = TRUE,
  labels = c("Low", "High")
)

fit_rna_km <- survfit(uromol_survival_object_full ~ rna_group, data = uromol)

ggsurvplot(
  fit_rna_km,
  data = uromol,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "RNA-score \n risk group",
  legend.labs = c("Low", "High"),
  palette = c("green", "red"),
  xlab = "Time (months)",
  ylab = "Recurrence-free survival probability",
  ggtheme = theme_minimal(),
  title = "UROMOL RNA-score Kaplan-Meier Curve"
)

# BUILD RNASEQ + CLINICAL FEATURE MODEL
clin_vars <- intersect(
  sub("\\.pvalue$", "", clin_var_kept_full),
  colnames(uromol)
)

formula_integrated <- as.formula(
  paste("uromol_survival_object_full ~ rna_score +",
        paste(clin_vars, collapse = " + "))
)

cox_integrated <- coxph(formula_integrated, data = uromol)
summary(cox_integrated)

# EVALUATE RNASEQ + CLINICAL FEATURES MODEL AND COMPARE RNASEQ TO RNASEQ + CLINICAL FEATURES
# C-INDEX
# Identify rows used in the integrated Cox model
rows_in_integrated <- rownames(model.frame(cox_integrated))

# Linear predictor (already aligned to these rows)
lp_integrated <- predict(cox_integrated, type = "lp")

cindex_integrated <- concordance.index(
  lp_integrated,
  surv.time  = uromol[rows_in_integrated, "RFS_time"],
  surv.event = uromol[rows_in_integrated, "Recurrence"]
)$c.index
cindex_integrated

# INTERNAL VALIDATION
dd <- datadist(uromol)
options(datadist = "dd")

fit_val_rna <- cph(
  Surv(RFS_time, Recurrence) ~ rna_score,
  data = uromol,
  x = TRUE, y = TRUE, surv = TRUE
)

val_rna <- validate(fit_val_rna, B = 100)

cindex_rna_internal <- (val_rna["Dxy", "index.corrected"] + 1) / 2
cindex_rna_internal

fit_val_integrated <- cph(
  Surv(RFS_time, Recurrence) ~ rna_score +
    Sex + Smoking + Concomitant.CIS +
    Incident.tumor + EAU.risk + BCG +
    UROMOL2021.classification,
  data = uromol,
  x = TRUE, y = TRUE, surv = TRUE
)

val_integrated <- validate(fit_val_integrated, B = 100)

cindex_integrated_internal <- (val_integrated["Dxy", "index.corrected"] + 1) / 2
cindex_integrated_internal

if (cindex_integrated_internal > cindex_rna_internal) {
  "Integrated model has a higher C-index"
} else if (cindex_rna_internal > cindex_integrated_internal) {
  "RNA-only model has a higher C-index"
} else {
  "Both models have the same C-index"
}
# C-INDEX FOR RNA ONLY MODEL HAS A HIGHER C-INDEX, I WILL PROCEED WITH THIS MODEL

# EXTERNAL VALIDATION
# EXTRACT MEAN AND SD FROM UROMOL COHORT RNASEQ Z-SCORING AND APPLY TO KNOWLES MICROARRAY DATA
# TRAINING COHORT SCALING PARAMETERS
gene_means <- rowMeans(expr_high_var)
gene_sds   <- apply(expr_high_var, 1, sd)

# IDENTIFY GENES PRESENT IN BOTH DATASETS
common_genes <- intersect(genes_rna, colnames(knowles$exprs))

# EXTRACT MATCHING GENES FROM KNOWLES (genes are columns)
knowles_expr <- knowles$exprs[, common_genes, drop = FALSE]

# Z-SCORE KNOWLES USING TRAINING (UROMOL) STATISTICS
knowles_expr_z <- sweep(t(knowles_expr), 1, gene_means[common_genes], "-")
knowles_expr_z <- sweep(knowles_expr_z, 1, gene_sds[common_genes], "/")

# RNA RISK SCORE (NO REFITTING)
rna_score_knowles <- t(knowles_expr_z) %*% coef_rna[common_genes, ]
knowles$rna_score <- as.numeric(rna_score_knowles)

# C-INDEX (EXTERNAL VALIDATION)
cindex_knowles <- concordance.index(
  knowles$rna_score,
  surv.time  = knowles$RFS_time,
  surv.event = knowles$Recurrence
)$c.index
cindex_knowles

# KM PLOT (MEDIAN SPLIT)
knowles$rna_group <- cut(
  knowles$rna_score,
  breaks = quantile(knowles$rna_score, c(0, 0.5, 1)),
  include.lowest = TRUE,
  labels = c("Low", "High")
)

fit_km_knowles <- survfit(
  Surv(RFS_time, Recurrence) ~ rna_group,
  data = knowles
)

ggsurvplot(
  fit_km_knowles,
  data = knowles,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  legend.title = "RNA-score \n risk group",
  legend.labs = c("Low", "High"),
  palette = c("green", "red"),
  xlab = "Time (months)",
  ylab = "Recurrence-free survival probability",
  ggtheme = theme_minimal(),
  title = "Knowles RNA-score Kaplan-Meier Curve"
)

# UROMOL TIME-DEPENDENT AUC

uromol_auc <- timeROC(
  T = uromol$RFS_time,
  delta = uromol$Recurrence,
  marker = uromol$rna_score,
  cause = 1,
  times = c(12, 24, 36),
  iid = TRUE
)

print(uromol_auc)

# PLOT
plot(
  uromol_auc,
  time = 12,
  col = "red",
  title = FALSE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "UROMOL TIME-DEPENDENT AUC (RNA-score)"
)
plot(
  uromol_auc,
  time = 24,
  add = TRUE,
  col = "blue"
)
plot(
  uromol_auc,
  time = 36,
  add = TRUE,
  col = "green"
)
title(main = "UROMOL Time-Dependent AUC")
legend(
  "bottomright",
  legend = c(
    paste0("12 months (AUC = ", round(uromol_auc$AUC[1], 3), ")"),
    paste0("24 months (AUC = ", round(uromol_auc$AUC[2], 3), ")"),
    paste0("36 months (AUC = ", round(uromol_auc$AUC[3], 3), ")")
  ),
  col = c("red", "blue", "green"),
  lwd = 2
)

# KNOWLES TIME-DEPENDENT AUC
knowles_auc <- timeROC(
  T = knowles$RFS_time,
  delta = knowles$Recurrence,
  marker = knowles$rna_score,
  cause = 1,
  times = c(12, 24, 36),
  iid = TRUE
)

print(knowles_auc)

# PLOT
plot(
  knowles_auc,
  time = 12,
  col = "red",
  title = FALSE,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "KNOWLES TIME-DEPENDENT AUC (RNA-score)"
)
plot(
  knowles_auc,
  time = 24,
  add = TRUE,
  col = "blue"
)
plot(
  knowles_auc,
  time = 36,
  add = TRUE,
  col = "green"
)
title(main = "Knowles Time-Dependent AUC")
legend(
  "bottomright",
  legend = c(
    paste0("12 months (AUC = ", round(knowles_auc$AUC[1], 3), ")"),
    paste0("24 months (AUC = ", round(knowles_auc$AUC[2], 3), ")"),
    paste0("36 months (AUC = ", round(knowles_auc$AUC[3], 3), ")")
  ),
  col = c("red", "blue", "green"),
  lwd = 2
)
