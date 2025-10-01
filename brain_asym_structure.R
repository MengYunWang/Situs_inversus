# This code is to do the linear regression on the brain structure asym
# 1. 
# 2. 

# WANG. 19-Jun-2025

# Remove all objects created before to prevent clash
rm(list=ls())

library(openxlsx)
library(readxl)
library(dplyr)
library(tidyverse)
library(car)
library(emmeans)
library(effectsize)
library(heplots)
library(patchwork) 
# Set the working directory to the path where your files are located
# setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Situs_inversus/")
setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Situs_inversus/")

############################################################
## 1.  Loading the data   ##################################
############################################################

data2read <- read_excel("brain_asym_SI.xlsx", sheet = "both") %>%
  select(-1:-5, -EHI, -years_education)


df <- data2read |>
  mutate(Solved = relevel(factor(Solved), ref = "0"))

#############################################################################
## 2.  UNIVARIATE ANCOVAs for bending and petalia; transverse_sinus##########
#############################################################################
dvs <- c("Frontal_bending", "Frontal_petalia",
         "Occipital_bending", "Occipital_petalia")

dvs <- c("Left_transverse_sinus", "Right_transverse_sinus")

df[dvs] <- lapply(df[dvs], \(x) as.numeric(as.character(x)))

ancova_results <- map(dvs, \(dv) {
  fmla <- as.formula(paste0(dv, " ~ Solved + Age + Sex"))
  fit  <- lm(fmla, data = df)
  # Type-III SS for covariance-adjusted F
  aov_tab <- car::Anova(fit, type = 3)
  # Partial η²
  eff     <- effectsize::eta_squared(fit, partial = TRUE)
  list(model = fit, anova = aov_tab, es = eff)
}) |> set_names(dvs)

# Pull raw p values for the Group term
p_raw <- map_dbl(ancova_results, \(lst) lst$anova["Solved", "Pr(>F)"])
p_adj <- p.adjust(p_raw, method = "holm")

tidy_univ <- tibble(
  DV     = dvs,
  F      = map_dbl(ancova_results, \(lst) lst$anova["Solved", "F value"]),
  df1    = map_dbl(ancova_results, \(lst) lst$anova["Solved", "Df"]),
  df2    = map_dbl(ancova_results, \(lst) lst$anova["Residuals", "Df"]),
  p_raw,
  p_holm = p_adj,
  partial_eta2 = map_dbl(ancova_results, \(lst)
                         lst$es |> filter(Parameter == "Solved") |> pull(Eta2_partial))
)

print(tidy_univ)

############################################################
## 3.  PAIRWISE COMPARISONS (Tukey, covariate-adjusted) #####
############################################################

library(emmeans)

# funtion to compute effect size and get pairwise group contrasts from an ANCOVA model
hedges_g_pairs <- function(fit, factor_name = "Solved") {
  emm  <- emmeans(fit, factor_name)
  cons <- contrast(emm, "pairwise")  # pairwise mean diffs
  eff_size(cons,
           sigma        = sigma(fit),          # model residual SD
           edf          = df.residual(fit),    # residual df for CIs
           bias.adjust  = TRUE,
           method       = "identity")                # Hedges' g (vs Cohen's d)
}

#--------------petalia and bending
# 1) Occipital_bending
fit_ob <- ancova_results$Occipital_bending$model
pairs_ob <- contrast(emmeans(fit_ob, "Solved"), "pairwise", adjust = "tukey")
summary(pairs_ob)
g_ob <- hedges_g_pairs(fit_ob)
summary(g_ob)

# 2) Occipital_petalia
fit_op <- ancova_results$Occipital_petalia$model
pairs_op <- contrast(emmeans(fit_op, "Solved"), "pairwise", adjust = "tukey")
summary(pairs_op)
g_op <- hedges_g_pairs(fit_op)
summary(g_op)

# 3) Frontal_petalia
fit_fp <- ancova_results$Frontal_petalia$model
pairs_fp <- contrast(emmeans(fit_fp, "Solved"), "pairwise", adjust = "tukey")
summary(pairs_fp)
g_fp <- hedges_g_pairs(fit_fp)
summary(g_fp)

#------------ transvers sinus
# 1) Left_transverse_sinus
fit_ls <- ancova_results$Left_transverse_sinus$model
pairs_ls <- contrast(emmeans(fit_ls, "Solved"), "pairwise", adjust = "tukey")
summary(pairs_ls)
g_ls <- hedges_g_pairs(fit_ls)
summary(g_ls)

# 2) Right_transverse_sinus
fit_rs <- ancova_results$Right_transverse_sinus$model
pairs_rs <- contrast(emmeans(fit_rs, "Solved"), "pairwise", adjust = "tukey")
summary(pairs_rs)
g_rs <- hedges_g_pairs(fit_rs)
summary(g_rs)


############################################################
## 4. (optional)  DIAGNOSTIC PLOTS  ########################
############################################################
# Q-Q plot & Scale-Location for each model
plot(ancova_results$Occipital_bending$model)  # repeat for others

# Homogeneity of regression slopes (Age × Group example)
lm_slopes <- lm(Occipital_bending ~ Solved * Age + Sex, data = df)
car::Anova(lm_slopes, type = 3)




############################################################
## 5.  variation test across three groups. #################
############################################################

result <- lapply(dvs, \(v) {
  fml <- reformulate("Solved", response = v)   # <-- v on the left, Solved on the right
  leveneTest(fml, data = df, center = median)[1, ]
})

out <- do.call(rbind, result)
rownames(out) <- dvs
print(out)

# adjust p-values (e.g., Holm)
out$p_adj <- p.adjust(out[ , "Pr(>F)"], method = "holm")
out



####################################################
## 6. Plot the figures
##################################################
df <- df %>%
  mutate(
    Group2 = case_when(
      Solved == "0" ~ "Control",
      Solved %in% c("1", "2") ~ "Case"
    ),
    Group2 = factor(Group2, levels = c("Control", "Case")),
    Solved = factor(Solved, levels = c("0", "1", "2"),
                    labels = c("Control", "Solved", "Unsolved"))
  )


custom_theme <- theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),               # No inside grids
    axis.line = element_line(size = 1.5, color = "black"),  # X and Y axes visible and thick
    axis.title = element_text(face = "bold", size = 16),    # Bold axis labels
    axis.text = element_text(size = 14),        # Axis tick labels size
    legend.position = "none"
  )

for (var in dvs) {
  
  # p1: Control vs Case
  p1 <- ggplot(df, aes(x = Group2, y = .data[[var]], fill = Group2)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#999999", "Case" = "#009E73")) +
    labs(x = "Group", y = var) +
    custom_theme+ 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
  
  # p2: Control, Solved, Unsolved
  p2 <- ggplot(df, aes(x = Solved, y = .data[[var]], fill = Solved)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#999999", "Solved" = "#E69F00", "Unsolved" = "#56B4E9")) +
    labs(x = "Group", y = var) +
    custom_theme
  
  # Combine and save
  combined_plot <- p1 + p2
  
  ggsave(
    filename = paste0("violin_plot_", var, ".pdf"),
    plot = combined_plot,
    device = cairo_pdf,
    width = 10,
    height = 7,
    units = "in"
  )
}


