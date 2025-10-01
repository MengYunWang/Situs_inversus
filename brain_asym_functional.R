# This code is to do the linear regression on the brain functional asym
# 2. use ANOVA to test group difference on four measures
# 3. test variants across three groups
# 4. levene test on atypical patterns across three groups
# 5.

# WANG. 19-Jun-2025

# Remove all objects created before to prevent clash
rm(list=ls())

# load relevant libraries
library(openxlsx)
library(readxl)
library(dplyr)
library(tidyverse)
library(car)
library(emmeans)
library(effectsize)
library(heplots)
library(nlme)

# Set the working directory to the path
setwd("/Users/wang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Situs_inversus/")
# setwd("/Users/joeywang/Library/CloudStorage/OneDrive-RadboudUniversiteit/Research_Project/Situs_inversus/")

############################################################
## 1.  Loading the data   ##################################
############################################################

data2read <- read_excel("brain_asym_SI.xlsx", sheet = "both") %>%
  select(-1:-5, -EHI, -years_education)


df <- data2read |>
  mutate(Solved = relevel(factor(Solved), ref = "0"))

############################################################
## 2.  Test the group difference between four measures######
############################################################
dvs <- c("WGEN", "Praxis",
         "Spatatt", "Face")

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

##################################################################
## 3.  Test the group variance difference between four measures###
##################################################################

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


############################################################
# 4.  Test the group varaince difference for N_typical
############################################################
# # Levene (mean-centred)            – classic version
# leveneTest(N_atypical ~ Solved, data = df)

# Brown–Forsythe (median-centred)  – more robust
leveneTest(N_atypical ~ Solved, data = df, center = median)

# # Fligner–Killeen (rank-based)     – fully non-parametric
# fligner.test(N_atypical ~ Solved, data = df)

# Quick descriptive variances & SDs per group
tapply(df$N_atypical, df$Solved, var)
tapply(df$N_atypical, df$Solved, sd)


###########################
## 5.  Test handedness  ###
###########################
handedness <- matrix(c(5, 13,
                       6, 13),
              nrow = 2,
              byrow = TRUE)

colnames(handedness) <- c("left_hand", "right_hand")
rownames(handedness) <- c("solved", "unsolved")

groups <- rownames(handedness)

fisher.test(handedness)

###########################
## 6.  Plot the figure ###
###########################
# Create new grouping variable: Control vs Case (1+2)
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





# ############################################################
# # 5.  Heteroskedastic ANCOVA via GLS for N_typical----------
# ############################################################
# gls_mod <- gls(
#   N_atypical ~ Solved + Age + Sex,
#   data     = df,
#   weights  = varIdent(form = ~ 1 | Solved)   # separate σ² per group
# )
# 
# anova(gls_mod)      # omnibus F tests
# summary(gls_mod)    # coefficients, group-specific σ̂
# 
# # Residual df for effect-size CIs
# df_res <- gls_mod$dims$N - length(coef(gls_mod))
# 
# ############################################################
# # 6.  Adjusted means & pairwise contrasts (Solved) ---------
# ############################################################
# emm_solved   <- emmeans(gls_mod, "Solved", data=df)
# summary(emm_solved)                          # adjusted means ± SE
# 
# pairs_solved <- contrast(emm_solved, "pairwise", adjust = "tukey")
# summary(pairs_solved)                        # Tukey-style pairwise p-values
# 
# # Hedges g (bias-corrected Cohen d) for those contrasts
# g_solved <- eff_size(
#   pairs_solved,
#   sigma        = sigma(gls_mod),
#   edf          = df_res,
#   bias.adjust  = TRUE
# )
# g_solved
# 
# # Quick plot of the marginal means with 95 % CIs
# emmip(gls_mod, Solved ~ 1, data=df, CIs = TRUE)
# 
# ############################################################
# # 7.  Follow-up on significant Handedness effect ----------
# ############################################################
# emm_hand <- emmeans(gls_mod, "Handedness", data=df)
# 
# # Pairwise contrasts (Tukey correction)
# contrast(emm_hand, "pairwise", adjust = "tukey")
# 
# # Hedges g for Handedness contrasts
# g_hand <- eff_size(
#   contrast(emm_hand, "pairwise"),
#   sigma        = sigma(gls_mod),
#   edf          = df_res,
#   bias.adjust  = TRUE
# )
# g_hand