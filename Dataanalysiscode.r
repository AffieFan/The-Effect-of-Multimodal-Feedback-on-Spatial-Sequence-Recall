library(MASS)
library(car)

# ── Load Data ─────────────────────────────────────────────────
data <- read.csv("/Users/affiefan/Desktop/Data Analysis/DataAnalysisFull.csv")
colnames(data) <- c("Participant", "FinalLevel", "TimeTotal", "Clicks", "TimePerClick", "Type", "Age", "Gender")
data$Gender <- trimws(data$Gender)

visual     <- subset(data, Type == "Visual")
multimodal <- subset(data, Type == "Multimodal")
acc_vis  <- visual$FinalLevel;    acc_mm   <- multimodal$FinalLevel
prec_vis <- visual$TimePerClick;  prec_mm  <- multimodal$TimePerClick

# ── Descriptive Statistics ────────────────────────────────────
cat("=== Accuracy (Final Level Reached) ===\n")
cat(sprintf("Visual     — Mean: %.4f  SD: %.4f  Var: %.4f\n", mean(acc_vis), sd(acc_vis), var(acc_vis)))
cat(sprintf("Multimodal — Mean: %.4f  SD: %.4f  Var: %.4f\n", mean(acc_mm),  sd(acc_mm),  var(acc_mm)))
cat("\n=== Precision (Time Per Click, s) ===\n")
cat(sprintf("Visual     — Mean: %.4f  SD: %.4f  Var: %.4f\n", mean(prec_vis), sd(prec_vis), var(prec_vis)))
cat(sprintf("Multimodal — Mean: %.4f  SD: %.4f  Var: %.4f\n", mean(prec_mm),  sd(prec_mm),  var(prec_mm)))

# ── Variable A: Poisson Verification ─────────────────────────
cat("\n=== Variable A: Poisson Equidispersion Check ===\n")
cat(sprintf("Visual     — Mean: %.4f  Var: %.4f\n", mean(acc_vis), var(acc_vis)))
cat(sprintf("Multimodal — Mean: %.4f  Var: %.4f\n", mean(acc_mm),  var(acc_mm)))

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
qqplot(qpois(ppoints(length(acc_vis)), lambda = mean(acc_vis)), acc_vis,
       main = "Poisson QQ : Visual Accuracy",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", pch = 16)
abline(0, 1, col = "red", lwd = 2)
qqplot(qpois(ppoints(length(acc_mm)), lambda = mean(acc_mm)), acc_mm,
       main = "Poisson QQ : Multimodal Accuracy",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", pch = 16)
abline(0, 1, col = "red", lwd = 2)

poisson_gof <- function(x, label) {
  lam <- mean(x); brk <- min(x):max(x)
  obs <- as.numeric(table(factor(x, levels = brk)))
  exp <- dpois(brk, lam) * length(x); keep <- exp >= 1
  if (sum(keep) < 2) { cat(label, ": Too few bins — interpret QQ plot instead.\n"); return(invisible(NULL)) }
  chi2 <- sum((obs[keep] - exp[keep])^2 / exp[keep]); df <- max(sum(keep) - 2, 1)
  cat(sprintf("%s — chi-sq: %.3f  df: %d  p-value: %.4f\n(p > 0.05 = no significant departure from Poisson)\n\n", label, chi2, df, pchisq(chi2, df, lower.tail = FALSE)))
}
cat("\n=== Variable A: Chi-Squared GoF for Poisson ===\n")
poisson_gof(acc_vis, "Visual    "); poisson_gof(acc_mm, "Multimodal")

# ── Variable B: Normal Verification ──────────────────────────
cat("=== Variable B: Shapiro-Wilk Normality Test ===\n")
sw_vis <- shapiro.test(prec_vis); sw_mm <- shapiro.test(prec_mm)
cat(sprintf("Visual     — W: %.4f  p-value: %.4f\n", sw_vis$statistic, sw_vis$p.value))
cat(sprintf("Multimodal — W: %.4f  p-value: %.4f\n", sw_mm$statistic,  sw_mm$p.value))
cat("(p > 0.05 supports normality)\n\n")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
qqnorm(prec_vis, main = "Normal QQ — Visual Precision",     pch = 16); qqline(prec_vis, col = "red", lwd = 2)
qqnorm(prec_mm,  main = "Normal QQ — Multimodal Precision", pch = 16); qqline(prec_mm,  col = "red", lwd = 2)

# ── Confidence Intervals ──────────────────────────────────────
ci_poisson <- function(x, label, conf = 0.95) {
  n <- length(x); total <- sum(x); alpha <- 1 - conf
  cat(sprintf("%s  lambda_hat: %.4f  95%% CI: [%.4f, %.4f]\n", label, mean(x),
              qgamma(alpha/2, total, n), qgamma(1 - alpha/2, total + 1, n)))
}
ci_normal <- function(x, label, conf = 0.95) {
  n <- length(x); m <- mean(x); se <- sd(x)/sqrt(n); t_cv <- qt((1+conf)/2, df = n-1)
  cat(sprintf("%s  Mean: %.4f  95%% CI: [%.4f, %.4f]\n", label, m, m - t_cv*se, m + t_cv*se))
}
cat("=== Variable A: 95% CIs (Exact Poisson) ===\n")
ci_poisson(acc_vis, "Visual:    "); ci_poisson(acc_mm, "Multimodal:")
cat("\n=== Variable B: 95% CIs (t-based) ===\n")
ci_normal(prec_vis, "Visual:    "); ci_normal(prec_mm, "Multimodal:")

# ── Inferential Statistics ────────────────────────────────────
cat("\n=== Variable A: Two-Sample Poisson Rate Test ===\n")
cat("H0: lambda_visual = lambda_multimodal\nH1: lambda_visual ≠ lambda_multimodal\n\n")
print(poisson.test(c(sum(acc_vis), sum(acc_mm)), c(length(acc_vis), length(acc_mm)), alternative = "two.sided"))

cat("\n=== Variable B: Welch's t-test ===\n")
cat("H0: mu_visual = mu_multimodal\nH1: mu_visual ≠ mu_multimodal\n\n")
print(t.test(prec_vis, prec_mm, alternative = "two.sided"))

# ── Speed-Accuracy Trade-off ──────────────────────────────────
cat("\n=== Speed-Accuracy Trade-off: Spearman Correlation ===\n")
cat("H0: no correlation between Time Per Click and Final Level Reached\n")
cat("H1: there is a correlation between Time Per Click and Final Level Reached\n\n")
print(cor.test(data$TimePerClick, data$FinalLevel, method = "spearman"))

# ── Assumption Check: Levene's Test ──────────────────────────
cat("\n=== Levene's Test: Homogeneity of Variance (Variable B) ===\n")
cat("H0: Variance of Visual = Variance of Multimodal\nH1: Variance of Visual ≠ Variance of Multimodal\n\n")
leveneTest(TimePerClick ~ Type, data = data)

# ── Demographic Analyses ──────────────────────────────────────
cat("\n=== Gender t-test: Final Level Reached ===\n")
cat("H0: mu_female = mu_male\nH1: mu_female ≠ mu_male\n\n")
print(t.test(FinalLevel ~ Gender, data = data, alternative = "two.sided"))

cat("\n=== Gender t-test: Time Per Click ===\n")
cat("H0: mu_female = mu_male\nH1: mu_female ≠ mu_male\n\n")
print(t.test(TimePerClick ~ Gender, data = data, alternative = "two.sided"))

cat("\n=== Age vs Final Level Reached: Spearman Correlation ===\n")
print(cor.test(data$Age, data$FinalLevel, method = "spearman"))
cat("\n=== Age vs Time Per Click: Spearman Correlation ===\n")
print(cor.test(data$Age, data$TimePerClick, method = "spearman"))

# ── Plots ─────────────────────────────────────────────────────
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
boxplot(FinalLevel ~ Type, data = data, main = "Final Level Reached by Feedback Type",
        xlab = "Feedback Type", ylab = "Final Level Reached", col = c("lightblue", "lightcoral"), pch = 16)
boxplot(TimePerClick ~ Type, data = data, main = "Time Per Click by Feedback Type",
        xlab = "Feedback Type", ylab = "Time Per Click (s)", col = c("lightblue", "lightcoral"), pch = 16)
plot(data$TimePerClick, data$FinalLevel, main = "Speed-Accuracy Trade-off",
     xlab = "Time Per Click (s)", ylab = "Final Level Reached", pch = 16,
     col = ifelse(data$Type == "Visual", "lightblue", "lightcoral"))
abline(lm(FinalLevel ~ TimePerClick, data = data), col = "black", lwd = 2)
legend("topright", legend = c("Visual", "Multimodal"), col = c("lightblue", "lightcoral"), pch = 16)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
boxplot(FinalLevel ~ Gender, data = data, main = "Final Level Reached by Gender",
        xlab = "Gender", ylab = "Final Level Reached", col = c("lightpink", "lightsteelblue"), pch = 16)
boxplot(TimePerClick ~ Gender, data = data, main = "Time Per Click by Gender",
        xlab = "Gender", ylab = "Time Per Click (s)", col = c("lightpink", "lightsteelblue"), pch = 16)

age_vis <- data$Age[data$Type == "Visual"]; age_mm <- data$Age[data$Type == "Multimodal"]
barplot(rbind(table(factor(age_vis, levels = 19:22)), table(factor(age_mm, levels = 19:22))),
        beside = TRUE, col = c("lightblue", "lightcoral"), names.arg = c("19", "20", "21", "22"),
        main = "Age Distribution by Feedback Type", xlab = "Age", ylab = "Number of Participants",
        legend.text = c("Visual", "Multimodal"), args.legend = list(x = "topright"))

gender_vis <- table(data$Gender[data$Type == "Visual"]); gender_mm <- table(data$Gender[data$Type == "Multimodal"])
barplot(t(rbind(Visual = c(Female = gender_vis["Female"], Male = gender_vis["Male"]),
                Multimodal = c(Female = gender_mm["Female"], Male = gender_mm["Male"]))),
        beside = TRUE, col = c("lightpink", "lightsteelblue"), names.arg = c("Visual", "Multimodal"),
        main = "Gender Distribution by Feedback Type", xlab = "Feedback Type", ylab = "Number of Participants",
        legend.text = c("Female", "Male"), args.legend = list(x = "topright"))

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(data$Age, data$FinalLevel, main = "Age vs Final Level Reached",
     xlab = "Age", ylab = "Final Level Reached", pch = 16,
     col = ifelse(data$Gender == "Female", "lightpink", "lightsteelblue"))
abline(lm(FinalLevel ~ Age, data = data), col = "black", lwd = 2)
legend("topright", legend = c("Female", "Male"), col = c("lightpink", "lightsteelblue"), pch = 16)
plot(data$Age, data$TimePerClick, main = "Age vs Time Per Click",
     xlab = "Age", ylab = "Time Per Click (s)", pch = 16,
     col = ifelse(data$Gender == "Female", "lightpink", "lightsteelblue"))
abline(lm(TimePerClick ~ Age, data = data), col = "black", lwd = 2)
legend("topright", legend = c("Female", "Male"), col = c("lightpink", "lightsteelblue"), pch = 16)