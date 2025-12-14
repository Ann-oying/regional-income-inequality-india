colnames(data_ecotrix)
colnames(data_ecotrix) <- c("state", 
                            "nsdp_pc", 
                            "gva_industry", 
                            "social_exp", 
                            "power_pc", 
                            "credit_agri")
##Descriptive statistics

summary(data_ecotrix)

#regress raw data
model <- lm(nsdp_pc ~ gva_industry + social_exp + power_pc + credit_agri, data = data_ecotrix)
summary(model)


###LINEARITY

##Component + Residual Plots and partial regression plots
library(car)
crPlots(model)
###individual_plots
##comp+residual plot
crPlot(model, variable = "gva_industry")
crPlot(model, variable = "social_exp")
crPlot(model, variable = "power_pc")
crPlot(model, variable = "credit_agri")

##Added variable plot
avPlots(model)

###Transformed model
data_ecotrix$ln_nsdp_pc<- log(data_ecotrix$nsdp_pc)
data_ecotrix$ln_gva_industry<- log(data_ecotrix$gva_industry)
data_ecotrix$ln_social_exp<- log(data_ecotrix$social_exp)
data_ecotrix$ln_power_pc<- log(data_ecotrix$power_pc)
data_ecotrix$ln_credit_agri<- log(data_ecotrix$credit_agri)
data_ecotrix$credit_agri_sqrt<- sqrt(data_ecotrix$credit_agri)

model2<- lm(nsdp_pc ~ ln_gva_industry + ln_social_exp + ln_power_pc + credit_agri_sqrt , data = data_ecotrix)
summary(model2)

crPlots(model2)
crPlot(model2, variable = "ln_gva_industry")
crPlot(model2, variable = "ln_social_exp")
crPlot(model2, variable = "ln_power_pc")
crPlot(model2, variable = "credit_agri_sqrt")
avPlots(model2)

###HETEROSCEDASTICITY

##model2
# Extract residuals and fitted values 
resid<- residuals(model2)
std_resid<- rstandard(model2)
std_resid_sq<- std_resid^2
y_hat<- fitted(model2)
     
##plot std_e_sq vs yhat 
plot(y_hat, std_resid_sq,
     main = "Std Residuals² vs Fitted Values (yhat)",
     xlab = "Fitted Values (yhat)",
     ylab = "Standardised Residuals²",
     pch = 19, col = "darkgrey")

# Add LOESS curve
lines(lowess(y_hat, std_resid_sq), col = "blue", lwd = 2)

#BP test
library(lmtest)
bptest(model2)

##model3- using log transfomed nsdp_pc

model3<- lm(ln_nsdp_pc ~ ln_gva_industry + ln_social_exp + ln_power_pc + credit_agri_sqrt , data = data_ecotrix)
summary(model3)

crPlots(model3)
crPlot(model3, variable = "ln_gva_industry")
crPlot(model3, variable = "ln_social_exp")
crPlot(model3, variable = "ln_power_pc")
crPlot(model3, variable = "credit_agri_sqrt")
avPlots(model3)

# Extract residuals and fitted values 
e<- residuals(model3)
std_e<- rstandard(model3)
std_e_sq<- std_e^2
yhat<- fitted(model3)

##plot std_e_sq vs yhat
plot(yhat, std_e_sq,
     main = "Std Residuals² vs Fitted Values (yhat)",
     xlab = "Fitted Values (yhat)",
     ylab = "Standardised Residuals²",
     pch = 19, col = "darkgrey")

# Add LOESS curve
lines(lowess(yhat, std_e_sq), col = "blue", lwd = 2)

#BP test
library(lmtest)
bptest(model3)
# Goldfeld–Quandt test
gqtest(model3)
gqtest(model3, order.by = ~ ln_gva_industry, data = data_ecotrix)
gqtest(model3, order.by = ~ ln_social_exp, data = data_ecotrix)
gqtest(model3, order.by = ~ ln_power_pc, data = data_ecotrix)
gqtest(model3, order.by = ~ credit_agri_sqrt, data = data_ecotrix)

###NORMALITY

#Shapiro-wilk Test
shapiro.test(data_ecotrix$ln_nsdp_pc)

#Q-Q Plot

qqnorm(data_ecotrix$ln_nsdp_pc, main = "Q-Q Plot: ln_nsdp_pc")
qqline(data_ecotrix$ln_nsdp_pc, col = "red")
x <- data_ecotrix$ln_nsdp_pc

qqnorm(x,
       main = "Q-Q Plot: ln_nsdp_pc (Zoomed Out)",
       xlim = c(-3, 3),        # expand/zoom as needed
       ylim = c(min(x) - 1, max(x) + 1))

qqline(x, col = "red", lwd = 2)

#Anderson Darling Test
install.packages("nortest")
library(nortest)
ad.test(data_ecotrix$ln_nsdp_pc)

###MODEL SELECTION

##Mallows' Cp
install.pmodelAinstall.packages("leaps")
library(leaps)
df <- data_ecotrix[, c("ln_nsdp_pc", 
                       "ln_gva_industry", 
                       "ln_social_exp", 
                       "ln_power_pc", 
                       "credit_agri_sqrt")]
cp_results <- regsubsets(
  ln_nsdp_pc ~ ln_gva_industry + ln_social_exp + ln_power_pc + credit_agri_sqrt,
  data = df,
  nbest = 1,         # best model of each size
  nvmax = 4          # maximum number of predictors
)

summary_cp <- summary(cp_results)
summary_cp$cp
plot(summary_cp$cp, type = "b",
     xlab = "Number of predictors",
     ylab = "Mallows' Cp",
     main = "Mallows' Cp for subset models")
abline(h = 5, col = "red", lty = 2, lwd = 2)

### MULTICOLLINEARITY

vif(model3)
# Select only the variables you want
vars <- data_ecotrix[, c("ln_gva_industry", "ln_social_exp", "ln_power_pc", "credit_agri_sqrt")]

# Correlation matrix
cor_mat <- cor(vars, use = "complete.obs")
cor_mat
# Eigenvalues
eigen(cor_mat)$values

eigen_vals <- eigen(cor_mat)$values
condition_number <- sqrt(max(eigen_vals) / min(eigen_vals))
condition_number

# 1. Extract model matrix
X <- model.matrix(model3)

# 2. Standardize X (excluding intercept)
Z <- scale(X[, -1], center = TRUE, scale = TRUE)

# 3. Compute cross-product matrix
ZtZ <- t(Z) %*% Z
ZtZ
# 4. Eigenvalues
eigen_vals <- eigen(ZtZ)$values

# 5. Condition indices
condition_indices <- sqrt(max(eigen_vals) / eigen_vals)

# 6. Attach variable names
names(condition_indices) <- colnames(Z)

# 7. Print nicely
print(condition_indices)


###Influence analysis

# leverage values (hat-values)
leverages <- hatvalues(model3)

# view first few
head(leverages)

# full summary
summary(leverages)

# number of predictors (excluding intercept)
k <- length(coefficients(model3)) - 1

# sample size
n <- nrow(data_ecotrix)

# leverage cutoff rule-of-thumb
cutoff <- 2*(k+1)/n
cutoff

# identify high-leverage observations
high_lev <- which(leverages > cutoff)
leverages[high_lev]

plot(leverages, 
     pch = 20, 
     main = "Leverage Values",
     ylab = "Leverage (hat values)")

abline(h = cutoff, col = "red", lty = 2)

##outliers

# Studentized residuals
stud_res <- rstudent(model3)

# Identify outliers (cutoff = |2|)
outliers_resid <- which(abs(stud_res) > 2)
outliers_resid

# Plot residuals
plot(stud_res, pch = 19, col = "darkblue")
abline(h = c(-2, 2), col = "red", lty = 2)


##dfbeta
# DFBETAS for each observation and each coefficient
dfb <- dfbetas(model3)

# View first few rows
head(dfb)

# Find observations with |DFBETAS| > 2/sqrt(n)
n <- nrow(data_ecotrix)
cutoff_dfb <- 2 / sqrt(n)

which(apply(abs(dfb), 1, function(x) any(x > cutoff_dfb)))


# DFFITS values
dff <- dffits(model3)


# View first few
head(dff)
dff

# Identify influential observations
cutoff_dff <- 2 * sqrt((length(coef(model3)) - 1) / n)

which(abs(dff) > cutoff_dff)

##cook's D
cook <- cooks.distance(model3)

# rule of thumb: > 4/n is influential
cutoff_cook <- 4 / nrow(data_ecotrix)
which(cook > cutoff_cook)

# plot Cook's distance
plot(cook, type = "h", col = "darkgreen")
abline(h = cutoff_cook, col = "red", lty = 2)

data_ecotrix[c(30, 31), ]
model4 <- lm(ln_nsdp_pc ~ ln_gva_industry + ln_social_exp + 
                     ln_power_pc + credit_agri_sqrt,
                   data = data_ecotrix[-c(30,31), ])
summary(model4)
summary(model3)

# cov ratio
cov_ratio_values <- covratio(model3)

# Print the values
print(cov_ratio_values)

n <- nrow(data_ecotrix)
k <- length(coefficients(model3))
n; k

##cook's D vs residuals
plot(rstudent(model3), 
     cooks.distance(model3),
     xlab = "Studentized Residuals",
     ylab = "Cook's Distance",
     main = "Cook's Distance vs Studentized Residuals",
     pch = 19, col = "blue")

abline(h = 4/(nrow(data_ecotrix)-length(coef(model3))), col = "red", lty = 2)
text(rstudent(model3), cooks.distance(model3), labels = 1:nrow(data_ecotrix), pos = 4, cex = 0.7)

library(car)
influencePlot(model3, 
              main = "Influence Plot for model3")


## cook's D vs hii/(1-hii)

# Leverage (hii)
hii <- hatvalues(model3)

length(hii)
length(cook)
# Transformation h/(1-h)
h_trans <- hii / (1 - hii)

# Plot
plot(h_trans, cook,
     xlab = "hii / (1 - hii)",
     ylab = "Cook's Distance",
     main = "Cook's D vs hii / (1 - hii)",
     pch = 19, col = "darkblue")

# Add cutoff line for Cook's D
abline(h = 4/(nrow(data_ecotrix) - length(coef(model3))), col = "red", lty = 2)

# Label points
text(h_trans, cook, labels = 1:length(hii), pos = 4, cex = 0.7)

# Perform ANOVA (Type I)
anova(model4)


