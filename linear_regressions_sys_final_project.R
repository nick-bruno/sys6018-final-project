# Final Sys Project #

# Upload libraries
library(car)
library(glmnet)
library(lubridate)
library(tidyverse)
library(leaps)

# Load in data
setwd('/Users/nickbruno/STAT6021')
data <- read.csv('stat_train.csv', header=T) # per-game stats (qbs who played more than 10 games)
head(data)
names(data)

##### Variable Investigation #####
summary(data$win_count)
qqnorm(data$win_count)
qqline(data$win_count)
  # qqplot shows a non-linear distribution in the win_count variable


logs <- log(data$win_count)
qqnorm(logs)
qqline(logs)
  # Taking the logs of the response variable "win_count" leads to a much more normal looking
  # qqplot

# Take a further look at "win_count"
summary(data$win_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    8.00   19.00   34.14   54.00  199.00 
boxplot(data$win_count, main='Box Plot of Wins per Quarterback')
  # Very uneven distribution

##### Starting linear models #####
### Linear model with all variables that do not have to do with win_count
full_lin_model <- lm(win_count ~ .-name - id - team_win - start_win, data=data)
summary(full_lin_model)
  # Only significant variable was Fmb, and Sk was almost significant at 10%
vif(full_lin_model) # serious issues
# Cmp        Att        Yds         TD        Int         Sk        Fmb       Cmp.        Y.A 
# 178.589916 100.772584 164.524789  47.121116  33.371322   1.669225   1.129366  41.888207  58.649848 
# TD_pct    Int_pct 
# 45.989280  38.203976 
  # vif values larger than 10 are problematic, so multicollinearity is a huge problem in this initial 
  # linear model


# Looking at influential points
hat_vals <- hatvalues(full_lin_model)
res <- resid(full_lin_model)
PRESS.res <- res / (1-hat_vals)
press_stat <- sum((PRESS.res)^2) # 362762.6

aic_lin <- AIC(full_lin_model) # 3079.781
bic_lin <- BIC(full_lin_model)
adj_lin <- summary(full_lin_model)$adj.r.squared # 0.06755128
mse_lin <- anova(full_lin_model)[12,3] # 1110.885
  # We will use these statistics as benchmarks to compare our future models

# Find the quantity of outliers
length(which(abs(rstandard(full_lin_model)) > 3)) # 4 observations
length(which(abs(rstandard(full_lin_model)) > 2)) # 10 observations

# Analyze the outliers
outliers_lin <- data[c(33,43,68,91,95,103,163,173,228,282),]
mean(outliers_lin$win_count) # 138.5
  # The outliers are caused by large win_counts, averaging a win_count of 138.5 from the 10 outlier
  # observations, much larger than the mean of quarterback wins of 34.14

# Examine residual plot
fitted_full <- fitted(full_lin_model)
residual_full <- resid(full_lin_model)
plot(fitted_full, residual_full)
abline(h=0, lty=1, lwd=3)
  # Terrible reiduals, obviously not a normal distribution

### Subjective model (chose regressors we felt would help predict win_count)
subj_model <- lm(win_count ~ Y.A + Int + TD, data=data)
summary(subj_model)
  # Only Y.A. is significant

# Plot residuals
fitted_subj <- fitted(subj_model)
residual_subj <- resid(subj_model)

plot(fitted_subj, residual_subj)
abline(h=0, lty=1, lwd=3)
  # slightly better, but still bad residuals

##### Try variable selection methods to the regular linear model #####
### Elimination suggestions

# Define variable selection functions (created by Professor Spitzner)
forward.step <- function(curr.var.in, alpha.in, resp.name, reg.names, data.name) {
  curr.var.out.idx <- which(!curr.var.in)
  enter.idx <- NA
  if (length(curr.var.out.idx) > 0) {
    k <- length(reg.names)
    pval.seq <- rep(x=Inf, times=k)
    for (iVAR in curr.var.out.idx) {
      cand.var.in <- curr.var.in
      cand.var.in[iVAR] <- TRUE
      cand.model.str <- get.model.str(cand.var.in, resp.name, reg.names)
      cand.model.lm <- eval.lm(cand.model.str, data.name)
      iROW <- which(row.names(summary(cand.model.lm)$coefficients) == reg.names[iVAR])
      pval.seq[iVAR] <- summary(cand.model.lm)$coefficients[iROW,4]
    }
    enter.idx <- which.min(pval.seq)
    if (pval.seq[enter.idx] < alpha.in) {
      print(paste("Variable ", reg.names[enter.idx], " enters the model (pval=", sprintf("%6.4f", pval.seq[enter.idx]), ")", sep=""))
    } else {
      print("No variables enter the model")
      enter.idx <- NA
    }
  } else {
    print("No variables available to enter the model")
  }
  return(enter.idx)
}

backward.step <- function(curr.var.in, alpha.out, resp.name, reg.names, data.name) {
  curr.var.in.idx <- which(curr.var.in)
  leave.idx <- NA
  if (length(curr.var.in.idx) > 0) {
    k <- length(reg.names)
    pval.seq <- rep(x=-Inf, times=k)
    curr.model.str <- get.model.str(curr.var.in, resp.name, reg.names)
    curr.model.lm <- eval.lm(curr.model.str, data.name)
    for (iVAR in curr.var.in.idx) {
      iROW <- which(row.names(summary(curr.model.lm)$coefficients) == reg.names[iVAR])
      pval.seq[iVAR] <- summary(curr.model.lm)$coefficients[iROW,4]
    }
    leave.idx <- which.max(pval.seq)
    if (pval.seq[leave.idx] >= alpha.out) {
      print(paste("Variable ", reg.names[leave.idx], " leaves the model (pval=", sprintf("%6.4f", pval.seq[leave.idx]), ")", sep=""))
    } else {
      print("No variables leave the model")
      leave.idx <- NA
    }
  } else {
    print("No variables available to leave the model")
  }
  return(leave.idx)
}

forward.selection <- function(alpha.in, resp.name, reg.names, data.name) {
  k <- length(reg.names)
  curr.var.in <- rep(x=FALSE, times=k)
  stop <- FALSE
  while(!stop) {
    enter.idx <- forward.step(curr.var.in, alpha.in, resp.name, reg.names, data.name)
    if (is.na(enter.idx)) {
      stop <- TRUE
    } else {
      curr.var.in[enter.idx] <- TRUE
    }
  }
  curr.model.str <- get.model.str(curr.var.in, resp.name, reg.names)
  print(paste("Final model: ", curr.model.str, sep=""))
  curr.model.lm <- eval.lm(curr.model.str, data.name)
  return(curr.model.lm)
}

backward.elimination <- function(alpha.out, resp.name, reg.names, data.name) {
  k <- length(reg.names)
  curr.var.in <- rep(x=TRUE, times=k)
  stop <- FALSE
  while(!stop) {
    leave.idx <- backward.step(curr.var.in, alpha.out, resp.name, reg.names, data.name)
    if (is.na(leave.idx)) {
      stop <- TRUE
    } else {
      curr.var.in[leave.idx] <- FALSE
    }
  }
  curr.model.str <- get.model.str(curr.var.in, resp.name, reg.names)
  print(paste("Final model: ", curr.model.str, sep=""))
  curr.model.lm <- eval.lm(curr.model.str, data.name)
  return(curr.model.lm)
}

stepwise.selection <- function(alpha.in, alpha.out, resp.name, reg.names, data.name) {
  k <- length(reg.names)
  curr.var.in <- rep(x=FALSE, times=k)
  stop <- FALSE
  while(!stop) {
    enter.idx <- forward.step(curr.var.in, alpha.in, resp.name, reg.names, data.name)
    if (is.na(enter.idx)) {
      stop <- TRUE
    } else {
      curr.var.in[enter.idx] <- TRUE
      leave.idx <- backward.step(curr.var.in, alpha.out, resp.name, reg.names, data.name)
      if (!is.na(leave.idx)) {
        curr.var.in[leave.idx] <- FALSE
        if (leave.idx == enter.idx) {
          stop <- TRUE
        }
      }
    }
  }
  curr.model.str <- get.model.str(curr.var.in, resp.name, reg.names)
  print(paste("Final model: ", curr.model.str, sep=""))
  curr.model.lm <- eval.lm(curr.model.str, data.name)
  return(curr.model.lm)
}

# Use these functions to find ideal models
names(data)
resp.name <- "win_count"
reg.names <- c("Att","Cmp", "Yds", "TD", "Int","Sk","Fmb","Cmp.","Y.A","TD_pct","Int_pct")
data.name <- "data"
alpha.out <- 0.10
alpha.in <- 0.25

backward_elim <- backward.elimination(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Cmp + Yds + Fmb + Cmp."
forward_select <- forward.selection(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Fmb + Y.A"
stepwise_select <- stepwise.selection(alpha.in, alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Fmb + Y.A"

##### Examining the suggest models #####
### Backward elimination model
backward_lin <- lm(win_count ~ Cmp + Yds + Fmb + Cmp., data=data)
summary(backward_lin) # adj. r-squared of 0.06942
vif(backward_lin)
  # Cmp      Yds      Fmb     Cmp. 
  # 6.433951 4.269884 1.094783 2.436830 
  # Much better in terms of multicollinearity

# Plot residuals
fitted_backward_lin <- fitted(backward_lin)
resids_backward_lin <- resid(backward_lin)
plot(fitted_backward_lin, resids_backward_lin, main='Linear Backward Elimination Model Residuals')
abline(h=0, lty=1, lwd=3)
  # Still very poor residuals

hat_vals_back_lin <- hatvalues(backward_lin)
res_back_lin <- resid(backward_lin)
PRESS.res_back_lin <- res_back_lin / (1-hat_vals_back_lin)
press_stat_back_lin <- sum((PRESS.res_back_lin)^2) # 351159.7

aic_back_lin <- AIC(backward_lin) # 3079.781
bic_back_lin <- BIC(backward_lin)
adj_back_lin <- summary(backward_lin)$adj.r.squared # 0.06941528
mse_back_lin <- anova(backward_lin)[5,3]
  # These statistics seem to improve, but not greatly.

### Forward and Stepwise selection model
forward_lin <- lm(win_count ~ Fmb + Y.A, data=data)
summary(forward_lin) # adjusted r-squared of 0.06523

# Analyze residuals
fitted_forward_lin <- fitted(forward_lin)
resids_forward_lin <- resid(forward_lin)
plot(fitted_forward_lin, resids_forward_lin, main='Linear Forward and Stepwise Selection Model Residuals')
abline(h=0, lty=1, lwd=3) # poor residuals

hat_vals_forward_lin <- hatvalues(forward_lin)
res_forward_lin <- resid(forward_lin)
PRESS.res_forward_lin <- res_forward_lin / (1-hat_vals_forward_lin)
press_stat_forward_lin <- sum((PRESS.res_forward_lin)^2) # 348975.2

aic_forward_lin <- AIC(forward_lin) # 3079.209
bic_forward_lin <- BIC(forward_lin)
adj_forward_lin <- summary(forward_lin)$adj.r.squared # 0.06522607
mse_forward_lin <- anova(forward_lin)[3,3] # 1113.655
  # Still poor statistics

##### Logarithmic functions #####
  # Does not initially work because two observations had 0 wins, so I will remove these two observations
rows_of_zero_wins <- as.numeric(rownames(data[which(data$win_count == 0), ]))
new_data <- data[-rows_of_zero_wins, ]

data <- new_data # reassigns the new_data to the 'data' title we have been using throughout this analysis

# Analyze a logarithmic regression with each variable
full_log_lin_model <- lm(log(win_count) ~ .-name - id - team_win - start_win, data=data)
summary(full_log_lin_model)
  # Statistically significant variables at 5%: Sk, Fmb, 
  # Statistically significant variables at 10%: Cmp, Yds, Cmp.

# Add log(wins) to the dataframe
data['log_wins'] <- log(data$win_count)
qqnorm(data$log_wins)
qqline(data$log_wins)
  # This qqplot is much more normal

### Elimination log models (including Cmp)
resp.name <- "log_wins"
reg.names <- c("Att", "Cmp", "Yds", "TD", "Int","Sk","Fmb","Cmp.","Y.A","TD_pct","Int_pct")
data.name <- "data"
alpha.out <- 0.10
alpha.in <- 0.25
backward_elim_log <- backward.elimination(alpha.out, resp.name, reg.names, data.name)
  # Final model: Cmp + Yds + Int + Sk + Fmb + Cmp. + Y.A"
forward_select_log <- forward.selection(alpha.out, resp.name, reg.names, data.name)
  # "Final model: log_wins ~ Sk + Fmb + Y.A"
stepwise_select <- stepwise.selection(alpha.in, alpha.out, resp.name, reg.names, data.name)
  # "Final model: log_wins ~ Sk + Fmb + Y.A"

### Backward elimination model on log(wins)
backward_lin_log <- lm(log_wins ~ Cmp + Yds + Int + Sk + Fmb + Cmp., data=data)
summary(backward_lin_log) # adj. r-squared = 0.1154, huge improvement

# Analyze residual plots
fitted_back_log <- fitted(backward_lin_log)
resid_back_log <- resid(backward_lin_log)
plot(fitted_back_log, resid_back_log, main='Logarithmic Regression Backward Elimination Model Residuals')
abline(h=0, lty=1, lwd=3)
  # RESIDS LOOK GREAT
back_subset <- new_data[,c(3,5,7:9,12)]
cor(back_subset)
# Cmp      Yds      Int       Sk      Fmb     Cmp. 
# 7.130904 4.373840 1.110906 1.550481 1.100481 2.470663
cor(new_data$Cmp, new_data$Cmp.) # 0.7616852, not great

studentized_resids <- rstandard(backward_lin_log)
major_outlier_rows <- which(abs(studentized_resids) > 3) # only one observation
new_data[21,]
outlier_rows <- which(abs(studentized_resids) > 2)
outliers <- new_data[c(21,33,54,66,91,209,228,239,276,292,304),] # outliers
mean(outliers$win_count)

hat_vals_back <- hatvalues(backward_lin_log)
res_back <- resid(backward_lin_log)
PRESS.res_back <- res_back / (1-hat_vals_back)
press_stat_back <- sum((PRESS.res_back)^2) # 397.779

aic_back_log <- AIC(backward_lin_log) # 958.3052
bic_back_log <- BIC(backward_lin_log)
adj_back_log <- summary(backward_lin_log)$adj.r.squared 
mse_back_log <- anova(backward_lin_log)[7,3]

summary(backward_lin_log)$coefficients
# Estimate  Std. Error   t value     Pr(>|t|)
# (Intercept)  1.22924206 0.541935785  2.268243 2.401797e-02
# Cmp         -0.17875896 0.046189218 -3.870145 1.332269e-04
# Yds          0.01686014 0.003418863  4.931506 1.347134e-06
# Int         -0.19465137 0.137948055 -1.411048 1.592563e-01
# Sk          -0.12665682 0.058178523 -2.177037 3.024939e-02
# Fmb         -1.09879792 0.351821688 -3.123167 1.962085e-03
# Cmp.         3.36111803 1.300834314  2.583817 1.023945e-02

confint(backward_lin_log, level=0.95)
# 2.5 %      97.5 %
#   (Intercept)  0.16280776  2.29567636
# Cmp         -0.26965121 -0.08786670
# Yds          0.01013242  0.02358787
# Int         -0.46610888  0.07680613
# Sk          -0.24114192 -0.01217172
# Fmb         -1.79112111 -0.40647472
# Cmp.         0.80130492  5.92093114



### Forward and stepwise selection model on log(wins)
forward_lin_log <- lm(log_wins ~ Sk + Fmb + Y.A, data=new_data)
summary(forward_lin_log)

# Analyze residual plots
fitted_forward_log <- fitted(forward_lin_log)
resids_forward_log <- resid(forward_lin_log)
plot(fitted_forward_log, resids_forward_log, main='Logarithmic Regression Forward Selection Model Residuals')
abline(h=0, lty=1, lwd=3)
  # Resids still look good

hat_vals_forward_log <- hatvalues(forward_lin_log)
res_forward_log <- resid(forward_lin_log)
PRESS.res_forward_log <- res_forward_log / (1-hat_vals_forward_log)
press_stat_forward_log <- sum((PRESS.res_forward_log)^2) # 397.779

aic_forward_log <- AIC(forward_lin_log) # 958.3052
bic_forward_log <- BIC(forward_lin_log)
adj_forward_log <- summary(forward_lin_log)$adj.r.squared 
mse_forward_log <- anova(forward_lin_log)[4,3]


## Conclusion: Still low R-squared values, but residuals are much better with a logarithmic transformation 
## However, f-statistic of each regression is large and the p-value is significant. ##

### Log functions with quadratic terms
quad_forward <- lm(log_wins ~ poly(Sk,2) + poly(Fmb,2) + poly(Y.A, 2), data=new_data)
summary(quad_forward)
  # Increases adj. R-squared but decreasese the F-statistic, so we will not use this in our analysis
quad_forward_interaction <- lm(log_wins ~ poly(Sk,2) + poly(Fmb,2) + poly(Y.A, 2) + Sk*Fmb + Sk*Y.A + Fmb*Y.A, data=new_data)
summary(quad_forward_interaction)
  # Interaction terms do not help

##### Ridge Regression #####
# All regressors model 
new_data 
y.vect <- new_data$log_wins
X0.mat <- as.matrix(new_data[,c(3:9, 12:15)])

ridge_all_regressors <- glmnet(X0.mat, y.vect, alpha=0)
summary(ridge_all_regressors)
plot(ridge_all_regressors, xvar="lambda", label=TRUE) # lambda around -1.75

ridge_all_regressors_fix <- glmnet(X0.mat, y.vect, alpha=0, lambda=exp(-1.75))
coefficients(ridge_all_regressors_fix)
# s0
# (Intercept)  1.704998557
# Cmp         -0.020147471
# Att         -0.012740333
# Yds          0.003880808
# TD           0.179671221
# Int         -0.043091488
# Sk          -0.103711274
# Fmb         -0.982848581
# Cmp.         0.124643714
# Y.A          0.212130911
# TD_pct      -2.460937787
# Int_pct     -2.190829206

ridge_all_regressors_fix$dev.ratio # 0.1217926 (better than the linear version)

# Backward elimination regressors applied to the ridge regression
y.vect <- new_data$log_wins
X0.mat.new <- as.matrix(new_data[,c(3,5, 7:9,12)])
ridge_back <- glmnet(X0.mat.new, y.vect, alpha=0)
summary(ridge_back)
plot(ridge_back, xvar="lambda", label=TRUE) # lambda=-2

ridge_back_fix <- glmnet(X0.mat.new, y.vect, alpha=0, lambda=exp(-2))
coefficients(ridge_back_fix)
#                 s0
# (Intercept)  1.888349425
# Cmp         -0.072119668
# Yds          0.008959173
# Int         -0.139149522
# Sk          -0.129215899
# Fmb         -1.046139197
# Cmp.         1.787391504

ridge_back_fix$dev.ratio # 0.1139675 (not very good)
  # compared to 0.06523 from the linear model with the same variables


# Logistic Ridge Regression
# Sk + Fmb + Y.A
X0.log.mat <- as.matrix(new_data[,c(8:9,13)])
y.vect <- new_data$log_wins

ridge_log_reg <- glmnet(X0.log.mat, y.vect, alpha=0)
plot(ridge_log_reg, xvar="lambda", label=TRUE) # lambda=0.75

ridge_log_reg_fix <- glmnet(X0.log.mat, y.vect, alpha=0, lambda=exp(0.75))
coefficients(ridge_log_reg_fix)

ridge_log_reg_fix$dev.ratio #  0.06653291 (not very good)

ridge_log_reg_fix_lower_lambda <- glmnet(X0.log.mat, y.vect, alpha=0, lambda=exp(-1.5))
coefficients(ridge_log_reg_fix_lower_lambda)
ridge_log_reg_fix_lower_lambda$dev.ratio # 0.1091631 (MUCH BETTER THAN THE RIDGE ABOVE)


##### Regressions with win percentage as the response variable
# Upload new dataset
setwd('/Users/nickbruno/STAT6021')
stat_data <- read.csv('new_stat_train.csv', header=T)
max(stat_data$win_percent)
names(stat_data)

qqnorm(stat_data$win_percent)
qqline(stat_data$win_percent)

resp.name <- "win_percent"
reg.names <- c("Att", "Yds", "TD", "Int","Sk","Fmb","Cmp.","Y.A","TD_pct", "Int_pct")
# "team_win"
data.name <- "stat_data"
alpha.out <- 0.10
alpha.in <- 0.25

backward_new_reg <- backward.elimination(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_percent ~ Att + Yds + Sk + Fmb + Int_pct"

forward_new_reg <- forward.selection(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_percent ~ Int + Sk + Fmb + Y.A"

stepwise_new_reg <- stepwise.selection(alpha.in, alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_percent ~ Int + Sk + Fmb + Y.A"

# Backward linear model
backward_lin_percent <- lm(win_percent ~ Att + Yds + Sk + Fmb + Int_pct, data=stat_data)
summary(backward_lin_percent) # adj. r-squared = 0.1421

fitted_backward_lin_percent <- fitted(backward_lin_percent)
resids_backward_lin_percent <- resid(backward_lin_percent)
plot(fitted_backward_lin_percent, resids_backward_lin_percent)
abline(h=0, lty=1, lwd=3)
  # Residuals look good
vif(backward_lin_percent)
# Att      Yds       Sk      Fmb  Int_pct 
# 3.293566 2.967618 1.524564 1.086949 1.298948 

# Forward and stepwise linear model
forward_lin_percent <- lm(win_percent ~ Y.A + Sk + Fmb + Int, data=stat_data)
summary(forward_lin_percent) # adj. r-squared = 0.1323

# Plot residuals
fitted_forward_lin_percent <- fitted(forward_lin_percent)
resids_forward_lin_percent <- resid(forward_lin_percent)
plot(fitted_forward_lin_percent, resids_forward_lin_percent)
abline(h=0, lty=1, lwd=3)
  # Also looks fine, nonconstant variance

vif(forward_lin_percent)
# Y.A       Sk      Fmb      Int 
# 1.004220 1.098664 1.057258 1.062665 
  # No multicollinearity present


##### Summary Statistics for 5 main models #####
adj.r.squared <- c(adj_lin,adj_back_lin,adj_forward_lin,adj_back_log,adj_forward_log)
press.stats <- c(press_stat, press_stat_back_lin, press_stat_forward_lin,press_stat_back ,press_stat_forward_log)
mse <- c(mse_lin, mse_back_lin, mse_forward_lin, mse_back_log, mse_forward_log)
aic <- c(aic_lin, aic_back_lin, aic_forward_lin, aic_back_log, aic_forward_log)
bic <- c(bic_lin, bic_back_lin, bic_forward_lin, bic_back_log, bic_forward_log)

# Matrix will be 5x5
summary_matrix <- rbind(adj.r.squared, press.stats, mse, aic, bic)
colnames(summary_matrix) <- c('Full Lin Mod','Backwards Lin Mod','Forward Lin Mod',
                           'Backwards Log Mod','Forward Log Mod')
new_sum_matrix <- cbind(adj.r.squared, press.stats, mse, aic, bic)
rownames(new_sum_matrix) <- c('Full Lin Mod','Backwards Lin Mod','Forward Lin Mod',
                              'Backwards Log Mod','Forward Log Mod')

##### Split the summary tables between logarithmic response and normal win_count response variable #####
adj.r.squared.lin <- c(adj_lin, adj_back_lin, adj_forward_lin)
press.stat.lin <- c(press_stat, press_stat_back_lin, press_stat_forward_lin)
mse.lin <- c(mse_lin, mse_back_lin, mse_forward_lin)
aic.lin <- c(aic_lin, aic_back_lin, aic_forward_lin)
bic.lin <- c(bic_lin, bic_back_lin, bic_forward_lin)


lin_sum_matrix <- cbind(adj.r.squared.lin, press.stat.lin, mse.lin, aic.lin, bic.lin)
rownames(lin_sum_matrix) <- c('Full Linear Mod','Backwards Linear Mod','Forward Linear Mod')
colnames(lin_sum_matrix) <- c('Adj.r.squared','press.statistic','mse','aic','bic')

adj.r.squared.log <- c(adj_back_log,adj_forward_log)
press.stat.log <- c(press_stat_back ,press_stat_forward_log)
mse.log <- c(mse_back_log, mse_forward_log)
aic.log <- c(aic_back_log, aic_forward_log)
bic.log <- c(bic_back_log, bic_forward_log)

log_sum_matrix <- cbind(adj.r.squared.log, press.stat.log, mse.log, aic.log, bic.log)
rownames(log_sum_matrix) <- c('Backwards Logarithmic Mod','Foward Linear Mod')
colnames(log_sum_matrix) <- c('Adj.r.squared','press.statistic','mse','aic','bic')
