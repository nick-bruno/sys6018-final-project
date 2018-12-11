# Final Stat Project #
# By: Nick Bruno (nhb3zf) #
library(car)
library(glmnet)
# Load in data
data <- read.csv('stat_train.csv', header=T) # per-game stats (qbs who played more than 10 games)
head(data)

data['games_played'] <- data$win_count/data$team_win # add new column
data$games_played[which(!is.finite(data$games_played))] <- 0

##### Starting linear models #####
### Linear model with all variables that do not have to do with win_count
full_lin_model <- lm(win_count ~ .-name - id - team_win - start_win - games_played, data=data)
summary(full_lin_model)
  # Only significant variable was Fmb, and Sk was almost significant at 10%

fitted_full <- fitted(full_lin_model)
residual_full <- resid(full_lin_model)
plot(fitted_full, residual_full)
abline(h=0, lty=1, lwd=3)
  # Terrible reiduals

### Subjective model 
subj_model <- lm(win_count ~ Y.A + Int + TD, data=data)
summary(subj_model)

fitted_subj <- fitted(subj_model)
residual_subj <- resid(subj_model)

plot(fitted_subj, residual_subj)
abline(h=0, lty=1, lwd=3)
  # slightly better, but still bad residuals
  # Yards per attempt is the only statistically significant model

##### Try variable selection methods to the regular linear model #####
### Elimination suggestions
library(leaps)
resp.name <- "win_count"
reg.names <- c("Att", "Yds", "TD", "Int","Sk","Fmb","Cmp.","Y.A","TD_pct","Int_pct")
data.name <- "data"
alpha.out <- 0.10
alpha.in <- 0.25

backward_elim <- backward.elimination(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Att + Yds + Fmb"
forward_select <- forward.selection(alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Fmb + Y.A"
stepwise_select <- stepwise.selection(alpha.in, alpha.out, resp.name, reg.names, data.name)
  # "Final model: win_count ~ Fmb + Y.A"

### Backward elimination model
backward_lin <- lm(win_count ~ Att + Yds + Fmb, data=data)
summary(backward_lin)
fitted_backward_lin <- fitted(backward_lin)
resids_backward_lin <- resid(backward_lin)
plot(fitted_backward_lin, resids_backward_lin)
abline(h=0, lty=1, lwd=3) # bad resids


### Forward and Stepwise selection model
forward_lin <- lm(win_count ~ Fmb + Y.A, data=data)
summary(forward_lin)
fitted_forward_lin <- fitted(forward_lin)
resids_forward_lin <- resid(forward_lin)
plot(fitted_forward_lin, resids_forward_lin)
abline(h=0, lty=1, lwd=3) # bad resids


##### Logarithmic functions #####
  # Does not work because two observations had 0 wins
rows_of_zero_wins <- as.numeric(rownames(data[which(data$win_count == 0), ]))
new_data <- data[-rows_of_zero_wins, ]
nrow(new_data)
full_log_lin_model <- lm(log(win_count) ~ .-name - id - team_win - start_win - games_played, data=new_data)
summary(full_log_lin_model)
  # Statistically significant variables at 5%: Sk, Fmb, 
  # Statistically significant variables at 10%: Cmp, Yds, Cmp.

# Add log(wins) to the dataframe
new_data['log_wins'] <- log(new_data$win_count)

### Elimination log models
resp.name <- "log_wins"
reg.names <- c("Att", "Yds", "TD", "Int","Sk","Fmb","Cmp.","Y.A","TD_pct","Int_pct")
data.name <- "new_data"
alpha.out <- 0.10
alpha.in <- 0.25
backward_elim_log <- backward.elimination(alpha.out, resp.name, reg.names, data.name)
  # Final model: log_wins ~ Att + Yds + Sk + Fmb"
forward_select_log <- forward.selection(alpha.out, resp.name, reg.names, data.name)
  # "Final model: log_wins ~ Sk + Fmb + Y.A"
stepwise_select <- stepwise.selection(alpha.in, alpha.out, resp.name, reg.names, data.name)
  # "Final model: log_wins ~ Sk + Fmb + Y.A"

### Backward elimination model on log(wins)
backward_lin_log <- lm(log_wins ~ Att + Yds + Sk + Fmb, data=new_data)
summary(backward_lin_log)
fitted_back_log <- fitted(backward_lin_log)
resid_back_log <- resid(backward_lin_log)
plot(fitted_back_log, resid_back_log)
abline(h=0, lty=1, lwd=3)
  # RESIDS LOOK GREAT
vif(backward_lin_log)
# Att      Yds       Sk      Fmb 
# 3.209680 2.972284 1.402497 1.069391 
  # SEEMS OK

### Forward and stepwise selection model on log(wins)
forward_lin_log <- lm(log_wins ~ Sk + Fmb + Y.A, data=new_data)
summary(forward_lin_log)
fitted_forward_log <- fitted(forward_lin_log)
resids_forward_log <- resid(forward_lin_log)
plot(fitted_forward_log, resids_forward_log)
abline(h=0, lty=1, lwd=3)
  # Resids still look good

## Conclusion: Still low R-squared values, but residuals are much better with a logarithmic transformation 
## However, f-statistic of each regression is large and the p-value is significant. ##

### Log functions with quadratic terms
quad_forward <- lm(log_wins ~ poly(Sk,2) + poly(Fmb,2) + poly(Y.A, 2), data=new_data)
summary(quad_forward)
  # Increases adj. R-squared but decreasese the F-statistic
quad_forward_interaction <- lm(log_wins ~ poly(Sk,2) + poly(Fmb,2) + poly(Y.A, 2) + Sk*Fmb + Sk*Y.A + Fmb*Y.A, data=new_data)
summary(quad_forward_interaction)
  # Interaction terms do not help

##### Ridge Regression #####
# All regressors model 
new_data 
y.vect <- new_data$win_count
X0.mat <- as.matrix(new_data[,c(3:9, 12:15)])

ridge_all_regressors <- glmnet(X0.mat, y.vect, alpha=0)
summary(ridge_all_regressors)
plot(ridge_all_regressors, xvar="lambda", label=TRUE) # lambda around 1

ridge_all_regressors_fix <- glmnet(X0.mat, y.vect, alpha=0, lambda=exp(1))
coefficients(ridge_all_regressors_fix)

ridge_all_regressors_fix$dev.ratio # 0.08412482 (better than the linear version)


# Forward selection regressors applied to the ridge regression
# Fmb + Y.A
X0.mat.new <- as.matrix(new_data[,c(9,13)])
ridge_forward_selection <- glmnet(X0.mat.new, y.vect, alpha=0)
summary(ridge_forward_selection)
plot(ridge_forward_selection, xvar="lambda", label=TRUE) # lambda=0.75

ridge_forward_selection_fix <- glmnet(X0.mat.new, y.vect, alpha=0, lambda=exp(0.75))
coefficients(ridge_forward_selection_fix)
  #                     s0
  # (Intercept)  -7.075024
  # Fmb         -30.909343
  # Y.A           6.506742

ridge_forward_selection_fix$dev.ratio # 0.06642173 (not very good)
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

