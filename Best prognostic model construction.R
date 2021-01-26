# Multivariate Cox regression with three penalities including least absolute shrinkage and selection operator(Lasso), Adaptive lasso and Elastic net algorithms for informative prognostic-related genes selection.
# we need two Package in R 
library(glmnet)
library(survival)
dat<-read.csv("working directroy", row.names = 1)
set.seed(1234)
# Prognostic gene with Elastic net algorithm selection 
y <- Surv(dat$os, dat$event) # os-overall survival, event-overall survival status
x<-model.matrix(y~., dat[,c(-1:-2)])
cv.fit <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 0.5)
plot(cv.fit)
fit <- glmnet(x,y, family = "cox", alpha = 1)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
coef(cv.fit, s = "lambda.min") # extracting selecting prognostic genes
# Prognostic gene with Lasso algorithm selection 
y <- Surv(dat$os, dat$event) 
x<-model.matrix(y~., dat[,c(-1:-2)])
cv.fit <- cv.glmnet(x,y, family="cox",nfold = 10, alpha= 1)
plot(cv.fit)
fit <- glmnet(x,y, family = "cox", alpha = 1)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
coef(cv.fit, s = "lambda.min") # extracting selecting prognostic genes

#addaptive lasso
cv.fit <- cv.glmnet(x,y, family="cox", nfold = 10, alpha= 0) # rigid Cox Model
plot(cv.fit)
fit <- glmnet(x,y, family = "cox",
              nfold = 10, alpha= 0)
plot(fit)
cv.fit$lambda.min

Coefficients <- coef(fit, s = fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
coef(cv.fit, s = "lambda.min")
# addaptive lasso
## Extract coefficients at the error-minimizing lambda
cv.fit$lambda.min
coef(cv.fit, s = "lambda.min")
## The intercept estimate should be dropped.
best_ridge_coef <- as.numeric(coef(cv.fit, s = "lambda.min"))
best_ridge_coef
## Perform adaptive LASSO
alasso1 <- glmnet(x,y, family = "cox", alpha = 1,
                  penalty.factor = 1 / abs(best_ridge_coef))
plot(alasso1, xvar = "lambda")

## Perform adaptive LASSO with 10-fold CV
alasso1_cv <- cv.glmnet(x,y, family = "cox",
                        ## type.measure: loss to use for cross-validation.
                        ## K = 10 is the default.
                        nfold = 10,
                        ## 'alpha = 1' is the lasso penalty, and 'alpha = 0' the ridge penalty.
                        alpha = 1,
                        penalty.factor = 1 / abs(best_ridge_coef),
                        ## prevalidated array is returned
                        keep = TRUE)
## Penalty vs CV MSE plot
plot(alasso1_cv)

## Extract coefficients at the error-minimizing lambda
alasso1_cv$lambda.min
## s: Value(s) of the penalty parameter 'lambda' at which
##    predictions are required. Default is the entire sequence used
##    to create the model.
best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
best_alasso_coef1 # extracting prognostic gene selected by Adaptive lasso 
# Best subset Cox regression model construction 
# Combine all gene selected by three algorithms and use glumti package to identify which combination of genes would produce optimal prognostic signture 

#best subset regression analysis 
library(survival)
library(glmulti)
dat2<-read.csv("working directory path", row.names = 1)
dat <- within(dat, {
  survival.vector    <- Surv(dat$os, dat$event)})
#*************************************************************
glmulti.coxph.out <-glmulti(survival.vector ~., data = dat,
                            level = 1,               # No interaction considered
                            method = "h",            # Exhaustive approach
                            crit = "aic",            # AIC as criteria
                            confsetsize = 5,         # Keep 5 best models
                            plotty = F, report = F,  # No plot or interim reports
                            fitfunction = "coxph")   # coxph function

## Show result for the best model
summary(glmulti.coxph.out@objects[[1]])
## Show 5 best models (Use @ instead of $ for an S4 object)
glmulti.coxph.out@formulas
plot(res)
plot(res, type="s")
summary(res@objects[[1]])
print(res)
top <- weightable(res)
top

# Risk score model construction based Best subset  prognostic genes 
model<- coxph(Surv(os,event) ~. , data =dat)
summary(model)
# 
dat<-read.csv("wd path")
fit<- survfit(Surv(dat$os,dat$event) ~dat$RS, data = dat)# Rs-Risk score
fit<- survfit(Surv(dat$os,dat$event) ~dat$RS, data = dat)
fit
ggsurvplot(fit, data = dat, pval = T,
           legend.labs = c("High-risk", "Low-risk"),
           risk.table = TRUE,
           xlab = "time(months)", 
           ylab = "survival probablity",
           legend.title = " Risk",
           risk.table.y.text.col = T,
           risk.table.y.text = T)
# Prepare dataset having OS time and status and Risk socre each HCC patients
# dat<-read.csv("wd path", row.names = 1)
library(timeROC)
library(survival)
ROC.DSST<-timeROC(T=dat$os,delta=dat$event,
                  marker=dat$RS,cause=1,
                  weighting="cox",
                  times=c(12,24,36),ROC=TRUE) # 1st year, 2nd year and 3rd year risk classification performance 

plot(ROC.DSST,time=12,title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col="dodgerblue4",add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=60,col="chartreuse",add=TRUE,title=FALSE,lwd=2)
legend("bottomright",text.font = 2, c("1-year(AUC = o.865 )", "2-year (AUC = 0.854)", "3-year (AUC = 0.779)"),col=c("red","dodgerblue4", "Green" "),lwd=2, cex = 0.75)
#******************************************************

