# mix-effect logistic regression for prediction 
# data should be a data frame, the rownames are the sample names, and the colnames are the feature names, 
# The data should contains a "target variable", which is the classification, "0" for control, "1" for case
# which should be noticed is, the matrix is the raw value, not the scale one, if you want to scale 
# if something is wrong when random_effects_var is 'factor', change it to 'character'


# 'data' should be a data frame, the rownames are the sample names, and the colnames are the feature names
# 'target_var' is the name of classification variable in the 'data'
# 'num_iterations' is the number of iteration that how many times you want the model to repeat
# 'random_effects_var' is the name of 'random effect', at most of time, it should be the 'batch'. For example, you merge data from two cohorts, than you should adjust the batch effect, a variable with information of 'cohort1' or 'cohort2' should be assign here
# 'family' depends on your data type. For binary prediction, you should use "binomial". More information in 'lme4::glmer()'

library(lme4)
library(pROC)
library(Matrix)
library(glmnet)
# *****************************************
# ******** Use whole data to train *******
# *****************************************
lasso_logistic_regression_for_select_noncolinear_predictor <- function(data, 
                                                                       target_var, 
                                                                       scale_or_not = TRUE,
                                                                       num_iterations = 100, 
                                                                       family = "binomial",
                                                                       typemeasure = "class",
                                                                       nfolds = 5,
                                                                       custom_lambda = 10^seq(3, -3, by = -0.1)    # Define your own sequence of lambda value (not really necessary)
                                                                       ){

  results <- data.frame()
  DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(data)))
  cv_fit_list = list()

  for (k in 1:num_iterations) {
    set.seed(k)
    #target_var <- "severity"
    predictors <- setdiff(names(data), c(target_var))

    # Separate predictors and response
    X <- model.matrix(as.formula(paste(target_var, "~ .")), data = data)[, -1]
    if (scale_or_not == TRUE) {
      X <- apply(X, 2, scale)
    } else {
      X <- X
    }

    y <- data[[target_var]]

    # Fit the logistic regression model with Lasso penalty using cross-validation
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = nfolds, lambda = custom_lambda, type.measure= "mse")
    cv_fit_list[[k]] <- cv_fit
    # Get the best lambda
    best_lambda <- cv_fit$lambda.min
    # Fit the final model with the best lambda
    fit <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)


    # Store coefficients
    coef_values <- coef(fit)[,]
    coef_values <- coef_values[coef_values != 0]
    coef_names <- names(coef_values)
    coef_names <- c(coef_names, rep(NA, nrow(DMR_predictor) - length(coef_names)))
    coef_values_attach <- c(coef_values, rep(NA, nrow(DMR_predictor) - length(coef_values)))
    iteration_coef <- data.frame(coef_names, coef_values_attach)
    DMR_predictor <- cbind(DMR_predictor, iteration_coef)

    # Prediction and ROC analysis
    train_pre <- predict(fit, newx = X, type = 'response', re.form = NA)
    train_modelroc <- roc(y, as.numeric(train_pre)) #计算ROC绘图数据
    roc_train = train_modelroc$auc
    accuracy_train <- sum(diag(table(y, ifelse(train_pre < 0.5, 0, 1)))) / length(y)


    result <- data.frame(accuracy_train, roc_train)
    results <- rbind(results, result)

    print(k)
    print(cv_fit)
    print(result)
    print(coef_values)
  }
  return(list(last_train_modelroc = roc_train, results, DMR_predictor, cv_fit_list)) 
}

# *****************************************
# ******** Train set are different *******
# *****************************************

# lasso_logistic_regression_for_select_noncolinear_predictor <- function(data, 
#                                                                        target_var, 
#                                                                        scale_or_not = TRUE,
#                                                                        num_iterations = 100, 
#                                                                        family = "binomial",
#                                                                        typemeasure = "class"){

#   results <- data.frame()
#   DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(data)))

#   for (k in 1:num_iterations) {
#     set.seed(k)
#     level1 <- subset(data, data[[target_var]] == 0)
#     level2 <- subset(data, data[[target_var]] == 1)
#     index_level1 <- sample(nrow(level1), 0.7*nrow(level1)) # this will affect the number of final non-colinear variables
#     index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
#     train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
#     test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))

#     train_x <- subset(train, select = setdiff(names(train), c(target_var)))
#     train_y <- as.numeric(train[[target_var]])
#     test_x <- subset(test, select = setdiff(names(test), c(target_var)) )
#     test_y <- as.numeric(test[[target_var]])

#     #scale data values
#     # scaled_train_x <- scale(train_x)
#     # scaled_test_x <- scale(test_x)
#     if (scale_or_not == TRUE) {
#       scaled_train_x <- apply(train_x, 2, scale)
#       scaled_test_x <- apply(test_x, 2, scale)
#     } else {
#       scaled_train_x <- train_x
#       scaled_test_x <- test_x
#     }

#     # assemble data frames
#     train <- cbind(data.frame(scaled_train_x), 
#                    subset(train, select = c(target_var)))
#     test <- cbind(data.frame(scaled_test_x), 
#                   subset(test, select = c(target_var)))

#     # Fit the logistic regression model with Lasso penalty using cross-validation
#     cv_fit <- cv.glmnet(as.matrix(scaled_train_x), train_y, family = "binomial", alpha = 1, type.measure="class") # a logistic model with lasso penalty 
#     # Get the best lambda
#     best_lambda <- cv_fit$lambda.min
#     # Fit the final model with the best lambda
#     fit <- glmnet(as.matrix(scaled_train_x), train_y, family = "binomial", alpha = 1, lambda = best_lambda)


#     # Store coefficients
#     coef_values <- coef(fit)[,]
#     coef_values <- coef_values[coef_values != 0]
#     coef_names <- names(coef_values)
#     coef_names <- c(coef_names, rep(NA, nrow(DMR_predictor) - length(coef_names)))
#     coef_values_attach <- c(coef_values, rep(NA, nrow(DMR_predictor) - length(coef_values)))
#     iteration_coef <- data.frame(coef_names, coef_values_attach)
#     DMR_predictor <- cbind(DMR_predictor, iteration_coef)

#     # Prediction and ROC analysis
#     train_pre <- predict(fit, newx = as.matrix(scaled_train_x), type = 'response', re.form = NA)
#     train_modelroc <- roc(train[[target_var]], as.numeric(train_pre)) #计算ROC绘图数据
#     roc_train = train_modelroc$auc
#     accuracy_train <- sum(diag(table(train_y, ifelse(train_pre < 0.5, 0, 1)))) / length(train_y)

#     test_pre <- predict(fit, newx = as.matrix(scaled_test_x), type = 'response', re.form = NA)
#     test_modelroc <- roc(test_y, as.numeric(test_pre)) # Calculate ROC for test data
#     roc_test <- test_modelroc$auc
#     accuracy_test <- sum(diag(table(test_y, ifelse(test_pre < 0.5, 0, 1)))) / length(test_y)

#     result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
#     results <- rbind(results, result)

#     print(k)
#     print(result)
#     print(coef_values)
#   }
#   return(list(last_test_modelroc = roc_test, results, DMR_predictor)) # roc_test will be from the last iteration
# }





logistic_regression_prediction_by_lme4 <- function(data, target_var, random_effects_var,scale_or_not = TRUE,
                                           num_iterations = 100, family = "binomial") {

  results <- data.frame()
  DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(data)))

  for (k in 1:num_iterations) {
    set.seed(k)
    level1 <- subset(data, data[[target_var]] == 0)
    level2 <- subset(data, data[[target_var]] == 1)
    index_level1 <- sample(nrow(level1), 0.7*nrow(level1)) # this will affect the number of final non-colinear variables
    index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
    test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))

    train_x <- subset(train, select = names(train)[-which(names(train) %in% c(target_var, random_effects_var))])
    train_y <- as.numeric(train[[target_var]])
    test_x <- subset(test, select = names(test)[-which(names(test) %in% c(target_var, random_effects_var))])
    test_y <- as.numeric(test[[target_var]])

    #scale data values
    # scaled_train_x <- scale(train_x)
    # scaled_test_x <- scale(test_x)
    if (scale_or_not == TRUE) {
      scaled_train_x <- apply(train_x, 2, scale)
      scaled_test_x <- apply(test_x, 2, scale)
    } else {
      scaled_train_x <- train_x
      scaled_test_x <- test_x
    }

    # assemble data frames
    train <- cbind(data.frame(scaled_train_x), 
                   subset(train, select = c(target_var, random_effects_var)))
    test <- cbind(data.frame(scaled_test_x), 
                  subset(test, select = c(target_var, random_effects_var)))

    predictors <- colnames(data)[-which(names(data) %in% c(target_var, random_effects_var))]
    formula_str <- paste(target_var, "~", paste(predictors, collapse = " + "), "+ (1 |", random_effects_var, ")")
    model_formula <- as.formula(formula_str)

    fit <- glmer(model_formula, data = train, family = family)

    # Store coefficients
    coef_values <- fixef(fit)
    coef_names <- names(coef_values)
    coef_names <- c(coef_names, rep(NA, nrow(DMR_predictor) - length(coef_names)))
    coef_values <- c(coef_values, rep(NA, nrow(DMR_predictor) - length(coef_values)))
    iteration_coef <- data.frame(coef_names, coef_values)
    DMR_predictor <- cbind(DMR_predictor, iteration_coef)

    # Prediction and ROC analysis
    train_pre <- predict(fit, newdata = train, type = 'response', re.form = NA)
    train_modelroc <- roc(train[[target_var]], as.numeric(train_pre)) #计算ROC绘图数据
    roc_train = train_modelroc$auc
    accuracy_train <- sum(diag(table(train_y, ifelse(train_pre < 0.5, 0, 1)))) / length(train_y)

    test_pre <- predict(fit, newdata = test, type = 'response', re.form = NA)
    test_modelroc <- roc(test_y, as.numeric(test_pre)) # Calculate ROC for test data
    roc_test <- test_modelroc$auc
    accuracy_test <- sum(diag(table(test_y, ifelse(test_pre < 0.5, 0, 1)))) / length(test_y)

    result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
    results <- rbind(results, result)

    print(k)
    print(result)

  }

  return(list(last_test_modelroc = roc_test, results, DMR_predictor)) # roc_test will be from the last iteration
}






library(data.table)
library(glmmTMB)
library(pROC)
library(randomForest)


logistic_regression_prediction_by_glmmTMB <- function(data, target_var, random_effects_var, scale_or_not = TRUE,
                                           num_iterations = 100, family = "binomial"){

  results <- data.frame()
  DMR_predictor <- data.frame(matrix(ncol = 1, nrow = ncol(data)))

  for (k in 1:num_iterations) {
    set.seed(k+500)
    level1 <- subset(data, data[[target_var]] == 0)
    level2 <- subset(data, data[[target_var]] == 1)
    index_level1 <- sample(nrow(level1), 0.7*nrow(level1))
    index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
    test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))

    train_x <- subset(train, select = names(train)[-which(names(train) %in% c(target_var, random_effects_var))])
    train_y <- as.numeric(train[[target_var]])
    test_x <- subset(test, select = names(test)[-which(names(test) %in% c(target_var, random_effects_var))])
    test_y <- as.numeric(test[[target_var]])

    #scale data values
    # scaled_train_x <- scale(train_x)
    # scaled_test_x <- scale(test_x)
    if (scale_or_not == TRUE) {
      scaled_train_x <- apply(train_x, 2, scale)
      scaled_test_x <- apply(test_x, 2, scale)
    } else {
      scaled_train_x <- train_x
      scaled_test_x <- test_x
    }
    # assemble data frames
    train <- cbind(data.frame(scaled_train_x), 
                   subset(train, select = c(target_var, random_effects_var)))
    test <- cbind(data.frame(scaled_test_x), 
                  subset(test, select = c(target_var, random_effects_var)))

    predictors <- colnames(data)[-which(names(data) %in% c(target_var, random_effects_var))]
    formula_str <- paste0(target_var, "~", paste(predictors, collapse = " + "), "+ (1 | ", random_effects_var, ")")
    model_formula <- as.formula(formula_str)

    fit <- glmmTMB(model_formula, data = train, family=family)

    # Store results
    coef_values <- fixef(fit)[["cond"]][!is.na(fixef(fit)[["cond"]])]
    coef_names <- names(coef_values)
    coef_names <- c(coef_names, rep(NA, nrow(DMR_predictor) - length(coef_names)))
    coef_values <- c(coef_values, rep(NA, nrow(DMR_predictor) - length(coef_values)))
    iteration_coef <- data.frame(coef_names, coef_values)
    DMR_predictor <- cbind(DMR_predictor, iteration_coef)

    # Prediction ----
    train_pre <- predict(fit, newdata = train, type='response', re.form=NA)
    train_modelroc <- roc(train[[target_var]], as.numeric(train_pre)) #计算ROC绘图数据
    roc_train = train_modelroc$auc
    tab_train <- table(train_y, ifelse(train_pre < 0.5, 0, 1), dnn = c("true", "pre"))
    accuracy_train <- sum(diag(tab_train)) / length(train_y)

    test_pre <- predict(fit, newdata = test, type='response', re.form=NA)
    test_modelroc <- roc(test_y, as.numeric(test_pre)) # Calculate ROC for test data
    roc_test <- test_modelroc$auc
    tab_test <- table(test_y, ifelse(test_pre < 0.5, 0, 1), dnn = c("true", "pre"))
    accuracy_test <- sum(diag(tab_test)) / length(test_y)

    # pdf("Logistic_ROC.pdf", width = 5, height = 5)
    # plot(test_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
    #      grid.col=c("green", "red"), max.auc.polygon=TRUE,
    #      auc.polygon.col="skyblue", print.thres=TRUE) #绘制ROC曲线图
    # plot(ci(test_modelroc, of = "thresholds", thresholds = "best")) #添加最佳截断值位置
    # dev.off()

    result <- data.frame(accuracy_train, accuracy_test, roc_train, roc_test)
    results <- rbind(results, result)

    print(result)
    print(k)
  }
  return(list(test_modelroc, results, DMR_predictor)) #test_modelroc will be from the last iteration
}

