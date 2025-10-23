# Function for elastic net regression
# This is a script for using ENR to reduce the seature size or you can say feature selection
# The data should contains a "target variable", which is the classification, "0" for control, "1" for case
# data should be a data frame, the rownames are the sample names, and the colnames are the feature names

library(data.table)
library(glmnet)

elastic_net_regression <- function(data, target_var, num_iterations = 50, target_num_vars = 100, potential_lambda) {
  sites_selection <- data.frame(matrix(ncol = 1, nrow = 5000))
  MSE <- c()
  
  for (j in 1:num_iterations) {
    print(paste("**** Iteration:", j))
    set.seed(j)
    level1 <- subset(data, data[[target_var]] == 0)
    level2 <- subset(data, data[[target_var]] == 1)
    index_level1 <- sample(nrow(level1), 0.7*nrow(level1))
    index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ]))
    test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))
    
    train_x <- subset(train, select = names(train)[-which(names(train) == target_var)])
    train_y <- train[[target_var]]
    test_x <- subset(test, select = names(test)[-which(names(test) == target_var)])
    test_y <- test[[target_var]]

    
    print("# -- Fit the model")
    fit_ENR <- glmnet(x = as.matrix(train_x), y = as.matrix(train_y),
                      lambda = potential_lambda, family = "binomial",
                      alpha = 0.5, standardize = FALSE, trace.it = 1)
    
    print("# -- save a variable_selection_flow figure")
    if (j ==1) {
      pdf("variable_selection_flow.pdf", width = 5, height = 4)
      plot(fit_ENR,xvar="lambda", lwd=2, label=F)
      dev.off()
      }
    #pdf(paste("variable_selection_flow", j, ".pdf", sep = ""), width = 7, height = 7)
    #plot(fit_ENR,xvar="lambda", lwd=2, label=F)
    #dev.off()

    print("# -- Find the lambda value based on the target number of variables")
    coef_matrix <- coef(fit_ENR)
    num_selected <- colSums(coef_matrix != 0)
    lambda_target <- which(abs(num_selected - target_num_vars) == min(abs(num_selected - target_num_vars))) # find the lambda leading to 'target_num_vars' of variable 
    lambda_values <- fit_ENR$lambda[lambda_target]
    
    print("# -- According to the lambda find the target sites and coefficients")
    coef_ENR <- predict(fit_ENR, s = lambda_values, type = "coefficients")
    site_name <- names(train)[1 + coef_ENR@i]
    print(paste("There are", length(site_name), "sites selected"))
    
    # Store results
    site_name <- c(site_name, rep(NA, nrow(sites_selection) - length(site_name)))
    coef <- c(coef_ENR@x, rep(NA, nrow(sites_selection) - length(coef_ENR@x)))
    site_coef <- data.frame(cbind(site_name, coef))
    sites_selection <- cbind(sites_selection, site_coef)
    
    print("# -- Check the accuracy of elastic net regression")
    predicted <- predict(fit_ENR, s = lambda_values, newx = as.matrix(test_x), type = 'response')
    predicted <- ifelse(predicted > 0.5, 1, 0) #should match with the factor number on target variable
    prediction <- data.frame(cbind(test_y, predicted))
    #print(prediction)
    names(prediction) <- c("sample_group", "predicted")
    MSE[j] <- sum(prediction$sample_group == prediction$predicted) / nrow(prediction)
    print(MSE[j])
    sites_selection[, 1] <- c(MSE, rep(NA, nrow(sites_selection) - length(MSE)))
  }
  
  return(sites_selection)
}

