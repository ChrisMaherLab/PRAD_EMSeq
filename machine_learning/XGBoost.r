# XGBoost to select the important features
# data should be a data frame, the rownames are the sample names, and the colnames are the feature names, 
# The data should contains a "target variable", which is the classification, "0" for control, "1" for case

# 'data' should be a data frame, the rownames are the sample names, and the colnames are the feature names
# 'target_var' is the name of classification variable in the 'data'
# 'num_iterations' is the number of iteration that how many times you want the model to repeat
# 'max_depth_values', 'eta_values', 'nround_values', 'gamma_values', 'min_child_weight_values' are the model parameter. More information in xgboost::xgboost()

library(data.table)
library(Matrix)
library(xgboost)



xgboost_for_feature_selection <- function(data, target_var, 
                                          num_iterations = 1,
                                          max_depth_values = 2:5,
                                          eta_values = 10^(-3:2),
                                          nround_values = seq(10, 300, by = 10),
                                          gamma_values = 10^(-3:2),
                                          min_child_weight_values = 10^(-3:1)) {
  
  # Initialize vectors to store results
  results <- data.frame()
  Xgboost_import <- data.frame(matrix(ncol = 1, nrow = ncol(data)))
  feature_name <- names(data)[-which(names(data) == target_var)]
  progress = 0
  N = num_iterations * length(max_depth_values) * length(eta_values)* length(nround_values)* length(gamma_values)* length(min_child_weight_values)
                                       
  for (seed in 1:num_iterations){
    set.seed(progress)
    # Split data into train and test sets
    level1 <- subset(data, data[[target_var]] == 0)
    level2 <- subset(data, data[[target_var]] == 1)
    index_level1 <- sample(nrow(level1), 0.7*nrow(level1))
    index_level2 <- sample(nrow(level2), 0.7*nrow(level2))
    train <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ])) # some sites might not be fitted in the model all the time
    test <- data.frame(rbind(level1[-index_level1, ], level2[-index_level2, ]))
    #train <- data # make sure all site will be through into algotism for selection.
    #test <- data.frame(rbind(level1[index_level1, ], level2[index_level2, ])) # Random 70% subset will be predicted, but this does not really matter, since the train set is the whole dataset. Just to have this data
    

    # Prepare data matrices
    train_x <- subset(train, select = names(train)[-which(names(train) == target_var)])
    train_y <- as.numeric(train[[target_var]]) 
    test_x <- subset(test, select = names(test)[-which(names(test) == target_var)])
    test_y <- as.numeric(test[[target_var]]) 
    
    # Scale features
    scaled_train_x <- scale(train_x)
    scaled_test_x <- scale(test_x)
    
    # merge features and results
    scaled_train <- cbind(sample_group = train_y, scaled_train_x)
    scaled_train <- data.frame(scaled_train)
    scaled_train$sample_group <- as.factor(scaled_train$sample_group)
    
    scaled_test <- cbind(sample_group = test_y, scaled_test_x)
    scaled_test <- data.frame(scaled_test)
    scaled_test$sample_group <- as.factor(scaled_test$sample_group)
    
    # Create sparse matrices
    train_matrix <- sparse.model.matrix(sample_group~ . - 1, data = scaled_train)
    test_matrix <- sparse.model.matrix(sample_group~ . - 1, data = scaled_test)
    
    # Define parameter grid
    max_depth_values <- max_depth_values
    eta_values <- eta_values
    nround_values <- nround_values
    gamma_values <- gamma_values
    min_child_weight_values <- min_child_weight_values
    
    for (max_depth_value in max_depth_values) {
      for (eta_value in eta_values) {
        for (nround_value in nround_values) {
          for (gamma_value in gamma_values) {
            for (min_child_weight_value in min_child_weight_values) {
              # Fit XGBoost model
              fit_xgb <- xgboost(data = train_matrix, label = train_y, nrounds = nround_value, 
                                 eta = eta_value, max_depth = max_depth_value, gamma = gamma_value, 
                                 min_child_weight = min_child_weight_value, verbose = 0, 
                                 objective = "binary:logistic", eval_metric = "error")
              
              # -- interpretation 
              importance <- xgb.importance(feature_names = feature_name, model = fit_xgb)
              importance <- c(importance$Feature, rep(NA, nrow(Xgboost_import) - length(importance$Feature)))
              Xgboost_import <- cbind(Xgboost_import, importance)
              
              # Prediction
              fit_xgb_predict_train <- predict(fit_xgb, newdata = train_matrix)
              fit_xgb_predict_test <- predict(fit_xgb, newdata = test_matrix)
              
              # Metrics
              MSE_train <- mean((fit_xgb_predict_train - train_y)^2)
              MSE_test <- mean((fit_xgb_predict_test - test_y)^2)
              accuracy_train <- sum(round(fit_xgb_predict_train) == train_y) / length(train_y)
              accuracy_test <- sum(round(fit_xgb_predict_test) == test_y) / length(test_y)
              
              # Store results
              result <- data.frame(max_depth = max_depth_value,
                                   eta = eta_value,
                                   nround = nround_value,
                                   gamma = gamma_value,
                                   min_child_weight = min_child_weight_value,
                                   MSE_train = MSE_train,
                                   MSE_test = MSE_test,
                                   accuracy_train = accuracy_train,
                                   accuracy_test = accuracy_test)
              
              results <- rbind(results, result)
              progress = progress + 1
              print(paste(progress, "/", N))
              print(result)
            }
          }
        }
      }
    }
  }
  return(list(Xgboost_import, results))
}


