# Pacman is a package management tool 
install.packages("pacman")

library(pacman)

# p_load automatically installs packages if needed
p_load(xgboost, ParBayesianOptimization, mlbench, dplyr, skimr, recipes, resample)

# Load up some data
data("BostonHousing2")

# Data summary
skim(BostonHousing2)

# Predicting median house prices
rec <- recipe(cmedv ~ ., data = BostonHousing2) %>%
  
  # Collapse categories where population is < 3%
  step_other(town, chas, threshold = .03, other = "Other") %>% 
  
  # Create dummy variables for all factor variables 
  step_dummy(all_nominal_predictors())

# Train the recipe on the data set
prep <- prep(rec, training = BostonHousing2)

# Create the final model matrix
model_df <- bake(prep, new_data = BostonHousing2)

# All levels have been one hot encoded and separate columns have been appended to the model matrix
colnames(model_df)
##  [1] "tract"                  "lon"                    "lat"                   
##  [4] "medv"                   "crim"                   "zn"                    
##  [7] "indus"                  "nox"                    "rm"                    
## [10] "age"                    "dis"                    "rad"                   
## [13] "tax"                    "ptratio"                "b"                     
## [16] "lstat"                  "cmedv"                  "town_Boston.Savin.Hill"
## [19] "town_Cambridge"         "town_Lynn"              "town_Newton"           
## [22] "town_Other"             "chas_X1"

splits <- rsample::initial_split(model_df, prop = 0.7)

# Training set
train_df <- rsample::training(splits)

# Test set
test_df <- rsample::testing(splits)

dim(train_df)
## [1] 354  23

dim(test_df)
## [1] 152  23

# The xgboost interface accepts matrices 
X <- train_df %>%
  # Remove the target variable
  select(!medv, !cmedv) %>%
  as.matrix()

# Get the target variable
y <- train_df %>%
  pull(cmedv)

# Cross validation folds
folds <- list(fold1 = as.integer(seq(1, nrow(X), by = 5)),
              fold2 = as.integer(seq(2, nrow(X), by = 5)))

# Function must take the hyper-parameters as inputs
obj_func <- function(eta, max_depth, min_child_weight, subsample, lambda, alpha) {
  
  param <- list(
    
    # Hyter parameters 
    eta = eta,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    lambda = lambda,
    alpha = alpha,
    
    # Tree model 
    booster = "gbtree",
    
    # Regression problem 
    objective = "reg:squarederror",
    
    # Use the Mean Absolute Percentage Error
    eval_metric = "mape")
  
  xgbcv <- xgb.cv(params = param,
                  data = X,
                  label = y,
                  nround = 50,
                  folds = folds,
                  prediction = TRUE,
                  early_stopping_rounds = 5,
                  verbose = 0,
                  maximize = F)
  
  lst <- list(
    
    # First argument must be named as "Score"
    # Function finds maxima so inverting the output
    Score = -min(xgbcv$evaluation_log$test_mape_mean),
    
    # Get number of trees for the best performing model
    nrounds = xgbcv$best_iteration
  )
  
  return(lst)
}

bounds <- list(eta = c(0.001, 0.2),
               max_depth = c(1L, 10L),
               min_child_weight = c(1, 50),
               subsample = c(0.1, 1),
               lambda = c(1, 10),
               alpha = c(1, 10))
set.seed(1234)
bayes_out <- bayesOpt(FUN = obj_func, bounds = bounds, initPoints = length(bounds) + 2, iters.n = 10)

# Show relevant columns from the summary object 
bayes_out$scoreSummary[1:5, c(3:8, 13)]
##           eta max_depth min_child_weight subsample   lambda    alpha      Score
## 1: 0.13392137         8         4.913332 0.2105925 4.721124 3.887629 -0.0902970
## 2: 0.19400811         2        25.454160 0.9594105 9.329695 3.173695 -0.1402720
## 3: 0.16079775         2        14.035652 0.5118349 1.229953 5.093530 -0.1475580
## 4: 0.08957707         4        12.534842 0.3844404 4.358837 1.788342 -0.1410245
## 5: 0.02876388         4        36.586761 0.8107181 6.137100 6.039125 -0.3061535

# Get best parameters
data.frame(getBestPars(bayes_out))
##         eta max_depth min_child_weight subsample lambda alpha
## 1 0.1905414         8         1.541476 0.8729207      1     1

# Combine best params with base params
opt_params <- append(list(booster = "gbtree", 
                          objective = "reg:squarederror", 
                          eval_metric = "mae"), 
                     getBestPars(bayes_out))

# Run cross validation 
xgbcv <- xgb.cv(params = opt_params,
                data = X,
                label = y,
                nround = 100,
                folds = folds,
                prediction = TRUE,
                early_stopping_rounds = 5,
                verbose = 0,
                maximize = F)

# Get optimal number of rounds
nrounds = xgbcv$best_iteration

# Fit a xgb model
mdl <- xgboost(data = X, label = y, 
               params = opt_params, 
               maximize = F, 
               early_stopping_rounds = 5, 
               nrounds = nrounds, 
               verbose = 0)

# Evaluate performance 
actuals <- test_df$cmedv
predicted <- test_df %>%
  select_at(mdl$feature_names) %>%
  as.matrix %>%
  predict(mdl, newdata = .)

# Compute MAPE
mean(abs(actuals - predicted)/actuals)
## [1] 0.008391282