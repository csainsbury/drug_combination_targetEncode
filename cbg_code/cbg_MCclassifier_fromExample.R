## original example from : https://rpubs.com/mharris/multiclass_xgboost

library("xgboost")  # the main algorithm
library("archdata") # for the sample dataset
library("caret")    # for the confusionmatrix() function (also needs e1071 package)
library("dplyr")    # for some data preperation
#library("Ckmeans.1d.dp") # for xgb.ggplot.importance

export1 <- fread(paste0('Documents/data/CBGdata/abstract_exports/MULTICLASS_export_admissionDuration_4_days_hypothresh_NAs_included_3.csv'))

      x <- export1

      no_event_ids <- unique(x[label==0]$id)
      hypo_ids <- unique(x[label==1]$id)
      hyper_ids <- unique(x[label==2]$id)
      both_ids <- unique(x[label==3]$id)
      
      ratio = 2
      hr <- length(hypo_ids)
      
      all_ids <- unique(x$id)
      
      if (length(no_event_ids) > hr * ratio) {
        sample_no_event_ids <- no_event_ids[sample(length(no_event_ids), hr *ratio)]
      }
      if (length(hyper_ids) > hr * ratio) {
        sample_hyper_ids <- hyper_ids[sample(length(hyper_ids), hr *ratio)]
      } else {
        sample_hyper_ids <- hyper_ids
      }
      if (length(both_ids) > hr * ratio) {
        sample_both_ids <- both_ids[sample(length(both_ids), hr *ratio)]
      } else {
        sample_both_ids <- both_ids
      }
      
      id_sample <- c(hypo_ids, sample_no_event_ids, sample_hyper_ids, sample_both_ids)
      
      sample_x <- x[id %in% id_sample]

      subsample <- sample_x

# max = nrow(sample_x)
# subsample <- export1[sample(nrow(export1), max), ]
# set random seed
set.seed(42)

# data(RBGlass1)
dat <- subsample 
dat$Site <- as.numeric(dat$label)
dat <- dat %>% select(-id, -label, -max_by_day_4, -min_by_day_4, -median_by_day_4, -day_N_4, -iqr_by_day_4)
# dat_add <- dat[which(dat$Site == 0),] %>%
#   rowwise() %>%
#   mutate_each(funs(./10 + rnorm(1,.,.*0.1))) %>%
#   mutate_each(funs(round(.,2))) %>%
#   mutate(Site = 3)
# 
# dat <- rbind(dat, dat_add) %>%
#   mutate(Site = Site - 1)

summary(dat)

##
# Make split index
train_index <- sample(1:nrow(dat), nrow(dat)*0.75)
# Full data set
data_variables <- as.matrix(dat[,-18])
data_label <- dat[,"Site"]
data_matrix <- xgb.DMatrix(data = as.matrix(dat), label = data_label$Site)
# split train data and make xgb.DMatrix
train_data   <- data_variables[train_index,]
train_label  <- data_label[train_index]
train_matrix <- xgb.DMatrix(data = train_data, label = train_label$Site)
# split test data and make xgb.DMatrix
test_data  <- data_variables[-train_index,]
test_label <- data_label[-train_index]
test_matrix <- xgb.DMatrix(data = test_data, label = test_label$Site)

##
numberOfClasses <- length(unique(dat$Site))
xgb_params <- list(max.depth= 2,
                   "objective" = "multi:softprob",
                   "eval_metric" = "mlogloss",
                   "num_class" = numberOfClasses)
nround    <- 85 # number of XGBoost rounds
cv.nfold  <- 5

# Fit cv.nfold * cv.nround XGB models and save OOF predictions
cv_model <- xgb.cv(params = xgb_params,
                   data = train_matrix, 
                   nrounds = nround,
                   nfold = cv.nfold,
                   verbose = 1,
                   prediction = TRUE)

##
OOF_prediction <- data.frame(cv_model$pred) %>%
  mutate(max_prob = max.col(., ties.method = "last"),
         label = train_label + 1)
head(OOF_prediction)

##
# confusion matrix
confusionMatrix(factor(as.matrix(OOF_prediction$max_prob)),
                factor(as.matrix(OOF_prediction$label)),
                mode = "everything")

##
bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = nround)

# Predict hold-out test set
test_pred <- predict(bst_model, newdata = test_matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test_label + 1,
         max_prob = max.col(., "last"))
# confusion matrix of test set
confusionMatrix(factor(as.matrix(test_prediction$max_prob)),
              factor(as.matrix(test_prediction$label)),
                mode = "everything")

##
# get the feature real names
names <-  colnames(dat[,-1])
# compute feature importance matrix
importance_matrix = xgb.importance(feature_names = names, model = bst_model)
head(importance_matrix)

gp = xgb.plot.importance(importance_matrix)
print(gp) 
