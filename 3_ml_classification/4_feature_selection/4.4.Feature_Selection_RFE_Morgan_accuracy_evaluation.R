#######################
#  R でもやってみる   #
# 　　　精度検証　　  #
#######################
library(randomForest)
library(caret)
library(jsonlite)
library(ggplot2)
library(reshape2)

### 1. load ############################
X_train <- read.csv("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.csv", header = F)
Y_train <- read.csv("./path/to/intermediate/data/4.Feature_Selection/Y_train_filtered_Morgan.csv", header = F)
X_test <- read.csv("./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.csv", header = F)
Y_test <- read.csv("./path/to/intermediate/data/4.Feature_Selection/Y_test_filtered_Morgan.csv", header = F)
Feature_names <- read.csv("./path/to/intermediate/data/4.Feature_Selection/feature_names_selected_Morgan.csv", header = F)

### 2. data 再構築 ########################
data_train           <- X_train
data_train$label     <- as.factor(Y_train$V1)
colnames(data_train)[1:64] <- gsub(" ", "_", Feature_names$V1)

data_test           <- X_test
data_test$label     <- as.factor(Y_test$V1)
colnames(data_test)[1:64] <- gsub(" ", "_", Feature_names$V1)

### 3. データ準備 #########################
# RFEコントロール設定
rf_model <- randomForest(label ~ ., data = data_train[,c("Feature_1398", "label")], ntree = 100, importance = TRUE)

# 予測値を取得（クラスラベル）
predicted <- predict(rf_model, newdata = data_test[,c("Feature_1398", "label")])

# 実際のラベル（正解）
true_labels <- data_test$label

# 正解率
mean(predicted == true_labels)

# 混合行列
confusionMatrix(predicted, true_labels)

　
### 3. データ準備 #########################
# RFEコントロール設定
rf_model <- randomForest(label ~ ., data = data_train[,c("Feature_650", "label")], ntree = 100, importance = TRUE)

# 予測値を取得（クラスラベル）
predicted <- predict(rf_model, newdata = data_test[,c("Feature_650", "label")])

# 実際のラベル（正解）
true_labels <- data_test$label

# 正解率
mean(predicted == true_labels)

# 混合行列
confusionMatrix(predicted, true_labels)

### 3. データ準備 #########################
# RFEコントロール設定
rf_model <- randomForest(label ~ ., data = data_train[,c("Feature_650", "Feature_1088", "label")], ntree = 100, importance = TRUE)

# 予測値を取得（クラスラベル）
predicted <- predict(rf_model, newdata = data_test[,c("Feature_650", "Feature_1088", "label")])

# 実際のラベル（正解）
true_labels <- data_test$label

# 正解率
mean(predicted == true_labels)

# 混合行列
confusionMatrix(predicted, true_labels)

### 3. データ準備 #########################
# RFEコントロール設定
rf_model <- randomForest(label ~ ., data = data_train[,c("Feature_650", "Feature_1088", "Feature_1398", "label")], ntree = 100, importance = TRUE)

# 予測値を取得（クラスラベル）
predicted <- predict(rf_model, newdata = data_test[,c("Feature_650", "Feature_1088", "Feature_1398", "label")])

# 実際のラベル（正解）
true_labels <- data_test$label

# 正解率
mean(predicted == true_labels)

# 混合行列
confusionMatrix(predicted, true_labels)
