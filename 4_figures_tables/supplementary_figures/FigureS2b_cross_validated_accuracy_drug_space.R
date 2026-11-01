#######################
#  R でもやってみる   #
# 　 Randomforest　　 #
#######################
library(randomForest)
library(caret)
library(jsonlite)
library(ggplot2)
library(reshape2)

### 1. load ################################
X_train <- read.csv("./path/to/intermediate/data/4.Feature_Selection/X_train_filtered_Morgan.csv", header = F)
Y_train <- read.csv("./path/to/intermediate/data/4.Feature_Selection/Y_train_filtered_Morgan.csv", header = F)

X_test <- read.csv("./path/to/intermediate/data/4.Feature_Selection/X_test_filtered_Morgan.csv", header = F)
Y_test <- read.csv("./path/to/intermediate/data/4.Feature_Selection/Y_test_filtered_Morgan.csv", header = F)

Feature_names <- read.csv("./path/to/intermediate/data/4.Feature_Selection/feature_names_selected_Morgan.csv", header = F)

rfe_results <- fromJSON("./path/to/intermediate/data/4.Feature_Selection/rfe_results_Morgan.json")
print(rfe_results)

### 2. data 再構築 ########################
data_train           <- X_train
data_train$label     <- as.factor(Y_train$V1)
colnames(data_train)[1:64] <- Feature_names$V1

### 3.  モデル構築 ########################
# クロスバリデーションの設定
ctrl <- trainControl(method = "cv", number = 5)

list_rf <-
  lapply(rfe_results, function(x){
    tmp_data       <- as.data.frame(data_train[,x])
    tmp_data$label <- data_train$label
    colnames(tmp_data) <- c(x, "label")
  # 学習 & 評価（randomForestを内部的に使用）
  model_cv <- train(label ~ ., data = tmp_data, method = "rf", trControl = ctrl)
})

# 結果表示
print(list_rf[[1]])

### 3.  結果表示 ########################
df <- melt(sapply(list_rf, function(x){(x$resample$Accuracy)}))

g <- ggplot(df, aes(x = as.factor(Var2), y = value)) + 
  geom_boxplot(fill="#00B0F0") + xlab("Number of Selected Features") + 
  ylab("Cross-Validated Accuracy") + theme_bw() + theme(text = element_text(size = 24))
plot(g)　#やはり、1か3かがバランスがいいようだ

ggsave("./path/to/intermediate/data/4.Feature_Selection/Cross_validated_Accuracy.pdf", width = 10.5, height = 7.8)
