#the used packages
library(caret)
library(readr)

#upload data
processed_squamous <- read_table2("processed_squamous.txt")
View(processed_squamous)

#Data splicing
intrain <- createDataPartition(y = processed_squamous$label, p= 0.7, list = F)
training <- processed_squamous[intrain,]
testing <- processed_squamous[-intrain,]


#Training data
svm_Linear <- train(label ~., data = training, method = "svmLinear",
                    trControl= trainControl(method = "repeatedcv", number = 10, repeats = 3))
                    #preProcess = c("center", "scale"),
                    #tuneLength = 10)
svm_Linear


#Testing data
test_pred <- predict(svm_Linear, newdata = testing)
test_pred
confusionMatrix(table(test_pred, testing$label))
