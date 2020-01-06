#the used packages
library(caret)
library(readr)

#upload data
processed_squamous_done <- read_table2("processed_squamous_done.txt")
View(processed_squamous_done)
x= processed_squamous_done

#Data splicing
intrain <- createDataPartition(y = x$label, p= 0.6, list = F)
training <- x[intrain,]
testing <- x[-intrain,]


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

sensitivity = 0.89474         
specificity = 0.99500 
F1_score = 2*((sensitivity*specificity)/(sensitivity+specificity))
