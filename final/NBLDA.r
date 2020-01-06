#needed packages  
library(S4Vectors)
library(DESeq2)
library(MLSeq)

#uploading the data and manipulate the file
exp <- read.table("C:/Users/QUEEN/Desktop/ml/ml/exp_res_degs_squamous[done].txt", quote="\"", comment.char="")
#exp=exp[-917]
#exp=t(exp)
#exp=as.data.frame(exp)
#row.names(exp) = exp[,1]
#exp=exp[-1]

#create a class label information in order to apply classification models
class <- DataFrame(condition = factor(rep(c("N","T"), c(49, 502))))

#convert character data frame to numeric 
genes=rownames(exp)
exp=apply(exp,2,as.integer)
row.names(exp)=genes

#splitting the data 60 training and 40 test because we have a small number of samples
nTest <- ceiling(ncol(exp) * 0.4)
ind <- sample(ncol(exp), nTest, FALSE)
data.train <- as.matrix(exp[ ,-ind] + 1)
data.test <- as.matrix(exp[ ,ind] + 1)
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])

#The training and test sets are stored in a DESeqDataSet using related functions from DESeq2. This object is then used as input for MLSeq.
data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr, design = formula(~condition))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts, design = formula(~condition))


#choosing the nblda algorithm to classify the data and its tuning parameters
#Negative-Binomial linear discriminant analysis (NBLDA) one of the most suitble 
#models for rna seq data because it is a discrete data and follow the negative binomial distrbution 
#choosing the deseq as a normlization method (estimateSizeFactors) because it is suitble for NBLDA
set.seed(2128)
ctrl.discrete <- discreteControl(method = "repeatedcv", number = 5, repeats = 10, tuneLength = 10)
fit.nblda <- classify(data = data.trainS4, method = "NBLDA", normalize = "deseq", ref = "T", control = ctrl.discrete)

#fit.svm <- classify(data = data.trainS4, method = "svmRadial", preProcessing = "deseq-vst", ref = "T", tuneLength = 10, control = trainControl(method = "repeatedcv", number = 5, repeats = 10, classProbs = TRUE))
#show(fit.svm)
#trained(fit.svm)
#plot(fit.svm)
#pred.svm<- predict(fit.svm, data.testS4)
#show(pred.svm)

#show the model results ana the tuning results are obtained using setter function trained
show(fit.nblda)
trained(fit.nblda)

#Predicting the class labels of test samples
pred.nblda<- predict(fit.nblda, data.testS4)
show(pred.nblda)
plot(pred.nblda)

#the model performance for the prediction is summarized as below using confusionMatrix from caret.
pred.nblda <- relevel(pred.nblda, ref = "T")
actual <- relevel(classts$condition, ref = "T")
tbl <- table(Predicted = pred.nblda, Actual = actual)
confusionMatrix(tbl, positive = "T")

#calculate f1 score
Sensitivity <- 0.9750         
Specificity <- 1.0000 
f1= 2*((Specificity*Sensitivity)/(Specificity+Sensitivity))
