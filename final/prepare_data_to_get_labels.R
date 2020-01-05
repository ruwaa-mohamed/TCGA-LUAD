#upload packages

install.packages('tidyverse')
library(caret)
library(readr)
library(permute)
library(tidyverse)

#upload data
exp_res_degs_squamous_done_ <- read_table2("exp_res_degs_squamous[done].txt")
View(exp_res_degs_squamous_done_)
x = exp_res_degs_squamous_done_ 

#table(grepl( "01A" , names( expression_alldata_551 ) ))
colnames(x)

#call 01A as tumor, 11A as normal 
colnames(x)[2:50] <- "Normal"
colnames(x)[51:dim(x)[2]] <- "tumor"
colnames(x)

# remove genes column and set it as rows and adjust labels of normal and tumor
x2 <- x[,-1]
rownames(x2) <- x$Genes
colnames(x2)
x2=t(x2)
names <- rownames(x2)
rownames(x2) <- NULL
data <- cbind(x2, names)

colnames(data)[colnames(data) == "names"] <- 'label'
data= as.data.frame(data)
summary(data$label)
write.table(data, file = "processed_squamous_done.txt", row.names = F, quote = F)

