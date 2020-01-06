#required packages
library(readr)
library(devtools)
library(ggbiplot)

#upload data
processed_all_data_551 <- read_table2("processed_all data_551.txt")
View(processed_all_data_551)
x = processed_all_data_551

table(x$label)
x.label <- as.factor(x$label)

#determine the class to ensure that they are factor
class(x.label)
x.abundance = x[ ,1:25554]
x.pca <- prcomp(x.abundance)

#PCA for label
g <- ggbiplot(x.pca, obs.scale = 1, var.scale = 1, var.axes=FALSE,groups = x.label, ellipse = F,circle = T)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

print(g)

