library(readr)
library(devtools)
library(ggbiplot)
processed_squamous <- read_table2("processed_squamous.txt")
View(processed_squamous)
x = processed_squamous
table(x$label)
x.label <- as.factor(x$label)
#determine the class of both to ensure that they are factor
class(x.label)
x.abundance = x[ ,1:916]
x.pca <- prcomp(x.abundance)

#PCA for IRIS(g)
g <- ggbiplot(x.pca, obs.scale = 1, var.scale = 1, var.axes=FALSE,groups = x.label, ellipse = F,circle = T)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

print(g)

