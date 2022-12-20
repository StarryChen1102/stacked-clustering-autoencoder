library(dplyr)
library(ggplot2)
library(reshape2)

scvi <- read.table("pbmc_scvi.csv", sep=",", header=TRUE)
mymodel <- read.table("pbmc_model.csv", sep=",", header=TRUE)
id <- read.table("pbmc_idents.csv", sep="\t", header=TRUE)
id.order <- order(id$x)

scvi <- scvi[id.order,]
mymodel <- mymodel[id.order,]

scvi.dist <- dist(scvi)
scvi.dist <- as.matrix(scvi.dist)
scvi.dist <- as.data.frame(scvi.dist)
colnames(scvi.dist) <- c(1:ncol(scvi.dist))
rownames(scvi.dist) <- c(1:nrow(scvi.dist))
scvi.plot <- scvi.dist %>% mutate(B=colnames(scvi.dist)) %>% melt()
scvi.plot$B <- as.integer(scvi.plot$B)
scvi.plot$variable <- as.integer(scvi.plot$variable)
scvi.plot$value <- 100*scvi.plot$value/max(scvi.plot$value)

p <- ggplot(data = scvi.plot, aes(x=B,y=variable)) + geom_tile(aes(fill = value))
p <- p + scale_fill_continuous(low = "blue", high = "yellow")
p <- p + xlab("Cells") + ylab("Cells")
p

mymodel.dist <- dist(mymodel)
mymodel.dist <- as.matrix(mymodel.dist)
mymodel.dist <- as.data.frame(mymodel.dist)
colnames(mymodel.dist) <- c(1:ncol(mymodel.dist))
rownames(mymodel.dist) <- c(1:nrow(mymodel.dist))
mymodel.plot <- mymodel.dist %>% mutate(B=colnames(mymodel.dist)) %>% melt()
mymodel.plot$B <- as.integer(mymodel.plot$B)
mymodel.plot$variable <- as.integer(mymodel.plot$variable)
mymodel.plot$value <- 100*mymodel.plot$value/max(mymodel.plot$value)

p <- ggplot(data = mymodel.plot, aes(x=B,y=variable)) + geom_tile(aes(fill = value))
p <- p + scale_fill_continuous(low = "blue", high = "yellow")
p <- p + xlab("Cells") + ylab("Cells")
p
