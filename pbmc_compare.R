library(dplyr)
library(ggplot2)
library(reshape2)

d.ori <- read.table("pbmc_latent/pbmc_00000.csv", sep=",", header=TRUE)
d.after <- read.table("pbmc_latent/pbmc_11111_2.csv", sep=",", header=TRUE)

d1 <- read.table("pbmc_latent/pbmc_no_clustering.csv", sep=",", header=TRUE)
d2 <- read.table("pbmc_latent/pbmc_no_kmeans.csv", sep=",", header=TRUE)
d3 <- read.table("pbmc_latent/pbmc_no_cell_cell.csv", sep=",", header=TRUE)
d4 <- read.table("pbmc_latent/pbmc_no_consistency.csv", sep=",", header=TRUE)

#d1 <- read.table("pbmc_latent/pbmc_1011.csv", sep=",", header=TRUE)
d2 <- read.table("pbmc_latent/pbmc_11011_2.csv", sep=",", header=TRUE)
d3 <- read.table("pbmc_latent/pbmc_11101_2.csv", sep=",", header=TRUE)
d4 <- read.table("pbmc_latent/pbmc_11110_2.csv", sep=",", header=TRUE)

d.ori <- 100*d.ori/max(d.ori)
d.after <- 100*d.after/max(d.after)
#d1 <- 100*d1/max(d1)
d2 <- 100*d2/max(d2)
d3 <- 100*d3/max(d3)
d4 <- 100*d4/max(d4)

id <- read.table("pbmc_latent/pbmc_idents.csv", sep="\t", header=TRUE)
id.order <- order(id$idents)
id2 <- id[id.order,]

d <- d.ori
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Original"
d.all <- d.plot

d <- d1
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Clustering"
d.all <- rbind(d.all, d.plot)

d <- d2
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Center"
d.all <- rbind(d.all, d.plot)

d <- d3
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "KL"
d.all <- rbind(d.all, d.plot)

d <- d4
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Consistency"
d.all <- rbind(d.all, d.plot)

p <- ggplot(data = d.all, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 250)
p

d.between <- d.all[id2$idents[d.all$B] != id2$idents[d.all$variable],]
p1 <- ggplot(data = d.between, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 250)
p1

d.within <- d.all[id2$idents[d.all$B] == id2$idents[d.all$variable],]
p2 <- ggplot(data = d.within, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 250)
p2

d <- d.ori
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "Before"
d.all <- d.plot

d <- d.after
d <- d[id.order,]
d.dist <- dist(d)
d.dist <- as.matrix(d.dist)
d.dist <- as.data.frame(d.dist)
colnames(d.dist) <- c(1:ncol(d.dist))
rownames(d.dist) <- c(1:nrow(d.dist))
d.plot <- d.dist %>% mutate(B=colnames(d.dist)) %>% melt()
d.plot$B <- as.integer(d.plot$B)
d.plot$variable <- as.integer(d.plot$variable)
#d.plot$value <- 100*d.plot$value/max(d.plot$value)
p <- ggplot(data = d.plot, aes(x=value)) + geom_freqpoly(bins = 250)
p
d.plot$type <- "After"
d.all <- rbind(d.all, d.plot)

p <- ggplot(data = d.all, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 250)
p

d.between <- d.all[id2$idents[d.all$B] != id2$idents[d.all$variable],]
p1 <- ggplot(data = d.between, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 100)
p1

d.within <- d.all[id2$idents[d.all$B] == id2$idents[d.all$variable],]
p2 <- ggplot(data = d.within, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 100)
p2

var(d.between$value[d.between$type=="KL"])
var(d.within$value[d.within$type=="KL"])
var(d.between$value[d.between$type=="Original"])
var(d.within$value[d.within$type=="Original"])
var(d.between$value[d.between$type=="After"])
var(d.within$value[d.within$type=="After"])

mean(d.between$value[d.between$type=="KL"])
mean(d.within$value[d.within$type=="KL"])
mean(d.between$value[d.between$type=="Original"])
mean(d.within$value[d.within$type=="Original"])
mean(d.between$value[d.between$type=="After"])
mean(d.within$value[d.within$type=="After"])

p <- ggplot(data = d.all[d.all$type=="After",], aes(x=B,y=variable)) + geom_tile(aes(fill = value))
p <- p + scale_fill_continuous(low = "blue", high = "yellow")
p <- p + xlab("Cells") + ylab("Cells")
p

d.all$group <- "within"
d.all[id2$idents[d.all$B] != id2$idents[d.all$variable],5] <- "between"

p <- ggplot(data = d.all[d.all$type != "Consistency",], mapping = aes(x = value, y=..density.., color = type, linetype = group)) + geom_freqpoly(bins = 250)
p

d.test <- d.all[d.all$type=="After" | d.all$type=="KL",]
p <- ggplot(data = d.test, mapping = aes(x = value, y=..density.., color = group, linetype = type)) + geom_freqpoly(bins = 250)
p

d.test <- d.all[d.all$type=="After" | d.all$type=="Clustering",]
p <- ggplot(data = d.test, mapping = aes(x = value, y=..density.., color = group, linetype = type)) + geom_freqpoly(bins = 250)
p

d.test <- d.all[d.all$type=="After" | d.all$type=="Center",]
p <- ggplot(data = d.test, mapping = aes(x = value, y=..density.., color = group, linetype = type)) + geom_freqpoly(bins = 250)
p

d.test <- d.all[d.all$type=="After" | d.all$type=="Consistency",]
p <- ggplot(data = d.test, mapping = aes(x = value, y=..density.., color = group, linetype = type)) + geom_freqpoly(bins = 250)
p

d.test <- d.all[d.all$type=="Original" | d.all$type=="After",]
p <- ggplot(data = d.test, mapping = aes(x = value, y=..density.., color = group, linetype = type)) + geom_freqpoly(bins = 250)
p

d.select <- d.all[id2$idents[d.all$B]=="NK" & id2$idents[d.all$variable]=="Memory CD4 T",]
p <- ggplot(data = d.select, mapping = aes(x = value, y=..density.., color = type)) + geom_freqpoly(bins = 100)
p














