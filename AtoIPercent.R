library(readr)
library(tidyverse)
library(ggplot2)
library(IdeoViz)

sample_1_AG <- read_delim("sample_1_AG.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
sample_1_AG <- dplyr::filter(sample_1_AG, sample_1_AG$X9 < .90 )
mean_sample_1 <- mean(sample_1_AG$X9)

myList <- c(mean_sample_1, mean_sample_2)
myListnames <- c("sample_1","sample_2")

means <- as.data.frame(myList, row.names = myListnames)
boxplot(means$myList)
