library(tidyverse)
library(readr)

E1_A_filtered <- read_csv("E1_A_edited.txt", col_names = FALSE)
E1_A <- read_delim("E1_A.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

E1_AllA <- sum(E1_A$X5)
E1_Gs <- sum(E1_A_edited$X3)
E1_percentEdited <- E1_Gs/E1_AllA
E1_AllA
E1_Gs
E1_percentEdited
