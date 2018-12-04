#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ape)
library(stringdist)
library(data.tree)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

a = 	args[1]
b = 	args[2]
c = 	args[3]
d	=   args[4]
r = paste(replicate(nchar(a), "0"),collapse = "")

my_cells  <- c(a,b,c,d)
my_labels <- c("1","2","3","4")
my_rand <- sample(x=c(1:4), size = 4,replace = F)

## create dist matrix
#mydist <- stringdistmatrix(c(args[1:4],r), useNames = F)
mydist <- stringdistmatrix(c(my_cells[my_rand],r), useNames = F)
mydist=as.matrix(mydist, labels=TRUE)
dimnames(mydist) <- list( c(my_labels[my_rand],"5"), c(my_labels[my_rand],"5"))

## calculate NJ
mytree <-  nj(mydist)
## reroot
r5 <- root(mytree, resolve.root = T, outgroup =  5 )
# write tree to output
write.tree(r5, file = paste(args[5], ".nw", sep =""), digits = 2)
