#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ape)
library(stringdist)
library(data.tree)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
r = paste(replicate(nchar(args[1]), "0"),collapse = "")

## create dist matrix
mydist <- stringdistmatrix(c(args[1:4],r), useNames = F)
## calculate NJ
mytree <-  nj(mydist)
## reroot
r5 <- root(mytree, resolve.root = T, outgroup =  5 )
# write tree to output
write.tree(r5, file = paste(args[5], ".nw", sep =""), digits = 2)
