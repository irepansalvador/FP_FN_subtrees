### this will read all the results from the SEQuoia (No Intertarget) terminal 4cells subtrees analysis,
### of the Unweighted option, at different depths of the trees.
### The total number of trees is 1,000, and for each trees there are 10 subsamples of 1000 cells

## This analysis takes cells directly from the simulations to reconstruct 4cell-subtrees. I take
## 4 cells (taxa 1-4) that I know come from a common ancestor and reconstruct the tree only with those
## and a "root" taxa that is the unmutated construct (taxon 5)
## I compare the inferred 4-cell tree with the reference (or real) tree, which is always ((1,2),(3,4),5)

## In cases when a tree has less than 1,000 cells (divisions 2-9) I extract all the 4cell-trees.
## When there are more than 1,000 cells (divisions 10-16) I extract 250 random 4cell trees (total cells = 1000)

## The result of each 4cell tree is in one file for each simulation. This script opens each file and store the
## mean accuracy of the 4cell trees in a matrix (15x1000, for cell divisions and cells). The final output 
## should be a boxplot with the accuracy of the 4cell-trees at different time points.

# Packages required
require(lattice)
require(reshape2)

## Set working directory
setwd("/home/irepan/Desktop/Github/FP_FN_subtrees/R_scripts/")

## declare important variables
divs <- 2:16
trees <- 1:1000

############# Open output

all_results = matrix(data = NA, nrow = length(trees), ncol = length(divs))
colnames(all_results) <- paste(divs,"div",sep = "")
rownames(all_results) <- paste("tree",trees,sep = "")

for (d in divs)
{
  for (i in trees)  ## for each output file (which contains 50 results)
    {
    t = sprintf("%.4d", i)
    my.result = read.table(file = paste("../GESTALT/4cell_trees/GESTALT_",
                                        sprintf("%.2d", d), "div_10_targets_60_states_",
                                        t,"rep-4cells_Acc",sep = ""),
                           colClasses = c("numeric") )
    if (d == 2)   {all_results[i,d-1]  = my.result[,1] }
    if (d > 2)    {all_results[i,d-1]  = mean(my.result[,1])}
    }
}

my_means = apply(all_results, 2, function(x) mean(x,na.rm= T))


svg(filename = "GESTALT_4cell_analysis_1Ksims.svg",width = 5,height = 4)

boxplot(all_results, outline = F, boxwex = 0.6, ylim = c(0,1),
        las = 3, main = "4cell-trees accuracy\n(GESTALT)",
        xlab = "Time point", ylab = "Accuracy (RF)")
points(my_means, pch = 5, col = "red", lwd= 2)
lines(my_means, pch = 5, col = "red", lwd = 2)
grid()

dev.off()
