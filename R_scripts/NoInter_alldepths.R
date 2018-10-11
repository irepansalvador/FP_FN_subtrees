### this will read all the results from the SEQuoia (No Intertarget) False positives estimation,
### of the Unweighted option, at different depths of the trees.
### The total number of trees is 1,000, and for each trees there are 10 subsamples of 1000 cells

# Packages required
require(lattice)
require(reshape2)

## Set working directory
setwd("/home/irepan/Desktop/Github/FP_FN_subtrees/R_scripts/")

## declare important variables
depth <- 1:9
samples <- 1:10
trees <- 1:1000

############# Open output

all_results = NULL
reps = NULL

for (i in trees)  ## for each output file (which contains 50 results)
{
  for (ii in samples)  ## for each output file (which contains 50 results)
  {
    my_tree = ((i-1)*10)+ii 
    
    my.result = read.table(file = paste("../No_InterTarget/ALL_results/FP_t",my_tree,sep = ""),
                           colClasses = c("numeric") )
    if (ii == 1) {reps  = my.result[,1] }
    if (ii > 1)  {reps <- rbind(reps,my.result[,1])}
  }

  if (i == 1) {all_results  = apply(reps, MARGIN = 2, FUN = mean) }
  if (i > 1)  {all_results <- rbind(all_results,apply(reps, MARGIN = 2, FUN = mean))}
}

head(all_results)
boxplot(all_results)

## format the dataframe
d = depth
d = paste("depth",d,sep = "" )
colnames(all_results) = d
rownames(all_results)= c(1:1000)
my_means = apply(all_results, MARGIN = 2, FUN = mean)


## plot

svg(filename = "NoInter_FP_alldepths_1Ksims_10samples.svg",width = 4,height = 4)

boxplot(all_results[,9:1], ylim= c(0,0.5), 
        main = "No InterTarget", ylab= "False Positives",
        las = 3, boxwex = 0.5, outline = F, notch = T)
lines(x = my_means[9:1], col = "red", lwd=2)
points(x = my_means[9:1], col = "red", lwd=2)
grid()

dev.off()
