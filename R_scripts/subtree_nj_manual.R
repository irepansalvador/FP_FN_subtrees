library(ape)
library(stringdist)
library(data.tree)

setwd(dir = "/home/irepan/Desktop/Github/FP_FN_subtrees/R_scripts/")

# a = 	"Oxx0O0v000xut0x00N_Oxxxx_Qx0P0v0"
# b = 	"Oxx0O00000x0t00x0N00xxxx_Qx0P000"
# c = 	"O0x0O0x000x0t0w00N00xxxx_Qx0Px00"
# d	=   "O0x0O0000dx0tw000N00xxxx_Qx0P000"
# r = paste(replicate(32, "0"),collapse = "")

a = 	"1111"
b = 	"1111"
c = 	"1111"
d	=   "1111"
my_cells  <- c(a,b,c,d)
my_labels <- c("1","2","3","4")

r = paste(replicate(4, "0"),collapse = "")
my_rand <- sample(x=c(1:4), size = 4,replace = F)

## create dist matrix
#mydist <- stringdistmatrix(c(a,b,c,d,r), useNames = F)
mydist <- stringdistmatrix(c(my_cells[my_rand],r), useNames = F)
mydist=as.matrix(mydist, labels=TRUE)
dimnames(mydist) <- list( c(my_labels[my_rand],"5"), c(my_labels[my_rand],"5"))

mydist
## calculate NJ
mytree <-  nj(mydist)
## reroot
r5 <- root(mytree, resolve.root = T, outgroup =  5 )
plot(r5)
# write tree to output
write.tree(r5, file = "test.nw")

