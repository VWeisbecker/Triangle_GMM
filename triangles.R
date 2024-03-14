library(geomorph)
library(ape)

################################################################################
#data <- as.matrix(read.table("triangles3.txt"))
#data_mat <- arrayspecs(data,p=3,k=2)
#rand <- runif(32,0.5,2)
#data_mat[3,1,] <- rand
#plotAllSpecimens(data_mat)
#save(data_mat, file="triangle_coords.rda")
################################################################################

load("triangle_coords.rda")
gpa <- gpagen(data_mat)
plotOutliers(gpa$coords)
allom <- procD.lm(gpa$coords~gpa$Csize, iter=999); summary(allom)

plotAllSpecimens(data_mat)
plotAllSpecimens(gpa$coords)

tree <- read.tree("triangles_phy")
tree$tip.label
match(tree$tip.label, dimnames(gpa$coords)[[3]])
tree$edge.length[11] <- 1 # error making the tree on my part :/

### First test - random placement ###
phy1 <- procD.pgls(gpa$coords~gpa$Csize, phy=tree, iter=999); summary(phy1)
physignal(gpa$Csize, tree, iter=999)
### Almost no difference

library(phytools)
size_phy1 <- contMap(tree, gpa$Csize, legend=FALSE)

### Second test - ordered by centroid size ###
size_reordered <- order(gpa$Csize, decreasing=F)
size2 <- gpa$Csize[size_reordered]
tree2 <- tree
tree2$tip.label <- names(size2)
shape2 <- gpa$coords[,,size_reordered]
dimnames(shape2)[[3]]
phy2 <- procD.pgls(shape2~size2, phy=tree2, iter=999); summary(phy2)
physignal(size2, tree2, iter=999)
size_phy2 <- contMap(tree2, size2, legend=FALSE)

### third test - ordered by middle ###
size_reordered2 <- c(28,20,7,6,19,31,14,16,4,11,17,13,23,10,27,8,25,9,24,3,29,18,21,5,26,1,12,30,22,15,32,2)
size3 <- gpa$Csize[size_reordered2]
tree3 <- tree
tree3$tip.label <- names(size3)
shape3 <- gpa$coords[,,size_reordered2]
dimnames(shape3)[[3]]
phy3 <- procD.pgls(shape3~size3, phy=tree3, iter=999); summary(phy3)
physignal(shape3, tree3, iter=999)
size_phy3 <- contMap(tree3, size3, legend=FALSE)



### New test by creating an aberration ###
data2 <- data_mat
data2[2,2,2] <- 2
gpa <- gpagen(data2)
plotOutliers(gpa$coords)
allom <- procD.lm(gpa$coords~gpa$Csize, iter=999); summary(allom) # Rsq = 0.67

plotAllSpecimens(data2)
plotAllSpecimens(gpa$coords)

### First test - random placement ###
phy1 <- procD.pgls(gpa$coords~gpa$Csize, phy=tree, iter=999); summary(phy1)
physignal(gpa$Csize, tree, iter=999) # non-significant now
### Almost no difference

size_phy1 <- contMap(tree, gpa$Csize, legend=FALSE)

### Second test - ordered by centroid size ###
size_reordered <- order(gpa$Csize, decreasing=F)
size2 <- gpa$Csize[size_reordered]
tree2 <- tree
tree2$tip.label <- names(size2)
shape2 <- gpa$coords[,,size_reordered]
dimnames(shape2)[[3]]
phy2 <- procD.pgls(shape2~size2, phy=tree2, iter=999); summary(phy2) # non-significant Rsq = 0.08
physignal(size2, tree2, iter=999)
size_phy2 <- contMap(tree2, size2, legend=FALSE)

### third test - ordered by middle ###
size_reordered2 <- c(28,20,7,6,2,19,31,14,4,11,17,13,16,23,10,27,25,9,24,3,8,29,18,21,26,1,12,30,5,22,15,32)
size3 <- gpa$Csize[size_reordered2]
tree3 <- tree
tree3$tip.label <- names(size3)
shape3 <- gpa$coords[,,size_reordered2]
dimnames(shape3)[[3]]
phy3 <- procD.pgls(shape3~size3, phy=tree3, iter=999); summary(phy3)
match(dimnames(shape3)[[3]], tree3$tip.label)
physignal(shape3, tree3, iter=999)
size_phy3 <- contMap(tree3, size3, legend=FALSE)
