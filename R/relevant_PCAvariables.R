#Identifies relevant variables based on PCA projection
#uses  reference line as in fviz_contrib from library factoextra
#calculated as 100/length(contrib)
relevant_PCAvariables <- function(res.pca) {
  var_coord <- function(loadings, comp.sdev) {
    loadings * comp.sdev
  }
  nPCs <- length(res.pca$sdev[res.pca$sdev >= 1])
  loadings <- res.pca$rotation
  sdev <- res.pca$sdev
  var.coord <- t(apply(loadings, 1, var_coord, sdev))
  var.cos2 <- var.coord ^ 2
  comp.cos2 <- apply(var.cos2, 2, sum)
  contrib <- function(var.cos2, comp.cos2) {
    var.cos2 * 100 / comp.cos2
  }
  var.contrib <- data.frame(t(apply(var.cos2, 1, contrib, comp.cos2)))
  LowerLimit <- 100 / nrow(var.contrib)
  vars.contrib.important <- (apply(var.contrib[1:nPCs], 2, function(x) which(x > LowerLimit)))
  if (nPCs > 1) {
    vars.contrib.important.all <- unique(unlist(lapply(vars.contrib.important, names)))
  } else {
    vars.contrib.important.all <- unique(rownames(vars.contrib.important))
  }
  return(vars.contrib.important.all)
}
