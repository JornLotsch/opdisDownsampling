#' Identify Relevant Variables Based on PCA Projection
#'
#' This function identifies the relevant variables based on the contribution of each
#' variable to the principal components. It uses a reference line calculated as
#' 100/length(contrib), similar to the `fviz_contrib` function from the `factoextra`
#' package.
#'
#' @param res.pca The output of the `prcomp` function, containing the principal
#'   component analysis results.
#'
#' @return A vector of the names of the relevant variables.
#'
relevant_PCAvariables <- function(res.pca) {
  # Define a helper function to calculate the variable coordinates
  var_coord <- function(loadings, comp.sdev) {
    loadings * comp.sdev
  }
  
  # Determine the number of principal components to consider
  nPCs <- sum(res.pca$sdev >= 1)
  
  # Calculate the variable coordinates, cosine squared, and contributions
  loadings <- res.pca$rotation
  sdev <- res.pca$sdev
  var.coord <- t(apply(loadings, 1, var_coord, sdev))
  var.cos2 <- var.coord^2
  comp.cos2 <- apply(var.cos2, 2, sum)
  var.contrib <- data.frame(t(apply(var.cos2, 1, function(x, y) x * 100 / y, comp.cos2)))
  
  # Determine the reference line for important variables
  LowerLimit <- 100 / nrow(var.contrib)
  
  # Identify the important variables for each principal component
  vars.contrib.important <- apply(var.contrib[1:nPCs], 2, function(x) which(x > LowerLimit))
  
  # Combine the important variables across all principal components
  if (nPCs > 1) {
    vars.contrib.important.all <- unique(unlist(lapply(vars.contrib.important, names)))
  } else {
    vars.contrib.important.all <- unique(rownames(var.contrib))
  }
  
  return(vars.contrib.important.all)
}
