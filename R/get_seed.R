get_seed <- function(){
    #Retrieve the current random number seed the RNGkind() for reproducibility. 
    if(!exists('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)) runif(1L)
    seed=get('.Random.seed', envir=globalenv(), mode='numeric', inherits=FALSE)
    attr(seed, 'RNGkind')=RNGkind()
    seed
  }
