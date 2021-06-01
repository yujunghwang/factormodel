#' @title makeDummy
#'
#' @description This function is to make dummy variables using a discrete variable.
#'
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#'
#' @param tZ An input vector
#'
#' @return Returns dZ, a matrix of size length(tZ)-by-card(tZ) : \describe{
#'      The ij-th element in dZ is 1 if tZ\[i\] is equal to the j-th largest value of tZ.
#'      And the ij-th element in DZ is 0 otherwise.
#'      The row sum of dZ must be 1 by construction.
#'      }
#'
#' @examples
#' makeDummy(c(1,2,3))
#'
#' @export
makeDummy <- function(tZ){

  if (!is.vector(tZ) & !is.matrix(tZ) & !is.data.frame(tZ) & !is.array(tZ)){
    stop("Please provide the right input for makeDummy.")
  }

  tZ <- as.vector(tZ)

  # length of tZ
  n1 <- length(tZ)
  # zval is a vector for tZ values
  zval <- sort(unique(tZ))

  nz <- length(zval)
  dZ <- rep(0,n1*nz)
  dim(dZ) <- c(n1,nz)
  for (i1 in 1:nz){
    dZ[tZ==zval[i1],i1] <- 1
  }
  dZ[is.na(tZ)==1,] <- NA
  return(dZ)
}
