RK4_matrix <- function(a, b, n, A){
  g <- lapply(1:(n+1), function(x){matrix(0, nrow = dim(as.matrix(A(a)))[1], ncol = dim(as.matrix(A(a)))[2])})
  g[[1]] <- diag(dim(as.matrix(A(a)))[1])

  if(a!=b){
    h <- (b-a)/n

    if(h>0){
      for(i in 1:n){

        G1 <- g[[i]]%*%A(a + h * (i-1))*h
        G2 <- (g[[i]]+1/2*G1)%*%A(a + h * (i-1) + 1/2*h)*h
        G3 <- (g[[i]]+1/2*G2)%*%A(a + h * (i-1) + 1/2*h)*h
        G4 <- (g[[i]]+G3)%*%A(a + h*i)*h

        g[[i+1]] <- g[[i]] + 1/6*(G1+2*G2+2*G3+G4)
      }
    }
    if(h<0){
      for(i in 1:n){

        G1 <- A(a + h * (i-1))%*%g[[i]]*h
        G2 <- A(a + h * (i-1) + 1/2*h)%*%(g[[i]]+1/2*G1)*h
        G3 <- A(a + h * (i-1) + 1/2*h)%*%(g[[i]]+1/2*G2)*h
        G4 <- A(a + h*i)%*%(g[[i]]+G3)*h

        g[[i+1]] <- g[[i]] - 1/6*(G1+2*G2+2*G3+G4)
      }
    }
  }
  return(g)
}


#' Calculate the product integral of a matrix function
#'
#' @param start Start time.
#' @param end End time.
#' @param step_size Step size of the grid.
#' @param Lambda A given matrix function.
#'
#' @return The product integral of the given matrix function.
#'
#' @export
prodint <- function(start, end, step_size, Lambda){
  RK4_matrix(start, end, end/step_size, Lambda)}
