#' Compute the conditional Aalen-Johansen estimator.
#' 
#' @param x A numeric value for conditioning.
#' @param data A list of trajectory data for each individual.
#' @param X A vector of covariates.
#' @param p An integer representing the number of states. The absorbing state is last.
#' @param a A bandwidth. Default uses an asymmetric version using alpha.
#' @param alpha A probability around the point x, for asymmetric sub-sampling.
#' 
#' @return A list containing the Aalen-Johansen estimator, the Nelson-Aalen estimator, and related quantities.
#' @export
#' 
#' @examples
aalen_johansen <- function(x, data, X = numeric(0), a = NA, p = NA, alpha = 0.05){
  
  # Get the relevant data by filtering for rows where (x - X)/a is within -1/2 and 1/2
  n <- length(data)
  if(length(X) != n){relevant_data <- data}else{
    if(is.na(a)){
      pvl <- ecdf(X)(x)
      upper <- quantile(X, pvl + alpha/2)
      lower <- quantile(X, pvl - alpha/2)
      relevant_data <- which((X<=upper)&(X>=lower))
    }else{
      relevant_data <- which(((x - X)/a<=1/2)&((x - X)/a>=-1/2))
    }
  }
  data_x <- data[relevant_data]
  prop <- length(relevant_data)/n
  n_x <- length(data_x)
  if(is.na(p)) p <- max(unique(unlist(lapply(data_x, function(Z) Z$states)))) - 1

  # Initialize output lists
  out <- out2 <- list()
  
  # Extract sojourn times, cumulative sojourn times, and state transition information
  R_times <- unlist(lapply(data_x, FUN = function(Z) tail(Z$times,1)))
  t_pool <- unlist(lapply(data_x, FUN = function(Z) Z$times[-1]))
  
  #R_times <- unlist(lapply(data_x, FUN = function(Z) sum(Z$sojourns)))
  #t_pool <- unlist(lapply(data_x, FUN = function(Z) cumsum(Z$sojourns)))
  individuals <- c()
  for(i in 1:n_x){
    individuals <- c(individuals,rep(i,length(data_x[[i]]$times[-1])))
  }
  jumps_pool <- matrix(NA,0,2)
  for(i in 1:n_x){
    v <- data_x[[i]]$states
    jumps_pool <- rbind(jumps_pool,cbind(rev(rev(v)[-1]),v[-1]))
  }
  
  # Sort the times, individuals, and transitions by time
  order_of_times <- order(t_pool)
  ordered_times <- c(0,t_pool[order_of_times])
  ordered_individuals <- c(NA,individuals[order_of_times])
  ordered_jumps <- jumps_pool[order_of_times,]
  ordered_N <- rep(list(matrix(0,p+1,p+1)),nrow(ordered_jumps))
  for(i in 1:nrow(ordered_jumps)){
    ordered_N[[i]][ordered_jumps[i,][1],ordered_jumps[i,][2]] <- 1
  }
  ordered_N <- lapply(ordered_N, FUN = function(Z) Z - diag(diag(Z)))
  colsums_of_N <- lapply(ordered_N,function(N)colSums(N-t(N)))
  
  # Compute out and out2
  out[[1]] <- ordered_N[[1]]*n^{-1}/prop
  for(tm in 2:(length(ordered_times)-1)){
    out[[tm]] <- out[[tm - 1]] + ordered_N[[tm]]*n^{-1}/prop
  }
  #
  decisions <- ordered_times %in% R_times
  out2[[1]] <- colsums_of_N[[1]] * n^{-1}/prop
  for(tm in 2:(length(ordered_times)-1)){
    out2[[tm]] <- out2[[tm - 1]] + colsums_of_N[[tm]] * n^{-1}/prop
    if(decisions[tm]){
      wch <- ordered_individuals[tm]
      end_state <- as.numeric(1:(p+1) == tail(data_x[[wch]]$states,1))
      out2[[tm]] <- out2[[tm]] - end_state * n^{-1}/prop
    }
  }
  
  # Extract initial status for each individual
  I_initial <-  lapply(data_x, FUN = function(Z) as.numeric(1:(p+1) == head(Z$states,1)))
  
  # Compute initial rate for all individuals
  I0 <- Reduce("+", I_initial) * n^{-1}/prop
  
  # Compute rates over time
  It <- lapply(out2, FUN = function(N) return(I0 + N) )
  
  # Initialize list for increments
  increments <- list()
  
  # Compute increments for each time point
  increments[[1]] <- out[[1]]
  for(i in 2:length(out)){
    increments[[i]] <- out[[i]] - out[[i - 1]]
  }
  
  # Compute contribution of each time point
  contribution_first <- increments[[1]]/I0
  contribution_first[is.nan(contribution_first)] <- 0
  contributions <- mapply(FUN = function(a,b){
    res <- b/a
    res[is.nan(res)] <- 0
    return(res)
  }, It[-length(It)], increments[-1],SIMPLIFY = FALSE)
  
  # Compute cumulative sum of contributions
  cumsums <- list()
  cumsums[[1]] <- contribution_first
  for(i in 2:length(out)){
    cumsums[[i]] <- contributions[[i-1]] + cumsums[[i - 1]]
  }

  # Compute the Aalen-Johansen estimator using difference equations
  aj <- list()
  aj[[1]] <- I0
  Delta <- contribution_first
  aj[[2]] <- aj[[1]] + as.vector(aj[[1]] %*% Delta) - aj[[1]] * rowSums(Delta)
  for(i in 2:length(out)){
    Delta <- contributions[[i-1]]
    aj[[i + 1]] <- aj[[i]] + as.vector(aj[[i]] %*% Delta) - aj[[i]] * rowSums((Delta))
  }
  
  # Final touch: adding the diagonal to the Nelson-Aalen estimator
  cumsums <- append(list(matrix(0,p+1,p+1)),lapply(cumsums, FUN = function(M){M_out <- M; diag(M_out) <- -rowSums(M); M_out}))

  # Return output as a list
  return(list(p = aj, Lambda = cumsums, N = out, I0 = I0, It = It, t = ordered_times))
}