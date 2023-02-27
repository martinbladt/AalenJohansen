sim_jump <- function(i, t, u = 0, tn, rate, dist, b = NA){
  if(is.na(b)){
    b <- -optimize(function(x){-rate(t+x, u+x)}, c(0, tn - t))$objective #not reliable; better set b
  }
  r <- rexp(1, rate = b)
  s <- t + r
  v <- u + r
  y <- runif(1)
  while(y > rate(s, v)/b){
    r <- rexp(1, rate = b)
    s <- s + r
    v <- v + r
    y <- runif(1)
  }
  pr <- dist(s, v)
  j <- sample(length(pr), 1, prob = pr)
  return(list(time = s, mark = j))
}

#' Simulate the path of a time-inhomogeneous (semi-)Markov process until a maximal time
#' 
#' @param i The initial state, integer.
#' @param t The initial time, numeric.
#' @param u The initial duration (since the last transition), numeric. By default equal to zero
#' @param tn The maximal time, numeric. By default equal to inifinity
#' @param rates The total transition rates out of states, a function with arguments state (integer), time (numeric), and duration (numeric) returning a rate (numeric).
#' @param dists The distribution of marks, a function with arguments state (integer), time (numeric), and duration (numeric) returning a probability vector.
#' @param abs Vector indicating which states are absorbing. By default the last state is absorbing.
#' @param bs Vector of upper bounds on the total transition rates. By default the bounds are determined using optimize, which might only identify a local maximum.
#'
#' @return A list concerning jump times and states, with the first time being the initial time t and state and the last time being tn (if not absorbed)
#' 
#' @export 
#'  
#' @examples
#'  
#' jump_rate <- function(i, t, u){if(i == 1){3*t} else if(i == 2){5*t} else{0}}
#' mark_dist <- function(i, s, v){if(i == 1){c(0, 1/3, 2/3)} else if(i == 2){c(1/5, 0, 4/5)} else{0}}
#' sim <- sim_path(sample(1:2, 1), t = 0, tn = 2, rates = jump_rate, dists = mark_dist)
#' sim

sim_path <- function(i, rates, dists, t = 0, u = 0, tn = Inf, abs = numeric(0), bs = NA){
  times <- t
  marks <- i
  if(length(abs) == 0){
    abs <- c(rep(FALSE, length(dists(i, t, u)) - 1), TRUE)
  }
  while(!abs[tail(marks, 1)]){
    z <- sim_jump(tail(marks, 1), tail(times, 1), u, tn, function(s, v){rates(tail(marks, 1), s, v)}, function(s, v){dists(tail(marks, 1), s, v)}, bs[tail(marks, 1)])
    if(z$time > tn){
      break
    }
    times <- c(times, z$time)
    marks <- c(marks, z$mark)
    u <- 0
  }
  if(!abs[tail(marks, 1)]){
    times <- c(times, tn)
    marks <- c(marks, tail(marks, 1))
  }
  return(list(times = times, states = marks))
}