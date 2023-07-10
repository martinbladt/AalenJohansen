library(AalenJohansen)

set.seed(2)

jump_rate <- function(i, t, u){
  if(i == 1){
    2 / (1+1/2*t)
  } else if(i == 2){
    3 / (1+1/2*t)
  } else{
    0
  }
}

mark_dist <- function(i, s, v){
  if(i == 1){
    c(0, 1/2, 1/2)
  } else if(i == 2){
    c(2/3, 0, 1/3)
  } else{
    0
  }
}

lambda <- function(t){
  A <- matrix(c(2/(1+1/2*t)*mark_dist(1, t, 0), 3/(1+1/2*t)*mark_dist(2, t, 0), rep(0, 3)),
              nrow = 3, ncol = 3, byrow = TRUE)
  diag(A) <- -rowSums(A)
  A
}

#############################################

n <- 10000

R <- runif(n, 0, 5)

sim <- list()
for(i in 1:n){
  sim[[i]] <- sim_path(sample(1:2, 1), rates = jump_rate, dists = mark_dist,
                       tn = R[i], bs = c(2*R[i], 3*R[i], 0))
}

## now we left-truncate
L <- runif(n, 0, R * rbinom(n, 1, 0.5))

sim_trunc <- mapply(FUN = function(Z,l){
  if(l>max(Z$times)) return()
  o <- list()
  o$times <- c(l, Z$times[Z$times > l])
  o$states <- c(Z$states[max(which(Z$times <= l))], Z$states[Z$times > l])
  o
},sim,L,SIMPLIFY = F)

sim_trunc
#sim[[1]]
L[1]
sim_trunc <- sim_trunc[!sapply(sim_trunc,is.null)]


## now we do an augmented model
#p <- max(unlist(lapply(sim_trunc, FUN = function(Z) Z$states)))

sim_aug <- mapply(FUN = function(Z){
  o <- Z
  o$states <- o$states + 1
  if(head(o$times,1)>0){o$times <- c(0,o$times); o$states <- c(1,o$states)}
  o
},sim_trunc, SIMPLIFY = F)

## now the fitting
# Extract initial status for each individual
fit <- aalen_johansen(sim_aug,L_tr = TRUE)

fit$p <- lapply(fit$p, FUN = function(Z) Z[-1])
fit$Lambda <- lapply(fit$Lambda, FUN = function(Z) Z[-1,-1])

## now plot
v1 <- unlist(lapply(fit$Lambda, FUN = function(L) L[2,1]))
v0 <- fit$t
p <- unlist(lapply(fit$p, FUN = function(L) L[2]))
P <- unlist(lapply(prodint(0, 5, 0.01, lambda), FUN = function(L) (c(1/2, 1/2, 0) %*% L)[2]))

par(mfrow = c(1, 2))
par(mar = c(2.5, 2.5, 1.5, 1.5))

plot(v0, v1, type = "l", lty = 2, xlab = "", ylab = "", main = "Hazard")
lines(v0, 4*log(1+1/2*v0))
plot(v0, p, type = "l", lty = 2, xlab = "", ylab = "", main = "Probability")
lines(seq(0, 5, 0.01), P)

