################### READING DATA ####################

#Libraries
library(dplyr)
library(Rfast)
library(neuralnet)
library(pso)

#Work directory
setwd("C:/Users/gonza/OneDrive/Documentos/lipschitz")

#Data frame
df <- read.csv("~/lipschitz/ciudades.csv", dec=",")
n <- nrow(df)
x <- df[,2:4]
y <- df[,5]

################### CROSS-VALIDATION ####################

N <- 20
e.standard <- NULL
e.phistandard <- NULL
e.mw <- NULL
e.phimw <- NULL
e.nn <- NULL
e.lf <- NULL
t.standard <- 0
t.phistandard <- 0
t.mw <- 0
t.phimw <- 0
t.nn <- 0
t.lf <- 0

for (iter in 1:N) {
  
  #Training and test
  train <- sort(sample(1:n,round(0.7*n)))
  test <- setdiff(1:n,train)
  
  #Scaling
  maxs <- as.vector(apply(x[train,], 2, max))
  mins <- as.vector(apply(x[train,], 2, min))
  x_scaled_train <- data.frame(scale(x[train,], mins, maxs-mins))
  x_scaled_test <- data.frame(scale(x[test,], mins, maxs-mins))
  y_min <- min(y[train])
  y_max <- max(y[train])
  y_scaled_train <- (y[train]-y_min)/(y_max-y_min)
  
  #Matrices of distances
  t0 <- Sys.time()
  
  D <- Dist(x_scaled_train)
  Dind <- Dist(y_scaled_train)
  
  t1 <- Sys.time()
  t.dist <- t1-t0
  
  #Optimization of phi
  t0 <- Sys.time()
  
  A <- matrix(0, length(train), length(train))
  for (i in 1:length(train)) {
    for (j in i:length(train)) {
      A[i,j] <- y_scaled_train[train[i]] + y_scaled_train[train[j]]
    }
  }
  A <- A + t(A)
  
  J <- function(p) {
    
    phi <- function(x) {
      x <- sqrt(x)
      func <- p[1]^2*x + p[2]^2*log(1+x) + p[3]^2*atan(x) + p[4]^2*x/(1+x)
      return(func)
    }
    
    K <- max(Dind/phi(D), na.rm = TRUE)
    Q <- max(phi(D)/A, na.rm = TRUE)
    return(K*Q)
  }
  optimum <- psoptim(c(1,0,0,0), J, 
                     control = list(maxit = 100, maxf = 200))
  
  p <- optimum$par
  phi <- function(x) {
    x <- sqrt(x)
    func <- p[1]^2*x + p[2]^2*log(1+x) + p[3]^2*atan(x) + p[4]^2*x/(1+x)
    return(func)
  }
  
  L <- max(Dind/D, na.rm=TRUE)
  K <- max(Dind/phi(D), na.rm = TRUE)
  a0 <- which.min(y_scaled_train)
  
  t1 <- Sys.time()
  t.opt <- t1-t0
  
  #Standard Lipschitz
  t0 <- Sys.time()
  
  standard <- L*dista(x_scaled_train[a0,],x_scaled_test)
  standard <- standard*(y_max-y_min)+y_min
  e.standard[iter] <- sqrt(sum((standard-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.standard <- t.standard + t1-t0 + t.dist
  
  #Standard phi-Lipschitz
  t0 <- Sys.time()
  
  standard <- K*phi(dista(x_scaled_train[a0,],x_scaled_test))
  standard <- standard*(y_max-y_min)+y_min
  e.phistandard[iter] <- sqrt(sum((standard-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.phistandard <- t.phistandard + t1-t0 + t.dist + t.opt
  
  #McShane-Whitney Lipschitz
  t0 <- Sys.time()
  
  M <- dista(x_scaled_test,x_scaled_train)
  v <- t(replicate(nrow(M),y_scaled_train))
  Whitney <- apply(v + L*M, 1, min)*(y_max-y_min)+y_min
  McShane <- apply(v - L*M, 1, max)*(y_max-y_min)+y_min
  dif <- Whitney - McShane
  alpha <- sum(dif*(Whitney-y[test]))/sum(dif^2)
  I <- (1-alpha)*Whitney + alpha*McShane
  e.mw[iter] <- sqrt(sum((I-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.mw <- t.mw + t1-t0 + t.dist
  
  #McShane-Whitney phi-Lipschitz
  t0 <- Sys.time()
  
  M <- dista(x_scaled_test,x_scaled_train)
  v <- t(replicate(nrow(M),y_scaled_train))
  Whitney <- apply(v + K*phi(M), 1, min)*(y_max-y_min)+y_min
  McShane <- apply(v - K*phi(M), 1, max)*(y_max-y_min)+y_min
  dif <- Whitney - McShane
  alpha <- sum(dif*(Whitney-y[test]))/sum(dif^2)
  I <- (1-alpha)*Whitney + alpha*McShane
  e.phimw[iter] <- sqrt(sum((I-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.phimw <- t.phimw + t1-t0 + t.dist + t.opt
  
  #Neural Net
  t0 <- Sys.time()
  
  nn <- neuralnet(y_scaled_train ~ walk + transit + bike, 
                  data.frame(x_scaled_train, y_scaled_train))
  pr.nn <- compute(nn, x_scaled_test)
  pr.nn <- pr.nn$net.result*(y_max-y_min)+y_min
  e.nn[iter] <- sqrt(sum((pr.nn-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.nn <- t.nn + t1-t0

  #Linear fit
  t0 <- Sys.time()
  
  fit <- lm(y_scaled_train ~ walk + transit + bike,
            data.frame(x_scaled_train, y_scaled_train))
  pr.lf <- apply(fit$coefficients[1]+fit$coefficients[2:4]*x_scaled_test,1,sum)
  pr.lf <- pr.lf*(y_max-y_min)+y_min
  e.lf[iter] <- sqrt(sum((pr.lf-y[test])^2)/length(test))
  
  t1 <- Sys.time()
  t.lf <- t.lf + t1-t0
  
  
}

#Summary
summary(e.standard)
summary(e.phistandard)
summary(e.mw)
summary(e.phimw)
summary(e.nn)
summary(e.lf)

#Times
t.standard/N
t.phistandard/N
t.mw/N
t.phimw/N
t.nn/N
t.lf/N

################### CANADA RANKING ####################

#Canada data frame
df.can <- read.csv("~/lipschitz/ciudadescan.csv", dec=",")
scaled.usa <- (df[2:4]-min(df[2:4]))/(max(df[2:4])-min(df[2:4]))
index <- df[,5]
index.scaled <- (index-min(index))/(max(index)-min(index))
scaled.can <- (df.can[2:4]-min(df[2:4]))/(max(df[2:4])-min(df[2:4]))

D <- Dist(scaled.usa)
Dind <- Dist(index.scaled)

A <- matrix(0, nrow(df), nrow(df))
for (i in 1:nrow(df)) {
  for (j in i:nrow(df)) {
    A[i,j] <- index.scaled[i] + index.scaled[j]
  }
}
A <- A + t(A)

J <- function(p) {
  
  phi <- function(x) {
    x <- sqrt(x)
    func <- p[1]^2*x + p[2]^2*log(1+x) + p[3]^2*atan(x) + p[4]^2*x/(1+x)
    return(func)
  }
  
  K <- max(Dind/phi(D), na.rm = TRUE)
  Q <- max(phi(D)/A, na.rm = TRUE)
  return(K*Q)
}
optimum <- psoptim(c(1,0,0,0), J, 
                   control = list(maxit = 100, maxf = 200))

p <- optimum$par
phi <- function(x) {
  x <- sqrt(x)
  func <- p[1]^2*x + p[2]^2*log(1+x) + p[3]^2*atan(x) + p[4]^2*x/(1+x)
  return(func)
}

L <- max(Dind/D, na.rm=TRUE)
K <- max(Dind/phi(D), na.rm = TRUE)
a0 <- which.min(index.scaled)


#Standard phi-Lipschitz
pr.standard <- K*phi(dista(scaled.usa[a0,],scaled.can))
pr.standard <- pr.standard*(max(index)-min(index))+min(index)

#McShane-Whitney phi-Lipschitz
M <- dista(scaled.can,scaled.usa)
v <- t(replicate(nrow(M),index.scaled))
Whitney <- apply(v + K*phi(M), 1, min)*(max(index)-min(index))+min(index)
McShane <- apply(v - K*phi(M), 1, max)*(max(index)-min(index))+min(index)
pr.mw <- 0.5*Whitney + 0.5*McShane


#Neural net
nn <- neuralnet(index.scaled ~ walk + transit + bike, 
                data.frame(scaled.usa, index.scaled))
pr.nn <- compute(nn, scaled.can)
pr.nn <- pr.nn$net.result*(max(index)-min(index))+min(index)

#Linear fit
fit <- lm(index.scaled ~ walk + transit + bike,
          data.frame(scaled.usa, index.scaled))
pr.lf <- apply(fit$coefficients[1]+fit$coefficients[2:4]*scaled.can,1,sum)
pr.lf <- pr.lf*(max(index)-min(index))+min(index)


#Rankings
ranking <- data.frame(df.can$city)
ranking[,2] <- t(pr.standard)
ranking[,3] <- df.can$city
ranking[,4] <- pr.mw
ranking[,5] <- df.can$city
ranking[,6] <- pr.nn
ranking[,7] <- df.can$city
ranking[,8] <- pr.lf

colnames(ranking) <- c("City","Index","City","Index","City","Index","City","Index")
write.csv(ranking, "C:/Users/gonza/OneDrive/Documentos/lipschitz/ranking.csv", 
          row.names=FALSE)