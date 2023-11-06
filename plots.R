library(plotly)

d <- function(u,v) {
  distancia <- (u-v)^2
  return(sqrt(distancia))
}

d0 <- function(u,v) {
  distancia <- d(u,0) + d(v,0)
  return(distancia)
}

# \phi función
phi <- function(t) {
  return(log(1+t))
}

dphi <- function(u,v) {
  distancia <- phi(d(u,v))
  return(distancia)
}

d0phi <- function(u,v) {
  distancia <- phi(d(u,0)) + phi(d(v,0))
  return(distancia)
}

y <- x <- seq(-5,5,0.1)
z1 <- outer(x, y, d)
z2 <- outer(x, y, d0)
z3 <- outer(x, y, dphi)
z4 <- outer(x, y, d0phi)

# Triangular inequality of d
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~z1, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)")))
fig <- fig %>% add_surface(z = ~z2, colorscale = list(c(0,1),c("rgb(107,184,255)","rgb(0,90,124)")))
fig

# Triangular inequality of d_\phi
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~z3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)")))
fig <- fig %>% add_surface(z = ~z4, colorscale = list(c(0,1),c("rgb(107,184,255)","rgb(0,90,124)")))
fig

# d and d_phi
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~z1, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)")))
fig <- fig %>% add_surface(z = ~z3, colorscale = list(c(0,1),c("rgb(107,184,255)","rgb(0,90,124)")))
fig