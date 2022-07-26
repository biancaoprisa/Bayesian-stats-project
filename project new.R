require(osmar)
require(geosphere)
require(ggmap)
require(mvtnorm)
require(coda)
require(reshape2)
require(MASS)
require(ggplot2)
require(osmdata)
require(magrittr)
require(sf)


register_google(key="AIzaSyAnGAw53RJtWLiJ76HKtnxarZK2VR4GL5s",write=TRUE)
map <- get_map(c(-0.70956,52.04623),zoom=11,maptype="road")
ggmap(map)

landmarks<-data.frame(lon=c(-0.7449436,-1.0074497, -0.6300633),lat=c(51.99919, 52.1318862,52.0694635))
alpha <- bearing(c(-0.70956,52.04623), c(landmarks[1,1], landmarks[2,1] ))
beta <- bearing(c(-0.70956,52.04623), c(landmarks[2,1], landmarks[2,2] ))
gamma <- bearing(c(-0.70956,52.04623), c(landmarks[3,1], landmarks[2,3] ))

d <- seq(0,3,0.0001) # Length of the line
line1 <- data.frame(lon=landmarks[1,1] + d*sin(alpha*pi/180+pi),
                    lat=landmarks[1,2] + d*cos(alpha*pi/180+pi))
line2 <- data.frame(lon=landmarks[2,1] + d*sin(beta*pi/180+pi),
                    lat=landmarks[2,2] + d*cos(beta*pi/180+pi))
line3 <- data.frame(lon=landmarks[3,1] + d*sin(gamma*pi/180+pi),
                    lat=landmarks[3,2] + d*cos(gamma*pi/180+pi))

map <- get_map(c(-0.70956,52.04623),zoom=11,maptype="watercolor")
mapPlot <- ggmap(map)+
  geom_point(aes(x = lon, y = lat), size = 1, data = landmarks, alpha = .5) +
  geom_line(aes(x=lon,y=lat),data=line1) +
  geom_line(aes(x=lon,y=lat),data=line2) +
  geom_line(aes(x=lon,y=lat),data=line3)
mapPlot

intersectBearings <- function(p1,b1,p2,b2) {
  x1 <- p1[1]
  x2 <- p1[1] + 0.1*sin(b1*pi/180)
  x3 <- p2[1]
  x4 <- p2[1] + 0.1*sin(b2*pi/180)
  y1 <- p1[2]
  y2 <- p1[2] + 0.1*cos(b1*pi/180)
  y3 <- p2[2]
  y4 <- p2[2] + 0.1*cos(b2*pi/180)
  x <- ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
  y <- ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
  return(as.numeric(c(x,y)))
}


intersectBearings(landmarks[1,],alpha,landmarks[2,],beta)

intersection <- intersectBearings(landmarks[1,],alpha,landmarks[2,],beta)


bearings <- function(lambda1, phi1, lambda2, phi2) {
  (360 + atan2(phi2 - phi1, lambda2 - lambda1) * 180/pi) %% 360
}

likelihood <- function(theta, alpha, beta, gamma) {
  if (theta[3] <= 0){
    return (log(0))
  }
  else {
    return (dnorm(alpha, bearings(theta[1], theta[2],landmarks[1,1], landmarks[1,2]), sqrt(theta[3]), log=T) +
              dnorm(beta, bearings(theta[1], theta[2],landmarks[2,1], landmarks[2,2]), sqrt(theta[3]), log=T) + 
              dnorm(gamma, bearings(theta[1], theta[2],landmarks[3,1], landmarks[3,2]), sqrt(theta[3]), log=T))
  }
}

prior <- function(theta) {
  sum(dunif(theta[1], -1.2, -0.2, log = T) + 
        dunif(theta[2], 51.5, 52.5, log=T) + 
        dexp(theta[3], 20, log=T))
}


draws <- array(0,dim=c(4000,3,3))
draws[4000,1,] <- runif(3,intersection[1] - 0.01, intersection[1] + 0.01)
draws[4000,2,] <- runif(3,intersection[2] - 0.01, intersection[2] + 0.01)
draws[4000,3,] <- rexp(3,20)

prop.cov <- c(1e-8,1e-8,1e-4)*diag(3)

converged <- FALSE
while (!converged) {
  draws[1,,] <- draws[4000,,]
  accepted <- 1
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed <- rmvnorm(1,draws[step-1,,chain],prop.cov)
      r <- likelihood(proposed, alpha, beta, gamma) + prior(proposed) -
        (likelihood(draws[step-1,,chain], alpha, beta, gamma ))-prior(draws[step-1,,chain])
      a <- min(0,r)
      u <- runif(1)
      if (log(u) < a) {
        draws[step,,chain] <- proposed
        accepted <- accepted + 1
      } else {
        draws[step,,chain] <- draws[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted/120))
  chainlist <- mcmc.list(Chain1=mcmc(draws[,,1]),
                         Chain2=mcmc(draws[,,2]),
                         Chain3=mcmc(draws[,,3]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
  if (!converged) {
    if (accepted/120 <= 20) {
      prop.cov <- prop.cov * 0.8
    }
    if (accepted/120 >= 60) {
      prop.cov <- prop.cov * 1.5
    }
  }
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}

s1 <- rbind(draws[,,1], draws[,,2], draws[,,3])
vc <- cov(s1)
vc

converged <- FALSE
while (!converged) {
  draws[1,,] <- draws[4000,,]
  accepted <- 1
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed <- rmvnorm(1,draws[step-1,,chain],vc)
      r <- likelihood(proposed, alpha, beta, gamma) + prior(proposed) -
        (likelihood(draws[step-1,,chain], alpha, beta, gamma ))-prior(draws[step-1,,chain])
      a <- min(0,r)
      u <- runif(1)
      if (log(u) < a) {
        draws[step,,chain] <- proposed
        accepted <- accepted + 1
      } else {
        draws[step,,chain] <- draws[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted/120))
  chainlist <- mcmc.list(Chain1=mcmc(draws[,,1]),
                         Chain2=mcmc(draws[,,2]),
                         Chain3=mcmc(draws[,,3]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
  if (!converged) {
    if (accepted/120 <= 20) {
      vc <- vc * 0.8
    }
    if (accepted/120 >= 60) {
      vc <- vc * 1.2
    }
  }
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}

#since convergence is achieved we begin collecting the final sample from the posterior

sample <- array(0,dim=c(10000,3,3))
sample[1,,] <- draws[4000,,]
accepted <- 1

for (step in 2:10000) {
  for (chain in 1:3) {
    proposed <- rmvnorm(1,sample[step-1,,chain],vc)
    r <- likelihood(proposed, alpha, beta, gamma) + prior(proposed) -
      (likelihood(sample[step-1,,chain], alpha, beta, gamma ))-prior(sample[step-1,,chain])
    a <- min(0,r)
    u <- runif(1)
    if (log(u) < a) {
      sample[step,,chain] <- proposed
      accepted <- accepted + 1
    } else {
      sample[step,,chain] <- sample[step-1,,chain]
    }
  }
}
print(sprintf("Acceptance rate: %f",accepted/300))

#distance to a road

mean(sample[,2,])
mean(sample[,1,])
mean(sample[,3,])

D <- kde2d(as.vector(sample[,1,]),as.vector(sample[,2,]),
           h=c(sd(sample[,1,]),sd(sample[,2,])),
           n=1024,
           lims=c(-1.5, 2.5,51, 55)) # Enough to cover map
z <- melt(D$z)
z$Var1<-D$x[z$Var1]
z$Var2<-D$y[z$Var2]
map <- get_map(c(mean(sample[,1,]),mean(sample[,2,])),zoom=13,maptype="road")
mapPoints <- ggmap(map)+
  geom_point(aes(x = lon, y = lat), size = 1, data = landmarks, alpha = .5) +
  geom_raster(data=z,aes(x=Var1,y=Var2,fill=value))+
  guides(fill=FALSE,alpha=FALSE)+
  scale_fill_gradientn(colours=c("#0000FF00","#0000FFFF"))+
  coord_cartesian() +
  geom_line(aes(x=lon,y=lat),data=line1) +
  geom_line(aes(x=lon,y=lat),data=line2) +
  geom_line(aes(x=lon,y=lat),data=line3)+
  geom_point(aes(x=lon,y=lat),
             data=data.frame(lon=mean(sample[,1,]),lat=mean(sample[,2,])),
             size=0.5,colour="#FF0000")
mapPoints


