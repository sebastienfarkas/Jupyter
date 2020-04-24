### TP de mod??lisation stochastique - s??ance 1 ---------------------------

# PARTIE 1: Simulation de variables al??atoires uniformes -----------------

#(b)

randu = function(n,seed=1,m=2^31-1,a=65539,b=0){
  #Les valeurs seed=1,m=2^31-1,a=65539,b=0 sont des valeurs par d??faut
  #La fonction peut ainsi ??tre appel??e uniquement avec le param??tre n
  
  #Cr??ation d'un vecteur
  u = double(length = n)
  #Affectation de la valeur initiale (Y1)
  u[1]=seed
  #Calcul des n-1 termes suivant de la suite (Yn) par r??currence
  for(i in 1:(n-1)){
    u[i+1]<-(a*u[i]) %% m
  }
  #Division pour m pour obtenir les n premiers termes de la suite (Xn)
  u=u/m
  return(u)
}

#(c)

par(mfrow=c(1,3))
#Simulation des r??alisations avec le g??n??rateur RANDU
n <- 10^2
u <- randu(n)

#V??rification de l'hypoth??se (identiquement distribu??)
hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur RANDU \n avec n=100")
abline(h=1,col="red")

n <- 10^3
u <- randu(n)

#V??rification de l'hypoth??se (identiquement distribu??)
hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur RANDU\n avec n=1000")
abline(h=1,col="red")

n <- 10^4
u <- randu(n)

#V??rification de l'hypoth??se (identiquement distribu??)
hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur RANDU\n avec n=10000")
abline(h=1,col="red")

#V??rification de l'hypoth??se (ind??pendant)
library(rgl)
n <- 10^4
u <- randu(n+2)
plot3d(cbind(u[1:n], u[2:(n+1)], u[3:(n+2)]))

#(d)

#Simulation des r??alisations avec le g??n??rateur RUNIF
n <- 10^2
u <- runif(n)
#V??rification de l'hypoth??se (identiquement distribu??)

hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur runif\n avec n=100")
abline(h=1,col="red")

#Simulation des r??alisations avec le g??n??rateur RUNIF
n <- 10^3
u <- runif(n)
#V??rification de l'hypoth??se (identiquement distribu??)

hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur runif\n avec n=1000")
abline(h=1,col="red")

#Simulation des r??alisations avec le g??n??rateur RUNIF
n <- 10^4
u <- runif(n)
#V??rification de l'hypoth??se (identiquement distribu??)

hist(u,
     freq=FALSE,
     xlab="x",
     breaks=20,
     ylim=c(0,1.5),
     ylab="Frequence d'observations",
     main="Histogramme des realisations uniformes\n sur [0;1] selon le generateur runif\n avec n=10000")
abline(h=1,col="red")

#V??rification de l'hypoth??se (ind??pendant)
n <- 10^4
u <- runif(n+2)
plot3d(cbind(u[1:n], u[2:(n+1)], u[3:(n+2)]))


# PARTIE 2: M??thodes g??n??rales de simulation de variables al??atoires -----

exponentielle = function(n,lambda){
  exp = -log(runif(n))/lambda
  return(exp)
}

par(mfrow=c(3,2))

#Simulation des r??alisations
n <- 10^2
exp <- exponentielle(n,lambda=1)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 5 via le generateur runif\n avec n=100")
lines(seq(from = 0,to = 4,length.out = 500),
      dexp(seq(from = 0,to = 4,length.out = 500),
           rate = 1),col="red")

#Simulation des r??alisations

exp <- exponentielle(n,lambda=10)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 10 via le generateur runif\n avec n=100")
lines(seq(from = 0,to = 4,length.out = 500),dexp(seq(from = 0,to = 4,length.out = 500),rate = 10),col="red")


#Simulation des r??alisations
n <- 10^3
exp <- exponentielle(n,lambda=1)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 5 via le generateur runif\n avec n=1000")
lines(seq(from = 0,to = 4,length.out = 500),
      dexp(seq(from = 0,to = 4,length.out = 500),
           rate = 1),col="red")

#Simulation des r??alisations

exp <- exponentielle(n,lambda=10)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 10 via le generateur runif\n avec n=1000")
lines(seq(from = 0,to = 4,length.out = 500),dexp(seq(from = 0,to = 4,length.out = 500),rate = 10),col="red")


#Simulation des r??alisations
n <- 10^4
exp <- exponentielle(n,lambda=1)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 5 via le generateur runif\n avec n=10000")
lines(seq(from = 0,to = 4,length.out = 500),
      dexp(seq(from = 0,to = 4,length.out = 500),
           rate = 1),col="red")

#Simulation des r??alisations

exp <- exponentielle(n,lambda=10)

hist(exp,
     freq=FALSE,
     breaks=20,
     xlim=c(0,4),
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations exponentielle\n de parametre 10 via le generateur runif\n avec n=10000")
lines(seq(from = 0,to = 4,length.out = 500),dexp(seq(from = 0,to = 4,length.out = 500),rate = 10),col="red")


#(b)

laplace=function(n){
  E=exponentielle(n = n,lambda = 1)
  U=runif(n=n)
  E[U<1/2]=-E[U<1/2]
  return(E)
}

dlaplace=function(x){
  return(exp(-abs(x))/2)
}

par(mfrow=c(1,3))
hist(x = laplace(10^2),
     xlim=c(-10,10),
     freq = FALSE,
     breaks = 50,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi de Laplace\n via la méthode des mélanges de modèles\n avec n=100")
lines(seq(from = -10,to = 10,length.out = 200),dlaplace(seq(from = -10,to = 10,length.out = 200)),col="red")

hist(x = laplace(10^3),
     xlim=c(-10,10),
     freq = FALSE,
     breaks = 50,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi de Laplace\n via la méthode des mélanges de modèles\n avec n=1000")
lines(seq(from = -10,to = 10,length.out = 200),dlaplace(seq(from = -10,to = 10,length.out = 200)),col="red")

hist(x = laplace(10^4),
     xlim=c(-10,10),
     freq = FALSE,
     breaks = 50,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi de Laplace\n via la méthode des mélanges de modèles\n avec n=10000")
lines(seq(from = -10,to = 10,length.out = 200),dlaplace(seq(from = -10,to = 10,length.out = 200)),col="red")



#Fonction rejet prenant en argument
#n le nombre de simulations d??sir??es
#f la densit?? selon laquelle nous voulons simuler
#gsim la fonction permettant de g??n??rer une r??alisation suivant la loi instrumentale
#g la fonction densit?? de la loi instrumentale
#M la constante multiplicative

reject <- function(n,f,gsim,g,M){
  i <- 1 
  ans <- double(n) 
  while (i <= n){
    U <- runif(n=1)
    Y <- gsim(n=1)
    if (U <= (f(Y) / (M * g(Y)))){
      ans[i] <- Y
      i <- i + 1
    }
  }
  return(ans)
}


aire_reject <- function(n,f,gsim,g,M){
  i <- 1
  j=1
  ans <- matrix(NA, n, 2)
  while (i <= n){
    U <- runif(n=1)
    Y <- gsim(n=1)
    if (U <= (f(Y) / (M * g(Y)))){
      ans[i,] <- c(Y,U*M*g(Y))
      i <- i + 1
    }
    j=j+1
  }
  print(i/j)
  return(ans)
}

perf_reject <- function(n,f,gsim,g,M){
  i <- 1
  j=1
  ans <- matrix(NA, n, 2)
  while (i <= n){
    U <- runif(n=1)
    Y <- gsim(n=1)
    if (U <= (f(Y) / (M * g(Y)))){
      ans[i,1]=i
      ans[i,2]=i/j
      i <- i + 1
    }
    j=j+1
  }
  return(ans)
}


dnorm01=function(x){
  return((1/sqrt(2*pi))*exp(-(x^2)/2))
}

Mlaplacenormal=sqrt(2*exp(1))/sqrt(pi)

par(mfrow=c(2,3))

n = 10^2
rnorm_reject=reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
hist(rnorm_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du rejet\n avec n=100")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^3
rnorm_reject=reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
hist(rnorm_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du rejet\n avec n=1000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^4
rnorm_reject=reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
hist(rnorm_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du rejet\n avec n=10000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")


n = 10^2
airenorm_reject=aire_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(-10,10),
     ylim=c(0,0.7),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du rejet\n avec n=100")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^3
airenorm_reject=aire_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(-10,10),
     ylim=c(0,0.7),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du rejet\n avec n=1000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^4
airenorm_reject=aire_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal)
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(-10,10),
     ylim=c(0,0.7),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du rejet\n avec n=10000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^5
par(mfrow=c(1,1))
plot(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
     type='l',
     ylim=c(1/Mlaplacenormal-0.1,1/Mlaplacenormal+0.1),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode du rejet vers 1/M = 0.76")
abline(h=1/Mlaplacenormal,col="red")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
     type='l',
     col="purple")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="orange")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="brown")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="blue")

#Fonction ratio prenant en argument
#n le nombre de simulations d??sir??es
#a la borne sup??rieure de Ch suivant la premi??re coordonn??e
#bm la borne inf??rieure de Ch suivant la deuxi??me coordonn??e
#bp la borne sup??rieure de Ch suivant la deuxi??me coordonn??e

ratio_reject <- function(n,h,a,bm,bp){
  ans <- rep(NA, n)
  for (i in 1:n){
    while (TRUE){
      U <- runif(n = 1,min =0,max = a)
      V <- runif(n = 1,min = bm,max = bp)
      if (U <= sqrt(do.call(what = h,args = list(x=V/U)))){
        ans[i] <- V / U
        break
      }
    }
  }
  return(ans)
}

aire_ratio_reject <- function(n,h,a,bm,bp){
  ans <- matrix(NA, n, 2)
  i=0
  j=0
  while(i<n){
    reject=TRUE
    while(reject){
      U <- runif(1, 0, a)
      V <- runif(1, bm, bp)
      if (U <= sqrt(do.call(what = h,args = list(x=V/U)))){
        ans[i,] <- c(U,V)
        reject=FALSE
      }
      j=j+1
    }
    i=i+1
  }
  print(i/j)
  return(ans)
}

perf_ratio_reject <- function(n,h,a,bm,bp){
  ans <- matrix(NA, n, 2)
  i=0
  j=0
  while(i<n){
    reject=TRUE
    while(reject){
      U <- runif(1, 0, a)
      V <- runif(1, bm, bp)
      if (U <= sqrt(do.call(what = h,args = list(x=V/U)))){
        ans[i,1] <- i
        ans[i,2] <- i/j
        reject=FALSE
      }
      j=j+1
    }
    i=i+1
  }
  return(ans)
}

#(d)

h=function(x){return(exp(-(x^2)/2))}

par(mfrow=c(2,3))

n=10^2
rnorm_ratio_reject=ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
hist(rnorm_ratio_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du ratio rejet\n avec n=100")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")

n=10^3
rnorm_ratio_reject=ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
hist(rnorm_ratio_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du ratio rejet\n avec n=1000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")


n=10^4
rnorm_ratio_reject=ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
hist(rnorm_ratio_reject,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode du ratio rejet\n avec n=10000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")

n = 10^2
airenorm_reject=aire_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(0, 1),
     ylim=c(-sqrt(2/exp(1)), sqrt(2/exp(1))),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du ratio rejet\n avec n=100")
abline(v=0,col="red")
abline(v=1,col="red")
abline(h=-sqrt(2/exp(1)),col="red")
abline(h=sqrt(2/exp(1)),col="red")

n = 10^3
airenorm_reject=aire_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(0, 1),
     ylim=c(-sqrt(2/exp(1)), sqrt(2/exp(1))),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du ratio rejet\n avec n=1000")
abline(v=0,col="red")
abline(v=1,col="red")
abline(h=-sqrt(2/exp(1)),col="red")
abline(h=sqrt(2/exp(1)),col="red")

n = 10^4
airenorm_reject=aire_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1)))
plot(airenorm_reject[,1],
     airenorm_reject[,2],
     xlim=c(0, 1),
     ylim=c(-sqrt(2/exp(1)), sqrt(2/exp(1))),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode du ratio rejet\n avec n=10000")
abline(v=0,col="red")
abline(v=1,col="red")
abline(h=-sqrt(2/exp(1)),col="red")
abline(h=sqrt(2/exp(1)),col="red")

n = 10^5
par(mfrow=c(1,1))
plot(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
     type='l',
     ylim=c(0.73-0.1,0.73+0.1),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode du ratio rejet")
abline(h=sqrt(pi*exp(1))/4,col="red")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="purple")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="orange")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="brown")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="blue")


boxmuller <-function(n){
  res <- matrix(data = NA, nrow = n, ncol = 2)
  i=0
  j=0
  while(i<n){
    reject=TRUE
    while(reject){
      u=runif(n = 1,min = 0,max = 1)
      v=runif(n = 1,min = 0,max = 1)
      u=2*u-1
      v=2*v-1
      s <- u^2+v^2
      if(s<=1){
        res[i+1,1] = u * sqrt(-2*log(s)/s)
        res[i+1,2] = v * sqrt(-2*log(s)/s)
        reject=FALSE
      }
      j=j+1
    }
    i=i+1
  }
  return(res)
}

aire_boxmuller <-function(n){
  res <- matrix(data = NA, nrow = n, ncol = 2)
  i=0
  j=0
  while(i<n){
    reject=TRUE
    while(reject){
      u=runif(n = 1,min = 0,max = 1)
      v=runif(n = 1,min = 0,max = 1)
      u=2*u-1
      v=2*v-1
      s <- u^2+v^2
      if(s<=1){
        res[i+1,1] = u 
        res[i+1,2] = v
        reject=FALSE
      }
      j=j+1
    }
    i=i+1
  }
  print(i/j)
  return(res)
}

perf_boxmuller <-function(n){
  res <- matrix(data = NA, nrow = n, ncol = 2)
  i=0
  j=0
  while(i<n){
    reject=TRUE
    while(reject){
      u=runif(n = 1,min = 0,max = 1)
      v=runif(n = 1,min = 0,max = 1)
      u=2*u-1
      v=2*v-1
      s <- u^2+v^2
      if(s<=1){
        res[i+1,1] = i
        res[i+1,2] = i/j
        reject=FALSE
      }
      j=j+1
    }
    i=i+1
  }
  return(res)
}


par(mfrow=c(2,3))

n = 10^2
rnorm_boxmuller=boxmuller(n = n)
hist(rnorm_boxmuller,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode boxmuller\n avec n=100")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^3
rnorm_boxmuller=boxmuller(n = n)
hist(rnorm_boxmuller,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode boxmuller\n avec n=1000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")

n = 10^4
rnorm_boxmuller=boxmuller(n = n)
hist(rnorm_boxmuller,
     xlim=c(-10,10),
     ylim=c(0,0.7),
     freq = FALSE,
     breaks = 30,
     xlab="x",
     ylab="Frequence d'observations",
     main="Histogramme des realisations d'une loi Normale\n via la méthode boxmuller\n avec n=10000")
lines(x=seq(-10,10,length.out=1000),y=dnorm01(seq(-10,10,length.out=1000)),col="red")
lines(x=seq(-10,10,length.out=1000),y=Mlaplacenormal*dlaplace(seq(-10,10,length.out=1000)),col="blue")


n = 10^2
airenorm_boxmuller=aire_boxmuller(n = n)
plot(airenorm_boxmuller[,1],
     airenorm_boxmuller[,2],
     xlim=c(-1,1),
     ylim=c(-1,1),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode boxmuller\n avec n=100")

n = 10^3
airenorm_boxmuller=aire_boxmuller(n = n)
plot(airenorm_boxmuller[,1],
     airenorm_boxmuller[,2],
     xlim=c(-1,1),
     ylim=c(-1,1),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode boxmuller\n avec n=1000")

n = 10^4
airenorm_boxmuller=aire_boxmuller(n = n)
plot(airenorm_boxmuller[,1],
     airenorm_boxmuller[,2],
     xlim=c(-1,1),
     ylim=c(-1,1),
     xlab="x",
     ylab="Frequence d'observations",
     main="Représentation des simulations acceptées\n via la méthode boxmuller\n avec n=10000")

n = 10^5
par(mfrow=c(1,1))
plot(perf_boxmuller(n = n),
     type='l',
     ylim=c(0.75,0.85),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode boxmuller")
abline(h=pi/4,col="red")
lines(perf_boxmuller(n = n),
      type='l',
      col="purple")
lines(perf_boxmuller(n = n),
      type='l',
      col="orange")
lines(perf_boxmuller(n = n),
      type='l',
      col="brown")
lines(perf_boxmuller(n = n),
      type='l',
      col="blue")

n = 10^5
par(mfrow=c(1,3))


plot(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
     type='l',
     ylim=c(0.7,0.85),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode du rejet")
abline(h=1/Mlaplacenormal,col="red")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="purple")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="orange")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="brown")
lines(perf_reject(n = n,f = dnorm01,gsim = laplace,g = dlaplace,M = Mlaplacenormal),
      type='l',
      col="blue")


plot(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
     type='l',
     ylim=c(0.7,0.85),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode du ratio rejet")
abline(h=sqrt(pi*exp(1))/4,col="red")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="purple")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="orange")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="brown")
lines(perf_ratio_reject(n=n,h=h,a=1,bm=-sqrt(2/exp(1)),bp=sqrt(2/exp(1))),
      type='l',
      col="blue")


plot(perf_boxmuller(n = n),
     type='l',
     ylim=c(0.7,0.85),
     xlab="Nombre de réalisations acceptées",
     ylab="Taux d'acceptation",
     main="Convergence du taux d'acceptation de la méthode boxmuller")
abline(h=pi/4,col="red")
lines(perf_boxmuller(n = n),
      type='l',
      col="purple")
lines(perf_boxmuller(n = n),
      type='l',
      col="orange")
lines(perf_boxmuller(n = n),
      type='l',
      col="brown")
lines(perf_boxmuller(n = n),
      type='l',
      col="blue")