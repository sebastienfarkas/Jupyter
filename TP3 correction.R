# Chaines de Markov -------------------------------------------------------

#A voir : définition périodique / irréductible 
#A changer : la manière dont les questions sont posées


#a) Introduction

transition <-function(x,Q){
  k <- nrow(Q) 
  C2 <- cumsum(Q[x,])
  C1 <- cumsum(c(0,Q[x,1:(k-1)]))
  t <- runif(1)
  y  <- (t > C1) & (t<=C2)
  state <- which(y ==1)
  return(state)
}

P=matrix(c(1/3,1/3,1/3,1/2,0,1/2,0,1/2,1/2),
         nrow=3,byrow = T)

n=20
X1=rep(NA,length.out=n)
X1[1]=1
for (i in 2:n){
  X1[i]=transition(X1[i-1],P)
}
X2=rep(NA,length.out=n)
X2[1]=1	
for (i in 2:n){
  X2[i]=transition(X2[i-1],P)
}
X3=rep(NA,length.out=n)
X3[1]=1	
for (i in 2:n){
  X3[i]=transition(X3[i-1],P)
}

par(mfrow = c(3,1)) 
plot(1:20,X1,col='red',type='s',
     main="Une réalisation de chaîne de Markov définie par \n la matrice P avec 1 comme état initial")
plot(X2,col='blue',type='s',
     main="Une réalisation de chaîne de Markov définie par \n la matrice P avec 1 comme état initial")
plot(X3,col='green',type='s',
     main="Une réalisation de chaîne de Markov définie par \n la matrice P avec 1 comme état initial")


#b) Irreductibilite

R = matrix(c(1/2,1/2,0,0,1/2,1/2,0,0,0,0,2/3,1/3,0,0,1/4,3/4),nrow=4,byrow = T)

n=20
XR1=rep(NA,length.out=n)
XR1[1]=1
for (i in 2:n){
  XR1[i]=transition(XR1[i-1],R)
}
XR2=rep(NA,length.out=n)
XR2[1]=2
for (i in 2:n){
  XR2[i]=transition(XR2[i-1],R)
}
XR3=rep(NA,length.out=n)
XR3[1]=3
for (i in 2:n){
  XR3[i]=transition(XR3[i-1],R)
}
XR4=rep(NA,length.out=n)
XR4[1]=4
for (i in 2:n){
  XR4[i]=transition(XR4[i-1],R)
}

par(mfrow = c(2,2))
plot(XR1,ylim=c(1,4),col='red',type='s',main="Une réalisation de chaîne de Markov définie par \n la matrice R avec 1 comme état initial")
plot(XR2,ylim=c(1,4),col='blue',type='s',main="Une réalisation de chaîne de Markov définie par \n la matrice R avec 2 comme état initial")
plot(XR3,ylim=c(1,4),col='green',type='s',main="Une réalisation de chaîne de Markov définie par \n la matrice R avec 3 comme état initial")
plot(XR4,ylim=c(1,4),col='orange',type='s',main="Une réalisation de chaîne de Markov définie par \n la matrice R avec 4 comme état initial")


#c) Probabilité invariante

diag <- eigen(t(P))
M <- diag$vectors
mu <- M[,1]
t(P) %*% mu - mu
pstat = mu / sum(mu)

pstat

#d) Loi forte des grands nombres


n=1000
indic=matrix(rep(0,3*n),ncol=3)
XP=1
indic[1,1] = 1
for (i in 2:n){
  XP[i]=transition(XP[i-1],P)
  indic[i,XP[i]] =1
}	

emp.mean<- apply(indic,2,cumsum)

par(mfrow=c(1,1))
plot(1:n,emp.mean[,1]/1:n,type="l",ylim=c(0,1))
lines(1:n,emp.mean[,2]/1:n,col="blue")
lines(1:n,emp.mean[,3]/1:n,col="red")
lines(1:n,rep(pstat[1],n))
lines(1:n,rep(pstat[2],n),col="blue")
lines(1:n,rep(pstat[3],n),col="red")


# Study of the distribution of Xn

Emp.distr <- function(n,ntraj){
  Xntraj=rep(0,ntraj)
  for (m in 1:ntraj){  
    XP=rep(1,n)	
    for (i in 2:n){
      XP[i]=transition(XP[i-1],P)
    }
    Xntraj[m] <- XP[n]	  
  }	
  return(table(Xntraj)/ntraj)  
}	


Mat <- matrix(c(Emp.distr(2,1000),Emp.distr(3,1000),Emp.distr(4,1000),
                Emp.distr(5,1000),Emp.distr(10,1000),pstat), byrow = FALSE
              ,ncol = 6)
colnames(Mat) <- c('2','3','4','5','10','pstat')

par(mfrow=c(1,1))
barplot(Mat,beside = TRUE)

# Metropolis Hasting ------------------------------------------------------

#Partie notée

#a) Définition

hastingstraj=function(M,ptinit,alpha){
  traj=matrix(rep(0,2*(M+1)),nrow=M+1)
  traj[1,]=ptinit
  j=0
  for (i in 1:M) {
    y <- traj[i,]+alpha*c(runif(n = 1,min = -1,max = 1),runif(n = 1,min = -1,max = 1))
    rhotilde <- (exp(-(1/2)*(t(y)%*%(y))))/(exp(-(1/2)*(t(traj[i,])%*%(traj[i,]))))
    u <- runif(1)
    if(u<=rhotilde){
      traj[i+1,] <- y
      j=j+1
    }else{
      traj[i+1,] <- traj[i,]
    }
  }
  #print(i/j)
  return(traj)
}

hastingsstep=function(M,ptinit,alpha){
  traj=matrix(rep(0,2*(M+1)),nrow=M+1)
  traj[1,]=ptinit
  step=rep(0,length.out=(M+1))
  step[1]=0
  j=0
  for (i in 1:M) {
    y <- traj[i,]+alpha*c(runif(1,-1,1),runif(1,-1,1))
    rhotilde <- (exp(-(1/2)*(t(y)%*%(y))))/(exp(-(1/2)*(t(traj[i,])%*%(traj[i,]))))
    u <- runif(1)
    if(u<=rhotilde){
      traj[i+1,] <- y
      step[i+1] = step[i]+1
      j=j+1
      
    }else{
      traj[i+1,] <- traj[i,]
      step[i+1] = step[i]
    }
  }
  return(list(traj,step,j/M))
}


#b) Exemple

#Trajectoire des couples our deux états initiaux différents

par(mfrow = c(2,4))	

lev <- seq(0,0.16,length.out = 500)
couleurs <- rainbow(length(lev) -1,start = 0.7,end = 0.95)
grid.dens <- dnorm(seq(-4,4,0.1)) %*% t(dnorm(seq(-4,4,0.1)))
M=2000

alpha=0.01
T1=hastingstraj(M,c(3,3),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 0.01 \n X_0=(3,3)')
points(T1,pch=20,cex = 0.5)

alpha=0.1
T1=hastingstraj(M,c(3,3),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 0.1 \n X_0=(3,3)')
points(T1,pch=20,cex = 0.5)

alpha=1
T2=hastingstraj(M,c(3,3),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 1 \n X_0=(3,3)')
points(T2,pch=20,cex = 0.5)	

alpha=10
T3=hastingstraj(M,c(3,3),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 10 \n X_0=(3,3)')
points(T3,pch=20,cex = 0.5)	

alpha=0.01
T1=hastingstraj(M,c(1,1),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 0.01 \n X_0=(1,1)')
points(T1,pch=20,cex = 0.5)

alpha=0.1
T1=hastingstraj(M,c(1,1),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 0.1 \n X_0=(1,1)')
points(T1,pch=20,cex = 0.5)

alpha=1
T2=hastingstraj(M,c(1,1),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 1 \n X_0=(1,1)')
points(T2,pch=20,cex = 0.5)	

alpha=10
T3=hastingstraj(M,c(1,1),alpha)
image(seq(-4,4,0.1),seq(-4,4,0.1), grid.dens,breaks = lev, col = couleurs,xlab = '',ylab = '',main = 'alpha = 10 \n X_0=(1,1)')
points(T3,pch=20,cex = 0.5)	


#Illustration de la convergence du moment d'ordre 2 vers 1

par(mfrow=c(1,2))

set.seed(1)

N_simu=1000
T0 <- hastingstraj(N_simu,c(3,3),0.01)
T1 <- hastingstraj(N_simu,c(3,3),0.1)
T2 <- hastingstraj(N_simu,c(3,3),1)
T3 <- hastingstraj(N_simu,c(3,3),10)
plot(cumsum(T0[,1]^2)/(1:(N_simu+1)),ylim = c(0,10),type = 'l',ylab = "Moment d'ordre 2 de X1",
     xlab = '',main ="Convergence du moment \n d'ordre 2 vers 1")
lines(cumsum(T1[,1]^2)/(1:(N_simu+1)),col = 'orange')
lines(cumsum(T2[,1]^2)/(1:(N_simu+1)),col = 'blue')
lines(cumsum(T3[,1]^2)/(1:(N_simu+1)),col = 'green')
abline(h=1,lty = 2)
#legend('topright',c('alpha = 0.01','alpha = 0.1','alpha = 1','alpha = 10'), lty = c(1,1,1,1),col = c('black','orange','blue',"green"))


N_simu=1000
T0 <- hastingsstep(N_simu,c(3,3),0.01)
T1 <- hastingsstep(N_simu,c(3,3),0.1)
T2 <- hastingsstep(N_simu,c(3,3),1)
T3 <- hastingsstep(N_simu,c(3,3),10)
plot(x=0:N_simu,T0[[2]],ylim = c(0,N_simu),type = 's',ylab = "Nombre d'acceptations",
     xlab = '',main ="Courbe escalier des sauts d'acceptation \n du couple candidat dans Metropolis Hasting")
lines(x=0:N_simu,T1[[2]],col = 'orange',type = 's')
lines(x=0:N_simu,T2[[2]],col = 'blue',type = 's')
lines(x=0:N_simu,T3[[2]],col = 'green',type = 's')
#legend('topright',c('alpha = 0.01','alpha = 0.1','alpha = 1','alpha = 10'), lty = c(1,1,1,1),col = c('black','orange','blue',"green"))



#+ajouter un graph au dessus qui fait un saut dès qu'on change de valeur

#c) Convergence

MC.simus <- function(N_simu,N_exp,liste_alpha){ 
  MCsimus <- matrix(NA,N_exp,length(liste_alpha))  
  for (i in 1:N_exp)
  {  
    for(j in 1: length(liste_alpha)){
      MCsimus[i,j]<-mean(hastingstraj(N_simu,c(2,2),liste_alpha[j])^2)
    }
  }
  return(MCsimus)
}

MCtest<-MC.simus(N_simu=100,N_exp=10000, liste_alpha=c(0.01,0.1,1,10))
par(mfrow=c(1,1))
boxplot(MCtest,
        names = c('alpha = 0.01','alpha = 0.1','alpha = 1','alpha = 10'),
        ylab = "Estimation du moment d'ordre 2 de X")
abline(h=1,lty = 2)

MC<-MC.simus(N_simu=100,N_exp=10000, liste_alpha=c(0.1,1,2,3,4,5,6,7,8,9,10))
par(mfrow=c(1,1))
colnames(MC)=c(0.1,1,2,3,4,5,6,7,8,9,10)
boxplot(MC)
abline(h=1,lty = 2)
MCbis<-MC.simus(N_simu=100,N_exp=10000, liste_alpha=c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4))
par(mfrow=c(1,1))
colnames(MCbis)=c(1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4)
boxplot(MCbis)
abline(h=1,lty = 2)

#L'alpha optimal est environ égal à 2.25