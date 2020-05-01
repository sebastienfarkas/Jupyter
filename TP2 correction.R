# Illustration des th??or??mes de convergence asymptotique ------------------
# Et introduction au principe de simulations de siulation

#(a)  La loi des grands nombres
#appliqu?? ?? la loi exponentielle 

n <- 10000
X <- rexp(n = n,rate = 2)
Y <- cumsum(X)/1:n
plot(1:n,Y,
     type = 'l',
     ylim = c(0,1),
     xlab='Nombre de simulations',
     ylab='Moyenne empirique',
     main="Illustration de la loi des grands nombres pour une variable exponentielle")
abline(h=0.5,col="red")

for (i in 2:10){
  X <- rexp(n,2)
  Y <- cumsum(X)/1:n
  lines(1:n,Y,col = i)
}

#(b)  Le th??or??me central limite
#appliqu?? ?? la loi exponentielle 

cltexpo <- function(n_simus,n_exp){
  empmean<-rep(0,n_exp)  
  for (i in 1:n_exp){  
    empmean[i]=mean(rexp(n_simus,2)) 
  }
  empmean
}

par(mfrow = c(2,3)) 

n_exp=100
for (i in 1:3)
{
  n_simus=10^(i+1)
  x<-cltexpo(n_simus=n_simus,n_exp = n_exp)
  hist(x,
       freq=FALSE,
       breaks = seq(0.20,0.80,length.out = 200),
       xlim=c(0.45,0.55),
       main=paste("Adéquation des moyennes empiriques \n avec la loi limite, pour ",n_simus," simulations \n et ",n_exp," expériences."),
       xlab="x",
       ylab="y")
  #Ecart type théorique
  sd=sqrt(1/n_simus*1/4)
  curve(dnorm(x,mean = 0.5,sd = sd),xlim=c(0.45,0.55),col=i+1,add=T,n=100)
}

n_exp=1000
for (i in 1:3)
{
  n_simus=10^(i+1)
  x<-cltexpo(n_simus=n_simus,n_exp = n_exp)
  hist(x,
       freq=FALSE,
       breaks = seq(0.20,0.80,length.out = 200),
       xlim=c(0.45,0.55),
       main=paste("Adéquation des moyennes empiriques \n avec la loi limite,  pour ",n_simus," simulations \n et ",n_exp," expériences."),
       xlab="x",
       ylab="y")
  #Ecart type théorique
  sd=sqrt(1/n_simus*1/4)
  curve(dnorm(x,0.5,sd),xlim=c(0.35,0.65),col=i+1,add=T,n=100)
}

par(mfrow=c(1,1))
boxplot(cltexpo(n_simus = 100,n_exp = 10),cltexpo(n_simus = 1000,n_exp = 10),cltexpo(n_simus = 10000,n_exp = 10))

#(c)  Contre exemple de le loi de Cauchy 
#pour la loi des grands nombres

par(mfrow = c(1,1)) 

n <- 100000
X <- rcauchy(n,location = 0,scale=10)
Y <- cumsum(X)/1:n
plot(1:n,Y,
     type = 'l',ylim = c(0,50),
     xlab='Nombre de simulations',
     ylab='Moyenne empirique',
     main="Divergences des moyennes empiriques de réalisations de loi de Cauchy")


for (i in 2:5){
  X <- rcauchy(n,location = 0,scale=10)
  Y <- cumsum(X)/1:n
  lines(1:n,Y,col = i)
}


cltcauchy<-function(N){
  empmean=rep(NA,100)
  for (i in 1:100){
    empmean[i]=mean(rcauchy(N,2))
  }
  empmean}
par(mfrow=c(1,2))
boxplot(main="Evolution de la variabilité des moyennes empiriques en fonction du \n nombre de simulation pour une loi exponentielle",
        cltexpo(n_simus = 100,n_exp = 10),cltexpo(n_simus = 1000,n_exp = 10),cltexpo(n_simus = 10000,n_exp = 10),cltexpo(n_simus = 100000,n_exp = 10))
boxplot(main="Evolution de la variabilité des moyennes empiriques en fonction du \n nombre de simulation pour une loi de Cauchy",
        cltcauchy(100),cltcauchy(1000),cltcauchy(10000),cltcauchy(100000),outline=T,ylim = c(-10,10))

# R??duction de la variance ------------------

#(a)  Crude MC

crudeMC=function(n,f,a=0,b=1){
  return((b-a)*mean(f(runif(n = n,min = a, max = b))))
}

#(b)  Variables antith??tiques

antitMC=function(n,f,a=0,b=1){
  x=runif(n = floor(n/2),min = a, max = b)
  return((b-a)*(mean(f(x)+f(1-x))/2))
}

#Application

inverse=function(x){
  return(1/(1+x))
}

par(mfrow=c(1,1))

Nexp=1000
n=100
ansMC=double(Nexp)
ansantitMC=double(Nexp)
for(i in 1:Nexp){
  ansMC[i]=crudeMC(n =100,f = inverse)
  ansantitMC[i]=antitMC(n =100,f = inverse)
}
abline(h=log(2),col="red")

par(mfrow=c(1,1))
boxplot(ansMC,ansantitMC, names = c('Crude MC','Antithetic'),
        ylab = 'estimation of I')

par(mfrow=c(1,2))

Nexp=100
x=seq(from=1,to=10^4,length.out = Nexp)
ans=double(Nexp)
for(i in 1:Nexp){
  ans[i]=crudeMC(n = x[i],f = inverse)
}

plot(x,ans,type="l",ylim=c(0.65,0.75),main = "Crude MC")
abline(h=log(2),col="red")

for (j in 2:9){
  Nexp=100
  x=seq(from=1,to=10^4,length.out = Nexp)
  ans=double(Nexp)
  for(i in 1:Nexp){
    ans[i]=crudeMC(n = x[i],f = inverse)
  }
  lines(x,ans,col = j)
}


Nexp=100
x=seq(from=1,to=10^4,length.out = Nexp)
ans=double(Nexp)
for(i in 1:Nexp){
  ans[i]=antitMC(n = x[i],f = inverse)
}

plot(x,ans,type="l",ylim=c(0.65,0.75),main = "Variables antithétiques MC")
abline(h=log(2),col="red")

for (j in 2:9){
  Nexp=100
  x=seq(from=1,to=10^4,length.out = Nexp)
  ans=double(Nexp)
  for(i in 1:Nexp){
    ans[i]=antitMC(n = x[i],f = inverse)
  }
  lines(x,ans,col = j)
}

#(c)  Variable de contr??le

controlMC=function(n,f,g,int_g,a=0,b=1){
  x=runif(n = n,min = a, max = b)
  return((b-a)*mean(f(x)-g(x))+int_g)
}


csensibilite <- function(n,C){
  
  var_control_inverse=function(x){
    return(C*(1+x))
  }
  
  int_control_inverse=C*3/2
  
  return(controlMC(n = n, f = inverse,g = var_control_inverse, int_g = int_control_inverse))
  
} 

Nexp=1000
n=1000
const <- seq(-1, 1, by = 0.05) 
A <- matrix(NA, length(const),1000) 
for(j in 1:length(const))
{
  for(i in 1:Nexp){
    A[j,i]=csensibilite(n, const[j])
  }
} 

rownames(A)=const
par(mfrow=c(1,1))
boxplot(t(A))

#Theorical value 
#(1-log(2)*3/2)*12 approx -0.48

var_control_inverse=function(x){
  return((1-log(2)*3/2)*12*(1+x))
}

int_control_inverse=(1-log(2)*3/2)*12*3/2



Nexp=100
n=1000
ansMC=double(Nexp)
ansantitMC=double(Nexp)
anscontrolMC=double(Nexp)
for(i in 1:Nexp){
  ansMC[i]=crudeMC(n =100,f = inverse)
  ansantitMC[i]=antitMC(n =100,f = inverse)
  anscontrolMC[i]=controlMC(n =100,f = inverse,g = var_control_inverse, int_g = int_control_inverse)
}

par(mfrow=c(1,2))
boxplot(ansMC,ansantitMC,anscontrolMC , names = c('Crude MC','Antithetic','Control'),
        ylab = 'estimation of I')


eps=0.005
proportionMC=sum(abs(ansMC-log(2))<=eps)/Nexp 
proportionantitMC=sum(abs(ansantitMC-log(2))<=eps)/Nexp
proportioncontrolMC=sum(abs(anscontrolMC-log(2))<=eps)/Nexp

plot(x=c(1,2,3),
     y=c(proportionMC,proportionantitMC,proportioncontrolMC),
     xlab="Methodes",
     ylab="Part des estimations précises à epsilon près")


#(c)  ??chantillonage pr??f??rentiel

#Les rd sont appel??s w en cours

preferentielMC=function(n,f,rd,rsimu,a,b){
  x=rsimu(n=n)
  x=x[x>a && x<b]
  return((b-a)*mean(f(x)*rd(x)))
}


h = function(x){
  return(10*exp(-2*abs(x-5)))
}


simupreference = function(n){
  return(rnorm(n,mean=5,sd=1))
}

rdunif_norm = function(x){
  return(dunif(x, 0, 10)/dnorm(x, mean=5, sd=1))
}

Nexp=1000
n=100
ansMC=double(Nexp)
ansimportanceMC=double(Nexp)
for(i in 1:Nexp){
  ansMC[i]=crudeMC(n = n,f = h,a = 0,b = 10)
  ansimportanceMC[i]=preferentielMC(a = 0,b = 10,n = n,f = h,rd = rdunif_norm,rsimu = simupreference)
}

par(mfrow=c(1,1))
boxplot(ansMC,ansimportanceMC , names = c('Crude MC','Importance'),
        ylab = 'estimation of I')
abline(h=10,col="red")

sigmasensibilite <- function(n,C){
  
  simupreference_sigma = function(n){
    return(rnorm(n,mean=5,sd=C))
  }
  
  rdunif_norm_sigma = function(x){
    return(dunif(x, 0, 10)/dnorm(x, mean=5, sd=C))
  }
  
  return(preferentielMC(n = n,f = h,rd = rdunif_norm_sigma,rsimu = simupreference_sigma,a=0,b=10))
  
} 

Nexp=1000
n=1000
const <- seq(0.05, 5, by=0.05) 
A <- matrix(NA, length(const),Nexp) 
for(j in 1:length(const))
{
  for(i in 1:Nexp){
    A[j,i]=sigmasensibilite(n, const[j])
  }
}

rownames(A)=const

par(mfrow=c(1,2))
plot(const, A[,1],type='p',
     xlab="Constante C",
     ylab="Estimation",
     main="Evolution de l'estimation \n en fonction de C",
     ylim=c(9.5,10.5))

par(mfrow=c(1,1))
boxplot(t(A),ylim=c(9,11),
        xlab="Constante C",
        ylab="Variabilité des estimations",
        main="Evolution de la variabilité des estimations \n en fonction de C",)
abline(h=10,col="red")

#Choix optimal : approximativement 0,75