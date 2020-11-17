# SS2_HW1

```{r}
theta0 = 3

n = 10
n = 100
#n = 500

Xsamples = rnorm( n, theta0, 1 ) # draw n samples 

thetaBar = mean(Xsamples) # calculate the MLE 

MTilde = prod( dnorm(Xsamples,thetaBar,1) ) # the estimate of \tilde{M} 

Nsim = 2500 

# This is \pi(\theta|x) (faster)
#   1) for each value of theta evaluate the pdf of the n Gaussians
#   2) multiply by the Cauchy density 
fFunc = function(theta,Xsamples){
  
  XS = outer( Xsamples, rep(1,length(theta)) )
  TS = outer( rep(1,length(Xsamples)), theta )

  E = exp( -0.5*(XS - TS)^2 )/sqrt(2*pi)
  apply(E, 2, prod)/( pi*(1+theta^2) )

}

x=NULL
while( length(x)<Nsim ){
  # generate samples from the candidate density
  y = rcauchy( Nsim, location=0, scale=1. ) # these are proposed samples from theta ... 
  x = c(x, y[ runif(Nsim)*MTilde < fFunc(y,Xsamples)/dcauchy(y,location=0.,scale=1.) ])
}
x = x[1:Nsim]

hist(x,freq=F)

print(sprintf("n=%4d; mean(theta samples)= %10.6f; var(theta samples)= %10.6f", n, mean(x), var(x)))

```

HW2

```{r}
y1=sample(c(-1,1),10^3,rep=T)*rexp(10^3)
w1=exp(-sqrt(abs(y1)))*sin(y1)^2*(y1>0)/.5*dexp(y1)
y2=rcauchy(10^3)*2
w2=exp(-sqrt(abs(y2)))*sin(y2)^2*(y2>0)/dcauchy(y2/2)
y3=rnorm(10^3)
w3=exp(-sqrt(abs(y3)))*sin(y3)^2*(y3>0)/dnorm(y3)
boxplot(as.data.frame(cbind(w1,w2,w3)))

4*10^6*var(y1*w1/sum(w1))/mean(y1*w1/sum(w1))^2
4*10^6*var(y2*w2/sum(w2))/mean(y2*w2/sum(w2))^2
4*10^6*var(y3*w3/sum(w3))/mean(y3*w3/sum(w3))^2
```

```{r}
h=function(y){ exp(-y)*sqrt(y)/gamma(3/2)}
y = rexp(10^5,1) + 12.5
I=exp(-12.5)*sqrt(y)/gamma(3/2)
estint=cumsum(I)/(1:10^5)
esterr=sqrt(cumsum((I-estint)^2))/(1:10^5)
plot(estint,xlab="Iterations",ty="l",lwd=2,
ylim=mean(I)+20*c(-esterr[10^5],esterr[10^5]),ylab="")
integrate(h,10,Inf)
pchisq(25,3,low=F)

```



```{r}
x=c(125,18,20,34)
n=sum(x)
start=EM=cur=diff=.1
while (diff>.001)
  {
EM=c(EM,((cur*x[1]/(2+cur))+x[4])/((cur*x[1]/(2+cur))+x[2]+x[3]+x[4]))
diff=abs(cur-EM[length(EM)])
cur=EM[length(EM)]
}
M=500
MCEM=matrix(start,ncol=length(EM),nrow=500)
for (i in 2:length(EM)){
MCEM[,i]=1/(1+(x[2]+x[3])/(x[4]+rbinom(500,M*x[1],
prob=1/(1+2/MCEM[,i-1]))/M))
}
plot(EM,type="l",xlab="iterations",ylab="MCEM sequences")
upp=apply(MCEM,2,max);dow=apply(MCEM,2,min)
polygon(c(1:length(EM),length(EM):1),c(upp,rev(dow)),col="red")
lines(EM,col="blue",lty=2,lwd=2)
```

```{r}
mixnormal=function(x,theta0)
{x = runif(25)
part1=(1-theta0[5])*dnorm(x,theta0[1],theta0[3])
part2=(1-theta0[5])*dnorm(x,theta0[2],theta0[4])
gam=part2/(part1+part2)
denom1=sum(1-gam)
denom2=sum(gam)
mu1=sum((1-gam)*x)/denom1
sig1=sqrt(sum((1-gam)*((x-mu1)^2))/denom1)
mu2=sum(gam*((x-2)/2))/denom2
sig2=sqrt((sum(gam)*((((x-2)/2))-mu2)^2)/denom2)
p=mean(gam)
mixnormal=c(mu1,mu2,sig1,sig2,p)
mixnormal
}
Nsim=25
theta0=c(0,0,1^2,1^2,0.3)
mixnormal(x,theta0)
```

HW3

```{r}
# Part 1.3
Nsim=10^5
G=rt(Nsim, df=4)
for (i in 2:Nsim){
F=rnorm(1,0,1) 
rho=dt(F,1)*dnorm(G[i-1])/(dt(G[i-1],1)*dnorm(F))
G[i]=G[i-1] + (F-G[i-1])*(runif(1)<rho)
}
mean(G)
plot(cumsum(G<3)/(1:Nsim),lwd=2,ty="l",ylim=c(0.80,1), xlab = "iterations", ylab = "Pr(X < 3)", main = " Cumulative Converage Plot")
```


```{r}
mcmc=Braking1(nsim = 501)
data=cars
speed=cars[,1]
BrDis=cars[,2]

N=length(speed)
speed2=speed^2

MSR=10
SSR=100

beta=coefficients(lm(BrDis~speed+speed2))

loglike=function(a,b,c,sig)
{-(N/2)*log(sig)
-sum((BrDis-a-b*speed-c*speed2)^2)/(2*sig)
  
dcandidate=function(a,b,c,sig)
{dnorm(a,beta[1],sd=15,log=TRUE)
+dnorm(b,beta[2],sd=2,log=TRUE)
+dnorm(c,beta[3],sd=.06,log=TRUE)
-(N/2)*log(sig)-SSR/(2*sig)}

beta1=array(beta[1],dim=c(nsim,1));
beta2=array(beta[2],dim=c(nsim,1))
beta3=array(beta[3],dim=c(nsim,1));
sigma=array(MSR,dim=c(nsim,1));
for (i in 2:nsim)  
{bcandidate=c(rnorm(1,mean=beta[1],sd=15),rnorm(1,mean=beta[2],sd=2),rnorm(1,mean=beta[3],sd=.06),
1/rgamma(1,N/2,rate=SSR/2))

test=min(exp(loglike(bcandidate[1],bcandidate[2],bcandidate[3],bcandidate[4])
-loglike(beta1[i-1],beta2[i-1],beta3[i-1],sigma[i-1])
+dcandidate(beta1[i-1],beta2[i-1],beta3[i-1],sigma[i-1])
-dcandidate(bcandidate[1],bcandidate[2],bcandidate[3],bcandidate[4])),1);

rho<-(runif(1)<test)
beta1[i]=bcandidate[1]*rho+beta1[i-1]*(1-rho);
beta2[i]=bcandidate[2]*rho+beta2[i-1]*(1-rho);
beta3[i]=bcandidate[3]*rho+beta3[i-1]*(1-rho);
sigma[i]=bcandidate[4]*rho+sigma[i-1]*(1-rho);}

plot(speed,beta1[nsim]+beta2[nsim]*speed+beta3[nsim]*speed2,type="l",col="red",xlab="",ylab="",ylim=c(0,120),lwd=2)

for (i in (nsim=1):nsim){
lines(speed,beta1[i]+beta2[i]*speed+beta3[i]*speed2,col="grey",lwd=2)}
lines(speed,beta[1]+beta[2]*speed+beta[3]*speed2,col="purple",lwd=2)
points(speed,BrDis,pch=19)}

par(mfrow=c(1,3),mar=c(4,4,2,1))
plot(mcmc$a,type="l",xlab="",ylab="a");
plot(mcmc$b,type="l",xlab="",ylab="b");
plot(mcmc$c,type="l",xlab="",ylab="c");

acf(mcmc$a)
acf(mcmc$b)
acf(mcmc$c)

hist(mcmc$a,prob=T,main="",yla="",xla="a",col="blue")
hist(mcmc$b,prob=T,main="",yla="",xla="b",col="blue")
hist(mcmc$c,prob=T,main="",yla="",xla="c",col="blue")

quantile(mcmc$a,c(.025,.975))
quantile(mcmc$b,c(.025,.975))
quantile(mcmc$c,c(.025,.975))
```


```{r}
#Milad mn852
# Question 6.13
#part 1
library(mcsm)
data(challenger)
Fail=challenger[,1]
Tp=challenger[,2]
summary(glm(Fail~Tp, family = binomial))
ch=summary(glm(Fail~Tp, family = binomial))
B=as.vector(ch$coef[,1])
ch$cov.unscaled
plot(Tp,100*Fail,pch=10,col="blue",xlab="Temperature [F]",ylab="Failure[%]",main = "Failure Probability of O-ring over temperature")
curve(100/(1+exp(-B[1]-B[2]*x)),add=TRUE,col="red",lwd=2)
# Part 2
n=10^5
x=Tp
y=Fail
sig1=5
sig2=5/sd(x)
LP=function(a,b){
sum(y*(a+b*x)-log(1+exp(a+b*x)))+
dnorm(a,sd=sig1,log=TRUE)+dnorm(b,sd=sig2,log=TRUE)
}
a=b=rep(0,n)
a[1]=B[1]
b[1]=B[2]
Sc1=sqrt(ch$cov.un[1,1])
Sc2=sqrt(ch$cov.un[2,2])
for (t in 2:n)
  {
Pr1=a[t-1]+sample(c(-1,1),1)*rexp(1)*Sc1

if (log(runif(1))<LP(Pr1,b[t-1])-LP(a[t-1],b[t-1])) 
  a[t]=Pr1
else 
  a[t]=a[t-1]
Pr2=b[t-1]+sample(c(-1,1),1)*rexp(1)*Sc2

if (log(runif(1))<LP(a[t],Pr2)-LP(a[t],b[t-1])) 
  b[t]=Pr2
else 
  b[t]=b[t-1]
}
#acceptance rate
length(unique(a))/n
length(unique(b))/n
# Part 3
hist(a,prob=TRUE,col="blue",xlab=expression(al),main="Graph alpha")
hist(b,prob=TRUE,col="blue",xlab=expression(B),main="Graph beta")

acf(b,ylab=expression(B),main="alpha lag")
acf(a,ylab=expression(al),main="beta lag")

plot(b,type="l",xlab="Sample",ylab=expression(B),main="beta sample")
plot(a,type="l",xlab="Sample",ylab=expression(al),main="alpha sample")
plot(a,b,type="l",xlab=expression(al),ylab=expression(B),main="alpha & beta")
plot(Tp,100*Fail,pch=10,col="red",xlab="Temperature [F]",ylab="Failure[%]",main="Temp/Fail")

for (t in seq(100,n,le=100)) curve(100/(1+exp(-a[t]-b[t]*x)),
add=TRUE,col="green",lwd=2)
curve(100/(1+exp(-mean(a)-mean(b)*x)),add=TRUE,col="black",lwd=2.5)
postal=rep(0,1000);i=1


for (t in seq(100,n,le=1000)){ postal[i]=LP(a[t],b[t]);i=i+1}
plot(seq(100,n,le=1000),postal,type="l",xlab="Sample",ylab="Log-posterior",main="posterior")
abline(h=LP(a[1],b[1]),col="sienna",lty=2,main="posterior")
# part 4
# Failure at 50F
M1=mean(1/(1+exp(-a-b*50)))
sqrt(var(1/(1+exp(-a-b*50)),na.rm=TRUE)/length(na.omit(x)))
# Failure at 60F
mean(1/(1+exp(-a-b*60)))
sqrt(var(1/(1+exp(-a-b*60)),na.rm=TRUE)/length(na.omit(x)))
# Failure at 70F
mean(1/(1+exp(-a-b*70)))
sqrt(var(1/(1+exp(-a-b*70)),na.rm=TRUE)/length(na.omit(x)))
```


HW4


```{r}
T=500 ;p=5 ;r=0.25
X=cur=rnorm(p)
for (t in 1 :T){
for (j in 1 :p){
m=sum(cur[-j])/(p-1)
cur[j]=rnorm(1,(p-1)*r*m/(1+(p-2)*r),sqrt((1+(p-2)*r-(p-1)*r^2)/(1+(p-2)*r)))
}
X=cbind(X,cur)
}
par(mfrow=c(1,5))
for (i in 1:p){
hist(X[i,],prob=TRUE,col="wheat2",xlab="",main="")
curve(dnorm(x),add=TRUE,col="sienna",lwd=2)}


```
```{r}
T=500 ;p=5 ;r=0.25
#X=cur=rnorm(p)
#for (t in 1 :T){
#for (j in 1 :p){
#m=sum(cur[-j])/(p-1)
#cur[j]=rnorm(1,(p-1)*r*m/(1+(p-2)*r),sqrt((1+(p-2)*r-(p-1)*r^2)/(1+(p-2)*r)))
#}
#X=cbind(X,cur)
#}
#par(mfrow=c(1,5))
#for (i in 1:p){
#hist(X[i,],prob=TRUE,col="wheat2",xlab="",main="")
#curve(dnorm(x),add=TRUE,col="sienna",lwd=2)}

J=matrix(1,ncol=5,nrow=5)
I=diag(c(1,1,1,1,1))
s=(1-r)*I+r*J
rmnorm(500,s)
par(mfrow=c(1,5))
for (i in 1:p){
hist(X[i,],prob=TRUE,col="wheat2",xlab="",main="")
curve(dnorm(x),add=TRUE,col="sienna",lwd=2)}
```
```{r}
T=500 ;p=5 ;r=0.25
X=cur=rnorm(p)
for (t in 1 :T){
for (j in 1 :p){
m=sum(cur[-j])/(p-1)
cur[j]=rnorm(1,(p-1)*r*m/(1+(p-2)*r),sqrt((1+(p-2)*r-(p-1)*r^2)/(1+(p-2)*r)))
}
X=cbind(X,cur)
}




for (j in 1:m){
mea=sum(cur[-j])/(p-1)
prop=rnorm(1,(p-1)*r*mea/(1+(p-2)*r),
sqrt((1+(p-2)*r-(p-1)*r^2)/(1+(p-2)*r)))
if (sum(cur[(1:m)[-j]]^2+prop^2)<sum(cur[(m+1):p]^2))
cur[j]=prop
}
for (j in (m+1):p){
mea=sum(cur[-j])/(p-1)
prop=rnorm(1,(p-1)*r*mea/(1+(p-2)*r),
sqrt((1+(p-2)*r-(p-1)*r^2)/(1+(p-2)*r)))
if (sum(cur[(1:m)]^2)<sum(cur[((m+1):p)[-j]]^2+prop^2))
cur[j]=prop
}


par(mfrow=c(1,5))
for (i in 1:p){
hist(X[i,],prob=TRUE,col="wheat2",xlab="",main="")
curve(dnorm(x),add=TRUE,col="sienna",lwd=2)}
```

```{r}
xdata=c(3.64,2.78,2.91,2.85,2.54,2.62,3.16,2.21,4.05,2.19,
2.97,4.32,3.56,3.39,3.59,4.13,4.21,1.68,3.88,4.33)
m=length(xdata)
n=30;a=3.5 #1/3 missing data
nsim=10^4
xbar=mean(xdata)
that=array(xbar,dim=c(nsim,1))
zbar=array(a,dim=c(nsim,1))
for (i in 2:nsim){
temp=runif(n-m,min=pnorm(a,mean=that[i-1],sd=1),max=1)
zbar[i]=mean(qnorm(temp,mean=that[i-1],sd=1))
that[i]=rnorm(1,mean=(m*xbar+(n-m)*zbar[i])/n,
sd=sqrt(1/n))
}
par(mfrow=c(1,2),mar=c(5,5,2,1))
hist(that[500:nsim],col="grey",breaks=25,
xlab=expression(theta),main="",freq=FALSE)
curve(dnorm(x,mean(that),sd=sd(that)),add=T,lwd=2)
hist(zbar[500:nsim],col="grey",breaks=25,
main="",xlab= expression(bar(Z)),freq=FALSE)
curve(dnorm(x,mean(zbar),sd=sd(zbar)),add=T,lwd=2)
```



