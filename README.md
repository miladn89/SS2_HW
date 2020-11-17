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

