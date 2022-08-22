
library('FME')
amp    <- 6
period <- 5
phase  <- 0.5
x <- runif(20)*13
y <- amp*sin(2*pi*x/period+phase) + rnorm(20, mean = 0, sd = 0.05)
plot(x, y, pch = 16)
cost <- function(par)
  sum((par[1] * sin(2*pi*x/par[2]+par[3])-y)^2)
p1 <- optim(par = c(amplitude = 1, phase = 1, period = 1), fn = cost)
p2 <- optim(par = c(amplitude = 1, phase = 1, period = 1), fn = cost,
            method = "SANN")
p3 <- pseudoOptim(p = c(amplitude = 1, phase = 1, period = 1),
                  lower = c(0, 1e-8, 0), upper = c(100, 2*pi, 100),
                  f = cost, control = c(numiter = 3000, verbose = TRUE))
curve(p1$par[1]*sin(2*pi*x/p1$par[2]+p1$par[3]), lty = 2, add = TRUE)
curve(p2$par[1]*sin(2*pi*x/p2$par[2]+p2$par[3]), lty = 3, add = TRUE)
curve(p3$par[1]*sin(2*pi*x/p3$par[2]+p3$par[3]), lty = 1, add = TRUE)
legend ("bottomright", lty = c(1, 2, 3),
        c("Price", "Mathematical", "Simulated annealing"))











library(ggplot2) #library for plotting
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm)

#load concentration data
df=read.table("~/Downloads/ABC_data.dat")
names(df)=c("time","ca","cb","cc")

# plot data
tmp=melt(df,id.vars=c("time"),variable.name="species",value.name="conc")
ggplot(data=tmp,aes(x=time,y=conc,color=species))+geom_point(size=3)

# rate function
rxnrate=function(t,c,parms){

  # rate constant passed through a list called parms
  k1=parms$k1
  k2=parms$k2

  # c is the concentration of species

  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["A"] #dcA/dt
  r[2]=k1*c["A"]-k2*c["B"] #dcB/dt
  r[3]=k2*c["B"] #dcC/dt

  # the computed derivatives are returned as a list
  # order of derivatives needs to be the same as the order of species in c
  return(list(r))

}


# predicted concentration for a given parameter set
cinit=c(A=1,B=0,C=0)
t=df$time
parms=list(k1=2,k2=1)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
head(out)




ssq=function(parms){

  # inital concentration
  cinit=c(A=1,B=0,C=0)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(0,5,0.1),df$time)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  k1=parms[1]
  k2=parms[2]
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rxnrate,parms=list(k1=k1,k2=k2))

  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% df$time,]
  # Evaluate predicted vs experimental residual
  preddf=melt(outdf,id.var="time",variable.name="species",value.name="conc")
  expdf=melt(df,id.var="time",variable.name="species",value.name="conc")
  ssqres=preddf$conc-expdf$conc

  # return predicted vs experimental residual
  return(ssqres)

}




# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
parms=c(k1=0.5,k2=0.5)
# fitting
fitval=nls.lm(par=parms,fn=ssq)



# Summary of fit
summary(fitval)

# Estimated parameter
parest=as.list(coef(fitval))
parest

# degrees of freedom: # data points - # parameters
dof=3*nrow(df)-2
dof

# mean error
ms=sqrt(deviance(fitval)/dof)
ms

# variance Covariance Matrix
S=vcov(fitval)
S



# plot of predicted vs experimental data

# simulated predicted profile at estimated parameter values
cinit=c(A=1,B=0,C=0)
t=seq(0,5,0.2)
parms=as.list(parest)
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)
outdf=data.frame(out)
names(outdf)=c("time","ca_pred","cb_pred","cc_pred")

# Overlay predicted profile with experimental data
tmppred=melt(outdf,id.var=c("time"),variable.name="species",value.name="conc")
tmpexp=melt(df,id.var=c("time"),variable.name="species",value.name="conc")
p=ggplot(data=tmppred,aes(x=time,y=conc,color=species,linetype=species))+geom_line()
p=p+geom_line(data=tmpexp,aes(x=time,y=conc,color=species,linetype=species))
p=p+geom_point(data=tmpexp,aes(x=time,y=conc,color=species))
p=p+scale_linetype_manual(values=c(0,1,0,1,0,1))
p=p+scale_color_manual(values=rep(c("red","blue","green"),each=2))+theme_bw()
print(p)



Sinv=solve(S)

# draw the confidence region
# get points for a circle with radius r
r=sqrt(qf(0.95,2,58)*2)
theta=seq(0,2*pi,length.out=100)
z=cbind(r*cos(theta),r*sin(theta))
# transform points of circle into points of ellipse using
# svd of inverse covariance matrix
Sinv_svd=svd(Sinv) # inverse of covariance matrix
xt=t(Sinv_svd$v)%*%diag(1/sqrt(Sinv_svd$d))%*%t(z) # transform from circle to ellispse
x=t(xt)
# translate the ellipse so that center is the estimated parameter value
x=x+matrix(rep(as.numeric(parest),100),nrow=100,byrow=T)

plot(x[,1],x[,2],type="l",xlab="k1",ylab="k2",lwd=2)
points(parest$k1,parest$k2,pch=20,col="blue",cex=2)



# Simulation based estimation of uncertainty

# store original experimental data in a separate dataframe
dforig=df

# conc profile based on estimated k1 and k2
cinit=c(A=1,B=0,C=0)
t=dforig$time
parms=parest
out=ode(y=cinit,times=t,func=rxnrate,parms=parms)

outsim=matrix(0,nrow=nrow(dforig),ncol=4)
outsim[,1]=out[,1]

# number of simulations
nsim=1000

parsim=matrix(0,nrow=nsim,ncol=2)
colnames(parsim)=c("k1","k2")

for (i in 1:nsim){

  # Simulate data set by adding normal random variable with mean 0 and stdev from fit

  outsim[,2:4]=out[,2:4]+matrix(rnorm(3*nrow(dforig)),nrow=nrow(dforig),ncol=3)*ms
  df=data.frame(outsim)
  names(df)=c("time","ca","cb","cc")

  # get parameter estimate for the simulated dataset
  parms=as.numeric(parest)
  fitsim=nls.lm(par=parms,fn=ssq)
  # store estimated parameters in the ith row
  parsim[i,]=coef(fitsim)


}

# plot the parameter estimates from the 1000 simulations
plot(parsim[,1],parsim[,2],xlab="k1",ylab="k2")
# overlay the 95% ellipse computed previously
lines(x[,1],x[,2],col="blue",lwd=2)






# percentage of parameters from simulation within the 95% ellipse
tmp=rep(0,length.out=nsim)
for(i in 1:nsim){
  tmp[i]=(parsim[i,]-as.numeric(parest))%*%Sinv%*%(parsim[i,]-as.numeric(parest))
}
sum(tmp <= qf(0.95,2,58)*2)/nsim














library("deSolve")

learning <- function(t, Knowledge, p) {
  with (as.list(p), {
    # knowledge increase
    dKnowledge <- rate*Knowledge*(1-Knowledge/maxKnowledge)
    list(dKnowledge)
  })
}
parms <- c (
  rate         = 0.1,    # rate at which new knowledge is acquired
  maxKnowledge = 1)      # maximal knowledge that can be acquired

# student 1
y1 <- c(Knowledge = 0.01)  # initial knowledge about the subject
times <- 1:100
out   <- ode(y = y1, times = times, parms = parms, func = learning)

# student 2: more initial knowledge
y2 <- c(Knowledge = 0.1)
out2  <- ode(y = y2, times = times, parms = parms, func = learning)

# student 3: good basic knowledge but too many other things
y3 <- c(Knowledge = 0.1)
parms2 <- parms; parms2["rate"] <- 0.05
out3  <- ode(y = y3, times = times, parms = parms2, func = learning)

plot(out, out2, out3)





library("growthrates")
library("lattice")
data(bactgrowth)
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))

## Get table from single experiment
dat <- splitted.data[["D:0:1"]]

## (1) Spline fit
fit1 <- fit_spline(dat$time, dat$value, spar=0.5)

## derive start parameters from spline fit
p <- coef(fit1)

## (2) exponential model using the first 10 data points
first10 <-  dat[1:10, ]
fit2 <- fit_growthmodel(grow_exponential, p=p, time=first10$time, y=first10$value)

## (3) Logistic model for all data
p <- c(coef(fit1), K = max(dat$value))
fit3 <- fit_growthmodel(grow_logistic, p=p, time=dat$time, y=dat$value, transform="log")

plot(fit1)
lines(fit2, col="green")
lines(fit3, col="red")












x <- 1.3
y <- "hello"
a <- x
a

x -> b
b

a
x <- a <- b

x = a

x <<- 2
x

hist(2)

plot(1:4, c(3,4,3,6), type="l", col="red")

x <- 1:20
x
y <- matrix(x, nrow=5, ncol=4)
y
as.vector


# row-wise creation of a matrix
x <- matrix(1:20, nrow=5, ncol=4, byrow=TRUE)
x

# transpose of a matrix
x <- t(x)
x


L1 <- list(a=1:10, b=c(1,2,3), x="hello")
L1

L2 <- list(a=5:7, b=L1)
str(L2)

x <- c(a=1.2, b=2.3, c=6)
x

L <- list(a=1:3, b="hello")
L

names(L)

names(L) <- c("numbers", "text")

x <- 1:5
names(x) <- letters[1:5]
x

# Select and reorder data frame columns
x <- matrix(1:16, nrow=4)
df <- as.data.frame(x)
df
names(df) <- c("N", "P", "O2", "C")
df
df2 <- df[c("C", "N", "P")]
df2


plot.default

?plot.default

?par


r  <- 0.5
N0 <- 10
dt <- 0.1
time <- seq(0, 10, dt)
# analytical solution
N <- N0 * exp(r * time)
plot(time, N, type="l")



N <- numeric(length(time))
N[1] <- N0
for (i in 2:length(time)) {
  N[i] <- N[i-1] + r * N[i-1] * dt
}
plot(time, N, type = "l")



library(deSolve)
model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <-   r * N
    list(c(dN))
  })
}
y0     <- c(N = 0.1)
parms  <- c(r = 0.1, K = 10)
times  <- seq(0, 100, 1)
out <- ode(y0, times, model, parms) # <---
plot(out)





library(deSolve)
model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <-   r * N * (1 - N / K)
    list(c(dN))
  })
}
y0     <- c(N = 0.1)
parms  <- c(r = 0.1, K = 10)
times  <- seq(0, 100, 1)
out <- ode(y0, times, model, parms) # <---
plot(out)








library("growthrates")
data(bactgrowth)
dat <- multisplit(value ~ time | strain + conc + replicate, data = bactgrowth)[["D:0:1"]]

# determine mumax for single data set from smoothing spline
fit0 <- fit_spline(dat$time, dat$value)

# all splines
fit1 <- all_splines(value ~ time | strain + conc + replicate, data = bactgrowth, spar = 0.5)

# extract results
results1 <- results(fit1)

# initial parameters
p <- c(coef(fit0), K = max(dat$value))

# avoid negative parameters
lower = c(y0 = 0, mumax = 0, K = 0)

# fit all models
fit2 <- all_growthmodels(value ~ time | strain + conc + replicate,
                         data = bactgrowth, FUN=grow_logistic, p = p, lower = lower, ncores = 2)



results2 <- results(fit2)
pdf("all_plots.pdf", width=12, height=64)
par(mfrow=c(24,3))
plot(fit2)
dev.off()







x <- seq(5, 100, 5)
y <- c(0.1, 2.2, 3.1, 1.5, 8.9, 8, 8.4, 9.8, 9.3, 10.6, 12, 13.6,
       13.1, 13.3, 11.6, 14.7, 12.6, 11.9, 10.9, 9.4)
plot(x,y)


# differential equations (ODE)
ode_K_exp <- function (time, init, parms, ...) {
  with(as.list(c(parms, init)), {
    dy <- mumax * y * (1 - y/K)
    dK <- d_K * K
    list(c(dy, dK), log_y = unname(log(y)))
  })
}

# numerical solution
grow_K_exp <- function(time, parms, ...) {
  init    <- parms[c("y0", "K")]            # initial values
  names(init) <- c("y", "K")                # force names of state variables
  odeparms <- parms[c("mumax", "d_K")]      # the parms of the ODE model
  out <- ode(init, time, ode_K_exp, parms = odeparms)
  out
}

## optional: make it a "growthmodel" for the package
grow_K_exp <- growthmodel(grow_K_exp, pnames = c("y0", "K", "mumax", "d_K"))


head(grow_K_exp(time = 1:10, c(y0 = .1, K = 1, mumax = 0.1, d_K = -0.02)))


fit <- fit_growthmodel(grow_K_exp,
                       p = c(y0 = 0.1, mumax = 0.1, K = 20, d_K = -0.01), time = x, y = y)
growthrates::summary(fit)

plot(fit)

