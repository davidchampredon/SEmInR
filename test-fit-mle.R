
### Example code using the functions 
### defined in SEmInR: 
### - Make synthetic data
### - fit SEmInR
### - forecast future incidence with SEmInR

source("./SEmInR_deterministic.R")
source("./SEmInR_deterministic_fit_mle.R")

set.seed(1234)
t1 <- as.numeric(Sys.time())

# Simple simulation:
prm.to.fit <- c(
	infectious_mean = 3,
	popSize         = 1E5,
	R0              = 2.0)

prm.fxd <-  c(
	latent_mean     = 2,
	horizon = 200,
	nE = 5,
	nI = 5,
	init_I1 = 1,
	n.time.steps = 1000,
	per.capita = FALSE
)

# Synthetic data
sim <- simul.SEmInR.det(prm.to.fit, prm.fxd)
df  <- sim$ts
true.inc <- df$inc
tt <- df$time
t.obs <- seq(5,40,by=5)
round.inc <- ceiling(true.inc[tt %in% t.obs])
inc.obs <- rpois(n=length(t.obs), lambda = round.inc)

mydebug <- TRUE

if(mydebug){
	inc.obs <- c(2,5,8,10,12,17,28,25,39,45) ; plot(inc.obs,log='y')
	t.obs <- 1:length(inc.obs)
}

# Test likelihood function:
llk <- llk.pois(prm.to.fit,
				prm.fxd, 
				t.obs,
				inc.obs = inc.obs)
print(llk)

# Fit SEmInR model to incidence data:
boundaries <- matrix(data=c(0.5,5,
							0.6,4,
							0.88,7),
					 ncol = 2, byrow = TRUE)

init.prm.fit <-  c(
	infectious_mean = 3,
	popSize         = mean(inc.obs[1:3])/0.01,
	R0              = 2.0)

# Does the minimization is on the log of parameters?
# (can be more stable, but only for POSITIVE parameters)
logparam <- TRUE
if(logparam) init.prm.fit <- log(init.prm.fit)


FIT <- fit.mle.SEmInR(prm.to.fit = init.prm.fit, 
					  prm.fxd, 
					  t.obs, 
					  inc.obs = inc.obs,
					  logparam = logparam,
					  # lower=boundaries[,1], 
					  # upper=boundaries[,2], 
					  method = "Nelder-Mead",#"SANN",# "CG",# "Nelder-Mead",#'SANN',#"L-BFGS-B",
					  maxit = 80)
prm.fitted <- FIT[['prm.fitted']]
llkmin     <- FIT[['llkmin']]    

# Try to approximate CI value:
cival <- CI.llk.sample(CIlevel = 0.95, 
					   nsample = 100, 
					   prm.fitted, 
					   llkmin, prm.fxd, t.obs,
					   inc.obs = inc.obs,
					   prop.search = 0.5)

M <- matrix(unlist(cival), ncol=length(prm.fitted),byrow = T)
if("R0" %in% names(prm.fitted)){
	j <- which(names(prm.fitted)=="R0")
	R0.lo <- min(M[,j])
	R0.hi <- max(M[,j])
}

# Simulate forward with fitted data (point estimate):
sim.fit <- simul.SEmInR.det(prm.fitted, prm.fxd)
df.fit  <- sim.fit$ts

# Simulate forward with fitted data (CI envelop):
inc.CI <- list()
for(s in 1:length(cival)){
	if(s%%10==0) print(paste('simulating CI incidence',s,'/',length(cival)))
	sim.CI <- simul.SEmInR.det(cival[[s]],prm.fxd)
	inc.CI[[s]] <- sim.CI$ts$inc
}

m.CI <- matrix(unlist(inc.CI),ncol=length(inc.CI[[1]]),byrow = TRUE)
inc.ci.lo <- apply(m.CI,MARGIN = 2,FUN=min)
inc.ci.hi <- apply(m.CI,MARGIN = 2,FUN=max)

# Likelihood surfaces (for each parameter to fit pairs):
if(FALSE){  # <-- take some time and not used, just for example.
	lsurf <- llk.surf(grid.size = 20, prm.to.fit ,boundaries,
					  prm.fxd ,t.obs ,inc.obs, prm.fitted = prm.fitted )
}


#  - - - Plots - - - #

inc.best <- df.fit$inc
idxplot <- inc.best>0
par(mfrow = c(1,1))
plot(df.fit$time[idxplot],
	 inc.best[idxplot],
	 typ = 'l',
	 lwd = 2,
	 ylim = 1+range(inc.best,inc.ci.lo,inc.ci.hi),
	 main = 'Maximum Likelihood fit of determinisitic SEmInR', 
	 xlab = 'time',
	 ylab = 'incidence',
	 las = 1,
	 log = 'y')
lines(tt,true.inc, 
	  lwd = 1,
	  lty = 2)
for(s in 1:length(cival)){
	lines(df.fit$time,inc.CI[[s]], col=rgb(0,0,0,0.1))
}
lines(df.fit$time, inc.ci.lo,col='blue')
lines(df.fit$time, inc.ci.hi,col='blue')
points(t.obs,inc.obs,pch=16,col='orange',cex=2)
abline(v=max(t.obs), lwd=2, lty=2, col='orange')
grid()

print("  Fitted parameters:")
print(prm.fitted)


t2 <- as.numeric(Sys.time())
message(paste("completed in",round((t2-t1)/60,1),'minutes'))
