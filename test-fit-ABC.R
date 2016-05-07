source("SEmInR_deterministic_fit_ABC.R")
source("SEmInR_deterministic.R")
set.seed(1234)

# TRUE PARAMTERS:
prm.to.fit <- c(
	latent_mean     = 2,
	infectious_mean = 2,
	popSize         = 1E5,
	R0              = 2.3)

priors <- list(infectious_mean = list('unif', prm=list(1,5)),
			   latent_mean     = list('unif', prm=list(1,5)),
			   popSize         = list('lnorm', prm=list(log(1e4),1.1)),
			   R0              = list('lnorm',prm=list(log(2),0.8))
)
priors.prm.to.fit <- priors
prm.fxd <-  c(
	horizon = 200,
	nE = 2,
	nI = 2,
	init_I1 = 1,
	n.time.steps = 500,
	per.capita = FALSE
)


# Synthetic data
sim.true <- simul.SEmInR.det(prm.to.fit, prm.fxd)
df  <- sim.true$ts
true.inc <- df$inc
tt <- df$time
t.obs <- seq(5,20,by=2)

dt  <- tt[2]-tt[1]
idx <- vector()
for(k in 1:length(t.obs)) { idx[k] <- which(abs(tt-t.obs[k]) < dt/1.99)[1] }
t <- tt[idx]
round.inc <- ceiling(true.inc[idx])

inc.obs <- rpois(n=length(t.obs), lambda = round.inc)
#inc.obs <- c(4,6,6,8,9,11,12,15) ; plot(inc.obs,log='y')
# inc.obs <- c(1,3,1,8,9,15,19,23)
# plot(true.inc)


nABC <- 1E3
post.prop <- 0.01
nABC*post.prop  
CI <- 0.95

par(mfrow=c(5,5))
system.time(
	fit <- fit.ABC.SEmInR(t.obs, inc.obs, 
						  prm.fxd, priors.prm.to.fit, 
						  nABC, post.prop)
)

Mpost <- fit[['posteriors']]

par(mfrow=c(2,2))
for(j in 1:(ncol(Mpost)-1)){
	trueval <- prm.to.fit[[colnames(Mpost)[j]]]
	hist(Mpost[,j],breaks = min(30,nABC),
		 col = 'lightgrey',border = NA,
		 main = colnames(Mpost)[j], 
		 xlab='',ylab='', yaxt='n',
		 xlim = range(Mpost[,j],trueval))
	abline(v=median(Mpost[,j]),lwd=3)
	abline(v=quantile(Mpost[,j], probs = 0.5+CI/2), lwd=1)
	abline(v=quantile(Mpost[,j], probs = 0.5-CI/2), lwd=1)
	abline(v=trueval,lwd=3,lty=2,col='red')
}

# Sample from posteriors:

postinc <- post.incidence(Mpost, t=tt, CI)
inc.md <- postinc[['inc.md']]
inc.lo <- postinc[['inc.lo']]
inc.hi <- postinc[['inc.hi']]


myplot <- function(dolog,showall=FALSE){
	plot(tt, true.inc,typ='l',lty=2,log=dolog,ylim=1+range(inc.hi,true.inc))
	points(t.obs, inc.obs,pch=16)
	lines(tt,inc.md,col=rgb(1,0,0),lwd=2)
	lines(tt,inc.lo,col=rgb(1,0,0,0.5))
	lines(tt,inc.hi,col=rgb(1,0,0,0.5))
	if(showall){
		for(i in 1:ncol(inc.post)) lines(tt,inc.post[,i],col=rgb(1,0,0,0.1))
	}
}

par(mfrow=c(2,1))
myplot('y', showall=TRUE)
myplot('')