# ===================================================================
# Calculate generation interval distribution
# for a SEmInR model. 
#
# Method explained in:
# Champredon D, Dushoff J. 
# Intrinsic and realized generation intervals in infectious-disease transmission. 
# Proceedings of the Royal Society B: Biological Sciences 2015; 282: 20152026.
# ===================================================================

source("GI_SEmInR.R")
source("SEmInR_deterministic.R")

### SEmInR parameters:
latent_mean = 3
infectious_mean = 4
horizon = 30
nE = 5
nI = 2
n.points.GI.crv = 20

g <- calc.GI.SEmInR(latent_mean,
					infectious_mean,
					nE,nI,
					n.points.GI.crv,
					horizon)
t <- g[['times']]
g.int <- g[['intrinsic']]

plot(t, g.int,
	 typ='l', 
	 main='Intrinsic generation interval distribution')
