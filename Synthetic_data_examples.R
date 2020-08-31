# Load required packages
library(mvnfast)
library(raster)
library(rasterVis)

# Load required functions
source("Supporting_files/MCMC.R")

# Global simulation settings
J <- 3
n.total <- 900
n.train <- 200
n.test <- 200
p.true <- 0.5 
beta.true <- 0
itr <- 10^5
burn.in <- 10^4

set.seed(3681)
i.sample <- sample(1:n.total,n.train+n.test)
i.train <- i.sample[1:n.train]
i.test <-  i.sample[-c(1:n.train)]

X.train <- matrix(1,n.train,1)
X.test <- matrix(1,n.test,1)
X.all <- matrix(1,n.total,1)

rl.S <- raster(,xmn=0, xmx=1, ymn=0, ymx=1, nrows=n.total^0.5,ncols=n.total^0.5)
sf.S <- rasterToPolygons(rl.S)
S.all <- coordinates(rl.S)
S.train <- S.all[i.train,] 
S.test <- S.all[i.test,] 

##################################################################################################
# Scenario 1 
##################################################################################################

# Generate synthetic data
s1.f.true <- disaggregate(0.05+raster(matrix(c(-1,-1,-1,-1,-2.5,-2.5,
                         -1,-1,-1,-1,-2.5,-2.5,
                         2,2,2,2,2,2,
                         2,2,2,2,2,2,
                         -1.4,-1.4,-1.4,0.1,0.1,0.1,
                         -1.4,-1.4,-1.4,0.1,0.1,0.1),6,6,byrow=TRUE)),fact=5)[] 

s1.psi.true <- pnorm(X.all%*%beta.true + s1.f.true)  
s1.rl.psi.true <- rl.S
s1.rl.psi.true[] <- s1.psi.true 

set.seed(4245)
s1.z.all <- rbinom(n.total,1,s1.psi.true)

set.seed(3923)
s1.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s1.z.all,n.total,J)
s1.Y.train <- s1.Y.all[i.train,]
s1.Y.test <- s1.Y.all[i.test,]

# Fit occupancy models
s1.mcmc.regtree <- MCMC(Y.train = s1.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s1.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 4502)

s1.mcmc.svr <- MCMC(Y.train = s1.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s1.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                        method = "svr",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 4928)

s1.mcmc.lrgp <- MCMC(Y.train = s1.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s1.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                    method = "lrgp",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed = 8929)

s1.mcmc.gmrf <- MCMC(Y.train = s1.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s1.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1919)

# Examine trace plots
par(mfrow=c(1,3))
plot(s1.mcmc.regtree$beta,typ="l")
plot(s1.mcmc.regtree$p,typ="l")
site <- 1
plot(s1.mcmc.regtree$psi[,site],typ="l")

plot(s1.mcmc.svr$beta,typ="l")
plot(s1.mcmc.svr$p,typ="l")
site <- 1
plot(s1.mcmc.svr$psi[,site],typ="l")

plot(s1.mcmc.lrgp$beta,typ="l")
plot(s1.mcmc.lrgp$p,typ="l")
site <- 1
plot(s1.mcmc.lrgp$psi[,site],typ="l")

plot(s1.mcmc.gmrf$beta,typ="l")
plot(s1.mcmc.gmrf$p,typ="l")
site <- 1
plot(s1.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s1.rl.E.psi.regtree <- rl.S
s1.rl.E.psi.svr <- rl.S
s1.rl.E.psi.lrgp <- rl.S
s1.rl.E.psi.gmrf <- rl.S
s1.rl.E.psi.regtree[] <- apply(s1.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s1.rl.E.psi.svr[] <- apply(s1.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s1.rl.E.psi.lrgp[] <- apply(s1.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s1.rl.E.psi.gmrf[] <- apply(s1.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure 3
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_3.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s1.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s1.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s1.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s1.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s1.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()

##################################################################################################
# Scenario 2
##################################################################################################

# Generate synthetic data
s2.f.true <- ((S.all[,1]-0.5)^2+(S.all[,2]-0.5)^2 - 0.1665)/0.12

s2.psi.true <- pnorm(X.all%*%beta.true + s2.f.true)  
s2.rl.psi.true <- rl.S
s2.rl.psi.true[] <- s2.psi.true 

set.seed(2220)
s2.z.all <- rbinom(n.total,1,s2.psi.true)

set.seed(4765)
s2.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s2.z.all,n.total,J)
s2.Y.train <- s2.Y.all[i.train,]
s2.Y.test <- s2.Y.all[i.test,]

# Fit occupancy models
s2.mcmc.regtree <- MCMC(Y.train = s2.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s2.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 1753)

s2.mcmc.svr <- MCMC(Y.train = s2.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s2.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                    method = "svr",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed = 2448)

s2.mcmc.lrgp <- MCMC(Y.train = s2.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s2.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "lrgp",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1999)

s2.mcmc.gmrf <- MCMC(Y.train = s2.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s2.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 8222)

# Examine trace plots
par(mfrow=c(1,3))
plot(s2.mcmc.regtree$beta,typ="l")
plot(s2.mcmc.regtree$p,typ="l")
site <- 1
plot(s2.mcmc.regtree$psi[,site],typ="l")

plot(s2.mcmc.svr$beta,typ="l")
plot(s2.mcmc.svr$p,typ="l")
site <- 1
plot(s2.mcmc.svr$psi[,site],typ="l")

plot(s2.mcmc.lrgp$beta,typ="l")
plot(s2.mcmc.lrgp$p,typ="l")
site <- 1
plot(s2.mcmc.lrgp$psi[,site],typ="l")

plot(s2.mcmc.gmrf$beta,typ="l")
plot(s2.mcmc.gmrf$p,typ="l")
site <- 1
plot(s2.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s2.rl.E.psi.regtree <- rl.S
s2.rl.E.psi.svr <- rl.S
s2.rl.E.psi.lrgp <- rl.S
s2.rl.E.psi.gmrf <- rl.S
s2.rl.E.psi.regtree[] <- apply(s2.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s2.rl.E.psi.svr[] <- apply(s2.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s2.rl.E.psi.lrgp[] <- apply(s2.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s2.rl.E.psi.gmrf[] <- apply(s2.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure s1
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_s1.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s2.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s2.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s2.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s2.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s2.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()


##################################################################################################
# Scenario 3
##################################################################################################

# Generate synthetic data
s3.f.true <- cos((2*pi*S.all[,1])+(2*pi*S.all[,2]))

s3.psi.true <- pnorm(X.all%*%beta.true + s3.f.true)  
s3.rl.psi.true <- rl.S
s3.rl.psi.true[] <- s3.psi.true 

set.seed(6320)
s3.z.all <- rbinom(n.total,1,s3.psi.true)

set.seed(1623)
s3.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s3.z.all,n.total,J)
s3.Y.train <- s3.Y.all[i.train,]
s3.Y.test <- s3.Y.all[i.test,]

# Fit occupancy models
s3.mcmc.regtree <- MCMC(Y.train = s3.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s3.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 4222)

s3.mcmc.svr <- MCMC(Y.train = s3.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s3.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                    method = "svr",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed = 2998)

s3.mcmc.lrgp <- MCMC(Y.train = s3.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s3.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "lrgp",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 3923)

s3.mcmc.gmrf <- MCMC(Y.train = s3.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s3.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1119)

# Examine trace plots
par(mfrow=c(1,3))
plot(s3.mcmc.regtree$beta,typ="l")
plot(s3.mcmc.regtree$p,typ="l")
site <- 1
plot(s3.mcmc.regtree$psi[,site],typ="l")

plot(s3.mcmc.svr$beta,typ="l")
plot(s3.mcmc.svr$p,typ="l")
site <- 1
plot(s3.mcmc.svr$psi[,site],typ="l")

plot(s3.mcmc.lrgp$beta,typ="l")
plot(s3.mcmc.lrgp$p,typ="l")
site <- 1
plot(s3.mcmc.lrgp$psi[,site],typ="l")

plot(s3.mcmc.gmrf$beta,typ="l")
plot(s3.mcmc.gmrf$p,typ="l")
site <- 1
plot(s3.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s3.rl.E.psi.regtree <- rl.S
s3.rl.E.psi.svr <- rl.S
s3.rl.E.psi.lrgp <- rl.S
s3.rl.E.psi.gmrf <- rl.S
s3.rl.E.psi.regtree[] <- apply(s3.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s3.rl.E.psi.svr[] <- apply(s3.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s3.rl.E.psi.lrgp[] <- apply(s3.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s3.rl.E.psi.gmrf[] <- apply(s3.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure s2
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_s2.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s3.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s3.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s3.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s3.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s3.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()

##################################################################################################
# Scenario 4
##################################################################################################

# Generate synthetic data
Q <- icar.Q(S.all,0.04,0.99,FALSE)
set.seed(3991)
s4.f.true <-rmvn(1,rep(0,n.total),solve(Q))

s4.psi.true <- pnorm(X.all%*%beta.true + s4.f.true)  
s4.rl.psi.true <- rl.S
s4.rl.psi.true[] <- s4.psi.true 

set.seed(3409)
s4.z.all <- rbinom(n.total,1,s4.psi.true)

set.seed(2223)
s4.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s4.z.all,n.total,J)
s4.Y.train <- s4.Y.all[i.train,]
s4.Y.test <- s4.Y.all[i.test,]

# Fit occupancy models
s4.mcmc.regtree <- MCMC(Y.train = s4.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s4.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 5252)

s4.mcmc.svr <- MCMC(Y.train = s4.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s4.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                    method = "svr",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed =2008)

s4.mcmc.lrgp <- MCMC(Y.train = s4.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s4.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "lrgp",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1213)

s4.mcmc.gmrf <- MCMC(Y.train = s4.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s4.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10^6,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 9197)

# Examine trace plots
par(mfrow=c(1,3))
plot(s4.mcmc.regtree$beta,typ="l")
plot(s4.mcmc.regtree$p,typ="l")
site <- 1
plot(s4.mcmc.regtree$psi[,site],typ="l")

plot(s4.mcmc.svr$beta,typ="l")
plot(s4.mcmc.svr$p,typ="l")
site <- 1
plot(s4.mcmc.svr$psi[,site],typ="l")

plot(s4.mcmc.lrgp$beta,typ="l")
plot(s4.mcmc.lrgp$p,typ="l")
site <- 1
plot(s4.mcmc.lrgp$psi[,site],typ="l")

plot(s4.mcmc.gmrf$beta,typ="l")
plot(s4.mcmc.gmrf$p,typ="l")
site <- 1
plot(s4.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s4.rl.E.psi.regtree <- rl.S
s4.rl.E.psi.svr <- rl.S
s4.rl.E.psi.lrgp <- rl.S
s4.rl.E.psi.gmrf <- rl.S
s4.rl.E.psi.regtree[] <- apply(s4.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s4.rl.E.psi.svr[] <- apply(s4.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s4.rl.E.psi.lrgp[] <- apply(s4.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s4.rl.E.psi.gmrf[] <- apply(s4.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure s3
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_s3.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s4.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s4.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s4.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s4.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s4.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()


##################################################################################################
# Scenario 5
##################################################################################################

# Generate synthetic data
D <- rdist(S.all)
C <- exp(-D/0.2)
set.seed(1312)
s5.f.true <- rmvn(1,rep(0,n.total),C)

s5.psi.true <- pnorm(X.all%*%beta.true + s5.f.true)  
s5.rl.psi.true <- rl.S
s5.rl.psi.true[] <- s5.psi.true 

set.seed(1400)
s5.z.all <- rbinom(n.total,1,s5.psi.true)

set.seed(2423)
s5.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s5.z.all,n.total,J)
s5.Y.train <- s5.Y.all[i.train,]
s5.Y.test <- s5.Y.all[i.test,]

# Fit occupancy models
s5.mcmc.regtree <- MCMC(Y.train = s5.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s5.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 5252)

s5.mcmc.svr <- MCMC(Y.train = s5.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s5.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                    method = "svr",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed =2008)

s5.mcmc.lrgp <- MCMC(Y.train = s5.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s5.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "lrgp",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1213)

s5.mcmc.gmrf <- MCMC(Y.train = s5.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s5.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 9197)

# Examine trace plots
par(mfrow=c(1,3))
plot(s5.mcmc.regtree$beta,typ="l")
plot(s5.mcmc.regtree$p,typ="l")
site <- 1
plot(s5.mcmc.regtree$psi[,site],typ="l")

plot(s5.mcmc.svr$beta,typ="l")
plot(s5.mcmc.svr$p,typ="l")
site <- 1
plot(s5.mcmc.svr$psi[,site],typ="l")

plot(s5.mcmc.lrgp$beta,typ="l")
plot(s5.mcmc.lrgp$p,typ="l")
site <- 1
plot(s5.mcmc.lrgp$psi[,site],typ="l")

plot(s5.mcmc.gmrf$beta,typ="l")
plot(s5.mcmc.gmrf$p,typ="l")
site <- 1
plot(s5.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s5.rl.E.psi.regtree <- rl.S
s5.rl.E.psi.svr <- rl.S
s5.rl.E.psi.lrgp <- rl.S
s5.rl.E.psi.gmrf <- rl.S
s5.rl.E.psi.regtree[] <- apply(s5.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s5.rl.E.psi.svr[] <- apply(s5.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s5.rl.E.psi.lrgp[] <- apply(s5.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s5.rl.E.psi.gmrf[] <- apply(s5.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure s4
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_s4.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s5.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s5.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s5.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s5.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s5.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()

##################################################################################################
# Scenario 6
##################################################################################################

# Generate synthetic data
D <- rdist(S.all)
C <- exp(-(D/0.2)^2)
set.seed(9532)
s6.f.true <- rmvn(1,rep(0,n.total),0.5*C)

s6.psi.true <- pnorm(X.all%*%beta.true + s6.f.true)  
s6.rl.psi.true <- rl.S
s6.rl.psi.true[] <- s6.psi.true 

set.seed(3409)
s6.z.all <- rbinom(n.total,1,s6.psi.true)

set.seed(2223)
s6.Y.all <- matrix(rbinom(n.total*J,1,p.true),n.total,J)*matrix(s6.z.all,n.total,J)
s6.Y.train <- s6.Y.all[i.train,]
s6.Y.test <- s6.Y.all[i.test,]

# Fit occupancy models
s6.mcmc.regtree <- MCMC(Y.train = s6.Y.train,X.train = X.train,S.train = S.train,
                        Y.test = s6.Y.test,X.test = X.test,S.test = S.test,
                        X.pred = X.all,S.pred = S.all,
                        method = "rpart",
                        sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                        beta.start = 0,p.start = 0.5,iter = itr,seed = 5252)

s6.mcmc.svr <- MCMC(Y.train = s6.Y.train,X.train = X.train,S.train = S.train,
                    Y.test = s6.Y.test,X.test = X.test,S.test = S.test,
                    X.pred = X.all,S.pred = S.all,
                    method = "svr",
                    sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                    beta.start = 0,p.start = 0.5,iter = itr,seed =2008)

s6.mcmc.lrgp <- MCMC(Y.train = s6.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s6.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "lrgp",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 1213)

s6.mcmc.gmrf <- MCMC(Y.train = s6.Y.train,X.train = X.train,S.train = S.train,
                     Y.test = s6.Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.all,S.pred = S.all,
                     method = "gmrf",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = itr,seed = 9197)

# Examine trace plots
par(mfrow=c(1,3))
plot(s6.mcmc.regtree$beta,typ="l")
plot(s6.mcmc.regtree$p,typ="l")
site <- 1
plot(s6.mcmc.regtree$psi[,site],typ="l")

plot(s6.mcmc.svr$beta,typ="l")
plot(s6.mcmc.svr$p,typ="l")
site <- 1
plot(s6.mcmc.svr$psi[,site],typ="l")

plot(s6.mcmc.lrgp$beta,typ="l")
plot(s6.mcmc.lrgp$p,typ="l")
site <- 1
plot(s6.mcmc.lrgp$psi[,site],typ="l")

plot(s6.mcmc.gmrf$beta,typ="l")
plot(s6.mcmc.gmrf$p,typ="l")
site <- 1
plot(s6.mcmc.gmrf$psi[,site],typ="l")

# Make raster layers to map E(psi|Y)
s6.rl.E.psi.regtree <- rl.S
s6.rl.E.psi.svr <- rl.S
s6.rl.E.psi.lrgp <- rl.S
s6.rl.E.psi.gmrf <- rl.S
s6.rl.E.psi.regtree[] <- apply(s6.mcmc.regtree$psi[-c(1:burn.in),],2,mean)
s6.rl.E.psi.svr[] <- apply(s6.mcmc.svr$psi[-c(1:burn.in),],2,mean)
s6.rl.E.psi.lrgp[] <- apply(s6.mcmc.lrgp$psi[-c(1:burn.in),],2,mean)
s6.rl.E.psi.gmrf[] <- apply(s6.mcmc.gmrf$psi[-c(1:burn.in),],2,mean)

# Construct figure s5
colramp <- rev(terrain.colors(30))
m1 <- c(0.2, 0.2) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- (8.5-(m2[1]+m2[3]+m2[3]+m2[3]+m2[3]+m2[3]+0.2+m2[2]))/4 # plotting panel height
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]
h<-m1[1]+d+m1[2]
pdf("fig_s5.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 13), 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6,0,
                      rep(0, 13)), 3, 13, byrow = TRUE), 
             heights = c(m1[1], d, m1[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], d, m2[3], d, m2[3],0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s6.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
plot(sf.S[i.train,],add=TRUE,border="gray47")
box()
image(s6.rl.E.psi.regtree,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s6.rl.E.psi.svr,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s6.rl.E.psi.lrgp,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
image(s6.rl.E.psi.gmrf,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.3, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0, 1, 0.25),cex.axis=1.2, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4,cex=0.9)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
dev.off()


##################################################################################################
# Synthetic data LPPD and absolute difference calculations
##################################################################################################

# Scenario 1
# Calculate -2*LPPD
-2*sum(log(apply(s1.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s1.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s1.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s1.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s1.rl.psi.true[] - s1.rl.E.psi.regtree[]))
mean(abs(s1.rl.psi.true[] - s1.rl.E.psi.svr[]))
mean(abs(s1.rl.psi.true[] - s1.rl.E.psi.lrgp[]))
mean(abs(s1.rl.psi.true[] - s1.rl.E.psi.gmrf[]))


# Scenario 2
# Calculate -2*LPPD
-2*sum(log(apply(s2.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s2.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s2.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s2.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s2.rl.psi.true[] - s2.rl.E.psi.regtree[]))
mean(abs(s2.rl.psi.true[] - s2.rl.E.psi.svr[]))
mean(abs(s2.rl.psi.true[] - s2.rl.E.psi.lrgp[]))
mean(abs(s2.rl.psi.true[] - s2.rl.E.psi.gmrf[]))

# Scenario 3
# Calculate -2*LPPD
-2*sum(log(apply(s3.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s3.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s3.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s3.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s3.rl.psi.true[] - s3.rl.E.psi.regtree[]))
mean(abs(s3.rl.psi.true[] - s3.rl.E.psi.svr[]))
mean(abs(s3.rl.psi.true[] - s3.rl.E.psi.lrgp[]))
mean(abs(s3.rl.psi.true[] - s3.rl.E.psi.gmrf[]))


# Scenario 4
# Calculate -2*LPPD
-2*sum(log(apply(s4.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s4.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s4.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s4.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s4.rl.psi.true[] - s4.rl.E.psi.regtree[]))
mean(abs(s4.rl.psi.true[] - s4.rl.E.psi.svr[]))
mean(abs(s4.rl.psi.true[] - s4.rl.E.psi.lrgp[]))
mean(abs(s4.rl.psi.true[] - s4.rl.E.psi.gmrf[]))

# Scenario 5
# Calculate -2*LPPD
-2*sum(log(apply(s5.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s5.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s5.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s5.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s5.rl.psi.true[] - s5.rl.E.psi.regtree[]))
mean(abs(s5.rl.psi.true[] - s5.rl.E.psi.svr[]))
mean(abs(s5.rl.psi.true[] - s5.rl.E.psi.lrgp[]))
mean(abs(s5.rl.psi.true[] - s5.rl.E.psi.gmrf[]))


# Scenario 6
# Calculate -2*LPPD
-2*sum(log(apply(s6.mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s6.mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s6.mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(s6.mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))

# Calculate average absolute distance between psi and E(psi|Y)
mean(abs(s6.rl.psi.true[] - s6.rl.E.psi.regtree[]))
mean(abs(s6.rl.psi.true[] - s6.rl.E.psi.svr[]))
mean(abs(s6.rl.psi.true[] - s6.rl.E.psi.lrgp[]))
mean(abs(s6.rl.psi.true[] - s6.rl.E.psi.gmrf[]))


##################################################################################################
# Figure 1
##################################################################################################


# Set color ramp palette
colramp <- rev(terrain.colors(30))

# Setup plotting device
m1 <- c(0.1, 0.1) # top and bottom margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d <- 2 # plotting panel height
h<-m1[1]+d+m1[2]+d+m1[1]
w<-m2[1]+d+m2[3]+d+m2[3]+d+m2[3]+0.2+m2[2]

pdf("fig_1.pdf",width=w,height=h, paper='special')

mp <- layout(matrix(c(rep(0, 9), 0, 1, 0, 2, 0, 3, 0, 4, 0, 
                      rep(0, 9), 0, 5, 0, 6, 0, 7, 0, 4, 0,
                      rep(0,9)), 5, 9, byrow = TRUE), 
             heights = c(m1[1], d, m1[2],d, m1[1]), 
             #widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], 0.2, m2[2]), 
             widths = c(m2[1], d, m2[3], d, m2[3], d, m2[3], 0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(s1.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
image(s2.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
image(s3.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
#legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 1.5,
#      adj = c(1.25, 0))
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xaxs = "i",yaxs = "i")
axis(4, at = seq(0, 1, 0.25), cex.axis=1.5,las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4.75, cex=1.3)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0,
            0, 1, 1)
box()
image(s4.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
image(s5.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "e.)",  bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
image(s6.rl.psi.true,col = colramp,zlim = c(0, 1),asp=TRUE,axes=FALSE)
legend("topleft", legend = "f.)",  bty = "n", text.col = 1, cex = 2,
       adj = c(1.2, 0))
box()
dev.off()


