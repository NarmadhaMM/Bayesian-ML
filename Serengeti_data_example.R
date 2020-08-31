# Load required packages
library(rgdal)
library(raster)
library(ggplot2)
library(ggmap)
library(ggsn)

# Load required functions
source("MCMC.R")

# Load Thomson's gazelle occupancy data
df.occ <- read.csv("Supporting_files/Occupancy data.csv")

# Load shapefile of sampling grid
sf.S <- readOGR("Supporting_files/Sampling grid.shp")

# Preliminary steps
pt.sampled.sites <- SpatialPoints(df.occ[,1:2],proj4string=crs(sf.S)) 
pt.sampled.sites <- spTransform(pt.sampled.sites,crs("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"))
sf.S <- spTransform(sf.S,crs("+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"))

# Construct Figure 2
Api_key<-"AIzaSyDkSJA0fWMBzDbe3Hsg5FyJgjOF9MT-bEc"
register_google(key = Api_key,write = TRUE)
rtp<-spTransform(sf.S,crs("+proj=longlat +datum=WGS84"))
rtp@data$id <- 1:nrow(rtp@data)   # add id column for join
rtpFort <- fortify(rtp, data = rtp@data)
rtpFortMer <- merge(rtpFort, rtp@data, by.x = 'id', by.y = 'id')  # join data

scale_x_longitude <- function(xmin, xmax, step) {
  xbreaks <- seq(xmin,xmax,step)
  xlabels<-ifelse(xbreaks < 0, parse(text=paste0(xbreaks,"^o", "*W")), ifelse(xbreaks > 0, parse(text=paste0(xbreaks,"^o", "*E")),x))
  return(scale_x_continuous("Longitude", breaks = xbreaks,labels = xlabels,expand = c(0, 0)))
}
scale_y_latitude <- function(ymin, ymax, step) {
  ybreaks <- seq(ymin,ymax,step)
  ylabels<-ifelse(ybreaks < 0, parse(text=paste0(ybreaks,"^o", "*S")), ifelse(ybreaks > 0, parse(text=paste0(ybreaks,"^o", "*N")),x))
  scale_y_continuous("Latitude", breaks = ybreaks,labels = ylabels,expand = c(0, 0))
}

g <- get_map(location="Serengeti",zoom=9)
gg <- ggmap(g, extent="normal") +
  scale_x_longitude(xmin=34, xmax=35.5, step=.5) +
  scale_y_latitude(ymin=-2.5, ymax=-1, step=.5) +
  theme(axis.text = element_text(size=12)) +
  geom_polygon(data = rtpFortMer,
               aes(x = long, y = lat, group=group),
               alpha = 0.3,
               size = 0.2,
               color='black')+
  scalebar(x.min = 34.1, x.max = 35.5,
           y.min = -2.63, y.max = -1.5,
           location = "bottomleft", st.dist = .04,st.size = 4,
           dist = 20, dist_unit = 'km',transform = TRUE, model = 'WGS84')
pdf(file = "Fig 2.pdf",width = 6, height = 6)
gg
dev.off()

# Make training data set
set.seed(1299)
i.train <- sample(1:195,100)
Y.train <- df.occ[i.train,-c(1:2)]
X.train <- model.matrix(~1,data=df.occ[i.train,])
S.train <- coordinates(pt.sampled.sites[i.train,])
colnames(S.train) <- c("x","y")
S.train <- S.train/max(c(extent(sf.S)[2]-extent(sf.S)[1],extent(sf.S)[4]-extent(sf.S)[3]))

# Make test data set
Y.test <- df.occ[-c(i.train),-c(1:2)]
X.test <- model.matrix(~1,data=df.occ[-c(i.train),])
S.test <- coordinates(pt.sampled.sites[-c(i.train),])
colnames(S.test) <- c("x","y")
S.test <- S.test/max(c(extent(sf.S)[2]-extent(sf.S)[1],extent(sf.S)[4]-extent(sf.S)[3]))

# Make variables needed for prediction (i.e., spatial mapping)
X.pred <- matrix(1,dim(sf.S)[1],1)
S.pred <- coordinates(sf.S)
colnames(S.pred) <- c("x","y")
S.pred <- S.pred/max(c(extent(sf.S)[2]-extent(sf.S)[1],extent(sf.S)[4]-extent(sf.S)[3]))


mcmc.regtree <- MCMC(Y.train = Y.train,X.train = X.train,S.train = S.train,
                     Y.test = Y.test,X.test = X.test,S.test = S.test,
                     X.pred = X.pred,S.pred = S.pred,
                     method = "rpart",
                     sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                     beta.start = 0,p.start = 0.5,iter = 10^5, seed = 9902)

mcmc.svr <- MCMC(Y.train = Y.train,X.train = X.train,S.train = S.train,
                 Y.test = Y.test,X.test = X.test,S.test = S.test,
                 X.pred = X.pred,S.pred = S.pred,
                 method = "svr",
                 sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                 beta.start = 0,p.start = 0.5,iter = 10^5,seed = 4999)

mcmc.lrgp <- MCMC(Y.train = Y.train,X.train = X.train,S.train = S.train,
                  Y.test = Y.test,X.test = X.test,S.test = S.test,
                  X.pred = X.pred,S.pred = S.pred,
                  method = "lrgp",
                  sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                  beta.start = 0,p.start = 0.5,iter = 10^5,seed = 4429)

mcmc.gmrf <- MCMC(Y.train = Y.train,X.train = X.train,S.train = S.train,
                  Y.test = Y.test,X.test = X.test,S.test = S.test,
                  X.pred = X.pred,S.pred = S.pred,
                  method = "gmrf",
                  sigma2.beta = 1,sigma2.f = 10,f.hat.start.var = 4,
                  beta.start = 0,p.start = 0.5,iter = 10^5,seed = 1313)

# Examine trace plots
par(mfrow=c(1,3))
plot(mcmc.regtree$beta,typ="l")
plot(mcmc.regtree$p,typ="l")
site <- 1
plot(mcmc.regtree$psi[,site],typ="l")

plot(mcmc.svr$beta,typ="l")
plot(mcmc.svr$p,typ="l")
site <- 1
plot(mcmc.svr$psi[,site],typ="l")

plot(mcmc.lrgp$beta,typ="l")
plot(mcmc.lrgp$p,typ="l")
site <- 1
plot(mcmc.lrgp$psi[,site],typ="l")

plot(mcmc.gmrf$beta,typ="l")
plot(mcmc.gmrf$p,typ="l")
site <- 1
plot(mcmc.gmrf$psi[,site],typ="l")


# Set burn-in
burn.in <- 10^4

# Calculate -2*LPPD
-2*sum(log(apply(mcmc.regtree$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(mcmc.svr$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(mcmc.lrgp$ppd[-c(1:burn.in),],2,mean)))
-2*sum(log(apply(mcmc.gmrf$ppd[-c(1:burn.in),],2,mean)))


# Construct figure 4
rl.E.psi.regtree <- rasterFromXYZ(cbind(coordinates(sf.S),colMeans(mcmc.regtree$psi[-c(1:burn.in),])),digits=0,crs="+init=epsg:32736",res=2236)
rl.E.psi.svr <- rasterFromXYZ(cbind(coordinates(sf.S),colMeans(mcmc.svr$psi[-c(1:burn.in),])),digits=0,crs="+init=epsg:32736",res=2236)
rl.E.psi.lrgp <- rasterFromXYZ(cbind(coordinates(sf.S),colMeans(mcmc.lrgp$psi[-c(1:burn.in),])),digits=0,crs="+init=epsg:32736",res=2236)
rl.E.psi.gmrf <- rasterFromXYZ(cbind(coordinates(sf.S),colMeans(mcmc.gmrf$psi[-c(1:burn.in),])),digits=0,crs="+init=epsg:32736",res=2236)


colramp <- rev(terrain.colors(30))
m1 <- c(0.1, 0.1) # top, bottom and vertical margins
m2 <- c(0.1, 0.8, 0.1) # left, right, and horizontal margins
d1 <- (8.5-(m2[1]+m2[3]+m2[3]+0.2+m2[2]))/2 # plotting panel width
d2 <- (d1/21*13)+0.3 # plotting panel height
h<-m1[1]+d2+m1[2]+d2+m1[1]
w<-m2[1]+d1+m2[3]+d1+m2[3]+0.2+m2[2]
pdf("fig_4.pdf",width=w,height=h, paper='special')
mp <- layout(matrix(c(rep(0, 7), 0, 1, 0, 2, 0, 3,0,
                      rep(0, 7), 0, 4, 0, 5, 0, 3,0,
                      rep(0,7)), 5, 7, byrow = TRUE), 
             heights = c(m1[1], d2, m1[2],d2, m1[1]), 
             widths = c(m2[1], d1, m2[3], d1, m2[3], 0.2, m2[2]), 
             respect = TRUE)

par(mar = c(0,0,0,0))
image(rl.E.psi.regtree,col = colramp,zlim = c(0.4, 1),asp=TRUE,axes=FALSE)
plot(sf.S,add=TRUE,border="black",lwd=.2)
legend("topleft", legend = "a.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.55, 0))
image(rl.E.psi.svr,col = colramp,zlim = c(0.4, 1),asp=TRUE,axes=FALSE)
plot(sf.S,add=TRUE,border="black",lwd=.2)
legend("topleft", legend = "b.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.55, 0))
plot(c(0, 1), c(0.4, 1), type = "n", axes = FALSE, xaxs = "i",
     yaxs = "i")
axis(4, at = seq(0.4, 1, 0.2),cex.axis=1.5, las = 1)
mtext(expression(paste("Probability of occupancy (",psi[italic(i)],")")),4, line = 4.75,cex=1.3)
rasterImage(as.raster(matrix(rev(colramp)), ncol = 1), 0, 0.40, 1, 1)
box()
image(rl.E.psi.lrgp,col = colramp,zlim = c(0.4, 1),asp=TRUE,axes=FALSE)
plot(sf.S,add=TRUE,border="black",lwd=.2)
legend("topleft", legend = "c.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.55, 0))
image(rl.E.psi.gmrf,col = colramp,zlim = c(0.4, 1),asp=TRUE,axes=FALSE)
plot(sf.S,add=TRUE,border="black",lwd=.2)
legend("topleft", legend = "d.)", bty = "n", text.col = 1, cex = 2,
       adj = c(1.55, 0))
dev.off()
