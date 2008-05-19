cat("K. Tabelow, J. Polzehl, V. Spokoiny, and H.U. Voss,\n Diffusion Tensor Imaging: Structural Adaptive Smoothing,\n Neuroimage, 39(4), 1763--1773 (2008)\n used linear tensor estimation for their examples.\n The package also contains non-linear tensor estimation which is the new default.")
a <- readline("Use non-linear (Y, default) or linear (N) tensor estimates ?")
if (toupper(a) == "N") {
  method <- "linear"
  cat("Note: Contrary to the paper above, due to numeric issues in this demo,\n there will be some additional non-positive definite tensors in the phantoms!\n")
lambda <- 47 
# this value was used in the Neuroimage paper
} else {
  method <- "nonlinear"
lambda <- 32
}
cat("---> using",method,"tensor estimation!\n")

a <- readline("Mask small non-diffusion weighted values ? (y/n)?")
if (toupper(a) == "N") {
  mins0value <- 0
} else {
  mins0value <- 100
}




# define some constants
sigma <- 1600
rho <- 1
ddim <- c(64,64,26)
ngrad <- 25
factor <- 2.5

# read the gradient data, these are 25 gradient directions + one non-diffusion weighted
bvec <- read.table(system.file("dat/b-directions.txt",package="dti"))
bvec <- t(bvec)
btb <- matrix(0,6,ngrad+1)
btb[1,] <- bvec[1,]*bvec[1,]
btb[4,] <- bvec[2,]*bvec[2,]
btb[6,] <- bvec[3,]*bvec[3,]
btb[2,] <- 2*bvec[1,]*bvec[2,]
btb[3,] <- 2*bvec[1,]*bvec[3,]
btb[5,] <- 2*bvec[2,]*bvec[3,]

# a useful function to create a tensor of specified anisotropy
eta <- function(ai){
  aindex <- function(tensor){
    values <- eigen(matrix(tensor[c(1,2,3,2,4,5,3,5,6)],3,3))$values
    sqrt(3/2*sum((values-mean(values))^2)/sum(values^2))
  }
  risk <- function(par,ai) (ai-aindex((1-par)*c(1,0,0,0,0,0)+par*c(1,0,0,1,0,1)))^2
  optimize(f = risk, interval = c(0,1),ai=ai)$minimum
}

#
# create Phantom data
#
etas <- numeric(1001)
for(i in 1:1001) etas[i] <- eta((i-1)/1000)

ind <- array(0,ddim)
dtiso <- array(c(1,0,0,1,0,1),dim=c(6,ddim))
phi <- (1:1000)*2*pi/1000 
rad1 <- 7
rad2 <- 8
for(rad in seq(rad1,rad2,.25)){
  x <- as.integer(rad*sin(phi)+32.5)
  y <- as.integer(rad*cos(phi)+32.5)
  for(i in 1:125) ind[x[i],y[i],2:25] <- 0
  for(i in 126:250) ind[x[i],y[i],2:25] <- .6
  for(i in 251:375) ind[x[i],y[i],2:25] <- .2
  for(i in 376:500) ind[x[i],y[i],2:25] <- .8
  for(i in 501:625) ind[x[i],y[i],2:25] <- .4
  for(i in 626:750) ind[x[i],y[i],2:25] <- .7
  for(i in 751:875) ind[x[i],y[i],2:25] <- .3
  for(i in 876:1000) ind[x[i],y[i],2:25] <- .9
  for(i in 1:1000){
    etai <- etas[ind[x[i],y[i],12]*1000+1]
    dtiso[,x[i],y[i],] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
  }
}

rad1 <- 13
rad2 <- 14
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:4] <- .6
  for(i in 1:1000) ind[x[i],y[i],5:7] <- .3
  for(i in 1:1000) ind[x[i],y[i],8:10] <- .9
  for(i in 1:1000) ind[x[i],y[i],11:13] <- 0
  for(i in 1:1000) ind[x[i],y[i],14:16] <- .5
  for(i in 1:1000) ind[x[i],y[i],17:19] <- .2
  for(i in 1:1000) ind[x[i],y[i],20:22] <- .8
  for(i in 1:1000) ind[x[i],y[i],23:25] <- .4 
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}
rad1 <- 19
rad2 <- 20
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:3] <- .6
  for(i in 1:1000) ind[x[i],y[i],4:5] <- .2
  for(i in 1:1000) ind[x[i],y[i],6:7] <- .8
  for(i in 1:1000) ind[x[i],y[i],8:9] <- 0
  for(i in 1:1000) ind[x[i],y[i],10:11] <- .5
  for(i in 1:1000) ind[x[i],y[i],12:13] <- .9
  for(i in 1:1000) ind[x[i],y[i],14:15] <- .2
  for(i in 1:1000) ind[x[i],y[i],16:17] <- .6 
  for(i in 1:1000) ind[x[i],y[i],18:19] <- 0 
  for(i in 1:1000) ind[x[i],y[i],20:21] <- .9 
  for(i in 1:1000) ind[x[i],y[i],22:23] <- .3 
  for(i in 1:1000) ind[x[i],y[i],24:25] <- .7
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}

#
#   to slow
#
rad1 <- 25
rad2 <- 26
for(rad in seq(rad1,rad2,.25)){
   x <- as.integer(rad*sin(phi)+32.5)
   y <- as.integer(rad*cos(phi)+32.5)
   sphi <- sin(phi)
   cphi <- cos(phi)
   for(i in 1:1000) ind[x[i],y[i],3:5] <- (1+sin(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],6:8] <- (1+cos(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],9:11] <- 0
   for(i in 1:1000) ind[x[i],y[i],12:13] <- (1+sin(phi[i]))/2
   for(i in 1:1000) ind[x[i],y[i],14:15] <- (1+cos(phi[i]))/2
   for(i in 1:1000) ind[x[i],y[i],16:18] <- 0
   for(i in 1:1000) ind[x[i],y[i],19:21] <- (2+sin(phi[i]))/3
   for(i in 1:1000) ind[x[i],y[i],22:24] <- (2+cos(phi[i]))/3
   for(j in 1:26) {
      for(i in 1:1000){
         etai <- etas[ind[x[i],y[i],j]*1000+1]
         dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
      }
     cat(".")
   }
}

for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) {
    dtiso[,i,j,]<-0
  }
}

#  yields a maximum eigenvalue of 2.5 within a cylinder of radius 26, zero tensor outside this cylinder
dtiso <- factor*dtiso

# show images of slice 14 for all tensor components
par(mfrow=c(2,3),mar=c(1,1,1,.1))
for(i in 1:6) image(dtiso[i,,,14],col=grey((0:255)/255))

# now we want to view the projection of the zylinders onto a plane
project.cylinder <- function(obj,radius,phi=(1:1000)*2*pi/1000){
  ddim <- dim(obj)
  img <- matrix(0,length(phi),ddim[3])
  x <- as.integer(radius*sin(phi)+32.5)
  y <- as.integer(radius*cos(phi)+32.5)
  for ( i in 1:length(phi)) img[i,] <- obj[x[i],y[i],]
  img
}

# that's it!
par(mfrow=c(2,3),mar=c(1,1,1,1))
image(ind[,,22],col=grey((0:255)/255))
image(project.cylinder(ind,7.5),col=grey((0:255)/255))
image(project.cylinder(ind,13),col=grey((0:255)/255))
image(project.cylinder(ind,19),col=grey((0:255)/255))
image(project.cylinder(ind,25),col=grey((0:255)/255))

# reset S0 image
s0offa <- read.table(system.file("dat/S0ofFA.txt",package="dti"))
s0 <- s0offa[as.integer(as.vector(ind)*500+1),2]
dim(s0) <- dim(ind)
for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) s0[i,j,]<-0
}

# create noisy data
createdata.dti <- function(file,dtensor,btb,s0,sigma,level=250){
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  dim(si)<-c(ddim,ngrad)
  for (j in 1:ngrad) {
    for (i in 1:ddim[3]) {
      si[,,i,j] <- abs(fft(fft(si[,,i,j])+complex(real=rnorm(ddim[1]*ddim[2],0,sigma),imaginary=rnorm(ddim[1]*ddim[2],0,sigma)),inverse=TRUE))/ddim[1]/ddim[2]
    }
  }
  con <- file(file,"wb")
  writeBin(as.integer(si),con,2)
  close(con)
}


#   create phantom - object
createdata.dti("S_all",dtiso,btb,s0,0)
dt0obj <- dtiData(bvec,paste("S_all",sep=""),mins0value=mins0value,ddim)
dt0 <- dtiTensor(dt0obj, method=method)
dt0aniso <- dtiIndices(dt0)
file.remove("S_all")

# create noisy data
set.seed(1)
createdata.dti("S_noise_all",dtiso,btb,s0,sigma)

#  Bias of FA can not be avoided since the Expected S_0 S_b  and therefore the Expected Tensor 
#  differ from the "True" quantities 
#  We may be better off if we compare the FA with the FA of the Expected tensor computet from E S_0 and E S_b

# Read noisy data 
dtobj <- dtiData(bvec,paste("S_noise_all",sep=""),mins0value=mins0value,ddim)
dthat1 <- dtiTensor(dtobj, method=method)
dthat1aniso <- dtiIndices(dthat1)
file.remove("S_noise_all")

# adaptive smoothing
dthat4 <- dti.smooth(dtobj,hmax=4,graph=TRUE,lambda=lambda,minanindex=0,slice=15,rho=rho,lseq=NULL,method=method)
dthat4aniso <- dtiIndices(dthat4)

# plot the color-coded directional maps, phantom, noisy, smoothed
par(mfrow=c(1,3))
plot(dt0aniso,slice=15)
plot(dthat1aniso,slice=15)
plot(dthat4aniso,slice=15)

# illustrate what print() does
print(dtobj)
print(dthat1)
print(dthat4)
print(dthat4aniso)

# illustrate what summary() does
summary(dtobj)
summary(dthat1)
summary(dthat4)
summary(dthat4aniso)

# write tensor to a NIFTY-file
tensor2medinria(dthat4, "dti_art")
# read tensor from  NIFTY-file
dthat4b <- medinria2tensor("dti_art")
# plot the resulting object
plot(dthat4b,slice=15)

z <- readline("End of demo, remove created objects (Y/N) :")

graphics.off()
if(toupper(z)!="N"){
file.remove("dti_art.nii")
rm(a,btb,bvec,cphi,createdata.dti,ddim,dt0,dt0aniso,dt0obj,dthat1,dthat1aniso,
dthat4,dthat4aniso,dthat4b,dtiso,dtobj,eta,etai,etas,factor,i,ind,j,lambda,method,
mins0value,ngrad,phi,project.cylinder,rad,rad1,rad2,rho,s0,s0offa,sigma,sphi,x,y,z)

}


