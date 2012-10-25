getsphericalharmonicseven <- function(order,theta,phi){
#
#   compute spherical harmonics 
#  (symmetric modified SH Basis used by Descoteaux (2008))
#
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
if(length(theta)!=length(phi)) stop("need same length of theta and phi")
kseq <- seq(0,order,2)
n <- length(phi)
values <- matrix(0,(order+1)*(order+2)/2,n)
if(require(gsl,quietly = TRUE,warn.conflicts = FALSE)){
for(k in kseq){
mseq <- seq(-k,k,1)
for(m in mseq){
ind <- (k^2+k+2)/2+m
z <- legendre_sphPlm(k,abs(m),cos(theta))
if(m < 0){
z <- sqrt(2)*z*cos(m*phi)
} 
if(m > 0){
z <- sqrt(2)*z*sin(m*phi)
}
values[ind,] <- z
}
}
} else {
warning("gsl package not available \n returning zeros instead of estimates")
}
values
}
getsphericalharmonicsevenR <- function(order,theta,phi){
#
#   compute spherical harmonics
#
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
if(length(theta)!=length(phi)) stop("need same length of theta and phi")
kseq <- seq(0,order,2)
n <- length(phi)
values <- matrix(0,(order+1)*(order+2)/2,n)
if(require(gsl,quietly = TRUE,warn.conflicts = FALSE)){
for(k in kseq){
mseq <- seq(-k,k,1)
for(m in mseq){
ind <- (k^2+k+2)/2+m
z <- legendre_sphPlm(k,abs(m),cos(theta))*cos(m*phi)
if(m < 0) z <- (-1)^abs(m)*z
values[ind,] <- z
}
}
} else {
warning("gsl package not available \n returning zeros instead of estimates")
}
values
}
getsphericalharmonicsevenI <- function(order,theta,phi){
#
#   compute spherical harmonics
#
order <- as.integer(max(0,order))
if(order%%2==1){
warning("maximum order needs to be even, increase order by one")
order <- order+1
} 
if(length(theta)!=length(phi)) stop("need same length of theta and phi")
kseq <- seq(0,order,2)
n <- length(phi)
values <- matrix(0,(order+1)*(order+2)/2,n)
if(require(gsl,quietly = TRUE,warn.conflicts = FALSE)){
for(k in kseq){
mseq <- seq(-k,k,1)
for(m in mseq){
ind <- (k^2+k+2)/2+m
z <- legendre_sphPlm(k,abs(m),cos(theta))*sin(m*phi)
if(m < 0) z <- (-1)^abs(m+1)*z
values[ind,] <- z
}
}
} else {
warning("gsl package not available \n returning zeros instead of estimates")
}
values
}

getsphericalharmonicsall <- function(order,theta,phi){
#
#   compute spherical harmonics
#  (symmetric modified SH Basis used by Descoteaux (2008))
#
order <- as.integer(max(0,order))
if(length(theta)!=length(phi)) stop("need same length of theta and phi")
kseq <- 0:order
n <- length(phi)
values <- matrix(0,(order+1)^2,n)
if(require(gsl)){
l <- 1
for(k in kseq){
mseq <- (-k):k
for(m in mseq){
z <- legendre_sphPlm(k,abs(m),cos(theta))
if(m < 0){
z <- sqrt(2)*z*cos(m*phi)
} 
if(m > 0){
z <- sqrt(2)*z*sin(m*phi)
}
values[l,] <- z
l <- l+1
}
}
} else {
warning("gsl package not available \n returning zeros instead of estimates")
}
values
}

