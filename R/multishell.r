sphtrarea <- function(g1,g2,g3){
##  Compute area of sherical triangle spanned by vectors 
##  g1,g2,g3 on unit sphere
##  use absolute values to identify opposite directions with each other
c12 <- abs(g1%*%g2)
c13 <- abs(g1%*%g3)
c23 <- abs(g2%*%g3)
s12 <- sqrt(1-c12^2)
s13 <- sqrt(1-c13^2)
s23 <- sqrt(1-c23^2)
b1 <- (c23-c12*c13)/s12/s13 
b2 <- (c13-c12*c23)/s12/s23 
b3 <- (c12-c23*c13)/s23/s13 
acos(b1)+acos(b2)+acos(b3)-pi
}


getsphwghts <- function(g,g1,g2,g3){
##
##   compute weights for linear interpolation in g using g1,g2,g3
##
w <- numeric(3)
w[1] <- sphtrarea(g,g2,g3)
w[2] <- sphtrarea(g,g1,g3)
w[3] <- sphtrarea(g,g1,g2)
w/sum(w) 
}

getnext3g <- function(grad,bv){
##
##  calculate next neighbors on the sphere and weights for interpolation
##
grad <- grad[,bv>0]
bv <- bv[bv>0]
ubv <- unique(bv[bv>max(bv)/50])
nbv <- length(ubv)
ng <- dim(grad)[2]
ind <- array(0,c(nbv,ng,3))
w <- array(0,c(nbv,ng,3))
for(i in 1:nbv){
   indb <- (1:ng)[bv==ubv[i]]
   ind[i,indb,1] <- indb
   ind[i,indb,2] <- indb
   ind[i,indb,3] <- indb
   w[i,indb,1] <- 1
   w[i,indb,2] <- 0
   w[i,indb,3] <- 0
   for(j in (1:nbv)[-i])
      indbk <- (1:ng)[bv==ubv[j]]
      for(k in indb){
         d <- abs(t(grad[,k])%*%grad[,indbk])
         ijk <- indbk[order(d,decreasing = TRUE)[1:3]]
         ind[j,k,] <- ijk
         if(max(d)>1-1e-8){
            w[j,k,] <- c(1,0,0)
         } else {
            w[j,k,] <- getsphwghts(grad[,k],grad[,ijk[1]],grad[,ijk[2]],grad[,ijk[3]])
         }
   }
}   
#   spheres identified by bvalues in ubv
#   ind[j,k,] contains indices of gradients to be used in spherical interpolation 
#             on shell i for gradient k
#   w[j,k,]   contains the corresponding weights
list(ind=ind, w=w, ubv = ubv, nbv = nbv, bv=bv)
}

interpolatesphere <- function(theta,n3g){
##  interpolate estimated thetas to get values on all spheres
##  n3g  generated by function  getnext3g
##  dim(theta) = c(n1,n2,n3,ngrad)
dtheta <- dim(theta)
mstheta <- array(0,c(n3g$nbv,dtheta))
dim(theta) <- c(prod(dtheta[1:3]),dtheta[4])
for(i in 1:n3g$nbv){
   for(j in 1:dim(theta)[2]){
      mstheta[i,,,,j] <- theta[,n3g$ind[i,j,]]%*%n3g$w[i,j,] 
   }
}
mstheta
}


lkfullse3msh <- function(h,kappa,gradstats,vext,n){
    nbv <- gradstats$nbv
    dist <- gradstats$dist
    h <- vr <- matrix(0,ngrad,kstar)
    bvind <- gradstats$bvind
    ind <- matrix(0,5,n)
    nn <- 0
    for(i in 1:nbv){
      gshell <- list(k456=gradstats$k456[[i]],dist=dist) 
      z <- lkfullse3(h[bvind],kappa[bvind],gshell,vext,n)
      ind[1:3,nn+1:z$nind] <- z$ind[1:3,]
      ind[4:5,nn+1:z$nind] <- bvind[z$ind[4:5,]] 
#
#  convert indeces on selected shell to total indeces
#
      w[nn+1:z$nind] <- z$w
      nn <- nn + z$nind
    }  
list(h=h,kappa=kappa,ind=ind,w=w,nind=nn)
}

gethseqfullse3msh <-
function (kstar, gradstats, kappa=NULL, vext = c(1, 1)) 
{
#
#  generate information on local bandwidths and variance reduction
#  for smoothing on multiple shells
#
    nbv <- gradstats$nbv
    dist <- gradstats$dist
    h <- vr <- matrix(0,ngrad,kstar)
    n <- 0
    for(i in 1:nbv){
       gshell <- list(k456=gradstats$k456[[i]],dist=dist)
       z <- gethseqfullse3(kstar, gshell, kappa=kappa, vext=vext)
       h[gradstats$bvind[[i]],] <- z$h
       vr[gradstats$bvind[[i]],] <- z$vred
       n <- n+z$n
    }
        cat("\n total number of positive weights:",n,"mean maximal bandwidth",signif(mean(h[,kstar]),3), "\n")
    list(h=h,kappa=kappa,vred=vr,n=n)
}


getkappasmsh <- function(grad, msstructure, trace = 0, dist = 1){
ngrad <- dim(grad)[2]
nbv <- msstructure$nbv
bv <- msstructure$bv
ubv <- msstructure$ubv
bvind <- k456 <- bghat <- list(NULL)
for(i in 1:nbv){
#
#   collect information for spherical distances on each schell separately
#
ind <- (1:ngrad)[bv==ubv[i]]
z <- getkappas(grad[,ind], trace = trace, dist = dist)
bvind[[i]] <- ind
k456[[i]] <- z$k456
bghat[[i]] <- z$bghat
}
list(k456 = k456, bghat = bghat, bvind = bvind, dist=dist, nbv = nbv, ngrad = ngrad)
}

