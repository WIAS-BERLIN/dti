
plmatrix <- function(x, FUN, ..., mc.cores = setCores(,reprt=FALSE)){
  dx <- dim(x)[2]
  if(mc.cores>dx) mc.cores <- dx
  n <- trunc((dx-1)/mc.cores)+1
  lx <- list(NULL)
  for(i in 1:(mc.cores-1)) lx[[i]] <- x[,(i-1)*n+1:n]
  lx[[mc.cores]] <- x[,((mc.cores-1)*n+1):dx]
  cl <- makeCluster(mc <- mc.cores)
  lz <- parLapply(cl, lx, FUN , ...)
  stopCluster(cl)
  z <- matrix(0,length(lz[[1]])/n, dx)
  for(i in 1:(mc.cores-1)) z[,(i-1)*n+1:n] <- lz[[i]]
  z[,((mc.cores-1)*n+1):dx] <- lz[[mc.cores]]
  z
}

ptensnl <- function(x,ngrad,btb,sdcoef,maxit=1000,reltol=1e-7){
  nvox <- dim(x)[2]
  matrix(.C(C_dtens,
            as.integer(nvox),
            param=as.double(x[1:7,]),
            as.double(x[-(1:7),]),
            as.integer(ngrad),
            as.double(btb),
            as.double(sdcoef),
            as.double(rep(0,ngrad)),#si
            as.double(rep(1,ngrad)),#var
            as.integer(maxit),#maxit
            as.double(reltol))$param,7,nvox)
}

pmixtens <- function(x,ngrad0,maxcomp,maxit,pen,grad,reltol,th,penIC,vert){
  nvox <- length(x)/(ngrad0+3+maxcomp)
  dim(x) <- c(ngrad0+3+maxcomp,nvox)
  z <- .C(C_mixture,
          as.integer(nvox),
          as.integer(x[(ngrad0+2):(ngrad0+3+maxcomp),]),
          as.integer(ngrad0),
          as.integer(maxcomp),
          as.integer(maxit),
          as.double(pen),
          as.double(t(grad)),
          as.double(reltol),
          as.double(th),
          as.double(penIC),
          as.double(x[ngrad0+1,]),
          as.double(vert),
          as.double(x[1:ngrad0,]),
          sigma2  = double(nvox),# error variance
          orient  = double(2*maxcomp*nvox), # phi/theta for all mixture tensors
          order   = integer(nvox),   # selected order of mixture
          lev     = double(2*nvox),         # logarithmic eigenvalues
          mix     = double(maxcomp*nvox))[c("sigma2","orient","order","lev","mix")]
  rbind(z$order,z$sigma2,matrix(z$lev,2,nvox),matrix(z$mix,maxcomp,nvox),matrix(z$orient,2*maxcomp,nvox))
}
pdkiQP <- function(x,TA,Dmat,Amat){
  ##
  ##  dkiTensor CLLS-QP parallel version
  ##
  #require(quadprog)
  nvox <- dim(x)[2]
  param <- matrix(0,22,nvox)
  for(i in 1:nvox){
    dvec <- -as.vector(t(TA) %*% x[,i])
    cterm <- sum(x[,i]^2)
    resqp <- quadprog::solve.QP(Dmat, dvec, Amat,factorized=FALSE)
    resQPsolution <- resqp$solution
    param[22, i] <- cterm+2*resqp$value
    param[1:6, i] <- resQPsolution[1:6]
    param[7:21, i] <- resQPsolution[7:21] / mean(resQPsolution[1:3])^2
  }
  param
}
pgetsii30 <- function(x,maxcomp,dgrad,th,isample0,nsi,nth,nvico,nguess){
  nvox <- length(x)/(nsi+2)
  dim(x) <- c(nsi+2,nvox)
  z <- .Fortran(C_pgtsii30,
                as.double(x[1:nsi,]),
                as.double(x[nsi+1,]),#sigma2
                as.integer(nsi),
                as.integer(nvox),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(nvico),
                as.double(th),
                as.integer(nth),
                as.integer(x[nsi+2,]),#indth
                double(nsi*nvico),
                as.integer(isample0),
                as.integer(nguess),
                double(nsi),
                double(nsi*(maxcomp+2)),
                siind=integer((maxcomp+2)*nvox),
                krit=double(nvox),
                as.integer(maxcomp+2))[c("siind","krit")]
  dim(z$siind) <- c(maxcomp+2,nvox)
  rbind(z$krit,z$siind)
}
pgetsii31 <- function(x,maxcomp,dgrad,th,isample1,nsi,nth,nvico,nguess,dgradi,maxc){
  nvox <- length(x)/(nsi+3)
  dim(x) <- c(nsi+3,nvox)
  z <- .Fortran(C_pgtsii31,
                as.double(x[1:nsi,]),
                as.double(x[nsi+1,]),#sigma2
                as.integer(nsi),
                as.integer(nvox),
                as.integer(maxcomp),
                as.double(dgrad),
                as.integer(nvico),
                as.integer(x[nsi+3,]),#iandir
                as.double(th),
                as.integer(nth),
                as.integer(x[nsi+2,]),#indth
                double(nsi*nvico),
                as.integer(isample1),
                as.integer(nguess),
                double(nsi),
                double(nsi*(maxcomp+2)),
                siind=integer((maxcomp+2)*nvox),
                krit=double(nvox),
                as.integer(maxcomp+2),
                as.double(dgradi),
                as.double(maxc))[c("siind","krit")]
  dim(z$siind) <- c(maxcomp+2,nvox)
  rbind(z$krit,z$siind)
}
