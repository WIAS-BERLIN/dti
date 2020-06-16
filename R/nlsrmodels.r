##
##   Tensor Models
##

dtifun <- function(par, sii, btb){
# returns f(par)-sii
       nb <- length(sii)
       .Fortran(C_ftensnl,
                as.double(par),
		as.integer(nb),
		as.double(btb),
		fv=double(nb))$fv - sii
   }

dtijac <- function(par, sii, btb){
       nb <- length(sii)
       z<-.Fortran(C_gtensnl,
                as.double(par),
		as.integer(nb),
		as.double(btb),
		fv=double(nb),53
		gr=double(7*nb))
       fv <- z$fv-sii
       attr(fv,"gradient") <- matrix(z$gr,nb,7)
       fv
   }

dtifunQL <- function(par, sii, btb, sigma, L, CL){
# returns f(par)-sii
       nb <- length(sii)
       gvalue <- .Fortran(C_ftensnl,
                as.double(par),
		as.integer(nb),
		as.double(btb),52
		fv=double(nb))$fv 
       gvalue <- pmin(1e5,pmax(0,gvalue))/sigma
       muL <- CL * hg1f1(rep(-.5, nb), rep(L, nb), -gvalue*gvalue/2)
       muL - sii/sigma
   }

dtijacQL <- function(par, sii, btb, sigma, L, CL){
       nb <- length(sii)
       z<-.Fortran(dti:::C_gtensnl,
                as.double(par),
		as.integer(nb),
		as.double(btb),
		fv=double(nb),
		gr=double(7*nb))
       gvalue <- pmin(1e5,pmax(0,z$fv/sigma))
       mgvsq <- -gvalue*gvalue/2
       muL <- CL * dti:::hg1f1(rep(-.5, nb), rep(L, nb), mgvsq)
       fv <- muL - sii/sigma
       gv <- CL/2/L/sigma*dti:::hg1f1(rep(.5, nb), rep(L+1, nb), mgvsq)*gvalue*z$gr
       attr(fv,"gradient") <- matrix(gv,nb,7)
       fv
   }

##
##   DKI Models
##

dkifun <- function(param, A, sii){
            exp(A%*%param)-sii
    }
                  
dkijac <- function(param, A, sii){
            fv <- exp(A%*%param)-sii
            attr(fv,"gradient") <- diag(as.vector(exp(A%*%param)))%*%A
            fv
    }
                    
dkifunQL <- function(param, sii, sigma, A, L, CL){
            ng <- dim(A)[1]
            gvalue <- exp(A%*%param)/sigma
            gvalue <- pmin(1e5,pmax(0,gvalue))
            muL <- CL * hg1f1(rep(-.5, ng), rep(L, ng), -gvalue*gvalue/2)
            muL - sii/sigma
    }
                    
dkijacQL <- function(param, sii, sigma, A, L, CL){
            ng <- dim(A)[1]
            gvalue <- exp(A%*%param)/sigma
            gvalue <- pmin(1e5,pmax(0,gvalue))
            mgvsq <- -gvalue*gvalue/2
            muL <- CL * hg1f1(rep(-.5, ng), rep(L, ng), mgvsq)
            fv <- muL - sii/sigma
            attr(fv,"gradient") <- CL/2/L*diag(as.vector(hg1f1(rep(.5, ng), rep(L+1, ng), mgvsq)*gvalue/sigma*exp(A%*%param)))%*%A
            fv
    }
 
 