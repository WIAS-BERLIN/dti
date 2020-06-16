##
##   DKI 
##

##  Diffusion Kurtosis model with nonlinear regression

dkiModel <- function(param, si, A){
    gvalue <- exp(A%*%param[1:21])
    sum((si - gvalue)^2)
    }

dkigrad <- function(param, si, A){
    gvalue <- exp(A%*%param[1:21])
    2*as.vector((gvalue-si)*gvalue)%*%A
    }

##
##  Diffusion Kurtosis model with
##  Gauss-approximation for noncentral chi
##   si/sigma is assumed to follow a noncentral chi_{2L} distribution
##   sigma should be of length ng here

dkiModelQL <- function(param, si, sigma, A, L, CL){
    ng <- dim(A)[1]
    gvalue <- exp(A%*%param)/sigma
    gvalue <- pmin(1e5,pmax(0,gvalue))
    muL <- CL * hg1f1(rep(-.5, ng), rep(L, ng), -gvalue*gvalue/2)
    sum((si/sigma - muL)^2)
    }

dkigradQL <- function(param, si, sigma, A, L, CL){
    ng <- dim(A)[1]
    gvalue <- exp(A%*%param)/sigma
    gvalue <- pmin(1e5,pmax(0,gvalue))
    mgvsq <- -gvalue*gvalue/2
    muL <- CL * hg1f1(rep(-.5, ng), rep(L, ng), mgvsq)
    CL/2/L*as.vector((muL-si/sigma)*hg1f1(rep(.5, ng), rep(L+1, ng), mgvsq)*gvalue/sigma*exp(A%*%param))%*%A
}
