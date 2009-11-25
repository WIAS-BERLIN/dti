setClass("dwi",
         representation(.Data = "list",
                        call = "list",
                        gradient = "matrix",
                        btb    = "matrix",
                        ngrad  = "integer", # = dim(btb)[2] = dim(gradient)[2]
                        s0ind  = "integer", # indices of s0 images
                        replind = "integer", # replications in gradient design
                        ddim   = "integer",
                        ddim0  = "integer",
                        xind   = "integer",
                        yind   = "integer",
                        zind   = "integer",
                        voxelext = "numeric",
                        level  = "numeric",
                        orientation = "integer",
                        source = "character"),
         )

dwi <- function(object,  ...) cat("This object has class",class(object),"\n")
setGeneric("dwi", function(object,  ...) 
standardGeneric("dwi"))

setClass("dtiData",
         representation(si   = "array",
                        sdcoef = "numeric"),
         contains=c("list","dwi"),
         validity=function(object){
          if (any(dim(object@si)!=c(object@ddim,object@ngrad))) {
            cat("incorrect dimension of image data",dim(object@si),"\n")
            return(invisible(FALSE))
          }
          if (length(object@s0ind)<1) {
            cat("no S_0 images, parameters not identifiable \n")
             return(invisible(FALSE))
         }
          if (length(object@orientation)!=3) {
            cat("invalid orientation \n")
             return(invisible(FALSE))
         }
          if (any(sort((object@orientation)%/%2) != 0:2)) {
            cat("invalid orientation \n")
             return(invisible(FALSE))
         }
          if (object@level < 0) {
            cat("invalid level \n")
             return(invisible(FALSE))
         }
          if (length(object@sdcoef) != 4) {
             cat("invalid model for error standard deviation \n")
             return(invisible(FALSE))
         }
          if (any(object@sdcoef<0)||object@sdcoef[3]>object@sdcoef[4]) {
             cat("illegal interval of linearity in model for  error standard deviation \n")
             return(invisible(FALSE))
         }
         }
         )
setClass("dtiTensor",
         representation(method = "character",
                        D      = "array",
                        th0    = "array",
                        sigma  = "array",
                        scorr  = "array",
                        bw     = "numeric",
                        mask   = "array",
                        hmax   = "numeric",
                        outlier = "numeric",
                        scale  = "numeric"),
         contains=c("list","dwi"),
         validity=function(object){
          if (any(dim(object@D)!=c(6,object@ddim))) {
            cat("invalid dimension of tensor array D \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@th0)!=object@ddim)) {
            cat("invalid dimension of array th0\n")
            return(invisible(FALSE))
          }
          if (object@method=="linear"&any(dim(object@sigma)!=object@ddim)) {
            cat("invalid dimension of array sigma\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mask)!=object@ddim)) {
            cat("dimension of mask:",dim(object@mask),"\n")
            cat("should be:",object@ddim,"\n")
            cat("invalid dimension of array mask\n")
            return(invisible(FALSE))
          }
          if (!is.logical(object@mask)) {
            cat("invalid type of array mask, should be logical\n")
            return(invisible(FALSE))
          }
          if (length(dim(object@scorr))!=3) {
            cat("invalid dimension of scorr\n")
            return(invisible(FALSE))
          }
          if (length(object@bw)!=3) {
            cat("invalid length of bw\n")
            return(invisible(FALSE))
          }
          if (!(object@method %in% c("linear","nonlinear","unknown"))) {
            cat("method should specify either linear or nonlinear or unknown\n")
            return(invisible(FALSE))
          }
         }
         )

setClass("dtiIndices",
         representation(method = "character",
                        fa     = "array",
                        ga     = "array",
                        md     = "array",
                        andir  = "array",
                        bary   = "array"),
         contains=c("list","dwi"),
          validity=function(object){
          if (any(dim(object@fa)!=object@ddim)) {
            cat("invalid dimension of array fa\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@ga)!=object@ddim)) {
            cat("invalid dimension of array ga\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@md)!=object@ddim)) {
            cat("invalid dimension of array ra\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@andir)!=c(3,object@ddim))) {
            cat("invalid dimension of array andir\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@bary)!=c(3,object@ddim))) {
            cat("invalid dimension of array bary\n")
            return(invisible(FALSE))
          }
         }
        )

setClass("dwiQball",
         representation(what = "character",
                        order  = "integer",
                        lambda = "numeric",
                        sphcoef = "array",
                        varsphcoef = "array",
                        th0 = "array",
                        sigma  = "array",
                        scorr  = "array",
                        bw     = "numeric",
                        mask   = "array",
                        hmax   = "numeric",
                        outlier = "numeric",
                        scale  = "numeric"),
         contains=c("list","dwi"),
         validity=function(object){
          if (object@order%%2!=0) {
            cat("invalid order of spherical harmonics \n")
            return(invisible(FALSE))
          }
          if(object@what%in%c("ODF","wODF","aODF")&any(dim(object@sphcoef)!=c((object@order+1)*(object@order+2)/2,object@ddim))) {
            cat("invalid dimension of ceofficient array \n")
            return(invisible(FALSE))
          }
          if (object@lambda<0) {
            cat("invalid regularization parameter \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mask)!=object@ddim)) {
            cat("dimension of mask:",dim(object@mask),"\n")
            cat("should be:",object@ddim,"\n")
            cat("invalid dimension of array mask\n")
            return(invisible(FALSE))
          }
          if (!is.logical(object@mask)) {
            cat("invalid type of array mask, should be logical\n")
            return(invisible(FALSE))
          }
          if (length(dim(object@scorr))!=3) {
            cat("invalid dimension of scorr\n")
            return(invisible(FALSE))
          }
          if (length(object@bw)!=3) {
            cat("invalid length of bw\n")
            return(invisible(FALSE))
          }
          if (!(object@what %in% c("ODF","wODF","aODF","ADC"))) {
            cat("what should specify ODF, wODF, aODF or ADC\n")
            return(invisible(FALSE))
          }
         }
         )
setClass("dwiFiber",
          representation(call = "list",
                         fibers = "matrix",
                         startind = "integer",
                         roix  = "integer",
                         roiy = "integer",
                         roiz = "integer",
                         method = "character",
                         minanindex = "numeric",
                         maxangle   = "numeric"),
         contains=c("list","dwi"),
         validity=function(object){
            if(dim(object@fibers)[2]!=6) {
            cat("invalid dimension of fibers matrix \n")
            return(invisible(FALSE))
            }
         }
        )
setClass("dwiMixtensor",
         representation(method = "character",
                        ev     = "array",#length 2 (eigenvalues)
                        mix    = "array",
                        orient = "array",
                        order  = "array",
                        p      = "numeric", # p in "method"=="Jian"
                        th0    = "array",
                        sigma  = "array",
                        scorr  = "array",
                        bw     = "numeric",
                        mask   = "array",
                        hmax   = "numeric",
                        outlier = "numeric",
                        scale  = "numeric"),
         contains=c("list","dwi"),
         validity=function(object){
          if (any(dim(object@ev)!=c(2,object@ddim))) {
            cat("invalid dimension of eigenvalue array ev \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mix)[-1]!=object@ddim)) {
            cat("invalid dimension of array of mixture weights \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@orient)[-2]!=c(2,object@ddim))) {
            cat("invalid dimension of orientations array orient \n")
            return(invisible(FALSE))
          }
          if (dim(object@orient)[2]!=dim(object@mix)[1]) {
            cat("dimension of orientations array orient incompatible with
                 number of mixtures \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@order)!=object@ddim)) {
            cat("  \n")
            return(invisible(FALSE))
          }
          if (any(object@mix<0)) {
            cat("negative mixture coefficients \n")
            return(invisible(FALSE))
          }
          if (any(dim(object@th0)!=object@ddim)) {
            cat("invalid dimension of array th0\n")
            return(invisible(FALSE))
          }
          if (any(dim(object@mask)!=object@ddim)) {
            cat("dimension of mask:",dim(object@mask),"\n")
            cat("should be:",object@ddim,"\n")
            cat("invalid dimension of array mask\n")
            return(invisible(FALSE))
          }
          if (!is.logical(object@mask)) {
            cat("invalid type of array mask, should be logical\n")
            return(invisible(FALSE))
          }
          if (length(dim(object@scorr))!=3) {
            cat("invalid dimension of scorr\n")
            return(invisible(FALSE))
          }
          if (length(object@bw)!=3) {
            cat("invalid length of bw\n")
            return(invisible(FALSE))
          }
          if (!(object@method %in% c("mixtensor","Jian"))) {
            cat("method should specify either mixtensor or Jian \n")
            return(invisible(FALSE))
          }
         }
         )
