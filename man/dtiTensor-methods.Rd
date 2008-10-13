\name{dtiTensor-methods}
\docType{methods}
\alias{dtiTensor-methods}
\alias{dtiTensor,ANY-method}
\alias{dtiTensor,dtiData-method}
\title{ Methods for Function `dtiTensor' in Package `dti' }
\description{ The method estimates, in each voxel, the diffusion tensor from
the DWI data contained in an object of class \code{dtiData}.
}
\usage{
  \S4method{dtiTensor}{dtiData}(object, method="nonlinear", varmethod="replicates", varmodel="local")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{dtiData}}
  \item{method}{ if \code{model=="linear"} estimates are obtained using a linearization of the tensor model. This was the estimate used in Tabelow et.al. (2008). \code{model=="nonlinear"} uses a nonlinear regression model with reparametrization that ensures the tensor to be  positive semidefinite, see  Koay et.al. (2006).}
  \item{varmethod}{specifies the method for estimating the error  variance. If \code{varmethod=="replicates"} the error variance is estimated from replicated gradient directions if possible. Otherwise an estimate is obtained from the residual sum of squares.}
  \item{varmodel}{if \code{varmodel=="global"} a homogeneous variance estimate is assumed and estimated as the median of the local variance estimates.}
  \item{\dots}{ one or several  of the parameters above }
}
\value{ an object of class \code{dtiTensor}.
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning. }

\item{obj = "dtiData"}{ Estimate diffusion tensor from data in each voxel with the different options for the regression type and model for variance estimation. }

}
}
\references{  K. Tabelow, J. Polzehl, H.U. Voss, and V. Spokoiny. 
              \emph{Diffusion Tensor Imaging: Structural adaptive smoothing}, 
              NeuroImage 39(4), 1763-1773 (2008).

              C.G. Koay, J.D. Carew, A.L. Alexander, P.J. Basser and M.E. Meyerand.
              \emph{ Investigation of Anomalous Estimates of Tensor-Derived Quantities in Diffusion Tensor Imaging}, 
              Magnetic Resonance in medicine, 2006, 55, 930-936.

              \url{ http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{Karsten Tabelow tabelow@wias-berlin.de, Joerg Polzehl polzehl@wias-berlin.de}
\seealso{ \code{\link{dtiData}}, \code{\link{readDWIdata}}, \code{\link{dtiIndices-methods}},  \code{\link{medinria2tensor}}, \code{\link{tensor2medinria}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiIndices}}}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ methods }
\keyword{ manip }
\keyword{ models }