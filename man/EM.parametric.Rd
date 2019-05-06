
\name{Pembedded.EM.P}
\alias{Pembedded.EM.P}
\title{EM algorithm}
\description{
   Parametric EM algorithm                          
}
\usage{ Pembedded.EM.P(data, time.var="X", trt.var="trt", status.var="event_status",
           effect_p=0.6, t1=1, probResponder=NULL, 
           stopTol=1e-5, maxiter=10000, print=0)
}
\arguments{
  \item{data}{Data frame or matrix containing a time-to-event variable (\code{time.var}),
              a treatment variable (\code{trt.var}), 
              and a censoring variable (\code{status.var}).}
  \item{time.var}{Observational time variable name in \code{data} (months). The default is "X".}
  \item{trt.var}{Binary treatment assignment indicator name in \code{data} coded as 
                 0 for controls and 1 for treated subjects.}
  \item{status.var}{Name of the binary censoring variable in \code{data} coded
                   as 0 for censored subjects and 1 for subjects that
                   experienced an event.}
  \item{effect_p}{Proportion of responders among the treated subjects. The default is 0.6.}
  \item{t1}{Delayed duration. The default is 1 (month).}
  \item{probResponder}{NULL or vector of initial probabilities of a treated subject being a responder.
                       The default is NULL so that the initial probability is 0.5 for
                       treated subjects.}
  \item{stopTol}{Stopping tolerance. The default is 1e-5.}
  \item{maxiter}{Maximum number of iterations. The default is 10000.}
  \item{print}{0-2 to print information. Larger values will print more information.
               The default is 0.}
}
\details{
 
}
\value{ A list containing the objects:
  \tabular{ll}{
    Name          \tab Description \cr
    converged     \tab TRUE if EM algorithm converged \cr
    lambda        \tab estimated hazard rate \cr
    baseline      \tab estimated baseline hazard \cr
    probResponder \tab estimated probability of a treated subject being a responder
  }
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="Immunotherapy.Design")
  ret <- Pembedded.EM.P(data)
  ret$lambda
}
\keyword{ EM, survival } 
\seealso{ \code{\link{Pembedded.EM.NP}} }
