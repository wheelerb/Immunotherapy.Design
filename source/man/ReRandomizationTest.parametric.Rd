
\name{Pembedded.ReRandomizationTest.P}
\alias{Pembedded.ReRandomizationTest.P}
\title{Randomization test}
\description{
   Compute a randomization test p-value where test statisitic is calculated based on a parametric model.                     
}
\usage{ Pembedded.ReRandomizationTest.P(data, time.var="X", trt.var="trt", status.var="event_status",
           effect_p=0.6, t1=1, stopTol=1e-5, maxiter=10000, print=0, num_rand=10000)
}
\arguments{
  \item{data}{Data frame or matrix containing a time-to-event variable (\code{time.var}),
              a treatment variable (\code{trt.var}), 
              and a censoring variable (\code{status.var}).}
  \item{time.var}{Observational variable name in \code{data}.}
  \item{trt.var}{Name of treatment assignment indicator in \code{data} coded as 
                 0 for control subjects and 1 for treated subjects.}
  \item{status.var}{Name of the binary censoring variable in \code{data} coded
                   as 0 for censored subjects and 1 for subjects that experienced an event.}
  \item{effect_p}{Proportion of responders among the treated subjects. The default is 0.6.}
  \item{t1}{Delayed duration. The default is 1 (month).}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-5.}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 10000.}
  \item{print}{0-2 to print information. Larger values will print more information.
               The default is 0.}
  \item{num_rand}{The number of replications in the re-randomization test. The default is 10000.}
}
\details{
  In each randomization, the treatment label is resampled and then the EM algorithm is called.
  The final p-value is based on all randomizations in which the EM algorithm converged.
}
\value{ A list containing the objects:
  \tabular{ll}{
    Name          \tab Description \cr
    p.val.rerand  \tab re-randomization test p-value \cr
    baseline      \tab estimated baseline hazard from observed data \cr
    lambda        \tab estimated hazard rate from observed data 
  }
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="Immunotherapy.Design")
  set.seed(1)
  ret <- Pembedded.ReRandomizationTest.P(data)
  ret$p.val.rerand
}
\keyword{ EM, survival } 
\seealso{ \code{\link{Pembedded.ReRandomizationTest.NP}} }
