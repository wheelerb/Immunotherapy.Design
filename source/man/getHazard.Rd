
\name{getHazard}
\alias{getHazard}
\title{Compute initial estimates for the baseline hazard}
\description{
   Calls the coxph function to compute initial estimates for the baseline hazard                        
}
\usage{ getHazard(time, treatment, event_status)
}
\arguments{
  \item{time}{Vector of times.}
  \item{treatment}{Binary vector of treatments (1=subject received treatment).}
  \item{event_status}{Binary vector of event status (1=subject experienced an event).}
}
\details{
 
}
\value{ Vector of baseline hazards ordered by the event times.
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="Immunotherapy.Design")
  lambda0 <- getHazard(data[, "X"], data[, "trt"], data[, "event_status"])
  lambda0[1:10]
}
\keyword{ EM, survival } 
\seealso{ \code{\link{Pembedded.EM.NP}} }
