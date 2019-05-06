
\name{Pow.Pembedded.P}
\alias{Pow.Pembedded.P}
\title{Power}
\description{
   Compute the power for the parametric model                   
}
\usage{ 
Pow.Pembedded.P(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0)
}
\arguments{
  \item{nmax}{Sample size}
  \item{rand_ratio}{Probability of assignment to treatment arm}
  \item{effect_p}{Proportion of responders in the treatment arm}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Hazard ratio}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration in months}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 10000.}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-4.}
  \item{alpha}{Significance level. The default is 0.05.}
  \item{num_rand}{The number of replications in the re-randomization test. The default is 1000.}
  \item{nsim}{The number of simulations. The default is 1000.}
  \item{print}{0 or 1 to print information. The default is 0.}
}
\details{
  For each simulation, a simulated data set is created from the 
  \code{\link{generate_data}} function and then an estimated p-value is computed 
  by calling \code{\link{Pembedded.ReRandomizationTest.P}}.
  The power is calculated as the proportion of iterations whose estimated p-value
  was less than or equal to \code{alpha}.
}
\value{ A list containing the power and the number of simulated datasets used
        in the calculation.
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
%\examples{}
\keyword{ power, EM, survival } 
\seealso{ \code{\link{N.Pembedded.P}} }
