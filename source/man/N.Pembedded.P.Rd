
\name{N.Pembedded.P}
\alias{N.Pembedded.P}
\title{Sample Size}
\description{
   Compute the sample size for a given power for the parametric model                   
}
\usage{ 
N.Pembedded.P(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
              lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
              alpha=0.05, num_rand=1000, nsim=1000, min.N=100, max.N=700,
              tol.power=0.01, tol.N=1, print=1)
}
\arguments{
  \item{power}{The desired power. The default is 0.8.}
  \item{rand_ratio}{Allocation ratio}
  \item{effect_p}{Proportion of responders in the treatment arm}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Hazard ratio}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 1000.}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-4.}
  \item{alpha}{Significance level. The default is 0.05.}
  \item{num_rand}{Number of replications in the re-randomization test. The default is 1000.}
  \item{nsim}{Number of simulations in computing power. The default is 1000.}
  \item{min.N}{Lower bound for the sample size. The default is 100.}
  \item{max.N}{Upper bound for the sample size. The default is 700.}
  \item{tol.power}{Stopping tolerance for the power. The default is 0.01.}
  \item{tol.N}{Stopping tolerance for the sample size. The default is 1.}
  \item{print}{0 or 1 to print information. The default is 1.}
}
\details{
  This uses a bisection method to estimate the sample size. At each iteration,
  the estimated power \code{power_est} is computed using \code{\link{Pow.Pembedded.P}} 
  for a given sample size holding all other parameters fixed.
  The algorithm terminates when \code{abs(power - power_est) <= tol.power} or
  when the length of the estimated interval containing the sample size is
  less than or equal to \code{tol.N}.
}
\value{ A list containing the sample size and power.
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
%\examples{}
\keyword{ power, EM, survival } 
\seealso{ \code{\link{Pow.Pembedded.P}} }
