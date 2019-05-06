
\name{generate_data}
\alias{generate_data}
\title{Simulated data}
\description{
   Generate simulated data                          
}
\usage{ 
generate_data(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
              lambda1=0.117, HR=0.5, tau=12*5, t1=1)
}
\arguments{
  \item{nmax}{Sample size}
  \item{rand_ratio}{Allocation ratio}
  \item{effect_p}{Proportion of responders in the treatment arm}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Hazard ratio}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration}
}
\details{
 
}
\value{ A data frame with columns:
  \tabular{ll}{
    Name          \tab Description                                \cr
    id            \tab id variable                                \cr
    trt           \tab treatment allocation: 1 = treatment arm    \cr
    Z             \tab patient's response status                  \cr
    tau           \tab total study duration                       \cr
    enroll_time   \tab patients' enrollment times                 \cr
    time_to_event \tab patients' event times                      \cr
    event_status  \tab censoring indicator                        \cr
    X             \tab observational time                         \cr
    t1            \tab delayed duration
  }
}
\references{
}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data <- generate_data()
  data[1:5, ]
}
\keyword{ data, survival } 
%\seealso{ \code{\link{getHazard}}, \code{\link{EM.nonparametric}} }
