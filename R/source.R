generate_data <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                          lambda1=0.117, HR=0.5, tau=12*5, t1=1) {
  
  # output data
  dat <- as.data.frame(matrix(0, nrow=nmax, ncol=9))
  names(dat) <- c("id", "trt", "Z", "tau", "enroll_time", "time_to_event", "event_status", "X", "t1")
  
  dat$id = 1:nmax
  
  # treatment allocation: 1 = experimental arm
  # sort by treatment assignment status to make the responders only arise among treated subjects 
  # set.seed(seed)
  dat$trt <- sort(rbinom(n = nmax, size = 1, prob=rand_ratio), decreasing = TRUE)
  n_trt = sum(dat$trt == 1)
  n_cnt = sum(dat$trt == 0)
  
  # simulate the response statue among the experimental arm with proportion effect_P: constraint: proportion of responders = effect_p
  # set.seed(seed)
  Z.trt = rbinom(n = n_trt, size = 1, prob=effect_p)
  Z.cnt = rep(0, n_cnt)
  dat$Z = c(Z.trt, Z.cnt)
  
  # total study duration 
  dat$tau = tau
  
  # enrollment follows expected trajectory from Poisson process with memoryless property
  ## waiting time between two consective subjects follows an Exponential distribution with enrollment rate
  ## arrival time for each subject is converted from the waitig by taking the cumulative sum of the waiting times
  waiting_times <- rexp(n = nmax, rate=enroll_rate)
  dat$enroll_time <- cumsum(waiting_times)
  
  # t*
  dat$t1 <- t1
  
  # time to event
  # set.seed(seed)
  n_resp  <- sum(dat$Z == 1) # number of responders among the experimental and control groups
  n_nresp <- sum(dat$Z == 0) # number of non-responders among the experimental and control groups
  # dat$time_to_event[dat$Z == 1] <- rexp(n_resp, rate=lambda1*HR)
  dat$time_to_event[dat$Z == 1] <- rpexp(n_resp, rate=c(lambda1, HR*lambda1), t=c(0, dat$t1[1]))
  dat$time_to_event[dat$Z == 0] <- rexp(n_nresp, rate=lambda1)
  
  # build in the administrative censoring 
  dat$event_status <- ifelse(dat$time_to_event <= tau - dat$enroll_time, 1, 0)
  
  # observational time
  dat$X <- ifelse(dat$time_to_event <= tau - dat$enroll_time, dat$time_to_event, tau - dat$enroll_time)
  if (any(dat$X <= 0)) warning("Some observational times are <= 0. Try increasing tau or changing other parameters.")  

  row.names(dat) <- NULL
  
  dat
}

