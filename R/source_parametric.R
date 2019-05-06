Pembedded.EM.P <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, probResponder=NULL, 
               stopTol=1e-5, maxiter=10000, print=0) {

  if ( (!is.data.frame(data)) && (!is.matrix(data)) ) stop("ERROR: data must be a data frame or matrix")
  treatment    <- as.numeric(data[, trt.var])
  event_status <- as.numeric(data[, status.var])
  time         <- as.numeric(data[, time.var])
  n            <- length(time)
  tmp <- (!is.finite(time)) | (time < 0)
  tmp[is.na(tmp)] <- TRUE
  if (any(tmp)) stop("ERROR: with time vector")
  if (!all(treatment %in% 0:1)) stop("ERROR: with treatment vector")
  if (!all(event_status %in% 0:1)) stop("ERROR: with event_status vector")
  if ((effect_p <= 0) || (effect_p >= 1)) stop("ERROR: effect_p must be greater than 0 and less than 1")
  if (length(t1) != 1) stop("ERROR: t1 must have length 1")
  if (t1 <= 0) stop("ERROR: t1 must be positive")
  if (maxiter < 1) stop("ERROR: maxiter must be greater than 0")
  if (stopTol <= 0) stop("ERROR: stopTol must be positive")
  if (is.null(probResponder)) probResponder <- ifelse(treatment == 1, 0.5, 0)
  if (length(probResponder) != n) stop("ERROR: length(probResponder) != length(time)")

  ret <- EM.P.main(time, treatment, event_status, probResponder, effect_p, t1,  
                      stopTol=stopTol, maxiter=maxiter, print=print)

  ret

} # END: Pembedded.EM.P

################################################################
# log-likelihood function 
################################################################
loglik.parametric = function(X, Zt, exst1, nevents, xt1, xt12, 
               trt1LogEffect, trt1Log1mEffect, lambda, h.nr) {
                     
  r1     = sum(Zt*exst1)
  r2     = nevents
  r3     = sum(Zt*xt1)
  r4     = sum(Zt*xt12)
  r5     = sum((1-Zt)*X)
  r6     = sum(Zt*trt1LogEffect)
  r7     = sum((1-Zt)*trt1Log1mEffect)
  loglik = (r1*log(lambda) + r2*log(h.nr) - h.nr*(r3 + lambda*r4 + r5) + r6 + r7) 
 
  loglik
}

################################################################
# EM Algorithm
# num_rep     = number of replications in the bias assessment
# nmax        = maximum sample size
# rand_ratio  = probability of assignment to experimental arm
# effect_P    = Proportion of responders in the treatment arm
# HR          = hazard ratio
# enroll_rate = enrollment rate in subjects per month
# tau         = total study duration
################################################################

EM.P.main <- function(X, trt, event_status, Zt, effect_p, t1, 
                      stopTol=0.00001, maxiter=1000, print=0) {

  N <- length(trt)
  if (is.null(Zt)) {
    Zt           <- rep(0, N)
    Zt[trt == 1] <- 0.5
  }

  ret_conv   <- as.integer(0)
  ret_lambda <- as.numeric(-9999)
  ret_h_nr   <- as.numeric(-9999)

  tmp <- .C("C_EM_mainP", as.integer(N), as.numeric(X), as.integer(trt), 
             as.integer(event_status), as.numeric(effect_p), as.numeric(t1), 
             as.numeric(stopTol), as.integer(maxiter), as.integer(print), 
             ret_conv=ret_conv, ret_lambda=ret_lambda, ret_h_nr=ret_h_nr, 
             Zt=as.numeric(Zt), PACKAGE="Immunotherapy.Design")
  
  result <- list(converged=as.logical(tmp$ret_conv), lambda=tmp$ret_lambda, 
                 baseline=tmp$ret_h_nr, probResponder=tmp$Zt)

  return(result)

} # END: EM.P.main

EM.P.main0 <- function(X, trt, event_status, Zt, effect_p, t1, 
                      stopTol=0.00001, maxiter=1000,  print=0) {
  N               <- length(trt)
  trtEq1          <- trt == 1
  trtEq0          <- trt == 0
  Xminust1        <- X - t1
  tmp             <- X > t1
  exst1           <- event_status*tmp
  nevents         <- sum(event_status) 
  X_LT_t1         <- X < t1
  xt1             <- rep(t1, N)
  xt1[X_LT_t1]    <- X[X_LT_t1]
  xt12            <- Xminust1*tmp
  eventEq1        <- event_status == 1
  eventEq0        <- event_status == 0
  trt1LogEffect   <- trtEq1*log(effect_p)
  trt1Log1mEffect <- trtEq1*log(1-effect_p) 
  loglik0         <- 0
  iter            <- 0
  conv            <- FALSE
  if (is.null(Zt)) {
    Zt         <- rep(0, N)
    Zt[trtEq1] <- 0.5
  }

  while (iter < maxiter)
  {
    iter <- iter + 1
    #########################################
    # M-step: estimate baseline hazard and HR
    #########################################
    r1     = sum(Zt*exst1)
    r2     = nevents
    r3     = sum(Zt*xt1)
    r4     = sum(Zt*xt12)
    r5     = sum((1-Zt)*X)
    
    h.nr   = (r2 - r1)/(r3 + r5)
    lambda = r1/(r4*h.nr)
    
    #########################################
    # E-step: update Zt responder status
    #########################################

    Zt0           <- Zt
    tmp           <- lambda*h.nr  
    vec           <- rep(tmp, N)
    vec[X_LT_t1]  <- h.nr 
    vec[eventEq0] <- 1
    vec2          <- -h.nr*t1 - tmp*Xminust1
    vec2[X_LT_t1] <- -h.nr*(X[X_LT_t1])
    pdf.r         <- vec*exp(vec2)
    vec           <- exp(-h.nr*X)
    vec[eventEq1] <- vec[eventEq1]*h.nr
    pdf.nr        <- vec
    vec           <- effect_p*pdf.r
    Zt            <- vec/(vec + (1-effect_p)*pdf.nr)
    Zt[trtEq0]    <- Zt0[trtEq0]

    #########################################
    # Update loglik with the updated parameters
    #########################################

    loglik <- loglik.parametric(X, Zt, exst1, nevents, xt1, xt12, 
               trt1LogEffect, trt1Log1mEffect, lambda, h.nr)
    if (iter > 1) {
      diff <- abs(loglik - loglik0)
      if (diff <= stopTol) {  
        conv <- TRUE
        break
      }
    }
    loglik0 <- loglik
    if (print > 1) {
      str <- paste("Iter=", iter, ", loglike=", round(loglik, digits=4), sep="")
      if (iter > 1) str <- paste(str, ", diff=", round(diff, digits=8), sep="")
      str <- paste(str, "\n", sep="")
      cat(str)
    }    
  }

  if (print) {
    if (conv) {
      cat(paste("EM algoritm converged in ", iter, " iteration(s)\n", sep=""))
    } else {
      cat("EM algoritm did not converge\n")
    }
  }
  
  result <- list(converged=conv, lambda=lambda, baseline=h.nr, probResponder=Zt)

  return(result)

} # END: EM.P.main

Pembedded.ReRandomizationTest.P <- function(data, time.var="X", trt.var="trt", 
               status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-5, 
               maxiter=10000, print=0, num_rand=10000) {

  # EM estimate parameters of interest and Zt
  EM.est <- Pembedded.EM.P(data, time.var=time.var, trt.var=trt.var, status.var=status.var,
               effect_p=effect_p, t1=t1, probResponder=NULL, 
               stopTol=stopTol, maxiter=maxiter, print=print)
  lambda.obs <- EM.est$lambda

  if (print) cat(paste("Observed lambda = ", lambda.obs, "\n", sep=""))
   
  # shuffle the treatment labels in the observed data and obtain the 
  # re-randomization distribution of lambda.hat on the shuffled data
  nr        <- nrow(data)
  ret_p     <- as.numeric(-1)
  ret_nrand <- as.integer(0)
  prt       <- 0

  tmp <- .C("C_ReRandP", as.integer(num_rand), as.integer(nr), as.numeric(data[, time.var]), 
            as.integer(data[, trt.var]), as.integer(data[, status.var]), as.numeric(effect_p), 
            as.numeric(t1), as.numeric(stopTol), as.integer(maxiter),
            as.integer(prt), as.numeric(lambda.obs), ret_nrand=ret_nrand, ret_p=ret_p,
            PACKAGE="Immunotherapy.Design")
  m   <- tmp$ret_nrand
  if (print) {  
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  if (!m) stop("ERROR: p-value could not be estimated")
 
  result <- list(p.val.rerand=tmp$ret_p, baseline=EM.est$baseline, lambda=lambda.obs)

  return(result)


} # END: Pembedded.ReRandomizationTest.P

Pembedded.ReRandomizationTest.P0 <- function(data, time.var="X", trt.var="trt", 
               status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-5, 
               maxiter=10000, print=0, num_rand=1000) {

  # EM estimate parameters of interest and Zt
  EM.est <- Pembedded.EM.P(data, time.var=time.var, trt.var=trt.var, status.var=status.var,
               effect_p=effect_p, t1=t1, probResponder=NULL, 
               stopTol=stopTol, maxiter=maxiter, print=print)
  lambda.obs <- EM.est$lambda

  if (print) cat(paste("Observed lambda = ", lambda.obs, "\n", sep=""))
   
  # shuffle the treatment labels in the observed data and obtain the 
  # re-randomization distribution of lambda.hat on the shuffled data
  lambda.rand.all <- rep(NA, num_rand)
  nr              <- nrow(data)
  trt.obs         <- as.numeric(data[, trt.var])
  time.obs        <- as.numeric(data[, time.var])
  status.obs      <- as.numeric(data[, status.var])

  for (j in 1:num_rand) {
    if (print > 1) cat(paste("Randomization ", j, "\n", sep=""))
    trtVec <- trt.obs[sample(nr)]
    tmp    <- EM.P.main(time.obs, trtVec, status.obs, NULL, effect_p, t1, 
                      stopTol=stopTol, maxiter=maxiter,  print=0)
    if (tmp$converged) lambda.rand.all[j] = tmp$lambda
  }
  if (print) {
    m <- sum(!is.na(lambda.rand.all))
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  p.val.rand <- 1 - mean(lambda.rand.all > lambda.obs, na.rm=TRUE)
  result     <- list(p.val.rerand=p.val.rand, baseline=EM.est$baseline,
                     lambda=lambda.obs)

  return(result)

} # END: Pembedded.ReRandomizationTest.P0

Pow.Pembedded.P <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0) {

  p.val.all <- rep(NA, nsim) 
  for (i in 1:nsim)
  {  
    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    data.o <- generate_data(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
                 enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1)

    ret <- try(Pembedded.ReRandomizationTest.P(data.o, effect_p=effect_p, t1=t1, 
              stopTol=stopTol, maxiter=maxiter, print=0, num_rand=num_rand), silent=FALSE)
    if (!("try-error" %in% class(ret))) p.val.all[i] <- ret$p.val.rerand
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  power <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=power, n=n)

} # END: Pow.Pembedded.P

N.Pembedded.P <- function(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, min.N=100, max.N=700, 
                       tol.power=0.01, tol.N=1, print=1) {

  if (min.N < 10) stop("ERROR: min.N must be at least 10")
  #if (max.N > 1e3) stop("ERROR: max.N cannot be larger than 1e3")
  if (max.N < min.N) stop("ERROR: max.N cannot be smaller than min.N")
  if ((power < 0) || (power > 1)) stop("ERROR: power must be between 0 and 1")
  if (tol.N < 1) stop("ERROR: tol.N must be at least 1")

  # Test endpoints
  pwr1 <- Pow.Pembedded.P(nmax=min.N, rand_ratio=rand_ratio, effect_p=effect_p, 
            enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
            maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim)$power
  if (print) cat(paste("min.N = ", min.N, " power = ", pwr1, "\n", sep=""))
  if (pwr1 >= power) return(list(sampleSize=min.N, power=pwr1)) 
  pwr2 <- Pow.Pembedded.P(nmax=max.N, rand_ratio=rand_ratio, effect_p=effect_p, 
            enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
            maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim)$power
  if (print) cat(paste("max.N = ", max.N, " power = ", pwr2, "\n", sep=""))
  if (pwr2 <= power) return(list(sampleSize=max.N, power=pwr2))

  diff1    <- abs(pwr1 - power)
  diff2    <- abs(pwr2 - power)
  step     <- 100
  N0       <- min.N
  N1       <- max.N
  if (diff1 < diff2) {
    N      <- min(N0 + step, N1)
  } else {
    N      <- floor((N0 + N1)/2)
  }
  iter     <- 0
  min.diff <- 1e100
  min.N    <- 0
  min.pwr  <- -1
  while (1) {
    pwr  <- Pow.Pembedded.P(nmax=N, rand_ratio=rand_ratio, effect_p=effect_p, 
              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
              maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, 
              nsim=nsim)$power
    if (print) cat(paste("N = ", N, " power = ", pwr, "\n", sep=""))
    diff <- abs(power - pwr)
    if (diff < min.diff) {
      min.diff <- diff
      min.N    <- N
      min.pwr  <- pwr
    }
    if (diff <= tol.power) break
    if (pwr < power) {
      N0 <- N
    } else {
      N1 <- N
    }
    N  <- floor((N0 + N1)/2) 
    if (abs(N1 - N0) <= tol.N) {
      min.N   <- N
      min.pwr <- Pow.Pembedded.P(nmax=N, rand_ratio=rand_ratio, effect_p=effect_p, 
                  enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, 
                  maxiter=maxiter, stopTol=stopTol, alpha=alpha, num_rand=num_rand, 
                  nsim=nsim)$power
      break
    }
  }
  
  list(sampleSize=min.N, power=min.pwr)

} # END: N.Pembedded.P


