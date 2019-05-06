
Delta.fn <- function(t, t.fail.o, N, n.failure)
{
   tmp <- matrix(rep(t.fail.o, each=N), nrow=n.failure, ncol=N, byrow=TRUE)
   ret <- t(t(tmp) == t)
   ret
}

M.fn <- function(t, t.all, N, n.failure)
{
  rc  <- -1
  M   <- as.numeric(rep(0, n.failure*N))
  tmp <- .C("C_M_fn", as.numeric(t), as.numeric(t.all), as.integer(N), 
           as.integer(n.failure), retCode=as.integer(rc), retVec=M,
           PACKAGE="Immunotherapy.Design")
  M   <- matrix(tmp$retVec, nrow=n.failure, ncol=N, byrow=TRUE)

  M
}

M1.fn <- function (lambda, t, t.all, N, n.failure, Mfn=NULL)
{ 
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)
  as.vector(t(lambda)%*%Mfn)
}

Mstar.fn <- function(tstar, t.all, n.failure)
{
  tvec   <- t.all[1:n.failure]
  tmp1   <- tvec <= tstar
  tmp2   <- pmin(rep(tstar, n.failure), t.all[-1]) - tvec
  M      <- tmp1*tmp2
  dim(M) <- c(n.failure, 1)
  
  M
}

Mstar1.fn <- function(lambda, tstar, t.all, n.failure, Mstarfn=NULL)
{
  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  t(lambda)%*%Mstarfn
}

M2.fn <- function(lambda, t, tstar, t.all, N, n.failure, M1fn=NULL, Mstar1fn=NULL)
{
  if (is.null(M1fn)) M1fn = M1.fn(lambda, t, t.all, N, n.failure)
  if (is.null(Mstar1fn)) Mstar1fn = Mstar1.fn(lambda, tstar, t.all, n.failure)
  M1fn - rep(Mstar1fn, N)
}

A.fn <- function(Zt, tstar, t.all, t, N, n.failure, Mstarfn=NULL, Mfn=NULL)
{
  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)

  tmp0     <- Zt*(tstar <= t)
  tmp      <- tmp0
  dim(tmp) <- c(N, 1)
  mat      <- Mfn%*%(1-tmp)
  tmp2     <- tmp0*matrix(rep(Mstarfn, each=N), nrow=N, ncol=n.failure, byrow=FALSE)
  M        <- colSums(tmp2) + as.vector(mat)

  M
}

B.fn <- function(beta, t, t.all, tstar, Zt, N, n.failure, Mstarfn=NULL, Mfn=NULL)
{
  if (is.null(Mstarfn)) Mstarfn = Mstar.fn(tstar, t.all, n.failure)
  if (is.null(Mfn)) Mfn = M.fn(t, t.all, N, n.failure)
 
  tmp      <- Zt*(tstar <= t)
  dim(tmp) <- c(N, 1)
  mat      <- Mfn - matrix(rep(Mstarfn, each=N), nrow=n.failure, ncol=N, byrow=TRUE)
  M        <- exp(beta)*as.vector((mat %*% tmp))

  M
}

C.fn <- function(Zt, t, tstar, event_status)
{
  sum(Zt*(t>=tstar)*event_status)
}

D.fn <- function(Zt, t, tstar, lambda, t.all, N, n.failure, M2fn=NULL)
{
  if (is.null(M2fn)) M2fn = M2.fn(lambda, t, tstar, t.all, N, n.failure)
  ret <- sum(Zt*(t>=tstar)*M2fn)
  ret
}

loglik.fn = function(X, trt, event_status,
                     tstar,
                     t.all,
                     beta,
                     lambda, 
                     Zt,
                     effect_p, 
                     N, n.failure)
{

  Mstar1fn = as.vector(Mstar1.fn(lambda, tstar, t.all, n.failure, Mstarfn=NULL))
  M1fn     = M1.fn(lambda, X, t.all, N, n.failure, Mfn=NULL)
  trteq1   <- trt == 1
  Xge      <- X >= tstar
  tmp      <- Zt*Xge

  r1     = sum(log(lambda))
  r2     = -sum(tmp)*Mstar1fn
  r3     = -t(M1fn)%*%(1-tmp)
  r4     = beta*C.fn(Zt, t=X, tstar, event_status = event_status)
  r5     = -(exp(beta))*D.fn(Zt, t=X, tstar, lambda, t.all, N, n.failure, M2fn=NULL)
  r6     = sum(Zt*trteq1*log(effect_p))
  r7     = sum((1-Zt)*trteq1*log(1-effect_p))
  loglik = r1 + r2 + r3 + r4 + r5 + r6 + r7  
  loglik
}

pdf.r.fn <- function(lambda, t, Deltafn, event_status, beta, tstar, N, n.failure, t.all)
{

  vec      <- 1:N
  Mstar1fn = as.vector(Mstar1.fn(lambda, tstar, t.all, n.failure, Mstarfn=NULL))
  M2fn     = M2.fn(lambda,t,tstar,t.all,N, n.failure, M1fn=NULL, Mstar1fn=NULL)[vec]

  r1    <- apply(lambda^Deltafn, 2, prod)[vec]
  tmp1  <- (t >= tstar)[vec]
  r2    <- exp(beta*event_status[vec]*tmp1)
  m1fn  <- M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL)
  r3    <- exp(-(1-tmp1)*m1fn)
  r4    <- exp((-Mstar1fn-exp(beta)*M2fn)*tmp1)
  pdf.r <- (r1*r2*r3*r4)

  pdf.r
}

pdf.nr.fn <- function(lambda, t, Deltafn, N, n.failure, t.all)
{

  r1     <- apply(lambda^Deltafn, 2, prod)[1:N]
  r5     <- exp(-M1.fn(lambda,t,t.all,N, n.failure, Mfn=NULL))
  pdf.nr <- r1*r5

  pdf.nr

}

get_tfailo <- function(X, event_status) {

  t.failure     <- X[event_status == 1]
  tmp           <- order(t.failure)
  t.fail.o      <- t.failure[tmp]

  t.fail.o

} # END: get_tfailo

####################################################
# !!! Data must be ordered befor calling EM.main !!!
####################################################
EM.main <- function(X, trt, event_status, Zt, effect_p, t1,
               lambda0, stopTol=1e-4, maxiter=10000, print=FALSE) {
  # X is time (t)
  N             <- length(trt)
  trtEq1        <- trt == 1
  n_trt         <- sum(trtEq1)
  trtEq0        <- trt == 0
  ll0           <- NA
  iter          <- 0
  conv          <- FALSE
  t.fail.o      <- get_tfailo(X, event_status)
  n.failure     <- length(t.fail.o)
  t.all         <- c(0, t.fail.o)
  if (is.null(Zt)) Zt <- ifelse(trtEq1, 0.5, 0)
  if (is.null(lambda0)) lambda0 <- getHazard(X, trt, event_status)
  lambda        <- lambda0
  if (n.failure != length(lambda0)) stop("ERROR: lambda0 has the wrong length")

  while (iter < maxiter) {
    iter <- iter + 1
    
    #########################################
    # M-step: update baseline hazard and log(HR)
    #########################################
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## update beta=log(HR)                               ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

    cfn  <- C.fn(Zt, X, t1, event_status)
    dfn  <- D.fn(Zt, X, t1, lambda, t.all, N, n.failure, M2fn=NULL)
    beta <- log(cfn/dfn)
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## Solve for gamma and update lambda                 ##
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    Al <- A.fn(Zt, t1, t.all, X, N, n.failure, Mstarfn=NULL, Mfn=NULL)
    Bl <- B.fn(beta, X, t.all, t1, Zt, N, n.failure, Mstarfn=NULL, Mfn=NULL)

    lambda <- (1/(Al+Bl))[1:n.failure]

    #########################################
    # E-step: update Zt responder status
    #########################################
    Zt.trt  <- Zt[trtEq1]
    Deltafn <- Delta.fn(X, t.fail.o, N, n.failure)
    pdf.r   <- pdf.r.fn(lambda, X, Deltafn, event_status, beta, t1, n_trt, n.failure, t.all)
    pdf.nr  <- pdf.nr.fn(lambda, X, Deltafn, n_trt, n.failure, t.all)
    Zt.trt  <- (effect_p*pdf.r/(effect_p*pdf.r + (1-effect_p)*pdf.nr))[1:n_trt]
    Zt      <- c(Zt.trt, Zt[trtEq0])
    
    #########################################
    # Update loglik with the updated parameters
    #########################################
    loglik <- loglik.fn(X, trt, event_status, t1, t.all, beta, lambda, 
                        Zt, effect_p, N, n.failure)
  
    if (iter > 1) {
      diff <- abs(loglik - ll0)
      if (diff <= stopTol) {  
        conv <- TRUE
        break
      }
    }
    ll0 <- loglik
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
  
  tmp <- cbind(t.fail.o, lambda)
  colnames(tmp) <- c("EventTime", "Lambda")
  pr <- effect_p*pdf.r/(effect_p*pdf.r + (1-effect_p)*pdf.nr)
  pr[trtEq0] <- Zt[trtEq0]
  result <- list(converged=conv, logHR=beta, baseline=tmp, probResponder=pr)

  result
}

getHazard <- function(time, treatment, event_status) {

  t.fail.o   <- get_tfailo(time, event_status)
  n.failure  <- length(t.fail.o)
  t.all      <- c(0, t.fail.o)
  t.diff     <- t.all[-1]-t.all[-length(t.all)]
  fit        <- coxph(Surv(time, event_status) ~ treatment)
  ss         <- survfit(fit)
  cum.hazard <- unique(-log(ss$surv))
  if (length(cum.hazard) != n.failure) {
    tmp <- cum.hazard != 0
    cum.hazard <- cum.hazard[tmp]
    if (length(cum.hazard) != n.failure) stop("ERROR computing hazard")
  }

  hazard.dis   = c(cum.hazard[1], cum.hazard[-1] - cum.hazard[-length(cum.hazard)])
  hazard       = hazard.dis / t.diff

  hazard

} # END: getLambda

Pembedded.EM.NP <- function(data, time.var="X", trt.var="trt", status.var="event_status",
               effect_p=0.6, t1=1, lambda0=NULL, probResponder=NULL, 
               stopTol=1e-4, maxiter=10000, print=0) {

  if ( (!is.data.frame(data)) && (!is.matrix(data)) ) stop("ERROR: data must be a data frame or matrix")
  # data must be ordered by treatment
  data         <- data[order(data[, trt.var], decreasing = TRUE), , drop=FALSE]
  treatment    <- data[, trt.var]
  event_status <- data[, status.var]
  time         <- data[, time.var]
  n            <- length(time)
  if (length(treatment) != n) stop("ERROR: length(treatment) != length(time)")
  if (length(event_status) != n) stop("ERROR: length(event_status) != length(time)")
  tmp <- (!is.finite(time)) | (time < 0)
  tmp[is.na(tmp)] <- TRUE
  if (any(tmp)) stop("ERROR: with time vector")
  if (!all(treatment %in% 0:1)) stop("ERROR: with treatment vector")
  if (!all(event_status %in% 0:1)) stop("ERROR: with event_status vector")
  if ((effect_p <= 0) || (effect_p >= 1)) stop("ERROR: effect_p must be greater than 0 and less than 1")
  if (t1 <= 0) stop("ERROR: t1 must be positive")
  if (maxiter < 1) stop("ERROR: maxiter must be greater than 0")
  if (stopTol <= 0) stop("ERROR: stopTol must be positive")
  if (is.null(lambda0)) lambda0 <- getHazard(time, treatment, event_status)
  if (is.null(probResponder)) probResponder <- ifelse(treatment == 1, 0.5, 0)
  if (length(probResponder) != n) stop("ERROR: length(probResponder) != length(time)")

  ret <- EM.main(time, treatment, event_status, probResponder, effect_p, t1,
               lambda0, stopTol=stopTol, maxiter=maxiter, print=print)

  ret

} # END: Pembedded.EM.NP

Pembedded.ReRandomizationTest.NP <- function(data, time.var="X", trt.var="trt", 
              status.var="event_status", effect_p=0.6, t1=1, stopTol=1e-4, 
              maxiter=10000, print=0, num_rand=100) {

  # EM estimate parameters of interest and Zt
  EM.est    <- Pembedded.EM.NP(data, time.var=time.var, trt.var=trt.var, status.var=status.var,
                  effect_p=effect_p, t1=t1, lambda0=NULL, probResponder=NULL,
                  stopTol=stopTol, maxiter=maxiter, print=print)
  logHR.obs <- EM.est$logHR
  if (print) cat(paste("Observed log(HR) = ", logHR.obs, "\n", sep=""))
   
  # shuffle the treatment labels in the observed data and obtain the 
  # re-randomization distribution of lambda.hat on the shuffled data
  logHR.rand.all <- rep(NA, num_rand)
  nr             <- nrow(data)
  trt.obs        <- as.vector(data[, trt.var])
  time.obs       <- as.vector(data[, time.var])
  status.obs     <- as.vector(data[, status.var])

  for (j in 1:num_rand) {
    if (print > 1) cat(paste("Randomization ", j, "\n", sep=""))
    trtVec <- trt.obs[sample(nr)]
    ord    <- order(trtVec, decreasing=TRUE)
    tmp    <- EM.main(time.obs[ord], trtVec[ord], status.obs[ord], NULL, effect_p, t1,
                   NULL, stopTol=stopTol, maxiter=maxiter, print=0)
    if (tmp$converged) logHR.rand.all[j] = tmp$logHR
  }
  m <- sum(!is.na(logHR.rand.all))
  if (print) {  
    if (m < num_rand) cat("NOTE: EM algorithm did not converge for all randomizations\n", sep="")
    cat(paste("P-value is based on ", m, " randomizations\n", sep=""))
  }
  if (!m) stop("ERROR: p-value could not be estimated")
  p.val.rand <- 1 - mean(logHR.rand.all > logHR.obs, na.rm=TRUE)
  result     <- list(p.val.rerand=p.val.rand, baseline=EM.est$baseline, 
                     logHR=logHR.obs)

  return(result)

} # END: Pembedded.ReRandomizationTest.NP

Pow.Pembedded.NP <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0) {

  p.val.all <- rep(NA, nsim) 
  for (i in 1:nsim)
  {  
    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    data.o <- generate_data(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
                 enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1)

    ret <- try(Pembedded.ReRandomizationTest.NP(data.o, effect_p=effect_p, t1=t1, 
              stopTol=stopTol, maxiter=maxiter, print=0, num_rand=num_rand), silent=FALSE)
    if (!("try-error" %in% class(ret))) p.val.all[i] <- ret$p.val.rerand
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  power <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=power, n=n)

} # END: Pow.Pembedded.NP

