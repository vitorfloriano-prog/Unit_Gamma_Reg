# ------------------------------
# Unit.gamma.reg - Final Version (Corrected)
ugamma.fit <- function(formula = NULL, data = NULL,
                       Y = NULL, X = NULL,
                       intercepto = TRUE,
                       q = 1, B = 1000) {
  # Preparation: accept formula + data (like lm/glm), or Y and X directly
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) stop("Argument 'formula' must be an R formula.")
    if (is.null(data)) {
      mf <- model.frame(formula, envir = parent.frame())
      Xmat <- model.matrix(formula, data = parent.frame())
    } else {
      mf <- model.frame(formula, data = data)
      Xmat <- model.matrix(formula, data = data)
    }
    Yvec <- model.response(mf)
    Y <- matrix(Yvec, ncol = 1) # <--- CORRECTION: Ensures the dimension n x 1
    X <- as.matrix(Xmat)
    if (!intercepto && "(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  } else if (!is.null(Y) && !is.null(X)) {
    Y <- as.matrix(Y)
    if (is.data.frame(X)) X <- as.matrix(X)
    X <- as.matrix(X)
    if (intercepto) {
      # adds intercept if not existing (assumes that if the first column = 1 then it already exists)
      if (!(all(X[,1] == 1))) X <- cbind(1, X)
    }
  } else {
    stop("Provide 'formula + data' OR 'Y' and 'X' directly.")
  }
  
  # basic checks
  if (!is.numeric(Y)) stop("Error: Y must be numeric.")
  if (nrow(Y) != nrow(X)) stop("Number of observations in Y and X must match.")
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(X)
  if (q >= p) stop("'q' must be less than the number of parameters p.")
  
  ## Vector psi.0
  psi.0 <- rep(0,q)
  
  set.seed(100)
  
  ## Function to generate the sample Y ####
  vY <- function(n,alpha,phi)
  {
    ##### Generation of random numbers
    # Generating the Gamma distribution
    W <- rgamma(n,shape = phi, scale = 1/alpha)
    
    # Generating the Unit Gamma distribution
    Y <- exp(-W)
    
    return(Y)
  }
  
  ## Unrestricted Log-Likelihood Function to be maximized ####
  log.ver <- function(theta,Y)
  {
    beta <- theta[1:p]
    phi <- theta[p+1]
    #tau
    
    eta <- X%*%beta
    
    mi <- exp(eta)/(exp(eta)+1)
    # Protection against mi too close to 1, which would make the denominator 0
    # mi[mi >= 1] <- 1 - .Machine$double.eps # I added this to protect the calculation of alpha
    alpha <- ( mi^(1/phi) ) / ( 1 - mi^(1/phi) )
    
    # Protection against alpha too close to 0 or Inf
    # alpha[alpha <= 0] <- .Machine$double.eps
    
    l.beta.phi <- sum( phi*log(alpha) - log(gamma(phi)) + (alpha-1)*log(Y) +
                         (phi-1)*log(-log(Y)) )
    
    # Adding check to avoid problems in optimization
    if (!is.finite(l.beta.phi)) return(-1e10) # Returns a very low value if Inf or NA
    
    return(l.beta.phi)
  }
  
  ## Restricted Log-Likelihood Function to be maximized ####
  log.ver.H0 <- function(theta,Y)
  {
    beta <- theta[1:(p-q)]
    phi <- theta[p-q+1]
    eta <- matrix(X[,-(2:(q+1))],n,p-q)%*%beta
    mi <- exp(eta)/(exp(eta)+1)
    alpha <- ( mi^(1/phi) ) / ( 1 - mi^(1/phi) )
    l.beta.phi <- sum( phi*log(alpha) - log(gamma(phi)) + (alpha-1)*log(Y) +
                         (phi-1)*log(-log(Y)) )
    if (!is.finite(l.beta.phi)) return(-1e10)
    return(l.beta.phi)
  }
  
  ## Function to obtain the score vector ####
  vU <- function(theta,Y)
  {
    beta <- theta[1:p]
    phi <- theta[p+1]
    eta <- X%*%beta
    mi <- exp(eta)/(exp(eta)+1)
    alpha <- ( mi^(1/phi) ) / ( 1 - mi^(1/phi) )
    h.linha.eta <- exp(eta)/((exp(eta)+1)^2)
    mi.star <- alpha/(mi^((1/phi)+1))
    y.star <- ((alpha^2)*log(Y))/(phi*(mi^((1/phi)+1)))
    mT<- diag(c(h.linha.eta))
    s<- mi.star + y.star
    U.beta <- t(X)%*%mT%*%s
    u<- log(-log(Y)) - digamma(phi) - log((mi^(1/phi))/(alpha)) -
      (1/phi)*alpha*log(mi)*(1 + ((alpha*log(Y))/(phi*mi^(1/phi))))
    U.phi <- sum(u)
    U <- rbind(U.beta,U.phi)
    return(U)
  }
  
  ## Function to obtain the Fisher Information Matrix ####
  mIF <- function(theta,Y)
  {
    beta <- theta[1:p]
    phi <- theta[p+1]
    eta <- X%*%beta
    mi <- exp(eta)/(exp(eta)+1)
    alpha <- ( mi^(1/phi) ) / ( 1 - mi^(1/phi) )
    h.linha.eta <- exp(eta)/((exp(eta)+1)^2)
    mi.star <- alpha/(mi^((1/phi)+1))
    y.star <- ((alpha^2)*log(Y))/(phi*(mi^((1/phi)+1)))
    
    
    c <- (alpha*log(mi))/(mi^(1/phi))
    W.beta.beta <- diag(c(((mi.star^2)/phi)*(h.linha.eta^2)))
    W.beta.phi <- diag(c(((-mi.star)/phi)*((c/phi)+1)*h.linha.eta))
    W.phi.phi<- diag(c(trigamma(phi) + ((2*c)/(phi^2)) + ((c^2)/(phi^3))))
    K.beta.beta <- t(X)%*%W.beta.beta%*%X
    K.beta.phi <- t(X)%*%W.beta.phi%*%matrix(rep(1,n),n,1)
    K.phi.phi<- matrix(rep(1,n),1,n)%*%W.phi.phi%*%matrix(rep(1,n),n,1)
    K <- rbind(cbind(K.beta.beta,K.beta.phi),cbind(t(K.beta.phi),K.phi.phi))
    return(K)
  }
  
  ## Initial guess for theta ####
  beta.inic<- solve(t(X)%*%X)%*%t(X)%*%Y
  phi.inic<- 20
  theta.inic <- matrix(c(beta.inic,phi.inic),p+1,1)
  
  ## Log-likelihood maximization - Unrestricted ####
  Estimacao <- optim(theta.inic, log.ver, Y=Y, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T)

  theta.hat <- Estimacao$par
  psi.hat<- theta.hat[2:(q+1)]
  Verossim <- Estimacao$value
  
  ## Inverse of the Fisher information
  Inv.I.Fisher.hat <- solve(mIF(theta.hat, Y))
  
  ## Standard Error for m.l.e.
  EP <- matrix((c(sqrt(diag(Inv.I.Fisher.hat)))),p+1,1)
  
  ## Log-likelihood maximization - Restricted to H0 ####
  theta.inic.til <- c(theta.inic[1],theta.inic[(q+2):(p+1)])
  Estimacao.H0<- optim(theta.inic.til, log.ver.H0, Y=Y, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T)
  # the estimation reported that log(phi) generated NA's, is this normal?
  theta.til<- c(Estimacao.H0$par[1],psi.0,Estimacao.H0$par[2:(p-q+1)])
  Verossim.H0 <- log.ver(theta.til,Y=Y)
  
  ## Score vector evaluated at theta.til
  U.til.psi <- vU(theta.til,Y=Y)[2:(q+1)]
  
  ## Likelihood Ratio Statistic and p-value - Original
  RV<- 2*(Verossim-Verossim.H0)
  Valor.p.RV <- pchisq(RV,q,lower=F)
  
  ## Gradient Statistic and p-value - Original
  Gradiente <- c(U.til.psi%*%t(t(psi.hat-psi.0)))
  Valor.p.G <- c(pchisq(Gradiente,q,lower=F))
  
  ## Fisher Information Matrix
  IF.hat <- mIF(theta.hat,Y=Y)
  IF.til <- mIF(theta.til,Y=Y)
  
  ## Inverse of the Fisher Information Matrix
  Inv.IF.hat <- try(solve(IF.hat),T)
  Inv.IF.til <- try(solve(IF.til),T)
  Inv.IF.hat.psi.psi <- Inv.IF.hat[2:(q+1),2:(q+1)]
  Inv.IF.til.psi.psi <- Inv.IF.til[2:(q+1),2:(q+1)]
  Inv.Inv.IF.hat.psi.psi <- try(solve(Inv.IF.hat.psi.psi),T)
  
  ## Score Statistic and p-value - Original
  Escore <- U.til.psi%*%Inv.IF.til.psi.psi%*%t(t(U.til.psi))
  Valor.p.E <- c(pchisq(Escore,q,lower=F))
  
  ## Wald Statistic and p-value - Original
  Wald<- t(psi.hat-psi.0)%*%Inv.Inv.IF.hat.psi.psi%*%t(t(psi.hat-psi.0)) # WHAT IS WALD FOR?
  Valor.p.W <- c(pchisq(Wald,q,lower=F))
  
  RV.b <- NULL
  Gradiente.b <- NULL
  Escore.b <- NULL
  Wald.b <- NULL
  
  theta.hat.star <- matrix(NA,p+1,B) # Changed to NA to help with exclusion
  cont.erro.vies <- 0
  cont.erro.B <- 0
  
  cont.RV.b<- 0
  cont.Grad.b <- 0
  cont.Esc.b<- 0
  cont.Wald.b <- 0
  
  # Preparation for the loop
  beta.hat<- theta.hat[1:p]
  phi.hat<- theta.hat[p+1]
  eta.hat<- X%*%beta.hat
  mi.hat<- exp(eta.hat)/(exp(eta.hat)+1)
  alpha.hat <- (mi.hat^(1/phi.hat))/( 1 - mi.hat^(1/phi.hat))
  
  beta.til <- theta.til[1:p]
  phi.til<- theta.til[p+1]
  eta.til<- X%*%beta.til
  mi.til <- exp(eta.til)/(exp(eta.til)+1)
  alpha.til <- (mi.til^(1/phi.til))/( 1 - mi.til^(1/phi.til))
  
  ## Start of the loop - Bootstrap
  for (b in 1:B)
  {
    
    # 1. ESTIMATE FOR BIAS CORRECTION
    Y.star.hat.b<- vY(n,alpha.hat,phi.hat)
    Estimacao.b.vies <- try(optim(theta.inic,log.ver, Y = Y.star.hat.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    if(inherits(Estimacao.b.vies, "try-error") || is.null(Estimacao.b.vies$par) )
    {
      cont.erro.vies = cont.erro.vies + 1
      # Try a second time with new seed and sample
      set.seed(2021*b)
      Y.star.hat.b <- vY(n,alpha.hat,phi.hat)
      Estimacao.b.vies <- try(optim(theta.inic,log.ver, Y = Y.star.hat.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    }
    
    if(!inherits(Estimacao.b.vies, "try-error") && !is.null(Estimacao.b.vies$par)) {
      theta.hat.star[,b] <- Estimacao.b.vies$par
    } else {
      # If it still fails, ignore this iteration for bias correction
      next 
    }
    
    # 2. SAMPLING UNDER H0 AND CALCULATION OF STATISTICS
    set.seed(1625*b) # Different seed for the H0 sample
    Y.star.b<- vY(n,alpha.til,phi.til)
    
    ## Initial guess for theta - Bootstrap ####
    beta.inic.b <- solve(t(X)%*%X)%*%t(X)%*%Y.star.b
    phi.inic.b<- 20
    theta.inic.b <- matrix(c(beta.inic.b, phi.inic.b),p+1,1)
    
    ## MLE - Unrestricted - Bootstrap ####
    # Using Estimacao.b.vies$par as initial guess if successful
    initial_par_b <- if(!inherits(Estimacao.b.vies, "try-error") && !is.null(Estimacao.b.vies$par)) {
      Estimacao.b.vies$par
    } else {
      theta.inic.b
    }
    Estimacao.b <- try(optim(initial_par_b, log.ver, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    ## MLE - Restricted under H0 - Bootstrap ####
    theta.inic.b.H0 <- c(theta.inic.b[1],theta.inic.b[(q+2):(p+1)])
    Estimacao.b.H0 <- try(optim(theta.inic.b.H0, log.ver.H0, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    # MAIN CHECK TO CALCULATE STATISTICS
    if((!inherits(Estimacao.b, "try-error") && !is.null(Estimacao.b$par)) && 
       (!inherits(Estimacao.b.H0, "try-error") && !is.null(Estimacao.b.H0$par)))
    {
      theta.hat.b <- Estimacao.b$par
      psi.hat.b<- theta.hat.b[2:(q+1)]
      theta.til.b <- c(Estimacao.b.H0$par[1],psi.0,Estimacao.b.H0$par[2:(p-q+1)])
      
      Verossim.b <- log.ver(theta.hat.b, Y = Y.star.b)
      Verossim.b.H0 <- log.ver(theta.til.b, Y = Y.star.b)
      
      # LR Statistic - Bootstrap
      rv.b <- 2*(Verossim.b-Verossim.b.H0)
      
      ## Score vector evaluated at theta.til
      U.til.psi.b <- vU(theta.til.b, Y = Y.star.b)[2:(q+1)]
      
      # Gradient Statistic - Bootstrap
      Grad.b <- U.til.psi.b%*%t(t(psi.hat.b-psi.0))
      
      ## Fisher Information Matrix - Bootstrap
      IF.hat.b <- mIF(theta.hat.b, Y = Y.star.b)
      IF.til.b <- mIF(theta.til.b, Y = Y.star.b)
      
      ## Inverse of the Fisher Information Matrix - Bootstrap
      Inv.IF.hat.b <- try(solve(IF.hat.b),T)
      Inv.IF.til.b <- try(solve(IF.til.b),T)
      
      if((!inherits(Inv.IF.hat.b, "try-error"))&&(!inherits(Inv.IF.til.b, "try-error")))
      {
        Inv.IF.hat.psi.psi.b<- Inv.IF.hat.b[2:(q+1),2:(q+1), drop=FALSE] # drop=FALSE to maintain dimension
        Inv.IF.til.psi.psi.b<- Inv.IF.til.b[2:(q+1),2:(q+1), drop=FALSE] # drop=FALSE to maintain dimension
        Inv.Inv.IF.hat.psi.psi.b <- try(solve(Inv.IF.hat.psi.psi.b),T)
        
        if(!inherits(Inv.Inv.IF.hat.psi.psi.b, "try-error"))
        {
          # Score Statistic - Bootstrap
          Esc.b <- U.til.psi.b%*%Inv.IF.til.psi.psi.b%*%t(t(U.til.psi.b))
          
          # Wald Statistic - Bootstrap
          Wa.b <- t(psi.hat.b-psi.0)%*%Inv.Inv.IF.hat.psi.psi.b%*%t(t(psi.hat.b-psi.0))
          
          # Storage and Count
          RV.b <- c(RV.b, max(0, rv.b))
          Gradiente.b <- c(Gradiente.b, max(0, Grad.b))
          Escore.b <- c(Escore.b, max(0, Esc.b))
          Wald.b <- c(Wald.b, max(0, Wa.b))
          
          cont.RV.b<- cont.RV.b+ (rv.b>RV)
          cont.Grad.b <- cont.Grad.b + (Grad.b>Gradiente)
          cont.Esc.b <- cont.Esc.b + (Esc.b>Escore)
          cont.Wald.b <- cont.Wald.b + (Wa.b>Wald)
          
        } else {
          # Problem in IF inversion for Wald
          cont.erro.B = cont.erro.B + 1
        }
      } else {
        # Problem in IF inversion (IF.hat.b or IF.til.b)
        cont.erro.B = cont.erro.B + 1
      }
    }
    else{
      # Problem in optimization (Estimacao.b or Estimacao.b.H0)
      cont.erro.B = cont.erro.B + 1
      # print(paste("Error in loop", b)) # Uncomment for debug
    }
  }#End of the Bootstrap loop
  
  # Filter NA columns of theta.hat.star (failed iterations)
  theta.hat.star <- theta.hat.star[, !is.na(theta.hat.star[1,])]
  
  # Recount B for p-values and bias correction
  B_efetivo_teste <- length(RV.b)
  B_efetivo_vies <- ncol(theta.hat.star)
  
  ## p-value - Bootstrap
  Valor.p.RV.b <- cont.RV.b/B_efetivo_teste
  Valor.p.G.b <- cont.Grad.b/B_efetivo_teste
  Valor.p.E.b <- cont.Esc.b/B_efetivo_teste
  Valor.p.W.b <- cont.Wald.b/B_efetivo_teste
  
  # Bias correction via Bootstrap
  Media.theta.hat.star<- apply(theta.hat.star,1,mean)
  vies.hat <- Media.theta.hat.star - theta.hat
  theta.corr.boot<- theta.hat - vies.hat
  Inv.I.Fisher.hat.boot <- try(solve(mIF(theta.corr.boot, Y)),T) ## Inverse of the Fisher information
  EP.boot<- matrix((c(sqrt(diag(Inv.I.Fisher.hat.boot)))),p+1,1) ## Standard Error for bootstrap m.l.e.
  
  
  
  ## Unrestricted MV estimates for the parameters ####
  EMV <- theta.hat
  
  ## Restricted MV estimates to H0 for the parameters ####
  EMV.H0 <- theta.til
  
  ## Bootstrap estimates for the parameters ###
  theta.corr.b <- theta.corr.boot
  
  coefficients <- data.frame(Estimate=EMV, Std.Error= EP)
  row.names(coefficients) <- c(colnames(X), "phi")
  
  bootstrap_corrected_coefficients <- data.frame(Estimate = theta.corr.b, Std.Error.boot = EP.boot)
  row.names(bootstrap_corrected_coefficients) <- c(colnames(X),"phi")
  Tabela <- data.frame(c("LR","G", "S", "W"),
                       c(round(Valor.p.RV,3),
                         round(Valor.p.G,3),
                         round(Valor.p.E,3),
                         round(Valor.p.W,3)))
  names(Tabela) <- c("Statistic", "p-value" )
  
  Tabela_boot= data.frame(c("LR","G", "S", "W"),
                          c(round(Valor.p.RV.b,3),
                            round(Valor.p.G.b,3),
                            round(Valor.p.E.b,3),
                            round(Valor.p.W.b,3)))
  names(Tabela_boot) <- c("Statistic", "p-value.b")
  
  
  res <- list(
    n = n,
    coefficients = coefficients ,
    bootstrap_corrected_coefficients=bootstrap_corrected_coefficients,
    Tabela = Tabela,
    Tabela_boot=Tabela_boot,
    RV.b = RV.b,
    Grad.b = Grad.b,
    Escore.b = Escore.b,
    Wald.b = Wald.b,
    AIC.Gama = -2 * Verossim + 2 * p,
    BIC.Gama = -2 * Verossim + log(n) * p,
    cont.erro.B = cont.erro.B
  )
  class(res) <- "Unit.gamma"
  return(res)
}
# --------------------
# summary for Unit.gamma
# --------------------
summary.Unit.gamma <- function(object, digits = 3, ...) {
  if (!inherits(object, "Unit.gamma")) stop("Use an object returned by Unit.gamma.reg")
  cat("===========================================\n")
  cat("Model Summary - Unit Gamma (Ugamma.fit)\n")
  cat("===========================================\n")
  cat("sample size:", object$n, "\n\n")
  
  cat("coefficients:\n")
  print(round(object$coefficients, digits))
  
  cat("\ncoefficients corrected by Bootstrap:\n")
  print(round(object$bootstrap_corrected_coefficients, digits))
  
  cat("\nHypothesis Testing Table (LR, G, S, W):\n")
  print(object$Tabela, row.names = FALSE)
  
  cat("\nBootstrap hypothesis testing table (LR, G, S, W):\n")
  print(object$Tabela_boot, row.names = FALSE)
  
  
  cat("\nAIC.Gama:", round(object$AIC.Gama, 4), " BIC.Gama:", round(object$BIC.Gama, 4), "\n")
  cat("Number of errors in bootstrap:", object$cont.erro.B, "\n")
  cat("===========================================\n")
  invisible(object)
  
}
