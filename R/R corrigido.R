# ------------------------------
# Unit.gamma.reg - Final Version - CORRIGIDA PARA PORTABILIDADE
# ------------------------------
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
    Y <- as.matrix(Yvec)
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
    eta <- X%*%beta
    
    mi <- exp(eta)/(exp(eta)+1)
    alpha <- ( mi^(1/phi) ) / ( 1 - mi^(1/phi) )
    l.beta.phi <- sum( phi*log(alpha) - log(gamma(phi)) + (alpha-1)*log(Y) +
                         (phi-1)*log(-log(Y)) )
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
    U.beta<- t(X)%*%mT%*%s
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
    W.beta.phi<- diag(c(((-mi.star)/phi)*((c/phi)+1)*h.linha.eta))
    W.phi.phi<- diag(c(trigamma(phi) + ((2*c)/(phi^2)) + ((c^2)/(phi^3))))
    K.beta.beta<- t(X)%*%W.beta.beta%*%X
    K.beta.phi<- t(X)%*%W.beta.phi%*%matrix(rep(1,n),n,1)
    K.phi.phi<- matrix(rep(1,n),1,n)%*%W.phi.phi%*%matrix(rep(1,n),n,1)
    K <- rbind(cbind(K.beta.beta,K.beta.phi),cbind(t(K.beta.phi),K.phi.phi))
    return(K)
  }
  
  
  ## Initial guess for theta ####
  beta.inic <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  z = log(Y/(1-Y))
  beta.inic <- solve(t(X)%*%X)%*%t(X)%*%z
  e.chapeu = z - X%*%solve(t(X)%*%X)%*%t(X)%*%z
  mi.chapeu = NULL
  glinha= NULL
  sigma.chapeu= NULL
  for(i in 1:length(Y))
  {
    mi.chapeu[i] = (exp(X[i,]%*%solve(t(X)%*%X)%*%t(X)%*%z))/(1+exp(X[i,]%*%solve(t(X)%*%X)%*%t(X)%*%z))
    glinha[i] = 1/(mi.chapeu[i]*(1-mi.chapeu[i]))
    sigma.chapeu[i] = (t(e.chapeu)%*%e.chapeu)/((length(Y)-p)*(glinha[i]^2))
  }
  
  phi.inic = mean((mi.chapeu*(1-mi.chapeu))/sigma.chapeu)-1
  
  theta.inic <- c(as.numeric(beta.inic), as.numeric(phi.inic))
  
  
  
  ## Maximization of the log-likelihood - Unrestricted ####
  Estimacao <- optim(theta.inic, log.ver, Y=Y, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T)
  # The estimation reported that log(alpha) generated NA's. Is this normal?
  theta.hat <- Estimacao$par
  psi.hat<- theta.hat[2:(q+1)]
  Verossim<- Estimacao$value
  
  ## AIC and BIC
  AIC.Gama <- -2*Verossim + 2*p
  BIC.Gama <- -2*Verossim + log(n)*p
  
  
  ## Inverse of the Fisher Information
  Inv.I.Fisher.hat <- solve(mIF(theta.hat, Y))
  
  
  ## Standard error for e.m.v. (Maximum Likelihood Estimator)
  EP <- matrix((c(sqrt(diag(Inv.I.Fisher.hat)))),p+1,1)
  
  
  
  
  ## Maximization of the log-likelihood - Restricted to H0 ####
  theta.inic.til <- c(theta.inic[1],theta.inic[(q+2):(p+1)])
  Estimacao.H0<- optim(theta.inic.til, log.ver.H0, Y=Y, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T)
  # The estimation reported that log(phi) generated NA's. Is this normal?
  theta.til<- c(Estimacao.H0$par[1],psi.0,Estimacao.H0$par[2:(p-q+1)])
  Verossim.H0<- log.ver(theta.til,Y=Y)
  
  
  ## Score vector evaluated at theta.til
  U.til.psi <- vU(theta.til,Y=Y)[2:(q+1)]
  
  ## Likelihood Ratio Statistic and p-value - Original
  RV<- 2*(Verossim-Verossim.H0)
  Valor.p.RV <- pchisq(RV,q,lower=F)
  
  
  # WHAT IS THE GRADIENT FOR AND WHY IS HAT USED (UNRESTRICTED) AND TIL (RESTRICTED)?
  ## Gradient Statistic and p-value - Original
  Gradiente <- c(U.til.psi%*%t(t(psi.hat-psi.0)))
  Valor.p.G <- c(pchisq(Gradiente,q,lower=F))
  
  # Why use so much information from both the restricted (til) and the unrestricted (hat)?
  
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
  Wald <- t(psi.hat-psi.0)%*%Inv.Inv.IF.hat.psi.psi%*%t(t(psi.hat-psi.0)) # WHAT IS WALD FOR?
  Valor.p.W <- c(pchisq(Wald,q,lower=F))
  
  
  RV.b<- NULL
  Gradiente.b<- NULL
  Escore.b<- NULL
  Wald.b<- NULL
  
  theta.hat.star <- matrix(0,p+1,B)
  cont.erro.vies <- 0
  cont.erro.B <- 0
  
  cont.RV.b<- 0
  cont.Grad.b<- 0
  cont.Esc.b<- 0
  cont.Wald.b<- 0
  
  ## Start of the loop - Bootstrap
  for (b in 1:B)
  {
    ## Generate the bootstrap sample (vector yb*) from theta.hat
    beta.hat<- theta.hat[1:p]
    phi.hat<- theta.hat[p+1]
    eta.hat<- X%*%beta.hat
    mi.hat<- exp(eta.hat)/(exp(eta.hat)+1)
    alpha.hat<- (mi.hat^(1/phi.hat))/( 1 - mi.hat^(1/phi.hat))
    
    Y.star.hat.b <- vY(n,alpha.hat,phi.hat)
    
    ## MLE - Unrestricted - Bootstrap ####
    Estimacao.b.vies <- try(optim(theta.inic,log.ver, Y = Y.star.hat.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    
    if((class(Estimacao.b.vies) != "try-error"))
    {theta.hat.star[,b] <- Estimacao.b.vies$par}else{
      cont.erro.vies = cont.erro.vies + 1
      set.seed(2021*b)
      Y.star.hat.b <- vY(n,alpha.hat,phi.hat)
      Estimacao.b.vies <- try(optim(theta.inic,log.ver, Y = Y.star.hat.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
      theta.hat.star[,b] <- Estimacao.b.vies$par
    }
    
    ## Generate the bootstrap sample (vector yb*)
    beta.til<- theta.til[1:p]
    phi.til<- theta.til[p+1]
    eta.til<- X%*%beta.til
    mi.til<- exp(eta.til)/(exp(eta.til)+1)
    alpha.til <- (mi.til^(1/phi.til))/( 1 - mi.til^(1/phi.til))
    
    Y.star.b <- vY(n,alpha.til,phi.til)
    
    ## Initial guess for theta - Bootstrap ####
    beta.inic.b <- solve(t(X)%*%X)%*%t(X)%*%Y.star.b
    phi.inic.b<- phi.inic
    theta.inic.b <- matrix(c(beta.inic.b, phi.inic.b),p+1,1)
    
    ## MLE - Unrestricted - Bootstrap ####
    Estimacao.b <- try(optim(theta.hat.star[,b],log.ver, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    ## MLE - Restricted under H0 - Bootstrap ####
    theta.inic.b.H0 <- c(theta.inic.b[1],theta.inic.b[(q+2):(p+1)])
    Estimacao.b.H0 <- try(optim(theta.inic.b.H0, log.ver.H0, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
    
    # CORREÇÃO 1: Garante que o IF funcione com IF, prevenindo o erro 'length > 1'
    if(all((class(Estimacao.b) != "try-error") & (class(Estimacao.b.H0) != "try-error")))
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
      Inv.IF.hat.psi.psi.b<- Inv.IF.hat.b[2:(q+1),2:(q+1)]
      Inv.IF.til.psi.psi.b<- Inv.IF.til.b[2:(q+1),2:(q+1)]
      Inv.Inv.IF.hat.psi.psi.b <- try(solve(Inv.IF.hat.psi.psi.b),T)
      
      # CORREÇÃO 2: Garante que o IF funcione com IF, prevenindo o erro 'length > 1'
      if(all((class(Inv.IF.hat.b) != "try-error") & (class(Inv.IF.til.b) != "try-error") & (class(Inv.Inv.IF.hat.psi.psi.b) != "try-error")))
      {
        # Score Statistic - Bootstrap
        Esc.b <- U.til.psi.b%*%Inv.IF.til.psi.psi.b%*%t(t(U.til.psi.b))
        
        # Wald Statistic - Bootstrap
        Wa.b <- t(psi.hat.b-psi.0)%*%Inv.Inv.IF.hat.psi.psi.b%*%t(t(psi.hat.b-psi.0))
        
        RV.b[b] <- max(0, rv.b)
        Gradiente.b[b] <- max(0, Grad.b)
        Escore.b[b] <- max(0, Esc.b)
        Wald.b[b] <- max(0, Wa.b)
        
        cont.RV.b<- cont.RV.b+ (rv.b>RV)
        cont.Grad.b <- cont.Grad.b + (Grad.b>Gradiente)
        cont.Esc.b<- cont.Esc.b + (Esc.b>Escore)
        cont.Wald.b <- cont.Wald.b + (Wa.b>Wald)
        
      }#End of the if to check for problem in IF inversion
    }
    else{
      cont.erro.B = cont.erro.B + 1
      print(Estimacao.b)
      print(b)
      print(cont.erro.B)
      
      ## Generate the bootstrap sample (vector yb*)
      beta.til <- theta.til[1:p]
      phi.til<- theta.til[p+1]
      eta.til<- X%*%beta.til
      mi.til <- exp(eta.til)/(exp(eta.til)+1)
      alpha.til <- (mi.til^(1/phi.til))/( 1 - mi.til^(1/phi.til))
      
      set.seed(1625*b)
      Y.star.b <- vY(n,alpha.til,phi.til)
      
      ## MLE - Unrestricted - Bootstrap ####
      Estimacao.b <- try(optim(theta.hat.star[,b],log.ver, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
      
      
      ## MLE - Restricted under H0 - Bootstrap ####
      Estimacao.b.H0 <- try(optim(theta.inic.til, log.ver.H0, Y = Y.star.b, NULL, method = "BFGS", control=list(fnscale=-1), hessian=T),T)
      
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
      Inv.IF.hat.psi.psi.b<- Inv.IF.hat.b[2:(q+1),2:(q+1)]
      Inv.IF.til.psi.psi.b<- Inv.IF.til.b[2:(q+1),2:(q+1)]
      Inv.Inv.IF.hat.psi.psi.b <- try(solve(Inv.IF.hat.psi.psi.b),T)
      
      # Score Statistic - Bootstrap
      Esc.b <- U.til.psi.b%*%Inv.IF.til.psi.psi.b%*%t(t(U.til.psi.b))
      
      # Wald Statistic - Bootstrap
      Wa.b <- t(psi.hat.b-psi.0)%*%Inv.Inv.IF.hat.psi.psi.b%*%t(t(psi.hat.b-psi.0))
      
      RV.b[b] <- max(0, rv.b)
      Gradiente.b[b] <- max(0, Grad.b)
      Escore.b[b] <- max(0, Esc.b)
      Wald.b[b]<- max(0, Wa.b)
      
      cont.RV.b<- cont.RV.b+ (rv.b>RV)
      cont.Grad.b <- cont.Grad.b + (Grad.b>Gradiente)
      cont.Esc.b <- cont.Esc.b + (Esc.b>Escore)
      cont.Wald.b <- cont.Wald.b + (Wa.b>Wald)
      
      
    } #End of the if to check for problem in estimation
  }#End of Bootstrap loop
  
  ## p-value - Bootstrap
  Valor.p.RV.b <- cont.RV.b/B
  Valor.p.G.b <- cont.Grad.b/B
  Valor.p.E.b <- cont.Esc.b/B
  Valor.p.W.b <- cont.Wald.b/B
  
  # Bias correction via Bootstrap
  Media.theta.hat.star <- apply(theta.hat.star,1,mean)
  vies.hat <- Media.theta.hat.star - theta.hat
  theta.corr.boot <- theta.hat - vies.hat
  Inv.I.Fisher.hat.boot <- solve(mIF(theta.corr.boot, Y)) ## Inverse of the Fisher Information
  EP.boot<- matrix((c(sqrt(diag(Inv.I.Fisher.hat.boot)))),p+1,1) ## Standard error for bootstrap e.m.v.
  
  
  
  ## Unrestricted MLE estimates for the parameters ####
  EMV <- theta.hat
  
  
  
  ## Restricted MLE estimates to H0 for the parameters ####
  EMV.H0 <- theta.til
  
  ## Bootstrap Estimates for the parameters ###
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