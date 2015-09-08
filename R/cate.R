#' High dimensional factor analysis and confounder adjusted testing and estimation (CATE)
#' 
#' @description Provides several methods for factor analysis in high dimension (both n,p >> 1) and methods to adjust for possible confounders in multiple hypothesis testing.
#' 
#' @seealso \code{\link{factor.analysis}}, \code{\link{cate}}
#'
#' @docType package
#' @name cate-package
NULL

#' The main function for confounder adjusted testing
#'
#' @param X treatment variable, n*d0 matrix
#' @param Y outcome, n*p matrix
#' @param X.nuis nuisance covarites, n*d1 matrix
#' @param fa.method factor analysis method
#' @param r number of latent factors, can be estimated using the function \code{est.confounder.num}
#' @param adj.method adjustment method
#' @param nc position of the negative controls,
#'           if d0 > 1, this should be a matrix with 2 columns
#' @param psi derivative of the loss function in robust regression
#' @param nc.var.correction correct asymptotic variance based on our formula
#' @param calibrate if TRUE, use the Median and the Mean Absolute Deviation(MAD) to calibrate the test statistics
#' 
#' @return a list of objects
#' \describe{
#' \item{alpha}{estimated alpha}
#' \item{alpha.p.value}{asymptotic p-value for the global chi squared test of alpha.}
#' \item{beta}{estimated beta}
#' \item{beta.cov.row}{estimated row covariance of \code{beta}}
#' \item{beta.cov.col}{estimated column covariance of \code{beta}}
#' \item{beta.t}{asymptotic z statistics for \code{beta}}
#' \item{beta.p.value}{asymptotic p-values for beta, based on \code{beta.t}}
#' \item{Y.tilde}{the transformed outcome matrix}
#' \item{Gamma}{estimated factor loadings}
#' \item{Z}{estimated latent factors}
#' \item{Sigma}{estimated noise variance matrix}
#' }
#' 
#' @references {
#' J. Wang, Q. Zhao, T. Hastie, and A. B. Owen (2015). Confounder adjustment in multiple testing. \emph{arXiv:1508.04178}.
#' }
#'
#' @details
#' Ideally \code{nc} can either be a vector of numbers between 1 and p, if d0 = 1 or the negative controls are the same for every treatment variable, or a 2-column matrix specifying which positions of beta are known to be zero. But this is yet implemented.
#' 
#' @examples 
#' ## simulate a dataset with 100 observations, 1000 variables and 5 confounders
#' data <- gen.sim.data(n = 100, p = 1000, r = 5)
#' ## linear regression without any adjustment
#' output.naive <- cate(X = data$X1, Y = data$Y, X.nuis = NULL, r = 0, adj.method = "naive")
#' ## confounder adjusted linear regression
#' output <- cate(X = data$X1, Y = data$Y, X.nuis = NULL, r = 5)
#' ## plot the histograms of unadjusted and adjusted regression statistics
#' par(mfrow = c(1, 2))
#' hist(output.naive$beta.t)
#' hist(output$beta.t)
#' 
#' @seealso \code{\link{wrapper}} for wrapper functions of some existing methods.
#' 
#' @import MASS
#' @export
#' 
cate <- function(X,
                 Y,
                 X.nuis = NULL,
                 r,
                 fa.method = c("ml", "pc", "esa"),
                 adj.method = c("rr", "nc", "lqs", "naive"),
                 psi = psi.huber,
                 nc = NULL,
                 nc.var.correction = TRUE,
                 calibrate = TRUE) {

    # dimension check
    n <- nrow(X)
    d1 <- ncol(X)
    p <- ncol(Y)
    X <- as.matrix(X)
    if (nrow(Y) != nrow(X)) {
        stop("Dimensions of X and Y don't match!")
    }
    if (!is.null(X.nuis)) {
        d0 <- ncol(X.nuis)
    } else {
        d0 <- 0
    }
    ### comment: the check here should be if d + r >= n
    d <- d0 + d1
    if (d + r >= n) {
        stop("Observations too few to perform factor analysis!")
    }

    # match argumennts
    fa.method <- match.arg(fa.method, c("ml", "pc", "esa"))
    adj.method <- match.arg(adj.method, c("rr", "nc", "lqs", "naive"))

    # argument check
    if ((adj.method == "nc") && (is.null(nc))) {
        stop("Using negative control adjustment but no negative control is given!")
    }

    # naive tests
    if (adj.method == "naive" | r == 0) {
        naive.lm <- lm(Y ~ cbind(X.nuis, X) - 1)
		output <- list();
        s <- summary(naive.lm)
        output$beta <- sapply(1:p, function(j) {s[[j]]$coefficients[(d0+1):d,1]})
        output$beta.t <- sapply(1:p, function(j) {s[[j]]$coefficients[(d0+1):d,3]})
        if (calibrate) {
            output$beta.t <- (output$beta.t - median(output$beta.t)) / mad(output$beta.t - median(output$beta.t))
        }
        output$beta.p.value <- 2 * (1 - pnorm(abs(output$beta.t)))

        return(output)
    }

    # QR transformation
    O <- t(qr.Q(qr(cbind(X.nuis, X)), complete = TRUE))
    U <- O[1:d,] %*% cbind(X.nuis, X)
    Y.tilde <- O %*% Y
	  U11 <- U[(d0+1):d, (d0+1):d, drop = FALSE]
    Y1.tilde.scale <- solve(U11) %*% Y.tilde[(d0+1):d, , drop = FALSE]
	  X.full <- cbind(X.nuis, X)
	  Sigma.X <- t(X.full) %*% X.full / n
    
    # factor analysis
    fa.output <- factor.analysis(Y.tilde[-(1:d),],
                                 r,
                                 fa.method)

    # adjustment after factor analysis
    output <- adjust.latent(t(Y1.tilde.scale),
    						n,
							Sigma.X,
                            fa.output$Gamma,
                            fa.output$Sigma,
                            adj.method,
                            psi,
                            nc,
                            nc.var.correction)

    output$beta <- output$beta

    # compute effect size and p-value
    output$beta.t <- t(t(output$beta / sqrt(output$beta.cov.row))) / sqrt(diag(output$beta.cov.col)) *
					 sqrt(n)
    if (calibrate) { # calibrate the effect size
        output$beta.t <- (output$beta.t - median(output$beta.t)) / mad(output$beta.t - median(output$beta.t))
    }
    output$beta.p.value <- 2 * (1 - pnorm(abs(output$beta.t)))
    
    Z.tilde1 <- U[1:d, (d0+1):d, drop = F] %*% t(output$alpha)
    Z.tilde <- rbind(Z.tilde1, fa.output$Z)
    fa.output$Z <- t(O) %*% Z.tilde

    output$Y.tilde <- Y.tilde
    
    output$alpha.p.value <- pchisq(n * sum(output$alpha^2), r, lower.tail = F)

    return(c(fa.output, output))

}

#' Estimate the number of confounders
#' 
#' @inheritParams cate
#' @param method method to estimate the number of factors. There are currently two choices, 
#' "ed" is the eigenvalue difference method proposed by Onatski (2010) and "bcv" is the 
#' bi-cross-validation method proposed by Owen and Wang (2015). "bcv" tends to estimate more
#' weak factors and takes longer time
#' @param rmax the maximum number of factors to consider. If the estimated number of factors is rmax, 
#'   then users are encouraged to increase rmax and run again. Default is 20.
#' @param nRepeat the number of repeats of bi-cross-validation. A larger nRepeat will result in a 
#' more accurate estimate of the bcv error, but will need longer time to run.
#' @param bcv.plot whether to plot the relative bcv error versus the number of estimated 
#' ranks. The relative bcv error is the entrywise mean square error devided by the average of 
#' the estimated noise variance.
#' @param log if \code{log = "y"}, then the y-axis of the bcv plot is in log scale.
#'
#' @return if \code{method} is "ed", then return the estimated number of confounders/factors.
#' If \code{method} is "bcv", then return the a list of objects
#' \describe{
#' \item{r}{estimated number of confounders/factors}
#' \item{errors}{the relative bcv errors of length \code{1 + rmax}}
#' }
#'
#' @references {
#' A. B. Owen and J. Wang (2015), Bi-cross-validation for factor analysis. \emph{arXiv:1503.03515}.
#'
#' A. Onatski (2010), Determining the number of factors from empirical distribution of eigenvalues. 
#' \emph{The Review of Economics and Statistics} 92(4).
#' }
#' 
#' @examples
#' 
#' ## example for est.confounder.num
#' data <- gen.sim.data(n = 50, p = 100, r = 5)
#' est.confounder.num(data$X1, data$Y, data$X0, method = "ed")
#' est.confounder.num(data$X1, data$Y, data$X0, method = "bcv")
#' 
#' @export
#' 
est.confounder.num <- function(X,
                               Y,
                               X.nuis = NULL,
                               method = c("bcv", "ed"),
                               rmax = 20,
                               nRepeat = 12,
                               bcv.plot = TRUE, log = "") {
  
  method <- match.arg(method, c("bcv", "ed"))
  
  ## dimension check
  n <- nrow(X)
  d1 <- ncol(X)
  p <- ncol(Y)
  X <- as.matrix(X)
  if (nrow(Y) != nrow(X)) {
    stop("Dimensions of X and Y don't match!")
  }
  if (!is.null(X.nuis)) {
    d0 <- ncol(X.nuis)
  } else {
    d0 <- 0
  }
  d <- d0 + d1
  if (d >= n) {
    stop("Observations too few to perform factor analysis!")
  }
  O <- t(qr.Q(qr(cbind(X.nuis, X)), complete = TRUE))
  Y.tilde <- O %*% Y
  
 
  est.factor.num(Y.tilde[-(1:d), ], 
                 method = method, rmax = rmax, nRepeat = nRepeat, bcv.plot = bcv.plot)
}
