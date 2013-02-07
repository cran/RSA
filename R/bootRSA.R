#' @title Compute bootstrap replications for model parameteres
#' @aliases CI.boot 
#' @description
#' Compute bootstrap replications for model parameteres
#'
#' @details
#' None so far.
#'
#' @export
#' @param x RSA object
#' @param model A string specifying the model; defaults to "full"
#' @param ... Additional parameters passed to the bootstrapLavaan function
#'
#' @seealso \code{\link{RSA}}
#'
#' @examples
#'
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	sqdiff <- (x-y)^2
#' 	z.sq <- sqdiff + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df)
#' \dontrun{
#' b1 <- bootRSA(r1, model="SSD", R=5000, parallel="multicore", ncpus=2)
#' r1.boot.CI <- CI.boot(b1)	# compute percentile confidence intervals and percentile p-value of the bootstrapped values
#'}

bootRSA <- function(x, model="full", ...) {
	data.frame(bootstrapLavaan(x$models[[model]], FUN=function(x) return(coef(x, type="all")), ...))
}


#' @export
CI.boot <- function(r.boot) {
	CIs <- apply(r.boot, 2, function(x) {
		p <- sum(x<0)/length(x)
		qu <- quantile(x, probs=c(.025, .975))
		res <- c(qu, p.value=min(p, 1-p)*2)	# *2 to make p-values two-sided
		return(res)
	})
	
	CIs <- t(CIs)
	return(CIs)
}