#' @title Retrieves several variables from an RSA object
#'
#' @description
#' Retrieves several variables from an RSA object
#'
#' @details
#' None so far.
#'
#' @export
#' @param x RSA object
#' @param type One of: "syntax", "coef", "R2"
#' @param model A string specifying the model; defaults to "full"
#' @param ... Additional parameters passed to the extraction function
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
#' getPar(r1, "syntax")
#' getPar(r1, "R2")
#' getPar(r1, "coef")


getPar <- function(x, type="syntax", model="full", ...) {
	type <- tolower(type)
	if (type=="syntax") {
		return(x$models[[model]]@Options$model)
	}
	if (type=="coef") {
		return(parameterEstimates(x$models[[model]], ...))
	}
	if (type %in% c("r2", "rsquared", "r.squared")) {
		return(inspect(x$models[[model]], "R2", ...))
	}
	if (type %in% c("summary")) {
		return(summary(x$models[[model]], ...))
	}
	warning(paste0("Type '", type, "' not recognized!"))
}
