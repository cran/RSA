#' @title Compare RSA models
#'
#' @description
#' Compare several fit indexes of all models computed from the RSA function
#'
#' @details
#' No details so far.
#'
#' @export
#' @param x An RSA object
#' @param verbose Should the summary be printed?




compare <- function(x, verbose=TRUE) {
	library(plyr)
	
	with(x$models, {
	
	res <- data.frame()
	
	if (!is.null(full)) {
		free.max <- getFreeParameters(full)
	
		if (verbose==TRUE) {
			cat("-------------------------------------------------------------------------\n")	
			cat("Standard polynomial models:\n")
			cat("-------------------------------------------------------------------------\n")
		}
	
		aL1 <- anovaList(list(cubic=cubic, full=full, IA=IA, additive=additive, diff=diff, null=null))
		if (aL1$n.mods > 1) {
			if (verbose==TRUE) {
				cat("Testing directed difference models: Interaction, additive main effects, difference model :\n")
				cat("-------------------------------------------------------------------------\n")
			}
			a1 <- cbind(aL1$ANOVA, ldply(aL1$models, function(X) {
				F <- fitmeasures(X)
				R <- inspect(X, "r2")
				names(R) <- "R2"
				n <- nobs(X)
				k <- free.max - F["df"]
				
				R2.p <- ifelse(k==0,
					NA,
					pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
				names(R2.p) <- "R2.p"
				return(c(F[c("cfi", "tli", "rmsea", "srmr")], R, R2.p))

			}))
			a1 <- a1[, !grepl(".id", colnames(a1))]
			a1$k <- free.max - a1$Df
			a1$R2.adj <- 1 - ((1-a1$R2))*((nobs(full)-1)/(nobs(full)-a1$k-1))
			a1$delta.R2 <- c(NA, a1$R2[1:(nrow(a1)-1)] - a1$R2[2:(nrow(a1))])
			if (verbose==TRUE) print(round(a1, 3))
			res <- a1
		}
	
		aL2 <- anovaList(list(cubic=cubic, full=full, SRSD=SRSD, SSD=SSD, sqdiff=sqdiff, null=null))
		if (aL2$n.mods > 1) {
			if (verbose==TRUE) {
				cat("\n\nTesting 'flat ridge' discrepancy models against full polynomial model:\n")
				cat("-------------------------------------------------------------------------\n")
			}
			a2 <- cbind(aL2$ANOVA, ldply(aL2$models, function(X) {
				F <- fitmeasures(X)
				R <- inspect(X, "r2")
				names(R) <- "R2"
				n <- nobs(X)
				k <- free.max - F["df"]
				R2.p <- ifelse(k==0,
					NA,
					pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
				names(R2.p) <- "R2.p"
				return(c(F[c("cfi", "tli", "rmsea", "srmr")], R, R2.p))
			}))
			a2 <- a2[, !grepl(".id", colnames(a2))]
			a2$k <- free.max - a2$Df
			a2$R2.adj <- 1 - ((1-a2$R2))*((nobs(full)-1)/(nobs(full)-a2$k-1))
			a2$delta.R2 <- c(NA, a2$R2[1:(nrow(a2)-1)] - a2$R2[2:(nrow(a2))])
			if (verbose==TRUE) print(round(a2, 3))
			res <- rbind(res, a2)
		}
		
		aL2b <- anovaList(list(cubic=cubic, full=full, RR=RR, sqdiff=sqdiff, null=null))
		if (aL2b$n.mods > 1) {
			if (verbose==TRUE) {
				cat("\n\nTesting 'rising ridge' against full polynomial model:\n")
				cat("-------------------------------------------------------------------------\n")
			}
			a2 <- cbind(aL2b$ANOVA, ldply(aL2b$models, function(X) {
				F <- fitmeasures(X)
				R <- inspect(X, "r2")
				names(R) <- "R2"
				n <- nobs(X)
				k <- free.max - F["df"]
				R2.p <- ifelse(k==0,
					NA,
					pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE))
				names(R2.p) <- "R2.p"
				return(c(F[c("cfi", "tli", "rmsea", "srmr")], R, R2.p))
			}))
			a2 <- a2[, !grepl(".id", colnames(a2))]
			a2$k <- free.max - a2$Df
			a2$R2.adj <- 1 - ((1-a2$R2))*((nobs(full)-1)/(nobs(full)-a2$k-1))
			a2$delta.R2 <- c(NA, a2$R2[1:(nrow(a2)-1)] - a2$R2[2:(nrow(a2))])
			if (verbose==TRUE) print(round(a2, 3))
			res <- rbind(res, a2)
		}
	}
	
	
	aL3 <- anovaList(list(absunc=absunc, absdiff=absdiff))
	if (aL3$n.mods > 1) {
		if (verbose==TRUE) {
			cat("\n\n-------------------------------------------------------------------------\n")	
			cat("Piecewise regression: absolute difference vs. unrestricted difference model\n")
			cat("-------------------------------------------------------------------------\n")
		}
		free.max2 <- getFreeParameters(absunc)
		a3 <- cbind(aL3$ANOVA, ldply(aL3$models, function(X) {
			F <- fitmeasures(X)
			R <- inspect(X, "r2")
			names(R) <- "R2"
			n <- nobs(X)
			k <- free.max2 - F["df"]
			R2.p <- pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE)
			names(R2.p) <- "R2.p"
			return(c(F[c("cfi", "tli", "rmsea", "srmr")], R, R2.p))
		}))
		a3 <- a3[, !grepl(".id", colnames(a3))]
		a3$k <- free.max2 - a3$Df
		a3$R2.adj <- 1 - ((1-a3$R2))*((nobs(absunc)-1)/(nobs(absunc)-a3$k-1))
		a3$delta.R2 <- c(NA, a3$R2[1:(nrow(a3)-1)] - a3$R2[2:(nrow(a3))])
		if (verbose==TRUE) print(round(a3, 3))
		res <- rbind(res, a3)
	}
	
	invisible(res)
	})
}




# compare CFI of intercept-only null model

# CFI2 <- function(x, m1="full", m0="null") {
# 	cfi2 <- 1-(inspect(x$models[[m1]], "fit")['chisq']-inspect(x$models[[m1]], "fit")['df']) / (inspect(x$models[[m0]], "fit")['chisq']-inspect(x$models[[m0]], "fit")['df'])
# 	names(cfi2) <- "CFI2"
# 	return(cfi2)
# }
# 
# TLI2 <- function(x, m1="full", m0="null") {
# 	tli2 <- ((inspect(x$models[[m0]], "fit")['chisq'] / inspect(x$models[[m0]], "fit")['df']) - (inspect(x$models[[m1]], "fit")['chisq'] / inspect(x$models[[m1]], "fit")['df'])) / ((inspect(x$models[[m0]], "fit")['chisq'] / inspect(x$models[[m0]], "fit")['df']) - 1)
# 	names(tli2) <- "TLI2"
# 	return(tli2)
# }