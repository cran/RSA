#' @title Performs several RSA model tests on a data set with two predictors
#'
#' @description
#' Performs several RSA model tests on a data set with two predictors
#'
#' @details
#' No details so far.
#'
#' @export
#' @param formula A formula in the form \code{z ~ x*y}, specifying the variable names used from the data frame, where z is the name of the response variable, and x and y are the names of the predictor variables.
#' @param data A data frame with the variables
#' @param sample.cor A matrix with the sample correlation matrix - you have to provide either raw data in the data argument, or a sample.cor and a sample.nobs
#' @param sample.nobs Number of observations in the sample (goes along with sample.cor)
#' @param center Should predictor variables be centered on the sample mean before analyses?
#' @param scale Should predictor variables be scaled to SD = 1 before analyses?
#' @param na.rm Remove missings before proceeding?
#' @param out.rm Should outliers according to Bollen & Jackman (1980) criteria be excluded from analyses?
#' @param breakline Should the breakline in the unconstrained absolute difference model be allowed (the breakline is possible from the model formulation, but empirically rather unrealistic ...)
#' @param verbose Should additional information during the computation process be printed?
#' @param models A vector with names of all models that should be computed. Should be any from c("absdiff", "absunc", "diff", "additive", "IA", "sqdiff", "SSD", "SRSD", "full"). For \code{models="all"}, all models are computed, for \code{models="default"} all models but absolute difference models are computed.
#' @param ... Additional parameters passed to the lavaan sem function. For example: \code{se="boot"}
#'
#'
#' @seealso \code{\link{demoRSA}}, \code{\link{plotRSA}}, \code{\link{RSA.ST}}
#'
#' @examples
#' # Compute response surface from a fake data set
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 15
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	sqdiff <- (x-y)^2
#' 	z.diff <- diff + rnorm(n, 0, err)
#' 	z.abs <- absdiff + rnorm(n, 0, err)
#' 	z.sq <- sqdiff + rnorm(n, 0, err)
#' 	z.add <- diff + 0.4*x + rnorm(n, 0, err)
#' 	z.complex <- 0.4*x + - 0.2*x*y + + 0.1*x^2 - 0.03*y^2 + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df)
#' print(r1)
#' compare(r1)
#' plot(r1)
#' plot(r1, model="SRSD")
#' plot(r1, model="full", type="c")
#' getPar(r1, "coef")	# print model parameters including SE and CI
#' RSA.ST(r1)	# get surface parameters
#'
#' # Motive congruency example
#' data(motcon)
#' r.m <- RSA(negAct~EM*IM, motcon)


# formula <- z.sq~x*y; data <- df; center=FALSE; scale=FALSE; out.rm=TRUE; breakline=FALSE; verbose=TRUE

RSA <- function(formula, data=NULL, sample.cor=NULL, sample.nobs=NULL, center=FALSE, scale=FALSE, na.rm=FALSE, out.rm=TRUE, breakline=FALSE, models="default", verbose=TRUE, ...) {

	if (length(models)==1 & models[1]=="all") {models <- c("absdiff", "absunc", "diff", "additive", "IA", "sqdiff", "SSD", "SRSD", "full")}
	if (length(models)==1 & models[1]=="default") {models <- c("diff", "additive", "IA", "sqdiff", "SSD", "SRSD", "full")}
	if (any(!models %in% c("absdiff", "absunc", "diff", "additive", "IA", "sqdiff", "SSD", "SSD", "SRSD", "SRSD", "full"))) {
		stop("Unknown model name provided in parameter 'models'.")
	}
	
	if (is.null(data) & is.null(sample.cor) & is.null(sample.nobs)) stop("Please provide either a data frame or a correlation matrix!")
		
	if (!is.null(data)) {
		mode <- "data"
	} else {
		mode <- "cor"
		if (any (is.null(c(sample.cor, sample.nobs)))) stop("You have to provide both sample.cor and sample.nobs")
	}
	
	# set all result objects to NULL as default
	s.full <- s.IA <- s.diff <- s.absdiff <- s.additive <- s.sqdiff <- s.sq.shift <- s.sq.rot <- s.absunc <- NULL
	sq.shape <- ""
	
	DV <- all.vars(formula)[1]
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]

	## Step 0a: Standardize values and calculate higher order terms
	if (mode=="data") {
		df <- data
		df[, IV1] <- scale(df[, IV1], center=center, scale=scale)
		df[, IV2] <- scale(df[, IV2], center=center, scale=scale)
		
		df <- add.variables(formula, data.frame(data.matrix(df)))
	
		# omit warnings if the zero point is outside of data range
		if (0 < min(df[, IV1], na.rm=TRUE) | 0 > max(df[, IV1], na.rm=TRUE)) {warning(paste("The numerical zero point is outside of the range of variable", IV1, ". Please consider a re-centering of the variable."))}
		if (0 < min(df[, IV2], na.rm=TRUE) | 0 > max(df[, IV2], na.rm=TRUE)) {warning(paste("The numerical zero point is outside of the range of variable", IV2, ". Please consider a re-centering of the variable."))}
	}
	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV_IA <- paste0(IV1, "_", IV2)
	W_IV1 <- paste0("W_", IV1)
	W_IV2 <- paste0("W_", IV2)


	## Run polynomial regression as a linear model

	if (mode=="data") {
		f <- paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA, IV22, sep=" + "))
		rs <- lm(f, df)
	} else {
		rs <- NA
	}
	
	
	if (out.rm == TRUE & mode=="data") {
		# get outliers and influential cases according to Bollen & Jackman, 1980
	
		inf <- influence.measures(rs)
		outs <- which(apply(inf$is.inf[, c("dffit", "cook.d", "hat")], 1, sum) == 3)
		if (verbose==TRUE) {
			print(paste("Removed", length(outs), "case(s) according to Bollen & Jackman (1980) criteria."))
		}
		if (length(outs)>0) {
			df <- df[-outs, ]
		}
	}
	

	## Test lower order models
	
	library(lavaan)
	
	#COV <- cov(df[, c(DV, IV1, IV2, IV_IA, IV12, IV22, "W", W_IV1, W_IV2)], use="complete.obs")
	
	poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22)
	
	if ("additive" %in% models) {
		if (verbose==TRUE) print("Computing additive model ...")
		m.additive <-  paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		sep="\n")
		if (mode=="data") {
			s.additive <- sem(m.additive, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.additive <- sem(m.additive, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		#summary(s.additive, fit.measures=TRUE)
	}

	if ("diff" %in% models) {
		if (verbose==TRUE) print("Computing difference model ...")
		m.diff <- paste(poly,
			"b3==0",
			"b4==0",
			"b5==0",
			"b1 == -b2",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			sep="\n")
			if (mode=="data") {
				s.diff <- sem(m.diff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
			} else {
				s.diff <- sem(m.diff, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
			}
		
		#summary(s.diff, fit.measures=TRUE)
	}

	if ("IA" %in% models) {
		if (verbose==TRUE) print("Computing interaction model ...")
		m.IA <- paste(poly,
			"b3==0",
			"b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
		sep="\n")
		if (mode=="data") {
			s.IA <- sem(m.IA, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.IA <- sem(m.IA, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.IA, fit.measures=TRUE)
	}
	
	if ("sqdiff" %in% models) {
		if (verbose==TRUE) print("Computing squared difference model ...")
		m.sqdiff <- paste(poly,
			"b1==0",
			"b2==0",
			"b3==b5",
			"b3+b4+b5==0",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			sep="\n")			
		if (mode=="data") {
			s.sqdiff <- sem(m.sqdiff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.sqdiff <- sem(m.sqdiff, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.sqdiff, fit.measures=TRUE)
	}
	
	if (any(models %in% c("SSD"))) {
		if (verbose==TRUE) print("Computing shifted squared difference model ...")
		m.sq.shift <- paste(poly,
			"b1==-b2",
			"b3==b5",
			"b3+b4+b5==0",
			"C := b1/(2*b3)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			
			sep="\n")
		if (mode=="data") {
			s.sq.shift <- sem(m.sq.shift, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.sq.shift <- sem(m.sq.shift, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.sq.shift, fit.measures=TRUE)
	}
	
	if (any(models %in% c("SRSD"))) {
		if (verbose==TRUE) print("Computing rotated squared difference model ...")
		m.sq.rot <- paste(paste(poly, " + start(0.05)*", IV22),
			"b5^2 > 0.00000001",
			"b1 == (b2*b4)/(2*b5)",
			"b3 == (b4*b4)/(4*b5)",
			"C := -.5*(b2/b5)",
			"S := -(b1/b2)",
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			sep="\n")
		if (mode=="data") {
			s.sq.rot <- sem(m.sq.rot, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.sq.rot <- sem(m.sq.rot, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
	
		if (coef(s.sq.rot)[3] >= 0) {
			sq.shape <- "up"
			} else {
				sq.shape <- "down"
			}
		#summary(s.sq.rot, fit.measures=TRUE)
	}
	
	
	if ("full" %in% models) {
		if (verbose==TRUE) print("Computing polynomial model ...")
		m.full <-  paste(poly,
			"a1 := b1+b2",
			"a2 := b3+b4+b5",
			"a3 := b1-b2",
			"a4 := b3-b4+b5",
			"X0 := (b2*b4 - 2*b1*b5) / (4*b3*b5 - b4^2)",
			"Y0 := (b1*b4 - 2*b2*b3) / (4*b3*b5 - b4^2)",
			"p11 := (b5 - b3 + sqrt(((b3 - b5)^2) + (b4^2))) / b4",
			"p10 := Y0 - p11*X0",
			"p21 :=  (b5 - b3 - sqrt((b3 - b5)^2 + b4^2)) / b4", 
			"p20 := Y0 - p21*X0",
			"as1X := b1 + p11*b2 + b4*p10 + 2*b5*p10*p11",
			"as2X := b3 + b4*p11 + (p11^2)*b5",
			"as1Y := b1/p11 + b2 - (2*b3*p10)/p11^2 - (b4*p10)/p11",
			"as2Y := b3/p11^2 + b4/p11 + b5",
			"as3X := b1 + p21*b2 + b4*p20 + 2*b5*p20*p21",
			"as4X := b3 + b4*p21 + (p21^2)*b5",
			"as3Y := b1/p21 + b2 - (2*b3*p20)/p21^2 - (b4*p20)/p21",
			"as4Y := b3/p21^2 + b4/p21 + b5",
			# eigenvalues
			"l1 := (b3 + b5 + sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2", 
			"l2 := (b3 + b5 - sqrt((b3+b5)^2 - 4*b3*b5 + b4^2))/2",
			sep="\n"
		)
		if (mode=="data") {
			s.full <- sem(m.full, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.full <- sem(m.full, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.full, fit.measures=TRUE)
	}
	
	#m.absdiff.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + 0*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	"b1 == -b2",
	#	"b7 == -b8",
	#	"b7 == -2*b1",
	#	sep="\n")
	#s.absdiff.JRE <-  sem(m.absdiff.JRE, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	#summary(s.absdiff.JRE, fit.measures=TRUE)
	
	# the unconstrained absolute difference model - Edwards (2002) formula
	#m.absunc.JRE <-  paste(
	#	paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22, " + b6*W.JRE + b7*W.JRE_", IV1, " + b8*W.JRE_", IV2),
	#	sep="\n")
	#s.absunc.JRE <-  sem(m.absunc.JRE, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
	#summary(s.absunc.JRE, fit.measures=TRUE)
	
	
	if ("absdiff" %in% models) {
		if (verbose==TRUE) print("Computing constrained absolute difference model ...")
		m.absdiff <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			"b1 == 0",
			"b2 == 0",
			"b6 == 0",
			"b7 == -b8",
			sep="\n")
		if (mode=="data") {
			s.absdiff <- sem(m.absdiff, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.absdiff <- sem(m.absdiff, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.absdiff.JRE, fit.measures=TRUE)
	}
	
	if ("absunc" %in% models) {
		# the unconstrained absolute difference model - new formula
		if (verbose==TRUE) print("Computing unconstrained absolute difference model ...")
		m.absunc <-  paste(
			paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b6*W + b7*W_", IV1, " + b8*W_", IV2),
			ifelse(breakline==TRUE, "b6==0", ""),
			sep="\n")
		if (mode=="data") {
			s.absunc <- sem(m.absunc, data=df, fixed.x=TRUE, meanstructure=TRUE, ...)
		} else {
			s.absunc <- sem(m.absunc, sample.cov=sample.cor, sample.nobs=sample.nobs, fixed.x=TRUE, meanstructure=TRUE, ...)
		}
		
		#summary(s.absunc.JRE, fit.measures=TRUE)
	}
	
	
	## Build results object
	res <- list(models = list(full=s.full, IA=s.IA, diff=s.diff, absdiff=s.absdiff, additive=s.additive, sqdiff=s.sqdiff, SSD=s.sq.shift, SRSD=s.sq.rot, absunc=s.absunc), sq.shape = sq.shape, LM=rs, formula=formula, data=df, DV=DV, IV1=IV1, IV2=IV2, IV12=IV12, IV22=IV22, IV_IA=IV_IA, W_IV1=W_IV1, W_IV2=W_IV2, r.squared = summary(rs)$r.squared)
	
	attr(res, "class") <- "RSA"
	return(res)
}



getfit <- function(x) {
	library(plyr)
	res <- laply(x$models, function(m) {fitMeasures(m)[c("aic", "bic", "cfi", "tli", "chisq", "df", "pvalue", "rmsea", "srmr")]})
	rownames(res) <- names(x$models)
	res <- data.frame(res)
	res$rel.aic <- (res$aic - min(res$aic))/(max(res$aic) - min(res$aic))
	res$rel.bic <- (res$bic - min(res$bic))/(max(res$bic) - min(res$bic))
	return(res)
}


bestmodel <- function(x) {
	
	F <- summary(x$LM)$fstatistic
	p.model <- 1-pf(F[1], F[2], F[3])
	if (p.model > .05) {
		return("Overall model is not significant.")
	}
	
	f <- getfit(x)
	
	# Choose family based on AIC and BIC
	
	# if only one model is selected by AIC: that's it!
	if (length(which(f$rel.aic < .01)) == 1) {
		m <- rownames(f)[which(f$rel.aic < .01)]
	} else 	
	# AIC is ambiguos: let BIC decide!
	if (length(which(f$rel.aic < .01)) > 1 & length(which(f$rel.bic < .01)) == 1) {
		m <- rownames(f)[which(f$rel.bic < .01)]
	} else 	
	# If AIC and BIC are ambiguos: let AIC decide
	if (length(which(f$rel.aic < .01)) > 1 & length(which(f$rel.bic < .01)) > 1) {
		m <- rownames(f)[which.min(f$rel.aic)]
	} else {
		warning("Could not determine best model! err1")
		return("err1")
	}
	
	
	## Let ALWAYS BIC decide
	#m <- rownames(f)[which.min(f$rel.bic)]
	
	if (m %in% c("full", "IA", "additive", "diff")) {
		# find the last model that is not significantly different from full model
		a1 <- with(x$models, {anova(diff, additive, IA, full)})
		M <- which(a1[, "Pr(>Chisq)"] < .05)
		if (length(M) > 0) {
			return(rownames(a1)[min(which(a1[, "Pr(>Chisq)"] < .05))-1])
			} else {
				return("diff")
			}
	} else 
	if (m %in% c("sqdiff")) {
		return("sqdiff")
	} else 
	if (m %in% c("absdiff", "absunc")) {
		a1 <- with(x$models, {anova(absdiff, absunc)})
		return(ifelse(a1[2, "Pr(>Chisq)"] < .05, "absunc", "absdiff"))
	} else {
		warning("Could not determine best model! err2")
		return("err2")
	}
}

