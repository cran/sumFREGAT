# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

check.method <- function(method) {
	if (method == 'Davies') method <- 'davies'
	if (method == 'Kuonen') method <- 'kuonen'
	if (method == 'Hybrid') method <- 'hybrid'

	method <- match.arg(method, c('davies', 'kuonen', 'hybrid'))
	return(method)
}

check.weights <- function(weights, beta.par) {
	if (is.null(weights)) {
#		fweights <- function(maf, beta.par) ifelse(maf > 0, dbeta(maf, beta.par[1], beta.par[2]), 0)
		fweights <- function(maf, a = as.numeric(beta.par[1]), b = as.numeric(beta.par[2])) ifelse(maf > 0, dbeta(maf, a, b), 0)
	} else {
		if (is.function(weights)) {
			fweights <- weights
		} else {
			stop("'weights.function' should be a function of MAF")
		}
	}
	return(fweights)
}

check.basis <- function(value, name, base = 'none', order) {
	if (base == 'bspline') {
		if (value < order) {
			value <- order
			warning(paste("bspline basis cannot be less than order, set to", value))
		}
	}
	if (value < 1) {
		value <- 1
		warning(paste(name, "cannot be less than 1, set to", value))
	}
	if (base == 'fourier' & (value + 1) %% 2 != 0) {
		if (ceiling(value) %% 2 != 0) { value <- ceiling(value)
		} else if (floor(value) %% 2 != 0) { value <- floor(value)
		} else value <- value - 1
		warning(paste("fourier basis should be an odd integer, set to", value))
	}
	if (value %% 1 != 0) {
		value <- round(value)
		warning(paste(name, "should be an integer number, set to", value))
	}
	return(value)
}

check.rho <- function(rho) {
	if (length(rho) == 1 & is.logical(rho) & rho) {
		#rho <- (0:10)/10
		rho <- c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2,0.5^2, 0.5, 1)
	} else {
		for (i in 1:length(rho)) {
			if (rho[i] < 0 || rho[i] > 1) {
				stop("rho should be >= 0 and <= 1")
			}
		}
	}
	return(rho)
}

check.spec.famSKAT <- function() {

	method <- check.method(method)
	if (length(rho) > 1 | (length(rho) == 1 & rho)) {
		rhos <- check.rho(rho)
		rho <- TRUE
		sapply(c('method', 'rhos', 'rho'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	} else {
		rho <- FALSE
		sapply(c('method', 'rho'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	}
}

check.spec.famFLM <- function() {

	BSF <- match.arg(BSF, c('fourier', 'bspline'))

	g <- TRUE

#	if (missing(GVF)) g <- FALSE
	if (is.null(GVF)) g <- FALSE
	if (is.logical(GVF)){
		if (!GVF) g <- FALSE else GVF <- 'bspline'
	}
	if (g) GVF <- match.arg(GVF, c('fourier', 'bspline'))

	kb0 <- kg0 <- order0 <- FALSE

	if (BSF == 'bspline') {
		order0 <- check.basis(order, 'order')
		kb0 <- check.basis(kb, 'basis', 'bspline', order0)
	} else { kb0 <- check.basis(kb, 'basis', 'fourier') }

	if (g) {
		if (GVF == 'bspline') {
			if (!order0) order0 <- check.basis(order, 'order')
			kg0 <- check.basis(kg, 'basis', 'bspline', order0)
		} else { kg0 <- check.basis(kg, 'basis', 'fourier') }
	}

	# base condition is kg >= kb

	if (g) {
		if (kg0 < kb0) {
			kb0 <- kg0
			warning (paste('kb cannot exceed kg, kb set to', kg0))
		}
		if (kg0 == kb0) {
			g <- FALSE
			if (GVF == 'bspline' & BSF == 'fourier') BSF <- 'bspline'
			if (GVF == 'fourier' & BSF == 'bspline') BSF <- 'fourier'
			warning ("GVF has no effect under kg = kb, the equivalent model with '", BSF, "' for BSF will be used")
		}
	}

	# if g then kg > kb
	# if kg == kb then !g

	if (BSF == 'bspline') betabasis0 = create.bspline.basis(norder = order0, nbasis = kb0)
	else betabasis0 <- create.fourier.basis(c(0, 1), nbasis = kb0)

	if (g) {
		if (GVF == 'bspline') { genobasis0 <- create.bspline.basis(norder = order0, nbasis = kg0)
		} else { genobasis0 <- create.fourier.basis(c(0, 1), nbasis = kg0) }
		J0 <- inprod(genobasis0, betabasis0)
		model0 <- paste(GVF, kg0, '-', BSF, kb0, sep = '')
	}
	else model0 <- paste(0, '-', BSF, kb0, sep = '')
	# all fourier bases are odd

	if (g) sapply(c('BSF', 'g', 'GVF', 'kb0', 'kg0', 'order0', 'model0', 'genobasis0', 'betabasis0', 'J0'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)
	else sapply(c('BSF', 'g', 'kb0', 'order0', 'model0', 'betabasis0'),
			function(x) get(x), simplify = FALSE, USE.NAMES = TRUE)

}
