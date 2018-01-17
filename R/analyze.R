# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

analyze.sumstat <- function() {

	Z <- get.sumstat()
	m0 <- Z$m0
	if (is.null(Z$Z)) {
		warning("No variants to analyze in region ", r, ", skipped")
		return(c(r, 0, 0, rep(NA, lgt - 3)))
	}
	if (length(Z$Z) == 1) {
		if (Z$Z == 'not in range') {
			warning("MAFs are not in range 0..1 in region ", r, ", skipped")
			return(c(r, m0, 1, rep(NA, lgt - 3)))
		}
		if (Z$Z == 'no corfile') {
		warning("No correlation file for gene ", r, ", skipped")
		return(c(r, rep(NA, lgt - 1)))
		}
		if (Z$Z == 'no variants in corfile') {
		warning("No variants in correlation file for gene ", r, ", skipped")
		return(c(r, 0, 0, rep(NA, lgt - 3)))
		}
		if (Z$Z == 'no variants with nonzero scores') {
		warning("No variants with nonzero scores for gene ", r, ", skipped")
		return(c(r, m0, 0, rep(NA, lgt - 3)))
		}
	}
	m1 <- length(Z$Z)
	if (m1 == 1) {
		if (test == 'famBT') return(c(r, m0, m1, rep(NA, lgt - 3)))
		return(c(r, m0, m1, pnorm(abs(Z$Z), lower.tail = F) * 2, rep(NA, lgt - 4)))
	}
	if (test == 'famFLM' & is.null(Z$pos)) {
		warning("Some positions are not assigned in region ", r, ", skipped")
		return(c(r, m0, m1, rep(NA, lgt - 3)))
	}
	if (test == 'famSKAT' & rho) {
		return(c(r, m0, m1, sumstat.famSKATO(Z)))
	}
	return(c(r, m0, m1, sumstat.region(Z)))

}

