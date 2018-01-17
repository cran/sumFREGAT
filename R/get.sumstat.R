# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

get.sumstat <- function() {
	corfile <- paste(cor.path, r, '.cor', sep = '')
	if (!file.exists(corfile)) {
		return(list(m0 = NULL, Z = 'no corfile', w = NULL, pos = NULL))
	}
	w <- pos <- NULL
	if (length(annoType) == 1) {
		Z <- read.vcf.info(scoreFile, geneFile, r, annoType)
		if (is.null(Z[[1]])) return(list(m0 = NULL, Z = NULL, w = NULL, pos = NULL))
		if (test == 'famFLM') pos <- read.vcf.pos(scoreFile, geneFile, r, annoType)
		ID <- Z$ID
		AF <- Z$EAF
		se.beta <- Z$SE.Beta
		Z <- Z$Z
	} else { 
		Z <- ID <- AF <- se.beta <- c()
		for (anno in annoType) {
			Z0 <- read.vcf.info(scoreFile, geneFile, r, anno)
			if (is.null(Z0[[2]])) next
			Z <- c(Z, Z0$Z)
			ID <- c(ID, Z0$ID)
			AF <- c(AF, Z0$EAF)
			se.beta <- c(se.beta, Z0$SE.Beta)
			if (test == 'famFLM') {
				pos0 <- read.vcf.pos(scoreFile, geneFile, r, anno)
				pos <- c(pos, pos0)
			}
		}
	}
	if (is.null(Z)) return(list(m0 = NULL, Z = NULL, w = NULL, pos = NULL))
	U <- read.table(corfile)
	v <- match(ID, rownames(U))
	m0 <- length(v)
	if (sum(!is.na(v)) == 0) return(list(m0 = m0, Z = 'no variants in corfile', w = NULL, pos = NULL))
	index <- !is.na(v) & se.beta != 'NaN'
	if (sum(index) == 0) return(list(m0 = m0, Z = 'no variants with nonzero scores', w = NULL, pos = NULL))
	AF <- as.numeric(AF[index])
	se.beta <- as.numeric(se.beta[index])
	Z <- as.numeric(Z[index])
	if (test == 'famFLM') pos <- as.numeric(pos[index])
	v <- v[index]
	U <- as.matrix(U[v, v])
	
	if (any(is.na(Z))) {
		if (all(is.na(Z))) return(list(m0 = m0, Z = NULL, w = w, pos = pos))
		v <- !is.na(Z)
		Z <- Z[v]
		if (length(Z) == 1) return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
		AF <- AF[v]
		U <- as.matrix(U[v, v])
		se.beta <- se.beta[v]
		if (test == 'famFLM') pos <- as.vector(pos[v])
	}
	
	if (any(is.na(AF))) {
		if (all(is.na(AF))) return(list(m0 = m0, Z = NULL, w = w, pos = pos))
		v <- !is.na(AF)
		Z <- Z[v]
		if (length(Z) == 1) return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
		AF <- AF[v]
		U <- as.matrix(U[v, v])
		se.beta <- se.beta[v]
		if (test == 'famFLM') pos <- as.vector(pos[v])
	}
	
	v <- sapply(1:length(AF), function(x) all(AF[x] >= 0, na.rm = T) & all(AF[x] <= 1, na.rm = T))
	if (sum(v) == 0) return(list(m0 = m0, Z = 'not in range', w = NULL, pos = NULL))
	Z <- Z[v]
	if (length(Z) == 1) return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
	AF <- AF[v]
	U <- as.matrix(U[v, v])
	se.beta <- se.beta[v]

	if (any(AF > .5) & test %in% c('famBT', 'famSKAT', 'MLR')) {
		v <- which(AF > .5)
		Z[v] <- -Z[v]
		AF[v] <- 1 - AF[v]
		for (i in v) { # reverse correlation for reverse coding
			U[i, ] <- -U[i, ]
			U[, i] <- -U[, i]
		}
	}

	if (!is.null(fweights)) w <- fweights(AF) / se.beta

	if (test == 'famFLM') {
		if (any(is.na(pos))) {
			pos <- NULL
		} else {
			if (Fan) {
				u <- Indep.Variants(U)
				Z <- as.vector(Z[u])
				if (length(Z) == 1) return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
				U <- U[u, u]
				w <- as.vector(w[u])
				pos <- as.vector(pos[u])
				AF <- AF[u]
			}
			v <- order(pos)
			Z <- as.vector(Z[v])
			w <- as.vector(w[v])
			pos <- pos[v]
			U <- as.matrix(U[v, v])
			AF <- as.vector(AF[v])
			if (sum(duplicated(pos)) > 0) {
				for (i in 1:length(pos)) {
					pos[duplicated(pos)] <- pos[duplicated(pos)] + 1
					if (sum(duplicated(pos)) == 0) break()
				}
			}
			if (max(pos) > min(pos)) pos <- (pos - min(pos)) / (max(pos) - min(pos))
			if (max(pos) == min(pos)) pos <- .5
			if (flip.genotypes == TRUE) {
				flip <- sumstat.flipper(U, Z, AF)
				U <- flip$U
				Z <- flip$Z
			}
		}
	}

	if (test == 'MLR') {
		u <- Indep.Variants(U)
		Z <- as.vector(Z[u])
		if (length(Z) == 1) return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
		U <- U[u, u]
	}

	U[U == 1] <- 0.999
	U[U == -1] <- -0.999
	diag(U) <- 1

	return(list(m0 = m0, Z = Z, U = U, w = w, pos = pos))
}

read.vcf.info <- function(data, geneFile, r, annoType) {
	if (getRversion() >= "3.2.0") {
		invisible(capture.output(invisible(capture.output(info <- seqminer::readVCFToListByGene(data, geneFile, r, annoType, 'ID', c('EAF', 'Z', 'SE.Beta'), ''), type = 'output')),  type = 'message'))
	} else { info <- seqminer::readVCFToListByGene(data, geneFile, r, annoType, 'ID', c('EAF', 'Z', 'SE.Beta'), '') }
	return(info)
}

read.vcf.pos <- function(genodata, geneFile, r, annoType) {
	if (getRversion() >= "3.2.0") {
		invisible(capture.output(invisible(capture.output(pos <- seqminer::readVCFToListByGene(genodata, geneFile, r, annoType, 'POS', '', '')$POS, type = 'output')),  type = 'message'))
	} else { pos <- seqminer::readVCFToListByGene(genodata, geneFile, r, annoType, 'POS', '', '')$POS }
	return(pos)
}

Indep.Variants <- function(U) {
	q <- qr(U, tol = 1.e-7)
	q$pivot[seq(q$rank)]
}

sumstat.flipper <- function(U, Z, AF) {
D1 <- sqrt(2. * AF*(1.-AF))
	for (j in 1:(length(Z)-1)) {
		A_B <- U[j,j+1]*D1[j]*D1[j+1] + (1.-2.*AF[j])*(1.-2.*AF[j+1])
		if (A_B < 0) {
			U[,j+1] <- -1.*U[,j+1]
			U[j+1,] <- -1.*U[j+1,]
			Z[j+1]  <- -1.*Z[j+1]
			AF[j+1] <-  1.-AF[j+1]
		}
	}
	return(list(U = U, Z = Z))#, AF = AF))
}
