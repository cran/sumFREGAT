# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'do.sumstat' <- function (scoreFile, geneFile, regions = 'all', cor.path,
annoType = '', n = NULL, beta.par = NULL, weights.function = NULL,
method = 'kuonen', acc = 1e-8, lim = 1e+6, GVF = FALSE, BSF = 'fourier', kg = 30,
kb = 25, order = 4, flip.genotypes = FALSE, Fan = TRUE, rho = FALSE,
var.fraction = .85, write.file = FALSE, test) {

############ COMMON CHECKS

if (missing(scoreFile)) stop("'scoreFile' is missing, with no default")
if (!file.exists(scoreFile)) {
	scoreFile1 <- paste(scoreFile, 'anno.vcf.gz', sep = '.')
	if (file.exists(scoreFile1)) {
		scoreFile <- scoreFile1
	} else { stop(paste(scoreFile, '- No such file or directory')) }
}
if (missing(cor.path)) {
	cor.path <- ''
} else {
	if (substr(cor.path, nchar(cor.path), nchar(cor.path)) != '/') cor.path <- paste(cor.path, '/', sep = '')
}

if (missing(geneFile)) geneFile <- system.file("testfiles/refFlat_hg19_6col.txt.gz", package = "sumFREGAT")
if (missing(regions)) regions <- read.table(geneFile)$V1

regions <- unique(regions)
nreg <- length(regions)

############ SPECIFIC CHECKS

fweights <- NULL # for get.sumstat()

if (test != 'MLR') {
	if (missing(beta.par)) { # ne komilfo
		if (test %in% c('famFLM', 'PCA')) beta.par <- c(1, 1) else beta.par <- c(1, 25)
	}
	fweights <- check.weights(weights.function, beta.par)
#	if (!is.null(fweights)) weights <- NULL
}
if (test %in% c('MLR', 'famFLM', 'PCA')) {
	if (is.null(n)) stop("'n' should be set for this type of test")
	n <- n - 1
}
if (test %in% c('famSKAT', 'famFLM')) {
	fun <- as.function(get(paste('check.spec', test, sep = '.')))
	environment(fun) <- environment()
	tmp <- fun()
	for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])
}

############ ANALYSIS

lgt <- 4 # n.columns in output
if (test == 'famFLM') {
	lgt <- 5
} else if (test %in% c('PCA', 'famBT')) lgt <- 6

environment(analyze.sumstat) <- environment()
environment(get.sumstat) <- environment()
sumstat.region <- as.function(get(paste('sumstat', test, sep = '.')))
environment(sumstat.region) <- environment()
if (test == 'famSKAT' & rho) environment(sumstat.famSKATO) <- environment()
if (test == 'famFLM') environment(sumstat.MLR) <- environment()

out <- as.data.frame(matrix(NA, nrow = nreg, ncol = lgt))
c <- 0
#t1 <- c()

for (i in 1:nreg) {
	#t0 <- proc.time()
	r <- as.character(regions[i])
	out[i, ] <- analyze.sumstat()
	if (write.file != FALSE) {
		if (c == 0) {
		cn <- c('region', 'markers', 'filtered.markers', 'pvalue')
		if (test == 'famFLM') {
			cn <- c(cn, 'model')
		} else if (test == 'famBT') { cn <- c(cn, 'beta', 'se.beta')
		} else if (test == 'PCA') cn <- c(cn, 'ncomponents', 'explained.variance.fraction')
		write.table(t(cn), file = write.file, col.names = F, row.names = F, quote = F)
		c <- 1
		}
		if (test == 'famFLM') {
			out[i, 5] <- gsub("fourier", "F", out[i, 5])
			out[i, 5] <- gsub("bspline", "B", out[i, 5])
		}
		write.table(out[i, ], file = write.file, append = T, col.names = F, row.names = F, quote = F)
	}
	#t1 <- rbind(t1, (proc.time() - t0))
}

#### FINAL MAKE UP

out[, 2:4] <- sapply(2:4, function(x) as.numeric(as.character(out[, x])))
if (test %in% c('PCA', 'famBT')) out[, 5:6] <- sapply(5:6, function(x) as.numeric(as.character(out[, x])))

colnames(out) [1:4] <- c('region', 'markers', 'filtered.markers', 'pvalue')
if (test == 'famFLM') {
	colnames(out) [5] <- 'model'
	out[, 5] <- gsub("fourier", "F", out[, 5])
	out[, 5] <- gsub("bspline", "B", out[, 5])
	out[, 2:4] <- sapply(2:4, function(x) as.numeric(as.character(out[, x])))
} else if (test == 'famBT') {
	colnames(out) [5:6] <- c('beta', 'se.beta')
} else if (test == 'PCA') colnames(out) [5:6] <- c('ncomponents', 'explained.variance.fraction')

############ FINISH

out
#cbind(out, t1)

}