\name{gene-based test functions}
\alias{SKAT}
\alias{SKATO}
\alias{sumchi}
\alias{BT}
\alias{PCA}
\alias{FLM}
\alias{simpleM}
\alias{minp}
\alias{ACAT}
\title{Gene-based tests on summary statistics}
\description{
A set of tests for gene-based association analysis on GWAS summary statistics:
sequence kernel association tests ('SKAT', 'SKATO'), sum of chi-squares ('sumchi'),
aggregated Cauchy association test ('ACAT'), burden test ('BT'), principal component
analysis-based test ('PCA'), functional linear model-based test ('FLM'),
Bonferroni correction test ('simpleM'), minimum P value test ('minp').
}
\usage{
SKAT(score.file, gene.file, genes = "all", cor.path = "cor/",
approximation = TRUE, anno.type = "", beta.par = c(1, 25),
weights.function = NULL, user.weights = FALSE, gen.var.weights = "se.beta",
method = "kuonen", acc = 1e-8, lim = 1e+6, rho = FALSE, p.threshold = 0.8,
write.file = FALSE, quiet = FALSE)

SKATO(score.file, gene.file, genes = "all", cor.path = "cor/",
approximation = TRUE, anno.type = "", beta.par = c(1, 25),
weights.function = NULL, user.weights = FALSE, method = "kuonen",
acc = 1e-8, lim = 1e+6, rho = TRUE, p.threshold = 0.8, write.file = FALSE,
quiet = FALSE)

sumchi(score.file, gene.file, genes = "all", cor.path = "cor/",
approximation = TRUE, anno.type = "", method = "kuonen", acc = 1e-8,
lim = 1e+6, write.file = FALSE, quiet = FALSE)

ACAT(score.file, gene.file, genes = "all", anno.type = "",
beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE,
gen.var.weights = "none", mac.threshold = NA, n = NA,
write.file = FALSE, quiet = FALSE)

BT(score.file, gene.file, genes = "all", cor.path = "cor/",
anno.type = "", beta.par = c(1, 25), weights.function = NULL,
user.weights = FALSE, write.file = FALSE, quiet = FALSE)

PCA(score.file, gene.file, genes = "all", cor.path = "cor/",
approximation = TRUE, anno.type = "", n, beta.par = c(1, 1),
weights.function = NULL, user.weights = FALSE, reference.matrix = TRUE,
fun = "LH", var.fraction = 0.85, write.file = FALSE, quiet = FALSE)

FLM(score.file, gene.file, genes = "all", cor.path = "cor/",
approximation = TRUE, anno.type = "", n, beta.par = c(1, 1),
weights.function = NULL, user.weights = FALSE, basis.function = "fourier",
k = 25, order = 4, flip.genotypes = FALSE, Fan = TRUE,
reference.matrix = TRUE, fun = "LH", write.file = FALSE, quiet = FALSE)

simpleM(score.file, gene.file, genes = "all", cor.path = "cor/",
anno.type = "", var.fraction = .995, write.file = FALSE, quiet = FALSE)

minp(score.file, gene.file, genes = "all", cor.path = "cor/",
anno.type = "", write.file = FALSE, quiet = FALSE)
}

\arguments{
	\item{score.file}{name of data file generated by \code{prep.score.files()}.}

	\item{gene.file}{can be a name of a custom text file listing genes in refFlat format. Values "hg19" and "hg38"
	can be set to use default gene files containing protein-coding genes with start and stop positions from corresponding build.
	Mind that the same build should be used in \code{score.file} and \code{gene.file}.}

	\item{genes}{character vector of gene names to be analyzed. Can be "chr1", "chr2" etc. to analyze
	all genes on a corresponding chromosome. If not set, function will
	attempt to analyze all genes listed in \code{gene.file}.}

	\item{cor.path}{path to a folder with correlation matrix files (one file per each gene to be analyzed).
	Correlation matrices in text format are allowed, though ".RData" is preferable as computationally efficient.
	Each file should contain a square matrix with correlation coefficients (r) between genetic variants
	of a gene. An example of correlation file format:\cr
	"snpname1" "snpname2" "snpname3" ...\cr
	"snpname1" 1 0.018 -0.003 ...\cr
	"snpname2" 0.018 1 0.081 ...\cr
	"snpname3" -0.003 0.081 1 ...\cr
	...\cr
	If genotypes are available, matrices can be generated as follows:\cr
	\code{cor.matrix <- cor(g)}\cr
	\code{save(cor.matrix, file = paste0(geneName, ".RData"))}\cr
	where \code{g} is a genotype matrix (nsample x nvariants) for a given gene with genotypes coded as 0, 1, 2
	(coding should be exactly the same that was used to generate GWAS statistics).\cr
	If genotypes are not available, correlations can be approximated through reference samples.
	Reference matrices from 1KG European sample can be downloaded at http://mga.bionet.nsc.ru/sumFREGAT/.\cr
	Names of correlation files should be constructed as "geneName.RData" (e.g. "ABCG1.RData", "ADAMTS1.RData", etc.)
	for ".RData" format or "geneName.txt" for text format.\cr
	Example corfiles can be found as:\cr
	system.file("testfiles/CFH.cor", package = "sumFREGAT")\cr
	system.file("testfiles/CFH.RData", package = "sumFREGAT")\cr
	}

	\item{approximation}{a logical value indicating whether approximation for large genes (>= 500 SNPs) should be used.
	Applicable for SKAT, SKATO, sumchi, PCA, FLM (default = TRUE for these methods).}

	\item{anno.type}{given (functional) annotations provided by user (see) \code{prep.score.files()},
	a character (or character vector) indicating annotation types to be used.}

	\item{n}{size of the sample on which summary statistics were obtained.}

	\item{beta.par}{two positive numeric shape parameters in the beta distribution to assign weights 
	for each genetic variant as a function of minor allele frequency (MAF) in the default weights function (see Details).}

	\item{weights.function}{a function of MAF to assign weights
	for each genetic variant. By default, the weights will be calculated using the beta distribution (see Details).}

	\item{user.weights}{a logical value indicating whether weights from the input file should be applied. Default = FALSE.}

	\item{gen.var.weights}{a character indicating whether scores should be weighted by the variance of genotypes:
	"none" - no weights applied, resulting in a sum chi-square test.
	"se.beta" - scores weighted by variance of genotypes estimated from P values and effect sizes (regression coefficients)
	if provided by user.
	"af" - scores weighted by variance of genotypes calculated as AF * (1 - AF), where AF is allele frequency. 
	In \code{SKAT()}, "se.beta" weighting results in a test analogous to a standard SKAT, while "af" is a way to approximate
	the standard test in case when betas are not available. Unweighted sum chi-square test gives more weights to common variants
	compared with SKAT-like tests.
	In \code{ACAT()}, the simplest default unweighted test combines P values (Z scores) without any additional info.
	Scores can be weighted by via gen.var.weights = 'se.beta' or 'af', similarly to SKAT.
	}

	\item{mac.threshold}{ an integer number. In ACAT, scores with MAC <= 10 will be combined using Burden test. MACs are calculated
	from MAFs, \code{n} must be set to apply \code{mac.threshold}.}

	\item{method}{the method for computing P value in kernel-based tests. Available methods are "kuonen", "davies" and "hybrid"
	(see Details). Default = "kuonen".}

	\item{acc}{accuracy parameter for "davies" method.}

	\item{lim}{limit parameter for "davies" method.}

	\item{rho}{if TRUE the optimal test (SKAT-O) is performed [Lee et al., 2012]. \code{rho} can be a vector of grid
	values from 0 to 1. The default grid is c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1).}

	\item{p.threshold}{the largest P value that will be considered as important when performing computational optimization in SKAT-O.
	All P values larger than p.threshold will be processed via burden test.}

	\item{basis.function}{a basis function type for beta-smooth. Can be set to "bspline" (B-spline basis) or
	"fourier" (Fourier basis, default).}

	\item{k}{the number of basis functions to be used for beta-smooth (default = 25).}

	\item{order}{a polynomial order to be used in "bspline". Default = 4 corresponds to the cubic B-splines.
	as no effect if only Fourier bases are used.}

	\item{flip.genotypes}{a logical value indicating whether the genotypes of some genetic variants should be
	flipped (relabelled) for their better functional representation [Vsevolozhskaya et al., 2014]. Default = FALSE.}

	\item{Fan}{if TRUE (default) then linearly dependent genetic variants will be omitted, as it
	was done in the original realization of FLM test by Fan et al. (2013).}

	\item{reference.matrix}{logical indicating whether the correlation matrices were generated using reference matrix.
	If TRUE, regularization algorithms will be applied in order to ensure invertibility and numerical stability of the matrices.
	Use 'reference.matrix = FALSE' ONLY IF you are sure that correlation matrices were generated using the same genotype data
	as for GWAS summary statistics in input.}

	\item{fun}{one of two regularization algorithms, 'LH' (default) or 'derivLH'. Currently both give similar results.}

	\item{var.fraction}{minimal proportion of genetic variance within region that should be explained by principal
	components used (see Details for more info).}

	\item{write.file}{output file name. If specified, output (as it proceeds) will be written 
	to the file.}

	\item{quiet}{quiet = TRUE suppresses excessive output from reading vcf file genewise.}

}
\details{
	'SKAT' uses the linear weighted kernel function to set the inter-individual
	similarity matrix \eqn{K = GWWG^T} for SKAT and \eqn{K = GW(I\rho + (1 - \rho)ee^T)WG^T} for SKAT-O,
	where \eqn{\mathit{G}} is the \eqn{\mathit{n\times p}}
	genotype matrix for \eqn{\mathit{n}} individuals and \eqn{\mathit{p}} genetic variants in the region, 
	\eqn{\mathit{W}} is the \eqn{\mathit{p\times p}} diagonal weight matrix, \eqn{I} is the \eqn{\mathit{p\times p}}
	identity matrix, \eqn{\rho} is pairwise correlation
	coefficient between genetic effects (which will be adaptively selected from given \code{rho}),
	and \eqn{e} is the \eqn{\mathit{p\times 1}} vector of ones. Given the shape parameters
	of the beta function, \code{beta.par = c(a, b)}, 
	the weights are defined using probability density function of the beta distribution:\cr
	\cr
	\eqn{W_{i}=(B(a,b))^{^{-1}}MAF_{i}^{a-1}(1-MAF_{i})^{b-1} },\cr
	\cr
	where \eqn{MAF_{i}} is a minor allelic frequency for the \eqn{i^{th}} genetic variant in the region,
	which is estimated from genotypes, and \eqn{B(a,b)} is the beta function. This way of defining weights
	is the same as in original SKAT (see [Wu et al., 2011] for details). \code{beta.par = c(1, 1)} corresponds
	to the unweighted SKAT. The same weighting principle can be applied in 'BT' (burden test, default weights \code{beta.par = c(1, 25)}),
	as well as in 'PCA' and 'FLM' tests (unweighted by default).\cr
	
	Depending on the method option chosen, either Kuonen or Davies method is used to calculate P values from the
	score statistics in SKAT. Both an Applied Statistics algorithm that inverts the characteristic function
	of the mixture chisq [Davies, 1980] and a saddlepoint approximation [Kuonen, 1999] are nearly exact,
	with the latter usually being a bit faster.\cr
	
	A hybrid approach has been proposed by Wu et al. [2016]. It uses the Davies' method with high accuracy,
	and then switches to the saddlepoint approximation method when the Davies' method fails to converge.
	This approach yields more accurate results in terms of type I errors, especially for small significance levels.
	However, 'hybrid' method runs several times slower than the saddlepoint approximation method itself (i.e. 'kuonen' method).
	We therefore recommend using the hybrid approach only for those regions that show significant (or nearly significant)
	P values to ensure their accuracy.\cr

	'ACAT' is an adaptation of a set-based aggregated Cauchy association test (ACAT-V). It has been recently proposed
	by Liu, Y. et al. (2019). This in the only gene-based test statistic that can be calculated from P values (Z scores) only,
	and does not require SNP-SNP correlation info.

	Burden test ('BT', collapsing technique) suggests that the effects of causal genetic variants within a region have
	the same direction and the majority of variants are causal. If this is not the case, other regional tests (SKAT and FLM)
	are shown to have higher power compared to burden test [Svishcheva et al., 2015]. By default, 'BT' assigns weights
	calculated using the beta distribution with shape parameters \code{beta.par = c(1, 25)}.\cr

	'PCA' test is based on the spectral decomposition of
	correlation matrix among genetic variants. The number of top principal components will be chosen
	in such a way that >= \code{var.fraction} of region variance can be explained by these PCs.
	By default, \code{var.fraction} = 0.85, i.e. PCs explain >= 85\% of region variance.\cr
	
	A similar principle is used in 'simpleM' to calculate the effective number of independent tests.\cr
	
	'FLM' test assumes that the effects of multiple genetic variants
	can be described as a continuous function, which can be modelled through B-spline
	or Fourier basis functions. When the number of basis functions
	(set by \eqn{k}) is less than the number of variants
	within the region, the FLM test has an advantage of using
	less degrees of freedom [Svishcheva, et al., 2015].\cr

	For genes with \eqn{m <= k}, functional linear models are equivalent to a
	standard multiple linear regression, and the latter is used for these cases.
	The ultimate model name is returned in output (the "model" column).\cr
	

}
\value{
	A data frame with map info, P values, numbers of variants
	and filtered variants for each of analyzed genes.
	\cr
	In addition:
	\cr
	- BT() returns gene-level estimates of effect sizes (betas) and their standard errors.
	\cr
	- FLM() returns the names of the functional models used for each region.
	Names shortly describe the functional basis and the number of basis functions used.
	E.g., "F25" means 25 Fourier basis functions, "B15" means 15 B-spline basis functions.
	\cr
	- PCA() returns the number of the principal components
	used for each region and the proportion of genetic variance
	they make up.

}
\references{
	Davies R.B. (1980) Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables, Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 29, N 3, P. 323-333.\cr
	Kuonen D. (1999) Saddlepoint Approximations for Distributions of Quadratic Forms in Normal Variables. Biometrika, Vol. 86, No. 4, P. 929-935.\cr
	Wu M.C., et al. (2011) Rare-variant association testing for sequencing data with the sequence kernel association test. Am. J. Hum. Genet., Vol. 89, P. 82-93.\cr
	Lee S., et al. (2012) Optimal unified approach for rare variant association testing with application to small sample case-control whole-exome sequencing studies. American Journal of Human Genetics, 91, 224-237.\cr
	Liu Y. et al. (2019) ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. Am. J. Hum. Genet. 104, 410-421.
	Svishcheva G.R., Belonogova N.M. and Axenovich T.I. (2015) Region-based association test for familial data under functional linear models. PLoS ONE 10(6): e0128999.\cr
	Wu B., et al. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association testing of variant set. Ann Hum Genet, 80(2): 123-135.\cr
	Vsevolozhskaya O.A., et al. (2014) Functional Analysis of Variance for Association Studies. PLoS ONE 9(9): e105074.\cr
	Fan R, Wang Y, Mills JL, Wilson AF, Bailey-Wilson JE, et al. (2013) Functional linear models for association analysis of quantitative traits. Genet Epidemiol 37: 726-42.
}
\examples{

# Using example score files "CFH.vcf.gz" and "CFH.full.vcf.gz"
# generated by prep.score.files() function (see examples for prep.score.files())

# Test available with minimal input (P values only):

score.file <- system.file("testfiles/CFH.vcf.gz",
	package = "sumFREGAT")
ACAT(score.file, gene.file = "hg19", genes = "CFH")

# Tests that require P values, rsID, and SNP-SNP correlation info:

score.file <- system.file("testfiles/CFH.vcf.gz",
	package = "sumFREGAT")
cor.path <- system.file("testfiles/", package = "sumFREGAT")

sumchi(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)
simpleM(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)
minp(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)

# Tests available with full input including P values, rsID,
# betas, effect allele, allele frequencies (for default weighting),
# and correlation matrices:

score.file <- system.file("testfiles/CFH.full.vcf.gz",
	package = "sumFREGAT")
cor.path <- system.file("testfiles/", package = "sumFREGAT")

BT(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)
SKAT(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)
SKATO(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path)

# Tests that require sample size to be provided as "n" argument:

PCA(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path, n = 85)
FLM(score.file, gene.file = "hg19", genes = "CFH", cor.path = cor.path, n = 85)

}
