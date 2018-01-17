\name{BT}
\alias{BT}
\title{Family Burden Test}
\description{
Burden test on summary statistics
}
\usage{

BT(scoreFile, geneFile, regions, cor.path = "", annoType = "",
beta.par = c(1, 25), weights.function = ifelse(maf > 0,
dbeta(maf, beta.par[1], beta.par[2]), 0), write.file = FALSE)
}
\arguments{
	\item{scoreFile}{name of data file generated by \code{prep.score.files()}.}

	\item{geneFile}{name of a text file listing genes in refFlat format. If not set, hg19 file
	will be used (see Examples below).}

	\item{regions}{character vector of gene names to be analysed. If not set, function will
	attempt to analyse all genes listed in \code{geneFile}.}

	\item{cor.path}{path to a folder with correlation files (one file per each gene to be analysed).
	Names of correlation files should be constructed as "geneName.cor" (e.g. "ABCG1.cor", "ADAMTS1.cor", etc.)
	Each file should contain a square matrix with correlation coefficients (r) between genetic variants
	of a gene. An example of correlation file format:\cr
	"snpname1" "snpname2" "snpname3" ...\cr
	"snpname1" 1 0.018 -0.003 ...\cr
	"snpname2" 0.018 1 0.081 ...\cr
	"snpname3" -0.003 0.081 1 ...\cr
	...\cr
	One way to generate such file from original genotypes is:\cr
	\code{write.table(cor(g), file = paste0(geneName, ".cor"))}\cr
	where \code{g} is a genotype matrix (nsample x nvariants) for a given gene with genotypes coded as 0, 1, 2
	(exactly the same coding that was used to generate betas).
	}

	\item{annoType}{for files annotated with the \code{seqminer} package, a character (or character vector) indicating annotation types to be used (e.g. 
	"Nonsynonymous", "Start_Loss", "Stop_loss", "Essential_Splice_Site")}

	\item{beta.par}{two positive numeric shape parameters in the beta distribution to assign weights 
	for each genetic variant as a function of MAF in the default weights function (see Details). Default = c(1, 25).}

	\item{weights.function}{a function of minor allele frequency (MAF) to assign weights
	for each genetic variant. By default, the weights will be calculated using the beta distribution (see Details).}

	\item{write.file}{output file name. If specified, output (as it proceeds) will be written 
	to the file.}
}
\details{
	Burden test (collapsing technique) suggests that the effects of causal genetic variants within a region have the same direction. If this is not the case, other regional tests (SKAT and FLM) are shown to have higher power compared to burden test [Svishcheva, et al., 2015].\cr\cr
	By default, BT assigns weights calculated using the beta distribution. Given the shape parameters of the beta function, \code{beta.par = c(a, b)}, 
	the weights are defined using probability density function of the beta distribution:\cr
	\cr
	\eqn{W_{i}=(B(a,b))^{^{-1}}MAF_{i}^{a-1}(1-MAF_{i})^{b-1} },\cr
	\cr
	where \eqn{MAF_{i}} is a minor allelic frequency for the \eqn{i^{th}} genetic variant in region, which is estimated from genotypes, and \eqn{B(a,b)} is the beta function.\cr\cr
	\code{beta.par = c(1, 1)} corresponds to the unweighted burden test.
}
\value{
	A data frame containing P values, estimates of betas and their s.e., numbers of variants and filtered variants for each of analyzed regions.
}
\references{
	Svishcheva G.R., Belonogova N.M. and Axenovich T.I. (2015) Region-based association test for familial data under functional linear models. PLoS ONE 10(6): e0128999.\cr
	}
\examples{

## Run BT with example files:
VCFfileName <- system.file("testfiles/CFH.scores.anno.vcf.gz",
	package = "sumFREGAT")
cor.path <- system.file("testfiles/", package = "sumFREGAT")
out <- BT(VCFfileName, region = 'CFH', cor.path = cor.path)

}