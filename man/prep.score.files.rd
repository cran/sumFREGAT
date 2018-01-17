\name{prep.score.files}
\alias{prep.score.files}
\title{Prepare score files}
\description{
Calculates Z scores from P values and beta input
}
\usage{
prep.score.files(input.file, output.file.prefix)
}

\arguments{
	\item{input.file}{a file with the following columns:\cr\cr
	"CHROM": chromosome\cr
	"POS": positions\cr
	"ID": names of genetic variants, same as in files with genetic correlations\cr
	"REF": reference allele\cr
	"ALT": alternative allele\cr
	"P": p value\cr
	"BETA": effect size (betas and genetic correlations should be calculated for the same genotype coding)\cr
	"EAF": effect allele freequency\cr\cr
	For example:\cr\cr
	CHROM POS ID REF ALT pvalue beta EAF\cr
	1 196632134 1:196632134 C T 0.80675 0.22946 0.00588\cr
	1 196632386 1:196632386 G A 0.48694 0.65208 0.00588\cr
	1 196632470 1:196632470 A G 0.25594 -0.19280 0.19412\cr
	Avoid rounding of betas and pvalues as this can affect the precision of regional tests.\cr
	}

	\item{output.file.prefix}{if not set, the input file name will be used as output prefix.}

}
\value{
	
	does not return any value, writes output files with Z scores to be used in any type of
	gene-based analysis: \code{SKAT()}, \code{BT()}, \code{MLR()}, \code{FLM()}, \code{PCA()}.

}
\examples{

input.file <- system.file("testfiles/CFH.pvalues.dat",
	package = "sumFREGAT")
prep.score.files(input.file, "CFH.scores")

}
