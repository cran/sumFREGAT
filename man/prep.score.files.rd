\name{prep.score.files}
\alias{prep.score.files}
\title{Prepare score files}
\description{
Calculates Z scores from P values and beta input
}
\usage{
prep.score.files(input.file, reference.file = "", output.file.prefix)
}

\arguments{
	\item{input.file}{a file with two mandatory columns (case-insensitive header):\cr\cr
	"ID": names of genetic variants (we suggest to provide rsIDs when possible)\cr
	"P": P value\cr
	Additional columns that can be present in input file:
	"CHROM": chromosome\cr
	"POS": positions for the same build as in \code{gene.file} (see \code{gene-based test functions}) and \code{reference.file} (37.3 with default files)\cr
	"EA": effect allele\cr
	"BETA": effect size (betas and genetic correlations should be calculated for the same genotype coding)\cr
	"EAF": effect allele frequency\cr
	"ANNO": functional annotations\cr\cr
	For example:\cr\cr
	CHROM POS ID EA P BETA EAF\cr
	1 196632134 1:196632134 T 0.80675 0.22946 0.00588\cr
	1 196632386 1:196632386 A 0.48694 0.65208 0.00588\cr
	1 196632470 1:196632470 G 0.25594 -0.19280 0.19412\cr\cr
	Avoid rounding of betas and P values as this can affect the precision of regional tests.\cr\cr
	The more data (columns) is present in input file, the more gene-based tests are available to run. Minimal input (rsIDs and P values)
	together with correlation matrices (reference matrices calculated from 1000G data are available at http://mga.bionet.nsc.ru/sumFREGAT/) allow to run \code{minp()},
	\code{simpleM}, and \code{sumchi()} tests. Adding info on effect allele ("EA") and effect size ("BETA") enables essentially all sumFREGAT tests.
	Adding allele frequencies enables standard weighting via beta distribution (see \code{gene-based test functions} for details).
	}

	\item{reference.file}{path to a reference file with additional data. Reference file from 1000G is available at http://mga.bionet.nsc.ru/sumFREGAT/.}

	\item{output.file.prefix}{if not set, the input file name will be used as output prefix.}

}
\value{
	
	does not return any value, writes output files with Z scores to be used in any type of
	gene-based analysis in sumFREGAT (see 'gene-based test functions').

}
\examples{

input.file <- system.file("testfiles/CFH.full.input.dat", package = "sumFREGAT")
prep.score.files(input.file, output.file.prefix = "CFH.scores.full")

\dontrun{

# requires reference file "ref1KG.MAC5.EUR_AF.RData" (can be downloaded
# at http://mga.bionet.nsc.ru/sumFREGAT/)

input.file <- system.file("testfiles/CFH.dat", package = "sumFREGAT")
prep.score.files(input.file, reference = "ref1KG.MAC5.EUR_AF.RData",
	output.file.prefix = "CFH.scores")

input.file <- system.file("testfiles/CFH.full.input.dat", package = "sumFREGAT")
prep.score.files(input.file, reference = "ref1KG.MAC5.EUR_AF.RData",
	output.file.prefix = "CFH.scores.full.ref")

}
}
