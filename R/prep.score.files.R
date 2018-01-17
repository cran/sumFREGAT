# sumFREGAT (2017) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

prep.score.files <- function(input.file, output.file.prefix) {#, annotate = FALSE, path.to.reference, annotation.param) { # 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'P', 'BETA', 'EAF'
	df <- read.table(input.file, header = TRUE, as.is = TRUE)
	colnames(df) <- toupper(colnames(df))
	if (any(colnames(df) %in% c('CHR', 'CHROMOSOME'))) {
		colnames(df)[colnames(df) %in% c('CHR', 'CHROMOSOME')] <- 'CHROM'}
	if (any(colnames(df) %in% c('POSITION', 'POSITIONS', 'MAP'))) {
		colnames(df)[colnames(df) %in% c('POSITION', 'POSITIONS', 'MAP')] <- 'POS'}
	if (any(colnames(df) %in% c('PVALUE', 'PV', 'PVAL', 'P.VALUE', 'P_VALUE'))) {
		colnames(df)[colnames(df) %in% c('PVALUE', 'PV', 'PVAL', 'P.VALUE', 'P_VALUE')] <- 'P'}
	ColNames <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'P', 'BETA', 'EAF')
	v <- !ColNames %in% colnames(df)
	if (sum(v)) stop(paste("Some columns are missing:", paste(ColNames[v], collapse = ', ')))
	df$Z <- qnorm(df$P / 2, lower.tail = FALSE) * sign(df$BETA)
	df$SE.BETA <- df$BETA / df$Z
	#df$Z <- df$BETA / df$SE.BETA
	df <- df[order(df[, 'POS']), ]
	df <- df[order(df[, 'CHROM']), ]
	vcf <- df[, c('CHROM', 'POS', 'ID', 'REF', 'ALT')]
	colnames(vcf)[1] <- '#CHROM'
	vcf$POS <- format(vcf$POS, scientific = FALSE)
	vcf$POS <- gsub(' ', '', vcf$POS)
	vcf <- cbind(vcf, QUAL = '.', FILTER = '.')
	vcf$INFO <- paste('EAF=', df$EAF, ';Z=', df$Z, ';SE.Beta=', df$SE.BETA, sep = '')
	title <- c('##INFO=<ID=EAF,Number=1,Type=Float,Description="Effect allele frequency">', '##INFO=<ID=Z,Number=1,Type=Float,Description="Z statistics">', '##INFO=<ID=SE.Beta,Number=1,Type=Float,Description="SE Beta">')

	if (!missing(output.file.prefix)) {
		fn <- paste(output.file.prefix, 'vcf', sep = '.')
	} else {
		fn <- paste(input.file, 'vcf', sep = '.')
	}
	write.table(title, fn, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
	suppressWarnings(write.table(vcf, fn, row.names = FALSE, quote = FALSE, append = TRUE, sep = '\t'))

	fn.gz <- paste(fn, 'gz', sep = '.')
	if (file.exists(fn.gz)) system(paste('rm', fn.gz))
	system(paste('bgzip', fn))
	system(paste('tabix -p vcf', fn.gz))
	print(paste('File', fn.gz, 'has been created'))

}

# annotate <- function (input.file, output.file.prefix, path.to.reference, annotation.param) {
	# if (!missing(output.file.prefix)) {
		# fn.anno <- paste(output.file.prefix, 'anno.vcf', sep = '.')
	# } else {
		# fn.anno <- paste(input.file, 'anno.vcf', sep = '.')
	# }
	# if (missing(path.to.reference)) {
		# path.to.reference <- ''
	# } else {
		# if (substr(path.to.reference, nchar(path.to.reference), nchar(path.to.reference)) != '/') path.to.reference <- paste0(path.to.reference, '/')
	# }
	# if (requireNamespace("seqminer", quietly = TRUE)) {

		# if (missing(annotation.param)) annotation.param <- seqminer::makeAnnotationParameter(list(
			# reference = paste0(path.to.reference, "hs37d5.fa"),
			# geneFile = paste0(path.to.reference, "refFlat_hg19.txt.gz"),
			# codonFile = paste0(path.to.reference, "codon.txt"),
			# priorityFile = paste0(path.to.reference, "priority.txt")))
##		annotation.param -> param
##		seqminer::annotateVcf(fn, fn.anno, params = param)
		# annotateVcfMy(fn, fn.anno, annotation.param)
		# fn.anno.gz <- paste(fn.anno, 'gz', sep = '.')
		# if (file.exists(fn.anno.gz)) system(paste('rm', fn.anno.gz))
		# system(paste('bgzip', fn.anno))
		# system(paste('tabix -p vcf', fn.anno.gz))
		# print(paste('Annotated file', fn.anno.gz, 'has been created'))
		# } else { stop(paste("Please install 'seqminer' package to process VCF file", sep = '')) }
# }

# annotateVcfMy <- function (inVcf, outVcf, params)
##a function from seqminer, with changes
# {
   # params$inputFormat = "vcf"
##  param <- makeAnnotationParameter(param)
##  res <- validateAnnotationParameter(param)
   # params <- seqminer::makeAnnotationParameter(params)
   # res <- seqminer::validateAnnotationParameter(params)
   # if (!res[[1]]) {
       # cat(paste(res[[2]], collapse = "\\n"))
       # stop("Stop due to critical error")
   # }
   # seqminer:::verifyFilename(inVcf, outVcf)
   # storage.mode(inVcf) <- "character"
   # storage.mode(outVcf) <- "character"
   # .Call("anno", inVcf, outVcf, params, PACKAGE = 'seqminer')
   # invisible(NULL)
# }

# download.annotation.resource <- function(outputDirectory)
## a function from 'seqminer' package
# {
 # outDir = outputDirectory
 # prepare a writable dir
 # if (!seqminer::dir.exists(outDir)) {
   # message(gettextf("Create output directory: %s", outDir))
   # seqminer::dir.create(outDir, recursive = TRUE)
 # }
 # if (!seqminer:::isDirWritable(outDir)) {
   # stop(gettextf("Unable to write to directory: %s", outDir))
 # }

 # download function
 # download <- function(url) {
   # fn <- basename(url)
    # destfile <- file.path(outDir, fn)
    # if (file.exists(destfile)) {
      # warning(gettextf("Overwriting %s", fn))
    # }
    # utils::download.file(url, destfile)
  # }
  # download resources
  # message("Begin download TabAnno resource files (human hg19)...")


  # message("Download reference file and its index:")
  # download("https://qbrc.swmed.edu/zhanxw/software/anno/resources/hs37d5.fa")
  # download("https://qbrc.swmed.edu/zhanxw/software/anno/resources/hs37d5.fa.fai")

  # message("Download gene definition:")
  # download("https://qbrc.swmed.edu/zhanxw/software/anno/resources/refFlat_hg19.txt.gz")

  # message("Download TabAnno codon definition and annotation priority files:")
  # download("https://qbrc.swmed.edu/zhanxw/software/anno/codon.txt")
  # download("https://qbrc.swmed.edu/zhanxw/software/anno/priority.txt")

  # message("Download completed")
  # message(gettextf("You can begin to use it:"))
  # message(gettextf(" param <- makeAnnotationParameter(list(reference = \"%s\", geneFile = \"%s\", codonFile = \"%s\", priorityFile = \"%s\" ))",
                   # file.path(outDir, "hs37d5.fa"),
                   # file.path(outDir, "refFlat_hg19.txt.gz"),
                   # file.path(outDir, "codon.txt"),
                   # file.path(outDir, "priority.txt")))
  # invisible(NULL)
# }
