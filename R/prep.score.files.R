# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

prep.score.files <- function(input.file, reference.file = '', output.file.prefix) {#, annotate = FALSE, path.to.reference, annotation.param) { # 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'P', 'BETA', 'EAF'
	df <- read.table(input.file, header = TRUE, as.is = TRUE)
	colnames(df) <- toupper(colnames(df))
	v <- which(colnames(df) %in% c('CHR', 'CHROMOSOME'))
	if (length(v) == 1) colnames(df)[v] <- 'CHROM'
	v <- which(colnames(df) %in% c('POSITION', 'POSITIONS', 'MAP'))
	if (length(v) == 1) colnames(df)[v] <- 'POS'
	v <- which(colnames(df) %in% c('PVALUE', 'PV', 'PVAL', 'P.VALUE', 'P_VALUE'))
	if (length(v) == 1) colnames(df)[v] <- 'P'
	v <- which(colnames(df) %in% c('RSID', 'RS.ID', 'RS_ID', 'SNP.ID', 'SNP_ID'))
	if (length(v) == 1) colnames(df)[v] <- 'ID'
	v <- which(colnames(df) == 'EA')
	if (length(v) == 1) colnames(df)[v] <- 'EFFECT.ALLELE'
	
	# ID and PVAL mandatory
	# others from user file or 1KG
	
	ColNames <- c('ID', 'P')
	v <- !ColNames %in% colnames(df)
	if (sum(v)) stop(paste("Mandatory column(s) missing:", paste(ColNames[v], collapse = ', ')))
	
	
	ColNames <- c('CHROM', 'POS', 'EAF')
	v <- !ColNames %in% colnames(df)
	take <- ColNames[v]
	if (sum(v)) print(paste("Columns that are missing and will be looked for in reference data:", paste(take, collapse = ', ')))
	take[take == 'EAF'] <- 'AF'

	if ('BETA' %in% colnames(df)) {
		if ('EFFECT.ALLELE' %in% colnames(df)) {
			if(!'REF' %in% colnames(df)) take <- c(take, 'REF', 'ALT')
		} else {
			print("Effect allele column not found, effect sizes cannot be linked")
		}
	} else {
		print("Effect sizes (beta) column not found")
	}
	if (length(take) > 0) {
		if (file.exists(reference.file)) {
			print('Reading reference file...')
			ref <- get(load(reference.file))
			if ('CHROM' %in% take & !'CHROM' %in% colnames(ref)) stop ("No CHROM column in reference data")
			if ('POS' %in% take & !'POS' %in% colnames(ref)) stop ("No POS column in reference data")
			v <- match(df$ID, ref$ID)
			
			if (!sum(v, na.rm = TRUE)) {
				if (all(c('CHROM', 'POS') %in% colnames(df))) {
					df$ind <- paste(df$CHROM, df$POS, sep = ':')
					print('No IDs matching, trying to link through map data...')
					ref$ind <- paste(ref$CHROM, ref$POS, sep = ':')
					v <- match(df$ind, ref$ind)
					if (sum(!is.na(v)) < (length(v) / 2)) {
						print("Too few variants match between input file and reference data")
						v <- NA
					}
				}
			}
			if (sum(v, na.rm = TRUE)) {
				print(paste(sum(!is.na(v)), "of", length(v), "variants found in reference"))
				vv <- take %in% colnames(ref)
				if (sum(!vv)) {
					print(paste("Columns that are missing in reference data:", paste(take[!vv], collapse = ', ')))
					if ('REF' %in% take & !'REF' %in% colnames(ref)) print ("Reference alleles not found, effect sizes cannot be linked")
					if ('AF' %in% take & !'AF' %in% colnames(ref)) print ("Allele frequencies not found, some weighted tests will be unavailable")
				}
				df <- cbind(df, ref[v, take[vv]])
			}
		} else {
			if (reference.file != '') print ("Reference file not found")
			if (any(c('CHROM', 'POS') %in% take)) stop ("Cannot find map data (chromosome, position)")
		}
	}

	if ('REF' %in% colnames(df)) {
		v <- df$EFFECT.ALLELE == df$REF
		df[is.na(v), 'BETA'] <- NA
		if ('EAF' %in% colnames(df)) df[is.na(v), 'EAF'] <- NA
		if ('ALT' %in% colnames(df)) {
			vv <- !v & df$EFFECT.ALLELE != df$ALT
			if (sum(vv, na.rm = T)) {
				print(paste("Effect alleles do not match reference/alternative alleles for", sum(vv), "variant(s)"))
				df[vv, 'BETA'] <- NA
				if ('EAF' %in% colnames(df)) df[vv, 'EAF'] <- NA
			}
		}
		#here we go
		v <- which(v)
		df$BETA[v] <- -df$BETA[v]
		if ('EAF' %in% colnames(df)) {
			df$EAF[v] <- 1 - df$EAF[v]
			colnames(df)[colnames(df) == 'EAF'] <- 'AF'
		}
		print(paste('Effect sizes recoded for', length(v), 'variant(s)'))
	}

	df$Z <- qnorm(df$P / 2, lower.tail = FALSE)
	if ('BETA' %in% colnames(df)) {
		df$Z <- df$Z * sign(df$BETA)
		df$SE.BETA <- df$BETA / df$Z
	}

	df <- df[order(df[, 'POS']), ]
	df <- df[order(df[, 'CHROM']), ]
	if (!'ALT' %in% colnames(df)) df$ALT <- NA
	if (!'REF' %in% colnames(df)) df$REF <- NA
	vcf <- df[, c('CHROM', 'POS', 'ID', 'REF', 'ALT')]
	colnames(vcf)[1] <- '#CHROM'
	vcf$POS <- format(vcf$POS, scientific = FALSE)
	vcf$POS <- gsub(' ', '', vcf$POS)
	vcf <- cbind(vcf, QUAL = '.', FILTER = '.')
	vcf$INFO <- paste0('Z=', df$Z)
	title <- c('##INFO=<ID=Z,Number=1,Type=Float,Description="Z statistics">')

	if ('BETA' %in% colnames(df)) {
		vcf$INFO <- paste0(vcf$INFO, ';SE.Beta=', df$SE.BETA)
		title <- c(title, '##INFO=<ID=SE.Beta,Number=1,Type=Float,Description="SE Beta">')
	}

	if ('EAF' %in% colnames(df)) colnames(df)[colnames(df) == 'EAF'] <- 'AF'
	if ('AF' %in% colnames(df)) {
		vcf$INFO <- paste0(vcf$INFO, ';AF=', df$AF)
		title <- c(title, '##INFO=<ID=AF,Number=1,Type=Float,Description="Frequency of alternative allele">')
		print(paste0('Allele frequencies found and linked'))
	}

	a <- grep('ANNO', colnames(df))
	if (length(a) == 1) {
		vcf$INFO <- paste0(vcf$INFO, ';ANNO=', df[, a])
		title <- c(title, '##INFO=<ID=ANNO,Number=1,Type=String,Description="Variants annotations">')
		print(paste0('Annotations ("', colnames(df)[a], '") found and linked'))
	}

	a <- grep('\\bW', colnames(df))
	if (length(a) == 1) {
		vcf$INFO <- paste0(vcf$INFO, ';W=', df[, a])
		title <- c(title, '##INFO=<ID=W,Number=1,Type=Float,Description="Weights">')
		print(paste0('User weights ("', colnames(df)[a], '") found and linked'))
	}

	a <- grep('\\bN', colnames(df))
	if (length(a) == 1) {
		vcf$INFO <- paste0(vcf$INFO, ';N=', df[, a])
		title <- c(title, '##INFO=<ID=N,Number=1,Type=Integer,Description="Sample size">')
		print(paste0('Sample size info ("', colnames(df)[a], '") found and linked'))
	}

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

