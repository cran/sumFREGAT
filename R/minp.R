# a wrapper to call functions of the 'GBJ' package

sumstat.minp <- function(obj) {
	#p <- minP_norestriction(test_stats = obj$df$Z, cor_mat = obj$U)$minP_pvalue
	p <- GBJ::minP(test_stats = obj$df$Z, cor_mat = obj$U)$minP_pvalue
	return(p)
}

'minp' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', write.file = FALSE, quiet = FALSE) {

############ COMMON CHECKS

tmp <- check.input(score.file, cor.path, gene.file, genes)
for (i in 1:length(tmp)) assign(names(tmp)[i], tmp[[i]])

############ SPECIFIC CHECKS

check.list <- get.check.list('minP', score.file, anno.type)

############ ANALYSIS

genewise(score.file, gene.file, gf, anno.type, cor.path, cor.file.ext, check.list, write.file, NULL, ncl = 3, quiet = quiet, test = 'minp')

}

# a function from 'GBJ', with changes

# minP_norestriction <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL)
# {
  # # Parse inputs, do some error checking.
  # param_list <- parse_input_norestriction(test_stats=test_stats, cor_mat=cor_mat,
                            # pairwise_cors=pairwise_cors)
  # t_vec <- param_list$t_vec
  # pairwise_cors <- param_list$pairwise_cors
  # d <- length(t_vec)

	# # minP bounds
  # minP_stat = 1-pchisq(t_vec[1]^2, df=1)
	# minP_p_bounds <- rep(minP_stat, d)
	# minP_z_bounds <- qnorm(1-minP_p_bounds/2)
	# minP_z_bounds <- sort(minP_z_bounds, decreasing=F)

	# minP_z_bounds[which(minP_z_bounds > 8.2)]= 8.2

	# # Send it to the C++.
	# if (sum(abs(pairwise_cors)) == 0) {
		# # For the independence flag in the c++, just have to send a number < -1.
		# minP_corp <- GBJ:::ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=rep(-999,2))
	# } else {
		# minP_corp <- GBJ:::ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=pairwise_cors)
	# }

	# return ( list(minP=minP_stat, minP_pvalue=minP_corp) )
# }

# # a function from 'GBJ', with changes

# parse_input_norestriction <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL) {

	# # Ensure that the thresholds are sorted in descending order, largest first.
	# t_vec <- sort(abs(test_stats), decreasing=TRUE)
	# d <- length(t_vec)

	# # Sometimes test stats are too big for R's precision
	# too_big <- which(t_vec > 8.2)
	# if (length(too_big) > 0) {t_vec[too_big] <- 8.2}

	# # Did they specify correlations?
	# if (is.null(cor_mat) & is.null(pairwise_cors)) {
		# stop("You must specify either cor_mat or pairwise_cors!")
	# }

	# # cor_mat specification gets priority
	# if (!is.null(cor_mat)) {
		# if(!isSymmetric(cor_mat)) {
			# stop("You did not specify a symmetric correlation matrix")
		# }

	  # # Put the correlation matrix into pairwise_cor
	  # pairwise_cors <- cor_mat[upper.tri(cor_mat)]
	# }

	# # Correct number of pairwise correlations?
	# if (length(pairwise_cors) != d*(d-1)/2) {
		# stop("Your pairwise correlation matrix/vector is of the wrong size!")
	# }

	# return( list(t_vec=t_vec, pairwise_cors=pairwise_cors))
# }
