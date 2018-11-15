# sumFREGAT (2017-2018) Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

'sumchi' <- function (score.file, gene.file, genes = 'all', cor.path = 'cor/',
anno.type = '', method = 'kuonen', acc = 1e-8, lim = 1e+6, write.file = FALSE) {

	SKAT(score.file, gene.file, genes, cor.path, anno.type, beta.par = c(1, 1), weights.function = NULL,
	user.weights = FALSE, gen.var.weights = 'none', method, acc, lim, write.file = write.file)

}