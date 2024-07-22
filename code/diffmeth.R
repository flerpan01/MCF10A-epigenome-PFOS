# Differential methylation analysis
# dose = 10µM PFOS
# DMR type = cpg, tile 100bp. Arguments are sent in via R CMD BATCH '--args' script.R

# Datasets are divided into 2 methylBaseDB objects based on:
# (1) CpG-site resolution (1bp) [DM cpgsite], coverage of 67%
# (4 / 6 samples with coverage)
# (2) region resolution (tilesize) [DMRs]
# (3) dose exposure, [low, high], low vs control & high vs control
# calculateDiffMeth(obj, overdispersion = "MN", padjust = "SLIM", test = "Chisq")

# module load R_packages/4.1.1

library(dplyr)
library(methylKit)
library(data.table) # fread, fwrite

library(ggplot2)
library(cowplot)
library(ggrepel)

# Code ran on HPC
projdir <- file.path("~/proj/PFOS/em_seq")
covfilesdir <- file.path(projdir, "results/bismark_methylation_calls/methylation_coverage")
reportdir <- file.path(projdir, "reports")

chromosomes <- paste0("chr", c(1:22, "X"))
nsamples <- 0.6

# = functions ================================================================ #

get_arguments <- function(){
	args <- commandArgs(TRUE)
	
	for (i in seq_along(args)){
		# add quotes to args or the parse will not work
		VAR <- sub("=", "='", args[i]) 	# behind =
		VAR <- sub("$", "'", VAR)				# last char
		
		eval(parse(text = VAR), envir = .GlobalEnv)
	}
	
	# quickfix. quotation problems taking arguments from CLI into R
	dose <<- as.numeric(dose)
	tile <<- if(tile == "NULL") NULL else as.numeric(tile)

	print(dose)
	print(tile)
}

pca_plot <- function(meth){
	res <- PCASamples(meth, obj.return = TRUE)

	data <- cbind(res$x[,1:2], treatment = obj@treatment) %>% data.frame

	data$PC1 <- as.numeric(data$PC1) / res$sdev[1] * sqrt(nrow(data))
	data$PC2 <- as.numeric(data$PC2) / res$sdev[2] * sqrt(nrow(data))
	data$treatment <- as.factor(data$treatment)
	data$id <- rownames(data)

	p <- ggplot(data, aes(PC1, PC2, fill = treatment)) +
	  geom_point(shape = 21, size = 3) +
	  geom_label_repel(aes(label = id, fill = treatment), size = 3, show.legend = FALSE) +
	    coord_fixed() +
	    labs(
	    	#title = title,
	    	#x = paste("PC1:", var_perc[1], "% var."),
	      #y = paste("PC2:", var_perc[2], "% var")
	    	fill = "") +
	    theme_linedraw(base_size = 20)
	  #stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))), data = df[df$V3 == "1" | df$V3 == "2",], size = 1) 

	return(p)
}

hist_plot <- function(std){
	ggplot() +
		geom_histogram(aes(std), color = "black", bins = 50) + 
		theme_linedraw(base_size = 20)
}

# REWRITE function! sample info is hardcoded
build_metadata <- function() {
  files <- list.files(covfilesdir, pattern = "*.gz", full.names = TRUE)

  metadata <- data.frame(
  	sample_id = c(paste0("C", 1:6), paste0("low", 1:6), paste0("high", 1:6)),
  	exposure = c(rep(0, 6), rep(1, 6), rep(2, 6)),
  	dose = c(rep(0, 6), rep(1, 6), rep(10, 6)),
  	files = files
  )

  return(metadata)
}

# import, call common cpg pos, filter, normalise
build_methylObj <- function(metadata, tile = NULL, chromosomes, cores = 10) {
	filelist <- split(
		as.character(metadata$files),
		seq_along(metadata$files)
	)

	samplelist <- split(
		as.character(metadata$sample_id),
		seq_along(metadata$sample_id)
	)

	obj <- methRead(
		location = filelist,
		sample.id = samplelist,
		treatment = metadata$exposure,
		header = FALSE,
		assembly = "GRCh38",

		# "amp", "bismark","bismarkCoverage", "bismarkCytosineReport"
		pipeline = "bismarkCoverage"
	)
	
	# cpg_ref <- pull_cpg_ref() # 29.4M CG-motifs in human genome hg38
	cpg_ref <- fread(file.path(projdir, "GRCh38/cg_pos_CRGh38.csv.gz"))

	# Filter CpG-site position that match the ref genome
  # iterate over each sample
  for (i in seq_along(obj)) {
    print(paste("Filtering sample:", i))

    cpg_id <- paste0(obj[[i]]$chr, ":", obj[[i]]$start)
    rows <- cpg_id %in% cpg_ref$pos # cpg_ref$id?

    obj[[i]] <- obj[[i]][rows, ]
  }
	
	if (!is.null(tile)){
		print(paste("calculating tiles w/ size:", tile))
		
		tile <- as.numeric(tile)

		obj <- tileMethylCounts(
			obj,
			win.size = tile,
			step.size = tile,
			cov.bases = 2,
			mc.cores = cores
		)
	}
	
	# Filter reads
	obj <- filterByCoverage(
		obj, 
		lo.count = 10,	# < 10
		hi.perc = 99.9 	# the top 99th percentile (PCR duplicates)
	) 

	# Normalisation, scaling factor between samples based on differences
	# between median of coverage distribution
	obj <- normalizeCoverage(
	  obj,
	  method = "median"
	)

	return(obj)
}

# unite, std filter > 2
make_methdata <- function(obj, tile = NULL, nsamples = 1L, save_methylationMatrix = TRUE, make_plots = FALSE, cores = 10) {
	dmr_type <- ifelse(is.null(tile), "cpg", paste0("tile", tile))

	print(paste("DMR type:", dmr_type))

	if (!is.integer(nsamples)) nsamples <- as.integer(ceiling(table(metadata$exposure) %>% mean * nsamples))

	print(paste("Minimum number of samples per group:", nsamples))

	meth <- methylKit::unite(
		obj, 
		min.per.group = nsamples,
		mc.cores = cores
	)
	
	print(paste("meth nrow =", nrow(meth)))


	## FIX HERE, 
	## move the dmr_id to keep as csv to complementary_files.R
	## then import
	if (!is.null(tile)){

		d <- data.frame(
			chr = meth$chr, 
			start = meth$start, 
			end = meth$end
		)

		d$dmr_id <- paste0(d$chr, ":", d$start)

		# remove mitochondrial chromosome
		d <- subset(d, chr != "chrM")

    d$num_cpg <- apply(d, 1, function(pos){
    	#chrom <- if (pos[1] == "chrMT") "chrM" else pos[1]
    	chrom <- pos[1]

    	# add 1 bp before and after as CG-motif are 2 bp long
    	start <- as.numeric(pos[2]) - 1
    	end <- as.numeric(pos[3]) + 1

    	countPattern("CG", genome[[chrom]][start:end])
    })

    if (any(d$num_cpg < 2)){
    	id <- subset(d, num_cpg < 2)$dmr_id
	    tiles <- paste0("chr", getData(meth)$chr, ":", getData(meth)$start)
	    rows <- tiles %in% id

	    meth <- meth[!rows, ]
    }
	}

	# remove CpG-sites with low variation, std < 2
	mat <- percMethylation(meth)
	std <- matrixStats::rowSds(mat, na.rm = TRUE)

	if (make_plots) {
		pca_before <- pca_plot(meth)
		hist_before <- hist_plot(std)
	}
	
	meth <- meth[std > 2]

	print(paste("# of CpG-sites removed (due to low variation):", sum(std < 2)))
	print(paste("# of CpG-sites passed the filtering:", nrow(meth)))

	if (make_plots) {
		pca_after <- pca_plot(meth)
		hist_after <- hist_plot(std[std > 2])

		# The plot 2 x 2
		title1 <- ggdraw() +
			draw_label("PCA and hist of std", x = 0.05, hjust = 0, vjust = 1)

		title2 <- ggdraw() +
			draw_label("std < 2 CpG-sites removed", x = 0.05, hjust = 0, vjust = 1)

		top <- plot_grid(pca_before, hist_before, nrow = 1, labels = c("A", "B"))
		bottom <- plot_grid(pca_after, hist_after, nrow = 1, labels = c("C", "D"))

		p <- plot_grid(
			title1,
			top,
			title2,
			bottom,
			ncol = 1,
			rel_heights = c(0.1, 1, 0.1, 1))


		filename <- file.path(
			reportdir,
			paste0("meth_stats_", dose, "_", dmr_type, ".pdf"))

		pdf(filename, width = 12)
		print(p)

		# clusterplot
		clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
		clusterSamples(meth, dist = "euclidean", method = "ward", plot = TRUE)

		dev.off()
	}

	if (save_methylationMatrix){
		print("Saving beta values matrix")
		dmr_id <- paste0(getData(meth)$chr, ":", getData(meth)$start)

		percMat <- percMethylation(meth)
		rownames(percMat) <- dmr_id

		dmr_type <- ifelse(is.null(tile), "cpg", paste0("tile", tile))

		name <- paste0("PFOS_MCF-10A_betavalues_matrix_", dmr_type, ".Rds")
		filename <- file.path(projdir, "data", name)

		saveRDS(percMat, file = filename)
	}

	return(meth)
}

meth_diff_analysis <- function(meth, tile = NULL, cores = 10) {
  dat <- calculateDiffMeth(
  	meth,
    adjust = "SLIM", 
    #effect = "predicted",
    overdispersion = "MN", 
    test = "Chisq",
    mc.cores = cores
  )

  print(paste("calc diff meth done"))

  # convert to dataframe
  data <- getData(dat)

  # add type and cpg_id
  data$type <- ifelse(data$meth.diff > 0, "hyper", "hypo")
  #data$chr <- paste0("chr", data$chr)
  data$dmr_id <- paste0(data$chr, ":", data$start)
  data$feature <- ifelse(is.null(tile), "cpg", paste0("tile", tile))
  data$dose <- dose
    
  if (is.null(tile)) data$num_cpg <- 1
  

  ## FIX THIS!
  ## already calculated this in the function make_methdata, above
  # remove mitochondrial chromosome
	#data <- subset(data, chr != "chrM")

  if (!is.null(tile)){
  	data$num_cpg <- apply(data, 1, function(pos){
  		#chrom <- if (pos[1] == "chrMT") "chrM" else pos[1]
  		chrom <- pos[1]

	  	# add 1 bp before and after as CG-motif are 2 bp long
	  	start <- as.numeric(pos[2]) - 1
	  	end <- as.numeric(pos[3]) + 1

	  	countPattern("CG", genome[[chrom]][start:end])
  	})
	}

  return(data)
}

# = code ===================================================================== #

get_arguments()

if (!is.null(tile)){
	# load human genome
	library("BSgenome.Hsapiens.UCSC.hg38")
	genome <- BSgenome.Hsapiens.UCSC.hg38
}

metadata <- build_metadata()
metadata <- subset(metadata, dose %in% c(0, 1))

obj <- build_methylObj(metadata, tile, chromosomes)
meth <- make_methdata(obj, tile, nsamples, make_plots = TRUE)
data <- meth_diff_analysis(meth, tile)

# = save output ============================================================== #

dmr_type <- ifelse(is.null(tile), "cpg", paste0("tile", tile))

name <- paste0("diffmeth_", dose, "_", dmr_type, ".csv.gz")
filename <- file.path(projdir, "dump", name)
fwrite(data, file = filename)

# stats
out <- c(
	dose = dose,
	dmr = dmr_type, 
	nrow = nrow(data),
	sign_0.05 = nrow(subset(data, qvalue <= 0.05)),
	sign_0.01 = nrow(subset(data, qvalue <= 0.01)),
	sign_0.05_5 = nrow(subset(data, qvalue <= 0.05 & abs(meth.diff) >= 5)), 	# for tile
	sign_0.01_5 = nrow(subset(data, qvalue <= 0.01 & abs(meth.diff) >= 5)), 	# for tile
	sign_0.05_15 = nrow(subset(data, qvalue <= 0.05 & abs(meth.diff) >= 15)),	# for cpg
	sign_0.01_15 = nrow(subset(data, qvalue <= 0.01 & abs(meth.diff) >= 15))	# for cpg
)	

print(out)