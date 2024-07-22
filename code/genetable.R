# generate a meth table based on genes and CGI rather then DMRs

# (1) DMG-table should contain info about which genomic regions DMRs are present in
# (2) CGI table should contain overlapping DMRs, these DMRs can be pulled from 
# DMR-table to identify genes

library(data.table)
library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(genomation)

# = functions ================================================================ #

sort_chrom <- function(d) {
  d$chr <- sub("chr", "", d$chr)

  # Sort table after chromosomes > start pos
  d$chr <- sub("MT", 100, d$chr)
  d$chr <- sub("M", 100, d$chr)
  d$chr <- sub("X", 101, d$chr)
  d$chr <- sub("Y", 102, d$chr)
  d <- d[order(as.numeric(d$chr), d$start), ]
  d$chr <- sub(100, "MT", d$chr)
  d$chr <- sub(101, "X", d$chr)
  d$chr <- sub(102, "Y", d$chr)
  rownames(d) <- NULL
  d$chr <- sub("^", "chr", d$chr)
  
  return(d)
}

build_genetable <- function(gene){
  print(paste("Testing", gene, "| # = ", which(genes %in% gene)))
  
  gene_info <- ens[ens$gene_id %in% gene, ]
  meth_info <- meth[meth$chr %in% gene_info$chr, ]

  gr_gene <- as(gene_info[, 1:4], "GRanges")
  gr_meth <- as(meth_info[, 1:3], "GRanges")

  rows <- findOverlaps(gr_gene, gr_meth) %>% data.frame()
  
  d <- meth_info[rows$subjectHits, ]

  if ( any(d$sign) ){
    tmp <- d[d$sign == TRUE, ]
    
    gr_tmp <- as(tmp[, 1:3], "GRanges")

    # introns
    rows <- findOverlaps(annotations$intron, gr_tmp) %>% data.frame()
    tmp$region <- "intragenic"
    tmp$region[rows$subjectHits] <- "intron"

    # exon
    rows <- findOverlaps(annotations$exon, gr_tmp) %>% data.frame()
    tmp$region[rows$subjectHits] <- "exon"

    # promoter
    rows <- findOverlaps(annotations$promoter, gr_tmp) %>% data.frame()
    tmp$region[rows$subjectHits] <- "promoter"

    # cgi
    rows <- findOverlaps(annotations$cgi, gr_tmp) %>% data.frame()
    tmp$cgi_id <- NA
    tmp$cgi_id[rows$subjectHits] <- annotations$cgi$id[rows$queryHits]
    tmp$cgi_id <- ifelse(!is.na(tmp$cgi_id), paste0("cgi_id:", tmp$cgi_id), tmp$cgi_id)

    # output table
    out <- gene_info
    
    out$DMG <- TRUE
    out$num_cpg <- sum(tmp$feature == "cpg")
    out$num_tile <- sum(tmp$feature == "tile100")

    out$dmr_in_cgi <- sum(!is.na(tmp$cgi_id))
    out$dmr_in_promoter <- sum(tmp$region == "promoter")
    out$dmr_in_exon <- sum(tmp$region == "exon")
    out$dmr_in_intron <- sum(tmp$region == "intron")
    out$dmr_in_intragenic <- sum(tmp$region == "intragenic")
    
    out$dmr_id_cgi <- list(tmp$dmr_id[!is.na(tmp$cgi_id)])
    out$dmr_id_promoter <- list(tmp$dmr_id[tmp$region == "promoter"])
    out$dmr_id_exon <- list(tmp$dmr_id[tmp$region == "exon"])
    out$dmr_id_intron <- list(tmp$dmr_id[tmp$region == "intron"])
    out$dmr_id_intragenic <- list(tmp$dmr_id[tmp$region == "intragenic"])

  }else{

    out <- cbind(
      gene_info,
      data.frame(FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA,NA)
    )

    colnames(out)[ (ncol(gene_info) + 1) : ncol(out) ] <- c(
      "DMG", "num_cpg", "num_tile", 
      "dmr_in_cgi", "dmr_in_promoter", "dmr_in_exon", "dmr_in_intron", "dmr_in_intragenic", 
      "dmr_id_cgi", "dmr_id_promoter", "dmr_id_exon", "dmr_id_intron", "dmr_id_intragenic"
    )
  }

  return(out)
}

# = code ===================================================================== #

get_args()

# table from DMR persective
projdir <- file.path("~/proj/PFOS/em_seq")
dumpdir <- file.path(projdir, "dump")

# import curated ensembl dataset
ens <- fread("GRCh38/ensembl_dataset_GRCh38.csv.gz")
ens$strand <- ifelse(ens$strand > 0, "+", "-")

meth <- readRDS(file.path(projdir, "data", "PFOS_MCF-10A_DMR.Rds"))

# import genomic features

# Genes table, each row = gene w/ info about overlapping DMRs
# RefSeq
file <- file.path("GRCh38", "refseq_UCSC_GRCh38.bed")
annotations <- readTranscriptFeatures(file)

# CpG-islands table, each row = CGI w/ info about overlapping DMRs
file <- file.path("GRCh38", "cpgislands_GRCh38.bed")
CGI <- readGeneric(file, keep.all.metadata = TRUE)

# only use well annotated chromsomes
chrom <- c(1:22, "X", "Y")
seqlevels(CGI, pruning.mode="coarse") <- seqlevels(CGI)[seqlevels(CGI) %in% paste0("chr", chrom)]
values(CGI) <- DataFrame(id = 1:length(CGI), values(CGI)) # assign each cgi an ID, numeric

annotations$cgi <- CGI

genes <- unique(ens$gene_id)

l1 <- lapply(genes, build_genetable)
genes <- Reduce(function(x,y) rbind(x,y), l1)
genes <- sort_chrom(genes)


# = save results ============================================================= #

name <- "PFOS_MCF-10A_DMG.Rds"
filename <- file.path(projdir, "data", name)
saveRDS(genes, file = filename)