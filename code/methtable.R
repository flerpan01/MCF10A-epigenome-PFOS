# import all differentally methylated tables from the folder data
# merge into one file
# DMR-table should contain overlapping gene IDs

library(data.table)
library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(genomation)

# = functions ================================================================ #

# FIX: chromosort finns i en "snyggare" funktion
sort_chrom <- function(d){
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

build_methtable <- function(dmr_type, ens){
  print(paste("building dmr table =", dmr_type))

  files <- list.files(
    dumpdir,
    pattern = dmr_type,
    full.names = TRUE
  )

  files <- data.frame(
    files = files,
    stringsAsFactors = FALSE
  )

  files$name <- sapply(files$files, function(x) {
    x <- strsplit(x, "[/]")[[1]]
    x <- x[length(x)]
    sub(".csv.gz" ,"", x)
  })

  files$dose <- sapply(files$name, function(x) {
    strsplit(x, "_")[[1]][2]
  })

  # import csv
  meth <- lapply(files$files, fread)
  meth <- Reduce(function(x,y) rbind(x,y), meth)

  # add sign threshold
  if (dmr_type == "cpg") meth$sign <- ifelse(meth$qvalue <= 0.05 & abs(meth$meth.diff) >= 15, TRUE, FALSE)
  if (dmr_type == "tile100") meth$sign <- ifelse(meth$qvalue <= 0.05 & abs(meth$meth.diff) >= 5, TRUE, FALSE)

  meth <- meth[meth$sign == TRUE, ] # reduce computation time
  #meth <- meth[ (nrow(meth) - 9) : nrow(meth),]

  # overlap meth data and ensembl gene data
  gr_ens <- as(ens[, 1:4], "GRanges")
  
  meth$gene_id <- apply(meth, 1, function(x){
    #show(x[9])
    gr_meth <- as( as.data.frame(t(x[1:3])) , "GRanges")
    rows <- findOverlaps(gr_ens, gr_meth) %>% data.frame()
    #show(nrow(rows))
    
    if (nrow(rows) == 0) out <- NA
    if (nrow(rows) > 0) out <- paste(ens$gene_id[rows$queryHits], collapse = ",")
    
    return(out)
  })

  # CpG-islands
  file <- file.path("~/genomes/GRCh38", "cpgislands_GRCh38.bed")
  cgi <- readGeneric(file, keep.all.metadata = TRUE)
  
  # only use well annotated chromsomes
  chrom <- c(1:22, "X", "Y")
  seqlevels(cgi, pruning.mode = "coarse") <- seqlevels(cgi)[seqlevels(cgi) %in% paste0("chr", chrom)]
  cgi$id <- paste0("cgi:", 1:length(cgi))

  # add cpgisland data, include as ID
  gr_meth <- as(meth[, 1:3], "GRanges")
  rows <- findOverlaps(cgi, gr_meth) %>% data.frame()
  meth$cgi_id <- NA
  meth$cgi_id[rows$subjectHits] <- cgi$id[rows$queryHits]

  # intergenic
  meth$region <- ifelse(is.na(meth$gene_id), "intergenic", "intragenic")

  # RefSeq
  file <- file.path("~/genomes/GRCh38", "refseq_UCSC_GRCh38.bed")
  annotations <- readTranscriptFeatures(file)

  # introns
  rows <- findOverlaps(annotations$intron, gr_meth) %>% data.frame()
  meth$region[rows$subjectHits] <- "intron"

  # exon
  rows <- findOverlaps(annotations$exon, gr_meth) %>% data.frame()
  meth$region[rows$subjectHits] <- "exon"

  # promoter
  rows <- findOverlaps(annotations$promoter, gr_meth) %>% data.frame()
  meth$region[rows$subjectHits] <- "promoter"

  # FIX, if no gene id set to intergenic
  meth$region <- ifelse(is.na(meth$gene_id), "intergenic", meth$region)

  return(meth)
}


# = code ===================================================================== #

# table from DMR persective
projdir <- file.path("~/proj/PFOS/em_seq")
dumpdir <- file.path(projdir, "dump")

# import curated ensembl dataset
ens <- fread("GRCh38/ensembl_dataset_GRCh38.csv.gz")
ens$strand <- ifelse(ens$strand > 0, "+", "-") # FIX!

print("Running analysis on DMR")

# Methylation table, each row = DMR
# divide between tile and cpg, make 2 meth tables and merge
dmr_types <- c("cpg", "tile100")

l1 <- lapply(dmr_types, function(dmr_type) build_methtable(dmr_type, ens))
meth <- Reduce(function(x,y) rbind(x,y), l1)

meth <- sort_chrom(meth)

# = save results ============================================================= #

name <- paste0("PFOS_MCF-10A_DMR.Rds")
filename <- file.path(projdir, "data", name)
saveRDS(meth, file = filename)