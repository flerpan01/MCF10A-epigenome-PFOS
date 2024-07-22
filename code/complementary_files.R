# R/4.1.1

# Generate files for subsequent analysis
# 1. ensembl dataset
# 2. table of CpG-positions
# 3. sequencing depth table

# = variables ================================================================ #

# select organism
organism_dataset <- "hsapiens_gene_ensembl"

# select chromosomes
chromosomes <- c(1:22, "MT", "X", "Y")

projdir <- file.path("~/projdir/PFOS/em_seq")

# = libraries ================================================================ #

library(biomaRt)
library(data.table)
library(dplyr)

# namespace = Hsapiens / BSgenome.Hsapiens.UCSC.hg38
library(BSgenome.Hsapiens.UCSC.hg38) 

# = functions ================================================================ #

sort_chrom <- function(d){
  d$chr <- sub("chr", "", d$chr)

  # Sort table after chromosomes > start pos
  d$chr <- sub("MT", 100, d$chr)
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

# add extra column which is concatinated version of gene types
get_ensembl <- function(chromosomes, organism_dataset = NULL){
  # listEnsembl() # list available datasets
  # mart <- useEnsembl(biomart="genes") # download all gene lists
  # searchDatasets(mart=mart, pattern="mus") # identify house mouse dataset
  # mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  # listAttributes(mart) # list of available attributes

  if (is.null(organism_dataset)) organism_dataset <- "hsapiens_gene_ensembl"

  print(paste("Loading dataset:", organism_dataset, "..."))

  # get mart
  mart <- useEnsembl(
    #mirror = "www", # values: uswest, useast, asia, www
    biomart = "genes",
    dataset = organism_dataset
  ) 

  # columns to import
  #attr <- listAttributes(mart)$name # holds all columns
  cols <- c(
    "chromosome_name", 
    "start_position", 
    "end_position", 
    "strand", 
    "ensembl_gene_id", 
    "external_gene_name",
    #"entrezgene_id",
    "description", 
    "gene_biotype")

  # import columns
  ens <- getBM(
    mart = mart, 
    attributes = cols
  )
  ens <- ens[, cols]
  #ens <- data.table(ens)

  print(paste("Generating ensembl table..."))

  colnames(ens) <- c(
    "chr", 
    "start", 
    "end", 
    "strand", 
    "gene_id",     
    "gene_name",
    #"entrez_id",  
    "gene_info", 
    "gene_type")

  ens <- subset(ens, chr %in% chromosomes)
  ens <- sort_chrom(ens)

  # remove un-needed info in gene_info column
  ens$gene_info <- sapply(ens$gene_info, function(x){
    gsub(" \\[.*\\]", "", x)
  })

  # reduce gene types, collapse pseudogenes
  ens$gene_type2 <- ens$gene_type
  rows <- grep("pseudo", ens$gene_type2)
  ens$gene_type2[rows] <- "pseudogene"
  rows <- ens$gene_type %in% c(
    "snoRNA", "misc_RNA", "sRNA", "scaRNA", "snoRNA",
    "snRNA", "scRNA", "miRNA")
  ens$gene_type2[rows] <- "ncRNA"

  # 400+ IG / TR genes hid in protein coding
  rows <- grep("_gene", ens$gene_type2)
  ens$gene_type2[rows] <- "protein_coding"

  ens$size <- ens$end - ens$start

  print(paste("Total number of genes imported =", nrow(ens)))
  
  return(ens)
}

# = code ===================================================================== #

if (!dir.exists("GRCh38")) dir.create("GRCh38")

# 1. Pull the latest gene list from the ensembl database
ens <- get_ensembl(chromosomes, organism_dataset)

#filename <- file.path(projdir, "data/ensembl_dataset_GRCm39.csv.gz")
filename <- file.path("GRCh38/ensembl_dataset_GRCh38.csv.gz")
fwrite(ens, file = filename)

# 2. Calculate all CpG positions on the reference genome
genome <- BSgenome.Hsapiens.UCSC.hg38

cpg <- list() # empty list to fill
# iterate over chosen chromosomes
for (chrom in chromosomes){
  chrom <- paste0("chr", chrom)
  if (chrom == "chrMT") chrom <- "chrM"
  
  print(chrom)
  
  out <- matchPattern("CG", genome[[chrom]]) %>%
    as.data.frame()

  out$chr <- chrom
  cpg[[chrom]] <- out # Save in a list
}
cpg <- Reduce(function(x, y) rbind(x, y), cpg) # convert list to data frame
cpg <- cpg[, c("chr", "start", "end", "seq")]
cpg$id <- paste0(cpg$chr, ":", cpg$start) # cpg_id

#filename <- file.path(projdir, "data/cg_pos_CRGh38.csv.gz")
filename <- file.path("GRCh38/cg_pos_CRGh38.csv.gz")
fwrite(cpg, file = filename)

# 3. Calculate the sequencing depth
#d <- fread("reports/seqkit_stats.txt")

# sum_len = total num reads * read lengths
#d$sum_len <- as.numeric(gsub(",", "", d$sum_len))
#d$num_seqs <- as.numeric(gsub(",", "", d$num_seqs))

# depth = sum_len / genome size
#d$depth <- d$sum_len / 3.1e9 

# dose names
#d$dose <- sub("Sample_UA-2815-", "", d$file)
#d$dose <- sub("_R.*", "", d$dose)
#d$dose <- ifelse(grepl("C", d$dose), "control", d$dose)
#d$dose <- ifelse(grepl("PFOS-10", d$dose), "10", d$dose)
#d$dose <- ifelse(grepl("PFOS-1", d$dose), "1", d$dose)

#out <- d %>%
#  group_by(dose) %>%
#  summarise(depth_avg = mean(depth), 
#            depth_median = median(depth),
#            num_seq_avg = mean(num_seqs)) %>%
#  mutate(across(where(is.numeric), round, 2))

#fwrite(out, file = "reports/seq_depth.txt", sep = "\t")