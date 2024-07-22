# run with R 4.3 for clusterProfiler to work
# module load R_packages/4.3.1 # on UPPMAX

# Gene ontology enrichment analysis
# Divided between subgroups: DMRs in (1) promoter, (2) exon, (3) intron,
# (4) CGI

projdir <- file.path("~/proj/PFOS/em_seq")

genes <- readRDS(file.path(projdir, "data", "PFOS_MCF-10A_DMG.Rds"))
meth <- readRDS(file.path(projdir, "data", "PFOS_MCF-10A_DMR.Rds"))

library(clusterProfiler) # doi: 10.1089/omi.2011.0118, 10.1016/j.xinn.2021.100141
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(cowplot)

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

go_analysis <- function(genelist, genes){
  # print(print(deparse(substitute(genelist)))) # print the name of the list?

  print("Running KEGG analysis...")
  go_kegg <- enrichKEGG(
    gene = genelist$entrez_id,
    organism = "hsa"
  ) %>% data.frame

  print("Running GO analysis...")
  go <- enrichGO(
    gene = genelist$gene_id,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    universe = genes$gene_id,
    readable = TRUE
  ) %>% data.frame

  print("Merging GO results...")
  # check if any sign GO list was found
  if (nrow(go) > 0 | nrow(go_kegg) > 0) {
    if (nrow(go_kegg) > 0) {
      go_kegg$ONTOLOGY <- "KEGG"
      
      # translate entrez_id to gene name
      go_kegg$geneID <- sapply(go_kegg$geneID, function(gene){
        genes <- strsplit(gene, "/")[[1]]
        genes <- subset(ens, entrez_id %in% genes)$gene_name
        paste(genes, collapse = "/")
      })
      
      go <- rbind(go, go_kegg)
    }
  }

  print("Sorting GO results...")
  # sort on qvalue & split the gene IDs
  go <- go[order(go$qvalue), ]
  go$ONTOLOGY <- factor(go$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))
  go$geneID <- sapply(go$geneID, function(gene){
    strsplit(gene, "/")[[1]]
  })
  go$gene_id <- sapply(go$geneID, function(gene){
    ens[ens$gene_name %in% gene, "gene_id"]
  })
  
  print("Calculating gene ranks...")
  # Rank genes by occurrence in the GO terms
  go_genes <- go$geneID %>% unlist %>% unique
  terms <- go$ONTOLOGY %>% unique

  d <- sapply(terms, function(x){
    tab <- subset(go, ONTOLOGY %in% x)$geneID %>% unlist %>% table
    rows <- match(names(tab), go_genes)
    
    dat <- vector(mode = "numeric", length = length(go_genes))
    names(dat) <- go_genes
    dat[rows] <- tab
    
    dat
  })
  colnames(d) <- terms

  # sort by freq. of genes in go terms by rank
  ranks <- apply(d, 2, rank) %>% data.frame
  ranks$tot <- rowSums(ranks)
  rows <- order(-ranks$tot)

  generanks <- data.frame(gene_name = go_genes, d)
  generanks <- generanks[rows,]
  rownames(generanks) <- NULL

  res <- list()
  res$go <- go
  res$generanks <- generanks

  return(res)
}

# import curated ensembl dataset
ens <- fread("GRCh38/ensembl_dataset_GRCh38.csv.gz")

geneset <- list()
geneset$promoter$gene_id <- subset(genes, dmr_in_promoter > 0)$gene_id %>% na.omit
geneset$promoter$entrez_id <- subset(genes, dmr_in_promoter > 0)$entrez_id %>% na.omit
geneset$exon$gene_id <- subset(genes, dmr_in_exon > 0)$gene_id %>% na.omit
geneset$exon$entrez_id <- subset(genes, dmr_in_exon > 0)$entrez_id %>% na.omit
geneset$intron$gene_id <- subset(genes, dmr_in_intron > 0)$gene_id %>% na.omit
geneset$intron$entrez_id <- subset(genes, dmr_in_intron > 0)$entrez_id %>% na.omit

cgi <- subset(meth, !is.na(cgi_id))$gene_id %>% na.omit
cgi <- sapply(cgi, function(x) strsplit(x, ",")[[1]]) %>% unlist %>% unname %>% unique

geneset$cgi$gene_id <- cgi
geneset$cgi$entrez_id <- ens[ens$gene_id %in% cgi, "entrez_id"] %>% na.omit

go_res <- lapply( geneset, function(x) go_analysis(x, genes) )

saveRDS(go_res, file = "data/PFOS_MCF-10A_GO.Rds")