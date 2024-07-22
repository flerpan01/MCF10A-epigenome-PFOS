The code in the repository was used to for computational analysis for the project: [Perfluorooctanesulfonic acid (PFOS) induced cancer related methylome alterations in human breast cells](doi)

# Project summary

+ Dataset is composed of Perflourooctanesulfonic acid (PFOS) exposed (72h) normal human breast cell line (MCF-10A). 12 samples, 6 control, 6 exposed.
+ The sequencing data is built upon the study: [PFOS induces proliferation, cell-cycle progression, and malignant phenotype in human breast epithelial cells](https://doi.org/10.1007/s00204-017-2077-8). From this study only the control and 1µM cells cultures were used.
+ As PFOS has potentially [endocrine disrupting estrogen activity](https://doi.org/10.1111/j.1365-2605.2008.00870.x) and has been [linked to breast cancer](https://doi.org/10.1186/1476-069X-10-88), the hypothesis is to identify DM regions which overlap genes related to the phenotypes listed above.
+ DNA methylation levels were mapped with Enzymatic Methyl sequencing (EM-seq)

>**The Aim** was to elucidate methylation alteration induced by PFOS exposure, aka "methylome fingerprint"

# Analysis workflow

Project file structure

```
project/
├── code
├── data
├── dump
├── merged
├── results
└── README.md
```

## Mapping & coverage

Run [nf-core methylseq pipeline](https://nf-co.re/methylseq) to align the sequencing reads to the reference genome and generate the methylation coverage files.

```sh
# Ran on HPC (UPPMAX). Use login node or pipeline is killed when node is killed

# setup env
PROJECT=""
EMAIL=""
PROJDIR=/home/$USER/proj/PFOS/em_seq/
cd $PROJDIR

# Load modules
ml bioinfo-tools 
ml Nextflow
nextflow pull nf-core/methylseq

# Nextflow parameters
export VERSION=1.6.1
export NXF_HOME=${PROJDIR}
export PATH=${NXF_HOME}:${PATH}
export NXF_TEMP=$SNIC_TMP
export NXF_LAUNCHER=$SNIC_TMP
export NXF_OPTS='-Xms1g -Xmx4g'

nextflow run nf-core/methylseq -r $VERSION \
  -profile uppmax \
  --project $PROJECT \
  --genome GRCh38 \
  --em_seq \
  --input 'merged/*_R{1,2}.fastq.gz' \
  --outdir 'results' \
  --aligner bismark \
  --email $EMAIL \
  -resume
```

## Complementary files

To elucidate the relevance of identified DMRs information about overlapping genomic features is needed (such as promoters, exon, intron, CpG-island). The database hosted at [University of California Santa Cruz (UCSC) Genomics Institute](https://genome-euro.ucsc.edu/index.html) holds genomic features for various species. Annotations can be exported from the [Table Browser tool](https://genome-euro.ucsc.edu/cgi-bin/hgTables). To download the annotations set the parameters as follow:

+ CpG-islands annotations, save as `cpgislands_GRCh38.bed`

```
clade = "Mammals" 
genome = "Human"
assembly = "Dec. 2013 (GRCh38/hg38)"
group = "Regulation"
track = "CpG-islands"
table = "cpgIslandExt"
output format = "BED"
output filename = "cpgislands_GRCh38.bed"

<click> "get output" 
<click> "get BED"
```

+ [Refseq](https://en.wikipedia.org/wiki/RefSeq) annotations, save as `refseq_UCSC_GRCh38.bed`

```
clade = "Mammals" 
genome = "Human"
assembly = "Dec. 2013 (GRCh38/hg38)"
group = "Genes and Gene Predictions"
track = "NCBI RefSeq"
table = "UCSC RefSeq (refGene)"
output format = "BED"
output filename = "refseq_UCSC_GRCh38.bed"

<click> "get output" 
<click> "get BED"
```

+ ensembl database and reference genome CpG-site positions, run `code/complementary_files.R`, will generate `ensembl_dataset_GRCm39.csv.gz` and `cg_pos_CRGh38.csv.gz`

```sh
Rscript code/complementary_files.R
```

The `GRCh38/` folder should contain the following: 

```
GRCh38/
├── cg_pos_CRGh38.csv.gz
├── cpgislands_GRCh38.bed
├── ensembl_dataset_GRCh38.csv.gz
└── refseq_UCSC_GRCh38.bed
```

## Differental methylation analysis

1. Differentially methylated regions (DMRs) were divided into 2 resolutions, (1) CpG and (2) tile of 100 bp. CpG-sites with low coverage (< 10 reads) and the top 99th percentile (PCR duplicates) were removed. Normalisation was done with scaling factor between samples based on differences between median of coverage distribution. The tiles were only considored if 2 or more CpG-sites where present. Finally, on a group level (control and exposed) CpG-sites were considored if 66% of the samples (4 out of 6) had coverage. Standard deviation (SD) filtering was applied where CpG-site with < 2 SD (little to no variation) were removed as they would not contribute information for downstream analysis.

```sh
# Ran at HPC (UPPMAX)
sbatch code/diffmeth.sh

## will start code/diffmeth.R with different arguments for CpG resolution
```

The `dump/` and `data/` folder should contain the following:

```
dump/
├── diffmeth_1_cpg.csv.gz
└── diffmeth_1_tile100.csv.gz

data/
├── PFOS_MCF-10A_betavalues_matrix_cpg.Rds
└── PFOS_MCF-10A_betavalues_matrix_tile100.Rds
```

2. Generate 3 tables: DMR and DMG, (1) DMR = each row is a dmr_id, (2) DMG = each row is gene with info about DMRs within it, genomic regions, dmr_id, hyper/hypo etc, (3) CGI = each row CpG-island. Significance threshold for DMRs were set to qvalue < 0.05 and meth.diff > ±15 and ±5, CpG-sites and 100 bp tiles, respectively. 

```sh
Rscript code/methtable.R
Rscript code/genetable.R
Rscript code/go_analysis.R
```

The `data/` folder should contain the following:

```
data/
├── PFOS_MCF-10A_DMG.Rds
├── PFOS_MCF-10A_DMR.Rds
└── PFOS_MCF-10A_GO.Rds
```
