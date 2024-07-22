# README

The code in the repository was used to for computational analysis for the project: [Perfluorooctanesulfonic acid (PFOS) induced cancer related methylome alterations in human breast cells](doi)

## Project summary

+ Dataset is composed of Perflourooctanesulfonic acid (PFOS) exposed (72h) normal human breast cell line (MCF-10A). 12 samples, 6 control, 6 exposed.
+ The sequencing data is built upon the study: [PFOS induces proliferation, cell-cycle progression, and malignant phenotype in human breast epithelial cells](https://doi.org/10.1007/s00204-017-2077-8). From this study only the control and 1µM cells cultures were used.
+ As PFOS has potentially [endocrine disrupting estrogen activity](https://doi.org/10.1111/j.1365-2605.2008.00870.x) and has been [linked to breast cancer](https://doi.org/10.1186/1476-069X-10-88), the hypothesis is to identify DM regions which overlap genes related to the phenotypes listed above.
+ DNA methylation levels were mapped with Enzymatic Methyl sequencing (EM-seq)

**Aim** to elucidate methylation alteration induced by PFOS exposure, aka "methylome fingerprint"

## Workflow

### File structure

```
project/
├── code
├── data
├── dump
├── merged
├── results
└── README.md
```

### Mapping & coverage

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

### Complementary files

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
├── ensembl_dataset_GRCm39.csv.gz
└── refseq_UCSC_GRCh38.bed
```

### Differental methylation analysis
