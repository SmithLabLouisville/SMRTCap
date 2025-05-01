# SMRTCap
SMRTCap analysis pipeline

## Requirements
samtools    (tested on v1.16.1)
seqtk       (tested on v1.3-r117-dirty)
minimap2    (tested on 2.24-v1122)
python      (tested on 2.7.5)
bedtools    (tested on 2.30.0)
ncbi-blast  (tested on 2.10.0+)
perl        (tested on v5.16.3)

## Required Perl libraries
Getopt::long

## Required Python libraries
pysam
sys
re
os.listdir
os.path.isfile
os.path.join
xlsxwriter
warnings

## Required R libraries
GenomicRanges
rtracklayer
TxDb.Hsapiens.UCSC.hg38.knownGene
ensembldb
EnsDb.Hsapiens.v86
annotatr
AnnotationHub
ggplot2
dplyr
stringr
readxml
ggforce
stringi
chromoMap
randomcoloR
data.table
patchwork
S4Vectors
tidyr
