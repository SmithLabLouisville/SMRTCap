# SMRTcap
SMRTcap is a PacBio long read sequencing protocol for analyzing viral integrations into a host genome. See: 


Sadri G, Nadakal ST, Sachs D, Lauer W, Kos J, Singh P, Elliott EM, Kaiser CM, Ford EE, Richardson N, Hudson E, Linden N, Powell J, Warburton P, Soto J, Emery M, Deikus G, Lee GQ, Lamers S, Reynolds SJ, Galiwango RM, Prodger JL, Tomusange S, Kityamuweesi T, Han T, Jones RB, Tobian AAR, Engelman A, Sebra R, Morgello S, Redd AD, Rouchka E, Smith ML. (2025) **Development and validation of HIV SMRTcap for the characterization of HIV reservoirs across tissues and subtypes**. (under review)

# SMRTcap analysis pipeline
The SMRTcap pipeline is a combination of bash, R, perl, and python scripts designed to analyze long sequencing reads for viral integrations into host genomes.  There are three main driver programs: 1. hiv_insert.sh which finds integration sites and flanking host sequences; 2. BMS.Insertion.v3.3.ECR.R which annotates the integrated virus according to where in the host genome (intergenic or specific genic regions); and 3. findHIVSIVGeneREgionsV9.5 which annotates the completeness in terms of each of the viral genes for the integration.

## Requirements
- samtools    (tested on v1.16.1)
- seqtk       (tested on v1.3-r117-dirty)
- minimap2    (tested on 2.24-v1122)
- python      (tested on 2.7.5)
- bedtools    (tested on 2.30.0)
- ncbi-blast  (tested on 2.10.0+)
- perl        (tested on v5.16.3)

## Required Perl libraries
- Getopt::long

## Required Python libraries
- pysam
- sys
- re
- os.listdir
- os.path.isfile
- os.path.join
- xlsxwriter
- warnings

## Required R libraries
- GenomicRanges
- rtracklayer
- TxDb.Hsapiens.UCSC.hg38.knownGene
- ensembldb
- EnsDb.Hsapiens.v86
- annotatr
- AnnotationHub
- ggplot2
- dplyr
- stringr
- readxml
- ggforce
- stringi
- chromoMap
- randomcoloR
- data.table
- patchwork
- S4Vectors
- tidyr
