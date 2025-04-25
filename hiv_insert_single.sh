###################################################################
## PROGRAM: hiv_insert_single.sh                                 ##
## AUTHOR:  David Sachs, Ichan School of Medicine at Mount Sinai ##
##          Modified by Eric Rouchka, University of Louisville   ##
##          Copyright, University of Louisville                  ##
## DATE:    4/23/2025                                            ##
## LAST MODIFIIED: 4/23/2025                                     ##
## GITHUB:                                                       ##
## LICENSE:                                                      ##
###################################################################

##################################################################################
## USAGE:                                                                       ##
##    hiv_insert_single.sh <INPUT_FILE> <OUTPUT_DIR> <VIRAL_REFERENCE> <HOST_REFERENCE> ##
##                                                                              ##
## In order to set up your analysis, you should place your sequences of         ##
## interest in the input directory as either bam files or fastq files           ##
## ending with a fastq or bam suffix.  All files ending in either fastq         ##
## or bam in the specified input directory will be used in the analysis.        ##
## The HIV/SIV/viral reference should be in fasta format, but can contain       ##
## multiple references.  The HOST reference should be in fasta format as        ##
## well.                                                                        ##
##################################################################################

################################
## Test for required programs ##
################################
if !(command zcat --version > /dev/null); then
   echo "ERROR: zcat is required, but is not found in the path"
   exit 5
fi
if !(command samtools --version > /dev/null); then
   echo "ERROR: samtools is required, but is not found in the path"
   exit 5
fi
if !(which seqtk > /dev/null); then
   echo "ERROR: seqtk is required, but is not found in the path"
   exit 5
fi
if !(command minimap2 --version > /dev/null); then
   echo "ERROR: minimap2 is required, but is not found in the path"
   exit 5
fi
if !(command python --version > /dev/null); then
   echo "ERROR: python is required, but is not found in the path"
   exit 5
fi
if !(command bedtools --version > /dev/null); then
   echo "ERROR: bedtools is required, but is not found in the path"
   exit 5
fi
exit 0



############################################################
## Set the global variables based on the input parameters ##
############################################################

if [ "$#" -ne 4 ]; then
   echo "USAGE: hiv_insert_single.sh <INPUT_FILE> <OUTPUT_DIRECTORY> <HIV_REFERENCE> <HOST_REFERENCE>"
   exit 1
fi
INITIAL_FILE=$1
OUT_DIR=$2
HIV_REF=$3
HOST_REF=$4

#######################################
## Check that input directory exists ##
#######################################
if [ ! -f "$INITIAL_FILE" ]; then
   echo "Input file $INITIAL_FILE not found"
   exit 2
fi

########################################
## Check that output directory exists ##
########################################
if [ ! -d "$OUT_DIR" ]; then
   echo "Output directory $OUT_DIR not found"
   exit 2
fi

##########################################
## Check that HIV reference file exists ##
##########################################
if [ ! -f "$HIV_REF" ]; then
   echo "HIV/SIV/viral Reference file $HIV_REF does not exist -- reference should be a fasta file (may contain multiple references)"
   exit 3
fi

###########################################
## Check that Host reference file exists ##
###########################################
if [ ! -f "$HOST_REF" ]; then
   echo "Host Reference file $HOST_REF does not exist -- reference should be a fasta file (may contain multiple references)"
   exit 4
fi


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILE_EXT=$(echo $INITIAL_FILE | tr "." "\n" | tail -1)

######################################################
## Convert the input file to an unzipped fastq file ##
## which may require:                               ##
##  (1) unzipping the file and                      ##
##  (2) converting bam file to fastq using samtools ##
######################################################
if [ $FILE_EXT == "gz" ]; then
    FASTQ=$OUT_DIR/$(basename $INITIAL_FILE | sed 's/.gz$//g')
    echo $FASTQ
    zcat $INITIAL_FILE > $FASTQ
elif [ $FILE_EXT == "bam" ]; then
    FASTQ=$OUT_DIR/$(basename $INITIAL_FILE | sed 's/.bam$/.fastq/g')
    echo "samtools fastq $INITIAL_FILE > $FASTQ"
    samtools fastq $INITIAL_FILE > $FASTQ	
else
    FASTQ=$INITIAL_FILE
fi

PREFIX=$OUT_DIR/$(basename $FASTQ ).hiv
echo $FASTQ 

#########################################################
## Trim adapters -- this is likely an unnecessary step ##
## and is currently set to trim 0 from beginning and   ##
## o from the end                                      ##
#########################################################
seqtk trimfq -b 0 -e 0 $FASTQ > $PREFIX.trimmed.fastq

################################################################
## Step 1: align the raw fastq files to the viral reference   ##
##         genome using 16 threads (-t 16) a minimal chaining ##
##         score of 0 to allow for more matches (-m 0), soft  ##
##         clipping for supplemental alignments (-Y) and      ##
##         alignment of PacBio reads (-ax map-pb)             ##
################################################################
minimap2 -t  16 -m 0 -Y -ax map-pb $HIV_REF $PREFIX.trimmed.fastq > $PREFIX.initial.1.sam

##############################################################################
## Loop HIV alignment until no more HIV sequences can be found in the reads ##
## this will include the header (-h) use a sam file (-S), filter out        ##
## secondary and unmapped alignments (-F 260) and remove  secondary         ##
## alignments, unmapped sequences, and supplementary alignments (0x904)     ##
##############################################################################
NUM=1
samtools view -h -S -F 260 $PREFIX.initial.$NUM.sam > $PREFIX.$NUM.sam
samtools view -h -S -F 0x904 $PREFIX.initial.1.sam > $PREFIX.1.sam

while [[ $(samtools view -c -F 4 $PREFIX.$NUM.sam) -ne 0 ]]
do
  python $SCRIPT_DIR/mask.py $PREFIX.$NUM.sam > $PREFIX.$NUM.fa
  minimap2 -t  16 -Y -p 0 -N 10000 -ax map-pb $HIV_REF $PREFIX.$NUM.fa > $PREFIX.initial.$((NUM+1)).sam
  python $SCRIPT_DIR/pick_reads.py $PREFIX.$NUM.sam $PREFIX.initial.$((NUM+1)).sam $PREFIX.$((NUM+1)).sam
  NUM=$((NUM+1))
  echo $NUM
  if [[ $(samtools view -c -F 4 $PREFIX.initial.$NUM.sam) = 0 ]]; then echo "Done"; fi;
done

python $SCRIPT_DIR/mask.py $PREFIX.$NUM.sam > $PREFIX.$NUM.fa

#Reverse mask to provide sequences containing only HIV aligned segments
python $SCRIPT_DIR/unmask.py $PREFIX.1.sam $PREFIX.$NUM.fa > $PREFIX.unmasked.fa

#Extract flanks from sequences and map to human genome
python $SCRIPT_DIR/get_flanks.py $PREFIX.$((NUM-1)).fa > $PREFIX.flanks.fa
minimap2 -t  16 -Y -p 0 -N 10000 -ax map-pb $HOST_REF $PREFIX.flanks.fa > $PREFIX.flanks.sam

#Confirm alignments by mapping original reads to human genome
samtools view -bS $PREFIX.1.sam > $PREFIX.1.bam
bedtools bamtofastq -i $PREFIX.1.bam -fq $PREFIX.1.fastq
samtools fastq $PREFIX.1.sam > $PREFIX.1.fastq 
minimap2 -t  16 -Y -ax map-pb $HOST_REF $PREFIX.1.fastq > $PREFIX.human.sam
samtools view -h -S -F 0x900 $PREFIX.human.sam > $PREFIX.human.filtered.sam
#python $SCRIPT_DIR/mask.py $OUT_DIR
python $SCRIPT_DIR/mask.py $PREFIX.human.filtered.sam

#Generate output files
#python $SCRIPT_DIR/combine_hiv.py $PREFIX > $PREFIX.combine.log
python $SCRIPT_DIR/combine_hiv_V2.py $PREFIX > $PREFIX.combine.log
