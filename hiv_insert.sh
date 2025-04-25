###################################################################
## PROGRAM: hiv_insert.sh                                        ##
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
##    hiv_insert.sh <INPUT_DIR> <OUTPUT_DIR> <VIRAL_REFERENCE> <HOST_REFERENCE> ##
##                                                                              ##
## In order to set up your analysis, you should place your sequences of         ##
## interest in the input directory as either bam files or fastq files           ##
## ending with a fastq or bam suffix.  All files ending in either fastq         ##
## or bam in the specified input directory will be used in the analysis.        ##
## The HIV/SIV/viral reference should be in fasta format, but can contain       ##
## multiple references.  The HOST reference should be in fasta format as        ##
## well.                                                                        ##
################################################################################## 

############################################################
## Set the global variables based on the input parameters ##
############################################################
if [ "$#" -ne 4 ]; then
   echo "USAGE: hiv_insert.sh <INPUT_DIRECTORY> <OUTPUT_DIRECTORY> <HIV_REFERENCE> <HOST_REFERENCE>"
   exit 1
fi
IN_DIR=$1
OUT_DIR=$2
HIV_REF=$3
HOST_REF=$4

#######################################
## Check that input directory exists ##
#######################################
if [ ! -d "$IN_DIR" ]; then
   echo "Input directory $IN_DIR not found"
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

##################################################
## Find all of the files ending in fastq or bam ##
##################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILES=( $(find "$IN_DIR" | egrep "fastq$|bam$" ) )
if [ "${#FILES[@]}" -eq 0 ]; then
   echo "No bam or fastq files found in input directory $IN_DIR -- files must end with fastq or bam extension"
   exit 5
fi

#####################################################
## All parameters seem ok -- proceed with analysis ##
#####################################################
mkdir -p $OUT_DIR

#########################################################
## Process each of the fastq or bam files individually ##
#########################################################
for FILE in ${FILES[@]}
do
    echo $FILE
    bash $SCRIPT_DIR/hiv_insert_single.sh $FILE $OUT_DIR $HIV_REF $HOST_REF
done

