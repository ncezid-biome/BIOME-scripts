#!/bin/bash -l

#-------------------------------------------------------------------------------------------------
# This is a qsub array job script whichs assembles a large batch of raw read data with Shovill.
# The ordered columns of the required CSV import are: READ1_NAME,READ2_NAME,STRAIN_ID,MIN_COVERAGE.
# Output files will have the named prefix structure of STRAIN_ID.
#-------------------------------------------------------------------------------------------------

VERSION="0.2";

# Read ARGV
INPUT_CSV=$1
SLOTS_PER_JOB=1 # manually change this as needed
PARENTDIR=$2
READDIR=$3

# If there is not enough arguments, exit and display usage
if [ "$READDIR" == "" ]; then
  echo "run shovill - trims, assembles, runs read correction; formats names, sets cov cutoff."
  echo "Set the cov cutoff to coverage calculated by CG pipeline divided by ten whichever is smaller since shovill downsamples to 100X."
  echo "The ordered columns of the required CSV import are: READ1_NAME,READ2_NAME,STRAIN_ID,MIN_COVERAGE."
  echo "The [READDIR] should contain both R1 and R2 fastq.gz files of STRAIN_IDs in CSV."
  echo "Usage: $0 [STRAIN-LIST].csv [OUTPUT-DIR] [READDIR]"
  exit 1;
fi

mkdir $PARENTDIR

# Create working temporary directory
TMP=$(mktemp --tmpdir='.' --directory qsubShovill.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# Set up array control structure
CTRL_FILE="$TMP/array.txt"
cp $INPUT_CSV $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"


qsub -N Shovill -q all.q -o $TMP/log -j y -pe smp $SLOTS_PER_JOB -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE","PARENTDIR=$PARENTDIR","READDIR=$READDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  source /etc/profile.d/modules.sh
  module load bwa/0.7.17 flash kmc lighter
  module load megahit pilon Skesa/2.3.0 SPAdes/3.14.0
  module load trimmomatic velvet/1.2.10 seqtk/1.3

  which shovill
  shovill --check

  # Retrieve column values from INPUT_CSV
  STRAIN=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' '{print $3}')
  MINCOV=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' '{print $4}')
  READ1=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' '{print $1}')
  READ2=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' '{print $2}')
  tmpdir=/scratch/$USER
  mkdir -p $tmpdir
  outdir=$(mktemp --tmpdir=$tmpdir --directory Shovill.XXXXXX);
  trap "rm -rf $outdir" EXIT

  echo "Shovill will be run on $STRAIN with minimum coverage of $MINCOV under $(hostname)"
  echo "Working directory is $outdir"

  shovill --outdir "$outdir"/shovill --R1 $READDIR/$READ1 --R2 $READDIR/$READ2 --mincov $MINCOV --trim --namefmt "$STRAIN"_contig%05d 

  # Rename Shovill output files to include strain name
  mv $outdir/shovill $PARENTDIR/$STRAIN
  mv $PARENTDIR/$STRAIN/contigs.fa $PARENTDIR/$STRAIN/$STRAIN.contigs.fa;
  mv $PARENTDIR/$STRAIN/contigs.gfa $PARENTDIR/$STRAIN/$STRAIN.contigs.gfa;
  mv $PARENTDIR/$STRAIN/shovill.corrections $PARENTDIR/$STRAIN/$STRAIN.shovill.corrections;
  mv $PARENTDIR/$STRAIN/shovill.log $PARENTDIR/$STRAIN/$STRAIN.shovill.log;
  mv $PARENTDIR/$STRAIN/spades.fasta $PARENTDIR/$STRAIN/$STRAIN.spades.fasta
  
  
END_OF_SCRIPT

