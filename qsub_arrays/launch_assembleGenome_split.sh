#!/bin/bash -l

#-------------------------------------------------------------------------------------------------
# This is a qsub array job script whichs downloads a large batch of SRA data.
# The ordered columns of the required CSV import are: Key,WGS_id,NCBI_ACCESSION,SRR_id. The second
# and fourth columns are the only columns which require content; the others can be blank.
# Output files will have the named prefix structure of WGSID_SRRID.
# Please adjust the user id of line number #69 to match your user id.
#-------------------------------------------------------------------------------------------------

# Read ARGV
INPUT_CSV=$1
SLOTS_PER_JOB=1 # manually change this as needed
PARENTDIR=$2

if [ "INPUT_CSV" == "" ]; then
  echo "run shovill - trims, assembles, runs read correction; formats names, sets cov cutoff."
  echo "Set the cov cutoff to coverage calculated by CG pipeline divided by ten whichever is smaller since shovill downsamples to 100X."
  echo "1 = Read 1 , 2 = Read 2, 3 = Strain, 4 = min coverage"
  echo "Usage: $0 [STRAIN-LIST].csv [OUTPUT-DIR]"
  exit 1;
fi

mkdir $PARENTDIR

TMP=$(mktemp --tmpdir='.' --directory qsubShovill.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

CTRL_FILE="$TMP/array.txt"
cp $INPUT_CSV $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"


qsub -N Shovill -q all.q -o $TMP/log -j y -pe smp $SLOTS_PER_JOB -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE","PARENTDIR=$PARENTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  source /etc/profile.d/modules.sh

  which shovill

  # Set up filenames
  STRAIN=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' {print $3}')
  MINCOV=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' {print $4}')
  READ1=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' {print $1}')
  READ2=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk -F',' {print $2}')
  tmpdir=/scratch/$USER
  mkdir -p $tmpdir
  outdir=$(mktemp --tmpdir=$tmpdir --directory Shovill.XXXXXX);
  trap "rm -rf $outdir" EXIT

  echo "shovill will be run on $STRAIN with minimum coverage of $MINCOV under $(hostname)"
  echo "Working directory is $outdir"

  echo shovill --outdir "$outdir"/shovill --R1 $READ1 --R2 $READ2 --mincov $MINCOV --trim --namefmt "$STRAIN"_contig%05d
  exit 0  

  #rename files to include strain name
  mv $outdir/shovill $PARENTDIR/$STRAIN
  #mv $PARENTDIR/$STRAIN/contigs.fa $PARENTDIR/$STRAIN/$STRAIN.contigs.fa; mv $PARENTDIR/$STRAIN/contigs.gfa $PARENTDIR/$STRAIN/$STRAIN.contigs.gfa;
  #mv $PARENTDIR/$STRAIN/shovill.corrections $PARENTDIR/$STRAIN/$STRAIN.shovill.corrections; mv $PARENTDIR/$STRAIN/shovill.log $PARENTDIR/$STRAIN/$STRAIN.shovill.log; mv $PARENTDIR/$STRAIN/spades.fasta $PARENTDIR/$STRAIN/$STRAIN.spades.fasta
  
  
END_OF_SCRIPT

