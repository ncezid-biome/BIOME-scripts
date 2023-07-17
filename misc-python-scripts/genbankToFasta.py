# Joseph S. Wirth
# July 2023

import getopt, os, sys
from io import TextIOWrapper
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


def __parseArgs() -> tuple[str,TextIOWrapper,TextIOWrapper]:
    """parses command line arguments

    Raises:
        ValueError: check allowed out formats
        FileNotFoundError: input file must exist
        FileExistsError: output file must not exist
        SyntaxError: input file and out format are required

    Returns:
        tuple[str,TextIOWrapper,TextIOWrapper]: out format, input filehandle, output filehandle
    """
    # constants
    SHORT_OPTS = ("i:f:o:h")
    LONG_OPTS = ("in=",
                 "format=",
                 "out=",
                 "help")
    FMT_FLAGS = ("-f", "--format")
    IN_FLAGS = ("-i", "--in")
    OUT_FLAGS = ("-o", "--out")
    HELP_FLAGS = ("-h", "--help")
    ALLOWED_FORMATS = ("nuc", "prot")
    ERR_MSG_1 = "invalid format specified"
    ERR_MSG_2 = "input genbank file not found"
    ERR_MSG_3 = "output file already exists"
    ERR_MSG_4 = "please provide an input file and a output format"
    IGNORE_MSG = "ignoring unused argument "

    # initialize variables
    outfmt = None
    inFH = None
    outFH = None
    
    # print the help message if requested then return
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv:
        __helpMessage()
        return outfmt, inFH, outFH

    # parse each command line argument
    opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
    for opt, arg in opts:
        # format
        if opt in FMT_FLAGS:
            if arg not in ALLOWED_FORMATS:
                raise ValueError(ERR_MSG_1)
            outfmt = arg
        
        # input file
        elif opt in IN_FLAGS:
            if not os.path.exists(arg):
                raise FileNotFoundError(ERR_MSG_2)
            inFH = open(arg, 'r')
    
        # output file
        elif opt in OUT_FLAGS:
            if os.path.exists(arg):
                raise FileExistsError(ERR_MSG_3)
            outFH = open(arg, 'w')
    
        # ignore any additional arguments
        else:
            print(IGNORE_MSG + opt)
        
    # make sure all required arguments have been provided
    if None in (inFH, outfmt):
        raise SyntaxError(ERR_MSG_4)
    
    # get the default out filehandle if one was not provided
    if outFH is None:
        outFH = __getOutFileHandle(inFH.name, outfmt)
    
    return outfmt, inFH, outFH


def __helpMessage() -> None:
    """prints the help message to the screen
    """
    GAP = " "*4
    HELP_MSG = "\nConverts a genbank file to a nucleotide (fna) or amino acid (faa) fasta\n" + \
                GAP + "Joseph S. Wirth, 2023\n\n" + \
                "Usage:\n" + GAP + "genbankToFasta.py [-ifoh]\n\n" + \
                "Required arguments:\n" + \
                GAP + "-i,--in        input file in genbank file format\n" + \
                GAP + "-f,--format    output file format {nuc,prot}\n\n" + \
                "Optional arguments:\n" + \
                GAP + "-o,--out       output filename (default: input filename with new extension)\n" + \
                GAP + "-h,--help      print this message\n"
    
    print(HELP_MSG)


def __getOutFileHandle(inFN:str, outfmt:str) -> TextIOWrapper:
    """gets a filehandle to an output file if one was not provided

    Args:
        inFN (str): the input filename
        outfmt (str): the output format

    Raises:
        FileExistsError: output file cannot already exist

    Returns:
         TextIOWrapper: filehandle for the output file
    """
    # constants
    FNA_FMT = "nuc"
    FNA_EXT = ".fna"
    FAA_FMT = "prot"
    FAA_EXT = ".faa"
    OUT_MSG = "file will be written to "
    ERR_MSG = "output file already exists: "

    # get basename
    outFN = os.path.splitext(inFN)[0]
    
    # get the nucleotide filename
    if outfmt == FNA_FMT:
        outFN += FNA_EXT

    # get the amino acid filename
    elif outfmt == FAA_FMT:
        outFN += FAA_EXT
    
    # print status
    print(OUT_MSG + outFN)
    
    # make sure file doesn't already exist
    if os.path.exists(outFN):
        raise FileExistsError(ERR_MSG + outFN)

    return open(outFN, 'w')


def __gbk2fna(gbkFH:TextIOWrapper, outFH:TextIOWrapper) -> None:
    """converts a genbank file to a nucleotide fasta

    Args:
        gbkFH (TextIOWrapper): input genbank filehandle
        outFH (TextIOWrapper): output fasta filehandle

    Raises:
        RuntimeError: input file needs to be properly formatted
    """
    # constants
    IN_FMT = "genbank"
    OUT_FMT = "fasta"
    ERR_MSG = "input file is not properly formatted"

    # parse the genbank file into a list of SeqRecords
    records = list(SeqIO.parse(gbkFH, IN_FMT))

    # make sure the file could be parsed
    if len(records) == 0:
        raise RuntimeError(ERR_MSG)

    # write records to file
    SeqIO.write(records, outFH, OUT_FMT)


def __gbk2faa(gbkFH:TextIOWrapper, outFH:TextIOWrapper) -> None:
    """converts a genbank file to an amino acid fasta

    Args:
        gbkFH (TextIOWrapper): input genbank filehandle
        outFH (TextIOWrapper): output fasta filehandle

    Raises:
        RuntimeError: input file needs to be properly formatted
    """
    # constants
    OUT_FMT = "fasta"
    ERR_MSG = "could not extract CDS records from input file"

    # retrieve a list of records; one for each CDS
    recordsL = __getCdsRecords(gbkFH)

    # check that records were successfully retrieved
    if len(recordsL) == 0:
        raise RuntimeError(ERR_MSG)

    # write records to file
    SeqIO.write(recordsL, outFH, OUT_FMT)


def __getCdsRecords(gbkFH:TextIOWrapper) -> list[SeqRecord]:
    """retrieves a list of CDS records from a genbank filehandle

    Args:
        gbkFH (TextIOWrapper): filehandle to an input genbank

    Returns:
        list[SeqRecord]: a list of CDS as SeqRecord objects
    """
    # constants
    FORMAT = "genbank"
    CDS = "CDS"
    TAG = 'locus_tag'
    GENE = "gene"
    ANNOTATION = "product"
    SEQ = 'translation'
    SEP = "|"

    # parse the genbank file into a list of SeqRecords
    parsed = SeqIO.parse(gbkFH, FORMAT)

    # build a list of CDS SeqRecords
    rec:SeqRecord
    feat:SeqFeature
    records = list()
    for rec in parsed:
        for feat in rec.features:
            # only process a feature if it is a CDS
            if feat.type == CDS:
                # get the locus tag
                tag = feat.qualifiers[TAG][0]

                # get the gene name
                if GENE in feat.qualifiers.keys():
                    gene = feat.qualifiers[GENE][0]
                else:
                    gene = ""
                
                # get the annotation
                ann = feat.qualifiers[ANNOTATION][0]
                
                # attempte to get the amino acid sequence
                if SEQ in feat.qualifiers.keys():
                    seq = Seq(feat.qualifiers[SEQ][0])
                else:
                    try:
                        seq = Seq(feat.translate(rec.seq))

                    # skip the record if AA sequence could not be extracted
                    except: continue

                # create a new SeqRecord object and add it to the list
                newRec = SeqRecord(seq)
                newRec.id = tag + SEP + gene + SEP + ann
                newRec.name = ""
                newRec.description = ""
                records.append(newRec)

    return records


def __main() -> None:
    """main runner function. parses command line arguments; converts genbank to fasta
    """
    # constants
    FNA = "nuc"
    FAA = "prot"

    # no arguments ==> help message
    if len(sys.argv) == 1:
        __helpMessage()
        return
    
    # parse command line arguments
    outfmt,inFH,outFH = __parseArgs()
    
    # if outfmt is None, then help was requested (nothing to do)
    if outfmt is not None:
        if outfmt == FNA:
            __gbk2fna(inFH, outFH)
        
        elif outfmt == FAA:
            __gbk2faa(inFH, outFH)


if __name__ == "__main__":
    __main()
