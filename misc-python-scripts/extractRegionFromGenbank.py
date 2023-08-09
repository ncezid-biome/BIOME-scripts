# Joseph S. Wirth
# July 2023

import getopt, os, sys
from io import TextIOWrapper
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def __extractRecord(inFN:str, contigs:list, start:int, end:int) -> list[SeqRecord]:
    """ retrieves a region of an input genbank and creates a new SeqRecord

    Args:
        inFN (str): input file (genbank file format)
        contigs (list): the names of the desired contigs
        start (int): the start of the desired region
        end (int): the end of the desired region

    Raises:
        ValueError: the input genbank is not readable
        ValueError: the contig does not exist

    Returns:
        list[SeqRecord]: the requested contigs
    """
    # constants
    GBK = "genbank"
    ERR_MSG_1 = "could not read the input file"
    ERR_MSG_2 = "could not find one or more of the specified contigs"
    
    # parse the input file
    parsed = SeqIO.parse(inFN, GBK)
    
    # go through each record until the specified contig if found
    num = 0
    rec:SeqRecord
    outL = list()
    for rec in parsed:
        num += 1
        # return the desired region
        if rec.id in contigs:
            # return contig to the very end if no end specified
            if end is None:
                outL.append(rec[start-1:])

            # otherwise return the specified region
            else:
                outL.append(rec[start-1:end])
    
    # make sure the file was properly formatted
    if num == 0:
        raise ValueError(ERR_MSG_1)
    
    # the contig was not found if at this line
    if len(contigs) != len(outL):
        raise ValueError(ERR_MSG_2)

    return outL


def __writeRecord(rec:list[SeqRecord], fh:TextIOWrapper, format:str) -> None:
    """writes a SeqRecord to file

    Args:
        rec (list): the SeqRecord to write to file
        fh (TextIOWrapper): a filehandle to where the record should be written
        format (str): the desired output format
    """
    SeqIO.write(rec, fh, format)


def __parseArgs() -> tuple[str,list[str],int,int,TextIOWrapper,str,bool]:
    """ parses command line arguments

    Raises:
        FileNotFoundError: input file needs to exist
        ValueError: output format must be genbank or fasta

    Returns:
        tuple[str,str,int,int,TextIOWrapper,str]: returns the parsed arguments:
            input file, contig, start position, end position, out filehandle, out format
    """
    # flags
    INPUT_FLAGS = ("-i","--in")
    CONTIG_FLAGS = ("-c","--contig")
    START_FLAGS = ("-s","--start")
    END_FLAGS = ("-e","--end")
    OUT_FLAGS = ("-o","--out")
    FORMAT_FLAGS = ("-f","--format")
    HELP_FLAGS = ("-h","--help")
    SHORT_OPTS = INPUT_FLAGS[0][-1] + ":" + \
                 CONTIG_FLAGS[0][-1] + ":" + \
                 START_FLAGS[0][-1] + ":" + \
                 END_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 FORMAT_FLAGS[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (INPUT_FLAGS[1][2:] + "=",
                 CONTIG_FLAGS[1][2:] + "=",
                 START_FLAGS[1][2:] + "=",
                 END_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 FORMAT_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:])
    
    # default values
    DEFAULT_START = 1
    DEFAULT_END = None
    DEFAULT_OUT = sys.stdout
    DEFAULT_FORMAT = "genbank"
    ALLOWED_FORMATS = ("fasta", "genbank")
    
    # messages
    IGNORE_MSG = "ignoring unused argument "
    NO_COORDS_ALLOWED = "ignoring start/end coordinates; incompatible with multiple contigs"
    ERR_MSG_1 = "could not find specified input file"
    ERR_MSG_2 = "output file already exists"
    ERR_MSG_3 = "could not open output file"
    ERR_MSG_4 = "unrecognized out format"
    ERR_MSG_5 = "input file not specified"
    ERR_MSG_6 = "contig name(s) not specified"
    
    # helper function for printing help messages
    def helpMessage() -> None:
        """ prints the help message
        """
        GAP = 4*" "
        EOL = "\n"
        SEP = ", "
        MSG = EOL + "Extracts a region from a genbank file" + EOL + \
             GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
             "usage:" +  EOL + \
             GAP + "python3 extractRegionFromGenbank.py [-icseofh]" + EOL*2 + \
             "required arguments:" + EOL + \
             GAP + f'{INPUT_FLAGS[0] + SEP + INPUT_FLAGS[1]:<16}{"[str] the filename of an input genbank file"}' + EOL +\
             GAP + f'{CONTIG_FLAGS[0] + SEP + CONTIG_FLAGS[1]:<16}{"[str] the name of the contig(s) (comma-separated) within the input file"}' + EOL*2 + \
             "optional arguments:" + EOL + \
             GAP + f'{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<16}{"[str] the output filename (default: stdout)"}' + EOL + \
             GAP + f'{START_FLAGS[0] + SEP + START_FLAGS[1]:<16}{"[int] the start position (default: contig start)"}' + EOL + \
             GAP + f'{END_FLAGS[0] + SEP + END_FLAGS[1]:<16}{"[int] the end position (default: contig end)"}' + EOL + \
             GAP + f'{FORMAT_FLAGS[0] + SEP + FORMAT_FLAGS[1]:<16}{"[str] the output format [fasta|genbank] (default: genbank)"}' + EOL + \
             GAP + f'{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<16}{"print this message"}' + EOL
            
        print(MSG)
    
    # set default values
    inFN = None
    contigs = None
    start = DEFAULT_START
    end = DEFAULT_END
    outFH = DEFAULT_OUT
    outfmt = DEFAULT_FORMAT
    helpRequested = False
    
    # print the help message then exit if help was requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpMessage()
        helpRequested = True
    
    else:
        # extract command line arguments
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        
        # parse each command line argument
        for opt,arg in opts:
            # get the input filename
            if opt in INPUT_FLAGS:
                if not os.path.exists(arg):
                    raise FileNotFoundError(ERR_MSG_1)
                inFN = arg
            
            # get the contig name(s)
            elif opt in CONTIG_FLAGS:
                contigs = arg.split(",")
            
            # get the start position
            elif opt in START_FLAGS:
                start = int(arg)
            
            # get the end position
            elif opt in END_FLAGS:
                end = int(arg)
            
            # get the output filename
            elif opt in OUT_FLAGS:
                if os.path.exists(arg):
                    raise FileExistsError(ERR_MSG_2)
                try:
                    outFH = open(arg, 'w')
                except:
                    raise IOError(ERR_MSG_3)
            
            # get the output file format
            elif opt in FORMAT_FLAGS:
                if arg not in ALLOWED_FORMATS:
                    raise ValueError(ERR_MSG_4)
                outfmt = arg
        
            # ignore any other arguments
            else:
                print(IGNORE_MSG + opt)
        
        # make sure all required arguments have been provided
        if inFN is None:
            raise ValueError(ERR_MSG_5)
        if contigs is None:
            raise ValueError(ERR_MSG_6)
        
        # check if the positional coordinates have been modified
        coordsModified = start != DEFAULT_START or end != DEFAULT_END
        
        # cannot use positional coordinates if extracting multiple contigs
        if len(contigs) > 1 and coordsModified:
            start = DEFAULT_START
            end = DEFAULT_END
            print(NO_COORDS_ALLOWED)
    
    return inFN, contigs, start, end, outFH, outfmt, helpRequested


def __main():
    """ main runner function
    """
    # parse the arguments
    inFN, contigs, start, end, outFH, outfmt, helpRequested = __parseArgs()
    
    # inFN will be None if a help flag was specified
    if not helpRequested:
        recs = __extractRecord(inFN, contigs, start, end)
        __writeRecord(recs, outFH, outfmt)
        if outFH != sys.stdout:
            outFH.close()


if __name__ == "__main__":
    __main()
