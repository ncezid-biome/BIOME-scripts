#!/usr/bin/env python3
__author__ = "Joseph S. Wirth"

import getopt,os,sys
from Bio import SeqIO
from Bio.Seq import Seq
from io import TextIOWrapper
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def __extractTags(tags:list[str], gbkFn:str, mol:str) -> list[SeqRecord]:
    """extracts locus tags from a genbank file

    Args:
        tags (list[str]): a list of locus tags to extract
        gbkFn (str): the genbank file to extract locus tags from
        mol (str): the molecule type (nucl|prot)

    Returns:
        list[SeqRecord]: a list of locus tags as SeqRecord objects
    """
    # constants
    FORMAT = 'genbank'
    TAG = "locus_tag"
    NUCL = "nucl"
    PROT = "prot"
    
    # initialize variables
    out = list()
    found = set()
    rec:SeqRecord
    feat:SeqFeature
    
    # go through each feature in the genbank file
    for rec in SeqIO.parse(gbkFn, FORMAT):
        for feat in rec.features:
            # check if the locus tag was requested
            if TAG in feat.qualifiers.keys():
                if feat.qualifiers[TAG][0] in tags and feat.qualifiers[TAG][0] not in found:
                    # get either the nucleotide or amino acid sequence
                    if mol == NUCL:
                        newRec = feat.extract(rec)
                    elif mol == PROT:
                        try:
                            newRec = SeqRecord(Seq(feat.qualifiers['translation'][0]))
                        except:
                            newRec = feat.translate(rec)
                            if type(newRec) == Seq:
                                newRec = SeqRecord(newRec)
                
                    # update the record data and save it
                    newRec.id = feat.qualifiers[TAG][0]
                    newRec.name = ''
                    newRec.description = ''
                    found.add(feat.qualifiers[TAG][0])
                    out.append(newRec)
    
    return out


def __parseArgs() -> tuple[str,list[str],str,TextIOWrapper,bool]:
    """parses command line arguments

    Raises:
        FileNotFoundError: input genbank is not a file
        ValueError: invalid molecule type specified
        RuntimeError: required argument(s) missing

    Returns:
        tuple[str,list[str],str,TextIOWrapper,bool]: gbk filename, locus tags, molecule type, filehandle, helpRequested
    """
    # flags
    GBK_FLAGS = ("-g", "--gbk")
    TAG_FLAGS = ("-t", "--tag")
    MOL_FLAGS = ("-m", "--mol")
    OUT_FLAGS = ("-o", "--out")
    HELP_FLAGS = ("-h", "--help")
    SHORT_OPTS = GBK_FLAGS[0][-1] + ":" + \
                 TAG_FLAGS[0][-1] + ":" + \
                 MOL_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (GBK_FLAGS[1][2:] + "=",
                 TAG_FLAGS[1][2:] + "=",
                 MOL_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:])
    
    # constants
    SEP = ","
    ALLOWED_MOL = ("nucl", "prot")
    
    # messages
    ERR_MSG_1 = "input file not found"
    ERR_MSG_2 = "invalid molecule type specified"
    ERR_MSG_3 = "must specify all required arguments"
    IGNORE_MSG = "ignoring unused argument: "
    
    # default values
    DEFAULT_MOL = "nucl"
    DEFAULT_OUT = sys.stdout
    
    # define a helper function for printing the help message
    def help() -> None:
        GAP = 4*" "
        EOL = "\n"
        MSG = EOL + "Extracts locus tag sequences from a genbank file" + EOL + \
              GAP + __author__ + ", 2023" + EOL*2 + \
              "usage:" + EOL + \
              GAP + "python extractLocusTagsFromGenbank.py [-" + GBK_FLAGS[0][-1] + \
                TAG_FLAGS[0][-1] + MOL_FLAGS[0][-1] + OUT_FLAGS[0][-1] + \
                HELP_FLAGS[0][-1] + "]" + EOL*2 + \
              "required arguments:" + EOL + \
              GAP + f"{GBK_FLAGS[0] + SEP + GBK_FLAGS[1]:<15}{'a genbank file':<}" + EOL + \
              GAP + f"{TAG_FLAGS[0] + SEP + TAG_FLAGS[1]:<15}{'a comma-separated list of locus tag(s)':<}" + EOL*2 + \
              "optional arguments:" + EOL + \
              GAP + f"{MOL_FLAGS[0] + SEP + MOL_FLAGS[1]:<15}{'the molecule type for the output file [nucl|prot] (default: ' + DEFAULT_MOL + ')':<}" + EOL + \
              GAP + f"{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<15}{'an output file (default: stdout)':<}" + EOL + \
              GAP + f"{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<15}{'print this message'}" + EOL*2
        
        print(MSG)

    # set default values
    gbk = None
    tag = None
    mol = DEFAULT_MOL
    out = DEFAULT_OUT
        
    # determine if help was requested
    helpRequested = HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1
    
    if helpRequested:
        help()
    
    else:
        # parse command line arguments
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        
        # extract each command line argument
        for opt,arg in opts:
            # get the genbank filename
            if opt in GBK_FLAGS:
                if not os.path.isfile(arg):
                    raise FileNotFoundError(ERR_MSG_1)
                gbk = arg
            
            # get the locus tag(s)
            elif opt in TAG_FLAGS:
                tag = arg.split(",")
        
            # get the molecule type
            elif opt in MOL_FLAGS:
                if arg not in ALLOWED_MOL:
                    raise ValueError(ERR_MSG_2)
                mol = arg
    
            # get the output filehandle
            elif opt in OUT_FLAGS:
                out = open(arg, 'w')
        
            # ignore all other arguments
            else:
                print(IGNORE_MSG + opt)
        
        # check for required arguments
        if None in (gbk, tag):
            raise RuntimeError(ERR_MSG_3)
    
    return gbk, tag, mol, out, helpRequested


def __main():
    """main runner function

    Raises:
        RuntimeError: failed to extract all specified tags
    """
    # constants
    FORMAT = "fasta"
    ERR_MSG = "could not extract all specified locus tags"
    
    # parse command line arguments
    gbk,tags,mol,fh,helpRequested = __parseArgs()
    
    
    if not helpRequested:
        # extract the locus tags as SeqRecords
        recs = __extractTags(tags, gbk, mol)
        
        # make sure all tags were extracted
        if len(recs) != len(tags):
            raise RuntimeError(ERR_MSG)
        
        # write records to file
        SeqIO.write(recs, fh, FORMAT)
        
        # close the file if not stdout
        if fh != sys.stdout:
            fh.close()


if __name__ == "__main__":
    __main()
