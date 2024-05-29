#!/usr/bin/env python3

import ftplib, getopt, glob, gzip, os, re, shutil, sys, time
from Bio import Entrez
from Bio.Entrez import Parser

__author__ = "Joseph S. Wirth"


def __parseArgs() -> tuple[str,str,str,str,bool,bool]:
    """parses command line arguments

    Raises:
        ValueError: missing or invalid input file
        ValueError: specified output is not a directory
        ValueError: invalid email address
        ValueError: missing one or more required arguments

    Returns:
        tuple[str,str,str,str,bool,bool]: email, in file, format, out directory, rename files?, help requested?
    """
    # flags
    IN_FLAGS = ("-i", "--in")
    EMAIL_FLAGS = ("-e", "--email")
    FORMAT_FLAGS = ("-f", "--format")
    OUT_FLAGS = ("-o", "--out")
    RENAME_FLAGS = ("-r", "--rename")
    HELP_FLAGS = ("-h", "--help")
    SHORT_OPTS = IN_FLAGS[0][-1] + ":" + \
                 EMAIL_FLAGS[0][-1] + ":" + \
                 FORMAT_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 RENAME_FLAGS[0][-1] + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (IN_FLAGS[1][2:] + "=",
                 EMAIL_FLAGS[1][2:] + "=",
                 FORMAT_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 RENAME_FLAGS[1][2:],
                 HELP_FLAGS[1][2:])
    
    # messages
    ERR_MSG_1 = "input file is invalid or missing"
    ERR_MSG_8 = "specified output is not a directory"
    ERR_MSG_2 = "specified email address is invalid"
    ERR_MSG_3 = "invalid file format specified: "
    ERR_MSG_4 = "must specify all required arguments"
    IGNORE_MSG = "ignoring unused argument: "
    
    # constant
    ALLOWED_FORMATS = ("genbank", "fasta")
    
    # helper function to print the help message
    def printHelpMsg() -> None:
        """ prints a help message if requested
        """
        GAP = " "*4
        EOL = "\n"
        SEP = ", "
        HELP_MSG = EOL + "Downloads gbff files from a file containing NCBI Assembly Accession numbers" + EOL + \
                   GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                   "Usage:" + EOL + \
                   GAP + "python3 downloadAssemblies.py [-ieoh]" + EOL*2 + \
                   "Required arguments:" + EOL + \
                   GAP + f'{IN_FLAGS[0] + SEP + IN_FLAGS[1]:<16}{"input file containing accession numbers (one per line)"}' + EOL + \
                   GAP + f'{EMAIL_FLAGS[0] + SEP + EMAIL_FLAGS[1]:<16}{"email address (for communicating with NCBI; not stored)"}' + EOL*2 + \
                   "Optional arguments:" + EOL + \
                   GAP + f'{FORMAT_FLAGS[0] + SEP + FORMAT_FLAGS[1]:<16}{"file format to download [fasta|genbank] (default: genbank)"}' + EOL + \
                   GAP + f'{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<16}{"output directory where files will be downloaded (default: cwd)"}' + EOL + \
                   GAP + f'{RENAME_FLAGS[0] + SEP + RENAME_FLAGS[1]:<16}{"rename downloaded files to match the input file (default: no)"}' + EOL + \
                   GAP + f'{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<16}{"print this message"}' + EOL
                
        print(HELP_MSG)

    # set default values
    fn = None
    email = None
    frmt = ALLOWED_FORMATS[0]
    outdir = os.getcwd()
    rename = False
    helpRequested = False
    
    # print help if requested or if no arguments provided
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequested = True
        printHelpMsg()

    else:
        # parse then command line arguments
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get the input file
            if opt in IN_FLAGS:
                if not os.path.isfile(arg):
                    raise ValueError(ERR_MSG_1)
                fn = arg
            
            # get the email address
            elif opt in EMAIL_FLAGS:
                if not __validEmailAddress(arg):
                    raise ValueError(ERR_MSG_2)
                email = arg
            
            # get the file format
            elif opt in FORMAT_FLAGS:
                if not arg in ALLOWED_FORMATS:
                    raise ValueError(f"{ERR_MSG_3}'{arg}'")
                frmt = arg
            
            # get the out directory
            elif opt in OUT_FLAGS:
                if not os.path.exists(arg):
                    os.mkdir(arg)
                elif not os.path.isdir(arg):
                    raise ValueError(ERR_MSG_8)
                outdir = arg
            
            # determine if renaming will occur
            elif opt in RENAME_FLAGS:
                rename = True
            
            # ignore all other arguments
            else:
                print(IGNORE_MSG + opt + " " + arg)
        
        # must provide an email and an input file    
        if None in (email,fn):
            raise ValueError(ERR_MSG_4)
    
    return email,fn,frmt,outdir,rename,helpRequested


def _parseAccessionFile(accnFN:str) -> list[str]:
    """ reads a file and extracts the accession numbers from it

    Args:
        accnFN (str): filename containing accession numbers (1 per line)

    Returns:
        list[str]: list of accession numbers
    """
    # initialize the list
    accnL = list()

    # extract each accession number; 1 per line
    fh = open(accnFN, 'r')
    for line in fh.readlines():
        # truncate any newline characters from end of line before appending
        if line[-1] == '\n':
            line = line[:-1]
        if line != "":
            accnL.append(line)

    return accnL


def __validEmailAddress(email:str) -> bool:
    """ checks if an email address is validly formatted

    Args:
        email (str): email address to evaluate

    Returns:
        bool: indicating if the address is valid
    """
    GREP_FIND_1 = r'\s'
    GREP_FIND_2 = r'^[^@]+@[^@]+\.[^@]+$'
    
    # if whitespace present, then invalid
    if re.search(GREP_FIND_1, email):
        return False
    
    # one '@' followed by at least one '.'
    elif re.search(GREP_FIND_2, email):
        return True
    
    # if here, then false
    return False


def _makeAssemblySearchString(accnL:list) -> str:
    """ constructs a search string for the assembly database

    Args:
        accnL (list): list of accession numbers (str)

    Returns:
        str: search string
    """
    # constants
    OR = " OR "

    # initialize the string then add each accession number to it
    searchStr = ""
    for accn in accnL:
        searchStr += accn +  OR
    
    # truncate the trailing OR
    searchStr = searchStr[:-len(OR)]

    return searchStr


def _assemblyIdsFromSearchTerm(searchStr:str, retmax:int) -> Parser.ListElement[str]:
    """ retrieves uids from the assembly database for a given search term

    Args:
        searchStr (str): the search term
        retmax (int): the max number of uids to retrieve

    Returns:
        Parser.ListElement (str): a list (Bio.Entrez format) of uids
    """
    # constants
    DATABASE = "assembly"
    MAX_TRIES = 5
    PAUSE = 1
    
    # initialize variables
    numTries = 0
    result = None
    
    # keep trying until success
    while result is None and numTries < MAX_TRIES:
        # retrieve result from NCBI
        try:
            handle = Entrez.esearch(db=DATABASE, term=searchStr, retmax=retmax)
            result = Entrez.read(handle)
            handle.close()
        
        # reset on a failure; increment number of tries
        except:
            result = None
            numTries += 1
            time.sleep(PAUSE)
    
    if result is None:
        raise RuntimeError(f"failed to retrieve ids from {DATABASE} database for '{searchStr}'")

    # return the list of ids that were found
    return result['IdList']


def __esummary(uid:str, database:str, validate:bool) -> Parser.DictionaryElement:
    """performs Entrez esummary queries to NCBI

    Args:
        uid (str): the uid to search
        database (str): the database to search
        validate (bool): indicates if a result should be validated

    Raises:
        RuntimeError: failure to retrieve a summary

    Returns:
        Parser.DictionaryElement: the esummary result
    """
    # constants
    MAX_TRIES = 5
    PAUSE = 1
    
    # initialize variables
    result = None
    numTries = 0
    
    # keep trying until a success
    while result is None and numTries < MAX_TRIES:
        # retrieve data from NCBI
        try:
            handle = Entrez.esummary(id=uid, db=database)
            result = Entrez.read(handle, validate=validate)
            handle.close()
        
        # reset on a failure; increment tries
        except:
            result = None
            numTries += 1
            time.sleep(PAUSE)
    
    if result is None:
        raise RuntimeError(f"failed to retrieve summary from {database} database for '{uid}'")
    
    return result


def __assemblySummaryFromIdList(idList:list) -> Parser.DictionaryElement:
    """ retrieves the summaries from the assembly database for a list of uids

    Args:
        idList (list): a list of uids (Parser.ListElement ok)

    Raises:
        ValueError: invalid input type (expects a list of ids or a single id)
        ValueError: invalid id present (all ids should be str(int))

    Returns:
        Parser.DictionaryElement: Bio.Entrez version of a dictionary
    """
    # constants
    DATABASE = "assembly"
    SEP_CHAR  = ', '
    SEP_LEN   = len(SEP_CHAR)
    RET_MAX   = 10000    # this is NCBI's default; cannot be exceeded
    RESULT_K1 = 'DocumentSummarySet'
    RESULT_K2 = 'DocumentSummary'
    VALIDATE  = False   # Bio.Entrez may throw a tantrum if modified

    # if not a list, assume it is an int or a string
    if type(idList) is str:
        try:
            int(idList)
        except ValueError:
            raise ValueError("invalid input: expected a list of ids or a single id")

    if type(idList) in [int, str]:
        queryL = [str(idList)]

    # otherwise assume it is a list (should work with stupid Entrez formats, too)
    else:
        # convert the list into a comma-separated string
        queryL = list()
        idStr  = str()
        count  = 0
        for id in idList:
            # make sure the ids are valid
            try:
                int(id)
            except ValueError:
                raise ValueError("invalid id present in the input list: " + str(id))

            # separate each id with a comma
            idStr += str(id) + SEP_CHAR
            count += 1

            # make sure each string contains no more than max allowed ids
            if count % RET_MAX == 0:
                # remove the trailing SEP_CHAR from the string
                idStr = idStr[:-SEP_LEN]

                # add the query to the list and begin making the next query
                queryL.append(idStr)
                idStr = str()
            
            # once all ids have been added to a query string
            elif count == len(idList):
                # remove the trailing SEP_CHAR from the string
                idStr = idStr[:-SEP_LEN]
                # add the query to the list
                queryL.append(idStr)

    # search NCBI with the first query
    result = __esummary(queryL[0], DATABASE, VALIDATE)

    # if there are multiple queries, then keep querying NCBI until done
    for i in range(1,len(queryL)):
        nextResult = __esummary(queryL[i], DATABASE, VALIDATE)
        result[RESULT_K1][RESULT_K2].append(nextResult[RESULT_K1][RESULT_K2])

    return result


def _getAssemblySummary(assId) -> list[Parser.DictionaryElement]:
    """ retrieves assembly summaries for a list of uids or a single uid

    Args:
        assId (list or str): a list of uids or a single uid

    Returns:
        list[Parser.DictionaryElement]: list of summaries
    """
    # constants
    RES_F1   = 'DocumentSummarySet'
    RES_F2   = 'DocumentSummary'

    # query NCBI
    result = __assemblySummaryFromIdList(assId)

    # parse the result
    return result[RES_F1][RES_F2]


def _getFtpPathFromAssSummary(assSummary:Parser.DictionaryElement, frmt:str) -> str:
    """ extracts the FTP path from an assembly's summary

    Args:
        assSummary (Parser.DictionaryElement): an assembly summary

    Returns:
        str: the FTP path (.gz)
    """
    # constants
    SLASH = '/'
    GBK = "genbank"
    FNA = "fasta"
    FTP_REF_SEQ = 'FtpPath_RefSeq'
    FTP_GENBANK = 'FtpPath_GenBank'
    FTP_SUFFIX = '_genomic'
    GBK_SUFFIX = '.gbff.gz'
    FNA_SUFFIX = '.fna.gz'
    GREP_FIND = r'^.+/([^/]+)/$'
    GREP_REPLACE = r'\1'

    # extract the ftp path from the entry (folder)
    ftpPath = assSummary[FTP_REF_SEQ]

    # if refseq doesn't have it, then get it from genbank
    if ftpPath == "":
        ftpPath = assSummary[FTP_GENBANK]

    # make sure it ends the correct character
    if ftpPath[-1] != SLASH:
        ftpPath += SLASH
    
    # add the file name of the gbff.gz to the path and return
    ftpPath += re.sub(GREP_FIND, GREP_REPLACE, ftpPath)
    ftpPath += FTP_SUFFIX
    
    if frmt == GBK:
        ftpPath += GBK_SUFFIX
    elif frmt == FNA:
        ftpPath += FNA_SUFFIX

    return ftpPath


def __downloadFileFromFTP(ftpHost:str, ftpFilePath:str, filename:str) -> None:
    """ downloads a file for a given FTP address

    Args:
        ftpHost (str): the host portion of the FTP address
        ftpFilePath (str): the path portion of the FTP address
        filename (str): the local filename where it should be saved
    """
    # initialize the ftp object, then connect and login to the server
    ftp = ftplib.FTP()
    ftp.connect(ftpHost)
    ftp.login()

    # construct the ftp retrieve command
    retrCmd = "RETR " + ftpFilePath
    
    # download the file and write it locally
    with open(filename, 'wb') as fileHandle:
        ftp.retrbinary(retrCmd, fileHandle.write)


def __decompressGZ(gzFileName:str) -> None:
    """ decompresses a gz file

    Args:
        gzFileName (str): the filename of the gz file
    """
    # get gbff file name
    gbffFileName = os.path.splitext(gzFileName)[0]

    # read gz file
    gzFH = gzip.open(gzFileName, "rt") # read text
    
    # create empty gbff file
    gbffFH = open(gbffFileName, 'w')

    # copy contents of decompressed file
    shutil.copyfileobj(gzFH, gbffFH)

    # close files
    gzFH.close()
    gbffFH.close()

    # remove the gzip file
    os.remove(gzFileName)


def _downloadGbff(ftpPath:str, gbffDir:str) -> str:
    """ a wrapper function to download gbff from NCBI using FTP protocol

    Args:
        ftpPath (str): the FTP address
        gbffDir (str): the directory where the file should be saved
    
    Returns:
        str: the filename for the downloaded gbff file
    """
    # constants
    NCBI_FTP  = "ftp.ncbi.nlm.nih.gov"
    LEAD_STR  = "ftp://" + NCBI_FTP

    # remove the leading string from the ftp path
    ftpPath = ftpPath[len(LEAD_STR):]

    # extract the filename from the ftp path
    gzFileName = os.path.basename(ftpPath)

    # add the folder info to the path
    gzFileName = os.path.join(gbffDir, gzFileName)

    # download the file
    __downloadFileFromFTP(NCBI_FTP, ftpPath, gzFileName)

    # get the filename for the unpacked gbff
    gbffFileName = os.path.splitext(gzFileName)[0]

    # unpack the gz
    __decompressGZ(gzFileName)
    
    return gbffFileName


def _renameFile(fn:str, base:str) -> None:
    """renames a file to match the accession number

    Args:
        fn (str): file to be renamed
        base (str): the basename for the renamed file
    """
    # get the directory and extension of the input file
    directory = os.path.dirname(fn)
    extension = os.path.splitext(fn)[1]
    
    # create the new filename and then rename the old file
    newFN = os.path.join(directory, base + extension)
    os.rename(fn, newFN)


def _cleanup() -> None:
    """ removes the unused dtd files from the current working directory
    """
    PATTERN = "*dtd"
    
    unneededFiles = glob.glob(os.path.join(os.curdir, PATTERN))
    for fn in unneededFiles:
        os.remove(fn)


def _runner(email:str, accessions:list[str], frmt:str, outdir:str, rename:bool) -> None:
    """downloads assemblies

    Args:
        email (str): email address
        accessions (list[str]): a list of accession numbers to download
        frmt (str): format [fasta|genbank]
        outdir (str): output directory
        rename (bool): should the sequence files be renamed to match the input file?
    """
    # set the email address
    Entrez.email = email
    
    # print status
    print(f'downloading {frmt}s ... ', flush=True)
    
    # if renaming the files
    if rename:
        # get the uid for each file
        for accn in accessions:
            try:
                uids = _assemblyIdsFromSearchTerm(accn, 1)
            
                # use the uids to get assembly summaries
                summaries = _getAssemblySummary(uids)
            
                # download one file and rename it
                ftp = _getFtpPathFromAssSummary(summaries.pop(), frmt)
                gbffFN = _downloadGbff(ftp, outdir)

                _renameFile(gbffFN, accn)
            
            except:
                print(f'{accn} failed')
            
            time.sleep(2)
    
    # if not renaming files then can process in batch
    else:
        # get all the assembly summaries for all accessions
        search = _makeAssemblySearchString(accessions)
        uids = _assemblyIdsFromSearchTerm(search, len(accessions))
        summaries = _getAssemblySummary(uids)
        
        # download each file
        for summary in summaries:
            ftp = _getFtpPathFromAssSummary(summary, frmt)
            _downloadGbff(ftp, outdir)

    # remove unnecessary files
    _cleanup()
    
    print('done')


def __main() -> None:
    """ main wrapper function
    """
    # parse the arguments
    email,fn,frmt,outdir,rename,helpRequested = __parseArgs()

    # run the program if help wasn't requested
    if not helpRequested:
        # get a list of uids for the assembly database
        accessions = _parseAccessionFile(fn)
        
        # run the rest of the file
        _runner(email, accessions, frmt, outdir, rename)

###############################################################################
###############################################################################


if __name__ == "__main__":
    __main()
