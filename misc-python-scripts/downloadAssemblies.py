# Joseph S. Wirth
# May 2023

import ftplib, glob, gzip, os, re, shutil, sys
from Bio import Entrez
from Bio.Entrez import Parser


def __parseArgs() -> tuple[str,list[str],str]:
    """ parses command line arguments

    Raises:
        BaseException: improper syntax (wrong num arguments)
        BaseException: email address not specified
        ValueError: email address invalid
        BaseException: email address not specified
        BaseException: accession filename not specified
        BaseException: accession filename not specified
        BaseException: specified output dir flag without a value
        ValueError: specified output directory is not a directory

    Returns:
        tuple[str,list[str],str]: email address; list of accessions; outdir
    """
    # constants
    BASIC_ERR = "invalid usage. use -h flag for help."
    EMAIL_FLAG = "-e"
    EMAIL_ERR = "must specify valid email address"
    INPUT_FLAG = "-i"
    INPUT_ERR = "must specify valid filename containing accession numbers (one per line)"
    OUTDIR_FLAG = "-o"
    OUTDIR_ERR = "specified output is not a directory"

    # check for the expected number of arguments
    if len(sys.argv) not in [5,7]:
        raise BaseException(BASIC_ERR)

    # get the email address
    if EMAIL_FLAG in sys.argv:
        # should be one index beyond the flag
        idx = sys.argv.index(EMAIL_FLAG)
        try:
            email = sys.argv[idx+1]
        except:
            raise BaseException(EMAIL_ERR)
    
        # make sure the email address is valid
        if not __validEmailAddress(email):
            raise ValueError(EMAIL_ERR)
    
    # email address is a required parameter
    else:
        raise BaseException(EMAIL_ERR)
    
    # get the accession numbers
    if INPUT_FLAG in sys.argv:
        # the accession filename should be one index beyond the flag
        idx = sys.argv.index(INPUT_FLAG)
        try:
            accnFN = sys.argv[idx+1]
        except:
            raise BaseException(INPUT_ERR)
        
        # extract the accession numbers from the file
        accnL = __parseAccessionFile(accnFN)
    
    # accession file is a required parameter
    else:
        raise BaseException(INPUT_ERR)
    
    # get the output directory
    if OUTDIR_FLAG in sys.argv:
        # directory should be one index beyond the flag
        idx = sys.argv.index(OUTDIR_FLAG)
        try:
            outdir = sys.argv[idx+1]

        except:
            raise BaseException(BASIC_ERR)
        
        # make the directory if it does not exist
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        # make sure it is a directory if it does exist
        else:
            if not os.path.isdir(outdir):
                raise ValueError(OUTDIR_ERR)
        
        
    
    # outdir is optional; default to current working directory if absent
    else:
        outdir = os.getcwd()
    
    return email,accnL,outdir


def __printHelpMsg() -> None:
    """ prints a help message if requested
    """
    GAP = "    "
    HELP_MSG = "\nDownloads gbff files from a file containing NCBI Assembly Accession numbers\n" + \
               GAP + "Joseph S. Wirth, 2023\n\n" + \
               "Usage:\n" + GAP + "python3 downloadAssemblies.py [-eioh]\n\n" + \
               "Required arguments:\n" + \
               GAP + "-e" + GAP*2 + "email address (for communicating with NCBI; not stored)\n" + \
               GAP + "-i" + GAP*2 + "input file containing accession numbers (one per line)\n\n" + \
               "Optional arguments:\n" + \
               GAP + "-o" + GAP*2 + "output directory where files will be downloaded\n" + \
               GAP + "-h" + GAP*2 + "print this help message; incompatible with all other flags\n"
               
    print(HELP_MSG)


def __parseAccessionFile(accnFN:str) -> list[str]:
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


def __makeAssemblySearchString(accnL:list) -> str:
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


def __assemblyIdsFromSearchTerm(searchStr:str, retmax:int) -> Parser.ListElement[str]:
    """ retrieves uids from the assembly database for a given search term

    Args:
        searchStr (str): the search term
        retmax (int): the max number of uids to retrieve

    Returns:
        Parser.ListElement (str): a list (Bio.Entrez format) of uids
    """
    # constants
    DATABASE = "assembly"

    # retrieve result from NCBI
    handle = Entrez.esearch(db=DATABASE, term=searchStr, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()

    # return the list of ids that were found
    return result['IdList']


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
    handle = Entrez.esummary(id=queryL[0], db=DATABASE)
    result = Entrez.read(handle, validate=VALIDATE)
    handle.close()

    # if there are multiple queries, then keep querying NCBI until done
    for i in range(1,len(queryL)):
        query = queryL[i]

        handle = Entrez.esummary(id=query, db=DATABASE)
        nextResult = Entrez.read(handle, validate=VALIDATE)
        handle.close()

        # something I don't wuote understand is happening here.
        # sometimes, Entrez returns a list, other times its a dict
        # (technically, they are Entrez's stupid format, and not list or dict)
        # the problem may be database-specific:
        ### I think taxonomy returns list and assembly returns dict
        ### it's not worth resolving since the try-except works
        try:
            # append the results if they are returned as lists
            result.extend(nextResult)
        except:
            # append the results if they are returned as dicts
            result[RESULT_K1][RESULT_K2] += nextResult[RESULT_K1][RESULT_K2]

    return result


def __getAssemblySummary(assId) -> list[Parser.DictionaryElement]:
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


def __getFtpPathFromAssSummary(assSummary:Parser.DictionaryElement) -> str:
    """ extracts the FTP path from an assembly's summary

    Args:
        assSummary (Parser.DictionaryElement): an assembly summary

    Returns:
        str: the FTP path (.gz)
    """
    # constants
    SLASH = '/'
    FTP_REF_SEQ = 'FtpPath_RefSeq'
    FTP_GENBANK = 'FtpPath_GenBank'
    FTP_SUFFIX = '_genomic.gbff.gz'
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


def __downloadGbff(ftpPath:str, gbffDir:str) -> None:
    """ a wrapper function to download gbff from NCBI using FTP protocol

    Args:
        ftpPath (str): the FTP address
        gbffDir (str): the directory where the file should be saved
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


def __cleanup() -> None:
    """ removes the unused dtd files from the current working directory
    """
    PATTERN = "*dtd"
    
    unneededFiles = glob.glob(os.path.join(os.curdir, PATTERN))
    for fn in unneededFiles:
        os.remove(fn)


def __main() -> None:
    """ main wrapper function

    Raises:
        BaseException: invalid usage
    """
    # constants
    HELP_FLAG = "-h"
    BASIC_ERR = "invalid usage. use -h flag for help."

    if len(sys.argv) == 1:
      __printHelpMsg() 

    elif len(sys.argv) == 2:
        if sys.argv[1] == HELP_FLAG:
            __printHelpMsg()
            return
        else:
            raise BaseException(BASIC_ERR)
    else:
        # parse the arguments
        email,accnL,outdir = __parseArgs()

        # set the email address
        Entrez.email = email

        # get a list of uids for the assembly database
        searchStr = __makeAssemblySearchString(accnL)
        uidsL = __assemblyIdsFromSearchTerm(searchStr, len(accnL))

        # use the uids to get assembly summaries
        sumL = __getAssemblySummary(uidsL)

        # download each of the files and extract them
        for sum in sumL:
            ftp = __getFtpPathFromAssSummary(sum)
            __downloadGbff(ftp, outdir)

        # remove unnecessary files
        __cleanup()

###############################################################################
###############################################################################


if __name__ == "__main__":
    __main()
