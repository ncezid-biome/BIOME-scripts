# Joseph S. Wirth
# August 2023

import getopt, glob, multiprocessing, os, subprocess, sys
from io import TextIOWrapper


def __downloadOneRead(srr:str, outdir:str, numThreads:int) -> None:
    """downloads a read by calling `fasterq-dump`

    Args:
        srr (str): the SRR number of the desired read
        outdir (str): the directory where the reads should be downloaded
        numThreads (int): the number of threads to use for this call
    """
    # constants
    CMD_PREFIX = ('fasterq-dump', '--split-files', '--threads')
    CMD_SUFFIX = '--outdir'
    MAX_ATTEMPTS = 3
    EOL = "\n"
    FAIL_MSG_1 = 'failed to download reads for '
    FAIL_MSG_2 = "    Error message:  "
    
    # build the command
    cmd = list(CMD_PREFIX)
    cmd.extend([str(numThreads), CMD_SUFFIX, outdir, srr])

    attempts = 0
    while attempts < MAX_ATTEMPTS:
        # run the command; exit loop if successful
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            failed = False
            break

        # try again if the run failed; mark it as a failure
        except subprocess.CalledProcessError as e:
            attempts += 1
            failed = True

    # print error if it failed
    if failed:
        sys.stderr.write(FAIL_MSG_1 + srr + EOL)
        sys.stderr.write(FAIL_MSG_2 + str(e) + EOL*2)


def __downloadReads(srrL:list, outdir:str, maxThreads:int, threadsPerCall:int) -> None:
    """downloads reads in parallel from SRA

    Args:
        srrL (list): a list of SRR numbers to download
        outdir (str): the directory where the reads should be saved
        maxThreads (int): the total number of threads that are available
        threadsPerCall (int): the number of threads that should be allocated to each call
    """
    # constants
    PATTERN = "*fastq"
    
    # make sure the total number of threads does not exceed the number specified
    numParallel = int(maxThreads // threadsPerCall)
    
    # get a set of SRR numbers that have already been downloaded
    existingFiles = glob.glob(os.path.join(outdir, PATTERN))
    existingFiles = {os.path.splitext(os.path.basename(x))[0][:-2] for x in existingFiles} # `<path>/SRR123456_1.fastq` -> `SRR123456`

    # build a list of arguments to be passed to __downloadOneRead
    argsL = list()
    for srr in srrL:
        # skip any existing files
        if srr not in existingFiles:
            argsL.append((srr, outdir, threadsPerCall))
    
    # map the arguments with the available processors
    pool = multiprocessing.Pool(processes=numParallel)
    pool.starmap(__downloadOneRead, argsL)
    pool.close()
    pool.join()


def _getSrrList(srrFH:TextIOWrapper) -> list[str]:
    """gets a list of SRR numbers to download

    Args:
        srrFH (TextIOWrapper): filehandle to a file containing SRR numbers; one per line

    Returns:
        list[str]: list of SRR numbers
    """
    # constant
    EOL = "\n"
    
    # add each SRR number (one per line) to the list
    outL = list()
    for line in srrFH:
        # chomp eol characters
        if line[-1] == EOL:
            line = line[:-1]
        
        # ignore empty lines
        if line != '':
            outL.append(line)
    
    return outL


def __helpMessage() -> None:
    """ prints the help message
    """
    GAP = 4*" "
    EOL = "\n"
    MSG = "\nDownloads SRA files in parallel" + EOL + \
          GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
          "Usage:" + EOL + \
          GAP + "downloadSRA.py [-idnph]" + EOL*2 + \
          "required arguments:" + EOL + \
          GAP + f'{"-i, --input":<30}{"[str] a file containing one SRR ID per line":<}' + EOL*2 + \
          "optional arguments:\n" + \
          GAP + f'{"-d, --dir":<30}{"[str] the directory where reads will be saved":<54}{"(default: current wd)":<}' + EOL + \
          GAP + f'{"-t, --max_threads":<30}{"[int] the maximum allowed number of threads":<54}{"(default: 1)":<}' + EOL + \
          GAP + f'{"-p, --threads_per_download":<30}{"[int] the number of threads to use for each download":<54}{"(default: 1)":<}' + EOL + \
          GAP + f'{"-h, --help":<30}{"print this help message":<}' + EOL
        
    sys.stdout.write(MSG)


def __parseArgs() -> tuple[str,str,int,int,bool]:
    """parses command line arguments

    Raises:
        FileNotFoundError: input file is invalid
        NotADirectoryError: output directory is invalid
        ValueError: invalid number of threads
        ValueError: invalid number of threads per download
        ValueError: must specify an input file
        ValueError: number of threads per download cannot exceed total number of threads

    Returns:
        tuple[str,str,int,int,bool]: input filename, output directory, max number of threads, threads per call, help requested?
    """
    # flags
    INPUT_FLAGS = ("-i", "--in")
    DIR_FLAGS = ("-d", "--dir")
    THREADS_FLAGS_1 = ("-t", "--max_threads")
    THREADS_FLAGS_2 = ("-p", "--threads_per_download")
    HELP_FLAGS = ("-h", "--help")
    SHORT_OPTS = INPUT_FLAGS[0][-1] + ":" + \
                 DIR_FLAGS[0][-1] + ":" + \
                 THREADS_FLAGS_1[0][-1] + ":" + \
                 THREADS_FLAGS_2[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (INPUT_FLAGS[1][2:] + "=",
                 DIR_FLAGS[1][2:] + "=",
                 THREADS_FLAGS_1[1][2:] + "=",
                 THREADS_FLAGS_2[1][2:] + "=",
                 HELP_FLAGS[1][2:])
    
    # default values
    DEFAULT_DIR = os.path.curdir
    DEFAULT_PARALLEL_CALLS = 1
    DEFAULT_THREADS_PER_CALL = 1
    
    # messages
    IGNORE_MSG = "ignoring unused argument "
    ERR_MSG_1 = "invalid input file"
    ERR_MSG_2 = "invalid output directory"
    ERR_MSG_3 = "invalid number of threads"
    ERR_MSG_4 = "invalid number of threads per download"
    ERR_MSG_5 = "must specify an input file (-i;--in)"
    ERR_MSG_6 = "number of threads per download cannot exceed the maximum number of threads"
    
    # set default values
    inFN = None
    outDir = DEFAULT_DIR
    numThreads = DEFAULT_PARALLEL_CALLS
    threadsPerCall = DEFAULT_THREADS_PER_CALL
    
    # determine if help was requested (or no arguments passed)
    helpRequested = HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1
    
    # print the help message then exit (if help was requested or no arguments)
    if helpRequested:
        __helpMessage()
        
    else:
        # extract command line arguments
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        
        # parse each command line argument
        for opt,arg in opts:
            # get the input filename; must be a file
            if opt in INPUT_FLAGS:
                if not os.path.isfile(arg):
                    raise FileNotFoundError(ERR_MSG_1)
                inFN = os.path.abspath(arg)
            
            # get the output directory; must be a directory
            elif opt in DIR_FLAGS:
                # make the directory if it does not exist; or ensure it is an existing directory
                if not os.path.exists(arg):
                    os.mkdir(arg)
                elif not os.path.isdir(arg):
                    raise NotADirectoryError(ERR_MSG_2)
                outDir = os.path.abspath(arg)
            
            # get the total number of threads
            elif opt in THREADS_FLAGS_1:
                try:
                    numThreads = int(arg)
                except:
                    raise ValueError(ERR_MSG_3)
            
            # get the number of threads per call
            elif opt in THREADS_FLAGS_2:
                try:
                    threadsPerCall = int(arg)
                except:
                    raise ValueError(ERR_MSG_4)
            
            # ignore all other arguments
            else:
                print(IGNORE_MSG + opt)

        # make sure all required arguments have been provided
        if inFN is None:
            raise ValueError(ERR_MSG_5)
        
        # make sure the number of threads per call does not exceed the number of threads
        if threadsPerCall > numThreads:
            raise ValueError(ERR_MSG_6)
    
    return inFN,outDir,numThreads,threadsPerCall,helpRequested
    

def __main() -> None:
    """main runner function
        * parses command line arguments
        * imports the SRR IDs to download
        * downloads the reads for each SRR ID in parallel
    """
    # parse the arguments
    inFN,outDir,numThreads,threadsPerCall,helpRequested = __parseArgs()
    
    # no work to do if help requested
    if not helpRequested:
        # get the SRR numbers to be downloaded
        with open(inFN, 'r') as fh:
            srrL = _getSrrList(fh)
        
        # download the reads in parallel
        __downloadReads(srrL, outDir, numThreads, threadsPerCall)


if __name__ == "__main__":
    __main()
