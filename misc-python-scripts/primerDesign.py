# Joseph S. Wirth, 2023

from __future__ import annotations
import getopt, multiprocessing, os, sys
from io import TextIOWrapper
from multiprocessing.managers import ListProxy, DictProxy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp


class Primer:
    def __init__(self, seq:Seq, contig:str, start:int, length:int) -> Primer:
        # initialize the attributes
        self.seq:Seq = None
        self.start:int = start
        self.end:int = start + length
        self.contig:str = contig
        self.Tm:float = None
        self.gcPer:float = None
        self.__length:int = length
        
        # run import methods
        self.__importSeq(seq)
        self.__calcPerGc()
        self.__calculateTm()
    
    # overloads
    def __hash__(self) -> int:
        return self.seq.__hash__()
    
    def __len__(self) -> int:
        return self.__length
    
    def __str__(self) -> str:
        return self.seq.__str__()
    
    def __repr__(self) -> str:
        return str(self)
    
    def __eq__(self, other:Primer) -> bool:
        return self.seq == other.seq
    
    def __ne__(self, other:Primer) -> bool:
        return not self.seq == other.seq
    
    # private methods
    def __importSeq(self, seq:Seq) -> None:
        # make sure that the sequence is upper-case
        self.seq = Seq(str(seq).upper())
    
    def __calcPerGc(self) -> None:
        # constant
        GC = ("G", "C")
        
        # count number of GCs
        numGc = 0
        for base in self.seq:
            if base in GC:
                numGc += 1
        
        # calculate and store the percentage of GC
        self.gcPer = numGc / len(self) * 100

    def __calculateTm(self) -> None:
        # calculate and store the melting temp
        self.Tm = MeltingTemp.Tm_Wallace(self.seq)
    
    # public methods
    def reverseComplement(self) -> Primer:
        """reverse complements the calling object

        Returns:
            Primer: the reverse complement of the calling object
        """
        new = Primer(self.seq.reverse_complement(), self.contig, self.end, len(self))
        new.end = self.start
        
        return new

    def getMinimizer(self, lmerSize:int) -> Seq:
        """finds the minimizer for the calling object

        Args:
            lmerSize (int): the length of the l-mer

        Raises:
            ValueError: primer length must exceed l-mer length

        Returns:
            Seq: the minimizer sequence
        """
        # constants
        ERR_MSG = "Window size should be less than or equal to the k-mer length."
        
        # make sure the lmer is smaller than the primer
        if len(self.seq) < lmerSize:
            raise ValueError(ERR_MSG)

        # Initialize the minimizer and its position
        minimizer = self.seq[:lmerSize]

        # Iterate through the k-mer with the sliding window
        for i in range(1, len(self.seq) - lmerSize + 1):
            current_window = self.seq[i:i + lmerSize]

            # Update the minimizer if the current window is lexicographically smaller
            if current_window < minimizer:
                minimizer = current_window

        return minimizer


def __kmpSearch(text:str, pattern:str) -> bool:
    """implementation of the KMP string search algorithm: O(n+m)

    Args:
        text (str): text to search
        pattern (str): pattern to search for

    Returns:
        bool: indicates if the pattern is found in the text
    """

    # helper function
    def computeLpsArray(patternLen:int) -> list[int]:
        """precomputes the longest prefix that is also a suffix array

        Args:
            patternLen (int): the length of the pattern

        Returns:
            list[int]: the precomputed lps array
        """
        # initialize the lps array
        lps = [0] * patternLen
        
        # initialize the indices; lps[0] always == 0 so start idx at 1
        matchLen = 0 
        idx = 1
        
        # go through each position in the pattern 
        while idx < patternLen:
            # if a match occurs update the match length, store in the array, and move to the next position
            if pattern[idx] == pattern[matchLen]:
                matchLen += 1
                lps[idx] = matchLen
                idx += 1
            
            # if a mismatch
            else:
                # if at previous string also did not match, then no lps; move to next position
                if matchLen == 0:
                    lps[idx] = 0
                    idx += 1 
                
                # if the previous string matched, we need to start at the beginning of the last match
                else:
                    matchLen = lps[matchLen - 1]
        
        return lps

    # get the length of the text and the pattern
    textLen = len(text)
    patternLen = len(pattern)
    
    lps = computeLpsArray(patternLen)
    
    # initialize indices
    idx = 0
    jdx = 0
    
    # keep evaluating the text until at the end
    while idx < textLen:
        # keep comparing strings while they match
        if text[idx] == pattern[jdx]:
            idx += 1
            jdx += 1
        
        # if a mismatch occurs
        else:
            # if at the start of the pattern, cannot reset; move to next text character
            if jdx == 0:
                idx += 1
            
            # otherwise, move to the previous character as defined by the lps array
            else:
                jdx = lps[jdx -1]
    
        # match found if at the end of the pattern
        if jdx == patternLen:
            return True
    
    return False


def __getUniqueKmers(inFN:str, frmt:str, minLen:int, maxLen:int) -> dict:   
    # helper functions to clarify boolean expressions
    def isOneEndGc(seq:Seq) -> bool:
        """ evaluates if one end of a sequence is a GC
        """
        # constants
        GC = {"G", "C", 'g', 'c'}
        
        return seq[-1] in GC or seq[0] in GC

    # initialize variables for looping
    seqsD = dict()
    bad = set()
    rec:SeqRecord
    
    # go through each contig
    for rec in SeqIO.parse(inFN, frmt):
        # go through each primer length
        for primerLen in range(minLen, maxLen+1):
            # get every possible primer start position
            for start in range(len(rec)- (primerLen - 1)):
                # extract the primer sequence
                seq:Seq = rec.seq[start:start+primerLen]
                
                # mark duplicate primers for removal
                if seq in seqsD.keys() or seq.reverse_complement() in seqsD.keys():
                    bad.add(seq)
                
                # only save primers that have GC at one end, no long repeats, and no complements
                if isOneEndGc(seq):
                    seqsD[seq] = (rec.id, start, primerLen)

    # discard any sequences that were duplicated
    for seq in bad:
        seqsD.pop(seq)
    
    return seqsD


def __evaluatePrimersAtOnePosition(contig:str, start:int, posL:list[tuple[Seq,int]], minGc:float, maxGc:float, minTm:float, maxTm:float, shareL:ListProxy) -> None:
    """evaluates all the primers at a single position in the genome; designed for parallel calls

    Args:
        contig (str): the name of the contig
        start (int): the start position in the sequence
        posL (list[tuple[Seq,int]]): a list of tuples; primer sequence and primer length
        minGc (float): the minimum percent GC allowed
        maxGc (float): the maximum percent GC allowed
        minTm (float): the minimum melting temperature allowed
        maxTm (float): the maximum melting temperature allowed
        shareL (ListProxy): a shared list for parallel calling;
    
    Returns:
        Does not return.
        Primer sequences passing the boolean checks are added to the shared list
    """
    # define helper functions to make booleans below more readable
    def isGcWithinRange(primer:Primer) -> bool:
        """is the percent GC within the acceptable range?"""
        return primer.gcPer >= minGc and primer.gcPer <= maxGc

    def isTmWithinRange(primer:Primer) -> bool:
        """is the Tm within the acceptable range?"""
        return primer.Tm >= minTm and primer.Tm <= maxTm
    
    def noLongRepeats(primer:Primer) -> bool:
        """verifies that a primer does not have long repeats
            O(1)
        """
        # constants
        MAX_LEN = 4
        REPEATS = ("A"*MAX_LEN, "T"*MAX_LEN, "C"*MAX_LEN, "G"*MAX_LEN)

        # check each repeat in the primer
        for repeat in REPEATS:
            if __kmpSearch(primer.seq, repeat):
                return False
        return True

    def noIntraPrimerComplements(primer:Primer) -> bool:
        """verifies that the primer does not have hairpin potential
            O(len(primer))!
        """
        # constants
        MAX_LEN = 3
        LEN_TO_CHECK = MAX_LEN + 1
        
        # go through each frame of the primer up to the max length
        for idx in range(len(primer)-MAX_LEN):
            # get the fragment and see if the reverse complement exists downstream
            revComp = primer.seq.reverse_complement()
            if __kmpSearch(primer.seq, revComp[idx:idx+LEN_TO_CHECK]):
                return False
        
        return True
    
    # initialize values for the while loop
    found = False
    idx = 0
    
    # continue to iterate through each primer in the list until a primer is found 
    while idx < len(posL) and not found:
        # extract data from the list
        seq,length = posL[idx]
        
        # create a Primer object
        primer = Primer(seq, contig, start, length)
        
        # evaluate the primer's percent GC, Tm, and homology; save if found
        if isGcWithinRange(primer) and isTmWithinRange(primer): # O(1)
            if noLongRepeats(primer): # O(1)
                if noIntraPrimerComplements(primer): # this runtime is the worst O(len(primer)); evaluate last
                    shareL.append(primer)
                    found = True
                
        # move to the next item in the list
        idx += 1


def __getAllCandidatePrimers(inFN:str, frmt:str, minLen:int, maxLen:int, minGc:float, maxGc:float, minTm:float, maxTm:float, numThreads:int) -> dict[str:list[Primer]]:
    """retrieves a list of all viable primer kmers from an input genome

    Args:
        inFN (str): the input filename
        frmt (str): the format of the file (genbank, fasta)
        minLen (int): the minimum primer length
        maxLen (int): the maximum primer length
        minGc (float): the minimum percent GC
        maxGc (float): the maximum percent GC
        minTm (float): the minimum melting temperature
        maxTm (float): the maximum meltding temperature
        numThreads (int): the number of threads to use for parallel processing

    Returns:
        dict[str:list[Primer]]: key = contig name; val = list of primers sorted by start position
    """
    # messages
    MSG_1 = 'finding all unique kmers'
    MSG_2 = 'evaluating candidate kmers at each position in parallel'
    END = ' ... '
    DONE = 'done.'
    
    # print status
    print(MSG_1, end=END, flush=True)
    
    # get a dictionary of all the non-duplicated kmers in the genome
    seqsD = __getUniqueKmers(inFN, frmt, minLen, maxLen)
    
    # reorganize data by each unique start positions
    positionD = dict()
    for seq in seqsD.keys():
        # extract data from the dictionary
        contig, start, length = seqsD[seq]
        
        # contig = top level key; start position = second level key; val = list
        positionD[contig] = positionD.get(contig, dict())
        positionD[contig][start] = positionD[contig].get(start, list())
        
        # add the sequence and its length to the list
        positionD[contig][start].append((seq, length))
    
    # print status
    print(DONE)
    print(MSG_2, end=END, flush=True)
    
    # initialize a shared list and a list of arguments
    primerL = multiprocessing.Manager().list()
    argsL = list()
    
    # each contig needs to be evalutated
    for contig in positionD.keys():
        # each start position within the contig needs to be evaluated
        for start in positionD[contig].keys():
            # save arguments to pass in parallel
            argsL.append((contig, start, positionD[contig][start], minGc, maxGc, minTm, maxTm, primerL))

    # parallelize primer evaluations
    pool = multiprocessing.Pool(processes=numThreads)
    pool.starmap(__evaluatePrimersAtOnePosition, argsL)
    pool.close()
    pool.join()
    
    # collapse the shared list to a list
    primerL:list[Primer] = list(primerL)
    
    # create a dictionary whose keys are contigs and values are the primers
    outD = dict()
    for primer in primerL:
        outD[primer.contig] = outD.get(primer.contig, list())
        outD[primer.contig].append(primer)
    
    # sort lists of primers by start position
    for contig in outD.keys():
        outD[contig] = sorted(outD[contig], key=lambda x: x.start)
    
    # print status
    print(DONE)
    
    return outD


def __evaluateOnePair(p1:Primer, p2:Primer, length:int, queue) -> None:
    """evaluates if a pair of primers is suitable using an alignment-based method
       designed to run in parallel

    Args:
        p1 (Primer): forward primer
        p2 (Primer): reverse primer
        length (int): PCR product length
        queue: a multiprocessing.Manager().Queue() object to store results for writing
    """
    # define a helper function to evaluate primer homology
    def noPrimerDimer() -> bool:
        """verifies that the primer pair will not form primer dimers
        """
        # constant
        MAX_PID = 0.9
        
        # create the aligner; cost-free end gaps; no internal gaps
        from Bio.Align import PairwiseAligner, Alignment
        aligner = PairwiseAligner(mode='global',
                                  end_gap_score=0,
                                  end_open_gap_score=0,
                                  internal_open_gap_score=-float('inf'),
                                  match_score=2,
                                  mismatch_score=-1)
        
        # go through the alignments and check for high percent identities
        aln:Alignment
        for aln in aligner.align(p1.seq, p2.seq):
            matches = aln.counts().identities
            pid1 = matches / len(p1)
            pid2 = matches / len(p2)
            
            if max(pid1,pid2) > MAX_PID:
                return False
        
        return True
    
    # the primer pair is good if it doesn't form dimers
    if noPrimerDimer():
        queue.put((p1,p2.reverseComplement(),length))


def __writePairsToFile(fn:str, terminator:str, queue) -> None:
    """writes primer pairs to file in parallel

    Args:
        fn (str): the filename for writing the results
        terminator (str): the string that will tell this process to terminiate
        queue: a multiprocessing.Manager().Queue() object to store results for writing
    """
    # constants
    SEP = '\t'
    EOL = '\n'
    NUM_DEC = 1
    HEADER = 'contig' + SEP + 'fwd_seq' + SEP + 'fwd_Tm' + SEP + 'fwd_GC' + SEP + \
             'rev_seq' + SEP +'rev_Tm' + SEP + 'rev_GC' + SEP + 'product_len' + EOL
    
    with open(fn, 'w') as fh:
        # write the header
        fh.write(HEADER)
        
        # keep checking the queue
        while True:
            # extract the value from the queue
            val = queue.get()
            
            # stop loooping when the terminate signal is received
            if val == terminator:
                break
            
            # extract the primers and PCR product length from the queue
            p1:Primer
            p2:Primer
            length:int
            p1,p2,length = queue.get()
            
            # write the row to file
            fh.write(p1.contig + SEP)
            fh.write(str(p1.seq) + SEP)
            fh.write(str(round(p1.Tm, NUM_DEC)) + SEP)
            fh.write(str(round(p1.gcPer, NUM_DEC)) + SEP)
            fh.write(str(p2.seq) + SEP)
            fh.write(str(round(p2.Tm, NUM_DEC)) + SEP)
            fh.write(str(round(p2.gcPer, NUM_DEC)) + SEP)
            fh.write(str(length) + EOL)
            fh.flush()


def __getPrimerPairs(primersD:dict, minLen:int, maxLen:int, maxTmDiff:float, outFN:str, numThreads:int) -> None:
    """picks primer pairs from a collection of candidate primers and write them to file

    Args:
        primersD (dict): the dictionary produced by __getAllCandidatePrimers
        minLen (int): min PCR product length
        maxLen (int): max PCR product length
        maxTmDiff (float): maximum difference in melting temps between primers
        outFN (str): filename to write the results
        numThreads (int): number of threads for parallel processing
    """
    # constants
    FWD = 'forward'
    REV = 'reverse'
    TERMINATOR = 'END'
    
    # messages
    MSG_1 = "evaluating 3' ends and primer pair Tm differences"
    MSG_2 = "screening primer pairs for dimer formation in parallel"
    END = ' ... '
    DONE = 'done.'
    
    # helper functions for evaluating primers
    def isThreePrimeGc(primer:Primer, direction:str=FWD) -> bool:
        """checks if a primer has a G or C at its 3' end
        """
        # constant
        GC = {"G", "C"}
        
        # different check depending if forward or reverse primer
        if direction == FWD:
            return primer.seq[-1] in GC
        if direction == REV:
            return primer.seq[0] in GC
    
    def isTmDiffWithinRange(p1:Primer, p2:Primer) -> bool:
        """ checks if the difference between melting temps is within the specified range
        """
        return abs(p1.Tm - p2.Tm) <= maxTmDiff

    # make a queue so results can be written to the file in parallel
    queue = multiprocessing.Manager().Queue()

    # intitialize the pool and start the writer process
    pool = multiprocessing.Pool(processes=numThreads)
    pool.apply_async(__writePairsToFile, (outFN, TERMINATOR, queue))

    # print status
    print(MSG_1, end=END, flush=True)

    # initialize a list of arguments for __evaluateOnePair
    argsL = list()
    
    # process each contig separately
    for contig in primersD.keys():
        # extract the list of primers for this contig
        primersL:list[Primer] = primersD[contig]
        
        # get all pairwise combinations of primers
        for idx in range(len(primersL)-1):
            primer1 = primersL[idx]
            
            # forward primer must have GC at 3' end
            if isThreePrimeGc(primer1): 
                for primer2 in primersL[idx+1:]:
                    # product length must be within the specified range
                    productLen = primer2.end - primer1.start
                    if minLen <= productLen <= maxLen:
                        # reverse primer must have GC at 3' end
                        if isThreePrimeGc(primer2, REV):
                            # tm diff must be acceptable
                            if isTmDiffWithinRange(primer1, primer2):
                                argsL.append((primer1, primer2, productLen, queue))
                        
                    # move on to next primer1 once length is too long
                    elif productLen > maxLen:
                        break

    # print status
    print(DONE)
    print(MSG_2, end=END, flush=True)

    # evaluate candidate pairs and write them in parallel
    pool.starmap(__evaluateOnePair, argsL)
    
    # stop the writer then close the pool
    queue.put(TERMINATOR)
    pool.close()
    pool.join()
    
    # print status
    print(DONE)


def __parseArgs() -> tuple[str,str,str,int,int,float,float,float,float,int,int,float,int,bool]:
    """parses command line arguments

    Raises:
        ValueError: invalid input file
        ValueError: invalid file format
        ValueError: invalid primer length
        ValueError: primer lengths not integers
        ValueError: specify range of GC
        ValueError: GC values non-numeric
        ValueError: specify range of Tm
        ValueError: Tm values non-numeric
        ValueError: invalid PCR product length
        ValueError: PCR product lengths not integers
        ValueError: Tm difference is non-numeric
        ValueError: num threads not an integer
        ValueError: must specify an input file

    Returns:
        tuple[str,str,str,int,int,float,float,float,float,int,int,float,int,bool]:
            inFN,outFN,format,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,helpRequested
    """
    # constants
    ALLOWED_FORMATS = ('genbank', 'fasta')
    SEP = ","
    
    # flags
    IN_FLAGS = ('-i', '--in')
    OUT_FLAGS = ('-o', '--out')
    FMT_FLAGS = ('-f', '--format')
    PRIMER_LEN_FLAGS = ('-p', '--primer_len')
    GC_FLAGS = ('-g', '--gc_range')
    TM_FLAGS = ('-t', '--tm_range')
    THREADS_FLAGS = ('-n', '--num_threads')
    PCR_LEN_FLAGS = ('-r', '--pcr_prod_len')
    TM_DIFF_FLAGS = ('-d', '--tm_diff')
    HELP_FLAGS = ('-h', '--help')
    SHORT_OPTS = IN_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 FMT_FLAGS[0][-1] + ":" + \
                 PRIMER_LEN_FLAGS[0][-1] + ":" + \
                 GC_FLAGS[0][-1] + ":" + \
                 TM_FLAGS[0][-1] + ":" + \
                 PCR_LEN_FLAGS[0][-1] + ":" + \
                 TM_DIFF_FLAGS[0][-1] + ":" + \
                 THREADS_FLAGS[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (IN_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 FMT_FLAGS[1][2:] + "=",
                 PRIMER_LEN_FLAGS[1][2:] + "=",
                 GC_FLAGS[1][2:] + "=",
                 TM_FLAGS[1][2:] + "=",
                 PCR_LEN_FLAGS[1][2:] + "=",
                 TM_DIFF_FLAGS[1][2:] + "=",
                 THREADS_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:])

    # default values
    DEF_FRMT = ALLOWED_FORMATS[0]
    DEF_MIN_LEN = 16
    DEF_MAX_LEN = 20
    DEF_MIN_GC = 40.0
    DEF_MAX_GC = 60.0
    DEF_MIN_TM = 55.0
    DEF_MAX_TM = 68.0
    DEF_MIN_PCR = 120
    DEF_MAX_PCR = 2400
    DEF_MAX_TM_DIFF = 5.0
    DEF_NUM_THREADS = 1

    # messages
    IGNORE_MSG = 'ignoring unused argument: '
    ERR_MSG_1  = 'invalid or missing input file'
    ERR_MSG_2  = 'invalid format'
    ERR_MSG_3  = 'can only specify one primer length or a range (min,max)'
    ERR_MSG_4  = 'primer lengths are not integers'
    ERR_MSG_5  = 'must specify a range of GC values (min,max)'
    ERR_MSG_6  = 'gc values are not numeric'
    ERR_MSG_7  = 'must specify a range of Tm values (min, max)'
    ERR_MSG_8  = 'Tm values are not numeric'
    ERR_MSG_9  = 'can only specify one PCR product length or a range (min,max)'
    ERR_MSG_10 = 'PCR product lengths are not integers'
    ERR_MSG_11 = 'max Tm difference is not numeric'
    ERR_MSG_12 = 'num threads is not an integer'
    ERR_MSG_13 = 'must specify an input file'
    ERR_MSG_14 = 'must specify an output file'

    def printHelp():
        GAP = " "*4
        EOL = "\n"
        SEP_1 = ", "
        SEP_2 = "|"
        SEP_3 = ","
        DEF_OPEN = ' (default: '
        CLOSE = ')'
        HELP_MSG = EOL + "Finds pairs of primers suitable for an input genome." + EOL + \
                   GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                   "usage:" + EOL + \
                   GAP + "python3 primerDesign.py [-iofpgtnrdh]" + EOL*2 + \
                   "required arguments:" + EOL + \
                   GAP + f"{IN_FLAGS[0] + SEP_1 + IN_FLAGS[1]:<22}{'[file] input filename'}" + EOL + \
                   GAP + f"{OUT_FLAGS[0] + SEP_1 + OUT_FLAGS[1]:<22}{'[file] output filename'}" + EOL*2 + \
                   "optional arguments:" + EOL + \
                   GAP + f"{FMT_FLAGS[0] + SEP_1 + FMT_FLAGS[1]:<22}{'[str] input file format '}{{{ALLOWED_FORMATS[0] + SEP_2 + ALLOWED_FORMATS[1]}}}{DEF_OPEN + DEF_FRMT + CLOSE}" + EOL + \
                   GAP + f"{PRIMER_LEN_FLAGS[0] + SEP_1 + PRIMER_LEN_FLAGS[1]:<22}{'[int(s)] a single primer length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_LEN) + SEP_3 + str(DEF_MAX_LEN) + CLOSE}" + EOL + \
                   GAP + f"{GC_FLAGS[0] + SEP_1 + GC_FLAGS[1]:<22}{'[float,float] a min and max percent GC specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_GC) + SEP_3 + str(DEF_MAX_GC) + CLOSE}" + EOL + \
                   GAP + f"{TM_FLAGS[0] + SEP_1 + TM_FLAGS[1]:<22}{'[float,float] a min and max melting temp (Tm) specified as a comma separated list' + DEF_OPEN + str(DEF_MIN_TM) + SEP_3 + str(DEF_MAX_TM) + CLOSE}" + EOL + \
                   GAP + f"{PCR_LEN_FLAGS[0] + SEP_1 + PCR_LEN_FLAGS[1]:<22}{'[int(s)] a single PCR product length or a range specified as '}" + "'min,max'" + f"{DEF_OPEN + str(DEF_MIN_PCR) + SEP_3 + str(DEF_MAX_PCR) + CLOSE}" + EOL + \
                   GAP + f"{TM_DIFF_FLAGS[0] + SEP_1 + TM_DIFF_FLAGS[1]:<22}{'[float] the maximum allowable Tm difference between a pair of primers' + DEF_OPEN + str(DEF_MAX_TM_DIFF) + CLOSE}" + EOL + \
                   GAP + f"{THREADS_FLAGS[0] + SEP_1 + THREADS_FLAGS[1]:<22}{'[int] the number of threads for parallel processing' + DEF_OPEN + str(DEF_NUM_THREADS) + CLOSE}" + EOL + \
                   GAP + f"{HELP_FLAGS[0] + SEP_1 + HELP_FLAGS[1]:<22}{'print this message'}" + EOL*2
        
        print(HELP_MSG)
        
    # set default values
    inFN = None
    outFN = None
    frmt = DEF_FRMT
    minLen = DEF_MIN_LEN
    maxLen = DEF_MAX_LEN
    minGc = DEF_MIN_GC
    maxGc = DEF_MAX_GC
    minTm = DEF_MIN_TM
    maxTm = DEF_MAX_TM
    minPcr = DEF_MIN_PCR
    maxPcr = DEF_MAX_PCR
    maxTmDiff = DEF_MAX_TM_DIFF
    numThreads = DEF_NUM_THREADS
    helpRequested = False
    
    # give help if requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequested = True
        printHelp()
    
    # parse command line arguments
    else:
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get the input file
            if opt in IN_FLAGS:
                if not os.path.isfile(arg):
                    raise ValueError(ERR_MSG_1)
                inFN = arg
            
            # get output filehandle
            elif opt in OUT_FLAGS:
                outFN = arg
            
            # get the input file format
            elif opt in FMT_FLAGS:
                if arg not in ALLOWED_FORMATS:
                    raise ValueError(ERR_MSG_2)
                frmt = arg
            
            # get the primer lengths
            elif opt in PRIMER_LEN_FLAGS:
                # split comma-separated list
                primerRange = arg.split(SEP)
                
                # make sure at one or two primers specified
                if len(primerRange) not in {1,2}:
                    raise ValueError(ERR_MSG_3)
                
                # coerce to lengths to ints
                try:
                    primerRange = [int(x) for x in primerRange]
                except:
                    raise ValueError(ERR_MSG_4)
                
                # save values
                minLen = min(primerRange)
                maxLen = max(primerRange)
            
            # get the allowed GC range
            elif opt in GC_FLAGS:
                # expecting two values separated by a comma
                gcRange = arg.split(SEP)
                if len(gcRange) != 2:
                    raise ValueError(ERR_MSG_5)
                
                # make sure the values are numeric
                try:
                    gcRange = [float(x) for x in gcRange]
                except:
                    raise ValueError(ERR_MSG_6)
            
                # save values
                minGc = min(gcRange)
                maxGc = max(gcRange)
            
            # get the allowed Tm range
            elif opt in TM_FLAGS:
                # expecting two values separated by a comma
                tmRange = arg.split(SEP)
                if len(tmRange) != 2:
                    raise ValueError(ERR_MSG_7)
            
                # make sure the values are numeric
                try:
                    tmRange = [float(x) for x in tmRange]
                except:
                    raise ValueError(ERR_MSG_8)
            
                # save values
                minTm = min(tmRange)
                maxTm = max(tmRange)
        
            # get the allowed PCR lengths
            elif opt in PCR_LEN_FLAGS:
                # expecting one or two values
                pcrRange = arg.split(SEP)
                if len(pcrRange) not in {1,2}:
                    raise ValueError(ERR_MSG_9)
                
                # coerce to integers
                try:
                    pcrRange = [int(x) for x in pcrRange]
                except:
                    raise ValueError(ERR_MSG_10)
            
                # save values
                minPcr = min(pcrRange)
                maxPcr = max(pcrRange)
            
            # get the allowed Tm difference between primer pairs
            elif opt in TM_DIFF_FLAGS:
                # make sure input is numeric
                try:
                    maxTmDiff = float(arg)
                except:
                    raise ValueError(ERR_MSG_11)
            
            # get the number of threads to use
            elif opt in THREADS_FLAGS:
                # make sure input is an integer
                try:
                    numThreads = int(arg)
                except:
                    raise ValueError(ERR_MSG_12)
            
            else:
                print(IGNORE_MSG + opt + " " + arg)
        
        # make sure an input file was specified
        if inFN is None:
            raise ValueError(ERR_MSG_13)
        
        # make sure an output file was specified
        if outFN is None:
            raise ValueError(ERR_MSG_14)
    
    return inFN,outFN,frmt,minLen,maxLen,minGc,maxGc,minTm,maxTm,minPcr,maxPcr,maxTmDiff,numThreads,helpRequested


def __main():
    """main runner function
        * parses command line arguments
        * retrieves all the candidate primers
        * finds suitable primer pairs
        * writes results to file
    """
    # parse command line arguments
    inFN,outFN,frmt,minPrimerLen,maxPrimerLen,minGc,maxGc,minTm,maxTm,minPcrLen,maxPcrLen,maxTmDiff,numThreads,helpRequested = __parseArgs()
    
    # no work to do if help was requested
    if not helpRequested:
        # get the candidate primer sequences
        candidatesD = __getAllCandidatePrimers(inFN,
                                            frmt,
                                            minPrimerLen,
                                            maxPrimerLen,
                                            minGc,
                                            maxGc,
                                            minTm,
                                            maxTm,
                                            numThreads)

        # evaluate primer pairs and write them to file in parallel
        __getPrimerPairs(candidatesD,
                         minPcrLen,
                         maxPcrLen,
                         maxTmDiff,
                         outFN,
                         numThreads)


if __name__ == "__main__":
    __main()
