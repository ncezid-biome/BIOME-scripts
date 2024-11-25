#!/usr/bin/env python3

from __future__ import annotations
from Bio import SeqIO
from typing import Generator
from khmer import Countgraph
from itertools import combinations
import os, pickle, shutil, sys, time
from multiprocessing import Event, Manager, Pool, Process, Queue

# global constants
__author__ = "Joseph S. Wirth"
__version__ = "0.0.1"
__JUNCTION_CHAR = "~"
__PICKLE_DIR = os.path.join(os.curdir, ".pickles")
__PICKLE_EXT = ".pkl"


# classes
class Wheel:
    """provides a wheel spinning
    """
    # public attribute
    PAUSE = 0.12
    
    # priave attribute
    __CHARS = "-\|/"
    __EVENT = Event()
        
    # overloads
    def __init__(self) -> Wheel:
        """constructor

        Returns:
            Wheel: a Wheel object
        """
        # initialize attributes
        self.__msg:str = ''
        self.__process:Process

    # private methods
    def __spin(self) -> None:
        """spinning function
        """
        # keep spinning until it is time to stop
        while not self.__EVENT.is_set():
            # print each character in the wheel then pause
            for char in Wheel.__CHARS:
                sys.stdout.write('\r' + self.__msg + char)
                sys.stdout.flush()
                time.sleep(Wheel.PAUSE)
    
    # public methods
    def start(self, msg:str) -> None:
        """starts the wheel spinning

        Args:
            msg (str): the message to preceed the wheel
        """
        # initialize attributes
        self.__msg = msg
        self.__process = Process(target=self.__spin)
        self.__EVENT.clear()
        
        # start spinning the wheel
        self.__process.start()
    
    def stop(self) -> None:
        """stops spinning the wheel
        """
        # this will raise AttributeError if the wheel isn't spinning
        try:
            # stop spinning the wheel
            if self.__process.is_alive():
                self.__EVENT.set()
                self.__process.join()
                
                # remove the wheel character
                sys.stdout.write('\r' + self.__msg)
                sys.stdout.flush()
        
        except AttributeError:
            pass


class Clock:
    """ this class is used to easily track durations in a pretty format
    """
    # private attributes; share a single wheel (easier to kill)
    __WHEEL = Wheel()
    
    # overload
    def __init__(self) -> Clock:
        """ constructor. accepts no inputs
            saves the start time and initializes member variables
        """
        # initialize attributes
        self.__startTime:float = 0.0
        self.__duration:float = 0.0
        self.__spin:bool = True
        
        # start the clock
        self.__start()
    
    # private methods
    def __start(self) -> None:
        """ accepts no inputs
            saves the start time
        """
        self.__startTime = time.time()
    
    def __end(self) -> None:
        """ accepts no inputs
            saves the elapsed duration
        """
        self.__duration = time.time() - self.__startTime
    
    def __parseDuration(self,digi:int) -> tuple[int,int,float]:
        """converts duration from seconds to hours, minutes, and seconds

        Args:
            digi (int): the number of decimal points for the seconds

        Returns:
            tuple[int,int,float]: hours, minutes, seconds
        """
        # constants
        TO_HRS = 3600
        TO_MIN = 60
        
        # parse the time; round seconds to specified digits
        hours = int(self.__duration // TO_HRS)
        minutes = int(self.__duration % TO_HRS // TO_MIN)
        seconds = round(self.__duration - minutes * TO_MIN - hours * TO_HRS, digi)
        
        return hours, minutes, seconds
    
    def __getDurationString(self, digi:int) -> str:
        """converts the duration to a formatted string

        Args:
            digi (int): the number of decimal points for the seconds
        
        Returns:
            str: the duration as a string in hh:mm:ss.ms format
        """
        # constant
        ZERO = "0"
        
        # parse the duration
        hours,minutes,seconds = self.__parseDuration(digi)
        
        # separate seconds from their decimals
        decimal = str(int(10**digi*seconds))
        seconds = str(int(seconds))
        
        # convert hours and minutes to strings
        hours = str(hours)
        minutes = str(minutes)
        
        # make sure hours are at least 2 chars long
        if len(hours) == 1:
            hours = ZERO + hours
        
        # make sure minutes are at least 2 chars long
        if len(minutes) == 1:
            minutes = ZERO + minutes
        
        # make sure seconds (int) are at least 2 chars long
        if len(seconds) == 1:
            seconds = ZERO + seconds
        
        # ensure that the decimal is exactly `digi` chars long
        if len(decimal) > digi:
            decimal = decimal[-digi:]
        elif len(decimal) < digi:
            decimal = ZERO * (digi - len(decimal)) + decimal
        
        # reconstruct the seconds to include decimals unless no decimals were requested
        if digi != 0:
            seconds = seconds + "." + decimal
        
        return ":".join([hours,minutes,seconds])

    # public methods
    def getTime(self, decimals:int=2) -> tuple[int,int,float]:
        """gets the elapsed time as a tuple

        Args:
            decimals (int, optional): the number of deimal points for the seconds. Defaults to 2.

        Returns:
            tuple[int,int,float]: hours, minutes, seconds
        """
        self.__end()
        return self.__parseDuration(decimals)

    def getTimeString(self, decimals:int=2) -> str:
        """gets the elapsed time as a string in hh:mm:ss.ms format

        Args:
            decimals (int, optional): the number of decimal points for the seconds. Defaults to 2.

        Returns:
            str: hh:mm:ss.ms
        """
        self.__end()
        return self.__getDurationString(decimals)
      
    def printTime(self, decimals:int=2) -> None:
        """prints the current time in hh:mm::ss.ms format

        Args:
            decimals (int, optional): the number of decimal points for the seconds. Defaults to 2.
        """
        self.__end()
        print(self.__getDurationString(decimals))

    def restart(self) -> None:
        """ restarts the clock object
        """
        self.__start()
    
    def printStart(self, msg:str, prefix:str='', end:str=' ... ', spin:bool=True,) -> None:
        """prints the start message and restarts the clock

        Args:
            msg (str): the message to print
            prefix (str, optional): the beginning of the printed message. Defaults to ''.
            end (str, optional): the end of the printed message. Defaults to ' ... '.
            spin (bool, optional): indicates if the wheel should spin. Defaults to True.
        """
        # save the spin status
        if not "CI" in os.environ:
            self.__spin = spin
        
        # do not spin if we are doing CI
        else:
            self.__spin = False
        
        # spin the wheel if requested; wheel handles printing
        if self.__spin:
            self.__WHEEL.start(prefix + msg + end)
        
        # otherwise print the message 
        else:
            print(prefix + msg, end=end)
            sys.stdout.flush()
        
        # restart the clock
        self.restart()
    
    def printDone(self) -> None:
        """prints the end message and the duration
        """
        # stop spinning the wheel if necessary
        Clock._killWheel()
        
        # print the done string
        print(f"done {self.getTimeString()}")
        sys.stdout.flush()
        
        Clock.__WHEEL.__msg = ''
    
    def _killWheel() -> None:
        """kills any spinning clocks
        """
        Clock.__WHEEL.stop()


# functions
def __printHelp() -> None:
    """prints the help message
    """
    # build the message
    GAP = " "*4
    EOL = "\n"
    MSG = f"{EOL}Calculates the Jaccard index for a set of fasta files{EOL}" + \
          f"{GAP}{__author__}, 2024{EOL*2}" + \
          f"Usage:{EOL}" + \
          f"{GAP}{os.path.basename(__file__)} fastas out klen cpus{EOL*2}" + \
          f"Positional arguments:{EOL}" + \
          f"{GAP}fastas{GAP}a list of fasta files (space separated) or a pattern (eg. path/*.fasta){EOL}" + \
          f"{GAP}out   {GAP}the output filename{EOL}" + \
          f"{GAP}klen  {GAP}the kmer length as an integer{EOL}" + \
          f"{GAP}cpus  {GAP}the number of cpus as an integer{EOL}"
    
    # print the message
    print(MSG)


def __pickleKmers(name:str, kmers:set[str]) -> None:
    """pickles a set of kmers

    Args:
        name (str): the name of the genome
        kmers (set[str]): the set of kmers
    """
    # determine the filename for the pickle
    fn = os.path.join(__PICKLE_DIR, name + __PICKLE_EXT)
    
    # dump the kmers to file
    with open(fn, 'wb') as fh:
        pickle.dump(kmers, fh)
        

def __getKmers(seq:str, k:int, name:str) -> None:
    """gets a set of kmers from a given sequence and pickles them

    Args:
        seq (str): the genomic sequence
        k (int): the kmer length
        name (str): the name of the sequence
    """
    # create a countgraph and use it to retrieve a set of kmers
    kh = Countgraph(k, 1e7, 1)
    
    # save the kmers
    __pickleKmers(name, {x for x in kh.get_kmers(seq) if not __JUNCTION_CHAR in x})


def __importSequences(files:list[str], k:int, cpus:int) -> list[str]:
    """imports kmer sequences from files in parallel

    Args:
        files (list[str]): a list of sequence files in fasta format
        k (int): the length of the kmers
        cpus (int): the number of cpus for parallel processing

    Returns:
        dict[str,set[str]]: key=genome name; val=set of kmers
    """
    # constant
    FORMAT = 'fasta'
    
    # helper function
    def genArgs(allNames:list) -> Generator[tuple[str,int,str],None,None]:
        """generates arguments for __getKmers
        """
        # determine the junction
        junction = __JUNCTION_CHAR * k
        
        # for each file
        for fn in files:
            # combine the contig sequences
            seq = junction.join(map(str, (x.seq for x in SeqIO.parse(fn, FORMAT))))
            
            # get the name and save it
            name = os.path.splitext(os.path.basename(fn))[0]
            allNames.append(name)
            
            # generate arguments
            yield (seq, k, name)
    
    # initialize output
    out = list()
    
    # make sure pickle directory exists
    if not os.path.isdir(__PICKLE_DIR):
        os.mkdir(__PICKLE_DIR)
    
    # get the kmers in parallel
    pool = Pool(cpus)
    pool.starmap(__getKmers, genArgs(out))
    pool.close()
    pool.join()
    
    return out


def __unpickleKmers(name:str) -> set[str]:
    """unpickles kmer file for a genome name

    Args:
        name (str): the name of the genome

    Returns:
        set[str]: the set of kmers
    """
    # get the name of the pickle file
    fn = os.path.join(__PICKLE_DIR, name + __PICKLE_EXT)
    
    # unpickle the file
    with open(fn, 'rb') as fh:
        out = pickle.load(fh)
    
    return out


def __writeTable(fn:str, terminator:str, queue:Queue) -> None:
    """writes jaccard index values in parallel

    Args:
        fn (str): the filename to write
        terminator (str): the kill signal
        queue (Queue): a queue to communicate with child processes
    """
    # constants
    SEP = "\t"
    EOL = "\n"
    
    # open the file
    with open(fn, 'w') as fh:
        # keep writing until the kill signal is received
        while True:
            # extract data from the queue
            val = queue.get()
            
            # check for kill signal
            if val == terminator:
                break
            
            # write the data from the queue to file
            fh.write(SEP.join(map(str, val)) + EOL)
            fh.flush()


def __calculateOneJaccard(name1:str, name2:str, queue:Queue) -> None:
    """calculates the jaccard index for a pair of genomes

    Args:
        name1 (str): one genome
        name2 (str): a second genome
        queue (Queue): a queue to communicate with the writer process
    """
    # unpickle the kmers for these two genomes
    set1 = __unpickleKmers(name1)
    set2 = __unpickleKmers(name2)
    
    # determine the intersection and union
    inter = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    # save the names and the jaccard index for these genomes
    queue.put((name1, name2, (inter / union)))


def __calculateAllJaccards(names:list[str], fn:str, cpus:int) -> None:
    """calculates all pairwise jaccard indices

    Args:
        names (list[str]): a list of genome names
        fn (str): the filename where jaccard indices will be written
        cpus (int): the number of cpus for parallel processing
    """
    # constant
    TERMINATOR = "kill"
    
    # helper function
    def genArgs(q) -> Generator[tuple[str,str,set[str],set[str]],None,None]:
        """generates arguments for __calculateOneJaccard
        """
        for name1,name2 in combinations(names, 2):
            yield name1, name2, q
    
    # create a queue and a pool
    queue = Manager().Queue()
    pool = Pool(cpus)
    
    # start the writer process
    pool.apply_async(__writeTable, (fn, TERMINATOR, queue))
    
    # calculate jaccard indices in parallel
    pool.starmap(__calculateOneJaccard, genArgs(queue))
    
    # send the kill signal to the writer process and close the pool
    queue.put(TERMINATOR)
    pool.close()
    pool.join()


def __importTableAsMatrix(fn:str) -> dict[tuple[str,str],float]:
    """imports the jaccard index table as a matrix

    Args:
        fn (str): the filename of the jaccard index table

    Returns:
        dict[tuple[str,str],float]: key=pairs of genomes; val=jaccard index
    """
    # constants
    SEP = "\t"
    SELF_JI = 1.0
    
    # initialize output
    out = dict()
    
    # go through the file line by line
    with open(fn, 'r') as fh:
        for line in fh:
            n1,n2,ji = line.rstrip().split(SEP)
            
            # create a square matrix for this line
            out[(n1,n2)] = ji
            out[(n2,n1)] = ji
    
    
    # set the self-vs-self values
    for n in {x[0] for x in out.keys()}:
        out[(n,n)] = SELF_JI
    
    return out


def __writeJaccardMatrix(jiMat:dict[tuple[str,str],float], fn:str) -> None:
    """writes a square matrix of percent identities to file

    Args:
        jiMat (dict[tuple[str,str],float]): the dictionary produced by __calculateAllJaccards
        fn (str): the filename where the data will be written
    """
    # constants
    SEP = "\t"
    EOL = "\n"
    
    # determine the order that the names will be written
    names = sorted({x[0] for x in jiMat.keys()})
    
    # open the file and write the header
    with open(fn, 'w') as fh:
        fh.write(SEP + SEP.join(names) + EOL)
        
        # for each genome
        for n1 in names:
            # write the row name
            fh.write(n1)
            
            # write the values for the comparison of this genome versus the others
            for n2 in names:
                fh.write(SEP + str(jiMat[(n1,n2)]))
            
            # end the current row
            fh.write(EOL)


def __main() -> None:
    """main runner function
    """
    # constant
    TEMP_FN = os.path.join(os.curdir, '.tmp')
    
    # handle no arguments
    if len(sys.argv) == 1:
        __printHelp()
    
    # handle all other conditions
    else:
        try:
            fnaFiles = sys.argv[1:-3]
            outFn  = sys.argv[-3]
            klen = int(sys.argv[-2])
            cpus = int(sys.argv[-1])
        
        # print the help and quit if anything fails
        except:
            __printHelp()
            return

        # create clocks
        clock = Clock()
        total = Clock()
        
        # import kmers
        clock.printStart('importing sequence data', spin=False)
        names = __importSequences(fnaFiles, klen, cpus)
        clock.printDone()

        # calculate the jaccard indices for all pairwise comparisons of genomes
        clock.printStart('calculating jaccard indices', spin=False)
        __calculateAllJaccards(names, TEMP_FN, cpus)
        clock.printDone()
        
        # convert the table to a matrix and write it
        clock.printStart('converting jaccard index table to matrix', spin=False)
        matrix = __importTableAsMatrix(TEMP_FN)
        __writeJaccardMatrix(matrix, outFn)
        clock.printDone()
        
        # delete temporary files
        clock.printStart('deleting intermediate files', spin=False)
        shutil.rmtree(__PICKLE_DIR)
        os.remove(TEMP_FN)
        clock.printDone()
        
        # print the total runtime
        print(f'\ntotal time: {total.getTimeString()}\n')


# entrypoint
if __name__ == "__main__":
    __main()
