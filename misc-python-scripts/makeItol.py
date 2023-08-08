# Joseph S. Wirth, 2023
import getopt, os, rpy2.robjects, sys

# global variables
JOB_0 = "help"
JOB_1 = "colorstrip"
JOB_2 = "binary"

def __main():
    """main runner function
    """
    # parse command line arguments
    task,inFN,outFN,label,delim,helpRequested = __parseArgs()
    
    # do nothing if help was requested
    if not helpRequested:
        # color strip
        if task == JOB_1:
            dataD = __parseCsv(inFN, JOB_1)
            __writeColorStrip(outFN, dataD, label, delim)
            
        # binary
        elif task == JOB_2:
            dataD = __parseCsv(inFN, JOB_2)
            __writeBinary(outFN, dataD, label, delim)


def __parseArgs() -> tuple[str,str,str,str,str,bool]:
    """parses command line arguments

    Raises:
        ValueError: invalid task specified
        ValueError: invalid or missing input file
        ValueError: invalid delimiter
        ValueError: missing one or more required arguments

    Returns:
        tuple[str,str,str,str,str,bool]: task; infile; outfile; data label; delim char; help requested?
    """
    # constants
    ALLOWED_DELIM = ("COMMA", "TAB", "SPACE")
    DELIM_D = {"COMMA": ",",
               "SPACE": " ",
               "TAB": "\t"}
    
    # flags
    IN_FLAGS = ("-i", "--in")
    OUT_FLAGS = ("-o", "--out")
    LABEL_FLAGS = ("-L", "--label")
    DELIM_FLAGS = ("-d", "--delim")
    HELP_FLAGS = ("-h", "--help")
    SHORT_OPTS = IN_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 LABEL_FLAGS[0][-1] + ":" + \
                 DELIM_FLAGS[0][-1] + ":" + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (IN_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + '=',
                 LABEL_FLAGS[1][2:] + "=",
                 DELIM_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:])
    
    # messages
    IGNORE_MSG = "ignoring unused argument: "
    ERR_MSG_0 = "invalid task: "
    ERR_MSG_1 = "invalid (or missing) input file: "
    ERR_MSG_2 = "invalid delimiter specified: "
    ERR_MSG_3 = "must specify all required arguments"
    
    # helper function to retrieve the appropriate help message
    def getHelpMessage(task:str) -> str:
        GAP = 4*' '
        EOL = "\n"
        SEP = ", "
        HEADER = EOL + "makes an iToL formatted dataset from an input text file" + EOL + \
                 GAP + "Joseph S. Wirth, 2023" + EOL*2 + "usage:" + EOL
        HELP_0 = GAP + "python3.11 makeItol.py TASK [OPTIONS]" + EOL * 2 + \
                 "Available tasks:" + EOL + \
                 GAP + f"{JOB_1:<14}{'makes a color strip dataset iToL file':<}" + EOL + \
                 GAP + f"{JOB_2:<14}{'makes a binary dataset iToL file':<}" + EOL + \
                 GAP + f"{JOB_0:<14}{'print this message':<}" + EOL*2 + \
                 "Run 'python3.11 makeItol.py TASK --help' for help with a specific task." + EOL
        HELP_1 = GAP + "python3.11 makeItol.py " + JOB_1 + " [-ioLdh]" + EOL*2 + \
                 "required arguments:" + EOL + \
                 GAP + f"{IN_FLAGS[0] + SEP + IN_FLAGS[1]:<14}{'[file] a csv containing a header and exactly two columns (name, value)'}" + EOL + \
                 GAP + f"{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<14}{'[file] the output filename'}" + EOL + \
                 GAP + f"{LABEL_FLAGS[0] + SEP + LABEL_FLAGS[1]:<14}{'[str] a label for the dataset'}" + EOL*2 + \
                 "optional arguments:" + EOL + \
                 GAP + f"{DELIM_FLAGS[0] + SEP + DELIM_FLAGS[1]:<14}{'[str] the delimiter for the output file [COMMA, SPACE, TAB] (default: COMMA)'}" + EOL + \
                 GAP + f"{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<14}{'print this message'}" + EOL
        HELP_2 = GAP + "python3.11 makeItol.py " + JOB_2 + " [-ioLdh]" + EOL*2 + \
                 "required arguments:" + EOL + \
                 GAP + f"{IN_FLAGS[0] + SEP + IN_FLAGS[1]:<14}{'[file] a csv containing a header and two or more columns (name, val1, val2, etc.)'}" + EOL + \
                 GAP + f"{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<14}{'[file] the output filename'}" + EOL + \
                 GAP + f"{LABEL_FLAGS[0] + SEP + LABEL_FLAGS[1]:<14}{'[str] a label for the dataset'}" + EOL*2 + \
                 "optional arguments:" + EOL + \
                 GAP + f"{DELIM_FLAGS[0] + SEP + DELIM_FLAGS[1]:<14}{'[str] the delimiter for the output file [COMMA, SPACE, TAB] (default: COMMA)'}" + EOL + \
                 GAP + f"{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<14}{'print this message'}" + EOL
        
        if task == JOB_0:
            return HEADER + HELP_0
    
        elif task == JOB_1:
            return HEADER + HELP_1
        
        elif task == JOB_2:
            return HEADER + HELP_2
    
    # set default values
    task = None
    inFN = None
    outFN = None
    label = None
    delim = ","
    helpRequested = False

    # handle no task
    if len(sys.argv) == 1:
        helpRequested = True
        print(getHelpMessage(JOB_0))
    
    else:
        # extract the task; make sure it is valid
        task = sys.argv[1]
        
        if task not in (JOB_0,JOB_1,JOB_2):
            raise ValueError(ERR_MSG_0 + task)
    
        if task == JOB_0:
            helpRequested = True
            print(getHelpMessage(task))
        
        # handle task provided; no arguments
        elif len(sys.argv) == 2:
            helpRequested = True
            print(getHelpMessage(task))
        
        else:
            # parse the command line arguments
            opts,args = getopt.getopt(sys.argv[2:], SHORT_OPTS, LONG_OPTS)
            
            for opt,arg in opts:
                # extract the input file
                if opt in IN_FLAGS:
                    if not os.path.isfile(arg):
                        raise ValueError(ERR_MSG_1 + arg)
                    inFN = os.path.abspath(arg)
                
                # extract the output file
                elif opt in OUT_FLAGS:
                    outFN = os.path.abspath(arg)
                
                # get the label for the dataset
                elif opt in LABEL_FLAGS:
                    label = arg
                
                # get the delim for the output file
                elif opt in DELIM_FLAGS:
                    if arg not in ALLOWED_DELIM:
                        raise ValueError(ERR_MSG_2 + arg)
                    delim = DELIM_D[arg]
            
                elif opt in HELP_FLAGS:
                    helpRequested = True
                    print(getHelpMessage(task))
            
                # ignore all other flags
                else:
                    print(IGNORE_MSG + opt + " " + arg)
            
            # make sure all required arguments were specified
            if None in (inFN, outFN, label) and not helpRequested:
                raise ValueError(ERR_MSG_3)
    
    return task,inFN,outFN,label,delim,helpRequested


def __parseCsv(fn:str, task:str) -> dict:
    """parses a CSV file for input to write functions

    Args:
        fn (str): filename to be parsed
        task (str): the task (colorstrip or binary)

    Returns:
        dict: key = name
                if colorstrip; val = value
                if binary; val = dict: key = column name; val = value
    """
    # constants
    SEP = ","
    EOL = "\n"
    
    # initialize looping vars
    outD = dict()
    header = True
    with open(fn, 'r') as fh:
        for line in fh:
            # chomp new line characters
            if line[-1] == EOL:
                line = line[:-1]
            
            # split the line into a list of columns
            row = line.split(SEP)
            
            # link indices to header names (used in binary)
            if header:
                indexD = dict()
                for idx in range(1, len(row)):
                    indexD[idx] = row[idx]
                header = False
            
            # if not the header
            else:
                if line != '':
                    # extract the tip name
                    name = row[0]
                    
                    # make the colorstrip dictionary
                    if task == JOB_1:
                        outD[name] = row[1]
                    
                    # make the binary dictionary
                    elif task == JOB_2:
                        outD[name] = dict()
                        for idx in range(1, len(row)):
                            outD[name][indexD[idx]] = str(int(row[idx]))
    
    return outD
  

def __writeColorStrip(fn:str, dataD:dict, label:str, delim:str) -> None:
    """writes the color strip file

    Args:
        fn (str): filename to be written
        dataD (dict): dictionary produced by __parseCsv (with colorstrip task)
        label (str): the label for the file
        delim (str): the delimiter character
    """
    # constants
    EOL = "\n"
    DELIM_D = {",": "COMMA",
               "\t": "TAB",
               " ": "SPACE"}
    
    # header strings
    HEADER_1 = "DATASET_COLORSTRIP" + EOL + "SEPARATOR "
    HEADER_2 = EOL + "DATASET_LABEL"
    HEADER_3 = EOL + "COLOR"
    HEADER_4 = EOL + "COLOR_BRANCHES"
    HEADER_5 = "0" + EOL + "DATA" + EOL
    
    # get a dictionary of values with their respective colors
    colorD = __getColorStripColors(dataD)
    
    # pick a random color as the dataset color
    for col in colorD.values():
        break
    
    with open(fn, 'w') as fh:
        # write the header string
        fh.write(HEADER_1)
        fh.write(DELIM_D[delim])
        fh.write(HEADER_2 + delim + label)
        fh.write(HEADER_3 + delim + col)
        fh.write(HEADER_4 + delim)
        fh.write(HEADER_5)
        
        # write the data
        for name in dataD.keys():
            fh.write(name + delim)
            fh.write(colorD[dataD[name]] + delim)
            fh.write(dataD[name] + EOL)


def __writeBinary(fn:str, dataD:dict, label:str, delim:str) -> None:
    """writes the binary file

    Args:
        fn (str): filename to be written
        dataD (dict): dictionary produced by __parseCsv (binary task)
        label (str): dataset label
        delim (str): delimiter character
    """
    # constants
    EOL = "\n"
    DELIM_D = {",": "COMMA",
               "\t": "TAB",
               " ": "SPACE"}
    
    # header strings
    HEADER_1 = "DATASET_BINARY" + EOL + "SEPARATOR "
    HEADER_2 = EOL + "DATASET_LABEL"
    HEADER_3 = EOL + "COLOR"
    HEADER_4 = EOL + "FIELD_SHAPES"
    HEADER_5 = EOL + "FIELD_LABELS"
    HEADER_6 = EOL + "FIELD_COLORS"
    HEADER_7 = EOL + "DATA" + EOL + "#name"
    
    # get a dictionary of fields with their respective colors
    colorD = __getColorBinary(dataD)
    
    # get the field colors, field labels, field shapes, and the field order
    colorStr,fieldStr,shapeStr,fieldOrder = __getBinaryFieldData(colorD, delim)

    # pick a random color as the dataset color
    for col in colorD.values():
        break
    
    with open(fn, 'w') as fh:
        # write the header string
        fh.write(HEADER_1 + DELIM_D[delim])
        fh.write(HEADER_2 + delim + label)
        fh.write(HEADER_3 + delim + col)
        fh.write(HEADER_4 + shapeStr)
        fh.write(HEADER_5 + fieldStr)
        fh.write(HEADER_6 + colorStr)
        fh.write(HEADER_7 + fieldStr + EOL)
        
        # write the data
        for name in dataD.keys():
            fh.write(name)
            for field in fieldOrder:
                fh.write(delim + dataD[name][field])
            fh.write(EOL)


def __getBinaryFieldData(colorD:dict, delim:str) -> tuple[str,str,str,list[str]]:
    """retrieves the field data for binary file

    Args:
        colorD (dict): the dictionary produced by __getColorBinary
        delim (str): the delimiter character

    Returns:
        tuple[str,str,str,list[str]]: color string, field string, shape string, field order (list)
    """
    # constant
    SHAPE = "1"
    
    # initialize outputs
    fieldOrder = list()
    colorStr = ""
    shapeStr = ""
    fieldStr = ""
    
    # go through each field in the color dictionary
    for field in colorD.keys():
        # keep track of the order
        fieldOrder.append(field)
        
        # make the color, shape, and field strings
        colorStr += delim + colorD[field]
        shapeStr += delim + SHAPE
        fieldStr += delim + field
    
    return colorStr,fieldStr,shapeStr,fieldOrder


def __getColorStripColors(dataD:dict) -> dict:
    """gets a dictionary of colors for a color strip

    Args:
        dataD (dict): the dictionary produced by __parseCsv (colorstrip task)

    Returns:
        dict: key = field; val = color
    """
    # get all the unique values and an equal number of colors
    vals = set(dataD.values())
    cols = __pickColors(len(vals))
    
    # key each value with a color
    outD = dict()
    for val in vals:
        outD[val] = cols.pop()
    
    return outD


def __getColorBinary(dataD:dict) -> dict:
    """gets a dictionary of colors for a binary

    Args:
        dataD (dict): the dictionary produced by __parseCsv (binary task)

    Returns:
        dict: key = field; val = hex color
    """
    # get all unique values and an equal number of colors
    vals = {y for x in dataD.keys() for y in dataD[x].keys()}
    cols = __pickColors(len(vals))
    
    # key each value with a color
    outD = dict()
    for val in vals:
        outD[val] = cols.pop()

    return outD


def __pickColors(num:int) -> set[str]:
    """gets a set of num colors

    Args:
        num (int): the number of colors to pick

    Returns:
        set[str]: a set of hex colors
    """
    # use R's `rainbow` function to pick colors
    rainbow = rpy2.robjects.r['rainbow']
    return set(rainbow(num))


if __name__ == "__main__":
    __main()