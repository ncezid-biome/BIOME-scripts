# BIOME-scripts
random scripts

## misc-python-scripts
### creating and activating a `conda` environment to run these scripts
First create the `conda` environment:

    conda env create -f BIOME-scripts/misc-python-scripts/environment.yml

Then activate the `conda` environment to ensure that all the dependencies are available:

    conda activate misc-python

### downloadAssemblies.py
#### Joseph S. Wirth, 2023
Downloads assemblies from a text file containing one accession number for NCBI Assembly per line. Can optionally specify an output directory (default: current directory) or have the files be renamed to match the accession numbers in the text file.

### downloadSRA.py
#### Joseph S. Wirth, 2023
Downloads reads from NCBI **_in parallel_** using `sratoolkit version 2.11.3`. Takes in a text file with one SRR id per line. Can optionally specify an output directory (default: current directory). To use the parallel feature, specify a maximum number of threads to use (`-t`, `--max_threads`) and the number of threads to use for each parallel process (`-p`, `threads_per_download`). An example for using this program with `qsub` can be found in the `qsub_arrays` folder

### extractLocusTagsFromGenbank.py
#### Joseph S. Wirth, 2023
Extracts a locus tag(s) from a genbank file and writes a new fasta file to stdout. Must provide a genbank file and one or more locus tags. If specifying multiple locus tags, then provide them as a comma-separated list. Can optionally specify an output file as the destination. Can optionally specify a molecule type (default: nucleotide) for the output file.

### extractRegionFromGenbank.py
#### Joseph S. Wirth, 2023
Extracts a region(s) from a genbank file and writes a new file to stdout. Must input a genbank file and one or more contigs. If specifying multiple contigs, then provide them as a comma-separated list. Can optionally specify an output file as the destination. Can optionally specify start and/or end coordinates for a contig; this only works when extracting a single contig. Can optionally generate a file in fasta format instead of genbank format.

### genbankToFasta.py
#### Joseph S. Wirth, 2023
Converts a genbank file to a nucleotide (fna) fasta for an amino acid (faa) fasta file. Requires a genbank file as input as well as the desired output format (nucleotide or protein). If nucleotide (`nuc`) is selected, the fasta will consist of the sequences for all the contigs in the input file. If amino acid (`prot`) is selected, the fasta will consist of the protein sequences for all annotated CDS in the input file. Can optionally specify an output file (by default the output is the same filename with a different extension).

### makeItol.py
#### Joseph S. Wirth, 2023
Generates text files formatted for annotating trees using the [interactive Tree of Life (iToL)](https://itol.embl.de/). Can generate files for the `color strip` and `binary` data sets. Requires an output filename and a label for the dataset. If making a `color strip`, requires a csv with a header and two columns: column 1 = tip name and column 2 = value. If making a `binary`, requires a csv with a header and two or more columns: column 1 = tip name and subsequent columns = field names; the values for the subsequent columns must be `0` (absent), `1` (present), or `-1` (do not plot). Can optionally specify a delimiter character; by default it is `COMMA`.

### primerDesign.py
#### Joseph S. Wirth, 2023
Finds pairs of primers that can be used based with a given input genome independent of a PCR product. First, it extracts all kmers of a given size range that appear only once in the genome, have appropriate GC% and Tm values, do not contain long (>3bp) repeats, do not form hairpins, and have a G or C at the 3' end. Next, it compares primer pairs that can generate a specified PCR product size and saves only those primers who do not form primer dimers and whose GC% differences and Tm differences fall within the specified ranges. Writes the pairs to specified output file along with information about the position and PCR product lengths.
