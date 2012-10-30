/*
PRICE was written by Graham Ruby.
This file is part of PRICE.

PRICE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PRICE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PRICE.  If not, see <http://www.gnu.org/licenses/>.
*/

/* This is the text-based interface for the PRICE assembler.
 *
 */


#include "Assembler.h"

#include <set>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <string.h>
#include <bitset>
#include <cstddef>
#include <omp.h>
using namespace::std;
using ::size_t;

string helpString = string(
"PRICE v0.18\n"
"\n"
"These are the command options for the PRICE assembler. \n"
"For more details, see the accompanying README.txt file. \n"
"Contact: price@derisilab.ucsf.edu \n"
"\n"
"Usage: ./PriceTI [args] \n"
"\n"
"INPUT FILES: \n"
" accepted formats are fasta (appended .fa, .fasta, .fna, .ffn, .frn), fastq (.fq, .fastq, or _sequence.txt), or priceq (.pq or .priceq) \n"
"  NOTE ABOUT FASTQ ENCODING: multiple encodings are currently used for fastq quality scores.  The traditional encoding is Phred+33,\n"
"                             and PRICE will interpret scores from any .fq or .fastq file according to that encoding.  The Phred+64 \n"
"                             encoding has been used extensively by Illumina, and so it is applied to Illumina's commonly-used _sequence.txt\n"
"                             file append.  Please make sure that your encoding matches your file append.\n"
"INPUT READ FILES: \n"
"  NOTE: these flags can be used multiple times in the same command to include multiple read datasets. \n"
"  (default % ID = 90%) \n"
" PAIRED-END FILES (reads are 3p of one another on opposite strands, i.e. pointing towards one another)\n"
" -fp a b c [d e [f]]: (a,b)input file pair, (c)amplicon insert size (including read) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" -fpp a b c d [e f [g]]: (a,b)input file pair, (c)amplicon insert size (including read), (d)required % identity for match (25-100 allowed) \n"
"                  (e,f,g) are optional; (e)the num. cycles to be skipped before this file is used;\n"
"                  if (g) is provided, then the file will alternate between being used for (f) cycles and not used for (g) cycles;\n"
"                  otherwise, the file will be used for (f) cycles then will not be used again. \n"
" -fs a b [c d [e]]: (a)input paired-end file (alternating entries are paired-end reads), (b)amplicon insert size (including read) \n"
"                  (c,d,e) are optional; (c)the num. cycles to be skipped before this file is used;\n"
"                  if (e) is provided, then the file will alternate between being used for (d) cycles and not used for (e) cycles;\n"
"                  otherwise, the file will be used for (d) cycles then will not be used again. \n"
" -fsp a b c [d e [f]]: (a)input paired-end file (alternating entries are paired-end reads), (b)amplicon insert size (including read),\n"
"                  (c)required % identity for match (25-100 allowed) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" MATE-PAIR FILES (reads are 5p of one another on opposite strands, i.e. pointing away from one another)\n"
" -mp a b c [d e [f]]: like -fp above, but with reads in the opposite orientation.\n"
" -mpp a b c d [e f [g]]: like -fpp above, but with reads in the opposite orientation.\n"
" -ms a b [c d [e]]: like -fs above, but with reads in the opposite orientation.\n"
" -msp a b c [d e [f]]: like -fsp above, but with reads in the opposite orientation.\n"
" FALSE PAIRED-END FILES (unpaired reads are split into paired ends, with the scores of double-use nuceotides halved)\n"
" -spf a b c [d e [f]]: (a)input file, (b)the length of the 'reads' that will be taken from each side of the input reads,\n"
"                  (c)amplicon insert size (including read) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" -spfp a b c d [e f [g]]: (a)input file, (b)the length of the 'reads' that will be taken from each side of the input reads,\n"
"                  (c)amplicon insert size (including read), (d)required % identity for match (25-100 allowed) \n"
"                  (e,f,g) are optional; (e)the num. cycles to be skipped before this file is used;\n"
"                  if (g) is provided, then the file will alternate between being used for (f) cycles and not used for (g) cycles;\n"
"                  otherwise, the file will be used for (f) cycles then will not be used again. \n"
"INPUT INITIAL CONTIG FILES: \n"
"  NOTE: these flags can be used multiple times in the same command to include multiple initial contig datasets. \n"
" -icf a b c d: (a)initial contig file, (b)number of addition steps, (c)number of cycles per step,\n"
"               (d)const by which to multiply quality scores \n"
" -picf a b c d e: (a)num of initial contigs from this file, (b)initial contig file, (c)num addition steps,\n"
"                  (d)num cycles per step (e)const by which to multiply quality scores \n"
" -icfNt/-picfNt: same as -icf/-picf, but if target mode is invoked, contigs with matches to these input sequences will not necessarily\n"
"                 be retained \n"
"OUTPUT FILES: \n"
" accepted formats are fasta (.fa or .fasta) or priceq (.pq or .priceq) \n"
" -o a: (a)output file name (.fasta or .priceq) \n"
" -nco a: (a)num. cycles that pass in between output files being written (default=1) \n"
"OTHER PARAMS: \n"
" -nc a: (a)num. of cycles \n"
" -link a: (a)max. number of contigs that are allowed to replace a read in a mini-assembly (default=2)\n"
" -mol a: (a)minimum overlap length for mini-assembly (default=35) \n"
" -tol a: (a)threshold seq num for scaling overlap for overhang assemblies (default=20)\n"
" NOTE: -mol and -tol do not affect the parameters for de-Bruijn-graph-based assembly.\n"
" -mpi a: (a)minimum % identity for mini-assembly (default=85) \n"
" -tpi a: (a)threshold seq num for scaling % ID for mini-assemblies (default=20)\n"
" -MPI, -TPI : same as above, but for meta-assembly (-MPI default=85, -TPI default=1000)\n"
" NOTE: there is no minimum overlap value for meta-assembly\n"
" -dbmax a: (a) the maximum length sequence that will be fed into de Bruijn assembly \n"
"           (default=100; recommended: max paired-end read length)\n"
" -dbk a: (a) the k-mer size for de Bruijn assembly (default=20; keep less than the read length)\n"
" -dbms a: (a) the minimum number of sequences to which de Bruijn assembly will be applied (default=3)\n"
" -r a: (a) alignment score reward for a nucleotide match; should be a positive integer (default=1)\n"
" -q a: (a) alignment score penalty for a nucleotide mismatch; should be a negative integer (default=-2)\n"
" -G a: (a) alignment score penalty for opening a gap; should be a negative integer (default=-5)\n"
" -E a: (a) alignment score penalty for extending a gap; should be a negative integer (default=-2)\n"
"FILTERING READS: \n"
" -rqf a b [c d]: filters pairs of reads if either has an unaccptably high number of low-quality nucleotides, as defined\n"
"                 by the provided quality scores (only applies to files whose formats include quality score information). \n"
"                 (a) the percentage of nucleotides in a read that must be high-quality; (b) the minimum allowed probability\n"
"                 of a nucleotide being correct (must be between 0 and 1, and will usually be a decimal value close to 1);\n"
"                 (c) and (d) optionally constrain this filter to use after (c) cycles have passed, to run for (d) cycles.\n"
"                 This flag may be called multiple times to generate variable behavior across a PRICE run.\n"
" -rnf a [b c]: filters pairs of reads if either has an unaccptably high number of uncalled nucleotides (Ns or other ambiguous\n"
"                 IUPAC codes).  Like -rqf, but will also filter fasta-format data.  (a) the percentage of nucleotides in a read\n"
"                 that must be called; (b) and (c) optionally constrain this filter to use after (b) cycles have passed, to run\n"
"                 for (c) cycles.  This flag may be called multiple times to generate variable behavior across a PRICE run.\n"
" -maxHp a: filters out a pair of reads if either read has a homo-polymer track >(a) nucleotides in length.\n"
" -maxDi a: filters out a pair of reads if either read has a repeating di-nucleotide track >(a) nucleotides in length.\n"
"           NOTE: this will also catch mono-nucleotide repeats of the specified length (a string of A's is also a string\n"
"           of AA's), so calling -maxHp in addition to -maxDi is superfluous unless -maxHp is given a smaller max value.\n"
" -badf a b: prevents reads with a match of at least (b)% identity to a sequence in file (a) from being mapped to contigs.\n"
" -repmask a b c d e f [g]: uses coverage levels of constructed and/or input contigs to find repetitive elements and mask them \n"
"                           as if they were sequences input using -badf.\n"
"                           (a) = cycle number (1-indexed) at which repeats will be detected.\n"
"                           (b) = 's' if repeats will be sought at the start of the cycle or 'f' if they will be sought at the finish.\n"
"                           (c) = the min. number of variance units above the median that will be counted as high-coverage.\n"
"                           (d) = the min. fold increase in coverage above the median that will be counted as high-coverage.\n"
"                           (e) = the min. size in nt for a detected repeat. \n"
"                           (f) = reads with a match of at least this % identity to a repeat will not be mapped to contigs.\n"
"                           (g) = an optional output file (.fasta or .priceq) to which the detected repeats will be written.\n"
" -reset a [b c d...]: re-introduces contigs that were previously not generating assembly jobs of their own\n"
"                      (a) is the one-indexed cycle where the contigs will be reset.  Same with b, c, d. \n"
"                      Any number of args may be added.\n"
"FILTERING INITIAL CONTIGS: \n"
" -icbf a b [c]: prevents input sequences with a match of at least (b)% identity to a sequence in file (a) from being used.\n"
"                This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
" -icmHp a [b]: filters out an initial contig if it has a homo-polymer track >(a) nucleotides in length.\n"
"               This filter is optionally not applied to sequences of length greater than (b) nucleotides.\n"
" -icmDi a [b]: filters out an initial contig if it has a repeating di-nucleotide track >(a) nucleotides in length.\n"
"               NOTE: this will also catch mono-nucleotide repeats of the specified length (a string of A's is also a string\n"
"               of AA's), so calling -icmHp in addition to -icmDi is superfluous unless -icmHp is given a smaller max value.\n"
"               This filter is optionally not applied to sequences of length greater than (b) nucleotides.\n"
" -icqf a b [c]: filters out an initial contig if it has an unaccptably high number of low-quality nucleotides, as defined\n"
"                by the provided quality scores (only applies to files whose formats include quality score information). \n"
"                (a) the percentage of nucleotides in a read that must be high-quality; (b) the minimum allowed probability\n"
"                of a nucleotide being correct (must be between 0 and 1, and will usually be a decimal value close to 1);\n"
"                This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
" -icnf a [b]: filters pairs of reads if either has an unaccptably high number of uncalled nucleotides (Ns or other ambiguous\n"
"              IUPAC codes).  Like -icqf, but will also filter fasta-format data.  (a) the percentage of nucleotides in a read\n"
"              that must be called. This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
"FILTERING/PROCESSING ASSEMBLED CONTIGS: \n"
" -lenf a b: filters out contigs shorter than (a) nt at the end of every cycle, after skipping (b) cycles.\n"
"            NOTE: multiple -lenf commands can be entered; for any cycle, the most recently-initiated filter is used.\n"
"            Example: -lenf 50 2 -lenf 300 4 -lenf 200 6 => no filter for the first two cycles, then a 50nt filter for cycles\n"
"                     3 & 4, then a 300nt filter for cycles 5 & 6, then a 200nt filter for cycles 7 onwards.\n"
" -trim a b [c]: at the end of the (a)th cycle (indexed from 1), trim off the edges of conigs until reaching the minimum coverage\n"
"                level (b), optionally deleting contigs shorter than (c) after trimming; this flag may be used repeatedly.\n"
" -trimB a b [c]: basal trim; after skipping (a) cycles, trim off the edges of conigs until reaching the minimum coverage\n"
"                 level (b) at the end of EVERY cycle, optionally deleting contigs shorter than (c) after trimming.\n"
"                 -trimB may be called many times, and multiple calls will interact in the same way as multiple -lenf calls\n"
"                 (explained above). A call to -trim will override the basal trim values for that specified cycle only.\n"
" -trimI a [b]: initial trim; input initial contigs are trimmed before being used by PRICE to seed assemblies. Contigs are\n"
"               trimmed from their outside edges until reaching the minimum coverage level (a), optionally deleting contigs\n"
"               shorter than (b) after trimming. This flag is most appropriate for .priceq input, can be appropriate for\n"
"               .fastq input, and is inappropriate for .fasta input.  It will be equally applied to ALL input contigs.\n"
" -target a b [c d]: limit output contigs to those with matches to input initial contigs at the end of each cycle.\n"
"                    (a) % identity to an input initial contig to count as a match (ungapped); (b)num cycles to skip\n"
"                    before applying this filter.  [c and d are optional, but must both be provided if either is]\n"
"                    After target filtering has begin, target-filtered/-unfiltered cycles will alternate with (c)\n"
"                    filtered cycles followed by (d) unfiltered cycles.\n"
" -targetF a b [c d]: the same as -target, but now matches to all reads in the input set will be specified, not just\n"
"                    the ones that have been introduced up to that point (this is FullFile mode).\n"
"COMPUTATIONAL EFFICIENCY: \n"
" -a x: (x)num threads to use (default=1) \n"
" -mtpf a: (a)max threads per file (default=1) \n"
"USER INTERFACE: \n"
" -log a: determines the type of outputmakes the output verbose (lots of time stamp tags) \n"
"         (a) = c: concise stdout (default)\n"
"         (a) = n: no stdout \n"
"         (a) = v: verbose stdout \n"
" -logf a: (a)the name of an output file for verbose log info to be written (doesn't change stdout format) \n"
" -, -h, or --help: user interface info. \n"
);

int main(int argc, char *argv[]){

  int initialInput = 0;
  bool willDoAssembly = true; // allows full set of errors to be printed
  bool distributeInitialContigs = false;
  int cyclesPerStep;
  int numSteps;

  // default is 1 core
  omp_set_num_threads( 1 );
  int maxThreadsPerFile = 1;

  bool listenerIsNull = false;
  bool listenerIsVerbose = false;
  bool includeListenerFile = false;
  string listenerFilename = "";

  // first, get the total number of cycles and create the cycle manager
  int numOfCyclesToDo;
  bool numCyclesSpecified = false;
  int argn = 1;
  while ( argn < argc ){
    string arg = string(argv[argn]);
    if (arg == "-h" or arg == "-" or arg == "--help"){
      cerr << helpString << endl;
      willDoAssembly = false;
      argn++;
    } else if (arg == "-nc"){
      numOfCyclesToDo = atoi( argv[argn+1] );
      numCyclesSpecified = true;
      argn += 2;
    } else { argn++; }
  }

  if (! numCyclesSpecified){
    cerr << "number of cycles must be specified with -nc flag" << endl;
    willDoAssembly = false;
  }


  if (willDoAssembly){

    Assembler* assembler = new Assembler(numOfCyclesToDo);
    assembler->setMetaFractIdContigThreshold(1000);

    // now go through again and get the command-line args organized
    argn = 1;
    while ( argn < argc ){
      string arg = string(argv[argn]);
      if (arg == "-fp"){
	string fileA = argv[argn+1];
	string fileB = argv[argn+2];
	long ampSize = atol(argv[argn+3]);
	if ( argn + 4 >= argc or argv[argn + 4][0] == '-'){
	  assembler->addReadFile(fileA,fileB,ampSize);
	  argn += 4;
	} else if ( argn + 6 >= argc or argv[argn + 6][0] == '-'){
	  int startCycle = atoi(argv[argn+4]);
	  int numCycles = atoi(argv[argn+5]);
	  assembler->addReadFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize);
	  argn += 6;
	} else {
	  int startCycle = atoi(argv[argn+4]);
	  int onCycles = atoi(argv[argn+5]);
	  int skipCycles = atoi(argv[argn+6]);
	  assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize);
	  argn += 7;
	}
      } else if (arg == "-fpp"){
	string fileA = argv[argn+1];
	string fileB = argv[argn+2];
	long ampSize = atol(argv[argn+3]);

	// helps deal with the optional last two args
	float fractId = float(atoi(argv[argn+4])) / float(100);
	if ( argn + 5 >= argc or argv[argn + 5][0] == '-'){
	  assembler->addReadFile(fileA,fileB,ampSize,fractId);
	  argn += 5;
	} else if ( argn + 7 >= argc or argv[argn + 7][0] == '-'){
	  int startCycle = atoi(argv[argn+5]);
	  int numCycles = atoi(argv[argn+6]);
	  assembler->addReadFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize);
	  argn += 7;
	} else {
	  int startCycle = atoi(argv[argn+5]);
	  int onCycles = atoi(argv[argn+6]);
	  int skipCycles = atoi(argv[argn+7]);
	  assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize);
	  argn += 8;
	}
      } else if (arg == "-fs"){
	string file = argv[argn+1];
	long ampSize = atol(argv[argn+2]);
	if ( argn + 3 >= argc or argv[argn + 3][0] == '-'){
	  assembler->addReadFile(file,ampSize);
	  argn += 3;
	} else if ( argn + 5 >= argc or argv[argn + 5][0] == '-'){
	  int startCycle = atoi(argv[argn+3]);
	  int numCycles = atoi(argv[argn+4]);
	  assembler->addReadFileLimitUse(startCycle,numCycles,file,ampSize);
	  argn += 5;
	} else {
	  int startCycle = atoi(argv[argn+3]);
	  int onCycles = atoi(argv[argn+4]);
	  int skipCycles = atoi(argv[argn+5]);
	  assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize);
	  argn += 6;
	}
      } else if (arg == "-fsp"){
	string file = argv[argn+1];
	long ampSize = atol(argv[argn+2]);

	// helps deal with the optional last two args
	float fractId = float(atoi(argv[argn+3])) / float(100);
	if ( argn + 4 >= argc or argv[argn + 4][0] == '-'){
	  assembler->addReadFile(file,ampSize,fractId);
	  argn += 4;
	} else if ( argn + 6 >= argc or argv[argn + 6][0] == '-'){
	  int startCycle = atoi(argv[argn+4]);
	  int numCycles = atoi(argv[argn+5]);
	  assembler->addReadFileLimitUse(startCycle,numCycles,file,ampSize);
	  argn += 6;
	} else {
	  int startCycle = atoi(argv[argn+4]);
	  int onCycles = atoi(argv[argn+5]);
	  int skipCycles = atoi(argv[argn+6]);
	  assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize);
	  argn += 7;
	}
      } else if (arg == "-mp"){
	string fileA = argv[argn+1];
	string fileB = argv[argn+2];
	long ampSize = atol(argv[argn+3]);
	if ( argn + 4 >= argc or argv[argn + 4][0] == '-'){
	  assembler->addMatePairFile(fileA,fileB,ampSize);
	  argn += 4;
	} else if ( argn + 6 >= argc or argv[argn + 6][0] == '-'){
	  int startCycle = atoi(argv[argn+4]);
	  int numCycles = atoi(argv[argn+5]);
	  assembler->addMatePairFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize);
	  argn += 6;
	} else {
	  int startCycle = atoi(argv[argn+4]);
	  int onCycles = atoi(argv[argn+5]);
	  int skipCycles = atoi(argv[argn+6]);
	  assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize);
	  argn += 7;
	}
      } else if (arg == "-mpp"){
	string fileA = argv[argn+1];
	string fileB = argv[argn+2];
	long ampSize = atol(argv[argn+3]);

	// helps deal with the optional last two args
	float fractId = float(atoi(argv[argn+4])) / float(100);
	if ( argn + 5 >= argc or argv[argn + 5][0] == '-'){
	  assembler->addMatePairFile(fileA,fileB,ampSize,fractId);
	  argn += 5;
	} else if ( argn + 7 >= argc or argv[argn + 7][0] == '-'){
	  int startCycle = atoi(argv[argn+5]);
	  int numCycles = atoi(argv[argn+6]);
	  assembler->addMatePairFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize);
	  argn += 7;
	} else {
	  int startCycle = atoi(argv[argn+5]);
	  int onCycles = atoi(argv[argn+6]);
	  int skipCycles = atoi(argv[argn+7]);
	  assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize);
	  argn += 8;
	}
      } else if (arg == "-ms"){
	string file = argv[argn+1];
	long ampSize = atol(argv[argn+2]);
	if ( argn + 3 >= argc or argv[argn + 3][0] == '-'){
	  assembler->addMatePairFile(file,ampSize);
	  argn += 3;
	} else if ( argn + 5 >= argc or argv[argn + 5][0] == '-'){
	  int startCycle = atoi(argv[argn+3]);
	  int numCycles = atoi(argv[argn+4]);
	  assembler->addMatePairFileLimitUse(startCycle,numCycles,file,ampSize);
	  argn += 5;
	} else {
	  int startCycle = atoi(argv[argn+3]);
	  int onCycles = atoi(argv[argn+4]);
	  int skipCycles = atoi(argv[argn+5]);
	  assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize);
	  argn += 6;
	}
      } else if (arg == "-msp"){
	string file = argv[argn+1];
	long ampSize = atol(argv[argn+2]);

	// helps deal with the optional last two args
	float fractId = float(atoi(argv[argn+3])) / float(100);
	if ( argn + 4 >= argc or argv[argn + 4][0] == '-'){
	  assembler->addMatePairFile(file,ampSize,fractId);
	  argn += 4;
	} else if ( argn + 6 >= argc or argv[argn + 6][0] == '-'){
	  int startCycle = atoi(argv[argn+4]);
	  int numCycles = atoi(argv[argn+5]);
	  assembler->addMatePairFileLimitUse(startCycle,numCycles,file,ampSize);
	  argn += 6;
	} else {
	  int startCycle = atoi(argv[argn+4]);
	  int onCycles = atoi(argv[argn+5]);
	  int skipCycles = atoi(argv[argn+6]);
	  assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize);
	  argn += 7;
	}
      } else if (arg == "-spf"){
	string file = argv[argn+1];
	long readSize = atol(argv[argn+2]);
	long ampSize = atol(argv[argn+3]);
	if ( argn + 4 >= argc or argv[argn + 4][0] == '-'){
	  assembler->addFalsePairFile(file,readSize,ampSize);
	  argn += 4;
	} else if ( argn + 6 >= argc or argv[argn + 6][0] == '-'){
	  int startCycle = atoi(argv[argn+4]);
	  int numCycles = atoi(argv[argn+5]);
	  assembler->addFalsePairFileLimitUse(startCycle,numCycles,file,readSize,ampSize);
	  argn += 6;
	} else {
	  int startCycle = atoi(argv[argn+4]);
	  int onCycles = atoi(argv[argn+5]);
	  int skipCycles = atoi(argv[argn+6]);
	  assembler->addFalsePairFileLimitUse(startCycle,onCycles,skipCycles,file,readSize,ampSize);
	  argn += 7;
	}
      } else if (arg == "-spfp"){
	string file = argv[argn+1];
	long readSize = atol(argv[argn+2]);
	long ampSize = atol(argv[argn+3]);

	// helps deal with the optional last two args
	float fractId = float(atoi(argv[argn+4])) / float(100);
	if ( argn + 5 >= argc or argv[argn + 5][0] == '-'){
	  assembler->addFalsePairFile(file,readSize,ampSize,fractId);
	  argn += 5;
	} else if ( argn + 7 >= argc or argv[argn + 7][0] == '-'){
	  int startCycle = atoi(argv[argn+5]);
	  int numCycles = atoi(argv[argn+6]);
	  assembler->addFalsePairFileLimitUse(startCycle,numCycles,file,readSize,ampSize);
	  argn += 7;
	} else {
	  int startCycle = atoi(argv[argn+5]);
	  int onCycles = atoi(argv[argn+6]);
	  int skipCycles = atoi(argv[argn+7]);
	  assembler->addFalsePairFileLimitUse(startCycle,onCycles,skipCycles,file,readSize,ampSize);
	  argn += 8;
	}
      } else if (arg == "-icf"){
	string filename = argv[argn+1];
	float countFactor = atof( argv[argn+4] );
	int numSteps = atoi( argv[argn+2] );
	int cyclesPerStep = atoi( argv[argn+3] );
	assembler->addInitialFile(filename,numSteps,cyclesPerStep,countFactor);
	argn += 5;
      } else if (arg == "-picf"){
	long initialInput = atol( argv[argn+1] );
	string filename = argv[argn+2];
	float countFactor = atof( argv[argn+5] );
	int numSteps = atoi( argv[argn+3] );
	int cyclesPerStep = atoi( argv[argn+4] );
	assembler->addInitialFile(initialInput,filename,numSteps,cyclesPerStep,countFactor);
	argn += 6;
      } else if (arg == "-icfNt"){
	string filename = argv[argn+1];
	float countFactor = atof( argv[argn+4] );
	int numSteps = atoi( argv[argn+2] );
	int cyclesPerStep = atoi( argv[argn+3] );
	assembler->addInitialFileNoTarget(filename,numSteps,cyclesPerStep,countFactor);
	argn += 5;
      } else if (arg == "-picfNt"){
	long initialInput = atol( argv[argn+1] );
	string filename = argv[argn+2];
	float countFactor = atof( argv[argn+5] );
	int numSteps = atoi( argv[argn+3] );
	int cyclesPerStep = atoi( argv[argn+4] );
	assembler->addInitialFileNoTarget(initialInput,filename,numSteps,cyclesPerStep,countFactor);
	argn += 6;
      } else if (arg == "-nc"){
	// nothing happens here because this info was already collected
	argn += 2;
      } else if (arg == "-repmask"){
	int cycleNum = atoi(argv[argn+1]);
	string whenInCycle = argv[argn+2];
	float minStdDev = atof(argv[argn+3]);
	float minFoldUp = atof(argv[argn+4]);
	long minSize = atol(argv[argn+5]);
	float minFractId = atof(argv[argn+6]) / float(100);
	if (minFractId < 0.25){
	  cerr << "-repmask min % ID must be at least 25%" << endl;
	  willDoAssembly = false;
	}
	argn += 7;
	if (argn < argc and argv[argn][0] != '-'){
	  if (whenInCycle[0] == 's'){
	    assembler->findRepsBeforeCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId, argv[argn]);
	  } else if (whenInCycle[0] == 'f'){
	    assembler->findRepsAfterCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId, argv[argn]);
	  } else {
	    cerr << "-repmask second arg must be 's' or 'f', not '" << whenInCycle << "'." << endl;
	    willDoAssembly = false;
	  }
	  ++argn;
	} else {
	  if (whenInCycle[0] == 's'){
	    assembler->findRepsBeforeCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId);
	  } else if (whenInCycle[0] == 'f'){
	    assembler->findRepsAfterCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId);
	  } else {
	    cerr << "-repmask second arg must be 's' or 'f', not '" << whenInCycle << "'." << endl;
	    willDoAssembly = false;
	  }
	}
      } else if (arg == "-target" or arg == "-targetF"){
	bool fullFileMode = (arg == "-targetF");
	float fractId = atof( argv[argn+1] ) / float(100.0);
	int cyclesToSkip = atoi( argv[argn+2] );
	if ( argn + 3 >= argc or argv[argn+3][0] == '-'){
	  assembler->addTargetMode(fractId,cyclesToSkip,fullFileMode);
	  argn += 3;
	} else {
	  int numFilterCycles = atoi( argv[argn+3] );
	  int numSkipCycles = atoi( argv[argn+4] );
	  assembler->addTargetMode(fractId,cyclesToSkip,numFilterCycles,numSkipCycles,fullFileMode);
	  argn += 5;
	}
      } else if (arg == "-r"){
	assembler->setNucMatchScore( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-q"){
	assembler->setNucMismatchPenalty( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-G"){
	assembler->setOpenGapPenalty( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-E"){
	assembler->setExtendGapPenalty( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-badf"){
	string filename = argv[argn+1];
	float fractId = atof( argv[argn+2] ) / float(100.0);
	assembler->addBadSequenceFilter(filename, fractId);
	argn += 3;
      } else if (arg == "-maxHp"){
	assembler->addHomopolymerFilter( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-maxDi"){
	assembler->addDinucRepeatFilter( atol(argv[argn+1]) );
	argn += 2;
      } else if (arg == "-rqf"){
	float minFractGood = atof(argv[argn+1]) / float(100);
	float minProbCorrect = atof(argv[argn+2]);
	argn += 3;
	if (argn < argc and argv[argn][0] != '-'){
	  int numSkipCycles = atoi(argv[argn]);
	  int numRunCycles = atoi(argv[argn+1]);
	  assembler->addReadQualityFilter(minFractGood, minProbCorrect, numSkipCycles, numRunCycles);
	  argn += 2;
	} else {
	  assembler->addReadQualityFilter(minFractGood, minProbCorrect);
	}
      } else if (arg == "-rnf"){
	float minFractGood = atof(argv[argn+1]) / float(100);
	argn += 2;
	if (argn < argc and argv[argn][0] != '-'){
	  int numSkipCycles = atoi(argv[argn]);
	  int numRunCycles = atoi(argv[argn+1]);
	  assembler->addReadCalledBasesFilter(minFractGood, numSkipCycles, numRunCycles);
	  argn += 2;
	} else {
	  assembler->addReadCalledBasesFilter(minFractGood);
	}
      } else if (arg == "-icbf"){
	string filename = argv[argn+1];
	float fractId = atof( argv[argn+2] ) / float(100.0);
	argn += 3;
	if (argn < argc and argv[argn][0] != '-'){
	  long maxSeqLen = atol( argv[argn] );
	  ++argn;
	  assembler->icBadSequenceFilter(filename, fractId, maxSeqLen);
	} else {
	  assembler->icBadSequenceFilter(filename, fractId);
	}
      } else if (arg == "-icmHp"){
	long maxHp = atol(argv[argn+1]);
	argn += 2;
	if (argn < argc and argv[argn][0] != '-'){
	  long maxSeqLen = atol( argv[argn] );
	  ++argn;
	  assembler->icHomopolymerFilter( maxHp, maxSeqLen );
	} else {
	  assembler->icHomopolymerFilter( maxHp );
	}
      } else if (arg == "-icmDi"){
	long maxDi = atol(argv[argn+1]);
	argn += 2;
	if (argn < argc and argv[argn][0] != '-'){
	  long maxSeqLen = atol( argv[argn] );
	  ++argn;
	  assembler->icDinucRepeatFilter( maxDi, maxSeqLen );
	} else {
	  assembler->icDinucRepeatFilter( maxDi );
	}
      } else if (arg == "-icqf"){
	float minFractGood = atof(argv[argn+1]) / float(100);
	float minProbCorrect = atof(argv[argn+2]);
	argn += 3;
	if (argn < argc and argv[argn][0] != '-'){
	  long maxSeqLen = atol(argv[argn]);
	  ++argn;
	  assembler->icReadQualityFilter(minFractGood, minProbCorrect, maxSeqLen);
	} else {
	  assembler->icReadQualityFilter(minFractGood, minProbCorrect);
	}
      } else if (arg == "-icnf"){
	float minFractGood = atof(argv[argn+1]) / float(100);
	argn += 2;
	if (argn < argc and argv[argn][0] != '-'){
	  long maxSeqLen = atol(argv[argn]);
	  ++argn;
	  assembler->icReadCalledBasesFilter(minFractGood, maxSeqLen);
	} else {
	  assembler->icReadCalledBasesFilter(minFractGood);
	}
      } else if (arg == "-lenf"){
	long minLength = atol( argv[argn+1] );
	int cyclesToSkip = atoi( argv[argn+2] );
	assembler->addLengthFilter(minLength, cyclesToSkip);
	argn += 3;
      } else if (arg == "-link"){
	long maxNumLinks = atol( argv[argn+1] );
	assembler->setLinkMax(maxNumLinks);
	argn += 2;
      } else if (arg == "-o"){
	string filename = argv[argn+1];
	assembler->addOutputFile(filename);
	argn += 2;
      } else if (arg == "-nco"){
	assembler->setOutputInterval(atoi( argv[argn+1] ));
	argn += 2;
      } else if (arg == "-a"){
	omp_set_num_threads( atoi( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-dbk"){
	assembler->setDeBruijnKmerSize( atoi( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-dbmax"){
	assembler->setDeBruijnMaxSeqLength( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-dbms"){
	assembler->setDeBruijnMinSeqNum( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-mtpf"){
	maxThreadsPerFile = atoi( argv[argn+1] );
	argn += 2;
      } else if (arg == "-trim"){
	int cycleNum = atoi( argv[argn+1] );
	float minCoverage = atof( argv[argn+2] );
	if ( argn + 3 < argc and argv[argn+3][0] != '-' ) {
	  long minLength = atol( argv[argn+3] );
	  assembler->addTrim(cycleNum, minCoverage, minLength);
	  argn += 4;
	} else {
	  assembler->addTrim(cycleNum, minCoverage);
	  argn += 3;
	}
      } else if (arg == "-trimB"){
	int skipCycles = atoi( argv[argn+1] );
	float minCoverage = atof( argv[argn+2] );
	if ( argn + 3 < argc and argv[argn+3][0] != '-' ) {
	  long minLength = atol( argv[argn+3] );
	  assembler->addBasalTrim(skipCycles, minCoverage, minLength);
	  argn += 4;
	} else {
	  assembler->addBasalTrim(skipCycles, minCoverage);
	  argn += 3;
	}
      } else if (arg == "-trimI"){
	float minCoverage = atof( argv[argn+1] );
	if ( argn + 3 < argc and argv[argn+2][0] != '-' ) {
	  long minLength = atol( argv[argn+2] );
	  assembler->addInitialTrim(minCoverage, minLength);
	  argn += 3;
	} else {
	  assembler->addInitialTrim(minCoverage);
	  argn += 2;
	}
      } else if (arg == "-reset"){
	int numResets = 0;
	argn += 1;
	while (argn < argc and argv[argn][0] != '-'){
	  assembler->addResetCycle( atoi( argv[argn] ) );
	  argn += 1;
	  numResets++;
	}
	if (numResets == 0 ){
	  cerr << "-reset should have at least one cycle number after it." << endl;
	  willDoAssembly = false;
	}
      } else if (arg == "-mol"){ // min overlap for overhang assemblies
	assembler->setMiniMinOverlap( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-tol"){ // threshold for scaling overlap for overhang assemblies
	assembler->setMiniMinOvlContigThreshold( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-mpi"){ // min percent ID for overhang assemblies
	int pctId = atoi(argv[argn+1]);
	if (pctId < 25 or pctId > 100){
	  cerr << "% ID required for mini-assembly is out of range 25-100" << endl;
	  willDoAssembly = false;
	}
	float fractId = float(pctId) / 100.0;
	assembler->setMiniMinFractId(fractId);
	argn += 2;
      } else if (arg == "-tpi"){ // threshold for scaling percent ID for overhang assemblies
	assembler->setMiniFractIdContigThreshold( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-MPI"){ // min percent ID for meta assemblies
	int pctId = atoi(argv[argn+1]);
	if (pctId < 25 or pctId > 100){
	  cerr << "% ID required for mini-assembly is out of range 25-100" << endl;
	  willDoAssembly = false;
	}
	float fractId = float(pctId) / 100.0;
	assembler->setMetaMinFractId(fractId);
	argn += 2;
      } else if (arg == "-TPI"){ // threshold for scaling percent ID for meta assemblies
	assembler->setMetaFractIdContigThreshold( atol( argv[argn+1] ) );
	argn += 2;
      } else if (arg == "-log"){
	string type = string(argv[argn+1]);
	if (type == "n"){ listenerIsNull = true; }
	else if (type == "c"){ }
	else if (type == "v"){ listenerIsVerbose = true; }
	else {
	  cerr << "'" << type << "' is not a valid stdout type ('n', 'c', or 'v')." << endl;
	  willDoAssembly = false;
	}
	argn += 2;
      } else if (arg == "-logf"){
	includeListenerFile = true;
	listenerFilename = argv[argn+1];
	argn += 2;
      } else {
	cerr << "'" << arg << "' is an invalid argument." << endl;
	willDoAssembly = false;
	argn++;
      }
      if (argn != argc and argv[argn][0] != '-'){
	ostringstream outMessage;
	cerr << arg << " was given an incorrect number of args." << endl;
	willDoAssembly = false;
      }
    }

    // check again!
    if ( willDoAssembly ){

      // TEMPORARY set-up for returning the command that launched the job
      ostringstream launchCommand;
      launchCommand << "Launch command:";
      for (int n = 0; n < argc; ++n){ launchCommand << " " << argv[n]; }
      cout << launchCommand.str().c_str() << endl;


      if (listenerIsNull){ assembler->makeLogNull(); }
      else if (listenerIsVerbose){ assembler->makeLogVerbose(); }
      if (includeListenerFile == true){ assembler->verboseLogFile(listenerFilename); }

      if (maxThreadsPerFile != 1){ assembler->setMaxThreadsPerFile(maxThreadsPerFile); }

      assembler->runAssembly();

    }
    delete assembler;
  }

  if (willDoAssembly){ return 0; }
  else { return 1; }
}



