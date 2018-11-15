# calSOV

## Description
calculate the Segment Overlap Score (SOV) score for the predicted secondary
structures

## Author
Nanjiang Shu

Email: nanjiang.shu@scilifelab.se

## Reference

Zemla A, Venclovas C, Fidelis K, Rost B., "A modified definition of Sov, a
segment-based measure for protein secondary structure prediction assessment.",
Proteins, 1999, 34(2), pp. 220-3

Note: when the observed secondary structure contains unknown structures, the
calculation might be different compared to e.g. SOV\_refine.

## Installation
This installation has been tested on Ubuntu and CentOS
Requirement: g++

		make

		make install

If need to be installed to a customized location

		make "BINPATH=/your/path" install

## Usage

```
Usage: calSOV [options] file1 file2 ...
OPTIONS:
  -o, --out <outfile>  : output the result to outfile, (default = stdout)
  -l, --list listfile  : set the file containing a list of predicted files
  -f, --format 0|1|2|3 : set the format of the input file, (default = 1)
  -q3, --isOutQ3 yes|no: wheter output q3, (default = yes)
  -m, --method  0|1|2|3: set the method to get predicted secondary structure state from Res file, (default = 1)
  -statq3              : anylyze the Q3 for residues predicted at different confidences 
  -chkfstdir           : folder to which Res from checkfirst program are stored, this need to be supplied when -statq3 is enabled 
  -selconf     0|1     : selection of the confidence, (default = 0) 
                       : 0 -- raw confidence, 1 -- normalized confidence 
  -binwidth <real>     : set halfbinwidth, (default = 1.0) 
  -max, --max-length   : set the max length of the sequence, (default = 200000)
  -proof, --proof      : enable proof reading
  -proofmethod int     : set the method for proof reading, (default = 0)
  -polypara a1 a2 a3   : input three parameters for polynormail function which normalize the confidence
                       : default a1=-0.007100 , a2 = 1.808400 , a3 = -11.400000 
  -h, --help           : print this help message and exit

Created on 2009-07-23, updated 2016-11-03, Nanjiang Shu

Format description of input file:
    Format 0: prediction Res* file
    Format 1: AA OSEC PSEC NUM
    Format 2: Fasta format, >AA, >OSEC, >PSEC, order not important

Examples:
    ./calSOV -f 2 test/test1_f2.txt  # the same example as in the reference
    ./calSOV -f 2 test/test2_f2.txt
    ./calSOV -f 1 test/16VPA.psipred.format_1.txt
```
