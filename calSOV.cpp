/*
 * =====================================================================================
 *
 *       Filename:  calSOV.cpp
 *
 *    Description:  calculate the SOV score and Q3 score given a prediction 
 *
 *        Version:  1.0
 *        Created:  07/23/2009 03:49:05 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nanjiang Shu (Shu), nanjiang.shu@dbb.su.se
 *        Company:  Department of Biochemistry and Biophysics, Stockholm Univesity
 *
 * =====================================================================================
 */

/* ChangeLog 2009-10-26
 *  The confidence normalization function, 
 *  a1, a2 and a3 can be modified from argument
 * */


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <set>
#include <list>
#include <vector>
#include <utility>
#include "array.h"
#include "mytemplate.h"
/*#include "myfunc.h"*/
using namespace std;

#if defined(_Windows) || defined(__WINDOWS__) || \
    defined(__WIN32__) || defined(WIN32) || \
defined(__WINNT__) || defined(__NT__)
#   ifndef WINDOWS
#       define WINDOWS
#   endif
#endif

#define NUM_HSR_STATE 3


int method_HSR  = 1 ;
bool isStatQ3Conf = false; /*whether analyze the Q3 for residues at each confidence bin*/
int typeConfidence = 0; /*default using the raw confidence*/
char HSRalphabet[] = "HSR-";
double binWidth = 1.0; /*halfBinWidth*/
int maxSeqLength = MAX_SEQ_LENGTH +1 ;
bool isProofReading = false; /*2009-08-03*/
int method_ProofReading = 0;

/*parameters for the function 
 * normConf = a1 * rawConf*rawConf + a2 * rawConf + a3 
 * 2009-10-26*/
double poly_a1 = -0.0071;
double poly_a2 = 1.8084;
double poly_a3 = -11.4;

void PrintHelp(char** argv)
{
    fprintf(stdout,"Usage: calSOV [options] file1 file2 ...\n");
    fprintf(stdout,"OPTIONS:\n");
    fprintf(stdout,"  -o, --out <outfile>  : output the result to outfile, (default = stdout)\n");
    fprintf(stdout,"  -l, --list listfile  : set the file containing a list of predicted files\n");
    fprintf(stdout,"  -f, --format 0|1|2|3 : set the format of the input file, (default = 1)\n");
    fprintf(stdout,"  -q3, --isOutQ3 yes|no: wheter output q3, (default = yes)\n");
    fprintf(stdout,"  -m, --method  0|1|2|3: set the method to get predicted secondary structure state from Res file, (default = 1)\n");
    fprintf(stdout,"  -statq3              : anylyze the Q3 for residues predicted at different confidences \n");
    fprintf(stdout,"  -chkfstdir           : folder to which Res from checkfirst program are stored, this need to be supplied when -statq3 is enabled \n");
    fprintf(stdout,"  -selconf     0|1     : selection of the confidence, (default = 0) \n");
    fprintf(stdout,"                       : 0 -- raw confidence, 1 -- normalized confidence \n");
    fprintf(stdout,"  -binwidth <real>     : set halfbinwidth, (default = 1.0) \n");
    fprintf(stdout,"  -max, --max-length   : set the max length of the sequence, (default = %d)\n", maxSeqLength-1);
    fprintf(stdout,"  -proof, --proof      : enable proof reading\n");
    fprintf(stdout,"  -proofmethod int     : set the method for proof reading, (default = 0)\n");
    fprintf(stdout,"  -polypara a1 a2 a3   : input three parameters for polynormail function which normalize the confidence\n");
    fprintf(stdout,"                       : default a1=%8.6lf , a2 = %8.6lf , a3 = %8.6lf \n", poly_a1, poly_a2 ,poly_a3);
    fprintf(stdout,"  -h, --help           : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2009-07-23, updated 2016-11-03, Nanjiang Shu\n");
    fprintf(stdout,"\nFormat description of input file:\n");
    fprintf(stdout,"    Format 0: prediction Res* file\n");
    fprintf(stdout,"    Format 1: AA OSEC PSEC NUM\n");
    fprintf(stdout,"    Format 2: Fasta format, >AA, >OSEC, >PSEC, order not important\n");
    fprintf(stdout,"\nExamples:\n");
    fprintf(stdout,"    %s -f 1 test/16VPA.psipred.format_1.txt\n", argv[0]);
    fprintf(stdout,"    %s -f 2 test/16VPA.psipred.format_2.txt\n", argv[0]);
}
void PrintVerboseHelp() { }

int coverage(int a1, int b1, int a2, int b2)
    /*
    return the coverage of two intervals
    a1, b1, a2, b2 are integers
    when the return value <=0, it means there is no coverage
    */
{
    return (min(b1,b2)-max(a1,a2));
}
int my_strcpy(char* to, const char* from, int max)/*{{{*/
/******************************************************************************
 * my modification of strcpy
 * copy max characters from "from" to "to", add NULL terminator automatically
 * updated 2008-04-23, memcpy win in speed when the copying string is long
 * updated 2011-10-27:
 *   since calling strlen will be time consuming for very large from string,
 *   e.g. when "from" is the buffer of a whole trunk of file. Therefore, strlen
 *   part is removed.
 *****************************************************************************/
{
    if(max < 200) {
        strncpy(to,from,max);
    } else {
        memcpy(to,from,max);
    }
    to[max] = '\0';
    return max;
}/*}}}*/
char *rootname(const char* filename, char* rtname, int max_rtname /*= MAX_PATH*/)/*{{{*/
/*****************************************************************************
 * rootname
 * given the file name, 
 * return the rootname of the filename
 ****************************************************************************/
{
    const char *pch;
    char *pstr;
    if((pch = strrchr(filename,'/')) != NULL)
        pstr = (char*) pch+1;
    else
        pstr = (char*) filename;

    if((pch = strrchr(pstr,'.')) != NULL)
        my_strcpy(rtname,pstr, min((int)(pch - pstr), max_rtname));
    else
        rtname = pstr;
    return rtname;
}
/*}}}*/
int option_parser_filename(int argc, char **argv, int beg, char *filename)/*{{{*/
/*****************************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--out"
 ****************************************************************************/
{
    int i ; 
    bool isNonOptionArg = false;
    bool isFileNameSet = false;

    for(i = beg +1 ; i < argc ; i++)
    {
        if(argv[i][0] == '-' && strcmp(argv[i], "--") != 0 && !isNonOptionArg)
        {
            fprintf(stderr,"option '%s' must be followed by a filename, not option\n", argv[beg]);
            return -1;
        }
        else if(strcmp(argv[i], "--") == 0 && !isNonOptionArg)
        {
            isNonOptionArg = true;
        }
        else
        {
            my_strcpy(filename, argv[i], MAX_PATH-1);
            isFileNameSet = true;
            break;
        }
    }

    if(!isFileNameSet)
    {
        fprintf(stderr,"option '%s' must be followed by a filename\n", argv[beg]);
        return -1;
    }
    else 
    {
        return i+1;
    }
}
/*}}}*/
template <class T> int option_parser_numeric(int argc, char **argv, int beg, T &x, bool isRangeSet /*= false*/, T min /*= MIN_INT*/, T max /*= MAX_INT*/)/*{{{*/
/*****************************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--value"
 ****************************************************************************/
{
    int i ; 
    bool isValueSet = false;
    double tmp;
    i = beg +1;

    if (i < argc)
    {
        if(IsNumeric(argv[i]))
        {
            tmp = atof(argv[i]);
            if(isRangeSet)
            {
                if(tmp < min || tmp > max)
                {
                    fprintf(stderr,"Invalid value! Value after option '%s' must be in the range of [%g %g]\n", argv[beg], double(min), double(max));
                    return -1;
                }
            }
            x = T(tmp);
            isValueSet = true;
        }
    }

    if(!isValueSet)
    {
        fprintf(stderr,"option '%s' must be followed by a numerical value\n", argv[beg]);
        return -1;
    }
    else 
    {
        return i+1;
    }
}
template int option_parser_numeric<int>   (int argc, char **argv, int beg, int &x   , bool isRangeSet /*= false*/, int min /*= MIN_INT*/      , int max/* = MAX_INT*/);
template int option_parser_numeric<float> (int argc, char **argv, int beg, float &x , bool isRangeSet /*= false*/, float min /*= MIN_FLOAT*/  , float max /*= MAX_FLOAT*/);
template int option_parser_numeric<double>(int argc, char **argv, int beg, double &x, bool isRangeSet /*= false*/, double min/* = MIN_DOUBLE*/, double max /*= MAX_DOUBLE*/);
template int option_parser_numeric<int8>(int argc, char **argv, int beg, int8 &x, bool isRangeSet /*= false*/, int8 min/* = MIN_DOUBLE*/, int8 max /*= MAX_DOUBLE*/);
/*}}}*/
int ReadNextSeq_FASTA(FILE *fp, char* seq, int *pSeq_type /*= NULL*/, int maxlength /*= LONGEST_SEQ*/, char* annotationLine/*=NULL*/, int maxSizeAnnotationLine /*=50*/)/*{{{*/
    /****************************************************************************
     * ReadNextSeq_FASTA()
     * read in the fasta format sequence from the file stream 
     * return the length of the sequence if successful
     * return 0 or minus value if no more sequence can be read from the file stream
     * The leading white spaces are ignored
     * check the type (DNA or AA) of fasta sequence file based the annotation line 
     * Last modified, 2007-02-12, Nanjiang Shu
     * Updated 2010-04-22: the annotation line (without the leading ">") can be
     * read in, note that the maxSizeAnnotationLine must be supplied
     ***************************************************************************/
{
    int   c;
    int   i;

    do{  /* ignore the leading white spaces */ 
        c = getc(fp);
    }while(isspace(c));

    if(pSeq_type != NULL) 
        *pSeq_type = UNKNOWN_SEQ_TYPE; /* initializing sequence type*/

    if(c  == '>') 
    { 
        Array1D <char> line_1darray(maxSizeAnnotationLine +1);
        char *line = line_1darray.array1D;
        fgetline(fp,line,maxSizeAnnotationLine);
        if( pSeq_type != NULL)
        {
            if(strstr(line, "protein") != NULL)
                *pSeq_type = AA_SEQ;
            else if( strstr(line,"nucleic") != NULL)
                *pSeq_type = DNA_SEQ;
            else if( strstr(line,"shape") != NULL)
                *pSeq_type = SHAPE_SEQ;
            else
                *pSeq_type = UNKNOWN_SEQ_TYPE;
        }
        if (annotationLine != NULL)
        {
            my_strcpy(annotationLine,line, maxSizeAnnotationLine-1);  /*read in the annotation line*/ 
        }
    }
    else  
    {
        fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream if there is no annotation line*/ 
    }

    i = 0 ;
    while((c = getc(fp)) != EOF) 
    {
        if(c == '>') 
        {
            fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream*/
            break;  /* read in the first sequence if there are multiple sequences in the file*/ 
        }
        if(! isspace(c)) /* neglect white spaces and return characters in sequence region*/ 
        {
            seq[i] = c ; i ++ ;
            if(i >= maxlength)
            {
                fprintf(stderr,"Error, sequence longer then maxlength = %d\n", maxlength);
                return 1;

            }
        }
    }
    seq[i] = '\0' ;
    if(c == EOF && i == 0)
        return EOF;
    else
        return i; /* return the length of sequence*/ 
}/*}}}*/

bool IsNumeric(const char* str)/*{{{*/
//**********************************************************************
//IsNumeric(char*)
//check if a string is a numeric number
{
	int i = 0;
    int cntdot = 0;
	while(str[i] != '\0')
	{
        if( i == 0)
        {
            if(!isdigit(str[i]) && str[i] != '+' && str[i] != '-' && str[i] != '.' )
                return false;
            if(str[i] == '.')
                cntdot ++;
        }
        else
        {
            if(! isdigit(str[i]) && str[i] != '.')
                return false;
            if(str[i] == '.')
                cntdot ++;
        }
        i ++;
	}
    if(cntdot >= 2)
        return false;
    else 
        return true;
}
/*}}}*/
bool IsInCharSet(const char ch, const char *charSet, int n /*= 0 */)/*{{{*/
/*****************************************************************************
 * check if the character "ch" is in charSet,
 ****************************************************************************/
{
	if(n == 0)
        n = strlen(charSet);
    int i;
	for(i = 0 ;i < n ; i ++)
	{
		if(ch == charSet[i])
			return true;
	}
	return false;
}/*}}}*/
int fgetline(FILE* fp, char* line, int max/* = 0x7FFFFFFF*/)/*{{{*/
/*****************************************************************************
 * Read one line from fp, copying it to line array (but no more than max
 * chars). Does not place terminating \n in line array.  
 * it can be called without "max" flag, but should make sure the allocated
 * memory for "line" is larger than the longest line.
 *
 * Returns: line length, or 0 for empty line, or "EOF" for end-of-file.
 * 
 * since fgetc(), getc(), getchar(), returns int value,always use an int
 * variable to store the result of the fgetc(), getc() and getchar().  
 * getc() is faster than fgetc()
 *   LOG: 2006-04-26 16:31:29 Wednesday  Week 17 <nanjiang@shu>
 *   bug fixed for '\n' return keys, 
 *   now it is valid for both dos and unix
 ****************************************************************************/
{
    int nch = 0; /* record number of characters actually read */
    int c;
    max = max - 1;			/* leave room for '\0' */
    while((c = getc(fp)) != EOF)
    {
        if(c == 0x0d ) continue; /* in unix, '\n'= 0x0a, thus 0x0d will be cheated as another character*/ 
        if(c == '\n') break; /* in dos, '\n' is also equal to 0x0a, but the preceding 0x0d will not be read*/ 

        if(nch < max)
        {
            line[nch] = c;
            nch = nch + 1;
        }
    }
    line[nch] = '\0';

    if(c == EOF && nch == 0) return EOF;
    else return nch;
}/*}}}*/
int checkfilestream(FILE *fp, const char* filename, const char *mode, bool isAssert /*= false*/)/*{{{*/
{
    if( fp == NULL)
    {
        fprintf(stderr,"Can not open file '%s' with mode '%s'\n", filename,mode);
        if(isAssert)
        {
            assert(fp != NULL);
        }
        return -1;
    }
    else
        return 0;
}
/*}}}*/
double GetRawHSRConfidence(int probH, int probS, int probR)/*{{{*/
/*Determine raw confidence of the prediction based on the probability on H, S
 * and R
 * 2009-07-21, Nanjiang*/
{
    double sumHSR =  double(probH+probS+probR) + 1e-6; /*plus 1e-6 to avoid division by zero*/
    double maxHSR =  double (max (max(probH, probS), probR ));
    double rawConf = maxHSR/sumHSR;
    return rawConf;
}/*}}}*/
float GetNormHSRConfidence(int probH, int probS, int probR)/*{{{*/
/*Determine normalized confidence of the prediction based on the probability on H, S
 * and R, parameters were obtained by polynormial regression (order 2) of the
 * raw confidence and two Q3 at different raw confidence bins
 * 2009-08-02 
 * */
{
    double sumHSR =  double(probH+probS+probR) + 1e-6; /*plus 0.00001 to avoid division by zero*/
    double maxHSR =  double (max (max(probH, probS), probR ));
    double rawConf =(maxHSR/sumHSR)*100.0;

    double normConf = poly_a1 * rawConf*rawConf + poly_a2 * rawConf + poly_a3;

    /*method 2, three order regression, 2009-08-02 */
   //y = -0.0001x3 + 0.0225x2 - 0.1356x + 28.889 
    //double a1 =-0.0001;
    //double a2 = 0.0225;
    //double a3 = -0.1356;
    //double a4 =28.889; 
    //double normConf =  a1 * rawConf*rawConf*rawConf + a2 * rawConf*rawConf + a3*rawConf + a4;

    normConf /= 100.0;
    if (normConf < 0.0)
    {
        normConf = 0.0;
    }
    else if (normConf> 1.0)
    {
        normConf = 1.0;
    }
    return normConf;

}/*}}}*/
int GetHSRState(int probH, int probS, int probR, int method_HSR /*= 1*/)/*{{{*/
/*Determine the state of secondary structure based on the probability on H, S
 * and R
 * the return value is 0 or 1 or 2, representing H, S and R respectively
 * 2009-07-13, Nanjiang*/
{
    int hsrState = 2;
    if (method_HSR == 0)
    {
        if  (    (probS>=probH) && (probS>=probR) ) { hsrState = 1; }
        else if ((probH>=probS) && (probH>=probR) ) { hsrState = 0; }
        else                                        { hsrState = 2; }
    }
    else if (method_HSR == 1)
    {
        if  (    (probR>=probH) && (probR>=probS) ) { hsrState = 2; }
        else if ((probS>=probH) && (probS>=probR) ) { hsrState = 1; }
        else { hsrState = 0; }
    }
    else 
    {
        if  (    (probR>=probH) && (probR>=probS) ) { hsrState = 2; }
        else if ((probH>=probS) && (probH>=probR) ) { hsrState = 0; }
        else { hsrState = 1; }
    }

    return hsrState;
}/*}}}*/
int ReadInSecPredFile(const char *infile, char *aaSeq, char *obsSec, char *predSec, double *predConf, int typeConfidence, int fileFormat)/*{{{*/
{
    int linesize;
    int maxline = 400;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char aa;
    char obs_sec;
    char pred_sec;
    int status_sscanf = 0;

    int seqLength = 0;
    FILE *fp = fopen(infile, "r");
    checkfilestream(fp, infile, "r", true);

    if (fileFormat == 0) /*Res file format*//*{{{*/
    {
        int num ; 
        char shape;
        int probH;
        int probS;
        int probR;
        int predHSRstate = 0; /*predicted HSR state, 0 -- H, 1 -- S, 2 -- R*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            status_sscanf = sscanf(line,"%d %c %c %c %d %d %d" , &num,&aa,&shape, &obs_sec,&probH,&probS,&probR);
            if (status_sscanf != 7)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 7);
            }
            predHSRstate = GetHSRState(probH, probS, probR, method_HSR);

            aaSeq[cntRes] = aa;
            obsSec[cntRes] = obs_sec;
            predSec[cntRes] = HSRalphabet[predHSRstate];
            if (typeConfidence == 0) /*Raw confidence*/
            {
                predConf[cntRes] = GetRawHSRConfidence(probH, probS, probR);
            }
            else 
            {
                predConf[cntRes] = GetNormHSRConfidence(probH, probS, probR);
            }

            cntRes ++;
        } 
        seqLength = cntRes;
    }/*}}}*/
    else if (fileFormat == 1) /*AA OSEC PSEC NUM*//*{{{*/
    {
        /*AA OSEC PSEC NUM*/
        linesize=fgetline(fp,line,maxline); /*neglect the header line*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            status_sscanf = sscanf(line," %c %c %c" , &aa,&obs_sec,&pred_sec);
            if (status_sscanf != 3)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 3);
            }
            if (obs_sec == 'E' ) { obs_sec = 'S'; }
            else if (obs_sec == 'C' || obs_sec == 'L') { obs_sec = 'R'; }

            if (pred_sec == 'E' ) { pred_sec = 'S'; }
            else if (pred_sec == 'C' || pred_sec == 'L') { pred_sec = 'R'; }

            aaSeq[cntRes] = aa;
            obsSec[cntRes] = obs_sec;
            predSec[cntRes] = pred_sec;
#ifdef DEBUG
            fprintf(stdout, "line=%s\n", line);
            fprintf(stdout, "cntRes=%d\n", cntRes);
#endif
            cntRes ++;
            if (cntRes >= maxSeqLength){
                fprintf(stdout,"Number of residues in the input file exceeds the maxSeqLength (%d). Please check your input file. Exit.\n", maxSeqLength);
                assert(cntRes<maxSeqLength);
            }
        } 
        seqLength = cntRes;
    }/*}}}*/
    else if (fileFormat == 2) /*Fasta format >AA, >OSEC, >PSEC*//*{{{*/
    {
        int maxSizeAnno = 100;
        Array1D <char> anno_1darray(maxSizeAnno+1);
        Array1D <char> tmpSeq_1darray(maxSeqLength+1);
        char *tmpSeq = tmpSeq_1darray.array1D;
        char *anno = anno_1darray.array1D;
        int seqlength = 0;
        /*>AA >OSEC >PSEC*/
        int cntSeq = 0;
        while((seqlength = ReadNextSeq_FASTA(fp, tmpSeq, NULL, maxSeqLength, anno, maxSizeAnno))!= EOF)
        {
            if(strstr(anno, "AA") != NULL){
                my_strcpy(aaSeq,tmpSeq, maxSeqLength-1);
            }else if(strstr(anno, "OSEC") != NULL){
                my_strcpy(obsSec,tmpSeq, maxSeqLength-1);
            }else if(strstr(anno, "PSEC") != NULL){
                my_strcpy(predSec,tmpSeq, maxSeqLength-1);
            }
            cntSeq ++;
            if (cntSeq>3){
                fprintf(stderr, "You input file contains more than 3 sequences. Please check the format of your input file\n");
                assert(cntSeq<=3);
            }
        }
        int i = 0;
        seqlength = strlen(obsSec);
        for (i = 0; i<seqlength; i++){
            if (obsSec[i] == 'E' ) { obsSec[i] = 'S'; }
            else if (obsSec[i] == 'C' || obsSec[i] == 'L') { obsSec[i] = 'R'; }
        }

        seqlength = strlen(predSec);
        for (i = 0; i<seqlength; i++){
            if (predSec[i]  == 'E' ) { predSec[i] = 'S'; }
            else if (predSec[i]  == 'C' || predSec[i]  == 'L') { predSec[i]  = 'R'; }
        }
        seqLength = strlen(aaSeq);
    }/*}}}*/
    if (fp != NULL) { fclose(fp); }
    return seqLength;
}/*}}}*/
int ProofReading(char *predSec, int seqLength, int method_ProofReading = 0)/*{{{*/
{
    //proof reading for predicted result, for 1 long strand (sheet)
    int Part1;
    int Part2;
    int ik;
    int ij;
    if (method_ProofReading == 0)/*{{{*/
    {
        int NShe = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'S'  )
            {
                NShe++;
            }
            else
            {
                if (  NShe == 1  )
                {
                    //previous
                    Part1 = ik-3;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='S')  )
                    {
                        predSec[Part1+1] = 'S';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='S')  )
                    {
                        predSec[Part2-1] = 'S';
                    }
                    else
                    {
                        for (ij=ik-NShe; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NShe = 0;
            }
        }


        //proof reading for predicted result, for 1 or 2 long helix
        int NHel = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'H'  )
            {
                NHel++;
            }
            else
            {
                if  (   (NHel>=1) &&  (NHel<=2)   )
                {
                    //previous
                    Part1 = ik-NHel-2;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='H')  )
                    {
                        predSec[Part1+1] = 'H';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='H')  )
                    {
                        predSec[Part2-1] = 'H';
                    }
                    else
                    {
                        for (ij=ik-NHel; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NHel = 0;
            }
        } 
    }/*}}}*/
    else if (method_ProofReading == 1)/*{{{*/
    {  /*do not proof reading single residue sheet, when using the DSSP8to3 scheme BE --> Sheet, 2009-10-29*/
        //proof reading for predicted result, for 1 or 2 long helix
        int NHel = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'H'  )
            {
                NHel++;
            }
            else
            {
                if  (   (NHel>=1) &&  (NHel<=2)   )
                {
                    //previous
                    Part1 = ik-NHel-2;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='H')  )
                    {
                        predSec[Part1+1] = 'H';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='H')  )
                    {
                        predSec[Part2-1] = 'H';
                    }
                    else
                    {
                        for (ij=ik-NHel; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NHel = 0;
            }
        } 
    }/*}}}*/
    return seqLength;
}/*}}}*/
int ReadConf(const char *infile, double *predConf, int typeConfidence)/*{{{*/
{
    int linesize;
    int maxline = 400;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char aa;
    char obs_sec;
    int status_sscanf = 0;

    int seqLength = 0;
    FILE *fp = fopen(infile, "r");
    checkfilestream(fp, infile, "r", true);
    int num ; 
    char shape;
    int probH;
    int probS;
    int probR;
    int cntRes = 0;
    while((linesize=fgetline(fp, line, maxline)) != EOF)
    { 
        status_sscanf = sscanf(line,"%d %c %c %c %d %d %d" , &num,&aa,&shape, &obs_sec,&probH,&probS,&probR);
        if (status_sscanf != 7)
        {
            fprintf(stderr,"sscanf error!");
            fprintf(stderr,"File: %s\n", infile);
            fprintf(stderr,"line: %s\n", line);
            assert (status_sscanf == 7);
        }

        if (typeConfidence == 0) /*Raw confidence*/
        {
            predConf[cntRes] = GetRawHSRConfidence(probH, probS, probR);
        }
        else 
        {
            predConf[cntRes] = GetNormHSRConfidence(probH, probS, probR);
        }

        cntRes ++;
    } 
    seqLength = cntRes;
    if (fp != NULL) { fclose(fp); }
    return seqLength;
}/*}}}*/
list < pair <int, int > > GetSecSeg(char *seq, char state, int seqLength)/*{{{*/
    /*Get segment list for a certain secondary structure state*/
{
    list < pair <int,int> > segList;
    pair <int, int> seg;
    int i,j;
    i = 0;
    while (i < seqLength) 
    {
        int b, e;
        if (seq[i] == state)
        {
            b = i;
            j = i;
            while (j < seqLength && seq[j] == state)
            {
                j ++;
            }
            e = j;
            i = j;
            seg = make_pair(b,e);
            segList.push_back(seg);
        }
        else
        {
            i ++;
        }
    }
    return segList;
}
/*}}}*/
int CalSOV(char *obsSec, char *predSec, int seqLength, double *numCorrectHSR, int *numHSR, double *sovHSR, int numState = 3)/*{{{*/
{
    list <list <pair <int, int> > > segment_obs;
    list <list <pair <int, int> > > segment_pred;
    list <list <pair <int, int> > >::iterator it_sl_1;
    list <list <pair <int, int> > >::iterator it_sl_2;
    list <pair <int, int> > ::iterator it_seg_1;
    list <pair <int, int> > ::iterator it_seg_2;
    int i = 0;

    for (i = 0; i < 3; i ++) /*iterator for secondary structure states*/
    {
        list <pair <int, int> > segList;
        segList = GetSecSeg(obsSec, HSRalphabet[i], seqLength);
        segment_obs.push_back(segList);
        segList = GetSecSeg(predSec, HSRalphabet[i], seqLength);
        segment_pred.push_back(segList);
    }
    vector <int> N(3,0); /*Normalization value*/
    vector <double> sum_un_norm(3,0.0); /*sum of un-normalized score*/
    int total_N = 0;
    double total_sum_un_norm = 0.0;

    int b1, e1, b2, e2, minov, maxov, delta, len_seg1, len_seg2;
    double score;
    for (i = 0; i < 3; i ++) /*iterator for secondary structure states*/
    {
        /*For all segments of the observed secondary structure*/
        it_sl_1 = segment_obs.begin();
        advance( it_sl_1, i );

        it_sl_2 = segment_pred.begin();
        advance( it_sl_2, i );

        for (it_seg_1 = (*it_sl_1).begin(); it_seg_1 != (*it_sl_1).end(); ++it_seg_1)
        {
            b1 = (*it_seg_1).first;
            e1 = (*it_seg_1).second;
            len_seg1 = e1-b1;
            bool isOverlap = false;
            for (it_seg_2 = (*it_sl_2).begin(); it_seg_2 != (*it_sl_2).end(); ++it_seg_2)
            {
                b2 = (*it_seg_2).first;
                e2 = (*it_seg_2).second;
                len_seg2 = e2-b2;
                minov = coverage(b1, e1, b2, e2);
                if (minov > 0) {
                    maxov = max(e1,e2) - min(b1,b2);
                } else {
                    maxov = 0;
                }
                if (minov > 0)
                {
                    isOverlap = true;
                    double min1 = min(maxov - minov, minov);
                    double min2 = min(int(len_seg1/2), int(len_seg2/2));
                    delta = min(min1, min2);
                    score = (minov + delta)/double(maxov)*len_seg1;
                    sum_un_norm[i] += score;
                    N[i] += len_seg1;
                }
            }
            if (isOverlap == false) {
                N[i] += len_seg1;
            }
        }
        sovHSR[i] = sum_un_norm[i]/N[i];
        total_sum_un_norm += sum_un_norm[i];
        total_N += N[i];
        numHSR[i] = N[i];
    }
    sovHSR[numState] = total_sum_un_norm/total_N;
    numHSR[numState] = total_N;

   return numHSR[numState];

}/*}}}*/

int CalSOV_obsolete(char *obsSec, char *predSec, int seqLength, double *numCorrectHSR, int *numHSR, double *sovHSR, int numState = 3)/*{{{*/
    /*=======================================================================
     * numHSR[4] Helix, Sheet, Coil, All
     * sovHSR[4]
     * numState = 3 by default
     *====================================================================*/
{
    int Min_overlap[500] ;
    int Max_extend[500]  ;
    int LenSeg_pred[500];
    int All_HSRLen[3];
    double All_Sum[3];

    int iterHSRState = 0; /*iterator for the HSR state, H, S and R*/
    char thisHSRstate = ' '; /*current HSR state*/
    int i = 0;
    int j = 0;

    for (i=0;i<3;i++)
    {
        All_HSRLen[i] = 0;
        All_Sum[i] = 0.0; 
    }

    int Length = seqLength +1 ;  /*Here Length is 1 residue longer than seqLength*/
    for (iterHSRState=0; iterHSRState<3; iterHSRState++)//corrsponding to H, S, R
    {
        //for helix
        int HSRLen = 0;/*iterator for secondary structure elements*/
        thisHSRstate = HSRalphabet[iterHSRState];
        for(i=0; i<Length; i++)
        { 
            if (  obsSec[i] == thisHSRstate  )
            {
                HSRLen++;
            }
            else
            { 
                if (  HSRLen >= 1  )
                { 
                    int Num_Segment = 0;//the number of segments in an observed element
                    int SubLen = 0;  /*length of the sub segment*/
                    //Min_overlap[]--the length of overlapping
                    //Mam_extend[] -- the extended length of two segments to be compared
                    for (j=i-HSRLen; j<i; j++)
                    { 
                        if (  predSec[j] == thisHSRstate  )
                        { 
                            SubLen++;
                        } 
                        else
                        { 
                            if (  SubLen >= 1  )
                            {
                                Min_overlap[Num_Segment] = SubLen;
                                Num_Segment++;
                            } 
                            SubLen = 0;
                        } 
                    }
                    if (  SubLen >= 1  )
                    {
                        Min_overlap[Num_Segment] = SubLen;
                        Num_Segment++;
                    }
                    double sum = 0;
                    int beg = 0;
                    int Nextend = 0;
                    int Num_HSRLen = HSRLen;
                    if  (  Num_Segment >= 1  )
                    {
                        for (j=0; j<Num_Segment; j++)
                        {
                            Max_extend[j] = HSRLen;
                            LenSeg_pred[j] = Min_overlap[j];
                        }
                        //for first segment
                        if (  predSec[i-HSRLen] == thisHSRstate  )
                        { 
                            beg = i - HSRLen - 1;
                            Nextend = 0;
                            while (   (beg>=0) && (predSec[beg]==thisHSRstate) && (obsSec[beg]!=thisHSRstate)   )
                            {
                                Nextend++;
                                beg--;
                            }
                            Max_extend[0] = Max_extend[0] + Nextend;
                            LenSeg_pred[0] = LenSeg_pred[0] + Nextend;
                        }
                        //for last segment
                        if (  predSec[i-1] == thisHSRstate  )
                        { 
                            beg = i;
                            Nextend = 0;
                            while (   (beg<Length) && (predSec[beg]==thisHSRstate) && (obsSec[beg]!=thisHSRstate)   )
                            {
                                Nextend++;
                                beg++;
                            }
                            Max_extend[Num_Segment-1] = Max_extend[Num_Segment-1] + Nextend;
                            LenSeg_pred[Num_Segment-1] = LenSeg_pred[Num_Segment-1] + Nextend;
                        }
                        for (j=0; j<Num_Segment; j++)
                        {
                            int Nmin1 = min( (Max_extend[j]-Min_overlap[j]), Min_overlap[j]);
                            int Nmin2 = min( HSRLen/2, LenSeg_pred[j]/2);
                            int Nmin = min(Nmin1, Nmin2);
                            double score_seg =  (Min_overlap[j]+Nmin)*1.0/double(Max_extend[j]) ;
                            sum = sum + score_seg;
                        }
                        sum = sum*HSRLen;
                        Num_HSRLen = HSRLen*Num_Segment;
                    }

                    All_HSRLen[iterHSRState] = All_HSRLen[iterHSRState] + Num_HSRLen;
                    All_Sum[iterHSRState] = All_Sum[iterHSRState] + sum;
                } 
                HSRLen = 0;
            }  
        } 
    }

    double sumOverlap = 0;
    int sumAllRes = 0;
    for (i=0; i<numState; i++)
    {
        numHSR[i] = All_HSRLen[i];
        numCorrectHSR[i] = All_Sum[i];
        sovHSR[i] = All_Sum[i]/double(All_HSRLen[i]+1e-6);
        sumAllRes  += All_HSRLen[i];
        sumOverlap += All_Sum[i];
    }
    numHSR[numState] = sumAllRes;
    numCorrectHSR[numState] = sumOverlap;
    sovHSR[numState] = sumOverlap / double(sumAllRes+1e-6);

    return numHSR[numState];
}/*}}}*/
int CalQ3(char *obsSec, char *predSec, int seqLength, int *numCorrectHSR, int *numHSR, double *q3HSR, int numState = 3)/*{{{*/
    /*=======================================================================
     * numCorrectHSR[4]
     * numHSR[4] Helix, Sheet, Coil, All
     * q3HSR[4]
     * numState = 3 by default
     *====================================================================*/
{

    int iterHSRState = 0; /*iterator for the HSR state, H, S and R*/
    char thisHSRstate = ' '; /*current HSR state*/
    int i = 0;

    for (i=0;i<=numState; i++)
    {
        numCorrectHSR[i] = 0;
        numHSR[i] = 0;
    }

    for (iterHSRState=0; iterHSRState<3; iterHSRState++)//corrsponding to H, S, R
    {
        thisHSRstate = HSRalphabet[iterHSRState];
        for(i=0; i<seqLength; i++)
        { 
            if (  obsSec[i] == thisHSRstate  )
            {
                numHSR[iterHSRState] ++;
                if (obsSec[i] == predSec[i])
                {
                    numCorrectHSR[iterHSRState] ++;
                }
            }
        } 
    }

    for (i=0; i<numState; i++)
    {   
        q3HSR[i] = numCorrectHSR[i] / double (numHSR[i]+1e-6);
        numCorrectHSR[numState]  += numCorrectHSR[i];
        numHSR[numState]  += numHSR[i];
    }
    q3HSR[numState] = numCorrectHSR[numState] / double (numHSR[numState]+1e-6); 

    return numHSR[numState];
}/*}}}*/
int CalQ3ConfBin(char *obsSec, char *predSec, double *predConf, int seqLength, int **numCorrectHSR, int **numHSR, double **q3HSR, int numState , double halfBinWidth)/*{{{*/
/*=======================================================================
 * calculate Q3 for residues lying in each confidence bin
 * numCorrectHSR[N][4]: N is the number of bins, and [4] is for Helix, Sheet, Coil and All
 * numHSR[N][4] Helix, Sheet, Coil, All
 * q3HSR[N][4]
 * numState = 3 by default
 *
 * predConf is in the range of [0,100]
 *
 * bins are setting like this. 
 * for numCorrectHSR[50][j], residues with predConf >= 50 - halfBinWidth and
 * predConf < 50 + halfBinWidth are included
 * predConf has two format: 
 *  1: rawConfidence: raw confidence for secondary structure is calculated as
 *  maxProb / sumProb
 *  2: nomalized confidence:  the normalization function is generated by
 *  second order polynomial regression 
 *  2009-07-31 
 *
 *====================================================================*/
{

    int iterHSRState = 0; /*iterator for the HSR state, H, S and R*/
    char thisHSRstate = ' '; /*current HSR state*/
    int i = 0;

    int iBin;         /*iterator of the bin*/
    for (i=0;i<=numState; i++)
    {
        for (iBin = 0; iBin <=100; iBin++)
        {
            numCorrectHSR[iBin][i] = 0;
            numHSR[iBin][i] = 0;
        }
    }

    int iBinBeg = 0;
    int iBinEnd = 0;
    for (iterHSRState=0; iterHSRState<3; iterHSRState++)//corrsponding to H, S, R
    {
        thisHSRstate = HSRalphabet[iterHSRState];
        for(i=0; i<seqLength; i++)
        { 
            iBinBeg = int(ceil(predConf[i]-halfBinWidth));
            iBinEnd = int(floor(predConf[i]+halfBinWidth));
            /* so if predConf = 35.1, halfBinWidth = 0.8
             * iBinBeg = ceil(35.1-0.8) = ceil(34.3) = 35
             * iBinEnd = floor(35.1+0.8) = floor(35.9) = 35
             * then this residue lies only in the bin [35]
             * */
            
            for (iBin = iBinBeg ; iBin <= iBinEnd; iBin++)
            { 
                if (  obsSec[i] == thisHSRstate  )
                {
                    numHSR[iBin][iterHSRState] ++;
                    if (obsSec[i] == predSec[i])
                    {
                        numCorrectHSR[iBin][iterHSRState] ++;
                    }
                }
            }
        } 
    }

    for (iBin = 0; iBin <= 100; iBin++)
    {
        for (i=0; i<numState; i++)
        {   
            q3HSR[iBin][i] = numCorrectHSR[iBin][i] / double (numHSR[iBin][i]+1e-6);
            numCorrectHSR[iBin][numState]  += numCorrectHSR[iBin][i];
            numHSR[iBin][numState]  += numHSR[iBin][i];
        }
        q3HSR[iBin][numState] = numCorrectHSR[iBin][numState] / double (numHSR[iBin][numState]+1e-6); 
    }

    return EXIT_SUCCESS;
}/*}}}*/

int CalQ3ConfAccumulate(char *obsSec, char *predSec, double *predConf, int seqLength, int **numCorrectHSR, int **numHSR, double **q3HSR, int numState , double halfBinWidth)/*{{{*/
/*=======================================================================
 * calculate Q3 for residues lying above a certain confidence
 * numCorrectHSR[N][4]: N is the number of bins, and [4] is for Helix, Sheet, Coil and All
 * numHSR[N][4] Helix, Sheet, Coil, All
 * q3HSR[N][4]
 * numState = 3 by default
 *
 * predConf is in the range of [0,100]
 *
 * bins are setting like this. 
 * for numCorrectHSR[50][j], residues with predConf >= 50 are included
 * predConf has two format: 
 *  1: rawConfidence: raw confidence for secondary structure is calculated as
 *  maxProb / sumProb
 *  2: nomalized confidence:  the normalization function is generated by
 *  second order polynomial regression 
 *  2009-07-31 
 *
 *====================================================================*/
{

    int iterHSRState = 0; /*iterator for the HSR state, H, S and R*/
    char thisHSRstate = ' '; /*current HSR state*/
    int i = 0;

    int iBin;         /*iterator of the bin*/
    for (i=0;i<=numState; i++)
    {
        for (iBin = 0; iBin <=100; iBin++)
        {
            numCorrectHSR[iBin][i] = 0;
            numHSR[iBin][i] = 0;
        }
    }

    int iBinBeg = 0;
    int iBinEnd = 0;
    for (iterHSRState=0; iterHSRState<3; iterHSRState++)//corrsponding to H, S, R
    {
        thisHSRstate = HSRalphabet[iterHSRState];
        for(i=0; i<seqLength; i++)
        { 
            iBinBeg = 0;
            iBinEnd = int(floor(predConf[i]));;
            /* so if predConf = 35.1, 
             * */
            
            for (iBin = iBinBeg ; iBin <= iBinEnd; iBin++)
            { 
                if (  obsSec[i] == thisHSRstate  )
                {
                    numHSR[iBin][iterHSRState] ++;
                    if (obsSec[i] == predSec[i])
                    {
                        numCorrectHSR[iBin][iterHSRState] ++;
                    }
                }
            }
        } 
    }

    for (iBin = 0; iBin <= 100; iBin++)
    {
        for (i=0; i<numState; i++)
        {   
            q3HSR[iBin][i] = numCorrectHSR[iBin][i] / double (numHSR[iBin][i]+1e-6);
            numCorrectHSR[iBin][numState]  += numCorrectHSR[iBin][i];
            numHSR[iBin][numState]  += numHSR[iBin][i];
        }
        q3HSR[iBin][numState] = numCorrectHSR[iBin][numState] / double (numHSR[iBin][numState]+1e-6); 
    }

    return EXIT_SUCCESS;
}/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;

    if(argc < 2) 
    {
        PrintHelp(argv);
        return 0;
    }
    int i,j;
    char outfile[MAX_PATH+1] = "";
    char fileListFile[MAX_PATH+1] ="";
    const char control_option[] = ""; //options which control the program, and does not take parameters
    //bool isAll = false;
    //bool isQuiet = false;
    //bool isSingle = false;
    bool isOutputQ3 = true;
    int fileFormat = 1;
    char chkFstFolder[MAX_PATH+1] ="";
    
    set <string> ::iterator iss;
    set <string> fileList_set;

    char infile [MAX_PATH+1] = "";
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < int(strlen(argv[i])); j++)
                {
                    switch (argv[i][j])
                    {
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[i][j]); return -1;
                    }
                }
                i ++;
            }
            else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(argv); 
                return 0;
            }
            else if(strcmp(argv[i],"-H") == 0 )
            {
                PrintVerboseHelp();
                return 0;
            }
            else if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i], "--out") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-l") == 0) || (strcmp(argv[i], "--list") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, fileListFile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-chkfstdir") == 0) || (strcmp(argv[i], "--chkfstdir") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, chkFstFolder)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-q3") == 0) || (strcmp(argv[i], "--isOutQ3") == 0))  
            {
                char tmpstr[MAX_PATH+1] ="";
                if( ( i = option_parser_filename(argc, argv, i, tmpstr)) == -1)
                    return -1;
                if (strncasecmp(tmpstr,"yes", 1) == 0)
                {
                    isOutputQ3 = true;
                }
                else
                {
                    isOutputQ3 = false;
                }
            }
            else if( (strcmp(argv[i],"-polypara") == 0) || (strcmp(argv[i], "--polypara") == 0))  
            {
                if (!IsNumeric(argv[i+1]) || !IsNumeric(argv[i+2]) || !IsNumeric(argv[i+3]) )
                {
                    fprintf(stderr,"Error! Option -polypara should be followed by three floating values\n");
                }
                else
                {
                    poly_a1 = atof(argv[i+1]);
                    poly_a2 = atof(argv[i+2]);
                    poly_a3 = atof(argv[i+3]);
                }
                i = i + 4;
            }
            else if( (strcmp(argv[i], "-f") == 0) || strcmp(argv[i], "--format") == 0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, fileFormat, true, 0, 4)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-max") == 0) || strcmp(argv[i], "--max-length") == 0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, maxSeqLength, true, 0, MAX_INT)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-statq3") == 0) )  
            {
                isStatQ3Conf = true;
                i ++;
            }
            else if( (strcmp(argv[i], "-proof") == 0) )  
            {
                isProofReading = true;
                i ++;
            }
            else if( (strcmp(argv[i], "-proofmethod") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, method_ProofReading, true, 0, 10)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-selconf") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, typeConfidence, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-binwidth") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, binWidth, true, 0.0, 100.0)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-m") == 0) || strcmp(argv[i], "--method") == 0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, method_HSR, true, 0, 4)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else //non-option argument
        {
            fileList_set.insert(argv[i]);
            i ++;
        }
    }/*}}}*/
    if (isStatQ3Conf)
    {
        if (strcmp(chkFstFolder, "") == 0)
        {
            fprintf (stderr,"chkFstFolder not set, since -q3stat is enabled. Exit!\n");
            return -1;
        }
    }

    double halfBinWidth = binWidth/2;

    if (strcmp(fileListFile , "") != 0)
    {
        FILE *fpList=fopen(fileListFile, "r");
        checkfilestream(fpList, fileListFile,"r", false);
        if (fpList != NULL)
        {
            while((linesize=fgetline(fpList, line, maxline)) != EOF)
            {
                fileList_set.insert(line);
            }
            fclose(fpList);
        }
    }
    FILE *fpout = NULL;
    if(strcmp(outfile,"") == 0 || strcasecmp(outfile, "stdout") == 0)
    {
        fpout = stdout;
    }
    else
    {
        fpout = fopen(outfile, "w");
        checkfilestream(fpout, outfile,"w");
    }
    /*add code here*/

    /*============================= main procedure goes here ================ */


    if(fileList_set.size()> 0 )
    {
        /*print out the header annotation*/
        fprintf(fpout,"%-10s %8s %6s %6s %6s %7s %7s %6s %6s %6s ", "Item", "All(SOV)", "H(SOV)", "S(SOV)", "R(SOV)", "numAll", "numH", "numS", "numR", "seqLen");
        if(isOutputQ3)
        {
            fprintf(fpout,"%7s %6s %6s %6s %7s %6s %6s %6s", "All(Q3)", "H(Q3)", "S(Q3)", "R(Q3)", "numAll", "numH", "numS", "numR");
        }
        fprintf(fpout,"\n");
    }
    Array1D <char> aaSeq_1darray(maxSeqLength);
    Array1D <char> obsSec_1darray(maxSeqLength);
    Array1D <char> predSec_1darray(maxSeqLength);
    aaSeq_1darray.Init('-');
    obsSec_1darray.Init('-');
    predSec_1darray.Init('-');
    char *aaSeq = aaSeq_1darray.array1D;
    char *obsSec = obsSec_1darray.array1D;
    char *predSec = predSec_1darray.array1D;


    Array1D <double> q3HSR_1darray(NUM_HSR_STATE+1);
    Array1D <double> sovHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> numQ3CorrectHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> numQ3HSR_1darray(NUM_HSR_STATE+1);
    Array1D <double> numSOVCorrectHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> numSOVHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> totalNumQ3CorrectHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> totalNumQ3HSR_1darray(NUM_HSR_STATE+1);
    Array1D <double> totalQ3HSR_1darray(NUM_HSR_STATE+1);
    Array1D <double> totalNumSOVCorrectHSR_1darray(NUM_HSR_STATE+1);
    Array1D <int> totalNumSOVHSR_1darray(NUM_HSR_STATE+1);
    Array1D <double> totalSOVHSR_1darray(NUM_HSR_STATE+1);

    q3HSR_1darray.Init(0);
    sovHSR_1darray.Init(0);
    numQ3CorrectHSR_1darray.Init(0);
    numQ3HSR_1darray.Init(0);
    numSOVCorrectHSR_1darray.Init(0);
    numSOVHSR_1darray.Init(0);
    totalNumQ3CorrectHSR_1darray.Init(0);
    totalNumQ3HSR_1darray.Init(0);
    totalQ3HSR_1darray.Init(0);
    totalNumSOVCorrectHSR_1darray.Init(0);
    totalNumSOVHSR_1darray.Init(0);
    totalSOVHSR_1darray.Init(0);


    double *q3HSR = q3HSR_1darray.array1D;
    double *sovHSR = sovHSR_1darray.array1D;

    int *numQ3CorrectHSR = numQ3CorrectHSR_1darray.array1D;
    int *numQ3HSR= numQ3HSR_1darray.array1D;   /*number of residues for each state when for analysis of Q3*/  
    double *numSOVCorrectHSR=numSOVCorrectHSR_1darray.array1D;
    int *numSOVHSR = numSOVHSR_1darray.array1D; /*number of residues for each state when for analysis of SOV*/

    int *totalNumQ3CorrectHSR = totalNumQ3CorrectHSR_1darray.array1D;
    int *totalNumQ3HSR = totalNumQ3HSR_1darray.array1D;
    double *totalQ3HSR = totalQ3HSR_1darray.array1D;

    double *totalNumSOVCorrectHSR = totalNumSOVCorrectHSR_1darray.array1D;
    int *totalNumSOVHSR = totalNumSOVHSR_1darray.array1D;
    double *totalSOVHSR = totalSOVHSR_1darray.array1D;

    /*variable for analysis of the q3 within each confidence bin*/
    Array1D <double> predConf_1darray(maxSeqLength);

    int numBin = 101;

    Array2D <int> statNumQ3CorrectHSR_2darray(numBin ,NUM_HSR_STATE+1);
    Array2D <int> statNumQ3HSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <double> statQ3HSR_2darray(numBin,NUM_HSR_STATE+1);

    Array2D <int> statTotalNumQ3CorrectHSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <int> statTotalNumQ3HSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <double> statTotalQ3HSR_2darray(numBin, NUM_HSR_STATE+1);

    Array2D <int> statAccumulateNumQ3CorrectHSR_2darray(numBin ,NUM_HSR_STATE+1);
    Array2D <int> statAccumulateNumQ3HSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <double> statAccumulateQ3HSR_2darray(numBin,NUM_HSR_STATE+1);

    Array2D <int> statTotalAccumulateNumQ3CorrectHSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <int> statTotalAccumulateNumQ3HSR_2darray(numBin, NUM_HSR_STATE+1);
    Array2D <double> statTotalAccumulateQ3HSR_2darray(numBin, NUM_HSR_STATE+1);


    predConf_1darray.Init(0.0);

    statNumQ3CorrectHSR_2darray.Init(0);
    statNumQ3HSR_2darray.Init(0);
    statQ3HSR_2darray.Init(0.0);

    statTotalNumQ3CorrectHSR_2darray.Init(0);
    statTotalNumQ3HSR_2darray.Init(0);
    statTotalQ3HSR_2darray.Init(0);

    statAccumulateNumQ3CorrectHSR_2darray.Init(0);
    statAccumulateNumQ3HSR_2darray.Init(0);
    statAccumulateQ3HSR_2darray.Init(0.0);

    statTotalAccumulateNumQ3CorrectHSR_2darray.Init(0);
    statTotalAccumulateNumQ3HSR_2darray.Init(0);
    statTotalAccumulateQ3HSR_2darray.Init(0);

    double *predConf =predConf_1darray.array1D;

    int **statNumQ3CorrectHSR = statNumQ3CorrectHSR_2darray.array2D;
    int **statNumQ3HSR = statNumQ3HSR_2darray.array2D;
    double **statQ3HSR = statQ3HSR_2darray.array2D;

    int **statTotalNumQ3CorrectHSR = statTotalNumQ3CorrectHSR_2darray.array2D;
    int **statTotalNumQ3HSR = statTotalNumQ3HSR_2darray.array2D;
    double **statTotalQ3HSR = statTotalQ3HSR_2darray.array2D;

    int **statAccumulateNumQ3CorrectHSR = statAccumulateNumQ3CorrectHSR_2darray.array2D;
    int **statAccumulateNumQ3HSR = statAccumulateNumQ3HSR_2darray.array2D;
    double **statAccumulateQ3HSR = statAccumulateQ3HSR_2darray.array2D;

    int **statTotalAccumulateNumQ3CorrectHSR = statTotalAccumulateNumQ3CorrectHSR_2darray.array2D;
    int **statTotalAccumulateNumQ3HSR = statTotalAccumulateNumQ3HSR_2darray.array2D;
    double **statTotalAccumulateQ3HSR = statTotalAccumulateQ3HSR_2darray.array2D;


    int iBin = 0;
    int seqLength = 0;
    int totalLength = 0;     /*total number of residues for multiple files, even with those without structure definition*/

    char rtname [MAX_PATH+1] = "";
    char chkFstResFile [MAX_PATH+1] = "";

    for(iss = fileList_set.begin(); iss != fileList_set.end(); iss ++)
    {
        my_strcpy( infile, (*iss).c_str(), MAX_PATH-1);
        rootname(infile, rtname);
        seqLength = ReadInSecPredFile(infile, aaSeq, obsSec, predSec,predConf, typeConfidence, fileFormat);

        if (isProofReading)
        {
            ProofReading(predSec, seqLength, method_ProofReading);
        }
        /*multiply predConf by 100.0, to be percentages*/

        if (isStatQ3Conf)
        {
            sprintf(chkFstResFile, "%s/%s.txt", chkFstFolder, rtname);
            ReadConf(chkFstResFile, predConf, typeConfidence);
            for(i = 0; i < seqLength; i ++) { predConf[i] *= 100.0; }
        }

        CalSOV(obsSec, predSec, seqLength, numSOVCorrectHSR, numSOVHSR, sovHSR, NUM_HSR_STATE);
        fprintf(fpout,"%-10.10s %8.2f %6.2f %6.2f %6.2f %7d %7d %6d %6d %6d ", rtname,  sovHSR[3]*100.0,  sovHSR[0]*100.0, sovHSR[1]*100.0, sovHSR[2]*100.0, numSOVHSR[3],numSOVHSR[0], numSOVHSR[1],numSOVHSR[2], seqLength);
        /*print out the result*/
        if (isOutputQ3)
        {
            CalQ3(obsSec, predSec, seqLength, numQ3CorrectHSR, numQ3HSR, q3HSR, NUM_HSR_STATE);
            fprintf(fpout,"%7.2f %6.2f %6.2f %6.2f %7d %6d %6d %6d", q3HSR[3]*100.0,  q3HSR[0]*100.0,  q3HSR[1]*100.0, q3HSR[2]*100.0, numQ3HSR[3], numQ3HSR[0],numQ3HSR[1],numQ3HSR[2]);
        }
        fprintf(fpout,"\n");

        if (isStatQ3Conf) /*analyze Q3 at different confidence; 2009-08-02*/
        {
            CalQ3ConfBin(obsSec, predSec, predConf, seqLength, statNumQ3CorrectHSR, statNumQ3HSR, statQ3HSR, NUM_HSR_STATE, halfBinWidth);
            CalQ3ConfAccumulate(obsSec, predSec, predConf, seqLength, statAccumulateNumQ3CorrectHSR, statAccumulateNumQ3HSR, statAccumulateQ3HSR, NUM_HSR_STATE, halfBinWidth);
            for (i = 0; i< NUM_HSR_STATE+1; i ++)
            {
                for (iBin = 0; iBin <= 100; iBin ++)
                {
                    statTotalNumQ3CorrectHSR[iBin][i] += statNumQ3CorrectHSR[iBin][i];
                    statTotalNumQ3HSR[iBin][i] += statNumQ3HSR[iBin][i];

                    statTotalAccumulateNumQ3CorrectHSR[iBin][i] += statAccumulateNumQ3CorrectHSR[iBin][i];
                    statTotalAccumulateNumQ3HSR[iBin][i] += statAccumulateNumQ3HSR[iBin][i];
                }
            }

        }


        for (i = 0; i< NUM_HSR_STATE+1; i++)
        {
            totalNumSOVCorrectHSR[i] += numSOVCorrectHSR[i];
            totalNumSOVHSR[i] += numSOVHSR[i];

            totalNumQ3CorrectHSR[i] += numQ3CorrectHSR[i];
            totalNumQ3HSR[i] += numQ3HSR[i];

        }
        totalLength += seqLength;
        
    }
    /*output the summary results*/
    for (i = 0; i< NUM_HSR_STATE+1; i++)
    {
        totalSOVHSR[i] = totalNumSOVCorrectHSR[i] / double (totalNumSOVHSR[i]+1e-6);
        totalQ3HSR[i] = totalNumQ3CorrectHSR[i] / double (totalNumQ3HSR[i]+1e-6);
    }
    if (isStatQ3Conf)
    {
        for (i = 0; i< NUM_HSR_STATE+1; i++)
        {
            for (iBin = 0; iBin <= 100; iBin ++)
            {
                statTotalQ3HSR[iBin][i] = statTotalNumQ3CorrectHSR[iBin][i] / double (statTotalNumQ3HSR[iBin][i]+1e-6);
                statTotalAccumulateQ3HSR[iBin][i] = statTotalAccumulateNumQ3CorrectHSR[iBin][i] / double (statTotalAccumulateNumQ3HSR[iBin][i]+1e-6);
            }
        }
        /*print out*/
        fprintf(stdout,"# BinWidth = %f\n", binWidth);
        fprintf(stdout,"#%28s    %20s\n", "non-accumulate" , "accumulate");
        fprintf(stdout,"#%-7s %6s %7s %7s   %8s %7s %7s\n", "Conf" , "Q3", "numRes", "%numRes", "Q3", "numRes", "%numRes");

        double percentNumResBin =0.0;
        double percentNumResAccumulate = 0.0;
        for(iBin = 0; iBin <= 100; iBin ++)
        {
            percentNumResBin = statTotalNumQ3HSR[iBin][3]/(statTotalAccumulateNumQ3HSR[0][3]+1e-6)*100.0;
            percentNumResAccumulate = statTotalAccumulateNumQ3HSR[iBin][3]/(statTotalAccumulateNumQ3HSR[0][3]+1e-6)*100.0;
            
            fprintf(stdout,"%-8d %6.2f %7d %6.2f   %8.2f %7d %6.2f\n", iBin, statTotalQ3HSR[iBin][3]*100.0, statTotalNumQ3HSR[iBin][3],percentNumResBin,  statTotalAccumulateQ3HSR[iBin][3]*100.0, statTotalAccumulateNumQ3HSR[iBin][3], percentNumResAccumulate);
        }

    }

    fprintf(fpout,"%-10s %8.2f %6.2f %6.2f %6.2f %7d %7d %6d %6d %6d ", "Sum",  totalSOVHSR[3]*100.0,  totalSOVHSR[0]*100.0, totalSOVHSR[1]*100.0, totalSOVHSR[2]*100.0, totalNumSOVHSR[3],totalNumSOVHSR[0], totalNumSOVHSR[1],totalNumSOVHSR[2], totalLength);
    if (isOutputQ3)
    {
        fprintf(fpout,"%7.2f %6.2f %6.2f %6.2f %7d %6d %6d %6d", totalQ3HSR[3]*100.0,  totalQ3HSR[0]*100.0,  totalQ3HSR[1]*100.0, totalQ3HSR[2]*100.0, totalNumQ3HSR[3], totalNumQ3HSR[0],totalNumQ3HSR[1],totalNumQ3HSR[2]);
    }
    fprintf(fpout,"\n");

    /*============================= main procedure ends here ================ */



    if(fpout != NULL && fpout != stdout) fclose(fpout);

    return 0;
}
/*}}}*/

