#ifndef HAS_MYFUNC_H
#define HAS_MYFUNC_H


#include <cassert>
#include <ctime>
#include <iostream>
#include <set>

#include <regex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Constant.h"
#include "DataType.h"

/*using namespace std;*/


//ChangeLog and Note/*{{{*/
/*****************************************************************************
 * Note 2008-04-23
 * function like macros are extremely dangerous, it will bring unexpected
 * consequences. Avoid using function like macros where arguments appear more
 * than once.
 *
 ****************************************************************************//*}}}*/

/*#if defined(linux) || defined (cygwin)*/
/*#define HAS_STRICMP*/
/*#define stricmp strcasecmp   */
/*#define _stricmp strcasecmp*/
/*#define _my_strupr my_strupr*/
/*#define _strdup strdup*/
/*#define _getcwd getcwd*/
/*#define _chdir chdir*/
/*#define __min MIN2*/
/*#define __max MAX2*/
/*#endif*/

// type of modm file
#ifndef MODM_FILE_TYPE
#define MODM_FILE_TYPE
#define MODM_LOG 0
#define MODM_PER 1
#endif

#ifndef INIT_CONSV
#define INIT_CONSV -99
#endif

#ifndef INIT_RESSEQ
#define INIT_RESSEQ -999
#endif

#ifndef INIT_AASEQINDEX
#define INIT_AASEQINDEX -1
#endif

#ifndef INIT_ICODE
#define INIT_ICODE '!'
#endif

#ifndef INIT_SHAPE
#define INIT_SHAPE '*'
#endif

#ifndef INIT_WATERACC
#define INIT_WATERACC  -1
#endif
// define the sequence type DNA or AA
#ifndef SEQ_TYPE
#define SEQ_TYPE
#define DNA_SEQ 0 // dna sequence
#define AA_SEQ  1 // amino acid sequence
#define SHAPE_SEQ 2 // shape string sequence
#define UNKNOWN_SEQ_TYPE -1
#endif


/*avoid using function like macros where variables appear twice 2008-04-23*/
/*#ifndef HAVE_MIN_MAX*/
/*#define HAVE_MIN_MAX*/
/*#define MIN2(x,y)	((x) < (y) ? (x) : (y))*/
/*#define MIN3(x,y,z) (MIN2((x),(y))<(z) ? (MIN2((x),(y))) : (z))*/
/*#define MAX2(x,y)   ((x)<(y) ? (y) : (x))*/
/*#define MAX3(x,y,z) (MAX2((x),(y))<(z) ? (z) : MAX2((x),(y)))*/
/*#endif*/

#define IS_EQUAL(x,y) (((x) == (y)) ? true : false)

/*#ifndef CheckValue*/
/*#define CheckValue*/
/*#define CheckBoolValue(x)    printf("%s = %s\n", #(x), ((x) == true) ?  "true" : "false");*/
/*#define CheckStringValue(x)  printf("%s = %s\n", #(x), (x));*/
/*#define CheckNumericValue(x) printf("%s = %g\n", #(x), double(x));*/
/*#endif [>!CheckValue<]*/

#ifndef SORT_ORDER
#define SORT_ORDER
#define ASCENDING 0
#define DESCENDING 1
#endif


#ifdef __cplusplus
// Do not define log2
#else
#define log2(x) (log(x) / M_LOG2_E)
#endif

/*#ifndef DEBUG_MESSAGE*/
/*#define DEBUG_MESSAGE*/
/*#define DEBUG_STDERR(x, fmt) fprintf(stderr,"%s:%u: %s=" fmt, __FILE__, __LINE__, #(x), (x))*/
/*#endif [>!DEBUG_MESSAGE<]*/


/*#ifndef BLOSUM_ALPHABET*/
/*#define BLOSUM_ALPHABET*/
/*char BLOSUM_alphabet_3L[24][4]={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","ASX","GLX","UNK"};*/
/*char BLOSUM_alphabet_1L[24]="ARNDCQEGHILKMFPSTWYVBZX"; // unknown <-> 'X'*/
/*// amino acids out of BLOSUM_alphabet_3L treat as unknown	*/
/*#endif*/


/*_io_stream*/
int checkfilestream(FILE *fp, const char* filename, const char *mode, bool isAssert = false);

/*file management*/
char *rootname(const char*filename, char*rtname, int max_rtname = MAX_PATH);
char *getfilepath(const char* filename, char *path, int max_path = MAX_PATH);
char *getfileext(const char* filename, char *ext, int max_ext = MAX_PATH);

void VerifyFolder(const char* folder);
int GetDataDir(char  *datadir);
int GetTMPDir(char   *tmpdir);
int GetWorkDir(char  *workdir);

char* GetPDBFilePath(const char* pdbid, char* pdbfilepath, const char pdbpath[] = "", const char pdbobsoletepath[] = "");
char* GetMODMFilePath(const char* id, char* MODMfilepath, const char modmpath[] = "");
char* GetPSSMFilePath(const char* id, char* pssmfilepath, const char pssmpath[] = "");
char* GetDSSPFilePath(const char* pdbid, char* dsspfilepath, const char dssppath[] = "");
char* GetPDBAAFilePath(const char* id, char* pdbaafilepath, const char pdbaapath[] = "");
char* GetSCOPAAFilePath(const char* id, char* scopaafilepath, const char scopaapath[] = "");
char* GetSEQMAPFilePath(const char* id, char* seqmapfilepath, const char seqmappath[] = "");
char* GetShapeStringFilePath(const char* id, char* shapestringfilepath, const char shapestringpath[] = "");

/*_io_ reading */
int   fgetlinesize(FILE* fp);
int   fgetline(FILE* fp, char* line, int max = 0x7FFFFFFF);
int   fgetdelim(FILE* fp, char* str, const char *delim = WHITE_SPACE, int max = 0x7FFFFFFF);
int   fgetlinecnt(const char* filename, int &maxline, bool is_count_blank_line = true);
int   fgetlinecnt(FILE* fp, bool is_count_black_line = true);
int   fgetlinecnt(const char* filename, bool is_count_black_line = true);
void  f_neglect_comment(FILE* fp, const char comment_char = '#');

/*string operation*/
int   my_strcpy(char *to, const char* from, int max);
char *my_strupr(char* str, int beg = 0, int end = 0x7FFFFFFF);
char *my_strlwr(char* str, int beg = 0, int end = 0x7FFFFFFF);
char *strrtrim(char *str, const char *trim = " \t\n\r");
char *strltrim(char *str, const char *trim = " \t\n\r");
char *strtrim (char *str, const char *trim = " \t\n\r");
char *ssubstitute(char* str, char chForSub, char chAfterSub,int start = 0, int end = 0x7FFFFFFF);
char *sreplace(char* to, char* from, int start = 0);
char *SpanExcluding(const char* strForSpan, char* strAfterSpan, const char* charSet = WHITE_SPACE);
char *strchomp(char *str);


/*binary checking*/
bool IsBlankLine(const char *buf);
bool IsInteger(double x);
bool IsZero(Complex *cplx, int dim);
bool IsNumeric(const char* str);
bool IsDigit(const char* str);
bool IsDigit(const char ch);
bool IsLower(const char c);
bool IsUpper(const char c);
bool IsInCharSet(const char ch, const char* charSet, int n = 0);

/* mathematical functions */
void FFT(Complex *cplx, int dim, int pow2, bool isInv);
int  BitSwap(int i, int pow2);
int  Integer(double x);
int GetNumDigit(int num);
void SmoothImage(unsigned short **image, int imageWidth, int imageHeight, int sizeN = 1);
void GetAmpPha(double re, double im, double& amp, double& pha);
double uniform_random(double min = 0.0, double max = 1.0);
double SigmoidScore(double a, double b, double x);
double GatingScore(double x, double y);


/*array operation, searching, sorting and shuffling*/
void  Shuffle(int* array, int n, unsigned int rand_seed = time(NULL));

template <class T> int   BinarySearch_String(T keyStr, T* strs, int n);
template <class T> int   LinearSearch_String(T keyStr, T* strs, int n) ;

void  BubbleSort(int* a, int n);
void QuickSort_String(int *idx, char** strs, int low, int high);

int Grouping(char **strs, int numStrs, char **strGroup, int *subTotal,int SIZE_STRGROUP);

/*bioinformatics, programming on protein or nucleic database,*/

int Char2Digit(char aa, const char* alphabet, int n = 0);
int Charcase2Digit(char aa, const char* alphabet, int n = 0);
int Digit2Char(int dc, const char* alphabet, int n = 0);

double Compute_ROC_score(int* label, int n);
double Compute_ROC50_score(int* label, int n);
double Compute_medianRFP_score(int* label, double* score, int n);
double Compute_medianRFP50_score(int* label, double* score, int n);

/*regular expression*/                                                                           
/* fint all matchs of the regular expression pattern in a string */
int reg_findall(const char* string, const char* pattern,  regmatch_t * pmatch ,bool isOverlap = false);

// functions for argument parser
int option_parser_filename(int argc, char **argv, int beg, char *filename);
template <class T> int option_parser_numeric(int argc, char **argv, int beg, T &x, bool isRangeSet = false, T min = MIN_INT, T max = MAX_INT); 

#endif //HAS_MYFUNC_H
