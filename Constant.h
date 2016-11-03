/*
 * =====================================================================================
 *       Filename:  Constant.h
 *    Description:  Define commonly used constant
 *        Version:  1.0
 *       Compiler:  g++
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

//use preprocessor directives in that file to prevent it from being inadvertently included more than once.
#ifndef HAS_CONSTANT_H
#define HAS_CONSTANT_H


#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef NULL
#define NULL 0
#endif

#ifndef MAX_PATH
#define MAX_PATH 500
#endif

#ifndef MAX_COMMAND_LINE
#define MAX_COMMAND_LINE 501
#endif

#ifndef WHITE_SPACE
#define WHITE_SPACE " \t\r\n"
#endif


#ifndef NULLCHAR
#define NULLCHAR '\0' //define NULLCHAR as '\0'
#endif

#ifdef _UNICODE
#define _TCHAR wchar_t
#else
#define _TCHAR char
#endif

#ifndef ERROR_CODE
#define READ_FILE_ERROR  1
#define WRITE_FILE_ERROR 2
//#define READ_OVERFLOW    5
//#define WRITE_OVERFLOW   5
//#define UNDERFLOW        10
#endif /*ERROR_CODE*/

// define some MAX and MIN value for different data types
#define MAX_INT         0x7FFFFFFF
#define MIN_INT         0x80000000

#define MAX_FLOAT       3.402823466E+38
#define MIN_FLOAT      -3.402823466E+38

#define MAX_DOUBLE      1.7976931348623157E308
#define MIN_DOUBLE     -1.7976931348623157E308


#define INIT_INT	      123456789
#define INIT_FLOAT	float(123456789.0)
#define INIT_DOUBLE	      123456789.0



// define some universal constant
#ifndef PI
#define PI 3.1415926535897932384626433832795028841972 
#endif

#ifndef E
#define E 2.7182818284590452353602874713526624977573 
#endif

#define _CONST_PI  3.1415926535897932384626433832795028841972 
#define _CONST_E   2.7182818284590452353602874713526624977573 
#define _CONST_h   6.626068e-34 // planck's constant,       unit: m2 kg / s
#define _CONST_G   6.6742e-11   // gravitation constant,    unit: m3·kg-1·s-2 
#define _CONST_c   299792458.0  // speed of light in vacumn unit: m/s 
#define _CONST_atm 101325       // standard atmosphere      unit: Pa 
#define _CONST_Na  6.022e23     // avogadro's constant      unit: no unit


//constant for SCOP file
#define SIZE_SCOP_ID  7              // length of the scop domain id, e.g. d9icwa1
#define SIZE_DOMAIN_DEF_RECORD   40  // maximal length of the scop domain definition, e.g. "A:205-288"
#define NUM_CHAIN_PER_DOMAIN     5   // maximal number of chain can be included in on scop domain
#define NUM_DOMAIN_PDB           200 // max number of domains for one pdbid	

// constant for PDB file 
#define SIZE_LINE_PDB            80
#define SIZE_PDBID               4

/*****************************************************************************
 * The SIZE_CHAIN_ID is set to a larger value, for user's input.
 ****************************************************************************/
#define SIZE_CHAIN_ID            30

#define SIZE_RECORD_ID           6
#define SIZE_RES_NAME            3
#define SIZE_METAL_ATOM_NAME     2
#define SIZE_METAL_ATOM_RES_NAME 4
#define MAX_ATOM_SERIAL          99999

//constant for positions of item on PDB ATOM record line, starting from 0
//record name can be retrieved by line[POS_RECORD_BEG..POS_RECORD_END-1]
#define POS_RECORD_BEG     0
#define POS_RECORD_END     6

#define POS_SERIAL_BEG     6
#define POS_SERIAL_END     11

#define POS_ATOMNAME_BEG   12
#define POS_ATOMNAME_END   16

#define POS_ALTLOC         16

#define POS_RESNAME_BEG    17
#define POS_RESNNAME_END   20

#define POS_CHAIN_ID       21

#define POS_RESSEQ_BEG     22
#define POS_RESSEQ_END     26

#define POS_ICODE          26

#define POS_X_COOR_BEG     30
#define POS_X_COOR_END     38
#define POS_Y_COOR_BEG     38
#define POS_Y_COOR_END     46
#define POS_Z_COOR_BEG     46
#define POS_Z_COOR_END     54

#define POS_OCCUPANCY_BEG  54
#define POS_OCCUPANCY_END  60
#define POS_TEMPFACTOR_BEG 60
#define POS_TEMPFACTOR_END 66
#define POS_SEGID_BEG      72
#define POS_SEGID_END      76
#define POS_ELEMENT_BEG    76
#define POS_ELEMENT_END    78
#define POS_CHARGE_BEG     78
#define POS_CHARGE_END     80


// constant for some max values
#define MAX_SSBOND_CHAIN     100  // maximal number of ssbonds per chain
#define MAX_SSBOND_PRO       900  // maximal number of ssbonded proteins

#define NUM_METAL_ELEMENT    100
#define MAX_ATOM_PER_RES     25   // maximal number of atoms per residue, only one location from alternative location will be recorded
#define MAX_ATOM_PER_RES_WITH_ALTATOM 100 // maximal number of atoms per residue, considering alternative locations
#define MAX_ALTATOM          20   // maximal number of alternative atoms in PDB, when altLoc != ' '

#define MAX_NUM_PDB          30000// maximal number of entries of PDB
#define MAX_NUM_NRPDB        5000 // maximal number of entries of nrPDB
#define MAX_NRPDB_LIST       9000 // maximal number of chains in nrPDB list

#define LONGEST_SEQ          200000 // longest sequence
#define LONGEST_SHAPE        LONGEST_SEQ // longest shape string
#define MAX_SEQ_LENGTH       LONGEST_SEQ // maximal length of amino acid sequence

#define MAX_BONDED_RES_CHAIN 1000 // maximal number of bonded residues per chain
#define MAX_METAL_CHAIN      1000 // maximal number of metal atoms per chain
#define MAX_CLOSE_RESIDUE    100  // maximal number of residues close to a metal atom

//Atom record in pdb file
#define SIZE_TITLE         SIZE_RECORD_ID
#define SIZE_ATOM_NAME     4

#define SIZE_ATOM_ORIGNAME 4
#define SIZE_ATOM_SEGID    4
#define SIZE_ATOM_ELEMENT  2
#define SIZE_ATOM_CHARGE   2

#endif  //HAS_CONSTANT_H
