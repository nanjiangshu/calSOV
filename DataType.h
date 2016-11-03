/*
 * =====================================================================================
 *       Filename:  DataType.h
 *    Description:  Define commonly used datatype
 *        Version:  1.0
 *       Compiler:  g++
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */
// typedef signed char __int8;
// typedef signed long long __int64;


//use preprocessor directives in that file to prevent it from being inadvertently included more than once.
#ifndef HAS_DATATYPE_H
#define HAS_DATATYPE_H

typedef int                 BOOL;
typedef unsigned char       BYTE;

typedef signed char         int8;
typedef unsigned char       unit8;
typedef signed short        int16;
typedef unsigned short      unit16;
typedef signed int          int32;
typedef unsigned int        unit32;
typedef signed long long    int64;
typedef unsigned long long  unit64;

typedef int8       __int8;
typedef int8        _int8;
typedef unit8      __unit8;
typedef unit8       _unit8;
typedef int16      __int16;  
typedef int16       _int16;
typedef unit16     __unit16;
typedef unit16      _unit16;
typedef int32      __int32;
typedef int32       _int32;
typedef unit32     __unit32;
typedef unit32      _unit32;
typedef int64      __int64;
typedef int64       _int64;
typedef unit64     __unit64;
typedef unit64      _unit64;

typedef struct 
{
	double x;
	double y;
	double module;
	double pha;
}Complex;

#endif
