/*
 * =====================================================================================
 * 
 *       Filename:  mytemplate.h
 * 
 *    Description:  header file for mytemplate.cpp
 * 
 *        Version:  1.0
 *        Created:  11/21/2006 12:52:58 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */

#ifndef HAS_MYTEMPLATE_H
#define HAS_MYTEMPLATE_H

#include "myfunc.h"

template <class T> double Average(T *array, int low, int high);
template <class T> void Swap(T* x1, T* x2);
template <class T> void Swap(T& x1, T& x2);
template <class T> T Coverage(T a1, T b1, T a2, T b2);
template <class T> T max_element(T *array, int lo, int hi);
template <class T> int max_element_index(T *array, int lo, int hi);

template <class T> int binarysearch(T key, T *array, int n);
template <class T> void QuickSort(T *sort, int low, int high, int order = ASCENDING) ;
template <class T> void QuickSort_index(int *idx, T *sort, int low, int high, int order  = ASCENDING);
template <class T, class InputIterator> void Set2Array(InputIterator __First, InputIterator __Last, T *array);

template <class T> int locate_range(T key, T *array, int lo, int hi);

template <class T> void Swap(T* x1, T* x2)/*{{{*/
{// call by pointers
	T tmp;
	tmp = *x1;
	*x1 = *x2;
	*x2 = tmp;
}/*}}}*/
template <class T> void Swap(T& x1, T& x2)/*{{{*/
{ //call by references
	T tmp;
	tmp = x1;
	x1 = x2;
	x2 = tmp;
}/*}}}*/
template <class T> T Coverage(T a1, T b1, T a2, T b2)/*{{{*/
/*****************************************************************************
 * Coverage()
 * return the length of coverage for 2 subset. for example, coverage of [2,5]
 * and [3,8] returns 2
 ****************************************************************************/
{
	return (min(b1,b2) - max(a1,a2));
}
/*}}}*/
template <class T> double Average(T *array, int low, int high)/*{{{*/
///*******************************************************
// return the average of array[low..high]
///******************************************************
{
    double sum = 0.0;
    int n = high - low +1;
    int i;
    for( i = low; i <= high ; i++)
        sum += double(array[i]);
    return sum / double(n);
}/*}}}*/
template <class T> T max_element(T *array, int lo, int hi)/*{{{*/
/*****************************************************************************
 * return the maximum element of a range of an array
 * both array[lo] and array[hi] are accessed. so 
 * max_element(array, 0, n-1) means returning the max_element of the whole
 * array
 * Example: max = max_element(array, 0, n-1);
 ****************************************************************************/
{
    int i ;
    T max_ele = array[0];
    for(i = lo; i <= hi ; i ++ )
    {
        if(array[i] > max_ele)
        {
            max_ele = array[i];
        }
    }
    return max_ele;
}/*}}}*/
template <class T> int max_element_index(T *array, int lo, int hi)/*{{{*/
/*****************************************************************************
 * return the index of the maximum element of a range of an array
 * both array[lo] and array[hi] are accessed. so 
 * max_element(array, 0, n-1) means returning the max_element of the whole
 * array
 * Example: index_max_ele = max_element(array, 0, n-1);
 ****************************************************************************/
{
    int i ;
    int index_max_ele = lo;
    T max_ele = array[0];
    for(i = lo; i <= hi ; i ++ )
    {
        if(array[i] > max_ele)
        {
            max_ele = array[i];
            index_max_ele = i;
        }
    }
    return index_max_ele;
}/*}}}*/
template <class T> T Sum(T *array, int low, int high)/*{{{*/
///*******************************************************
// return the sum of array[low..high]
///******************************************************
{
    T sum = T(0);
    int i;
    for( i = low; i <= high ; i++)
        sum += array[i];
    return sum ;
}/*}}}*/

template <class T> int binarysearch(T key, T *array, int n)/*{{{*/
///*******************************************************
// BinarySearch()
// search the key in array in which the element can be compared directly.
// n: the number of items of strs
// return the index of key in array if successful
// else return -1
///******************************************************
{
    int lo = 0;
    int hi = n - 1;
    int mid;
    /* Repeat while there are elements in range */
    while (lo <= hi) 
    {
        mid = (lo + hi) / 2; /* Compute midpoint */
        if (array[mid] == key)            /* found element */
        { return(mid); } 
        else if (array[mid] > key)      /* target is in first half */
        { hi = mid - 1; } 
        else /* array[mid] < target , target is in second half */
        { lo = mid + 1; }
    }
    return(-1); /* Nothing left in range; failure */
}/*}}}*/
template <class T>  void QuickSort(T *sort, int low, int high, int order /* = ASCENDING*/)/*{{{*/
/*****************************************************************************
 * QuickSort()
 * sort the numeric array in range [low,high] in asccending order
 * using binary sort algorithm
 * content of array 'sort' is going to be changed
 * note: sort[high] is included, so if sizeof(sort) = n, low = 0 , high = n-1
 *       for sorting the whole array
 ****************************************************************************/
{
    T pivot;
    int m;
    int i;
    if(low < high)
    {
		Swap(&sort[low], &sort[(high+low)/2]);
	    pivot = sort[low];
		m = low;
	    for (i = low + 1; i <= high; i++)
		{
			if(((order == ASCENDING) && (sort[i] <  pivot))       
                   || (order == DESCENDING && (sort[i] > pivot)))
            /*if(sort[i] < pivot)*/
			{
				m++;
                Swap(&sort[m], &sort[i]);
			}
	   }
	   Swap(&sort[low], &sort[m]);
	   QuickSort(sort, low, m - 1, order);
	   QuickSort(sort, m + 1, high, order);
	}
}/*}}}*/
template <class T>  void QuickSort_index(int *idx, T *sort, int low, int high, int order /*= ASCENDING*/)/*{{{*/
/*****************************************************************************
 * QuickSort_index()
 * sort the numeric array in range [low,high] in asccending order
 * using binary sort algorithm, the original index will be kept
 * in the input, the idx value is 0,1,2,3,4,5...
 * content of array 'sort' will not be changed.
 ****************************************************************************/
{
    T pivot;
    int m;
    int i;
    if(low < high)
    {
		Swap(&idx[low], &idx[(high+low)/2]);
	    pivot = sort[idx[low]];
		m = low;
	    for (i = low + 1; i <= high; i++)
		{
			if(((order == ASCENDING) && (sort[idx[i]] <  pivot))
                    || (order == DESCENDING && (sort[idx[i]] > pivot))) 
			{
				m++;
                Swap(&idx[m], &idx[i]);
			}
	   }
	   Swap(&idx[low], &idx[m]);
	   QuickSort_index(idx, sort, low, m - 1, order);
	   QuickSort_index(idx, sort, m + 1, high, order);
	}
}/*}}}*/
template <class T, class InputIterator> void Set2Array(InputIterator __First, InputIterator __Last, T *array)/*{{{*/
// convert array stored in set to a C-like normal array
// Example: Set2Array(intset.begin(), intset.end(), array);
{
    InputIterator is;
    for(is = __First ; is != __Last; is++)
        *array++ = *is;
}
/*}}}*/

template <class T> int locate_range(T key, T *array, int lo, int hi)/*{{{*/
/*****************************************************************************
 * locate the range of key in the array which is in accending order,  
 * e.g., array[5] = 0.8, array[6] = 0.9, for key= 0.85, retern 5      
 * if key < array[low], return low -1;                                
 * if key > array[high] , return high;                                
 * Example: index = locate_range(key, array, 0, n-1);
 ****************************************************************************/
{
    int mid;
    /* Repeat while there are elements in range */
    if(key < array[lo]) return (lo-1);
    if(key >= array[hi]) return hi;
    while (lo < hi)  
    {
        mid = (lo + hi) / 2; /* Compute midpoint */
        if (array[mid] <= key && array[mid+1] > key)            /* found element */
        { return(mid); } 
        else if (array[mid] > key)      /* target is in first half */
        { hi = mid ; }  
        else /* array[mid] < target , target is in second half */
        { lo = mid ; }
    }
    return(-1); /* Nothing left in range; failure */
}/*}}}*/

#endif /*HAS_MYTEMPLATE_H*/
