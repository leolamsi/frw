/***************************************************************************
                          FileIO.h  -  description
                             -------------------
    begin                : Mon Aug 5 2002
    copyright            : (C) 2002 by Ricky
    email                : ricky@y-trace.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#define MESSAGE(x) cerr << x << endl;
//#define DIE(x) {MESSAGE(x); exit(1); }
#define SYSTEM(s) {MESSAGE(s); system(s);};

#define WARN(s) MESSAGE("Warning at __FILE__(" << __LINE__ << "): " s); 
#define DIE(s) {MESSAGE("Error at " __FILE__ "(" << __LINE__ << "): " s); exit(1); }
#define CHKPTR(s) {if (s==NULL){DIE(); }};
#define CHK(s) {if (s){DIE(s); }};

#define ENTER_MESSAGE(s)  cout << "Entering " #s << endl;  
#define EXIT_MESSAGE(s)   cout << "Exiting " #s << endl;

#define DUMP(var) cout << #var "=" << var << "\t";
#define LDUMP(var) cout << #var "=" << var << endl;
#define DDUMP(var) cout << #var "=" << (double)var << "\t";
#define DUMPARRAY(arr,arr_size,ncol) {cout << #arr << endl;		\
    for (int dumpi=0; dumpi<arr_size; ++dumpi){ cout << setw(10) << (arr)[dumpi] << ' '; if ((dumpi&(ncol-1))==(ncol-1)) cout << endl; } \
    cout << endl;}
  // print one array of size arr_size in ncol-column format

#define DUMPARRAY2D(arr,sizex, sizey) {cout << #arr << endl;		\
    for (int i=0; i<sizey; ++i){    for (int j=0; j<sizex; ++j) cout << setw(10) << (arr)[i][j] << ' ';  cout << endl; } \
    cout << endl;}

#define LOGPLOT(array,size) { \
 system("rm plot.tmp"); \
 FileIO::AppendArrs("plot.tmp", size, 1, array); \
 system("gp l plot.tmp"); \
}
#define PLOT(array,size) { \
 system("rm plot.tmp"); \
 FileIO::AppendArrs("plot.tmp", size, 1, array); \
 system("gp plot.tmp"); \
}

using namespace std;

#ifndef FILEIO_H
#define FILEIO_H

#ifdef _MSC_VER
	#pragma warning( disable: 4786 ) // identifier was truncated to '255' characters in the debug information
#endif	// _MSC_VER

#include <fstream>
#include <stdarg.h>
#include <assert.h>


// Macros for multiple arguments (own non-potable version)
#define my_va_start(VA_LIST,LAST_NAMED_ARG) \
	((VA_LIST) = ( (char *)(&(LAST_NAMED_ARG)) + sizeof(LAST_NAMED_ARG) ))
#define my_va_arg(VA_LIST, SOME_TYPE) \
	((VA_LIST) = (char *)(VA_LIST) + sizeof(SOME_TYPE), \
	*(SOME_TYPE*)((char *)(VA_LIST) - sizeof(SOME_TYPE)))
#define my_va_end(VA_LIST)

#define varptr2list(n_array,varptr,arrlist)	\
  va_list lptr; \
  va_start(lptr, varptr); \
  T** arrlist= new T*[n_array]; \
  arrlist[0]=varptr; \
  for (int j=1; j< n_array; ++j) arrlist[j]=va_arg(lptr, T*); \
  va_end(lptr);

#define var2list(n_array,var,arrlist)	\
  va_list lptr; \
  va_start(lptr, var); \
  T** arrlist= new T*[n_array]; \
  T*  valuelist= new T[n_array]; \
  arrlist[0]=&valuelist[0]; valuelist[0]=var;                            \
  for (int j=1; j< n_array; ++j){ arrlist[j]=&valuelist[j]; valuelist[j]=va_arg(lptr, T);} \
  va_end(lptr);

namespace FileIO {

//std::ostream::openmode appmode = std::ostream::app;
  std::ostream::openmode appmode = (std::ostream::app | std::ostream::binary) ; // binary forces \n to be 0a for all OS

static void SetPrecision(int prec);

static void setnblank(unsigned int blank);

template <class T>
void SaveArrs(	const char* filename,		// The file name
				std::ostream::openmode mode,	// File opening option, [ ofstream::out | app ] etc...
				int arr_size,				// Size of the array
				int n_array,				// Number of array
				T* arrptr, ... );			// Variable number of array to be saved

namespace PRIVATE_NS {	// Private namespace, functions for internal use are put here
template <class T>
void SaveArrs_(	const char* filename,		// The file name
				std::ostream::openmode mode,// File opening option, [ ofstream::out | app ] etc...
				int arr_size,				// Size of the array
				int n_array,				// Number of array
				va_list prev_list,			// The start pointer of the previous function's argument
				T* arrptr, ... );			// Variable number of array to be saved
};

template <class T>
int ReadArrs(	const char* filename,	// The file name
				int n_array,			// Size of the array
				T* arrptr, ... );		// Variable number of array to be saved

template <class T>
void ReadArrs(	const char* filename,	// The file name
				int arr_size,			// Size of the array
				va_list prev_list,		// The start pointer of the previous function's argument
				T* arrptr, ... );		// Variable number of array to be saved

template <class T>
void AppendArrs(	const char* filename,	// The file name
					int arr_size,			// Size of the array
					int n_array,			// Number of array
					T* arrptr, ... );		// Variable number of array to be saved

template <class T>							// Warning: This function is compiler dependent
void AppendArrs(	const char* filename,	// The file name
					int arr_size,			// Size of the array
					T* arrptr, ... );		// Variable number of array to be saved

template <class T>
void AppendVars(	const char* filename,	// The file name
					int n_var,				// Number of variables
					T var, ... );		    // Variable number of variables to be saved

template <class T>
void AppendArr2d(	const char* filename,	// The file name
				 	int maxi,				// Dimension of the array
				 	int maxj,				//
				 	T* p	);				// Array's pointer

template <class T>
void ReadArr2d(	const char* filename,	// The file name
				int maxi,				// Get dimension of the array
				int maxj,				//
				T* p	);				// Array's pointer

class NumIfstream : public std::ifstream{
 public:

  bool eoln(){
    char chkchr;
    while (peek()==' ' || peek()=='\t'){
      get(chkchr); // skip blanks
    }
    bool eolnflg=(peek()=='\n' || peek()=='\r' || eof())?1:0;
    //std::cout << "eoln=" << eolnflg << " peek=" << peek() << line_separator;
    return eolnflg;
  }

  void skipline(){
  char chkchr=' ';
  do{
    get(chkchr);
  }while(chkchr!='\n' && (!eof()));
  }

};


static std::ofstream OutFile;
//static std::ifstream InFile;
static NumIfstream InFile;

static char col_separator = '\t';
static char line_separator = '\n';
static int nblank = 1;
static int i0 = 0;

//void skipline();
static int GetArrCount(const char* filename);

/****************************************************
 * Program implementation
 */


/* Save a numbers of one dimension arrays
 * into a file 
 */

// Set the output precision
static void SetPrecision(int prec)
{
	OutFile.precision(prec);
}

static void setnblank(unsigned int blank)
{
	FileIO::nblank = blank;
}

static void seti0(unsigned int i0)
{
	FileIO::i0 = i0;
}

template <class T>
void SaveArrs(	const char* filename,		// The file name
				std::ostream::openmode mode,// File opening option, [ ofstream::out | app ] etc...
				int arr_size,				// Size of the array
				int n_array,				// Number of array
				T* arrptr, ... )			// Variable number of array to be saved
{
	va_list lptr;
	va_start(lptr, arrptr);
	PRIVATE_NS::SaveArrs_(filename, mode, arr_size, n_array, lptr, arrptr);
	va_end(lptr);
}


namespace PRIVATE_NS {
template <class T>
void SaveArrs_(	const char* filename,		// The file name
				std::ostream::openmode mode,// File opening option, [ ofstream::out | app ] etc...
				int arr_size,				// Size of the array
				int n_array,				// Number of array
				va_list prev_list,			// The start pointer of the previous function's argument
				T* arrptr, ... )			// Variable number of array to be saved
{
  assert(arr_size > 0);
  assert(n_array > 0);
  assert(arrptr);
  
  int i, j;						// Loop variable
  //va_list lptr;
  T* p;
  OutFile.open(filename, mode);
  
  //va_start(lptr, arrptr);
  va_start(prev_list, arrptr);
  for(i=0; i<arr_size; i++) {
    j = n_array;
    //Modified 05-12-06 to avoid: incompatible types in assignment of `__va_list_tag*' to `...
    //if(prev_list == NULL)
    //va_start(lptr, arrptr);	// Initialize the list of arguments
    //else
    //lptr = prev_list;		// Initialize the list from previous function
    p = arrptr;					// Initialize p with the first argument
    while(j-- > 0) {
      OutFile << p[i];
      if(j > 0)				// If not the last column, print out the col_separator
	OutFile << col_separator;
      p = va_arg(prev_list, T*);	// Get the next argument
    }
    OutFile << line_separator;
  }
  
  for(i = 0; i < FileIO::nblank; ++i)
    OutFile << line_separator;
  
  OutFile.close();
  va_end(prev_list);				// End of all arguments
}

template <class T>
void SaveArrs_List(	const char* filename,		// The file name
				std::ostream::openmode mode,// File opening option, [ ofstream::out | app ] etc...
				int arr_size,				// Size of the array
				int n_array,				// Number of array
			T** arrlistptr)			// Variable number of array to be saved
// save arrays based on array list
// 07-02-01: previous version cannot handle 2 or more arrays in Suse 10.2
{
  assert(arr_size > 0);
  assert(n_array > 0);
  assert(arrlistptr);
  
  OutFile.open(filename, mode);
  
  for(int i=i0; i<arr_size+i0; i++) {
    for (int j=0; j< n_array; ++j){
      OutFile << arrlistptr[j][i];
      if(j < n_array-1)	// If not the last column, print col_separator
	OutFile << col_separator;
    }
    OutFile << line_separator;
  }
  
  for(int i = 0; i < FileIO::nblank; ++i)
    OutFile << line_separator;
  
  OutFile.close();
}


};


template <class T>
void ReadArrs_(	const char* filename,	// The file name
				int arr_size,			// Size of the array
				T* arrptr, ... )		// Variable number of array to be saved
{
	va_list lptr;
	va_start(lptr, arrptr);
	ReadArrs(filename, arr_size, lptr, arrptr);
	va_end(lptr);
}

#if 0
template <class T> 
int ReadArrs(	const char* filename,	// The file name
		int n_array,
		T* arrptr, ... )		// Variable number of array to be read
// added 070605 
{
	int line,column;
	assert(arrptr);
	assert(n_array > 0);

	varptr2list(n_array,arrptr,arrlist);
	InFile.open(filename);
	line=0; column=99;
	while((!InFile.eof()) && column>0) {
	  column=0;
	  while((!InFile.eoln()) && (column<n_array)) {
	    InFile >> arrlist[column][line];
	    //cout << arrlist[0][line] << ' ' << line << line_separator;
	    ++column;
	  }
	  //InFile.get(); // skip newline
	  skipline();
	  if (column>0) ++line;
	}
	InFile.close();
	return line;
}

template <class T>
int ReadArrs(	const char* filename,	// The file name
		T* arrptr, ... )		// Variable number of array to be read
// obsolete
{
	int line,column;
	va_list lptr ;
	T* p;
	int arr_size=0;
	assert(arrptr);
	InFile.open(filename);
	va_start(lptr, arrptr);
	line=0; column=99;
	while((!InFile.eof()) && column>0) {
	  va_start(lptr, arrptr);	// Initialize the list of arguments
	  p = arrptr;		// Initialize p with the first argument
	  column=0;
	  while((!InFile.eoln()) && (p!=0)) {
	    ++column;
	    InFile >> p[line];
	    p = va_arg(lptr, T*);	// Get the next argument
	  }
	  InFile.get(); // skip newline
	  if (column>0) ++line;
	}
	InFile.close();
	va_end(lptr);				// End of all arguments
	return line;
}
#endif 

template <class T>
int ReadLastLine(	const char* filename,	// The file name
		T* arrptr, ... )		// Variable number of array to be read
{
	int line,column;
	va_list lptr ;
	T* p;
	int arr_size=0;
	assert(arrptr);
	InFile.open(filename);
	va_start(lptr, arrptr);
	line=0; column=99;
	while((!InFile.eof()) && column>0) {
	  va_start(lptr, arrptr);	// Initialize the list of arguments
	  p = arrptr;		// Initialize p with the first argument
	  column=0;
	  while((!InFile.eoln()) && (p!=0)) {
	    ++column;
	    InFile >> p[0];
	    p = va_arg(lptr, T*);	// Get the next argument
	  }
	  InFile.get(); // skip newline
	  if (column>0) ++line;
	}
	InFile.close();
	va_end(lptr);				// End of all arguments
	return line;
}

template <class T>
void ReadArrs___(	const char* filename,	// The file name
		int arr_size,			// Size of the array
		//va_list prev_list,		// The start pointer of the previous function's argument
		T* arrptr, ... )		// Variable number of array to be saved
{
	assert(arr_size > 0);
	assert(arrptr);

	int n_array = GetArrCount(filename);	// Get the column count
	assert(n_array > 0);			// Should have some data inside the file
	int i, j;						// Loop variable
	va_list lptr ;
	T* p;
	InFile.open(filename);

	va_start(lptr, arrptr);
	for(i=0; i<arr_size; i++) {

		j = n_array;
		//		if(prev_list == NULL)
			va_start(lptr, arrptr);	// Initialize the list of arguments
			//else
			//lptr = prev_list;		// Initialize the list from previous function
		p = arrptr;					// Initialize p with the first argument
		while(j-- > 0) {
			InFile >> p[i];
			p = va_arg(lptr, T*);	// Get the next argument
		}
	}

	InFile.close();
	va_end(lptr);				// End of all arguments
}


template <class T>
void AppendArrs( const char* filename,	// The file name
		 int arr_size,			// Size of the array
		 int n_array,			// Number of array
		 T* arrptr, ... )		// Variable number of array to be saved
{
  varptr2list(n_array,arrptr,arrlist);
  PRIVATE_NS::SaveArrs_List<T>(filename, appmode, arr_size, n_array, arrlist);
}

template <class T>
void WriteArrs( const char* filename,	// The file name
		 int arr_size,			// Size of the array
		 int n_array,			// Number of array
		 T* arrptr, ... )		// Variable number of array to be saved
{
  varptr2list(n_array,arrptr,arrlist);
  PRIVATE_NS::SaveArrs_List<T>(filename, std::ofstream::out, arr_size, n_array, arrlist);
}

template <class T>
void AppendArrs(	const char* filename,	// The file name
					int arr_size,			// Size of the array
					T* arrptr, ... )		// Variable number of array to be saved
  // counting variables does not work for Suse 10.2
{
	va_list lptr ;
	T* p;	T* p2;
	va_start(lptr, arrptr);
	p = p2 = arrptr;
	int count = 0;				// Count for number of arguments

	while(p) { 
		count ++;
		p = va_arg(lptr, T*);	// Get the next argument
		if(p == p2)		// Check for arguments termination, either NULL or same as the last
			p = NULL;	//
		else			//
			p2 = p;		//
	}

	va_end(lptr);
	va_start(lptr, arrptr);
	PRIVATE_NS::SaveArrs_<T>(filename, appmode, arr_size, count, lptr, arrptr);
	va_end(lptr);
}


template <class T>
void AppendVars(	const char* filename,	// The file name
					int n_var,				// Number of variables
					T var, ... )		    // Variable number of variables to be saved
{
  int nblanksave=nblank; setnblank(0);
  int i0save=i0; seti0(0);
  var2list(n_var,var,arrlist);
  PRIVATE_NS::SaveArrs_List<T>(filename, appmode, 1, n_var, arrlist);
  setnblank(nblanksave);
  seti0(i0save);
}


template <class T>
void AppendArr2d(	const char* filename,	// The file name
				 	int maxi,				// Dimension of the array
				 	int maxj,				//
				 	T* p	)				// Array's pointer
{
	int i;
	for(i=0; i<maxi; i++) {
		SaveArrs(filename, appmode, maxj, 1, p[i]);
	}
	OutFile.open(filename, appmode);
	OutFile << line_separator;
	OutFile.close();
}


template <class T>  // old version
void ReadArr2d(	const char* filename,	// The file name
				int maxi,				// Get dimension of the array
				int maxj,				//
				T* p	)				// Array's pointer
{
	assert(GetArrCount(filename) == 1);		// The file should hold only one column
	int i, j;
	InFile.open(filename);

	for(i=0; i<maxi; i++) {
		for(j=i0; j<i0+maxj; j++) {
			InFile >> p[i][j];
		}
	}

	InFile.close();
}

//////////////////////////////////////////////////////////////////////////

template <class T>
void AppendFunc(const char* filename,		// The file name
		   T (*func)(T),
		   T xa, T xb,
		   int npt)				// Size of the array
{
  T *X = new T[npt];
  T *Y = new T[npt];
  T dx=(xb-xa)/(npt-1);
  for (int i=0; i<npt; ++i){
    T x=xa+i*dx;
    X[i]=x;
    Y[i]=(*func)(x);
  }
  int i0save=i0; seti0(0);
  AppendArrs(filename, npt, 2, X, Y);
  seti0(i0save);
}

//////////////////////////////////////////////////////////////////////////

// Get the number of array stored int the file (column count)

static int GetArrCount(const char* filename)
{
    InFile.open(filename, std::fstream::binary);
	char chkchr=' ';
	bool LastIsSep = true;		// Is the last character the separator?
	int colcount = 0;

	// Check for count of non-continuous non-col_separator in the first line of the file
	while(!InFile.eof() && chkchr != line_separator) {
		InFile.get(chkchr);
		if(chkchr != col_separator && chkchr != line_separator && LastIsSep == true) {
			colcount++;
			LastIsSep = false;
		}
		else if(chkchr == col_separator) {
			LastIsSep = true;
		}
	}; 
	InFile.close();
    LDUMP(colcount);

	return colcount;
}

/////////////////////////////////////////////////////////////////////////
// 071218 new read array procedure
/////////////////////////////////////////////////////////////////////////

class ArrIfstream : public std::ifstream{

 public:

  void skipline(){
    char chkchr=' ';
    while((!eof()) && (chkchr!='\n')){
      get(chkchr);
    }
  }

  bool eoln(){
    char chkchr;
    while (peek()==' ' || peek()=='\t'){
      get(chkchr); // skip blanks
    }
    bool eolnflg=(peek()=='\n' || peek()=='\r' || eof())?1:0;
    //std::cout << "eoln=" << eolnflg << " peek=" << peek() << line_separator;
    return eolnflg;
  }

  template <class T> 
    int readarrs( int n_array, T* arrptr, ... )
    {
      int line,column;
      assert(arrptr);
      assert(n_array > 0);
      
      varptr2list(n_array,arrptr,arrlist);
      line=0; column=99;
      while((!eof()) && column>0) {
        column=0;
        while((!eoln()) && (column<n_array)) {
          *this >> arrlist[column][line];
          //DUMP(line); LDUMP(arrlist[column][line]);
          ++column;
        }
        skipline();
        if (column>0) ++line;
      }
      return line;
    }

  template <class T>  // T is typically double 
    int readarr(T* arrptr)
    {
      return readarrs(1, arrptr);
    }

  template <class T> // T is typically double []
    bool readarr2d( int* maxiptr,	// Get dimension of the array
		    int* maxjptr,	//
		    T* p	)	// Array's pointer
    {
      //assert(GetArrCount(filename) == 1);		// The file should hold only one column
      int i=0, j=0, maxj=0, linemax=0;
      bool rectflg=true;
      while(!eof()){
        int line=readarr(p[i]);
        //LDUMP(line);
        if (line>0){
          ++i;
          if (linemax==0)
            linemax=line;
          else if (line!=linemax){
            DIE("Array not in rectangular format"); 
            rectflg=false;
          }
        }
      }		
      *maxiptr=i;
      *maxjptr=linemax;
      return rectflg;
    }
};

template <class T>
  int ReadArrs(	const char* fname,	// The file name
				int n_array,			// Size of the array
				T* arrptr, ... )		// Variable number of array to be saved
  {
    ArrIfstream ifile;
    ifile.open(fname);
    int line=ifile.readarrs(n_array, arrptr);
    ifile.close();
    return line;
  }

template <class T>
void ReadArr2d(	const char* fname,	// The file name
				int* maxiptr,				// Get dimension of the array
				int* maxjptr,				//
				T* p	)				// Array's pointer
{
  ArrIfstream ifile;
  ifile.open(fname);
  if ( !ifile ) DIE("Cannot Open: " << fname);

  ifile.readarr2d(maxiptr, maxjptr, p);
  ifile.close();
}



} // End namespace FileIO

#endif // FILEIO_H
