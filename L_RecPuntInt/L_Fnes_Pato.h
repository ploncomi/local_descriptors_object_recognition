#ifndef __L_FNES_PATO_H_
#define __L_FNES_PATO_H_

// Use dynamic_bitset instead of vector<bool>

// Includes interfaces for ATL, OpenCV, FLTK, SVS by using some #defines
// Ansi C++ 99, Ansi C 89
// Author: Patricio Loncomilla (ploncomi@gmail.com)

// Code properties:
// * All objects manage their own allocated memory, all destroy their allocated memory when destructed
// * swaps between objects are used when needed


// Copyright (c) 2011, Patricio Loncomilla
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. All advertising materials mentioning features or use of this software
//    must display the following acknowledgement:
//    This product includes software developed by Patricio Loncomilla.
// 4. Neither the name of Patricio Loncomilla nor the
//    names of its contributors may be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY Patricio Loncomilla ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Patricio Loncomilla BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS






// Operator = for Borland compiler for POD objects
#define __OP_ASSIGN(X) X &operator=(const X &other) {memcpy(this, &other, sizeof(X)); return *this;}

// Redirectable assert
#define throw_L_ArgException_if(X,S)  ( (X) ? (throw_L_ArgException_fn_NULL(S)) : 0 )  // build un null pointer exception
//#define throw_L_ArgException_if(X,S)  ( (X) ? (throw_L_ArgException_fn(S)) : 0 )  // Lanza una excepcion y un mensaje, se puede separar de la siguiente instruccion con una coma
//#define throw_L_ArgException_if(X,S)  ( (X) ? (printing_L_ArgException_fn(S)) : 0 )  // Lanza una excepcion y un mensaje, se puede separar de la siguiente instruccion con una coma
//#define throw_L_ArgException_if(X,S)  ( (X) ? (throw_L_ArgException_fn(NULL)) : 0 )  // Es piola, se puede atrapar con try-catch sin que se note
//#define throw_L_ArgException_if(X,S) (0)  // No hace nada, rapido y mortal

template<bool> class L_StaticAssertClass {public: L_StaticAssertClass(...) {}};
template<> class L_StaticAssertClass<false>{};

#define L_StaticAssert(test, errormsg)          do {struct ERROR_##errormsg {}; typedef L_StaticAssertClass< (test) != 0 > tmplimpl; tmplimpl aTemp = tmplimpl(ERROR_##errormsg());sizeof(aTemp);} while (0)


#define DEBUG_BREAK_IF(X,S)  ( (X) ? (throw_L_ArgException_fn_NULL(S)) : 0 )

//#define USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED
//#define DEFINE_NONEFFICIENT_MATRIX_OPERATORS
//#define USE_NONEFFICIENT_MATRIX_OPERATORS


// FLAGS PARA DEBUG IMPORTANTES

#define L_ARRAY_ACCESS_DEBUG
//#define L_STRING_ACCESS_DEBUG
#define L_MATRIX_ACCESS_DEBUG
//#define L_IMAGE_ACCESS_DEBUG

//#define L_BRACKET_LEVEL_DEBUG   // Slow, associated to L_ptrM, and L_ptrME, only for debugging very strange cases...



#define L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE 1      // Tamano extra NULL para las matrices con doble puntero
#define L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_VALUE (NULL)  // Un lugar de la memoria que no sea valido, podria ser NULL...



#define L_ARREGLO_ALLOC_SIZE_FACTOR 10  // 10 .. Cantidad en que se aumentan los arreglos por defecto

#if defined(_MSC_VER) || defined(__GNUC__) // En visual studio tiene un name raro...
#define L_restrict __restrict
#else
#define L_restrict
#endif

// Memory hints for Visual Studio:
// 0xCDCDCDCD	Allocated in heap, but not initialized
// 0xDDDDDDDD	Released heap memory.
// 0xFDFDFDFD	"NoMansLand" fences automatically placed at boundary of heap memory. Should never be overwritten. If you do overwrite one, you're probably walking off the end of an array.
// 0xCCCCCCCC	Allocated on stack, but not initialized

#define nothrows() throw()

//#define L_DEFINE_PRIME_NUMBERS   // Non included for zero overhead
#ifdef L_DEFINE_PRIME_NUMBERS
size_t L_PRIMOS[100];
#endif

//#define L_USE_SMOOTH_SIGMA_STORING_TREE  // This can be included for storing sigmas used in smoothings
//#define L_BASIC_DEBUG  // Ensure greater-than-zero allocations


#ifndef __TURBOC__ // Support for Turbo C :)
#include <cstdio> // ANSI C++ 98
#include <cmath> // ANSI C++ 98
#include <cstdlib> // ANSI C++ 98
#include <cstdarg> // ANSI C++ 98
#include <cstring> // ANSI C++ 98
#include <ctime> // ANSI C++ 98
#else
#include <stdio.h> // ANSI C++ 98
#include <math.h> // ANSI C++ 98
#include <stdlib.h> // ANSI C++ 98
#include <stdarg.h> // ANSI C++ 98
#include <string.h> // ANSI C++ 98
#include <time.h> // ANSI C++ 98
#endif

#include <exception> // ANSI C++ 98
#include <new> // ANSI C++ 98
#include <vector> // ANSI C++ 98
#include <list> // ANSI C++ 98
#include <algorithm>  // ANSI C++ 98
#include <iostream>  // ANSI C++ 98
#include <fstream>   // ANSI C++ 98

/////////////////////
/////////////////////   #ifdefs para bibliotecas
/////////////////////


#ifdef __COMPAT_EIGEN__  // Optional interface to Eigen (faster matrix operations)
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/LU>
#endif

#ifdef __COMPAT_ATLIMAGE__ // Optional interface for ATL images (Windows)
#include <atlimage.h>
typedef ATL::CImage ATL_CImage; // For enabling ATL: select Project -> X Properties -> General -> Use Of ATL
#endif

#ifdef __COMPAT_IPLIMAGE__ // Optional interface for OpenCV
#include "cv.h"
#include "highgui.h"
#endif

#ifdef __COMPAT_CVMAT__
#include <opencv2/opencv.hpp>
#endif


#ifdef __COMPAT_TBITMAP__ // Optional interface for Borland
#include <graphics.hpp>
#include <jpeg.hpp>
#endif

#ifdef __COMPAT_SVS__ // Opcional interface for SVS (stereo camera Videre)
#include <csignal> // For managing internal segmentation fault in SVS library
#include <csetjmp> // For managing internal segmentation fault in SVS library
#include "svsclass.h"
#endif

#ifdef __COMPAT_OPENNI__  // Opcional interface for OpenNI (Kinect library)
#include <XnCppWrapper.h>
#include <XnTypes.h>
#endif // __COMPAT_OPENNI__

#ifdef __COMPAT_KINECTSDK__
#include <MSR_NuiApi.h>
#endif // __COMPAT_KINECTSDK__

#ifdef __COMPAT_FLTK__ // Opcional interface for FLTK (portable GUI / window system)
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/fl_draw.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <GL/glu.h>
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_UBLAS__
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp>
#endif

#ifdef __COMPAT_OPENMP__
#include "omp.h"
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef SQRT2
#define SQRT2 1.41421356237309504880168872421
#endif

#ifndef SQRT3
#define SQRT3 1.732050807568877293527446341505872366943
#endif

#if defined(L_PUSH_EXECUTING_FN) || defined(L_POP_EXECUTING_FN) || defined(L_DEBUGPRINTF)
	#error Confilctos de #defines
#endif

#define L_PUSH_EXECUTING_FN(x) //L_DebugPushEnterFuncion(x) // Activar las 2 funciones si se desea ver el stack
#define L_POP_EXECUTING_FN(x) //L_DebugPopExitFuncion(x)
#define L_DEBUGPRINTF(x,y) printf(x,y) //Para hacer debug "temporal", fácil de encontrar en el código

#if defined(L_ABS) || defined(L_SETABS) || defined(L_SETMIN) || defined(L_SETMAX) || defined(L_TIME) || defined(L_MIN) || defined(L_MAX) || defined(L_FLOOR) || defined(L_ROUND) || defined(L_XOR)
#error #define conflicts
#endif

#define L_ABS(x) ( ((x)>=0)?(x):-(x) )
#define L_SETABS(x,y) {(x) = (y); if ((x)<0) (x)=-(x);} // Preferible a L_ABS(x)
#define L_SETMIN(x,y,z) {if ((y)<(z)) (x)=(y); else (x)=(z);}
#define L_SETMAX(x,y,z) {if ((y)>(z)) (x)=(y); else (x)=(z);}
#define L_TIME() (clock()/(double)CLOCKS_PER_SEC)
#define L_MIN(x,y) (  ((x)<(y))?(x):(y)  )
#define L_MAX(x,y) (  ((x)>(y))?(x):(y)  )
#define L_FLOOR(x) (int)floor(x)
#define L_FLOOR1000(x) ((int)((x)+1000)-1000)
#define L_FLOOR10000(x) ((int)((x)+10000)-10000)
#define L_ROUND(x) L_FLOOR((x)+0.5)
#define L_ROUND1000(x) L_FLOOR1000((x)+0.5)
#define L_ROUND10000(x) L_FLOOR10000((x)+0.5)
#define L_XOR(a,b) ((a)^(b)) //( ((a)|(b))&(~(a)|~(b)) )
#define L_RANDOMIZE() (srand((unsigned int)time(NULL)))
#define L_RANDOM(a,b) ((a)+(((b)-(a))*1.0*rand())/(RAND_MAX+1.0))
template<class T> inline void L_randomPermutation(T *arr, int n) {int i, u; T t; for (i=0; i<n; i++) {u = rand()%(n-i); t=arr[u]; arr[u] = arr[i]; arr[i] = t;}}
inline double L_difang(double xMinusY)  {while (xMinusY>M_PI) xMinusY-=2*M_PI; while (xMinusY<-M_PI) xMinusY+=2*M_PI; return xMinusY;}
inline double L_promang(double x, double y) {return x+0.5*L_difang(y-x);}



#if defined __GNUC__
inline void *L_aligned_malloc_128(size_t size) {void *ptr; posix_memalign(&ptr, 16, size); return ptr;}
inline void L_aligned_free_128(void *mem) {free(mem);}
#elif defined _WIN32 
inline void *L_aligned_malloc_128(size_t size) {return _aligned_malloc(size * sizeof(float), 16);}
inline void L_aligned_free_128(void *mem) {_aligned_free(mem);}
#else
inline void *L_aligned_malloc_128(size_t size) {return malloc(size);}
inline void L_aligned_free_128(void *mem) {free(mem);}
#endif





inline int L_sgn(double x) {if (x>0) return 1; return -1;}

#if defined(L_LOG_2) || defined(L_log_base_2)
#error #define conflicts
#endif

#define L_LOG_2 (0.6931471806) // =log(2.0)
#define L_log_base_2(X) ( log(X) / L_LOG_2 )

#if defined(L_INTERP_MAX_PARABOL_ARG) || defined(L_INTERP_MAX_PARABOL_VAL)
#error #define conflicts
#endif

#define L_INTERP_MAX_PARABOL(a,b,c) (0.5 * ((a) - (c)) / ((a) - 2.0 * (b) + (c)))
#define L_INTERP_MAX_PARABOL_VAL(a,b,c,resArg,tmp1,tmp2,resVal) if ((tmp1=0.5*(a+c)-b)!=0) {tmp2=0.5*(c-a); resArg=-tmp2/(2.0*tmp1); resVal=(tmp1*arg+tmp2)*arg+b;} else {resArg=(a>c)?-1:1; resVal=(a>c)?a:c;}

#if defined _WIN32 || defined _WIN32_ || defined __WIN32__
#define FOLDER_SEPARATION_CHAR '/'
#define L_DEL_SYSTEM_CMD "del"
#else
#define FOLDER_SEPARATION_CHAR '/'
#define L_DEL_SYSTEM_CMD "rm"
#endif

void L_hard_shutdown_fn(const char *file, long line, const char *x, ...);
bool L_printError(const char *formato, ...); //! Función para lanzar excepción grave en el sistema

#if defined(L_hard_shutdown) || defined(L_hard_shutdown2) || defined(L_hard_shutdown3) || defined(L_hard_shutdown4) || defined(L_hard_shutdown5)
#error #define conflicts
#endif

#define L_hard_shutdown(X) L_hard_shutdown_fn(__FILE__, __LINE__, X)
#define L_hard_shutdown2(X,Y1) L_hard_shutdown_fn(__FILE__, __LINE__, X,Y1)
#define L_hard_shutdown3(X,Y1,Y2) L_hard_shutdown_fn(__FILE__, __LINE__, X,Y1,Y2)
#define L_hard_shutdown4(X,Y1,Y2,Y3) L_hard_shutdown_fn(__FILE__, __LINE__, X,Y1,Y2,Y3)
#define L_hard_shutdown5(X,Y1,Y2,Y3,Y4) L_hard_shutdown_fn(__FILE__, __LINE__, X,Y1,Y2,Y3,Y4)

#define L_REVHEAP _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF | _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG )) // Si se addFrom esto en Win se revisa el heap

/////// Tipos basicos
//
typedef unsigned char L_uchar; //! tipo basico, hay que revisar que los tamanos sean correctos
typedef char L_int8;
typedef short L_int16;
typedef long L_int32;
typedef unsigned char L_uint8;
typedef unsigned short L_uint16;
typedef unsigned long L_uint32;
union L_LongChars {long lo; char ch[4];};  // Se usa solo para guardar los tamanos

enum L_PrintingLanguage {L_espanol, L_english};


//////// Storing classes
class L_HuffmanTree; // For compressing data
class L_SignalDouble; //! 1D signal, double type
class L_ImageGrayDouble;
class L_ImageGrayUchar;
class L_ImageRGBDouble; // interleaved
class L_ImageRGBUchar; // interleaved
class L_ImageBGRAuchar_row; // interleaved
class L_ImageRGBAuchar_row; // interleaved
class L_ImagenRGBuchar_row; // interleaved
class L_Matrix; // matrix "in heap"

//// POD classes
class L_ComplexDouble;
class L_Quaternion;
class L_CoordsCart3D;
class L_PanTiltRoll3D;
class L_Pose3D;
class L_Pose3D_cuat;
class L_Rot3DMatrix; // 3x3 rotation matrix "in stack"
class L_HomogeneousMatrix;   // 4d rotation & translation matrix "in stack"
class L_EssentialMatrix;  // q'Eq = 0 for correspondences between two images

/////// Global functions and defined
// Gaussianas
inline int L_GaussianHalfWidth(double s) {return (int)(3*s+0.4);} //(int)(3*s+0.4)  //! Half of width for a Gaussian
inline int L_GaussianWidth(double s) {return 2*L_GaussianHalfWidth(s)+1;} //! Width for a Gaussian

// fourier
#define L_FFT1D(dir,fouR,fouI,log2n) {double *x96 = fouR; double *y96 = fouI; long m96 = log2n; {long n96, i96,i1_96,j96,k96,i2_96,l_96,l1_96,l2_96; double c1_96,c2_96,tx_96,ty_96,t1_96,t2_96,u1_96,u2_96,z_96;n96 = 1; for (i96=0;i96<m96;i96++) n96 *= 2; i2_96 = n96 >> 1; j96 = 0; for (i96=0;i96<n96-1;i96++)\
	{ if (i96 < j96) { tx_96 = x96[i96]; ty_96 = y96[i96]; x96[i96] = x96[j96]; y96[i96] = y96[j96]; x96[j96] = tx_96; y96[j96] = ty_96; } k96 = i2_96; while (k96 <= j96) { j96 -= k96; k96 >>= 1; } j96 += k96; } c1_96 = -1.0; c2_96 = 0.0; l2_96 = 1; for (l_96=0;l_96<m96;l_96++) {l1_96 = l2_96; l2_96 <<= 1; u1_96 = 1.0; u2_96 = 0.0; for (j96=0;j96<l1_96;j96++) \
	{for (i96=j96;i96<n96;i96+=l2_96) { i1_96 = i96 + l1_96; t1_96 = u1_96 * x96[i1_96] - u2_96 * y96[i1_96]; t2_96 = u1_96 * y96[i1_96] + u2_96 * x96[i1_96]; x96[i1_96] = x96[i96] - t1_96; y96[i1_96] = y96[i96] - t2_96; x96[i96] += t1_96; y96[i96] += t2_96;}z_96 =  u1_96 * c1_96 - u2_96 * c2_96;u2_96 = u1_96 * c2_96 + u2_96 * c1_96; u1_96 = z_96;}\
	c2_96 = sqrt((1.0 - c1_96) / 2.0); if (dir == 1) c2_96 = -c2_96; c1_96 = sqrt((1.0 + c1_96) / 2.0);}if (dir == 1) {for (i96=0;i96<n96;i96++) { x96[i96] /= n96; y96[i96] /= n96;}}}}
#define L_FFT1D_2(dir,fouR,fouI,log2n,N2) {double **x97 = fouR; double **y97 = fouI; long m97 = log2n; long n2_97 = N2; {long n_97, i97,i1_97,j97,k97,i2_97,l_97,l1_97,l2_97; double c1_97,c2_97,tx_97,ty_97,t1_97,t2_97,u1_97,u2_97,z_97;n_97 = 1; for (i97=0;i97<m97;i97++) n_97 *= 2; i2_97 = n_97 >> 1; j97 = 0; for (i97=0;i97<n_97-1;i97++) \
	{ if (i97 < j97) { tx_97 = x97[i97][n2_97]; ty_97 = y97[i97][n2_97]; x97[i97][n2_97] = x97[j97][n2_97]; y97[i97][n2_97] = y97[j97][n2_97]; x97[j97][n2_97] = tx_97; y97[j97][n2_97] = ty_97; } k97 = i2_97; while (k97 <= j97) { j97 -= k97; k97 >>= 1; } j97 += k97; } c1_97 = -1.0; c2_97 = 0.0; l2_97 = 1; for (l_97=0;l_97<m97;l_97++) \
	{l1_97 = l2_97; l2_97 <<= 1; u1_97 = 1.0; u2_97 = 0.0; for (j97=0;j97<l1_97;j97++) {for (i97=j97;i97<n_97;i97+=l2_97) { i1_97 = i97 + l1_97; t1_97 = u1_97 * x97[i1_97][n2_97] - u2_97 * y97[i1_97][n2_97]; t2_97 = u1_97 * y97[i1_97][n2_97] + u2_97 * x97[i1_97][n2_97]; x97[i1_97][n2_97] = x97[i97][n2_97] - t1_97; y97[i1_97][n2_97] = y97[i97][n2_97] - t2_97;\
	x97[i97][n2_97] += t1_97; y97[i97][n2_97] += t2_97;}z_97 =  u1_97 * c1_97 - u2_97 * c2_97;u2_97 = u1_97 * c2_97 + u2_97 * c1_97; u1_97 = z_97;} c2_97 = sqrt((1.0 - c1_97) / 2.0); if (dir == 1) c2_97 = -c2_97; c1_97 = sqrt((1.0 + c1_97) / 2.0);}if (dir == 1) {for (i97=0;i97<n_97;i97++) { x97[i97][n2_97] /= n_97; y97[i97][n2_97] /= n_97;}}}}


#if defined(L_CrossProduct2D)
	#error #define conflicts
#endif
#define L_CrossProduct2D(xA, yA, xB, yB, tipo) (  (xA) * (tipo)(yB) - (xB) * (tipo)(yA)  )
#define L_AreaTri(xA, yA, xB, yB, xC, yC, tipo)  ( 0.5*L_CrossProduct2D((xB)-(xA),(yB)-(yA),(xC)-(xA),(yC)-(yA),tipo) )

// Combinatorials
double NR_gammln(double xx);
double NR_betai(double a, double b, double x);
double NR_betacf(double a, double b, double x);
double NR_gammp(double a, double x);
double NR_gammq(double a, double x);
void NR_gcf(double* gammcf, double a, double x, double* gln);
void NR_gser(double* gamser, double a, double x, double* gln);
double NR_ran2(long *idum); // Initial seed must be negative
void L_pruebaRand();
inline double L_binom(double n, double k, double p) {return NR_betai(k, n-k+1, p);}
// Distributions
double L_GaussianNoise();
inline double L_facchi2(double n) {return 1/( pow(2,n/2) * exp(NR_gammln(n/2)) );}
inline double L_chi2(double x, double n) {return pow(x,n/2-1)*exp(-x/2)*(x>0);} // El resultado es facchi2(n)*chi2(x,n)
inline double L_chi2Accumulated(double x, double n) {return NR_gammp(n/2,x/2);}
// Strings
void L_div_name_file2(const char *dirnameext, char *dir, char *name, char *ext);
// Verify sizes of types
void L_verifTypeSizes();
// Debug for function enter and exit
void L_DebugPushEnterFuncion(const char *s);
void L_DebugPopExitFuncion(const char *s);
// Language
inline L_PrintingLanguage L_compilationLanguage() {return L_espanol;}


#if defined(_MSC_VER) // En visual studio tiene un name raro...
#define L_snprintf sprintf_s
#else
#define L_snprintf snprintf  // Si no esta va a haber un error de compilacion, lo cual esta bien...
#endif

#define L_memcpy memcpy

char *L_gets(char *buf, int sizebuf);
#define gets dont_use_gets_as_it_is_deprecated(;  // use L_gets or gets_s()

class L_Int
{
public:
	static int cmp(const void *inta, const void *intb);
	static int cmpinv(const void *inta, const void *intb);
	static void print(int val);
};

class L_Double
{
public:
	static int cmp(const void *inta, const void *intb);
	static int cmpinv(const void *inta, const void *intb);
	static void print(double val);
};

class L_IndOrd
{
public:
	int i;
	double d;
	static int cmp(const void *indOrdA, const void *indOrdB);
	int operator < (const L_IndOrd &b) const {return d < b.d;}
	int operator <= (const L_IndOrd &b) const {return d <= b.d;}
	int operator > (const L_IndOrd &b) const {return d > b.d;}
	int operator >= (const L_IndOrd &b) const {return d >= b.d;}
};

class L_IndOrdi
{
public:
	int i;
	int val;
	static int cmp(const void *indOrdiA, const void *indOrdiB);
	int operator < (const L_IndOrdi &b) const {return val < b.val;}
	int operator <= (const L_IndOrdi &b) const {return val <= b.val;}
	int operator > (const L_IndOrdi &b) const {return val > b.val;}
	int operator >= (const L_IndOrdi &b) const {return val >= b.val;}
};

class L_IndOrdL
{
public:
	int i;
	long val;
	static int cmp(const void *indOrdLA, const void *indOrdLB);
	int operator < (const L_IndOrdL &b) const {return val < b.val;}
	int operator <= (const L_IndOrdL &b) const {return val <= b.val;}
	int operator > (const L_IndOrdL &b) const {return val > b.val;}
	int operator >= (const L_IndOrdL &b) const {return val >= b.val;}
};




class L_ArgException:public std::exception //! For logical exceptions, internal errors in code
{
public:
    L_ArgException() nothrows() {}
    L_ArgException(const L_ArgException& other) nothrows() {}
    L_ArgException& operator =(const L_ArgException& other) nothrows() {return *this;}
	~L_ArgException () nothrows() {}
	const char* what () const nothrows() {return "L_ArgException: argumento no valido.";};
};

int throw_L_ArgException_fn_NULL(const char *s);
int throw_L_ArgException_fn(const char *s);
int printing_L_ArgException_fn(const char *s);

class L_HardErrorException:public std::exception //! Exception by zero size allocation
{
public:
	char *message;
    int messageLen;
    char *fileName;
    int fileNameLen;
    long lineOfCode;

    L_HardErrorException() nothrows() {message=NULL; messageLen=0; fileName=NULL; fileNameLen=0; lineOfCode=-1;}

	// Las llaves están puestas abreviadas para comprimir el código
	L_HardErrorException(const char *file, long line, const char *message, va_list ppp) nothrows();
    L_HardErrorException(const L_HardErrorException& other) nothrows();
    L_HardErrorException& operator =(L_HardErrorException& other) nothrows();
	void print(FILE *fp=stdout)
	{
		if (L_compilationLanguage() == L_espanol)
			fprintf(fp, "L_HardErrorException() en el archivo %s, linea de codeMapping %ld.\n",fileName,lineOfCode);
		else
			fprintf(fp, "L_HardErrorException() at fileName %s, lineOfCode %ld.\n",fileName,lineOfCode);
		fprintf(fp, "message: %s\n",message);
	}
	void pause()
	{
		if (L_compilationLanguage() == L_espanol)
			printf("Presione ENTER para salir\n");
		else
			printf("Press ENTER to exit\n");
		getchar();
	}

	~L_HardErrorException () nothrows()
    {
    	if (message!=NULL)
			delete[] message;
		if (fileName!=NULL)
			delete[] fileName;
    }
	const char* what () const nothrows() {return "L_HardErrorException";};
};

class L_ZeroDivisionException:public std::exception //! Division by zero
{
public:
    L_ZeroDivisionException() nothrows() {printf("L_ZeroDivisionException\n");}
    L_ZeroDivisionException(const L_ZeroDivisionException& other) nothrows() {}
    L_ZeroDivisionException& operator =(L_ZeroDivisionException& other) nothrows() {return *this;}
	~L_ZeroDivisionException () nothrows() {}
	const char* what () const nothrows() {return "L_ZeroDivisionException";};
};

class L_NonImplementedException:public std::exception //! When calling non-implemented code
{
public:
	L_NonImplementedException() nothrows() {printf("L_NonImplementedException\n");}
    L_NonImplementedException(const L_NonImplementedException& other) nothrows() {}
    L_NonImplementedException& operator =(const L_NonImplementedException& other) nothrows() {return *this;}
	~L_NonImplementedException () nothrows() {}
	const char* what () const nothrows() {return "L_NonImplementedException";};
};

#define NRR_NTAB 32
class L_Rand
{
public:
	long idum2;                            // initialized below
	long iy;
	long iv[NRR_NTAB];
	long seed;
	long seed0;
	long numCalls;
	long numCalls2;

private:
	L_Rand(const L_Rand &dontUseThatConstructor);
public:

	L_Rand() {seed = 10; seed0 = 10; numCalls = -1;  numCalls2 = -1; iy = 0;}
	void init(long semillaNeg) {throw_L_ArgException_if(semillaNeg > 0, "L_Rand() : inicializacion con seed positiva"); this->seed = semillaNeg; this->seed0 = semillaNeg; numCalls = 0; numCalls2 = 0; iy = 0;}
	double _ran2();
	double random(double a, double b) {return a + _ran2()*(b-a);}
	void iterate_until_state(long seed0, int numCalls, int numCalls2) {this->seed = seed0; this->numCalls = 0; this->numCalls2 = 0; iy = 0; while (this->numCalls < numCalls || this->numCalls2 < numCalls2) _ran2();}
};


template <class T, int N>
class L_CommaListParser
{
public:
	T *p;
	L_CommaListParser(T *pe) {p = pe;}
	L_CommaListParser<T,N-1> operator = (T d) {*p++ = d; return L_CommaListParser<T,N-1>(p);}
	L_CommaListParser<T,N-1> operator , (T d) {*p++ = d; return L_CommaListParser<T,N-1>(p);}
};

template <class T>
class L_CommaListParser<T, -1>
{
};


// Templates for overloading operator[] for debugging

template <class T>
class L_ptr
{
	int v110;
	T *el;
	int v111;
	long offset;
	int v112;
	L_ptr(const L_ptr &dontCallMe);
public:
	L_ptr(int &n) : el(NULL), v110(110), v111(111), v112(112) {offset = (long)((char *)this - (char *)&n);}
	L_ptr &operator=(const L_ptr &other) {el=other.el; return *this;}
	T &operator[](char i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[i];}
	T &operator[](unsigned char i) {if (i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[i];}
	T &operator[](int i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[i];}
	T &operator[](unsigned int i) {int ia = i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ia];}
	T &operator[](long i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[i];}
	T &operator[](unsigned long i) {long ia = i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ia];}

	T &_ref() {return el;}
	T* operator=(T* elem) {el=elem; return el;}
	operator T*() {return el;}
	operator const T*() const {return el;}
};


template <class T>
class L_ptrMint
{
	int v110;
	T **el;
	int v111;
	int offset;
	int v112;
	mutable int ulti;
	int v113;

	L_ptrMint(); // Don't call this
	L_ptrMint(const L_ptrMint &dontCallMe);
public:
	L_ptrMint(int &lj) : el(NULL), v110(110), v111(111), v112(112), v113(113) {offset = (long)((char *)this - (char *)&lj);}
	L_ptrMint<T> &operator=(const L_ptrMint<T> &other) {el=other.el; return *this;}
	T &operator[](char j) {if (j<0 || j>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][j];}
	T &operator[](unsigned char j) {if (j>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][j];}
	T &operator[](int j) {if (j<0 || j>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][j];}
	T &operator[](unsigned int j) {int ja = (long) j; if (ja<0 || ja>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][ja];}
	T &operator[](long j) {if (j<0 || j>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][j];}
	T &operator[](unsigned long j) {long ja = (long) j; if (ja<0 || ja>=(*(int *)((char *)this-offset))) *(char *)NULL=5; return el[ulti][ja];}
	T** operator=(T** elem) {el=elem; return el;}
	operator T*() {return el[ulti];}
	operator const T*() const {return el[ulti];}
	template <class U> friend class L_ptrM;
};

template <class T>
class L_ptrM
{
	int v110;
	L_ptrMint<T> memo;
	int v111;
	int offset;
	int v112;

	L_ptrM(); // Don't call this
	L_ptrM(const L_ptrM &dontCallMe);
public:
	L_ptrM(int &li, int &lj) : memo(lj), v110(110), v111(111), v112(112) {offset = (long)((char *)this - (char *)&li);}
	L_ptrM<T> &operator=(const L_ptrM<T> &other) {memo = other.memo; return *this;}
	L_ptrMint<T> &operator[](char i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	L_ptrMint<T> &operator[](unsigned char i) {if (i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	L_ptrMint<T> &operator[](int i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	L_ptrMint<T> &operator[](unsigned int i) {int ia = (long) i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=ia; return memo;}
	L_ptrMint<T> &operator[](long i) {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	L_ptrMint<T> &operator[](unsigned long i) {long ia = (long) i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=ia; return memo;}
	const L_ptrMint<T> &operator[](char i) const {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	const L_ptrMint<T> &operator[](unsigned char i) const {if (i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	const L_ptrMint<T> &operator[](int i) const {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	const L_ptrMint<T> &operator[](unsigned int i) const {int ia = (long) i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=ia; return memo;}
	const L_ptrMint<T> &operator[](long i) const {if (i<0 || i>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=i; return memo;}
	const L_ptrMint<T> &operator[](unsigned long i) const {long ia = (long) i; if (ia<0 || ia>=(*(int *)((char *)this-offset))) *(char *)NULL=5; memo.ulti=ia; return memo;}
	T** operator=(T** elem) {memo.el=elem; return memo.el;}
	operator T**() {return memo.el;}
	operator T**() const {return memo.el;}
};


template <int LI, int LJ> 
class L_ptrMEint
{
public:
	int v110;
	double el[LI][LJ];
	int v111;
	int _lj;
	int v112;
	mutable int ulti;
	int v113;
	L_ptrMEint() :_lj(LJ), v110(110), v111(111), v112(112) {}
	double &operator[](char j) {char *u=NULL; if (j<0 || j>=_lj) *u=5; return el[ulti][j];}
	double &operator[](unsigned char j) {char *u=NULL; if (j>=_lj) *u=5; return el[ulti][j];}
	double &operator[](int j) {char *u=NULL; if (j<0 || j>=_lj) *u=5; return el[ulti][j];}
	double &operator[](unsigned int j) {int ja = (long) j; char *u=NULL; if (ja<0 || ja>=_lj) *u=5; return el[ulti][ja];}
	double &operator[](long j) {char *u=NULL; if (j<0 || j>=_lj) *u=5; return el[ulti][j];}
	double &operator[](unsigned long j) {long ja = (long) j; char *u=NULL; if (ja<0 || ja>=_lj) *u=5; return el[ulti][ja];}
	operator double*() {return el[ulti];}
	operator const double*() const {return el[ulti];}
};

template <int LI, int LJ> 
class L_ptrME
{
	int v110;
	L_ptrMEint<LI,LJ> memo;
	int v111;
	int _li;
	int v112;

public:
	L_ptrME() : _li(LI) ,v110(110), v111(111), v112(112) {}
	L_ptrMEint<LI,LJ> &operator[](char i) {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	L_ptrMEint<LI,LJ> &operator[](unsigned char i) {char *u=NULL; if (i>=_li) *u=5; memo.ulti=i; return memo;}
	L_ptrMEint<LI,LJ> &operator[](int i) {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	L_ptrMEint<LI,LJ> &operator[](unsigned int i) {int ia = (long) i; char *u=NULL; if (ia<0 || ia>=_li) *u=5; memo.ulti=ia; return memo;}
	L_ptrMEint<LI,LJ> &operator[](long i) {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	L_ptrMEint<LI,LJ> &operator[](unsigned long i) {long ia = (long) i; char *u=NULL; if (ia<0 || ia>=_li) *u=5; memo.ulti=ia; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](char i) const {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](unsigned char i) const {char *u=NULL; if (i>=_li) *u=5; memo.ulti=i; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](int i) const {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](unsigned int i) const {int ia = (long) i; char *u=NULL; if (ia<0 || ia>=_li) *u=5; memo.ulti=ia; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](long i) const {char *u=NULL; if (i<0 || i>=_li) *u=5; memo.ulti=i; return memo;}
	const L_ptrMEint<LI,LJ> &operator[](unsigned long i) const {long ia = (long) i; char *u=NULL; if (ia<0 || ia>=_li) *u=5; memo.ulti=ia; return memo;}
	operator double*() {return &memo.el[0][0];}
	operator double*() const {return &memo.el[0][0];}
};


template <class T> class L_ArrayInverseIterator
{
	T* p;
public:
	L_ArrayInverseIterator(T* x) :p(x) {}
	L_ArrayInverseIterator& operator++() {--p;return *this;}
	L_ArrayInverseIterator& operator--() {++p;return *this;}
	bool operator==(const L_ArrayInverseIterator& rhs) {return p==rhs.p;}
	bool operator!=(const L_ArrayInverseIterator& rhs) {return p!=rhs.p;}
	T& operator*() {return *p;}
};


// Memory growth model for L_Array
#define L_Array_calc_allocated_max(arr) ((arr).growthMem == 0  ? (int)(arr).size()  :  (arr).growthMem + 4*(int)(arr).size())
#define L_Array_calc_allocated(arr)   ((arr).growthMem == 0 ? (int)(arr).size()  :  (arr).growthMem + 2*(int)(arr).size())

template <class T>
class L_Array
{
public:
	typedef T* iterator;
	typedef const T* const_iterator;
	typedef L_ArrayInverseIterator<T> inverse_iterator;
	typedef const L_ArrayInverseIterator<T> const_inverse_iterator;
	
	typedef T value_type;
	typedef int size_type;

	T *data() {return __elem;}
	const T *data() const {return __elem;}
	void setData(T *buf, int nel) {__elem = buf; __n = nel;}
	iterator begin() {return &(__elem[0]);}
	iterator end() {return &(__elem[0])+__n;}
	const_iterator begin() const {return &(__elem[0]);}
	const_iterator end() const {return &(__elem[0])+__n;}
	size_type size() const {return (size_type)__n;}
	void resize(size_type __n) {this->__n = __n; adjustMem();}
	void resize_swapping(size_type __n) {this->__n = __n; adjustMem_swap();}
	void clear() {destroy();}

#ifdef L_ARRAY_ACCESS_DEBUG
	T& operator[](size_type i) {throw_L_ArgException_if(i<0 || i>=__n || i>=nMem, "Out of bounds in L_Array"); return __elem[i];}
	const T& operator[](size_type i) const {throw_L_ArgException_if(i<0 || i>=__n || i>=nMem, "Out of bounds in L_Array"); return __elem[i];}
#else
	T& operator[](size_type i) {return __elem[i];}
	const T& operator[](size_type i) const {return __elem[i];}
#endif

	template <size_type N> L_CommaListParser<T,N> defineCommas() {resize(N); return L_CommaListParser<T,N>(data());}

#ifdef L_BRACKET_LEVEL_DEBUG
	L_Array(size_type __n, size_type growthMem = L_ARREGLO_ALLOC_SIZE_FACTOR) : __elem(this->nMem) {this->__elem = NULL; this->__n = __n; this->nMem = 0; this->growthMem = growthMem; if (this->__n > 0) adjustMemDown();}
	L_Array() : __elem(this->nMem) {this->__elem=NULL; this->__n=0; this->nMem=0; this->growthMem=L_ARREGLO_ALLOC_SIZE_FACTOR;}
	L_Array(const L_Array &other) : __elem(this->nMem) {printf("L_Array(const L_Array &other) : posible paso incorrecto a funcion\n"); __n=other.__n; nMem=other.nMem; growthMem=other.growthMem; if (other.__elem!=NULL) {__elem=new T[nMem]; size_type i; for (i=0; i<__n; i++) __elem[i]=other.__elem[i];}}
#else
	L_Array(size_type __n, size_type growthMem = L_ARREGLO_ALLOC_SIZE_FACTOR) {this->__elem = NULL; this->__n = __n; this->nMem = 0; this->growthMem = growthMem; if (this->__n > 0) adjustMemDown();}
	L_Array() {this->__elem=NULL; this->__n=0; this->nMem=0; this->growthMem=L_ARREGLO_ALLOC_SIZE_FACTOR;}
	L_Array(const L_Array &other) {printf("L_Array(const L_Array &other) : posible paso incorrecto a funcion\n"); __n=other.__n; nMem=other.nMem; growthMem=other.growthMem; if (other.__elem!=NULL) {__elem=new T[nMem]; size_type i; for (i=0; i<__n; i++) __elem[i]=other.__elem[i];}}
#endif

	void swap(L_Array &other)
	{
		T *elem_t;
		size_type n_t;
		size_type nMem_t;
		size_type growthMem_t;
		// tmp = este
		elem_t = __elem;
		n_t = __n;
		nMem_t = nMem;
		growthMem_t = growthMem;
		// este = other
		__elem = other.__elem;
		__n = other.__n;
		nMem = other.nMem;
		growthMem = other.growthMem;
		// other = tmp
		other.__elem = elem_t;
		other.__n = n_t;
		other.nMem = nMem_t;
		other.growthMem = growthMem_t;
	}

	void OP_assign(const L_Array<T> &other) { if (__elem==NULL) {if (other.__elem==NULL) return; __elem=new T[other.nMem];} else if (nMem != other.nMem) {delete[] (T*)__elem; __elem=new T[other.nMem];}; __n=other.__n; nMem=other.nMem; growthMem=other.growthMem; size_type i; for (i=0; i<__n; i++) __elem[i]=other.__elem[i];}
	L_Array &operator =(const L_Array &other) {OP_assign(other); return *this;}

	void reserve(size_type nelem)
	{
		size_type n0 = __n;
		if (nelem > __n)
			__n = nelem;
		adjustMem();
		__n = n0;
	}

	void reserve_swap(size_type nelem)
	{
		size_type n0 = __n;
		if (nelem > __n)
			__n = nelem;
		adjustMem_swap();
		__n = n0;
	}

	void adjustMemDown()
	{
		if (__n < nMem && nMem <= L_Array_calc_allocated_max(*this))
			return;
		size_type newAllocNumber = L_Array_calc_allocated(*this);
#ifdef L_BASIC_DEBUG
		if (__n<0)
			L_hard_shutdown("L_Array<T>::adjustMemDown() : number of elements is negative");
		if (newAllocNumber<__n)
			L_hard_shutdown("L_Array<T>::adjustMemDown() : internal error");
#endif
		intern_reserve(newAllocNumber);
	}

	void adjustMemDown_swap()
	{
		if (__n < nMem && nMem <= L_Array_calc_allocated_max(*this))
			return;
		size_type newAllocNumber = L_Array_calc_allocated(*this);
#ifdef L_BASIC_DEBUG
		if (__n<0)
			L_hard_shutdown("L_Array<T>::adjustMemDown_swap() : number of elements is negative");
		if (newAllocNumber<__n)
			L_hard_shutdown("L_Array<T>::adjustMemDown_swap() : internal error");
#endif
		intern_reserve_swap(newAllocNumber);
	}


	void getVector(std::vector<T> &vec) const
	{
		vec.reserve(nMem);
		vec.clear();
		size_type i;
		for (i=0; i<__n; i++)
			vec.push_back(__elem[i]);		
	}

	void setVector(const std::vector<T> &vec)
	{
		__n=vec.capacity();
		if (vec.capacity() > 2*vec.size())
			__n=(size_type)(vec.size()*1.2);
		if (__n==0)
			return;
		adjustMem();
		for (__n=0; __n<vec.size(); __n++)
			__elem[__n]=vec[__n];
	}

	void push_back(const T &nuevo) {__n++; if (__n > nMem) intern_reserve(L_Array_calc_allocated(*this)); __elem[__n-1]=nuevo;}
	void push_back_swapping(const T &nuevo) {__n++; if (__n > nMem) intern_reserve_swap(L_Array_calc_allocated(*this)); nuevo.swap(__elem[__n-1]);}

	void erase_preserving_order(size_type i)
	{
		for(; i<size()-1; i++)
			(*this)[i] = (*this)[i+1];
		resize(size()-1);
	}
	void erase_preserving_order_swapping(size_type i)
	{
		for(; i<size()-1; i++)
			(*this)[i].swap((*this)[i+1]);
		resize_swapping(size()-1);
	}
	void erase_fast_non_preserving_order(size_type i)
	{
		(*this)[i] = (*this)[size()-1];
		resize(size()-1);
	}
	void erase_fast_non_preserving_order_swapping(size_type i)
	{
		(*this)[size()-1].swap((*this)[i]);
		resize_swapping(size()-1);
	}
	void erase_preserving_order(const std::vector<int> &indexes)
	{
		size_type i, u=0;
		for (i=0; i<size()-(size_type)indexes.size(); i++)
		{
			while (u < (size_type)indexes.size() && indexes[u] == i+u)
				u++;
			if (u > 0)
				(*this)[i] = (*this)[i+u];
		}
		resize(size()-(size_type)indexes.size());
	}
	void erase_preserving_order(const L_Array<bool> &elim)
	{
		size_type i=0, u=0;
		while(true)
		{
			while (i<size() && elim[i] == true)
				i++;
			if (i==__n)
				break;
			while (i<size() && elim[i] == false)
			{
				(*this)[u] = (*this)[i];
				i++;
				u++;
			}
			if (i==size())
				break;
		}
		resize(u);
	}
	void erase_preserving_order_swapping(const std::vector<int> &ind)
	{
		size_type i, u=0;
		for (i=0; i<size()-(size_type)ind.size(); i++)
		{
			while (u < (size_type)ind.size() && ind[u] == i+u)
				u++;
			if (u > 0)
				(*this)[i+u].swap((*this)[i]);
		}
		resize_swapping(size()-(size_type)ind.size());
	}

	void fwriteArr_POD(FILE *fp)  // Only use for POD elements (no pointers, no allocations, no virtual functions)
	{
		fwrite(&__n,sizeof(__n),1,fp);
		fwrite(&nMem,sizeof(__n),1,fp);
		fwrite(&growthMem,sizeof(__n),1,fp);
		if (__n>0)
			fwrite(__elem, sizeof(__elem[0]), __n, fp);
	}
	bool freadArr_POD(FILE *fp)  // Only use for POD elements (no pointers, no allocations, no virtual functions)
	{
		size_type nMemmm;
		bool ret = true;
		ret |= fread(&__n,sizeof(__n),1,fp)>0;
		ret |= fread(&nMemmm,sizeof(__n),1,fp)>0;
		ret |= fread(&growthMem,sizeof(__n),1,fp)>0;
		intern_reserve(nMemmm);
		if (__n>0)
			ret |= fread(__elem, sizeof(__elem[0]), __n, fp) == (size_t)__n;  // Se leyeron los __n elementos
		return ret;
	}

	void sort(size_type (*qCmp)(const void *pT1, const void *pT2));
	size_type isolateRepeatedElementsAtEnd(int (*qCmp)(const void *pT1, const void *pT2));
	size_type sortIsolatingRepeatedAtEnd(int (*qCmp)(const void *pT1, const void *pT2)) {sort(qCmp); return isolateRepeatedElementsAtEnd(qCmp);}; // Deja los elementos indistinguibles al final :  // size_type qCmp(T *p1, T *p2)
	
	void randomPermutation() {L_randomPermutation(__elem, __n);}

	size_type find_in_sorted_array(const T& val, size_type (*qCmp)(const void *pT1, const void *pT2)) {return (size_type)( ((T*)bsearch(&val, __elem, __n, sizeof(T), qCmp)) - __elem );}

	void assign(int num, const T &zero) {size_type i; resize(num); for (i=0; i<num; i++) __elem[i] = zero;}
	void setValueFromRange(const T &inic, const T&delta, const T &final) // arr=inic:delta:final
		{T val; size_type i; __n=(size_type)((final-inic)/delta)+1; adjustMem(); val=inic; for (i=0; i<__n; i++) {__elem[i]=val; val+=delta;}}

	T sum() const {T s=__elem[0]; size_type i; for(i=1; i<__n; i++) s+=__elem[i]; return s;}
	T sumSquares() const {T s=__elem[0]; size_type i; for(i=1; i<__n; i++) s+=__elem[i]*__elem[i]; return s;}
	T mean() const {T s=__elem[0]; size_type i; for(i=1; i<__n; i++) s+=__elem[i]; return s/__n;} // Si es arreglo de size_type devuelve un size_type
	T variance() const {T s=__elem[0], s2=__elem[0]*__elem[0]; size_type i; for(i=1; i<__n; i++) {s+=__elem[i]; s2+=__elem[i]*__elem[i];} return s2/__n - s*s/(__n*__n);}
	void sum_to_each_element(const T &e) {size_type i; for (i=0; i<__n; i++) __elem[i] +=e;}
	void subtract_to_each_element(const T &e) {size_type i; for (i=0; i<__n; i++) __elem[i] -= e;}
	void multiply_to_each_element(double factor) {size_type i; for (i=0; i<__n; i++) __elem[i] *= factor;}

	void sumElement(const L_Array<T> &other, const T &e) {__n=other.__n; adjustMem(); size_type i; for (i=0; i<__n; i++) __elem[i] = other.__elem[i] + e;}
	void subtractElement(const L_Array<T> &other, const T &e) {__n=other.__n; adjustMem(); size_type i; for (i=0; i<__n; i++) __elem[i] = other.__elem[i] - e;}
	void sum(const L_Array<T> &first, const L_Array<T> &second) {throw_L_ArgException_if(first.__n!=second.__n, "sumaDe"); __n=first.__n; adjustMem(); size_type i; for (i=0; i<__n; i++) __elem[i] = first.__elem[i] + second.__elem[i];}
	void subtract(const L_Array<T> &first, const L_Array<T> &second) {throw_L_ArgException_if(first.__n!=second.__n, "restaDe"); __n=first.__n; adjustMem(); size_type i; for (i=0; i<__n; i++) __elem[i] = first.__elem[i] - second.__elem[i];}
	void mult(const L_Array<T> &first, const L_Array<T> &second) {throw_L_ArgException_if(first.__n!=second.__n,"multDe"); __n=first.__n; adjustMem(); size_type i; for (i=0; i<__n; i++) __elem[i] = first.__elem[i] * second.__elem[i];}
	void unionOf(const L_Array &arr1,const L_Array &arr2) {__n = arr1.__n+arr2.__n; adjustMem(); size_type i; for (i=0; i<arr1.__n; i++) __elem[i] = arr1.__elem[i]; for (i=0; i<arr2.__n; i++) __elem[i+arr1.__n] = arr2.__elem[i];}

	size_type countElementsDifferentTo(const T &val) const {size_type i, ene=0; for (i=0; i<__n; i++) ene+=(__elem[i]!=val); return ene;}
	size_type countElementsEqualTo(const T &val) const {size_type i, ene=0; for (i=0; i<__n; i++) ene+=(__elem[i]==val); return ene;}
	size_type countElementsLesserThan(const T &val) const {size_type i, ene=0; for (i=0; i<__n; i++) ene+=(__elem[i]<val); return ene;}
	size_type countElementsGreaterThan(const T &val) const {size_type i, ene=0; for (i=0; i<__n; i++) ene+=(__elem[i]>val); return ene;}
	T maxElement() const {size_type i; T max=__elem[0]; for(i=1; i<__n; i++) if (__elem > max) max = __elem; return max;}
	T minElement() const {size_type i;T min=__elem[0]; for(i=1; i<__n; i++) if (__elem < min) min = __elem; return min;}


	void print()
	{
		for (size_type i=0; i<__n; i++)
			std::cout << __elem[i] << "\n";
	}

	template <class Fn> void fill(Fn func_int) {size_type i; for (i=0; i<__n; i++) __elem[i] = func_int(i);}
	template <class Fn> void transformUsing(Fn func_T) {size_type i; for (i=0; i<__n; i++) __elem[i] = func_T(__elem[i]);}
	template <class Fn> void apply(Fn func_T) {size_type i; for (i=0; i<__n; i++) func_T(__elem[i]);}

	~L_Array() {if (__elem!=NULL) delete[] (T*)__elem; __elem=NULL;}

private:
	T &operator=(const T &dontCallMe);

	void intern_reserve(size_type newAllocNumber) // Allocating new space copying current elements to the new array
	{
		throw_L_ArgException_if(newAllocNumber<0, "L_Array::intern_reserve() : allocation of negative size");
		if (newAllocNumber == 0)
		{
			if (__elem!=NULL)
				delete[] (T*)__elem;
			__elem = NULL;
			nMem = 0;
			return;
		}
		if (__elem==NULL)
			__elem=new T[newAllocNumber];
		else
		{
			T *arr2=new T[newAllocNumber];
			size_type i;
			for (i=0; i<__n && i<nMem && i<newAllocNumber; i++)
				arr2[i]=__elem[i];
			delete[] (T*)__elem;
			__elem=arr2;
		}
		nMem = newAllocNumber;
	}

	void intern_reserve_swap(size_type newAllocNumber) // Allocating new space swapping current elements to the new array
	{
		throw_L_ArgException_if(newAllocNumber<0, "L_Array::intern_reserve() : peticion de cantidad de memoria negativa");
		if (newAllocNumber == 0)
		{
			if (__elem!=NULL)
				delete[] (T*)__elem;
			__elem = NULL;
			nMem = 0;
			return;
		}
		if (__elem==NULL)
			__elem=new T[newAllocNumber];
		else
		{
			T *arr2=new T[newAllocNumber];
			size_type i;
			for (i=0; i<__n && i<nMem && i<newAllocNumber; i++)
				__elem[i].swap(arr2[i]);
			delete[] (T*)__elem;
			__elem=arr2;
		}
		nMem = newAllocNumber;
	}

	void adjustMem()
	{
		if (__n < nMem) // Ya esta ok
			return;
		size_type newAllocNumber = L_Array_calc_allocated(*this);
#ifdef L_BASIC_DEBUG
#if L_ARREGLO_ALLOC_SIZE_FACTOR != 0
		if (newAllocNumber<=0)
			L_hard_shutdown("L_Array<T>::adjustMem : zero size allocation");
#endif
		if (newAllocNumber<__n)
			L_hard_shutdown("L_Array<T>::adjustMem() : internal error");
#endif
		intern_reserve(newAllocNumber);
	}

	void adjustMem_swap()
	{
		if (__n < nMem)
			return;
		size_type newAllocNumber = L_Array_calc_allocated(*this);
#ifdef L_BASIC_DEBUG
#if L_ARREGLO_ALLOC_SIZE_FACTOR != 0
		if (newAllocNumber<=0)
			L_hard_shutdown("L_Array<T>::adjustMem_swap : zero size allocation");
#endif
		if (newAllocNumber<__n)
			L_hard_shutdown("L_Array<T>::adjustMem_swap() : internal error");
#endif
		intern_reserve_swap(newAllocNumber);
	}

	void destroy() {if (__elem!=NULL) {delete[] (T*)__elem; __elem = NULL;}; __n=0; nMem=0;}

private:
#ifdef L_BRACKET_LEVEL_DEBUG
	L_ptr<T> __elem;
#else
	T *__elem;
#endif
	int __n; // Number of elements currently stored
	int nMem; // Number of elements currently allocated
	int growthMem; // Parameter for allocation size increment
};


inline void exampleOfUseOfArrays()
{
	L_Array<double> a, b;
	a.defineCommas<5>() = 1, 2, 3, 4, 5;
	a.print();
	b.resize(3);
	b[0] = 2;
	b[1] = 1;
	b[2] = 3;
	std::sort(b.begin(), b.end());
	b.print();

	// Selfdestruction...
}



template <class T>
int L_Array<T>::isolateRepeatedElementsAtEnd(int (*qCmp)(const void *pT1, const void *pT2)) //  // int qCmp(T *p1, T *p2)
{
	size_type i, j, ndie2, ndie = 0;

#ifdef L_BASIC_DEBUG
		if (size()<=0)
			L_hard_shutdown("L_Array<T>::eliminaRep : non-positive numbre of elements in array");
#endif

	bool *die=new bool[size()];

	// Array for indexing repeated elements (comparison function = 0)
	for (i=0; i<size(); i++)
		die[i] = false;

	for (i=0; i<size()-1; i++)
	{
		if (qCmp(&__elem[i],&__elem[i+1])==0)
		{
			die[i]=true;
			die[i+1]=true;
		}
	}

	for (i=0; i<size(); i++)
		if (die[i])
			ndie++;

	if (ndie == 0)  // No elements to isolate
	{
		delete[] die;
		return 0;
	}

	int *order = new int[size()];
	T *elem2 = new T[nMem];

	ndie2 = ndie;
	j=0;
	for (i=0; i<size(); i++)
	{
		if (die[i])
			order[i] = size() - ndie--; // Isolate repeated elements at the ending of the array
		else
			order[i] = j++;
	}

	for (i=0; i<size(); i++)
		elem2[order[i]]=__elem[i];
	delete[] die;
	delete[] order;
	delete[] (T*)__elem;
	__elem = elem2;
	return ndie2;
};

template <class T>
void L_Array<T>::sort(int (*qCmp)(const void *pT1, const void *pT2)) {qsort(__elem,__n,sizeof(T),qCmp);} // int qCmp(T *p1, T *p2)

//

template <class V>
class L_VectMap
{
private:
	L_VectMap();
public:
	V &v;
	typedef typename V::value_type T;
	typedef typename V::value_type value_type;
	typedef typename V::size_type size_type;

	L_VectMap(std::vector<T> &vect) : v(vect) { }

	void setValueElementWise(const T &zero) {typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = zero;}
	void setValueFromRange(const T &inic, const T&delta, const T &final) // arr=inic:delta:final
		{T val; typename V::size_type i; v.resize((int)((final-inic)/delta)+1); val=inic; for (i=0; i<v.size(); i++) {v[i]=val; val+=delta;}}

	void sum_to_each_element(const T &e) {typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] +=e;}
	void subtract_to_each_element(const T &e) {typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] -= e;}
	void multiply_to_each_element(double factor) {typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] *= factor;}

	void sumElement(const std::vector<T> &other, const T &e) {v.resize(other.size()); typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = other[i] + e;}
	void subtractElement(const std::vector<T> &other, const T &e) {v.resize(other.size()); typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = other[i] - e;}
	void sum(const std::vector<T> &first, const std::vector<T> &second) {throw_L_ArgException_if(first.size()!=second.size(), "sumaDe"); v.resize(first.size()); typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = first[i] + second[i];}
	void subtract(const std::vector<T> &first, const std::vector<T> &second) {throw_L_ArgException_if(first.size()!=second.size(), "restaDe"); v.resize(first.size()); typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = first[i] - second[i];}
	void mult(const std::vector<T> &first, const std::vector<T> &second) {throw_L_ArgException_if(first.size()!=second.size(),"multDe"); v.resize(first.size()); typename std::vector<T>::size_type i; for (i=0; i<v.size(); i++) v[i] = first[i] * second[i];}

	void fwriteArr_POD(FILE *fp)  // Only use for POD elements (no pointers, no allocations, no virtual functions)
	{
		size_type __n = (size_type)v.size();
		size_type __nMem = (size_type) v.capacity();
		size_type _growthMem = 10;
		fwrite(&__n,sizeof(__n),1,fp);
		fwrite(&__nMem,sizeof(__n),1,fp);
		fwrite(&_growthMem,sizeof(__n),1,fp);
		for (size_type i=0; i<v.size(); i++)
			fwrite(&v[i], sizeof(v[i]), 1, fp);
	}
	bool freadArr_POD(FILE *fp)  // Only use for POD elements (no pointers, no allocations, no virtual functions)
	{
		bool ret = true;
		size_type __n;
		size_type __nMem = (size_type)v.capacity();
		size_type _growthMem = 10;
		ret |= (fread(&__n,sizeof(__n),1,fp) > 0);
		ret |= (fread(&__nMem,sizeof(__n),1,fp) > 0);
		ret |= (fread(&_growthMem,sizeof(__n),1,fp) > 0);
		v.reserve(__nMem);
		v.resize(__n);
		for (typename std::vector<T>::size_type i=0; i<(size_type)v.size(); i++)
			ret |= (fread(&v[i], sizeof(v[i]), 1, fp) > 0);
		return ret;
	}
	void erase_preserving_order(const std::vector<int> &indexes)
	{
		typename std::vector<int>::size_type i, u=0;
		for (i=0; i<v.size()-indexes.size(); i++)
		{
			while (u < indexes.size() && indexes[u] == i+u)
				u++;
			if (u > 0)
				v[i] = v[i+u];
		}
		v.resize(v.size()-indexes.size());
	}
	void erase_preserving_order(const std::vector<bool> &elim)
	{
		typename std::vector<bool>::size_type i=0, u=0;
		while(true)
		{
			while (i<v.size() && elim[i] == true)
				i++;
			if (i==v.size())
				break;
			while (i<v.size() && elim[i] == false)
			{
				v[u] = v[i];
				i++;
				u++;
			}
			if (i==v.size())
				break;
		}
		v.resize(u);
	}
	void erase_preserving_order_swapping(const std::vector<int> &ind)
	{
		typename std::vector<int>::size_type i, u=0;
		for (i=0; i<v.size()-ind.size(); i++)
		{
			while (u < ind.size() && ind[u] == i+u)
				u++;
			if (u > 0)
				v[i+u].swap(v[i]);
		}
		v.resize(v.size()-ind.size());
	}
	void erase_fast_non_preserving_order(size_type i)
	{
		v[i] = v[v.size()-1];
		v.resize(v.size()-1);
	}
	void erase_fast_non_preserving_order_swapping(size_type i)
	{
		v[v.size()-1].swap(v[i]);
		v.resize(v.size()-1);
	}
	int isolateRepeatedElementsAtEnd(int (*qCmp)(const void *pT1, const void *pT2));
	size_type sortIsolatingRepeatedAtEnd(int (*qCmp)(const void *pT1, const void *pT2)) {sort(qCmp); return isolateRepeatedElementsAtEnd(qCmp);}; // Deja los elementos indistinguibles al final :  // size_type qCmp(T *p1, T *p2)
	// Funcion no segura
	void sort(int (*qCmp)(const void *pT1, const void *pT2)) {qsort(&v[0],v.size(),sizeof(v[0]),qCmp);} // int qCmp(T *p1, T *p2)
};

template <class V>  // Use L_VM() or L_VMC()
class L_VectMapC
{
private:
	L_VectMapC();
public:
	const V &v;
	typedef typename V::value_type T;

	L_VectMapC(const std::vector<T> &vect) : v(vect) { }

	T sum() {T s=v[0]; typename V::size_type i; for(i=1; i<v.size(); i++) s+=v[i]; return s;}
	T sumSquares() {T s=v[0]; typename V::size_type i; for(i=1; i<v.size(); i++) s+=v[i]*v[i]; return s;}
	T mean() {T s=v[0]; typename V::size_type i; for(i=1; i<v.size(); i++) s+=v[i]; return s/v.size();} // Si es arreglo de int devuelve un int
	T variance() {T s=v[0], s2=v[0]*v[0]; typename V::size_type i; for(i=1; i<v.size(); i++) {s+=v[i]; s2+=v[i]*v[i];} return s2/v.size() - s*s/(v.size()*v.size());}

	int countElementsDifferentTo(const T &val) {typename V::size_type i, ene=0; for (i=0; i<v.size(); i++) ene+=(v[i]!=val); return ene;}
	int countElementsEqualTo(const T &val) {typename V::size_type i, ene=0; for (i=0; i<v.size(); i++) ene+=(v[i]==val); return ene;}
	int countElementsLesserThan(const T &val) {typename V::size_type i, ene=0; for (i=0; i<v.size(); i++) ene+=(v[i]<val); return ene;}
	int countElementsGreaterThan(const T &val) {typename V::size_type i, ene=0; for (i=0; i<v.size(); i++) ene+=(v[i]>val); return ene;}
	T maxElement() {typename V::size_type i; T max=v[0]; for(i=1; i<v.size(); i++) if (v[i] > max) max = v[i]; return max;}
	T minElement() {typename V::size_type i;T min=v[0]; for(i=1; i<v.size(); i++) if (v[i] < min) min = v[i]; return min;}
};

template <class V> L_VectMapC<V> L_VMC(V &v) {return L_VectMapC<V>(v);}
template <class V> L_VectMap<V> L_VM(V &v) {return L_VectMap<V>(v);}


template <class V>
int L_VectMap<V>::isolateRepeatedElementsAtEnd(int (*qCmp)(const void *pT1, const void *pT2)) //  // int qCmp(T *p1, T *p2)
{
	size_type i, j, ndie2, ndie = 0;

#ifdef L_BASIC_DEBUG
		if (size()<=0)
			L_hard_shutdown("L_Array<T>::eliminaRep : non-positive numbre of elements in array");
#endif

	bool *die=new bool[v.size()];

	// Array for indexing repeated elements (comparison function = 0)
	for (i=0; i<v.size(); i++)
		die[i] = false;

	for (i=0; i<v.size()-1; i++)
	{
		if (qCmp(&v[i],&v[i+1])==0)
		{
			die[i]=true;
			die[i+1]=true;
		}
	}

	for (i=0; i<v.size(); i++)
		if (die[i])
			ndie++;

	if (ndie == 0)  // No elements to isolate
	{
		delete[] die;
		return 0;
	}

	int *order = new int[v.size()];
	V elem2(v.size());

	ndie2 = ndie;
	j=0;
	for (i=0; i<v.size(); i++)
	{
		if (die[i])
			order[i] = (int)v.size() - (int)ndie--; // Isolate repeated elements at the ending of the array
		else
			order[i] = (int)j++;
	}

	for (i=0; i<v.size(); i++)
		elem2[order[i]]=v[i];
	delete[] die;
	delete[] order;

	v.swap(elem2);
	return (int)ndie2;
};



template <class V1, class V2>
void L_complementOf(const V1 &orig, V2 &complem, int nTot)
{
	int i;
	int u=0, v=0;
	std::vector<int> arr;
	arr.resize(orig.size());
	for (i=0; i<(int)arr.size(); i++)
		arr[i] = orig[i];

	L_VM(arr).sort(L_Int::cmp); // arr[] = orig[] pero ordered
	complem.resize(nTot); // ACA DECIA orig.size()

	for (i=0; i<nTot; i++)
	{
		if (u < (int)arr.size() && arr[u] == i)
			while (u < (int)arr.size() && arr[u] == i) // Por si hay varios repetidos en arr[]
				u++;
		else
			complem[v++] = i;
	}
	complem.resize(v);
	throw_L_ArgException_if(u < (int)orig.size() && arr[u] >= nTot, "L_complementOf(): argumentos incompatibles (E u | orig[u] > nTot)");
	throw_L_ArgException_if(u != (int)orig.size(), "L_complementOf(): error interno");
}


int L_searchInSortedArray(const std::vector<double> valOrd, double val);
void L_randomSelection(std::vector<int> &seleccion, int nCota);
void L_shiftPermutation_right(std::vector<int> &v, std::vector<int> &stackTmp); // Va cambiando los elementos a la derecha primero
void L_shiftPermutation_left(std::vector<int> &v, std::vector<int> &stackTmp); // Va cambiando los elementos a la izquierda primero
void L_getOrderedIndexesMajorMinor(std::vector<int> &inds, const std::vector<double> &puntaje); // Devuelve el orden de v[] en inds[]
void L_getOrderedIndexesMinorMajor(std::vector<int> &inds, const std::vector<double> &puntaje); // Devuelve el orden de v[] en inds[]

void L_calcGausH(std::vector<double> &buf, double s, int width, int facpaso=1);    //! compute kernel gaussiano discreto
void L_calcGausHNorm(std::vector<double> &buf, double s, int width, int facpaso=1);    //! compute kernel gaussiano discreto normalizado (sum 1)
void L_calcGausHInt(std::vector<double> &buf, double s, int width, int facpaso=1); //! compute kernel gaussiano aproximado continuo

// Gabor
void L_calcGaborSeparable_cos_cos(std::vector<double> &h, double s, double f, double theta, int width, int facpaso=1);
void L_calcGaborSeparable_cos_sin(std::vector<double> &h, double s, double f, double theta, int width, int facpaso=1);
void L_calcGaborSeparable_sin_cos(std::vector<double> &h, double s, double f, double theta, int width, int facpaso=1);
void L_calcGaborSeparable_sin_sin(std::vector<double> &h, double s, double f, double theta, int width, int facpaso=1);


class L_String // Simple string class, make copies all the times.
{
public:
	L_String() {}
	L_String(const L_String &ret) {sv = ret.sv;}

	L_String(const char *txt)
	{
		throw_L_ArgException_if(txt==NULL, "L_String(const char *txt) : argumento NULL");
		resize((int)strlen(txt) + 1);
		strcpy(data(),txt);
	}
	L_String(short num, int numberWidthZeros=1) {setInteger(num, numberWidthZeros);}
	L_String(int num, int numberWidthZeros=1)  {setInteger(num, numberWidthZeros);}
	L_String(unsigned short num, int numberWidthZeros=1)  {setInteger(num, numberWidthZeros);}
	L_String(unsigned int num, int numberWidthZeros=1)  {setInteger(num, numberWidthZeros);}
	L_String(long num, int numberWidthZeros=1)  {setInteger(num, numberWidthZeros);}
	L_String(float num, int nDecimals=4)  {setDouble(num, nDecimals);}
	L_String(double num, int nDecimals=4)  {setDouble(num, nDecimals);}

	typedef std::string::size_type size_type;

#ifdef L_STRING_ACCESS_DEBUG
	char &operator[](size_type i) {throw_L_ArgException_if(i<0 || i>=sv.size(), "L_String: Out of bounds"); return sv[i];}
	const char &operator[](size_type i) const {throw_L_ArgException_if(i<0 || i>=sv.size(), "L_String: Out of bounds"); return sv[i];}
#else
	char &operator[](size_type i) {return sv[i];}
	const char &operator[](size_type i) const {return sv[i];}
#endif

	char *data() {if (sv.size() == 0) return NULL; return &(sv[0]);}
	const char *data() const {if (sv.size() == 0) return NULL; return &(sv[0]);}
	size_type size() const {return (int)sv.size();}
	void resize(size_type num) {sv.resize(num);}

	char *begin() {return &(sv[0]);}
	char *end() {return &(sv[0])+sv.size();}
	const char *begin() const {return &(sv[0]);}
	const char *end() const {return &(sv[0])+sv.size();}

	void clear() {sv.clear();}

	void push_back(const char &c) {resize(size()+1); data()[size()-2] = c; data()[size()-1] = 0;}

	void setInteger(long num, int numberWidthZeros=1);
	void setDouble(double num, int nDecimals = 4);
	bool find(const char *subString, int &pos);
	void copySubString(const char *original, size_type i1, size_type i2) {size_type i; if (i1>i2) {resize(1);} else {resize(i2-i1+1);} for ( i=0 ; ((*this)[i]=original[i+i1])!=0 && i+i1<=i2 ; i++ ); (*this)[i]=0; resize(i+1);}
	void uppercase() {size_type i; for (i=0; ((*this)[i]!=0 && i<size()-1); i++) if ((*this)[i] >= 'a' && (*this)[i] <='z') (*this)[i] = (*this)[i] + 'A' - 'a'; return;}
	void lowercase() {size_type i; for (i=0; ((*this)[i]!=0 && i<size()-1); i++) if ((*this)[i] >= 'A' && (*this)[i] <='Z') (*this)[i] = (*this)[i] + 'a' - 'A'; return;}
	void erase_spacing_characters_at_beginning_and_at_ending(); // Saca caracteres ' ', '\t', '\n' al inicio y al final del string
	int sprintf_cpy(size_type lenmax, const char *format, ...);
	int sprintf_cat(size_type lenmax, const char *format, ...);
	void encript(unsigned long semilla=0);
	void codifBase16(const char *txt);
	void decodifBase16(const char *txt);
	L_String& operator =(const L_String &s);
	friend L_String operator+ (const L_String &s1, const L_String &s2);
	L_String& operator +=(const L_String &s);

	void tokens(L_Array<L_String> &arrStr, const char *delimitadores);

	static bool saveStringToFile(const char *name, const char *string) {FILE *fp = fopen(name, "w"); if (fp == NULL) return false; fputs(string, fp); fclose(fp); return true;}

	//operator const char* () const { return elem; }
	const char *c_str() const {return data();}
	void swap(L_String &other) {sv.swap(other.sv);}

	static char cHexDigit(int num) {if (num<10) return (char)('0'+num); else return (char)('A'+num-10);}
	static int cHexDigitInv(char hexD) {if (hexD<'A') return (char)(hexD-'0'); else if (hexD<'a') return (char)(hexD+10-'A'); else return (char)(hexD+10-'a');}
private:
	std::vector<char> sv;
};




class L_ComplexDouble
{
public:
	double re;
	double im;
	L_ComplexDouble() {}
	L_ComplexDouble(double real, double imag) {re=real; im=imag;}

	inline void operator+=(const L_ComplexDouble &other)
	{
		re+=other.re;
		im+=other.im;
	}
	inline void operator-=(const L_ComplexDouble &other)
	{
		re-=other.re;
		im-=other.im;
	}
	inline void operator*=(const L_ComplexDouble &other)
	{
		double _re=re*other.re-im*other.im;
		im=re*other.im+im*other.re;
		re=_re;
	}
	inline void operator/=(const L_ComplexDouble &other)
	{
		double abs2=other.re*other.re+other.im+other.im;
		// a/b = a*b'/|b|^2 => other.im with flipped sign
		double _re=( re*other.re+im*other.im )/abs2;
		im=( -re*other.im+im*other.re )/abs2;
		re=_re;
	}

	inline void operator+=(double other) {re+=other;}
	inline void operator-=(double other) {re-=other;}
	inline void operator*=(double other) {re*=other; im*=other;}
	inline void operator/=(double other) {re/=other; im/=other;}

	inline void copyTo(L_ComplexDouble &other) const
	{
		other.re=re;
		other.im=im;
	}
	inline void define(double real, double imag) {re=real; im=imag;}
	inline void definePolar(double radio, double angRadianes) {re=radio*cos(angRadianes); im=radio*sin(angRadianes);}
	inline void calcExp()
	{
		double abs=exp(re);
		re=abs*cos(im);
		im=abs*sin(im);
	}
	inline void calcLog()
	{
		double re_=0.5*log(re*re+im*im);
		if (im==0 && re==0)
			im=0;
		else
			im=atan2(im,re);
		re=re_;
	}
	inline void calcCos()
	{
		double re_=cos(re)*cosh(im);
		im=-sin(re)*sinh(im);
		re=re_;
	}
	inline void calcSin()
	{
		double re_=cos(re)*cosh(im);
		im=cos(re)*sinh(im);
		re=re_;
	}
	inline void calcCosh()
	{
		double re_=cosh(re)*cos(im);
		im=sinh(re)*sin(im);
		re=re_;
	}
	inline void calcSinh()
	{
		double re_=sinh(re)*cos(im);
		im=cosh(re)*sin(im);
		re=re_;
	}
	inline void calcPow(const L_ComplexDouble &other)
	{
		calcLog();
		(*this)*=other;
		calcExp();
	}

	inline void calcSqrt() // two solutions: calcSqrt() and -calcSqrt()
	{
		double r, ang, r2, ang2;
		r=sqrt(re*re + im*im);
		ang=atan2(im, re);
		r2 = sqrt(r);
		ang2 = ang / 2;
		re = r2*cos(ang2);
		im = r2*sin(ang2);
	}

	inline void calcSqrt(double d) // two solutions: calcSqrt() and -calcSqrt()
	{
		if (d>0)
		{
			re = sqrt(d);
			im = 0;
		}
		else
		{
			re = 0;
			im = sqrt(d);
		}
	}

	inline void calcCubicRoot(L_ComplexDouble &res1, L_ComplexDouble &res2, L_ComplexDouble &res3)  // three solutions: calcCubicRoot() and calcCubicRoot()*(1<120°) and calcCubicRoot()*(1<40°)
	{
		double r, ang, r3, ang3;
		r=sqrt(re*re + im*im);
		ang=atan2(im, re);
		r3 = pow(r, 1.0/3);
		ang3 = ang / 3;
		res1.re = r3*cos(ang3);
		res1.im = r3*sin(ang3);
		res2.re = r3*cos(ang3 + 2*M_PI/3);
		res2.im = r3*sin(ang3 + 2*M_PI/3);
		res3.re = r3*cos(ang3 + 2 * 2*M_PI/3);
		res3.im = r3*sin(ang3 + 2 * 2*M_PI/3);
	}
	inline void calcCubicRoot() //   // three solutions: calcCubicRoot() and calcCubicRoot()*(1<120°) and calcCubicRoot()*(1<40°)
	{
		double r, ang, r3, ang3;
		r=sqrt(re*re + im*im);
		ang=atan2(im, re);
		r3 = pow(r, 1.0/3);
		ang3 = ang / 3;
		re = r3*cos(ang3);
		im = r3*sin(ang3);
	}
	double abs() {return sqrt(re*re+im*im);}
	double abs2() {return re*re+im*im;}
	double ang() {if (im==0 && re==0) return 0; else return atan2(im,re);}
	double dist(const L_ComplexDouble &other) {L_ComplexDouble x(this->re-other.re,this->im-other.im); return x.abs();}
	bool isReal(double err) {if (re==0) {if (im==0) return true; return false;} double d = im/re; if (d<0) d=-d; return d<err;}	
	bool isNaN() {return re!=re || im!=im;}

	static int solveQuadraticEquation(double a, double b, double c, L_ComplexDouble &x1, L_ComplexDouble &x2);
	static int solveCubicEquation(double a, double b, double c, double d, L_ComplexDouble &x1, L_ComplexDouble &x2, L_ComplexDouble &x3);
	static int solveQuarticEquation(double a, double b, double c, double d, double e, L_ComplexDouble &x1, L_ComplexDouble &x2, L_ComplexDouble &x3, L_ComplexDouble &x4);
	static L_ComplexDouble csqrt(const L_ComplexDouble &a) {L_ComplexDouble b(a); b.calcSqrt(); return b;}
	static L_ComplexDouble ccurt(const L_ComplexDouble &a) {L_ComplexDouble b(a); b.calcCubicRoot(); return b;}

	L_ComplexDouble operator +(const L_ComplexDouble &b) const {return L_ComplexDouble(re+b.re, im+b.im);}
	L_ComplexDouble operator -(const L_ComplexDouble &b) const {return L_ComplexDouble(re-b.re, im-b.im);}
	L_ComplexDouble operator *(const L_ComplexDouble &b) const {return L_ComplexDouble(re*b.re-im*b.im, re*b.im+im*b.re);}
	L_ComplexDouble operator /(const L_ComplexDouble &b) const {double m=b.re*b.re+b.im*b.im; return L_ComplexDouble((re*b.re+im*b.im)/m, (-re*b.im+im*b.re)/m);} // a/b = a*b'/|b|^2 => b.im con el signo cambiado
	L_ComplexDouble operator +(double b) const {return L_ComplexDouble(re+b, im);}
	L_ComplexDouble operator -(double b) const {return L_ComplexDouble(re-b, im);}
	L_ComplexDouble operator *(double b) const {return L_ComplexDouble(re*b, im*b);}
	L_ComplexDouble operator /(double b) const {return L_ComplexDouble(re/b, im/b);}
	L_ComplexDouble operator +() const {return *this;}
	L_ComplexDouble operator -() const {return L_ComplexDouble(-re, -im);}

};

// Slow, but nice to use
inline L_ComplexDouble operator + (double a, const L_ComplexDouble &b) {return L_ComplexDouble(a+b.re, b.im);}
inline L_ComplexDouble operator - (double a, const L_ComplexDouble &b) {return L_ComplexDouble(a-b.re, -b.im);}
inline L_ComplexDouble operator * (double a, const L_ComplexDouble &b) {return L_ComplexDouble(a*b.re, a*b.im);}
inline L_ComplexDouble operator / (double a, const L_ComplexDouble &b) {double m=b.re*b.re+b.im*b.im; return L_ComplexDouble(a*b.re/m, -a*b.im/m);} // a/b = a*b'/|b|^2

// Fractions
class L_Frac // Must have denominator>0 and must be simplified for being valid
{
public:
	long num;
	long den;
protected:
	double val;
public:
	L_Frac() {num=0; den=0; val=0;}
	L_Frac(long num, long den) {this->num=num; this->den=den; this->val=num/(double)den; if (den<0) {num*=-1;den*=-1;}}

	void setValue(long num, long den) {this->num=num; this->den=den; this->val=num/(double)den; if (den<0) {num*=-1;den*=-1;};}
	static long mcd(long num, long den);
	void simplify();
	inline void OP_add(L_Frac &other)
	{
		num=num*other.den+other.num*den;
		den*=other.den;
		simplify();
		val=num/(double)den;
		return;
	}
	inline void OP_subtract(L_Frac &other)
	{
		num=num*other.den-other.num*den;
		den*=other.den;
		simplify();
		val=num/(double)den;
		return;
	}
	inline void OP_mult(L_Frac &other)
	{
		num*=other.num;
		den*=other.den;
		simplify();
		val=num/(double)den;
		return;
	}
	inline void OP_div(L_Frac &other)
	{
		num*=other.den;
		den*=other.num;
		if (den<0)
		{
			num*=-1;
			den*=-1;
		}
		simplify();
		val=num/(double)den;
		return;
	}

	inline operator double() const {return val;}
	inline operator long() const {if (den==1) return num; else return (long)val;}
	inline L_Frac operator +=(L_Frac &other) {OP_add(other); return *this;}
	inline L_Frac operator -=(L_Frac &other) {OP_subtract(other); return *this;}
	inline L_Frac operator *=(L_Frac &other) {OP_mult(other); return *this;}
	inline L_Frac operator /=(L_Frac &other) {OP_div(other); return *this;}
	inline bool operator ==(L_Frac &other) {return (num*other.den-other.num*den)==0;}
	inline bool operator !=(L_Frac &other) {return (num*other.den-other.num*den)!=0;}
	inline bool operator >=(L_Frac &other) {return (num*other.den-other.num*den)>=0;}
	inline bool operator <=(L_Frac &other) {return (num*other.den-other.num*den)<=0;}
	inline bool operator >(L_Frac &other) {return (num*other.den-other.num*den)>0;}
	inline bool operator <(L_Frac &other) {return (num*other.den-other.num*den)<0;}
	friend L_Frac operator +(L_Frac &f1, L_Frac&f2)	{L_Frac ret(f1);ret.OP_add(f2);return ret;	}
	friend L_Frac operator -(L_Frac &f1, L_Frac&f2) {L_Frac ret(f1);ret.OP_subtract(f2);return ret; }
	friend L_Frac operator *(L_Frac &f1, L_Frac&f2) {L_Frac ret(f1);ret.OP_mult(f2);return ret; }
	friend L_Frac operator /(L_Frac &f1, L_Frac&f2) {L_Frac ret(f1);ret.OP_div(f2);return ret; }
};


// 2D and 3D allocators and deallocators

template <class T> T** L_new2d(int lx0, int ly0) //! 2D allocator
{
	int i;
    T ** im=NULL;
	try{
#ifdef L_BASIC_DEBUG
		if (lx0<=0 || ly0<=0)
			L_hard_shutdown("L_New2d : non-positive allocation size");
#endif

#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
		im=new T*[lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE]; // throws std::bad_alloc
#else
		im=new T*[lx0]; // throws std::bad_alloc
#endif
		if (im==NULL) throw std::bad_alloc();
		try
		{
			im[0]=new T[lx0*ly0];
			if (im[0]==NULL) throw std::bad_alloc();
			for (i=1; i<lx0; i++)
			{
				im[i]=im[0]+ly0*i;
			}
#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
			for (i=lx0; i<lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE; i++)
				im[i]=L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_VALUE;
#endif
		}
		catch (std::bad_alloc &)
		{
			delete[] im;
			throw;
		}
	}
	catch (std::bad_alloc &)
	{
		im=NULL;
		throw;
	}
	return im;
}

template <class T> void L_delete2d(T** im)
{
	delete[] im[0];
	im[0] = NULL;
	delete[] im;
}

template <class T> T** L_new2d_ext(int lx0, int ly0, T *buf, int step_in_chars) //! External memory embedding
{
	int i;
    T ** im=NULL;
	try{
#ifdef L_BASIC_DEBUG
		if (lx0<=0 || ly0<=0)
			L_hard_shutdown("L_New2d : non-positive allocation size");
#endif

#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
		im=new T*[lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE]; // throw std::bad_alloc
#else
		im=new T*[lx0]; // lanza std::bad_alloc
#endif
		if (im==NULL) throw std::bad_alloc();
		try
		{
			char *p = (char *)buf;
			if (p==NULL) throw std::bad_alloc();
			for (i=0; i<lx0; i++)
			{
				im[i]=(T*)((char*)p+step_in_chars*i);
			}
#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
			for (i=lx0; i<lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE; i++)
				im[i]=L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_VALUE;
#endif
		}
		catch (std::bad_alloc &)
		{
			delete[] im;
			throw;
		}
	}
	catch (std::bad_alloc &)
	{
		im=NULL;
		throw;
	}
	return im;
}

template <class T> void L_delete2d_ext(T **im) //! destroy 2D memo
{
	im[0] = NULL;
	delete[] im;
}



//#define __FORCE_ALIGN_FLOAT_DOUBLE__
#undef __FORCE_ALIGN_FLOAT_DOUBLE__
void L_test_for_aligned_malloc(size_t blo=sizeof(double[2])); // Para probar si se devuelven bloques ya alineados

// Funciones que devuelven memoria alineada, permiten usar SSE2
#ifdef __FORCE_ALIGN_FLOAT_DOUBLE__
void *L_malloc_aligned(int bytes, size_t alin);
void L_free_aligned(void *mem);
float *L_malloc_float(int n);
void L_free_float(float *mem);
double *L_malloc_double(int n);
void L_free_double(double *mem);
template <> double** L_new2d<double>(int lx0, int ly0);
template <> void L_delete2d<double>(double **im);
template <> float** L_new2d<float>(int lx0, int ly0);
template <> void L_delete2d<float>(float **im);
#endif // __FORCE_ALIGN_FLOAT_DOUBLE__


template <class T> T*** L_new3d(int lx0, int ly0, int lz0) //! 3D allocator, slow... don't use if possible
{
    int i;
	T ***im=NULL;
	try
	{
#ifdef L_BASIC_DEBUG
		if (lx0<=0 || ly0<=0 || lz0<=0)
			L_hard_shutdown("L_New3d : non-positive size allocation");
#endif

#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
		im=new T**[lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE];
#else
		im=new T**[lx0];
#endif
		if (im==NULL) throw std::bad_alloc();
		try
		{
#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
			im[0]=new T*[lx0*ly0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE];
#else
			im[0]=new T*[lx0*ly0];
#endif
			if (im[0]==NULL) throw std::bad_alloc();
			try
			{
				im[0][0]=new T[lx0*ly0*lz0];
				if (im[0][0]==NULL) throw std::bad_alloc();
				for (i=1; i<lx0*ly0; i++)
				{
					im[0][i]=&im[0][0][0]+lz0*i;
				}
#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
				for (i=lx0*ly0; i<lx0*ly0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE; i++)
					im[0][i]=L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_VALUE;
#endif
				for (i=1; i<lx0; i++)
					im[i]=&im[0][0]+i*ly0;
#ifdef L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE
				for (i=lx0; i<lx0+L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_SPACE; i++)
					im[i]=L_DOUBLE_POINTER_ALLOCATE_EXTRA_BLANK_VALUE;
#endif
			}
			catch(std::bad_alloc &)
			{
				delete[] im[0];
				throw;
			}
		}
		catch(std::bad_alloc &)
		{
			delete[] im;
			throw;
		}
	}
	catch(std::bad_alloc &)
	{
		im=NULL;
		throw;
	}
    return im;
}
template <class T> void L_delete3d(T ***im)
{
   	delete[] im[0][0];
	delete[] im[0];
	im[0] = NULL;
	delete[] im;
}


template <class T> T** L_new2dM(int li0, int lj0) //! 2D allocator, independant rows
{
	int i=0;
    T ** elem=NULL;

	try
	{
#ifdef L_BASIC_DEBUG
		if (li0<=0 || lj0<=0)
			L_hard_shutdown("L_New2d : non-positive allocation size");
#endif
		elem=new T *[li0];
		try
		{
			for (i=0; i<li0; i++)
			{
    			elem[i]=new T[lj0];
				if (elem[i]==NULL) throw;
			}
		}
		catch(std::bad_alloc &)
		{
			while(i>=0)
				delete[] elem[i--];
			delete[] elem;
			throw;
		}
	}
	catch(std::bad_alloc &)
	{
		elem=NULL;
		throw;
	}
    return elem;
}
template <class T> void L_delete2dM(T **&elem, long li) //! destroy memoria bidimensional por columnas
{
	int i;
	for (i=0; i<li; i++)
       	delete[] elem[i];
	delete[] elem;
	elem = NULL;
}

template <class T> T* L_vNewArray(int nelem, T elem1, va_list ppp)
{
	int i;
	T* ptr;
#ifdef L_BASIC_DEBUG
	if (nelem<=0)
			L_hard_shutdown("L_New2d : non-positive allocation size");
#endif
	ptr=new T[nelem];
	ptr[0]=elem1;
	for (i=1; i<nelem; i++)
		ptr[i]=va_arg(ppp,T);
	return ptr;
}

template <class T> T* L_newArreglo(int nelem, T elem1, ...)
{
	va_list ppp;
	T* ptr;
	va_start(ppp, elem1);
	ptr=L_vNewArray<T>(nelem, elem1, ppp);
	va_end(ppp);
	return ptr;
}

#if defined __TURBOC__
const L_ImageGrayDouble ** L_newArreglo(int nelem, const L_ImageGrayDouble *elem1, ...);
#endif


class L_ArgcArgv  // ° is used as argument separator in cmdList
{
public:
	int argc;
	char** argv;
	L_ArgcArgv():argc(0),argv(NULL) {}
	L_ArgcArgv(const char *str):argc(0),argv(NULL) {define(str);}
	void define(const char *cmdList);
	void ask_arguments_by_console(int *argc, char ***argv) {printf("Ingrese los argumentos del programa: ");char lin[1000]; if (L_gets(lin, 1000)==NULL) return; define(lin); *argc=this->argc; *argv=this->argv;}
	~L_ArgcArgv() {if (argv!=NULL) L_delete2d<char>(argv);}
};


/////// Template nodes+lists
template <class T> class L_MemPool;

// Basic node
template <class T>
class L_Node
{
public:
	T c;
	L_Node<T> *sig;
	L_Node() : c() {sig=NULL;}
	inline L_Node& operator =(const L_Node<T> &other) // Clones the object
	{
		c = other.c;
		sig = other.sig;
		return *this;
	}
	inline L_Node& operator =(const T &other) // Clones the object
	{
		c = other;
		return *this;
	}
	~L_Node() {}
};

template <class T>
class L_ListIterator
{
public:
	L_Node<T> *ptr;
	T &operator *() {return ptr->c;}
	void operator++() {ptr = ptr->sig;}
	void operator++(int) {ptr = ptr->sig;}
	bool operator == (const L_ListIterator &other) {return ptr == other.ptr;}
	bool operator != (const L_ListIterator &other) {return ptr != other.ptr;}
};

template <class T>
class L_List
{
private:
	L_List(const L_List<T> &other);
public:
	L_Node<T> *root;
	L_Node<T> **pult;
	int n;

	typedef L_ListIterator<T> iterator;
	iterator begin() {return reinterpret_cast<L_ListIterator<T> &>(root);}
	iterator end() {return reinterpret_cast<L_ListIterator<T> &>(*pult);}
	typedef size_t size_type;

	size_type size() const {return n;}

	L_List() {root=NULL; pult=&root;n=0;} // *pult se debe mantener = NULL  ,  pult = &(ultimo existente)->sig   o   pult=&root
	inline void recalcula_pult()                   {n=0; pult=&root; while(*pult!=NULL) {n++;pult=&(*pult)->sig;} }

	inline void push_back(const T &obj)	{*pult=new L_Node<T>; (*pult)->c=obj; pult=&(*pult)->sig; *pult=NULL; n++;}
	inline void push_front(const T &obj)  {L_Node<T> *nod=new L_Node<T>; (*nod).c=obj; nod->sig=root; root=nod; if (root->sig==NULL) pult=&root->sig; n++;}
	inline void push_back_swapping(T &obj)            {*pult=new L_Node<T>; obj.swap((**pult).c); pult=&(*pult)->sig; *pult=NULL; n++; }
	inline L_Node<T> *push_back_ret(T &obj)      {L_Node<T> *ret; *pult=new L_Node<T>; (**pult).c=obj; ret=*pult; pult=&(*pult)->sig; *pult=NULL; n++; return ret;}
	inline L_Node<T> *push_back_swapping_ret(T &obj)  { L_Node<T> *ret; *pult=new L_Node<T>; obj.swap((**pult).c); ret=*pult; pult=&(*pult)->sig; *pult=NULL; n++; return ret;}

	// Have care with these functions...
	inline void push_back_node(L_Node<T> *p)  {*pult=p; pult=&(*pult)->sig; *pult=NULL; n++;}
	inline void push_front_node(L_Node<T> *p) {p->sig=root; root=p; if (root->sig==NULL) pult=&root->sig; n++;}
	inline L_Node<T> *pop_front_ptr_ret()          {L_Node<T> *res; res=root; root=root->sig; if (root==NULL) pult = &root; n--; res.sig = NULL; return res;}

	// Reuse nodes from another list
	inline void push_back(const T &obj, L_List<T> *pool)
	{
		if (pool != NULL && pool->root != NULL)
		{
			*pult = pool->root;
			//
			pool->root=pool->root->sig;
			if (pool->root == NULL)
				pool->pult = &pool->root;
			pool->n--;
			//
			(*pult)->c = obj;
			pult=&(*pult)->sig;
			*pult = NULL;
			n++;
		}
		else
		{
			*pult=new L_Node<T>;
			(*pult)->c=obj;
			pult=&(*pult)->sig;
			*pult=NULL;
			n++;
		}
	}

	void getList(std::list<T> &li) const
	{
		L_Node<T> *nodo;
		li.clear();
		for (nodo=root; nodo!=NULL; nodo=nodo->sig);
			li.push_back(nodo->c);
	}

	void setList(const std::list<T> &li)
	{
		clear();
		typedef typename std::list<T>::const_iterator cIt;
		for (cIt iter = li.begin(); iter != li.end(); ++iter)
		{
			(*pult) = new L_Node<T>;
			(*pult)->c = *iter;
			pult = &(*pult)->sig;
		}
		*pult = NULL;
		n=li.size();
	}

	bool verifNumElem(void)
	{
		int i;
		L_Node<T> *ptr;
		for (i=0,ptr=root;ptr!=NULL;i++,ptr=ptr->sig);
		if (i!=n)
		{
			n=i;
			return false;
		}
		return true;
	}
	void sort(int delete_or_isolate_repeated, int (*qqcmp)(const void *ppT1, const void *ppT2), L_List<T> *pool = NULL);  // USA DOBLE PUNTERO, delete_or_isolate_repeated={0:no clean , 1:dejar 1 por grupo , 2: dejar 0 por grupo}.
    inline void clear()
    {
		int n0 = 0;
    	L_Node<T> *p;
		while(root!=NULL)
		{
			p=root->sig;
			delete root;
			root=p;
			n0++;
		}
		pult=&root;
		throw_L_ArgException_if(n0 != n , "L_List::clear() : numero incorrecto de nodos");
		n=0;
    }
	inline void pop_front()
	{
		L_Node<T> *p;
		if (root!=NULL)
		{
			p=root;
			root=root->sig;
			delete p;
			if (root==NULL)
				pult=&root;
			n--;
		}
	}

	inline void erase_preserving_order(bool (*funcElim)(const T &))
	{
		L_Node<T> *p;
		int n0 = n;
		int nElim = 0;
		n=0;
		for (pult = &root; *pult!=NULL; )
		{
			if ( (*funcElim)( (*pult)->c ) == true )
			{
				p = *pult;
				*pult = (*pult)->sig;
				delete p;
				nElim++;
			}
			else
			{
				pult = &(*pult)->sig;
				n++;
			}
		}
		throw_L_ArgException_if(nElim + n != n0 , "L_List::erase_preserving_order() : wrong number of nodes");
	}

	inline void erase_preserving_order(const L_Array<bool> &elim)
	{
		L_Node<T> *p;
		int n0 = n;
		int nElim = 0;
		int i;
		n=0;
		for (i=0,pult = &root; *pult!=NULL;)
		{
			if ( elim[i] == true )
			{
				p = *pult;
				*pult = (*pult)->sig;
				i++;
				delete p;
				nElim++;
			}
			else
			{
				pult = &(*pult)->sig;
				i++;
				n++;
			}
		}
		throw_L_ArgException_if(nElim + n != n0 , "L_List::erase_preserving_order() : wrong number of nodes");
	}

	inline void resize(int nMax)
	{
		int n0 = n;
		if (n <= nMax || root == NULL)
		{
			std::cout << "L_List::resize() : Case in which list must be enlarged is not implemented" << std::endl;
			return;
		}
		if (n==0)
		{
			clear();
			return;
		}
		L_List<T> tmpL;
		n=0;

		for (pult = &root; *pult!=NULL && n<nMax; n++, pult=&(*pult)->sig) {}

		throw_L_ArgException_if(*pult==NULL, "L_List::resize() : wrong number of nodes");
		tmpL.root = (*pult);
		tmpL.pult = pult;
		tmpL.n = n0-n;
		(*pult) = NULL;
		tmpL.clear();
	}

	inline void cloneListOn(L_List<T> &clone)
	{
		L_Node<T> *ptr;
		L_Node<T> nodo;
		int n0 = 0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
		{
			ptr->c.clonaObjetoEn(nodo);
			clone.push_back_swapping(nodo);
			n0++;
		}
		throw_L_ArgException_if (n0!=n, "L_List::cloneListOn() : wrong number of nodes");
	}
	inline bool copyListOn(L_List<T> &clone) const
	{
		L_Node<T> *ptr;
		int n0 = 0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
		{
			clone.push_back(ptr->c);
			n0++;
		}
		throw_L_ArgException_if (n0!=n, "L_List::copyListOn() : wrong number of nodes");
		if (n0 == n)
			return true;
		else
			{*const_cast<int *>(&n) = n0; return false;}
	}
	inline bool copyListOn(L_List<T> &clone, L_List<T> *pool) const
	{
		L_Node<T> *ptr;
		int n0 = 0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
		{
			clone.push_back(ptr->c, pool);
			n0++;
		}
		throw_L_ArgException_if (n0!=n, "L_List::copyListOn() : wrong number of nodes");
		if (n0 == n)
			return true;
		else
			{*const_cast<int *>(&n) = n0; return false;}
	}

	inline void copyListOn(L_Array<T> &clone) const // Clona nodos y los copia usando operator=()
	{
		L_Node<T> *ptr;
		int n0 = 0;
		int i;
		clone.resize(n);
		for (ptr=root,i=0; ptr!=NULL; ptr=ptr->sig,i++)
		{
			clone[i] = ptr->c;
			n0++;
		}
		throw_L_ArgException_if (n0!=n, "L_List::copyListOn() : wrong number of nodes");
	}

	inline void swap(L_List<T> &other)
	{
		L_Node<T> *raiz_t = other.root;
		L_Node<T> **pult_t = other.pult;
		int n_t = (int)other.size();
		other.root = root;
		other.pult = pult;
		other.n = n;
		root = raiz_t;
		pult = pult_t;
		n = n_t;
		if (root == NULL)
			pult = &root;
		if (other.root == NULL)
			other.pult = &other.root;
	}

	inline void moveListTo(L_List<T> &other)
	{
		if (root == NULL)
			return;
		*other.pult = root;
		other.pult = pult;
		other.n += n;
		root = NULL;
		pult = &root;
		n = 0;
	}

	~L_List() {clear();}
	L_List<T> &operator=(const L_List<T> &other)// Use cloneListOn() o copyListOn() instead of this
	{
		printf("warning: calling operator= for lists - replace by copyListOn()\n");
		other.copyListOn(*this);
		return *this;
	}
	inline void writeListTxt(FILE *fp)
	{
		L_Node<T> *ptr;
		n=0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
			n++;
		fprintf(fp,"%d\n",n);
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
			ptr->c.writeElementTxt(fp);
	}

	inline void readListTxt(FILE *fp)
	{
		int num=0;
		if (fscanf(fp,"%d",&num)!=1)
		{
			printf("Lista::readListTxt() : wrong file format\n");
			return;
		}
		while (num-- > 0)
		{
			*pult=new L_Node<T>;
			(*pult).c.readElementTxt(fp);
			pult=&(*pult)->sig;
			n++;
		}
		*pult=NULL;
	}
	inline void writeListBin(FILE *fp)
	{
		L_Node<T> *ptr;
		n=0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
			n++;
		fwrite(&n,sizeof(n),1,fp);
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
			ptr->c.writeElementBin(fp);
	}

	inline void readListBin(FILE *fp)
	{
		int num=0;
		if (fread(&num,sizeof(num),1,fp)!=1)
		{
			printf("Lista::readListBin() : wrong file format\n");
			return;
		}
		while (num-- > 0)
		{
			*pult=new L_Node<T>;
			(*pult)->c.readElementBin(fp);
			pult=&(*pult)->sig;
			n++;
		}
		*pult=NULL;
	}
	inline void readListBin(FILE *fp, ...)
	{
		va_list va;
		va_start(va, fp);
		int num=0;
		if (fread(&num,sizeof(num),1,fp)!=1)
		{
			printf("Lista::leeListaCopia() : wrong file format\n");
			return;
		}
		while (num-- > 0)
		{
			*pult=new L_Node<T>;
			(*pult)->c.readElementBin(fp, va);
			pult=&(*pult)->sig;
			n++;
		}
		*pult=NULL;
		va_end(va);
	}

	inline void writeListMem(std::vector<char> &buf)
	{
		L_Node<T> *ptr;
		int i, di;
		n=0;
		for (ptr=root; ptr!=NULL; ptr=ptr->sig)
			n++;
		i = (int)buf.size();
		buf.resize(buf.size() + sizeof(n));
		memcpy(&buf[i], &n, sizeof(n));
		i += sizeof(buf.size());
		di=root->c.writeElementMem(buf); // width de un elemento de la lista
		buf.reserve(buf.size() + n*di*2);
		if (root != NULL)
			for (ptr=root->sig; ptr!=NULL; ptr=ptr->sig)
					ptr->c.writeElementMem(buf);
	}

	inline void readListMem(const std::vector<char> &buf, int i0)
	{
		int num=0;
		int i = i0;
		memcpy(&num, &buf[i], sizeof(num));
		i += sizeof(num);
		while (num-- > 0)
		{
			*pult=new L_Node<T>;
			i+=(*pult)->c.readElementMem(buf,i);
			pult=&(*pult)->sig;
			n++;
		}
		*pult=NULL;
	}

	inline void readListMem(const std::vector<char> &buf, int i0, ...)
	{
		va_list va;
		va_start(va, i0);
		int num=0;
		int i = i0;
		memcpy(&num, &buf[i], sizeof(num));
		i += sizeof(num);
		while (num-- > 0)
		{
			*pult=new L_Node<T>;
			i+=(*pult)->c.readElementMem(buf,i,va);
			pult=&(*pult)->sig;
			n++;
		}
		*pult=NULL;
		va_end(va);
	}
};

// Sort list and isolate/delete repeated elements
//
// delete_or_isolate_repeated==0 => Only sort
// delete_or_isolate_repeated==1 => Sort and isolate repeated elements at end of list
// delete_or_isolate_repeated==2 => Sort and delete repeated elements

template <class T>
void L_List<T>::sort(int delete_or_isolate_repeated, int (*qqCmp)(const void *ppT1, const void *ppT2), L_List<T> *pool) // int qqCmp(T **ppT1, T **ppT2)
{
	typedef L_Node<T> * L_NPtr;
#ifdef L_BASIC_DEBUG
	if (n<=0)
		L_hard_shutdown("L_List<T>::sort : non-positive size allocation");
#endif
	L_NPtr *arr=new L_NPtr[n];
	bool *die=new bool[n];
	int i,n0;
	L_Node<T> *ptr;
	ptr=root;
	for (i=0; i<n; i++)
	{
		arr[i]=ptr;
		die[i]=false;
		ptr=ptr->sig;
		arr[i]->sig=NULL;
	}
	qsort(arr,n,sizeof(L_NPtr),qqCmp);
	if (delete_or_isolate_repeated==2)
	{
		for (i=0; i<n-1; i++)
		{
			if (qqCmp(&arr[i],&arr[i+1])==0) // Del grupo de iguales, no dejar ninguno
			{
				die[i]=true;
				die[i+1]=true;
			}
		}
	}
	else if (delete_or_isolate_repeated==1)
	{
		for (i=0; i<n-1; i++)
		{
			if (qqCmp(&arr[i],&arr[i+1])==0) // Del grupo de iguales, dejar 1
			{
				die[i+1]=true;
			}
		}
	}
	n0=n;
	n=0;
	root=NULL;
	pult=&root;
	if (pool == NULL)
	{
		for (i=0; i<n0; i++)
		{
			if (die[i]==false)
			{
				*pult=arr[i];
				n++;
				pult=&(*pult)->sig;
				*pult=NULL;
			}
			else
				delete arr[i];
		}
	}
	else
	{
		for (i=0; i<n0; i++)
		{
			if (die[i]==false)
			{
				*pult=arr[i];
				n++;
				pult=&(*pult)->sig;
				*pult=NULL;
			}
			else
				pool->push_back_node(arr[i]);
		}
	}
	delete[] arr;
	delete[] die;
	return;
};

// List-related classes

template <class T>
class L_NodePtr
{
public:
	T *c;
	L_NodePtr<T> *sig;
	L_NodePtr() {c=NULL;sig=NULL;}
	~L_NodePtr() {}
};

template <class T>
class L_ListPtr
{
public:
	L_NodePtr<T> *root;
	L_NodePtr<T> **pult;
	int n;

	L_ListPtr() {root=NULL; pult=&root;n=0;}

	size_t size() const {return n;}

	inline bool push_back_ptr(T *obj)
	{
		if (obj==NULL)
			return false;
		*pult=new L_NodePtr<T>;
		(*pult)->c=obj;
		pult=&(*pult)->sig;
		*pult=NULL;
		n++;
		return true;
	}
	inline bool push_back_swapping(T& obj)
	{
		return push_back_ptr(obj.clonaRegalaCalces());
	}
	inline bool push_back(T& obj)
	{
		return push_back_ptr(obj.clona());
	}
	inline void swap(L_ListPtr<T> &other)
	{
		L_NodePtr<T> *raiz_t = other.root;
		L_NodePtr<T> **pult_t = other.pult;
		int n_t = (int)other.size();
		other.root = root;
		other.pult = pult;
		other.n = n;
		root = raiz_t;
		pult = pult_t;
		n = n_t;
		if (root == NULL)
			pult = &root;
		if (other.root == NULL)
			other.pult = &other.root;
	}
	inline void moveListTo(L_ListPtr<T> &other)
	{
		if (root == NULL)
			return;
		*other.pult = root;
		other.pult = pult;
		other.size() += n;
		root = NULL;
		pult = &root;
		n = 0;
	}

	bool verifNumElem(void)
	{
		int i;
		L_NodePtr<T> *ptr;
		for (i=0,ptr=root;ptr!=NULL;i++,ptr=ptr->sig);
		if (i!=n)
		{
			n=i;
			return false;
		}
		return true;
	}
    inline void clear()
    {
    	L_NodePtr<T> *p;
		while(root!=NULL)
		{
			p=root->sig;
			delete root->c;
			delete root;
			root=p;
		}
		pult=&root;
		n=0;
    }
	inline bool empty() {return root==NULL;}
	~L_ListPtr() {clear();}
};

class L_CountingTreeNode
{
public:
	int nobj;
	int nvot;
	L_CountingTreeNode *right;
	L_CountingTreeNode *left;
	L_CountingTreeNode() {right=NULL; left=NULL; nvot=0;}
};

class L_CountingTree
{
protected:
	L_CountingTreeNode *root;
public:
	int nVotTot;
	int nCellsTot;

protected:
	void destroyRec(L_CountingTreeNode *ptr)
	{
		if (ptr==NULL)
			return;
		destroyRec(ptr->right);
		destroyRec(ptr->left);
		delete ptr;
	}
public:
	L_CountingTree() {root=NULL;nVotTot=0;nCellsTot=0;}
	void addVoteToFrequencies(int nCelda, int votos=1);
	int votosDe(int nCelda);
	void swap(L_CountingTree &other)
	{
		L_CountingTreeNode *raiz_t = other.root;
		int nVotTot_t = other.nVotTot;
		int nCeldasTot_t = other.nCellsTot;

		other.root=root;
		other.nVotTot=nVotTot;
		other.nCellsTot=nCellsTot;

		root=raiz_t;
		nVotTot=nVotTot_t;
		nCellsTot=nCeldasTot_t;
	}
	void destroy()	{destroyRec(root); root=NULL; nVotTot=0; nCellsTot=0;}
	~L_CountingTree() { destroy(); }
};

class L_HuffmanTreeNode
{
public:
	// Para estructura de tree
	L_HuffmanTreeNode *padre;
	L_HuffmanTreeNode *hijoIzq;
	L_HuffmanTreeNode *hijoDer;
	// Para estructura de lista de arboles
	L_HuffmanTreeNode *sig;
	long peso;
	unsigned char value;
	int prof;
	L_HuffmanTreeNode() {padre=hijoIzq=hijoDer=sig=NULL; prof=0;}
	void destroyRec() {  if (hijoIzq!=NULL) {hijoIzq->destroyRec(); delete hijoIzq;} if (hijoDer!=NULL) {hijoDer->destroyRec(); delete hijoDer;}  }
	int cuentaNodos(int nmax=1000) {if (this==NULL) return 0; if (nmax<=0) {printf("Overflow\n"); return 1;} if (sig!=NULL) return 1+sig->cuentaNodos(nmax-1); else return 1;}
	static int cmp(const void *p1, const void *p2);
	static int cmp2(const void *p1, const void *p2);
	int calcProf() {if (this==NULL) return 0; int a=1+hijoIzq->calcProf(); int b=1+hijoDer->calcProf(); if (a>b) return a; else return b;}
};

#define L_HuffmanTreeList_pop_front_ptr(ret_min,lis) {ret_min=(lis).root; (lis).root=(lis).root->sig; if ((lis).root==NULL) (lis).pult = &(lis).root; ret_min->sig = NULL; (lis).n--;}
#define L_HuffmanTreeList_pop_minor_element(ret_min,lis1,lis2) {if (lis1.root != NULL) { if (lis2.root != NULL) { if (lis1.root->peso <= lis2.root->peso) {L_HuffmanTreeList_pop_front_ptr(ret_min,lis1);} else {L_HuffmanTreeList_pop_front_ptr(ret_min,lis2);} } else {L_HuffmanTreeList_pop_front_ptr(ret_min,lis1);} } else {L_HuffmanTreeList_pop_front_ptr(ret_min,lis2);} }

class L_HuffmanTreeList
{
public:
	L_HuffmanTreeNode *root;
	L_HuffmanTreeNode **pult;
	int n;

	L_HuffmanTreeList() {root=NULL; pult=&root; n=0;} // *pult must be = NULL  ,  pult = &(last)->sig   or   pult=&root

	inline void push_back_node(L_HuffmanTreeNode *p)  {*pult=p; pult=&(*pult)->sig; *pult=NULL; n++;}
	inline L_HuffmanTreeNode *pop_front_ptr_ret()          {L_HuffmanTreeNode *res; res=root; root=root->sig; if (root==NULL) pult = &root; res->sig = NULL; n--; return res;}

	static L_HuffmanTreeNode *pop_minor_element(L_HuffmanTreeList &lis1, L_HuffmanTreeList &lis2)
	{
		L_HuffmanTreeNode *ret;
		if (lis1.root != NULL)
		{
			if (lis2.root != NULL)
			{
				if (lis1.root->peso <= lis2.root->peso)
					ret = lis1.pop_front_ptr_ret();
				else
					ret = lis2.pop_front_ptr_ret();
			}
			else
				ret = lis1.pop_front_ptr_ret(); // lis1!=NULL, lis2==NULL
		}
		else
			ret = lis2.pop_front_ptr_ret();  // lis1==NULL  lis2?
		return ret;
	}

	~L_HuffmanTreeList() {L_HuffmanTreeNode *p; while (root!=NULL) {p=root; root=root->sig; delete p;} }
};


class L_HuffmanEncoder;


#define new_L_HuffmanTreeNode() (pA[iult].padre=pA[iult].hijoIzq=pA[iult].hijoDer=pA[iult].sig=NULL, pA[iult].prof=0, &pA[iult++])
#define new_L_HuffmanTreeNode2() (tree.pA[tree.iult].padre=tree.pA[tree.iult].hijoIzq=tree.pA[tree.iult].hijoDer=tree.pA[tree.iult].sig=NULL, tree.pA[tree.iult].prof=0, &tree.pA[tree.iult++])


// Bit-wise functions
#define L_CH_getBit(buf,bitPos) (  (((long)(buf)[(bitPos)/8]) >> ((bitPos)%8)) & 1L   )
#define L_CH_Bit_0(buf,bitPos) ( (buf)[(bitPos)/8] &= ~(1L<<((bitPos)%8)) ) // Solo apaga el bit
#define L_CH_Bit_1(buf,bitPos) ( (buf)[(bitPos)/8] |= (1L<<((bitPos)%8)) ) // Solo enciende el bit
#define L_CH_setBit(buf,bitPos,val_bit) (L_CH_Bit_0(buf, bitPos), val_bit && L_CH_Bit_1(buf, bitPos))

class L_HuffmanTree // Huffman tree for compressing sequences of bytes
{
public:
	long pesos[256]; // Frequencies for characters. Index is unsigned char.
	L_HuffmanTreeNode *root; // tree
	L_HuffmanTreeNode *children[256]; // Pointers to leaves
	L_HuffmanTreeNode pA[256 * 4]; // Nodes of the tree preallocated here
	int iult;

	L_HuffmanTree() {  root=NULL; for (int i=0; i<=255; i++) {pesos[i]=0; children[i]=NULL;} iult = 0;}
	inline void destroyTree() {  if (root!=NULL && iult == 0) {root->destroyRec(); delete root;} root=NULL;  iult = 0;}
	inline void clear() {  destroyTree(); for (int i=0; i<=255; i++) {pesos[i]=0; children[i]=NULL;}  }
	inline void addVoteToFrequencies(unsigned char elem, int num=1) { pesos[elem]+=num; }
	inline void addVotesToFrequencies(const char *fuente, long l) { for (long i=0; i<l; i++) pesos[(unsigned char)fuente[i]]++; }
	void findPathInTree(unsigned char c, long &ruta, unsigned char &nBits);  // build la codificacion para el caracter c

	void build_tree_using_stored_frequeencies();
	bool search_for_incoherences_in_tree(const L_HuffmanTreeNode *inic=NULL) const;

	~L_HuffmanTree() { destroyTree(); }
};

class L_HuffmanPath
{
public:
	long trainOfBits; 
	unsigned char numBits;
	unsigned char value;

	L_String getString() {int j; L_String ret; ret.resize(numBits+1); for (j=0; j<numBits; j++) ret[numBits-j-1] = ((trainOfBits>>j)&1 ) ? '1' : '0'; ret[numBits] = 0; return ret;}

	static int canonicalOrdering(const void *huffRuta1, const void *huffRuta2);
	static int lexicalOrdering(const void *huffRuta1, const void *huffRuta2);
};


class L_HuffmanEncoder
{
public:
	typedef int size_type; // Could be increased for some applications
	// Caracter c encoded as a train of bits
	L_HuffmanPath codeMapping[256]; // Path from the leave to the root: (>>0)&1, (>>1)&1, ...
	L_HuffmanTree tree;
	bool hasCanonicalOrdering;

private:
	void buildPathsFromTree();
	bool buildTreeFromPaths();
	bool check_for_incoherences_beetween_tree_and_paths();
	void buildCanonicalOrdering(); // using only codeMapping[].numBits

	int encode(const char *fuente, long len, char *dest); // Requires both paths[] and rutasNumBits[] built, DANGEROUS
	int encode0(const char *fuente, long len, char *dest); // Requires both paths[] and rutasNumBits[] built

	static void testTree();
	static void testBits();
	static void testReadWriteHuffman();
public:
	std::vector<char> res; // Result from calling encode() or decode()

	// Funciones para read y guardar las paths
	void writeHuffmanTxt(FILE *fp);
	void readHuffmanTxt(FILE *fp);
	void writeHuffmanBin(FILE *fp);
	int readHuffmanBin(FILE *fp);
	int lengthOfHuffmanBin();
	void writeHuffman(std::vector<char> &buffer);
	int readHuffman(const std::vector<char> &buffer);
	void writeHuffmanCanonical(std::vector<char> &buffer);
	int readHuffmanCanonical(const std::vector<char> &buffer);

	void clear() {tree.clear();}

	//These functions require paths yet built
	void build(const char *fuente, long len);
	void encode(const char *fuente, long len); // Requires paths[] and rutasNumBits[]
	bool decode(const char *fuente, long len, long lenFinal); // Requires the tree

	///////////////////////
	// Simple interface
	///////////////////////

	static void encodeAll(const std::vector<char> &input, std::vector<char> &output);
	static bool decodeAll(const std::vector<char> &input, std::vector<char> &output);

	static double entropy_in_bits(const std::vector<char> &input);

	// For RLE and delta encoding before Huffman
	static void encodeRLE0(const std::vector<char> &input, std::vector<char> &salidaRLE);
	static void decodeRLE0(const std::vector<char> &entradaRLE, std::vector<char> &output);

	static void encodeDelta(const std::vector<char> &input, std::vector<char> &salidaD);
	static void decodeDelta(const std::vector<char> &entradaD, std::vector<char> &output);

	static void testEncoding();
	static void test(bool soloCompr = false);

	static void testEncoderRLE();
};


template <class T>
class L_Arr2D_POD
{
private:
	L_Arr2D_POD &operator =(const L_Arr2D_POD &other) {return *this;}
	L_Arr2D_POD(const L_Arr2D_POD<T> &other) {li=other.li; lj=other.lj; elem=L_new2d<T>(li,lj); memcpy(*elem, *other.elem, li*lj*sizeof(T));}
public:
	T **elem;
	int li, lj;
	L_Arr2D_POD(int li, int lj) {this->li=li; this->lj=lj; elem=L_new2d<T>(li, lj);}
	~L_Arr2D_POD() {if (elem!=NULL) L_delete2d<T>(elem);}
};

// Only for POD obejcts, it uses memcpy()
template <class T>
class L_Arr3D_POD
{
private:
	L_Arr3D_POD &operator =(const L_Arr3D_POD &other) {return *this;}
	L_Arr3D_POD(const L_Arr3D_POD<T> &other) {li=other.li; lj=other.lj; lk=other.lk; elem=L_new3d<T>(li,lj, lk); memcpy(**elem, **other.elem, li*lj*lk*sizeof(T));}
public:
	T ***elem;
	int li, lj, lk;
	L_Arr3D_POD(int li, int lj, int lk) {this->li=li; this->lj=lj; this->lk=lk; elem=L_new3d<T>(li, lj, lk);}
	~L_Arr3D_POD() {if (elem!=NULL) L_delete3d<T>(elem);}
};



/////// Clases

#if defined(L_FUNCTION_CALL_SIZE)
	#error #define conflicts
#endif
#define L_FUNCTION_CALL_SIZE 100
class L_FunctionCall
{
public:
	char name[L_FUNCTION_CALL_SIZE];
	L_FunctionCall & operator = (const L_FunctionCall &other) {memcpy(name, other.name, L_FUNCTION_CALL_SIZE); return *this;}
};

typedef L_Node<L_FunctionCall> L_FunctionCallNode;
typedef L_List<L_FunctionCall> L_FunctionCallList;


////////////////
// Parameter reading and writing
////////////////

// tipos de variables
enum L_DataType {L_ty_undefined, L_ty_bool, L_ty_char, L_ty_uchar, L_ty_schar, L_ty_short, L_ty_int, L_ty_long, L_ty_ushort, L_ty_uint, L_ty_ulong, L_ty_double, L_ty_float, L_ty_string};


#ifdef L_USE_NEW_PARAMETER_SERVER
union L_MultiTypeBasic
{
	bool b;
	char c;
	unsigned char uc;
	signed char sc;
	short s;
	int i;
	long l;
	unsigned short us;
	unsigned int ui;
	unsigned long ul;
	double d;
	float f;
};

#define L_REFDATA_SEPARATOR '\"' // Caracter para aislar cada string en el archivo de parametros

class L_RefData
{
private:
	// Cuando el tipo es L_ty_undefined, ptr siempre apunta a temporal (como double) o a strTemporal (como string)
	// Cuando el tipo es other, ptr puede apuntar a una variable externa, a temporal (como multitipo) o a strTemporal (como string)
	//
	L_MultiTypeBasic temporal; // Almacenamiento temporal para numeros, debe poder contener cualquier tipo basico
	L_String strTemporal;  // Almacenamiento temporal para strings
	L_DataType tipo;
	void *ptr; // Puede apuntar a una variable externa, a temporal o a strTemporal
public:

	static bool compat(L_DataType a, L_DataType b) {if (a==L_ty_undefined || b==L_ty_undefined || a==b) return true; return false;}
	bool coherent() const {return  (tipo==L_ty_undefined) <= (ptr==&temporal||ptr==&strTemporal);} // menor_igual_logico = implicancia
	bool hasString() const {return (tipo==L_ty_undefined && ptr==&strTemporal) ||  tipo==L_ty_string;}
	bool hasNum() const {return (tipo==L_ty_undefined && ptr==&temporal) ||  (tipo!=L_ty_undefined && tipo!=L_ty_string);}
	bool coherentWith(const L_RefData &other) const {return coherent() && other.coherent() && (tipo==L_ty_undefined || hasString()==other.hasString());}
	bool useTemporal() {return ptr==&temporal || ptr==&strTemporal;}

	L_RefData() {tipo=L_ty_undefined; ptr=&temporal; temporal.d=0;} // Para que parta inicializado
	template <class T> L_RefData(T *var) {tipo=typeOfData(var); ptr=var;}

	void copia(const L_RefData &other);

	L_RefData(const L_RefData &other) {tipo=L_ty_undefined; ptr=&temporal; copia(other);}
	const L_RefData& operator=(const L_RefData &other) {throw_L_ArgException_if(!coherentWith(other)); copia(other); return *this;}

	inline static L_DataType typeOfData(const bool *b) {return L_ty_bool;}

	inline static L_DataType typeOfData(const char *b) {return L_ty_char;}
	inline static L_DataType typeOfData(const unsigned char *b) {return L_ty_uchar;}
	inline static L_DataType typeOfData(const signed char *b) {return L_ty_schar;}

	inline static L_DataType typeOfData(const short *b) {return L_ty_short;}
	inline static L_DataType typeOfData(const int *b) {return L_ty_int;}
	inline static L_DataType typeOfData(const long *b) {return L_ty_long;}

	inline static L_DataType typeOfData(const unsigned short *b) {return L_ty_ushort;}
	inline static L_DataType typeOfData(const unsigned int *b) {return L_ty_uint;}
	inline static L_DataType typeOfData(const unsigned long *b) {return L_ty_ulong;}

	inline static L_DataType typeOfData(const double *b) {return L_ty_double;}
	inline static L_DataType typeOfData(const float *b) {return L_ty_float;}

	inline static L_DataType typeOfData(const L_String *b) {return L_ty_string;}

	inline void copyValueOn(void *var) const // Copia usando directamente el tipo de data especifico
	{
		throw_L_DimException_if(ptr==NULL, "L_RefData::copyValueOn");
		switch(tipo)
		{
		case L_ty_undefined: throw_L_ArgException_if(true,"L_RefData::copyValueOn");
		case L_ty_bool: *(bool *)var=*(bool *)ptr; break;
		case L_ty_char: *(char *)var=*(char *)ptr; break;
		case L_ty_uchar: *(unsigned char *)var=*(unsigned char *)ptr; break;
		case L_ty_schar: *(signed char *)var=*(signed char *)ptr; break;
		case L_ty_short: *(short *)var=*(short *)ptr; break;
		case L_ty_int: *(int *)var=*(int *)ptr; break;
		case L_ty_long: *(long *)var=*(long *)ptr; break;
		case L_ty_ushort: *(unsigned short *)var=*(unsigned short *)ptr; break;
		case L_ty_uint: *(unsigned int *)var=*(unsigned int *)ptr; break;
		case L_ty_ulong: *(unsigned long *)var=*(unsigned long *)ptr; break;
		case L_ty_double: *(double *)var=*(double *)ptr; break;
		case L_ty_float:  *(float *)var=*(float *)ptr; break;
		case L_ty_string: *(L_String *)var=*(L_String *)ptr; break;
		default: throw_L_ArgException_if(true,"L_RefData::copyValueOn");
		}
	}

	inline void setValue(double val) // Para copia usando double
	{
		throw_L_ArgException_if(ptr==NULL, "L_RefData::setValue");
		if (ptr==&strTemporal)
			ptr=&temporal;
		switch(tipo)
		{
		case L_ty_undefined: throw_L_ArgException_if(ptr!=&temporal, "L_RefData::setValue"); temporal.d=val; break;
		case L_ty_bool: *(bool *)ptr=(val>0.1 || val<=-0.1)?1:0; break;
		case L_ty_char: *(char *)ptr=(char)val; break;
		case L_ty_uchar: *(unsigned char *)ptr=(unsigned char)val; break;
		case L_ty_schar: *(signed char *)ptr=(signed char)val; break;
		case L_ty_short: *(short *)ptr=(short)val; break;
		case L_ty_int: *(int *)ptr=(int)val; break;
		case L_ty_long: *(long *)ptr=(long)val; break;
		case L_ty_ushort: *(unsigned short *)ptr=(unsigned short)val; break;
		case L_ty_uint: *(unsigned int *)ptr=(unsigned int)val; break;
		case L_ty_ulong: *(unsigned long *)ptr=(unsigned long)val; break;
		case L_ty_double: *(double *)ptr=(double)val; break;
		case L_ty_float:  *(float *)ptr=(float)val; break;
		default: throw_L_ArgException_if(ptr!=&temporal, "L_RefData::setValue");
		}
	}

	inline void setValue(L_String &val) // Para copia usando L_String
	{
		throw_L_ArgException_if(ptr==NULL, "L_RefData::setValue");
		if (ptr==&temporal)
			ptr=&strTemporal;
		switch(tipo)
		{
		case L_ty_undefined: throw_L_ArgException_if(ptr!=&strTemporal,"L_RefData::setValue"); strTemporal=val; break;
		case L_ty_string: *(L_String *)ptr=val; break;
		default: throw_L_ArgException_if(true,"L_RefData::setValue");
		}
	}

	inline double getValueDouble() const // Para copia usando double
	{
		double ret;
		throw_L_ArgException(ptr==NULL || ptr==&strTemporal,"L_RefData::getValueDouble");
		switch(tipo)
		{
		case L_ty_undefined: ret=temporal.d; break;
		case L_ty_bool: ret=(double)*(bool *)ptr; break;
		case L_ty_char: ret=(double)*(char *)ptr; break;
		case L_ty_uchar: ret=(double)*(unsigned char *)ptr; break;
		case L_ty_schar: ret=(double)*(signed char *)ptr; break;
		case L_ty_short: ret=(double)*(short *)ptr; break;
		case L_ty_int: ret=(double)*(int *)ptr; break;
		case L_ty_long: ret=(double)*(long *)ptr; break;
		case L_ty_ushort: ret=(double)*(unsigned short *)ptr; break;
		case L_ty_uint: ret=(double)*(unsigned int *)ptr; break;
		case L_ty_ulong: ret=(double)*(unsigned long *)ptr; break;
		case L_ty_double: ret=(double)*(double *)ptr; break;
		case L_ty_float:  ret=(double)*(float *)ptr; break;
		default: throw_L_ArgException(true,"L_RefData::getValueDouble");
		}
		return ret;
	}

	inline L_String getValueString() const // Para copia usando L_String
	{
		L_String ret;
		throw_L_ArgException_if(ptr==NULL || ptr==&temporal, "L_RefData::getValueString");
		switch(tipo)
		{
		case L_ty_undefined: ret=strTemporal; break;
		case L_ty_string: ret=*(L_String *)ptr; break;
		default: throw_L_ArgException(true,"L_RefData::getValueString");
		}
		return ret;
	}

	inline char *allToText(char *buf)
	{
		throw_L_ArgException_if(ptr==NULL, "L_RefData::allToText");
		switch(tipo)
		{
		case L_ty_undefined: if (ptr==&temporal) sprintf(buf,"(indef)%g",temporal.d); else sprintf(buf,"(indef)%c%s%c",L_REFDATA_SEPARATOR,strTemporal.str,L_REFDATA_SEPARATOR); break;
		case L_ty_bool: if (*(bool*)ptr==true) strcpy(buf,"(bool)true"); else strcpy(buf,"(bool)false"); break;
		case L_ty_char: sprintf(buf,"(char)%d",*(char*)ptr); break;
		case L_ty_uchar: sprintf(buf,"(unsigned char)%d",*(unsigned char*)ptr); break;
		case L_ty_schar: sprintf(buf,"(signed char)%d",*(signed char*)ptr); break;
		case L_ty_short: sprintf(buf,"(short)%d",*(short*)ptr); break;
		case L_ty_int: sprintf(buf,"(int)%d",*(int*)ptr); break;
		case L_ty_long: sprintf(buf,"(long)%ld",*(long*)ptr); break;
		case L_ty_ushort: sprintf(buf,"(unsigned short)%d",*(unsigned short*)ptr); break;
		case L_ty_uint: sprintf(buf,"(unsigned int)%d",*(unsigned int*)ptr); break;
		case L_ty_ulong: sprintf(buf,"(unsigned long)%ld",*(unsigned long*)ptr); break;
		case L_ty_double: sprintf(buf,"(double)%g",*(double*)ptr); break;
		case L_ty_float:  sprintf(buf,"(float)%g",*(float*)ptr); break;
		case L_ty_string: sprintf(buf,"(L_String)%c%s%c", L_REFDATA_SEPARATOR, strTemporal.str, L_REFDATA_SEPARATOR); break;
		default: throw_L_ArgException(true, "L_RefData::allToText");
		}
		return buf;
	}

	bool readValuesFromLine(const char *buf);

	template <class T> inline void set(T *var) {tipo=typeOfData(var); ptr=var;} // Para tipo basico (no objetos)

	void eraseType() {ptr=&temporal; tipo=L_ty_undefined;}
	void setType(L_DataType tipo) {throw_L_ArgException(this->tipo!=L_ty_undefined,"L_RefData::setType"); this->tipo=tipo; if (tipo==L_ty_string && ptr==&temporal) ptr=&strTemporal; else if (tipo!=L_ty_undefined && ptr==&strTemporal) ptr=&temporal;}
};

// parametro etiquetado
class L_ParamLabel
{
private:
	L_RefData data;
	char *name;
	char *class_name;

public:
	L_ParamLabel() {name=NULL; class_name=NULL;}
	L_ParamLabel(const L_ParamLabel &other) {set(other.name,other.data); setClass(other.class_name);}
	L_ParamLabel &operator=(const L_ParamLabel &other) {set(other.name,other.data); setClass(other.class_name); return *this;}
	void setName(const char *name);
	void setClass(const char *class_name);
	void setFromString(const char *value);
	template <class T> inline void set(const char *name, T *ptr) {setName(name); data.set(ptr);}
	inline void set(const char *name, L_RefData otherDato) {setName(name); data=otherDato;}
	void swap(L_ParamLabel &other) {L_ParamLabel t; t.data = other.data; t.name = other.name; t.class_name = other.class_name; other.data=data; other.name=name; other.class_name=class_name; data = t.data; name=t.name; class_name=t.class_name;}
	~L_ParamLabel() {if (name!=NULL) delete[] name; if (class_name!=NULL) delete[] class_name;}
	static int cmpInv(const void *a, const void *b);
	static int cmpNorm(const void *a, const void *b);
	friend class L_ParamBlock;
	friend class L_ParamLabelList;
};

typedef L_Node<L_ParamLabel> L_ParamLabelNode;

class L_ParamLabelList:public L_List<L_ParamLabel>
{
public:
	bool readFile(FILE *fp, bool cerrarArchivo=true);
	bool saveFile(FILE *fp, bool cerrarArchivo=true);
	bool readFile(const char *arch) {FILE *fp; if ( (fp=fopen(arch,"r"))!=NULL) return readFile(fp, true); return false;}
	bool saveFile(const char *arch) {FILE *fp; if ( (fp=fopen(arch,"w"))!=NULL) return saveFile(fp, true); return false;}
};

class L_ParamManagerLocal;

// Set of labeled parameters
class L_ParamBlock
{
private:
	L_ParamLabel *pa;
	int nPaTot; // Number of elements stored
	int nPaMem; // Number of elements allocated
	bool ordered;
	char *class_name;

	// Internal functions
	L_ParamBlock() {} // Must no exist
	void addFrom(const char *name, void *ptr) {} // Must no exist
	L_ParamBlock(const L_ParamBlock &other) {} // Must no exist
	L_ParamBlock &operator=(const L_ParamBlock &other) {return *this;} // Must no exist
	void copyOn_naive(L_ParamBlock &other) {other.pa=pa; other.nPaTot=nPaTot; other.nPaMem=nPaMem; other.ordered=ordered; other.class_name=class_name;}
	void undefine_naive() {pa=NULL; class_name=NULL;}
	// Funciones de uso interno mas generales
	inline void sort() {qsort(pa, nPaTot, sizeof(L_ParamLabel), &L_ParamLabel::cmpNorm); ordered=true;}
	L_RefData *findVariable(const char *name);
	void revMemo();
	bool belongTo(const L_ParamLabel &param) {return strcmp(class_name, param.class_name)==0;}
	bool updateValue(const L_ParamLabel &other);

public:
	L_ParamBlock(const char *class_name)
	{
        class_name=new char[strlen(class_name)+1];
		strcpy(class_name,class_name);
		pa=NULL;
		nPaTot=0;
		nPaMem=0;
		ordered=true;
		revMemo();
	}
	void swap(L_ParamBlock &other) {L_ParamBlock tmp; copyOn_naive(tmp); other.copyOn_naive(*this); tmp.copyOn_naive(other); tmp.undefine_naive();}
	template <class T> inline void addFrom(const char *name, T *ptr) {revMemo(); pa[nPaTot].set(name, ptr); pa[nPaTot].setClass(class_name); ordered=false; nPaTot++;}
	void updateAssociatedValues(L_ParamLabelList &lista, bool borrarValoresUsados);
	void buildListFromAssociatedValues(L_ParamLabelList &lista);

	~L_ParamBlock() {if (pa!=NULL) delete[] pa; if (class_name!=NULL) delete[] class_name;}
};

class L_ParamManagerLocal_POD
{
protected:
	L_ParamManagerLocal **children;
	int nChildrenTot;
	int nChildrenMem;
public:
	bool modifiedFlowOfInfo;
};

// Local parameter manager with children
class L_ParamManagerLocal : public L_ParamManagerLocal_POD
{
private:
	L_ParamBlock params;

public:
	L_ParamManagerLocal(const char *class_name, int nChildrenMem):params(class_name)
	{
		nChildrenTot=0;
		modifiedFlowOfInfo=true;
		this->nChildrenMem=nChildrenMem;
		if (nChildrenMem>0)
			children=new L_ParamManagerLocal*[nChildrenMem];
		else
			children=NULL;
	}
	void swap(L_ParamManagerLocal &other);
	bool addChildren(L_ParamManagerLocal *nuevoHijo);
	void updateValues_me_and_children(L_ParamLabelList &lista);
	void buildListOfValues_me_and_children(L_ParamLabelList &lista);
	void setValue_modifiedFlowOfInfo(bool modified) {modifiedFlowOfInfo=modified; for (int i=0; i<nChildrenTot; i++) children[i]->setValue_modifiedFlowOfInfo(alterar);}
	template <class T> inline void addFrom(const char *name, T *ptr) {params.addFrom(name, ptr);}
	~L_ParamManagerLocal() {if (children!=NULL) delete[] children;}
};

#else
class L_RefData
{
private:
	double temporal;
public:
	L_DataType tipo;
	void *ptr;

	L_RefData() {tipo=L_ty_undefined; ptr=&temporal;} // Para que parta inicializado
	template <class T> L_RefData(T *var) {tipo=typeOfData(var); ptr=var;}
	L_RefData(const L_RefData &other) {tipo=other.tipo; if (tipo==L_ty_undefined) {temporal=other.temporal; ptr=&temporal;} else ptr=other.ptr;}
	const L_RefData& operator=(const L_RefData &other) {tipo=other.tipo; if (tipo==L_ty_undefined) {temporal=other.temporal; ptr=&temporal;} else ptr=other.ptr; return *this;}

	inline static L_DataType typeOfData(const bool *b) {return L_ty_bool;}

	inline static L_DataType typeOfData(const char *b) {return L_ty_char;}
	inline static L_DataType typeOfData(const unsigned char *b) {return L_ty_uchar;}
	inline static L_DataType typeOfData(const signed char *b) {return L_ty_schar;}

	inline static L_DataType typeOfData(const short *b) {return L_ty_short;}
	inline static L_DataType typeOfData(const int *b) {return L_ty_int;}
	inline static L_DataType typeOfData(const long *b) {return L_ty_long;}

	inline static L_DataType typeOfData(const unsigned short *b) {return L_ty_ushort;}
	inline static L_DataType typeOfData(const unsigned int *b) {return L_ty_uint;}
	inline static L_DataType typeOfData(const unsigned long *b) {return L_ty_ulong;}

	inline static L_DataType typeOfData(const double *b) {return L_ty_double;}
	inline static L_DataType typeOfData(const float *b) {return L_ty_float;}

	inline static int width(L_DataType t)
	{
		switch(t)
		{
		case L_ty_undefined: return sizeof(double);
		case L_ty_bool: return sizeof(bool);
		case L_ty_char: return sizeof(char);
		case L_ty_uchar: return sizeof(unsigned char);
		case L_ty_schar: return sizeof(signed char);
		case L_ty_short: return sizeof(short);
		case L_ty_int: return sizeof(int);
		case L_ty_long: return sizeof(long);
		case L_ty_ushort: return sizeof(unsigned short);
		case L_ty_uint: return sizeof(unsigned int);
		case L_ty_ulong: return sizeof(unsigned long);
		case L_ty_double: return sizeof(double);
		case L_ty_float: return sizeof(float);
		default: throw_L_ArgException_if(true,"L_DataType::width");
		}
 		return 100;
	}
	inline void copyValueOn(void *var) const
	{
		throw_L_ArgException_if(ptr==NULL, "L_DataType::copyValueOn");
		switch(tipo)
		{
		case L_ty_undefined: throw_L_ArgException_if(true,"L_DataType::copyValueOn");
		case L_ty_bool: *(bool *)var=*(bool *)ptr; break;
		case L_ty_char: *(char *)var=*(char *)ptr; break;
		case L_ty_uchar: *(unsigned char *)var=*(unsigned char *)ptr; break;
		case L_ty_schar: *(signed char *)var=*(signed char *)ptr; break;
		case L_ty_short: *(short *)var=*(short *)ptr; break;
		case L_ty_int: *(int *)var=*(int *)ptr; break;
		case L_ty_long: *(long *)var=*(long *)ptr; break;
		case L_ty_ushort: *(unsigned short *)var=*(unsigned short *)ptr; break;
		case L_ty_uint: *(unsigned int *)var=*(unsigned int *)ptr; break;
		case L_ty_ulong: *(unsigned long *)var=*(unsigned long *)ptr; break;
		case L_ty_double: *(double *)var=*(double *)ptr; break;
		case L_ty_float:  *(float *)var=*(float *)ptr; break;
		default: throw_L_ArgException_if(true,"L_DataType::copyValueOn");
		}
	}

	inline void setValue(double val)
	{
		throw_L_ArgException_if(ptr==NULL, "L_DataType::setValue");
		switch(tipo)
		{
		case L_ty_undefined: *(double *)ptr=(double)val; break;
		case L_ty_bool: *(bool *)ptr=(val>0.1 || val<=-0.1)?1:0; break;
		case L_ty_char: *(char *)ptr=(char)val; break;
		case L_ty_uchar: *(unsigned char *)ptr=(unsigned char)val; break;
		case L_ty_schar: *(signed char *)ptr=(signed char)val; break;
		case L_ty_short: *(short *)ptr=(short)val; break;
		case L_ty_int: *(int *)ptr=(int)val; break;
		case L_ty_long: *(long *)ptr=(long)val; break;
		case L_ty_ushort: *(unsigned short *)ptr=(unsigned short)val; break;
		case L_ty_uint: *(unsigned int *)ptr=(unsigned int)val; break;
		case L_ty_ulong: *(unsigned long *)ptr=(unsigned long)val; break;
		case L_ty_double: *(double *)ptr=(double)val; break;
		case L_ty_float:  *(float *)ptr=(float)val; break;
		default: throw_L_ArgException_if(true, "L_DataType::setValue");
		}
	}

	inline char *allToText(char *buf)
	{
		throw_L_ArgException_if(ptr==NULL, "L_DataType::allToText");
		switch(tipo)
		{
		case L_ty_undefined: sprintf(buf,"(indef)%g",*(double*)ptr);break;
		case L_ty_bool: if (*(bool*)ptr==true) strcpy(buf,"(bool)true"); else strcpy(buf,"(bool)false"); break;
		case L_ty_char: sprintf(buf,"(char)%d",*(char*)ptr); break;
		case L_ty_uchar: sprintf(buf,"(unsigned char)%d",*(unsigned char*)ptr); break;
		case L_ty_schar: sprintf(buf,"(signed char)%d",*(signed char*)ptr); break;
		case L_ty_short: sprintf(buf,"(short)%d",*(short*)ptr); break;
		case L_ty_int: sprintf(buf,"(int)%d",*(int*)ptr); break;
		case L_ty_long: sprintf(buf,"(long)%ld",*(long*)ptr); break;
		case L_ty_ushort: sprintf(buf,"(unsigned short)%d",*(unsigned short*)ptr); break;
		case L_ty_uint: sprintf(buf,"(unsigned int)%d",*(unsigned int*)ptr); break;
		case L_ty_ulong: sprintf(buf,"(unsigned long)%ld",*(unsigned long*)ptr); break;
		case L_ty_double: sprintf(buf,"(double)%g",*(double*)ptr); break;
		case L_ty_float:  sprintf(buf,"(float)%g",*(float*)ptr); break;
		default: throw_L_ArgException_if(true, "L_DataType::allToText");
		}
		return buf;
	}

	bool readValuesFromLine(const char *buf);

	template <class T> inline void set(T *var) {tipo=typeOfData(var); ptr=var;} // Para tipo basico (no objetos)

	void clean() {ptr=&temporal; tipo=L_ty_undefined;}
};

// parametro etiquetado
class L_ParamLabel
{
private:
	L_RefData data;
	char *name;
	char *class_name;

public:
	L_ParamLabel() {name=NULL; class_name=NULL;}
	L_ParamLabel(const L_ParamLabel &other) {set(other.name,other.data); setClass(other.class_name);}
	L_ParamLabel &operator=(const L_ParamLabel &other) {set(other.name,other.data); setClass(other.class_name); return *this;}
	void setName(const char *name);
	void setClass(const char *class_name);
	void setFromString(const char *value);
	template <class T> inline void set(const char *name, T *ptr) {setName(name); data.set(ptr);}
	inline void set(const char *name, L_RefData otherDato) {setName(name); data=otherDato;}
	void swap(L_ParamLabel &other) {L_RefData dato_t = other.data; char *nombre_t = other.name; char *clase_t = other.class_name; other.data=data; other.name=name; other.class_name=class_name; data = dato_t; name = nombre_t; class_name = clase_t;}
	~L_ParamLabel() {if (name!=NULL) delete[] name; if (class_name!=NULL) delete[] class_name;}
	static int cmpInv(const void *a, const void *b);
	static int cmpNorm(const void *a, const void *b);
	friend class L_ParamBlock;
	friend class L_ParamLabelList;
};

typedef L_Node<L_ParamLabel> L_ParamLabelNode;

class L_ParamLabelList:public L_List<L_ParamLabel>
{
public:
	bool readFile(FILE *fp, bool cerrarArchivo=true);
	bool saveFile(FILE *fp, bool cerrarArchivo=true);
};

class L_ParamManagerLocal;

// Conjunto de parametros etiquetados
class L_ParamBlock
{
private:
	L_ParamLabel *pa;
	int nPaTot; // Numero de elementos de pa que tienen un value fijado
	int nPaMem; // Numero de elementos de pa para los que hay memoria pedida
	bool ordered;
	char *class_name;

	L_ParamBlock() {} // No debe existir
	void addFrom(const char *name, void *ptr); // No debe existir
	void copyOn_naive(L_ParamBlock &other) {other.pa=pa; other.nPaTot=nPaTot; other.nPaMem=nPaMem; other.ordered=ordered; other.class_name=class_name;}
	void undefine_naive() {pa=NULL; class_name=NULL;}

	inline void sort() {qsort(pa, nPaTot, sizeof(L_ParamLabel), &L_ParamLabel::cmpNorm); ordered=true;}
	L_RefData *findVariable(const char *name);
	void revMemo();
	bool belongTo(const L_ParamLabel &param) {return strcmp(class_name, param.class_name)==0;}
	bool updateValue(const L_ParamLabel &other);

	L_ParamBlock(const L_ParamBlock &other)
	{
		int i;
		class_name=new char[strlen(other.class_name)+1];
		strcpy(class_name,other.class_name);
		nPaTot = other.nPaTot;
		nPaMem = other.nPaMem;
		pa = new L_ParamLabel[nPaMem];
		for (i=0; i<nPaTot; i++)
			pa[i] = other.pa[i];
		ordered=other.ordered;
		revMemo();
	}
	L_ParamBlock &operator=(const L_ParamBlock &other)
	{
		int i;
		if (class_name != NULL)
			delete[] class_name;
		class_name=new char[strlen(other.class_name)+1];
		strcpy(class_name,other.class_name);
		nPaTot = other.nPaTot;
		nPaMem = other.nPaMem;
		if (pa != NULL)
			delete[] pa;
		pa = new L_ParamLabel[nPaMem];
		for (i=0; i<nPaTot; i++)
			pa[i] = other.pa[i];
		ordered=other.ordered;
		revMemo();
		return *this;
	}

public:
	L_ParamBlock(const char *className)
	{
        class_name=new char[strlen(className)+1];
		strcpy(class_name,className);
		pa=NULL;
		nPaTot=0;
		nPaMem=0;
		ordered=true;
		revMemo();
	}
	void cambiaNombreClase(const char *className)
	{
		int i;
		if (class_name!=NULL)
			delete[] class_name;
		class_name=new char[strlen(className)+1];
		strcpy(class_name,className);
		for (i=0; i<nPaTot; i++)
			pa[i].setClass(class_name);
	}

	void swap(L_ParamBlock &other)
	{
		L_ParamLabel *pa_t = other.pa;
		int nPaTot_t = other.nPaTot;
		int nPaMem_t = other.nPaMem;
		bool ordenado_t = other.ordered;
		char *nombreClase_t = other.class_name;

		other.pa = pa;
		other.nPaTot = nPaTot;
		other.nPaMem = nPaMem;
		other.ordered = ordered;
		other.class_name = class_name;

		pa = pa_t;
		nPaTot = nPaTot_t;
		nPaMem = nPaMem_t;
		ordered = ordenado_t;
		class_name = nombreClase_t;
	}

	template <class T> inline void addFrom(const char *name, T *ptr) {revMemo(); pa[nPaTot].set(name, ptr); pa[nPaTot].setClass(class_name); ordered=false; nPaTot++;}
	void updateAssociatedValues(L_ParamLabelList &lista, bool borrarValoresUsados);
	void buildListFromAssociatedValues(L_ParamLabelList &lista);

	~L_ParamBlock() {if (pa!=NULL) delete[] pa; if (class_name!=NULL) delete[] class_name;}
};


// Local parameter manager with children
class L_ParamManagerLocal
{
private:
	L_ParamBlock params;
	L_ParamManagerLocal **children;
	int nHijosTot;
	int nChildrenMem;
public:
	bool modifiedFlowOfInfo;

private:
	L_ParamManagerLocal(const L_ParamManagerLocal &other); // Es complicado copiar este object pq contiene punteros a variables externas
	L_ParamManagerLocal &operator =(const L_ParamManagerLocal &other); // Es complicado copiar este object pq contiene punteros a variables externas

public:
	L_ParamManagerLocal(const char *class_name, int nChildrenMem):params(class_name)
	{
		nHijosTot=0;
		modifiedFlowOfInfo=true;
		this->nChildrenMem=nChildrenMem;
		if (nChildrenMem>0)
			children=new L_ParamManagerLocal*[nChildrenMem];
		else
			children=NULL;
	}
	void cambiaNombreClase(const char *class_name)
	{
		params.cambiaNombreClase(class_name);
	}
	void swap(L_ParamManagerLocal &other);
	bool addChildren(L_ParamManagerLocal *nuevoHijo);
	void olvidaHijos();
	void updateValues_me_and_children(L_ParamLabelList &lista);
	void buildListOfValues_me_and_children(L_ParamLabelList &lista);
	void setValue_modifiedFlowOfInfo(bool alterar) {modifiedFlowOfInfo=alterar; for (int i=0; i<nHijosTot; i++) children[i]->setValue_modifiedFlowOfInfo(alterar);}
	template <class T> inline void addFrom(const char *name, T *ptr) {params.addFrom(name, ptr);}

	~L_ParamManagerLocal() {if (children!=NULL) delete[] children;}
};
#endif


// Abstract optimizer
class L_Optimizer
{
public:
	L_Optimizer *children;
	std::vector<double> vectorLimInf, vectorLimSup; // Son punteros a memo externa
	int lenVector;
	L_Optimizer() {children=NULL; lenVector=0;}
	virtual bool readyToWork(bool print) = 0;
	inline void set_lenVector_to_children() {if (children!=NULL) {children->lenVector=lenVector; children->set_lenVector_to_children();}}
	inline void set_vectorLim_to_children() {if (children!=NULL) {children->vectorLimInf=vectorLimInf; children->vectorLimSup=vectorLimSup; children->set_vectorLim_to_children();}}
	inline bool listoParaTrabajarBase() {if (lenVector==0) return false; if (children!=NULL) {if (children->lenVector!=lenVector) return false;} return true;}
	virtual double minimize(const void *object, std::vector<double> &vector, double fnToMinimize(const void *object, double *vector)) = 0;
	~L_Optimizer() { }
};

class L_LevenbergMarquardt: public L_Optimizer
{
public:
	long nIterationsMax;
	double errorToStop;
	double factorGrowth; // Growth/reduction factor for lambda
	double epsilon;
	double lambda; // LM dumping parameter
	bool useHessian;
	L_LevenbergMarquardt() {nIterationsMax=16777215L; errorToStop=-1e30; epsilon=1e-4; factorGrowth=10;lambda = 0.1; useHessian=true;}
	bool readyToWork(bool print);

	double minimize(const void *object, std::vector<double> &vector, double fnToMinimize(const void *object, double *vector));

	double minimize_vect(const void *object, L_Matrix &x, void (*fnToMinimize)(const void *object, const L_Matrix &x, L_Matrix &errVect), void (*jacob)(const void *object, const L_Matrix &x, L_Matrix &Jdedx)=NULL, void (*preIter)(const void *object, L_Matrix &x) = NULL);
};

class L_OptimizationPath
{
public:
	L_Array<L_Matrix> puntos;
	std::vector<double> errores;
};

class L_Random_plus_LM_optimizer: public L_Optimizer
{
public:
	L_LevenbergMarquardt LM;
	L_Array<L_OptimizationPath> paths;
	long nIterationsMax;
	bool usarLM;

	L_Random_plus_LM_optimizer() {nIterationsMax = 16777215L; usarLM = true;}
	bool readyToWork(bool print);
	double minimize(const void *object, std::vector<double> &vector, double fnToMinimize(const void *object, double *vector));
	double minimize_vect(const void *object, L_Matrix &x, void (*fnToMinimize)(const void *object, const L_Matrix &x, L_Matrix &errVect), void (*jacob)(const void *object, const L_Matrix &x, L_Matrix &Jdedx) = NULL, void (*preIter)(const void *object, L_Matrix &x) = NULL);
};

#define L_LinearQuantizer_nCluster(cuant,in) (L_FLOOR1000(cuant.m*(in)+cuant.n))
#define L_LinearQuantizer_nClusterLeft(cuant,in) (L_FLOOR1000(cuant.m*(in)+cuant.n-0.5))
#define L_LinearQuantizer_nClusterRight(cuant,in) (L_FLOOR1000(cuant.m*(in)+cuant.n+0.5))
#define L_LinearQuantizer_rCenCluster(cuant,nclu) (L_FLOOR1000(((nclu)-cuant.n+0.5)/cuant.m))

class L_LinearQuantizer
{
public:

	double m;
	double n;

	inline L_LinearQuantizer()	{m=0;n=0;}
	// Sirve si se sabe que un rango [rmin1,*[ corresponde a un cluster nclu1 y un rango [rmin2,*[ corresponde a un cluster nclu2
	inline void setFromTransitions(double rmin1, int nclu1, double rmin2, int nclu2)
	{
		m=(nclu2-nclu1)/(rmin2-rmin1);
		n=nclu1-m*rmin1;
	}
	// Sirve si se sabe que un rango [rmin,rmax[ corresponde a un cluster nclu
	inline void setFromCellRange(double rmin, double rmax, int nclu)
	{
		m=1/(rmax-rmin);
		n=nclu-m*rmin;
	}
	// Number of cluster associated to in
	inline int nCluster(double in)	{ return L_FLOOR(m*in+n);}
	// 	Number of cluster associated to in (floor)
	inline int nClusterLeft(double in) { return L_FLOOR(m*in+n-0.5); }
	// 	Number of cluster associated to in (ceil)
	inline int nClusterRight(double in) { return L_FLOOR(m*in+n+0.5); }

	inline int nCluster_Left_or_Right(double in, int selecc) { return L_FLOOR(m*in+n-0.5+selecc); }

	inline double rCenCluster(int nclu)	{ return (nclu-n+0.5)/m; }

	inline double rMinCluster(int nclu)	{ return (nclu-n)/m; }

	inline double rMaxCluster(int nclu) { return (nclu-n+1)/m; }
};

class L_Statistics1D
{
private:
	long n;
	double sx;
	double sxx;
public:
	L_Statistics1D() {clear();}
	void push(double x) {n++; sx+=x; sxx+=x*x;}
	void pop(double x) {n--; sx-=x; sxx-=x*x;}
	double getMean() {return sx/n;}
	double getVariance() {return sxx/n - (sx/n)*(sx/n);}
	int size() {return n;}
	void clear() {n=0; sx=0; sxx=0;}
};

#define L_Statistics2D_push(est,x1,x2) ((est).sx+=(x1), (est).sy+=(x2), (est).sxx+=(x1)*(x1), (est).sxy+=(x1)*(x2), (est).syy+=(x2)*(x2), (est).n++)
#define L_Statistics2D_clear(est) ((est).n = 0, (est).sx=0, (est).sy=0, (est).sxx=0, (est).sxy=0, (est).syy=0)

class L_Statistics2D
{
public:
	double sx;
	double sy;
	double sxx;
	double sxy;
	double syy;
	int n;
	L_Statistics2D() {sx=0; sy=0; sxx=0; sxy=0; syy=0; n=0;}
	void clear() {n=0; sx=0; sy=0; sxx=0; sxy=0; syy=0;}
	size_t size() const {return n;}
	void invertAxis() {double t; t=sx; sx=sy; sy=t;  t=sxx; sxx=syy; syy=t;}
	void mirrorFrom(const L_Statistics2D &other) {n=other.n; sx=other.sy; sy=other.sx; sxx = other.syy; sxy=other.sxy; syy=other.sxx;}
	void push_pair(double x1, double x2) {sx+=x1;sy+=x2;sxx+=x1*x1;sxy+=x1*x2;syy+=x2*x2; n++;}
	void pop_pair(double x1, double x2) {sx-=x1;sy-=x2;sxx-=x1*x1;sxy-=x1*x2;syy-=x2*x2; n--;}
	double getMeanX1() {return sx/n;}
	double getMeanX2() {return sy/n;}
	double getVarianceX1X1() {return sxx/n - (sx/n)*(sx/n);}
	double getVarianceX1X2() {return sxy/n - (sx/n)*(sy/n);}
	double getVarianceX2X2() {return syy/n - (sy/n)*(sy/n);}
	double getRadioFromVariance() {return sqrt(sxx/n - (sx/n)*(sx/n)) + sqrt(syy/n - (sy/n)*(sy/n)); } // r = sqrt(Vxx+Vyy);
	double getRadioFromVariance_X() {return 2*sqrt(sxx/n - (sx/n)*(sx/n)); } // r = sqrt(Vxx+Vyy);
	double getRadioFromVariance_Y() {return 2*sqrt(syy/n - (sy/n)*(sy/n)); } // r = sqrt(Vxx+Vyy);
	double getSideOfSquareFromVariance_X() {return 2*sqrt(3.0)*sqrt(sxx/n - (sx/n)*(sx/n));};
	double getSideOfSquareFromVariance_Y() {return 2*sqrt(3.0)*sqrt(syy/n - (sy/n)*(sy/n));};
	double getCorrelationCoefficient() {return (sxy/n-sx/n*sy/n) / sqrt((sxx/n-sx/n*sx/n)*(syy/n-sy/n*sy/n)) ;}
	void getLine(double &mRecta, double &nRecta, double vex=0) {double v11 = getVarianceX1X1(), v12 = getVarianceX1X2(); mRecta = (v12 / (v11 - vex)); nRecta = (sy/n) - mRecta*(sx/n);}
	double getLineOrientation(double vex=0) {double m, n; if (sxx > sxy) {getLine(m, n, vex); return atan2(m,1);} else {L_Statistics2D other; other.mirrorFrom(*this); other.getLine(m, n); return atan2(1,m);}}
	// Straight line considering errors in both axis: x=xideal+ex; y=yideal+ey; delta = var(ey)/var(ex); => y = m + n*xideal + ey; 
	void getLineDeming(double &mRecta, double &nRecta, double delta=1.0) {double vxx = getVarianceX1X1(), vxy = getVarianceX1X2(), vyy = getVarianceX2X2(); mRecta = (vyy-delta*vxx+sqrt((vyy-delta*vxx)*(vyy-delta*vxx)+4*delta*vxy*vxy)) / (2*vxy); nRecta = (sy/n) - mRecta*(sx/n);}
	double getLineOrientationDeming(double delta=1.0) {double m, n; if (sxx > sxy) {getLineDeming(m, n, delta); return atan2(m,1);} else {L_Statistics2D other; other.mirrorFrom(*this); other.getLineDeming(m, n, 1/delta); return atan2(1,m);}}
	void correctPointDeming(double mRecta, double nRecta, double delta, double x, double y, double &xc, double &yc) {xc = x + mRecta/(mRecta*mRecta+delta)*(y-nRecta-mRecta*x); yc = mRecta+xc+nRecta;}
	double correctX1Deming(double x1, double x2, double mRecta, double nRecta, double delta=1.0) {return x1+mRecta/(mRecta*mRecta+delta)*(x2-nRecta-mRecta*x1);}
};

class L_FileName
{
private:
	L_String _dirNameExt; // It is used only for returning a (const char *)
public:
	L_String dir;
	L_String name;
	L_String ext;

	L_FileName() {}
	L_FileName(const char *dir, const char *name, const char *ext) {this->dir=dir;this->name=name;this->ext=ext;}
	L_FileName(const char *dirnameext) {asign(dirnameext);}
	L_FileName(const L_FileName& other) {dir=other.dir;name=other.name;ext=other.ext;}
	L_FileName& operator=(const L_FileName& other) {dir=other.dir;name=other.name;ext=other.ext;return *this;}
	L_FileName& asign(const char *dirnameext)
	{
		char dir_[1024];
		char nom_[1024];
		char ext_[1024];
		L_div_name_file2(dirnameext, dir_, nom_, ext_);
		dir=dir_;
		name=nom_;

		ext=ext_;
		return *this;
	}

	const char *dirnameext()
	{
		char dirsep[2];
		dirsep[0]=FOLDER_SEPARATION_CHAR;
		dirsep[1]=0;

		_dirNameExt=name;
		if (*(const char *)dir.begin()!=0 && dir[0] != 0)
			_dirNameExt=dir+dirsep+_dirNameExt;
		if (*(const char *)ext.begin()!=0 && ext[0] != 0)
			_dirNameExt=_dirNameExt+"."+ext;
		return _dirNameExt.begin();
	}
	const char * c_str() { return dirnameext();}

	bool hasExtension(const char *str)
	{
		int i;
		char c1, c2;
		for (i=0; str[i]!=0 && ext[i]!=0 && i<(int)ext.size(); i++)
		{
			c1 = str[i];
			c1 += (c1>='A' && c1<='Z') * ('a'-'A');
			c2 = ext[i];
			c2 += (c2>='A' && c2<='Z') * ('a'-'A');
			if (c1 != c2)
				return false;
		}
		return true;
	}

	~L_FileName()
	{
	}
};

enum L_ShapeType
{
	L_Shape_undefined,
	L_Shape_line,
	L_Shape_circle,
	L_Shape_ellipse,
	L_Shape_dashed_ellipse
};

class L_Shape
{
public:
	int x1;
	int y1;
	int n1; // nº de imagen 1 (generalmente imagen de referencia)
	int x2;
	int y2;
	int n2; // nº de imagen 2 (generalmente imagen de test)
	double rx, ry, ang; // Parametros para elipse centrada en (x1,y1)
	L_uchar R;
	L_uchar G;
	L_uchar B;
	L_ShapeType tipo;

	L_Shape() {}
	L_Shape(int x1, int y1, int n1, int x2, int y2, int n2, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B, L_ShapeType tipo)
	{
		this->x1=x1; this->y1=y1; this->n1=n1;
		this->x2=x2; this->y2=y2; this->n2=n2;
		this->rx=rx; this->ry=ry; this->ang=ang;
		this->R=R; this->G=G; this->B=B;
		this->tipo = tipo;
	}
};

//typedef L_Node<L_Shape> L_LineaNodo;

class L_ShapeArray
{
public:
	std::vector<L_Shape> v;
	typedef int size_type;

	static void genRandomColor(L_uchar &R, L_uchar &G, L_uchar &B);
	void drawLine(int x1, int y1, int x2, int y2, L_uchar R, L_uchar G, L_uchar B)    {_drawShape(0, 0, x1, y1, x2, y2, R, G, B);}
	void drawCircle(int xc, int yc, double r, L_uchar R, L_uchar G, L_uchar B)        {_drawCircle(0,0,xc,yc,r,R,G,B);}
	void drawEllipse(int xc, int yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B)         {_drawEllipse(0,0,xc,yc,rx,ry,ang,R,G,B);}
	void drawDashedEllipse(int xc, int yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B)  {_drawDashedEllipse(0,0,xc,yc,rx,ry,ang,R,G,B);}
	void drawRectangle(int x1, int y1, int x2, int y2, L_uchar R, L_uchar G, L_uchar B)           {_drawRectangle(0, 0, x1, y1, x2, y2, R, G, B);}
	void drawRotatedSquare(int xc, int yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B)    {_drawRotatedSquare(0, 0, xc, yc, r, ang, R, G, B);}
	void drawArrow(int xc, int yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B)       {_dibFlecha(0, 0, xc, yc, r, ang, R, G, B);}
	void drawNumber(int xc, int yc, int num, int radio, L_uchar R, L_uchar G, L_uchar B);
	void changeSizeFactor(double sizeFactor);
	void drawTracks(L_Matrix &m, L_Array<bool> *v=NULL);
	void drawLineByParameters(double a, double b, double c, int lx, int ly, L_uchar R=255, L_uchar G=255, L_uchar B=255); // Dibuja recta a*x+b*y+c=0

	static void intersectionsWithScreen(double m, double n, int lx, int ly, int x1, int x2) {double w, v; int xa, xb; x1 = 0; x2 = lx-1; if (m!=0.0 && m!=-0.0) {w=1/m; v=-n/m; x1 = (int)(v); x2 = (int)(w*(lx-1) + v); if (xa > xb) {int t; t=xa; xa=xb; xb=t;} if (xa > x1) x1 = xa; if (xb < x2) x2 = xb;} if (x1 < 0) x1 = 0; if (x1 >= lx) x1 = lx-1; if (x2 < 0) x2 = 0; if (x2 >= lx) x2 = lx-1; if (m>0) {} if (m<0) {}}

	size_type size() const {return (size_type)v.size();}
	void resize(size_type n) {v.resize(n);}
	void swap(L_ShapeArray &other) {v.swap(other.v);}
	L_Shape &operator [] (size_type i) {return v[i];}
	const L_Shape &operator [] (size_type i) const {return v[i];}
	void push_back(const L_Shape &sh) {v.push_back(sh);}

	void setColorForAllShapes(L_uchar R, L_uchar G, L_uchar B)
	{
		for (size_type i=0; i<(size_type)v.size(); i++)
		{
			v[i].R=R;
			v[i].G=G;
			v[i].B=B;
		}
	}

	void drawLine(size_type x1, size_type y1, size_type x2, size_type y2)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		drawLine(x1, y1, x2, y2, R, G, B);
	}
	void drawCircle(size_type xc, size_type yc, double r)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		drawCircle(xc, yc, r, R, G, B);
	}
	void drawRectangle(size_type x1, size_type y1, size_type x2, size_type y2)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		drawRectangle(x1, y1, x2, y2, R, G, B);
	}
	void drawRotatedSquare(size_type xc, size_type yc, double r, double ang)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		drawRotatedSquare(xc, yc, r, ang, R, G, B);
	}
	void drawArrow(size_type xc, size_type yc, double r, double ang)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		drawArrow(xc, yc, r, ang, R, G, B);
	}

	void _drawShape(size_type im1, size_type im2, size_type x1, size_type y1, size_type x2, size_type y2, L_uchar R, L_uchar G, L_uchar B);
	void _drawShape_push_front(size_type im1, size_type im2, size_type x1, size_type y1, size_type x2, size_type y2, L_uchar R, L_uchar G, L_uchar B);

	void write(FILE *fp);
	void read(FILE *fp);

	void _drawCircle(size_type im1, size_type im2, size_type xc, size_type yc, double r, L_uchar R, L_uchar G, L_uchar B);
	void _drawEllipse(size_type im1, size_type im2, size_type xc, size_type yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B);
	void _drawDashedEllipse(size_type im1, size_type im2, size_type xc, size_type yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B);
	void _drawRectangle(size_type im1, size_type im2, size_type x1, size_type y1, size_type x2, size_type y2, L_uchar R, L_uchar G, L_uchar B);
	void _drawRotatedSquare(size_type im1, size_type im2, size_type xc, size_type yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B);
	void _dibFlecha(size_type im1, size_type im2, size_type xc, size_type yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B);

	void _drawShape(size_type im1, size_type im2, size_type x1, size_type y1, size_type x2, size_type y2)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		_drawShape(im1, im2, x1, y1, x2, y2, R, G, B);
	}
	void _drawCircle(size_type im1, size_type im2, size_type xc, size_type yc, double r)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		_drawCircle(im1, im2, xc, yc, r, R, G, B);
	}
	void _drawRectangle(size_type im1, size_type im2, size_type x1, size_type y1, size_type x2, size_type y2)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		_drawRectangle(im1, im2, x1, y1, x2, y2, R, G, B);
	}
	void _drawRotatedSquare(size_type im1, size_type im2, size_type xc, size_type yc, double r, double ang)
	{
		L_uchar R, G, B;
		genRandomColor(R, G, B);
		_drawRotatedSquare(im1, im2, xc, yc, r, ang, R, G, B);
	}

	void setImageNumber_ref(size_type numRef)
	{
		for (size_type i=0; i<(size_type)v.size(); i++)
			v[i].n1 = numRef;
	}

	void setImageNumber_cur(size_type numCur)
	{
		for (size_type i=0; i<(size_type)v.size(); i++)
			v[i].n2 = numCur;
	}

	static void testShapes();

public:
	~L_ShapeArray() {}
};

// Nested representation for mathematical expressions
class L_MatParsParenNode
{
public:
	int i1;
	int i2;
	L_Array<L_MatParsParenNode *> children;
	L_MatParsParenNode *father;
	L_MatParsParenNode() {father=NULL;}
	void destroyRec() {int i; for (i=0; i<(int)children.size(); i++) if (children[i]!=NULL) {children[i]->destroyRec(); delete children[i]; children[i]=NULL;}}
};

enum L_MatParsNodeType {LMP_undef, LMP_variable, LMP_constant, LMP_operator, LMP_function, LMP_preoperator, LMP_postoperator};

class L_MatParser;

// Nested representation for mathematical expressions
class L_MatParsNode
{
public:
	bool useful; // flag for useful content
	L_String str;
	int level;
	int preced;
	int pre_i1, pre_i2, pos_i1, pos_i2;
	L_MatParsNodeType mytype;
	L_MatParsNode *preArg; // Left children
	L_MatParsNode *posArg; // Right children
	L_MatParsNode() {useful=false; preArg=NULL; posArg=NULL; mytype=LMP_undef;}
	L_MatParsNode *clonaRec() {L_MatParsNode *ret=new L_MatParsNode(); *ret=*this; if (preArg!=NULL) ret->preArg=preArg->clonaRec(); if (posArg!=NULL) ret->posArg=posArg->clonaRec(); return ret;}
	void destroyRec() {if (preArg!=NULL) {preArg->destroyRec(); delete preArg; preArg=NULL;} if (posArg!=NULL) {posArg->destroyRec(); delete posArg; posArg=NULL;}}
	const char *tipoStr();
	const static char * strUndef;
	const static char * strVariable;
	const static char * strConstant;
	const static char * strOperator;
	const static char * strFunction;
	const static char * strPreoperator;
	const static char * strPostoperator;
	L_ComplexDouble evaluateComplexRec(const L_MatParser &mat);
	double evaluateRealRec(const L_MatParser &mat);
	long evaluateIntegerRec(const L_MatParser &mat);
};

class L_MatParsPlaceholder
{
public:
	L_MatParsPlaceholder() {}
	L_MatParsPlaceholder(const char *str, double re, double im=0) {varNom=str; varRe=re; varIm=im;}
	L_String varNom;
	double varRe;
	double varIm;
	L_MatParsPlaceholder &operator =(const L_MatParsPlaceholder &other) {varRe=other.varRe; varIm=other.varIm; varNom=other.varNom; return *this;}
};

// Math expression parser and evaluator
class L_MatParser
{
public:
	L_MatParsNode *root; // Math expression
	std::vector<L_MatParsPlaceholder> *reemplazos; // Variables in expression

	L_MatParser() {root=NULL; reemplazos = NULL;}
	~L_MatParser() {if (root!=NULL) {root->destroyRec(); delete root; root=NULL;}}

	// Analyze expression
	void buildFrom(const char *expresion, bool debug=false);

	// Evaluate expression
	L_ComplexDouble evaluateComplex();
	long evaluateInteger(); // Redondea todos los numeros involucrados a enteros
	double evaluateReal();

	void destroy() {if (root!=NULL) {root->destroyRec(); delete root; root=NULL;}}
	L_MatParser &operator=(const L_MatParser &other)
	{
		destroy();
		memcpy(this, &other, sizeof(*this));
		if (other.root!=NULL)
			root=other.root->clonaRec();
		return *this;
	} // No copia reemplazos_

	// Para debug
	void showStructureDebug(FILE *fp=stdout);
	void showStructure(FILE *fp=stdout);

	// Search for variable value from string
	int searchIndexOfPlaceholder(const char *str) const;
	bool searchVariableValue(const char *str, double &num) const {int i; i=searchIndexOfPlaceholder(str); if (i==-1) return false; num=(*reemplazos)[i].varRe; return true;}
	bool searchVariableValue(const char *str, L_ComplexDouble &num) const {int i; i=searchIndexOfPlaceholder(str); if (i==-1) return false; num.re=(*reemplazos)[i].varRe; num.im=(*reemplazos)[i].varIm; return true;}
	long searchVariableValue(const char *str, long &num) const {int i; i=searchIndexOfPlaceholder(str); if (i==-1) return false; num=(long)(*reemplazos)[i].varRe; return true;}

private:
	bool buildNodeRecursive(L_MatParsNode *raizPars, L_String &ultExpr);
	void genCanonicalExpr(const char *expresion, L_String &ultExpr);
	void analyzeParenthesisRecursive(L_MatParsNode *raizPars, L_Array<L_MatParsParenNode *> &mapaParen, const char *expr, int level);
	void showStructureDebugRecursive(FILE *fp, L_MatParsNode *nodo, int level);
	void showStructureRecursive(FILE *fp, L_MatParsNode *nodo, int level);
};

// For generating strings depending on a number, ex. "nombre00001.ext"
class L_StringWithNumber : public L_String
{
public:
	L_String prefix;
	int numDigits;
	L_String suffix;
	int numMax;

	L_StringWithNumber() {numDigits = -1; numMax = -1;}
	L_StringWithNumber(const L_StringWithNumber &other) {L_String::operator =(other); prefix = other.prefix; numDigits = other.numDigits; suffix = other.suffix; numMax = other.numMax;}
	L_StringWithNumber(const char *prefix, int numDigits, const char *suffix, int numMax) {this->prefix = prefix; this->numDigits = numDigits; this->suffix = suffix; this->numMax = numMax;}

	const L_StringWithNumber& operator=(const char *s) {L_String::operator =(s); numDigits = -1; numMax = -1; return *this;}
	const L_StringWithNumber& operator=(const L_String &other) {L_String::operator =(other); numDigits = -1; numMax = -1; return *this;}
	const L_StringWithNumber& operator=(const L_StringWithNumber &other) {L_String::operator =(other); prefix = other.prefix; numDigits = other.numDigits; suffix = other.suffix; numMax = other.numMax; return *this;}

	void define(const char *prefix, int numDigits, const char *suffix, int numMax) {this->prefix = prefix; this->numDigits = numDigits; this->suffix = suffix; this->numMax = numMax;}
	void define(const L_String &prefix, int numDigits, const char *suffix, int numMax) {define(prefix.c_str(), numDigits, suffix, numMax);}

	const char *generateString(int i);
	void write(FILE *fp);
	bool read(FILE *fp);

	static void test();
};


class L_BMPheaderInfo
{
public:
	//header (14B)
	L_uint8 _bmp_B; // 0
	L_uint8 _bmp_M; // 1
	L_uint32 _sizeOfFile; //2-5
	L_uint16 _reserv1; //6-7
	L_uint16 _reserv2; //8-9
	L_uint32 _offsetData;  //10-13 //Ah-Dh
	// info (40B)
	L_uint32 headerSize; //14-17 // Eh-11h
	L_int32 width; //18-21 //12h-15h
	L_int32 height;  //22-25 //16h-19h
	L_uint16 planes; //26-27 //1Ah-1Bh
	L_uint16 bits;   //28-29 //1Ch-1Dh
	L_uint32 compr;  //30-33 //1Eh-21h
	L_uint32 sizeOfImage; //34-37 //22h-25h
	L_int32 pixByEachMeterX; //38-41 //26h-29h
	L_int32 pixByEachMeterY; //42-45 //2Ah-2Dh
	L_uint32 nColors; //46-49 //30h-33h
	L_uint32 nColorsImp; //50-53 //34h-37h
	// extra
	long whiteSpace;
	L_BMPheaderInfo();
	bool setSize(L_uint32 width, L_uint32 Alto, L_uint16 Bits=24);
	void write(FILE *fp);
	bool read(FILE *fp);
	~L_BMPheaderInfo() {}
};


template <class T>  // x must have operators +, -
class L_Spline1DGen
{
public:
	std::vector<double> x;  // Values for independant variable
	std::vector<T> y;  // Values from the function
	std::vector<T> y2; // Internal values (deriatives)
	int n;

	L_Spline1DGen() : n(0) { }
	inline void addFrom(double xVal, T yVal) {n++; x.resize(n);y.resize(n); x[n-1] = xVal; y[n-1]=yVal;}
	inline void push_back(double xVal, T yVal) {n++; x.resize(n);y.resize(n); x[n-1] = xVal; y[n-1]=yVal;}

	void compute()
	{
		throw_L_ArgException_if (n<=0, "L_Spline1D::compute() : datos indefinidos");
		std::vector<double> aim(n), ai(n), aiM(n), bim(n), bi(n), biM(n);
		std::vector<T> c(n);
		T zero;
		int i;
		double fac;
		zero = y[0] - y[0];
		y2.resize(n);
		// a[i][i-1] = aim[i] ; a[i][i] = ai[i] ; a[i][i+1] = aiM[i]
		// b[i][i-1] = bim[i] ; b[i][i] = bi[i] ; b[i][i+1] = biM[i]
		for (i=1; i<n-1; i++)
		{
			aim[i] = (x[i] - x[i-1])/6;
			ai[i]  = (x[i+1]-x[i-1])/3;
			aiM[i] = (x[i+1]-x[ i ])/6;
		}
		ai[0] = 1; // Para spline suave
		aiM[0] = 0;
		aim[n-1] = 0;
		ai[n-1] = 1;
		for (i=1; i<n-1; i++)
			c[i] = (y[i+1]-y[i]) /(x[i+1]-x[i]) - (y[i]-y[i-1]) /(x[i]-x[i-1]);
		c[0] = zero;
		c[n-1] = zero;
		for (i=0; i<n; i++)
		{
			bim[i] = 0;  // La identity al principio
			bi[i] = 1;
			biM[i] = 0;
		}
		// [...                                                   | ...]
		// [...  aim[i-1]  ai[i-1]   aiM[i-1]                     | ...]
		// [...            aim[i]     ai[i]     aiM[i]            | ...]
		// [...                      aim[i-1]   ai[i]    aiM[i+1] | ...]
		// [...                                                   | ...]
		for (i=0; i<n-2; i++) // eliminar hacia abajo [A|B]
		{
			// Dejar la diagonal en 1
			// aim[i] = 0
			fac = 1.0/ai[i];
			ai[i] *= fac;
			aiM[i] *= fac;
			bim[i] *= fac;
			bi[i] *= fac;
			biM[i] *= fac;
			// Restar la fila de arriba a la de abajo
			fac = aim[i+1]/ai[i];
			aim[i+1] -= fac*ai[i];
			ai[i+1] -= fac*aiM[i];
			bim[i+1] -= fac*bi[i+1];
			bi[i+1] -= fac*biM[i+1];
		}
		fac = 1.0/ai[n-2];
		ai[n-2] *= fac;
		bim[n-2] *= fac;
		bi[n-2] *= fac;

		// [...                                                   | ...]
		// [...  aim[i-1]  ai[i-1]   aiM[i-1]                     | ...]
		// [...            aim[i]     ai[i]     aiM[i]            | ...]
		// [...                      aim[i-1]   ai[i]    aiM[i+1] | ...]
		// [...                                                   | ...]
		for (i=n-2; i>=1; i--)
		{
			// Los aim = 0
			fac = aiM[i-1]/ai[i];
			aiM[i-1] -= fac*ai[i];
			ai[i-1] -= fac*aim[i];
			biM[i-1] -= fac*bi[i];
			bi[i-1] -= fac*bim[i];
		}
		fac = 1.0/ai[0];
		ai[0] *= fac;
		bi[0] *= fac;
		biM[0] *= fac;
		for (i=1; i<n-1; i++)
			y2[i] = bim[i]*c[i-1] + bi[i]*c[i] + biM[i]*c[i+1];
		y2[0] = zero;
		y2[n-1] = zero;
	}

	T evaluate(double xVal)
	{
		double h, a, b;
		int k,k1,k2;
		T ret;
		k1=0;
		k2=n-1;
		while (k2-k1 > 1)
		{
			k=(k2+k1) >> 1;
			if (x[k] > xVal)
				k2=k;
			else
				k1=k;
		}
		if (xVal < x[0] || xVal > x[x.size()-1])
			printf("L_Spline1DGen::evaluate() : advertencia: evaluacion fuera de la spline\n");
		throw_L_ArgException_if(k1>1 && k2<n-1 && (x[k1] > xVal || x[k2] < xVal), "L_Spline::evaluate() : error de rango");
		h=x[k2]-x[k1];
		a=(x[k2]-xVal)/h;
		b=(xVal-x[k1])/h;
		ret = a*y[k1]+b*y[k2]+((a*a*a-a)*y2[k1]+(b*b*b-b)*y2[k2])*(h*h)/6.0;
		return ret;
	}
};


class L_Spline1D
{
public:
	std::vector<double> x;  // Values for independant variable
	std::vector<double> y;  // Values from the function
	std::vector<double> y2; // Internal values (deriatives)
	int n;
	bool xEntero; // Shows if B-spline; i.e., sampling is regular (space = 1)

	L_Spline1D() {n=0; xEntero=false;}
	inline void resize(int n)
	{
		if (this->n!=n || this->n==0)
		{
#ifdef L_BASIC_DEBUG
			if (n<=0)
				L_hard_shutdown("L_Spline1D::reallocate : non-positive size allocation");
#endif
			this->n = n; x.resize(n); y.resize(n); y2.resize(n);
		}
	}
	void swap(L_Spline1D &other) {int ene = n; n = other.n; other.resize(ene); x.swap(other.x); y.swap(other.y); y2.swap(other.y2);}
	void fillIndexes() {xEntero = true; for (int i=0; i<n; i++) x[i]=i;}
	void compute();
	double evaluate(double xVal);

	static void test();
	static void pruebaGen();
};

class L_Spline2DSlow
{
public:
	L_Array<L_Spline1D> splArr;
	int lx;
	int ly;
	L_Spline1D splAux;
	double xAux;
	L_Spline2DSlow() {xAux=-1e30;}
	inline void resize(int lx, int ly)
	{
		xAux=-1e30;
		if (this->lx!=lx || this->ly!=ly)
		{
			int i;
#ifdef L_BASIC_DEBUG
			if (lx<=0||ly<=0)
				L_hard_shutdown("L_Spline2DSlow<T>::reallocate : non-positive size allocation");
#endif
			this->lx=lx; this->ly=ly;
			splArr.resize_swapping(ly);
			for (i=0; i<ly; i++)
				splArr[i].resize(lx);
			splAux.resize(ly);
		}
	}
	inline void fillIndexes() { for (int j=0; j<ly; j++) splArr[j].fillIndexes(); }
	void compute() {int j; xAux=-1e30; for (j=0; j<ly; j++) splArr[j].compute();}
	double evaluate(double x, double y);
};

class L_SigmaConvGaus
{
public:
	double sigma;
	long codeMapping;
	static int cmp(const void *a, const void *b);
};

class L_SigmasConvGausCounter
{
public:
	L_SigmaConvGaus arr[100];
	int nSigmas;
	bool activo;
	L_SigmasConvGausCounter()
	{
		nSigmas=0;
		activo=false;
	}
	void agregaSigma(double s);
	void imprCodigo(FILE *fp);
};

class L_SignalDouble
{
public:
	std::vector<double> s;
	typedef int size_type;
	long fs; // Frecuencia de muestreo (muestras/seg), debe ser entera

	L_SignalDouble() : fs(-1) {}
	void destroy() {fs=-1; s.clear();}
	void swap(L_SignalDouble &other) {s.swap(other.s); long fs_t = other.fs; other.fs=fs; fs = fs_t;}

	double &operator[](int i) {return s[i];}
	const double &operator[](int i) const {return s[i];}
	int size() const {return (int)s.size();}
	void resize(int size) {s.resize(size);}

	bool readWAV(const char *name); // resize la senal en [-1,1]
	bool writeWAV(const char *name, int bitsPerSample=16); // fs sale de las senales, resize la senal en [-1,1]

	void resample_simple(const L_SignalDouble &other, long fsNuevo);

	void gaussian(double sigma, int width) {L_calcGausHNorm(this->s, sigma, width);}

	void fourier(int dir, L_SignalDouble &imag); // size of imag must be power of 2
	void fourierByWindowing(int dir, int width, L_SignalDouble &imag); // El length debe ser potencia de 2
	bool genGraphic(L_ImageRGBUchar &im, long i1=0, long i2=-1);
	bool scanLineFromImage(const L_ImageGrayDouble &im, double i1, double j1, double i2, double j2, bool bilineal=true);
	bool normalizeHistogram(int borde=0, double valMin = 0, double valMax = 1);
	bool countMaxMin(int &nMax, int &nMin, double umbralMax, double umbralMin, int semiAnchoVentana);
	void filter_IR_right(double a); //elem[i]=a*elem[i]+(1-a)*elem[i-1]
	void filter_IR_left(double a); //elem[i]=a*elem[i]+(1-a)*elem[i+1]
	L_SignalDouble& convSim(const L_SignalDouble &other, L_SignalDouble &h, int paso = 1);
	~L_SignalDouble() {destroy();}
};

class L_Arr_SignalDouble : public L_Array<L_SignalDouble>
{
public:
	bool readWAV(const char *name);  // resizes signal in range [-1,1]
	bool writeWAV(const char *name, int bitsPerSample=16); // resizes signal in range [-1,1]
};









//! Generic matrix class with internal pointer and width step

template <class T>
class L_MatrixBaseWidthStepNonAllocator
{
public:
	typedef T value_type;
	typedef int size_type;
	typedef long size_type2;
	typedef T* iterator;
	typedef const T* const_iterator;

	T * data() {return iniptr;}
	T * row_begin(size_type i) {return (T*)((char*)iniptr+i*ljStep);}
	T * row_end(size_type i) {return ((T*)((char*)iniptr+i*ljStep))+lj;}
	const T * data() const {return iniptr;}
	const T * row_begin(size_type i) const {return (T*)((char*)iniptr+i*ljStep);}
	const T * row_end(size_type i) const {return ((T*)((char*)iniptr+i*ljStep))+lj;}

	size_type rows() const {return li;}
	size_type cols() const {return lj;}
	size_type2 size() const {return rows()*(size_type2)ljStep;}

	L_MatrixBaseWidthStepNonAllocator() : iniptr(NULL), endptr(NULL), li(0), lj(0), ljStep(0) {}

#ifdef L_IMAGE_ACCESS_DEBUG
	T& operator() (size_type i, size_type j) {throw_L_ArgException_if(iniptr == NULL || i<0 || j<0 || i>=li || j>=lj, "Matrix index out of range"); return ((T*)((char*)iniptr + ljStep*i))[j];}
	const T& operator() (size_type i, size_type j) const {throw_L_ArgException_if(iniptr == NULL || i<0 || j<0 || i>=li || j>=lj, "Matrix index out of range"); return ((T*)((char*)iniptr + ljStep*i))[j];}
#else
	T& operator() (size_type i, size_type j) {return ((T*)((char*)iniptr + ljStep*i))[j];}
	const T& operator() (size_type i, size_type j) const {return ((T*)((char*)iniptr + ljStep*i))[j];}
#endif

protected:
	T *iniptr;
	T *endptr;
public:
	union {size_type li; size_type ly;};
	union {size_type lj; size_type lx;};
	union {size_type lxStep; size_type ljStep;};  // For internally allocated memory, ljStep = lj*sizeof(T)
};



// Matrix with allocator and width step

template <class T, bool POD = true>  // T must be a POD object
class L_MatrixBaseWidthStep : public L_MatrixBaseWidthStepNonAllocator<T>
{
public:
	typedef typename L_MatrixBaseWidthStepNonAllocator<T>::size_type size_type;
	using L_MatrixBaseWidthStepNonAllocator<T>::data;
	using L_MatrixBaseWidthStepNonAllocator<T>::iniptr;
	using L_MatrixBaseWidthStepNonAllocator<T>::endptr;
	using L_MatrixBaseWidthStepNonAllocator<T>::li;
	using L_MatrixBaseWidthStepNonAllocator<T>::ly;
	using L_MatrixBaseWidthStepNonAllocator<T>::lj;
	using L_MatrixBaseWidthStepNonAllocator<T>::lx;
	using L_MatrixBaseWidthStepNonAllocator<T>::lxStep;
	using L_MatrixBaseWidthStepNonAllocator<T>::ljStep;
	using L_MatrixBaseWidthStepNonAllocator<T>::operator();


private:
	L_MatrixBaseWidthStep(const L_MatrixBaseWidthStep &other);
	void operator=(const L_MatrixBaseWidthStep &other);
	bool mem_ext; // For externally allocated memory, "use with care"

public:
	L_MatrixBaseWidthStep() {mem_ext = false;}
	L_MatrixBaseWidthStep(size_type li, size_type lj, const T *arr)
	{
		L_StaticAssert(POD == true, L_MatrixBaseWidthStep__applies__only__in__POD__objects);
		mem_ext = false;
		reallocate_fc(li, lj);  // lj = ljStep
		memcpy(data(), arr, sizeof(T)*li*lj);
	}
	L_MatrixBaseWidthStep(size_type li, size_type lj) {mem_ext = false; reallocate_fc(li, lj);}

	// Siempre deja ljStep = lj
	void reallocate_fc(size_type li2, size_type lj2)
	{
		T *buf = NULL;

		throw_L_ArgException_if(li2 == 0 || lj2 == 0, "L_MatrixBaseWidthStep::reallocate_fc() : zero size allocation");

		int ljStep2 = lj2*sizeof(T);
		if (iniptr == NULL)
		{
			mem_ext = false;
			if (li2>0 && lj2 > 0)
				buf = (T*)L_aligned_malloc_128(li2*ljStep2);
			iniptr = buf;
			endptr = (T*)((char *)iniptr + li2*ljStep2);
		}
		else if (li!=li2 || lj!=lj2 || ljStep != ljStep2) // Error: decia ljStep != lj
		{
			throw_L_ArgException_if(mem_ext == true && (li != li2 || lj != lj2), "L_MatrixBaseWidthStep::reallocate_fc() : intento de redimensionar matriz con memoria externa");
			if (mem_ext == false)
			{
				mem_ext = false;
				if (li2>0 && lj2 > 0)
					buf = (T*)L_aligned_malloc_128(li2*ljStep2);
				destroy();
				iniptr = buf;
				endptr = (T*)((char *)iniptr + li2*ljStep2);
			}
		}
		li = li2;
		lj = lj2;
		ljStep = ljStep2;
		return;
	}

	void reallocate_fc_ext(size_type li2, size_type lj2, T *buf, size_type ljStep2)  // Siempre asign la memoria
	{
		if (iniptr == NULL)
		{
			mem_ext = true;
			iniptr = buf;
			endptr = (T*)((char *)iniptr + li2*ljStep2);
		}
		else if (li != li2 || lj != lj2 || (ljStep2 != 0 && ljStep2 != ljStep) || buf != iniptr)
		{
			if (mem_ext)
				destroy_ext();
			else
				destroy();
			mem_ext = true;
			iniptr = buf;
			endptr = (T*)((char *)iniptr + li2*ljStep2);
		}
		li = li2;
		lj = lj2;
		ljStep = ljStep2;
		return;
	}

	// Siempre deja lxStep = lx
	void reallocate_cf(size_type lx, size_type ly) {reallocate_fc(ly, lx);}
	void reallocate_cf_ext(size_type lx, size_type ly, T *buf, size_type lx_step) {reallocate_fc_ext(ly, lx, buf, lx_step);} // Siempre asign la memoria

#ifdef __COMPAT_IPLIMAGE__
	void reallocate_cf_ext_ipl(IplImage *img) {throw_L_ArgException_if(img == NULL, "reallocate_cf_ext_ipl : IplImage is NULL"); reallocate_cf_ext(img->width, img->height, img->imageData, img->widthStep);}
#endif

	void conservativeResize_fc(size_type li2, size_type lj2, T zero)
	{
		L_MatrixBaseWidthStep<T> m2;
		m2.reallocate_fc(li2, lj2);
		if (li2 > li || lj2 > lj)
			m2.setZero();
		size_type lix, ljx;
		lix = li;
		if (lix > li2)
			lix = li2;
		ljx = lj;
		if (ljx > lj2)
			ljx = lj2;
		for (size_type i=0; i<lix; i++)
			memcpy(&m2(i,0), &operator()(i,0), sizeof(T)*ljx);
		m2.swap(*this);
	}

	void conservativeResize_cf(size_type lx2, size_type ly2, T zero) {conservativeResize_fc(ly2, lx2, zero);}

	void setConstant(T val) {for (size_type i=0; i<li; i++) {T *p, *pbeg = &operator()(i,0), *pend = &operator()(i,0)+lj; for (p=pbeg; p<pend; ) *p++ = val;}}
	void setZero() {for (size_type i=0; i<li; i++) {T *p, *pbeg = &operator()(i,0), *pend = &operator()(i,0)+lj; for (p=pbeg; p<pend; ) *p++ = 0;}}

	void destroy()
	{
		throw_L_ArgException_if(mem_ext == true, "L_MatrixBaseWidthStep::destroy() : intento de destruir matriz con memoria externa");
		li = 0;
		lj = 0;
		if (iniptr != NULL)
		{
			L_aligned_free_128(iniptr);
			iniptr = NULL;
			endptr = NULL;
		}
	}

	void destroy_ext()
	{
		throw_L_ArgException_if(mem_ext == false, "L_MatrixBaseWidthStep::destroy_ext() : allocated memory is internal");
		li = 0;
		lj = 0;
		if (iniptr != NULL)
		{
			iniptr = NULL;
			endptr = NULL;
		}
		mem_ext = false;
	}

	size_type mem() {return li*ljStep;}

	void swap(L_MatrixBaseWidthStep<T> &other)
	{
		throw_L_ArgException_if (mem_ext == true || other.mem_ext == true, "L_MatrixBaseWidthStep::swap() on external memory");
		bool mem_ext_t; // External memory, use with care
		size_type li_t;
		size_type lj_t;
		size_type lxStep_t;
		T *iniptr_t;
		T *endptr_t;

		mem_ext_t = mem_ext;
		li_t = li;
		lj_t = lj;
		lxStep_t = lxStep;
		iniptr_t = iniptr;
		endptr_t = endptr;
		
		mem_ext = other.mem_ext;
		li = other.li;
		lj = other.lj;
		lxStep = other.lxStep;
		iniptr = other.iniptr;
		endptr = other.endptr;

		other.mem_ext = mem_ext_t;
		other.li = li_t;
		other.lj = lj_t;
		other.lxStep = lxStep_t;
		other.iniptr = iniptr_t;
		other.endptr = endptr_t;
	}

	void copyTo(L_MatrixBaseWidthStep<T> &other) const
	{
		other.reallocate_fc(li, lj);
		if (ljStep == other.ljStep)
			memcpy(other.data(), data(), li*ljStep);
		else
		{
			for (size_type i=0; i<li; i++)
				memcpy(&other(i,0), &operator()(i,0), sizeof(T)*lj);
		}
	}

	bool copyReferenceIn(L_MatrixBaseWidthStep<T> &other)
	{
		other.reallocate_fc_ext(li, lj, &this->_elem[0][0], lj);
		return true;
	}

	void OP_assign(const L_MatrixBaseWidthStep<T> &other)
	{
		reallocate_fc(other.li, other.lj);
		if (other.ljStep == ljStep)
			memcpy(data(), other.data(), li*ljStep);
		else
		{
			for (size_type i=0; i<li; i++)
				memcpy(&operator()(i,0), &other(i,0), sizeof(T)*lj);
		}
	}

	void copyToBuffer(T *buf, size_type lj_step = -1) const
	{
		if (lj_step == -1)
			lj_step = lj*sizeof(T);
		if (ljStep == lj_step)
			memcpy(buf, data(), li*lj_step);
		else
			for (size_type i=0; i<li; i++)
				memcpy((char*)buf+i*lj_step, &operator()(i,0), sizeof(T)*lj);
	}

	void copyFromBuffer(const T *buf, size_type lj_step = -1)
	{
		if (lj_step == -1)
			lj_step = lj*sizeof(T);
		if (lj_step == ljStep)
			memcpy(data(), buf, li*ljStep);
		else
			for (size_type i=0; i<li; i++)
				memcpy(&operator()(i,0), (char *)buf+i*ljStep, sizeof(T)*lj);
	}

	~L_MatrixBaseWidthStep() {destroy();}
};



// Matrix with allocator and width step

template <class T>
class L_MatrixBaseNonAllocator
{
public:
	typedef T value_type;
	typedef int size_type;
	typedef long size_type2;
	typedef T* iterator;
	typedef const T* const_iterator;

	T * data() {return iniptr;}
	T * begin() {return iniptr;}
	T * end() {return endptr;}
	const T * data() const {return iniptr;}
	const T * begin() const {return iniptr;}
	const T * end() const {return endptr;}

	size_type rows() const {return li;}
	size_type cols() const {return lj;}
	size_type2 size() const {return rows()*(size_type2)cols();}

	L_MatrixBaseNonAllocator() : iniptr(NULL), endptr(NULL), li(0), lj(0) {}

#ifdef L_MATRIX_ACCESS_DEBUG
	T& operator() (size_type i, size_type j) {throw_L_ArgException_if(iniptr == NULL || i<0 || j<0 || i>=li || j>=lj, "Matrix index out of range"); return (iniptr + lj*i)[j];}
	const T& operator() (size_type i, size_type j) const {throw_L_ArgException_if(iniptr == NULL || i<0 || j<0 || i>=li || j>=lj, "Matrix index out of range"); return (iniptr + lj*i)[j];}
#else
	T& operator() (size_type i, size_type j) {return (iniptr + lj*i)[j];}
	const T& operator() (size_type i, size_type j) const {return (iniptr + lj*i)[j];}
#endif

protected:
	T *iniptr;
	T *endptr;
public:
	size_type li;
	size_type lj;
};



// Matrix with allocator

template <class T, bool POD = true>  // T must be a POD object
class L_MatrixBase : public L_MatrixBaseNonAllocator<T>
{
public:
  typedef typename L_MatrixBaseNonAllocator<T>::value_type value_type;
  typedef typename L_MatrixBaseNonAllocator<T>::size_type size_type;
  typedef typename L_MatrixBaseNonAllocator<T>::size_type2 size_type2;
  typedef typename L_MatrixBaseNonAllocator<T>::iterator iterator;
  typedef typename L_MatrixBaseNonAllocator<T>::const_iterator const_iterator;
  using L_MatrixBaseNonAllocator<T>::data;
  using L_MatrixBaseNonAllocator<T>::begin;
  using L_MatrixBaseNonAllocator<T>::end;
  using L_MatrixBaseNonAllocator<T>::rows;
  using L_MatrixBaseNonAllocator<T>::cols;
  using L_MatrixBaseNonAllocator<T>::size;
  using L_MatrixBaseNonAllocator<T>::iniptr;
  using L_MatrixBaseNonAllocator<T>::endptr;
  using L_MatrixBaseNonAllocator<T>::li;
  using L_MatrixBaseNonAllocator<T>::lj;
  using L_MatrixBaseNonAllocator<T>::operator();

private:
	L_MatrixBase(const L_MatrixBase &other);
	void operator=(const L_MatrixBase &other);
	bool mem_ext; // For externally allocated memory, "use with care"

public:
	L_MatrixBase() {mem_ext = false;}
	L_MatrixBase(size_type li, size_type lj, const T *arr)
	{
		mem_ext = false;
		reallocate_fc(li, lj);
		memcpy(begin(), arr, sizeof(T)*li*lj);
	}
	L_MatrixBase(size_type li, size_type lj) {mem_ext = false; reallocate_fc(li, lj);}

	// Siempre deja ljStep = lj
	void reallocate_fc(size_type li2, size_type lj2)
	{
		T *buf = NULL;

		throw_L_ArgException_if(li2<0 || lj2<0, "L_MatrixBase::reallocate_fc() : tamano negativo");
		if (iniptr == NULL)
		{
			mem_ext = false;
			if (li2 > 0 && lj2 > 0)
			{
				if (POD)
					buf = (T*)L_aligned_malloc_128(sizeof(T)*li2*lj2);
				else
					buf = new T[li2*lj2];
			}
			iniptr = buf;
			endptr = (iniptr + li2*lj2);
		}
		else if (li!=li2 || lj!=lj2)
		{
			throw_L_ArgException_if(mem_ext == true && (li != li2 || lj != lj2), "L_MatrixBase::reallocate_fc() : intento de redimensionar matriz con memoria externa");
			if (mem_ext == false)
			{
				mem_ext = false;
				if (li2 > 0 && lj2 > 0)
				{
					if (POD)
						buf = (T*)L_aligned_malloc_128(sizeof(T)*li2*lj2);
					else
						buf = new T[li2*lj2];
				}
				destroy();
				iniptr = buf;
				endptr = (iniptr + li2*lj2);
			}
		}
		li = li2;
		lj = lj2;
		return;
	}

	void reallocate_fc_ext(size_type li2, size_type lj2, T *buf)  // Siempre asign la memoria
	{
		throw_L_ArgException_if(li2<0 || lj2<0, "L_MatrixBase::reallocate_fc_ext() : tamano negativo");
		if (iniptr == NULL)
		{
			mem_ext = true;
			iniptr = buf;
			endptr = (iniptr + li2*lj2);
		}
		else if (li != li2 || lj != lj2 || buf != iniptr)
		{
			if (mem_ext)
				destroy_ext();
			else
				destroy();
			mem_ext = true;
			iniptr = buf;
			endptr = (iniptr + li2*lj2);
		}
		li=li2;
		lj=lj2;
		return;
	}

	// Siempre deja lxStep = lx

	void conservativeResize_fc(size_type li2, size_type lj2, T zero)
	{
		L_MatrixBase<T> m2;
		m2.reallocate_fc(li2, lj2);
		if (li2 > li || lj2 > lj)
			m2.setZero();
		size_type lix, ljx;
		lix = li;
		if (lix > li2)
			lix = li2;
		ljx = lj;
		if (ljx > lj2)
			ljx = lj2;
		for (size_type i=0; i<lix; i++)
			memcpy(&m2(i,0), &operator()(i,0), sizeof(T)*ljx);
		m2.swap(*this);
	}

	void conservativeResize_cf(size_type lx2, size_type ly2, T zero) {conservativeResize_fc(ly2, lx2, zero);}

	void setConstant(T val) {for (size_type i=0; i<li; i++) {T *p, *pbeg = &operator()(i,0), *pend = &operator()(i,0)+lj; for (p=pbeg; p<pend; ) *p++ = val;}}
	void setZero() {for (size_type i=0; i<li; i++) {T *p, *pbeg = &operator()(i,0), *pend = &operator()(i,0)+lj; for (p=pbeg; p<pend; ) *p++ = 0;}}

	void destroy()
	{
		throw_L_ArgException_if(mem_ext == true, "L_MatrixBase::destroy() : intento de destruir matriz con memoria externa");
		li = 0;
		lj = 0;
		if (iniptr != NULL)
		{
			if (POD)
				L_aligned_free_128(iniptr);
			else
				delete[] iniptr;
			iniptr = NULL;
			endptr = NULL;
		}
	}

	void destroy_ext()
	{
		throw_L_ArgException_if(mem_ext == false, "L_MatrixBase::destroy_ext() : allocated memory is internal");
		li = 0;
		lj = 0;
		if (iniptr != NULL)
		{
			iniptr = NULL;
			endptr = NULL;
		}
		mem_ext = false;
	}

	size_type mem() {return li*lj*sizeof(T);}

	void swap(L_MatrixBase<T> &other)
	{
		throw_L_ArgException_if (mem_ext == true || other.mem_ext == true, "L_MatrixBase::swap() on external memory");
		bool mem_ext_t; // External memory, use with care
		size_type li_t;
		size_type lj_t;
		T *iniptr_t;
		T *endptr_t;

		mem_ext_t = mem_ext;
		li_t = li;
		lj_t = lj;
		iniptr_t = iniptr;
		endptr_t = endptr;
		
		mem_ext = other.mem_ext;
		li = other.li;
		lj = other.lj;
		iniptr = other.iniptr;
		endptr = other.endptr;

		other.mem_ext = mem_ext_t;
		other.li = li_t;
		other.lj = lj_t;
		other.iniptr = iniptr_t;
		other.endptr = endptr_t;
	}

	void copyTo(L_MatrixBase<T> &other) const
	{
		other.reallocate_fc(li, lj);
		if (this->ljStep == other.ljStep)
			memcpy(other.begin(), begin(), li*lj*sizeof(T));
		else
		{
			for (size_type i=0; i<li; i++)
				memcpy(&other(i,0), &operator()(i,0), sizeof(T)*lj);
		}
	}

	bool copyReferenceIn(L_MatrixBase<T> &other)
	{
		other.reallocate_fc_ext(li, lj, begin(), lj);
		return true;
	}

	void OP_assign(const L_MatrixBase<T> &other)
	{
		reallocate_fc(other.li, other.lj);
		if (other.ljStep == this->ljStep)
			memcpy(begin(), other.begin(), li*lj*sizeof(T));
		else
		{
			for (size_type i=0; i<li; i++)
				memcpy(&operator()(i,0), &other(i,0), sizeof(T)*lj);
		}
	}

	void copyToBuffer(T *buf) const
	{
		memcpy(buf, begin(), li*lj*sizeof(T));
	}

	void copyFromBuffer(const T *buf)
	{
		memcpy(begin(), buf, li*lj*sizeof(T));
	}

	~L_MatrixBase() {destroy();}
};

// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_ImageGrayFloat_isMax2D_c(I,cen,i,j) (((cen)>(I).pix((i)+1,(j)+1) && (cen)>(I).pix((i)+1,(j)) && (cen)>(I).pix((i)+1,(j)-1) &&(cen)>(I).pix((i),(j)+1) && (cen)>(I).pix((i),(j)) && (cen)>(I).pix((i),(j)-1) &&(cen)>(I).pix((i)-1,j+1) && (cen)>(I).pix(i-1,j) && (cen)>(I).pix(i-1,j-1)))
#define L_ImageGrayFloat_isMin2D_c(I,cen,i,j) (((cen)<(I).pix((i)+1,(j)+1) && (cen)<(I).pix((i)+1,(j)) && (cen)<(I).pix((i)+1,(j)-1) && (cen)<(I).pix((i),(j)+1) && (cen)<(I).pix((i),(j)) && (cen)<(I).pix((i),(j)-1) && (cen)<(I).pix((i)-1,j+1) && (cen)<(I).pix(i-1,j) && (cen)<(I).pix(i-1,j-1)))

#define L_ImageGrayFloat_isMax2D(I,i,j) ((((I).pix(i,j))>(I).pix((i)+1,(j)+1) && ((I).pix(i,j))>(I).pix((i)+1,(j)) && ((I).pix(i,j))>(I).pix((i)+1,(j)-1) && ((I).pix(i,j))>(I).pix((i),(j)+1) && ((I).pix(i,j))>(I).pix((i),(j)-1) &&((I).pix(i,j))>(I).pix((i)-1,j+1) && ((I).pix(i,j))>(I).pix(i-1,j) && ((I).pix(i,j))>(I).pix(i-1,j-1)))
#define L_ImageGrayFloat_isMin2D(I,i,j) ((((I).pix(i,j))<(I).pix((i)+1,(j)+1) && ((I).pix(i,j))<(I).pix((i)+1,(j)) && ((I).pix(i,j))<(I).pix((i)+1,(j)-1) && ((I).pix(i,j))<(I).pix((i),(j)+1) && ((I).pix(i,j))<(I).pix((i),(j)-1) && ((I).pix(i,j))<(I).pix((i)-1,j+1) && ((I).pix(i,j))<(I).pix(i-1,j) && ((I).pix(i,j))<(I).pix(i-1,j-1)))

#define L_ImageGrayFloat_rhoHough(theta,x,y) (x*cos(theta)+y*sin(theta))

// Estas funciones no piden memoria al sistema
#define L_ImageGrayFloat_OP_add(M, A, B) {int i63, j73; throw_L_ArgException_if((M).data()==NULL || (A).data()==NULL || (B).data()==NULL || (M).lx!=(A).lx || (M).ly!=(A).ly || (M).lx!=(B).lx || (M).ly!=(A).ly, "L_ImageGrayFloat_OP_add"); for (j73=0; j73<(M).ly; j73++) for (i63=0; i63<(M).lx; i63++) (M)(j73,i63) = (A)(j73,i63) + (B)(j73,i63);}
#define L_ImageGrayFloat_OP_subtract(M, A, B) {int i62, j72; throw_L_ArgException_if((M).data()==NULL || (A).data()==NULL || (B).data()==NULL || (M).lx!=(A).lx || (M).ly!=(A).ly || (M).lx!=(B).lx || (M).ly!=(A).ly, "L_ImageGrayFloat_OP_subtract"); for (j72=0; j72<(M).ly; j72++) for (i62=0; i62<(M).lx; i62++) (M)(j72,i62) = (A)(j72,i62) - (B)(j72,i62);}

class L_ImageGrayFloat : public L_MatrixBaseWidthStep<float>
{
public:
  typedef typename L_MatrixBaseWidthStep<float>::size_type size_type;
  using L_MatrixBaseWidthStep<float>::data;
  using L_MatrixBaseWidthStep<float>::iniptr;
  using L_MatrixBaseWidthStep<float>::endptr;
  using L_MatrixBaseWidthStep<float>::li;
  using L_MatrixBaseWidthStep<float>::ly;
  using L_MatrixBaseWidthStep<float>::lj;
  using L_MatrixBaseWidthStep<float>::lx;
  using L_MatrixBaseWidthStep<float>::lxStep;
  using L_MatrixBaseWidthStep<float>::ljStep;
  using L_MatrixBaseWidthStep<float>::operator();

private:
	L_ImageGrayFloat(const L_ImageGrayFloat &other);
public:
	L_ImageGrayFloat() : L_MatrixBaseWidthStep<float>() {}
	L_ImageGrayFloat(int lx, int ly) : L_MatrixBaseWidthStep<float>(ly, lx) {}
	L_ImageGrayFloat(int lx, int ly, const float *buf) : L_MatrixBaseWidthStep<float>(ly, lx, buf) {}

	inline float &pix(int c, int f) {return operator()(f,c);}
	inline const float &pix(int c, int f) const {return operator()(f,c);}

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}

	void operator=(const L_ImageGrayFloat &other) {other.copyTo(*this);}

	void operator+=(const L_ImageGrayFloat &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator+= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) += other(i,j);}
	void operator-=(const L_ImageGrayFloat &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator-= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) -= other(i,j);}
	void OP_add(const L_ImageGrayFloat &im1, const L_ImageGrayFloat &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayFloat::OP_add"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i<this->li; i++) for (j=0; j<this->lj; j++) operator()(i,j) = im1(i,j) + im2(i,j);}
	void OP_subtract(const L_ImageGrayFloat &im1, const L_ImageGrayFloat &im2) { throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayFloat::OP_subtract"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) - im2(i,j);}
	void OP_mult_elementwise(const L_ImageGrayFloat &im1, const L_ImageGrayFloat &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayFloat::OP_mult_elementwise"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) * im2(i,j);}
	void transpOf(const L_ImageGrayFloat &other) {int i, j; reallocate_fc(other.lj, other.li); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(j,i) = other(i,j);}

	static float channelY(float R, float G, float B) {return +(float)0.299*R +(float)0.587*G +(float)0.114*B;}
	static float channelU(float R, float G, float B) {return -(float)0.168*R -(float)0.330*G +(float)0.498*B;}
	static float channelV(float R, float G, float B) {return +(float)0.498*R -(float)0.417*G -(float)0.081*B;}

	static float channelR(float Y, float U, float V) {return (float)1.0*Y +(float)0.000*U +(float)1.4075*V;}
	static float channelG(float Y, float U, float V) {return (float)1.0*Y -(float)0.3455*U -(float)0.7169*V;}
	static float channelB(float Y, float U, float V) {return (float)1.0*Y +(float)1.779*U +(float)0.000*V;}
};

class L_ImageGrayDouble : public L_MatrixBaseWidthStep<double>
{
private:
	L_ImageGrayDouble(const L_ImageGrayDouble &other);
public:
	L_ImageGrayDouble() : L_MatrixBaseWidthStep<double>() {}
	L_ImageGrayDouble(int lx, int ly) : L_MatrixBaseWidthStep<double>(ly, lx) {}
	L_ImageGrayDouble(int lx, int ly, const double *buf) : L_MatrixBaseWidthStep<double>(ly, lx, buf) {}

	inline double &pix(int c, int f) {return operator()(f,c);}
	inline const double &pix(int c, int f) const {return operator()(f,c);}

	inline const double &pixc(int c, int f) const {return operator()(f,c);} // const

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}

	int readBMP(const char *name);
	bool writeBMP(const char *name, short nb);
	bool writeBMP(const char *name);
	bool writeBMP_RB(const char *name, int marg=5);

	bool readImage(const char *name);
	bool saveImage(const char *name);

#ifdef __COMPAT_ATLIMAGE__
	int readImageATL(const char *name);
	bool writeImageATL(const char *name);
#endif
    bool normalizeHistogram(int borde=0);

	L_ImageGrayDouble& genChannelY(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB);
	L_ImageGrayDouble& genChannelU(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB);
	L_ImageGrayDouble& genChannelV(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB);

	L_ImageGrayDouble& genChannelR(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV);
	L_ImageGrayDouble& genChannelG(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV);
	L_ImageGrayDouble& genChannelB(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV);


	void operator+=(const L_ImageGrayDouble &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator+= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) += other(i,j);}
	void operator-=(const L_ImageGrayDouble &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator-= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) -= other(i,j);}
	void OP_add(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayDouble::OP_add"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i<this->li; i++) for (j=0; j<this->lj; j++) operator()(i,j) = im1(i,j) + im2(i,j);}
	void OP_subtract(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2) { throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayDouble::OP_subtract"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) - im2(i,j);}
	void OP_mult_elementwise(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayDouble::OP_mult_elementwise"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) * im2(i,j);}
	void transpOf(const L_ImageGrayDouble &other) {int i, j; reallocate_fc(other.lj, other.li); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(j,i) = other(i,j);}
	void transposeMe() {L_ImageGrayDouble other; other.transpOf(*this); other.swap(*this);}

	void multiply_to_each_element(double factor) {int c, f; for (f=0; f<ly; f++) for (c=0; c<lx; c++) pix(c,f)*=(double)factor;}
	void sum_to_each_element(double ilum) {int c, f; for (f=0; f<ly; f++) for (c=0; c<lx; c++) pix(c,f)+=(double)ilum;}

	L_ImageGrayDouble& OP_mult_elementwise_sub(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2, int x0, int y0, int rx, int ry);
	L_ImageGrayDouble& OP_div_elementwise(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2); // im2 no debe tener pixeles 0

	static void smoothXSqrt2(const L_ImageGrayDouble &imEntr, L_ImageGrayDouble &res);
	static void smoothYSqrt2(const L_ImageGrayDouble &imEntr, L_ImageGrayDouble &res);
	
	L_ImageGrayDouble& smooth(const L_ImageGrayDouble &im1, double s, int paso=1);
	L_ImageGrayDouble& smooth_sub(L_ImageGrayDouble &imTmp, const L_ImageGrayDouble &im1, double s, int x0, int y0, int rx1, int ry1, int rx2, int ry2);
	L_ImageGrayDouble& smooth_changeSigma(const L_ImageGrayDouble &im1, double sInic, double sFinal, double sMinSubm=100.0, FILE *fp=NULL);
	L_ImageGrayDouble& smoothXYSqrt2(const L_ImageGrayDouble &im1) {L_ImageGrayDouble a; smoothXSqrt2(im1,a); smoothYSqrt2(a,*this); return *this;}
	L_ImageGrayDouble& z_smooth0Tr(const L_ImageGrayDouble &im1, double s);
	L_ImageGrayDouble& z_smooth1Tr(const L_ImageGrayDouble &im1, double s);
	L_ImageGrayDouble& z_smooth2Tr(const L_ImageGrayDouble &im1, double s);
	L_ImageGrayDouble& z_smooth0Tr_sub(L_ImageGrayDouble &imTmp, const L_ImageGrayDouble &im1, double s, int x0, int y0, int rx1, int ry1, int rx2, int ry2);

	L_ImageGrayDouble& filterGabor(const L_ImageGrayDouble &other, L_ImageGrayDouble *gaborIm, double s, double lambda, double ang, double diametroEnSigmas=-1.0);
	L_ImageGrayDouble& monogenicDecomposition_slow(int semiancho, L_ImageGrayDouble &I_Cos, L_ImageGrayDouble &I_SinCos, L_ImageGrayDouble &I_SinSin); // I_Cos es la miama señal original. I_Sin(fase)Cos(ori)
	void fourier(int dir, L_ImageGrayDouble &imag); // lx and ly must be power of 2
	void fftshift() {double d; int c, f, L2 = lx - lx/2; for (f=0; f<ly; f++) for (c=0; c<lx/2; c++) {d = pix(c+L2,(f+ly/2)%ly); pix(c+L2,(f+ly/2)%ly) = pix(c,f); pix(c,f) = d;}}

	void houghDetPuntosBordes(const L_ImageGrayDouble &imorig, L_ImageGrayUchar &binariz, std::vector<double> &px, std::vector<double> &py, std::vector<double> &ang, double umbral, bool supresionNoMaximo);
	void houghDetLineas(const L_ImageGrayDouble &imorig, std::vector<double> &px, std::vector<double> &py, std::vector<double> &ang, int nThetas, int nRhos, bool usaGrad);
	void houghDibLineas(const L_ImageGrayDouble &imorig, L_ShapeArray &lins, double umbral);

	L_ImageGrayDouble& conv2d(const L_ImageGrayDouble &other, const L_ImageGrayDouble &mask, int xCenMasc, int yCenMasc); // La mask se le pasa ya dada vuelta en 180°
	L_ImageGrayDouble& convXSim(const L_ImageGrayDouble &other, const std::vector<double> &h);
    L_ImageGrayDouble& convYSim(const L_ImageGrayDouble &other, const std::vector<double> &h);
    L_ImageGrayDouble& convYSimTr(const L_ImageGrayDouble &other, const std::vector<double> &h);

	L_ImageGrayDouble& convXSim_sub(const L_ImageGrayDouble &other, const std::vector<double> &h, int x0, int y0, int rx, int ry);
    L_ImageGrayDouble& convYSim_sub(const L_ImageGrayDouble &other, const std::vector<double> &h, int x0, int y0, int rx, int ry);

	L_ImageGrayDouble& z_smooth0Tr_step(const L_ImageGrayDouble &im1, double s, int paso);
    L_ImageGrayDouble& z_smoothTr_step(const L_ImageGrayDouble &im1, double s, int paso);
    L_ImageGrayDouble& z_smooth2Tr_step(const L_ImageGrayDouble &im1, double s, int paso);

	L_ImageGrayDouble& convXSim_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso);
    L_ImageGrayDouble& convYSim_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso);
    L_ImageGrayDouble& convYSimTr_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso);

    L_ImageGrayDouble& gradRight_3(const L_ImageGrayDouble &other);
    L_ImageGrayDouble& gradUp_3(const L_ImageGrayDouble &other);
    L_ImageGrayDouble& gradRight_2(const L_ImageGrayDouble &other);
    L_ImageGrayDouble& gradUp_2(const L_ImageGrayDouble &other);
    L_ImageGrayDouble& gradRightUp_2(L_ImageGrayDouble &gradArr, const L_ImageGrayDouble &other);
	L_ImageGrayDouble& gradR_2(L_ImageGrayDouble &gradAng, const L_ImageGrayDouble &other); // Gradiente en polares. gradAng=atan2(gradArr,gradDer);
	L_ImageGrayDouble& gradR_3(L_ImageGrayDouble &gradAng, const L_ImageGrayDouble &other); // Gradiente en polares. gradAng=atan2(gradArr,gradDer);
    L_ImageGrayDouble& gradRight_3_sub(const L_ImageGrayDouble &other, int x0, int y0, int rx, int ry);
    L_ImageGrayDouble& gradUp_3_sub(const L_ImageGrayDouble &other, int x0, int y0, int rx, int ry);
	void integralOf(const L_ImageGrayDouble &other);
	L_ImageGrayDouble& integralDebugOf(const L_ImageGrayDouble &other) {reallocate(other.lx, other.ly); int c, f, u, v; for (f=0; f<ly; f++) for (c=0; c<lx; c++) {pix(c,f) = 0; for (v=0; v<=f; v++) for (u=0; u<=c; u++) pix(c,f) += other.pix(u,v);} return *this;}
	double sum_rectangle_using_integral(int c1, int f1, int c2, int f2) {throw_L_ArgException_if(c1 > c2 || f1 > f2 || c1==0 || f1==0, "L_ImageGrayDouble::sum_rectangle_using_integral"); return pix(c2,f2) - pix(c1-1,f2) - pix(c2,f1-1) + pix(c1-1,f1-1);}
	double sum_rectangle_pixelwise(int c1, int f1, int c2, int f2) {double d=0; int c, f; for (f=f1; f<=f2; f++) for (c=c1; c<=c2; c++) d+=pix(c,f); return d;}


	L_ImageGrayDouble& laplac(const L_ImageGrayDouble &other);

    L_ImageGrayDouble& z_subm50(const L_ImageGrayDouble &other);
    L_ImageGrayDouble& z_subm67(const L_ImageGrayDouble &other);
	inline L_ImageGrayDouble& subm(const L_ImageGrayDouble &other)
	{
		return z_subm50(other);
	}
	inline L_ImageGrayDouble& subm(const L_ImageGrayDouble &other, L_Frac &f)
	{
		if (f.num==1 && f.den==2)
			return z_subm50(other);
		else if (f.num==2 && f.den==3)
			return z_subm67(other);
		else
			throw_L_ArgException_if(true,"L_ImageGrayDouble::subm");
		return *this;
	}
	void doubleResolution(const L_ImageGrayDouble &other);
    void doubleResolutionBilineal(const L_ImageGrayDouble &other);
	void halfResolution(const L_ImageGrayDouble &other);
	void halfResolution_noResample(const L_ImageGrayDouble &other);
	void subImageOf(const L_ImageGrayDouble &imOrig, int x0, int y0, int x1, int y1);

	void AplicaTransfAfinPixelado(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano=true); // Debe usarse antes un filtro antialiasing
	void AplicaTransfSemejanzaPixelado(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano=true); // Debe usarse antes un filtro antialiasing
	bool AplicaTransfProyectivaPixelado(const L_ImageGrayDouble &im, const double *coef_8); // Debe usarse antes un filtro antialiasing

	void AplicaTransfAfinBilineal(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano=true); // Debe usarse antes un filtro antialiasing
	void AplicaTransfSemejanzaBilineal(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano=true); // Debe usarse antes un filtro antialiasing

	void doubleResolution() {L_ImageGrayDouble im; im.doubleResolution(*this); im.swap(*this);}
	void halfResolution() {L_ImageGrayDouble im; im.halfResolution(*this); im.swap(*this);}
	void halfResolution_noResample() {L_ImageGrayDouble im; im.halfResolution_noResample(*this); im.swap(*this);}

	void calcSpline2D(L_Spline2DSlow &spline) const;
	void cambiaTamanoSpline(const L_ImageGrayDouble &im, double factor); // Debe usarse antes un filtro antialiasing
	void cambiaTamanoPixelado(const L_ImageGrayDouble &im, int nuevolx, int nuevoly); // Sin antialiasing
	void AplicaTransfAfinSpline(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano=true); // Debe usarse antes un filtro antialiasing
	void AplicaTransfSemejanzaSpline(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano=true); // Debe usarse antes un filtro antialiasing

	bool calculaHistGradiente(int lu, int lv, int lw, double ***hist);

	inline bool esMax2D(int i, int j)
	{
		double & cen=pix(i,j);
		return (cen>pix(i+1,j+1) && cen>pix(i+1,j) && cen>pix(i+1,j-1) &&
	    cen>pix(i,j+1) && cen>pix(i,j-1) &&
		cen>pix(i-1,j+1) && cen>pix(i-1,j) && cen>pix(i-1,j-1)
    	);
    }
	inline bool esMax2D(int i, int j, double porc) // idea: porc < 1
	{
		double cen=pix(i,j)*(double)porc;
		return (cen>pix(i+1,j+1) && cen>pix(i+1,j) && cen>pix(i+1,j-1) &&
	    cen>pix(i,j+1) && cen>pix(i,j-1) &&
		cen>pix(i-1,j+1) && cen>pix(i-1,j) && cen>pix(i-1,j-1)
    	);
    }
	bool esMax2D_anillo(double val, int i, int j, double r);
	inline bool esMin2D(int i, int j)
	{
		double & cen=pix(i,j);
		return (cen<pix(i+1,j+1) && cen<pix(i+1,j) && cen<pix(i+1,j-1) &&
		    cen<pix(i,j+1) && cen<pix(i,j-1) &&
			cen<pix(i-1,j+1) && cen<pix(i-1,j) && cen<pix(i-1,j-1)
	    );
	}
	inline bool esMin2D(int i, int j, double porc) // idea: porc < 1
	{
		double cen=pix(i,j)*(double)porc;
		return (cen<pix(i+1,j+1) && cen<pix(i+1,j) && cen<pix(i+1,j-1) &&
		    cen<pix(i,j+1) && cen<pix(i,j-1) &&
			cen<pix(i-1,j+1) && cen<pix(i-1,j) && cen<pix(i-1,j-1)
	    );
	}
	bool esMin2D_anillo(double val, int i, int j, double r);
	inline bool esMax2D(double &cen, int i, int j)
	{
		return (cen>pix(i+1,j+1) && cen>pix(i+1,j) && cen>pix(i+1,j-1) &&
		    cen>pix(i,j+1) && cen>pix(i,j) && cen>pix(i,j-1) &&
			cen>pix(i-1,j+1) && cen>pix(i-1,j) && cen>pix(i-1,j-1)
	    );
	}
	inline bool esMax2D(double &cen, int i, int j, double porc) // idea: porc < 1
	{
		double cen2=cen*(double)porc;
		return (cen2>pix(i+1,j+1) && cen2>pix(i+1,j) && cen2>pix(i+1,j-1) &&
		    cen2>pix(i,j+1) && cen2>pix(i,j) && cen2>pix(i,j-1) &&
			cen2>pix(i-1,j+1) && cen2>pix(i-1,j) && cen2>pix(i-1,j-1)
	    );
	}
	inline bool esMin2D(double &cen, int i, int j)
	{
		return (cen<pix(i+1,j+1) && cen<pix(i+1,j) && cen<pix(i+1,j-1) &&
		    cen<pix(i,j+1) && cen<pix(i,j) && cen<pix(i,j-1) &&
			cen<pix(i-1,j+1) && cen<pix(i-1,j) && cen<pix(i-1,j-1)
	    );
	}
	inline bool esMin2D(double &cen, int i, int j, double porc) // idea: porc < 1
	{
		double cen2=cen*(double)porc;
		return (cen2<pix(i+1,j+1) && cen2<pix(i+1,j) && cen2<pix(i+1,j-1) &&
		    cen2<pix(i,j+1) && cen2<pix(i,j) && cen2<pix(i,j-1) &&
			cen2<pix(i-1,j+1) && cen2<pix(i-1,j) && cen2<pix(i-1,j-1)
	    );
	}
	L_ImageGrayDouble& harris(const L_ImageGrayDouble &other, double sd=1.3, double si=2.0, double alfa=0.04); // Usar normalizaHist(30) para mostrarlo
	void harrisEliminaLineas(double cornernessMin, const L_ImageGrayDouble &other, double *x, double *y, double n, bool *eliminar, double sd=1.3, double si=2.0, double alfa=0.04);
	static double harrisIndividual(L_Array<L_ImageGrayDouble> &imArrTmp, const L_ImageGrayDouble &other, double x, double y, double sd=1.3, double si=2.0, double alfa=0.04);

	static void z_difumXGenCodigo(double s, int s3=-1, FILE *fp=stdout);
	static void z_difumYGenCodigo(double s, int s3=-1, FILE *fp=stdout);
	static inline void z_difumXGenCodigo(double s, FILE *fp) {z_difumXGenCodigo(s,-1,fp);}
	static inline void z_difumYGenCodigo(double s, FILE *fp) {z_difumYGenCodigo(s,-1,fp);}

	// Funciones de convolucion gaussian precalculada
	#if defined(L_defFnBNxfloat)
		#error #define conflicts
	#endif
	#define L_defFnBNxfloat(S) L_ImageGrayDouble& z_difumX_##S(const L_ImageGrayDouble &im1, int paso);L_ImageGrayDouble& z_difumY_##S(const L_ImageGrayDouble &im1, int paso)
	//
	#ifdef _s0_0_5__nio_2
	// Para scale space, s0=0.5, 2 imagenes por octava
	L_defFnBNxfloat(8660);  // solo para sDoG
	L_defFnBNxfloat(10000);
	L_defFnBNxfloat(14142);
	L_defFnBNxfloat(20000);
	L_defFnBNxfloat(28284);
	L_defFnBNxfloat(30000); // solo para Harris-Laplace
	L_defFnBNxfloat(42426); // solo para Harris-Laplace
	#endif
	#ifdef _s0_1_6__nio_2
	L_defFnBNxfloat(12489);
	L_defFnBNxfloat(16000);
	L_defFnBNxfloat(22627);
	L_defFnBNxfloat(32000);
	L_defFnBNxfloat(45254);
	L_defFnBNxfloat(48000); // solo para Harris-Laplace
	L_defFnBNxfloat(67882); // solo para Harris-Laplace
	#endif
	//
	#undef L_defFnBNxfloat

	L_ImageGrayDouble & operator=(const L_ImageGrayDouble &other) {OP_assign(other); return *this;}
	L_ImageGrayDouble & operator=(const L_ImageGrayUchar &other);
	L_ImageGrayDouble & operator=(const L_ImageRGBDouble &other);
	L_ImageGrayDouble & operator=(const L_ImageRGBUchar &other);
	#ifdef __COMPAT_ATLIMAGE__
	void operator=(const ATL_CImage &other);
	void copyTo(ATL_CImage &other);
	#endif

	#ifdef __COMPAT_IPLIMAGE__
	void operator=(const IplImage &other);
	void copyTo(IplImage **other); // Los IplImage deben: o ser NULL, o tener ya un object
	#endif

#ifdef __COMPAT_CVMAT__
	void operator=(const cv::Mat& other);
	void copyTo(cv::Mat& other);
#endif


#ifdef __COMPAT_EIGEN__
	Eigen::MatrixXd crearEigen() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenLenta() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenRef() const;  // createFrom una referencia a una L_Matrix
	void operator=(const Eigen::MatrixXd &other);
#endif

	static double channelY(double R, double G, double B) {return +(double)0.299*R +(double)0.587*G +(double)0.114*B;}
	static double channelU(double R, double G, double B) {return -(double)0.168*R -(double)0.330*G +(double)0.498*B;}
	static double channelV(double R, double G, double B) {return +(double)0.498*R -(double)0.417*G -(double)0.081*B;}

	static double channelR(double Y, double U, double V) {return (double)1.0*Y +(double)0.000*U +(double)1.4075*V;}
	static double channelG(double Y, double U, double V) {return (double)1.0*Y -(double)0.3455*U -(double)0.7169*V;}
	static double channelB(double Y, double U, double V) {return (double)1.0*Y +(double)1.779*U +(double)0.000*V;}

#ifdef L_USE_SMOOTH_SIGMA_STORING_TREE
	// overhead
	static L_SigmasConvGausCounter cuentaSigmas; // Para contar cuales sigmas se usan en las convoluciones gaussianas
#endif
};

class L_ImageGrayUchar : public L_MatrixBaseWidthStep<L_uchar>
{
private:
	L_ImageGrayUchar(const L_ImageGrayUchar &other);
public:
	L_ImageGrayUchar() : L_MatrixBaseWidthStep<L_uchar>() {}
	L_ImageGrayUchar(int lx, int ly) : L_MatrixBaseWidthStep<L_uchar>(ly, lx) {}
	L_ImageGrayUchar(int lx, int ly, const L_uchar *buf) : L_MatrixBaseWidthStep<L_uchar>(ly, lx, buf) {}

	inline L_uchar &pix(int c, int f) {return operator()(f,c);}
	inline const L_uchar &pix(int c, int f) const {return operator()(f,c);}

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}

	void fillRandomValues(int v1=0, int v2=256) {int c, f; for (f=0; f<ly; f++) for (c=0; c<lx; c++) pix(c,f) = (rand()+v1)%(v2-v1);}
	void multiply_to_each_element(double val) {int c, f; for (f=0; f<li; f++) for (c=0; c<lj; c++) {double res = operator()(f,c)*val; if (res < 0) res = 0; if (res > 255) res = 255; operator()(f,c) = (L_uchar)res;}}

	bool difusion_pinta(int i, int j, L_uchar lum);
	L_ImageGrayUchar& genDrawing(const L_ShapeArray& lins) {return  genDibujo_CenGrav(lins);} // Usa el canal R para elegir el tono de gris
	L_ImageGrayUchar& genDibujo_CenGrav(const L_ShapeArray& lins, double *cx=NULL, double *cy=NULL); // Usa el canal R para elegir el tono de gris

	void supresionNoMaximo(const L_ImageGrayDouble &gradR, const L_ImageGrayDouble &gradA);
	void morfol_eliminaPolvo(const L_ImageGrayUchar &orig);
	int morfol_hitAndMiss(const L_ImageGrayUchar &orig, const char *mascara3x3); // mask = 0, 1, 2 (omitir)
	int morfol_hitAndMiss_resta(const L_ImageGrayUchar &orig, const char *mascara3x3); // mask = 0, 1, 2 (omitir)
	int morfol_erosiona(const L_ImageGrayUchar &orig, const char *mascara3x3); // mask = 0, 1
	int morfol_dilata(const L_ImageGrayUchar &orig, const char *mascara3x3); // mask = 0, 1
	int morfol_adelgaza(const L_ImageGrayUchar &orig);
	void OP_max(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2);
	void OP_min(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2);
	void OP_max(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2, L_ImageGrayUchar &im3);
	void OP_min(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2, L_ImageGrayUchar &im3);
	void lbp(L_ImageGrayUchar &orig);
	void lbp(L_ImageGrayDouble &orig);


	void encodeToBuffer(std::vector<char> &buf) const;
	void decodeFromBuffer(const std::vector<char> &buf);
	void encodeDelta(std::vector<char> &bufD) const;
	void decodeDelta(std::vector<char> &bufD);

	void intercambiaObjetoCon(L_ImageGrayUchar &other) {swap(other);}

	int readBMP(const char *name);
    bool writeBMP(const char *name, short nb);
    bool writeBMP(const char *name);

	bool readImage(const char *name);
	bool saveImage(const char *name);

#ifdef __COMPAT_ATLIMAGE__
	int readImageATL(const char *name);
	bool writeImageATL(const char *name);
#endif

	bool leeB8(const char *name);
	bool grabaB8(const char *name);

	bool grabaPGM(const char *name);

	int grabaComprGris(FILE *fp, double *le=NULL, int forzMetodo=-1, double t=0, long ncam=0) const;
	int leeComprGris(FILE *fp, double *t=NULL, long *ncam=NULL);
	int leeComprGris_contar(FILE *fp, int *histMetodo = NULL);
	static void pruebaComprGris();

	bool normalizeHistogram();

	L_ImageGrayUchar& subm50(const L_ImageGrayUchar &other);

	void operator+=(const L_ImageGrayUchar &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator+= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) += other(i,j);}
	void operator-=(const L_ImageGrayUchar &other) {throw_L_ArgException_if(this->li != other.li || this->lj != other.lj, "operator-= : tamanos de imagenes diferentes"); for (int i=0; i < this->li; i++) for (int j=0; j < this->lj; j++) operator()(i,j) -= other(i,j);}
	void OP_add(const L_ImageGrayUchar &im1, const L_ImageGrayUchar &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayUchar::OP_add"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i<this->li; i++) for (j=0; j<this->lj; j++) operator()(i,j) = im1(i,j) + im2(i,j);}
	void OP_subtract(const L_ImageGrayUchar &im1, const L_ImageGrayUchar &im2) { throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayUchar::OP_subtract"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) - im2(i,j);}
	void OP_mult_elementwise(const L_ImageGrayUchar &im1, const L_ImageGrayUchar &im2) {throw_L_ArgException_if(im1.li!=im2.li || im1.lj!=im2.lj, "L_ImageGrayUchar::OP_mult_elementwise"); int i, j; reallocate_fc(im1.li, im2.lj); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(i,j) = im1(i,j) * im2(i,j);}
	void transpOf(const L_ImageGrayUchar &other) {int i, j; reallocate_fc(other.lj, other.li); for (i=0; i < this->li; i++) for (j=0; j < this->lj; j++) operator()(j,i) = other(i,j);}
	void transposeMe() {L_ImageGrayUchar other; other.transpOf(*this); other.swap(*this);}

	L_ImageGrayUchar & operator=(const L_ImageGrayDouble &other);
	L_ImageGrayUchar & operator=(const L_ImageGrayUchar &other) {OP_assign(other); return *this;}
	L_ImageGrayUchar & operator=(const L_ImageRGBDouble &other);
	L_ImageGrayUchar & operator=(const L_ImageRGBUchar &other);
	#ifdef __COMPAT_ATLIMAGE__
	void operator=(const ATL_CImage &other);
	void copyTo(ATL_CImage &other);
	#endif

	#ifdef __COMPAT_IPLIMAGE__
	void operator=(const IplImage &imagen1);
	void copyTo(IplImage **imagen2);
	#endif

#ifdef __COMPAT_CVMAT__
	void operator=(const cv::Mat& imagen1);
	void copyTo(cv::Mat& imagen2);
#endif
};

class L_ImageGrayUint16 : public L_MatrixBaseWidthStep<L_uint16>
{
private:
	L_ImageGrayUint16(const L_ImageGrayUint16 &other);
public:
	L_ImageGrayUint16() : L_MatrixBaseWidthStep<L_uint16>() {}
	L_ImageGrayUint16(int lx, int ly) : L_MatrixBaseWidthStep<L_uint16>(ly, lx) {}
	L_ImageGrayUint16(int lx, int ly, const L_uint16 *buf) : L_MatrixBaseWidthStep<L_uint16>(ly, lx, buf) {}

	inline L_uint16 &pix(int c, int f) {return operator()(f,c);}
	inline const L_uint16 &pix(int c, int f) const {return operator()(f,c);}

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}
	void intercambiaObjetoCon(L_ImageGrayUint16 &other) {swap(other);}

	void operator=(const L_ImageGrayUint16 &other) {other.copyTo(*this);}

	void doubleResolution(const L_ImageGrayUint16 &other);
    void doubleResolutionBilineal(const L_ImageGrayUint16 &other);
	void halfResolution(const L_ImageGrayUint16 &other);
	void halfResolution_noResample(const L_ImageGrayUint16 &other);

	void doubleResolution() {L_ImageGrayUint16 im; im.doubleResolution(*this); im.swap(*this);}
	void halfResolution() {L_ImageGrayUint16 im; im.halfResolution(*this); im.swap(*this);}
	void halfResolution_noResample() {L_ImageGrayUint16 im; im.halfResolution_noResample(*this); im.swap(*this);}

	bool saveImage_16_txt(const char *name);
	bool saveImage_16_pgm(const char *name);
	bool saveImage_16(const char *name);
	bool saveImage_8(const char *name, double factor) {L_ImageGrayUchar im; im.reallocate(lx, ly); for (int i=0; i<li; i++) for (int j=0; j<lj; j++) im(i,j) = (L_uchar)(operator()(i,j)*factor); return im.saveImage(name);}
	bool savePoints_xyzrgb(const char *name, const L_ImageRGBUchar &im, double dfx, double dfy, double factor);
	bool savePoints_xyz(const char *name, double dfx, double dfy, double factor);
	bool readImage_16_txt(const char *name);


#ifdef __COMPAT_IPLIMAGE__
	void copyOn_16(IplImage **imagen2);
	bool saveImageIPL_16(const char *name) {IplImage *imIPL; copyOn_16(&imIPL); cvSaveImage(name, imIPL); cvReleaseImage(&imIPL); return true;}
#endif // __COMPAT_IPLIMAGE__

	void fillRandomValues(int v1=0, int v2=256) {int c, f; for (f=0; f<ly; f++) for (c=0; c<lx; c++) pix(c,f) = (rand()+v1)%(v2-v1);}
};

class L_ImageRGBDouble : public L_MatrixBaseWidthStep<double[3]>
{
private:
	L_ImageRGBDouble(const L_ImageRGBDouble &other);
public:
	L_ImageRGBDouble() : L_MatrixBaseWidthStep<double[3]>() {}
	L_ImageRGBDouble(int lx, int ly) : L_MatrixBaseWidthStep<double[3]>(ly, lx) {}
	L_ImageRGBDouble(int lx, int ly, const double *buf) : L_MatrixBaseWidthStep<double[3]>(ly, lx, (const double(*)[3])buf) {}

	inline double &pix(int c, int f, int rgb) {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}
	inline const double &pix(int c, int f, int rgb) const {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}

	inline double &operator() (int f, int c, int rgb) {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}
	inline const double & operator() (int f, int c, int rgb) const {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}

	void fill(double num) {int i, j; for (i=0; i<ly; i++) for (j=0; j<lx; j++) {operator()(i,j,0) = operator()(i,j,1) = operator()(i,j,2) = num;}}
	void setZero() {memset(iniptr, 0, lxStep*ly);}

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}

	void copyFromBuffer(const double *buf) {L_MatrixBaseWidthStep<double[3]>::copyFromBuffer((const double(*)[3])buf);}
	void copyToBuffer(double *buf) const {L_MatrixBaseWidthStep<double[3]>::copyToBuffer((double(*)[3])buf);}

	int readBMP(const char *name); // devuelve nº bits por pixel, 0 si hay error
    bool writeBMP(const char *name, short nb=24);

	bool readImage(const char *name);
	bool saveImage(const char *name);

	L_ImageRGBDouble &conv2d(const L_ImageRGBDouble &other, const L_ImageRGBDouble &mask, int xCenMasc, int yCenMasc); // La mask se le pasa ya dada vuelta en 180°

#ifdef __COMPAT_ATLIMAGE__
	int readImageATL(const char *name);
	bool writeImageATL(const char *name);
#endif

	bool normalizeHistogram();

    L_ImageRGBDouble& transpOf(const L_ImageRGBDouble &im);
	L_ImageRGBDouble& transposeMe() {L_ImageRGBDouble other; other.transpOf(*this); other.swap(*this); return *this;}

	L_ImageGrayDouble& copyR_into(L_ImageGrayDouble &imR);
	L_ImageGrayDouble& copyG_into(L_ImageGrayDouble &imG);
	L_ImageGrayDouble& copyB_into(L_ImageGrayDouble &imB);

    L_ImageRGBDouble& composeFromChannels(const L_ImageGrayDouble& imR, const L_ImageGrayDouble& imG, const L_ImageGrayDouble& imB);

	L_ImageRGBDouble &RGB_to_YPbPr();
	L_ImageRGBDouble &YPbPr_to_RGB();

	L_ImageRGBDouble & operator=(const L_ImageGrayDouble &other);
	L_ImageRGBDouble & operator=(const L_ImageGrayUchar &other);
	L_ImageRGBDouble & operator=(const L_ImageRGBDouble &other) {OP_assign(other); return *this;}
	L_ImageRGBDouble & operator=(const L_ImageRGBUchar &other);
	#ifdef __COMPAT_ATLIMAGE__
	void operator=(const ATL_CImage &other);
	void copyTo(ATL_CImage &other);
	#endif

	#ifdef __COMPAT_IPLIMAGE__
	void operator=(const IplImage &imagen1);
	void copyTo(IplImage **imagen2);
	#endif

#ifdef __COMPAT_EIGEN__
	Eigen::MatrixXd crearEigenLentaR() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenLentaG() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenLentaB() const;  // createFrom una referencia a una L_Matrix
	void operator=(const Eigen::MatrixXd &other);
#endif
};




class L_ImageRGBUchar : public L_MatrixBaseWidthStep<L_uchar[3]>
{
private:
	L_ImageRGBUchar(const L_ImageRGBUchar &other);
public:
	L_ImageRGBUchar() : L_MatrixBaseWidthStep<L_uchar[3]>() {}
	L_ImageRGBUchar(int lx, int ly) : L_MatrixBaseWidthStep<L_uchar[3]>(ly, lx) {}
	L_ImageRGBUchar(int lx, int ly, const L_uchar *buf) : L_MatrixBaseWidthStep<L_uchar[3]>(ly, lx, (const L_uchar(*)[3])buf) {}

	void reallocate_cf_ext(int lx, int ly, char *buf, int lx_step) {L_MatrixBaseWidthStep<L_uchar[3]>::reallocate_cf_ext(lx, ly, (L_uchar(*)[3])buf, lx_step);}
	void reallocate_fc_ext(int li, int lj, char *buf, int lx_step) {L_MatrixBaseWidthStep<L_uchar[3]>::reallocate_fc_ext(li, lj, (L_uchar(*)[3])buf, lx_step);}

	inline L_uchar &pix(int c, int f, int rgb) {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}
	inline const L_uchar &pix(int c, int f, int rgb) const {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}

	inline L_uchar &operator()(int f, int c, int rgb) {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}
	inline const L_uchar &operator()(int f, int c, int rgb) const {return L_MatrixBaseWidthStep::operator()(f,c)[rgb];}

	void setZero() {memset(iniptr,0,lxStep*ly);}
	void fill(L_uchar num=0) {memset(iniptr,num,lxStep*ly);}
	void fill(L_uchar R, L_uchar G, L_uchar B);
	unsigned char *pedirPunteroInterno() {return (unsigned char *)&pix(0,0,0);}
	void fillRandomValues(int v1=0, int v2=256) {int i, j; for (j=0; j<ly; j++) for (i=0; i<lx; i++) {pix(i,j,0) = (rand()+v1)%(v2-v1); pix(i,j,1) = (rand()+v1)%(v2-v1); pix(i,j,2) = (rand()+v1)%(v2-v1);}}

	void reallocate(int lx, int ly) {reallocate_cf(lx,ly);}

	void copyFromBuffer(const L_uchar *buf) {L_MatrixBaseWidthStep<L_uchar[3]>::copyFromBuffer((const L_uchar(*)[3])buf);}
	void copyToBuffer(L_uchar *buf) const {L_MatrixBaseWidthStep<L_uchar[3]>::copyToBuffer((L_uchar(*)[3])buf);}

    int readBMP(const char *name); // devuelve nº bits por pixel, 0 si hay error
    bool writeBMP(const char *name, short nb=24);

	bool readImage(const char *name);
	bool saveImage(const char *name);

#ifdef __COMPAT_ATLIMAGE__
	bool readImageATL(const char *name);
	bool writeImageATL(const char *name);
#endif

#ifdef __COMPAT_IPLIMAGE__
	bool leeImagenIPL(const char *name);
	bool grabaImagenIPL(const char *name);
#endif

#ifdef __COMPAT_CVMAT__
	void copyTo(cv::Mat& other);
#endif

#ifdef __COMPAT_TBITMAP__
	bool leeImagenTBitmap(const char *name);
	bool grabaImagenTBitmap(const char *name);
#endif

	bool leeB24(const char *name); // read raw image: {int16 lx, int16 ly, char R, char G, char B, ...}
	bool grabaB24(const char *name);

	int writeComprRGB(FILE *fp, double *le = NULL, int forzMetodo = -1, double t=0, long ncam=0) const;
	int readComprRGB(FILE *fp, double *t=NULL, long *ncam=NULL);
	static int leeComprRGB_contar(FILE *fp, int *histMetodo = NULL);
	static void testComprRGB();
	static void pruebaVelocComprRGB();

    bool normalizeHistogram();
	L_ImageRGBUchar& doubleResolution(const L_ImageRGBUchar &other);
	L_ImageRGBUchar& doubleResolution() {L_ImageRGBUchar im; im.doubleResolution(*this); im.swap(*this); return *this;}

	L_ImageRGBUchar& halfResolution(const L_ImageRGBUchar &other);
	L_ImageRGBUchar& halfResolution() {L_ImageRGBUchar im; im.halfResolution(*this); im.swap(*this); return *this;}

	L_ImageRGBUchar& halfResolution_noResample(const L_ImageRGBUchar &other);
	L_ImageRGBUchar& halfResolution_noResample() {L_ImageRGBUchar im; im.halfResolution_noResample(*this); im.swap(*this); return *this;}

	void applyHomography(const L_Matrix &h8, const L_ImageRGBUchar &orig);
	void changeBrightness(double factor) {unsigned char *p = pedirPunteroInterno(), *pf = p+lxStep*ly; for ( ; p<pf; p++) *p = (unsigned char)(*p*factor);}

	L_ImageRGBUchar& transpOf(const L_ImageRGBUchar &im);
	L_ImageRGBUchar& transposeMe() {L_ImageRGBUchar other; other.transpOf(*this); other.swap(*this); return *this;}

	L_ImageGrayUchar& copyR_into(L_ImageGrayUchar &imR);
	L_ImageGrayUchar& copyG_into(L_ImageGrayUchar &imG);
	L_ImageGrayUchar& copyB_into(L_ImageGrayUchar &imB);

	void concatenateHorz(const L_ImageRGBUchar &im1, const L_ImageRGBUchar &im2);
	void divideHorz(L_ImageRGBUchar &im1, L_ImageRGBUchar &im2) const;
	L_ImageRGBUchar& composeFromChannels(const L_ImageGrayUchar& imR, const L_ImageGrayUchar& imG, const L_ImageGrayUchar& imB);
	L_ImageRGBUchar& transformaDaltonicoLonco(bool transfInversa);
	L_ImageRGBUchar& transformaDaltonicoLonco2(bool transfInversa);

	void encodeToBuffer(std::vector<char> &buf) const;
	void decodeFromBuffer(const std::vector<char> &buf);
	void encodeDelta(std::vector<char> &bufD) const;
	void decodeDelta(std::vector<char> &bufD);

	L_ImageRGBUchar& copyImageWithOffset(const L_ImageGrayDouble &im, int xc, int yc);
	L_ImageRGBUchar& copyImageWithOffset(const L_ImageRGBUchar &im, int xc, int yc);
	L_ImageRGBUchar& genDrawing(const L_ShapeArray& lins);
	L_ImageRGBUchar& genDrawingMatches(L_Array<const L_ImageGrayDouble *> &imarr, L_ShapeArray& lins);
	void convertToGray() {size_type c, f; for (f=0; f<ly; f++) for (c=0; c<lx; c++) pix(c,f,2)=pix(c,f,1)=pix(c,f,0)=(L_uchar)(((int)pix(c,f,0)+pix(c,f,1)+pix(c,f,2))/3);}

	void subImageOf(L_ImageRGBUchar &other, int x0, int y0, int nx, int ny) { reallocate(nx, ny); int x, y; for (y=0; y<ny; y++) for (x=0; x<nx; x++) {pix(x,y,0) = other.pix(x+x0,y+y0,0); pix(x,y,1) = other.pix(x+x0,y+y0,1); pix(x,y,2) = other.pix(x+x0,y+y0,2);}}

	bool searchInside(L_ImageRGBUchar &imGrande, int &x, int &y);
	bool searchInside_delta(L_ImageRGBUchar &imGrande, int delta, int &x, int &y);
	bool searchInside_transp(L_ImageRGBUchar &imGrande, L_uchar R, L_uchar G, L_uchar B, int x0, int y0, int x1, int y1, int &x, int &y);
	bool searchInside_transp_delta(L_ImageRGBUchar &imGrande, int delta, L_uchar R, L_uchar G, L_uchar B, int x0, int y0, int x1, int y1, int &x, int &y);

	void computeMeanShiftHistogram(double xCen, double yCen, std::vector<double> &hist, double radio_h_x, double radio_h_y, int nBinsPorCanal = 32);
	void computeMeanShiftTracking_private(double &xCen, double &yCen, std::vector<double> &histAntiguo, std::vector<double> &histNuevo, double radio_h_x, double radio_h_y, int nBinsPorCanal = 32);
	double computeMeanShiftTracking_private(double &xCen, double &yCen, std::vector<double> &hist, double radio_h_x, double radio_h_y, int nBinsPorCanal = 32);
	double computeMeanShiftTracking(double &xCen, double &yCen, std::vector<double> &hist, double &radio_h_x, double &radio_h_y, bool adaptTamano = false, int nBinsPorCanal = 32);

	L_ImageRGBUchar & operator=(const L_ImageGrayDouble &other);
	L_ImageRGBUchar & operator=(const L_ImageGrayUchar &other);
	L_ImageRGBUchar & operator=(const L_ImageRGBDouble &other);
	L_ImageRGBUchar & operator=(const L_ImageRGBUchar &other) {OP_assign(other); return *this;}
	void operator=(const L_ImageBGRAuchar_row &other);
	void operator=(const L_ImageRGBAuchar_row &other);
	void operator=(const L_ImagenRGBuchar_row &other);

	#ifdef __COMPAT_ATLIMAGE__
	void operator=(const ATL_CImage &other);
	void copyTo(ATL_CImage &other) const;
	bool copyFromsde_lento(const ATL_CImage &other); // Espero que este buena, es complicada
	#endif

	#ifdef __COMPAT_IPLIMAGE__
	void operator=(const IplImage &imagen1);
	void copyTo(IplImage **imagen2) const;
	#endif

	#ifdef __COMPAT_SVS__
	bool SVSPideIzq(const svsStereoImage *imageObject) {if (!SVSPideIzqColor(imageObject)) return SVSPideIzqBN(imageObject); return true;}
	bool SVSPideDer(const svsStereoImage *imageObject) {if (!SVSPideDerColor(imageObject)) return SVSPideDerBN(imageObject); return true;}
	bool SVSPideIzqBN(const svsStereoImage *imageObject);
	bool SVSPideDerBN(const svsStereoImage *imageObject);
	bool SVSPideIzqColor(const svsStereoImage *imageObject);
	bool SVSPideDerColor(const svsStereoImage *imageObject);
	#endif

	static void computePositionsForComposingImages(int n, int *imLx, int *imLy, int *imPosX, int *imPosY, int &tamTotalX, int &tamTotalY);
};

class L_ImageBGRAuchar_row
{
public:
	L_Array<L_uchar> buf;
	int lx, ly, length;
	void reallocate(int lx, int ly) {if ((L_uchar *)buf.begin()!=NULL && this->lx==lx && this->ly==ly) return; this->lx=lx; this->ly=ly; this->length=lx*ly*4; buf.resize(length);}
	L_ImageBGRAuchar_row() { }
	L_ImageBGRAuchar_row(const L_ImageBGRAuchar_row &other) {reallocate(other.lx, other.ly); memcpy((L_uchar *)buf.begin(), (const L_uchar *)other.buf.begin(), sizeof(L_uchar)*length);}
	void operator=(const L_ImageRGBUchar &im) {int i, j, k=0; reallocate(im.lx, im.ly); for (j=0; j<ly; j++) for (i=0; i<lx; i++) {buf[k++]=im.pix(i,j,2); buf[k++]=im.pix(i,j,1); buf[k++]=im.pix(i,j,0); buf[k++]=0;}}
	void operator=(L_ImageBGRAuchar_row &im) {reallocate(im.lx, im.ly); memcpy(buf.begin(), im.buf.begin(), sizeof(L_uchar)*length);}
	void copyFrom(const L_uchar *im, int lx, int ly) {reallocate(lx, ly); memcpy(buf.begin(), im, sizeof(L_uchar)*length);}
	~L_ImageBGRAuchar_row() { }
};

class L_ImageRGBAuchar_row
{
public:
	L_Array<L_uchar> buf;
	int lx, ly, length;
	void reallocate(int lx, int ly) {if ((const L_uchar*)buf.begin()!=NULL && this->lx==lx && this->ly==ly) return; this->lx=lx; this->ly=ly; this->length=lx*ly*4; buf.resize(length);}
	L_ImageRGBAuchar_row() { }
	L_ImageRGBAuchar_row(const L_ImageRGBAuchar_row &other) {reallocate(other.lx, other.ly); memcpy(buf.begin(), other.buf.begin(), sizeof(L_uchar)*length);}
	void operator=(const L_ImageRGBUchar &im) {reallocate(im.lx, im.ly); const L_uchar *im1 = &im.pix(0,0,0); L_uchar *im2 = &buf[0], *im2f = im2 + buf.size()-4; im2[0] = im1[0]; im2[1] = im1[1]; im2[2] = im1[2]; im2[3] = 0; im2 += 4-1; im1 += 3-1; while (im2 < im2f) {*(++im2) = *(++im1);*(++im2) = *(++im1);*(++im2) = *(++im1);*(++im2) = 0;}}
	void operator=(const L_ImageRGBAuchar_row &im) {reallocate(im.lx, im.ly); memcpy((L_uchar *)&buf[0], (L_uchar *)&im.buf[0], sizeof(L_uchar)*length);}
	void copyFrom(const L_uchar *im, int lx, int ly) {reallocate(lx, ly); memcpy(buf.begin(), im, sizeof(L_uchar)*length);}
	~L_ImageRGBAuchar_row() { }
};

class L_ImagenRGBuchar_row
{
public:
	L_Array<L_uchar> buf;
	int lx, ly, length;
	void reallocate(int lx, int ly) {if ((L_uchar *)buf.begin()!=NULL && this->lx==lx && this->ly==ly) return; this->lx=lx; this->ly=ly; this->length=lx*ly*3; buf.resize(length);}
	L_ImagenRGBuchar_row() {}
	L_ImagenRGBuchar_row(const L_ImagenRGBuchar_row &other) {reallocate(other.lx, other.ly); memcpy((L_uchar *)buf.begin(), (const L_uchar *)other.buf.begin(), sizeof(L_uchar)*length);}
	void operator=(const L_ImageRGBUchar &im) {int i, j, k=0; reallocate(im.lx, im.ly); for (j=0; j<ly; j++) for (i=0; i<lx; i++) {buf[k++]=im.pix(i,j,0); buf[k++]=im.pix(i,j,1); buf[k++]=im.pix(i,j,2);}}
	void operator=(L_ImagenRGBuchar_row &im) {reallocate(im.lx, im.ly); memcpy((L_uchar *)buf.begin(), (L_uchar *)im.buf.begin(), sizeof(L_uchar)*length);}
	void copyFrom(const L_uchar *im, int lx, int ly) {reallocate(lx, ly); memcpy(buf.begin(), im, sizeof(L_uchar)*length);}
	~L_ImagenRGBuchar_row() { }
};


#if defined(L_M_NUMERITO) || defined(L_M_NAN) || defined(L_M_NEG) || defined(L_M_ASI) || defined(L_M_DIA) || defined(L_M_CER) || defined(L_M_COV) || defined(L_M_RAN)
#error Nombres para define ocupados
#endif

#define L_M_NUMERITO (0.37586383729e5)



// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes

// Macros para cualquier tamano de matriz
#define L_StaticMatrix_checkIfDefined(M,LI,LJ) ( throw_L_ArgException_if((M).rows()=!(LI) || (M).cols()!=(LJ), "L_StaticMatrix_revDef()")) // build error si la matriz M no esta definida de tamano (li,lj)
#define L_StaticMatrix_clean(M) {int df_i925, df_j925; for (df_i925=0; df_i925<(M).rows(); df_i925++) for (df_j925=0; df_j925<(M).cols(); df_j925++) (M)(df_i925,df_j925) = 0.0;}
#define L_StaticMatrix_OP_assign(M,A) {int df_i94, df_j94; for (df_i94=0; df_i94<(M).rows(); df_i94++) for (df_j94=0; df_j94<(M).cols(); df_j94++) (M)(df_i94,df_j94) = (A)(df_i94,df_j94);}
#define L_StaticMatrix_OP_add(M,A,B) {int df_i923, df_j923; for (df_i923=0; df_i923<(M).rows(); df_i923++) for (df_j923=0; df_j923<(M).cols(); df_j923++) (M)(df_i923,df_j923) = (A)(df_i923,df_j923) + (B)(df_i923,df_j923);}
#define L_StaticMatrix_OP_subtract(M,A,B) {int df_i912, df_j912; for (df_i912=0; df_i912<(M).rows(); df_i912++) for (df_j912=0; df_j912<(M).cols(); df_j912++) (M)(df_i912,df_j912) = (A)(df_i912,df_j912) - (B)(df_i912,df_j912);}
#define L_StaticMatrix_OP_mult(M,A,B) {int df_i96, df_j96, df_u96; for (df_i96=0; df_i96<(A).rows(); df_i96++) for (df_j96=0; df_j96<(B).cols(); df_j96++) {(M)(df_i96,df_j96)=0; for (df_u96=0; df_u96<(A).cols(); df_u96++) (M)(df_i96,df_j96) += (A)(df_i96,df_u96) * (B)(df_u96,df_j96);}} // ESTABA MALO!!!
#define L_StaticMatrix_OP_mult_sub(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2) {int df_i98,df_j98,df_u98; throw_L_ArgException_if((lj1)!=(li2) || (M).begin()==NULL || (M).rows()<(i01)+(li1) || (M).cols()<(j02)+(lj2), "L_StaticMatrix_OP_mult_sub"); for (df_i98=0; df_i98<(li1); df_i98++) for (df_j98=0; df_j98<(lj2); df_j98++) { (M)(df_i98+(i0),df_j98+(j0))=0; for (df_u98=0; df_u98<(li2); df_u98++) (M)(df_i98+(i0),df_j98+(j0))+=(m1)(df_i98+(i01),df_u98+(j01))*(m2)(df_u98+i02,df_j98+j02);}}
#define L_StaticMatrix_OP_mult_sub_sum(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2) {int df_i99,df_j99,df_u99; throw_L_ArgException_if((lj1)!=(li2) || (M).begin()==NULL || (M).rows()<(i01)+(li1) || (M).cols()<(j02)+(lj2), "L_StaticMatrix_OP_mult_sub_sum"); for (df_i99=0; df_i99<(li1); df_i99++) for (df_j99=0; df_j99<(lj2); df_j99++) { for (df_u99=0; df_u99<(li2); df_u99++) (M)(df_i99+(i0),df_j99+(j0))+=(m1)(df_i99+(i01),df_u99+(j01))*(m2)(df_u99+i02,df_j99+j02);}}
#define L_StaticMatrix_OP_mult_sub1_sum(M,m1,m2,i01,j01,li1,lj1) {int m74i,m74j,m74u; (M).reallocate((m1).rows(), (m2).cols()); for (m74i=(i01); m74i<(i01)+(li1); m74i++) { for (m74j=0; m74j<(m2).cols(); m74j++) { for (m74u=(j01); m74u<(j01)+(lj1); m74u++) (M)(m74i,m74j)+=(m1)(m74i,m74u)*(m2)(m74u,m74j); } } }
#define L_StaticMatrix_OP_mult_sub2_sum(M,m1,m2,i02,j02,li2,lj2) {int m93i,m93j,m93u; (M).reallocate((m1).rows(), (m2).cols()); for (m93i=0; m93i<(m1).rows(); m93i++) { for (m93j=(j02); m93j<(j02)+(lj2); m93j++) { for (m93u=(i02); m93u<(i02)+(li2); m93u++) (M)(m93i,m93j)+=(m1)(m93i,m93u)*(m2)(m93u,m93j); } } }

#define L_StaticMatrix_OP_amplifica(M,val) {int df_i934, df_j934; for (df_i934=0; df_i934<(M).li; df_i934++) for (df_j934=0; df_j934<(M).lj; df_j934++) (M)(df_i934,df_j934)*=(val);}

// Macros especializadas para 2x2, 3x3, 4x4
#define L_MEme2(M,A,B,i,j) ((M)(i,j)=(A)(i,0)*(B)(0,j) + (A)(i,1)*(B)(1,j))
#define L_MEme3(M,A,B,i,j) ((M)(i,j)=(A)(i,0)*(B)(0,j) + (A)(i,1)*(B)(1,j) + (A)(i,2)*(B)(2,j))
#define L_MEme4(M,A,B,i,j) ((M)(i,j)=(A)(i,0)*(B)(0,j) + (A)(i,1)*(B)(1,j) + (A)(i,2)*(B)(2,j) + (A)(i,3)*(B)(3,j))
#define L_StaticMatrix_OP_mult2x2(M,A,B) (L_MEme2(M,A,B,0,0),L_MEme2(M,A,B,0,1),L_MEme2(M,A,B,1,0),L_MEme2(M,A,B,1,1))
#define L_StaticMatrix_OP_mult3x3(M,A,B) (L_MEme3(M,A,B,0,0),L_MEme3(M,A,B,0,1),L_MEme3(M,A,B,0,2),L_MEme3(M,A,B,1,0),L_MEme3(M,A,B,1,1),L_MEme3(M,A,B,1,2),L_MEme3(M,A,B,2,0),L_MEme3(M,A,B,2,1),L_MEme3(M,A,B,2,2))
#define L_StaticMatrix_OP_mult4x4(M,A,B) (L_MEme4(M,A,B,0,0),L_MEme4(M,A,B,0,1),L_MEme4(M,A,B,0,2),L_MEme4(M,A,B,0,3),L_MEme4(M,A,B,1,0),L_MEme4(M,A,B,1,1),L_MEme4(M,A,B,1,2),L_MEme4(M,A,B,1,3),L_MEme4(M,A,B,2,0),L_MEme4(M,A,B,2,1),L_MEme4(M,A,B,2,2),L_MEme4(M,A,B,2,3),L_MEme4(M,A,B,3,0),L_MEme4(M,A,B,3,1),L_MEme4(M,A,B,3,2),L_MEme4(M,A,B,3,3))
#define L_StaticMatrix_OP_mult4x4_0001(M,A,B) (L_MEme4(M,A,B,0,0),L_MEme4(M,A,B,0,1),L_MEme4(M,A,B,0,2),L_MEme4(M,A,B,0,3),L_MEme4(M,A,B,1,0),L_MEme4(M,A,B,1,1),L_MEme4(M,A,B,1,2),L_MEme4(M,A,B,1,3),L_MEme4(M,A,B,2,0),L_MEme4(M,A,B,2,1),L_MEme4(M,A,B,2,2),L_MEme4(M,A,B,2,3),(M)(3,0)=0,(M)(3,1)=0,(M)(3,2)=0,M(3,3)=1)

#define L_StaticMatrix_transpDe(M, A) {int df_i97, df_j97;for (df_i97=0; df_i97<(A).rows(); df_i97++) for (df_j97=0; df_j97<(A).cols(); df_j97++) (M)(df_j97,df_i97)=(A)(df_i97,df_j97);}

#define L_MEtre(M,A,i,j) ((M)(i,j)=(A)(j,i))
#define L_StaticMatrix_transpOf2x2(M,A) (L_MEtre(M,A,0,0),L_MEtre(M,A,0,1),L_MEtre(M,A,1,0),L_MEtre(M,A,1,1))
#define L_StaticMatrix_transpOf3x3(M,A) (L_MEtre(M,A,0,0),L_MEtre(M,A,0,1),L_MEtre(M,A,0,2),L_MEtre(M,A,1,0),L_MEtre(M,A,1,1),L_MEtre(M,A,1,2),L_MEtre(M,A,2,0),L_MEtre(M,A,2,1),L_MEtre(M,A,2,2))
#define L_StaticMatrix_transpOf4x4(M,A) (L_MEtre(M,A,0,0),L_MEtre(M,A,0,1),L_MEtre(M,A,0,2),L_MEtre(M,A,0,3),L_MEtre(M,A,1,0),L_MEtre(M,A,1,1),L_MEtre(M,A,1,2),L_MEtre(M,A,1,3),L_MEtre(M,A,2,0),L_MEtre(M,A,2,1),L_MEtre(M,A,2,2),L_MEtre(M,A,2,3),L_MEtre(M,A,3,0),L_MEtre(M,A,3,1),L_MEtre(M,A,3,2),L_MEtre(M,A,3,3))

#define L_StaticMatrix_identity(M) {int i_23,j_23; for (i_23=0; i_23<(M).rows(); i_23++) for (j_23=0; j_23<(M).cols(); j_23++) (M)(i_23,j_23)=(i_23==j_23);}
#define L_StaticMatrix_identity_val(M,a) {int i_23,j_23; for (i_23=0; i_23<(M).rows(); i_23++) for (j_23=0; j_23<(M).cols(); j_23++) (M)(i_23,j_23)=(i_23==j_23)*a;}
#define L_StaticMatrix_identity2x2(M) ((M)(0,0)=1, (M)(0,1)=0, (M)(1,0)=0, (M)(1,1)=1)
#define L_StaticMatrix_identity2x2_val(M,a) ((M)(0,0)=a, (M)(0,1)=0, (M)(1,0)=0, (M)(1,1)=a)
#define L_StaticMatrix_identity3x3(M) ((M)(0,0)=1, (M)(0,1)=0, (M)(0,2)=0, (M)(1,0)=0, (M)(1,1)=1, (M)(1,2)=0, (M)(2,0)=0, (M)(2,1)=0, (M)(2,2)=1)
#define L_StaticMatrix_identity3x3_val(M,a) ((M)(0,0)=a, (M)(0,1)=0, (M)(0,2)=0, (M)(1,0)=0, (M)(1,1)=a, (M)(1,2)=0, (M)(2,0)=0, (M)(2,1)=0, (M)(2,2)=a)
#define L_StaticMatrix_identity3x3s(M,i0,j0) ((M)(0+i0,0+j0)=1, (M)(0+i0,1+j0)=0, (M)(0+i0,2+j0)=0, (M)(1+i0,0+j0)=0, (M)(1+i0,1+j0)=1, (M)(1+i0,2+j0)=0, (M)(2+i0,0+j0)=0, (M)(2+i0,1+j0)=0, (M)(2+i0,2+i0)=1)

#define L_StaticMatrix_invert2x2(Mret,M) {double idet373 = 1/((M)(0,0)*(M)(1,1)-(M)(0,1)*(M)(1,0)); (Mret)(0,0) = (M)(1,1)*idet373; (Mret)(0,1)=-(M)(0,1)*idet373; (Mret)(1,0)=-(M)(1,0)*idet373; (Mret)(1,1)=M(0,0)*idet373;}
void L_invertMe(double *elem, int li, int lj);

// Contenedor para expresiones matematicas con matrices
// DEBE TENER TAMANO CERO, es decir, es solo una interfaz

template <class A>
class L_MatrizExpr {
public:
	operator const A&() const {return static_cast<const A&>(*this);} // Esta funcion permite recuperar una referencia a la expresion contenida
	const A& cont() const {return static_cast<const A&>(*this);} // Esta funcion permite recuperar una referencia a la expresion contenida
	const L_Matrix &e() const {return static_cast<const A&>(*this).e();} // Esta funcion evaluate la expresion en una matriz
};


// Aunque g++ es rapido con las funciones inline, Visual C++ es lento al evaluarlas
// Para llamadas rapidas en Visual C++, se reemplaza codeMapping tipo elemento-a-elemento
// El problema con los defines es que si se declara dentro una variable contenida en el argumento, no funcionan
#define L_Matrix_destroy(M) {(M).destroy()}
#define L_Matrix_reallocate(M, LI, LJ) {(M).reallocate(LI,LJ);}
#define L_Matrix_checkIfDefined(M, LI, LJ) ( throw_L_ArgException_if((M).rows()!=(LI) || (M).cols()!=(LJ) || M.begin()==NULL, "L_Matrix_revDef")) // build error si la matriz M no esta definida de tamano (li,lj)
#define L_Matrix_checkIfNoNULL(M) ( throw_L_ArgException_if((M).rows()==0 || (M).cols()==0 || (M).begin()==NULL, "L_Matrix_revNoNula") ) // build error si la matriz M no esta definida de tamano (li,lj)
#define L_Matrix_clean(M) {L_Matrix_checkIfNoNULL(M); L_StaticMatrix_clean(M);}
#define L_Matrix_OP_assign(M,A) {L_Matrix_reallocate(M,(A).rows(),(A).cols()); L_StaticMatrix_OP_assign(M,A);}
#define L_Matrix_OP_add(M,A,B) {L_Matrix_reallocate(M,(A).rows(),(A).cols()); L_StaticMatrix_OP_add(M,A,B);}
#define L_Matrix_OP_subtract(M,A,B) {L_Matrix_reallocate(M,(A).rows(),(A).cols()); L_StaticMatrix_OP_subtract(M,A,B);}
#define L_Matrix_OP_mult(M,A,B) {L_Matrix_reallocate(M,(A).rows(),(B).cols()); L_StaticMatrix_OP_mult(M,A,B);}
#define L_Matrix_transpOf(M,A) {L_Matrix_reallocate(M,(A).cols(),(A).rows()); L_StaticMatrix_transpDe(M,A);}
#define L_Matrix_OP_mult_sub(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2) {L_StaticMatrix_OP_mult_sub(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2);}
#define L_Matrix_OP_mult_sub_sum(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2) {L_StaticMatrix_OP_mult_sub_sum(M,i0,j0,m1,m2,i01,j01,li1,lj1,i02,j02,li2,lj2);}


#define L_Matrix_OP_amplify(M,val) {L_Matrix_checkIfNoNULL(M); L_StaticMatrix_OP_amplifica(M,val);}


template <int LI, int LJ>  // Esta tipo de matriz conviene para matrices pequenas, se evita la memoria dinamica
class L_StaticMatrix;

// Lo ideal seria usar __m128d para hacer los calculos...
// Los expression templates son muy complicados para usarlos aca...


class L_Matrix : public L_MatrixBase<double>
#ifdef USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED
: public L_MatrizExpr<L_Matrix>
#endif
{
public:
  typedef typename L_MatrixBase<double>::size_type size_type;
  using L_MatrixBase<double>::data;
  using L_MatrixBase<double>::iniptr;
  using L_MatrixBase<double>::endptr;
  using L_MatrixBase<double>::li;
  using L_MatrixBase<double>::lj;
  using L_MatrixBase<double>::operator();


private:
	L_Matrix(const L_Matrix &other);
public:

	template <int LI, int LJ> L_CommaListParser<double,LI*LJ> defineCommas() {reallocate(LI, LJ); return L_CommaListParser<double,LI*LJ>(iniptr);}

	const L_Matrix &expr() const {return *this;}

	L_Matrix() : L_MatrixBase<double>() { }

	L_Matrix(int li, int lj, const double *arr) : L_MatrixBase<double>()
	{
		int i, j, n=0;
		reallocate(li, lj);
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=arr[n++];
	}
	L_Matrix(int li, int lj) {reallocate(li, lj);}
	template <class A> L_Matrix(const L_MatrizExpr<A> &other) : L_MatrixBaseWidthStep<double>() {OP_assign(other.e());}

	double &elem_deb(int i, int j) {throw_L_ArgException_if(iniptr==NULL || i<0 || j<0 || i>=li || j>=lj , "L_Matrix::elem_deb() : indexes fuera de la matriz"); return operator()(i,j);}


#define OJO______USE_OPTIMIZED_MULTIPLICATION
#define OJO______USE_OPTIMIZED_EIGENDECOMPOSITION

#ifdef __COMPAT_UBLAS__
	boost::numeric::ublas::matrix<double> crearUblas() const;
	boost::numeric::ublas::matrix<double> crearUblasRef() const;  // createFrom una referencia, no copia la memoria
	L_Matrix &operator=(const boost::numeric::ublas::matrix<double> &other);
	void OP_mult_ublas(const L_Matrix &m1, const L_Matrix &m2); // Rapida x2
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_EIGEN__
	Eigen::MatrixXd crearEigen() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenLenta() const;  // createFrom una referencia a una L_Matrix
	Eigen::MatrixXd crearEigenRef() const;  // createFrom una referencia a una L_Matrix
	L_Matrix &operator=(const Eigen::MatrixXd &other);
	void OP_mult_eigen(const L_Matrix &m1, const L_Matrix &m2);   // MUY rapida x4
#endif

	void reallocate(int li0, int lj0) { reallocate_fc(li0, lj0); }
	void conservativeResize(int li0, int lj0) {conservativeResize_fc(li0, lj0, 0); }

	void identity() {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)=(i==j);}
	void identity(double diag) {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)=diag*(i==j);}
    void fill(const double *lista); // fill de izquierda a derecha por filas
	void setSubMatrixOffset(const L_Matrix subMatrix, int i0, int j0) {int i, j; for (i=0; i<subMatrix.li; i++) for (j=0; j<subMatrix.lj; j++) operator()(i+i0,j+j0) = subMatrix(i,j);}
	void concatenateSubMatrices(const L_Matrix &aa, const L_Matrix &ab, const L_Matrix &ba, const L_Matrix &bb); // Forma la matriz a partir de ni x nj submatrices
	void concatenateSubMatrices(int ni, int nj, const L_Matrix *m00, ...); // Forma la matriz a partir de ni x nj submatrices
	//void separarEnSubMatrices(int ni, int nj, L_Matrix *m00, ...); // Forma la matriz a partir de ni x nj submatrices
	//void copyVectorFrom(const std::vector<double> &vector) {reallocate(vector.size(), 1); for (int i=0; i<(int)vector.size(); i++) operator()(i,0) = vector[i];}
	void copyVectorFrom(const std::vector<double> &vector) {reallocate((size_type)vector.size(), 1); for (typename std::vector<double>::size_type i=0; i<vector.size(); i++) operator()((size_type)i,0) = vector[i];}
	void copyVectorTo(std::vector<double> &vector) {vector.resize(li); for (size_type i=0; i<(size_type)vector.size(); i++) vector[i] = operator()(i,0);}
	//void copyVectorTo(std::vector<double> &vector) {vector.resize(li); for (int i=0; i<(int)vector.size(); i++) vector[i] = operator()(i,0);}

	bool readFile(FILE *fp);
	bool readFile(const char *name) {FILE *fp; fp = fopen(name, "r"); if (fp == NULL) return false; if (readFile(fp) == false) return false; fclose(fp); return true;}
	bool saveFile(FILE *fp);
	bool saveFile(const char *name) {FILE *fp; fp = fopen(name, "w"); if (fp == NULL) return false; if (saveFile(fp) == false) return false; fclose(fp); return true;}

	double maxElem() {int i, j; double maxe = operator()(0,0); for (i=0; i<li; i++) for (j=0; j<lj; j++) maxe = (operator()(i,j) > maxe) ? operator()(i,j):maxe; return maxe;}
	double minElem() {int i, j; double mine = operator()(0,0); for (i=0; i<li; i++) for (j=0; j<lj; j++) mine = (operator()(i,j) < mine) ? operator()(i,j):mine; return mine;}
	double maxElemDiag() {int i; double maxe = operator()(0,0); for (i=0; i<li; i++) maxe = (operator()(i,i) > maxe) ? operator()(i,i):maxe; return maxe;}
	double minElemDiag() {int i; double mine = operator()(0,0); for (i=0; i<li; i++) mine = (operator()(i,i) < mine) ? operator()(i,i):mine; return mine;}
	double maxElemAbs() {int i, j; double maxe = operator()(0,0); for (i=0; i<li; i++) for (j=0; j<lj; j++) maxe = (L_ABS(operator()(i,j)) > maxe) ? L_ABS(operator()(i,j)):maxe; return maxe;}
	double minElemAbs() {int i, j; double mine = operator()(0,0); for (i=0; i<li; i++) for (j=0; j<lj; j++) mine = (L_ABS(operator()(i,j)) < mine) ? L_ABS(operator()(i,j)):mine; return mine;}
	double maxElemDiagAbs() {int i; double maxe = operator()(0,0); for (i=0; i<li; i++) maxe = (L_ABS(operator()(i,i)) > maxe) ? L_ABS(operator()(i,i)):maxe; return maxe;}
	double minElemDiagAbs() {int i; double mine = operator()(0,0); for (i=0; i<li; i++) mine = (L_ABS(operator()(i,i)) < mine) ? L_ABS(operator()(i,i)):mine; return mine;}


	// simetrica => vect propios ortogonales (PDPT)
	// pos def => D > 0       simetrica + pos def -> LLT con diag(L)>0
	bool isGoodCovarianceMatrix(double factorAsi=0.01, bool print = true, const char *name = NULL); // Test rapido, se debe agregar isPositiveDefinite()
	bool isGoodCovarianceMatrix(int &iErr, int &jErr, double factorAsi=0.01, bool print = true, const char *name = NULL); // Test rapido, se debe agregar isPositiveDefinite()
	bool isSymmetric(double precision) {int i, j; double d; for (i=0; i<li; i++) { for (j=0; j<lj; j++) { d=operator()(i,j)-operator()(j,i); if (d>precision || d<-precision) return false;}} return true;}
	bool isPositiveDefinite(); // Lenta...
	bool isNaN() {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) if (operator()(i,j) != operator()(i,j)) return true; return false;}
	double computeMaximumSubtraction(L_Matrix &C); // mayor lambda (aprox) tal que *this - lambda*C >= 0
	bool checkRangeValues(double limInf, double limSup, bool print=true, const char *name=NULL);
	void repairPositiveDefinite_valPr(); // slow
	//bool repararDefinidaPositivaSimetrica_LLT(); // Idea don't work
	//void repararDefinidaPositivaSimetrica_LLT_forzado(); // Idea don't work
	void repairSymmetry();
	void repairRangeValues(double limInf, double limMax);
	void repararNans(int i0=0, int j0=0, int li0=-1, int lj0=-1) {throw_L_ArgException_if(data()==NULL || i0+li0>li || j0+lj0>lj, "L_Matrix::repararNans() : matriz no creada"); int i, j; if (li0<0) li0=li; if(lj0<0) lj0=lj; for (i=0; i<li0; i++) for(j=0; j<lj0; j++) if (operator()(i+i0,j+j0) != operator()(i+i0,j+j0)) operator()(i+i0,j+j0) = 1.0e60;}

	// N = (1:n) - a;  P -> [P(a,a) P(a,N) ; P(N,a) P(N,N)] -> [J*P(a,a)*JT J*P(a,N) ; P(N,a)*JT P(N,N)] -> P
	void covarianceThroughJacobian(const std::vector<int> &a, const L_Matrix &J, const std::vector<int> &b);
	double covarianceThroughJacobianTest(const std::vector<int> &a, const std::vector<int> &b, const L_Matrix &J,  const L_Matrix *Ji = NULL, double umbral = -1, bool *aceptado=NULL, std::vector<L_StaticMatrix<3,3> > *cov3x3=NULL, L_Matrix *rot=NULL);
	void covarianceThroughJacobianTestNew(const std::vector<int> &a, const std::vector<int> &b, const L_Matrix &J, std::vector<L_StaticMatrix<3,3> > &cov3x3, bool forzar_maxima_cov = false);

	void fillRandomValues(double rango1, double rango2) {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)=rango1 + (rango2-rango1)*rand()/(double)RAND_MAX;}
	void integralOf(const L_Matrix &other);
	void integralSquaresOf(const L_Matrix &other);
	void integralDebugOf(const L_Matrix &other);
	void integralCuadradosDe_prueba(const L_Matrix &other);
	double sum_rectangle_using_integral(int i1, int j1, int i2, int j2) {throw_L_ArgException_if(i1 > i2 || j1 > j2, "L_Matrix::sum_rectangle_using_integral"); if (i1==0 && j1==0) {return operator()(i2,j2);} if (i1==0) {return operator()(i2,j2)-operator()(i2,j1-1);} if (j1==0) {return operator()(i2,j2)-operator()(i1-1,j2);} return operator()(i2,j2) - operator()(i1-1,j2) - operator()(i2,j1-1) + operator()(i1-1,j1-1);}
	double mean_rectangle_using_integral(int i1, int j1, int i2, int j2) {double N = (i2-i1+1)*(j2-j1+1), s; s = sum_rectangle_using_integral(i1,j1,i2,j2); return s/N;}
	double variance_rectangle_using_integral(int i1, int j1, int i2, int j2, L_Matrix &sumaCuad) {int N = (i2-i1+1)*(j2-j1+1); double s, s2; s = sum_rectangle_using_integral(i1,j1,i2,j2); s2 = sumaCuad.sum_rectangle_using_integral(i1,j1,i2,j2); return s2/N - (s/N)*(s/N);}

	void subMatrixOf(const L_Matrix &other, int i1, int i2, int j1, int j2) {int i, j; reallocate(i2-i1+1, j2-j1+1); for (i=i1; i<=i2; i++) for (j=j1; j<=j2; j++) operator()(i-i1,j-j1)=other(i,j);}  // MUY LENTO, es como escribir m(i1:i2,j1:j2) en matlab pero con indexes desde zero
	L_Matrix subMatrix(int i1, int i2, int j1, int j2) const {L_Matrix ret; ret.subMatrixOf(*this, i1, i2, j1, j2); return ret;}
	void concatenateHorz(L_Matrix &m1, L_Matrix &m2) {int i, j; reallocate(m1.li, m1.lj+m2.lj); for (i=0; i<li; i++) for (j=0; j<m1.lj; j++) operator()(i,j)=m1(i,j); for (i=0; i<li; i++) for (j=0; j<m2.lj; j++) operator()(i,m1.lj+j)=m2(i,j);}
	void concatenateVert(L_Matrix &m1, L_Matrix &m2) {int i, j; reallocate(m1.li+m2.li, m1.lj); for (i=0; i<m1.li; i++) for (j=0; j<lj; j++) operator()(i,j)=m1(i,j); for (i=0; i<m2.li; i++) for (j=0; j<lj; j++) operator()(m1.li+i,j)=m2(i,j);}
	void separateHorz(int lj1, L_Matrix &m1, L_Matrix &m2) {int i, j; m1.reallocate(li, lj1); m2.reallocate(li, lj-lj1); for (i=0; i<li; i++) for (j=0; j<lj1; j++) m1(i,j)=operator()(i,j); for (i=0; i<li; i++) for (j=0; j<lj-lj1; j++) m2(i,j)=operator()(i,j+lj1);}

	void setRotationComponents(const L_Matrix &matrRot) {int i, j; for (i=0; i<3; i++) for (j=0; j<3; j++) operator()(i,j)=matrRot(i,j);} // set componentes de rotacion
	void setTranslationComponents(const L_Matrix &vectTrans) {int i; for (i=0; i<3; i++) operator()(i,3)=vectTrans(i,0);} // set componentes de traslacion en la cuarta columna
	void getRotationComponents(L_Matrix &R) {R.reallocate(3,3); R(0,0)=operator()(0,0); R(0,1)=operator()(0,1); R(0,2)=operator()(0,2); R(1,0)=operator()(1,0); R(1,1)=operator()(1,1); R(1,2)=operator()(1,2); R(2,0)=operator()(2,0); R(2,1)=operator()(2,1); R(2,2)=operator()(2,2);}
	void getTranslationComponents(L_Matrix &t) {t.reallocate(3,1); t(0,0)=operator()(0,3); t(1,0)=operator()(1,3); t(2,0)=operator()(2,3);} // Pide domponentes de traslacion en la cuarta columna

	void moveObjectTo_NULL(L_Matrix &other)
    {
    	other.destroy();
		other.iniptr = iniptr; other.endptr = endptr; other.li = li; other.lj = lj;
		iniptr=NULL; endptr = NULL; li = 0; lj = 0;
    }

    L_Matrix& transpOf(const L_Matrix &im);
	L_Matrix& transposeMe() {L_Matrix other; other.transpOf(*this); other.swap(*this); return *this;}

	double trace();
	double traceRange(int i1, int i2);
	double det();
	double det_v0();
	double normL2() const {double sum=0; int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) sum+=operator()(i,j)*operator()(i,j); return sum;}
	double normL2between (const L_Matrix &other) const {throw_L_ArgException_if(li!=other.li||lj!=other.lj,"normL2between()");double sum=0; int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) sum+=operator()(i,j)*other(i,j); return sum;}
	void normalizeL2() {double factor=1/sqrt(normL2()); int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)*=factor;}
	double cosineDistanceTo(const L_Matrix &other) const {return normL2between(other) / sqrt(normL2()*other.normL2());}
	static double cosineDistanceBetweenArrays(double * L_restrict v1, double * L_restrict v2, int n) {int i; double m12=0, m11=0, m22=0; for(i=0; i<n; i++) {m12=v1[i]*v2[i]; m11=v1[i]*v1[i]; m22=v2[i]*v2[i];} return m12/sqrt(m11*m22);}

	void OP_changeSign(const L_Matrix &other);
	void subMatrDeleteRowCol_at(L_Matrix &ret, int i0, int j0);
	void insertRowsColumns(int i0, int ni, int j0, int nj, const L_Matrix &m); // Las columna i0 es "empujada" adelante en ni
	void deleteRowsColumns(int *indi, int ni, int *indj, int nj, const L_Matrix &m);
	void deleteRowsColumns(int i0, int ni, int j0, int nj, const L_Matrix &m);
	void copySubMatrixFrom(const L_Matrix &other, int i0, int j0) {int i, j; for (i=0; i<other.li; i++) for (j=0; j<other.lj; j++) operator()(i+i0,j+j0) = other(i,j);}

	void genDrawing(L_ImageRGBUchar &im); // Positivos en azul, negativos en verde

	void crossProductMatrix(double t1, double t2, double t3) {reallocate(3,3); operator()(0,0)=0; operator()(0,1)=-t3; operator()(0,2)=t2;  operator()(1,0)=t3; operator()(1,1)=0; operator()(1,2)=-t1;  operator()(2,0)=-t2; operator()(2,1)=t1; operator()(2,2)=0;} // [t)x
	void antiSymmetricElemental() {reallocate(3,3); operator()(0,0)=0; operator()(0,1)=1; operator()(0,2)=0;  operator()(1,0)=-1; operator()(1,1)=0; operator()(1,2)=0;  operator()(2,0)=0; operator()(2,1)=0; operator()(2,2)=1;}

	void OP_assign(const L_Matrix &other);
	void OP_assignInvertingSigns(const L_Matrix &other);
	void OP_add(const L_Matrix &other);
	void OP_subtract(const L_Matrix &other);
	void OP_mult(const L_Matrix &other);
	void OP_mult(const L_Matrix &other, double val);
	void OP_multByDiagonal(const L_Matrix &diagonal);
	void OP_div(const L_Matrix &other);
	void OP_add(const L_Matrix &m1, const L_Matrix &m2);
	void OP_subtract(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult(const L_Matrix &m1, const L_Matrix &m2);

	template <class A> void OP_subtract(const L_Matrix &m1, const L_MatrizExpr<A> &m2) {OP_subtract(m1, m2.e());}

	void OP_mult_private(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult_ABT(const L_Matrix &m1, const L_Matrix &m2); // Un poquito mas rapida
	void OP_mult_APAT(const L_Matrix &A, const L_Matrix &P) {L_Matrix AP; AP.OP_mult(A,P); OP_mult_ABT(AP,A);}
	void OP_mult_APAT_cholesky(const L_Matrix &A, const L_Matrix &P, bool reparar); // ahora que lo pienso, no es mucho mas lento que lo other
	void OP_mult_sorting(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult_longDouble(const L_Matrix &m1, const L_Matrix &m2);
	

	// LAS FUNCIONES sub TIENDEN A GENERAR ERRORES (pero son muy eficientes)
	void OP_mult_sub(const L_Matrix &m1, const L_Matrix &m2, int li1, int lj1, int li2, int lj2);
	void OP_mult_sub(int i0, int j0, const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1, int i02, int j02, int li2, int lj2);
	void OP_mult_sub_sum(int i0, int j0, const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1, int i02, int j02, int li2, int lj2);

	void OP_mult_sub1_sum(const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1);
	void OP_mult_sub2_sum(const L_Matrix &m1, const L_Matrix &m2, int i02, int j02, int li2, int lj2);

	void OP_mult_Abstr(const L_Matrix *m1, const L_Matrix *m2);
	void OP_multByDiagonal(const L_Matrix &normal, const L_Matrix &diagonal);
	void OP_div(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult(double val);
	void OP_amplify(double val); // Lo mismo
	void OP_div(double val);
	void OP_productElementwise(const L_Matrix &m1, const L_Matrix &m2);

	// Si algunas de las matrices son NULL, estas funciones la asumen = 0
	void OP_assign_NULL(const L_Matrix &other);
	void OP_assignInvertingSigns_NULL(const L_Matrix &other);
	void OP_add_NULL(const L_Matrix &other);
	void OP_subtract_NULL(const L_Matrix &other);
	void OP_mult_NULL(const L_Matrix &other) {L_Matrix aux; aux.OP_mult_NULL(*this, other); aux.swap(*this);}
	void OP_multByDiagonal_NULL(const L_Matrix &diagonal);
	void OP_div_NULL(const L_Matrix &other);
	void OP_add_NULL(const L_Matrix &m1, const L_Matrix &m2);
	void OP_subtract_NULL(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult_NULL(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult_Abstr_NULL(const L_Matrix *m1, const L_Matrix *m2);
	void OP_multByDiagonal_NULL(const L_Matrix &normal, const L_Matrix &diagonal);
	void OP_div_NULL(const L_Matrix &m1, const L_Matrix &m2);
	void OP_mult_NULL(double val);
	void OP_div_NULL(double val);
	void OP_productElementwise_NULL(const L_Matrix &m1, const L_Matrix &m2);

#ifdef __COMPAT_UBLAS__
	bool invertMe_ublas();
#endif
#ifdef __COMPAT_EIGEN__
	bool invertMe_eigen();
#endif
#ifdef __COMPAT_EIGEN__
	void vectPr_Eigen_sim(L_Matrix &valPr); // No logro hacer que funcione
#endif
#ifdef __COMPAT_UBLAS__
	bool cholesky_L_ublas_De(const L_Matrix &original); // Para matrices def positivas. Devuelve L, se tiene: original = L * LT
#endif
#ifdef __COMPAT_EIGEN__
	bool cholesky_L_eigen_De(const L_Matrix &original); // Para matrices def positivas. Devuelve L, se tiene: original = L * LT
#endif

	void copiaObjetoEn(L_Matrix &other) {other.OP_assign(*this);}

	L_Matrix& operator =(const L_Matrix &other) {OP_assign(other); return *this;}
    L_Matrix& operator +=(const L_Matrix &other){OP_add(other); return *this;}
    L_Matrix& operator -=(const L_Matrix &other){OP_subtract(other); return *this;}
    L_Matrix& operator *=(const L_Matrix &other){OP_mult(other); return *this;}
	L_Matrix& operator *=(double val) {OP_mult(val); return *this;}
	L_Matrix& operator /=(double val) {OP_div(val); return *this;}

	template <class A> L_Matrix &operator=(const L_MatrizExpr<A> &other) {OP_assign(other.e()); return *this;}
	template <class A> L_Matrix &operator+=(const L_MatrizExpr<A> &other) {OP_add(other.e()); return *this;}
	template <class A> L_Matrix &operator-=(const L_MatrizExpr<A> &other) {OP_subtract(other.e()); return *this;}
	template <class A> L_Matrix &operator*=(const L_MatrizExpr<A> &other) {OP_mult(other.e()); return *this;}

#if !defined(USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED) && defined(DEFINE_NONEFFICIENT_MATRIX_OPERATORS)
	L_Matrix operator +(const L_Matrix &otra2) const {L_Matrix ret; ret.OP_add(*this, otra2); return ret;}
	L_Matrix operator -(const L_Matrix &otra2) const {L_Matrix ret; ret.OP_subtract(*this, otra2); return ret;}
	L_Matrix operator *(const L_Matrix &otra2) const {L_Matrix ret; ret.OP_mult(*this, otra2); return ret;}
	L_Matrix operator /(const L_Matrix &otra2) const {L_Matrix ret; ret.OP_div(*this, otra2); return ret;}
	L_Matrix operator *(double val) const {L_Matrix ret=*this; ret.OP_mult(val);return ret;}
	L_Matrix operator /(double val) const {L_Matrix ret=*this; ret.OP_div(val);return ret;}
	L_Matrix operator -() const {int i, j; L_Matrix ret; ret.reallocate(li,lj); for (i=0; i<li; i++) for (j=0; j<lj; j++) ret[i][j]=-elem[i][j]; return ret;}
	L_Matrix operator ~() const {L_Matrix ret; ret.transpOf(*this); return ret;}
	inline friend L_Matrix operator *(double val, const L_Matrix &other) {L_Matrix ret=other; ret.OP_mult(val);return ret;}
#endif // USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED

	bool diagonalize_me() {if (li <= lj) return diagonalize_me_wide(); return diagonalizame_tall();}
	bool diagonalize_me_wide();
	bool diagonalizame_tall();
	void diagonalVectorOf(const L_Matrix &other) {int i; reallocate(other.li, 1); for (i=0; i<other.li; i++) operator()(i,0) = other(i,i);}
	void diagonalVectorOf(const double *d, int n) {int i; reallocate(n, 1); for (i=0; i<li; i++) operator()(i,0) = d[i];}
	void diagonalMatrixOf(const double *d, int n)  {int i, j; reallocate(n, n); for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j) = (i==j)*d[i];}
	void diagonalVectorAddMe(L_Matrix &diag) {int i; for (i=0; i<diag.li; i++) operator()(i,i) += diag(i,0);}
	bool solveLinearSystem_x_Ab_result_in_b(L_Matrix &A_destruible); // Resuelve Ax=b para un solo b. destroy A en el proceso. Usa metodo de Gauss. Inicialmente *this=b. Al final, *this=x.
	bool solveLinearSystem(const L_Matrix &A, const L_Matrix &b) {L_Matrix Adestr; Adestr=A; *this=b; return solveLinearSystem_x_Ab_result_in_b(Adestr);}
	bool invertMe(); // Invierte matriz usando Gauss-Jordan. Util para resolver Ax=b con muchos b distintos
	bool invertMe_private(); // Invierte matriz usando Gauss-Jordan. Util para resolver Ax=b con muchos b distintos
	bool invertMeDefPos(); // Invierte matriz definida positiva
	bool invertMeSimetricaDefPos(); // Depende de Cholesky...
	void invertMe_invertirValoresPropios(double minVal);
	bool invertMegj();  // Busca el mejor pivote, mas lenta, pero mas estable numericamente (para matrices raras)
	bool resolverSisLinealgj(const L_Matrix &A, const L_Matrix &b);
	bool inverseOf(const L_Matrix &other) {*this = other; return invertMe();}
	static void prueba_inversas_veloc();
	static void prueba_inversas();

	void espacioNuloPivoteo(const L_Matrix &original_ancha); // a=espacioNuloPivoteo(b); => b*a=0
	void gramSchmidtColumnas(); // Reflexiones de HouseHolder podrian funcionar mejor que esto

	bool pruebaIgualA(L_Matrix other, double tolerancia) {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) if (fabs(operator()(i,j)-other(i,j)) > tolerancia) return false; return true;}


	bool svd(L_Matrix &d_vector, L_Matrix &v); // Descompone la matriz (*this) en u*d*vT, con u de rotacion, d diagonal y v de rotacion. u queda almacenada en (*this).  Si esto demora 2, los vectores propios (simetrica def pos) demora 3.
	bool svd_0(L_Matrix &d_vector, L_Matrix &v); // Descompone la matriz (*this) en u*d*vT, con u de rotacion, d diagonal y v de rotacion. u queda almacenada en (*this).  Si esto demora 2, los vectores propios (simetrica def pos) demora 3.
	bool ordenaDiagonal(L_Matrix &Dvector, L_Matrix &V); // Para svd y vect propios
	bool svd_ordenado(L_Matrix &d_vector, L_Matrix &v) {if (svd(d_vector, v) == false) return false; ordenaDiagonal(d_vector, v); return true;}
	bool svd_reparaNegativos();
	bool svd_reparaEsencial();
	bool svd_reparaFundamental();
	bool svd_reparaOrtonormal();
	void prueba_svd(double error);
	void prueba_svd_ordenado(double error);
	void pruebasvd_cholesky_vectpr();
	bool cholesky_L_De(const L_Matrix &original); // Para matrices def positivas. Devuelve L, se tiene: original = L * LT
	void cholesky_L_forzado_De(const L_Matrix &original); // Forzado para matrices no definidas positivas
	bool unscented_desde_media_cov(const L_Matrix &media, const L_Matrix &cov, double alfa=1.0e-3, double kappa=0, double beta=2.0); // build 2*n+1 columnas que se deben propagar por f()
	bool unscented_hacia_media_cov(L_Matrix &fmedia, L_Matrix &fcov, double alfa=1.0e-3, double kappa=0, double beta=2.0); // Esto devuelve la propagacion de media y cov por f()
	bool inverseOf2x2(L_Matrix &M2x2) {throw_L_ArgException_if(M2x2.li!=2 || M2x2.lj!=2 || li!=2 || lj!=2, "L_Matrix::inv2x2De() : matrices no inicializadas"); double det = M2x2(0,0)*M2x2(1,1)-M2x2(0,1)*M2x2(1,0), idet; if (det==0) return false; idet=1/det; operator()(0,0) = M2x2(1,1)*idet; operator()(0,1) = -M2x2(0,1)*idet; operator()(1,0) = -M2x2(1,0)*idet; operator()(1,1) = M2x2(0,0)*idet; return true;} // Lo que hay que hacer para que el intelissense sea comodo... // ERROR: decia M2x2(2,2)
	bool svd2x2(double &ang1, double &ang2, double &e1, double &e2) const; // Hay que revisarla mas
	bool svd2x2(double &rot, double &eMax, double &eMin) const; // Hay que revisarla mas
	bool svd2x2_lento_a_mano(double &ang1, double &ang2, double &e1, double &e2, bool debuggg = false) const; // Hay que revisarla mas
	bool vectPr_2x2(L_Matrix &valPr, bool columna=true);
	void vectPr_sim(L_Matrix &valPr); // Los vectores propios (P) se reescriben en la misma matriz (A) -> A = P D P^-1 = P D PT, requiere P simetrica. Si svd demora 2, esto demora 3.
	void vectPr_jacobi_sim(L_Matrix &valPr); // Los vectores propios (P) se reescriben en la misma matriz (A) -> A = P D P^-1 = P D PT, requiere P simetrica. Si svd demora 2, esto demora 3.
	void PCA_ordenado(L_Matrix &P, L_Matrix &d_nx1) const;
	void PCA_proyectar(L_Matrix &P, int ndims=-1);

	bool interseccionRayos(L_Matrix &p1, L_Matrix &v1, L_Matrix &p2, L_Matrix &v2, double &dis1, double &dis2);

	// Estas funciones de FORTRAN fueron adaptadas al formato de C a mano (orden de filas y rango desde 0)
	static void blas_zrotg(double &CA, double &CB, double &C, double &S);
	static double blas_dnrm2(int N, const L_Matrix &X, int i0, int j0, int INCX);
	static double blas_zdotc(int N, const L_Matrix &ZX, int i0, int j0, int INCX, const L_Matrix &ZY, int i1, int j1, int INCY);
	static void linpack_zchud(L_Matrix &r, int ldr, int p, const L_Matrix &x, L_Matrix &z, int ldz, int nz, L_Matrix &y, L_Matrix &rho, L_Matrix &c, L_Matrix &s); // cholesky update de LINPACK, r triangular superior
	static void linpack_zchdd(L_Matrix &r, int ldr, int p, const L_Matrix &x, L_Matrix &z, int ldz, int nz, L_Matrix &y, L_Matrix &rho, L_Matrix &c, L_Matrix &s, int &info); // cholesky downdate de LINPACK, r triangular superior
	static void linpack_zchex(L_Matrix &r, int ldr, int p, int k, int l, L_Matrix &z, int ldz, int nz, L_Matrix &c, L_Matrix &s, int job);  // Intercambio de filas/columnas de la matriz original usando solo la representacion Cholesky, LINPACK, r triangular superior

	bool cholupdate_U(const L_Matrix &z, char op = '+');
	bool cholupdate_L(const L_Matrix &z, char op = '+');  // llama la funcion de Fortran pasada a C++

	// Funciones para restar covarianzas sin ensuciarlas, las primera es O(n^3),
	// la segunda es O(C.li^2) + O(n^2*Z.lj), la tercera es O(n^2*Z.lj)
	int restarCovarianzaRecomponer(L_Matrix &C, double factor = 1.0);  // P := P - C
	int restarCovarianzaRecomponer(L_Matrix &J, L_Matrix &C, double factor = 1.0);  // P := P - J*C*JT,  "J = [J|0]"
	int restarCovarianzaCholesky_U(L_Matrix &Z, double factor = 1.0);  // UT*U := UT*U - Z*ZT
	int restarCovarianzaCholesky_L(L_Matrix &Z, double factor = 1.0);  // L*LT := L*LT - Z*ZT

	

	void print(const char *name, FILE *fp=stdout) const;
	static void imprime_arreglo(double *v, int n, const char *name, FILE *fp=stdout);
	void printMatlabFormat(FILE *fp, const char *name) const;
	//void printMatlabFormat(FILE *fp, const char *name, int i) const {char str[100]; sprintf(str, "%s(%d)", name, i); printMatlabFormat(fp, name);}
	void printMatlabFormat(const char *arch, const char *name) const {FILE *fp=fopen(arch,"w"); if (fp!=NULL) {printMatlabFormat(fp, name); fclose(fp);}}
	double debugValue(int i, int j);
    ~L_Matrix() {destroy();}

	static void pruebaDefinesMatrices();
	static void pruebaSvd2x2();
	static void pruebaVectPr_2x2();
	static void prueba_det();
	static void benchmarks();


private:
	static bool svdcmp_0(double **A, int M, int N, int MP, int NP, double *W, double **V );
	static void Jacobi_Cyclic_Method(double *valpr, double *vectpr, double *A, int n);
};



template <int LI, int LJ>  // Esta tipo de matriz conviene para matrices pequenas, se evita la memoria dinamica
class L_StaticMatrix  // Sin embargo, para matrices muy grandes (>10x10) es poco eficiente
{
public:
#ifdef L_BRACKET_LEVEL_DEBUG
	L_ptrME<LI,LJ> _elem;
#else
	double _elem[LI][LJ];
#endif
	L_StaticMatrix() {} // No necesita constructor de copia ni operator= ni destructor

	//   Dado que estas matrices suelen ser pequenas, para maximizar
	// la velocidad se deben usar las macros listadas arriba.

	typedef double * iterator;
	typedef const double * const_iterator;
	typedef int size_type;
	typedef size_t size_type2;

	double *data() {return &_elem[0][0];}
	double *begin() {return &_elem[0][0];}
	double *end() {return &(_elem[0][0]) + rows()*cols();}

	size_type rows() const {return LI;}
	size_type cols() const {return LJ;}
	size_type2 size() const {return rows()*cols()*sizeof(double);}

	const double *data() const {return &_elem[0][0];}
	const double *begin() const {return &_elem[0][0];}
	const double *end() const {return &(_elem[0][0]) + rows()*cols();}

	double &operator()(int i, int j) {return _elem[i][j];}
	const double &operator()(int i, int j) const {return _elem[i][j];}

	L_CommaListParser<double,LI*LJ> defineCommas() {return L_CommaListParser<double,LI*LJ>(&(_elem[0][0]));}

	double &elem_deb(int i, int j) {throw_L_ArgException_if(i<0 || j<0 || i>=rows() || j>=cols() , "L_StaticMatrix::elem_deb() : indexes fuera de la matriz"); return _elem[i][j];}

	void swap(L_StaticMatrix<LI,LJ> &other) {L_StaticMatrix<LI,LJ> t = other; other = *this; *this = t;}

	bool inverseOf(const L_StaticMatrix &other) {*this = other; return invertMe();}

	void copyOn(L_Matrix &m) {m.conservativeResize(LI,LJ); for (int i=0; i<LI; i++) for (int j=0; j<LJ; j++) m(i,j) = (*this)(i,j);}

	bool invertMe()  // Tal vez la saque
	{
		int i, j, c;
		double factor;
		double temp;
		throw_L_ArgException_if(rows()!=cols(), "L_StaticMatrix::invertMe()");
		// Crear matriz unitaria other
		L_StaticMatrix<LI,LJ> other;
		for (i=0; i<rows(); i++)
			for (j=0; j<cols(); j++)
				other._elem[i][j] = (i==j);

		for (c=0; c<cols(); c++) // Limpiar columna c
		{
			if (_elem[c][c]==0) // Diagonal con zero, intercambiar por other fila mas abajo que c
			{
				for (i=c+1; i<rows(); i++) // Encontrar fila i con _elemento != 0 en la columna c
					if (_elem[i][c]!=0)
						break;
				if (i == rows()) // fila i solo contiene _elementos zero para j>=c
					return false;
				for (j=0; j<c; j++) // intercambiar fila c por fila i. Notar que _elem[i][j] = 0 para j<c
				{
					temp = other._elem[i][j];
					other._elem[i][j]=other._elem[c][j];
					other._elem[c][j] = temp;
				}
				for (; j<cols(); j++) // intercambiar
				{
					temp = other._elem[i][j];
					other._elem[i][j]=other._elem[c][j];
					other._elem[c][j] = temp;
					temp = _elem[i][j];
					_elem[i][j]=_elem[c][j];
					_elem[c][j] = temp;
				}
			}
			// Amplificar fila c para que quede con 1 en la diagonal
			factor = 1 / _elem[c][c];
			for (j=0; j<c; j++)
			{
				other._elem[c][j]*=factor;
			}
			for (; j<cols(); j++)
			{
				_elem[c][j]*=factor;
				other._elem[c][j]*=factor;
			}
			// Restarla a las otras filas para dejar la columna c "unitaria"
			for (i=0; i<cols(); i++)
			{
				if (i==c)
					continue;
				factor = -_elem[i][c];
				for (j=0; j<c; j++)
					other._elem[i][j]+=factor*other._elem[c][j];
				for (; j<cols(); j++)
				{
					_elem[i][j]+=factor*_elem[c][j];
					other._elem[i][j]+=factor*other._elem[c][j];
				}
			}
		}
		// La inversa esta en other
		for (i=0; i<rows(); i++)
			for (j=0; j<cols(); j++)
				_elem[i][j] = other._elem[i][j];
		return true;
	}

	void print(const char *name, FILE *fp=stdout) const
	{
		L_String str;
		L_String str2;
		int lar=0;
		int i, j;
		bool imprNom=false;

		if (name!=NULL)
		{
			lar=(int)( strlen(name)+3 );
			str.resize(lar+1);
			str2.resize(lar+1);
			strcpy(str.data(), name);
			strcat(str.data(), " = ");
			for (i=0; i<lar; i++)
				str2[i]=' ';
			str2[lar]=0;
		}

		for (i=0; i<rows(); i++)
		{
			if (str.size()!=0)
			{
				if (imprNom==false && i>=rows()/2)
				{
					fprintf(fp, "%s", str.c_str());
					imprNom=true;
				}
				else
					fprintf(fp, "%s", str2.c_str());
			}
			fprintf(fp,"[");
			for (j=0; j<cols(); j++)
				fprintf(fp, "%+.3g  ", _elem[i][j]);
			fprintf(fp,"]\n");
		}
		fprintf(fp,"\n");
	}

	void printMatlabFormat(FILE *fp, const char *name) const
	{
		int i, j;
		fprintf(fp, "%s = [", name);
		for (i=0; i<rows(); i++)
		{
			for (j=0; j<cols(); j++)
			{
				if (_elem[i][j] == _elem[i][j])
					fprintf(fp, "%.10g, ", _elem[i][j]);
				else
					fprintf(fp, "NaN, ");
			}
			fprintf(fp,";");
		}
		fprintf(fp,"];\n");
	}
};


// Expression templates
/////////////////////////////////////////////////////////////////////////////

#ifdef USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED
//----------------------------------------------
template <class A, class B>
class L_MatrizSuma : public L_MatrizExpr<L_MatrizSuma<A,B> > {
	const A &a_;
	const B &b_;
	mutable L_Matrix tmp;
public:
	L_MatrizSuma(const A a, const B b) : a_(a),b_(b) {}
	const L_Matrix &e() const {
		tmp.OP_add(a_.e(), b_.e()); return tmp;
	}
};

//----------------------------------------------
template <class A, class B>
inline L_MatrizSuma<A,B>
operator+(const L_MatrizExpr<A>& a, const L_MatrizExpr<B>& b){
	return L_MatrizSuma<A,B>(a,b);
}

//----------------------------------------------
template <class A, class B>
class L_MatrizResta : public L_MatrizExpr<L_MatrizResta<A,B> > {
	const A &a_;
	const B &b_;
	mutable L_Matrix tmp;
public:
	L_MatrizResta(const A& a, const B& b) : a_(a), b_(b){}
	const L_Matrix &e() const {
		tmp.OP_subtract(a_.e(), b_.e()); return tmp;
	}
};
//----------------------------------------------
template <class A, class B>
inline L_MatrizResta<A,B>
operator-(const L_MatrizExpr<A>& a, const L_MatrizExpr<B>& b){
	return L_MatrizResta<A,B>(a,b);
};
 
template <class A>
class L_MatrizInvSigno : public L_MatrizExpr<L_MatrizInvSigno<A> > {
	const A &a_;
	mutable L_Matrix tmp;
public:
	L_MatrizInvSigno(const A& a) : a_(a){}
	const L_Matrix &e() const {
		tmp.OP_changeSign(a_.e()); return tmp;
	}
	const A & no_aplicar() const {return a_;}
};
//----------------------------------------------
template <class A>
inline L_MatrizInvSigno<A>
operator- (const L_MatrizExpr<A>& a) {
	return L_MatrizInvSigno<A>(a);
}

// doble signo -
template <class A>
inline A
operator- (const L_MatrizInvSigno<A>& a) {
	return a.no_aplicar();
}
 //----------------------------------------------
template <class A>
class L_MatrizMult : public L_MatrizExpr<L_MatrizMult<A> > {
	const A &a_;
	double b_;
	mutable L_Matrix tmp;
public:
	L_MatrizMult(const A& a, const double& b)
	: a_(a), b_(b){}
	const L_Matrix &e() const {
		tmp.OP_mult(a_.e(),b_); return tmp;
	}
};
//----------------------------------------------
template <class A>
inline L_MatrizMult<A>
operator* (const L_MatrizExpr<A>& a, double b) {
	return L_MatrizMult<A>(a,b);
}

template <class A>
inline L_MatrizMult<A>
operator* (double b, const L_MatrizExpr<A>& a) {
	return L_MatrizMult<A>(a,b);
}

//----------------------------------------------
template <class A>
inline L_MatrizMult<A>
operator* (int b, const L_MatrizExpr<A>& a) {
	return L_MatrizMult<A>(a,b);
}
//----------------------------------------------
template <class A, class B>
class L_MatrizMultM : public L_MatrizExpr<L_MatrizMultM<A,B> > {
	const A &a_;
	const B &b_;
	mutable L_Matrix tmp;
public:
	L_MatrizMultM(const A& a, const B& b)
	: a_(a), b_(b){}
	const L_Matrix &e() const {
		tmp.OP_mult(a_.e(),b_.e()); return tmp;
	}
};
//----------------------------------------------
template <class A, class B>
inline L_MatrizMultM<A,B>
operator* (const L_MatrizExpr<A>& a, const L_MatrizExpr<B>& b) {
	return L_MatrizMultM<A,B>(a,b);
}


//----------------------------------------------
template <class A>
class L_MatrizTransp : public L_MatrizExpr<L_MatrizTransp<A> > {
	const A &a_;
	mutable L_Matrix tmp;
	public:
	L_MatrizTransp(const A& a, void *nada) : a_(a){}
	const L_Matrix &e() const {
		tmp.transpOf(a_.e()); return tmp;
	}
	const A & no_aplicar() const {return a_;}
};

template <class A>
inline L_MatrizTransp<A>
operator~ (const L_MatrizExpr<A>& a) {
	return L_MatrizTransp<A>(a, NULL);
}

inline L_MatrizTransp<L_Matrix>
operator~ (const L_Matrix& a) {
	return L_MatrizTransp<L_Matrix>(a, NULL);
}

// Doble transpuesta
template <class A>
inline A
operator~ (const L_MatrizTransp<A>& a) {
	return a.no_aplicar();
}

//----------------------------------------------
template <class A, class B>
class L_MatrizMultM_ABT : public L_MatrizExpr<L_MatrizMultM_ABT<A,B> > {
	const A &a_;
	const B &b_;
	mutable L_Matrix tmp;
public:
	L_MatrizMultM_ABT(const A& a, const B& b)
	: a_(a), b_(b){}
	const L_Matrix &e() const {
		tmp.OP_mult_ABT(a_.e(),b_.e()); return tmp;
	}
};
//----------------------------------------------
template <class A, class B>
inline L_MatrizMultM_ABT<A,B>
operator* (const L_MatrizExpr<A>& a, const L_MatrizTransp<B>& b) {
	return L_MatrizMultM_ABT<A,B>(a,b.no_aplicar());
}

//----------------------------------------------
#endif // USE_EXPRESSION_TEMPLATES_NON_IMPLEMENTED
/////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////


inline void ejemploUsoMatriz()
{
	L_Matrix m, m2, m3;
	m.defineCommas<3,2>() =  1, 2,
							3, 4,
							5, 6;
	m.print("m");
	m2.reallocate(2, 2);
	m2(0,0) = 4;
	m2(0,1) = 3;
	m2(1,0) = 1;
	m2(1,1) = 2;
	std::sort(m2.begin(), m2.end());
	m2.print("m2");
#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	m3 = m * ~m2;
#else
	L_Matrix m2t;
	m2t.transpOf(m2);
	m3.OP_mult(m, m2t);
#endif
	m3.print("m*transp(m2)");
}






// Esta class_name es solo un lector, la memoria se pide por fuera
class L_ArregloFortranInt
{
private:
	int *elem;
public:
#ifdef L_ARRAY_ACCESS_DEBUG
	int nmem;
	L_ArregloFortranInt() {elem=NULL; nmem=0;}
	int &operator() (int i) {throw_L_ArgException_if (i<1 || (i-1)>=nmem || elem==NULL, "L_ArregloFortranInt::operator()");return elem[i-1];}
	void fijarRef(std::vector<int> &arr) {elem = &(arr[0]); nmem=(int)arr.size();}
	void fijarRefNula() {elem=NULL; nmem=0;}
#else
	int &operator() (int i) {return elem[i-1];}
	void fijarRef(std::vector<int> &arr) {elem = &(arr[0]);}
	void fijarRefNula() {elem=NULL;}
#endif
};

//#define L_FORTRAN_INTERFAZ

// Esta class_name es solo un lector de matriz traspuesta que parte de (1,1)
// Sirve para arreglos y matrices double precision Fortran
class L_MatrizFortranTr
{
private:
	double *elem; // referencia, no maneja la memo
	int li;
public:
#ifdef L_MATRIX_ACCESS_DEBUG
	int nmem;
	L_MatrizFortranTr() {elem=NULL; li = 0; nmem=0;}
	double &operator() (int i) {throw_L_ArgException_if(i<1 || (i-1) >= nmem || elem==NULL, "L_MatrizFortranTr:operator()");return elem[i-1];}
	double &operator() (int i, int j) {throw_L_ArgException_if(i<1 || j<1 || (i-1)>=li || (j-1)*li+(i-1) >= nmem || elem==NULL, "L_MatrizFortranTr:operator()"); return elem[(j-1)*li+(i-1)];}
	L_MatrizFortranTr ref(int i) {L_MatrizFortranTr ret; ret.elem=&elem[i-1]; ret.li=li; ret.nmem = nmem-(i-1); return ret;}
	L_MatrizFortranTr ref(int i, int j) {L_MatrizFortranTr ret; ret.elem=&elem[(j-1)*li+(i-1)]; ret.li=li; ret.nmem = nmem-((j-1)*li+(i-1)); return ret;}
	void dimension(int li) {this->li = li;}
	void dimension(int li, int lj) {this->li = li; this->nmem = li*lj;}
	void fijarRefTraspuesta(L_Matrix &other) {li = other.lj; elem = other.data(); nmem=other.li*other.lj;}
	void fijarRefTraspuesta(const L_Matrix &other) {li = other.lj; elem = const_cast<double *>(other.data()); nmem=other.li*other.lj;} // Es lo que hay no mas
	void fijarRefNula() {elem=NULL; li=0; nmem=0;}
#else
	double &operator() (int i) {return elem[i-1];}
	double &operator() (int i, int j) {return elem[(j-1)*li+(i-1)];}
	L_MatrizFortranTr ref(int i) {L_MatrizFortranTr ret; ret.elem=&elem[i-1]; ret.li=li; return ret;}
	L_MatrizFortranTr ref(int i, int j) {L_MatrizFortranTr ret; ret.elem=&elem[(j-1)*li+(i-1)]; ret.li=li; return ret;}
	void dimension(int li) {this->li = li;}
	void dimension(int li, int lj) {this->li = li;}
	void fijarRefTraspuesta(L_Matrix &other) {li = other.lj; elem = other.data();}
	void fijarRefTraspuesta(const L_Matrix &other) {li = other.lj; elem = const_cast<double *>(other.data());} // Es lo que hay no mas
	void fijarRefNula() {elem=NULL; li=0;}
#endif

	// Estas funciones fueron traducidas a mano de fortran a ratfor, el cual el compilador C++ puede entender.
	// http://www.netlib.org/lapack/  las originales en fortran
	static double dsqrt(double d) {return sqrt(d);}
	static void zrotg(double &CA, double &CB, double &C, double &S);
	static double dnrm2(int N, L_MatrizFortranTr X, int INCX);
	static double zdotc(int N, L_MatrizFortranTr ZX, int INCX, L_MatrizFortranTr ZY, int INCY);
	static double ddotc(int N, L_MatrizFortranTr ZX, int INCX, L_MatrizFortranTr ZY, int INCY) {return zdotc(N,ZX,INCX,ZY,INCY);}
	static void dswap(int N, L_MatrizFortranTr DX, int INCX, L_MatrizFortranTr DY, int INCY);
	static void daxpy(int N, double DA, L_MatrizFortranTr DX, int INCX, L_MatrizFortranTr DY, int INCY);
	static void zchud(L_MatrizFortranTr r, int ldr, int p, L_MatrizFortranTr x, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr y, L_MatrizFortranTr rho, L_MatrizFortranTr c, L_MatrizFortranTr s); // cholesky update de LINPACK, r triangular superior
	static void zchdd(L_MatrizFortranTr r, int ldr, int p, L_MatrizFortranTr x, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr y, L_MatrizFortranTr rho, L_MatrizFortranTr c, L_MatrizFortranTr s, int &info); // cholesky downdate de LINPACK, r triangular superior
	static void zchex(L_MatrizFortranTr r, int ldr, int p, int k, int l, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr c, L_MatrizFortranTr s, int job);  // Intercambio de filas/columnas de la matriz original usando solo la representacion Cholesky, LINPACK, r triangular superior
	static void dchdc(L_MatrizFortranTr a, int lda, int p, L_MatrizFortranTr work, L_ArregloFortranInt &jpvt, int job, int &info);
	static int ilaenv(int ispec, char *name, char *opts, int n1, int n2, int n3, int n4);
	static void dgesvd(char *jobu, char *jobvt, int m, int n, L_MatrizFortranTr a, int lda, L_MatrizFortranTr s, L_MatrizFortranTr u, int ldu, L_MatrizFortranTr vt, int ldvt, L_MatrizFortranTr work, int lwork, int &info);
	static void dsyev(char *jobz, char *uplo, int n, L_MatrizFortranTr a, int lda, L_MatrizFortranTr w, L_MatrizFortranTr work, int lwork, int &info);

	static void prueba_daxpy();

};

class L_MatrixFloat : public L_MatrixBase<float>
{
public:
	inline void copyFrom(const L_Matrix &other) {int i, j; reallocate(other.li, other.lj); for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)=(float)other(i,j);}
	inline void copyTo(L_Matrix &other) const {int i, j; other.reallocate(other.li, other.lj); for (i=0; i<li; i++) for (j=0; j<lj; j++) other(i,j)=operator()(i,j);}

	void reallocate(int li, int lj) {reallocate_fc(li, lj);}

	void OP_mult(const L_MatrixFloat &m1, const L_MatrixFloat &m2);
	void OP_mult_private(const L_MatrixFloat &m1, const L_MatrixFloat &m2);
	void resize(int li2, int lj2) {conservativeResize_fc(li2,lj2,0);}
	bool invertMeDefPos(); // Invierte matriz definida positiva. No es notoriamente mas rapido que con double
	bool invertMeDefPosx4(); // Invierte matriz definida positiva (unrollx4). No es notoriamente mas rapida que la de arriba

#ifdef __COMPAT_UBLAS__
	boost::numeric::ublas::matrix<float> crearUblas() const;
	boost::numeric::ublas::matrix<float> crearUblasRef() const;  // createFrom una referencia, no copia la memoria
	L_MatrixFloat &operator=(const boost::numeric::ublas::matrix<float> &other);
	bool invertMe_ublas();
	void OP_mult_ublas(const L_MatrixFloat &m1, const L_MatrixFloat &m2); // MUY rapida
#endif // __COMPAT_UBLAS__

};



class L_MatrizLongDouble : public L_MatrixBase<long double>
{
public:
	inline void copyFrom(const L_Matrix &other) {int i, j; reallocate(other.li, other.lj); for (i=0; i<li; i++) for (j=0; j<lj; j++) operator()(i,j)=(float)other(i,j);}
	inline void copyTo(L_Matrix &other) const {int i, j; other.reallocate(other.li, other.lj); for (i=0; i<li; i++) for (j=0; j<lj; j++) other(i,j)=operator()(i,j);}

	void reallocate(int li, int lj) {reallocate_fc(li, lj);}

	void OP_mult(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2);
	void OP_mult_private(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2);
	void OP_mult_ABT(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2);
	void resize(int li2, int lj2) {conservativeResize_fc(li2,lj2,0);}
	bool invertMeDefPos(); // Invierte matriz definida positiva. No es notoriamente mas rapido que con double

#ifdef __COMPAT_UBLAS__
	boost::numeric::ublas::matrix<long double> crearUblas() const;
	boost::numeric::ublas::matrix<long double> crearUblasRef() const;  // createFrom una referencia, no copia la memoria
	L_MatrizLongDouble &operator=(const boost::numeric::ublas::matrix<long double> &other);
	bool invertMe_ublas();
	void OP_mult_ublas(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2); // MUY rapida
#endif // __COMPAT_UBLAS__

};



class L_Matrix_intercambia // Matriz en que las filas tienen memoria independiente y se pueden intercambiar los punteros
{
public:
	double **elem;
	int li, lj;
	L_Matrix_intercambia() {elem = NULL; li = 0; lj = 0;}
	void reallocate(int li0, int lj0);
	void destroy();
	double &operator() (int i, int j) {return elem[i][j];}
	const double &operator() (int i, int j) const {return elem[i][j];}
	void OP_assign(const L_Matrix &other);
	void copyTo(L_Matrix &m);
	void identity() {int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) elem[i][j] = (i==j);}
	bool gaussj(L_Matrix_intercambia &B);
	~L_Matrix_intercambia() {destroy();}
};

// Los mismos defines-macro de L_Matrix deberian funcionar aca, excepto los que incluyen L_Matrix_reallocate()

class L_MegaMatriz
{
public:
	L_Matrix **m;  // Si m[i][j].elem == NULL, se asume fill de ceros
	int ni, nj;
	bool tamanosFijados;

	typedef int size_type;

	L_MegaMatriz() {m=NULL; ni=0; nj=0; tamanosFijados = false;}

	void destroy();
	void reallocate(int ni, int nj); // hace que las m[][] existan, aunque inicialmente van a tener m[][].elem == NULL

	bool coherent(std::vector<int> *lis=NULL, std::vector<int> *ljs=NULL, int *numTotI=NULL, int *numTotJ=NULL) const; // devuelve true si la matriz no tiene submatrices de tamano zero
	bool pideTamanos(std::vector<int> &lis, std::vector<int> &ljs, int *numTotI=NULL, int *numTotJ=NULL);
	void reconstruyeFijaTamanos(const std::vector<int> &lis, const std::vector<int> &ljs);

	void pideMatriz(L_Matrix &resultado, std::vector<int> *lis=NULL, std::vector<int> *ljs=NULL) const;
	void fijaMatriz(const L_Matrix &origen, const std::vector<int> &lis, const std::vector<int> &ljs, double zero = 0.0);
	void fijaMatrizDeltas(const L_Matrix &origen, int deltai, int deltaj);

	void OP_assign(L_MegaMatriz &other) {reallocate(other.ni, other.nj); int i, j; for (i=0; i<ni; i++) for (j=0; j<nj; j++) m[i][j].OP_assign_NULL(other.m[i][j]); tamanosFijados = other.tamanosFijados;}
	void OP_add(const L_MegaMatriz &a, const L_MegaMatriz &b); // Supone que son compatibles
	void OP_subtract(const L_MegaMatriz &a, const L_MegaMatriz &b); // Supone que son compatibles
	void OP_mult(const L_MegaMatriz &a, const L_MegaMatriz &b); // Supone que son compatibles
	void inverseOf(const L_MegaMatriz &other);
	void invertMe();
	void traspuestaDe(const L_MegaMatriz &other);

	void swap(L_MegaMatriz &other)
	{
		L_Matrix **m_t = other.m;
		int ni_t = other.ni;
		int nj_t = other.nj;
		bool tamanosFijados_t = other.tamanosFijados;
		
		other.m = m;
		other.ni = ni;
		other.nj = nj;
		other.tamanosFijados = tamanosFijados;
		
		m = m_t;
		ni = ni_t;
		nj = nj_t;
		tamanosFijados = tamanosFijados_t;
	}
	void copiaObjetoEn(L_MegaMatriz &other) {other.reallocate(ni, nj); int i, j; for (i=0; i<ni; i++) for (j=0; j<nj; j++) m[i][j].copiaObjetoEn(other.m[i][j]); other.tamanosFijados = tamanosFijados;}

	void insertaFilaColumna_regala(int i0, int li, int j0, int lj, const L_MegaMatriz &other);
	void eliminaFilasColumnas_regala(int *indi, int ni, int *indj, int nj, const L_MegaMatriz &other);
	void eliminaFilasColumnas_regala(int i0, int ni, int j0, int nj, const L_MegaMatriz &other);


	~L_MegaMatriz() {destroy();}
};

class L_MatrizFundamental:public L_Matrix
{
public:
	void createFrom(L_Matrix &H);  // Interfaz para convencion de ejes de robotica (ade,izq,arr)
	void createFrom(L_Matrix &RDAA, L_Matrix &tDAA); // Interfaz para convencion de ejes de vision (der,aba,ade)

	double calcError(double der1, double aba1, double der2, double aba2);
	double calcErrorProm_qTEq(L_Matrix &u1v1_u2v2); // u = (der,aba)
	bool reparaFundamental() {return svd_reparaFundamental();}
	void dibLinsF(L_ShapeArray &lins, int lx, int ly);
};

class L_EssentialMatrix:public L_MatrizFundamental // Es de 3x3 con valores singulares {x,x,0} si es esencial, {x,y,0} si es fundamental
{
public:
	// Funciones nuevas confiables
	bool descompone_Rt(double tanIzq1, double tanArr1, double tanIzq2, double tanArr2, L_HomogeneousMatrix &H); // Interfaz para convencion de ejes de robotica
	bool descompone_Rt(double tanDer1, double tanAba1, double tanDer2, double tanAba2, L_Matrix &RFinalDAA, L_Matrix &tFinalDAA, double *c1x=NULL, double *c2x=NULL); // Interfaz para convencion de ejes de vision

	static bool triangula(L_Matrix &RaDAA, L_Matrix &tDAA, double tanDer1, double tanAba1, double tanDer2, double tanAba2, double &xDer, double &yAba, double &zAde);
	bool reparaEsencial() {return svd_reparaEsencial();}
	void dibLinsE(L_ShapeArray &lins, double fx, double fy, int lx, int ly);
};


//#define DEBUG_STURM

namespace Sturm
{
// codeMapping copiado para poder calcular raices de un polinomio
// Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
// by D.G. Hook and P.R. McAree
// from "Graphics Gems", Academic Press, 1990
#define	 STURM_MAX_ORDER 20
// maxElement orden para los polinomios
#define	 STURM_RELERROR	 1.0e-14
// El menor error relativo aceptable
#define	 STURM_MAXPOW	 32
// mayor potencia de 10 que se usa
#define	 STURM_MAXIT	 800
// Numero maxElement de iteraciones
#define	 STURM_SMALL_ENOUGH	 1.0e-12
// Umbral para dejar los coeficientes igual a zero

#ifdef DEBUG_STURM
	bool main_sturm(int ord, double *coefs, int &nroots, L_ptr<double> &roots); // ord <= STURM_MAX_ORDER, coefs de tamano ord+1; roots de tamano n, coefs[0] es la constante...
#else
	bool main_sturm(int ord, double *coefs, int &nroots, double *roots); // ord <= STURM_MAX_ORDER, coefs de tamano ord+1; roots de tamano n, coefs[0] es la constante...
#endif
};


class L_PolinGr11
{
public:
	std::vector<double> c;
	int gr;

	L_PolinGr11() : c(12) {gr=0; c[0]=0; c[1]=0; c[2]=0; c[3]=0; c[4]=0; c[5]=0; c[6]=0; c[7]=0; c[8]=0; c[9]=0; c[10]=0; c[11]=0;}
	L_PolinGr11(const L_PolinGr11 &other) : c(12) {gr=0; for (int i=0; i<12; i++) c[i] = other.c[i]; gr = other.gr;}
	L_PolinGr11 &operator=(const L_PolinGr11 &other) {gr=0; for (int i=0; i<12; i++) c[i] = other.c[i]; gr = other.gr; return *this;}

	friend L_PolinGr11 operator*(const L_PolinGr11 &p1, const L_PolinGr11 &p2);
	friend L_PolinGr11 operator+(const L_PolinGr11 &p1, const L_PolinGr11 &p2);
	friend L_PolinGr11 operator-(const L_PolinGr11 &p1, const L_PolinGr11 &p2);

	void print(FILE *fp = stdout) const {for (int i=gr; i>=1; i--) fprintf(fp, "%.3g*z^%d + ", c[i], i); fprintf(fp, "%.3g", c[0]);}

	L_PolinGr11 operator+=(const L_PolinGr11 &p) {if (p.gr > gr) gr=p.gr; for (int i=0; i<=gr; i++) c[i]+=p.c[i]; return *this;}

	L_PolinGr11 operator*(int d) const {L_PolinGr11 p; p.gr=gr; for (int i=0; i<=gr; i++) p.c[i]=c[i]*d; return p;}
	L_PolinGr11 operator*(double d) const {L_PolinGr11 p; p.gr=gr; for (int i=0; i<=gr; i++) p.c[i]=c[i]*d; return p;}

	void crear(double c4, double c3, double c2, double c1, double c0) {gr=4; c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3; c[4]=c4;} // Error: decia c[4] = 4
	void crear(double c3, double c2, double c1, double c0) {gr=3; c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3;}
	void crear(double c2, double c1, double c0) {gr=2; c[0]=c0; c[1]=c1; c[2]=c2;}
	void crear(double c1, double c0) {gr=1; c[0]=c0; c[1]=c1;}

	double evaluate(double variable) {double ret = 0; for (int i=gr; i>=0; i--) {ret = ret * variable + c[i];} return ret;}

#ifdef DEBUG_STURM
	bool calculaRaicesSturm(int &nraices, double *raices) {L_ptr<double> u(11); return Sturm::main_sturm(gr, c.elem, nraices, u);}
#else
	bool calculaRaicesSturm(int &nraices, double *raices) {return Sturm::main_sturm(gr, &(c[0]), nraices, raices);}
#endif

	double pruebaSturm(); // Devuelve el error mean al evaluar el polinomio en las raices obtenidas con sturm
};

class L_PolinGr11_matriz
{
public:
#ifdef L_BRACKET_LEVEL_DEBUG
	L_ptrM<L_PolinGr11> elem;
	L_PolinGr11_matriz() : elem(li,lj) {elem=NULL; li=0; lj=0;}
	L_PolinGr11_matriz(const L_PolinGr11_matriz &other) : elem(li,lj) {elem = NULL; reallocate(other.li, other.lj); int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) elem[i][j]=other.elem[i][j];}
#else
	L_PolinGr11 **elem;
	L_PolinGr11_matriz() {elem=NULL; li=0; lj=0;}
	L_PolinGr11_matriz(const L_PolinGr11_matriz &other) {elem = NULL; reallocate(other.li, other.lj); int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) elem[i][j]=other.elem[i][j];}
#endif
	int li, lj;
	L_PolinGr11_matriz &operator=(const L_PolinGr11_matriz &other) {reallocate(other.li, other.lj); int i, j; for (i=0; i<li; i++) for (j=0; j<lj; j++) elem[i][j]=other.elem[i][j]; return *this;}
	void reallocate(int li, int lj) {this->li=li; this->lj=lj; if (elem!=NULL) L_delete2d<L_PolinGr11>((L_PolinGr11 **)elem); elem=L_new2d<L_PolinGr11>(li, lj);}

	L_PolinGr11 &operator()(int i, int j) {return elem[i][j];}
	const L_PolinGr11 &operator()(int i, int j) const {return elem[i][j];}

	L_PolinGr11 det();
	void subMatrElimFilCol(L_PolinGr11_matriz &ret, int i0, int j0);
	L_PolinGr11_matriz subMatrix(int i1, int i2, int j1, int j2) {int i, j; L_PolinGr11_matriz r; r.reallocate(i2-i1+1, j2-j1+1); for (i=i1; i<=i2; i++) for (j=j1; j<=j2; j++) r.elem[i-i1][j-j1]=elem[i][j]; return r;} // es como escribir m(i1:i2,j1:j2) en matlab pero con indexes desde zero
	void evaluaEn(L_Matrix &m, double variable);

	void print(FILE *fp, const char *name) const;
	void printMatlabFormat(FILE *fp, const char *name) const;
	void printMatlabFormat(const char *arch, const char *name) const {FILE *fp=fopen(arch,"w"); if (fp!=NULL) {printMatlabFormat(fp, name); fclose(fp);}}

	void destroy() {if (elem!=NULL) {L_delete2d<L_PolinGr11>((L_PolinGr11 **)elem);} li=0; lj=0; elem=NULL; }
	~L_PolinGr11_matriz() {destroy();}
};

#define L_POLIN_XYZ_NCOL 29
class L_Polin_xyz // representa una sum de monomios k * x^a*y^b*z^c con a<=3, b<=3, c<=4
{
public:
	std::vector<double> c;
	L_Polin_xyz() : c(L_POLIN_XYZ_NCOL) {if (!inicializado) inicializa(); for (int i=0; i<L_POLIN_XYZ_NCOL; i++) c[i]=0;}
	L_Polin_xyz(const L_Polin_xyz &other) : c(L_POLIN_XYZ_NCOL) {if (!inicializado) inicializa(); for (int i=0; i<L_POLIN_XYZ_NCOL; i++) c[i]=other.c[i];}
	L_Polin_xyz& operator=(const L_Polin_xyz &other) {if (!inicializado) inicializa(); for (int i=0; i<L_POLIN_XYZ_NCOL; i++) c[i]=other.c[i]; return *this;}

	void sacaFila(const L_Matrix &m, int i) {int j; for (j=0; j<m.lj; j++) c[j]=m(i,j); for(;j<L_POLIN_XYZ_NCOL; j++) c[j]=0;}
	void multX(L_Polin_xyz &fil);
	void multY(L_Polin_xyz &fil);
	void multZ(L_Polin_xyz &fil);
	void factoriza_xy_x_y_1(L_PolinGr11 &p1, L_PolinGr11 &p2, L_PolinGr11 &p3, L_PolinGr11 &p4, int exp1, int exp2, int exp3, int exp4);
	void factorizaInv_xy_x_y_1(L_PolinGr11 &p1, L_PolinGr11 &p2, L_PolinGr11 &p3, L_PolinGr11 &p4, int exp1, int exp2, int exp3, int exp4);

	static void inicializa();
	static int col(int expx, int expy, int expz) {return grInd[expx + expy*4 + expz*16];}
	static void exps(int col, int &expx, int &expy, int &expz) {int w=grIndInv[col]; expx=w%4; expy=(w/4)%4; expz=w/16;}

	void imprimeMonomio(double c, int e1, int e2, int e3) {printf("%.3g", c); if (e1>0)	printf("*x^%d", e1); if (e2>0) printf("*y^%d", e2); if (e3>0) printf("*z^%d", e3);}
	void print() {int u,v,w; for (int i=0; i<L_POLIN_XYZ_NCOL; i++) {if (c[i]!=0) {u=i%4; v=(i/4)%4; w=i/16; imprimeMonomio(c[i], u, v, w); printf(" + ");}} printf(".\n");}


	friend L_Polin_xyz operator+(const L_Polin_xyz &p1, const L_Polin_xyz &p2) {L_Polin_xyz r; for (int i=0; i<L_POLIN_XYZ_NCOL; i++) r.c[i]=p1.c[i]+p2.c[i]; return r;}
	friend L_Polin_xyz operator-(const L_Polin_xyz &p1, const L_Polin_xyz &p2) {L_Polin_xyz r; for (int i=0; i<L_POLIN_XYZ_NCOL; i++) r.c[i]=p1.c[i]-p2.c[i]; return r;}
	friend L_Polin_xyz operator*(const L_Polin_xyz &p, double d) {L_Polin_xyz r; for (int i=0; i<L_POLIN_XYZ_NCOL; i++) r.c[i]=p.c[i]*d; return r;}
	friend L_Polin_xyz operator*(double d, const L_Polin_xyz &p) {L_Polin_xyz r; for (int i=0; i<L_POLIN_XYZ_NCOL; i++) r.c[i]=p.c[i]*d; return r;}

private:
	static bool inicializado;
	static std::vector<int> grInd; // LUT para saber la relacion entre numero de columna y tripleta de exponentes
	static std::vector<int> grIndInv; // LUT para saber la relacion entre numero de columna y tripleta de exponentes
};

// Esto funciona bien, pero no es facil aplicarlo
// Usado por la class_name L_RelEpipolarEsencial
class L_Algoritmo5Puntos
{
public:
	static bool algoritmo5Puntos(L_Matrix E[10], int &nSols, const L_Matrix &u1v1_u2v2); // Convencion de ejes der,aba,ade
	static bool algoritmo5_mas_1_puntos(L_Matrix &E, const L_Matrix &u1v1_u2v2); // Convencion de ejes der,aba,ade


	static bool ejemploAlgoritmo5Puntos();
	static bool ejemploAlgoritmo5Puntos_xAde();  // Incompleto
	static bool ejemploAlgoritmo6Puntos(); // 5+1 puntos
	static bool ejemploAlgoritmo6Puntos_xAde(); // 5+1 puntos   Incompleto
private:

	static bool despejar_XYZW_5p(L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W, const L_Matrix &QagrT); // QagrT de 5x9
	static bool despejar_XYZW_np(L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W, const L_Matrix &QagrT); // QagrT de nx9, n>5
	static bool calcMatrizCubicas9x20(L_Matrix &Rc, L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W); // ESTA BUENA
	static bool calcMatrizCubicas10x20(L_Matrix &Rc, L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W); // ESTA BUENA
	static bool calcMatricesBC(L_PolinGr11_matriz &B, L_PolinGr11_matriz &C, L_Matrix &Rc);
	static bool calcPolinGr10(L_PolinGr11 &p, L_PolinGr11_matriz &B, L_PolinGr11_matriz &C);
};



// Tipo de datos L_ImagenBNxfloatFnConv : "puntero a fn de convolucion en L_ImageGrayDouble"
typedef L_ImageGrayDouble& (L_ImageGrayDouble::*L_ImagenBNxfloatFnConv)(const L_ImageGrayDouble &im1, int paso);


class L_ImagenBNxfloatFnConvNodo
{
public:
	L_ImagenBNxfloatFnConv fX;
	L_ImagenBNxfloatFnConv fY;
	long sigN;
	static int cmp(const void *nodo1, const void *nodo2)
	{
		return ( (*(L_ImagenBNxfloatFnConvNodo *)nodo1).sigN>(*(L_ImagenBNxfloatFnConvNodo *)nodo2).sigN )
			-
			( (*(L_ImagenBNxfloatFnConvNodo *)nodo1).sigN<(*(L_ImagenBNxfloatFnConvNodo *)nodo2).sigN );
	}
	static long calcN(double sigma) {return (long)(sigma*10000);}
};

class L_ImagenBNxfloatFnConvArreglo
{
public:
	L_ImagenBNxfloatFnConvNodo *arr;
	int n;
	int nMem;
	bool activo;
	L_ImagenBNxfloatFnConvArreglo(bool activar)
	{
		arr=NULL;
		n=0;
		nMem=0;
		activo=activar;

		#if defined L_defAgregaFnBNxfloat
		#error Este name no debe estar definido externamente
		#endif

		// Agregar aqui todos los sigma que tienen funcion precalculada
		#if defined(L_defAgregaFnBNxfloat)
			#error #define conflicts
		#endif
		#define L_defAgregaFnBNxfloat(S) ( addFrom(&L_ImageGrayDouble::z_difumX_##S, &L_ImageGrayDouble::z_difumY_##S, S) )
		//
		// Para scale space, s0=0.5, s1=1.0
		#ifdef _s0_0_5__nio_2
		L_defAgregaFnBNxfloat(8660); // Solo para sDoG
		L_defAgregaFnBNxfloat(10000);
		L_defAgregaFnBNxfloat(14142);
		L_defAgregaFnBNxfloat(20000);
		L_defAgregaFnBNxfloat(28284);
		L_defAgregaFnBNxfloat(30000); // Solo para Harris Laplace
		L_defAgregaFnBNxfloat(42426); // Solo para Harris Laplace
		#endif
		#ifdef _s0_1_6__nio_2
		L_defAgregaFnBNxfloat(12489);
		L_defAgregaFnBNxfloat(16000);
		L_defAgregaFnBNxfloat(22627);
		L_defAgregaFnBNxfloat(32000);
		L_defAgregaFnBNxfloat(45254);
		L_defAgregaFnBNxfloat(48000); // Solo para Harris Laplace
		L_defAgregaFnBNxfloat(67882); // Solo para HArris Laplace
		#endif
		//
		#undef L_defAgregaFnBNxfloat

		sort();
	};
	void addFrom(L_ImagenBNxfloatFnConv fX, L_ImagenBNxfloatFnConv fY, int sigN)
	{
		if (n==nMem)
		{
			L_ImagenBNxfloatFnConvNodo *arr2;
			arr2=(L_ImagenBNxfloatFnConvNodo *)realloc(arr, (nMem+10) * sizeof(L_ImagenBNxfloatFnConvNodo));
			if (arr2==NULL)
				throw std::bad_alloc();
			else
				arr=arr2;
			nMem+=10;
		}
		arr[n].fX=fX;
		arr[n].fY=fY;
		arr[n].sigN=sigN;
		n++;
	}
	void sort()
	{
		qsort(arr,n,sizeof(L_ImagenBNxfloatFnConvNodo),L_ImagenBNxfloatFnConvNodo::cmp);
	}
	L_ImagenBNxfloatFnConvNodo *buscaFn(int sigN) const
	{
		L_ImagenBNxfloatFnConvNodo f0;
		if (!activo || n==0)
			return NULL;
		f0.fX=(L_ImagenBNxfloatFnConv)NULL;
		f0.fY=(L_ImagenBNxfloatFnConv)NULL;
		f0.sigN=sigN;
		return (L_ImagenBNxfloatFnConvNodo *)bsearch(&f0,arr,n,sizeof(L_ImagenBNxfloatFnConvNodo),L_ImagenBNxfloatFnConvNodo::cmp);
	}
	~L_ImagenBNxfloatFnConvArreglo()
	{
		if (arr!=NULL)
			free(arr);
	}
};

#define L_BLOB_BORDE 10
// Segmentador colores
class L_BlobColor
{
public:
	// Estadisticas y caracteristicas del blob
	int indiceColor;
	L_Statistics2D pos;
	L_Statistics1D colorR;
	L_Statistics1D colorG;
	L_Statistics1D colorB;
	int area;
	bool saleIzq;
	bool saleDer;
	bool saleArr;
	bool saleAba;
	int keyp[4][2]; // right, down, left, up
	L_BlobColor() : indiceColor(-1), pos(), colorR(), colorG(), colorB(), area(0), saleIzq(false), saleDer(false), saleArr(false), saleAba(false) {reseteaKeypoint();}
	void resetea() {indiceColor=-1; pos.clear(); colorR.clear(); colorG.clear(); colorB.clear(); area=0; saleIzq=false; saleDer=false; saleArr=false; saleAba=false; reseteaKeypoint();}
	// Funciones para detectar cuando el blob sale de la imagen
	void revSaleIm(int i, int j, int lx, int ly) {saleIzq|=i<L_BLOB_BORDE; saleDer|=i>lx-L_BLOB_BORDE-1; saleArr|=j<L_BLOB_BORDE; saleAba|=j>ly-L_BLOB_BORDE-1;}
	void limpiaSaleIm() {saleIzq=false; saleDer=false; saleArr=false; saleAba=false;}
	void calcSaleImUsaVar(int lx, int ly, bool modeloCuadrado0_circulo1);
	bool saleIm() {return saleIzq||saleDer||saleArr||saleAba;}
	// Estos agregaPixel no revisan que cada punto se salga de la pantalla (por rapidez)
	inline void agregaPixel(L_ImageRGBUchar &im, int c, int f) {pos.push_pair(c,f); colorR.push(im.pix(c,f,0)); colorG.push(im.pix(c,f,1)); colorB.push(im.pix(c,f,2)); area++;}
	inline void agregaPixel(L_uchar R, L_uchar G, L_uchar B, int i, int j) {pos.push_pair(i,j); colorR.push(R); colorG.push(G); colorB.push(B); area++;}
	inline void agregaPixel(int c, int f) {pos.push_pair(c,f); area++;}
	inline void reseteaKeypoint() {keyp[0][0] = -1; keyp[1][1] = -1; keyp[2][0] = 32767; keyp[3][1] = 32767;}
	inline void agregaKeypoint(int c, int f) {if (c > keyp[0][0]) {keyp[0][0] = c; keyp[0][1] = f;} if (f > keyp[1][1]) {keyp[1][0] = c; keyp[1][1] = f;} if (c < keyp[2][0]) {keyp[2][0] = c; keyp[2][1] = f;} if (f < keyp[3][1]) {keyp[3][0] = c; keyp[3][1] = f;}}
	// Agregar funciones tipo agregaTira(...) para RLE futuro
};

class L_SegmentadorAbstracto
{
public:
	std::vector<L_BlobColor> blobs;
	L_ImageGrayUchar imSegmentada;
	typedef int size_type;

	L_SegmentadorAbstracto():blobs(0) {blobs.reserve(500);}
	virtual bool segmentaImagen(const L_ImageRGBUchar &original)=0; // blobs debe limpiarse manualmente antes de llamar a esta funcion
	virtual bool pintaImagen(L_ImageRGBUchar &output)=0; // blobs debe limpiarse manualmente antes de llamar a esta funcion

	int buscaMayorBlob(int indiceColor) {size_type i, iMax=-1, area=-1; for (i=0; i<(size_type)blobs.size(); i++) {if (blobs[i].indiceColor==indiceColor && blobs[i].area>area) {area=blobs[i].area; iMax = i;}} return iMax;}
	virtual ~L_SegmentadorAbstracto() {}
};

class L_ColorIndice
{
public:
	double colorCentralR;
	double colorCentralG;
	double colorCentralB;
	double n;
	L_ColorIndice() : n(0), colorCentralR(0), colorCentralG(0), colorCentralB(0) {}
	void addFrom(L_uchar R, L_uchar G, L_uchar B, double peso=1) {colorCentralR=(n*colorCentralR+peso*R)/(n+peso);colorCentralG=(n*colorCentralG+peso*G)/(n+peso);colorCentralB=(n*colorCentralB+peso*B)/(n+peso);n+=peso;}
	void resetea() {n=0;}
};

class L_SegmentadorColoresLUTdif : public L_SegmentadorAbstracto
{
public:
	char ***LUT;  //  segm[i][j] = LUT(R/dR, g/dG, b/dB);
	int ****hist;
	long nPixEjemplos;
	int nIndices;
	int dR, dG, dB; // 32 tipicamente
	L_ImageGrayUchar recursion;
	L_ColorIndice *colorIndice;
	bool usarRLE;
	int areaMin;
	// Funciones no implementadas
private:
	L_SegmentadorColoresLUTdif(); // {L_hard_shutdown("Llamada a constructor no valido");}
	L_SegmentadorColoresLUTdif(const L_SegmentadorColoresLUTdif &other); // {L_hard_shutdown("Llamada a constructor no valido");}
	// Funciones implementadas
	void pideMemoHist();
public:
	L_SegmentadorColoresLUTdif(int numR, int numG, int numB, int nIndicesTot) : dR(256/numR), dG(256/numG), dB(256/numB), nIndices(nIndicesTot), nPixEjemplos(0), usarRLE(false), areaMin(16) {LUT = L_new3d<char>(numR, numG, numB); memset(LUT[0][0],0,numR*numG*numB); hist=NULL; colorIndice=new L_ColorIndice[nIndices];}
	void generaImagenSegmentada(const L_ImageRGBUchar &im) {imSegmentada.reallocate(im.lx, im.ly); for (int j=0; j<im.ly; j++) for (int i=0; i<im.lx; i++) imSegmentada.pix(i,j)=LUT[im.pix(i,j,0)/dR][im.pix(i,j,1)/dG][im.pix(i,j,2)/dB];}
	void entrenaHistograma(const L_ImageRGBUchar &im, const L_ImageGrayUchar &indexes); // value 255 en indexes = saltar pixel
	void entrenaHistograma(const L_ImageRGBUchar &im, const L_ShapeArray &polig);
	void fusionaHistograma(const int ***histExt, int nPixEjemplosExt);
	void borraHistograma(bool borrarLUT=true);
	bool generaLUT_mayorVotacion();
	bool segmentaImagen(const L_ImageRGBUchar &original);
	bool pintaImagen(L_ImageRGBUchar &output);

	bool segmentaDif(const L_ImageRGBUchar &original);
	bool segmentaRLE(const L_ImageRGBUchar &original);
	bool segmentaDif(const L_ImageGrayUchar &original);  // No usa LUT

	void grabarHist(FILE *fp);
	void leerHist(FILE *fp);

	void reallocate(int numR, int numG, int numB, int nIndicesTot);

	~L_SegmentadorColoresLUTdif() {L_delete3d<char>(LUT); if (hist!=NULL) {for (int i=0; i<nIndices; i++) L_delete3d<int>(hist[i]); delete[] hist;}; delete[] colorIndice;}
};

class L_Segmentador2ColoresAdaptivo : public L_SegmentadorAbstracto
{
private:
	double facR, facG, facB, fac1;
	bool inicializado;
	L_ImageGrayUchar recursion;
	int areaMin;
	int histogramaProy[600];

	L_Segmentador2ColoresAdaptivo(); // {L_hard_shutdown("Llamada a constructor no valido");}
	L_Segmentador2ColoresAdaptivo(const L_SegmentadorColoresLUTdif &other); // {L_hard_shutdown("Llamada a constructor no valido");}
	// Funciones implementadas
private:
	void init(L_uchar r1, L_uchar g1, L_uchar b1, L_uchar r0, L_uchar g0, L_uchar b0);
	void adaptar(const L_ImageRGBUchar &original);
	bool segmentaDif(const L_ImageRGBUchar &original);
public:
	L_Segmentador2ColoresAdaptivo(L_uchar r1, L_uchar g1, L_uchar b1, L_uchar r0, L_uchar g0, L_uchar b0) : areaMin(16) { init(r1, g1, b1, r0, g0, b0);}
	bool segmentaImagen(const L_ImageRGBUchar &original) { adaptar(original); return segmentaDif(original);}
	bool pintaImagen(L_ImageRGBUchar &output);
};




class L_CaracterizadorBlobs
{
public:
	static bool ajustaCirculoAlgebraico(double *x, double *y, bool *acept, int n, double &r, double &cx, double &cy);
	static bool ajustaCirculoRansac(double *x, double *y, int n, int nIntentos, bool *acept, int nMin, double maxErr, double &r, double &cx, double &cy);
	static bool scanLineFromImage(L_ImageGrayUchar &imSegm, L_uchar c, double x, double y, double dx, double dy, int toler, int &xf, int &yf);
	static bool calcCirculo(L_ImageGrayUchar &imSegm, L_uchar c, int x, int y, double maxErr, double &r, double &cx, double &cy);
};


/////
// Objetos geometricos y modelos de camaras
/////

enum L_SistemaEjes3D {L_AdeIzqArr, L_DerAbaAde, L_DerArrAtr, L_IzqArrAde};

class L_CoordsCart2D
{
public:
	double x;
	double y;
	L_CoordsCart2D() {}
	L_CoordsCart2D(double theta) {x=cos(theta); y=sin(theta);}
	L_CoordsCart2D(double xVal, double yVal) : x(xVal), y(yVal) {}
	L_CoordsCart2D &fijaXY(double xVal, double yVal) {x=xVal; y=yVal; return *this;}
	L_CoordsCart2D &fijaRTheta(double r, double theta) {x=r*cos(theta); y=r*sin(theta); return *this;}
	L_CoordsCart2D operator - () {return L_CoordsCart2D(-x, -y);}
	L_CoordsCart2D & operator+=(const L_CoordsCart2D &other) {x+=other.x; y+=other.y; return *this;}
	L_CoordsCart2D & operator-=(const L_CoordsCart2D &other) {x+=other.x; y+=other.y; return *this;}
	L_CoordsCart2D operator+(const L_CoordsCart2D &other) const {return L_CoordsCart2D(x+other.x, y+other.y);}
	L_CoordsCart2D operator-(const L_CoordsCart2D &other) const {return L_CoordsCart2D(x-other.x, y-other.y);}
	double punto (const L_CoordsCart2D &other) const {return x*other.x+y*other.y;}
	double cruz (const L_CoordsCart2D &other) const {return x*other.y;}
	L_CoordsCart2D & multiply_to_each_element(double factor) {x*=factor; y*=factor; return *this;}
	L_CoordsCart2D prodComplejo (const L_CoordsCart2D &other) const {return L_CoordsCart2D(x*other.x-y*other.y, x*other.y+y*other.x);}
	L_CoordsCart2D &rota(double ang) {double c=cos(ang),s=sin(ang), x0=x; x=x*c+y*s; y=-x0*s+y*c; return *this;}
	double pideR () const {return sqrt(x*x+y*y);}
	double pideTheta () const {return atan2(y,x);}
};

// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_CoordsCart3D_fijaXYZ(pun,X,Y,Z) ((pun).x=(X), (pun).y=(Y), (pun).z=(Z))
#define L_CoordsCart3D_fijaCero(pun) ((pun).x=(pun).y=(pun).z=0.0)
#define L_CoordsCart3D_OP_igual(pret,p) ((pret).x=(p).x,(pret).y=(p).y,(pret).z=(p).z)
#define L_CoordsCart3D_OP_suma(pret,p1,p2) ((pret).x=(p1).x+(p2).x, (pret).y=(p1).y+(p2).y, (pret).z=(p1).z+(p2).z)
#define L_CoordsCart3D_OP_resta(pret,p1,p2) ((pret).x=(p1).x-(p2).x, (pret).y=(p1).y-(p2).y, (pret).z=(p1).z-(p2).z)
#define L_CoordsCart3D_punto(p1,p2) ((p1).x*(p2).x + (p1).y*(p2).y + (p1).z*(p2).z)
#define L_CoordsCart3D_cruz(pret,p1,p2) ((pret).x = (p1).y*(p2).z - (p1).z*(p2).y, (pret).y = (p1).z*(p2).x - (p1).x*(p2).z, (pret).z = (p1).x*(p2).y - (p1).y*(p2).x)
#define L_CoordsCart3D_distanciaA(p1,p2) (sqrt(((p1).x-(p2).x)*((p1).x-(p2).x)+((p1).y-(p2).y)*((p1).y-(p2).y)+((p1).z-(p2).z)*((p1).z-(p2).z)))
#define L_CoordsCart3D_distancia2A(p1,p2) ((((p1).x-(p2).x)*((p1).x-(p2).x)+((p1).y-(p2).y)*((p1).y-(p2).y)+((p1).z-(p2).z)*((p1).z-(p2).z)))
#define L_CoordsCart3D_normaliza(pun) {double d=sqrt((pun).x*(pun).x+(pun).y*(pun).y+(pun).z*(pun).z); if (d!=0) {(pun).x/=d; (pun).y/=d; (pun).z/=d;} else {pun.x=0; pun.y=0; pun.z=0;}}

class L_CoordsCart3D
{
public:
	double x; // En general es el eje que va hacia adelante
	double y; // En general es el eje que va hacia la izquierda
	double z; // En general es el eje que va hacia arriba
	L_CoordsCart3D() {}
	L_CoordsCart3D(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}
	double &el(int i) {switch(i) {case 0:return x; case 1:return y; case 2:return z; default: throw_L_ArgException_if(true,"L_CoordsCart3D::el"); return x;}}
	const double &el(int i) const {switch(i) {case 0:return x; case 1:return y; case 2:return z; default: throw_L_ArgException_if(true,"L_CoordsCart3D::el"); return x;}}
	void define(double xVal, double yVal, double zVal) {x=xVal; y=yVal; z=zVal;}
	double *data() {return &x;}
	const double *data() const {return &x;}
	void swap(L_CoordsCart3D &other) {L_CoordsCart3D t = other; other = *this; *this = t;}
	L_CoordsCart3D &fijaXYZ(double xVal, double yVal, double zVal) {x=xVal; y=yVal; z=zVal; return *this;}
	L_CoordsCart3D &fijaRThetaPhi(double r, double theta, double phi) {double ct=cos(theta), st=sin(theta), cf=cos(phi), sf=sin(phi); x=r*ct*cf; y=r*st*cf; z=r*sf; return *this;}
	L_CoordsCart3D operator - () {return L_CoordsCart3D(-x, -y, -z);}
	L_CoordsCart3D & operator+=(const L_CoordsCart3D &other) {x+=other.x; y+=other.y; z+=other.z; return *this;}
	L_CoordsCart3D & operator-=(const L_CoordsCart3D &other) {x-=other.x; y-=other.y; z-=other.z; return *this;} // ESTABA MALO!!
	double punto (const L_CoordsCart3D &other) const {return x*other.x+y*other.y+z*other.z;}
	L_CoordsCart3D cruz (const L_CoordsCart3D &other) const {return L_CoordsCart3D(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);}
	L_CoordsCart3D &operator*=(double factor) {x*=factor; y*=factor; z*=factor; return *this;}
	L_CoordsCart3D &operator/=(double divisor) {x/=divisor; y/=divisor; z/=divisor; return *this;}
	L_CoordsCart3D &rotaEnEjeX(double ang) {double c=cos(ang),s=sin(ang), y0=y; y=y*c+z*s; z=-y0*s+z*c; return *this;}
	L_CoordsCart3D &rotaEnEjeY(double ang) {double c=cos(ang),s=sin(ang), x0=x; x=x*c+z*s; z=-x0*s+z*c; return *this;}
	L_CoordsCart3D &rotaEnEjeZ(double ang) {double c=cos(ang),s=sin(ang), x0=x; x=x*c+y*s; y=-x0*s+y*c; return *this;}

	L_CoordsCart3D rotarVectRodrigues(L_CoordsCart3D &p)
	{
		L_CoordsCart3D n;
		double ang, c, s;
		n=*this;
		ang=n.pideR();
		n/=ang;
		c=cos(ang);
		s=sin(ang);
		return p*c +
			n.cruz(p)*s +
			n*(n.punto(p))*(1-c);
	}
	L_CoordsCart3D componerVectorRot(L_CoordsCart3D &w2);

	void cambiaEjes_a_AdeIzqArr() {double X,Y,Z; X=z; Y=-x; Z=-y;  x=X; y=Y; z=Z;} // Convencion de ejes en robotica
	void cambiaEjes_a_DerAbaAde() {double X,Y,Z; X=-y; Y=-z; Z=x;  x=X; y=Y; z=Z;} // Convencion de ejes en vision
	void cambiaEjes(const L_CoordsCart3D &orig, L_SistemaEjes3D inicial, L_SistemaEjes3D final);

	double pideR () const {return sqrt(x*x+y*y+z*z);}
	double pideR2 () const {return x*x+y*y+z*z;}
	double pideTheta () const {return atan2(y,x);}
	double pidePhi () const {return atan2(z, sqrt(x*x+y*y));}
	double distanciaA(const L_CoordsCart3D &other) const {return sqrt((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z));}
	double distancia2A(const L_CoordsCart3D &other) const {return ((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z));}
	void normaliza() {double d=sqrt(x*x+y*y+z*z); if (d!=0) {x/=d; y/=d; z/=d;} else {printf("Problema en L_CoordsCart3D::normaliza()\n"); x=1; y=0; z=0;}}
	bool normaliza_ret();
	L_CoordsCart3D unitario() const {double d = sqrt(x*x+y*y+z*z); if (d==0) {printf("Problema en L_CoordsCart3D::unitario()\n"); return L_CoordsCart3D(1,0, 0);}; return L_CoordsCart3D(x/d, y/d, z/d);}
	double cosineDistanceTo(const L_CoordsCart3D &other) const {return (x*other.x+y*other.y+z*other.z)/sqrt((x*x+y*y+z*z)*(other.x*other.x+other.y*other.y+other.z*other.z));}
	void fijaCero() {x=y=z=0;}
	void fijaAzar(double dist = 10) {x=L_RANDOM(-dist,dist); y=L_RANDOM(-dist,dist); z=L_RANDOM(-dist,dist);}
	void fijaAzarSemilla(double dist, L_Rand &randF) {x=randF.random(-dist,dist); y=randF.random(-dist,dist); z=randF.random(-dist,dist);}
	void print(const char *name, FILE *fp=stdout) const {fprintf(fp,"%s = (%.2f %.2f %.2f)\n", name, x, y, z);}
	char *print(const char *name, char *buffer) const {sprintf(buffer,"%s = (%.2f %.2f %.2f)", name, x, y, z); return buffer;}

	L_CoordsCart3D operator+(const L_CoordsCart3D &b) const {return L_CoordsCart3D(x+b.x, y+b.y, z+b.z);}
	L_CoordsCart3D operator-(const L_CoordsCart3D &b) const {return L_CoordsCart3D(x-b.x, y-b.y, z-b.z);}
	L_CoordsCart3D operator*(double factor) const {return L_CoordsCart3D(x*factor,y*factor,z*factor);}
	L_CoordsCart3D operator/(double divisor) const {return L_CoordsCart3D(x/divisor,y/divisor,z/divisor);}

	static void test();
};

inline L_CoordsCart3D operator*(double factor, const L_CoordsCart3D &a) {return L_CoordsCart3D(a.x*factor,a.y*factor,a.z*factor);}


class L_CoordsPolar2D
{
public:
	double r;
	double theta;
	L_CoordsPolar2D() {}
	L_CoordsPolar2D(double rVal, double thetaVal) : r(rVal),theta(thetaVal) {}
	L_CoordsPolar2D &fijaXY(double xVal, double yVal) {r=sqrt(xVal*xVal+yVal*yVal); theta=atan2(yVal,xVal); return *this;}
	L_CoordsPolar2D &fijaRTheta(double rVal, double thetaVal) {r=rVal; theta=thetaVal; return *this;}
	L_CoordsPolar2D &fijaCoordsCart2D(const L_CoordsCart2D &other) {r=sqrt(other.x*other.x+other.y*other.y); theta=atan2(other.y,other.x); return *this;}
	L_CoordsPolar2D operator- () {return L_CoordsPolar2D(r,theta+M_PI);}
	L_CoordsPolar2D & operator+=(const L_CoordsPolar2D &other) {L_CoordsCart2D sum(r*cos(theta)+other.r*cos(other.theta), r*sin(theta)+other.r*cos(other.theta)); r=sum.pideR(); theta=sum.pideTheta(); return *this;}
	L_CoordsPolar2D & operator-=(const L_CoordsPolar2D &other) {L_CoordsCart2D sum(r*cos(theta)-other.r*cos(other.theta), r*sin(theta)-other.r*cos(other.theta)); r=sum.pideR(); theta=sum.pideTheta(); return *this;}
	L_CoordsPolar2D operator+(const L_CoordsPolar2D &other) const {L_CoordsCart2D sum(r*cos(theta)+other.r*cos(other.theta), r*sin(theta)+other.r*cos(other.theta)); return L_CoordsPolar2D(sum.pideR(), sum.pideTheta());}
	L_CoordsPolar2D operator-(const L_CoordsPolar2D &other) const {L_CoordsCart2D sum(r*cos(theta)-other.r*cos(other.theta), r*sin(theta)-other.r*cos(other.theta)); return L_CoordsPolar2D(sum.pideR(), sum.pideTheta());}
	L_CoordsPolar2D & multiply_to_each_element(double factor) {r*=factor; return *this;}
	double punto (const L_CoordsPolar2D &other) const {return r*other.r*cos(theta-other.theta);}
	double cruz (const L_CoordsPolar2D &other) const {return r*other.r*sin(theta-other.theta);}
	L_CoordsPolar2D prodComplejo (const L_CoordsPolar2D &other) const {return L_CoordsPolar2D(r*other.r, theta+other.theta);}
	L_CoordsPolar2D &rota(double ang) {r+=ang; return *this;}
	double pideX () const {return r*cos(theta);}
	double pideY () const {return r*sin(theta);}
};

class L_CoordsEsf3D
{
public:
	double r;
	double theta; // Es como el paralelo (en radianes)
	double phi;   // Es como el meridiano (en radianes)
	L_CoordsEsf3D() {}
	L_CoordsEsf3D(double rVal, double thetaVal, double phiVal) : r(rVal),theta(thetaVal), phi(phiVal) {}
	// MALA L_CoordsEsf3D &fijaXYZ(double xVal, double yVal, double zVal) {r=sqrt(xVal*xVal+yVal*yVal+zVal*zVal); theta=atan2(yVal,xVal); phi=atan2(zVal,r); return *this;}
	L_CoordsEsf3D &fijaXYZ(double xVal, double yVal, double zVal) {r=sqrt(xVal*xVal+yVal*yVal+zVal*zVal); theta=atan2(yVal,xVal); phi=atan2(zVal,sqrt(xVal*xVal+yVal*yVal)); return *this;}
	L_CoordsEsf3D &fijaRThetaPhi(double rVal, double thetaVal, double phiVal) {r=rVal; theta=thetaVal; phi=phiVal; return *this;}
	L_CoordsEsf3D &fijaCoordsCart3D(const L_CoordsCart3D &other) {r=sqrt(other.x*other.x+other.y*other.y+other.z*other.z); theta=atan2(other.y,other.x); phi=atan2(other.z,sqrt(other.x*other.x+other.y*other.y)); return *this;}
	L_CoordsCart3D pideCoordsCart3D() {L_CoordsCart3D ret; ret.x=r*cos(theta)*cos(phi); ret.y=r*sin(theta)*cos(phi); ret.z=r*sin(phi); return ret;}
	L_CoordsEsf3D operator- () {return L_CoordsEsf3D(r,theta+M_PI, -phi);}
	L_CoordsEsf3D & operator+=(const L_CoordsEsf3D &other) {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(r,theta,phi); p2.fijaRThetaPhi(other.r, other.theta, other.phi); p1+=p2; fijaXYZ(p1.x, p1.y, p1.z); return *this;}
	L_CoordsEsf3D & operator-=(const L_CoordsEsf3D &other) {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(r,theta,phi); p2.fijaRThetaPhi(other.r, other.theta, other.phi); p1-=p2; fijaXYZ(p1.x, p1.y, p1.z); return *this;}
	L_CoordsEsf3D operator+(const L_CoordsEsf3D &other) const {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(r,theta,phi); p2.fijaRThetaPhi(other.r, other.theta, other.phi); p1+=p2; L_CoordsEsf3D a; return a.fijaXYZ(p1.x, p1.y, p1.z);}
	L_CoordsEsf3D operator-(const L_CoordsEsf3D &other) const {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(r,theta,phi); p2.fijaRThetaPhi(other.r, other.theta, other.phi); p1-=p2; L_CoordsEsf3D a; return a.fijaXYZ(p1.x, p1.y, p1.z);}
	L_CoordsEsf3D & multiply_to_each_element(double factor) {r*=factor; return *this;}
	double punto (const L_CoordsEsf3D &other) const {double cf1=cos(phi), cf2=cos(other.phi); return r*other.r*( cos(theta)*cf1*sin(other.theta)*cf2+sin(theta)*cf1*sin(other.theta)*cf2+sin(phi)*sin(other.phi) );}
	L_CoordsEsf3D cruz (const L_CoordsEsf3D &other) const {L_CoordsCart3D p1, p2; p1.fijaRThetaPhi(r,theta,phi); p2.fijaRThetaPhi(other.r,other.theta,other.phi); L_CoordsEsf3D ret; return ret.fijaCoordsCart3D(p1.cruz(p2)); }
	L_CoordsEsf3D &rotaEnEjeZ(double ang) {theta+=ang; return *this;}
	L_CoordsEsf3D &rotaEnEjeY(double ang) {L_CoordsCart3D p1; p1.fijaRThetaPhi(r,theta,phi); fijaCoordsCart3D(p1.rotaEnEjeY(ang)); return *this;}
	L_CoordsEsf3D &rotaEnEjeX(double ang) {L_CoordsCart3D p1; p1.fijaRThetaPhi(r,theta,phi); fijaCoordsCart3D(p1.rotaEnEjeX(ang)); return *this;}
	double pideX () const {return r*cos(theta)*cos(phi);}
	double pideY () const {return r*sin(theta)*cos(phi);}
	double pideZ () const {return r*sin(phi);}
};


// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_CoordsEsfInv3D_pideCoordsCart3D(ret,cei3d) ((ret).x=1/(cei3d).invR*cos((cei3d).theta)*cos((cei3d).phi), (ret).y=1/(cei3d).invR*sin((cei3d).theta)*cos((cei3d).phi), (ret).z=1/(cei3d).invR*sin((cei3d).phi))
#define L_CoordsEsfInv3D_jacob_pideCoordsCart3D(J3x3,cei3d,i0,j0) {double rho_cij = 1/(cei3d).invR; double rho2_cij = rho_cij*rho_cij; double cos_t_cij = cos((cei3d).theta); double sin_t_cij = sin((cei3d).theta); double cos_p_cij = cos((cei3d).phi); double sin_phi = sin((cei3d).phi); (J3x3)((i0)+0,(j0)+0) =  -cos_t_cij*cos_p_cij*rho2_cij; (J3x3)((i0)+0,(j0)+1) =  -sin_t_cij*cos_p_cij*rho_cij; (J3x3)((i0)+0,(j0)+2) =  -cos_t_cij*sin_phi*rho_cij; (J3x3)((i0)+1,(j0)+0) =  -sin_t_cij*cos_p_cij*rho2_cij; (J3x3)((i0)+1,(j0)+1) =   cos_t_cij*cos_p_cij*rho_cij; (J3x3)((i0)+1,(j0)+2) =  -sin_t_cij*sin_phi*rho_cij; (J3x3)((i0)+2,(j0)+0)  = -sin_phi*rho2_cij; (J3x3)((i0)+2,(j0)+1)  =  0; (J3x3)((i0)+2,(j0)+2)  =  cos_p_cij*rho_cij;}


// Para inverse depth parametrization
class L_CoordsEsfInv3D
{
public:
	double invR;
	double theta; // Es como el paralelo (en radianes)
	double phi;   // Es como el meridiano (en radianes)

	double &el(int i) {switch(i) {case 0: return invR; case 1: return theta; case 2: return phi; default: printf("L_CoordsEsfInv3D::el() fuera de rango\n");} return invR;}

	L_CoordsEsfInv3D() {}
	L_CoordsEsfInv3D(double invRVal, double thetaVal, double phiVal) : invR(invRVal),theta(thetaVal), phi(phiVal) {}
	L_CoordsEsfInv3D &fijaXYZ(double xVal, double yVal, double zVal) {invR=1/sqrt(xVal*xVal+yVal*yVal+zVal*zVal); theta=atan2(yVal,xVal); phi=atan2(zVal,sqrt(xVal*xVal+yVal*yVal)); return *this;}
	L_CoordsEsfInv3D &fijaRThetaPhi(double rVal, double thetaVal, double phiVal) {invR=1/rVal; theta=thetaVal; phi=phiVal; return *this;}
	L_CoordsEsfInv3D &fijaInvRThetaPhi(double invRval, double thetaVal, double phiVal) {invR=invRval; theta=thetaVal; phi=phiVal; return *this;}
	L_CoordsEsfInv3D &fijaCoordsCart3D(const L_CoordsCart3D &other) {double r=sqrt(other.x*other.x+other.y*other.y+other.z*other.z); invR=1/r; theta=atan2(other.y,other.x); phi=atan2(other.z,sqrt(other.x*other.x+other.y*other.y)); return *this;}
	void fijaAzar() {L_CoordsCart3D p; p.fijaAzar(); fijaCoordsCart3D(p);}
	L_CoordsCart3D pideCoordsCart3D() const {L_CoordsCart3D ret; ret.x=1/invR*cos(theta)*cos(phi); ret.y=1/invR*sin(theta)*cos(phi); ret.z=1/invR*sin(phi); return ret;}
	void jacob_pideCoordsCart3D(L_Matrix &J3x3) const; // deriv = dX3D(invR,theta,phi)/d(invR,theta,phi)
	L_CoordsEsfInv3D operator- () {return L_CoordsEsfInv3D(invR,theta+M_PI, -phi);}
	L_CoordsEsfInv3D & operator+=(const L_CoordsEsfInv3D &other) {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(1/invR,theta,phi); p2.fijaRThetaPhi(1/other.invR, other.theta, other.phi); p1+=p2; fijaXYZ(p1.x, p1.y, p1.z); return *this;}
	L_CoordsEsfInv3D & operator-=(const L_CoordsEsfInv3D &other) {L_CoordsCart3D p1,p2; p1.fijaRThetaPhi(1/invR,theta,phi); p2.fijaRThetaPhi(1/other.invR, other.theta, other.phi); p1-=p2; fijaXYZ(p1.x, p1.y, p1.z); return *this;}
	L_CoordsEsfInv3D & multiply_to_each_element(double factor) {invR/=factor; return *this;}
	double punto (const L_CoordsEsfInv3D &other) const {double cf1=cos(phi), cf2=cos(other.phi); return 1/(invR*other.invR)*( cos(theta)*cf1*sin(other.theta)*cf2+sin(theta)*cf1*sin(other.theta)*cf2+sin(phi)*sin(other.phi) );}
	L_CoordsEsfInv3D cruz (const L_CoordsEsfInv3D &other) const {L_CoordsCart3D p1, p2; p1.fijaRThetaPhi(1/invR,theta,phi); p2.fijaRThetaPhi(1/other.invR,other.theta,other.phi); L_CoordsEsfInv3D ret; return ret.fijaCoordsCart3D(p1.cruz(p2)); }
	L_CoordsEsfInv3D &rotaEnEjeZ(double ang) {theta+=ang; return *this;}
	L_CoordsEsfInv3D &rotaEnEjeY(double ang) {L_CoordsCart3D p1; p1.fijaRThetaPhi(1/invR,theta,phi); fijaCoordsCart3D(p1.rotaEnEjeY(ang)); return *this;}
	L_CoordsEsfInv3D &rotaEnEjeX(double ang) {L_CoordsCart3D p1; p1.fijaRThetaPhi(1/invR,theta,phi); fijaCoordsCart3D(p1.rotaEnEjeX(ang)); return *this;}
	double pideX () const {return 1/invR*cos(theta)*cos(phi);}
	double pideY () const {return 1/invR*sin(theta)*cos(phi);}
	double pideZ () const {return 1/invR*sin(phi);}

	friend L_CoordsEsfInv3D operator+(const L_CoordsEsfInv3D &a, const L_CoordsEsfInv3D &b);
	friend L_CoordsEsfInv3D operator-(const L_CoordsEsfInv3D &a, const L_CoordsEsfInv3D &b);

	double cosineDistanceTo(L_CoordsEsfInv3D &other) const {return (invR*other.invR + theta*other.theta + phi*other.phi) / sqrt((invR*invR + theta*theta + phi*phi)*(other.invR*other.invR + other.theta*other.theta + other.phi*other.phi));}

	void print(const char *name, FILE *fp = stdout) {printf("%s = (1/%f  <  %f <  %f)\n", name, invR, theta, phi);}

	static void test();
};




// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_CoordsInv3D_pideCoords3D(ret,puni) ((ret).x=(puni).pos.x+1/(puni).dir.invR*cos((puni).dir.theta)*cos((puni).dir.phi) , (ret).y = (puni).pos.y+1/(puni).dir.invR*sin((puni).dir.theta)*cos((puni).dir.phi), (ret).z = (puni).pos.z + 1/(puni).dir.invR*sin((puni).dir.phi))
#define L_CoordsInv3D_jacob_pideCoords3D(vi,J3x6,i0,j0) {L_StaticMatrix_identity3x3(J3x6); L_CoordsEsfInv3D_jacob_pideCoordsCart3D(J3x6,(vi).dir,i0,j0+3);}


class L_CoordsInv3D
{
public:
	L_CoordsCart3D pos;
	L_CoordsEsfInv3D dir;
	double &el(int i) {switch(i) {case 0: return pos.x; case 1: return pos.y; case 2: return pos.z; case 3: return dir.invR; case 4: return dir.theta; case 5: return dir.phi; default: throw_L_ArgException_if(true,"L_CoordsInv3D::el()"); return pos.x;}}
	void define(double x, double y, double z, double dx, double dy, double dz) {pos.define(x,y,z); dir.fijaXYZ(dx,dy,dz);}
	L_CoordsCart3D pideCoords3D() const {return pos + dir.pideCoordsCart3D();}
	void jacob_pideCoords3D(L_Matrix &J3x6, int i0=0, int j0=0) const;

	void fijaAzar() {pos.fijaAzar(); dir.fijaAzar();}

	void swap(L_CoordsInv3D &other) {L_CoordsInv3D t = other; other = *this; *this = t;}

	void print(const char *name, FILE *fp=stdout) const {fprintf(fp,"%s = (%.2f %.2f %.2f) (%.2f %.2f %.2f)\n", name, pos.x, pos.y, pos.z, dir.invR, dir.theta, dir.phi);}
	static void test();
};

#define L_RayoCamara_proyeccionDe(r,p) ((r).tanIzq = (p).y/(p).x, (r).tanArr = (p).z/(p).x)
#define L_RayoCamara_jacob_proyeccionDe(p,J2x3,i0,j0) ((J2x3)(i0+0,j0+0) = -(p).y/((p).x*(p).x), (J2x3)(i0+0,j0+1) = 1/(p).x, (J2x3)(i0+0,j0+2) = 0, (J2x3)(i0+1,j0+0) = -(p).z/((p).x*(p).x), (J2x3)(i0+1,j0+1) = 0, (J2x3)(i0+1,j0+2)=1/(p).x)

class L_RayoCamara // Convencion de ejes de Robotica
{
public:
	double tanIzq; // tangente de un rayo que sale de la camara proyectado sobre el eje horizontal de la camara
	double tanArr; // tangente de un rayo que sale de la camara proyectado sobre el eje vertical de la camara

	double &el(int i) {switch(i){case 0: return tanIzq;  case 1: return tanArr;  default: throw_L_ArgException_if(true,"L_RayoCamara::el"); return tanIzq;}}
	L_RayoCamara() {}
	L_RayoCamara(double tanIzqVal, double tanArrVal) : tanIzq(tanIzqVal), tanArr(tanArrVal) {}
	L_RayoCamara(double xRel, double yRel, double zRel) : tanIzq(yRel/xRel), tanArr(zRel/xRel) {}
	L_RayoCamara(const L_CoordsCart3D &p) : tanIzq(p.y/p.x), tanArr(p.z/p.x) {}
	void set(double tanIzqVal, double tanArrVal) {tanIzq=tanIzqVal; tanArr=tanArrVal;}
	void define(double tanIzqVal, double tanArrVal) {tanIzq=tanIzqVal; tanArr=tanArrVal;}
	double distTanA(const L_RayoCamara &other) const {return sqrt( (tanIzq-other.tanIzq)*(tanIzq-other.tanIzq) + (tanArr-other.tanArr)*(tanArr-other.tanArr) );}
	double dist2TanA(const L_RayoCamara &other) const {return (tanIzq-other.tanIzq)*(tanIzq-other.tanIzq) + (tanArr-other.tanArr)*(tanArr-other.tanArr);}
	L_CoordsCart3D pideUnitario() const {double d = sqrt(1 + tanIzq*tanIzq + tanArr*tanArr); return L_CoordsCart3D(1/d, tanIzq/d, tanArr/d);}
	double cosineDistanceTo(L_RayoCamara &other) {return (tanIzq*other.tanIzq + tanArr*other.tanArr) / sqrt((tanIzq*tanIzq+tanArr*tanArr)*(other.tanIzq*other.tanIzq+other.tanArr*other.tanArr));}
	void print(const char *name, FILE *fp = stdout) {printf("%s = (u=%.3f v=%.3f)\n", name, tanIzq, tanArr);}

	void proyeccionDe(const L_CoordsCart3D &p) {tanIzq = p.y/p.x; tanArr = p.z/p.x;}
	static void jacob_proyeccionDe(const L_CoordsCart3D &p, L_Matrix &J2x3, int i0, int j0) {J2x3(i0+0,j0+0) = -p.y/(p.x*p.x); J2x3(i0+0,j0+1) = 1/p.x; J2x3(i0+0,j0+2) = 0; J2x3(i0+1,j0+0) = -p.z/(p.x*p.x); J2x3(i0+1,j0+1) = 0; J2x3(i0+1,j0+2)=1/p.x;}
};

class L_Rayo3D
{
public:
	L_CoordsCart3D origen;
	L_CoordsCart3D direcc;
	L_Rayo3D() {}
	L_CoordsCart3D intersectionWith(const L_Rayo3D &other, L_CoordsCart3D *errorInt=NULL); // Para que tenga sentido, deben estar en el mismo sistema de coordenadas
	L_CoordsCart3D nearestPointTo(const L_CoordsCart3D &other) {return origen + direcc*direcc.punto(other-origen)/direcc.pideR2();}
	L_Rayo3D &OP_assign(const L_RayoCamara &r) {double invR=1.0/sqrt(1.0+r.tanIzq*r.tanIzq+r.tanArr*r.tanArr); origen.x=0; origen.y=0; origen.z=0; direcc.x=invR; direcc.y=r.tanIzq*invR; direcc.z=r.tanArr*invR; return *this;}
	bool normaliza() {return direcc.normaliza_ret();}
};


// Clases que representan rotaciones
// class_name intermedia para transformar rotaciones de forma rapida, reduccion de L_HomogeneousMatrix
class L_Quaternion;
class L_PanTiltRoll;

class L_PanTiltRoll
{
public:
	double pan;  // = yaw, es como el meridiano correspondiente (en radianes)
	double tilt; // = pitch, es como el paralelo correspondiente (en radianes)
	double roll; // Rotacion en torno al propio eje (en radianes)
	L_PanTiltRoll() {}
	L_PanTiltRoll(double panVal, double tiltVal, double rollVal):pan(panVal), tilt(tiltVal), roll(rollVal) {}
	double &el(int i) {switch(i) {case 0: return pan; case 1: return tilt; case 2: return roll; default: throw_L_ArgException_if(true,"L_PanTiltRoll::el"); return pan;}}
	void fijaRotacionCero() {pan=0; tilt=0; roll=0;}
	// ejeRot apunta en el sentido del eje de rotacion y su tamano es el angulo de rotacion
	double errorRespectoA(L_PanTiltRoll &other);
};

// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_Cuaternion_rotarVector(res,q,v) {double dqr_f2 = 1.0/((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); double dqr_ab = (q).a*(q).b*dqr_f2; double dqr_ac = (q).a*(q).c*dqr_f2;double dqr_ad = (q).a*(q).d*dqr_f2;double dqr_bb = (q).b*(q).b*dqr_f2;double dqr_bc = (q).b*(q).c*dqr_f2;double dqr_bd = (q).b*(q).d*dqr_f2;double dqr_cc = (q).c*(q).c*dqr_f2;double dqr_cd = (q).c*(q).d*dqr_f2;double dqr_dd = (q).d*(q).d*dqr_f2;(res).x = 2*( (-dqr_cc -dqr_dd)*(v).x + (dqr_bc - dqr_ad)*(v).y + (dqr_ac + dqr_bd)*(v).z ) + (v).x;(res).y = 2*( (dqr_ad + dqr_bc)*(v).x + (-dqr_bb -dqr_dd)*(v).y + (dqr_cd - dqr_ab)*(v).z ) + (v).y;(res).z = 2*( (dqr_bd - dqr_ac)*(v).x + (dqr_ab + dqr_cd)*(v).y + (-dqr_bb -dqr_cc)*(v).z ) + (v).z;}
#define L_Cuaternion_rotarVectorInv(res,q,v) {double dqr_f2 = 1.0/((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); double dqr_ab = (q).a*(q).b*dqr_f2; double dqr_ac = (q).a*(q).c*dqr_f2; double dqr_ad = (q).a*(q).d*dqr_f2; double dqr_bb = (q).b*(q).b*dqr_f2; double dqr_bc = (q).b*(q).c*dqr_f2; double dqr_bd = (q).b*(q).d*dqr_f2; double dqr_cc = (q).c*(q).c*dqr_f2; double dqr_cd = (q).c*(q).d*dqr_f2; double dqr_dd = (q).d*(q).d*dqr_f2; (res).x = 2*( (-dqr_cc -dqr_dd)*(v).x + (dqr_bc +dqr_ad)*(v).y + (-dqr_ac + dqr_bd)*(v).z ) + (v).x; (res).y = 2*( (-dqr_ad + dqr_bc)*(v).x + (-dqr_bb -dqr_dd)*(v).y + (dqr_cd +dqr_ab)*(v).z ) + (v).y; (res).z = 2*( (dqr_bd +dqr_ac)*(v).x + (-dqr_ab + dqr_cd)*(v).y + (-dqr_bb -dqr_cc)*(v).z ) + (v).z;}

#define L_Cuaternion_OP_igual(res,q) ((res).a=(q).a,(res).b=(q).b,(res).c=(q).c,(res).d=(q).d)
#define L_Cuaternion_OP_mult(res,qizq,qder) ( (res).a=(qizq).a*(qder).a - (qizq).b*(qder).b - (qizq).c*(qder).c - (qizq).d*(qder).d, (res).b=(qizq).a*(qder).b + (qizq).b*(qder).a + (qizq).c*(qder).d - (qizq).d*(qder).c, (res).c=(qizq).a*(qder).c - (qizq).b*(qder).d + (qizq).c*(qder).a + (qizq).d*(qder).b, (res).d=(qizq).a*(qder).d + (qizq).b*(qder).c - (qizq).c*(qder).b + (qizq).d*(qder).a )
#define L_Cuaternion_inverseOf(res,q) {double ir_lcid=1.0/((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); (res).a=ir_lcid*(q).a; (res).b=-ir_lcid*(q).b; (res).c=-ir_lcid*(q).c; (res).d=-ir_lcid*(q).d;}
#define L_Cuaternion_normaliza(q) {double ir_lcn; ir_lcn=1.0/sqrt((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); (q).a*=ir_lcn; (q).b*=ir_lcn; (q).c*=ir_lcn; (q).d*=ir_lcn;}

#define L_Cuaternion_jacob_izq_OP_mult(qizq,qder,J4x4,i0,j0) ((J4x4)(0+(i0),0+(j0)) = (qder).a, (J4x4)(0+(i0),1+(j0)) = -(qder).b, (J4x4)(0+(i0),2+(j0)) = -(qder).c, (J4x4)(0+(i0),3+(j0)) = -(qder).d, (J4x4)(1+(i0),0+(j0)) = (qder).b, (J4x4)(1+(i0),1+(j0)) = (qder).a, (J4x4)(1+(i0),2+(j0)) = (qder).d, (J4x4)(1+(i0),3+(j0)) = -(qder).c, (J4x4)(2+(i0),0+(j0)) = (qder).c, (J4x4)(2+(i0),1+(j0)) = -(qder).d, (J4x4)(2+(i0),2+(j0)) = (qder).a, (J4x4)(2+(i0),3+(j0)) = (qder).b, (J4x4)(3+(i0),0+(j0)) = (qder).d, (J4x4)(3+(i0),1+(j0)) = (qder).c, (J4x4)(3+(i0),2+(j0)) = -(qder).b, (J4x4)(3+(i0),3+(j0)) = (qder).a)
#define L_Cuaternion_jacob_der_OP_mult(qizq,qder,J4x4,i0,j0) ((J4x4)(0+(i0),0+(j0)) = (qizq).a, (J4x4)(0+(i0),1+(j0)) = -(qizq).b, (J4x4)(0+(i0),2+(j0)) = -(qizq).c, (J4x4)(0+(i0),3+(j0)) = -(qizq).d, (J4x4)(1+(i0),0+(j0)) = (qizq).b, (J4x4)(1+(i0),1+(j0)) = (qizq).a, (J4x4)(1+(i0),2+(j0)) = -(qizq).d, (J4x4)(1+(i0),3+(j0)) = (qizq).c, (J4x4)(2+(i0),0+(j0)) = (qizq).c, (J4x4)(2+(i0),1+(j0)) = (qizq).d, (J4x4)(2+(i0),2+(j0)) = (qizq).a, (J4x4)(2+(i0),3+(j0)) = -(qizq).b, (J4x4)(3+(i0),0+(j0)) = (qizq).d, (J4x4)(3+(i0),1+(j0)) = -(qizq).c, (J4x4)(3+(i0),2+(j0)) = (qizq).b, (J4x4)(3+(i0),3+(j0)) = (qizq).a)

// Explicacion sobre los jacobianos precalculados
// rot :   jizq_rot -> pre[36]    ,  jder_rot -> J3x3
// rot :   jizq_roti -> prem[36]  ,  jder_roti -> J3x3s

#define L_Cuaternion_calc_pre(q,pre36) ( (q).jacob_izq_rotarVector_pre(pre36) )
#define L_Cuaternion_calc_prem(q,pre36m) ( (q).jacob_izq_rotarVectorInv_pre(pre36m) )

#define L_Cuaternion_jacob_izq_rotarVector(v,pre,J3x4,i0,j0) ((J3x4)(0+i0,0+j0) = (pre)[0] * (v).x + (pre)[1] * (v).y + (pre)[2] * (v).z, (J3x4)(0+i0,1+j0) = (pre)[3] * (v).x + (pre)[4] * (v).y + (pre)[5] * (v).z, (J3x4)(0+i0,2+j0) = (pre)[6] * (v).x + (pre)[7] * (v).y + (pre)[8] * (v).z, (J3x4)(0+i0,3+j0) = (pre)[9] * (v).x + (pre)[10] * (v).y + (pre)[11] * (v).z, (J3x4)(1+i0,0+j0) = (pre)[12] * (v).x + (pre)[13] * (v).y + (pre)[14] * (v).z, (J3x4)(1+i0,1+j0) = (pre)[15] * (v).x + (pre)[16] * (v).y + (pre)[17] * (v).z,\
	(J3x4)(1+i0,2+j0) = (pre)[18] * (v).x + (pre)[19] * (v).y + (pre)[20] * (v).z, (J3x4)(1+i0,3+j0) = (pre)[21] * (v).x + (pre)[22] * (v).y + (pre)[23] * (v).z, (J3x4)(2+i0,0+j0) = (pre)[24] * (v).x + (pre)[25] * (v).y + (pre)[26] * (v).z, (J3x4)(2+i0,1+j0) = (pre)[27] * (v).x + (pre)[28] * (v).y + (pre)[29] * (v).z, (J3x4)(2+i0,2+j0) = (pre)[30] * (v).x + (pre)[31] * (v).y + (pre)[32] * (v).z, (J3x4)(2+i0,3+j0) = (pre)[33] * (v).x + (pre)[34] * (v).y + (pre)[35] * (v).z)

#define L_Cuaternion_jacob_der_rotarVector(q,J3x3,i0,j0) {double f2_lcjdrv = 1.0/((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); double ab_lcjdrv = (q).a*(q).b*f2_lcjdrv; double ac_lcjdrv = (q).a*(q).c*f2_lcjdrv; double ad_lcjdrv = (q).a*(q).d*f2_lcjdrv; double bb_lcjdrv = (q).b*(q).b*f2_lcjdrv; double bc_lcjdrv = (q).b*(q).c*f2_lcjdrv; double bd_lcjdrv = (q).b*(q).d*f2_lcjdrv; double cc_lcjdrv = (q).c*(q).c*f2_lcjdrv; double cd_lcjdrv = (q).c*(q).d*f2_lcjdrv; double dd_lcjdrv = (q).d*(q).d*f2_lcjdrv;\
	(J3x3)(0+(i0),0+(j0)) = 2*(-cc_lcjdrv -dd_lcjdrv) + 1; (J3x3)(0+(i0),1+(j0)) = 2*(bc_lcjdrv - ad_lcjdrv); (J3x3)(0+(i0),2+(j0)) = 2*(ac_lcjdrv + bd_lcjdrv); (J3x3)(1+(i0),0+(j0)) = 2*(ad_lcjdrv + bc_lcjdrv); (J3x3)(1+(i0),1+(j0)) = 2*(-bb_lcjdrv -dd_lcjdrv) + 1; (J3x3)(1+(i0),2+(j0)) = 2*(cd_lcjdrv - ab_lcjdrv); (J3x3)(2+(i0),0+(j0)) = 2*(bd_lcjdrv - ac_lcjdrv); (J3x3)(2+(i0),1+(j0)) = 2*(ab_lcjdrv + cd_lcjdrv); (J3x3)(2+(i0),2+(j0)) = 2*(-bb_lcjdrv -cc_lcjdrv) + 1;}

#define L_Cuaternion_jacob_izq_rotarVectorInv(v,prem,J3x4,i0,j0) ( (J3x4)(0+i0,0+j0) = (prem)[0] * (v).x + (prem)[1] * (v).y + (prem)[2] * (v).z, (J3x4)(0+i0,1+j0) = (prem)[3] * (v).x + (prem)[4] * (v).y + (prem)[5] * (v).z, (J3x4)(0+i0,2+j0) = (prem)[6] * (v).x + (prem)[7] * (v).y + (prem)[8] * (v).z, (J3x4)(0+i0,3+j0) = (prem)[9] * (v).x + (prem)[10] * (v).y + (prem)[11] * (v).z, (J3x4)(1+i0,0+j0) = (prem)[12] * (v).x + (prem)[13] * (v).y + (prem)[14] * (v).z, (J3x4)(1+i0,1+j0) = (prem)[15] * (v).x + (prem)[16] * (v).y + (prem)[17] * (v).z , \
	(J3x4)(1+i0,2+j0) = (prem)[18] * (v).x + (prem)[19] * (v).y + (prem)[20] * (v).z, (J3x4)(1+i0,3+j0) = (prem)[21] * (v).x + (prem)[22] * (v).y + (prem)[23] * (v).z, (J3x4)(2+i0,0+j0) = (prem)[24] * (v).x + (prem)[25] * (v).y + (prem)[26] * (v).z, (J3x4)(2+i0,1+j0) = (prem)[27] * (v).x + (prem)[28] * (v).y + (prem)[29] * (v).z, (J3x4)(2+i0,2+j0) = (prem)[30] * (v).x + (prem)[31] * (v).y + (prem)[32] * (v).z , (J3x4)(2+i0,3+j0) = (prem)[33] * (v).x + (prem)[34] * (v).y + (prem)[35] * (v).z)

#define L_Cuaternion_jacob_der_rotarVectorInv(q,J3x3,i0,j0) { double f2_lcji = 1.0/((q).a*(q).a+(q).b*(q).b+(q).c*(q).c+(q).d*(q).d); double ab_lcji = (q).a*(q).b*f2_lcji; double ac_lcji = (q).a*(q).c*f2_lcji; double ad_lcji = (q).a*(q).d*f2_lcji; double bb_lcji = (q).b*(q).b*f2_lcji; double bc_lcji = (q).b*(q).c*f2_lcji; double bd_lcji = (q).b*(q).d*f2_lcji; double cc_lcji = (q).c*(q).c*f2_lcji; double cd_lcji = (q).c*(q).d*f2_lcji; double dd_lcji = (q).d*(q).d*f2_lcji; (J3x3)((i0)+0,(j0)+0) = -2*cc_lcji-2*dd_lcji+1; (J3x3)((i0)+0,(j0)+1) = 2*bc_lcji+2*ad_lcji; (J3x3)((i0)+0,(j0)+2) = -2*ac_lcji+2*bd_lcji; (J3x3)((i0)+1,(j0)+0) = -2*ad_lcji+2*bc_lcji;\
	(J3x3)((i0)+1,(j0)+1) = -2*bb_lcji-2*dd_lcji+1; (J3x3)((i0)+1,(j0)+2) = 2*cd_lcji+2*ab_lcji; (J3x3)((i0)+2,(j0)+0) = 2*bd_lcji+2*ac_lcji; (J3x3)((i0)+2,(j0)+1) = -2*ab_lcji+2*cd_lcji; (J3x3)((i0)+2,(j0)+2) = -2*bb_lcji-2*cc_lcji+1;}

// q.rot(p) = q*p*q^-1
class L_Quaternion
{
public:
	double a, b, c, d; // q = a + bi + cj + dk
	L_Quaternion() {}
	L_Quaternion(double aVal, double bVal, double cVal, double dVal) : a(aVal), b(bVal), c(cVal), d(dVal) {}

	double &el(int i) {switch(i) {case 0: return a; case 1: return b; case 2: return c; case 3: return d; default: throw_L_ArgException_if(true,"L_Quaternion::el"); return a;}}

	void setZero() {a=1; b=0; c=0; d=0;} // Sin rotacion
	void fijaRotacionCero() {a=1; b=0; c=0; d=0;}

	void define(double aVal, double bVal, double cVal, double dVal) {a = aVal; b = bVal; c = cVal; d = dVal;}
	void setRotationComponents(double theta, double x, double y, double z) {a=cos(theta/2); double s=sin(theta/2); b=x*s; c=y*s; d=z*s;}
	void normaliza() {double r; r=sqrt(a*a+b*b+c*c+d*d); if(r==0 || r!=r) {printf("Error en L_Quaternion::normalizar() : cuaternion == 0 => NaN\n"); *(char *)NULL = 0;} a/=r; b/=r; c/=r; d/=r;}
	bool normaliza_ret() {double r; r=sqrt(a*a+b*b+c*c+d*d); a/=r; b/=r; c/=r; d/=r; return r==r && r!=0;}
	double abs() {return sqrt(a*a+b*b+c*c+d*d);}
	double abs2() {return a*a+b*b+c*c+d*d;}
	bool isNaN() {return a!=a || b!=b || c!=c || d!=d;}

	L_CoordsCart3D pideImag() {L_CoordsCart3D p; p.x=b; p.y=c; p.z=d; return p;}
	void set(double r, L_CoordsCart3D imag) {a=r; b=imag.x; c=imag.y; d=imag.z;}
	void fijaAzar() {a=L_RANDOM(-1,1); b=L_RANDOM(-1,1); c=L_RANDOM(-1,1); d=L_RANDOM(-1,1); L_Cuaternion_normaliza(*this);}

	// rotarVector_rapido() solo funciona con el cuaternion ya normalizado, rotarVector() funciona siempre
	L_CoordsCart3D rotarVector(const L_CoordsCart3D &inicial) const;  // = q * p * q^-1
	L_CoordsCart3D rotarVector_lento(const L_CoordsCart3D &inicial) const {L_Quaternion p, q, qi, qp, qpqi; p.set(0,inicial); q=*this; qi.inverseOf(*this); qp.OP_mult(q,p);  qpqi.OP_mult(qp,qi); return qpqi.pideImag();}

	void jacob_izq_rotarVector(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0=0, int j0=0) const; // Requiere *this ya calculado
	void jacob_izq_rotarVector_pre(double *pre) const; // Requiere *this ya calculado, precalcula cosas
	static void jacob_izq_rotarVector(const L_CoordsCart3D &inicial, const double *pre, L_StaticMatrix<3,4> &J3x4); // Requiere *this ya calculado
	void jacob_der_rotarVector(L_Matrix &J3x3, int i0=0, int j0=0) const; // Requiere *this ya calculado

	L_CoordsCart3D rotarVector_rapido(const L_CoordsCart3D &inicial) const;
	void jacob_izq_rotarVector_rapido(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0=0, int j0=0) const; // Requiere *this ya calculado
	void jacob_der_rotarVector_rapido(L_Matrix &J3x3, int i0=0, int j0=0) const; // Requiere *this ya calculado

	L_CoordsCart3D rotarVectorInv(const L_CoordsCart3D &inicial) const; // inversa = q^1 * p * q
	L_CoordsCart3D rotarVectorInv_lento(const L_CoordsCart3D &inicial) const {L_Quaternion p, q, qi, qip, qipq; p.set(0,inicial); q=*this; qi.inverseOf(*this); qip.OP_mult(qi,p);  qipq.OP_mult(qip,q); return qipq.pideImag();}

	void jacob_izq_rotarVectorInv(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0=0, int j0=0) const; // Requiere cuaternión normalizado
	void jacob_izq_rotarVectorInv_pre(double *pre) const; // Requiere *this ya calculado, precalcula cosas
	static void jacob_izq_rotarVectorInv(const L_CoordsCart3D &inicial, const double *pre, L_StaticMatrix<3,4> &J3x4); // Requiere *this ya calculado
	void jacob_der_rotarVectorInv(L_Matrix &J3x3, int i0=0, int j0=0) const; // Requiere cuaternión normalizado

	L_CoordsCart3D rotarVectorInv_rapido(const L_CoordsCart3D &inicial) const; // Requiere cuaternión normalizado
	void jacob_izq_rotarVectorInv_rapido(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0=0, int j0=0) const; // Requiere cuaternión normalizado
	void jacob_der_rotarVectorInv_rapido(L_Matrix &J3x3, int i0=0, int j0=0) const; // Requiere cuaternión normalizado

	// R(q1*q2)=R(q1)*R(q2)
	void OP_add(const L_Quaternion &izq, const L_Quaternion &der) {a=izq.a+der.a; b=izq.b+der.b; c=izq.c+der.c; d=izq.d+der.d;}
	void OP_subtract(const L_Quaternion &izq, const L_Quaternion &der) {a=izq.a-der.a; b=izq.b-der.b; c=izq.c-der.c; d=izq.d-der.d;}
	void OP_mult(const L_Quaternion &izq, const L_Quaternion &der);
	void OP_div(const L_Quaternion &izq, const L_Quaternion &der);
	
	void copiaSignoDe(L_Quaternion &other) {if (a*other.a + b*other.b + c*other.c + d*other.d < 0) {a=-a; b=-b; c=-c; d=-d;}}
	void invSigno() {a=-a; b=-b; c=-c; d=-d;} // No altera la rotacion representada por el cuaternion
	static void jacob_izq_OP_mult(const L_Quaternion &izq, const L_Quaternion &der, L_Matrix &J4x4, int i0=0, int j0=0); // No resize J!
	static void jacob_der_OP_mult(const L_Quaternion &izq, const L_Quaternion &der, L_Matrix &J4x4, int i0=0, int j0=0); // No resize J!

	void jacob_normaliza(L_Matrix &J4x4);

	double cosineDistanceTo(L_Quaternion &other) const {L_Quaternion q = other; if (OP_mult_elementwise(q) < 0) q.invSigno(); return (a*q.a+b*q.b+c*q.c+d*q.d) / sqrt((a*a+b*b+c*c+d*d)*(q.a*q.a+q.b*q.b+q.c*q.c+q.d*q.d));}

	double OP_mult_elementwise(const L_Quaternion &other) const {return a*other.a + b*other.b + c*other.c + d*other.d;}

	void inverseOf(const L_Quaternion &other) {double ir=1/(other.a*other.a+other.b*other.b+other.c*other.c+other.d*other.d); a=ir*other.a; b=-ir*other.b; c=-ir*other.c; d=-ir*other.d;}
	void jacob_inverseOf(const L_Quaternion &other, L_Matrix &J4x4, int i0=0, int j0=0);
	void conjugadoDe(const L_Quaternion &other) {a=other.a; b=-other.b; c=-other.c; d=-other.d;}

	void fijaVectorRotacion(const L_CoordsCart3D &ejeRot) {double theta = ejeRot.pideR(); double mod = theta; if (theta == 0) {a=1; b=0; c=0; d=0;} else {double s = sin(theta/2); a = cos(theta/2); b = s*(ejeRot.x/mod); c = s*(ejeRot.y/mod); d=s*(ejeRot.z/mod);}}
	void fijaVectorRotacion(const L_CoordsCart3D &ejeRotUn, double theta) {if (theta == 0) {a=1; b=0; c=0; d=0;} else {double s = sin(theta/2); a = cos(theta/2); b = ejeRotUn.x*s; c= ejeRotUn.y*s; d=ejeRotUn.z*s;}} // Requiere ejeRot normalizado
	L_CoordsCart3D pideVectorRotacion() {L_CoordsCart3D p; double bcd = sqrt(b*b + c*c + d*d); if (bcd == 0) return L_CoordsCart3D(0,0,0); double ang = atan2(bcd, a); p.x = b/bcd*2*ang; p.y=c/bcd*2*ang; p.z=d/bcd*2*ang; return p;}

	void pow(L_Quaternion &other, double exp) {fijaVectorRotacion(exp*pideVectorRotacion());} // Para slerp

	bool fijaDosPuntos(const L_CoordsCart3D &fijo, L_CoordsCart3D &rotado); // Asegura que rotado = q'*fijo*q, pero nada mas
	double fijaCuatroPuntos(const L_CoordsCart3D &fijo1, const L_CoordsCart3D &fijo2, const L_CoordsCart3D &rotado1, const L_CoordsCart3D &rotado2);
	double fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, int niter = 3); // Mas ambicioso
	double fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido); // Requiere invertir matriz
	double jacob_fijo_fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int niter = 3, int i0=0, int j0=0); // Mas ambicioso
	double jacob_fijo_fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int i0=0, int j0=0); // Requiere invertir matriz
	double jacob_movido_fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int niter = 3, int i0=0, int j0=0); // Mas ambicioso
	double jacob_movido_fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int i0=0, int j0=0); // Requiere invertir matriz

	static void pruebaCuaternion();
	static void pruebaCuaternion_def();

	void print(const char *name, FILE *fp = stdout) {fprintf(fp, "%s = %g + %g*i + %g*j + %g*k\n", name, a, b, c, d);}
	void print(FILE *fp = stdout) {fprintf(fp, "%g + %g*i + %g*j + %g*k\n", a, b, c, d);}
};

inline L_Quaternion operator+(const L_Quaternion &izq, const L_Quaternion &der) {L_Quaternion ret; ret.OP_add(izq, der); return ret;} // Es asimetrico
inline L_Quaternion operator-(const L_Quaternion &izq, const L_Quaternion &der) {L_Quaternion ret; ret.OP_subtract(izq, der); return ret;} // Es asimetrico
inline L_Quaternion operator*(const L_Quaternion &izq, const L_Quaternion &der) {L_Quaternion ret; ret.OP_mult(izq, der); return ret;} // Es asimetrico
inline L_Quaternion operator/(const L_Quaternion &izq, const L_Quaternion &der) {L_Quaternion ret; ret.OP_div(izq, der); return ret;} // Es asimetrico
inline L_Quaternion operator-(const L_Quaternion &q) {L_Quaternion ret; ret.a=-q.a; ret.b=-q.b; ret.c=-q.c; ret.d=-q.d; return ret;} // Es asimetrico

class L_Pose2D
{
public:
	L_CoordsCart2D pos;
	double ori;
};

// Clases para trabajar con poses 3D

class L_Pose3D;
class L_Pose3D_cuat;
class L_HomogeneousMatrix;


class L_Pose3D
{
public:
	L_CoordsCart3D pos;
	L_PanTiltRoll ori;
	L_Pose3D() {}
	L_Pose3D(double x, double y, double z, double pan, double tilt, double roll) : pos(x,y,z),ori(pan,tilt,roll) {}
	L_Pose3D(L_CoordsCart3D &posic, L_PanTiltRoll &orien) : pos(posic), ori(orien) {}
	double &el(int i) {switch(i) {case 0: return pos.x; case 1: return pos.y; case 2: return pos.z; case 3: return ori.pan; case 4: return ori.tilt; case 5: return ori.roll; default:throw_L_ArgException_if(true, "L_Pose3D::el()");return pos.x;}}
	const double &el(int i) const {switch(i) {case 0: return pos.x; case 1: return pos.y; case 2: return pos.z; case 3: return ori.pan; case 4: return ori.tilt; case 5: return ori.roll; default:throw_L_ArgException_if(true, "L_Pose3D::el()");return pos.x;}}
	void fijaPose(double x, double y, double z, double pan, double tilt, double roll) {pos.x=x; pos.y=y; pos.z=z; ori.pan=pan; ori.tilt=tilt; ori.roll=roll;}
	void fijaAzar(double rMax=1000, double angMax=M_PI) {pos.x=L_RANDOM(-rMax,rMax); pos.y=L_RANDOM(-rMax,rMax); pos.z=L_RANDOM(-rMax,rMax); ori.pan=L_RANDOM(-angMax,angMax); ori.tilt=L_RANDOM(-angMax,angMax); ori.roll=L_RANDOM(-angMax,angMax);}
	L_Pose3D &fijaPose_abs(double x, double y, double z, double xAde, double yAde, double zAde, double xArr, double yArr, double zArr); // Como gluLookAt()
	L_Pose3D &fijaPose_rel(double x, double y, double z, double xDirAde, double yDirAde, double zDirAde, double xDirArr, double yDirArr, double zDirArr); // Como gluLookAt()
	void fijaPosOri(const L_CoordsCart3D &pos, const L_PanTiltRoll &ori) {this->pos=pos; this->ori=ori;}
	void fijaCero() {pos.fijaCero(); ori.fijaRotacionCero();}
	void llenaArreglo(double *a) {a[0]=pos.x; a[1]=pos.y; a[2]=pos.z; a[3]=ori.pan; a[4]=ori.tilt; a[5]=ori.roll;}
	void leeArreglo(const double *a) {pos.x=a[0]; pos.y=a[1]; pos.z=a[2]; ori.pan=a[3]; ori.tilt=a[4]; ori.roll=a[5];}
	void print(const char *name, FILE *fp=stdout) {fprintf(fp, "%s = [%.6g %.6g %.6g pan=%.6g tilt=%.6g roll=%.6g]\n\n", name, pos.x, pos.y, pos.z, ori.pan, ori.tilt, ori.roll);}
	void imprime_prec(const char *name, FILE *fp=stdout) {fprintf(fp, "%s = [%.12g %.12g %.12g pan=%.12g tilt=%.12g roll=%.12g]\n\n", name, pos.x, pos.y, pos.z, ori.pan, ori.tilt, ori.roll);}
	void imprime_gr(const char *name, FILE *fp=stdout) {fprintf(fp, "%s = [%.3g %.3g %.3g pan=%.3g[gr] tilt=%.3g[gr] roll=%.3g[gr]]\n\n", name, pos.x, pos.y, pos.z, ori.pan*180/M_PI, ori.tilt*180/M_PI, ori.roll*180/M_PI);}
	void agregaDeltaPose(const L_Pose3D &deltaPose);

	static void pruebaPose3D();
};

class L_Pose3D_debug
{
public:
	L_CoordsCart3D pos;
	double panGrados;
	double tiltGrados;
	double rollGrados;

	L_Pose3D_debug() {pos.x=666; pos.y=666; pos.z=666; panGrados = 666; tiltGrados = 666; rollGrados = 666;}
	L_Pose3D_debug(const L_Pose3D &p) {pos=p.pos; panGrados = p.ori.pan*180/M_PI; tiltGrados = p.ori.tilt*180/M_PI;rollGrados = p.ori.roll*180/M_PI;}
	void operator=(const L_Pose3D &p) {pos=p.pos; panGrados = p.ori.pan*180/M_PI; tiltGrados = p.ori.tilt*180/M_PI;rollGrados = p.ori.roll*180/M_PI;}

};

class L_Pose3D_cuat_nube_puntos
{
public:
	const std::vector<L_CoordsCart3D> *movido;
	const std::vector<L_CoordsCart3D> *fijo;
};

class L_Pose3D_cuat_nube_puntos_pesos
{
public:
	const std::vector<L_CoordsCart3D> *movido;
	const std::vector<L_CoordsCart3D> *fijo;
	const std::vector<double> *pesos;
};

class L_Pose3D_cuat_rayos_rmpf
{
public:
	const L_Array<const L_Pose3D_cuat *> *pCam;
	const std::vector<L_RayoCamara> *movido;
	const std::vector<L_CoordsCart3D> *fijo;
};

class L_Pose3D_cuat_rayos_pmrf
{
public:
	const L_Array<const L_Pose3D_cuat *> *pCam21;
	const std::vector<L_CoordsCart3D> *movido;
	const std::vector<L_RayoCamara> *fijo;
	int comp; // Por si se necesita calcular solo una componente
};

// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_Pose3D_cuat_movido_a_fijo(res,p3d,v) {L_Cuaternion_rotarVector(res,(p3d).ori,v); res.x+=(p3d).pos.x; res.y+=(p3d).pos.y; res.z+=(p3d).pos.z;} // Lo normal
#define L_Pose3D_cuat_fijo_a_movido(res,p3d,v) {L_CoordsCart3D resta_dp3dc; resta_dp3dc.x=(v).x-(p3d).pos.x; resta_dp3dc.y=(v).y-(p3d).pos.y; resta_dp3dc.z=(v).z-(p3d).pos.z; L_Cuaternion_rotarVectorInv(res,(p3d).ori,resta_dp3dc);} // RT(p-t)
#define L_Pose3D_cuat_inverseOf(res,p3d) {L_Cuaternion_inverseOf((res).ori,(p3d).ori); L_Cuaternion_rotarVector((res).pos,(res).ori,(p3d).pos); (res).pos.x=-(res).pos.x; (res).pos.y=-(res).pos.y; (res).pos.z=-(res).pos.z;}
#define L_Pose3D_cuat_normaliza(p3d) {L_Cuaternion_normaliza(p3d.ori);}
#define L_Pose3D_cuat_OP_igual(res,p3d) (L_CoordsCart3D_OP_igual((res).pos,(p3d).pos),L_Cuaternion_OP_igual((res).ori,(p3d).ori))
#define L_Pose3D_cuat_OP_mult(res,p3d1,p3d2) {L_Cuaternion_rotarVector((res).pos,(p3d1).ori,(p3d2).pos); (res).pos.x+=p3d1.pos.x; (res).pos.y+=p3d1.pos.y; (res).pos.z+=p3d1.pos.z; L_Cuaternion_OP_mult((res).ori,(p3d1).ori,(p3d2).ori);}
#define L_Pose3D_cuat_OP_mult_p1_por_p2inv(res,p3d1,p3d2) {L_Pose3D_cuat p2inv_lcomppi; L_Pose3D_cuat_inverseOf(p2inv_lcomppi,(p3d2)); L_Pose3D_cuat_OP_mult((res),(p3d1),p2inv_lcomppi);}
//
#define L_Pose3D_cuat_jacob_fijo_a_movido_pre(p3d,J3x3s,pre36m) {L_Cuaternion_jacob_der_rotarVectorInv((p3d).ori,(J3x3s),0,0); (p3d).ori.jacob_izq_rotarVectorInv_pre((pre36m)); int jiri_i; for(jiri_i=0; jiri_i<36; jiri_i++) pre36m[jiri_i]=-pre36m[jiri_i];}
#define L_Pose3D_cuat_jacob_izq_fijo_a_movido(p3d,v,prem,J3x3s,J3x7,i0,j0) {L_CoordsCart3D resta_dp3dc; resta_dp3dc.x=(p3d).pos.x-(v).x; resta_dp3dc.y=(p3d).pos.y-(v).y; resta_dp3dc.z=(p3d).pos.z-(v).z; J3x7((i0)+0,(j0)+0) = -J3x3s(0,0); J3x7((i0)+0,(j0)+1) = -J3x3s(0,1); J3x7((i0)+0,(j0)+2) = -J3x3s(0,2); J3x7((i0)+1,(j0)+0) = -J3x3s(1,0); J3x7((i0)+1,(j0)+1) = -J3x3s(1,1); J3x7((i0)+1,(j0)+2) = -J3x3s(1,2); J3x7((i0)+2,(j0)+0) = -J3x3s(2,0); J3x7((i0)+2,(j0)+1) = -J3x3s(2,1); J3x7((i0)+2,(j0)+2) = -J3x3s(2,2); L_Cuaternion_jacob_izq_rotarVectorInv(resta_dp3dc,prem,J3x7,(i0)+0,j0+3);}
#define L_Pose3D_cuat_jacob_der_fijo_a_movido(p3d,J3x3,i1,j1) {L_Cuaternion_jacob_der_rotarVectorInv((p3d).ori,(J3x3),i1,j1);}
// inutiles...
#define L_Pose3D_cuat_jacob_movido_a_fijo_pre(p3d,J3x3d,pre36) {L_Cuaternion_jacob_der_rotarVector((p3d).ori,(J3x3d),0,0); (p3d).ori.jacob_izq_rotarVector_pre((pre36));}
#define L_Pose3D_cuat_jacob_izq_movido_a_fijo(p3d,movido,pre,J3x7) {L_StaticMatrix_identity3x3(J3x7); L_Cuaternion_jacob_izq_rotarVector(movido,pre,J3x7,0,3);}
#define L_Pose3D_cuat_jacob_der_movido_a_fijo(p3d,J3x3d) {L_Cuaternion_jacob_der_rotarVector(p3d.ori,J3x3d,0,0);}

//f(x) = -rotarInv(px-x);
//df/dx = -drotarInv(px-x)/dX * -1 = dRotarInv(px-x)/dX

// Ahora las funciones de proyeccion para el SLAM
// EN EL SLAM LAS PROYECCIONES SON FIJO_A_MOVIDO

//#define L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_pre(p3d,J3x3s,pre36m) L_Pose3D_cuat_jacob_fijo_a_movido_pre(...)
//
#define L_Pose3D_cuat_proyectar_fijo_a_movido_v(res,p3d,v) {L_CoordsCart3D movi_lp3dcpmf; L_Pose3D_cuat_fijo_a_movido(movi_lp3dcpmf,(p3d),v); (res).tanIzq=movi_lp3dcpmf.y/movi_lp3dcpmf.x; (res).tanArr=movi_lp3dcpmf.z/movi_lp3dcpmf.x;}
#define L_Pose3D_cuat_proyectar_fijo_a_movido_vAde(res,p3d,v,xs) {L_CoordsCart3D movi_lp3dcpmf; L_Pose3D_cuat_fijo_a_movido(movi_lp3dcpmf,(p3d),(v)); xs=movi_lp3dcpmf.x; (res).tanIzq=movi_lp3dcpmf.y/movi_lp3dcpmf.x; (res).tanArr=movi_lp3dcpmf.z/movi_lp3dcpmf.x;}
#define L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v(p3d,v,J3x3s,pre36m,J2x7,i0,j0,J2x3,i1,j1) {L_StaticMatrix<2,3> J2x3s; L_StaticMatrix<3,7> J3x7; L_CoordsCart3D movi_xd; L_Pose3D_cuat_fijo_a_movido(movi_xd,(p3d),v); L_Pose3D_cuat_jacob_izq_fijo_a_movido((p3d),(v),(pre36m),J3x3s,J3x7,0,0); L_RayoCamara_jacob_proyeccionDe(movi_xd,J2x3s,0,0); L_StaticMatrix_OP_mult_sub((J2x7),(i0),(j0),J2x3s,J3x7,0,0,2,3,0,0,3,7); L_StaticMatrix_OP_mult_sub((J2x3),(i1),(j1),J2x3s,(J3x3s),0,0,2,3,0,0,3,3);}
//
#define L_Pose3D_cuat_proyectar_fijo_a_movido_vi(res,p3d,vi) {L_CoordsCart3D pos_lpmfi; L_CoordsInv3D_pideCoords3D(pos_lpmfi,vi); L_Pose3D_cuat_proyectar_fijo_a_movido_v(res,p3d,pos_lpmfi);}
#define L_Pose3D_cuat_proyectar_fijo_a_movido_viAde(res,p3d,vi,xs) {L_CoordsCart3D pos_lpmfi; L_CoordsInv3D_pideCoords3D(pos_lpmfi,vi); L_Pose3D_cuat_proyectar_fijo_a_movido_vAde(res,p3d,pos_lpmfi,xs);}
#define L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_vi(p3d,vi,J3x3s,pre36m,J2x7,i0,j0,J2x6,i1,j1)  {L_StaticMatrix<2,3> J2x3_ljpmfi; L_StaticMatrix<3,6> J3x6_ljpmfi; L_CoordsCart3D pos_ljpmfi; L_CoordsInv3D_pideCoords3D(pos_ljpmfi,(vi)); L_CoordsInv3D_jacob_pideCoords3D((vi),J3x6_ljpmfi,0,0); L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v((p3d),pos_ljpmfi,(J3x3s),(pre36m),(J2x7),(i0),(j0),J2x3_ljpmfi,0,0); L_StaticMatrix_OP_mult_sub((J2x6),(i1),(j1),J2x3_ljpmfi,J3x6_ljpmfi,0,0,2,3,0,0,3,6);}

// Funciones inutiles (en su mayoria)
#define L_Pose3D_cuat_proyectar_movido_a_fijo_v(res,p3d,v) {L_CoordsCart3D fijo_lp3dcpmf; L_Pose3D_cuat_movido_a_fijo(fijo_lp3dcpmf,p3d,v); res.tanIzq=fijo_lp3dcpmf.y/fijo_lp3dcpmf.x; res.tanArr=fijo_lp3dcpmf.z/fijo_lp3dcpmf.x;}
#define L_Pose3D_cuat_proyectar_movido_a_fijo_vAde(res,p3d,v,xs) {L_CoordsCart3D fijo_lp3dcpmf; L_Pose3D_cuat_movido_a_fijo(fijo_lp3dcpmf,(p3d),(v)); xs=fijo_lp3dcpmf.x; (res).tanIzq=fijo_lp3dcpmf.y/fijo_lp3dcpmf.x; (res).tanArr=fijo_lp3dcpmf.z/fijo_lp3dcpmf.x;}
#define L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_v(p3d,v,J3x3,pre36,J2x7,i0,j0,J2x3,i1,j1) {L_StaticMatrix<2,3> J2x3s; L_StaticMatrix<3,7> J3x7; L_CoordsCart3D fijo_xd; L_Pose3D_cuat_movido_a_fijo(fijo_xd,(p3d),v); L_Pose3D_cuat_jacob_izq_movido_a_fijo((p3d),v,(pre36),J3x7); L_RayoCamara_jacob_proyeccionDe(fijo_xd,J2x3s,0,0); L_StaticMatrix_OP_mult_sub((J2x7),(i0),(j0),J2x3s,J3x7,0,0,2,3,0,0,3,7); L_StaticMatrix_OP_mult_sub((J2x3),(i1),(j1),J2x3s,(J3x3),0,0,2,3,0,0,3,3);}
//
#define L_Pose3D_cuat_proyectar_movido_a_fijo_vi(res,p3d,vi) {L_CoordsCart3D pos_lpmfi; L_CoordsInv3D_pideCoords3D(pos_lpmfi,vi); L_Pose3D_cuat_proyectar_movido_a_fijo_v(res,p3d,pos_lpmfi);}
#define L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_vi(p3d,vi,J3x3,pre36,J2x7,i0,j0,J2x6,i1,j1)  {L_StaticMatrix<2,3> J2x3_ljpmfi; L_StaticMatrix<3,6> J3x6_ljpmfi; L_CoordsCart3D pos_ljpmfi; L_CoordsInv3D_pideCoords3D(pos_ljpmfi,(vi)); L_CoordsInv3D_jacob_pideCoords3D((vi),J3x6_ljpmfi,0,0); L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_v((p3d),pos_ljpmfi,(J3x3),(pre36),(J2x7),(i0),(j0),J2x3_ljpmfi,0,0); L_StaticMatrix_OP_mult_sub((J2x6),(i1),(j1),J2x3_ljpmfi,J3x6_ljpmfi,0,0,2,3,0,0,3,6);}


class L_Pose3D_cuat
{
public:
	L_CoordsCart3D pos;
	L_Quaternion ori;
	L_Pose3D_cuat() {}
	L_Pose3D_cuat(const L_CoordsCart3D &posic, const L_Quaternion &orient)  : pos(posic), ori(orient) {}
	L_Pose3D_cuat(double x, double y, double z, double a, double b, double c, double d) {pos.x=x; pos.y=y; pos.z=z; ori.a=a; ori.b=b; ori.c=c; ori.d=d;}

	double &el(int i) {switch(i) {case 0: return pos.x; case 1: return pos.y; case 2: return pos.z; case 3: return ori.a;  case 4: return ori.b;  case 5: return ori.c;  case 6: return ori.d; default: throw_L_ArgException_if(true,"L_Pose3D_cuat::el"); return pos.x;}}
	const double &el(int i) const {switch(i) {case 0: return pos.x; case 1: return pos.y; case 2: return pos.z; case 3: return ori.a;  case 4: return ori.b;  case 5: return ori.c;  case 6: return ori.d; default: throw_L_ArgException_if(true,"L_Pose3D_cuat::el"); return pos.x;}}
	int ncomp() {return 7;}

	void swap(L_Pose3D_cuat &other) {L_Pose3D_cuat t = other; other = *this; *this = t;}

	void arr_a_puntocuat(const double *arr) {pos.x = arr[0]; pos.y = arr[1]; pos.z = arr[2]; ori.a = arr[3]; ori.b = arr[4]; ori.c = arr[5]; ori.d = arr[6];}
	void matriz_a_puntocuat(const L_Matrix &m) {pos.x = m(0,0); pos.y = m(1,0); pos.z = m(2,0); ori.a = m(3,0); ori.b = m(4,0); ori.c = m(5,0); ori.d = m(6,0);}
	void puntocuat_a_arr(double *arr) {arr[0] = pos.x; arr[1] = pos.y; arr[2] = pos.z; arr[3] = ori.a; arr[4] = ori.b; arr[5] = ori.c; arr[6] = ori.d;}
	void puntocuat_a_matriz(L_Matrix &m) {m(0,0) = pos.x; m(1,0) = pos.y; m(2,0) = pos.z; m(3,0) = ori.a; m(4,0) = ori.b; m(5,0) = ori.c; m(6,0) = ori.d;}

	void fijaPosOri(const L_CoordsCart3D &pos, const L_Quaternion &ori) {this->pos=pos; this->ori=ori;}
	void mirarA(double xMirado, double yMirado, double zMirado, double roll=0);
	void fijaPoseMirarA(double x, double y, double z, double xMirado, double yMirado, double zMirado, double roll=0);
	void fijaCero() {pos.x=0; pos.y=0; pos.z=0; ori.a=1; ori.b=0; ori.c=0; ori.d=0;}
	void fijaAzar() {pos.fijaAzar(); ori.fijaAzar();}
	void define(double x, double y, double z, double a, double b, double c, double d) {pos.x=x; pos.y=y; pos.z=z; ori.a=a; ori.b=b; ori.c=c; ori.d=d;}
	void definePose3D(L_Pose3D &other) {definePosPanTiltRoll(other.pos.x, other.pos.y, other.pos.z, other.ori.pan, other.ori.tilt, other.ori.roll);}
	void definePosPanTiltRoll(double x, double y, double z, double pan, double tilt, double roll);
	void print(const char *name, FILE *fp=stdout) {fprintf(fp, "%s = [%.4g %.4g %.4g] ori=%.4g + %.4gi + %.4gj + %.4gk\n\n", name, pos.x, pos.y, pos.z, ori.a, ori.b, ori.c, ori.d);}
	void imprimeGrados(const char *name, FILE *fp=stdout);
	void normaliza() {ori.normaliza();}
	bool normaliza_ret() {return ori.normaliza_ret();}

	void copiaSignoDe(L_Pose3D_cuat &other) {ori.copiaSignoDe(other.ori);}

	bool isNaN() {return pos.x!=pos.x || pos.y!=pos.y || pos.z!=pos.z || ori.a!=ori.a || ori.b!=ori.b || ori.c!=ori.c || ori.d!=ori.d || ori.a*ori.a + ori.b*ori.b + ori.c*ori.c + ori.d*ori.d == 0;}

	void inverseOf(const L_Pose3D_cuat &eta);
	void jacob_inverseOf(const L_Pose3D_cuat &eta, L_Matrix &J7x7); // Calcula *this por la ambiguedad de signo

	// Cambiar el signo (o el módulo) del cuaternion no altera la rotacion que representa
	void cambiaSigno(L_Pose3D_cuat &other) {pos = other.pos; ori.a = -other.ori.a; ori.b=-other.ori.b; ori.c=-other.ori.c; ori.d=-other.ori.d;}
	void jacob_cambiaSigno(L_Pose3D_cuat &other, L_Matrix &J7x7) {J7x7.reallocate(7,7); J7x7.identity(); J7x7(3,3) = -1; J7x7(4,4) = -1; J7x7(5,5) = -1; J7x7(6,6) = -1;}
	void jacob_cambiaSigno_porM(L_Pose3D_cuat &other, L_Matrix &J7x7, L_Matrix &m7xn, L_Matrix &J7x7m7xn) {int i, j; J7x7.reallocate(7,7); J7x7.identity(); J7x7(3,3) = -1; J7x7(4,4) = -1; J7x7(5,5) = -1; J7x7(6,6) = -1; J7x7m7xn.reallocate(m7xn.li, m7xn.lj); for (i=0; i<7; i++) for (j=0; j<m7xn.lj; j++) {J7x7m7xn(i,j) = m7xn(i,j) * (1 - 2*(i>=3));}}

	void OP_suma_punto_a_punto(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2);
	void OP_resta_punto_a_punto(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2);

	// Funciones directas, trabajan con q
	// Dadas vuelta
	void OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2); // Composicion directa, H(eta1*eta2) = H(eta2)*H(eta1)
	void OP_div(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2); // Composicion directa
	static void jacob_izq_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0); // Jacobiano c/r a eta1
	static void jacob_der_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0); // Jacobiano c/r a eta2
	static void jacob_izq_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_StaticMatrix<7,7> &J7x7); // Jacobiano c/r a eta1
	static void jacob_der_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_StaticMatrix<7,7> &J7x7); // Jacobiano c/r a eta1
	static void jacob_izq_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, const double *pre_p2_36, L_StaticMatrix<7,7> &J7x7); // Jacobiano c/r a eta2
	//static void jacob_der_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, const double *pre_p2_36, L_StaticMatrix<7,7> &J7x7); // Jacobiano c/r a eta2
	//static void jacob_der_OP_mult_pre(const L_Pose3D_cuat &eta2, double *pre_p2_36) {L_Cuaternion_calc_pre(eta2.ori,pre_p2_36);}
	static void jacob_izq_OP_mult_pre(const L_Pose3D_cuat &eta2, double *pre_p2_36) {L_Cuaternion_calc_pre(eta2.ori,pre_p2_36);}

	// Estas no es necesario darlas vuelta
	void OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2);
	static void jacob_izq_OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0);
	static void jacob_der_OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0); // Para reutilizar JInvP2

	void OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2);
	static void jacob_izq_OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0);
	static void jacob_der_OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0=0, int j0=0);

	L_CoordsCart3D movido_a_fijo(const L_CoordsCart3D &movido) const {return ori.rotarVector(movido) + pos;} // Lo normal  // NOTA: aca pos es la traslacion
	void jacob_izq_movido_a_fijo(const L_CoordsCart3D &movido, L_Matrix &J3x7, int i0=0, int j0=0) const {int i, j; for (i=0; i<3; i++) for (j=0; j<3; j++) J3x7(i+i0,j+j0) = (i==j); ori.jacob_izq_rotarVector(movido, J3x7, i0, j0+3);}
	void jacob_der_movido_a_fijo(L_Matrix &J3x3, int i0=0, int j0=0) const {return ori.jacob_der_rotarVector(J3x3, i0, j0);} // Lo normal

	L_CoordsCart3D fijo_a_movido(const L_CoordsCart3D &fijo) const {L_CoordsCart3D resta = fijo - pos;	return ori.rotarVectorInv(resta);} // RT(p-t) //NOTA: aca pos es la traslacion :)

	L_Pose3D calcPose3D() const;

	void calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Pose3D_cuat *inic=NULL, L_LevenbergMarquardt *opt=NULL); // niter cerca de 30 demora en converger
	void calcPose_movido_a_fijo_Nrayos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial=NULL, L_Array<L_Pose3D_cuat> *candidatos=NULL, int nIntentos3pt = 100);
	void calcPose_movido_a_fijo_Nrayos_sqrt(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial=NULL, L_Array<L_Pose3D_cuat> *candidatos=NULL, int nIntentos3pt = 100); // Un M-estimador de value absoluto
	void calcPose_movido_a_fijo_Npuntos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_CoordsCart3D> &pf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial=NULL, L_Array<L_Pose3D_cuat> *candidatos=NULL, int nIntentos3pt = 100);
	void calcPose_movido_a_fijo_Npuntos_sqrt(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_CoordsCart3D> &pf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial=NULL, L_Array<L_Pose3D_cuat> *candidatos=NULL, int nIntentos3pt = 100); // Un M-estimador de value absoluto

	L_RayoCamara proyectar_movido_a_fijo(const L_CoordsCart3D &movido) const {L_CoordsCart3D fijo; fijo = movido_a_fijo(movido); return L_RayoCamara(fijo.y/fijo.x, fijo.z/fijo.x);}

	// Usado por RANSAC
	bool calcPose_movido_a_fijo_4rayos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, int i0 = 0); // algoritmo de los 3 puntos

private:
	// Funciones directas, trabajan con q, inutiles para el SLAM
	L_Rayo3D movido_a_fijo(const L_Rayo3D &movido) const {L_Rayo3D ret; ret.origen = movido_a_fijo(movido.origen); ret.direcc = movido_a_fijo(movido.origen+movido.direcc)-ret.origen; return ret;} // Lo normal
	void jacob_izq_proyectar_movido_a_fijo(const L_CoordsCart3D &movido, L_Matrix &J2x7, int i0, int j0) const {L_Matrix J2x3s(2,3), J3x7(3,7); L_CoordsCart3D fijo; fijo = movido_a_fijo(movido); jacob_izq_movido_a_fijo(movido, J3x7, 0, 0); L_RayoCamara::jacob_proyeccionDe(fijo, J2x3s, 0, 0); J2x7.OP_mult_sub(i0,j0,J2x3s,J3x7,0,0,2,3,0,0,3,7);}
	void jacob_der_proyectar_movido_a_fijo(const L_CoordsCart3D &movido, L_Matrix &J2x3, int i0, int j0) const {L_Matrix J2x3s(2,3), J3x3(3,3); L_CoordsCart3D fijo; fijo = movido_a_fijo(movido); jacob_der_movido_a_fijo(J3x3, 0, 0); L_RayoCamara::jacob_proyeccionDe(fijo, J2x3s, 0, 0); J2x3.OP_mult_sub(i0,j0,J2x3s,J3x3,0,0,2,3,0,0,3,3);}
	void jacob_proyectar_movido_a_fijo_v(const L_CoordsCart3D &movido, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const {L_Matrix J2x3s(2,3), J3x3(3,3), J3x7(3,7); L_CoordsCart3D fijo; fijo = movido_a_fijo(movido); jacob_izq_movido_a_fijo(movido, J3x7, 0, 0); jacob_der_movido_a_fijo(J3x3, 0, 0); L_RayoCamara::jacob_proyeccionDe(fijo, J2x3s, 0, 0); J2x7.OP_mult_sub(i0,j0,J2x3s,J3x7,0,0,2,3,0,0,3,7); J2x3.OP_mult_sub(i1,j1,J2x3s,J3x3,0,0,2,3,0,0,3,3);} // Error: decia (i1,j1) y (i1,j1)
	L_RayoCamara proyectar_movido_a_fijo(const L_CoordsInv3D &movido) const {L_CoordsCart3D pos; pos = movido.pideCoords3D(); return proyectar_movido_a_fijo(pos);}
	void jacob_izq_proyectar_movido_a_fijo(const L_CoordsInv3D &movido, L_Matrix &J2x7, int i0, int j0) const {L_CoordsCart3D pos; pos = movido.pideCoords3D(); jacob_izq_proyectar_movido_a_fijo(pos, J2x7,i0,j0);}
	void jacob_der_proyectar_movido_a_fijo(const L_CoordsInv3D &movido, L_Matrix &J2x6, int i0, int j0) const {L_Matrix J2x3(2,3), J3x6(3,6); L_CoordsCart3D pos; pos = movido.pideCoords3D(); movido.jacob_pideCoords3D(J3x6,0,0); jacob_der_proyectar_movido_a_fijo(pos, J2x3,0,0);  J2x6.OP_mult_sub(i0,j0,J2x3,J3x6,0,0,2,3,0,0,3,6);}
	/*muy mula*/void calcPose_movido_a_fijo_rapidoMuyImpreciso(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo); // Muy impreciso, funciona para desplazamientos pequenos
	void calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo); // Muy impreciso, funciona para desplazamientos pequenos
	/*muy mula*/void jacob_movido_calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0=0, int j0=0);
	void jacob_movido_calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Matrix &J7x3n, int i0=0, int j0=0, double delta=1.0e-3); // NO FUNCIONA MUY BIEN
	void jacob_movido_calcPose_movido_a_fijo_matriz_diffin(std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0=0, int j0=0); // Movido se modifica para hacerle variaciones
	/*algo mula*/void jacob_fijo_calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0=0, int j0=0);
	void jacob_fijo_calcPose_movido_a_fijo_matriz_diffin(const std::vector<L_CoordsCart3D> &movido, std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0=0, int j0=0); // Se modifica fijo para hacer variaciones
	void jacob_fijo_calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Matrix &J7x3n, int i0=0, int j0=0, double delta=1.0e-3); // NO FUNCIONA MUY BIEN
	bool calcpose_movido_a_fijo_3puntos(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo);
	// rayos movidos, puntos fijos
	bool _calcPose_movido_a_fijo_4rayos(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, int i0 = 0); // algoritmo de los 3 puntos
	bool _calcPose_movido_a_fijo_4rayos3D(const L_Array<L_Rayo3D> &rm, const std::vector<L_CoordsCart3D> &pf, int i0 = 0); // algoritmo de los 3 puntos generalizado
	void _calcPose_movido_a_fijo_Nrayos(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, int niter = 10, L_Pose3D_cuat *inicial=NULL); // usa el algoritmo de 3 puntos para obtener una estimacion inicial, causa confusion
	void _calcPose_movido_a_fijo_Nrayos3D(const L_Array<L_Rayo3D> &rm, const std::vector<L_CoordsCart3D> &pf, int niter = 10, L_Pose3D_cuat *inicial=NULL); // usa el algoritmo de 3 puntos para obtener una estimacion inicial
	/*muy mula*/void _jacob_rm_calcPose_movido_a_fijo_Nrayos_matriz_malo(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, L_Matrix &J7x2n, int i0=0, int j0=0); // Horrible, esta para mostrar que es penca
	void _jacob_rm_calcPose_movido_a_fijo_Nrayos_diffin(std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, L_Matrix &J7x2n, int i0=0, int j0=0); // Se modifica rm para hacerle variaciones
	// rayos fijos, puntos movidos, SIRVEN PARA EL LANDMARK-CONJUNTO
	bool calcPose_movido_a_fijo_4rayos3D(const std::vector<L_CoordsCart3D> &pm, const L_Array<L_Rayo3D> &rf, int i0 = 0) {L_Pose3D_cuat q; if (q._calcPose_movido_a_fijo_4rayos3D(rf, pm, i0)==false) return false; inverseOf(q); return true;} // algoritmo de los 3 puntos generalizado // No implementado
	void calcPose_movido_a_fijo_Nrayos3D(const std::vector<L_CoordsCart3D> &pm, const L_Array<L_Rayo3D> &rf, int niter = 10, L_Pose3D_cuat *inicial=NULL);
	//
	bool jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &Hess7x7inv, double delta = 1.0e-4);
	void jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, const L_Matrix &Hess7x7inv, L_Matrix &J7x2n, int i0=0, int j0=0, double delta = 1.0e-4);
	void jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, const L_Matrix &Hess7x7inv, L_Matrix &J7x3n, int i0=0, int j0=0, double delta = 1.0e-4);
	//
	void jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x2n, int i0=0, int j0=0, L_Pose3D_cuat *inic=NULL, L_Array<L_Pose3D_cuat> *candidatos=NULL, double delta = 1.0e-3, double epsilon = 1.0e-8, double errorToStop = 1.0e-6, int nIterationsMax = 10, int nIteracParar2 = 4); // Error: este calculaba respecto a pm y no rf
	void jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x3n, int i0=0, int j0=0, L_Pose3D_cuat *inic=NULL, double delta = 1.0e-3, double epsilon = 1.0e-8, double errorToStop = 1.0e-6, int nIterationsMax = 10, int nIteracParar2 = 4); //
	void jacob_rf_calcPose_movido_a_fijo_Nrayos_matriz_malo(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x2n, int i0=0, int j0=0); // Requiere *this ya calculado, no funciona muy bien

	// Funciones inversas, trabajan con q^-1
	// (q1*q2).fijo_a_movido(p) = q2.fijo_a_movido(q1.fijo_a_movido(p));
	void jacob_der_fijo_a_movido(L_Matrix &J3x3, int i0=0, int j0=0) const {return ori.jacob_der_rotarVectorInv(J3x3, i0, j0);}//L_CoordsCart3D resta = pos - movido;	return -ori.rotarVectorInv(resta);} // -RT(p-t)
	L_Rayo3D fijo_a_movido(const L_Rayo3D &fijo) const {L_Rayo3D ret; ret.origen = fijo_a_movido(fijo.origen); ret.direcc = fijo_a_movido(fijo.origen+fijo.direcc)-ret.origen; return ret;} // Lo normal
	L_RayoCamara proyectar_fijo_a_movido_v(const L_CoordsCart3D &fijo) const {L_RayoCamara res; L_CoordsCart3D movi_lp3dcpmf; L_Pose3D_cuat_fijo_a_movido(movi_lp3dcpmf,(*this),fijo); res.tanIzq=movi_lp3dcpmf.y/movi_lp3dcpmf.x; res.tanArr=movi_lp3dcpmf.z/movi_lp3dcpmf.x; return res;}
	L_RayoCamara proyectar_fijo_a_movido_vAde(const L_CoordsCart3D &fijo, double &xs) const {L_RayoCamara res; L_CoordsCart3D movi_lp3dcpmf; L_Pose3D_cuat_fijo_a_movido(movi_lp3dcpmf,(*this),fijo); xs=movi_lp3dcpmf.x; res.tanIzq=movi_lp3dcpmf.y/movi_lp3dcpmf.x; res.tanArr=movi_lp3dcpmf.z/movi_lp3dcpmf.x;return res;}
	void jacob_proyectar_fijo_a_movido_pre(L_StaticMatrix<3,3> &J3x3s, double *pre36m) const {L_Cuaternion_jacob_der_rotarVectorInv(ori,J3x3s,0,0); ori.jacob_izq_rotarVectorInv_pre(pre36m); int i; for(i=0; i<36; i++) pre36m[i]=-pre36m[i];}
	void jacob_proyectar_fijo_a_movido_v(const L_CoordsCart3D fijo, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const {L_StaticMatrix<2,3> J2x3s; L_StaticMatrix<3,7> J3x7; L_CoordsCart3D movi_xd; L_Pose3D_cuat_fijo_a_movido(movi_xd,((*this)),fijo); L_Pose3D_cuat_jacob_izq_fijo_a_movido(((*this)),fijo,(pre36m),J3x3s,J3x7,0,0); L_RayoCamara_jacob_proyeccionDe(movi_xd,J2x3s,0,0); L_StaticMatrix_OP_mult_sub((J2x7),(i0),(j0),J2x3s,J3x7,0,0,2,3,0,0,3,7); L_StaticMatrix_OP_mult_sub((J2x3),(i1),(j1),J2x3s,(J3x3s),0,0,2,3,0,0,3,3);}
	L_RayoCamara proyectar_fijo_a_movido_vi(const L_CoordsInv3D &fijo) const {L_RayoCamara res; L_CoordsCart3D pos_lpmfi; L_CoordsInv3D_pideCoords3D(pos_lpmfi,fijo); L_Pose3D_cuat_proyectar_fijo_a_movido_v(res,(*this),pos_lpmfi); return res;}
	L_RayoCamara proyectar_fijo_a_movido_viAde(const L_CoordsInv3D &fijo, double &xs) const {L_RayoCamara res; L_CoordsCart3D pos_lpmfi; L_CoordsInv3D_pideCoords3D(pos_lpmfi,fijo); L_Pose3D_cuat_proyectar_fijo_a_movido_vAde(res,(*this),pos_lpmfi,xs); return res;}
	void jacob_proyectar_fijo_a_movido_vi(const L_CoordsInv3D &fijo, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const {L_StaticMatrix<2,3> J2x3_ljpmfi; L_StaticMatrix<3,6> J3x6_ljpmfi; L_CoordsCart3D pos_ljpmfi; L_CoordsInv3D_pideCoords3D(pos_ljpmfi,(fijo)); L_CoordsInv3D_jacob_pideCoords3D((fijo),J3x6_ljpmfi,0,0); L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v(((*this)),pos_ljpmfi,(J3x3s),(pre36m),(J2x7),(i0),(j0),J2x3_ljpmfi,0,0); L_StaticMatrix_OP_mult_sub((J2x6),(i1),(j1),J2x3_ljpmfi,J3x6_ljpmfi,0,0,2,3,0,0,3,6);}
	// Las funciones que siguen son del tipo (pCM*pC2_C).proyfm(fijo)
	// Requieren *this = (pCM*pC2_C) ya calculado
	/*no funka*/void jacob_proyectar_fijo_a_movido_v_c2_ant(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, const double *pre36_pr, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const;
	/*no funka*/void jacob_proyectar_fijo_a_movido_v_c2(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &prJ3x3s, const L_StaticMatrix<3,3> &pc2J3x3s, const L_StaticMatrix<3,3> &prpc2J3x3s, const double *prpre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const;
	static void jacob_proyectar_fijo_a_movido_v_c2_(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1);
	/*no funka*/void jacob_proyectar_fijo_a_movido_vi_c2_ant(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, const double *pre36_pr, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const;
	/*no funka*/void jacob_proyectar_fijo_a_movido_vi_c2(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &prJ3x3s, const L_StaticMatrix<3,3> &pc2J3x3s, const L_StaticMatrix<3,3> &prpc2J3x3s, const double *prpre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const;
	static void jacob_proyectar_fijo_a_movido_vi_c2_(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1);

	void jacob_movido_calcPose_fijo_a_movido_matriz_v2(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0=0, int j0=0);
	bool calcpose_fijo_a_movido_3puntos(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido) {return calcpose_movido_a_fijo_3puntos(movido, fijo);}

	void calcPose_fijo_a_movido_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo) {return calcPose_movido_a_fijo_matriz(fijo, movido);} // Muy impreciso, funciona para desplazamientos pequenos
	void calcPose_fijo_a_movido_matriz_v2(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo);
	void calcPose_fijo_a_movido_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter) {return calcPose_movido_a_fijo_iter(fijo, movido, niter);}; // niter cerca de 30 demora en converger

	static void err_mf_vect(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	static void err_mf_vect_rapido(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	static void err_mf_vect_rapido_sqrt(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	//
	static void err_proy_fm_vect(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	static void err_proy_fm_vect_rapido(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	//
	static void err_proy_mf_vect(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	static void err_proy_mf_vect_rapido(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido
	static void err_proy_mf_vect_rapido_sqrt(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores); // Fijo a movido

public:

	double cosineDistanceTo(L_Pose3D_cuat &eta) const {return 0.5*pos.cosineDistanceTo(eta.pos) + 0.5*ori.cosineDistanceTo(eta.ori);}

	L_Pose3D_cuat operator*(const L_Pose3D_cuat &eta2) const {L_Pose3D_cuat eta; eta.OP_mult(*this, eta2);return eta;}
	L_Pose3D_cuat operator/(const L_Pose3D_cuat &eta2) const {L_Pose3D_cuat eta; eta.OP_div(*this, eta2);return eta;}
	void pow(const L_Pose3D_cuat &pose, double exp); // Para slerp con poses


	// Notacion usada en la tesis
	void mov(L_CoordsCart3D &res, const L_CoordsCart3D &v) const  {L_Pose3D_cuat_movido_a_fijo(res,(*this),v);}
	void movi(L_CoordsCart3D &res, const L_CoordsCart3D &v) const  {L_Pose3D_cuat_fijo_a_movido(res,(*this),v);}
	void jmovi_pre(L_StaticMatrix<3,3> &J3x3s, double *pre36m) const  {L_Pose3D_cuat_jacob_fijo_a_movido_pre((*this),J3x3s,pre36m);}
	void jmovi(L_CoordsCart3D &v, const L_StaticMatrix<3,3> &J3x3s, const double *prem, L_Matrix &J3x7, int i0, int j0, L_Matrix &J3x3, int i1, int j1) const {L_Pose3D_cuat_jacob_izq_fijo_a_movido((*this),v,prem,J3x3s,J3x7,i0,j0); L_Pose3D_cuat_jacob_der_fijo_a_movido((*this),J3x3,i1,j1);}
//
	void prmov_v(L_RayoCamara &res, const L_CoordsCart3D &v, double &xs) const {L_Pose3D_cuat_proyectar_movido_a_fijo_vAde(res,(*this),v,xs);}
	void jprmov_v(const L_CoordsCart3D &v, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const {jacob_proyectar_movido_a_fijo_v(v,J2x7,i0,j0,J2x3,i1,j1);}

	void prmovi_v(L_RayoCamara &res, const L_CoordsCart3D &v, double &xs) const {L_Pose3D_cuat_proyectar_fijo_a_movido_vAde(res,(*this),v,xs);}
	void jprmovi_pre(L_StaticMatrix<3,3> &J3x3s, double *pre36m) const {jacob_proyectar_fijo_a_movido_pre(J3x3s,pre36m);}
	void jprmovi_v(const L_CoordsCart3D &v, const L_StaticMatrix<3,3> &J3x3s, const double *prem, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const {jacob_proyectar_fijo_a_movido_v(v,J3x3s,prem,J2x7,i0,j0,J2x3,i1,j1);}
	static void jprmovi_c2_v(const L_CoordsCart3D &v, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) {jacob_proyectar_fijo_a_movido_v_c2_(v, pCM, pC2_C, J2x7, i0, j0, J2x3, i1, j1);}
	void prmovi_vi(L_RayoCamara &res, const L_CoordsInv3D &vi, double &xs) const {L_Pose3D_cuat_proyectar_fijo_a_movido_viAde(res,(*this),vi,xs);}
	void jprmovi_vi(const L_CoordsInv3D &vi, const L_StaticMatrix<3,3> &J3x3s, const double *prem, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const {jacob_proyectar_fijo_a_movido_vi(vi,J3x3s,prem,J2x7,i0,j0,J2x6,i1,j1);}
	static void jprmovi_c2_vi(const L_CoordsInv3D &vi, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) {jacob_proyectar_fijo_a_movido_vi_c2_(vi, pCM, pC2_C, J2x7, i0, j0, J2x6, i1, j1);}

	// Las funciones mov() mueven la media, y las funciones movCov() la covarianza
	void movCov(L_StaticMatrix<3,3> &res, const L_StaticMatrix<3,3> &orig);
	void moviCov(L_StaticMatrix<3,3> &res, const L_StaticMatrix<3,3> &orig);

	bool promedioIR(const L_Pose3D_cuat &other, double factor);

	static void pruebaPoseCuaternion(int nPrueba = -1); // Varios , 0:16
	static void pruebaPoseCuaternion_2(); // Funciones de proyeccion y sus jacobianos con fnes y defines
	static void pruebaPoseCuaternion_3(); // Funciones de proyeccion para la segunda camara
	static void pruebaPoseCuaternion_4(); 
	static void pruebaPoseCuaternion_d(); 
	static void pruebaPoseCuaternion_converg();  // Funcion chanta

	static void pruebaPoseCuaternion_veloc();

	friend class L_LandSolido;
	friend class L_LandSolidoKin;
	friend class L_RotEscTrasl3D;

#ifdef __TURBOC__
	__OP_ASSIGN(L_Pose3D_cuat)
#endif
};

class L_Pose3D_cuat_sumable : public L_Pose3D_cuat // sum no distributiva...
{
public:
	L_Pose3D_cuat_sumable() {}
	L_Pose3D_cuat_sumable(double x, double y, double z, double a, double b, double c, double d) {pos.x=x; pos.y=y; pos.z=z; ori.a=a; ori.b=b; ori.c=c; ori.d=d;}
	void OP_add(const L_Pose3D_cuat_sumable &a, const L_Pose3D_cuat_sumable &b) {pos.x = a.pos.x+b.pos.x; pos.y = a.pos.y+b.pos.y; pos.z = a.pos.z+b.pos.z; ori.a = a.ori.a+b.ori.a; ori.b = a.ori.b+b.ori.b; ori.c = a.ori.c+b.ori.c; ori.d = a.ori.d+b.ori.d;}
	void OP_subtract(const L_Pose3D_cuat_sumable &a, const L_Pose3D_cuat_sumable &b) {pos.x = a.pos.x-b.pos.x; pos.y = a.pos.y-b.pos.y; pos.z = a.pos.z-b.pos.z; ori.a = a.ori.a-b.ori.a; ori.b = a.ori.b-b.ori.b; ori.c = a.ori.c-b.ori.c; ori.d = a.ori.d-b.ori.d;}
	void OP_mult(const L_Pose3D_cuat_sumable &a, double fac) {pos.x = a.pos.x*fac; pos.y = a.pos.y*fac; pos.z = a.pos.z*fac; ori.a = a.ori.a*fac; ori.b = a.ori.b*fac; ori.c = a.ori.c*fac; ori.d = a.ori.d*fac;}
	void OP_div(const L_Pose3D_cuat_sumable &a, double fac) {pos.x = a.pos.x/fac; pos.y = a.pos.y/fac; pos.z = a.pos.z/fac; ori.a = a.ori.a/fac; ori.b = a.ori.b/fac; ori.c = a.ori.c/fac; ori.d = a.ori.d/fac;}
	L_Pose3D_cuat_sumable operator+(const L_Pose3D_cuat_sumable &eta2) const {L_Pose3D_cuat_sumable eta; eta.OP_add(*this, eta2);return eta;}
	L_Pose3D_cuat_sumable operator-(const L_Pose3D_cuat_sumable &eta2) const {L_Pose3D_cuat_sumable eta; eta.OP_subtract(*this, eta2);return eta;}
	L_Pose3D_cuat_sumable operator*(double fac) const {L_Pose3D_cuat_sumable eta; eta.OP_mult(*this, fac);return eta;}
	L_Pose3D_cuat_sumable operator/(double fac) const {L_Pose3D_cuat_sumable eta; eta.OP_div(*this, fac);return eta;}
	L_Pose3D_cuat_sumable operator*(const L_Pose3D_cuat_sumable &eta2) const {L_Pose3D_cuat_sumable eta; eta.L_Pose3D_cuat::OP_mult(*this, eta2);return eta;}
	L_Pose3D_cuat_sumable operator/(const L_Pose3D_cuat_sumable &eta2) const {L_Pose3D_cuat_sumable eta; eta.L_Pose3D_cuat::OP_div(*this, eta2);return eta;}
};
inline L_Pose3D_cuat_sumable operator*(double fac, const L_Pose3D_cuat_sumable &a) {L_Pose3D_cuat_sumable eta; eta.OP_mult(a, fac); return eta;}

class L_Pose3DSpline : public L_Spline1DGen<L_Pose3D_cuat_sumable>
{
public:
	bool retroceder_t;
	double t1, t2;
	L_Pose3DSpline() : retroceder_t(false) {}
	void agregaPose(double t, const L_Pose3D_cuat &offset, double x, double y, double z, double pan, double tilt, double roll);
	void agregaPoseMirarA(double t, const L_Pose3D_cuat &offset, double x, double y, double z, double xMirado, double yMirado, double zMirado, double roll=0);
	void fijarSaltoTiempoPeriodico(double desde, double hasta) {throw_L_ArgException_if(hasta > desde, "L_Pose3DSpline::fijarSaltoTiempoPeriodico() : salto hacia adelante en el tiempo\n");retroceder_t = true; this->t1 = hasta; this->t2 = desde;}

	bool leerPuntosControl(const char *nomArch);
	bool grabarPuntosControl(const char *nomArch);
	void generarRutaEjemplo(int numRuta, L_Rand &randF); // Semilla negativa inicialmente
	bool generarGrabarRutasAzar(L_StringWithNumber &nomArch, int nRutas, L_Rand &randF);

	void genDibujoFlaite(L_ImageRGBUchar &im);

	L_Pose3D_cuat evaluaPose(double t) {L_Pose3D_cuat_sumable ps; if (retroceder_t && t>t2) t=t1+fmod(t-t1, t2-t1); ps=evaluate(t); ps.normaliza(); return ps;}
};

extern double L_P3DCB_S; // en cm supuestamente

class L_Pose3D_esf
{
public:
	L_CoordsEsf3D pos;
	L_PanTiltRoll ori;
	void fijaPose3D_esf(double r, double theta, double phi, double pan, double tilt, double roll) {pos.fijaRThetaPhi(r, theta, phi); ori.pan=pan; ori.tilt=tilt; ori.roll=roll;}
	void fijaPose3D(const L_Pose3D &pose) {pos.fijaCoordsCart3D(pose.pos); ori=pose.ori;}
	L_Pose3D pidePose3D() {L_Pose3D p; p.pos=pos.pideCoordsCart3D(); p.ori=ori; return p;}
};

typedef L_Node<L_Pose3D> L_Pose3DNodo;
typedef L_List<L_Pose3D> L_Pose3DLista;

// R(t) = [cos(t) -sin(t) ; sin(t) cos(t)]
class L_Rot3DMatrix : public L_StaticMatrix<3,3>
{
public:
	L_Rot3DMatrix &fijarMatrizRotacion(L_Matrix &rot) {int i; for (i=0; i<9; i++) operator()(i/3,i%3)=rot(i/3,i%3); return *this;}
	//
	L_Rot3DMatrix &fijarCuaternion(const L_Quaternion &q);
	void jacob_fijarCuaternion(const L_Quaternion &q, L_Matrix &J4x9, int i0=0, int j0=0);
	L_Rot3DMatrix &fijarCuaternionInv(const L_Quaternion &q);
	void jacob_fijarCuaternionInv(const L_Quaternion &q, L_Matrix &J4x9, int i0=0, int j0=0);
	//
	L_Rot3DMatrix &fijarVectorRotacion(const L_CoordsCart3D &eje);
	L_Rot3DMatrix &fijarPanTiltRoll(const L_PanTiltRoll &ori);
	//
	L_Quaternion pideCuaternion();
	L_Quaternion jacob_pideCuaternion(L_Matrix &J4x9, int i0=0, int j0=0);
	L_Quaternion pideCuaternionInv();
	L_Quaternion jacob_pideCuaternionInv(L_Matrix &J4x9, int i0=0, int j0=0);
	//
	L_CoordsCart3D pideVectorRotacion();
	L_PanTiltRoll pidePanTiltRoll();
	void transpOf(const L_Rot3DMatrix &other);

	void invertMeTrasponiendo() {L_StaticMatrix<3,3> other; L_StaticMatrix_OP_assign(other,*this); L_StaticMatrix_transpDe(*this, other);}

	double normalizaSVD(double *eMax=NULL, double *eMed=NULL, double *eMin=NULL);
	double normalizaSVD_corrigeJ(L_Matrix &J9xN);

	void OP_mult(const L_Rot3DMatrix &m1, const L_Rot3DMatrix &m2);
	void OP_mult(const L_Matrix &m1, const L_Matrix &m2);
	L_CoordsCart3D OP_mult(const L_CoordsCart3D &pun) {return L_CoordsCart3D(operator()(0,0)*pun.x+operator()(0,1)*pun.y+operator()(0,2)*pun.z,operator()(1,0)*pun.x+operator()(1,1)*pun.y+operator()(1,2)*pun.z,operator()(2,0)*pun.x+operator()(2,1)*pun.y+operator()(2,2)*pun.z);}
};

// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_MatrizHomogenea_reconstruye(M) {if ((M).size() == 0) (M).L_Matrix::reallocate(4, 4); (M)(3,0)=0; (M)(3,1)=0; (M)(3,2)=0; (M)(3,3)=1;}
#define L_MatrizHomogenea_fijaCero(M) {L_MatrizHomogenea_reconstruye((M)); int i_32, j_32; for (i_32=0; i_32<3; i_32++) for (j_32=0; j_32<4; j_32++) (M)(i_32,j_32)=0;}
#define L_MatrizHomogenea_fijaIdent(M) {L_MatrizHomogenea_reconstruye((M)); int i_33, j_33; for (i_33=0; i_33<3; i_33++) for (j_33=0; j_33<4; j_33++) (M)(i_33,j_33)=(i_33==j_33);}

#define L_MatrizHomogenea_OP_mult_M(mRet,M1,M2) {int iOMM,jOMM,uOMM; for (iOMM=0; iOMM<3; iOMM++) { for (jOMM=0; jOMM<4; jOMM++) { (mRet)(iOMM,jOMM)=0; for (uOMM=0; uOMM<4; uOMM++) (mRet)(iOMM,jOMM)+=(M1)(iOMM,uOMM)*(M2)(uOMM,jOMM); } } (mRet)(3,0) = (mRet)(3,1) = (mRet)(3,2) = 0.0; (mRet)(3,3) = 1.0;}
#define L_MatrizHomogenea_OP_mult_P(pRet,M,pOtro) ((pRet).x=(pOtro).x*(M)(0,0)+(pOtro).y*(M)(0,1)+(pOtro).z*(M)(0,2)+(M)(0,3), (pRet).y=(pOtro).x*(M)(1,0)+(pOtro).y*(M)(1,1)+(pOtro).z*(M)(1,2)+(M)(1,3), (pRet).z=(pOtro).x*(M)(2,0)+(pOtro).y*(M)(2,1)+(pOtro).z*(M)(2,2)+(M)(2,3))
#define L_MatrizHomogenea_OP_rota(pRet,M,pOtro) ((pRet).x=(pOtro).x*(M)(0,0)+(pOtro).y*(M)(0,1)+(pOtro).z*(M)(0,2), (pRet).y=(pOtro).x*(M)(1,0)+(pOtro).y*(M)(1,1)+(pOtro).z*(M)(1,2), (pRet).z=(pOtro).x*(M)(2,0)+(pOtro).y*(M)(2,1)+(pOtro).z*(M)(2,2))
#define L_MatrizHomogenea_revisa(M) ((M)(3,0) > -0.1 && (M)(3,0) < 0.1 && (M)(3,1) > -0.1 && (M)(3,1) < 0.1 && (M)(3,2) > -0.1 && (M)(3,2) < 0.1 && (M)(3,3) > 0.9 && (M)(3,3) < 1.1)

#define L_MatrizHomogenea_fijaPose3D_cuat(res,poseMovido) {L_Pose3D_cuat pm_lmfp; double a2_lmfp, b2_lmfp, c2_lmfp, d2_lmfp, ab_lmfp, ac_lmfp, ad_lmfp, bc_lmfp, bd_lmfp, cd_lmfp; L_Pose3D_cuat_OP_igual(pm_lmfp,(poseMovido)); L_Pose3D_cuat_normaliza(pm_lmfp); a2_lmfp = pm_lmfp.ori.a*pm_lmfp.ori.a; b2_lmfp = pm_lmfp.ori.b*pm_lmfp.ori.b; c2_lmfp = pm_lmfp.ori.c*pm_lmfp.ori.c; d2_lmfp = pm_lmfp.ori.d*pm_lmfp.ori.d;\
	ab_lmfp = pm_lmfp.ori.a*pm_lmfp.ori.b; ac_lmfp = pm_lmfp.ori.a*pm_lmfp.ori.c; ad_lmfp = pm_lmfp.ori.a*pm_lmfp.ori.d; bc_lmfp = pm_lmfp.ori.b*pm_lmfp.ori.c; bd_lmfp = pm_lmfp.ori.b*pm_lmfp.ori.d; cd_lmfp = pm_lmfp.ori.c*pm_lmfp.ori.d; (res)(0,0) = a2_lmfp+b2_lmfp-c2_lmfp-d2_lmfp; (res)(0,1) = 2*bc_lmfp-2*ad_lmfp; (res)(0,2) = 2*bd_lmfp+2*ac_lmfp;  (res)(0,3)=pm_lmfp.pos.x;\
	(res)(1,0) = 2*bc_lmfp+2*ad_lmfp; (res)(1,1) = a2_lmfp-b2_lmfp+c2_lmfp-d2_lmfp; (res)(1,2) = 2*cd_lmfp - 2*ab_lmfp; (res)(1,3)=pm_lmfp.pos.y; (res)(2,0) = 2*bd_lmfp - 2*ac_lmfp; (res)(2,1) = 2*cd_lmfp + 2*ab_lmfp; (res)(2,2) = a2_lmfp-b2_lmfp-c2_lmfp+d2_lmfp; (res)(2,3)=pm_lmfp.pos.z; (res)(3,0)=0; (res)(3,1)=0; (res)(3,2)=0; (res)(3,3)=1;}

#define L_MatrizHomogenea_fijaPose3D_cuat_fijo_a_movido(res,poseMovido) {L_Pose3D_cuat pm_lmfpi; double a2_lmfpi, b2_lmfpi, c2_lmfpi, d2_lmfpi, ab_lmfpi, ac_lmfpi, ad_lmfpi, bc_lmfpi, bd_lmfpi, cd_lmfpi; L_Pose3D_cuat_OP_igual(pm_lmfpi,(poseMovido)); L_Pose3D_cuat_normaliza(pm_lmfpi); a2_lmfpi = pm_lmfpi.ori.a*pm_lmfpi.ori.a; b2_lmfpi = pm_lmfpi.ori.b*pm_lmfpi.ori.b; c2_lmfpi = pm_lmfpi.ori.c*pm_lmfpi.ori.c;\
	d2_lmfpi = pm_lmfpi.ori.d*pm_lmfpi.ori.d; ab_lmfpi = pm_lmfpi.ori.a*pm_lmfpi.ori.b; ac_lmfpi = pm_lmfpi.ori.a*pm_lmfpi.ori.c; ad_lmfpi = pm_lmfpi.ori.a*pm_lmfpi.ori.d; bc_lmfpi = pm_lmfpi.ori.b*pm_lmfpi.ori.c; bd_lmfpi = pm_lmfpi.ori.b*pm_lmfpi.ori.d; cd_lmfpi = pm_lmfpi.ori.c*pm_lmfpi.ori.d; (res)(0,0) = a2_lmfpi+b2_lmfpi-c2_lmfpi-d2_lmfpi; (res)(0,1) = 2*bc_lmfpi+2*ad_lmfpi; (res)(0,2) = 2*bd_lmfpi - 2*ac_lmfpi;\
	(res)(1,0) = 2*bc_lmfpi-2*ad_lmfpi; (res)(1,1) = a2_lmfpi-b2_lmfpi+c2_lmfpi-d2_lmfpi; (res)(1,2) = 2*cd_lmfpi + 2*ab_lmfpi; (res)(2,0) = 2*bd_lmfpi+2*ac_lmfpi; (res)(2,1) = 2*cd_lmfpi - 2*ab_lmfpi; (res)(2,2) = a2_lmfpi-b2_lmfpi-c2_lmfpi+d2_lmfpi; (res)(0,3) = -(res)(0,0)*(poseMovido).pos.x - (res)(0,1)*(poseMovido).pos.y - (res)(0,2)*(poseMovido).pos.z;\
	(res)(1,3) = -(res)(1,0)*(poseMovido).pos.x - (res)(1,1)*(poseMovido).pos.y - (res)(1,2)*(poseMovido).pos.z; (res)(2,3) = -(res)(2,0)*(poseMovido).pos.x - (res)(2,1)*(poseMovido).pos.y - (res)(2,2)*(poseMovido).pos.z; res(3,0)=0; res(3,1)=0; res(3,2)=0; res(3,3)=1;}


// La matriz de rotacion theta es: [cos(theta) -sin(theta); sin(theta) cos(theta)]
class L_HomogeneousMatrix : public L_StaticMatrix<4,4> // Matriz homogenea de 4x4 exclusivamente. Siempre la ultima fila es [0 0 0 1]
{                                         // Considera coordenadas CoordsCart3D del tipo AIA (x adelante, y izq, z arr)
public:
	L_HomogeneousMatrix() {operator()(3,0)=0; operator()(3,1)=0; operator()(3,2)=0; operator()(3,3)=1;}
	bool esValida() {return rows()==4 && cols()==4 && operator()(3,0)==0 && operator()(3,1)==0 && operator()(3,2)==0 && operator()(3,3)==1;} // Se debe cumplir siempre para este object

	// Funciones directas, trabajan con H
	L_CoordsCart3D operator *(const L_CoordsCart3D &other); // coords_fijo = matriz * coords_movido (apply el desplaz)
	L_Rayo3D operator *(const L_Rayo3D &other);
	L_Rayo3D operator *(const L_RayoCamara &other) {return (*this)*L_Rayo3D().OP_assign(other);}
	void OP_mult(const L_HomogeneousMatrix &m1, const L_HomogeneousMatrix &m2);
	L_CoordsCart3D movido_a_fijo(const L_CoordsCart3D &other);
	L_CoordsCart3D fijo_a_movido(const L_CoordsCart3D &other);
	L_CoordsCart3D OP_rota(const L_CoordsCart3D &other); // solo rotacion
	L_HomogeneousMatrix & operator*=(const L_HomogeneousMatrix &other);
	void fijaCeroTodo() {operator()(0,0)=0; operator()(0,1)=0;operator()(0,2)=0; operator()(0,3)=0; operator()(1,0)=0; operator()(1,1)=0; operator()(1,2)=0; operator()(1,3)=0; operator()(2,0)=0; operator()(2,1)=0;operator()(2,2)=0; operator()(2,3)=0; operator()(3,0)=0; operator()(3,1)=0; operator()(3,2)=0; operator()(3,3)=0;}
	void fijaIdent() {operator()(0,0)=1; operator()(0,1)=0;operator()(0,2)=0; operator()(0,3)=0; operator()(1,0)=0; operator()(1,1)=1; operator()(1,2)=0; operator()(1,3)=0; operator()(2,0)=0; operator()(2,1)=0;operator()(2,2)=1; operator()(2,3)=0; operator()(3,0)=0; operator()(3,1)=0; operator()(3,2)=0; operator()(3,3)=1;}

	void izqArrAde_2_adeIzqArr(L_HomogeneousMatrix &other); // Cambia a coordenadas (ade,izq,arr) (las normales)
	void adeIzqArr_2_izqArrAde(L_HomogeneousMatrix &other); // Cambia a coordenadas (izq,arr,ade) (para algunos papers)

	void inversaTrasponiendoDe(const L_HomogeneousMatrix& other); // Requiere rotacion de "other" normalizada
	void invertMeTrasponiendo();
	static void jacob_inversaTrasponiendoDe(const L_HomogeneousMatrix& other, L_Matrix &J12x12); // Aproximacion mala para la inversa, requiere rotacion de "other" normalizada
	static void jacob_inversaAproxDe(const L_HomogeneousMatrix& other, L_Matrix &J12x12); // Aproximacion mala para la inversa (0.3)

	// Funciones directas, trabajan con H
	void fijaTrX(double dist) {fijaIdent(); operator()(0,3)=dist;} // transformUsing coords de movido a fijo
	void fijaTrY(double dist) {fijaIdent(); operator()(1,3)=dist;} // transformUsing coords de movido a fijo
	void fijaTrZ(double dist) {fijaIdent(); operator()(2,3)=dist;} // transformUsing coords de movido a fijo
	void fijaTr(double distX, double distY, double distZ) {fijaIdent(); operator()(0,3)=distX; operator()(1,3)=distY; operator()(2,3)=distZ;} // transformUsing coords de movido a fijo

	void fijaRotEjeX(double theta) {fijaIdent(); double c=cos(theta), s=sin(theta); operator()(1,1)=c; operator()(1,2)=-s; operator()(2,1)=s; operator()(2,2)=c;} // transformUsing coords de movido a fijo
	void fijaRotEjeY(double theta) {fijaIdent(); double c=cos(theta), s=sin(theta); operator()(2,2)=c; operator()(2,0)=-s; operator()(0,2)=s; operator()(0,0)=c;} // transformUsing coords de movido a fijo
	void fijaRotEjeZ(double theta) {fijaIdent(); double c=cos(theta), s=sin(theta); operator()(0,0)=c; operator()(0,1)=-s; operator()(1,0)=s; operator()(1,1)=c;} // transformUsing coords de moviso a fijo

	void setRotationComponents(L_Matrix &matrRot) {int i, j; for (i=0; i<3; i++) for (j=0; j<3; j++) operator()(i,j)=matrRot(i,j);} // set componentes de rotacion
	void fijaTraslacion(L_Matrix &vectTras) {operator()(0,3)=vectTras(0,0); operator()(1,3)=vectTras(1,0); operator()(2,3)=vectTras(2,0);} // set componentes de traslacion en la cuarta columna
	void getRotationComponents(L_Matrix &R) {R.reallocate(3,3); R(0,0)=operator()(0,0); R(0,1)=operator()(0,1); R(0,2)=operator()(0,2); R(1,0)=operator()(1,0); R(1,1)=operator()(1,1); R(1,2)=operator()(1,2); R(2,0)=operator()(2,0); R(2,1)=operator()(2,1); R(2,2)=operator()(2,2);}
	void getTranslationComponents(L_Matrix &t) {t.reallocate(3,1); t(0,0)=operator()(0,3); t(1,0)=operator()(1,3); t(2,0)=operator()(2,3);} // Pide domponentes de traslacion en la cuarta columna
	void fija4columnas(L_CoordsCart3D &r1, L_CoordsCart3D &r2, L_CoordsCart3D &r3, L_CoordsCart3D &t) {operator()(0,0) = r1.x; operator()(0,1) = r2.x; operator()(0,2) = r3.x; operator()(0,3) = t.x; operator()(1,0) = r1.y; operator()(1,1) = r2.y; operator()(1,2) = r3.y; operator()(1,3) = t.y; operator()(2,0) = r1.z; operator()(2,1) = r2.z; operator()(2,2) = r3.z; operator()(2,3) = t.z;}
	void pide4columnas(L_CoordsCart3D &r1, L_CoordsCart3D &r2, L_CoordsCart3D &r3, L_CoordsCart3D &t) {r1.x=operator()(0,0); r2.x=operator()(0,1); r3.x=operator()(0,2); t.x=operator()(0,3); r1.y=operator()(1,0); r2.y=operator()(1,1); r3.y=operator()(1,2); t.y=operator()(1,3); r1.z=operator()(2,0); r2.z=operator()(2,1); r3.z=operator()(2,2); t.z=operator()(2,3);}

	// Funciones directas, trabajan con H
	void fijaPose3D(const L_Pose3D &poseMovido); // Uso directo, al usar esto queda coordsFijo=H*coordsMovido
	void fijaPose3D_lento(const L_Pose3D &poseMovido); // La definicion formal de la transformacion
	L_Pose3D calcPose3D() const; // calcPose3D() y fijaTrPose() deben ser compatibles. Devuelve la pose de movido.
	void calcEsencial(L_Matrix &E) const;
	bool fijaPose3D_cuat(const L_Pose3D_cuat &poseMovido);
	L_Pose3D_cuat calcPose3D_cuat() const;
	L_Pose3D_cuat jacob_calcPose3D_cuat(L_Matrix &J7x12, int i0=0, int j0=0);

	// Funciones inversas, trabajan con H^-1
	void fijaPose3D_fijo_a_movido(const L_Pose3D &poseMovido); // Uso inverso, al usar esto queda coordsMovido=H*coordsFijo
	void fijaPose3D_fijo_a_movido_lento(const L_Pose3D &poseMovido); // La definicion formal de la transformacion inversa
	L_Pose3D calcPose3D_fijo_a_movido() const; // calcPose3D() y fijaTrPose() deben ser compatibles. Devuelve la pose de movido.
	void fijaPose3D_cuat_fijo_a_movido(const L_Pose3D_cuat &poseMovido);
	L_Pose3D_cuat calcPose3D_cuat_fijo_a_movido() const;
	void jacob_der_fijaPose3D_cuat_fijo_a_movido_lento(L_Pose3D_cuat &poseMovido, L_Array<L_Matrix> &J); // Calcula la transformacion y el jacobiano. El jacobiano tiene 4x4x7 componentes

	void fijaRotacionOpenGL() {operator()(0,0) =  0; operator()(0,1) = -1; operator()(0,2) =  0; operator()(1,0) =  0;  operator()(1,1) =  0;  operator()(1,2) =  1; operator()(2,0) = -1;  operator()(2,1) =  0;  operator()(2,2) =  0; operator()(0,3) = 0; operator()(1,3) = 0; operator()(2,3) = 0;}

	void guardaInfoEn(L_Matrix &vector);
	void leeInfoDe(const L_Matrix &vector);
	double normalizaSVD(double *eMax=NULL, double *eMed=NULL, double *eMin=NULL); // Deja la rotacion como 3 vectores ortonormales. Devuelve un error de normalizacion
	double normalizaSVD_corrigeJ(L_Matrix &J12xN); // Deja la rotacion como 3 vectores ortonormales. Devuelve un error de normalizacion
	double normaliza_mincuad(); // setValueElementWise la matriz de rotacion de forma iterativa
	void normalizaTraslacion() {double e=sqrt(operator()(0,3)*operator()(0,3)+operator()(1,3)*operator()(1,3)+operator()(2,3)*operator()(2,3)); operator()(0,3)/=e; operator()(1,3)/=e; operator()(2,3)/=e;}
	double normaL2rot() {return operator()(0,0)*operator()(0,0) + operator()(0,1)*operator()(0,1) + operator()(0,2)*operator()(0,2) + operator()(1,0)*operator()(1,0) + operator()(1,1)*operator()(1,1) + operator()(1,2)*operator()(1,2) + operator()(2,0)*operator()(2,0) + operator()(2,1)*operator()(2,1) + operator()(2,2)*operator()(2,2);}

	void errorRot_prepara(void *matriz, std::vector<double> &panTiltRoll, bool prenormalizar); // La memoria de matriz y panTiltRoll debe estar ya pedida
	static double errorRot(const void *matriz, double *panTiltRoll);

	void copiarLargoTraslacionEn(L_HomogeneousMatrix &other);
	double pedirLargoTraslacion() {return sqrt(operator()(0,3)*operator()(0,3) + operator()(1,3)*operator()(1,3) + operator()(2,3)*operator()(2,3));}
	void fijarLargoTraslacion(double length) {double lant = pedirLargoTraslacion(); operator()(0,3)*=length/lant; operator()(1,3)*=length/lant; operator()(2,3)*=length/lant;}

	bool fijaOrigenAdeArr(const L_CoordsCart3D &ori, const L_CoordsCart3D &dirAde, const L_CoordsCart3D &dirArr);

	void cambiaEjes_a_AdeIzqArr(); // Convencion de ejes en robotica
	void cambiaEjes_a_DerAbaAde(); // Convension de ejes en vision
	void cambiaEjes_transf(const L_HomogeneousMatrix &H, L_SistemaEjes3D inicial, L_SistemaEjes3D final); // Cuando se cambia xyz=H*xyz -> uvw=H*uvw (las coordenadas a ambos lados)
	void cambiaEjes_pose(const L_HomogeneousMatrix &H, L_SistemaEjes3D inicial, L_SistemaEjes3D final); // Cuando se cambia xyz=H*xyz -> uvw=H*xyz (solo las coordenadas al lado izquierdo)

	static L_HomogeneousMatrix rotacionX(double angulo) { L_HomogeneousMatrix m; m.fijaRotEjeX(angulo); return m;}
	static L_HomogeneousMatrix rotacionY(double angulo) { L_HomogeneousMatrix m; m.fijaRotEjeY(angulo); return m;}
	static L_HomogeneousMatrix rotacionZ(double angulo) { L_HomogeneousMatrix m; m.fijaRotEjeZ(angulo); return m;}
	static L_HomogeneousMatrix traslacion(double x, double y, double z) { L_HomogeneousMatrix m; m.fijaTr(x,y,z); return m;}

	static void test();
};

inline L_HomogeneousMatrix operator*(const L_HomogeneousMatrix &uno, const L_HomogeneousMatrix &other) {L_HomogeneousMatrix a; a.OP_mult(uno,other); return a;}


class L_Ellipsoid3D
{
public:
	L_HomogeneousMatrix H;  // global (x) a elipsoide (u): u = H*x  (fijo a movido)
	double sx, sy, sz; // elipsoide (u) a circulo (v): v = S*u  (fijo a movido)
	double rx() {return 1/sx;}
	double ry() {return 1/sy;}
	double rz() {return 1/sz;}
	double mahalanobis(const L_CoordsCart3D &p) const {double dx = H(0,0)*p.x+H(0,1)*p.y+H(0,2)*p.z+H(0,3); double dy = H(1,0)*p.x+H(1,1)*p.y+H(1,2)*p.z+H(1,3); double dz = H(2,0)*p.x+H(2,1)*p.y+H(2,2)*p.z+H(2,3); return sx*sx*dx*dx+sy*sy*dy*dy+sz*sz*dz*dz;}
	void setCovarianceAndTranslation(const double *cov3x3, const double *t3x1, double chi2 = 12.838);
};



class L_Pose3DGrados // Para debug
{
public:
	double x, y, z, panG, tiltG, rollG;
	L_Pose3DGrados() {x=0; y=0; z=0; panG=0; tiltG=0; rollG=0;}
	L_Pose3DGrados(const L_Pose3D &q) {x=q.pos.x; y=q.pos.y; z=q.pos.z; panG=q.ori.pan*180/M_PI; tiltG=q.ori.tilt*180/M_PI; rollG=q.ori.roll*180/M_PI;}
	L_Pose3DGrados(const L_Pose3D_cuat &q) {*this = q;}
	L_Pose3DGrados& operator= (const L_Pose3D &q) {x=q.pos.x; y=q.pos.y; z=q.pos.z; panG=q.ori.pan*180/M_PI; tiltG=q.ori.tilt*180/M_PI; rollG=q.ori.roll*180/M_PI; return *this;}
	L_Pose3DGrados& operator= (const L_Pose3D_cuat &q) {L_HomogeneousMatrix H; H.fijaPose3D_cuat(q); *this = H.calcPose3D(); return *this;}
	void print(const char *name, FILE *fp=stdout) {fprintf(fp, "%s = [%.6g %.6g %.6g pan=%.6g[gr] tilt=%.6g[gr] roll=%.6g[gr]]\n\n", name, x, y, z, panG, tiltG, rollG);}
};

class L_ConjuntoRectas3D  // Modelo de un poligono para ajustar a los puntos
{
public:
	std::vector<L_CoordsCart3D> p; // Puntos que se desea ajustar al poligono
	std::vector<L_CoordsCart3D> vertices;
	std::vector<int> segmentos; // Datos para puntos con asociacion conocida, indica segmentos de los puntos que pertenecen a cada recta
	//
	std::vector<double> a, b, c; // Parametros de las rectas del poligono (contenidas en z=0): ax+by+c=0;
	std::vector<L_CoordsCart3D> pSegm; // Arreglo de puntos con asociacion conocida, calculado internamente
	std::vector<int> nSegm; // Arreglo de indexes de numero de recta para cada punto con asociacion conocida, calculado internamente

	void agregarVertices(double x1, double y1, double x2, double y2);
	void agregarElipseVertices(double lx, double ly, int nVertices, double theta0, double sentido);
	void agregarSegmento(int frame1, int frame2) {segmentos.resize(segmentos.size()+2); segmentos[segmentos.size()-2] = frame1; segmentos[segmentos.size()-1] = frame2;}
	void generarAsociacionesSegmentos();
	double calcularMejorRectaPorPunto(std::vector<L_CoordsCart3D> &puntosMov, std::vector<int> &mejorRecta);
	void resetear() {a.resize(0); b.resize(0); c.resize(0); p.resize(0); segmentos.resize(0); pSegm.resize(0); nSegm.resize(0);}
	void eliminar() {a.clear(); b.clear(); c.clear(); p.clear(); segmentos.clear(); pSegm.clear(); nSegm.clear();}

	void calcularErrorGroundTruth(const char *nomArchM);
	void listaVideos(int numeroVideo);

	static void calcError1(const void *obj, const L_Matrix &x, L_Matrix &err); // Para minimizacion inicial usando puntos con asociacion conocida
	static void calcError2(const void *obj, const L_Matrix &x, L_Matrix &err); // Para minimizacion final usando todos los puntos
};

// Requiere al menos 3 puntos para estimarse
class L_RotEscTrasl3D
{
public:
	// transf(pun) = esc*q*pun*q^-1 + t
	L_Pose3D_cuat pcuat;
	double escala;

	double &el(int i) {throw_L_ArgException_if(i<0 || i>7, "L_RotEscTrasl3D::el()");if (i<7) return pcuat.el(i); else return escala;}
	const double &el(int i) const {throw_L_ArgException_if(i<0 || i>7, "L_RotEscTrasl3D::el()");if (i<7) return pcuat.el(i); else return escala;}

	void fijaCero() {pcuat.fijaCero(); escala = 1;}

	L_CoordsCart3D transformar(const L_CoordsCart3D &orig) const {L_CoordsCart3D res = orig; res = pcuat.ori.rotarVector(res); res*=escala; res+=pcuat.pos; return res;}
	L_CoordsCart3D transformarInv(const L_CoordsCart3D &orig) const {L_CoordsCart3D res = orig; res -= pcuat.pos; res/=escala; res = pcuat.ori.rotarVectorInv(res); return res;}
	static void err_mf_vect(const void *pose3D_cuat_nube_pesos, const L_Matrix &pcuatesc, L_Matrix &errores);
	void encontrarTransfMinCuad(const std::vector<L_CoordsCart3D> &orig, const std::vector<L_CoordsCart3D> &res);
	double encontrarTransfMinCuadTrayectorias(const std::vector<L_CoordsCart3D> &orig, const std::vector<L_CoordsCart3D> &res);
	void calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo);

	void pideMatriz4x4(L_Matrix &m);

	void calcTransfPoligono_1(const L_ConjuntoRectas3D &rectas, int niter = 10);
	void calcTransfPoligono_2(const L_ConjuntoRectas3D &rectas, int niter = 10);
};




class L_CalculadorPoseTriangulo
{
public:
	// Entradas
	L_RayoCamara q1, q2, q3; // input: Direcciones en el espacio de los 3 puntos p1, p2 y p3
	double a, b, c; // input: Distancia entre p2-p3, p1-p3 y p1-p2
	// Salidas
	L_CoordsCart3D j1, j2, j3; // Vectores unitarios que apuntan a los puntos p1, p2 y p3
	double s1[4], s2[4], s3[4]; // output: distancias del foco de la cámara a los puntos p1, p2 y p3, hay 4 posibles triangulos
	bool valida[4]; // Indica cuales resultados son validos, usualmente hay 2 o 4 validos

	bool metodoFinsterwalder(); // Devuelve 4 triangulos que tienen lados a, b, c y que pasan por q1, q2 y q3, NO FUNKA
	bool metodoDirecto(); // Devuelve 4 triangulos que tienen lados a, b, c y que pasan por q1, q2 y q3. Si funka.
	bool pruebame();
};

class L_CalculadorPoseTrianguloGeneralizado
{
public:
	// Entradas
	L_Rayo3D r1,r2,r3;
	L_CoordsCart3D q1,q2,q3;
	L_HomogeneousMatrix H[8];

	bool metodoNister();
};


enum L_CamaraTipo {cam_indef, camPinhole, camRadTang, camDistRad, camDistRadAj};

class L_Camara // Modelo generico indeterminado para cualquier camara
{
public:
	L_CamaraTipo tipo;
private:
	L_Camara() : tipo(cam_indef) {}
public:
	L_Camara(L_CamaraTipo tipoCam) : tipo(tipoCam) {}

	virtual void calcRayo(L_RayoCamara &rayo, double xPix, double yPix) const = 0;
	virtual void calcPixel(const L_RayoCamara &rayo, double &xPix, double &yPix) const = 0;
	virtual double radioErrorTan() const = 0;
	virtual void calcRayos(L_Matrix &uv_nx2) const; // Sobreescribe el resultado en la misma matriz
	virtual void calcPixeles(L_Matrix &uv_nx2) const; // Sobreescribe el resultado en la misma matriz
	virtual double tanIzqMax() const = 0;
	virtual double tanArrMax() const = 0;
};

#define L_CamaraPinhole_calcPixel(cam,rayo,xPix,yPix) ((xPix) = (cam).cenX-(rayo).tanIzq*(cam).dfX, (yPix)=(cam).cenY-(rayo).tanArr*(cam).dfY)
#define L_CamaraPinhole_calcPixeli(cam,rayo,xPix,yPix) ((xPix) = (int)((cam).cenX-(rayo).tanIzq*(cam).dfX+0.5), (yPix)=(int)((cam).cenY-(rayo).tanArr*(cam).dfY+0.5))
#define L_CamaraPinhole_calcRayo(cam,rayo,xPix,yPix)  ((rayo).tanIzq=((cam).cenX-(xPix))/(cam).dfX, (rayo).tanArr=((cam).cenY-(yPix))/(cam).dfY)


class L_CamaraPinhole : public L_Camara
{
public:
	int resX; // Resolucion en el eje X de la imagen
	int resY; // Resolucion en el eje Y de la imagen
	double cenX; // X del Pixel que corresponde al centro de la camara
	double cenY; // Y del Pixel que corresponde al centro de la camara
    double fovX; // Campo visual en radianes para la direccion X de la imagen, no se usa
	double fovY; // Campo visual en radianes para la direccion Y de la imagen, no se usa
	double dfX; // Distancia focal en pixeles para la direccion X de la imagen
	double dfY; // Distancia focal en pixeles para la direccion Y de la imagen

	static inline double calcDistFocX(int resX, double fovX) {return (resX/2)/tan(fovX/2);}
	static inline double calcDistFocY(int resY, double fovY) {return (resY/2)/tan(fovY/2);}

public:
	L_CamaraPinhole() : L_Camara(camPinhole) {}

	// Funciones para fijar parametros
	void fijaParametros(int resX, int resY, double fovX, double fovY, double cenX=-1, double cenY=-1)
	{
		this->resX=resX;
		this->resY=resY;
		this->fovX=fovX;
		this->fovY=fovY;
		if (cenX>=0)
			this->cenX=cenX;
		else
			this->cenX=(resX-1.0)/2;
		if (cenY>=0)
			this->cenY=cenY;
		else
			this->cenY=(resY-1.0)/2;
		dfX=calcDistFocX(this->resX, this->fovX);
		dfY=calcDistFocY(this->resY, this->fovY);
	}

	void fijaParametrosUnitarios() {cenX=0; cenY=0; dfX=-1; dfY=-1;} // Para que tanIzq y xPix sean lo mismo

	void fijaParametrosPhillips(int resX=320, int resY=240) {fijaParametros(resX,resY,45*M_PI/180,37*M_PI/180);} // Resoluciones: 160x120, 320x240, 640x480, camara SPC900NC
	void fijaParametrosVidereNormales(int resX=320, int resY=240) {fijaParametros(resX,resY,30*M_PI/180,28*M_PI/180);} // Al ojo
	void fijaParametrosKinectVisible(int resX=640, int resY=480) {fijaParametros(resX,resY,62.72660441*M_PI/180,49.13434264*M_PI/180);}
	void fijaParametrosKinectIR(int resX=640, int resY=480) {fijaParametros(resX,resY,57.77316353*M_PI/180,44.95886878*M_PI/180);}

	void fijaCentro(double relCenX, double relCenY) // relCenX y relCenY entre 0 y 1
	{
		this->cenX=relCenX*this->resX;
		this->cenY=relCenY*this->resY;
	}

	void fijaMatrizCalibracion(const L_Matrix &K3x3); // (xPix*a yPix*a a) = K*(uv1), esta camara pinhole no tiene shear ni rotacion/traslacion
	void pedirMatrizCalibracion(L_Matrix &K3x3) const; // (xPix*a yPix*a a) = K*(uv1), esta camara pinhole no tiene shear ni rotacion/traslacion

	bool dentroDeLaImagen(double xPix, double yPix) {return xPix>=0 && yPix>=0 && xPix<=resX-1 && yPix<=resY-1;}
	inline double mX(int xPix) {return (cenX-xPix)/dfX;}
	inline double mY(int yPix) {return (cenY-yPix)/dfY;}
	inline double thetaX(int xPix) {if (cenX-xPix==0 && dfX==0) return 0; else return atan2(cenX-xPix, dfX);}
	inline double thetaY(int yPix) {if (cenY-yPix==0 && dfY==0) return 0; else return atan2(cenY-yPix, dfY);}
	inline double mTOT(int xPix, int yPix) {return sqrt(mX(xPix)*mX(xPix)+ mY(yPix)*mY(yPix));}
	inline double thetaTOT(int xPix, int yPix) {return atan(mTOT(xPix,yPix));}
	inline double vignettingFactor(double thetaTotal) {return pow(cos(thetaTotal),4);}
	inline double errorAng_0() const {return L_MAX(1/dfX,1/dfY);} // Aproximacion rapida, revisar si conviene algo mejor

	inline double pideResX() const {return resX;}
	inline double pideResy() const {return resX;}

	inline void calcRayo(L_RayoCamara &rayo, double xPixDer, double yPixAba) const {rayo.tanIzq=(cenX-xPixDer)/dfX; rayo.tanArr=(cenY-yPixAba)/dfY;}
	inline void calcPixel(const L_RayoCamara &rayo, double &xPixDer, double &yPixAba) const {xPixDer = cenX-rayo.tanIzq*dfX; yPixAba=cenY-rayo.tanArr*dfY;}
	inline void calcPixel(const L_RayoCamara &rayo, int &xPixDer, int &yPixAba) const {xPixDer = (int)(cenX-rayo.tanIzq*dfX+0.5), yPixAba=(int)(cenY-rayo.tanArr*dfY+0.5);}
	virtual double radioErrorTan() const {return L_MAX(1/dfX,1/dfY);} // Aproximacion rapida, revisar si conviene algo mejor
	virtual void calcRayos(L_Matrix &uv_nx2) const;
	virtual void calcPixeles(L_Matrix &uv_nx2) const;
	virtual double tanIzqMax() const {return (resX-cenX)/dfX;}
	virtual double tanArrMax() const {return (resY-cenY)/dfY;}
};

class L_CamaraDistorsionRadial : public L_Camara // Modelo de Tsai: directo de pixel a rayo
{
private:
	// transformUsing de pixel a rayo (Tsai)
	// ui = (cenX-xPix)/dfX, vi = (cenY-yPix)/dfY  // coord norm en la imagen = pixel
	// ue = y/x , ve = z/x;                        // coord norm en el espacio = rayo
	//
	// ri2 = ui*ui + vi*vi;
	// ue  = (1 + k1*ri2 + k2*ri2*ri2)*ui;, ve=...
	//
	// Correccion: (xPixImg,yPixImg) -> (ui,vi) -> (ue,ve) -> (xPixIdeal,yPixIdeal)
	//
	// Usualmente se desea copiar un (xPixImg,yPixImg) desconocido en un (xPixIdeal,yPixIdeal) conocido

	int resX; // Resolucion en el eje X de la imagen
	int resY; // Resolucion en el eje Y de la imagen
	double cenX; // X del Pixel que corresponde al centro de la camara
	double cenY; // Y del Pixel que corresponde al centro de la camara
	double k1, k2; // Constantes de distorsion radial directas: (tanIzqD,tanArrD)*( 1 + k1*rD^2 + k2*rD^4 + ...) = (tanIzq,tanArr)
	double k3, k4; // Constantes de distorsion cruzadas: dx = 2*k3*x*y+k4*(r*r+2*x*x) ; dy = k3*(r*r+2*y*y) + 2*k4*x*y
	double dfX, dfY;  // Distancia focal "central" [pix] o [pix/tan]: xPix = -dfX*tanIzqD + cenX , yPix = -dfY*tanArrD + cenY
	bool modif; // Para saber si la LUT es valida o hay que recalcularla
	int tipoModelo;  // 0 = modelo tsai;  1 = modelo Heikkil & Silven (MATLAB)

	L_ImageGrayDouble LUTx, LUTy; // Permiten pasar de (ui,vi) -> (uc,vc), necesarias para deshacer la distorsion

	double xim, yim; // Internos, para minimizar error

public:
	bool activo;  // Variable para uso externo

	L_CamaraDistorsionRadial() : L_Camara(camDistRad) {resX=320; resY=240; cenX=159.5; cenY=119.5; k1=0; k2=0; k3 = 0; k4 = 0; dfX=500; dfY=500; modif = true; tipoModelo = 0; activo = false;} // Valores por defecto

	int resXv() {return resX;}
	int resYv() {return resY;}
	double cenXv() {return cenX;}
	double cenYv() {return cenY;}
	double k1v() {return k1;}
	double k2v() {return k2;}
	double dfXv() {return dfX;}
	double dfYv() {return dfY;}

	void fijaParametros(int resX, int resY, double cenX, double cenY, double k1, double k2, double k3, double k4, double dfX, double dfY, int modelo) {this->resX=resX; this->resY=resY; this->cenX=cenX; this->cenY=cenY; this->k1=k1; this->k2=k2; this->k3 = k3; this->k4 = k4; this->dfX=dfX; this->dfY=dfY; this->tipoModelo = modelo; this->modif = true;}

	// Estas funciones no usan LUT, evitarlas en lo posible
	void distorsionaPixel_k1_Tsai(double xPix, double yPix, double &xPixDist, double &yPixDist) const;
	void distorsionaPixel_k1_Tsai(double xPix, double yPix, int &xPixDist, int &yPixDist) const {double u, v; distorsionaPixel_k1_Tsai(xPix, yPix, u, v); xPixDist=(int)u; yPixDist=(int)v;}
	void distorsionaPixel_k1k2k3k4_Heikkil(double xPix, double yPix, double &xPixDist, double &yPixDist) const;
	void distorsionaPixel_k1k2k3k4_Heikkil(double xPix, double yPix, int &xPixDist, int &yPixDist) const {double u, v; distorsionaPixel_k1k2k3k4_Heikkil(xPix, yPix, u, v); xPixDist=(int)u; yPixDist=(int)v;}
	
	static void error_Heikkil(const void *modelo, const L_Matrix &x, L_Matrix &e);

	void calcLUT_k1_Tsai(int lx, int ly) {if (LUTx.lx != lx || LUTx.ly != ly || LUTy.lx != lx || LUTy.ly != ly || modif) {modif = false; LUTx.reallocate(lx,ly); LUTy.reallocate(lx,ly); int i, j; for (i=0; i<lx; i++) for (j=0; j<ly; j++) distorsionaPixel_k1_Tsai(i,j,LUTx.pix(i,j),LUTy.pix(i,j));}}
	void calcLUT_k1k2k3k4_Heikkil(int lx, int ly) {if (LUTx.lx != lx || LUTx.ly != ly || LUTy.lx != lx || LUTy.ly != ly || modif) {modif = false; LUTx.reallocate(lx,ly); LUTy.reallocate(lx,ly); int i, j; for (i=0; i<lx; i++) for (j=0; j<ly; j++) distorsionaPixel_k1k2k3k4_Heikkil(i,j,LUTx.pix(i,j),LUTy.pix(i,j));}}
	void distorsionaPixel_LUT(int xPix, int yPix, double &xPixDist, double &yPixDist) const {xPixDist = LUTx.pix(xPix,yPix); yPixDist = LUTy.pix(xPix,yPix);}
	void distorsionaPixel_LUT(int xPix, int yPix, int &xPixDist, int &yPixDist) const {xPixDist = L_ROUND10000(LUTx.pix(xPix,yPix)); yPixDist = L_ROUND10000(LUTy.pix(xPix,yPix));}

	void resetearParametros() {modif = true; cenX=(resX-1)*0.5; cenY=(resY-1)*0.5; k1=0; k2=0; dfX=500; dfY=500; tipoModelo = 0;}

	virtual void calcRayo(L_RayoCamara &rayo, double xPix, double yPix) const; // Calcula rayo en el espacio
	virtual void calcPixel(const L_RayoCamara &rayo, double &xPix, double &yPix) const; // Calcula pixel proyectado en la imagen
	virtual void calcRayos(L_Matrix &uv_nx2) const;
	virtual void calcPixeles(L_Matrix &uv_nx2) const;
	virtual double tanIzqMax() const {double tanIzq = cenX/dfX; double r2 = tanIzq*tanIzq; return tanIzq * (1 + k1*r2 + k2*r2*r2);}
	virtual double tanArrMax() const {double tanArr = cenY/dfY; double r2 = tanArr*tanArr; return tanArr * (1 + k1*r2 + k2*r2*r2);}

	void fijaParametrosEstandar(int lx, int ly, double fX=500, double fY=500) {resX=lx; resY=ly; cenX=(resX+1.0)/2; cenY=(resY+1.0)/2; dfX=fX; dfY=fY; /*q1=0; q2=0;*/; k1=0; k2=0;}
	void fijaParametrosVidere_lentesOrig(int resX=320, int resY=240) {this->resX=resX; this->resY=resY; cenX=(resX+1.0)/2; cenY=(resY+1.0)/2; dfX=180.0*(resX/320.0); dfY=180.0*(resY/240.0); /*q1=0.402617; q2=-0.101390;*/ k1=0.50000; k2=0.00000;} // Camara con los lentes que deforman mas la imagen
	void fijaParametrosVidere_lentes2(int resX=320, int resY=240) {this->resX=resX; this->resY=resY; cenX=(resX+1.0)/2; cenY=(resY+1.0)/2; dfX=180.0*(resX/320.0); dfY=180.0*(resY/240.0); /*q1=0.402617; q2=-0.101390;*/ k1=0.50000; k2=0.00000;} // Camara con los lentes que deforman mas la imagen

	// Estas funciones permiten rectificar una imagen y obtener una camara pinhole equivalente para ella
	void corrigeImagenPixelado(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr); // Conserva el FOV central
	void corrigeImagenBilineal(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr); // Conserva el FOV central
	void corrigeImagenBilineal_k1(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr); // Conserva el FOV central
	void corrigeImagenPixelado(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr); // Conserva el FOV central
	void corrigeImagenBilineal(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr); // Conserva el FOV central
	void corrigeImagenBilineal_k1(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr); // Conserva el FOV central
	void muestraImagenError(const L_ImageRGBUchar &orig, L_ImageRGBUchar &final, int niter);
	void muestraImagenError(const L_ImageGrayDouble &orig, L_ImageGrayDouble &final, int niter);
	void distorsionaPixel(double xPix, double yPix, double &xPixDist, double &yPixDist); // Conserva el FOV central

	void generaCamaraPinhole(L_CamaraPinhole &cam) {cam.fijaParametros(resX, resY, 2*atan2(resX, 2*dfX), 2*atan2(resY, 2*dfY), cenX, cenY);}
	double fovXCorregida() {return 2.0*atan2(resX, 2*dfX);}
	double fovYCorregida() {return 2.0*atan2(resY, 2*dfY);}
	virtual double radioErrorTan() const {return L_MAX(1/dfX,1/dfY);} // Aproximacion rapida, revisar si conviene algo mejor

	bool fijaDistorsionK1_Tsai(double k1) {modif = true; this->k1=k1; tipoModelo = 0; return true;}
	bool fijaDistorsionK1K2K3K4_Heikkil(double k1, double k2, double k3, double k4) {modif = true; this->k1=k1; this->k2=k2; this->k3=k3; this->k4=k4; tipoModelo = 1; return true;}

	void guardarParametros(FILE *fp);
	bool leerParametros(FILE *fp);

//	static double fnErrorEntreQyK_k(void *object, double *vectorK);
//	static double fnErrorEntreQyK_q(void *object, double *vectorQ);
};

class L_CamaraDistorsionRadialInv_Ajustador // No implementado ahora
{
public:
	L_CamaraDistorsionRadial cam;
	L_Matrix m;  // xyzuv_nx5 : (x,y,z)comparten un mismo sistema global desconocido, (u,v) pixeles asociados
	double tx, ty, tz; // Traslaciones
	double theta, psi, gamma; //Rotaciones
	// Dimensionalidad del problema: 10 (3 traslaciones, 3 rotaciones, 4 parametros camara)

	L_CamaraDistorsionRadialInv_Ajustador() {}

	void ajuste_inicial(); // Ajuste inicial suponiendo que se esta mirando un plano de frente y camara algo ajustada
	static double fnObjetivo_L2(void *object, double *vector); // Error cuadratico
	static double fnObjetivo_L1(void *object, double *vector); // Error no cuadratico (inestable pero robusto a falsos positivos)
};

// Para  calibrar camaras
class L_PuntoCalibradorCamara
{
public:
	double derMM, abaMM, adeMM, xPix, yPix;
};

class L_CuadroCalibradorCamara
{
public:
	L_Array<L_PuntoCalibradorCamara> v;

	L_CuadroCalibradorCamara() {v.reserve(200);}
	bool guardar(FILE *fp);
	bool read(FILE *fp);
	bool guardar_TsaiCM(const char *nomarch);
	void agregar(double derMM, double abaMM, double adeMM, double xPix, double yPix) {L_PuntoCalibradorCamara e; e.derMM = derMM; e.abaMM=abaMM; e.adeMM=adeMM; e.xPix=xPix; e.yPix=yPix; v.push_back(e);}
};

class L_CuadroCalibradorCamaraArr
{
public:
	std::vector<L_CuadroCalibradorCamara> v;

	L_Statistics1D nPuntosEstad;

	size_t size() {return v.size();}
	L_CuadroCalibradorCamara &operator[](int i) {return v[i];}
	const L_CuadroCalibradorCamara &operator[](int i) const {return v[i];}
	bool guardar(FILE *fp);
	bool read(FILE *fp);
	void agregarCuadro() {v.resize(v.size()+1); nPuntosEstad.push(v.size());}
	void agregarPunto(double derMM, double abaMM, double adeMM, double xPix, double yPix) {if (v.size()>0) v[v.size()-1].agregar(derMM, abaMM, adeMM, xPix, yPix);}
	void resetear() {v.clear(); nPuntosEstad.clear();}
};



/////////
// Clases para define transformaciones geometricas

class L_CamaraPinhole;

class L_RelGeom2D // Relacion geometrica:  f(datos,parametros) = 0
{
public:
	virtual double &el(int i) = 0;
	virtual const double &el(int i) const = 0;
	virtual int nEl() = 0;
	virtual int gradosLibertad() = 0; // Num calces para num finito de soluciones
	virtual int gradosLibertadU() = 0; // Num calces para solucion unica
	virtual int numSolsGradosLibertadU() = 0;
	virtual void precalcula() {} // Se llama antes de errorRelGeom() despues de modificar los parametros
	virtual double errorRelGeom(double x1, double y1, double x2, double y2) = 0; // Aca da lo mismo la convencion de ejes
	virtual void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1) {err.resize(x1y1x2y2.li); int nMax = x1y1x2y2.li; if (nRev >= 0 && iIni + nRev < nMax) nMax = iIni + nRev; for (int i=iIni; i<nMax; i++) err[i] = errorRelGeom(x1y1x2y2(i,0), x1y1x2y2(i,1), x1y1x2y2(i,2), x1y1x2y2(i,3));}
	virtual int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt); // Funcion por defecto, lenta
	virtual int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes) = 0;
	virtual int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2) = 0;
	virtual bool calcExacto(const L_Matrix &x1y1x2y2) = 0;
	virtual bool calcMinCuad(const L_Matrix &x1y1x2y2) = 0;
	virtual bool invertMe() = 0;
	virtual bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) = 0;
	virtual bool normaliza() {return true;}
	virtual void OP_assign(L_RelGeom2D *other) = 0;
	virtual ~L_RelGeom2D() {}

	// calcRansac() devuelve la mejor solucion que encontro hasta el momento
	int calcRansac(const L_Matrix &x1y1x2y2, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c=-1, int d=-1);
	int calcProsacSimple(const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c=-1, int d=-1);
	int calcRansacPreemptivo(const L_Matrix &x1y1x2y2, L_Array<bool> &seleccionados, double umbralError, int nHipotesis, int nScorings, int c, int d); // La idea es que c y d sean pequenos y c < d ... nScorings debe ser mucho mayor al numero de datos... intenta ser "Preemptive RANSAC for Live Structure and Motion Estimation" de Nister 
	int calcProsacPreemptivo(const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nHipotesis, int nScorings, int c, int d); // La idea es que c y d sean pequenos y c < d ... nScorings debe ser mucho mayor al numero de datos... intenta ser "Preemptive RANSAC for Live Structure and Motion Estimation" de Nister 
	int numIntentosRansac(double probInlierIndividual, double probEncontrarSolucion) {return (int)(1+log(1-probEncontrarSolucion) / log(1-pow(probInlierIndividual,gradosLibertadU())));}
};

class L_RelEpipolarFundamental : public L_RelGeom2D
{
public:
	L_MatrizFundamental F;
	// El error es con signo, ojo
	double &el(int i) {return F(i/3,i%3);}
	const double &el(int i) const {return F(i/3,i%3);}
	int nEl() {return 9;}
	int gradosLibertad() {return 7;} // Se set con el algoritmo de los 7 puntos
	int gradosLibertadU() {return 8;} // Algoritmo de los 8 puntos
	int numSolsGradosLibertadU() {return 3;} // Algoritmo de los 7 puntos, 1 o 3 soluciones
	void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1);
	int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2);
	bool calcExacto(const L_Matrix &x1y1x2y2) {return L_RelEpipolarFundamental::calcHartley(F, x1y1x2y2);}
	bool calcMinCuad(const L_Matrix &x1y1x2y2) {return L_RelEpipolarFundamental::calcHartley(F, x1y1x2y2);}
	bool invertMe() {return F.invertMe();}
	bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) {return false;} // Es innecesario aca
	bool normaliza() {return F.reparaFundamental();}
	virtual void OP_assign(L_RelGeom2D *other) {*this = *dynamic_cast<L_RelEpipolarFundamental*>(other);}

	static bool calcHartley(L_Matrix &F, const L_Matrix &x1y1x2y2, bool escalamiento = true);
	static bool escalamientoIsotherpico(L_Matrix &uv, L_Matrix &T1, L_Matrix &T2);
	static void escalamientoIsotherpicoInv(L_Matrix &F, L_Matrix &T1, L_Matrix &T2);

	static bool algoritmo7puntos(L_Matrix F[3], int &nSols, const L_Matrix &x1y1x2y2);  // esta MALA, de ahi la arreglo
	static bool algoritmo7_mas_1_puntos(L_Matrix &F, const L_Matrix &x1y1x2y2);  // esta MALA, de ahi la arreglo
	static bool despejar_XY_7p(L_Matrix &X, L_Matrix &Y, const L_Matrix &QagrT);

	inline double errorRelGeom(double x1, double y1, double x2, double y2)
	{
		double ux, uy, uz, d, e1, e2;
		ux = F(0,0)*x1+F(0,1)*y1+F(0,2);
		uy = F(1,0)*x1+F(1,1)*y1+F(1,2);
		uz = F(2,0)*x1+F(2,1)*y1+F(2,2);
		d = sqrt(ux*ux+uy*uy);
		ux/=d;
		uy/=d;
		uz/=d;
		e1 = x2*ux+y2*uy+uz;
		ux = F(0,0)*x2+F(1,0)*y2+F(2,0);
		uy = F(0,1)*x2+F(1,1)*y2+F(2,1);
		uz = F(0,2)*x2+F(1,2)*y2+F(2,2);
		d = sqrt(ux*ux+uy*uy);
		ux/=d;
		uy/=d;
		uz/=d;
		e2 = x1*ux+y1*uy+uz;
		return sqrt(e1*e1+e2*e2);
	}
	static void probando_probando();
};

class L_RelEpipolarEsencial : public L_RelGeom2D
{
public:
	L_EssentialMatrix E;
	// El error es con signo, ojo
	double &el(int i) {return E(i/3,i%3);}
	const double &el(int i) const {return E(i/3,i%3);}
	int nEl() {return 9;}
	int gradosLibertad() {return 5;} // Se set con el algoritmo de los 5 puntos
	int gradosLibertadU() {return 6;} // Algoritmo de los 5+1 puntos
	int numSolsGradosLibertadU() {return 10;} // Algoritmo de los 5

	void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1);
	int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2);
	bool calcExacto(const L_Matrix &x1y1x2y2) {return L_Algoritmo5Puntos::algoritmo5_mas_1_puntos(E, x1y1x2y2);}
	bool calcMinCuad(const L_Matrix &x1y1x2y2) {return calcHartleyEsencial(E, x1y1x2y2);}
	bool invertMe() {return E.invertMe();}
	bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) {return false;} // Es innecesario aca
	bool normaliza() {return E.reparaEsencial();}
	virtual void OP_assign(L_RelGeom2D *other) {*this = *dynamic_cast<L_RelEpipolarEsencial*>(other);}

	static bool calcHartleyEsencial(L_Matrix &E, const L_Matrix &x1y1x2y2) {return L_RelEpipolarFundamental::calcHartley(E, x1y1x2y2) && E.svd_reparaEsencial();}

	inline double errorRelGeom(double x1, double y1, double x2, double y2)
	{
		double ux, uy, uz, fac, e1, e2;
		ux = E(0,0)*x1+E(0,1)*y1+E(0,2);
		uy = E(1,0)*x1+E(1,1)*y1+E(1,2);
		uz = E(2,0)*x1+E(2,1)*y1+E(2,2);
		fac = sqrt(ux*ux+uy*uy);
		ux/=fac;
		uy/=fac;
		uz/=fac;
		e1 = x2*ux+y2*uy+uz;
		ux = E(0,0)*x2+E(1,0)*y2+E(2,0);
		uy = E(0,1)*x2+E(1,1)*y2+E(2,1);
		uz = E(0,2)*x2+E(1,2)*y2+E(2,2);
		fac = sqrt(ux*ux+uy*uy);
		ux/=fac;
		uy/=fac;
		uz/=fac;
		e2 = x1*ux+y1*uy+uz;
		return sqrt(e1*e1+e2*e2);
	}

	static void probando_probando();
};




////////
//


class L_Transf2D : public L_RelGeom2D //  Transformacion geometrica: datosSalida=f(datosEntrada,parametros)
{
public:
	virtual void proyeccionDe(double &resx, double &resy, double x, double y) const = 0; // Aca da lo mismo la convencion de ejes
	double errorRelGeom(double x1, double y1, double x2, double y2) {double xP, yP, ex, ey; proyeccionDe(xP,yP,x1,y1); ex = x2-xP; ey = y2-yP; return sqrt(ex*ex+ey*ey);}
};

class L_TransfAfin2D : public L_Transf2D
{
public:
	double m11, m12, m21, m22, tx, ty;
	double &el(int i) {switch(i) {case 0: return m11; case 1: return m12; case 2: return m21; case 3: return m22; case 4: return tx; case 5: return ty; default: throw_L_ArgException_if(true, "L_TransfAfin2D::el() : fuera de rango"); return m11;}}
	const double &el(int i) const {switch(i) {case 0: return m11; case 1: return m12; case 2: return m21; case 3: return m22; case 4: return tx; case 5: return ty; default: throw_L_ArgException_if(true, "L_TransfAfin2D::el() : fuera de rango"); return m11;}}
	int nEl() {return 6;}
	int gradosLibertad(void) {return 3;}
	int gradosLibertadU(void) {return 3;}
	int numSolsGradosLibertadU(void) {return 3;}
	void proyeccionDe(double &resx, double &resy, double x, double y) const {resx=m11*x+m12*y+tx; resy=m21*x+m22*y+ty;}
	void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1) {int nMax = x1y1x2y2.li; if (nRev >= 0 && iIni + nRev < nMax) nMax = iIni + nRev; double resx, resy, dx, dy; err.resize(x1y1x2y2.li); for (int i=iIni; i<nMax; i++) {resx=m11*x1y1x2y2(i,0)+m12*m11*x1y1x2y2(i,1)+tx; resy=m21*m11*x1y1x2y2(i,0)+m22*m11*x1y1x2y2(i,1)+ty; dx = resx-x1y1x2y2(i,2); dy = resy-x1y1x2y2(i,3); err[i] = sqrt(dx*dx + dy*dy);}}
	int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2);
	bool calcExacto(const L_Matrix &x1y1x2y2);
	bool calcMinCuad(const L_Matrix &x1y1x2y2); 
	bool invertMe();
	bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) {return false;} // Es innecesario aca
	virtual void OP_assign(L_RelGeom2D *other) {*this = *dynamic_cast<L_TransfAfin2D*>(other);}

	void calcMatricesRayos(const L_CamaraPinhole &cR, const L_CamaraPinhole &cP, L_Matrix &M, L_Matrix &t);
	void descomponeSVD(double &escMax, double &escMin, double &rot, const L_CamaraPinhole *cR=NULL, const L_CamaraPinhole *cP=NULL);
};

class L_TransfSimil2D : public L_Transf2D
{
public:
	double esc, ang, tx, ty;
	double m11, m12, m21, m22; // Para calcular rapido
	double &el(int i) {switch(i) {case 0: return esc; case 1: return ang; case 2: return tx; case 3: return ty; default: throw_L_ArgException_if(true, "L_TransfProyectiva2D::el() : fuera de rango"); return esc;}}
	const double &el(int i) const {switch(i) {case 0: return esc; case 1: return ang; case 2: return tx; case 3: return ty; default: throw_L_ArgException_if(true, "L_TransfProyectiva2D::el() : fuera de rango"); return esc;}}
	int nEl() {return 4;}
	int gradosLibertad(void) {return 2;}
	int gradosLibertadU(void) {return 2;}
	int numSolsGradosLibertadU(void) {return 1;}
	void precalcula() {double c=cos(-ang), s=sin(-ang); m11=esc*c; m12=-esc*s; m21=esc*s; m22=esc*c;}
	void proyeccionDe(double &resx, double &resy, double x, double y) const {resx=m11*x+m12*y+tx; resy=m21*x+m22*y+ty;}
	void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1);
	int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2);
	bool calcExacto(const L_Matrix &x1y1x2y2);
	bool calcMinCuad(const L_Matrix &x1y1x2y2); 
	bool invertMe();
	bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) {return false;} // Es innecesario aca
	virtual void OP_assign(L_RelGeom2D *other) {*this = *dynamic_cast<L_TransfSimil2D*>(other);}
};

class L_TransfProyectiva2D : public L_Transf2D
{
public:
	double esc, ang;
	double m11, m12, m21, m22, m31, m32, tx, ty;
	double &el(int i) {switch(i) {case 0: return m11; case 1: return m12; case 2: return m21; case 3: return m22; case 4: return m31; case 5: return m32; case 6: return tx; case 7: return ty; default: throw_L_ArgException_if(true, "L_TransfProyectiva2D::el() : fuera de rango"); return m11;}}
	const double &el(int i) const {switch(i) {case 0: return m11; case 1: return m12; case 2: return m21; case 3: return m22; case 4: return m31; case 5: return m32; case 6: return tx; case 7: return ty; default: throw_L_ArgException_if(true, "L_TransfProyectiva2D::el() : fuera de rango"); return m11;}}
	int nEl() {return 8;}
	int gradosLibertad(void) {return 4;}
	int gradosLibertadU(void) {return 4;}
	int numSolsGradosLibertadU(void) {return 1;}
	void proyeccionDe(double &resx, double &resy, double x, double y) const {resx=m11*x+m12*y+tx; resy=m21*x+m22*y+ty; double resz=m31*x+m32*y+1; resx/=resz; resy/=resz;}
	void calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni=0, int nRev = -1) {int nMax = x1y1x2y2.li; if (nRev >= 0 && iIni + nRev < nMax) nMax = iIni + nRev; double resx, resy, resz, dx, dy; err.resize(x1y1x2y2.li); for (int i=iIni; i<nMax; i++) {resx=m11*x1y1x2y2(i,0)+m12*m11*x1y1x2y2(i,1)+tx; resy=m21*m11*x1y1x2y2(i,0)+m22*m11*x1y1x2y2(i,1)+ty; resz = m31*m11*x1y1x2y2(i,0)+m32*m11*x1y1x2y2(i,1) + 1; dx = (resx-x1y1x2y2(i,2))/resz; dy = (resy-x1y1x2y2(i,3))/resz; err[i] = sqrt(dx*dx + dy*dy);}}
	int test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeom2D *> &hip, const L_Matrix &x1y1x2y2);
	bool calcExacto(const L_Matrix &x1y1x2y2);
	bool calcMinCuad(const L_Matrix &x1y1x2y2);
	bool invertMe();
	bool componer(const L_RelGeom2D &prim, const L_RelGeom2D &seg) {return false;} // Es innecesario aca
	virtual void OP_assign(L_RelGeom2D *other) {*this = *dynamic_cast<L_TransfProyectiva2D*>(other);}
};


class L_RelGeomMuestra2D // Relacion geometrica para una muestra:  f(datos,parametros) = 0 solo dentro de un set de datos dado
{
public:
	// La class_name derivada ya debe contener el set de datos
	virtual double &el(int i) = 0;
	virtual const double &el(int i) const = 0;
	virtual int nEl() = 0;
	virtual int tamDatos() = 0;
	virtual int gradosLibertad() = 0; // Num calces para num finito de soluciones
	virtual int gradosLibertadU() = 0; // Num calces para solucion unica
	virtual int numSolsGradosLibertadU() = 0; // Num calces para solucion unica
	virtual void precalcula() {} // Se llama antes de errorRelGeom() despues de modificar los parametros
	virtual double errorRelGeom(int i) = 0; // Aca da lo mismo la convencion de ejes
	virtual void calcVectorError(std::vector<double> &err, int iIni=0, int nRev = -1) = 0;
	virtual int test_c_de_d(double errMax, int c, int d, int iIni, int &iUlt) = 0;
	virtual int generarHipotesis(L_Array<L_RelGeomMuestra2D *> &hip, const std::vector<double> &puntajes) = 0;
	virtual int generarHipotesis(L_Array<L_RelGeomMuestra2D *> &hip) = 0;
	virtual bool calcExacto(const std::vector<int> &selecc) = 0;
	virtual bool calcMinCuad() = 0;
	virtual bool normaliza() {return true;}
	virtual void OP_assign(L_RelGeomMuestra2D *other) = 0;
	virtual ~L_RelGeomMuestra2D() {}

	// calcRansac() devuelve la mejor solucion que encontro hasta el momento
	int calcRansac(L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c=-1, int d=-1);
	int calcProsacSimple(const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c=-1, int d=-1);
	int numIntentosRansac(double probInlierIndividual, double probEncontrarSolucion) {return (int)(1+log(1-probEncontrarSolucion) / log(1-pow(probInlierIndividual,gradosLibertadU())));}
};

class L_TransfPoseProy2D : public L_RelGeomMuestra2D
{
public:
	L_Pose3D_cuat pose;
	std::vector<L_CoordsCart3D> pun; // Los puntos pun se mueven segun pose y luego se proyectan
	std::vector<L_RayoCamara> ray;
	std::vector<L_CoordsCart3D> punSel;
	std::vector<L_RayoCamara> raySel;

	L_TransfPoseProy2D() : punSel(4), raySel(4) { }

	double &el(int i) {return pose.el(i);}
	const double &el(int i) const {return pose.el(i);}
	int nEl() {return 7;}
	int tamDatos() {return (int)pun.size();}
	int gradosLibertad(void) {return 3;} // Se set con el algoritmo de los 3 puntos
	int gradosLibertadU(void) {return 4;}  // 4+1 puntos
	int numSolsGradosLibertadU(void) {return 4;}  // algoritmo 3 puntos, 4 soluciones
	int test_c_de_d(double errMax, int c, int d, int iIni, int &iUlt);
	int generarHipotesis(L_Array<L_RelGeomMuestra2D *> &hip, const std::vector<double> &puntajes);
	int generarHipotesis(L_Array<L_RelGeomMuestra2D *> &hip);
	double errorRelGeom(int i) {L_RayoCamara r; double errx, erry, err2; double xs; pose.prmovi_v(r, pun[i], xs); errx = r.tanIzq - r.tanIzq; erry = r.tanArr - r.tanArr; err2 = errx*errx + erry*erry; return (xs>0) ? sqrt(err2) : err2*1e10*xs;}
	void calcVectorError(std::vector<double> &err, int iIni=0, int nRev = -1);
	bool calcExacto(const std::vector<int> &selecc);
	bool calcMinCuad(); 
	bool normaliza() {return pose.ori.normaliza_ret();} // Es para que lo llame Ransac
	virtual void OP_assign(L_RelGeomMuestra2D *other) {*this = *dynamic_cast<L_TransfPoseProy2D*>(other);}
};

class L_TransfCov_tresMatrices  // class_name para dividir covarianzas, se usa para el slam...
{
public:
	L_Matrix Ppp;
	L_Matrix Ptridiag;
	L_Matrix Pdiag;
	void init(const L_Matrix &Pppe);
	double error(double lambda1, double lambda2) const;
	static double errorGen(const void *object, double *vect);
	static void errorGen_vect(const void *object, const L_Matrix &vect, L_Matrix &err);
	void encontrarCovarianzasIndividuales(const L_Matrix &Ppp, L_Matrix &P_arr_3x3);
	void encontrarCovarianzasIndividuales_vect(const L_Matrix &Ppp, L_Matrix &P_arr_3x3);
	static void pruebame();
};


///////
// Capturadores de camaras
///////

class L_CapturadorImagenAbstracto
{
public:
	L_ImageRGBUchar im;
	// For stereo
	L_ImageRGBUchar imDer; // Por si es estereo...
	// For Kinect
	L_ImageGrayUint16 imKin;
	L_ImageGrayDouble imDepth;
	L_ImageGrayDouble imDepthRGB;
	bool activo;
	struct L_soporteEstereo {bool esEstereo; int usarIzq0Der1Ambas2;}
		estereo;
	struct L_soporteProf {bool tieneProf; bool capturarProf; bool calcularProf; bool calcularProfRGB;}
		prof;
	struct L_forzarResolucion {bool im_activo; bool im_resample; int im_lx; int im_ly; bool imDer_activo; bool imDer_resample; int imDer_lx; int imDer_ly; bool imProf_activo; bool imProf_resample; int imProf_lx; int imProf_ly;}
		redimensionar;

	L_CapturadorImagenAbstracto() : activo(false) {estereo.esEstereo = false; estereo.usarIzq0Der1Ambas2 = 0; prof.tieneProf = false; prof.capturarProf = false; prof.calcularProf = false; redimensionar.im_activo = false; redimensionar.imDer_activo = false; redimensionar.imProf_activo = false;}
	virtual bool crear() = 0;
	virtual bool fijarResolucion(int width, int height) {return false;}
	bool forzarResolucion(bool im_activo, int im_lx, int im_ly, bool imDer_activo, int imDer_lx, int imDer_ly, bool imProf_activo, int imProf_lx, int imProf_ly) {redimensionar.im_activo = im_activo; redimensionar.im_lx = im_lx; redimensionar.im_ly = im_ly; redimensionar.imDer_activo = imDer_activo; redimensionar.imDer_lx = imDer_lx; redimensionar.imDer_ly = imDer_ly; redimensionar.imProf_activo = imProf_activo; redimensionar.imProf_lx = imProf_lx; redimensionar.imProf_ly = imProf_ly; return true;}
	bool forzarResolucion_resample(bool im_resample, bool imDer_resample, bool imProf_resample) {redimensionar.im_resample = im_resample; redimensionar.imDer_resample = imDer_resample; redimensionar.imProf_resample = imProf_resample; return true;}

	virtual bool capturarImagen() = 0; // Debe abrir el dispositivo al llamarla si estaba cerrado	virtual bool cerrar()=0; // Se debe llamar automaticamente antes de destruir el object
	virtual bool cerrar() = 0;
	virtual double tiempoSeg() = 0; // Alguna medida de tiempo en segundos
	virtual bool esArchivo() = 0;
	virtual ~L_CapturadorImagenAbstracto() {} // Destructor virtual, se llama siempre al children y luego a este
};

#ifdef __COMPAT_IPLIMAGE__
class L_CapturadorImagenCV:public L_CapturadorImagenAbstracto
{
public:
	CvCapture *capture;
	IplImage *imCam; // Memoria de esto manejada por OpenCV
	L_String name; // name de archivo
	double dt;
	bool suponerAnchaEstereo;
	long nFrame;
	//
	L_CapturadorImagenCV() {capture=NULL; imCam=NULL;nFrame=0; suponerAnchaEstereo = true;}
	L_CapturadorImagenCV(const char *nomarch) {capture=NULL; imCam=NULL; name = nomarch; nFrame = 0;}
	virtual bool crear();
	virtual bool capturarImagen();
	virtual bool cerrar() {if (capture != NULL) cvReleaseCapture(&capture); capture = NULL; name.clear(); activo=false; return true;}
	virtual double tiempoSeg() {if (name.begin() == NULL) return L_TIME(); return nFrame*dt;} // Alguna medida de tiempo en segundos
	virtual bool esArchivo() {return name.begin() != NULL;}
	virtual ~L_CapturadorImagenCV() {if (activo) cerrar();}
};
#endif

#ifdef __COMPAT_SVS__
class L_CapturadorImagenSVS:public L_CapturadorImagenAbstracto
{
private:
	bool crearRec(int n);
public:
	svsStereoImage *imageObject; // Memoria de esto manejada por SVS
	svsVideoImages *videoObject; // Memoria de esto manejada por SVS
	//svsStereoProcess *processObject; // Se puede incluir para que procese las imagenes
	int lx, ly;
	L_CapturadorImagenSVS(int w=320, int h=240) {imageObject=NULL; videoObject=NULL; lx=w; ly=h; estereo.esEstereo=true; estereo.usarIzq0Der1Ambas2 = 2;}
	virtual bool crear() {return crearRec(2);} // La camara se niega a abrir como un 30% de las veces, con esto eso disminuye...
	virtual bool capturarImagen();
	virtual bool cerrar() {if (videoObject != NULL) videoObject->Close(); videoObject = NULL; activo=false; return true;}
	virtual double tiempoSeg() {return L_TIME();}
	virtual bool esArchivo() {return false;}
	virtual ~L_CapturadorImagenSVS() {cerrar();}
};
#endif

class L_CapturadorImagenBMZ:public L_CapturadorImagenAbstracto
{
public:
	L_String name;
	FILE *fp;
	double t;
	int lx, ly;
	L_CapturadorImagenBMZ(const char *name) {this->name = name; activo = true;fp = fopen(name, "rb"); if (fp==NULL) activo = false;}
	virtual bool fijarResolucion(int width, int height) { return false; }
	virtual bool crear() {return fp!=NULL;}
	virtual bool capturarImagen();
	virtual bool cerrar() {if (fp != NULL) {fclose(fp); fp = NULL;} activo=false; name.clear(); return true;}
	virtual double tiempoSeg() {return t;}
	virtual bool esArchivo() {return true;}
	virtual ~L_CapturadorImagenBMZ() {if (fp != NULL) fclose(fp);}
};

class L_CapturadorImagenesCarpeta : public L_CapturadorImagenAbstracto
{
public:
	L_StringWithNumber name;
	L_StringWithNumber nombreDer;
	L_StringWithNumber nombreKin;
	L_StringWithNumber nombreProf;

	int i;
	L_CapturadorImagenesCarpeta(L_StringWithNumber name, int iIni) {this->name = name; i = iIni;}
	L_CapturadorImagenesCarpeta(L_StringWithNumber name, L_StringWithNumber nombreKin, int iIni) {this->name = name; this->nombreKin = nombreKin; i = iIni;}
	virtual bool crear() {return name.prefix.begin() != NULL;}
	virtual bool capturarImagen();
	virtual bool cerrar() {return true;}
	virtual double tiempoSeg() {return L_TIME();} // Aca si que no se puede estimar
	virtual bool esArchivo() {return false;}
	virtual bool grabarInfoExtra(const char *nomarch_base) {return true;}
};

#ifdef __COMPAT_OPENNI__
class L_CapturadorImagenKinectNI : public L_CapturadorImagenAbstracto
{
public:
	xn::Context context;
	xn::DepthGenerator depthGen;
	xn::ImageGenerator imageGen;

	L_CapturadorImagenKinectNI() {prof.tieneProf = true; prof.capturarProf = true; prof.calcularProf = true; prof.calcularProfRGB = false;}
	virtual bool crear();
	virtual bool capturarImagen();
	virtual bool cerrar();
	virtual double tiempoSeg() {return L_TIME();}
	virtual bool esArchivo() {return false;}

	double fn_prof1(int depth) {return tan(depth/1024.0f+0.5)*33.825+5.7;}
	double fn_prof2(int depth) {return 1.0 / ((double)(depth) * -0.0030711016 + 3.3309495161);}
	double IR_FOV() {return 580.0;}  // @ 640x480
	double RGB_FOV() {return 525.0;} // @ 640x480

	static void test_me();

	virtual ~L_CapturadorImagenKinectNI() {if (activo) cerrar();}
};
#endif // __COMPAT_OPENNI__

#ifdef __COMPAT_KINECTSDK__
class L_CapturadorImagenKinectSDK : public L_CapturadorImagenAbstracto
{
public:
	HANDLE imageEvent;
	HANDLE imageHandle;
	HANDLE depthEvent;
	HANDLE depthHandle;

	double depthLUT[2048];
	int maxDepthLUT;

	L_CapturadorImagenKinectSDK() : imageEvent(NULL), imageHandle(NULL), depthEvent(NULL), depthHandle(NULL) {prof.tieneProf = true; prof.capturarProf = true; prof.calcularProf = true;}
	virtual bool crear();
	virtual bool capturarImagen();
	virtual bool capturarImagenOld();
	virtual bool cerrar();
	virtual double tiempoSeg() {return L_TIME();}
	virtual bool esArchivo() {return false;}

	double fn_prof1(int depth) {return tan(depth/1024.0f+0.5)*33.825+5.7;}
	double fn_prof2(int depth) {return 1.0 / ((double)(depth) * -0.0030711016 + 3.3309495161);}
	double IR_FOV() {return 580.0;}  // @ 640x480
	double RGB_FOV() {return 525.0;} // @ 640x480

	double fn_prof(int depth) {return fn_prof1(depth);}

	virtual ~L_CapturadorImagenKinectSDK() {if (activo) cerrar();}
};
#endif // __COMPAT_KINECTSDK__



// Capturador de imagen que busca lo que haya y lo usa
class L_CapturadorImagen:public L_CapturadorImagenAbstracto // El mas comodo, usa lo que hay...
{
public:
	L_Array <L_CapturadorImagenAbstracto *> capt;
	int numActivo;

	L_CapturadorImagen();
	L_CapturadorImagen(const char *name); // Para read video .bmz 1 o 2 cam
	L_CapturadorImagen(L_StringWithNumber name, int iIni); // Para read video .bmz 1 o 2 cam
	L_CapturadorImagen(L_StringWithNumber name, L_StringWithNumber nombreKin, int iIni); // Para read video .bmz 1 o 2 cam
	virtual bool fijarResolucion(int width, int height);
	virtual bool crear(); // La primera vez que se llama, elige cual camara usar
	virtual bool capturarImagen(); // La primera vez que se llama, elige cual camara usar
	virtual bool redimensionarImagenes();
	virtual bool cerrar();
	virtual double tiempoSeg();
	virtual bool esArchivo() {if (activo) return capt[numActivo]->esArchivo(); return false;}
	virtual ~L_CapturadorImagen();

	const char *nombreVideo();
};

#ifdef __COMPAT_FLTK__
double L_PideNumFLTK(const char *mensaje, double porDefecto);
#endif

#ifdef __COMPAT_IPLIMAGE__
void L_cvCallback(int event, int x, int y, int flags, void* param);
#endif

enum L_Mouse {L_MouseLeido=0, L_MousePresiona, L_MouseSuelta, L_MouseArrastra};

class L_VentanaImagenAbstracta
{
public:
	virtual int mouseX() = 0;
	virtual int mouseY() = 0;
	virtual L_Mouse &mouseE() = 0;
	virtual bool dibuja(const L_ImageRGBUchar &im) = 0;
	virtual bool dibujaRedimensiona(const L_ImageRGBUchar &im) = 0;
	virtual int w() = 0; // width
	virtual int h() = 0; // height
	virtual ~L_VentanaImagenAbstracta() {}
};

#ifdef __COMPAT_IPLIMAGE__
class L_VentanaImagenCV : public L_VentanaImagenAbstracta
{
private:
	L_VentanaImagenCV(void); // No valido
	L_VentanaImagenCV(const L_VentanaImagenCV &other); // No valido
	int _mouseX;
	int _mouseY;
	L_Mouse _mouseE;
public:
	char name[128];
	IplImage *imCV;
	int _w;
	int _h;

	friend void L_cvCallback(int event, int x, int y, int flags, void* param);

	L_VentanaImagenCV(int lx, int ly, const char *name) : _w(lx), _h(ly) {imCV = NULL; _mouseE = L_MouseLeido; L_snprintf(this->name, 127, "%s", name); cvNamedWindow(name); cvResizeWindow(name, lx, ly); cvSetMouseCallback(name, &L_cvCallback, this);}
	L_VentanaImagenCV(int x, int y, int lx, int ly, const char *name) :  _w(lx) , _h(ly) {imCV = NULL; _mouseE = L_MouseLeido; L_snprintf(this->name, 127, "%s", name); cvNamedWindow(name); cvMoveWindow(name, x, y); cvResizeWindow(name, lx, ly); cvSetMouseCallback(name, &L_cvCallback, this);}
	bool dibuja(const L_ImageRGBUchar &im) {im.copyTo(&imCV); cvShowImage(name, imCV); return true;}
	bool dibujaRedimensiona(const L_ImageRGBUchar &im) {im.copyTo(&imCV); cvResizeWindow(name, im.lx, im.ly); _w = im.lx; _h = im.ly; cvShowImage(name, imCV); return true;}
	bool resize(int w, int h) {cvResizeWindow(name, w, h); _w = w; _h = h; return true;}
	bool resize(int x, int y, int w, int h) {cvMoveWindow(name, x, y); cvResizeWindow(name, w, h); _w = w; _h = h; return true;}
	int w() {return _w;}
	int h() {return _h;}
	int mouseX() {return _mouseX;}
	int mouseY() {return _mouseY;}
	L_Mouse &mouseE() {return _mouseE;}
	static void check() {cvWaitKey(1);}
	virtual ~L_VentanaImagenCV() {cvDestroyWindow(name);}
};
#endif // __COMPAT_IPLIMAGE__

#ifdef __COMPAT_FLTK__
class L_Fl_Double_Window : public Fl_Double_Window
{
private:
	L_Fl_Double_Window(void); // No valido
	L_Fl_Double_Window(const L_Fl_Double_Window &other); // No valido
	int _mouseX;
	int _mouseY;
	L_Mouse _mouseE;
public:
	L_Fl_Double_Window(int lx, int ly, const char *name): Fl_Double_Window(lx, ly, name) , _mouseE(L_MouseLeido) {show(); im.reallocate(lx, ly);}
	L_Fl_Double_Window(int x, int y, int lx, int ly, const char *name): Fl_Double_Window(x, y, lx, ly, name) , _mouseE(L_MouseLeido) {show(); im.reallocate(lx, ly);}
public:
	L_ImageRGBAuchar_row im;
	int mouseX() {return _mouseX;}
	int mouseY() {return _mouseY;}
	L_Mouse &mouseE() {return _mouseE;}
	virtual void draw() {Fl_Double_Window::draw(); fl_draw_image(im.buf.data(), 0, 0, im.lx, im.ly, 4);}
	int w() {return Fl_Double_Window::w();}
	int h() {return Fl_Double_Window::h();}
	virtual int handle(int ev);
};
#endif

#ifdef __COMPAT_FLTK__
class L_VentanaImagenFLTK : public L_VentanaImagenAbstracta
{
private:
	L_VentanaImagenFLTK(void); // No valido
	L_VentanaImagenFLTK(const L_VentanaImagenFLTK &other); // No valido
public:
	L_Fl_Double_Window vent;

	L_VentanaImagenFLTK(int lx, int ly, const char *name): vent(lx, ly, name) {}
	L_VentanaImagenFLTK(int x, int y, int lx, int ly, const char *name): vent(x, y, lx, ly, name) {}
	bool dibuja(const L_ImageRGBUchar &im) {vent.im=im; vent.redraw(); return true;}
	bool dibujaRedimensiona(const L_ImageRGBUchar &im) {vent.im=im; vent.resize(vent.x(),vent.y(),im.lx,im.ly); vent.redraw(); return true;}
	bool resize(int w, int h) {vent.resize(vent.x(),vent.y(),w,h); return true;}
	bool resize(int x, int y, int w, int h) {vent.resize(x,y,w,h); return true;}
	int w() {return vent.w();}
	int h() {return vent.h();}
	int mouseX() {return vent.mouseX();}
	int mouseY() {return vent.mouseY();}
	L_Mouse &mouseE() {return vent.mouseE();}
	static void check() {Fl::check();}
	~L_VentanaImagenFLTK() {}
};
#endif //__COMPAT_FLTK__


#if defined (__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
#define __L_VENTANA__
#endif


class L_VentanaImagen : public L_VentanaImagenAbstracta  // Usa OpenCV o FLTK para las ventanas
{
#if defined (__COMPAT_IPLIMAGE__)
  private:  L_VentanaImagenCV vent;
  public:  static void check() {L_VentanaImagenCV::check();}
#elif defined (__COMPAT_FLTK__)
  private:  L_VentanaImagenFLTK vent;
  public:   static void check() {L_VentanaImagenFLTK::check();}
#endif

private:
	L_VentanaImagen(void); // No valido
	L_VentanaImagen(const L_VentanaImagen &other); // No valido

public:

#if defined __L_VENTANA__
	// En este caso, hay soporte para ventanas
	L_VentanaImagen(int lx, int ly, const char *name) : vent(lx, ly, name) {}
	L_VentanaImagen(int x, int y, int lx, int ly, const char *name) : vent(x,y,lx, ly, name) {}
	bool dibuja(const L_ImageRGBUchar &im) {return vent.dibuja(im);}
	bool dibujaRedimensiona(const L_ImageRGBUchar &im) {return vent.dibujaRedimensiona(im);}
	bool resize(int w, int h) {return vent.resize(w, h);}
	bool resize(int x, int y, int w, int h) {return vent.resize(x, y, w, h);}
	int w() {return vent.w();}
	int h() {return vent.h();}
	int mouseX() {return vent.mouseX();}
	int mouseY() {return vent.mouseY();}
	L_Mouse &mouseE() {return vent.mouseE();}
	virtual ~L_VentanaImagen() {}
#else
	// En este caso, no hay soporte para ventanas
	L_VentanaImagen(int lx, int ly, const char *name) {mouseE() = L_MouseLeido; printf("No hay soporte para ventanas: %s\n", name);}
	L_VentanaImagen(int x, int y, int lx, int ly, const char *name) {mouseE() = L_MouseLeido; printf("No hay soporte para ventanas: %s\n", name);}
	bool dibuja(const L_ImageRGBUchar &im) {return false;}
	bool dibujaRedimensiona(const L_ImageRGBUchar &im) {return false;}
	bool resize(int w, int h) {return false;}
	int mouseX() {return -1;}
	int mouseY() {return -1;}
	L_Mouse &mouseE() {return mouseE_int;}
	static void check() {}
	virtual ~L_VentanaImagen() {}
	L_Mouse mouseE_int;
#endif
};




#ifdef __COMPAT_FLTK__
class L_VentMapa3D : public Fl_Gl_Window  // OpenGL creo que no se puede llamar desde OpenCV
{
private:
	L_Pose3D poseCamara; // Oculta para asegurar que mGL este actualizado
	L_CoordsCart3D posLuz;
	L_Pose3D poseRobot;
	bool poseRobotVal;
	L_HomogeneousMatrix mGL; // transformUsing del sistema de la camara (xADE,yIZQ,zARR) al de openGL (xDER,yARR,zATR)
	double maxDistFrustrum;
public:
	enum L_TipoFigura {L_punto, L_cubo, L_esfera, L_sis_coords};

	std::vector<L_Pose3D_cuat> puntos;
	std::vector<long> colores;
	std::vector<L_TipoFigura> tipos;  // 0 = cubo, 1 = pose, 2 = punto
	std::vector<double> radios; // radio de los puntos o poses para dibujarlos
	std::vector<L_CoordsCart3D> lineas; // Tiene que ser par

	std::vector<L_Ellipsoid3D> elipsoides;
	std::vector<long> coloresElipsoides;

	double radioRobot;

	L_VentMapa3D(int W, int H, const char *l=0):Fl_Gl_Window(W,H,l), maxDistFrustrum(400) {mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE); poseRobotVal = false;}
	L_VentMapa3D(int X, int Y, int W, int H, const char *l=0):Fl_Gl_Window(X,Y,W,H,l), maxDistFrustrum(400) {mode(FL_RGB | FL_ALPHA | FL_DEPTH | FL_DOUBLE); poseRobotVal = false;}

	// El sis de ejes de OpenGL es (der, arr, ade). Esta mariz se debe recalcular al cambiar la pose de la camara observadora
	void recalc_mGL() {L_HomogeneousMatrix m1, m2; m1.fijaPose3D_fijo_a_movido(L_Pose3D(0,0,0,-M_PI/2,0,M_PI/2)); m2.fijaPose3D_fijo_a_movido(poseCamara); mGL=m1*m2;}

	void fijaPoseCamara(const L_Pose3D &poseCamara) {this->poseCamara=poseCamara; recalc_mGL();}
	void fijaPoseCamara_abs(const L_CoordsCart3D &cam, const L_CoordsCart3D &obj, const L_CoordsCart3D &arr) {poseCamara.fijaPose_abs(cam.x,cam.y,cam.z,obj.x,obj.y,obj.z,arr.x,arr.y,arr.z); recalc_mGL();}
	void fijaPoseCamara_rel(const L_CoordsCart3D &cam, const L_CoordsCart3D &dirObj, const L_CoordsCart3D &dirArr) {poseCamara.fijaPose_rel(cam.x,cam.y,cam.z,dirObj.x,dirObj.y,dirObj.z,dirArr.x,dirArr.y,dirArr.z); recalc_mGL();}
	L_Pose3D pidePoseCamara() {return poseCamara;}

	void fijaDistanciaMaxima(double dist) {maxDistFrustrum = dist;}

	void fijaPosLuz(const L_CoordsCart3D &luz) {posLuz = luz;}

	void fijaPoseRobot(const L_Pose3D &poseRobot, double radioRob = 5) {this->poseRobot=poseRobot; radioRobot = radioRob; poseRobotVal = true;}
	L_Pose3D pidePoseRobot() {return poseRobot;}

	void clear() {poseRobotVal = false; puntos.resize(0); colores.resize(0); tipos.resize(0); radios.resize(0); elipsoides.resize(0); coloresElipsoides.resize(0); lineas.resize(0);} // Llamar redraw() cuando este listo para redibujar
	void limpiaPuntos() {puntos.resize(0); colores.resize(0); tipos.resize(0); radios.resize(0); elipsoides.resize(0); coloresElipsoides.resize(0);} // Llamar redraw() cuando este listo para redibujar

	void agregaPunto3D_cubo(const L_CoordsCart3D &punto, long color=0xFFFFFFFFL, double radio = 2);
	void agregaPunto3D_esfera(const L_CoordsCart3D &punto, long color=0xFFFFFFFFL, double radio = 2);
	void agregaPunto3D_punto(const L_CoordsCart3D &punto, long color=0xFFFFFFFFL);
	void agregaElipsoide3D(const double *cov3x3, const double *tr3x1, long color=0xFFFFFFFFL);
	void agregaPose3D(const L_Pose3D_cuat &pose, long color=0xFFFFFFFFL, double radio = 20);
	void agregaLinea(const L_CoordsCart3D &ini, const L_CoordsCart3D &fin);
	void agregaLineas(const std::vector<L_CoordsCart3D> &lins); // Lineas pegadas

	void draw();

	void draw_fijaPoseGL(const L_Pose3D &pose); // Debe ser rapido, equivale a hacer un "glLoadIdentity(); glMultMatrixd(pose;poseCamara)" para dibujar con OpenGL
	void draw_fijaPoseGL(const L_Pose3D_cuat &pose); // Debe ser rapido, equivale a hacer un "glLoadIdentity(); glMultMatrixd(pose;poseCamara)" para dibujar con OpenGL
	void draw_fijaPoseGL(const L_HomogeneousMatrix &pose); // Debe ser rapido, equivale a hacer un "glLoadIdentity(); glMultMatrixd(pose;poseCamara)" para dibujar con OpenGL
	void draw_dibCubo(double radio); // Dibuja un cubo centrado en (0,0,0). Llamar a fijaPoseGL() para cambiar la pose en que se va a dibujar.
	void draw_dibParalelepipedo(double x0, double y0, double z0, double x1, double y1, double z1); //x0<x1, y0<y1; z0<z1
	void draw_dibEsfera(double radio, GLUquadricObj *quadratic = NULL); // Dibuja un cubo centrado en (0,0,0). Llamar a fijaPoseGL() para cambiar la pose en que se va a dibujar.
	void draw_dibElipsoide(double rx, double ry, double rz, GLUquadricObj *quadratic = NULL); // Dibuja elipsoide. Llamar a fijaPoseGL() para cambiar la pose en que se va a dibujar.
	void draw_dibLineas(); // No es una primitiva
	void draw_dibRobot(double radio); // Dibuja un robot centrado en (0,0,0). Llamar a fijaPoseGL() para cambiar la pose en que se va a dibujar.
	void draw_dibSisCoords(double radio, double grosor); // Dibuja un sistema de coordenadas centrado en (0,0,0). Llamar a fijaPoseGL() para cambiar la pose en que se va a dibujar.

	void check() {Fl::check();}

	void test();
};
#endif //__COMPAT_FLTK__


class L_Avi
{
	FILE *fp;
	int flujo; // 0 = lectura, 1 = escritura

	L_Avi() {fp = NULL; flujo = -1;}
	
	// Escritura
	bool escribirHeader() {throw L_NonImplementedException();}
	bool escribirCuadro(L_ImageRGBUchar &im) {throw L_NonImplementedException();}
	bool escribirSonido(L_SignalDouble &wav) {throw L_NonImplementedException();}

	// Lectura
	int tipoChunk() {throw L_NonImplementedException();} // Para ver que se viene al read.  0 = read imagen, 1 = read sonido
	bool leerCuadro(L_ImageRGBUchar &im) {throw L_NonImplementedException();}
};

///

class L_KdNodo3D
{
public:
	L_KdNodo3D *der;
	L_KdNodo3D *izq;
	double val[3]; // Punto 3D (para los nodos hoja solamente)
	int d; // Dimension del corte (para los nodos interiores solamente)
	double umbral; // value del corte (para los nodos interiores solamente)
	int i; // Indice del elemento apuntado en el arreglo que se le paso (no quiero usar punteros aca)

	L_KdNodo3D() {der=NULL; izq=NULL; d=-1; umbral = 0; val[0] = 0; val[1] = 0; val[2]= 0; i=-1;}

	static L_KdNodo3D *creaRec(L_KdNodo3D *arr, int nTot);
	static L_KdNodo3D *creaNoRec(L_KdNodo3D *arr, int nTot);
	void destroyRec();
	void destruyeNoRec();

	static int cmpRepetidos(const void *kdNodo3D_1, const void *kdNodo3D_2);
	static int cmpDimension(const void *kdNodo3D_1, const void *kdNodo3D_2);

	~L_KdNodo3D() {}
};

class L_KdTree3D
{
public:
	L_KdNodo3D *root;
	int nHojas;
	int nElim;

	bool usePercentage; // Usar porcRev o nCompMax como parametro principal
	int nCompMax; // Parametro: nº máximo de comparaciones a realizar
	double porcRev; // Parametro: porcentaje de la base de datos que se revisan

	L_KdTree3D() {root = NULL; nHojas=0; nElim=0; usePercentage=true;nCompMax=200;porcRev=0.25;}

	bool createFrom(const std::vector<L_CoordsCart3D> &arr);
	void buscaMasCercanos(L_CoordsCart3D &pun, int n, int *indexes, double *dist);
	int buscaMasCercano(L_CoordsCart3D &pun);
	void destroy() {if (root != NULL) {root->destruyeNoRec(); delete root; root=NULL;}}
	~L_KdTree3D() {destroy();}
};

class L_ArregloP3D
{
public:
	std::vector<L_CoordsCart3D> arr;

	bool createFrom(const std::vector<L_CoordsCart3D> &arr) {this->arr = arr; return true;}
	void buscaMasCercanos(L_CoordsCart3D &pun, int n, int *indexes, double *dist);
	int buscaMasCercano(L_CoordsCart3D &pun);
	void destroy() {arr.clear();}
};

#define MAXLASERPOINTS (361)

struct L_LaserDato
{
	L_uint32 tiempo;
	L_uint16 cx;
	L_uint16 rayos[MAXLASERPOINTS];
	L_uint16 crc;
	L_uint16 status;
};

class L_CurvaPuntos
{
public: // No tiene mas miembros, solo funciones extra
	std::vector<L_CoordsCart3D> v;

	L_CurvaPuntos() {}
	L_CurvaPuntos(int n) : v(n) {}

	L_CoordsCart3D &operator[](int n) {return v[n];}
	const L_CoordsCart3D &operator[](int n) const {return v[n];}

	size_t size() const {return v.size();}
	void resize(size_t s) {v.resize(s);}

	void swap(L_CurvaPuntos &other) {v.swap(other.v);}

	bool leerBarridoTxt(FILE *fp, int nPuntos=361);
	bool guardarBarridoTxt(FILE *fp);
	bool guardarBarridoTxt2col(FILE *fp);
	void suavizar(L_CurvaPuntos &resultado, double sigma);
	void calcularNormales2D(L_CurvaPuntos &normales, L_KdTree3D &arb); // Calcular normales cuando los datos son planes
	void OP_mult(L_HomogeneousMatrix &H, L_CurvaPuntos &curva);
};

class L_ArregloCurvaPuntos : public L_Array<L_CurvaPuntos>
{
public:
	L_ArregloCurvaPuntos() : L_Array<L_CurvaPuntos>() {}
	L_ArregloCurvaPuntos(int n, int growthMem = 50) : L_Array<L_CurvaPuntos>(n, growthMem) {}
	bool leerTxt(FILE *fp, int nPuntos=361);
	bool guardarTxt(FILE *fp);
	bool guardarTxt2col(FILE *fp);
	static bool transformar_dat_a_txt(const char *nombreDat);
};

// Interfaz para archivos de laser .dat
class L_ArregloLaserDato : public L_Array<L_LaserDato>
{
public:
	void leerDat(FILE *fp);
	void grabarDat(FILE *fp);
	void guardarPuntosEn(L_CurvaPuntos &arr, int nBarrido);
	void copiarPuntosEn(L_ArregloCurvaPuntos &arr);
};

struct L_ICP_Dato
{
	std::vector<L_CoordsCart3D> *arr1; // Datos del escaneo 1
	std::vector<L_CoordsCart3D> *arr2; // Datos del escaneo 2
	std::vector<L_CoordsCart3D> *per2; // Perpendiculares al escaneo 2
	L_KdTree3D *arb2; // Kdtree del escaneo 2 (para encontrar puntos cercanos)
	std::vector<int> *asoc;
	double distMax;
};

class L_ICP
{
public:
	// Funciones de error para metodos de optimizacion iterativa
	static void fnObjetivo_calcAsoc(const void *datos, L_Matrix &xytheta);
	static void fnObjetivo_gradiente(const void *datos, const L_Matrix &xytheta, L_Matrix &errVect);
	static void fnObjetivo_puntos(const void *datos, const L_Matrix &xytheta, L_Matrix &errVect);
	// Metodos iterativos sin minimizacion iterativa
	static double ICP_lineal_norm_2d(L_ICP_Dato *datos, L_Matrix &xytheta);
	static double ICP_lineal_norm_3d(L_ICP_Dato *datos, L_Matrix &xytheta); // Para extension futura
	static double ICP_lineal_pun_2d(L_ICP_Dato *datos, L_Matrix &xytheta);
	// La funcion que se llama desde fuera
	static double calcICP(L_CurvaPuntos &inic, L_CurvaPuntos &final, double dMax, L_Pose3D &pose, int nMetodo = 0, L_CurvaPuntos *normales=NULL); // Lo que hay

	static int main_icp(int argc, char *argv[]);
};

#define VERIF_INT // Verificar integridad del codeMapping de L_Fnes_Pato.cpp en memoria. Para que un stack overflow no se coma el codeMapping...



class L_ProgramasFnesPato
{
public:
	static int main_L_for(int argc, char *argv[], bool debug=false);
	static int main_L_encriptado(int argc, char *argv[], char *sal1=NULL, char *sal2=NULL);
	static int main_L_bmplz(int argc, char *argv[]); // Usa readComprRGB() y writeComprRGB() para manipular bitmaps .bmplz  (bmz requiere trasponer la imagen)
	static int main_pruebaMitos(int ntest = -1);
	static int main_capturaBmplz(int argc, char *argv[]); // Usa readComprRGB() y writeComprRGB() para manipular bitmaps .bmplz
	static int main_capturaPNG(int argc, char *argv[]);

	static int main_pruebaClases(int argc, char *argv[]);
	static int main_reemplazar(int argc, char *argv[]);
	static int main_imprimirArgv(int argc, char *argv[]);
	static int main_interpolarVideo(int argc, char *argv[]);
	static int main_rectificar(int argc, char *argv[]);

#ifdef __COMPAT_IPLIMAGE__
	static int main_L_GrabaFramesOpenCV(int argc, char *argv[]);
	static int main_L_bmplzavi(int argc, char *argv[]); // bmz a avi y viceversa. Requiere OpenCV 2.0+ para usar codecs
#endif

#if defined(__COMPAT_IPLIMAGE__)
	static int main_L_capturaCV(int argc, char *argv[]); // bmz a avi y viceversa. Requiere OpenCV 2.0+ para usar codecs
#endif

// Calibrador distorsion manual
#if defined (__COMPAT_FLTK__)
	static int main_L_CalibraDistorsionRadialManual(int argc, char *argv[]);
	static int main_dib3d(int argc, char *argv[]);
	static int main_muestra_plpts(int argc, char *argv[]);
#endif

	static int main_muestraBmz(int argc, char *argv[]);

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
	static int main_pruebaFourier(int argc, char *argv[]);
	static int main_pruebaHough(int argc, char *argv[]);
	static int main_pruebaGabor(int argc, char *argv[]);
	static int main_pruebaGaborUmbral(int argc, char *argv[]);
	static int main_pruebaMeanShift(int argc, char *argv[]);
	static int main_pruebaPiel(int argc, char *argv[]);
	static int main_fractal(int argc, char *argv[]);
#endif
};


class L_ProgramasFnesPato_prueba
{
public:

	static int main_infoMutua(int argc, char *argv[]);
	static int main_infoMutuaRotEsc(int argc, char *argv[]);
};


// Algo gracioso
// #define STRING char *
// #define IF if(
// #define THEN ){
// #define ELSE } else {
// #define FI ;}
// #define WHILE while {
// #define DO ){
// #define od ;}
// #define INT int
// #define BEGIN {
// #define END }

// Otra cosa curiosa:  #if 0  ... #endif  se puede usar para comentar codeMapping

#endif //__L_FNES_PATO_H_



