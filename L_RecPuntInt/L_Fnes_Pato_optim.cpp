#include "L_Fnes_Pato.h"
#include <cfloat>

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


#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif




#ifdef __USAR_L_LAPACK_H__
#include "L_Lapack.hpp"  // Es una buena biblioteca. Pide que lo incluyan en los paper: Grosse-Kunstleve RW, Terwilliger TC, Adams PD: Experience converting a large Fortran-77 program to C++. Newsletter of the IUCr Commission on Crystallographic Computing 2009, 10, 75-84. (Reprint)
#endif // __USAR_L_LAPACK_H__

//#include "fnes_difgaus.h"

// Incluye manipulacion de: numeros complejos, fracciones, listas, arboles, compresion Huffman,
// arreglos, strings, parametros (servidor de parametros), imagenes, matrices y modelos de camaras
// Incluye interfaces a ATL, OpenCV, SVS si se usan "defines" adecuados
// Nota: Las imagenes son representadas como [columna][fila], es decir, se almacenan como columnas consecutivas
// Compatibilidad: Ansi C++
// Ojo: Poner acentos en los comentarios es anti ANSI
// Autor: Patricio Loncomilla (ploncomi@gmail.com)




#ifdef __COMPAT_UBLAS__
namespace boost { namespace numeric { namespace ublas {
    template<class T>
    class L_readonly_array_adaptor:
        public storage_array<L_readonly_array_adaptor<T> > {

        typedef L_readonly_array_adaptor<T> self_type;
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;
        typedef T value_type;
        typedef const T &const_reference;
        typedef const T *const_pointer;
    public:

        // Construction and destruction
        BOOST_UBLAS_INLINE
        L_readonly_array_adaptor ():
            size_ (0), data_ (0) {
        }
        BOOST_UBLAS_INLINE
        L_readonly_array_adaptor (size_type size, const_pointer data):
            size_ (size), data_ (data) {
        }
        BOOST_UBLAS_INLINE
        ~L_readonly_array_adaptor () {
        }

        L_readonly_array_adaptor (const L_readonly_array_adaptor& rhs)
          : size_(rhs.size_), data_(rhs.data_)
        { }

        // Resizing
        BOOST_UBLAS_INLINE
        void resize (size_type size) {
            size_ = size;
        }
        BOOST_UBLAS_INLINE
        void resize (size_type size, const_pointer data) {
            size_ = size;
            data_ = data;
        }

        // Random Access Container
        BOOST_UBLAS_INLINE
        size_type max_size () const {
            return std::numeric_limits<size_type>::max ();
        }
        
        BOOST_UBLAS_INLINE
        bool empty () const {
            return size_ == 0;
        }
            
        BOOST_UBLAS_INLINE
        size_type size () const {
            return size_;
        }

        // Element access
        BOOST_UBLAS_INLINE
        const_reference operator [] (size_type i) const {
            BOOST_UBLAS_CHECK (i < size_, bad_index ());
            return data_ [i];
        }

        // Iterators simply are pointers.
        typedef const_pointer const_iterator;

        BOOST_UBLAS_INLINE
        const_iterator begin () const {
            return data_;
        }
        BOOST_UBLAS_INLINE
        const_iterator end () const {
            return data_ + size_;
        }

        // this typedef is used by vector and matrix classes
        typedef const_pointer iterator;

        // Reverse iterators
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;

        BOOST_UBLAS_INLINE
        const_reverse_iterator rbegin () const {
            return const_reverse_iterator (end ());
        }
        BOOST_UBLAS_INLINE
        const_reverse_iterator rend () const {
            return const_reverse_iterator (begin ());
        }

    private:
        size_type size_;
        const_pointer data_;
    };

    /// converts a chunk of memory into a (readonly) usable ublas vector.
    template <class T>
    vector<const T, L_readonly_array_adaptor<T> >
    L_make_vector_from_pointer(const size_t size, const T * data)
    {
        typedef L_readonly_array_adaptor<T> a_t;
        typedef vector<const T, a_t>      v_t;
        return v_t(size, a_t(size, data));
    }

    // converts a chunk of memory into a (readonly) usable dense matrix.
    template <class LAYOUT, class T>
    matrix<const T, LAYOUT, L_readonly_array_adaptor<T> >
    L_make_matrix_from_pointer(const size_t size1, const size_t size2, const T * data)
    {
        typedef L_readonly_array_adaptor<T> a_t;
        typedef matrix<const T, LAYOUT, a_t>      m_t;
        return m_t(size1, size2, a_t(size1*size2, data));
    }
    // default layout: row_major
    template <class T>
    matrix<const T, row_major, L_readonly_array_adaptor<T> >
    L_make_matrix_from_pointer(const size_t size1, const size_t size2, const T * data)
    {
	  return L_make_matrix_from_pointer<row_major>(size1, size2, data);
    }

    // brief converts a C-style 2D array into a (readonly) usable dense matrix.
    template <class T, size_t M, size_t N>
    matrix<const T, row_major, L_readonly_array_adaptor<T> >
    L_make_matrix_from_pointer(const T (&array)[M][N])
    {
        typedef L_readonly_array_adaptor<T> a_t;
        typedef matrix<const T, row_major, a_t>      m_t;
        return m_t(M, N, a_t(M*N, array[0]));
    }
    template <class T, size_t M, size_t N>
    matrix<const T, row_major, L_readonly_array_adaptor<T> >
    L_make_matrix_from_pointer(const T (*array)[M][N])
    {
        typedef L_readonly_array_adaptor<T> a_t;
        typedef matrix<const T, row_major, a_t>      m_t;
        return m_t(M, N, a_t(M*N, (*array)[0]));
    }

}}}  // namespace boost { namespace numeric { namespace ublas {
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
// Uses lu_factorize and lu_substitute in uBLAS to invert a matrix
namespace boost { namespace numeric { namespace ublas {
template<class T>
bool L_InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = (int)lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}



// \brief decompose the symmetric positive definit matrix A into product L L^T.
//
// \param MATRIX type of input matrix 
// \param TRIA type of lower triangular output matrix
// \param A square symmetric positive definite input matrix (only the lower triangle is accessed)
// \param L lower triangular output matrix 
// \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
//
template < class MATRIX, class TRIA >
size_t cholesky_decompose(const MATRIX& A, TRIA& L)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  assert( A.size1() == A.size2() );
  assert( A.size1() == L.size1() );
  assert( A.size2() == L.size2() );

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
        
    double qL_kk = A(k,k) - inner_prod( project( row(L, k), range(0, k) ),
                                        project( row(L, k), range(0, k) ) );
    
    if (qL_kk <= 0) {
      return 1 + k;
    } else {
      double L_kk = sqrt( qL_kk );
      L(k,k) = L_kk;
      
      matrix_column<TRIA> cLk(L, k);
      project( cLk, range(k+1, n) )
        = ( project( column(A, k), range(k+1, n) )
            - prod( project(L, range(k+1, n), range(0, k)), 
                    project(row(L, k), range(0, k) ) ) ) / L_kk;
    }
  }
  return 0;      
}


// \brief decompose the symmetric positive definit matrix A into product L L^T.
//
// \param MATRIX type of matrix A
// \param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
// \param A output: the lower triangle of A is replaced by the cholesky factor
// \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
//
template < class MATRIX >
size_t cholesky_decompose(MATRIX& A)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  const MATRIX& A_c(A);

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
        
    double qL_kk = A_c(k,k) - inner_prod( project( row(A_c, k), range(0, k) ),
                                          project( row(A_c, k), range(0, k) ) );
    
    if (qL_kk <= 0) {
      return 1 + k;
    } else {
      double L_kk = sqrt( qL_kk );
      
      matrix_column<MATRIX> cLk(A, k);
      project( cLk, range(k+1, n) )
        = ( project( column(A_c, k), range(k+1, n) )
            - prod( project(A_c, range(k+1, n), range(0, k)), 
                    project(row(A_c, k), range(0, k) ) ) ) / L_kk;
      A(k,k) = L_kk;
    }
  }
  return 0;      
}

// \brief decompose the symmetric positive definit matrix A into product L L^T.
//
// \param MATRIX type of matrix A
// \param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
// \param A output: the lower triangle of A is replaced by the cholesky factor
// \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
//
template < class MATRIX >
size_t incomplete_cholesky_decompose(MATRIX& A)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  // read access to a const matrix is faster
  const MATRIX& A_c(A);

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
    
    double qL_kk = A_c(k,k) - inner_prod( project( row( A_c, k ), range(0, k) ),
                                          project( row( A_c, k ), range(0, k) ) );
    
    if (qL_kk <= 0) {
      return 1 + k;
    } else {
      double L_kk = sqrt( qL_kk );

      // aktualisieren
      for (size_t i = k+1; i < A.size1(); ++i) {
        T* Aik = A.find_element(i, k);

        if (Aik != 0) {
          *Aik = ( *Aik - inner_prod( project( row( A_c, k ), range(0, k) ),
                                      project( row( A_c, i ), range(0, k) ) ) ) / L_kk;
        }
      }
        
      A(k,k) = L_kk;
    }
  }
        
  return 0;
}




}}} // namespace boost { namespace numeric { namespace ublas {
#endif // __COMPAT_UBLAS__















int L_ComplexDouble::solveCubicEquation(double a, double b, double c, double d, L_ComplexDouble &x1, L_ComplexDouble &x2, L_ComplexDouble &x3) // Resuelve ax^3 + x + d = 0
{
	double ainv, p, q, p3, q2, dM;
	L_ComplexDouble A, B, UUU[2], U[6], X[6];
	int i, iM, n;
	// Cambio a la forma x^3 + a*x^2 + b*x + c = 0
	ainv = 1/a;
	a = b*ainv;
	b = c*ainv;
	c = d*ainv;
	// Cambio a la forma t^3 + p*t + q = 0
	p = b - a*a/3;
	q = c + (2*a*a*a -9*a*b) / 27;
	p3 = p / 3;
	q2 = q / 2;
	p3 = p3*p3*p3;
	q2 = q2*q2;
	A.define(-q/2,0); // Cardano: u + v = t  y   3*u*v + p = 0  =>  permite despejar u
	B.define(q2+p3, 0); // Soluciones para u: raices_cubicas( A + raices_cuadradas(B) )
	B.calcSqrt(); // Las raices son B y -B
	UUU[0] = A+B;
	UUU[1] = A-B;
	UUU[0].calcCubicRoot(U[0],U[1],U[2]); // Las tres raices cubicas
	UUU[1].calcCubicRoot(U[3],U[4],U[5]); // Las tres raices cubicas
	for (n=0; n<6; n++)
		X[n] = U[n] - p / (3*U[n]) -a/3; // De estas 6 soluciones para x, solo 3 son distintas
	
	for (n=0; n<3; n++)
	{
		// Busco el mas cercano a 2*n -> dejarlo en 2*n+1
		iM = 2*n+1;
		dM = X[2*n].dist(X[2*n+1]);

		for (i=2*n+2; i<6; i++)
		{
			d = X[2*n].dist(X[i]);
			if (d < dM)
			{
				iM = i;
				dM = d;
			}
		}
		A = X[iM];
		X[iM] = X[2*n+1];
		X[2*n+1] = A;
	}
	x1 = X[0]; // Las tres soluciones distintas
	x2 = X[2];
	x3 = X[4];

	return 3;
}

int L_ComplexDouble::solveQuarticEquation(double a, double b, double c, double d, double e, L_ComplexDouble &x1, L_ComplexDouble &x2, L_ComplexDouble &x3, L_ComplexDouble &x4)
{
	if (a==0)
		return -1;
	double a2 = a*a;
	double b2 = b*b;
	double alfa = -3.0/8.0 * b2/(a2) + c/a;
	double beta = 1.0/8 * b2*b/(a2*a) - b*c/(2*a2) + d/a;
	double gamma = -3.0/256 * b2*b2/(a2*a2) + 1.0/16 * c*b2/(a2*a) - 1/4.0 * b*d/(a2) + e/a;
	double P = -alfa*alfa/12.0 - gamma;
	double Q = -alfa*alfa*alfa/108 + alfa*gamma/3.0 - beta*beta/8.0;
	L_ComplexDouble R = -Q/2 + L_ComplexDouble::csqrt(L_ComplexDouble(Q*Q/4.0 + P*P*P/27.0, 0));
	L_ComplexDouble U = L_ComplexDouble::ccurt(R);

	L_ComplexDouble y;
	if (U.re == 0 && U.im == 0)
		y = -5.0/6 * alfa - L_ComplexDouble::ccurt(L_ComplexDouble(Q,0));
	else
		y = -5.0/6 * alfa + U - P/(3*U);
	L_ComplexDouble W = L_ComplexDouble::csqrt(alfa + 2*y);

	if (W.re == 0 && W.im == 0)
		return -1;

	x1 = -b/(4*a) + 0.5 * ( +W - L_ComplexDouble::csqrt(-(3*alfa+2*y + 2*beta/W)) );
	x2 = -b/(4*a) + 0.5 * ( +W + L_ComplexDouble::csqrt(-(3*alfa+2*y + 2*beta/W)) );
	x3 = -b/(4*a) + 0.5 * ( -W - L_ComplexDouble::csqrt(-(3*alfa+2*y - 2*beta/W)) );
	x4 = -b/(4*a) + 0.5 * ( -W + L_ComplexDouble::csqrt(-(3*alfa+2*y - 2*beta/W)) );

	return 4;
}

void L_Spline1D::compute()
{
#ifdef L_BASIC_DEBUG
	if (n<=0)
		L_hard_shutdown("L_Spline1D::compute : non-positive size allocation");
#endif
	std::vector<double> u(n);
	double sig, p;
	int i;
	y2[0]=0;
	u[0]=0;
	for (i=1; i<n-1; i++)
	{
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=1/(sig*y2[i-1]+2);
		y2[i]=(sig-1.0)*p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])*p;
	}
	y2[n-1]=0;
	u[n-1]=0;
	for (i=n-2;i>=0;i--)
		y2[i]=y2[i]*y2[i+1]+u[i];
}


double L_Spline1D::evaluate(double xVal)
{
	double h, a, b;
	int k,k1,k2;
	if (xEntero)
	{
		k1=L_FLOOR(xVal+x[0]);
		if (k1<0)
			k1=0;
		if (k1>=n-2)
			k1=n-2;
		k2=k1+1;
	}
	else
	{
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
	}
	h=x[k2]-x[k1];
	a=(x[k2]-xVal)/h;
	b=(xVal-x[k1])/h;
	return a*y[k1]+b*y[k2]+((a*a*a-a)*y2[k1]+(b*b*b-b)*y2[k2])*(h*h)/6.0;
}



//-----------------------------------------
//-----------------------------------------
//  Funciones tomadas de "Recipes in C"
//-----------------------------------------
//-----------------------------------------
#if defined(MAXIT) || defined(EPS) || defined(FPMIN)
	#error #define conflicts
#endif
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double NR_gammln(double xx)
//! Recipes in C: Returns the value ln[gamma(xx)] for xx > 0.
{
	//Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
	//accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	if (xx<1.0e10*DBL_MIN) // Para evitar que se caiga
		xx=1.0e10*DBL_MIN;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double NR_betai(double a, double b, double x)
//! Recipes in C: Returns the incomplete beta function Ix(a, b).
{
	double NR_betacf(double a, double b, double x);
	double NR_gammln(double xx);
	double bt;
	throw_L_ArgException_if(x < 0.0 || x > 1.0, "betai"); //Bad x in routine betai
	if (x == 0.0 || x == 1.0) bt=0.0;
	else //Factors in front of the continued fraction.
		bt=exp(NR_gammln(a+b)-NR_gammln(a)-NR_gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) //Use continued fraction directly.
		return bt*NR_betacf(a,b,x)/a;
	else //Use continued fraction after making the symmetry transformation.
		return 1.0-bt*NR_betacf(b,a,1.0-x)/b;
}

double NR_betacf(double a, double b, double x)
//! Recipes in C: Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz’s method (§5.2).
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b; //These q’s will be used in factors that occur in the coefficients (6.4.6). qap=a+1.0;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; //First step of Lentz’s method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; //One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break; //Are we done?
	}
	throw_L_ArgException_if (m > MAXIT, "betacf"); // a or b too big, or MAXIT too small in betacf
	return h;
}

double NR_gammp(double a, double x)
{
	double gamser,gammcf,gln;

	throw_L_ArgException_if (x < 0.0 || a <= 0.0, "gammp() : argumentos no validos");
	if (x < (a+1.0))
	{
		NR_gser(&gamser,a,x,&gln);
		return gamser;
	}
	else
	{
		NR_gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

double NR_gammq(double a, double x)
{
	double gamser,gammcf,gln;
	
	throw_L_ArgException_if (x < 0.0 || a <= 0.0, "gammq() : argumentos no validos");
	if (x < (a+1.0))
	{
		NR_gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	}
	else
	{
		NR_gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

void NR_gcf(double* gammcf, double a, double x, double* gln)
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;

	*gln=NR_gammln(a);
	a1=x;
	for (n=1;n<=MAXIT;n++)
	{
		an=(double) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1)
		{
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS)
			{
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
	printf(" a too large, NUMREC_ITMAX too small in routine GCF\n");
}

void NR_gser(double* gamser, double a, double x, double* gln)
{
	int n;
	double sum,del,ap;

	*gln=NR_gammln(a);
	if (x <= 0.0)
	{
		if (x < 0.0)
		{
			printf(" x less than 0 in routine GSER\n");
			*gamser=0.0;
			return;
		}
	}

	else
	{
		ap=a;
		del=sum=1.0/a;

		for (n=1;n<=MAXIT;n++)
		{
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS)
			{
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf(" a too large, NUMREC_ITMAX too small in routine GSER\n");
		return;
	}
}

#undef MAXIT
#undef EPS
#undef FPMIN




#define NRR_IM1   2147483563
#define NRR_IM2   2147483399
#define NRR_AM    (1.0/NRR_IM1)
#define NRR_IMM1  (NRR_IM1-1)
#define NRR_IA1   40014
#define NRR_IA2   40692
#define NRR_IQ1   53668
#define NRR_IQ2   52774
#define NRR_IR1   12211
#define NRR_IR2   3791
#define NRR_NTAB  32
#define NRR_NDIV  (1+NRR_IMM1/NRR_NTAB)
#define NRR_EPS   1.2e-7
#define NRR_RNMX  (1.0 - NRR_EPS)

double NR_ran2(long *idum) // float
{
	int j;
	long k;
	static long idum2;                            // initialized below
	static long iy = 0;
	static long iv[NRR_NTAB];
	double temp; // float

	if (*idum <= 0)
	{                             // initialize
		if (-(*idum) < 1)                           // prevent idum == 0
			*idum = 1;
		else
			*idum = -(*idum);                         // make idum positive
		for (j = NRR_NTAB + 7; j >= 0; j--)
		{           // load the shuffle table
			k = (*idum) / NRR_IQ1;
			*idum = NRR_IA1 * (*idum - k*NRR_IQ1) - k*NRR_IR1;
			if (*idum < 0)
				*idum += NRR_IM1;
			if (j < NRR_NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
		idum2 = iv[NRR_NTAB/2];      // Added for urbcsp so that a negative
	}                          // idum always starts the same sequence.

	k = (*idum) / NRR_IQ1;
	*idum = NRR_IA1 * (*idum - k*NRR_IQ1) - k*NRR_IR1;
	if (*idum < 0)
		*idum += NRR_IM1;
	k = idum2/NRR_IQ2;
	idum2 = NRR_IA2 * (idum2 - k*NRR_IQ2) - k*NRR_IR2;
	if (idum2 < 0)
		idum2 += NRR_IM2;
	j = iy / NRR_NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if (iy < 1)
		iy += NRR_IMM1;
	if ((temp = NRR_AM * iy) > NRR_RNMX)
		return NRR_RNMX;                                // avoid endpoint
	else
	return temp;
}

double L_Rand::_ran2() // float
{
	int j;
	long k;
	double temp; // float

	throw_L_ArgException_if(numCalls == -1, "L_Rand::ran2() : non initialized (or overflow in state)");
	throw_L_ArgException_if(seed <= 0 && numCalls != 0, "L_Rand::ran2() : incorrect initialization");

	if (seed <= 0)
	{
		seed0 = seed;
		if (-(seed) < 1) 
			seed = 1;
		else
			seed = -(seed);
		for (j = NRR_NTAB + 7; j >= 0; j--)
		{
			k = (seed) / NRR_IQ1;
			seed = NRR_IA1 * (seed - k*NRR_IQ1) - k*NRR_IR1;
			if (seed < 0)
				seed += NRR_IM1;
			if (j < NRR_NTAB)
				iv[j] = seed;
		}
		iy = iv[0];
		idum2 = iv[NRR_NTAB/2];
	}

	k = (seed) / NRR_IQ1;
	seed = NRR_IA1 * (seed - k*NRR_IQ1) - k*NRR_IR1;
	if (seed < 0)
		seed += NRR_IM1;
	k = idum2/NRR_IQ2;
	idum2 = NRR_IA2 * (idum2 - k*NRR_IQ2) - k*NRR_IR2;
	if (idum2 < 0)
		idum2 += NRR_IM2;
	j = iy / NRR_NDIV;
	iy = iv[j] - idum2;
	iv[j] = seed;
	if (iy < 1)
		iy += NRR_IMM1;

	numCalls++;
	if (numCalls > 1000000)
	{
		numCalls = 0;
		numCalls2++;
	}

	if ((temp = NRR_AM * iy) > NRR_RNMX)
		return NRR_RNMX;
	else
	return temp;
}

#undef NRR_IM1
#undef NRR_IM2
#undef NRR_AM
#undef NRR_IMM1
#undef NRR_IA1
#undef NRR_IA2
#undef NRR_IQ1
#undef NRR_IQ2
#undef NRR_IR1
#undef NRR_IR2
#undef NRR_NTAB
#undef NRR_NDIV
#undef NRR_EPS
#undef NRR_RNMX


int L_Int::cmp(const void *inta, const void *intb)
{
	return ( (*(int*)inta) > (*(int*)intb) ) - ( (*(int*)inta) < (*(int*)intb) );
}

int L_Int::cmpinv(const void *inta, const void *intb)
{
	return ( (*(int*)inta) < (*(int*)intb) ) - ( (*(int*)inta) > (*(int*)intb) );
}

int L_Double::cmp(const void *inta, const void *intb)
{
	return ( (*(double*)inta) > (*(double*)intb) ) - ( (*(double*)inta) < (*(double*)intb) );
}

int L_Double::cmpinv(const void *inta, const void *intb)
{
	return ( (*(double*)inta) > (*(double*)intb) ) - ( (*(double*)inta) < (*(double*)intb) );
}


double L_Spline2DSlow::evaluate(double x, double y)
{
	if (x!=xAux)
	{
		int i;
		splAux.fillIndexes();
		for (i=0; i<ly; i++)
			splAux.y[i]=splArr[i].evaluate(x); // fill intensidad
		splAux.compute();
		xAux=x;
	}
	return splAux.evaluate(y);
}



int L_IndOrd::cmp(const void *indOrdA, const void *indOrdB)
{
	const L_IndOrd *a = (const L_IndOrd *)indOrdA;
	const L_IndOrd *b = (const L_IndOrd *)indOrdB;
	return (a->d > b->d) - (a->d < b->d);
}

int L_IndOrdi::cmp(const void *indOrdiA, const void *indOrdiB)
{
	const L_IndOrdi *a = (const L_IndOrdi *)indOrdiA;
	const L_IndOrdi *b = (const L_IndOrdi *)indOrdiB;
	return (a->val > b->val) - (a->val < b->val);
}

int L_IndOrdL::cmp(const void *indOrdiA, const void *indOrdiB)
{
	const L_IndOrdL *a = (const L_IndOrdL *)indOrdiA;
	const L_IndOrdL *b = (const L_IndOrdL *)indOrdiB;
	return (a->val > b->val) - (a->val < b->val);
}
/////









///

void L_CountingTree::addVoteToFrequencies(int nobj, int votos)
{
	L_CountingTreeNode **pptr;
	pptr=&root;
	while (*pptr!=NULL)
	{
		if (nobj==(*pptr)->nobj)
		{
			(*pptr)->nvot+=votos;
			nVotTot+=votos;
			return;
		}
		if (nobj<(*pptr)->nobj)
			pptr=&(*pptr)->left;
		else
			pptr=&(*pptr)->right;
	}
	*pptr=new L_CountingTreeNode;
	(*pptr)->nobj=nobj;
	(*pptr)->nvot=votos;
	nVotTot+=votos;
	nCellsTot++;
}

int L_CountingTree::votosDe(int nobj)
{
	L_CountingTreeNode *ptr;
	ptr=root;
	while (ptr!=NULL)
	{
		if (nobj==ptr->nobj)
		{
			return ptr->nvot;
		}
		if (nobj<ptr->nobj)
			ptr=ptr->left;
		else
			ptr=ptr->right;
	}
	return 0;
}



int L_HuffmanTreeNode::cmp(const void *p1, const void *p2)
{
#ifdef LAHN
	#error #define conflicts
#endif
#define LAHN(x) ( (L_HuffmanTreeNode  *)x )
	if (LAHN(p1)->peso > LAHN(p2)->peso)
		return 1;
	else if (LAHN(p1)->peso < LAHN(p2)->peso)
		return -1;
	else
		return 0;
#undef LAHN
}

int L_HuffmanTreeNode::cmp2(const void *p1, const void *p2)
{
#ifdef LAHN
	#error #define conflicts
#endif
#define LAHN(x) (* (L_HuffmanTreeNode  **)x )
	if (LAHN(p1)->peso > LAHN(p2)->peso)
		return 1;
	else if (LAHN(p1)->peso < LAHN(p2)->peso)
		return -1;
	else
		return 0;
#undef LAHN
}



void L_HuffmanTree::build_tree_using_stored_frequeencies()
{
	// El nuevo codeMapping es comprensible, el inicial es realmente horrible
	L_HuffmanTreeList lista1;
	L_HuffmanTreeList lista2;
	L_HuffmanTreeNode *ptr, *menor;
	int i;
	L_HuffmanTreeNode *arr[256];
	for (i=0; i<256; i++)
	{
		arr[i]=new_L_HuffmanTreeNode();
		arr[i]->peso=pesos[i];
		arr[i]->value=(unsigned char)i;
		children[i] = arr[i];
	}
	qsort(arr,256,sizeof(arr[0]),&L_HuffmanTreeNode::cmp2);
	// La idea es:
	// - Se tiene una lista que se mantiene siempre ordenada (el nodo con menor frecuencia al inicio)
	// - Se sacan los dos nodos A y B con menor frecuencia
	// - Se build un nuevo nodo C y se le agregan los children A y B
	// - La frecuencia del nodo C es la sum de las frecuencias de A y B
	// - Se addFrom el nodo C a la lista y se reordena
	//
	// Para hacerlo eficiente se usan dos listas (lista1 y lista2)
	// - lista1 contiene solo nodos aislados (inicialmente vacia)
	// - lista2 contiene solo nodos padres (inicialmente vacia)
	// - Se agregan los nodos iniciales (aislados) a lista1 y se ordenan
	// - Se sacan los dos menores nodos (revisando las dos listas, se prefiere lista1 para minimizar la variance)
	// - Se forma un nodo padre con esos 2 children
	// - se addFrom el nodo nuevo a la lista2
	// - De este modo, las dos listas se mantienen siemre ordenadas, es mas rapido asi
	// - Se termina cuando queda solo 1 nodo

	// Hay dos listas: lista1 y lista2

	// Agregar nodos a lista1
	for (i=0 ; i<256 ; i++)
		lista1.push_back_node(arr[i]); // Esta lista guarda los objetos sin modificar su contenido (salvo ->sig)

	// Mantener la generacion de padres mientras hasta que lista1.size()==0 y lista2.size()==1 (queda solo 1 nodo)
	while( lista1.n!=0 || lista2.n!=1 )
	{
		ptr=new_L_HuffmanTreeNode(); // El nuevo papa, parte con sub-punteros nulos

		// Sacar el menor de las dos listas
		L_HuffmanTreeList_pop_minor_element(menor, lista1, lista2); //menor = L_HuffmanTreeList::pop_minor_element(lista1, lista2);
		// Agregar el children izquierdo
		ptr->hijoIzq=menor;
		ptr->hijoIzq->padre=ptr;

		// Sacar el menor de las dos listas
		L_HuffmanTreeList_pop_minor_element(menor, lista1, lista2); //menor = L_HuffmanTreeList::pop_minor_element(lista1, lista2);
		// Agregar el children izquierdo
		ptr->hijoDer=menor;
		ptr->hijoDer->padre=ptr;

		ptr->prof = L_MAX(ptr->hijoIzq->prof, ptr->hijoDer->prof) + 1;

		// Otros campos por si acaso...
		ptr->padre = NULL;
		ptr->sig = NULL;
		ptr->value = 255;

		// sum de los pesos
		ptr->peso=ptr->hijoIzq->peso+ptr->hijoDer->peso;

		lista2.push_back_node(ptr);
	}
	// Robarse el tree
	root=lista2.root;
	lista2.root = NULL;
	// En teoria desde cualquier hoja se deberia poder subir hasta la root... :S
}

bool L_HuffmanTree::search_for_incoherences_in_tree(const L_HuffmanTreeNode *inic) const
{
	if (inic == NULL)
	{
		if (root == NULL)
			return false;
		inic = root;
	}
	if (inic->sig != NULL)
		return false;
	if (inic->hijoIzq != NULL || inic->hijoDer != NULL)
		if (inic->value != 255)
			return false;
	if (inic->hijoIzq != NULL)
	{
		if (inic->hijoIzq->padre != inic)
			return false;
		if (search_for_incoherences_in_tree(inic->hijoIzq) == false)
			return false;
	}
	if (inic->hijoDer != NULL)
	{
		if (inic->hijoDer->padre != inic)
			return false;
		if (search_for_incoherences_in_tree(inic->hijoIzq) == false)
			return false;
	}
	return true;
}


void L_HuffmanTree::findPathInTree(unsigned char c, long &ruta, unsigned char &nBits)
{
	L_HuffmanTreeNode *ptr;
	int bit;
	long lo;
	if (root==NULL)
		return;
	ptr=children[c];  // Se recorre de la hoja a la root -> al reves del formato en que se guarda
	if (ptr->value != c)
		throw L_ArgException();

	nBits=0;
	lo = 0;
	ptr=children[c];
	// Las paths indican como llegar de la root a la punta ... esta invertido respecto a lo deseado
	while(ptr!=root)
	{
		if (ptr->padre->hijoDer==ptr) 
			lo|=(1L<<nBits);        
		nBits++;
		ptr=ptr->padre;
	}
	if (nBits > 32)
		printf("Profundidad excesiva en el tree\n");

	// Se recorrio al reves -> dar vuelta el resultado para que vaya de la root al children
	ruta = 0;
	for (bit=0; bit<nBits; bit++)
		ruta|=((lo>>bit)&1L)<<(nBits-1 - bit);
	// El resultado indica como ir de la root a la hoja
}

int L_HuffmanPath::canonicalOrdering(const void *huffRuta1, const void *huffRuta2)
{
	const L_HuffmanPath *r1 = (const L_HuffmanPath *)huffRuta1;
	const L_HuffmanPath *r2 = (const L_HuffmanPath *)huffRuta2;
	if (r1->numBits > r2->numBits)
		return 1;
	else if (r1->numBits < r2->numBits)
		return -1;
	else if (r1->value > r2->value)
		return 1;
	else if (r1->value < r2->value)
		return -1;
	return 0;
}

int L_HuffmanPath::lexicalOrdering(const void *huffRuta1, const void *huffRuta2)
{
	const L_HuffmanPath *r1 = (const L_HuffmanPath *)huffRuta1;
	const L_HuffmanPath *r2 = (const L_HuffmanPath *)huffRuta2;
	if (r1->value > r2->value)
		return 1;
	else if (r1->value < r2->value)
		return -1;
	return 0;
}

void L_HuffmanEncoder::buildPathsFromTree()
{
	int i;
	hasCanonicalOrdering = false;
	throw_L_ArgException_if(tree.root == NULL, "L_HuffmanEncoder::buildPathsFromTree");
	for (i=0; i<256; i++)
	{
		codeMapping[i].value = (unsigned char)i;
		tree.findPathInTree((unsigned char)i, codeMapping[i].trainOfBits, codeMapping[i].numBits);
	}
}


bool L_HuffmanEncoder::buildTreeFromPaths()
{
	int i;
	int j;
	L_HuffmanTreeNode **p;
	L_HuffmanTreeNode *padre;
	tree.destroyTree();
	for (i=0; i<256; i++)
	{
		p=&tree.root; // Se parte con un **p desde la root
		padre=NULL;
		for (j=0; j<codeMapping[i].numBits; j++) // paths[] va de la root a las hojas igual que aca (1 = derecha)
		{
			// Si faltan nodos intermedios, crearlos
			if (*p==NULL)
			{
				*p=new_L_HuffmanTreeNode2();
				(*p)->padre=padre;
				(*p)->value = 255;
			}
			padre=*p;
			// Elegir lado hacia el cual ir
			if (    ( (codeMapping[i].trainOfBits>>j)&1 )  == 0    )
				p=&(*p)->hijoIzq;
			else
				p=&(*p)->hijoDer;
		}
		if (*p!=NULL)
		{
			printf("Error en buildTreeFromPaths() : cadena duplicada\n");
			return false;
		}
		*p=new_L_HuffmanTreeNode2();
		(*p)->padre=padre;
		(*p)->value=i;
		(*p)->peso=0; // No hay como saber esto, pero no importa mucho
		tree.children[i]=*p; // Puntero-referencia, no maneja la memoria
	}
	return true;
	// tree.root queda apuntando-guardando al tree
}


bool L_HuffmanEncoder::check_for_incoherences_beetween_tree_and_paths()
{
	int i, bits, lado;
	L_HuffmanTreeNode *ptr;
	for (i=0; i<256; i++)
	{
		bits = 0;
		ptr=tree.root; // Recorrer desde la root a las hojas guiandose por (codeMapping[i].trainOfBits>>bits)&1
		while(true)
		{
			if (ptr==NULL) // Esto es porque se siguio un camino que lleva fuera del tree
				return false;
			if (ptr->hijoIzq==NULL && ptr->hijoDer==NULL) // Hoja
				break;
			lado = (codeMapping[i].trainOfBits>>bits)&1;
			if (lado == 0)  // // Derecha -> addFrom un (1L<<nBits)
				ptr=ptr->hijoIzq;
			else
				ptr=ptr->hijoDer;
			bits++;
			if (bits > codeMapping[i].numBits)
				return false;
		}
		if (ptr->value != (unsigned char)i)
			return false;
		if (bits != codeMapping[i].numBits)
			return false;
	}
	return true;
	// tree.root queda apuntando-guardando al tree
}

void L_HuffmanEncoder::buildCanonicalOrdering()
{
	unsigned int i, bit, nBits;
	long lo;
	L_String str[256];
	tree.destroyTree(); // El tree queda obsoleto con esto
	qsort(codeMapping, 256, sizeof(L_HuffmanPath), &L_HuffmanPath::canonicalOrdering);
	// Ahora asignar codigos nuevos a partir de los largos de paths
	// EL primero es un tren de ceros
	codeMapping[0].trainOfBits = 0;
	codeMapping[0].numBits = codeMapping[0].numBits;
	// Copiado de la wikipedia, ojala que funke
	for (i=1; i<256; i++)
		codeMapping[i].trainOfBits = (codeMapping[i-1].trainOfBits+1) << (codeMapping[i].numBits - codeMapping[i-1].numBits);
	// Invertir la secuencia de las paths para que indiquen como ir del children a la root
	for (i=0; i<256; i++)
	{
		lo = codeMapping[i].trainOfBits;
		nBits = codeMapping[i].numBits;
		codeMapping[i].trainOfBits = 0;
		for (bit=0; bit<nBits; bit++)
			codeMapping[i].trainOfBits |= ((lo>>bit)&1L)<<(nBits-1 - bit);
	}
	qsort(codeMapping, 256, sizeof(L_HuffmanPath), &L_HuffmanPath::lexicalOrdering);
	hasCanonicalOrdering = true;
}

void L_HuffmanEncoder::build(const char *fuente, long len)
{
	tree.clear();
	tree.addVotesToFrequencies(fuente, len);
	tree.build_tree_using_stored_frequeencies();
	buildPathsFromTree();
}

void L_HuffmanEncoder::encode(const char *fuente, long len)
{
	long bitTotal=0; // Podria ser bastante grande
	long i;

	// Pasada para contar la cantidad de bytes, es penca presuponerlo

	bitTotal = 0;
	for (i=0; i<len; i++)
		bitTotal+=codeMapping[(unsigned char)fuente[i]].numBits;
	res.resize(bitTotal/8 + (bitTotal%8 > 0));

	res.reserve(res.size() + 8);

	encode(fuente, len, &(res[0]));
}


int L_HuffmanEncoder::encode0(const char *fuente, long len, char *dest)
{
	long bitTotal=0; // Podria ser bastante grande
	int anchoBits;
	long i, bit;
	long value;

	bitTotal = 0;  // Se van acumulando
	for (i=0; i<len; i++)
	{
		// Guardar fuente[i] en (bitTotal, byteTotal)
		value=codeMapping[(unsigned char)fuente[i]].trainOfBits;
		anchoBits=codeMapping[(unsigned char)fuente[i]].numBits;
		for (bit=0; bit<anchoBits; bit++)
			L_CH_setBit(dest,bitTotal+bit,(value>>bit)&1);
		bitTotal += anchoBits;
	}
	int lenFinal=bitTotal/8+(bitTotal%8 > 0);
	return lenFinal;
}



int L_HuffmanEncoder::encode(const char *fuente, long len, char *dest) // PELIGROSA
{
	long bitTotal=0; // Podria ser bastante grande
	long i;
	int maxCodigo = 0;
	//int res; // Para evitar confusiones
	//res = 0;

	for (i=0; i<256; i++)
	{
		if (codeMapping[i].numBits > maxCodigo)
			maxCodigo = codeMapping[i].numBits;
	}

	bitTotal = 0;

	// (buf)[(bitPos)/8] |= 1 << ((bitPos)%8)
	*(long *)&dest[0] = 0;
	for (i=0; i<len; i++)
	{
		// Visual studio cool.. x>>32 = x !!!. La segunda linea sirve ademas para limpiar el arreglo
		*(long *)&dest[bitTotal/8] |= ( codeMapping[(unsigned char)fuente[i]].trainOfBits << (bitTotal%8) );
		*(long *)&dest[bitTotal/8+sizeof(long)] = (bitTotal%8 > 0) * ( codeMapping[(unsigned char)fuente[i]].trainOfBits >> (8*sizeof(long) - (bitTotal%8)) );
		bitTotal += codeMapping[(unsigned char)fuente[i]].numBits;
	}

	int lenCompr = bitTotal/8+(bitTotal%8 > 0);
	return lenCompr;
}

bool L_HuffmanEncoder::decode(const char *fuente, long len, long lenFinal)
{
	res.resize(lenFinal); // Por si acaso
	long bitTotal=0;
	long n=0;
	int lado, bit;
	L_HuffmanTreeNode *ptr;
	throw_L_ArgException_if(tree.root == NULL, "L_HuffmanEncoder::decode() : tree nulo");
	while(true)
	{
		ptr=tree.root; // Recorrer desde la root a las hojas guiandose por fuente[byteTotal]&(1L<<bitTotal)
		bit = 0;
		while(true)
		{
			throw_L_ArgException_if(ptr==NULL, "L_HuffmanEncoder::decode() : tree incoherente"); // Camino fuera del tree
			if (ptr->hijoIzq==NULL && ptr->hijoDer==NULL) // Hoja
				break;
			lado = L_CH_getBit(fuente,bitTotal+bit);
			if (lado == 0)  // // Derecha -> addFrom un (1L<<nBits)
				ptr=ptr->hijoIzq;
			else
				ptr=ptr->hijoDer;
			if (bitTotal/8+(bitTotal%8 > 0) > len || n==lenFinal) // output de emergencia
				return false;
			bit++;
		}
		res[n++]=ptr->value;
		bitTotal += bit;
		if (bit != codeMapping[ptr->value].numBits)
			return false;
		if (bitTotal/8+(bitTotal%8 > 0) > len || n==lenFinal) // output de emergencia
			break;
	}
	int lenCompr = bitTotal/8+(bitTotal%8 > 0);
	if (lenCompr != len || n!=lenFinal)
		return false;
	return true;
}

void L_HuffmanEncoder::writeHuffmanTxt(FILE *fp)
{
	int i;
	for (i=0; i<256; i++)
		fprintf(fp, "%d %ld", codeMapping[i].numBits, codeMapping[i].trainOfBits);
	printf("\n");
}

void L_HuffmanEncoder::readHuffmanTxt(FILE *fp)
{
	int i, nb;
	clear();
	for (i=0; i<256; i++)
	{
		if(fscanf(fp, "%d%ld", &nb, &codeMapping[i].trainOfBits) != 2)
			{printf("L_HuffmanEncoder::readHuffmanTxt() : error en formato de archivo\n"); return;}
		codeMapping[i].numBits = nb;
		codeMapping[i].value = i;
	}
	buildTreeFromPaths();
	printf("\n");
}

int L_HuffmanEncoder::lengthOfHuffmanBin()
{
	int ret=256;
	int i;
	for (i=0; i<256; i++)
		ret+=codeMapping[i].numBits/8 + ((codeMapping[i].numBits%8)>0);
	return ret;
}

void L_HuffmanEncoder::writeHuffmanBin(FILE *fp)
{
	std::vector<char> buf;
	writeHuffman(buf); // Guarda las paths de Huffman
	fwrite(&(buf[0]), sizeof(char), buf.size(), fp);
}

int L_HuffmanEncoder::readHuffmanBin(FILE *fp) // read las paths de Huffman
{
	int i, j, bit, u=0, anchoBits, anchoBytes;
	int c;
	unsigned char buf[4];
	for (i=0; i<256; i++)
	{
		c = fgetc(fp);
		u++;
		if (c==EOF)
			return -1;
		codeMapping[i].numBits=c;
		codeMapping[i].trainOfBits=0;
		codeMapping[i].value = i;
		anchoBits = codeMapping[i].numBits;
		anchoBytes = anchoBits/8 + ((anchoBits%8)>0);
		for (j=0; j<anchoBytes; j++)
		{
			buf[j] = fgetc(fp);
			u++;
		}
		for (bit=0; bit<anchoBits; bit++)
			codeMapping[i].trainOfBits|=L_CH_getBit(buf, bit)<<bit;
	}
	buildTreeFromPaths();
	return u;
}


void L_HuffmanEncoder::writeHuffman(std::vector<char> &buffer)
{
	int i, anchoBits, bit=0, u=0;

	buffer.resize(lengthOfHuffmanBin());

	for (i=0; i<256; i++)
	{
		buffer[u++] = codeMapping[i].numBits;
		anchoBits = codeMapping[i].numBits;
		for (bit=0; bit<anchoBits; bit++)
			L_CH_setBit(&(buffer[0]), u*8+bit, (codeMapping[i].trainOfBits>>bit)&1);
		u += anchoBits/8 + (anchoBits%8 > 0);
	}
}

int L_HuffmanEncoder::readHuffman(const std::vector<char> &buffer)
{
	int i, anchoBits, bit=0, u=0;

	for (i=0; i<256; i++)
	{
		codeMapping[i].numBits=buffer[u++];
		codeMapping[i].trainOfBits=0;
		codeMapping[i].value = i;
		anchoBits = codeMapping[i].numBits;
		for (bit=0; bit<anchoBits; bit++)
			codeMapping[i].trainOfBits|=L_CH_getBit(&(buffer[0]), u*8+bit)<<bit;
		u += anchoBits/8 + (anchoBits%8 > 0);
	}
	buildTreeFromPaths();
	return u;
}

void L_HuffmanEncoder::writeHuffmanCanonical(std::vector<char> &buffer)
{
	// Guardar solo los tamanos de las paths => "HFL32" + 160 bytes
	int bit=0, i;
	buffer.resize(160 + 5);

	buffer[0] = 'H';
	buffer[1] = 'F';
	buffer[2] = 'L';
	buffer[3] = '3';
	buffer[4] = '2';
	bit += 40;

	for (i=0; i<256; i++)
	{
		L_CH_setBit(&(buffer[0]), bit, (codeMapping[i].numBits>>0)&1);
		bit++;
		L_CH_setBit(&(buffer[0]), bit, (codeMapping[i].numBits>>1)&1);
		bit++;
		L_CH_setBit(&(buffer[0]), bit, (codeMapping[i].numBits>>2)&1);
		bit++;
		L_CH_setBit(&(buffer[0]), bit, (codeMapping[i].numBits>>3)&1);
		bit++;
		L_CH_setBit(&(buffer[0]), bit, (codeMapping[i].numBits>>4)&1);
		bit++;
	}
}

int L_HuffmanEncoder::readHuffmanCanonical(const std::vector<char> &buffer)
{
	// Leer solo los tamanos de las paths
	// Guardar solo los tamanos de las paths => "HFL32" + 160 bytes
	int bit=0, i;
	throw_L_ArgException_if(buffer.size() < 165, "L_HuffmanEncoder::readHuffmanCanonical()");

	if (buffer[0] != 'H' || buffer[1] != 'F' || buffer[2] != 'L' || buffer[3] != '3' || buffer[4] != '2')
		return 0; // No es este tipo de archivo

	bit += 40;
	for (i=0; i<256; i++)
	{
		codeMapping[i].numBits = 0;
		codeMapping[i].numBits |= L_CH_getBit(&(buffer[0]), bit) << 0;
		bit++;
		codeMapping[i].numBits |= L_CH_getBit(&(buffer[0]), bit) << 1;
		bit++;
		codeMapping[i].numBits |= L_CH_getBit(&(buffer[0]), bit) << 2;
		bit++;
		codeMapping[i].numBits |= L_CH_getBit(&(buffer[0]), bit) << 3;
		bit++;
		codeMapping[i].numBits |= L_CH_getBit(&(buffer[0]), bit) << 4;
		bit++;
		codeMapping[i].value = i;
	}
	return (40 + 256*5) / 8; // 165
}

void L_HuffmanEncoder::encodeAll(const std::vector<char> &input, std::vector<char> &output)
{
	union L_LongChars nA, nD, nC;
	std::vector<char> datosArbol;
	L_HuffmanEncoder huffman;
	int i, u;
	// Construir tree y guardarlo en datosArbol
	huffman.build(&(input[0]), (int)input.size());
	huffman.buildCanonicalOrdering();
	huffman.writeHuffmanCanonical(datosArbol);
	// Construir datos codificados y guardarlos en datosCodif
	nA.lo = (long)datosArbol.size(); // datosArbol[]
	nD.lo = (long)input.size(); // input[]
	nC.lo = (long)input.size() * 3; // quizas cuanto va a ser...

	output.resize((int)(4+4+4 + nA.lo + 5 + nC.lo));

	u = 0; // write la output correlativamente
	// 12 bytes para los tamanos
	for (i=0; i<4; i++,u++) // u=0
		output[u] = nA.ch[i];
	for (i=0; i<4; i++,u++) // u=4
		output[u] = nD.ch[i];
	for (i=0; i<4; i++,u++) // u=8
		output[u] = nC.ch[i];
	// nA.lo bytes para la tabla Huffman (= 165 bytes)
	for (i=0; i<(int)datosArbol.size(); i++,u++) // u=12
		output[u] = datosArbol[i];
	// 5 bytes para encabezado
	output[u++] = 'H';  // u = 12 + datosArbol.size()
	output[u++] = 'F';
	output[u++] = 'D';
	output[u++] = 'A';
	output[u++] = 'T';
	// u = 12 + datosArbol.size() + 5
	// nC.lo bytes para los datos


	nC.lo = huffman.encode(&(input[0]), (long)input.size(), &(output[0])+u); //  res[]... PELIGROSO pero eficiente
	//nC.lo = huffman.encode0(input.data(), input.size(), output.data()+u); //  res[]... PELIGROSO pero eficiente

	// Reescribir el verdadero value de nC (recien ahora se conoce)
	u = 8;
	for (i=0; i<4; i++,u++)
		output[u] = nC.ch[i];
	output.resize((int)(12 + nA.lo + 5 + nC.lo));
	// Fin de la grabacion de la output
}

bool L_HuffmanEncoder::decodeAll(const std::vector<char> &input, std::vector<char> &output)
{
	union L_LongChars nA, nD, nC;
	L_HuffmanEncoder huffman;
	std::vector<char> datosArbol, datosCompr;
	int i, u, n;

	u = 0; // read la input correlativamente
	// 12 bytes para los tamanos
	for (i=0; i<4; i++,u++)
		nA.ch[i] = input[u];
	for (i=0; i<4; i++,u++)
		nD.ch[i] = input[u];
	for (i=0; i<4; i++,u++)
		nC.ch[i] = input[u];

	datosArbol.resize(nA.lo);
	datosCompr.resize(nC.lo);

	// nA.lo bytes para la tabla Huffman (= 165 bytes)
	for (i=0; i<(int)datosArbol.size(); i++,u++)
		datosArbol[i] = input[u];
	// 5 bytes para encabezado
	if (input[u++] != 'H' || input[u++] != 'F' || input[u++] != 'D' || input[u++] != 'A' || input[u++] != 'T')
		return false;
	// nC.lo bytes para los datos
	for (i=0; i<(int)datosCompr.size(); i++,u++)
		datosCompr[i] = input[u];
	// Fin de la lectura de la input

	n = huffman.readHuffmanCanonical(datosArbol); // n = bytes del tree
	huffman.buildCanonicalOrdering();
	huffman.buildTreeFromPaths();
	if (n != nA.lo) // Esto no provoca problemas de memoria
		return false;
	if (huffman.decode(&(datosCompr[0]), (long)datosCompr.size(), nD.lo) == false)
		return false;
	if (huffman.res.size() != nD.lo)
		return false;
	output.resize(nD.lo);

	for (i=0; i<(int)output.size(); i++)
		output[i] = huffman.res[i];
	return true;
}

double L_HuffmanEncoder::entropy_in_bits(const std::vector<char> &input)
{
	double entrop = 0;
	double distr[256];
	int i;
	// Distribucion de probabilidad de los caracteres
	for (i=0; i<256; i++)
		distr[i] = 0;
	for (i=0; i<(size_type)input.size(); i++)
		distr[(unsigned char)input[i]]++;
	for (i=0; i<256; i++)
		distr[i] /= input.size();
	// Entropia de la distribucion (en bits)
	for (i=0; i<256; i++)
		if (distr[i] > 0)
			entrop += log(1/distr[i])/L_LOG_2*distr[i];
	return entrop;
}

void L_HuffmanEncoder::encodeRLE0(const std::vector<char> &input, std::vector<char> &salidaRLE)
{
	int i=1, u=1, nceros=0;
	salidaRLE.resize(2*input.size()); // Mejor ser previsor

	salidaRLE[0] = input[0];

	while (i<(int)input.size())
	{
		// El caso mas comun, un ciclo simple y rapido
		while (i < (int)input.size() && nceros == 0 && (input[i] != 0 || input[i-1] != 0))
			salidaRLE[u++] = input[i++];
		// Si llegamos aca es porque hay que empezar a acumular ceros
		if (i==input.size()) // Tres casos: termina con xx, con x0 o con 00 (en los 2 ultimos casos hay que agregar 2 ceros mas)
		{
			if (input[i] == 0) // Termina con x0: agregar solo un zero
			{
				salidaRLE[u++] = 0;
				if (input[i-1] == 0) // Termina con 00: agregar dos ceros mas
				{
					salidaRLE[u++] = 0;
					salidaRLE[u++] = 0;
				}
			}
			break; // termina con xx: no hacer nada
		}
		if (nceros != -1)
		{
			salidaRLE[u++] = 0; // Agregar el segundo zero pq no se guardo, pero evitarlo si hay un tren gigante de ceros
			i++; // Avanzar uno...
		}
		nceros = 0;
		// El siguiente caso mas comun, dentro de un RLE
		while (i < (int)input.size() && nceros < 255 && input[i] == 0)
		{
			i++;
			nceros++;
		}
		// output del RLE
		salidaRLE[u++] = nceros;
		if (nceros == 255)
			nceros = -1; // Se saturo, seguir guardando ceros
		else
			nceros = 0; // Se acabo el tren de ceros
		// Si hay 255 ceros, va a quedar guardado (0 0 255 ...), lo que es un patron reconocible
		// que indica que siguen ceros consecutivos, ej: si hay 255+20 ceros deberia quedar (0 0 255 20)
	}
	salidaRLE.resize(u);
}

void L_HuffmanEncoder::decodeRLE0(const std::vector<char> &entradaRLE, std::vector<char> &output)
{
	int i=1, u=1, k, nceros=0;
	if (output.size() < (int)entradaRLE.size())
	{
		output.resize(entradaRLE.size()*2);
	}
	output[0] = entradaRLE[0];
	while(i < (int)entradaRLE.size())
	{
		// El caso mas comun, un ciclo simple y rapido
		while (i < (int)entradaRLE.size() && u < (int)output.size() && nceros == 0 && (entradaRLE[i] != 0 || entradaRLE[i-1] != 0))
			output[u++] = entradaRLE[i++];
		if (i == entradaRLE.size())
			break;
		if (u == output.size()) // Se acabo la memoria
		{
			output.resize(2*output.size());
			continue;
		}
		output[u++] = 0; // Falto este
		i++;
		// Si llegamos aca es porque hay ceros acumulador
		nceros = 255; // Solo para comenzar el ciclo
		while (i < (int)entradaRLE.size() && nceros == 255)
		{
			nceros = (unsigned char)entradaRLE[i];
			// De acá hacia abajo no se modifica i...
			for (k=0; k < nceros && u < (int)output.size(); k++)
				output[u++] = 0;
			while (u == output.size()) // Hay que agrandar la output
			{
				output.resize(2*output.size());
				for (; k < nceros && u < (int)output.size(); k++)
					output[u++] = 0;
			}
			i++;
		}
		nceros = 0;
	}
	output.resize(u);
}


void L_HuffmanEncoder::encodeDelta(const std::vector<char> &input, std::vector<char> &salidaD)
{
	int i;
	salidaD.resize(input.size());
	salidaD[0] = input[0];
	for (i=1; i<(int)salidaD.size(); i++)
		salidaD[i] = input[i] - input[i-1];
}

void L_HuffmanEncoder::decodeDelta(const std::vector<char> &entradaD, std::vector<char> &output)
{
	int i;
	output.resize(entradaD.size());
	output[0] = entradaD[0];
	for (i=1; i<(int)output.size(); i++)
		output[i] = output[i-1] + entradaD[i];
}


double L_LevenbergMarquardt::minimize(const void *object, std::vector<double> &vector, double (*fnToMinimize)(const void *, double *))
{
	throw_L_ArgException_if(!readyToWork(true), "L_LevenbergMarquardt::minimize");
	double lambda=0.1;
	std::vector<double> grad(lenVector);
	L_Matrix hess;
	std::vector<double> vectorAnt(lenVector);
	double error = 666.0, errorTemp, errorAnt=1e30;
	bool vectorPegado=false, vectorPegado2=false;
	if (useHessian)
		hess.reallocate(lenVector,lenVector);
	long nIter;
	int i, j;
	for (nIter=0; nIter<nIterationsMax && errorAnt > errorToStop; nIter++)
	{
		if (nIter==0)
		{
			errorAnt=(*fnToMinimize)(object,&(vector[0]));
			error=errorAnt;
			for (i=0; i<lenVector; i++)
				vectorAnt[i]=vector[i];
			continue;
		}

		for (i=0; i<lenVector; i++)
		{
			vector[i]+=epsilon;
			errorTemp=(*fnToMinimize)(object,&(vector[0]));
			grad[i]=( errorTemp-errorAnt )/epsilon;
			vector[i]-=epsilon;
			if (useHessian)
			{
				for (j=0; j<lenVector; j++)
				{
					vector[j]+=epsilon;
					error=(*fnToMinimize)(object,&(vector[0]));
					vector[i]+=epsilon;
					errorTemp=(*fnToMinimize)(object,&(vector[0]));
					hess(i,j)= ( (errorTemp-error)/epsilon - grad[i] )/epsilon;
					vector[i]-=epsilon;
					vector[j]-=epsilon;
				}
			}
		}
		if (useHessian)
		{
			for (i=0; i<lenVector; i++)
				hess(i,i)+=lambda;
			if (hess.invertMe()==false)
				printf("Aca\n");
			for (i=0; i<lenVector; i++)
				for (j=0; j<lenVector; j++)
					vector[i]-=hess(i,j)*grad[j];
		}
		else
		{
			for (i=0; i<lenVector; i++)
				vector[i]-=grad[i]/lambda;
		}
		if (useHessian)
		{
			vectorPegado=true;
			for (i=0; i<lenVector; i++)
				if (vector[i]!=vectorAnt[i])
					vectorPegado=false;
			if (vectorPegado)
				useHessian=false;
		}
		else
		{
			vectorPegado2=true;
			for (i=0; i<lenVector; i++)
				if (vector[i]!=vectorAnt[i])
					vectorPegado2=false;
			if (vectorPegado2)
				break; // Ya no se puede optimizar mas
		}
		error=(*fnToMinimize)(object,&(vector[0]));
		if (error > 1e30) // La minimizacion no esta funcionando
			return 1e30;
		if (error < errorAnt)
		{
			if (lambda>1e-30)
				lambda/=factorGrowth;
			for (i=0; i<lenVector; i++)
				vectorAnt[i]=vector[i];
			errorAnt=error;
		}
		else
		{
			if (lambda<1e30)
				lambda*=factorGrowth;
			else
				printf("Wea\n");
			for (i=0; i<lenVector; i++)
				vector[i]=vectorAnt[i];
		}
	}
	if (vectorPegado==true)
		useHessian=true; // Devolverlo a su estado original
	return error;
}





double L_LevenbergMarquardt::minimize_vect(const void *object, L_Matrix &x, void (*fnToMinimize)(const void *, const L_Matrix &, L_Matrix &), void (*jacob)(const void *object, const L_Matrix &x, L_Matrix &Jdedx), void (*preIter)(const void *object, L_Matrix &x))
{
	throw_L_ArgException_if(!readyToWork(true),"L_LevenbergMarquardt::minimize_vect");
	int n=lenVector;
	// variables del object:
	//long nIterationsMax;
	//double errorToStop;
	//double factorGrowth; // lambda crece/disminuye en este factor para adaptarse a la funcion de error
	//double epsilon;
	//double lambda; // Factor de escala para pasar el gradiente a desplazamiento
	//bool usarInversion;
	//int lenVector;

	L_Matrix J, JT, h2, JTJ, JTJrep, JTy;
	L_Matrix dx, xSeg, y, ySeg, yMod, B, z;
	double err, errSeg; // err = norma de y
	int iter, i, j;
	double epsilon = this->epsilon; // Para proteger la variable interna pq se cambia
	int signo;
	bool calcJacobiano = true;

	// Preparar el object para la primera evaluacion
	if (preIter != NULL)
		(*preIter)(object, x);
	// Evaluar el error inicial
	(*fnToMinimize)(object, x, y);
	err = y.normL2();

	xSeg = x;
	ySeg = y;
	errSeg = err;

	J.reallocate(y.li, x.li);
	z.reallocate(y.li, 1);
	for (iter = 0; iter < nIterationsMax; iter++)
	{
		if (calcJacobiano) // Se necesita calcular un jacobiano en esta iteracion
		{
			if (jacob == NULL)
			{
				for (j=0; j<x.li; j++)
				{
					signo = 2*(((iter+1)*j)%2)-1; // Para no tener sesgo con el signo de epsilon
					x(j,0) += (epsilon*signo);
					(*fnToMinimize)(object, x, yMod);
					for (i=0; i<y.li; i++)
						J(i,j) = (yMod(i,0) - y(i,0)) / (epsilon*signo);
					x(j,0) -= (epsilon*signo);
				}
			}
			else
				(*jacob)(object, x, J);
			JT.transpOf(J);
			JTJ.OP_mult(JT,J);
		}
		JTJrep.reallocate(n, n);
		for (i=0; i<n; i++)
			for (j=0; j<n; j++)
				JTJrep(i,j) = JTJ(i,j) + (i==j)*lambda;//*JTJ(i,j)); // Solucion mixta
		// "B = -JT*y;"
		JTy.OP_mult(JT,y);
		B.OP_mult(JTy,-1);
		if (B.solveLinearSystem_x_Ab_result_in_b(JTJrep) == false) // Generar desplazamiento en Y
			printf("minimize_vect(): JTJ+lambda no invertible\n");
		dx = B;
		x+=dx;

		// Preparar el object para el cambio abrupto de punto
		if (preIter != NULL)
			(*preIter)(object, x);
		(*fnToMinimize)(object, x, y); // Evaluar el error en el nuevo punto
		err = y.normL2();

		calcJacobiano = true;

		if (err < errSeg)
		{
			xSeg = x; // Aca si sirvio, el error disminuyo
			ySeg = y;
			errSeg = err;
			lambda /= factorGrowth;
		}
		else if (err == errSeg) // El desplazamiento fue demasiado pequeno, el error no cambio
			epsilon*= factorGrowth; // Capaz que sea por un epsilon muy pequeno
		else
		{
			x = xSeg; // Devolverse de punto, el desplazamiento fue excesivo
			y = ySeg;
			err = errSeg;
			if (preIter != NULL)
				(*preIter)(object, x);
			lambda *= factorGrowth;
			if (lambda > 1e60)
				epsilon/= factorGrowth; // Capaz que sea por un epsilon muy grande
			calcJacobiano = false;
		}
		if (err == 0) // Que maravilloso
			break;
	}
	x = xSeg;
	err = errSeg;
	return err;
}






double L_Random_plus_LM_optimizer::minimize_vect(const void *object, L_Matrix &x, void (*fnToMinimize)(const void *, const L_Matrix &, L_Matrix &), void (*jacob)(const void *object, const L_Matrix &x, L_Matrix &Jdedx), void (*preIter)(const void *object, L_Matrix &x))
{
	int i, j;
	L_Matrix errVect;
	double errMin = 1.0e60;
	for (i=0; i<nIterationsMax; i++)
	{
		paths.resize(paths.size()+1);
		paths[i].puntos.resize(2);
		paths[i].errores.resize(2);

		for (j=0; j<x.li; j++)
			x(i,0) = L_RANDOM(vectorLimInf[j],  vectorLimSup[j]);

		(*fnToMinimize)(object, x, errVect);
		paths[i].puntos[0] = x;
		paths[i].errores[0] = errVect.normL2();
		if (paths[i].errores[0] < errMin)
			errMin = paths[i].errores[0];

		if (usarLM)
		{
			LM.minimize_vect(object, x, fnToMinimize, jacob, preIter);

			paths[i].puntos[1] = x;
			paths[i].errores[1] = errVect.normL2();
			if (paths[i].errores[0] < errMin)
				errMin = paths[i].errores[0];
		}
	}
	return errMin;
}



bool L_VerificaIntegridadCodigo();
bool L_VerificaIntegridadCodigoRand(double prob);

#ifdef VERIF_INT
extern char strIntegr_1[];
#endif

bool L_ImageGrayDouble::normalizeHistogram(int borde)
{
	int i, j;
	double max, min;
	double dif;
	bool ret = true;

	throw_L_ArgException_if(data()==NULL || 2*borde > lx || 2*borde > ly, "L_ImageGrayDouble::normalizeHistogram() : Llamada con imagen nula\n");

	max=min=pix(borde,borde);

	for (j=borde; j<ly-borde; j++)
	{
		for (i=borde; i<lx-borde; i++)
		{
			if (pix(i,j)<min)
				min=pix(i,j);
			else if (pix(i,j)>max)
				max=pix(i,j);
		}
	}
	dif=max-min;
	if (dif==0)
	{
		dif = 1;
		ret = false;
	}
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j)-=min;
			pix(i,j)/=dif;
		}
	}
	return ret;
}

L_ImageGrayDouble& L_ImageGrayDouble::genChannelY(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB)
{
	int i, j;
	throw_L_ArgException_if(imR.lx!=imG.lx || imG.lx!=imB.lx || imR.ly!=imG.ly || imG.ly!=imB.ly, "L_ImageGrayDouble::genChannelY");
	reallocate(imR.lx, imR.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelY(imR.pix(i,j),imG.pix(i,j),imB.pix(i,j));
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::genChannelU(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB)
{
	int i, j;
	throw_L_ArgException_if(imR.lx!=imG.lx || imG.lx!=imB.lx || imR.ly!=imG.ly || imG.ly!=imB.ly, "L_ImageGrayDouble::genChannelU");
	reallocate(imR.lx, imR.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelU(imR.pix(i,j),imG.pix(i,j),imB.pix(i,j));
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::genChannelV(L_ImageGrayDouble &imR, L_ImageGrayDouble &imG, L_ImageGrayDouble &imB)
{
	int i, j;
	throw_L_ArgException_if(imR.lx!=imG.lx || imG.lx!=imB.lx || imR.ly!=imG.ly || imG.ly!=imB.ly, "L_ImageGrayDouble::genChannelV");
	reallocate(imR.lx, imR.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelV(imR.pix(i,j),imG.pix(i,j),imB.pix(i,j));
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::genChannelR(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV)
{
	int i, j;
	throw_L_ArgException_if(imY.lx!=imU.lx || imU.lx!=imV.lx || imY.ly!=imU.ly || imU.ly!=imV.ly , "L_ImageGrayDouble::genChannelR");
	reallocate(imY.lx, imY.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelR(imY.pix(i,j),imU.pix(i,j),imV.pix(i,j));
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::genChannelG(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV)
{
	int i, j;
	throw_L_ArgException_if(imY.lx!=imU.lx || imU.lx!=imV.lx || imY.ly!=imU.ly || imU.ly!=imV.ly, "L_ImageGrayDouble::genChannelG");
	reallocate(imY.lx, imY.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelG(imY.pix(i,j),imU.pix(i,j),imV.pix(i,j));
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::genChannelB(L_ImageGrayDouble &imY, L_ImageGrayDouble &imU, L_ImageGrayDouble &imV)
{
	int i, j;
	throw_L_ArgException_if(imY.lx!=imU.lx || imU.lx!=imV.lx || imY.ly!=imU.ly || imU.ly!=imV.ly, "L_ImageGrayDouble::genChannelB");
	reallocate(imY.lx, imY.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(double)channelB(imY.pix(i,j),imU.pix(i,j),imV.pix(i,j));
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::OP_mult_elementwise_sub(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2, int x0, int y0, int rx, int ry)
{
	throw_L_ArgException_if(im1.lx!=im2.lx || im1.ly!=im2.ly, "L_ImageGrayDouble::OP_mult_elementwise");
	int i, j;
	int iIni, jIni;
	int iFin, jFin;

	reallocate(im1.lx, im2.ly);

	iIni = x0 - rx;
	iFin = x0 + rx + 1;
	jIni = y0 - ry;
	jFin = y0 + ry + 1;

	if (iIni < 0)
		iIni = 0;
	if (iFin > lx)
		iFin = lx;
	if (jIni < 0)
		jIni = 0;
	if (jFin > ly)
		jFin = ly;

	for (j=jIni; j<jFin; j++)
		for (i=iIni; i<iFin; i++)
			pix(i,j) = im1.pix(i,j) * im2.pix(i,j); // Mejor asi que con punteros, sino tiene que revisar *NULL muchas veces
	return *this;
}

//!< a.OPdiv(b,c) deja en "a" la division "b/c"
L_ImageGrayDouble& L_ImageGrayDouble::OP_div_elementwise(const L_ImageGrayDouble &im1, const L_ImageGrayDouble &im2)
{
	throw_L_ArgException_if(im1.lx!=im2.lx || im1.ly!=im2.ly, "L_ImageGrayDouble::OP_div_elementwise");
	int i, j;
	reallocate(im1.lx, im2.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j) = im1.pix(i,j) / im2.pix(i,j); // Mejor asi que con punteros, sino tiene que revisar *NULL muchas veces
	return *this;
}


void L_ImageGrayDouble::fourier(int dir, L_ImageGrayDouble &imag)
{
	int log2x=0, log2y=0;
	int lxx, lyy;
	int i, j;
	// Revisar que el tamano es potencia de 2
	lxx = lx;
	lyy = ly;
	while (lxx > 1)
	{
		lxx /= 2;
		log2x++;
	}
	while (lyy > 1)
	{
		lyy /= 2;
		log2y++;
	}
	lxx = 1;
	lyy = 1;
	for (i=0; i<log2x; i++)
		lxx *= 2;
	for (j=0; j<log2y; j++)
		lyy *= 2;

	throw_L_ArgException_if(lxx != lx || lyy != ly, "L_ImageGrayDouble::fourier");

	// lx, ly invertidos
	double **pp = L_new2d_ext<double>(li, lj, data(), ljStep);
	double **pi = L_new2d_ext<double>(imag.li, imag.lj, imag.data(), imag.ljStep);
	for (i=0; i<ly; i++)
		L_FFT1D(dir, pp[i], pi[i], log2x);
	for (j=0; j<lx; j++)
		L_FFT1D_2(dir, pp, pi, log2y, j);
	L_delete2d_ext<double>(pp);
	L_delete2d_ext<double>(pi);
}


void L_ImageGrayDouble::houghDetPuntosBordes(const L_ImageGrayDouble &imorig, L_ImageGrayUchar &binariz, std::vector<double> &px, std::vector<double> &py, std::vector<double> &ang, double umbral, bool supresionNoMaximo)
{
	L_ImageGrayDouble T, Gx, Gy, Gr2, imang;
	L_ImageGrayUchar imBN;
	int i, j;
	double umbral2 = umbral*umbral;
	T.reallocate(imorig.lx, imorig.ly);
	Gx.reallocate(imorig.lx, imorig.ly);
	Gy.reallocate(imorig.lx, imorig.ly);
	Gr2.reallocate(imorig.lx, imorig.ly);
	imang.reallocate(imorig.lx, imorig.ly);
	binariz.reallocate(imorig.lx, imorig.ly);
	T.setZero();
	px.resize(0);
	py.resize(0);
	ang.resize(0);
	// Limpiar bordes de los gradientes
	for (i=0; i<Gx.lx; i++) {Gx.pix(i,0) = Gx.pix(i,Gx.ly-1) = Gy.pix(i,0) = Gy.pix(i,Gy.ly-1) = 0.0;}
	for (j=0; i<Gy.ly; i++) {Gx.pix(0,j) = Gx.pix(Gx.lx-1,j) = Gy.pix(0,j) = Gy.pix(Gy.lx-1,j) = 0.0;}
	// Convoluciones para gradientes
	for (j=0; j<imorig.ly; j++)
		for (i=1; i<imorig.lx-1; i++)
			T.pix(i,j) = imorig.pix(i+1,j)-imorig.pix(i-1,j);
	for (j=1; j<imorig.ly-1; j++)
		for (i=1; i<imorig.lx-1; i++)
			Gx.pix(i,j) = (T.pix(i,j-1) + 2*T.pix(i,j) + T.pix(i,j+1)) * 1.0 / 4.0;
	for (j=1; j<imorig.ly-1; j++)
		for (i=0; i<imorig.lx; i++)
			T.pix(i,j) = imorig.pix(i,j+1)-imorig.pix(i,j-1);
	for (j=1; j<imorig.ly-1; j++)
		for (i=1; i<imorig.lx-1; i++)
			Gy.pix(i,j) = (T.pix(i-1,j) + 2*T.pix(i,j) + T.pix(i+1,j)) * 1.0 / 4.0;
	// Deteccion bordes
	for (j=0; j<imorig.ly; j++)
	{
		for (i=0; i<imorig.lx; i++)
		{
			Gr2.pix(i,j) = Gx.pix(i,j)*Gx.pix(i,j) + Gy.pix(i,j)*Gy.pix(i,j);
			binariz.pix(i,j) = (Gr2.pix(i,j) > umbral2) * 255;
			if (binariz.pix(i,j))
				imang.pix(i,j) = atan2(Gy.pix(i,j), Gx.pix(i,j));
		}
	}

	if (supresionNoMaximo)
		binariz.supresionNoMaximo(Gr2, imang);
	// La other opcion es usar adelgazar() pero es innecesario

	for (j=1; j<imorig.ly-1; j++)
	{
		for (i=1; i<imorig.lx-1; i++)
		{
			if (binariz.pix(i,j))
			{
				px.resize(px.size()+1);
				py.resize(py.size()+1);
				ang.resize(ang.size()+1);
				px[px.size()-1] = i;
				py[py.size()-1] = j;
				ang[ang.size()-1] = imang.pix(i,j);
				if (ang[ang.size()-1] < 0)
					ang[ang.size()-1] += M_PI;
			}
		}
	}
}

void L_ImageGrayDouble::houghDetLineas(const L_ImageGrayDouble &imorig, std::vector<double> &px, std::vector<double> &py, std::vector<double> &ang, int nThetas, int nRhos, bool usaGrad)
{
	reallocate(nThetas, nRhos);
	setZero();
	int i, r, t, t1, t2, t3, t4;
	double rho, theta, x, y;
	double errAng = 15*M_PI/180, diag = sqrt(imorig.lx*1.0*imorig.lx+imorig.ly*1.0*imorig.ly);
	double rangoRho = 2*diag;

	//    variable  rango-celda              rango-datos               0...2*diag
	//     rho      r = [1,...,nRho]     rho = [-diag,...,diag]    r = (rho+diag)/(2*diag)*(nRho-1)+1
	//     theta    t =[1,...,nTh]       theta=[0,...,pi]          t = theta/pi*(nTh-1)+1

	//        B-A          N
	//     x:[A,B] -> xE:[1,N]      (x-A)/(A+B) -> [0,1]     ....    (x-A)/(A+B) * (N-1) + 1     

	for (i=0; i<(int)px.size(); i++)
	{

		if (usaGrad)
		{
			t1 = (int)((ang[i]-errAng)/(M_PI)*nThetas);
			t2 = (int)((ang[i]+errAng)/(M_PI)*nThetas);
			t3 = 0;
			t4 = 0;
			// Usualmente ocurren uno de los dos casos
			if (t1 < 0) // El intervalo se sale por la izquierda
			{
				t1 = 0;
				t3 = t1 + nThetas;
				t4 = nThetas - 1;
			}
			if (t2 >= nThetas) // El intervalo se sale por la derecha
			{
				t2 = nThetas-1;
				t3 = 0;
				t4 = t2 - nThetas;
			}
		}
		else
		{
			t1 = 0;
			t2 = nThetas-1;
			t3 = 0;
			t4 = 0;
		}
		x = px[i];
		y = py[i];

		for (t=t1; t<t2; t++)
		{
			theta = t/(double)nThetas * M_PI;
			rho = x*cos(theta) + y*sin(theta);
			r = (int)((rho+diag) / rangoRho * nRhos);
			pix(t,r) = pix(t,r)+1;
		}
		for (t=t3; t<t4; t++)
		{
			theta = t/(double)nThetas * M_PI;
			rho = x*cos(theta) + y*sin(theta);
			r = (int)((rho+diag) / rangoRho * nRhos);
			pix(t,r) = pix(t,r)+1;
		}
	}
}

void L_ImageGrayDouble::houghDibLineas(const L_ImageGrayDouble &imorig, L_ShapeArray &lins, double umbral)
{
	int i, j;
	double nThetas, nRhos, m, n, minv, ninv, diag, rangoRho, ang, rho;
	nThetas = lx;
	nRhos = ly;
	diag = sqrt(imorig.lx*1.0*imorig.lx + imorig.ly*1.0*imorig.ly);
	rangoRho = 2*diag;
	lins.resize(0);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			if (pix(i,j) >= umbral)
			{
				ang = i/nThetas * M_PI; // lx = nRhos
				rho = j/nRhos * rangoRho - diag; // ly = nThetas
				if (ang < M_PI/2 || ang > 3/4*M_PI/2) // Recta horizontal
				{
					// rho = x cos(theta) + y sin(theta)
					// y = rho/sin(theta) - x*cos(theta)/sin(theta)
					m = -cos(ang)/sin(ang);
					n = rho / sin(ang);

					// for i=0...height
					//   j = m*i+b;
					//   if (j>=0 & j<width)
					//       im(i,j) = 255;
					//   end
					//  end

					lins.drawLine(0, (int)(m*0+n), imorig.lx, (int)(m*imorig.lx+n));
				}
				else
				{
					ang = M_PI-ang; // Cambiar los ejes : x = minv*y + ninv
					minv = -cos(ang)/sin(ang);
					ninv = rho / sin(ang);
					lins.drawLine((int)(minv*0+ninv), 0, (int)(minv*imorig.ly+ninv), imorig.ly);
				}
			}
}

L_ImageGrayDouble& L_ImageGrayDouble::conv2d(const L_ImageGrayDouble &other, const L_ImageGrayDouble &mask, int xCenMasc, int yCenMasc)
{
	int i, j, u, v;
	reallocate(other.lx, other.ly);
	for (i=0; i<lx; i++)
		for (j=0; j<ly; j++)
			pix(i,j)=0;

	//  0  1  ... xCenMasc ... ... lx-2  lx-1
	//  <--- xC+1 ---> <-----  lx-xC  ---->

	for (j=yCenMasc; j<ly-(mask.ly-yCenMasc); j++)
	{
		for (i=xCenMasc; i<lx-(mask.lx-xCenMasc); i++)
		{
			//pix(i,j)=0;
			for (u=0; u<mask.lx; u++)
				for (v=0; v<mask.ly; v++)
					pix(i,j)+=other.pix(i+u-xCenMasc,j+v-yCenMasc)*mask.pix(u,v);
		}
	}
	// De ahi veo los bordes...
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::convXSim(const L_ImageGrayDouble &other, const std::vector<double> &h)
{
	long i,j;
	int k;
	int n = (int)h.size();

	double sum;
	int u,v;
	int s3;

	reallocate(other.lx, other.ly);

	s3=n/2;

	for (j=0; j<ly; j++)
	{
		for (i=0; i<s3; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
 		}

		for (i=s3; i<lx-s3; i++)
		{
			sum=0;
			for (k=i-s3; k<=i+s3; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
		}
		for (i=lx-s3; i<lx; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::convYSim(const L_ImageGrayDouble &other, const std::vector<double> &h)
{
	long i,j;
	int k;

	double sum;
	int u,v;
	int s3;
	int n = (int)h.size();

	reallocate(other.lx, other.ly);

	s3=n/2;

	for (i=0; i<lx; i++)
	{
		for (j=0; j<s3; j++)
		{

			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (j=s3; j<ly-s3; j++)
		{
			sum=0;
			for (k=j-s3; k<=j+s3; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (j=ly-s3; j<ly; j++)
		{

			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}


L_ImageGrayDouble &L_ImageGrayDouble::convXSim_sub(const L_ImageGrayDouble &other, const std::vector<double> &h, int x0, int y0, int rx, int ry)
{
	long i,j;
	int k;
	int iIni, jIni;
	int iFin, jFin;

	double sum;
	int u,v;
	int s3;
	int width = (int)h.size();

	throw_L_ArgException_if(other.lx < width || other.ly < width, "L_ImageGrayDouble::convXSim_sub() : imagen muy pequena");

	reallocate(other.lx, other.ly);

	s3=width/2;

	iIni = x0 - rx;
	iFin = x0 + rx + 1;
	jIni = y0 - ry;
	jFin = y0 + ry + 1;

	if (iIni < 0)
		iIni = 0;
	if (iFin > lx)
		iFin = lx-1;
	if (jIni < 0)
		jIni = 0;
	if (jFin > ly)
		jFin = ly-1;

	for (j=jIni; j<jFin; j++)
	{
		for (i=iIni; i<s3; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
 		}

		for (i=s3; i<iFin-s3; i++)
		{
			sum=0;
			for (k=i-s3; k<=i+s3; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
		}
		i=iFin-s3;
		if (i<0)
			i=0;
		for ( ; i<iFin; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(k,j)*h[i-k+s3];
			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}

L_ImageGrayDouble &L_ImageGrayDouble::convYSim_sub(const L_ImageGrayDouble &other, const std::vector<double> &h, int x0, int y0, int rx, int ry)
{
	long i,j, iIni, jIni, iFin, jFin;
	int k;

	double sum;
	int u,v;
	int s3;

	int width = (int)h.size();

	throw_L_ArgException_if(other.lx < width || other.ly < width, "L_ImageGrayDouble::convYSim_sub() : imagen muy pequena");

	reallocate(other.lx, other.ly);

	s3=width/2;

	iIni = x0 - rx;
	iFin = x0 + rx + 1;
	jIni = y0 - ry;
	jFin = y0 + ry + 1;

	if (iIni < 0)
		iIni = 0;
	if (iFin > lx)
		iFin = lx-1;
	if (jIni < 0)
		jIni = 0;
	if (jFin > ly)
		jFin = ly-1;

	for (i=iIni; i<iFin; i++)
	{
		for (j=jIni; j<s3; j++)
		{

			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (j=s3; j<jFin-s3; j++)
		{
			sum=0;
			for (k=j-s3; k<=j+s3; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		j = jFin-s3;
		if (j<0)
			j=0;
		for ( ; j<jFin; j++)
		{
			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k++)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}

// static
void L_ImageGrayDouble::smoothYSqrt2(const L_ImageGrayDouble &imEntr, L_ImageGrayDouble &res)
{
	int lx, ly, i, j;
	const double **elem;
	double **relem, pc;

	res.reallocate(imEntr.lx, imEntr.ly);

	lx = imEntr.lx;
	ly = imEntr.ly;
	elem = L_new2d_ext<const double>(imEntr.li, imEntr.lj, imEntr.data(), imEntr.ljStep);
	relem = L_new2d_ext<double>(res.li, res.lj, res.data(), res.ljStep);;

	for (j = 0; j < ly; j++)
	{
		for (i = 3; i < lx - 3; i++)
		{
			relem[i][j] = 0.03001 * elem[i-3][j] + 0.10502 * elem[i-2][j] +
			0.22201 * elem[i-1][j] + 0.28602 * elem[i][j] + 0.22201 * elem[i+1][j] +
			0.10502 * elem[i+2][j] + 0.03001 * elem[i+3][j];
		}
		for (i = 0; i < 3; i++)
		{
			pc = (i < 1) ? elem[0][j] : elem[i-1][j];
			relem[i][j] = 0.135 * elem[0][j] + 0.22201 * pc +
			0.28602 * elem[i][j] + 0.22201 * elem[i+1][j] +
			0.10502 * elem[i+2][j] + 0.03001 * elem[i+3][j];
		}
		for (i = lx - 3; i < lx; i++)
		{
			pc = (i >= lx - 1) ? elem[lx-1][j] : elem[i+1][j];
			relem[i][j] = 0.03001 * elem[i-3][j] + 0.10502 * elem[i-2][j] +
			0.22201 * elem[i-1][j] + 0.28602 * elem[i][j] + 0.22201 * pc +
			0.13502 * elem[lx-1][j];
		}
	}
	L_delete2d_ext<const double>(elem);
	L_delete2d_ext<double>(relem);
}

// static
void L_ImageGrayDouble::smoothXSqrt2(const L_ImageGrayDouble &imEntr, L_ImageGrayDouble &res)
{
	int lx, ly, i, j;
	const double **elem;
	double **relem, pc;

	res.reallocate(imEntr.lx, imEntr.ly);

	lx = imEntr.lx;
	ly = imEntr.ly;
	elem = L_new2d_ext<const double>(imEntr.li, imEntr.lj, imEntr.data(), imEntr.ljStep);
	relem = L_new2d_ext<double>(res.li, res.lj, res.data(), res.ljStep);

	for (i = 0; i < lx; i++)
	{
		for (j = 3; j < ly - 3; j++)
		{
			relem[i][j] = 0.03001 * elem[i][j-3] + 0.10503 * elem[i][j-2] + 0.22202 * elem[i][j-1] +
			0.28601 * elem[i][j] + 0.22202 * elem[i][j+1] + 0.10503 * elem[i][j+2] + 0.03001 * elem[i][j+3];
		}
		for (j = 0; j < 3; j++)
		{
			pc = (j < 1) ? elem[i][0] : elem[i][j-1];
			relem[i][j] = 0.13501 * elem[i][0] + 0.22202 * pc +
			0.28601 * elem[i][j] + 0.22202 * elem[i][j+1] + 0.10501 * elem[i][j+2] + 0.03001 * elem[i][j+3];
		}
		for (j = ly - 3; j < ly; j++)
		{
			pc = (j >= ly - 1) ? elem[i][ly-1] : elem[i][j+1];
			relem[i][j] = 0.03001 * elem[i][j-3] + 0.10503 * elem[i][j-2] + 0.22202 * elem[i][j-1] +
			0.28601 * elem[i][j] + 0.22202 * pc + 0.13503 * elem[i][ly-1];
		}
	}
	L_delete2d_ext<const double>(elem);
	L_delete2d_ext<double>(relem);
}


L_ImageGrayDouble& L_ImageGrayDouble::convXSim_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso)
{
	long i,j;
	long k;

	double sum;
	long u,v;
	long s3;

	int n = (int)h.size();

	reallocate(other.lx, other.ly);

	s3=n/2;

	for (j=0; j<ly; j++)
	{
		for (i=0; i<s3; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k+=paso)
				sum+=other.pix(k,j)*h[i-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (i=s3; i<lx-s3; i++)
		{
			sum=0;
			for (k=i-s3; k<=i+s3; k+=paso)
				sum+=other.pix(k,j)*h[i-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (i=lx-s3; i<lx; i++)
		{
			u=i-s3;
			v=i+s3;
			if (u<0)
				u=0;
			if (v>=lx)
				v=lx-1;
			sum=0;
			for (k=u; k<=v; k+=paso)
				sum+=other.pix(k,j)*h[i-k+s3];
 			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::convYSim_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso)
{
	long i,j;
	long k;

	double sum;
	long u,v;
	long s3;

	int n = (int)h.size();

	reallocate(other.lx, other.ly);

	s3=n/2;

	for (i=0; i<lx; i++)
	{
		for (j=0; j<s3; j++)
		{
			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k+=paso)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (j=s3; j<ly-s3; j++)
		{
			sum=0;
			for (k=j-s3; k<=j+s3; k+=paso)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
		for (j=ly-s3; j<ly; j++)
		{
			u=j-s3;
			v=j+s3;
			if (u<0)
				u=0;
			if (v>=ly)
				v=ly-1;
			sum=0;
			for (k=u; k<=v; k+=paso)
				sum+=other.pix(i,k)*h[j-k+s3];
 			pix(i,j)=(double)sum;
 		}
	}
	return *this;
}



L_ImageGrayDouble& L_ImageGrayDouble::gradRight_3(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			pix(i,j)=(other.pix(i+1,j)-other.pix(i-1,j))/2;
		}
	}
   	for (j=0; j<ly; j++)
	{
		pix(0,j)=pix(1,j);
		pix(lx-1,j)=pix(lx-2,j);
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradUp_3(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=1; j<ly-1; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j)=(other.pix(i,j-1)-other.pix(i,j+1))/2;
		}
	}
   	for (i=0; i<lx; i++)
	{
		pix(i,0)=pix(i,1);
		pix(i,ly-1)=pix(i,ly-2);
	}
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::gradRight_3_sub(const L_ImageGrayDouble &other, int x0, int y0, int rx, int ry)
{
	int i, j;
	int iIni, jIni, iFin, jFin;
	reallocate(other.lx, other.ly);

	iIni = x0-rx;
	iFin = x0+rx+1;
	jIni = y0-ry;
	jFin = y0+ry+1;

	if (iIni < 1)
		iIni = 0;
	if (iFin > lx-1)
		iFin = lx-1;
	if (jIni < 0)
		jIni = 0;
	if (jFin > ly)
		jFin = ly;

	for (j=jIni; j<jFin; j++)
	{
		for (i=iIni+1; i<iFin-1; i++)
		{
			pix(i,j)=(other.pix(i+1,j)-other.pix(i-1,j))/2;
		}
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradUp_3_sub(const L_ImageGrayDouble &other, int x0, int y0, int rx, int ry)
{
	int i, j;
	int iIni, jIni, iFin, jFin;
	reallocate(other.lx, other.ly);

	iIni = x0-rx;
	iFin = x0+rx+1;
	jIni = y0-ry;
	jFin = y0+ry+1;

	if (iIni < 0)
		iIni = 0;
	if (iFin > lx)
		iFin = lx;
	if (jIni < 1)
		jIni = 1;
	if (jFin > ly-1)
		jFin = ly-1;

	for (j=jIni+1; j<jFin-1; j++)
	{
		for (i=iIni; i<iFin; i++)
		{
			pix(i,j)=(other.pix(i,j-1)-other.pix(i,j+1))/2;
		}
	}
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::gradRight_2(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);

	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx-1; i++)
		{
			pix(i,j)=other.pix(i+1,j)-other.pix(i,j);
		}
	}
   	for (j=0; j<ly; j++)
	{
		pix(lx-1,j)=pix(lx-2,j);
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradUp_2(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);

	for (j=0; j<ly-1; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j)=other.pix(i,j)-other.pix(i,j+1);
		}
	}
   	for (i=0; i<lx; i++)
	{
		pix(i,ly-1)=pix(i,ly-2);
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradRightUp_2(L_ImageGrayDouble &gradArr, const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	gradArr.reallocate(other.lx, other.ly);

	for (j=0; j<ly-1; j++)
	{
		for (i=0; i<lx-1; i++)
		{
			pix(i,j)=other.pix(i+1,j)-other.pix(i,j);
			gradArr.pix(i,j) = other.pix(i,j)-other.pix(i,j+1);
		}
	}
   	for (i=0; i<lx; i++)
	{
		pix(i,ly-1)=pix(i,ly-2);
		gradArr.pix(i,ly-1)=gradArr.pix(i,ly-2);
	}
   	for (j=0; j<ly; j++)
	{
		pix(lx-1,j)=pix(lx-2,j);
		gradArr.pix(lx-1,j)=gradArr.pix(lx-2,j);
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradR_2(L_ImageGrayDouble &gradAng, const L_ImageGrayDouble &other)
{
	register int i, j;
	register double gradDer, gradArr;
	reallocate(other.lx, other.ly);
	gradAng.reallocate(other.lx, other.ly);

	for (j=0; j<ly-1; j++)
	{
		for (i=0; i<lx-1; i++)
		{
			gradDer=other.pix(i+1,j)-other.pix(i,j);
			gradArr=other.pix(i,j)-other.pix(i,j+1);
			if (gradDer==0 && gradArr==0)
			{
				pix(i,j)=0;
				gradAng.pix(i,j)=0;
			}
			else
			{
				pix(i,j)=(double)sqrt(gradDer*gradDer+gradArr*gradArr);
				if (gradArr==0 && gradDer==0)
					gradAng.pix(i,j)=0;
				else
					gradAng.pix(i,j)=(double)atan2(gradArr,gradDer);
			}
		}
	}
	for (i=0; i<lx-1; i++)
	{
		pix(i,ly-1)=pix(i,ly-2);
		gradAng.pix(i,ly-1)=gradAng.pix(i,ly-2);
	}
	for (j=0; j<ly-1; j++)
	{
		pix(lx-1,j)=pix(lx-2,j);
		gradAng.pix(lx-1,j)=gradAng.pix(lx-2,j);
	}
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::gradR_3(L_ImageGrayDouble &gradAng, const L_ImageGrayDouble &other)
{
	int i, j;
	double gradDer, gradArr;
	reallocate(other.lx, other.ly);
	gradAng.reallocate(other.lx, other.ly);

	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			gradDer=(other.pix(i+1,j)-other.pix(i-1,j))/2;
			gradArr=(other.pix(i,j-1)-other.pix(i,j+1))/2;
			pix(i,j)=(double)sqrt(gradDer*gradDer+gradArr*gradArr);
			if (gradArr==0 && gradDer==0)
				gradAng.pix(i,j)=0;
			else
				gradAng.pix(i,j)=(double)atan2(gradArr,gradDer);
		}
	}
	for (i=0; i<lx; i++)
	{
		pix(i,0)=pix(i,1);
		gradAng.pix(i,0)=gradAng.pix(i,1);
		pix(i,ly-1)=pix(i,ly-2);
		gradAng.pix(i,ly-1)=gradAng.pix(i,ly-2);
	}
	for (j=0; j<ly; j++)
	{
		pix(0,j)=pix(1,j);
		gradAng.pix(0,j)=gradAng.pix(1,j);
		pix(lx-1,j)=pix(lx-2,j);
		gradAng.pix(lx-1,j)=gradAng.pix(lx-2,j);
	}
	return *this;
}

void L_ImageGrayDouble::integralOf(const L_ImageGrayDouble &other)
{
	int i, j;
	double sFila = 0;
	reallocate_fc(other.li,other.lj);
	operator()(0,0) = other(0,0);
	for (j=1; j<lj; j++)
		operator()(0,j)=operator()(0,j-1)+other(0,j);
	for (i=1; i<li; i++)
	{
		sFila = 0;
		for (j=0; j<lj; j++)
		{
			sFila += other(i,j);
			operator()(i,j) = operator()(i-1,j) + sFila;
		}
	}
}

L_ImageGrayDouble& L_ImageGrayDouble::laplac(const L_ImageGrayDouble &other)
{
	int i, j=1;
	reallocate(other.lx, other.ly);

	for (j=1; j<ly-1; j++)
		for (i=1; i<lx-1; i++)
			pix(i,j)=(double)0.5*other.pix(i-1,j)+(double)0.5*other.pix(i+1,j)+(double)0.5*other.pix(i,j-1)+(double)0.5*other.pix(i,j+1)-2*other.pix(i,j);
	for (i=1; i<lx-1; i++)
	{
		pix(i,0)=(double)0.5*other.pix(i-1,0)+(double)0.5*other.pix(i+1,j)-other.pix(i,0);
		pix(i,ly-1)=(double)0.5*other.pix(i-1,ly-1)+(double)0.5*other.pix(i+1,ly-1)-other.pix(i,ly-1);
	}
	for (j=1; j<ly-1; j++)
	{
		pix(0,j)=(double)0.5*other.pix(0,j-1)+(double)0.5*other.pix(0,j+1)-other.pix(0,j);
		pix(lx-1,j)=(double)0.5*other.pix(lx-1,j-1)+(double)0.5*other.pix(lx-1,j+1)-other.pix(lx-1,j);
	}
	pix(0,0)=0;
	pix(0,ly-1)=0;
	pix(lx-1,0)=0;
	pix(lx-1,ly-1)=0;
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::z_subm50(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx/2, other.ly/2);

	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=other.pix(2*i,2*j);
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::z_subm67(const L_ImageGrayDouble &other)
{
	register int u, v;
	register double _temp1, _temp2, _tempc;
	int lxB,lyB;
	lxB=other.lx/3;
	lyB=other.ly/3;
	reallocate((other.lx*2)/3, (other.ly*2)/3);
	for (v=0; v<lyB; v++)
	{
		for (u=0; u<lxB; u++)
		{
			pix(2*u,2*v)=other.pix(3*u,3*v)*(double)(4.0/9);
			pix(2*u+1,2*v)=other.pix(3*u+2,3*v)*(double)(4.0/9);
			pix(2*u,2*v+1)=other.pix(3*u,3*v+2)*(double)(4.0/9);
			pix(2*u+1,2*v+1)=other.pix(3*u+2,3*v+2)*(double)(4.0/9);
			_temp1=other.pix(3*u,3*v+1)*(double)(2.0/9);
			_temp2=other.pix(3*u+1,3*v)*(double)(2.0/9);
			_tempc=other.pix(3*u+1,3*v+1)*(double)(1.0/9);
			pix(2*u,2*v)+=_temp1+_temp2+_tempc;
			_temp1=other.pix(3*u+2,3*v+1)*(double)(2.0/9);
			pix(2*u+1,2*v)+=_temp1+_temp2+_tempc;
			_temp2=other.pix(3*u+1,3*v+2)*(double)(2.0/9);
			pix(2*u+1,2*v+1)+=_temp1+_temp2+_tempc;
			_temp1=other.pix(3*u,3*v+1)*(double)(2.0/9);
			pix(2*u,2*v+1)+=_temp1+_temp2+_tempc;
		}
	}
	if (other.lx%3==2)
	{
		u=lxB;
		for (v=0; v<lyB; v++)
		{
			_tempc=other.pix(3*u+1,3*v+1)*(double)(1.0/9)+other.pix(3*u,3*v+1)*(double)(2.0/9);
			pix(2*u,2*v)=other.pix(3*u,3*v)*(double)(4.0/9)+other.pix(3*u+1,3*v)*(double)(2.0/9)+_tempc;
			pix(2*u,2*v+1)=other.pix(3*u,3*v+2)*(double)(4.0/9)+other.pix(3*u+1,3*v+2)*(double)(2.0/9)+_tempc;
		}
	}
	if (other.ly%3==2)
	{
		v=lyB;
		for (u=0; u<lxB; u++)
		{
			_tempc=other.pix(3*u+1,3*v+1)*(double)(1.0/9)+other.pix(3*u+1,3*v)*(double)(2.0/9);
			pix(2*u,2*v)=other.pix(3*u,3*v)*(double)(4.0/9)+other.pix(3*u,3*v+1)*(double)(2.0/9)+_tempc;
			pix(2*u+1,2*v)=other.pix(3*u+2,3*v)*(double)(4.0/9)+other.pix(3*u+2,3*v+1)*(double)(2.0/9)+_tempc;
		}
	}
	if (other.lx%3==2 && other.ly%3==2)
	{
		u=lxB;
		v=lyB;
		pix(2*u,2*v)=other.pix(3*u,3*v)*(double)(4.0/9)+other.pix(2*u+1,2*v)*(double)(2.0/9)+other.pix(2*u,2*v+1)*(double)(2.0/9)+other.pix(2*u+1,2*v+1)*(double)(1.0/9);
	}
	return *this;
}
void L_ImageGrayDouble::doubleResolution(const L_ImageGrayDouble &other)
{
	int i, j;
	int u, v;
	reallocate(2*other.lx, 2*other.ly);
	for (j=0; j<other.ly; j++)
		for (i=0; i<other.lx; i++)
		{
			u=2*i;
			v=2*j;
			pix(u,v)=pix(u+1,v)=pix(u,v+1)=pix(u+1,v+1)=other.pix(i,j);
		}
}

void L_ImageGrayDouble::doubleResolutionBilineal(const L_ImageGrayDouble &other)
{
	int i, j;
	int u, v;

	reallocate(2*other.lx-1, 2*other.ly-1);
	for (i=0; i<other.lx-1; i++)
	{
		u=2*i;
		pix(u,0)=other.pix(i,0);
		pix(u+1,0)=(other.pix(i,0)+other.pix(i+1,0))/2;
		pix(u,1)=(other.pix(i,0)+other.pix(i,1))/2;
		for (j=1; j<other.ly-1; j++)
		{
			v=2*j;
			pix(u,v)=other.pix(i,j);
			pix(u+1,v)=(other.pix(i,j)+other.pix(i+1,j))/2;
			pix(u,v+1)=(other.pix(i,j)+other.pix(i,j+1))/2;
			pix(u+1,v-1)=(other.pix(u+1,v-2)+other.pix(u+1,v))/2;
		}
		pix(u,ly-1)=other.pix(i,other.ly-1);
		pix(u+1,ly-1)=(other.pix(i,other.ly-1)+other.pix(i+1,other.ly-1))/2;
		pix(u+1,ly-2)=(other.pix(u+1,ly-3)+other.pix(u+1,ly-1))/2;
	}
	for (j=0; j<other.lx-1; j++)
	{
		v=2*j;
		pix(lx-1,v)=other.pix(other.lx-1,j);
		pix(lx-1,v+1)=(other.pix(other.lx-1,j)+other.pix(other.lx-1,j+1))/2;
	}
}

void L_ImageGrayDouble::subImageOf(const L_ImageGrayDouble &imOrig, int x0, int y0, int x1, int y1)
{
	int i, j;
	int temp;
	if (x0 > x1)
	{
		temp=x0;
		x0=x1;
		x1=temp;
	}
	if (y0 > y1)
	{
		temp=y0;
		y0=y1;
		y1=temp;
	}
	if (x0 < 0)
	{
		x0 = 0;
		if (x1<0)
			x1=0;
	}
	if (y0 < 0)
	{
		y0 = 0;
		if (y1<0)
			y1=0;
	}
	if (x1 >= imOrig.lx)
	{
		x1 = imOrig.lx-1;
		if (x0 >= imOrig.lx)
			x0 = imOrig.lx-1;
	}
	if (y1 >= imOrig.ly)
	{
		y1 = imOrig.ly-1;
		if (y0 >= imOrig.ly)
			y0 = imOrig.ly-1;
	}
	reallocate(x1-x0+1, y1-y0+1);
	for (j=y0; j<=y1; j++)
		for (i=x0; i<=x1; i++)
			pix(i-x0,j-y0)=imOrig.pix(i,j);
}

bool L_ImageGrayDouble::calculaHistGradiente(int lu, int lv, int lw, double ***hist)
{
	if (data() == NULL || lx < 2 || ly < 2 || hist == NULL || lu < 0 || lv < 0 || lw < 0 || lu > lx/2 || lv > ly/2)
	{
		printf("L_ImageGrayDouble::calculaHistGradiente() : parametros inadecuados\n");
		return false;
	}
	L_ImageGrayDouble grR, grT; // Gradiente como radio y angulo
	double gx, gy;
	int u, v, w;
	double uFRAC, vFRAC, wFRAC;
	int u0,v0,w0;
	double du, dv, dw, r;
	int i, j;
	grR.reallocate(lx-1,ly-1);
	grT.reallocate(lx-1,ly-1);
	for (j=0; j<ly-1; j++)
	{
		for (i=0; i<lx-1; i++)
		{
			gx = pix(i+1,j)-pix(i,j);
			gy = pix(i,j+1)-pix(i,j);
			grR.pix(i,j)=sqrt(gx*gx+gy*gy);
			grT.pix(i,j)=atan2(gy,gx);
		}
	}
	for (u=0; u<lu; u++)
		for (v=0; v<lv; v++)
			for (w=0; w<lw; w++)
				hist[u][v][w]=0;
	for (j=0; j<ly-1; j++)
	{
		for (i=0; i<lx-1; i++)
		{
			uFRAC = i*(double)lu/(lx-1);
			vFRAC = j*(double)lv/(ly-1);
			wFRAC = grT.pix(i,j)*(double)lw/(2*M_PI);
			u0 = (int)uFRAC;  // uFRAC - u0 entre 0 y 1
			v0 = (int)vFRAC;
			w0 = (int)wFRAC;
			du = uFRAC-u0;
			dv = vFRAC-v0;
			dw = wFRAC-w0;
			r = grR.pix(i,j);
			hist[u0][v0][w0] += (1-du)*(1-dv)*(1-dw)*r;
			hist[u0+1][v0][w0] += (du)*(1-dv)*(1-dw)*r;
			hist[u0][v0+1][w0] += (1-du)*(dv)*(1-dw)*r;
			hist[u0+1][v0+1][w0] += (du)*(dv)*(1-dw)*r;
			hist[u0][v0][w0+1] += (1-du)*(1-dv)*(dw)*r;
			hist[u0+1][v0][w0+1] += (du)*(1-dv)*(dw)*r;
			hist[u0][v0+1][w0+1] += (1-du)*(dv)*(dw)*r;
			hist[u0+1][v0+1][w0+1] += (du)*(dv)*(dw)*r;
		}
	}
	return true;
}


bool L_ImageGrayDouble::esMax2D_anillo(double val, int i, int j, double r)
{
	int d, nDir;
	double rot_cos, rot_sin, x, y, t;
	nDir=(int)(8*r);
	rot_cos=cos(2*M_PI/nDir);
	rot_sin=sin(2*M_PI/nDir);
	x=r;
	y=0;
	if (i+r >= lx || i-r < 0 || j+r >= ly || j-r < 0)
		return false;
	// rotar (x,y) en cada paso
	for (d=0; d<nDir; d++)
	{
		if ( val <= pix((int)x+i,(int)y+j) )
			return false;
		t = x*rot_cos  + y*rot_sin;
		y = x*-rot_sin + y*rot_cos;
		x=t;
	}
	return true;
}

bool L_ImageGrayDouble::esMin2D_anillo(double val, int i, int j, double r)
{
	int d, nDir;
	double rot_cos, rot_sin, x, y, t;
	nDir=(int)(8*r);
	rot_cos=cos(2*M_PI/nDir);
	rot_sin=sin(2*M_PI/nDir);
	x=r;
	y=0;
	if (i+r >= lx || i-r < 0 || j+r >= ly || j-r < 0)
		return false;
	// Rotar (x,y) en cada paso
	for (d=0; d<nDir; d++)
	{
		if ( val >= pix((int)x+i,(int)y+j) )
			return false;
		t = x*rot_cos  + y*rot_sin;
		y = x*-rot_sin + y*rot_cos;
		x=t;
	}
	return true;
}

#define L_TAM_COPIA_MEMCPY 100

L_ImageGrayDouble & L_ImageGrayDouble::operator=(const L_ImageGrayUchar &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=other.pix(i,j)/(double)255.0;
	return *this;
}
L_ImageGrayDouble & L_ImageGrayDouble::operator=(const L_ImageRGBDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(other.pix(i,j,0)+other.pix(i,j,1)+other.pix(i,j,2))/3;
	return *this;
}
L_ImageGrayDouble & L_ImageGrayDouble::operator=(const L_ImageRGBUchar &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=(other.pix(i,j,0)*(double)1.0+other.pix(i,j,1)*(double)1.0+other.pix(i,j,2))/(3*255);
	return *this;
}

bool L_ImageGrayUchar::difusion_pinta(int i, int j, L_uchar lum)
{
	L_ImageGrayUchar recursion;
	L_uchar L;
	bool seguir=true;

	if (lx < 1 || ly < 1 || data()==NULL || i<0 || i>=lx || j<0 || j>=ly)
		return false;

	recursion.reallocate(lx, ly);
	recursion.setConstant(0);
	L=pix(i,j);

	while(seguir)
	{
		if (pix(i,j)!=L)
			recursion.pix(i,j)=(recursion.pix(i,j)&0xF0) | 8; // No hay que pintar aca
		switch(recursion.pix(i,j)&0x0F)
		{
		case 0:
			recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 1;
			if (i<recursion.lx-1 && recursion.pix(i+1,j)==0)
			{
				i++;
				recursion.pix(i,j) = 16;
			}
			else
				recursion.pix(i,j) = recursion.pix(i,j)&0xF0 + 1;
			break;
		case 1:
			recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 2;
			if (j>0 && recursion.pix(i,j-1)==0)
			{
				j--;
				recursion.pix(i,j) = 32;
			}
			else
				recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 2;
			break;
		case 2:
			recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 4;
			if (i>0 && recursion.pix(i-1,j)==0)
			{
				i--;
				recursion.pix(i,j) = 64;
			}
			else
				recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 4;
			break;
		case 4:
			recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 8;
			if (j<recursion.ly-1 && recursion.pix(i,j+1)==0)
			{
				j++;
				recursion.pix(i,j) = 128;
			}
			else
				recursion.pix(i,j) = (recursion.pix(i,j)&0xF0) | 8;
			break;
		case 8:
			pix(i,j)=lum;
			switch(recursion.pix(i,j)&0xF0)
			{
			case 0:
				seguir=false;
				break;
			case 16:
				i--;
				break;
			case 32:
				j++;
				break;
			case 64:
				i++;
				break;
			case 128:
				j--;
				break;
			}
			break;
		}
	}

	return true;
}


void L_ImageGrayUchar::supresionNoMaximo(const L_ImageGrayDouble &gradR, const L_ImageGrayDouble &gradA)
{
	int i, j;
	double ang;
	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			if (!pix(i,j))
				continue;
			ang = (gradA.pix(i,j)-M_PI/8) / (M_PI/4);
			while (ang < 0)
				ang += 4;
			while (ang >= 4)
				ang -= 4;
			switch ((int)ang)
			{
			case 0:
				pix(i,j) = 255 * (gradR.pix(i,j) > gradR.pix(i-1,j) && gradR.pix(i,j) > gradR.pix(i+1,j));
				break;
			case 1:
				pix(i,j) = 255 * (gradR.pix(i,j) > gradR.pix(i-1,j-1) && gradR.pix(i,j) > gradR.pix(i+1,j+1));
				break;
			case 2:
				pix(i,j) = 255 * (gradR.pix(i,j) > gradR.pix(i,j-1) && gradR.pix(i,j) > gradR.pix(i,j+1));
				break;
			case 3:
				pix(i,j) = 255 * (gradR.pix(i,j) > gradR.pix(i-1,j+1) && gradR.pix(i,j) > gradR.pix(i+1,j-1));
				break;
			}
		}
	}
}


void L_ImageGrayUchar::morfol_eliminaPolvo(const L_ImageGrayUchar &orig)
{
	int i, j;
	reallocate(orig.lx, orig.ly);
	for (j=1; j<orig.ly-1; j++)
		for (i=1; i<orig.lx-1; i++)
			pix(i,j) = (orig.pix(i,j) && (
				orig.pix(i-1,j-1) || orig.pix(i-1,j) || orig.pix(i-1,j+1) ||
				orig.pix(i,j-1)   ||                     orig.pix(i,j+1)   ||
				orig.pix(i+1,j-1) || orig.pix(i+1,j) || orig.pix(i+1,j+1)    )*255);
}

int L_ImageGrayUchar::morfol_hitAndMiss(const L_ImageGrayUchar &orig, const char *mascara3x3)
{
	int i, j;
	int n = 0;
	reallocate(orig.lx, orig.ly);
#define L_IBNx_rev(i,j,u,v)  (mascara3x3[u+1 + (v+1)*3]==2 || (!orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==0) || (orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==1))
	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			pix(i,j) =
			L_IBNx_rev(i,j,-1,-1) && L_IBNx_rev(i,j,0,-1) && L_IBNx_rev(i,j,1,-1) && 
			L_IBNx_rev(i,j,-1,0) && L_IBNx_rev(i,j,0,0) && L_IBNx_rev(i,j,1,0) &&
			L_IBNx_rev(i,j,-1,1) && L_IBNx_rev(i,j,0,1) && L_IBNx_rev(i,j,1,1);
			pix(i,j) *= 255;
			n += (!orig.pix(i,j) && pix(i,j)) || (orig.pix(i,j) && !pix(i,j)); // Cambios de estado
		}
	}
#undef L_IBNx_rev
	return n;
}

int L_ImageGrayUchar::morfol_hitAndMiss_resta(const L_ImageGrayUchar &orig, const char *mascara3x3)
{
	int i, j;
	int n = 0;
	reallocate(orig.lx, orig.ly);
#define L_IBNx_rev(i,j,u,v)  (mascara3x3[u+1 + (v+1)*3]==2 || (!orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==0) || (orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==1))
	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			pix(i,j) = orig.pix(i,j) && !(
			L_IBNx_rev(i,j,-1,-1) && L_IBNx_rev(i,j,0,-1) && L_IBNx_rev(i,j,1,-1) && 
			L_IBNx_rev(i,j,-1,0) && L_IBNx_rev(i,j,0,0) && L_IBNx_rev(i,j,1,0) &&
			L_IBNx_rev(i,j,-1,1) && L_IBNx_rev(i,j,0,1) && L_IBNx_rev(i,j,1,1)
			);
			pix(i,j) *= 255;
			n += (!orig.pix(i,j) && pix(i,j)) || (orig.pix(i,j) && !pix(i,j)); // Cambios de estado
		}
	}
#undef L_IBNx_rev
	return n;
}

int L_ImageGrayUchar::morfol_erosiona(const L_ImageGrayUchar &orig, const char *mascara3x3)
{
	int i, j;
	int n = 0;
	reallocate(orig.lx, orig.ly);
#define L_IBNx_rev(i,j,u,v)  (orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==1)
	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			n += pix(i,j) =
			L_IBNx_rev(i,j,-1,-1) && L_IBNx_rev(i,j,0,-1) && L_IBNx_rev(i,j,1,-1) && 
			L_IBNx_rev(i,j,-1,0) && L_IBNx_rev(i,j,0,0) && L_IBNx_rev(i,j,1,0) &&
			L_IBNx_rev(i,j,-1,1) && L_IBNx_rev(i,j,0,1) && L_IBNx_rev(i,j,1,1);
			pix(i,j) *= 255;
		}
	}
#undef L_IBNx_rev
	return n;
}


int L_ImageGrayUchar::morfol_dilata(const L_ImageGrayUchar &orig, const char *mascara3x3)
{
	int i, j;
	int n = 0;
	reallocate(orig.lx, orig.ly);
#define L_IBNx_rev(i,j,u,v)  (orig.pix(i+u,j+v) && mascara3x3[u+1 + (v+1)*3]==1)
	for (j=1; j<ly-1; j++)
	{
		for (i=1; i<lx-1; i++)
		{
			n += pix(i,j) =
			L_IBNx_rev(i,j,-1,-1) || L_IBNx_rev(i,j,0,-1) || L_IBNx_rev(i,j,1,-1) || 
			L_IBNx_rev(i,j,-1,0) || L_IBNx_rev(i,j,0,0) || L_IBNx_rev(i,j,1,0) ||
			L_IBNx_rev(i,j,-1,1) || L_IBNx_rev(i,j,0,1) || L_IBNx_rev(i,j,1,1);
			pix(i,j)*=255;
		}
	}
#undef L_IBNx_rev
	return n;
}

int L_ImageGrayUchar::morfol_adelgaza(const L_ImageGrayUchar &orig)
{
	int n = 0, dn;
	L_ImageGrayUchar im1;
	// [0 0 0]  [  0 0]
	// [  1  ]  [1 1 0]
	// [1 1 1]  [  1  ]
	const char m[8][9] = {
		{0,0,0,2,1,2,1,1,1},
		{2,0,0,1,1,0,2,1,2},
		{1,2,0,1,1,0,1,2,0},
		{2,1,2,1,1,0,2,0,0},
		{1,1,1,2,1,2,0,0,0},
		{2,1,2,0,1,1,0,0,2},
		{0,2,1,0,1,1,0,2,1},
		{0,0,2,0,1,1,2,1,2}};
	*this = orig;
	dn = 0;
	*this = orig;
	im1 = orig;
	while(dn < 8)
	{
		intercambiaObjetoCon(im1);
		if (morfol_hitAndMiss_resta(im1, m[n%8]) == 0)
			dn++;
		else
			dn = 0;
		n++;
	}
	return n;
}

void L_ImageGrayUchar::OP_max(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2)
{
	int i, j;
	throw_L_ArgException_if(im1.lx != im2.lx || im2.ly != im2.ly, "L_ImageGrayUchar::OP_max");
	reallocate(im1.lx, im1.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j) = (im1.pix(i,j) > im2.pix(i,j)) * im1.pix(i,j) + (im2.pix(i,j)-im1.pix(i,j));
}

void L_ImageGrayUchar::OP_min(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2)
{
	int i, j;
	throw_L_ArgException_if(im1.lx != im2.lx || im2.ly != im2.ly, "L_ImageGrayUchar::OP_min");
	reallocate(im1.lx, im1.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j) = (im1.pix(i,j) < im2.pix(i,j)) * im1.pix(i,j) + (im2.pix(i,j)-im1.pix(i,j));
}


void L_ImageGrayUchar::OP_max(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2, L_ImageGrayUchar &im3)
{
	int i, j;
	throw_L_ArgException_if(im1.lx != im2.lx || im2.ly != im2.ly, "L_ImageGrayUchar::OP_max");
	reallocate(im1.lx, im1.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			if (im1.pix(i,j) > im2.pix(i,j))
				if (im1.pix(i,j) > im3.pix(i,j)) // im1 > im2, im1 > im3
					pix(i,j) = im1.pix(i,j);
				else // im1 > im2, im3 > im1
					pix(i,j) = im3.pix(i,j);
			else // im2 > im1
				if (im2.pix(i,j) > im3.pix(i,j))
					pix(i,j) = im2.pix(i,j); // im2 > im1, im2 > im3
				else // im2 > im1, im3 > im2
					pix(i,j) = im3.pix(i,j);
}

void L_ImageGrayUchar::OP_min(L_ImageGrayUchar &im1, L_ImageGrayUchar &im2, L_ImageGrayUchar &im3)
{
	int i, j;
	throw_L_ArgException_if(im1.lx != im2.lx || im2.ly != im2.ly, "L_ImageGrayUchar::OP_min");
	reallocate(im1.lx, im1.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			if (im1.pix(i,j) < im2.pix(i,j))
				if (im1.pix(i,j) < im3.pix(i,j)) // im1 < im2, im1 < im3
					pix(i,j) = im1.pix(i,j);
				else // im1 < im2, im3 < im1
					pix(i,j) = im3.pix(i,j);
			else // im2 < im1
				if (im2.pix(i,j) < im3.pix(i,j))
					pix(i,j) = im2.pix(i,j); // im2 < im1, im2 < im3
				else // im2 < im1, im3 < im2
					pix(i,j) = im3.pix(i,j);
}


void L_ImageGrayUchar::encodeToBuffer(std::vector<char> &buf) const
{
	buf.resize(lx*ly);
	copyToBuffer((L_uchar *)&buf[0]);
}

void L_ImageGrayUchar::decodeFromBuffer(const std::vector<char> &buf)
{
	throw_L_ArgException_if(buf.size()!=lx*ly, "L_ImageGrayUchar::decodeFromBuffer() : tamano de imagen distinto al buffer");
	copyFromBuffer((L_uchar *)&buf[0]);
}

void L_ImageGrayUchar::encodeDelta(std::vector<char> &bufD) const
{
	int i, j;
	// invertido lx y ly
	bufD.resize(lx*ly);
	bufD[0] = pix(0,0);
	for (i=1; i<ly; i++)
		bufD[i*lxStep+0] = operator()(i,0) - operator()(i-1,0);
	for (i=0; i<ly; i++)
		for (j=1; j<lx; j++)
			bufD[i*lxStep+j+0] = operator()(i,j) - operator()(i,j-1);
}
void L_ImageGrayUchar::decodeDelta(std::vector<char> &bufD)
{
	int i, j;
	throw_L_ArgException_if(bufD.size() != lx*ly, "L_ImageRGBUchar::decodeDelta()");
	operator()(0,0) = bufD[0];
	for (i=1; i<ly; i++)
		operator()(i,0) = operator()(i-1,0) + bufD[i*lxStep+0];
	for (i=0; i<ly; i++)
		for (j=1; j<lx; j++)
			operator()(i,j) = operator()(i,j-1) + bufD[i*lxStep+j+0];
}


bool L_ImageGrayUchar::normalizeHistogram()
{
	int i, j;
	L_uchar max, min;
	L_uchar dif;
	bool ret=true;
	max=min=pix(0,0);

	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			if (pix(i,j)<min)
				min=pix(i,j);
			else if (pix(i,j)>max)
				max=pix(i,j);
		}
	}
	dif=max-min;
	if (dif==0)
	{
		dif = 1;
		ret = false;
	}
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j)-=min;
			pix(i,j)=(L_uchar)(pix(i,j)*1.0/dif);
		}
	}
	return ret;
}

L_ImageGrayUchar& L_ImageGrayUchar::subm50(const L_ImageGrayUchar &other)
{
	int i, j;
	reallocate(other.lx/2, other.ly/2);

	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			pix(i,j)=other.pix(2*i,2*j);
	return *this;
}


L_ImageGrayUchar & L_ImageGrayUchar::operator=(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			if (other.pix(i,j)<0)
				pix(i,j)=0;
			else if(other.pix(i,j)>=1)
				pix(i,j)=255;
			else
				pix(i,j)=(L_uchar)(other.pix(i,j)*255);
		}
	return *this;
}

L_ImageGrayUchar & L_ImageGrayUchar::operator=(const L_ImageRGBDouble &other)
{
	int i, j;
	double sumita;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			sumita=(other.pix(i,j,0)+other.pix(i,j,1)+other.pix(i,j,2))/3;
			if (sumita<0)
				pix(i,j)=0;
			else if(sumita>=1)
				pix(i,j)=255;
			else
				pix(i,j)=(L_uchar)(sumita*255);
		}
	return *this;
}
L_ImageGrayUchar & L_ImageGrayUchar::operator=(const L_ImageRGBUchar &other)
{
	int i, j;
	unsigned int sumita;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			sumita=((((unsigned int)other.pix(i,j,0))+other.pix(i,j,1))+other.pix(i,j,2))/3;
			pix(i,j)=(L_uchar)sumita;
		}
	return *this;
}

bool L_ImageRGBDouble::normalizeHistogram()
{
	int i, j, c;
	double max, min;
	double dif;
	bool ret = true;
	max=min=pix(0,0,0);

	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			for (c=0; c<3; c++)
			{
				if (pix(i,j,c)<min)
					min=pix(i,j,c);
				else if (pix(i,j,c)>max)
					max=pix(i,j,c);
			}
		}
	}
	dif=max-min;
	if (dif==0)
	{
		dif = 1;
		ret = false;
	}
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			for (c=0; c<3; c++)
			{
				pix(i,j,c)-=min;
				pix(i,j,c)/=dif;
			}
		}
	}
	return ret;
}

L_ImageRGBDouble& L_ImageRGBDouble::transpOf(const L_ImageRGBDouble &im)
{
	int i, j;
	reallocate(im.ly, im.lx);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=im.pix(j,i,0);
			pix(i,j,1)=im.pix(j,i,1);
			pix(i,j,2)=im.pix(j,i,2);
		}
	}
	return *this;
}



L_ImageGrayDouble& L_ImageRGBDouble::copyR_into(L_ImageGrayDouble &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,0);
	return imC;
}
L_ImageGrayDouble& L_ImageRGBDouble::copyG_into(L_ImageGrayDouble &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,1);
	return imC;
}
L_ImageGrayDouble& L_ImageRGBDouble::copyB_into(L_ImageGrayDouble &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,2);
	return imC;
}

L_ImageRGBDouble& L_ImageRGBDouble::composeFromChannels(const L_ImageGrayDouble& imR, const L_ImageGrayDouble& imG, const L_ImageGrayDouble& imB)
{
	int i, j;
	throw_L_ArgException_if(imR.lx!=imG.lx || imR.lx!=imB.lx || imR.ly!=imG.ly || imR.ly!=imB.ly, "L_ImageRGBDouble::composeFromChannels");
	reallocate(imR.lx, imR.ly);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=imR.pix(i,j);
			pix(i,j,1)=imG.pix(i,j);
			pix(i,j,2)=imB.pix(i,j);
		}
	}
	return *this;
}

L_ImageRGBDouble &L_ImageRGBDouble::RGB_to_YPbPr()
{
	double Y, Pb, Pr;
	double R,G,B;
	int i, j;
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			R=pix(i,j,0); G=pix(i,j,1); B=pix(i,j,2);
			Y  = +0.299   *R + 0.587   *G + 0.114   *B;
			Pb = -0.168736*R - 0.331264*G + 0.5     *B;
			Pr = +0.5     *R - 0.418688*G - 0.081312*B;
			pix(i,j,0)=(double)Y;
			pix(i,j,1)=(double)Pb;
			pix(i,j,2)=(double)Pr;
		}
	}
	return *this;
}

L_ImageRGBDouble &L_ImageRGBDouble::YPbPr_to_RGB()
{
	double Y, Pb, Pr;
	double R,G,B;
	int i, j;
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			Y=pix(i,j,0); Pb=pix(i,j,1); Pr=pix(i,j,2);
			R = Y + 1.402*Pr;
			G = Y -0.3441356781*Pb -0.7141361555*Pr;
			B = Y +1.772*Pb;
			pix(i,j,0)=(double)R;
			pix(i,j,1)=(double)G;
			pix(i,j,2)=(double)B;
		}
	}
	return *this;
}

L_ImageRGBDouble & L_ImageRGBDouble::operator=(const L_ImageGrayDouble &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.pix(i,j);
			pix(i,j,1)=other.pix(i,j);
			pix(i,j,2)=other.pix(i,j);
		}
	return *this;
}
L_ImageRGBDouble & L_ImageRGBDouble::operator=(const L_ImageGrayUchar &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.pix(i,j)/(double)255.0;
			pix(i,j,1)=pix(i,j,0);
			pix(i,j,2)=pix(i,j,1);
		}
	return *this;
}

L_ImageRGBDouble & L_ImageRGBDouble::operator=(const L_ImageRGBUchar &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.pix(i,j,0)/(double)255.0;
			pix(i,j,1)=other.pix(i,j,1)/(double)255.0;
			pix(i,j,2)=other.pix(i,j,2)/(double)255.0;
		}
	}
	return *this;
}

void L_ImageRGBUchar::fill(L_uchar R, L_uchar G, L_uchar B)
{
	int i, j;
	throw_L_ArgException_if(data()==NULL || lx<0 || ly<0, "L_ImageRGBUchar::fill() : imagen no valida");
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0) = R;
			pix(i,j,1) = G;
			pix(i,j,2) = B;
		}
	}
}


L_ImageRGBUchar& L_ImageRGBUchar::transpOf(const L_ImageRGBUchar &im)
{
	int i, j;
	reallocate(im.ly, im.lx);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=im.pix(j,i,0);
			pix(i,j,1)=im.pix(j,i,1);
			pix(i,j,2)=im.pix(j,i,2);
		}
	}
	return *this;
}

L_ImageRGBUchar& L_ImageRGBUchar::doubleResolution(const L_ImageRGBUchar &other)
{
	int i, j;
	int u, v;
	reallocate(2*other.lx, 2*other.ly);
	for (j=0; j<other.ly; j++)
		for (i=0; i<other.lx; i++)
		{
			u=2*i;
			v=2*j;
			pix(u,v,0)=pix(u+1,v,0)=pix(u,v+1,0)=pix(u+1,v+1,0)=other.pix(i,j,0);
			pix(u,v,1)=pix(u+1,v,1)=pix(u,v+1,1)=pix(u+1,v+1,1)=other.pix(i,j,1);
			pix(u,v,2)=pix(u+1,v,2)=pix(u,v+1,2)=pix(u+1,v+1,2)=other.pix(i,j,2);
		}
	return *this;
}

L_ImageRGBUchar& L_ImageRGBUchar::halfResolution(const L_ImageRGBUchar &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v,0)=(L_uchar)(0.25*(1.0*other.pix(i,j,0)+other.pix(i+1,j,0)+other.pix(i,j+1,0)+other.pix(i+1,j+1,0)));
			pix(u,v,1)=(L_uchar)(0.25*(1.0*other.pix(i,j,1)+other.pix(i+1,j,1)+other.pix(i,j+1,1)+other.pix(i+1,j+1,1)));
			pix(u,v,2)=(L_uchar)(0.25*(1.0*other.pix(i,j,2)+other.pix(i+1,j,2)+other.pix(i,j+1,2)+other.pix(i+1,j+1,2)));
		}
	return *this;
}

L_ImageRGBUchar& L_ImageRGBUchar::halfResolution_noResample(const L_ImageRGBUchar &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v,0)=other.pix(i,j,0);
			pix(u,v,1)=other.pix(i,j,1);
			pix(u,v,2)=other.pix(i,j,2);
		}
	return *this;
}

void L_ImageRGBUchar::applyHomography(const L_Matrix &h8, const L_ImageRGBUchar &orig)
{
	double xh, yh, zh;
	double xf, yf;
	double esqix[4], esqiy[4];
	double esqhx[4], esqhy[4];
	L_Matrix hi(8,1);
	L_Matrix m(3,3);
	double a, b;
	double u, v;
	int i, j, k;
	int tx=0, ty=0;

	esqix[0] = 0;
	esqiy[0] = 0;

	esqix[1] = orig.lx;
	esqiy[1] = 0;

	esqix[2] = orig.lx;
	esqiy[2] = orig.ly;

	esqix[3] = 0;
	esqiy[3] = orig.ly;

	for (i=0; i<4; i++)
	{
		xh = h8(0,0) * esqix[i] + h8(1,0) * esqiy[i] + h8(2,0);
		yh = h8(3,0) * esqix[i] + h8(4,0) * esqiy[i] + h8(5,0);
		zh = h8(6,0) * esqix[i] + h8(7,0) * esqiy[i] + 1.0;
		esqhx[i] = xh/zh;
		esqhy[i] = yh/zh;
		if (esqhx[i] > tx)
			tx = (int)(esqhx[i]+1);
		if (esqhy[i] > ty)
			ty = (int)(esqhy[i]+1);
	}

	reallocate(tx, ty);

	// Hay que aplicar la homografia inversa
	m(0,0) = h8(0,0);
	m(0,1) = h8(1,0);
	m(0,2) = h8(2,0);
	m(1,0) = h8(3,0);
	m(1,1) = h8(4,0);
	m(1,2) = h8(5,0);
	m(2,0) = h8(6,0);
	m(2,1) = h8(7,0);
	m(2,2) = 1.0;

	m.invertMe();

	hi(0,0) = m(0,0) / m(2,2);
	hi(1,0) = m(0,1) / m(2,2);
	hi(2,0) = m(0,2) / m(2,2);
	hi(3,0) = m(1,0) / m(2,2);
	hi(4,0) = m(1,1) / m(2,2);
	hi(5,0) = m(1,2) / m(2,2);
	hi(6,0) = m(2,0) / m(2,2);
	hi(7,0) = m(2,1) / m(2,2);


	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			xh = hi(0,0) * i + hi(1,0) * j + hi(2,0);
			yh = hi(3,0) * i + hi(4,0) * j + hi(5,0);
			zh = hi(6,0) * i + hi(7,0) * j + 1.0;
			xf = xh/zh;
			yf = yh/zh;
			
			if ((int)xf >= 0 && (int)yf >= 0 && (int)xf < orig.lx-1 && (int)yf < orig.ly-1)
			{
				a = xf - (int)xf;
				b = yf - (int)yf;
				for (k=0; k<3; k++)
				{
					u = orig.pix((int)xf,(int)yf,k)   * (1-a) + orig.pix((int)xf+1,(int)yf,k)   * a;
					v = orig.pix((int)xf,(int)yf+1,k) * (1-a) + orig.pix((int)xf+1,(int)yf+1,k) * a;
					pix(i,j,k) = (L_uchar)( u*(1-b) + v*b );
				}
			}
			else
			{
				pix(i,j,0) = 0;
				pix(i,j,1) = 0;
				pix(i,j,2) = 0;
			}
		}
	}
}


L_ImageGrayUchar& L_ImageRGBUchar::copyR_into(L_ImageGrayUchar &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,0);
	return imC;
}
L_ImageGrayUchar& L_ImageRGBUchar::copyG_into(L_ImageGrayUchar &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,1);
	return imC;
}
L_ImageGrayUchar& L_ImageRGBUchar::copyB_into(L_ImageGrayUchar &imC)
{
	int i, j;
	imC.reallocate(lx, ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
			imC.pix(i,j)=pix(i,j,2);
	return imC;
}

void L_ImageRGBUchar::concatenateHorz(const L_ImageRGBUchar &im1, const L_ImageRGBUchar &im2)
{
	int i, j;
	reallocate(im1.lx+im2.lx, L_MAX(im1.ly, im2.ly));
	for (j=0; j<im1.ly; j++)
	{
		for (i=0; i<im1.lx; i++)
		{
			pix(i,j,0) = im1.pix(i,j,0);
			pix(i,j,1) = im1.pix(i,j,1);
			pix(i,j,2) = im1.pix(i,j,2);
		}
	}
	for (; j<ly; j++) // Pintar negra la parte de abajo
	{
		for (i=0; i<im1.lx; i++)
		{
			pix(i,j,0) = 0;
			pix(i,j,1) = 0;
			pix(i,j,2) = 0;
		}
	}

	for (j=0; j<im2.ly; j++)
	{
		for (i=0; i<im2.lx; i++)
		{
			pix(im1.lx+i,j,0) = im2.pix(i,j,0);
			pix(im1.lx+i,j,1) = im2.pix(i,j,1);
			pix(im1.lx+i,j,2) = im2.pix(i,j,2);

		}
	}
	for (; j<ly; j++) // Pintar negra la parte de abajo
	{
		for (i=0; i<im2.lx; i++)
		{
			pix(im1.lx+i,j,0) = 0;
			pix(im1.lx+i,j,1) = 0;
			pix(im1.lx+i,j,2) = 0;
		}
	}
}

void L_ImageRGBUchar::divideHorz(L_ImageRGBUchar &im1, L_ImageRGBUchar &im2) const
{
	int i, j;
	im1.reallocate(lx/2, ly);
	im2.reallocate(lx - im1.lx, ly);
	for (j=0; j<im1.ly; j++)
	{
		for (i=0; i<im1.lx; i++)
		{
			im1.pix(i,j,0) = pix(i,j,0);
			im1.pix(i,j,1) = pix(i,j,1);
			im1.pix(i,j,2) = pix(i,j,2);
		}
	}
	for (j=0; j<im2.ly; j++)
	{
		for (i=0; i<im2.lx; i++)
		{
			im2.pix(i,j,0) = pix(im1.lx+i,j,0);
			im2.pix(i,j,1) = pix(im1.lx+i,j,1);
			im2.pix(i,j,2) = pix(im1.lx+i,j,2);
		}
	}
}


L_ImageRGBUchar& L_ImageRGBUchar::composeFromChannels(const L_ImageGrayUchar& imR, const L_ImageGrayUchar& imG, const L_ImageGrayUchar& imB)
{
	int i, j;
	throw_L_ArgException_if(imR.lx!=imG.lx || imR.lx!=imB.lx || imR.ly!=imG.ly || imR.ly!=imB.ly, "L_ImageRGBUchar::composeFromChannels");
	reallocate(imR.lx, imR.ly);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=imR.pix(i,j);
			pix(i,j,1)=imG.pix(i,j);
			pix(i,j,2)=imB.pix(i,j);
		}
	}
	return *this;
}




void L_ImageRGBUchar::encodeToBuffer(std::vector<char> &buf) const
{
	buf.resize(lx*ly*3);
	copyToBuffer((L_uchar *)&buf[0]);
}
void L_ImageRGBUchar::decodeFromBuffer(const std::vector<char> &buf)
{
	throw_L_ArgException_if(buf.size() != lx*ly*3, "L_ImageRGBUchar::decodeFromBuffer()");
	copyFromBuffer((L_uchar *)&buf[0]);
}
void L_ImageRGBUchar::encodeDelta(std::vector<char> &bufD) const
{
	int i, j;
	bufD.resize(3*lx*ly);
	bufD[0] = operator()(0,0,0);
	bufD[1] = operator()(0,0,1);
	bufD[2] = operator()(0,0,2);
	for (i=1; i<ly; i++)
	{
		bufD[3*i*lx+0] = operator()(i,0,0) - operator()(i-1,0,0);
		bufD[3*i*lx+1] = operator()(i,0,1) - operator()(i-1,0,1);
		bufD[3*i*lx+2] = operator()(i,0,2) - operator()(i-1,0,2);
	}
	for (i=0; i<ly; i++)
	{
		for (j=1; j<lx; j++)
		{
			bufD[3*i*lx+3*j+0] = operator()(i,j,0) - operator()(i,j-1,0);
			bufD[3*i*lx+3*j+1] = operator()(i,j,1) - operator()(i,j-1,1);
			bufD[3*i*lx+3*j+2] = operator()(i,j,2) - operator()(i,j-1,2);
		}
	}
}
void L_ImageRGBUchar::decodeDelta(std::vector<char> &bufD)
{
	int i, j;
	throw_L_ArgException_if(bufD.size() != lx*ly*3, "L_ImageRGBUchar::decodeDelta()");
	operator()(0,0,0) = bufD[0];
	operator()(0,0,1) = bufD[1];
	operator()(0,0,2) = bufD[2];
	for (i=1; i<ly; i++)
	{
		operator()(i,0,0) = operator()(i-1,0,0) + bufD[3*i*lx+0];
		operator()(i,0,1) = operator()(i-1,0,1) + bufD[3*i*lx+1];
		operator()(i,0,2) = operator()(i-1,0,2) + bufD[3*i*lx+2];
	}
	for (i=0; i<ly; i++)
	{
		for (j=1; j<lx; j++)
		{
			operator()(i,j,0) = operator()(i,j-1,0) + bufD[3*i*lx+3*j+0];
			operator()(i,j,1) = operator()(i,j-1,1) + bufD[3*i*lx+3*j+1];
			operator()(i,j,2) = operator()(i,j-1,2) + bufD[3*i*lx+3*j+2];
		}
	}
}


L_ImageRGBUchar & L_ImageRGBUchar::copyImageWithOffset(const L_ImageGrayDouble &im, int xc, int yc)
{
	int i, j;
	double val;
	L_uchar val2;
	for (j=0; j<im.ly; j++)
	{
		for (i=0; i<im.lx; i++)
		{
			val=im.pix(i,j)*255;
			if (val>255.0)
				val=255.0;
			else if (val<0)
				val=0;
			val2=(L_uchar)val;
			pix(i+xc,j+yc,0)=val2;
			pix(i+xc,j+yc,1)=val2;
			pix(i+xc,j+yc,2)=val2;
		}
	}
	return *this;
}

L_ImageRGBUchar & L_ImageRGBUchar::copyImageWithOffset(const L_ImageRGBUchar &im, int xc, int yc)
{
	int i, j;
	for (j=0; j<im.ly; j++)
	{
		for (i=0; i<im.lx; i++)
		{
			pix(i+xc,j+yc,0)=im.pix(i,j,0);
			pix(i+xc,j+yc,1)=im.pix(i,j,1);
			pix(i+xc,j+yc,2)=im.pix(i,j,2);
		}
	}
	return *this;
}


bool L_ImageRGBUchar::searchInside(L_ImageRGBUchar &imGrande, int &x, int &y)
{
	int i, j, u, v;
	bool valido, sale;
	for (j=0; j<imGrande.ly-ly; j++)
	{
		for (i=0; i<imGrande.lx-lx; i++)
		{
			valido = true;
			sale = false;
			for (v=0; v<ly; v++)
			{
				for (u=0; u<lx; u++)
				{
					if (pix(u,v,0) != imGrande.pix(i+u,j+v,0))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,1) != imGrande.pix(i+u,j+v,1))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,2) != imGrande.pix(i+u,j+v,2))
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (valido)
			{
				x = i;
				y = j;
				return true;
			}
		}

	}
	return false;
}


bool L_ImageRGBUchar::searchInside_delta(L_ImageRGBUchar &imGrande, int delta, int &x, int &y)
{
	int i, j, u, v;
	bool valido, sale;
	for (j=0; j<imGrande.ly-ly; j++)
	{
		for (i=0; i<imGrande.lx-lx; i++)
		{
			valido = true;
			sale = false;
			for (v=0; v<ly; v++)
			{
				for (u=0; u<lx; u++)
				{
					if (fabs((double)pix(u,v,0) - imGrande.pix(i+u,j+v,0)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (fabs((double)pix(u,v,1) - imGrande.pix(i+u,j+v,1)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (fabs((double)pix(u,v,2) - imGrande.pix(i+u,j+v,2)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (valido)
			{
				x = i;
				y = j;
				return true;
			}
		}

	}
	return false;
}

bool L_ImageRGBUchar::searchInside_transp(L_ImageRGBUchar &imGrande, L_uchar R, L_uchar G, L_uchar B, int x0, int y0, int x1, int y1, int &x, int &y)
{
	int i, j, u, v;
	bool valido, sale;
	for (j=0; j<imGrande.ly-ly; j++)
	{
		for (i=0; i<imGrande.lx-lx; i++)
		{
			valido = true;
			sale = false;
			for (v=y0; v<y1; v++)
			{
				for (u=x0; u<x1; u++)
				{
					if (pix(u,v,0) != R && pix(u,v,0) != imGrande.pix(i+u,j+v,0))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,1) != G && pix(u,v,1) != imGrande.pix(i+u,j+v,1))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,2) != B && pix(u,v,2) != imGrande.pix(i+u,j+v,2))
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (sale)
				continue;
			for (v=0; v<ly; v++)
			{
				for (u=0; u<lx; u++)
				{
					if (pix(u,v,0) != R && pix(u,v,0) != imGrande.pix(i+u,j+v,0))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,1) != G && pix(u,v,1) != imGrande.pix(i+u,j+v,1))
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,2) != B && pix(u,v,2) != imGrande.pix(i+u,j+v,2))
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (valido)
			{
				x = i;
				y = j;
				return true;
			}
		}

	}
	return false;
}

bool L_ImageRGBUchar::searchInside_transp_delta(L_ImageRGBUchar &imGrande, int delta, L_uchar R, L_uchar G, L_uchar B, int x0, int y0, int x1, int y1, int &x, int &y)
{
	int i, j, u, v;
	bool valido, sale;
	for (j=0; j<imGrande.ly-ly; j++)
	{
		for (i=0; i<imGrande.lx-lx; i++)
		{
			valido = true;
			sale = false;
			for (v=y0; v<y1; v++)
			{
				for (u=x0; u<x1; u++)
				{
					if (pix(u,v,0) != R && fabs(pix(u,v,0) - (double)imGrande.pix(i+u,j+v,0)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,1) != G && fabs(pix(u,v,1) - (double)imGrande.pix(i+u,j+v,1)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,2) != B && fabs(pix(u,v,2) - (double)imGrande.pix(i+u,j+v,2)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (sale)
				continue;
			for (v=0; v<ly; v++)
			{
				for (u=0; u<lx; u++)
				{
					if (pix(u,v,0) != R && fabs(pix(u,v,0) - (double)imGrande.pix(i+u,j+v,0)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,1) != G && fabs(pix(u,v,1) - (double)imGrande.pix(i+u,j+v,1)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
					if (pix(u,v,2) != B && fabs(pix(u,v,2) - (double)imGrande.pix(i+u,j+v,2)) > delta)
					{
						valido = false;
						sale = true;
						break;
					}
				}
				if (sale)
					break;
			}
			if (valido)
			{
				x = i;
				y = j;
				return true;
			}
		}

	}
	return false;
}

void L_ImageRGBUchar::computeMeanShiftHistogram(double xCen, double yCen, std::vector<double> &hist, double radio_h_x, double radio_h_y, int nBinsPorCanal)
{
	int i, j, i0, j0, i1, j1;
	int ind, anchoBin = 256 / nBinsPorCanal;
	double dist2, kernelVal, deni, den=0;
	double radX2i, radY2i;
	hist.resize(nBinsPorCanal * nBinsPorCanal * nBinsPorCanal);

	for (i=0; i<(int)hist.size(); i++)
		hist[i] = 0;

	i0 = (int)(xCen - radio_h_x);
	j0 = (int)(yCen - radio_h_y);
	i1 = (int)(xCen + radio_h_x);
	j1 = (int)(yCen + radio_h_y);
	
	if (i0 < 0)
		i0 = 0;
	if (j0 < 0)
		j0 = 0;
	if (i1 > lx-1)
		i1 = lx-1;
	if (j1 > ly-1)
		j1 = ly-1;

	radX2i = 1/(radio_h_x*radio_h_x);
	radY2i = 1/(radio_h_y*radio_h_y);

	for (j=j0; j<=j1; j++)
	{
		for (i=i0; i<=i1; i++)
		{
			ind = ((pix(i,j,0) / anchoBin)*nBinsPorCanal + (pix(i,j,1) / anchoBin))*nBinsPorCanal + (pix(i,j,2) / anchoBin);
			dist2 = (i-xCen)*(i-xCen)*radX2i + (j-yCen)*(j-yCen)*radY2i;
			kernelVal = (dist2 < 1) ? 1-dist2 : 0; // Opcion 1
			// Calculo de los histogramas = integral de la distribucion de probabilidad para cada color elem[ind]
			hist[ind] += kernelVal * dist2;
			den += kernelVal * dist2;
		}
	}
	deni = 1/den; // Normalizar el histograma
	for (i=0; i<(int)hist.size(); i++)
		hist[i] *= deni;
}

void L_ImageRGBUchar::computeMeanShiftTracking_private(double &xCen, double &yCen, std::vector<double> &histAntiguo, std::vector<double> &histNuevo, double radio_h_x, double radio_h_y, int nBinsPorCanal)
{

	int i, j, i0, j0, i1, j1;
	int ind, anchoBin = 256 / nBinsPorCanal;
	double dist2, kernelVal, w, fac;
	double radX2i, radY2i;
	double xFin=0, yFin=0, den = 0;

	i0 = (int)(xCen - radio_h_x);
	j0 = (int)(yCen - radio_h_y);
	i1 = (int)(xCen + radio_h_x);
	j1 = (int)(yCen + radio_h_y);
	
	if (i0 < 0)
		i0 = 0;
	if (j0 < 0)
		j0 = 0;
	if (i1 > lx-1)
		i1 = lx-1;
	if (j1 > ly-1)
		j1 = ly-1;

	radX2i = 1/(radio_h_x*radio_h_x);
	radY2i = 1/(radio_h_y*radio_h_y);

	for (j=j0; j<=j1; j++)
	{
		for (i=i0; i<=i1; i++)
		{
			ind = ((pix(i,j,0) / anchoBin)*nBinsPorCanal + (pix(i,j,1) / anchoBin))*nBinsPorCanal + (pix(i,j,2) / anchoBin);
			dist2 = (i-xCen)*(i-xCen)*radX2i + (j-yCen)*(j-yCen)*radY2i;
			kernelVal = (dist2 < 1) ? 1-dist2 : 0; // Opcion 1
			// Estimacion del movimiento de la distribucion de color elem[ind] sobre la imagen
			if (kernelVal > 0 && histNuevo[ind] != 0) // Esto ultimo NO DEBERIA PASAR, pero a veces pasa...
			{
				w = sqrt(histAntiguo[ind] / histNuevo[ind]);
				fac = w * kernelVal * dist2;
			}
			else
				fac = 0;
			xFin += (i-xCen) * fac;
			yFin += (j-yCen) * fac;
			den += fac;
		}
	}
	xFin /= den;
	yFin /= den;
	xCen += xFin;
	yCen += yFin;
}


L_ImageRGBUchar & L_ImageRGBUchar::operator=(const L_ImageGrayDouble &other)
{
	int i, j;
	double prom;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			prom=other.pix(i,j);
			if (prom<0) prom=0;
			if (prom>1) prom=1;
			pix(i,j,0)=(L_uchar)(prom*255);
			pix(i,j,1)=(L_uchar)(prom*255);
			pix(i,j,2)=(L_uchar)(prom*255);
		}
	return *this;
}
L_ImageRGBUchar & L_ImageRGBUchar::operator=(const L_ImageGrayUchar &other)
{
	int i, j;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.pix(i,j);
			pix(i,j,1)=other.pix(i,j);
			pix(i,j,2)=other.pix(i,j);
		}
	return *this;
}
L_ImageRGBUchar & L_ImageRGBUchar::operator=(const L_ImageRGBDouble &other)
{
	int i, j, c;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			for (c=0; c<3; c++)
			{
				if (other.pix(i,j,c)<0)
					pix(i,j,c)=0;
				else if (other.pix(i,j,c)>=1)
					pix(i,j,c)=255;
				else
					pix(i,j,c)=(L_uchar)(other.pix(i,j,c)*255);
			}
		}
	}
	return *this;
}

void L_ImageRGBUchar::operator=(const L_ImageBGRAuchar_row &other)
{
	int i, j, k=0;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,2)=other.buf[k++]; //B
			pix(i,j,1)=other.buf[k++]; //G
			pix(i,j,0)=other.buf[k++]; //R
			k++; //A
		}
}

void L_ImageRGBUchar::operator=(const L_ImageRGBAuchar_row &other)
{
	int i, j, k=0;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.buf[k++]; //B
			pix(i,j,1)=other.buf[k++]; //G
			pix(i,j,2)=other.buf[k++]; //R
			k++; //A
		}
}

void L_ImageRGBUchar::operator=(const L_ImagenRGBuchar_row &other)
{
	int i, j, k=0;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=other.buf[k++]; //R
			pix(i,j,1)=other.buf[k++]; //G
			pix(i,j,2)=other.buf[k++]; //B
		}
}

#ifdef __COMPAT_ATLIMAGE__
void L_ImageRGBUchar::operator=(const ATL_CImage &other)
{
	int i, j;
	const L_uchar *ptr;
	int nb;
	reallocate(other.GetWidth(), other.GetHeight());
	nb=other.GetBPP();
	if (nb==8)
	{
		for (j=0; j<ly; j++)
		{
			ptr=(const L_uchar *)other.GetPixelAddress(0,j);
			for (i=0; i<lx; i++)
			{
				pix(i,j,0)=ptr[i];
				pix(i,j,1)=ptr[i];
				pix(i,j,2)=ptr[i];
			}
		}
	}
	else if (nb==24)
	{
		for (j=0; j<ly; j++)
		{
			ptr=(const L_uchar *)other.GetPixelAddress(0,j);
			for (i=0; i<lx; i++)
			{
				// R: x&0xFF; G: (x>>8)&0xFF; B:(x>>16)&0xFF
				// Ver funciones GetRValue, ...
				pix(i,j,0)=ptr[i*3+ 2];
				pix(i,j,1)=ptr[i*3+ 1];
				pix(i,j,2)=ptr[i*3+ 0];
			}
		}
	}
	else if (nb==32)
	{
		for (j=0; j<ly; j++)
		{
			ptr=(const L_uchar *)other.GetPixelAddress(0,j);
			for (i=0; i<lx; i++)
			{
				pix(i,j,0)=ptr[i*4+ 2];
				pix(i,j,1)=ptr[i*4+ 1];
				pix(i,j,2)=ptr[i*4+ 0];
			}
		}
	}
}
#endif // __COMPAT_ATLIMAGE__

#ifdef __COMPAT_ATLIMAGE__
void L_ImageRGBUchar::copyTo(ATL_CImage &other) const
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	int i, j;
	L_uchar *ptr;
	if (!other.IsNull())
		other.Destroy();
	other.Create(lx, ly, 24);
	if (other.GetBPP()==24)
	{
		for (j=0; j<ly; j++)
		{
			ptr=(L_uchar *)other.GetPixelAddress(0,j);
			for (i=0; i<lx; i++)
			{
				ptr[i*3+ 2]=pix(i,j,0);
				ptr[i*3+ 1]=pix(i,j,1);
				ptr[i*3+ 0]=pix(i,j,2);
			}
		}
	}
	else if (other.GetBPP()==32)
	{
		for (j=0; j<ly; j++)
		{
			ptr=(L_uchar *)other.GetPixelAddress(0,j);
			for (i=0; i<lx; i++)
			{
				ptr[i*4+ 2]=pix(i,j,0);
				ptr[i*4+ 1]=pix(i,j,1);
				ptr[i*4+ 0]=pix(i,j,2);
			}
		}
	}
}
#endif // __COMPAT_ATLIMAGE__


#ifdef __COMPAT_IPLIMAGE__
void L_ImageRGBUchar::operator=(const IplImage &imagen1)
{
	unsigned char *pixel;
	int i, j;
	int origin = imagen1.origin;

#if defined(i_)
	#error #define conflicts
#endif
	#define i_ (ly-i-1)
	reallocate(imagen1.width,imagen1.height);
	if(imagen1.nChannels == 3)
	{   // RGB
		if (origin==IPL_ORIGIN_BL)
		{
			for (i=0; i<imagen1.height; i++)
			{
				for (j=0; j<imagen1.width; j++)
				{
					pixel   = &((unsigned char*)(imagen1.imageData + imagen1.widthStep*i))[j*3];
					pix(j,i_,0)=(L_uchar)pixel[2]; // 0
					pix(j,i_,1)=(L_uchar)pixel[1]; // 1
					pix(j,i_,2)=(L_uchar)pixel[0]; // 2
				}
			}
		}
		else
		{
			for (i=0; i<imagen1.height; i++)
			{
				for (j=0; j<imagen1.width; j++)
				{
					pixel   = &((unsigned char*)(imagen1.imageData + imagen1.widthStep*i))[j*3];
					pix(j,i,0)=(L_uchar)pixel[2]; // 0
					pix(j,i,1)=(L_uchar)pixel[1]; // 1
					pix(j,i,2)=(L_uchar)pixel[0]; // 2
				}
			}
		}
		return;
	}
	else if(imagen1.nChannels == 1)
	{  // gris
		if (origin==IPL_ORIGIN_BL)
		{
			for (i=0; i<imagen1.height; i++)
			{
				for (j=0; j<imagen1.width; j++)
				{
					pixel   = &((unsigned char*)(imagen1.imageData + imagen1.widthStep*i))[j];
					pix(j,i_,0)=(L_uchar)pixel[0];
					pix(j,i_,1)=(L_uchar)pixel[0];
					pix(j,i_,2)=(L_uchar)pixel[0];
				}
			}
		}
		else
		{
			for (i=0; i<imagen1.height; i++)
			{
				for (j=0; j<imagen1.width; j++)
				{
					pixel   = &((unsigned char*)(imagen1.imageData + imagen1.widthStep*i))[j];
					pix(j,i,0)=(L_uchar)pixel[0];
					pix(j,i,1)=(L_uchar)pixel[0];
					pix(j,i,2)=(L_uchar)pixel[0];
				}
			}
		}
		return;
	}
	#undef i_
	return;
}



void L_ImageRGBUchar::copyTo(IplImage **imagen2) const
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	unsigned char *pixel;
	int i, j;

	if ( (*imagen2) == NULL || (*imagen2)->width!=lx || (*imagen2)->height!=ly )
		(*imagen2) = cvCreateImage(cvSize(lx,ly),IPL_DEPTH_8U,3);

	int origin = (*imagen2)->origin;

	#if defined(i_)
		#error #define conflicts
	#endif
	#define i_ (ly-i-1)
	if((*imagen2)->nChannels == 3) // RGB
	{
		if (origin==IPL_ORIGIN_BL) // Eje y "invertido" respecto a L_ImageGrayDouble
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j*3];
					pixel[0]=(unsigned char)pix(j,i_,2); // 0
					pixel[1]=(unsigned char)pix(j,i_,1); // 1
					pixel[2]=(unsigned char)pix(j,i_,0); // 2
				}
			}
		}
		else
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j*3];
					pixel[0]=(unsigned char)pix(j,i,2); // 0
					pixel[1]=(unsigned char)pix(j,i,1); // 1
					pixel[2]=(unsigned char)pix(j,i,0); // 2
				}
			}
		}
		return;
	}
	else if((*imagen2)->nChannels == 1) // gris
	{
		if (origin==IPL_ORIGIN_BL)
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j];
					pixel[0]=(unsigned char)( ((double)pix(j,i_,0)+pix(j,i_,1)+(double)pix(j,i_,2))/3 );
			
				}
			}
		}
		else
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j];
					pixel[0]=(unsigned char)( ((double)pix(j,i,0)+pix(j,i,1)+(double)pix(j,i,2))/3 );
				}
			}
		}
		return;
	}
	#undef i_
	return;
}
#endif // __COMPAT_IPLIMAGE__


#ifdef __COMPAT_SVS__
bool L_ImageRGBUchar::SVSPideIzqBN(const svsStereoImage *imageObject)
{
	//int i, j;
	const L_uchar *bufBN;
	if (imageObject==NULL || imageObject->haveImages==false || imageObject->left == NULL)
		return false;
	reallocate(imageObject->sp.width, imageObject->sp.height);
	bufBN = imageObject->left;
	copyFromBuffer(bufBN);
	return true;
}
bool L_ImageRGBUchar::SVSPideDerBN(const svsStereoImage *imageObject)
{
	//int i, j;
	const L_uchar *bufBN;
	if (imageObject==NULL || imageObject->haveImages==false || imageObject->right == NULL)
		return false;
	reallocate(imageObject->sp.width, imageObject->sp.height);
	bufBN = imageObject->right;
	copyFromBuffer(bufBN);
	return true;
}
bool L_ImageRGBUchar::SVSPideIzqColor(const svsStereoImage *imageObject)
{
	if (imageObject==NULL || imageObject->haveColor==false || imageObject->color == NULL)
		return false;
	reallocate(imageObject->sp.width, imageObject->sp.height);
	L_uchar *buf = imageObject->color;
	L_uchar *p = &pix(0,0,0);
	L_uchar *pf = &pix(0,0,0) + lxStep*ly;
	while (p < pf)
	{
		p[0]=buf[2];
		p[1]=buf[1];
		p[2]=buf[0];
		buf+=4;
		p+=3;
	}
	return true;
}
bool L_ImageRGBUchar::SVSPideDerColor(const svsStereoImage *imageObject)
{
	if (imageObject==NULL || imageObject->haveColorRight==false || imageObject->color_right == NULL)
		return false;
	reallocate(imageObject->sp.width, imageObject->sp.height);
	L_uchar *buf = imageObject->color_right;
	L_uchar *p = &pix(0,0,0);
	L_uchar *pf = &pix(0,0,0) + lxStep*ly;
	while (p < pf)
	{
		p[0]=buf[2];
		p[1]=buf[1];
		p[2]=buf[0];
		buf+=4;
		p+=3;
	}
	return true;
}
#endif //__COMPAT_SVS__



void L_Matrix_intercambia::reallocate(int li0, int lj0)
{
	double **elem_;
	if (elem==NULL)
	{
		elem=L_new2dM<double>(li0, lj0); // Es ineficiente esto...
		li=li0;
		lj=lj0;
	}
	else if (li!=li0 || lj!=lj0)
	{
		elem_=L_new2dM<double>(li0, lj0);
		destroy();
		elem=elem_;
		li=li0;
		lj=lj0;
	}
}


void L_Matrix_intercambia::destroy()
{
	if (elem!=NULL)
	{
		L_delete2dM<double>(elem, li); // Es ineficiente esto...
		elem = NULL;
	}
	return;
}

void L_Matrix_intercambia::OP_assign(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.begin()==NULL, "L_Matrix_intercambia::OP_assign"); // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=other(i,j);
	return;
}

void L_Matrix_intercambia::copyTo(L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(elem==NULL, "L_Matrix_intercambia::copyTo"); // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	other.reallocate(li, lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			other(i,j)=operator()(i,j);
	return;
}


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
bool L_Matrix_intercambia::gaussj(L_Matrix_intercambia &B)
{
	double **a = elem;
	int n = li;
	double **b = B.elem;
	int m = B.lj;
	std::vector<int> indxc(n);
	std::vector<int> indxr(n);
	std::vector<int> ipiv(n);
	int i,icol=0,irow=0,j,k,l,ll;
	double big,dum,pivinv,temp;

	for (j=0;j<n;j++)
		ipiv[j]=0;

	for (i=0;i<n;i++)
	{
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++)
				{
					if (ipiv[k] == 0)
					{
						if (fabs(a[j][k]) >= big)
						{
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
					else if (ipiv[k] > 1)
						return false;
				}
		++(ipiv[icol]);
		if (irow != icol)
		{
			for (l=0;l<n;l++)
				SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++)
				SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0)
			return false;
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++)
			a[icol][l] *= pivinv;
		for (l=0;l<m;l++)
			b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol)
			{
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++)
					a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++)
					b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>1;l--)
	{
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	return true;
}
#undef SWAP

void L_MatrixFloat::OP_mult(const L_MatrixFloat &m1, const L_MatrixFloat &m2)
{
#ifdef __COMPAT_UBLAS__
	return OP_mult_ublas(m1, m2);
#else
	return OP_mult_private(m1, m2);
#endif
}


void L_MatrixFloat::OP_mult_private(const L_MatrixFloat &m1, const L_MatrixFloat &m2)
{
	int i,j,u;

	throw_L_ArgException_if(m1.lj!=m2.li || m1.begin() == begin() || m2.begin()==begin(), "L_Matrix::OP_mult()");

	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=0; j<m2.lj; j++)
		{
			operator()(i,j)=0;
			for (u=0; u<m1.lj; u++)
				operator()(i,j)+=m1(i,u)*m2(u,j);
		}
	}
	return;
}


bool L_MatrixFloat::invertMeDefPos()
{
    if (li == 0 || lj == 0)
		return false;
    if (li*lj == 1)
	{
		operator()(0,0) = 1/operator()(0,0);
		return true;
	}
	float *data = begin();

	for (int i=1; i < li; i++)
		data[i] /= data[0]; // normalize row 0
	for (int i=1; i < li; i++)
	{ 
		for (int j=i; j < li; j++)
		{ // do a column of L
			float sum = 0.0;
			for (int k = 0; k < i; k++)  
				sum += data[j*lj+k] * data[k*lj+i];
			data[j*lj+i] -= sum;
		}
		if (i == li-1)
			continue;
		for (int j=i+1; j < li; j++)
		{  // do a row of U
			float sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += data[i*lj+k]*data[k*lj+j];
			data[i*lj+j] = 
			(data[i*lj+j]-sum) / data[i*lj+i];
		}
	}
	for ( int i = 0; i < li; i++ )  // invert L
		for ( int j = i; j < li; j++ )
		{
			float x = 1.0;
			if ( i != j )
			{
				x = 0.0;
				for ( int k = i; k < j; k++ ) 
				x -= data[j*lj+k]*data[k*lj+i];
			}
			data[j*lj+i] = x / data[j*lj+j];
		}
	for ( int i = 0; i < li; i++ )   // invert U
		for ( int j = i; j < li; j++ )
		{
			if ( i == j )
				continue;
			float sum = 0.0;
			for ( int k = i; k < j; k++ )
				sum += data[k*lj+j]*( (i==k) ? 1.0f : data[i*lj+k] );
			data[i*lj+j] = -sum;
		}
	for ( int i = 0; i < li; i++ )   // final inversion
		for ( int j = 0; j < li; j++ )
		{
			float sum = 0.0;
			for ( int k = ((i>j)?i:j); k < li; k++ )  
				sum += ((j==k)?1.0f:data[j*lj+k])*data[k*lj+i];
			data[j*lj+i] = sum;
		}
	return true;
};


bool L_MatrixFloat::invertMeDefPosx4()
{
    if (li == 0 || lj == 0)
		return false;
    if (li*lj == 1)
	{
		operator()(0,0) = 1/operator()(0,0);
		return true;
	}
	int k;
	float *data = begin();

	for (int i=1; i < li; i++)
		data[i] /= data[0]; // normalize row 0
	for (int i=1; i < li; i++)
	{ 
		for (int j=i; j < li; j++)
		{ // do a column of L
			float sum[4] = {0.0f,0.0f,0.0f,0.0f};
			for ( k = 0; k < i; k+=4)
			{
				sum[0] += data[j*lj+k+0] * data[k*lj+i+0];
				sum[1] += data[j*lj+k+1] * data[k*lj+i+1];
				sum[2] += data[j*lj+k+2] * data[k*lj+i+2];
				sum[3] += data[j*lj+k+3] * data[k*lj+i+3];
			}
			k-=4; // Faltaba esto
			switch (i%4)
			{
			case 3: sum[0] += data[j*lj+k+0] * data[k*lj+i+0]; k++;
			case 2: sum[0] += data[j*lj+k+0] * data[k*lj+i+0]; k++;
			case 1: sum[0] += data[j*lj+k+0] * data[k*lj+i+0]; k++;
			}
			data[j*lj+i] -= sum[0]+sum[1]+sum[2]+sum[3];
		}
		if (i == li-1)
			continue;
		for (int j=i+1; j < li; j++)
		{  // do a row of U
			float sum[4] = {0.0f,0.0f,0.0f,0.0f};
			for ( k = 0; k < i; k+=4)
			{
				sum[0] += data[i*lj+k+0]*data[k*lj+j+0];
				sum[1] += data[i*lj+k+1]*data[k*lj+j+1];
				sum[2] += data[i*lj+k+2]*data[k*lj+j+2];
				sum[3] += data[i*lj+k+3]*data[k*lj+j+3];
			}
			k-=4; // Faltaba esto
			switch (k%4)
			{
			case 3: sum[0] += data[i*lj+k+0]*data[k*lj+j+0]; k++;
			case 2: sum[0] += data[i*lj+k+0]*data[k*lj+j+0]; k++;
			case 1: sum[0] += data[i*lj+k+0]*data[k*lj+j+0]; k++;
			}
			data[i*lj+j] = (data[i*lj+j]-(sum[0]+sum[1]+sum[2]+sum[3])) / data[i*lj+i];
		}
	}
	for ( int i = 0; i < li; i++ )  // invert L
		for ( int j = i; j < li; j++ )
		{
			float x[4] = {1.0f,1.0f,1.0f,1.0f};
			if ( i != j )
			{
				x[0] = 0.0; x[1]=0.0; x[2]=0.0; x[3]=0.0;
				for ( k = i; k < j; k+=4 ) 
				{
					x[0] -= data[j*lj+k+0]*data[k*lj+i+0];
					x[1] -= data[j*lj+k+1]*data[k*lj+i+1];
					x[2] -= data[j*lj+k+2]*data[k*lj+i+2];
					x[3] -= data[j*lj+k+3]*data[k*lj+i+3];
				}
			}
			k-=4; // Faltaba esto
			switch (k%4)
			{
			case 3: x[0] -= data[j*lj+k+0]*data[k*lj+i+0]; k++;
			case 2: x[0] -= data[j*lj+k+0]*data[k*lj+i+0]; k++;
			case 1: x[0] -= data[j*lj+k+0]*data[k*lj+i+0]; k++;
			}
			data[j*lj+i] = (x[0]+x[1]+x[2]+x[3]) / data[j*lj+j];
		}
	for ( int i = 0; i < li; i++ )   // invert U
		for ( int j = i; j < li; j++ )
		{
			if ( i == j )
				continue;
			float sum[4] = {0.0f,0.0f,0.0f,0.0f};
			for ( k = i; k < lj; k+=4 )
			{
				sum[0] += data[k*lj+j+0]*( (i==k+0) ? 1.0f : data[i*lj+k+0] );
				sum[1] += data[k*lj+j+1]*( (i==k+1) ? 1.0f : data[i*lj+k+1] );
				sum[2] += data[k*lj+j+2]*( (i==k+2) ? 1.0f : data[i*lj+k+2] );
				sum[3] += data[k*lj+j+3]*( (i==k+3) ? 1.0f : data[i*lj+k+3] );
			}
			k-=4; // Faltaba esto
			switch (k%4)
			{
			case 3: sum[0] += data[k*lj+j]*( (i==k) ? 1.0f : data[i*lj+k] ); k++;
			case 2: sum[0] += data[k*lj+j]*( (i==k) ? 1.0f : data[i*lj+k] ); k++;
			case 1: sum[0] += data[k*lj+j]*( (i==k) ? 1.0f : data[i*lj+k] ); k++;
			}
			data[i*lj+j] = -(sum[0]+sum[1]+sum[2]+sum[3]);
		}
	for ( int i = 0; i < li; i++ )   // final inversion
		for ( int j = 0; j < li; j++ )
		{
			float sum[4] = {0.0f,0.0f,0.0f,0.0f};
			for ( k = ((i>j)?i:j); k < li; k+=4 )
			{
				sum[0] += ((j==k+0)?1.0f:data[j*lj+k+0])*data[k*lj+i+0];
				sum[1] += ((j==k+1)?1.0f:data[j*lj+k+1])*data[k*lj+i+1];
				sum[2] += ((j==k+2)?1.0f:data[j*lj+k+2])*data[k*lj+i+2];
				sum[3] += ((j==k+3)?1.0f:data[j*lj+k+3])*data[k*lj+i+3];
			}
			k-=4; // Faltaba esto
			switch (k%4)
			{
				case 3: sum[0] += ((j==k)?1.0f:data[j*lj+k+0])*data[k*lj+i+0]; k++;
				case 2: sum[0] += ((j==k)?1.0f:data[j*lj+k+0])*data[k*lj+i+0]; k++;
				case 1: sum[0] += ((j==k)?1.0f:data[j*lj+k+0])*data[k*lj+i+0]; k++;
			}
			data[j*lj+i] = sum[0]+sum[1]+sum[2]+sum[3];
		}
	return true;
}



#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<float> L_MatrixFloat::crearUblas() const
{
	size_t lI = (size_t)li, lJ = (size_t)lj;
	boost::numeric::ublas::matrix<float> mb(lI,lJ);
	memcpy(&(mb.data()[0]), data(), sizeof(float)*li*lj);
	return mb;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<float> L_MatrixFloat::crearUblasRef() const
{
	size_t lI = (size_t)li, lJ = (size_t)lj;
	return boost::numeric::ublas::L_make_matrix_from_pointer(lI, lJ, data());
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
L_MatrixFloat &L_MatrixFloat::operator=(const boost::numeric::ublas::matrix<float> &other)
{
	reallocate((int)other.size1(), (int)other.size2());
	memcpy(data(), &(other.data()[0]), sizeof(float)*li*lj);
	return *this;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
void L_MatrixFloat::OP_mult_ublas(const L_MatrixFloat &m1, const L_MatrixFloat &m2)
{
	boost::numeric::ublas::matrix<float> ma = m1.crearUblasRef(), mb = m2.crearUblasRef(), mc;
	mc = boost::numeric::ublas::prod(ma, mb);
	*this = mc;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
bool L_MatrixFloat::invertMe_ublas()
{
	boost::numeric::ublas::matrix<float> input = crearUblasRef(), output(li, lj);
	bool ret;
	ret = boost::numeric::ublas::L_InvertMatrix(input, output);
	if (ret == false)
		return false;
	*this = output;
	return ret;
}
#endif // __COMPAT_UBLAS__


#ifdef __COMPAT_EIGEN__
void L_Matrix::OP_mult_eigen(const L_Matrix &m1, const L_Matrix &m2)
{
	Eigen::MatrixXd mat1 = m1.crearEigenRef(), mat2 = m2.crearEigenRef();
	*this = mat1*mat2;
}
#endif // __COMPAT_EIGEN__


#ifdef __COMPAT_EIGEN__
bool L_Matrix::invertMe_eigen()
{
	Eigen::MatrixXd mat = crearEigenRef();
#if EIGEN_VERSION_AT_LEAST(2,90,0)
	Eigen::FullPivLU<Eigen::MatrixXd> lu(mat);
#else
	Eigen::LU<Eigen::MatrixXd> lu(mat);
#endif

	if(lu.isInvertible() == false)
		return false;

	*this = lu.inverse(); // No avisa si no pudo invertirla, que mala onda
	return true;
}
#endif // __COMPAT_EIGEN__


#ifdef __COMPAT_EIGEN__
void L_Matrix::vectPr_Eigen_sim(L_Matrix &valPr)
{
	Eigen::MatrixXd M = crearEigenRef();
#if EIGEN_VERSION_AT_LEAST(2,90,0)
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> sol(M, Eigen::ComputeEigenvectors); // Para matrices Hermitian
#else
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> sol(M, true); // Para matrices Hermitian
#endif
	valPr = sol.eigenvalues().template cast<double>();  // QUE COSA MAS IMPRESIONANTE ESTO
	*this = sol.eigenvectors();
}
#endif // __COMPAT_EIGEN__


void L_MatrizLongDouble::OP_mult(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2)
{
#ifdef __COMPAT_UBLAS__
	return OP_mult_ublas(m1, m2);
#else
	return OP_mult_private(m1, m2);
#endif
}


void L_MatrizLongDouble::OP_mult_private(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2)
{
	int i,j,u;
	double r;

	throw_L_ArgException_if(m1.lj!=m2.li || m1.begin() == begin() || m2.begin()==begin(), "L_Matrix::OP_mult()");

	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=0; j<m2.lj; j++)
		{
			r=0;
			for (u=0; u<m1.lj; u++)
				r+=m1(i,u)*m2(u,j);
			operator()(i,j) = r;
		}
	}
	return;
}

void L_MatrizLongDouble::OP_mult_ABT(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2)
{
	int i,j,u;
	double r;

	throw_L_ArgException_if(m1.lj!=m2.li || m1.begin() == begin() || m2.begin() == begin(), "L_MatrizLongDouble::OP_mult()");

	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=0; j<m2.li; j++)
		{
			r=0;
			for (u=0; u<m1.lj; u++)
				r+=m1(i,u)*m2(j,u);
			operator()(i,j) = r;
		}
	}
	return;
}


bool L_MatrizLongDouble::invertMeDefPos()
{
    if (li == 0 || lj == 0)
		return false;
    if (li*lj == 1)
	{
		operator()(0,0) = 1/operator()(0,0);
		return true;
	}
	long double *data = begin();

	for (int i=1; i < li; i++)
		data[i] /= data[0]; // normalize row 0
	for (int i=1; i < li; i++)
	{ 
		for (int j=i; j < li; j++)
		{ // do a column of L
			long double sum = 0.0;
			for (int k = 0; k < i; k++)  
				sum += data[j*lj+k] * data[k*lj+i];
			data[j*lj+i] -= sum;
		}
		if (i == li-1)
			continue;
		for (int j=i+1; j < li; j++)
		{  // do a row of U
			long double sum = 0.0;
			for (int k = 0; k < i; k++)
				sum += data[i*lj+k]*data[k*lj+j];
			data[i*lj+j] = 
			(data[i*lj+j]-sum) / data[i*lj+i];
		}
	}
	for ( int i = 0; i < li; i++ )  // invert L
		for ( int j = i; j < li; j++ )
		{
			long double x = 1.0;
			if ( i != j )
			{
				x = 0.0;
				for ( int k = i; k < j; k++ ) 
				x -= data[j*lj+k]*data[k*lj+i];
			}
			data[j*lj+i] = x / data[j*lj+j];
		}
	for ( int i = 0; i < li; i++ )   // invert U
		for ( int j = i; j < li; j++ )
		{
			if ( i == j )
				continue;
			long double sum = 0.0;
			for ( int k = i; k < j; k++ )
				sum += data[k*lj+j]*( (i==k) ? 1.0f : data[i*lj+k] );
			data[i*lj+j] = -sum;
		}
	for ( int i = 0; i < li; i++ )   // final inversion
		for ( int j = 0; j < li; j++ )
		{
			long double sum = 0.0;
			for ( int k = ((i>j)?i:j); k < li; k++ )  
				sum += ((j==k)?1.0f:data[j*lj+k])*data[k*lj+i];
			data[j*lj+i] = sum;
		}
	return true;
};



#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<long double> L_MatrizLongDouble::crearUblas() const
{
	size_t lI = (size_t)li, lJ = (size_t)lj;
	boost::numeric::ublas::matrix<long double> mb(lI,lJ);
	memcpy(&(mb.data()[0]), data(), sizeof(long double)*li*lj);
	return mb;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<long double> L_MatrizLongDouble::crearUblasRef() const
{
	size_t lI = (size_t)li, lJ = (size_t)lj;
	return boost::numeric::ublas::L_make_matrix_from_pointer(lI, lJ, data());
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
L_MatrizLongDouble &L_MatrizLongDouble::operator=(const boost::numeric::ublas::matrix<long double> &other)
{
	reallocate((int)other.size1(), (int)other.size2());
	memcpy(data(), &(other.data()[0]), sizeof(long double)*li*lj);
	return *this;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
void L_MatrizLongDouble::OP_mult_ublas(const L_MatrizLongDouble &m1, const L_MatrizLongDouble &m2)
{
	boost::numeric::ublas::matrix<long double> ma = m1.crearUblasRef(), mb = m2.crearUblasRef(), mc;
	mc = boost::numeric::ublas::prod(ma, mb);
	*this = mc;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
bool L_MatrizLongDouble::invertMe_ublas()
{
	boost::numeric::ublas::matrix<long double> input = crearUblasRef(), output(li, lj);
	bool ret;
	ret = boost::numeric::ublas::L_InvertMatrix(input, output);
	if (ret == false)
		return false;
	*this = output;
	return ret;
}
#endif // __COMPAT_UBLAS__


void L_Matrix::OP_changeSign(const L_Matrix &other)
{
	int i, j;
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j) = -operator()(i,j);
}



bool L_Matrix::isGoodCovarianceMatrix(double factorAsi, bool print, const char *name)
{
	int i, j;
	const char *retr = "\b";
	int tipo = 0;
	double val, errMax=0;
	int iA, jA;
	double val1, val2;
	throw_L_ArgException_if(begin() == NULL || li!=lj, "L_Matrix::isGoodCovarianceMatrix() : matriz indefinida o no cuadrada");
	if (name == NULL)
		name = retr;
	for (i=0; i<li; i++)
	{
		// Diagonal > 0
		if (tipo <=3 && operator()(i,i) < 0)
		{
			if (tipo < 3)
				errMax = 0;
			tipo = 3;
			if (-operator()(i,i) > errMax)
			{
				iA = i; jA = i;
				errMax = -operator()(i,i);
			}
		}
		for (j=0; j<lj; j++)
		{
			// covarianzas propias > covarianzas cruzadas
			if (tipo <= 2 && operator()(i,j)*operator()(i,j) > operator()(i,i)*operator()(j,j))
			{
				if (tipo < 2)
					errMax = 0;
				tipo = 2;
				val = operator()(i,j)*operator()(i,j) - operator()(i,i)*operator()(j,j);
				if (val > errMax)
				{
					iA = i; jA = j;
					val1 = operator()(i,j)*operator()(i,j);
					val2 = operator()(i,i)*operator()(j,j);
					errMax = operator()(i,j)*operator()(i,j) - operator()(i,i)*operator()(j,j);
				}
			}
			// Simetria
			val = fabs(operator()(i,j) - operator()(j,i)) / (fabs(operator()(i,i)) + fabs(operator()(j,j)));
			if (tipo <= 1 && val > factorAsi)
			{
				if (tipo < 1)
					errMax = 0;
				tipo = 1;
				if (val > errMax)
				{
					iA = i; jA = j;
					errMax = val;
				}
			}
		}
	}
	if (print)
	{
		if (tipo == 1)
			L_printError("\noperator()(i,i) < 0 (%g < 0) en %s : pos(%d, %d) de (%d x %d)\n", val, name, iA, jA, li, lj);
		else if (tipo == 2)
			L_printError("\noperator()(i,i)*operator()(j,j) < operator()(i,j)*operator()(j,i) (%g < %g) en %s : pos(%d, %d) de (%d x %d)\n", val1, val2, name, iA, jA, li, lj);
		else if (tipo == 3)
			L_printError("\nAsimetria : elem-elem / elem+elem (%g) en %s : pos(%d, %d) de (%d x %d)\n", val, name, iA, jA, li, lj);
	}
	return tipo == 0;
}


bool L_Matrix::isGoodCovarianceMatrix(int &iErr, int &jErr, double factorAsi, bool print, const char *name)
{
	int i, j;
	const char *retr = "\b";
	double errMax = 0;
 	throw_L_ArgException_if(begin() == NULL || li!=lj, "L_Matrix::isGoodCovarianceMatrix() : matriz indefinida o no cuadrada");
	if (name == NULL)
		name = retr;
	for (i=0; i<li; i++)
	{
		// Diagonal > 0
		if (operator()(i,i) < 0)
		{
			if (print)
				L_printError("\noperator()(i,i) < 0 (%g < 0) en %s (%d, %d) de (%d x %d)\n", operator()(i,i), name, i, i, li, lj);
			iErr = i;
			jErr = j;
			return false;
		}
		for (j=0; j<lj; j++)
		{
			// covarianzas propias > covarianzas cruzadas
			if (operator()(i,j)*operator()(i,j) > operator()(i,i)*operator()(j,j))
			{
				if (operator()(i,j)*operator()(i,j) - operator()(i,i)*operator()(j,j) > errMax)
				{
					iErr = i;
					jErr = j;
					errMax = operator()(i,j)*operator()(i,j) - operator()(i,i)*operator()(j,j);
				}
			}
		}
	}
	if (errMax > 0)
		return false;
	return true;
}

bool L_Matrix::checkRangeValues(double limInf, double limSup, bool print, const char *name)
{
	int i, j;
	const char *retr = "\b";
	if (name == NULL)
		name = retr;

	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			if (operator()(i,j) != operator()(i,j))
			{
				if (print)
					L_printError("\nNaN en matriz %s de %d x %d en elemento (%d, %d)\n", name, li, lj, i, j);
				return false;
			}
			if (operator()(i,j) < limInf)
			{
				if (print)
					L_printError("\noperator()(i,j) < %f en matriz %s de %d x %d\n en elemento (%d, %d)\n", limInf, name, li, lj, i, j);
				return false;
			}
			if (operator()(i,j) > limSup)
			{
				if (print)
					L_printError("\noperator()(i,j) > %f en matriz %s de %d x %d\n en elemento (%d, %d)\n", limSup, name, li, lj, i, j);
				return false;
			}
		}
	}
	return true;
}

bool L_Matrix::isPositiveDefinite()
{
	L_Matrix yo, D;
	int i;

	yo = *this;
	yo.vectPr_sim(D);
	
	for (i=0; i<D.li; i++)
		if (D(i,0) <= 0)
			return false;
	return true;
}

double L_Matrix::computeMaximumSubtraction(L_Matrix &C)
{
	// ambas deben ser simetricas
	double lambdaDer = -1e60, lambdaIzq = 1e60, lambda1, lambda2, discr;
	double Mii, Mij, Mjj, Cii, Cij, Cjj;
	int i, j;
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			Mii = operator()(i,i);
			Mij = operator()(i,j);
			Mjj = operator()(j,j);
			Cii = C(i,i);
			Cij = C(i,j);
			Cjj = C(j,j);

			if (Cii == 0 && Cij == 0 && Cjj == 0)
				continue;
			if (Cii*Cjj == Cij*Cij)
			{
				discr = Mii*Cjj+Cii*Mjj+2*Mij*sqrt(Cii*Cjj);
				if (discr == 0)
					printf("No puede ser...\n");
				lambda1 = lambda2 = (Mii*Mjj-Mij*Mij)/(Mii*Cjj+Cii*Mjj+2*Mij*sqrt(Cii*Cjj));
			}
			else
			{
				double t1 = (Cii * Cjj);
				double t2 = (Cij * Cij);
				double t5 = (Mii * Cjj);
				double t6 = (Cii * Mjj);
				double t7 = (Cij * Mij);
				double t9 = (Mii * Mii);
				double t10 = (Cjj * Cjj);
				double t16 = (Cii * Cii);
				double t17 = (Mjj * Mjj);
				double t21 = (Mij * Mij);
				discr = (t9 * t10 - 2 * t5 * t6 + 4 * t5 * t7 + t16 * t17 + 4 * t6 * t7 + 4 * t1 * t21 + 4 * t2 * Mii * Mjj);
				if (discr < 0)
					printf("No puede ser...\n");
				double t28 = sqrt(discr);
				lambda1 = 0.1e1 / (t1 - t2) * (t5 + t6 + (2 * t7) - t28) / 0.2e1;
				lambda2 = 0.1e1 / (t1 - t2) * (t5 + t6 + (2 * t7) + t28) / 0.2e1;
			}
			if (lambda1 > lambdaDer)
				lambdaDer = lambda1;
			if (lambda2 > lambdaDer)
				lambdaDer = lambda2;
			if (lambda1 < lambdaIzq)
				lambdaIzq = lambda1;
			if (lambda2 < lambdaIzq)
				lambdaIzq = lambda2;
		}
	}
	// El rango (lambdaMin,lambdaMax) para el factor de la ganancia de Kalman
	// asegura que la matriz de covarianza resultante tiene sus submatrices de 2x2 def pos0+
	// Devolver el mayor lambda valido
	return lambdaDer;
}

void L_Matrix::concatenateSubMatrices(const L_Matrix &aa, const L_Matrix &ab, const L_Matrix &ba, const L_Matrix &bb)
{
	int i, j;
	throw_L_ArgException_if (aa.li != ab.li || aa.lj != ba.lj || ba.li != bb.li || ab.lj != bb.lj, "L_Matrix::concatenateSubMatrices(m,m,m,m)");
	reallocate(aa.li+ba.li, aa.lj+ab.lj);
	for (i=0; i<aa.li; i++)
		for (j=0; j<aa.lj; j++)
			operator()(i,j) = aa(i,j);
	for (i=0; i<ab.li; i++)
		for (j=0; j<ab.lj; j++)
			operator()(i,aa.lj+j) = ab(i,j);
	for (i=0; i<ba.li; i++)
		for (j=0; j<ba.lj; j++)
			operator()(i+aa.li,j) = ba(i,j);
	for (i=0; i<bb.li; i++)
		for (j=0; j<bb.lj; j++)
			operator()(i+aa.li,aa.lj+j) = bb(i,j);
}

void L_Matrix::concatenateSubMatrices(int ni, int nj, const L_Matrix *m00,...)
{
	va_list ppp;
	int ii, jj, i, j;
	int lii=0, ljj=0;
	L_Array<L_Array<const L_Matrix *> > matr;
	std::vector<int> lis, ljs;
	matr.resize(ni);
	for (ii=0; ii<ni; ii++)
	{
		matr[ii].resize(nj);
	}
	matr[0][0] = m00;
	va_start(ppp, m00);
	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			if (ii!=0 || jj!=0)
				matr[ii][jj] = va_arg(ppp, const L_Matrix *);
	va_end(ppp);

	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			throw_L_ArgException_if(matr[ii][jj]->li != matr[ii][0]->li || matr[ii][jj]->lj != matr[0][jj]->lj, "L_Matrix::concatenateSubMatrices");

	lis.resize(li);
	ljs.resize(lj);

	lis[0] = 0;
	lii = matr[0][0]->li;
	for (ii=1; ii<ni; ii++)
	{
		lis[ii] = lis[ii-1] + matr[ii-1][0]->li;
		lii = lii + matr[ii][0]->li;
	}
	ljs[0] = 0;
	ljj = matr[0][0]->lj;
	for (jj=1; jj<nj; jj++)
	{
		ljs[jj] = ljs[jj-1] + matr[0][jj-1]->lj;
		ljj = ljj + matr[0][jj]->lj;
	}

	reallocate(lii, ljj);

	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			for (i=0; i<matr[ii][jj]->li; i++)
				for (j=0; j<matr[ii][jj]->lj; j++)
					operator()(lis[ii]+i,ljs[jj]+j) = matr[ii][jj]->operator()(i,j);
}

/*
void L_Matrix::separarEnSubMatrices(int ni, int nj, L_Matrix *m00,...)
{
	va_list ppp;
	int ii, jj, i, j;
	int lii=0, ljj=0;
	L_Array<L_Array<const L_Matrix *> > matr;
	std::vector<int> lis, ljs;
	matr.resize(ni);
	for (ii=0; ii<ni; ii++)
	{
		matr[ii].resize(nj);
	}
	matr[0][0] = m00;
	va_start(ppp, m00);
	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			if (ii!=0 || jj!=0)
				matr[ii][jj] = va_arg(ppp, const L_Matrix *);
	va_end(ppp);

	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			throw_L_ArgException_if(matr[ii][jj]->li != matr[ii][0]->li || matr[ii][jj]->lj != matr[0][jj]->lj, "L_Matrix::separarEnSubMatrices");

	lis.resize(li);
	ljs.resize(lj);

	lis[0] = 0;
	lii = matr[0][0]->li;
	for (ii=1; ii<ni; ii++)
	{
		lis[ii] = lis[ii-1] + matr[ii-1][0]->li;
		lii = lii + matr[ii][0]->li;
	}
	ljs[0] = 0;
	ljj = matr[0][0]->lj;
	for (jj=1; jj<nj; jj++)
	{
		ljs[jj] = ljs[jj-1] + matr[0][jj-1]->lj;
		ljj = ljj + matr[0][jj]->lj;
	}

	throw_L_ArgException_if(lii != li || ljj != lj, "L_Matrix::separarEnSubMatrices");

	for (ii=0; ii<ni; ii++)
		for (jj=0; jj<nj; jj++)
			for (i=0; i<matr[ii][jj]->li; i++)
				for (j=0; j<matr[ii][jj]->lj; j++)
					matr[ii][jj]->elem[i][j] = elem[lis[ii]+i][ljs[jj]+j];
}
*/

void L_Matrix::covarianceThroughJacobian(const std::vector<int> &a, const L_Matrix &J, const std::vector<int> &b)
{
	L_Matrix Paa((int)a.size(),(int)a.size()), PaM((int)a.size(),li-(int)a.size()), PMa(li-(int)a.size(),(int)a.size()), PMM(li-(int)a.size(),li-(int)a.size());
	L_Matrix Pbb, PbN, PNb; // PNN = PMM
	L_Matrix JT, JPaa;
	std::vector<int> M, N;
	int i, j;

	L_complementOf(a, M, li);
	L_complementOf(b, N, (int)(li+b.size()-a.size()));

	for (i=0; i<(int)a.size(); i++)
		for (j=0; j<(int)a.size(); j++)
			Paa(i,j) = operator()(a[i],a[j]);
	for (i=0; i<(int)a.size(); i++)
		for (j=0; j<(int)M.size(); j++)
			PaM(i,j) = operator()(a[i],M[j]);
	for (i=0; i<(int)M.size(); i++)
		for (j=0; j<(int)a.size(); j++)
			PMa(i,j) = operator()(M[i],a[j]);
	for (i=0; i<(int)M.size(); i++)
		for (j=0; j<(int)N.size(); j++)
			PMM(i,j) = operator()(M[i],M[j]);
	
	JT.transpOf(J);

	JPaa.OP_mult(J,Paa);
	Pbb.OP_mult(JPaa,JT);
	PbN.OP_mult(J,PaM);
	PNb.OP_mult(PMa,JT);

	reallocate((int)(b.size() + N.size()), (int)(b.size() + N.size()));

	for (i=0; i<(int)b.size(); i++)
		for (j=0; j<(int)b.size(); j++)
			operator()(b[i],b[j]) = Pbb(i,j);
	for (i=0; i<(int)b.size(); i++)
		for (j=0; j<(int)N.size(); j++)
			operator()(b[i],N[j]) = PbN(i,j);
	for (i=0; i<(int)N.size(); i++)
		for (j=0; j<(int)b.size(); j++)
			operator()(N[i],b[j]) = PNb(i,j);
	for (i=0; i<(int)N.size(); i++)
		for (j=0; j<(int)N.size(); j++)
			operator()(N[i],N[j]) = PMM(i,j);

}


void L_Matrix::integralOf(const L_Matrix &other)
{
	int i, j;
	double sFila = 0;
	reallocate(other.li,other.lj);
	operator()(0,0)=other(0,0);
	for (j=1; j<lj; j++)
		operator()(0,j)=operator()(0,j-1)+other(0,j);
	for (i=1; i<li; i++)
	{
		sFila = 0;
		for (j=0; j<lj; j++)
		{
			sFila += other(i,j);
			operator()(i,j) = operator()(i-1,j) + sFila;
		}
	}
}


void L_Matrix::integralSquaresOf(const L_Matrix &other)
{
	int i, j;
	double sFila = 0;
	reallocate(other.li,other.lj);
	operator()(0,0)=other(0,0)*other(0,0);
	for (j=1; j<lj; j++)
		operator()(0,j)=operator()(0,j-1)+other(0,j)*other(0,j);
	for (i=1; i<li; i++)
	{
		sFila = 0;
		for (j=0; j<lj; j++)
		{
			sFila += other(i,j)*other(i,j);
			operator()(i,j) = operator()(i-1,j) + sFila;
		}
	}
}

L_Matrix& L_Matrix::transpOf(const L_Matrix &other)
{
	int i, j;
	reallocate(other.lj, other.li);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=other(j,i);
	return *this;
}


double L_Matrix::det()
{
	int i, j, k;
	L_Matrix m;  
	double det = 1, tmp;

	m.reallocate(li, li);

	for ( i = 0; i < li; i++ ) {
		for ( j = 0; j < li; j++ )
		m(i,j) = operator()(i,j);
	}

	for ( k = 0; k < li; k++ ) {
		if ( m(k,k) == 0 ) {
			bool ok = false;

			for ( j = k; j < li; j++ ) {
				if ( m(j,k) != 0 )
				{
					ok = true;  /// ERROR: faltaba el break, es cuatico esto
					break;
				}
			}

			if ( !ok )
			return 0;

			for ( i = k; i < li; i++ )
			{
				tmp = m(i,j);
				m(i,j) = m(i,k);
				m(i,k) = tmp;
			}
			det = -det;
		}

		det *= m(k,k);

		if ( k + 1 < li ) {
		for ( i = k + 1; i < li; i++ ) {
			for ( j = k + 1; j < li; j++ )
				m(i,j) = m(i,j) - m(i,k) * 
				m(k,j) / m(k,k);
			}
		}
	}
	return det;
}


void L_Matrix::subMatrDeleteRowCol_at(L_Matrix &ret, int i0, int j0)
{
	int i, j;
	ret.reallocate(li-1, lj-1);
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
		{
			ret(i,j)=operator()(i,j);
		}
		for (j=j0; j<lj-1; j++)
		{
			ret(i,j)=operator()(i,j+1);
		}
	}
	for (i=i0; i<li-1; i++)
	{
		for (j=0; j<j0; j++)
		{
			ret(i,j)=operator()(i+1,j);
		}
		for (j=j0; j<lj-1; j++)
		{
			ret(i,j)=operator()(i+1,j+1);
		}
	}
}

void L_Matrix::insertRowsColumns(int i0, int ni, int j0, int nj, const L_Matrix &m)
{
	int i, j;
	reallocate(m.li+ni, m.lj+nj);
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
			operator()(i,j) = m(i,j);
		for (j=j0; j<j0+nj; j++)
			operator()(i,j) = 0;
		for (j=j0+nj; j<lj; j++)
			operator()(i,j) = m(i,j-nj);
	}
	for (i=i0; i<i0+ni; i++)
		for (j=0; j<lj; j++)
			operator()(i,j) = 0;
	for (i=i0+ni; i<li; i++)
	{
		for (j=0; j<j0; j++)
			operator()(i,j) = m(i-ni,j);
		for (j=j0; j<j0+nj; j++)
			operator()(i,j) = 0;
		for (j=j0+nj; j<lj; j++)
			operator()(i,j) = m(i-ni,j-nj);
	}
}

void L_Matrix::deleteRowsColumns(int *indi, int ni, int *indj, int nj, const L_Matrix &m)
{
	reallocate(m.li-ni, m.lj-nj);

	int i, u=0, j, v;
	for (i=0; i<li; i++)
	{
		while (u < ni && indi[u] == i+u)
			u++; // No copiar
		v = 0;
		for (j=0; j<lj; j++)
		{
			while (v < nj && indj[v] == j+v)
				v++; // No copiar
			operator()(i,j) = m(i+u,j+v);
		}
	}
}

void L_Matrix::deleteRowsColumns(int i0, int ni, int j0, int nj, const L_Matrix &m)
{
	reallocate(m.li-ni, m.lj-nj);

	int i, j;
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
			operator()(i,j) = m(i,j);
		for (j=j0; j<m.lj-nj; j++)
			operator()(i,j) = m(i,j+nj);
	}
	for (i=i0; i<m.li-ni; i++)
	{
		for (j=0; j<j0; j++)
			operator()(i,j) = m(i+ni,j);
		for (j=j0; j<m.lj-nj; j++)
			operator()(i,j) = m(i+ni,j+nj);
	}
}


void L_Matrix::OP_assign(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.begin()==NULL, "L_Matrix::OP_assign"); // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	reallocate(other.li, other.lj);
	if (li*lj < L_TAM_COPIA_MEMCPY)
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=other(i,j);
	}
	else
		memcpy(begin(),other.begin(),sizeof(double)*li*lj);
	return;
}

void L_Matrix::OP_assignInvertingSigns(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.begin()==NULL, "L_Matrix::OP_assignInvertingSigns"); // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=-other(i,j);
	return;
}

void L_Matrix::OP_add(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.li!=li || other.lj!=lj, "L_Matrix::OP_add");
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)+=other(i,j);
	return;
}

void L_Matrix::OP_subtract(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(begin()==NULL || other.begin()==NULL || other.li!=li || other.lj!=lj, "L_Matrix::OP_subtract");
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)-=other(i,j);
	return;
}

void L_Matrix::OP_multByDiagonal(const L_Matrix &diagonal)
{
	int i, j;
	throw_L_ArgException_if(begin()==NULL || diagonal.begin()==NULL || li==0 || lj != diagonal.li, "L_Matrix::OP_multByDiagonal");
	if (diagonal.lj==1) // La matriz diagonal esta representada como vector (para ahorrar espacio)
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)*=diagonal(j,0);
	}
	else // La matriz diagonal esta representada realmente como matriz diagonal
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)*=diagonal(j,j);
	}
}

void L_Matrix::OP_multByDiagonal(const L_Matrix &normal, const L_Matrix &diagonal)
{
	int i, j;
	throw_L_ArgException_if(normal.begin()==NULL || diagonal.begin()==NULL || normal.lj!=diagonal.li || (diagonal.lj!=1 && diagonal.li!=diagonal.lj), "L_Matrix::OP_multByDiagonal_NULL");
	reallocate(normal.li,normal.lj);
	if (diagonal.lj==1) // La matriz diagonal esta representada como vector (para ahorrar espacio)
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=normal(i,j)*diagonal(j,0);
	}
	else // La matriz diagonal esta representada realmente como matriz diagonal
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=normal(i,j)*diagonal(j,j);
	}
}

void L_Matrix::OP_div(const L_Matrix &other)
{
	L_Matrix aux;
	aux=other;

	if (lj == 1)
	{
		if (!solveLinearSystem_x_Ab_result_in_b(aux))
			throw L_ZeroDivisionException();
		return;
	}
	if (!aux.invertMe())
		throw L_ZeroDivisionException();
	aux.swap(*this);
	return;
}



void L_Matrix::OP_add(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_add()");
	if (begin()==NULL || li!=m1.li || lj!=m1.lj)
		reallocate(m1.li, m1.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=m1(i,j)+m2(i,j);
	return;
}

void L_Matrix::OP_subtract(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_subtract()");
	if (begin()==NULL || li!=m1.li || lj!=m1.lj)
		reallocate(m1.li, m1.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=m1(i,j)-m2(i,j);
	return;
}

#ifdef OJO______USE_OPTIMIZED_MULTIPLICATION
void L_Matrix::OP_mult(const L_Matrix &m1, const L_Matrix &m2)
{
	throw_L_ArgException_if(m1.lj!=m2.li || m1.begin() == begin() || m2.begin()==begin(), "L_Matrix::OP_mult()");
#if defined (__COMPAT_EIGEN__)
	OP_mult_eigen(m1, m2);
	return;
#elif defined(__COMPAT_UBLAS__)
	OP_mult_ublas(m1, m2);
	return;
#endif // __COMPAT_UBLAS__
	if (m1.li< 200 || m1.lj < 200)
		OP_mult_private(m1, m2);
	else
	{
		L_Matrix m2T;
		m2T.transpOf(m2);
		OP_mult_ABT(m1, m2T);
		return;
	}
}
#else
void L_Matrix::OP_mult(const L_Matrix &m1, const L_Matrix &m2) {OP_mult_private(m1, m2);}
#endif // OJO______USE_OPTIMIZED_MULTIPLICATION

void L_Matrix::OP_mult_private(const L_Matrix &m1, const L_Matrix &m2)
{
	int i,j,u;
	double r;

	throw_L_ArgException_if(m1.lj!=m2.li, "L_Matrix::OP_mult()");
	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=0; j<m2.lj; j++)
		{
			r=0;
			for (u=0; u<m1.lj; u++)
				r+=m1(i,u)*m2(u,j);
			operator()(i,j) = r;
		}
	}
	return;
}


void L_Matrix::OP_mult_ABT(const L_Matrix &A, const L_Matrix &B)
{
	throw_L_ArgException_if(A.lj!=B.lj || A.begin()==NULL || B.begin()==NULL || A.begin() == begin() || B.begin()==begin(), "L_Matrix::OP_mult_ABT()");
#if defined(__COMPAT_EIGEN__)   &&   defined(OJO______USE_OPTIMIZED_MULTIPLICATION)
	Eigen::MatrixXd mat1 = A.crearEigenRef(), mat2 = B.crearEigenRef();
	*this = mat1*mat2.transpose();
	return;
#elif defined(__COMPAT_UBLAS__)   &&   defined(OJO______USE_OPTIMIZED_MULTIPLICATION)
	boost::numeric::ublas::matrix<double> ma = A.crearUblasRef(), mb = B.crearUblasRef(), mc;
	mc = boost::numeric::ublas::prod(ma, boost::numeric::ublas::trans(mb));
	*this = mc;
	return;
#else
	int i, j, k;
	double r;
	const double * L_restrict q1;
	const double * L_restrict q2;
	throw_L_ArgException_if(A.data() == 0 || B.data() == 0, "L_Matrix::OP_mult_ABT() : matriz argumento no inicializada");
	reallocate(A.li, B.li);
	for (i=0; i<A.li; i++)
	{
		q1 = &(A(i,0));
		for (j=0; j<B.li; j++)
		{
			q2 = &(B(j,0));
			r = 0;
			for (k=0; k<A.lj; k++)
				r += q1[k]*q2[k];
			operator()(i,j) = r;
		}
	}
#endif
}

// Permite calcular A*P*AT de forma numericamente simetrica pq usando double, (A*P)*AT != A*(P*AT)
void L_Matrix::OP_mult_APAT_cholesky(const L_Matrix &A, const L_Matrix &P, bool reparar)
{
	L_Matrix L, AL;
	L_Matrix VP, vp, AT;
	int i;
	if (L.cholesky_L_De(P) == true)
	{
		AL.OP_mult(A,L);
		OP_mult_ABT(AL,AL);
	}
	else if (reparar)
	{
		VP = P;
		VP.vectPr_sim(vp);
		for (i=0; i<P.li; i++)
		{
			if (vp(i,0) < 0)
				vp(i,0) = 1.0e-30;
			vp(i,0) = sqrt(vp(i,0));
		}
		AT.OP_mult(A,VP);
		AL.OP_mult(AT,vp);
		OP_mult_ABT(AL,AL);
	}
	else
	{
		AT.transpOf(A);
		VP.OP_mult(A,P);
		OP_mult(VP, AT);
	}
}


void L_Matrix::OP_mult_sorting(const L_Matrix &m1, const L_Matrix &m2)
{
	int i,j,u;
	double res;
	std::vector<double> sumitas(m1.lj);

	throw_L_ArgException_if(m1.lj!=m2.li || m1.begin() == begin() || m2.begin()==begin(), "L_Matrix::OP_mult_sorting()");
	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=0; j<m2.lj; j++)
		{
			for (u=0; u<m1.lj; u++)
				sumitas[u] = m1(i,u)*m2(u,j);
			std::sort(sumitas.begin(), sumitas.end()); //sumitas.sort(L_Int::cmp);
			res = 0;
			for (u=0; u<m1.lj; u++)
				res += sumitas[u];
			operator()(i,j) = res;
		}
	}
	return;
}

void L_Matrix::OP_mult_longDouble(const L_Matrix &m1, const L_Matrix &m2)
{
	L_MatrizLongDouble M1, M2, M3;
	M1.copyFrom(m1);
	M2.copyFrom(m2);
	M3.OP_mult(M1, M2);
	M3.copyTo(*this);
	return;
}


void L_Matrix::OP_mult_sub(const L_Matrix &m1, const L_Matrix &m2, int li1, int lj1, int li2, int lj2)
{
	int i,j,u;
	double res;

	throw_L_ArgException_if(lj1!=li2, "L_Matrix::OP_mult_sub() : indices incorrectos");
	if (begin()==NULL || li!=li1 || lj!=lj2)
		reallocate(li1, lj2);

	for (i=0; i<li1; i++) for (j=0; j<lj2; j++)
	{
		res=0;
		for (u=0; u<li2; u++)
			res+=m1(i,u)*m2(u,j);
		operator()(i,j) = res;
	}
	return;
}

void L_Matrix::OP_mult_sub(int i0, int j0, const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1, int i02, int j02, int li2, int lj2)
{
	int i,j,u;
	double res;
	throw_L_ArgException_if(lj1!=li2 || begin()==NULL || li<i01+li1 || lj<j02+lj2, "L_Matrix::OP_mult_sub() : indices incorrectos"); // No es la idea

	for (i=0; i<li1; i++) for (j=0; j<lj2; j++)
	{
		res = 0;
		for (u=0; u<li2; u++)
			res+=m1(i+i01,u+j01)*m2(u+i02,j+j02);
		operator()(i+i0,j+j0) = res;
	}
	return;
}

void L_Matrix::OP_mult_sub_sum(int i0, int j0, const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1, int i02, int j02, int li2, int lj2)
{
	int i,j,u;

	// Multiplicar m1*m2 desplazando el resultado

	throw_L_ArgException_if(lj1!=li2 || begin()==NULL || li<i01+li1 || lj<j02+lj2, "L_Matrix::OP_mult_sub_sum()");

	for (i=0; i<li1; i++) for (j=0; j<lj2; j++)
	{
		for (u=0; u<li2; u++)
			operator()(i+i0,j+j0)+=m1(i+i01,u+j01)*m2(u+i02,j+j02);
	}
	return;
}

void L_Matrix::OP_mult_sub1_sum(const L_Matrix &m1, const L_Matrix &m2, int i01, int j01, int li1, int lj1)
{
	int i,j,u;

	// Multiplicar m1*m2 sin desplazar el resultado
	// La restriccion de zona se indica para m2

	throw_L_ArgException_if (m1.lj != m2.li || m1.li < i01+li1 || m1.lj < j01+lj1 || i01 < 0 || j01 < 0 || li1 < 0 || lj1 < 0, "L_Matrix::OP_mult_sub1_sum() : indices incorrectos");

	reallocate(m1.li, m2.lj);

	for (i=i01; i<i01+li1; i++)
	{
		for (j=0; j<m2.lj; j++)
		{
			for (u=j01; u<j01+lj1; u++)
				operator()(i,j)+=m1(i,u)*m2(u,j);
		}
	}
	return;
}

void L_Matrix::OP_mult_sub2_sum(const L_Matrix &m1, const L_Matrix &m2, int i02, int j02, int li2, int lj2)
{
	int i,j,u;

	// Multiplicar m1*m2 sin desplazar el resultado
	// La restriccion de zona se indica para m2

	throw_L_ArgException_if (m1.lj != m2.li || m2.li < i02+li2 || m2.lj < j02+lj2 || i02 < 0 || j02 < 0 || li2 < 0 || lj2 < 0, "L_Matrix::OP_mult_sub1_sum() : indices incorrectos");

	reallocate(m1.li, m2.lj);

	for (i=0; i<m1.li; i++)
	{
		for (j=j02; j<j02+lj2; j++)
		{
			for (u=i02; u<i02+li2; u++)
				operator()(i,j)+=m1(i,u)*m2(u,j);
		}
	}
	return;
}


void L_Matrix::OP_mult(double val)
{
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)*=val;
	return;
}

void L_Matrix::OP_amplify(double val)
{
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)*=val;
	return;
}

void L_Matrix::OP_div(double val)
{
	int i, j;
	if (val==0)
		throw L_ZeroDivisionException();
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)/=val;
	return;
}

void L_Matrix::OP_productElementwise(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_productElementwise");
	if (li!=m1.li || lj!=m1.lj)
		reallocate(m1.li, m1.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=m1(i,j)*m2(i,j);
	return;
}



void L_Matrix::OP_assign_NULL(const L_Matrix & other)
{
	int i, j;
	if (other.begin()==NULL) // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	{
		destroy();
		li = other.li;
		lj = other.lj;
		return;
	}
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=other(i,j);
	return;
}

void L_Matrix::OP_assignInvertingSigns_NULL(const L_Matrix & other)
{
	int i, j;
	if (other.begin()==NULL) // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	{
		destroy();
		li = other.li;
		lj = other.lj;
		return;
	}
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=-other(i,j);
	return;
}

void L_Matrix::OP_add_NULL(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.li!=li || other.lj!=lj, "L_Matrix::OP_add_NULL");
	if (other.begin() == NULL)
		return;
	else if (begin() == NULL)
		OP_assign(other);
	else
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)+=other(i,j);
	return;
}

void L_Matrix::OP_subtract_NULL(const L_Matrix & other)
{
	int i, j;
	throw_L_ArgException_if(other.li!=li || other.lj!=lj, "L_Matrix::OP_subtract_NULL");
	if (other.begin() == NULL)
		return;
	else if (begin() == NULL)
		OP_assignInvertingSigns(other);
	else
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)-=other(i,j);
	return;
}

void L_Matrix::OP_multByDiagonal_NULL(const L_Matrix &diagonal)
{
	int i, j;
	throw_L_ArgException_if(lj != diagonal.li, "L_Matrix::OP_multByDiagonal_NULL");
	if (diagonal.begin() == NULL)
	{
		destroy();
		return;
	}
	else if (begin() == NULL)
		return;
	if (diagonal.lj==1) // La matriz diagonal esta representada como vector (para ahorrar espacio)
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)*=diagonal(j,0);
	}
	else // La matriz diagonal esta representada realmente como matriz diagonal
	{
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)*=diagonal(j,j);
	}
}

void L_Matrix::OP_multByDiagonal_NULL(const L_Matrix &normal, const L_Matrix &diagonal)
{
	int i, j;

	throw_L_ArgException_if(normal.lj!=diagonal.li || (diagonal.lj!=1 && diagonal.li!=diagonal.lj), "L_Matrix::OP_multByDiagonal_NULL");
	if (normal.begin() == NULL)
	{
		destroy();
		li = normal.li;
		lj = normal.lj;
		return;
	}
	else if (diagonal.begin() == NULL)
	{
		li = normal.li;
		lj = normal.lj;
		return;
	}

	reallocate(normal.li,normal.lj);
	if (diagonal.lj==1) // La matriz diagonal esta representada como vector (para ahorrar espacio)
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=normal(i,j)*diagonal(j,0);
	else // La matriz diagonal esta representada realmente como matriz diagonal
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j)=normal(i,j)*diagonal(j,j);
}


void L_Matrix::OP_add_NULL(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_add_NULL");
	if (m1.begin() == NULL)
	{
		if (m2.begin() == NULL)
		{
			destroy();
			li = m1.li;
			lj = m1.lj;
		}
		else
			OP_assign(m2);
	}
	else
	{
		if (m2.begin() == NULL)
			OP_assign(m1);
		else
		{
			reallocate(m1.li, m1.lj);
			for (i=0; i<li; i++)
				for (j=0; j<lj; j++)
					operator()(i,j)=m1(i,j)+m2(i,j);
		}
	}
	return;
}

void L_Matrix::OP_subtract_NULL(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_subtract_NULL");
	if (m1.begin() == NULL)
	{
		if (m2.begin() == NULL)
		{
			destroy();
			li = m1.li;
			lj = m1.lj;
		}
		else
			OP_assignInvertingSigns(m2);
	}
	else
	{
		if (m2.begin() == NULL)
			OP_assign(m1);
		else
		{
			reallocate(m1.li, m1.lj);
			for (i=0; i<li; i++)
				for (j=0; j<lj; j++)
					operator()(i,j)=m1(i,j)-m2(i,j);
		}
	}
	return;
}

void L_Matrix::OP_mult_NULL(const L_Matrix &m1, const L_Matrix &m2)
{
	int i_m1,j_m2,u;

	throw_L_ArgException_if(m1.lj!=m2.li, "L_Matrix::OP_mult_NULL");
	if (m1.begin() == NULL || m2.begin() == NULL)
	{
		destroy();
		li = m1.li;
		lj = m2.lj;
		return;
	}
	if (begin()==NULL || li!=m1.li || lj!=m2.lj)
		reallocate(m1.li, m2.lj);

	for (i_m1=0; i_m1<m1.li; i_m1++) for (j_m2=0; j_m2<m2.lj; j_m2++)
	{
		operator()(i_m1,j_m2)=0;
		for (u=0; u<m2.li; u++)
			operator()(i_m1,j_m2)+=m1(i_m1,u)*m2(u,j_m2);
	}
	return;
}


void L_Matrix::OP_mult_NULL(double val)
{
	if (begin() == NULL)
		return;
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)*=val;
	return;
}
void L_Matrix::OP_div_NULL(double val)
{
	int i, j;
	if (val==0)
		throw L_ZeroDivisionException();
	if (begin() == NULL)
		return;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)/=val;
	return;
}

void L_Matrix::OP_productElementwise_NULL(const L_Matrix &m1, const L_Matrix &m2)
{
	int i, j;
	throw_L_ArgException_if(m1.li!=m2.li || m1.lj!=m2.lj, "L_Matrix::OP_productElementwise_NULL");
	if (m1.begin() == NULL || m2.begin() == NULL)
	{
		destroy();
		li = m1.li;
		lj = m1.lj;
	}
	if (li!=m1.li || lj!=m1.lj)
		reallocate(m1.li, m1.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=m1(i,j)*m2(i,j);
	return;
}

bool L_Matrix::invertMe() {return invertMe_private();}  // Es mas rapida que las otras

bool L_Matrix::invertMe_private()
{
	int i, j, c;
	double factor;
	double temp;
	throw_L_ArgException_if(li!=lj, "L_Matrix::invertMe");
	// Crear matriz unitaria other
	L_Matrix other;
	other.reallocate(li, lj);
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			if (i==j)
				other(i,j)=1;
			else
				other(i,j)=0;
		}
	}
	for (c=0; c<lj; c++) // Limpiar columna c
	{
		if (operator()(c,c)==0) // Diagonal con zero, intercambiar por other fila mas abajo que c
		{
			for (i=c+1; i<li; i++) // Encontrar fila i con elemento != 0 en la columna c
				if (operator()(i,c)!=0)
					break;
			if (i == li) // fila i solo contiene elementos zero para j>=c
				return false;
			for (j=0; j<c; j++) // intercambiar fila c por fila i. Notar que operator()(i,j) = 0 para j<c
			{
				temp = other(i,j);
				other(i,j)=other(c,j);
				other(c,j) = temp;
			}
			for (; j<lj; j++) // intercambiar
			{
				temp = other(i,j);
				other(i,j)=other(c,j);
				other(c,j) = temp;
				temp = operator()(i,j);
				operator()(i,j)=operator()(c,j);
				operator()(c,j) = temp;
			}
		}
		// Amplificar fila c para que quede con 1 en la diagonal
		factor = 1 / operator()(c,c);
		for (j=0; j<c; j++)
		{
			other(c,j)*=factor;
		}
		for (; j<lj; j++)
		{
			operator()(c,j)*=factor;
			other(c,j)*=factor;
		}
		// Restarla a las otras filas para dejar la columna c "unitaria"
		for (i=0; i<lj; i++)
		{
			if (i==c)
				continue;
			factor = -operator()(i,c);
			for (j=0; j<c; j++)
				other(i,j)+=factor*other(c,j);
			for (; j<lj; j++)
			{
				operator()(i,j)+=factor*operator()(c,j);
				other(i,j)+=factor*other(c,j);
			}
		}
	}
	swap(other);
	return true;
}

// Probe paralelizar con #pragma omp parallel for shared(...) private(...), y es peor (por el cache supongo)
bool L_Matrix::invertMeDefPos()
{
    if (li == 0 || lj == 0)
	{
		printf("L_Matrix::invertMeDefPos() : matriz de 0x0\n");
		return false;
	}
    if (li*lj == 1)
	{
		operator()(0,0) = 1/operator()(0,0);
		return true;
	}
	double *data = begin();
	int i, j, k, jxlj, ixlj;
	double sum, x;

	for (i=1; i < li; i++)
		data[i] /= data[0]; // normalize row 0


	for (i=1; i < li; i++)
	{
		ixlj = i*lj;
		for (j=i; j < li; j++)
		{ // do a column of L
			jxlj = j*lj;
			sum = 0.0;
			for (k = 0; k < i; k++)  
				sum += data[jxlj+k] * data[k*lj+i];
			data[jxlj+i] -= sum;
		}
		if (i == li-1)
			continue;
		for (j=i+1; j < li; j++)
		{  // do a row of U
			sum = 0.0;
			for (k = 0; k < i; k++)
				sum += data[ixlj+k]*data[k*lj+j];
			data[ixlj+j] = 
			(data[ixlj+j]-sum) / data[ixlj+i];
		}
	}
	for ( i = 0; i < li; i++ )  // invert L
	{
		for ( j = i; j < li; j++ )
		{
			jxlj = j*lj;
			x = 1.0;
			if ( i != j )
			{
				x = 0.0;
				for ( k = i; k < j; k++ ) 
					x -= data[jxlj+k]*data[k*lj+i];
			}
			data[jxlj+i] = x / data[jxlj+j];
		}
	}

	for ( i = 0; i < li; i++ )   // invert U
	{
		ixlj = i*lj;
		for ( j = i+1; j < li; j++ )
		{
			sum = 0.0;
			for ( k = i; k < j; k++ )
				sum += data[k*lj+j]*( (i==k) ? 1.0 : data[ixlj+k] );
			data[ixlj+j] = -sum;
		}
	}
	for ( i = 0; i < li; i++ )   // final inversion
	{
		for ( j = 0; j < li; j++ )
		{
			jxlj = j*lj;
			sum = 0.0;
			for ( k = ((i>j)?i:j); k < li; k++ )  
				sum += ((j==k)?1.0:data[jxlj+k])*data[k*lj+i];
			data[jxlj+i] = sum;
		}
	}
	return true;
}


bool L_Matrix::invertMeSimetricaDefPos()
{
	L_Matrix L(li,lj), Li(li,lj), LiT(li,lj);
	int i, j, k;
	double factor;

	throw_L_ArgException_if(li != lj, "L_Matrix::invertMeSimetricaDefPos(): matriz no cuadrada");

	if (L.cholesky_L_De(*this) == false)
		return false;

	Li.identity();
	// Calcular la inversa Li destruyendo L, ambas son triangulares inferiores
	for (j=0; j<li; j++)
	{
		if (L(j,j) == 0)
			return false;
		factor = 1/L(j,j);
		for (k=0; k<=j; k++)
			L(j,k)*=factor;
		for (k=0; k<=j; k++)
			Li(j,k)*=factor;
		for (i=j+1; i<li; i++)
		{
			double factor = L(i,j);
			for (k=0; k<=j; k++)
				L(i,k) -= factor*L(j,k);
			for (k=0; k<=j; k++)
				Li(i,k) -= factor*Li(j,k);
		}
	}

	//Li.print("Li");

	// m  =  L*LT
	// m^-1 = (L*LT)^-1
	// m^-1 = LT^-1 * L^-1
	// m^-1 = L^-1T * L^-1

	LiT.transpOf(Li);

	// Multiplicacion de LiT * Li
	reallocate(li,lj);
	for (i=0; i<li; i++)
	{
		for (j=0; j<li; j++) // Decia j<=li
		{
			double sum = 0;
			for (k=i; k<li; k++)
				sum += LiT(i,k)*Li(k,j);
			operator()(i,j) = sum;
		}
	}

	return true;
}


void L_Matrix::invertMe_invertirValoresPropios(double minVal)
{
	L_Matrix P, D, PD, tmp;
	double maxVal;
	int i;
	swap(P);
	P.vectPr_sim(D);

	maxVal = 1/minVal;

	for (i=0; i<D.li; i++)
	{
		if (D(i,0) < minVal)
			D(i,0) = minVal;  // Para que no se destruya la matriz
		if (D(i,0) > maxVal)
			D(i,0) = maxVal;  // Para que no se destruya la matriz
		D(i,0) = 1/D(i,0);
		D(i,0) = sqrt(D(i,0));
	}
	PD.OP_multByDiagonal(P, D); // P * sqrt(D)
	OP_mult_ABT(PD,PD);
	return;
}



bool L_Matrix::diagonalize_me_wide()
{
	int i, j, c;
	int lmin;
	double factor;
	double temp;
	lmin = li;
	if (li > lj)
		lmin = lj; // No tiene sentido que deje puros ceros si tiene dimensiones inadecuadas
	for (c=0; c<lmin; c++) // Limpiar columna c para que quede triangular
	{
		if (operator()(c,c)==0) // Diagonal con zero, intercambiar por other fila mas abajo que c
		{
			for (i=c+1; i<li; i++) // Encontrar fila i con elemento != 0 en la columna c
				if (operator()(i,c)!=0)
					break;
			if (i == li) // fila i solo contiene elementos zero para j>=c. No se puede diagonalizar
				return false;
			for (j=c; j<lj; j++) // intercambiar. Notar que operator()(>c,<c) = 0
			{
				temp = operator()(i,j);
				operator()(i,j)=operator()(c,j);
				operator()(c,j) = temp;
			}
		}
		// Amplificar fila c para que quede con 1 en la diagonal
		factor = 1 / operator()(c,c);
		for (j=c; j<lj; j++)
		{
			operator()(c,j)*=factor;
		}
		// Ayudar a dejar triangular
		for (i=c+1; i<lmin; i++)
		{
			factor = -operator()(i,c);
			for (j=c; j<lj; j++)
				operator()(i,j)+=factor*operator()(c,j);
		}
	}
	for (c=lmin-1; c>0; c--) // Limpiar columna c para que quede triangular
	{
		for (i=0; i<c; i++)
		{
			factor = -operator()(i,c);
			for (j=0; j<lj; j++)
				operator()(i,j)+=factor*operator()(c,j);
		}
	}
	return true;
}


bool L_Matrix::diagonalizame_tall() // Es como diagonalize_me(), pero con [i,li) y [j,lj) invertidos
{
	int i, j, c;
	int lmin;
	double factor;
	double temp;
	lmin = lj;
	if (lj > li)
		lmin = li; // No tiene sentido que deje puros ceros si tiene dimensiones inadecuadas
	for (c=0; c<lmin; c++) // Limpiar fila c para que quede triangular
	{
		if (operator()(c,c)==0) // Diagonal con zero, intercambiar por other columna mas abajo que c
		{
			for (j=c+1; j<lj; j++) // Encontrar columna j con elemento != 0 en la fila c
				if (operator()(c,j)!=0)
					break;
			if (j == lj) // columna j solo contiene elementos zero para i>=c. No se puede diagonalizar
				return false;
			for (i=c; i<li; i++) // intercambiar. Notar que operator()(>c,<c) = 0
			{
				temp = operator()(i,j);
				operator()(i,j)=operator()(i,c);
				operator()(i,c) = temp;
			}
		}
		// Amplificar columna c para que quede con 1 en la diagonal
		factor = 1 / operator()(c,c);
		for (i=c; i<li; i++)
		{
			operator()(i,c)*=factor;
		}
		// Ayudar a dejar triangular
		for (j=c+1; j<lmin; j++)
		{
			factor = -operator()(c,j);
			for (i=c; i<li; i++)
				operator()(i,j)+=factor*operator()(i,c);
		}
	}
	for (c=lmin-1; c>0; c--) // Limpiar fila c para que quede triangular
	{
		for (j=0; j<c; j++)
		{
			factor = -operator()(c,j);
			for (i=0; i<li; i++)
				operator()(i,j)+=factor*operator()(i,c);
		}
	}
	return true;
}

bool L_Matrix::solveLinearSystem_x_Ab_result_in_b(L_Matrix &A_Mod)
{
	int i, j, c;
	double factor;
	double temp;
	throw_L_ArgException_if(A_Mod.li!=A_Mod.lj || li!=A_Mod.li || lj != 1, "L_Matrix::solveLinearSystem_x_Ab_result_in_b");
	for (c=0; c<A_Mod.lj; c++) // Limpiar columna c desde i=c hacia abajo
	{
		if (A_Mod(c,c)==0) // Diagonal con zero, intercambiar por other fila mas abajo que c
		{
			for (i=c+1; i<A_Mod.li; i++) // Encontrar fila i con elemento != 0 en la columna c
				if (A_Mod(i,c)!=0)
					break;
			if (i == A_Mod.li) // fila i solo contiene elementos zero para j>=c
				return false;
			for (j=c; j<A_Mod.lj; j++) // intercambiar fila c por fila i. Notar que operator()(i,j) = 0 para j<c
			{
				temp = A_Mod(i,j);
				A_Mod(i,j)=A_Mod(c,j);
				A_Mod(c,j) = temp;
			}
			temp=operator()(i,0); // apply la misma operacion sobre el vector "b", que esta a la derecha de la matriz "A"
			operator()(i,0)=operator()(c,0);
			operator()(c,0)=temp;
		}
		// Amplificar fila c para que quede con 1 en la diagonal
		factor = 1 / A_Mod(c,c);
		for (j=c; j<A_Mod.lj; j++)
			A_Mod(c,j)*=factor;
		operator()(c,0)*=factor; // apply la misma operacion sobre el vector "b", que esta a la derecha de la matriz "A"
		// Pivotear todas las filas para dejar la columna c triangular superior
		for (i=c+1; i<A_Mod.lj; i++)
		{
			factor = -A_Mod(i,c);
			for (j=c; j<A_Mod.lj; j++)
			{
				A_Mod(i,j)+=factor*A_Mod(c,j);
			}
			operator()(i,0)+=factor*operator()(c,0); // apply la misma operacion sobre el vector "b", que esta a la derecha de la matriz "A"
		}
	}
	// Hacer sustituciones para despejar x. operator()(u,0) corresponde a "b" para u<i, "x" para u>i
	for (i=A_Mod.li-1; i>=0; i--)
	{
		temp = 0;
		for (j=i+1; j<A_Mod.lj; j++)
			temp+=A_Mod(i,j)*operator()(j,0);
		operator()(i,0) = operator()(i,0) - temp; // Supuestamente A_Mod(i,i)=1
	}
	return true;
}


bool L_Matrix::svd(L_Matrix &d_vector, L_Matrix &v)
{
	return svd_0(d_vector, v);
}

bool L_Matrix::svd_0(L_Matrix &d_vector, L_Matrix &v)
{
	//void L_HomogeneousMatrix::svdcmp(double **a, int m, int n, double w[], double **v)
	// Given a matrix a[0..m-1][0__n-1], m>=n, this routine computes its singular value decomposition, A = U·W·VT.
	// The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[0__n-1].
	// The matrix V (not the transpose VT ) is output as v[0__n-1][0__n-1].
	
	if (li < lj)
	{
		L_Matrix vT, uT;
		vT.transpOf(*this);
		if (vT.svd(d_vector,uT) == false)  // thisT = vT*d*uT
			return false;
		v.transpOf(vT);
		(*this).transpOf(uT);
		return true;
		
	}
	int i;
	std::vector<double> w(lj+1);
	bool ret;
	v.reallocate(lj, lj);
	double **pelem, **velem;
	pelem = L_new2d_ext<double>(li, lj, begin(), lj*sizeof(double));
	velem = L_new2d_ext<double>(v.li, v.lj, v.begin(), v.lj*sizeof(double));
	ret=svdcmp_0(pelem, li, lj, li, lj, &(w[0]), velem);
	L_delete2d_ext<double>(pelem);
	L_delete2d_ext<double>(velem);
	d_vector.reallocate(lj, 1);
	for (i=0; i<lj; i++)
		d_vector(i,0)=w[i];
	return ret;
}

// Dada una matriz P, la inversa es como la traspuesta dividida por el determinante, es facil...
bool L_Matrix::vectPr_2x2(L_Matrix &valPr, bool columna)
{
	throw_L_ArgException_if( li != 2 ||  lj != 2 || valPr.li!= 2 || begin() == NULL, "L_Matrix::vectPr_2x2() : argumentos no validos");
	double a, b, c, d, discr, det, tr, val1, val2;
	L_StaticMatrix<2,2> m2;
	a = operator()(0,0);
	b = operator()(0,1);
	c = operator()(1,0);
	d = operator()(1,1);
	det = a*d-b*c;
	tr = a+d;
	discr = tr*tr/4-det;
	if (discr < 0)  // Valores propios complejos
		return false;
	discr = sqrt(discr);
	val1 = tr/2 + discr;
	val2 = tr/2 - discr;
	valPr(0,0) = val1;
	valPr(1,columna?0:1) = val2;
	// Calcular condicionamiento de los calculos, a veces dan resultados horribles
	double ca, cb, cc, cd;
	double v1,dd,v2,aa;
	bool malCond, calcm1=false, calcm2=false;
	// Calcular condicionamiento de las matrices
	aa = L_ABS(a);
	dd = L_ABS(d);
	v1 = L_ABS(val1);
	v2 = L_ABS(val2);
	ca = fabs(val1-d) / (v1+dd);
	cb = fabs(val2-d) / (v2+dd);
	cc = fabs(val1-a) / (v1+aa);
	cd = fabs(val2-a) / (v1+aa);
	ca = L_MIN(ca,cb);
	cc = L_MIN(cc,cd);
	malCond = (ca < 1e-8 && cc < 1e-8); // Mal condicionada, hay que evaluar todas las opciones y elegir la mejor

	if (c!=0 && (ca >= cc || b==0 || malCond))
	{
		operator()(0,0) = val1-d;
		operator()(1,0) = c;
		operator()(0,1) = val2-d;
		operator()(1,1) = c;
		calcm1 = true;
	}
	if (b!=0 && (ca < cc || c==0 || malCond))
	{
		m2(0,0) = b;
		m2(1,0) = val1-a;
		m2(0,1) = b;
		m2(1,1) = val2-a;
		calcm2 = true;
	}
	if (calcm1 == false && calcm2 == false)
	{
		// La matriz es diagonal, asi que no incluye rotaciones  =>  P = I
		operator()(0,0) = 1;
		operator()(0,1) = 0;
		operator()(1,0) = 0;
		operator()(1,1) = 1;
		return true;
	}
	if (calcm2 == true && calcm1 == false)
	{
		calcm1 = true;
		calcm2 = false;
		operator()(0,0) = m2(0,0);
		operator()(0,1) = m2(0,1);
		operator()(1,0) = m2(1,0);
		operator()(1,1) = m2(1,1);
	}
	// Arreglar la matriz de vectores propios para que sea unitaria (en lo posible)
	double r1 = sqrt(operator()(0,0)*operator()(0,0)+operator()(1,0)*operator()(1,0));
	double r2 = sqrt(operator()(0,1)*operator()(0,1)+operator()(1,1)*operator()(1,1));
	det = operator()(0,0) * operator()(1,1) - operator()(1,0) * operator()(0,1);
	if (det < 0)  // Para que dejarlo con det negativo...
		r2 = -r2;
	if (operator()(0,0) < 0)
	{
		r1 = -r1;
		r2 = -r2;
	}
	if (r1 != 0)
	{
		operator()(0,0) /= r1;
		operator()(1,0) /= r1;
	}

	if (r2 != 0)
	{
		operator()(0,1) /= r2;
		operator()(1,1) /= r2;
	}
	if (calcm2 == true) // Matriz mal condicionada, hay que ver las 2 opciones, es lento pero ocurre poco
	{
		// Arreglar la matriz de vectores propios para que sea unitaria (en lo posible)
		r1 = sqrt(m2(0,0)*m2(0,0)+m2(1,0)*m2(1,0));
		r2 = sqrt(m2(0,1)*m2(0,1)+m2(1,1)*m2(1,1));
		det = m2(0,0) * m2(1,1) - m2(1,0) * m2(0,1);
		if (det < 0)  // Para que dejarlo con det negativo...
			r2 = -r2;
		if (m2(0,0) < 0)
		{
			r1 = -r1;
			r2 = -r2;
		}
		if (r1 != 0)
		{
			m2(0,0) /= r1;
			m2(1,0) /= r1;
		}

		if (r2 != 0)
		{
			m2(0,1) /= r2;
			m2(1,1) /= r2;
		}
		double sima, simb, simc, simd;
		sima = operator()(0,0) - operator()(1,1);
		simb = operator()(0,1) + operator()(1,0);
		simc = m2(0,0) - m2(1,1);
		simd = m2(0,1) + m2(1,0);
		sima = L_ABS(sima);
		simb = L_ABS(simb);
		simc = L_ABS(simc);
		simd = L_ABS(simd);
		sima = (sima<simb) ? sima : simb;
		simc = (simc<simd) ? simc : simd;
		// Elegir la mejor simetria
		if (simc < sima) // Es mejor la matriz m2
		{
			operator()(0,0) = m2(0,0);
			operator()(0,1) = m2(0,1);
			operator()(1,0) = m2(1,0);
			operator()(1,1) = m2(1,1);
		}
		if (fabs(sima) > 0.1 || fabs(simb) > 0.1) // Matriz resultante horrible, asimetrica
			return false; // Se hizo lo que se pudo
	}
	return true;
}

bool L_Matrix::svd2x2(double &ang1, double &ang2, double &e1, double &e2) const
{
	return svd2x2_lento_a_mano(ang1, ang2, e1, e2);
	return true;
}

bool L_Matrix::svd2x2(double &rot, double &eMax, double &eMin) const
{
	L_Matrix U(2,2),D(2,1),V(2,2), VT(2,2), R(2,2);

	U(0,0) = operator()(0,0);
	U(0,1) = operator()(0,1);
	U(1,0) = operator()(1,0);
	U(1,1) = operator()(1,1);

	if (U.svd_ordenado(D,V) == false)
		printf("T_T\n");
	VT.transpOf(V);
	eMax = D(0,0);
	eMin = D(1,0);
	R.OP_mult(U,VT);
	rot = atan2(R(1,0), R(0,0));
	return true;
}

bool L_Matrix::svd2x2_lento_a_mano(double &ang1, double &ang2, double &e2, double &e1, bool debuggg) const
{
	// Wikipedia:
	// MMT = (UDVT)(VDTUT) -> MMT = U * D*DT * UT
	// MTM = (VDTUT)(UDVT) -> MTM = V * DT*D * VT
	const L_Matrix &M = *this;
	L_Matrix MT(2,2);
	L_Matrix U(2,2), D(2,1), V(2,2), VT(2,2), mtmp(2,2);
	double c, s;
	if (operator()(0,0) == operator()(1,1) && operator()(0,1) == -operator()(1,0)) // Las 2 escalas son iguales
	{
		e1 = e2 = sqrt(operator()(0,0) * operator()(1,1) - operator()(0,1) * operator()(1,0));
		ang1 = atan2(operator()(1,0), operator()(0,0));
		ang2 = 0;
		return true;
	}
	MT.transpOf(M);
	L_StaticMatrix_OP_mult2x2(U,M,MT);
	U.vectPr_2x2(D);  // Hace su mejor esfuerzo
	D(0,0)=sqrt(D(0,0));
	D(1,0)=sqrt(D(1,0));

	// Mejorar un poco la estabilidad numerica
	c = (U(0,0) + U(1,1))/2;
	s = (U(1,0) - U(0,1))/2;
	U(0,0) = c;
	U(0,1) = -s;
	U(1,0) = s;
	U(1,1) = c;
	e1 = D(0,0);
	e2 = D(1,0);

	// U = [c s ; -s c)
	ang1 = atan2(U(1,0), U(0,0));

	L_StaticMatrix_OP_mult2x2(V,MT,M);
	V.vectPr_2x2(D); // Hace su mejor esfuerzo
	VT.transpOf(V);
	D(0,0) = sqrt(D(0,0));
	D(1,0) = sqrt(D(1,0));

	// Mejorar un poco la estabilidad numerica
	c = (VT(0,0) + VT(1,1))/2;
	s = (VT(1,0) - VT(0,1))/2;
	VT(0,0) = c;
	VT(0,1) = -s;
	VT(1,0) = s;
	VT(1,1) = c;	// VT = [c s ; -s c)
	ang2 = atan2(VT(1,0), VT(0,0));

	// A ver si el angulo ang2 es el adecuado
	U(0,0) = cos(ang1);
	U(0,1) = -sin(ang1);
	U(1,0) = -U(0,1);
	U(1,1) = U(0,0);
	D(0,0) = e1;
	D(1,0) = e2;
	VT(0,0) = cos(ang2);
	VT(0,1) = -sin(ang2);
	VT(1,0) = -VT(0,1);
	VT(1,1) = VT(0,0);

	mtmp.OP_multByDiagonal(U,D);
	L_StaticMatrix_OP_mult2x2(MT,mtmp,VT);
	double punto = MT(0,0) * operator()(0,0) + MT(0,1) * operator()(0,1) + MT(1,0) * operator()(1,0) + MT(1,1) * operator()(1,1);
	double m2 = operator()(0,0)*operator()(0,0) + operator()(0,1)*operator()(0,1) + operator()(1,0)*operator()(1,0) + operator()(1,1)*operator()(1,1);
	double rel = punto / m2;
	if (debuggg)
	{
		printf("%f  %f\n", punto, m2);
		print("M");
		MT.print("m");
	}
	if (rel < 0)
		rel = -rel;
	if (punto < 0)
		ang2 = ang2 + M_PI;
	while (ang1 > M_PI)
		ang1 -= 2*M_PI;
	while (ang1 < -M_PI)
		ang1 += 2*M_PI;
	while (ang2 > M_PI)
		ang2 -= 2*M_PI;
	while (ang2 < -M_PI)
		ang2 += 2*M_PI;
	if (rel > 0.99999999 && rel < 1.000000001 && rel==rel)
		return true;
	return false;
}

bool L_Matrix::cholesky_L_De(const L_Matrix &original)
{
	throw_L_ArgException_if (original.li!=original.lj, "L_Matrix::cholesky_L_De() : matriz no cuadrada");

	reallocate(original.li, original.lj);

	for( int i=0; i<li; i++ )
	{
		const double *fila_i_this = &(operator()(i,0));
		const double *fila_i_original = &(original(i,0));

		for( int j=i; j<li; j++ )
		{
			double *fila_j_this = &operator()(j,0);
			double sum = fila_i_original[j];
			for( int k=i-1; k>=0; k-- )
				sum -= fila_i_this[k] * fila_j_this[k];
			if( i == j )
			{
				if( sum <= 0 )
					return false;

				if( begin() != NULL )
					operator()(j,i) = sqrt( sum );

			}
			else
			{
				if( begin() != NULL )
					operator()(j,i) = sum / operator()(i,i); 
			}		
			if( i<j )
				operator()(i,j) = 0;
		}
	}
	return true;
}

void L_Matrix::cholesky_L_forzado_De(const L_Matrix &original)
{
	throw_L_ArgException_if (original.li!=original.lj, "L_Matrix::cholesky_L_forzado_De() : matriz no cuadrada");

	reallocate(original.li, original.lj);

	for( int i=0; i<li; i++ )
	{
		double *fila_i_this = &(operator()(i,0));
		const double *fila_i_original = &(original(i,0));

		for( int j=i; j<li; j++ )
		{
			double *fila_j_this = &operator()(j,0);
			double sum = fila_i_original[j];
			for( int k=i-1; k>=0; k-- )
				sum -= fila_i_this[k] * fila_j_this[k];
			if( i == j )
			{
				if( sum <= 0 )
					sum = 1e-30;

				if( begin() != NULL )
					operator()(j,i) = sqrt( sum );

			}
			else
			{
				if( begin() != NULL )
					operator()(j,i) = sum / operator()(i,i); 
			}		
			if( i<j )
				operator()(i,j) = 0;
		}
	}
	return;
}

bool L_Matrix::unscented_desde_media_cov(const L_Matrix &media, const L_Matrix &cov, double alfa, double kappa, double beta)
{
	L_Matrix cov2, root;
	int i, j;
	double lambda = alfa*alfa*(cov.li+kappa) - cov.li;
	throw_L_ArgException_if(cov.li != cov.lj,  "L_Matrix::unscented_desde_media_cov");
	reallocate(cov.li, 2*cov.lj + 1);

	// "cov2 = (cov.li + lambda)*cov;"
	cov2.OP_mult(cov, cov.li+lambda);

	if (root.cholesky_L_De(cov2) == false) // La root cuadrada si es def pos
		return false;

	for (i=0; i<cov.li; i++)
		operator()(i,0) = media(i,0); // Punto central
	for (j=0; j<cov.lj; j++)
		for (i=0; i<cov.li; i++)
			operator()(i,j+1) = media(i,0) + root(i,j); // Puntos positivos
	for (j=0; j<cov.lj; j++)
		for (i=0; i<cov.li; i++)
			operator()(i,cov.lj+j+1) = media(i,0) - root(i,j); // puntos negativos
	return true;
}

bool L_Matrix::unscented_hacia_media_cov(L_Matrix &fmedia, L_Matrix &fcov, double alfa, double kappa, double beta)
{
	double lambda = alfa*alfa*(fcov.li+kappa) - fcov.li;
	double ws0, wc0, wi, dx, dxT;
	int i, j, k;
	throw_L_ArgException_if(2*li+1 != lj, "L_Matrix::unscented_hacia_media_cov");
	fcov.reallocate(li, li);

	ws0 = lambda / (fcov.li + lambda);
	wc0 = lambda / (fcov.li + lambda) + (1 - alfa*alfa + beta);
	wi =  0.5 / (fcov.li + lambda);

	fmedia.setZero();
	fcov.setZero();

	for (i=0; i<fcov.li; i++)
		fmedia(i,0) = ws0*operator()(i,0);
	for (k=1; k<lj; k++)
		for (i=0; i<fcov.li; i++)
			fmedia(i,0) += wi*operator()(i,k);

	for (i=0; i<fcov.li; i++)
	{
		for (j=0; j<fcov.li; j++)
		{
			dx = operator()(i,0) - fmedia(i,0);
			dxT = operator()(j,0) - fmedia(j,0); // (i,j) intercambiados
			fcov(i,j) = wc0 * dx * dxT;
		}
	}

	for (k=1; k<lj; k++)
	{
		for (i=0; i<fcov.li; i++)
		{
			for (j=0; j<fcov.li; j++)
			{
				dx = operator()(i,k) - fmedia(i,0);
				dxT = operator()(j,k) - fmedia(j,0); // (i,j) intercambiados
				fcov(i,j) += wi * dx * dxT;
			}
		}
	}
	return true;
}


bool L_Matrix::svdcmp_0( double **A, int M, int N,
    int MP, int NP, double *W, double **V )
{
    /*
       Give a matrix A, with logical dimensions M by N and physical
       dimensions MP by NP, this routine computes its singular value
       decomposition, A = U * W * transpose V. The matrix U replaces
       A on output. The diagonal matrix of singular values, W, is output
       as a vector W. The matrix V (not the transpose of V) is output as
       V. M must be greater or equal to N. If it is smaller then A should
       be filled up to square with zero rows.
    */

	std::vector<double> rv1(N);
	//double rv1[4000]; // OJO CON ESTA LINEA

    /* Householder reduction to bidiagonal form. */
    int NM=0;
    double C;
    double F;
    double G = 0.0;
    double H;
    double S;
    double X;
    double Y;
    double Z;
    double Scale = 0.0;
    double ANorm = 0.0;
    double tmp;
    int flag;
    int i;
    int its;
    int j;
    int jj;
    int k;
    int l=0;

    if( M < N ) {
        fprintf( stderr, "You must augment A with extra zero rows.\n" );
        return false;
    }

    for( i = 0; i < N; ++i ) {
        l = i + 1;
        rv1[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( i < M ) {
            for( k = i; k < M; ++k ) {
                Scale = Scale + fabs( A[k][i] );
            }
            if( Scale != 0.0 ) {
                for( k = i; k < M; ++k ) {
                    A[k][i] = A[k][i] / Scale;
                    S = S + A[k][i] * A[k][i];
                }
                F = A[i][i];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][i] = F - G;
                if( i != (N-1) ) {
                    for( j = l; j < N; ++j ) {
                        S = 0.0;
                        for( k = i; k < M; ++k ) {
                            S = S + A[k][i] * A[k][j];
                        }
                        F = S / H;
                        for( k = i; k < M; ++k ) {
                            A[k][j] = A[k][j] + F * A[k][i];
                        }
                    }
                }
                for( k = i; k < M; ++k ) {
                    A[k][i] = Scale * A[k][i];
                }
            }
        }

        W[i] = Scale * G;
        G = 0.0;
        S = 0.0;
        Scale = 0.0;
        if( (i < M) && (i != (N-1)) ) {
            for( k = l; k < N; ++k ) {
                Scale = Scale + fabs( A[i][k] );
            }
            if( Scale != 0.0 ) {
                for( k = l; k < N; ++k ) {
                    A[i][k] = A[i][k] / Scale;
                    S = S + A[i][k] * A[i][k];
                }
                F = A[i][l];
                G = sqrt(S);
                if( F > 0.0 ) {
                    G = -G;
                }
                H = F * G - S;
                A[i][l] = F - G;
                for( k = l; k < N; ++k ) {
                    rv1[k] = A[i][k] / H;
                }
                if( i != (M-1) ) {
                    for( j = l; j < M; ++j ) {
                        S = 0.0;
                        for( k = l; k < N; ++k ) {
                            S = S + A[j][k] * A[i][k];
                        }
                        for( k = l; k < N; ++k ) {
                            A[j][k] = A[j][k] + S * rv1[k];
                        }
                    }
                }
                for( k = l; k < N; ++k ) {
                    A[i][k] = Scale * A[i][k];
                }
            }
        }
        tmp = fabs( W[i] ) + fabs( rv1[i] );
        if( tmp > ANorm )
            ANorm = tmp;
    }

    /* Accumulation of right-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        if( i < (N-1) ) {
            if( G != 0.0 ) {
                for( j = l; j < N; ++j ) {
                    V[j][i] = (A[i][j] / A[i][l]) / G;
                }
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < N; ++k ) {
                        S = S + A[i][k] * V[k][j];
                    }
                    for( k = l; k < N; ++k ) {
                        V[k][j] = V[k][j] + S * V[k][i];
                    }
                }
            }
            for( j = l; j < N; ++j ) {
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        V[i][i] = 1.0;
        G = rv1[i];
        l = i;
    }

    /* Accumulation of left-hand transformations. */
    for( i = N-1; i >= 0; --i ) {
        l = i + 1;
        G = W[i];
        if( i < (N-1) ) {
            for( j = l; j < N; ++j ) {
                A[i][j] = 0.0;
            }
        }
        if( G != 0.0 ) {
            G = 1.0 / G;
            if( i != (N-1) ) {
                for( j = l; j < N; ++j ) {
                    S = 0.0;
                    for( k = l; k < M; ++k ) {
                        S = S + A[k][i] * A[k][j];
                    }
                    F = (S / A[i][i]) * G;
                    for( k = i; k < M; ++k ) {
                        A[k][j] = A[k][j] + F * A[k][i];
                    }
                }
            }
            for( j = i; j < M; ++j ) {
                A[j][i] = A[j][i] * G;
            }
        } else {
            for( j = i; j < M; ++j ) {
                A[j][i] = 0.0;
            }
        }
        A[i][i] = A[i][i] + 1.0;
    }

    /* Diagonalization of the bidiagonal form.
       Loop over singular values. */
    for( k = (N-1); k >= 0; --k ) {
        /* Loop over allowed iterations. */
        for( its = 1; its <= 30; ++its ) {
            /* Test for splitting.
               Note that rv1[0] is always zero. */
            flag = true;
            for( l = k; l >= 0; --l ) {
                NM = l - 1;
                if( (fabs(rv1[l]) + ANorm) == ANorm ) {
                    flag = false;
                    break;
                } else if( (fabs(W[NM]) + ANorm) == ANorm ) {
                    break;
                }
            }

            /* Cancellation of rv1[l], if l > 0; */
            if( flag ) {
                C = 0.0;
                S = 1.0;
                for( i = l; i <= k; ++i ) {
                    F = S * rv1[i];
                    if( (fabs(F) + ANorm) != ANorm ) {
                        G = W[i];
                        H = sqrt( F * F + G * G );
                        W[i] = H;
                        H = 1.0 / H;
                        C = ( G * H );
                        S = -( F * H );
                        for( j = 0; j < M; ++j ) {
                            Y = A[j][NM];
                            Z = A[j][i];
                            A[j][NM] = (Y * C) + (Z * S);
                            A[j][i] = -(Y * S) + (Z * C);
                        }
                    }
                }
            }
            Z = W[k];
            /* Convergence. */
            if( l == k ) {
                /* Singular value is made nonnegative. */
                if( Z < 0.0 ) {
                    W[k] = -Z;
                    for( j = 0; j < N; ++j ) {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }

            if( its >= 30 ) {
				fprintf( stderr, "svdcmp_0() : No convergence in 30 iterations.\n" );
                return false;
            }

            X = W[l];
            NM = k - 1;
            Y = W[NM];
            G = rv1[NM];
            H = rv1[k];
            F = ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / (2.0*H*Y);
            G = sqrt( F * F + 1.0 );
            tmp = G;
            if( F < 0.0 )
                tmp = -tmp;
            F = ((X-Z)*(X+Z) + H*((Y/(F+tmp))-H)) / X;

            /* Next QR transformation. */
            C = 1.0;
            S = 1.0;
            for( j = l; j <= NM; ++j ) {
                i = j + 1;
                G = rv1[i];
                Y = W[i];
                H = S * G;
                G = C * G;
                Z = sqrt( F * F + H * H );
                rv1[j] = Z;
                C = F / Z;
                S = H / Z;
                F = (X * C) + (G * S);
                G = -(X * S) + (G * C);
                H = Y * S;
                Y = Y * C;
                for( jj = 0; jj < N; ++jj ) {
                    X = V[jj][j];
                    Z = V[jj][i];
                    V[jj][j] = (X * C) + (Z * S);
                    V[jj][i] = -(X * S) + (Z * C);
                }
                Z = sqrt( F * F + H * H );
                W[j] = Z;

                /* Rotation can be arbitrary if Z = 0. */
                if( Z != 0.0 ) {
                    Z = 1.0 / Z;
                    C = F * Z;
                    S = H * Z;
                }
                F = (C * G) + (S * Y);
                X = -(S * G) + (C * Y);
                for( jj = 0; jj < M; ++jj ) {
                    Y = A[jj][j];
                    Z = A[jj][i];
                    A[jj][j] = (Y * C) + (Z * S);
                    A[jj][i] = -(Y * S) + (Z * C);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = F;
            W[k] = X;
        }
    }
    return true;
}



void L_Matrix::Jacobi_Cyclic_Method(double *valpr, double *vectpr, double *A, int n)
{
	int i, j, k, m;
	double *pAk, *pAm, *p_r, *p_e;
	double threshold_norm;
	double threshold;
	double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
	double sin_2phi, cos_2phi, cot_2phi;
	double dum1;
	double dum2;
	double dum3;
	double max;

  // Take care of trivial cases
	if ( n < 1) return;
	if ( n == 1) {
		valpr[0] = *A;
		*vectpr = 1.0;
		return;
	}

  // Initialize the valpr to the identity matrix.
	for (p_e = vectpr, i = 0; i < n; i++)
		for (j = 0; j < n; p_e++, j++)
			if (i == j) *p_e = 1.0; else *p_e = 0.0;
  
	// Calculate the threshold and threshold_norm.
	for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
		for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
	threshold = sqrt(threshold + threshold);
	threshold_norm = threshold * DBL_EPSILON;
	max = threshold + 1.0;
	while (threshold > threshold_norm)
	{
		threshold /= 10.0;
		if (max < threshold) continue;
		max = 0.0;
		for (pAk = A, k = 0; k < (n-1); pAk += n, k++)
		{
			for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++)
			{
				if ( fabs(*(pAk + m)) < threshold )
					continue;

                 // Calculate the sin and cos of the rotation angle which
                 // annihilates A[k][m].

				cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
				dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
				if (cot_2phi < 0.0) dum1 = -dum1;
				tan_phi = -cot_2phi + dum1;
				tan2_phi = tan_phi * tan_phi;
				sin2_phi = tan2_phi / (1.0 + tan2_phi);
				cos2_phi = 1.0 - sin2_phi;
				sin_phi = sqrt(sin2_phi);
				if (tan_phi < 0.0) sin_phi = - sin_phi;
				cos_phi = sqrt(cos2_phi); 
				sin_2phi = 2.0 * sin_phi * cos_phi;
				cos_2phi = cos2_phi - sin2_phi;

				// Rotate columns k and m for both the matrix A 
				//     and the matrix of vectpr.
				p_r = A;
				dum1 = *(pAk + k);
				dum2 = *(pAm + m);
				dum3 = *(pAk + m);
				*(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
				*(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
				*(pAk + m) = 0.0;
				*(pAm + k) = 0.0;
				for (i = 0; i < n; p_r += n, i++)
				{
					if ( (i == k) || (i == m) )
						continue;
					if ( i < k )
						dum1 = *(p_r + k);
					else
						dum1 = *(pAk + i);
					if ( i < m )
						dum2 = *(p_r + m);
					else
						dum2 = *(pAm + i);
					dum3 = dum1 * cos_phi + dum2 * sin_phi;
					if ( i < k )
						*(p_r + k) = dum3;
					else
						*(pAk + i) = dum3;
					dum3 = - dum1 * sin_phi + dum2 * cos_phi;
					if ( i < m )
						*(p_r + m) = dum3;
					else
						*(pAm + i) = dum3;
				}
				for (p_e = vectpr, i = 0; i < n; p_e += n, i++)
				{
					dum1 = *(p_e + k);
					dum2 = *(p_e + m);
					*(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
					*(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
				}
			}
			for (i = 0; i < n; i++)
				if ( i == k )
					continue;
			else if ( max < fabs(*(pAk + i)))
				max = fabs(*(pAk + i));
		}
	}
	for (pAk = A, k = 0; k < n; pAk += n, k++)
		valpr[k] = *(pAk + k); 
}




// http://www.netlib.org/blas/zrotg.f
void L_MatrizFortranTr::zrotg(double &CA, double &CB, double &C, double &S) // SUBROUTINE ZROTG(CA,CB,C,S)
{
	// DOUBLE COMPLEX CA,CB,S
	// DOUBLE PRECISION C
	//
	// *  Purpose
	// *  =======
	// *
	// *    determines a double complex Givens rotation.

	double ALPHA; //      DOUBLE COMPLEX ALPHA
	double NORM, SCALE; //      DOUBLE PRECISION NORM,SCALE
	// *    ..
	// *    .. Intrinsic Functions ..
	//      INTRINSIC CDABS,DCMPLX,DCONJG,DSQRT
	// *    ..
	if (CA == 0.0) // IF (CDABS(CA).NE.0.0d0) GO TO 10
	{
		C = 0.0;
		S = 1.0;
		CA = CB;
		//GO TO 20
	}//10 CONTINUE
	else
	{
		SCALE = fabs(CA) + fabs(CB); // SCALE = CDABS(CA) + CDABS(CB)
		double casc = CA/SCALE;
		double sbsc = CB/SCALE;
		NORM = SCALE * sqrt(casc*casc + sbsc*sbsc);// NORM = SCALE*DSQRT((CDABS(CA/DCMPLX(SCALE,0.0d0)))**2+(CDABS(CB/DCMPLX(SCALE,0.0d0)))**2)
		ALPHA = CA/fabs(CA); // ALPHA = CA/CDABS(CA)
		C = fabs(CA)/NORM; // CDABS(CA)/NORM
		S = ALPHA*CB/NORM; // ALPHA*DCONJG(CB)/NORM
		CA = ALPHA*NORM; //
	}  // 20 CONTINUE
	return;   //RETURN
}  //    END

// http://www.netlib.org/blas/dnrm2.f
double L_MatrizFortranTr::dnrm2(int N, L_MatrizFortranTr X, int INCX)
{
//       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
// *     .. Scalar Arguments ..
//       INTEGER INCX,N
// *     ..
// *     .. Array Arguments ..
//       DOUBLE PRECISION X(*)
// *     ..
// *
// *  Purpose
// *  =======
// *
// *  DNRM2 returns the euclidean norm of a vector via the function
// *  name, so that
// *
// *     DNRM2 := sqrt( x'*x )
// *
// *
// *  -- This version written on 25-October-1982.
// *     Modified on 14-October-1993 to inline the call to DLASSQ.
// *     Sven Hammarling, Nag Ltd.
// *
// *
// *     .. Parameters ..
//       DOUBLE PRECISION ONE,ZERO
//       PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
// *     ..
// *     .. Local Scalars ..
//       DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
//       INTEGER IX
// *     ..
// *     .. Intrinsic Functions ..
//       INTRINSIC ABS,SQRT
//
	double NORM,SCALE,SSQ,TEMP;
    int IX;

	if (N < 1 || INCX < 1)  //IF (N.LT.1 .OR. INCX.LT.1) THEN
		NORM = 0;
	else //  ELSE
	{
		SCALE = 0;
		SSQ = 1;
		//        The following loop is equivalent to this call to the LAPACK
		//        auxiliary routine:
		//        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
		//
		for (IX=1; IX<=1+(N-1)*INCX; IX+=INCX)//DO 10 IX = 1,1 + (N-1)*INCX,INCX
		{
			if (X(IX) != 0) //IF (DBLE(X(IX)).NE.ZERO) THEN
			{
				TEMP = fabs(X(IX));//TEMP = ABS(DBLE(X(IX)))
				if (SCALE < TEMP) //IF (SCALE.LT.TEMP) THEN
				{
					SSQ = 1 + SSQ* (SCALE/TEMP)*(SCALE/TEMP);
					SCALE = TEMP;
				}
				else //ELSE
					SSQ = SSQ + (TEMP/SCALE)*(TEMP/SCALE);
				//END IF
			}// END IF
			// IF (DIMAG(X(IX)).NE.ZERO) THEN ... ENDIF -> aca no hay imaginarios
		} // 10     CONTINUE
		NORM = SCALE*sqrt(SSQ);
	}//END IF

	return NORM;
}

// http://www.netlib.org/blas/zdotc.f
double L_MatrizFortranTr::zdotc(int N, L_MatrizFortranTr ZX, int INCX, L_MatrizFortranTr ZY, int INCY)
{
	int I, IX, IY;
	double ZTEMP = 0;
	if (N == 0)  // IF (N.LE.0) RETURN
		return 0;
	if (INCX < 0)
		IX = (-N+1)*INCX+1;
	if (INCY < 0)
		IY = (-N+1)*INCY+1;

	if (INCX != 1 || INCY != 1) // IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
	{
		for (I=1; I<=N; I++) // DO 10 I = 1,N
		{
			ZTEMP += ZX(IX) * ZY(IY);  // ZTEMP = ZTEMP + DCONJG(ZX(IX))*ZY(IY)
			IX += INCX;
			IY += INCY;
		} // 10 CONTINUE
	}
	else
	{
		for (I=1; I<=N; I++)  // 20   DO 30 I = 1,N
			ZTEMP += ZX(I)*ZY(I);  // ZTEMP = ZTEMP + DCONJG(ZX(I))*ZY(I)
	} // 30 CONTINUE
	return ZTEMP;  // ZDOTC = ZTEMP    RETURN      END
}

// http://www.netlib.org/blas/daxpy.f
void L_MatrizFortranTr::daxpy(int N, double DA, L_MatrizFortranTr DX, int INCX, L_MatrizFortranTr DY, int INCY)
{
	  int IX, IY, I, M, MP1;
      if (N == 0) return;
      if (DA == 0.0e0) return;
      if ( ! (INCX == 1 && INCY == 1) )// GO TO 20
	  {
// *
// *        code for unequal increments or equal increments
// *          not equal to 1
// *
      IX = 1;
      IY = 1;
      if (INCX < 0) IX = (-N+1)*INCX + 1;
      if (INCY < 0) IY = (-N+1)*INCY + 1;
      for (I=1; I<=N; I++) // DO 10 I = 1,N
	  {
          DY(IY) = DY(IY) + DA*DX(IX);
          IX = IX + INCX;
          IY = IY + INCY;
	  } // 10 CONTINUE
      return;
// *
// *        code for both increments equal to 1
// *
// *
// *        clean-up loop
// *
	  } // 20
	  M = N%4; // MOD(N,4)
      if (M != 0)// if (M.EQ.0) // GO TO 40
	  {
      for (I = 1; I<=M; I++)//DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I);
   //30 CONTINUE
      if (N < 4) return;
	  } //
	  // 40
	  MP1 = M + 1;
      for (I=MP1; I<=N; I+=4) // DO 50 I = MP1,N,4
	  {
          DY(I) = DY(I) + DA*DX(I);
          DY(I+1) = DY(I+1) + DA*DX(I+1);
          DY(I+2) = DY(I+2) + DA*DX(I+2);
          DY(I+3) = DY(I+3) + DA*DX(I+3);
	  } //50 CONTINUE
      return;
      //END
}


void L_MatrizFortranTr::dswap(int N, L_MatrizFortranTr DX, int INCX, L_MatrizFortranTr DY, int INCY)
{
// ***END PROLOGUE  DSWAP
      double DTEMP1, DTEMP2, DTEMP3;
	  int IX, IY, I, M, MP1, NS;
// ***FIRST EXECUTABLE STATEMENT  DSWAP
     if (N <= 0) return;

	 if (INCX != INCY || INCX < 0 || INCY < 0)
	 {
//       Code for unequal or nonpositive increments.

      IX = 1;
      IY = 1;
      if (INCX < 0) IX = (-N+1)*INCX + 1;
      if (INCY < 0) IY = (-N+1)*INCY + 1;
      for (I=1; I<=N; I++) //DO 10 I = 1,N
	  {
        DTEMP1 = DX(IX);
        DX(IX) = DY(IY);
        DY(IY) = DTEMP1;
        IX = IX + INCX;
        IY = IY + INCY;
	  } // 10 CONTINUE
      return;
	 }
	 else if (INCX == INCY && INCX == 1)
	 {

//       Code for both increments equal to 1.

//       Clean-up loop so remaining vector length is a multiple of 3.

      M=N%3; // 20  M = MOD(N,3)
      if (M != 0) // if (M .EQ. 0) GO TO 40
	  {
      for (I=1; I<=M; I++) // DO 30 I = 1,M
	  {
        DTEMP1 = DX(I);
        DX(I) = DY(I);
        DY(I) = DTEMP1;
	  } // 30 CONTINUE
      if (N < 3) return;
	  } // 40
	  MP1 = M + 1;
      for (I=MP1; I<=N; I+=3) //DO 50 I = MP1,N,3
	  {
        DTEMP1 = DX(I);
        DTEMP2 = DX(I+1);
        DTEMP3 = DX(I+2);
        DX(I) = DY(I);
        DX(I+1) = DY(I+1);
        DX(I+2) = DY(I+2);
        DY(I) = DTEMP1;
        DY(I+1) = DTEMP2;
        DY(I+2) = DTEMP3;
	  } //50 CONTINUE
      return;
	 }
	 else
	 {
//       Code for equal, positive, non-unit increments.

      NS = N*INCX; // 60
      for (I=1; I<=NS; I+=INCX) // DO 70 I = 1,NS,INCX
	  {
        DTEMP1 = DX(I);
        DX(I) = DY(I);
        DY(I) = DTEMP1;
	  } // 70 CONTINUE
	 }
	 return;
      //END
}

// http://www.netlib.org/linpack/zchud.f
void L_MatrizFortranTr::zchud(L_MatrizFortranTr r, int ldr, int p, L_MatrizFortranTr x, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr y, L_MatrizFortranTr rho, L_MatrizFortranTr c, L_MatrizFortranTr s)
{
	int i, j, jm1;
	double azeta, scale;
	double t, xj, zeta;  // complex*16

	// c   update r.

	for (j=1; j<=p; j++)//do 30 j = 1, p
	{
		xj = x(j);// xj = x(j)

		// c      apply the previous rotations.

		jm1 = j - 1;
		if (jm1 >= 1)  //   if (jm1 .lt. 1) go to 20
		{
			for (i=1; i<=jm1; i++)// do 10 i = 1, jm1
			{
				t = c(i)*r(i,j) + s(i)*xj; // t = c(i)*r(i,j) + s(i)*xj
				xj = c(i)*xj - s(i)*r(i,j); // xj = c(i)*xj - dconjg(s(i))*r(i,j), todo es real aca
				r(i,j) = t; // r(i,j) = t
			} //	10    continue
		} //20    continue

	// c      compute the next rotation.

		zrotg(r(j,j),xj,c(j),s(j));   // call zrotg(r(j,j),xj,c(j),s(j))
	} // 30 continue

	// c   if required, update z and rho.

	if (nz >= 1) //if (nz .lt. 1) go to 70
	{
		for (j=1; j<=nz; j++)//do 60 j = 1, nz
		{
			zeta = y(j); // zeta = y(j)
			for (i=1; i<=p; i++)// do 40 i = 1, p
			{
				t = c(i)*z(i,j) + s(i)*zeta; // t = c(i)*z(i,j) + s(i)*zeta
				zeta = c(i)*zeta - s(i)*z(i,j); // zeta = c(i)*zeta - dconjg(s(i))*z(i,j)
				z(i,j) = t;  // z(i,j) = t
			} // 40    continue
			azeta = fabs(zeta);  // azeta = cdabs(zeta)
			if (azeta != 0.0 && rho(j) >= 0)  // //if (azeta .eq. 0.0d0 .or. rho(j) .lt. 0.0d0) go to 50
			{
				scale = azeta + rho(j); // scale = azeta + rho(j)
				double azs = azeta/scale;
				double rhs = rho(j)/scale;
				rho(j) = scale*sqrt(azs*azs+rhs*rhs); // rho(j) = scale*dsqrt((azeta/scale)**2+(rho(j)/scale)**2)
			} // 50    continue
		}//	60 continue
	} // 70 continue
	return;
} // end


void L_MatrizFortranTr::zchdd(L_MatrizFortranTr r, int ldr, int p, L_MatrizFortranTr x, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr y, L_MatrizFortranTr rho, L_MatrizFortranTr c, L_MatrizFortranTr s, int &info)
{
	int i,ii,j, iter = 0;
	double a,alpha,azeta,norm;
	double t,zeta,b,xx; // complex*16
	double scale; // No se por que no estaba declarada
	//complex*16 zdumr,zdumi   //complex*16
	//dreal(zdumr) = zdumr
	//dimag(zdumi) = (0.0d0,-1.0d0)*zdumi

// c   solve the system ctrans(r)*a = x, placing the result
// c   in the array s.

	info = 0;
	s(1) = x(1) / r(1,1); //s(1) = dconjg(x(1))/dconjg(r(1,1))
	if (p >= 2) //if (p .lt. 2) go to 20
	{
		for (j=2; j<=p; j++)//do 10 j = 2, p
		{
			s(j) = x(j) - zdotc(j-1, r.ref(1,j), 1, s, 1); //s(j) = dconjg(x(j)) - zdotc(j-1,r(1,j),1,s,1)
			s(j) = s(j)/r(j,j); // s(j) = s(j)/dconjg(r(j,j))
		} //10 continue
	} //20 continue
	norm = dnrm2(p, s, 1); // dznrm2(p,s,1)

	if (norm >= 1.0e0)   // if (norm .lt. 1.0d0) go to 30
	{
		info = -1;
		return;  // go to 120
	}

	// 30 continue
	alpha = sqrt(1.0e0 - norm*norm);// dsqrt(1.0d0-norm**2)

	// c      determine the transformations.

	for (ii=1; ii<=p ; ii++)//do 40 ii = 1, p
	{
		i = p - ii + 1;
		scale = alpha + fabs(s(i)); // scale = alpha + cdabs(s(i))
		a = alpha/scale;
		b = s(i)/scale;// b = s(i)/scale
		norm = sqrt(a*a+b*b); // norm = dsqrt(a**2+dreal(b)**2+dimag(b)**2)
		c(i) = a/norm; // c(i) = a/norm
		s(i) = b/norm; //s(i) = dconjg(b)/norm
		alpha = scale*norm;
	} // 40    continue

	// c      apply the transformations to r.

	for (j=1; j<=p; j++)//do 60 j = 1, p
	{
		xx = 0.0; //(0.0d0,0.0d0)
		for (ii=1; ii<=j; ii++)  //do 50 ii = 1, j
		{
			i = j - ii + 1;
			t = c(i)*xx + s(i)*r(i,j); // t = c(i)*xx + s(i)*r(i,j)
			r(i,j) = c(i) * r(i,j) - s(i)*xx; // r(i,j) = c(i)*r(i,j) - dconjg(s(i))*xx   // Error: decia * en vez de -
			xx = t;
		}  // 50       continue
	} //60    continue

	// c      if required, downdate z and rho.

	if (nz >= 1) // if (nz .lt. 1) go to 110
	{
		for (j=1; j<=nz; j++) //do 100 j = 1, nz
		{
			zeta = y(j); //zeta = y(j)
			for (i=1; i<=p; i++)  // do 70 i = 1, p
			{
				z(i,j) = (z(i,j) - s(i)*zeta) / c(i); //z(i,j) = (z(i,j) - dconjg(s(i))*zeta)/c(i)  // Error: faltaba *zeta  // Error: decia c[i][0]
				zeta = c(i)*zeta - s(i)*z(i,j); // zeta = c(i)*zeta - s(i)*z(i,j)
			} // 70       continue
			azeta = fabs(zeta); //azeta = cdabs(zeta)
			if (azeta >= rho(j))// if (azeta .le. rho(j)) go to 80    // Error: puse <
			{
				// OJO ACA
				double adentro = 1.0e0-(azeta/rho(j))*(azeta/rho(j));
				if (adentro < 0)
					adentro = 0;
				rho(j) = rho(j)*sqrt(adentro); // rho(j) = rho(j)*dsqrt(1.0d0-(azeta/rho(j))**2)
			}
			else // go to 90
			{// 80       continue
				info = 1;
				rho(j) = -1.0e0; // rho(j) = -1.0d0
			}  // 90       continue
		} //100    continue
	}// 110  continue
	//  120 continue
	return;
	//      end
}



void L_MatrizFortranTr::zchex(L_MatrizFortranTr r, int ldr, int p, int k, int l, L_MatrizFortranTr z, int ldz, int nz, L_MatrizFortranTr c, L_MatrizFortranTr s, int job)
{
      int i,ii,il,iu,j,jj,km1,kp1,lmk,lm1;
      double t;  // complex*16 rjp1j, ...
      //double precision dreal,dimag;
      //complex*16 zdumr,zdumi
      //dreal(zdumr) = zdumr
      //dimag(zdumi) = (0.0d0,-1.0d0)*zdumi

//c   initialize

      km1 = k - 1;
      kp1 = k + 1;
      lmk = l - k;
      lm1 = l - 1;

//c   perform the appropriate task.

      
	  if (job == 0) // go to (10,130), job
	  {
		  //10 continue

//c   right circular shift.

   //c      reorder the columns.

         for (i=1; i<=l; i++)  //do 20 i = 1, l
		 {
            ii = l - i + 1;
			s(i) = r(ii,l); //s(i) = r(ii,l)
		 }//20    continue
         for (jj=k; jj<=lm1; jj++) //do 40 jj = k, lm1
		 {
            j = lm1 - jj + k;
            for (i=1; i<=j; i++)  //do 30 i = 1, j
				r(i,j+1) = r(i,j); //r(i,j+1) = r(i,j)
            //30       continue
			r(j+1,j+1) = 0; //r(j+1,j+1) = (0.0d0,0.0d0)
		 }//40    continue
         if (k != 1) //if (k .eq. 1) go to 60
		 {
            for (i=1; i<=km1; i++) //do 50 i = 1, km1
			{
               ii = l - i + 1;
			   r(i,k) = s(ii); //r(i,k) = s(ii)
		    } // 50       continue
		 } //60    continue

//c      calculate the rotations.

		 t = s(1); //t = s(1)
         for (i=1; i<=lmk; i++) //do 70 i = 1, lmk
		 {
			 zrotg(s(i+1), t, c(i), s(i)); //call zrotg(s(i+1),t,c(i),s(i))
			 t = s(i+1); //t = s(i+1)
		 } // 70    continue
		 r(k,k) = t;  //r(k,k) = t
         for (j=kp1; j<=p; j++) //do 90 j = kp1, p
		 {
            il = L_MAX(1,l-j+1); // max0(1,l-j+1)
            for (ii=il; ii<=lmk; ii++) //do 80 ii = il, lmk
			{
               i = l - ii;
			   t = c(ii)*r(i,j) + s(ii)*r(i+1,j); // t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
			   r(i+1,j) = c(ii)*r(i+1,j) - s(ii)*r(i,j); // r(i+1,j) = c(ii)*r(i+1,j) - dconjg(s(ii))*r(i,j)
			   r(i,j) = t; // r(i,j) = t
			} //80       continue
		 } //90    continue

//c      if required, apply the transformations to z.

         if (nz >= 1) //if (nz .lt. 1) go to 120
		 {
         for (j=1; j<=nz; j++) //do 110 j = 1, nz
		 {
            for (ii=1; ii<=lmk; ii++) //do 100 ii = 1, lmk
			{
               i = l - ii;
			   t = c(ii)*z(i,j) + s(ii)*z(i+1,j); // t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
			   z(i+1,j) = c(ii)*z(i+1,j) - s(ii)*z(i,j); //  z(i+1,j) = c(ii)*z(i+1,j) - dconjg(s(ii))*z(i,j)
			   z(i,j) = t; // z(i,j) = t
			}// 100       continue
		 } //110    continue
		 } // 120    continue
      //go to 260
  	  } // 130 continue
	  else  // job
	  {
//c   left circular shift
//c      reorder the columns

         for (i=1; i<=k; i++) // do 140 i = 1, k
		 {
            ii = lmk + i;
			s(ii) = r(i,k); // s(ii) = r(i,k)
  
		 } // 140    continue
         for (j=k; j<=lm1; j++) // do 160 j = k, lm1
		 {
            for (i=1; i<=j; i++) // do 150 i = 1, j
				r(i,j) = r(i,j+1); // r(i,j) = r(i,j+1)
			// 150       continue
            jj = j - km1;
			s(jj) = r(j+1,j+1);// s(jj) = r(j+1,j+1)
		 } // 160    continue
         for (i=1; i<=k; i++) // do 170 i = 1, k
		 {
            ii = lmk + i;
			r(i,l) = s(ii); // r(i,l) = s(ii)
		 } // 170    continue
         for (i=kp1; i<=l; i++) //do 180 i = kp1, l
		 {
			 r(i,l) = 0; //r(i,l) = (0.0d0,0.0d0)
		 } // 180    continue

//c      reduction loop.

         for (j=k; j<=p; j++) // do 220 j = k, p
		 {
            if (j != k) // if (j .eq. k) go to 200
			{

//c            apply the rotations.

               iu = L_MIN(j-1,l-1); // iu = min0(j-1,l-1)
               for (i=k; i<=iu; i++) // do 190 i = k, iu
			   {
                  ii = i - k + 1;
				  t = c(ii)*r(i,j) + s(ii)*r(i+1,j); // t = c(ii)*r(i,j) + s(ii)*r(i+1,j)
				  r(i+1,j) = c(ii)*r(i+1,j) - s(ii)*r(i,j); // r(i+1,j) = c(ii)*r(i+1,j) - dconjg(s(ii))*r(i,j)
				  r(i,j) = t; // r(i,j) = t
			   } // 190          continue
			} // 200       continue
            if (j < l) // if (j .ge. l) go to 210
			{
			   jj = j - k + 1;
			   t = s(jj); // t = s(jj)
			   zrotg(r(j,j),t,c(jj),s(jj)); // call zrotg(r(j,j),t,c(jj),s(jj))
			} // 210       continue
		 } // 220    continue

  //c      apply the rotations to z.

         if (nz > 1) // if (nz .lt. 1) go to 250
		 {
         for (j=1; j<=nz; j++) // do 240 j = 1, nz
		 {
            for (i=k; i<=lm1; i++) // do 230 i = k, lm1
			{
               ii = i - km1;
			   t = c(ii)*z(i,j) + s(ii)*z(i+1,j); // t = c(ii)*z(i,j) + s(ii)*z(i+1,j)
			   z(i+1,j) = c(ii)*z(i+1,j) - s(ii)*z(i,j); // z(i+1,j) = c(ii)*z(i+1,j) - dconjg(s(ii))*z(i,j)
			   z(i,j) = t; // z(i,j) = t
			} // 230       continue
		 } // 240    continue
		 } // 250    continue
      } // 260 continue, job
      return;
      //end
}

void L_MatrizFortranTr::dchdc(L_MatrizFortranTr a, int lda, int p, L_MatrizFortranTr work, L_ArregloFortranInt &jpvt, int job, int &info)
{
      int pu,pl,plp1,j,jp,jt,k,kb,km1,kp1,l,maxl;
      double temp;
      double maxdia;
      bool swapk,negk;

	  pl = 1;
      pu = 0;
      info = p;
      if (job != 0)// if (job .eq. 0) go to 160
	  {

//          pivoting has been requested. rearrange the
//          the elements according to jpvt.

         for (k=1; k<=p; k++) // do 70 k = 1, p
		 {
            swapk = jpvt(k) > 0;
            negk = jpvt(k) < 0;
            jpvt(k) = k;
            if (negk) jpvt(k) = -jpvt(k);
            if (swapk == true) // if (.not.swapk) go to 60
			{
               if (k != pl) // if (k .eq. pl) go to 50
			   {
                  dswap(pl-1,a.ref(1,k),1,a.ref(1,pl),1);
                  temp = a(k,k);
                  a(k,k) = a(pl,pl);
                  a(pl,pl) = temp;
                  plp1 = pl + 1;
                  if (p >= plp1) // if (p .lt. plp1) go to 40
				  {
                  for (j = plp1; j<=p; j++) // do 30 j = plp1, p
				  {
                     if (j < k) //if (j .ge. k) go to 10
					 {
                        temp = a(pl,j);
                        a(pl,j) = a(j,k);
                        a(j,k) = temp;
                        continue; //go to 20
					 } // 10                continue
                     if (j != k) //if (j .eq. k) go to 20
					 {
                        temp = a(k,j);
                        a(k,j) = a(pl,j);
                        a(pl,j) = temp;
					 } // 20                continue
				  } // 30             continue
				  } // 40             continue
                  jpvt(k) = jpvt(pl);
                  jpvt(pl) = k;
			   } // 50          continue
               pl = pl + 1;
			} //60       continue
		 } //70    continue
         pu = p;
         if (p >= pl) // if (p .lt. pl) go to 150
		 {
         for (kb = pl; kb<=p; kb++) // do 140 kb = pl, p
		 {
            k = p - kb + pl;
            if (jpvt(k) < 0) // if (jpvt(k) .ge. 0) go to 130
			{
               jpvt(k) = -jpvt(k);
               if (pu != k) //if (pu .eq. k) go to 120
			   {
                  dswap(k-1,a.ref(1,k),1,a.ref(1,pu),1);
                  temp = a(k,k);
                  a(k,k) = a(pu,pu);
                  a(pu,pu) = temp;
                  kp1 = k + 1;
                  if (p >= kp1) // if (p .lt. kp1) go to 110
				  {
                  for (j=kp1; j<=p; j++) //do 100 j = kp1, p
				  {
                     if (j < pu) // if (j .ge. pu) go to 80
					 {
                        temp = a(k,j);
                        a(k,j) = a(j,pu);
                        a(j,pu) = temp;
                        continue; //go to 90
					 } // 80                continue
                     if (j != pu) // if (j .eq. pu) go to 90
					 {
                        temp = a(k,j);
                        a(k,j) = a(pu,j);
                        a(pu,j) = temp;
					 } // 90                continue
				  } // 100             continue
				  } // 110             continue
                  jt = jpvt(k);
                  jpvt(k) = jpvt(pu);
                  jpvt(pu) = jt;
			   } // 120          continue
               pu = pu - 1;
			} // 130       continue
		 } //140    continue
		 } // 150    continue
	  } // 160 continue
      for (k=1; k<=p; k++) // do 270 k = 1, p
	  {

//          reduction loop.

         maxdia = a(k,k);
         kp1 = k + 1;
         maxl = k;

//          determine the pivot element.

         if (!  (k < pl || k >= pu)  )//if (k .lt. pl .or. k .ge. pu) go to 190
		 {
            for (l=kp1; l<=pu; l++) //do 180 l = kp1, pu
			{
               if (a(l,l) > maxdia) // if (a(l,l) .le. maxdia) go to 170
			   {
                  maxdia = a(l,l);
                  maxl = l;
			   } // 170          continue
			} // 180       continue
		 } // 190    continue

//          quit if the pivot element is not positive.

         if (maxdia <= 0.0e0)// if (maxdia .gt. 0.0d0) go to 200
		 {
            info = k - 1;
//       ......exit
            return; // go to 280
		 } // 200    continue
         if (k != maxl)// if (k .eq. maxl) go to 210
		 {

//             start the pivoting and update jpvt.

            km1 = k - 1;
            dswap(km1,a.ref(1,k),1,a.ref(1,maxl),1);
            a(maxl,maxl) = a(k,k);
            a(k,k) = maxdia;
            jp = jpvt(maxl);
            jpvt(maxl) = jpvt(k);
            jpvt(k) = jp;
		 } // 210    continue

//          reduction step. pivoting is contained across the rows.

         work(k) = dsqrt(a(k,k));
         a(k,k) = work(k);
         if (p >= kp1) // if (p .lt. kp1) go to 260
		 {
         for (j=kp1; j<=p; j++) // do 250 j = kp1, p
		 {
            if (k != maxl)// if (k .eq. maxl) go to 240
			{
               if (j < maxl)// if (j .ge. maxl) go to 220
			   {
                  temp = a(k,j);
                  a(k,j) = a(j,maxl);
                  a(j,maxl) = temp;
                  continue; //go to 230
			   } // 220          continue
               if (j != maxl)// if (j .eq. maxl) go to 230
			   {
                  temp = a(k,j);
                  a(k,j) = a(maxl,j);
                  a(maxl,j) = temp;
			   } // 230          continue
			} // 240       continue
            a(k,j) = a(k,j)/work(k);
            work(j) = a(k,j);
            temp = -a(k,j);
            daxpy(j-k,temp,work.ref(kp1),1,a.ref(kp1,j),1);
		 } // 250    continue
		 } // 260    continue
	  } // 270 continue
  //280 continue
      return;
	 // end
}

int L_MatrizFortranTr::ilaenv(int ispec, char *name, char *opts, int n1, int n2, int n3, int n4)
{
#ifdef __USAR_L_LAPACK_H__
	fem::str_cref name_(name, (int)strlen(name)), opts_(opts, (int)strlen(opts));
	return lapack_fem::ilaenv(ispec, name_, opts_, n1, n2, n3, n4);
#else
	return -1;
#endif
}

void L_MatrizFortranTr::dsyev(char *jobz, char *uplo, int n, L_MatrizFortranTr a, int lda, L_MatrizFortranTr w, L_MatrizFortranTr work, int lwork, int &info)
{
#ifdef __USAR_L_LAPACK_H__
	lapack_fem::common cmn;

	fem::str_cref JOBZ(jobz, (int)strlen(jobz));
	fem::str_cref UPLO(uplo, (int)strlen(uplo));
	fem::arr_ref<double, 2> A(a[0]);
	fem::arr_ref<double> W(w[0]);
	fem::arr_ref<double> WORK(work[0]);

	lapack_fem::dsyev(cmn, JOBZ, UPLO, n, A, lda, W, WORK, lwork, info);
#else
	printf("dsyev no dispoinble\n");
	info = -1;
        return;
#endif
}

void L_MatrizFortranTr::dgesvd(char *jobu, char *jobvt, int m, int n, L_MatrizFortranTr a, int lda, L_MatrizFortranTr s, L_MatrizFortranTr u, int ldu, L_MatrizFortranTr vt, int ldvt, L_MatrizFortranTr work, int lwork, int &info)
{
#ifdef __USAR_L_LAPACK_H__
	lapack_fem::common cmn;
	fem::str_ref JOBU(jobu,(int)strlen(jobu));
	fem::str_ref JOBVT(jobvt, (int)strlen(jobvt));
	fem::arr_ref<double, 2> A(a[0]);
	fem::arr_ref<double> S(s[0]);
	fem::arr_ref<double, 2> U(u[0]);
	fem::arr_ref<double, 2> VT(vt[0]);
	fem::arr_ref<double> WORK(work[0]);

	lapack_fem::dgesvd(cmn, JOBU, JOBVT, m, n, A, lda, S, U, ldu, VT, ldvt, WORK, lwork, info);
#else
	printf("dgesvd no disponible\n");
	info = -1;
        return;
#endif
}


void L_MatrizFortranTr::prueba_daxpy()
{
#ifdef __USAR_L_LAPACK_H__
	L_Matrix dx, dy0, dy1, dy2;
	dx.reallocate(10,1);
	dx.fillRandomValues(-100,100);
	dy0.reallocate(10,1);
	dy0.fillRandomValues(-100,100);

	dy1 = dy0;
	dy2 = dy0;

	L_MatrizFortranTr dx1_, dy1_;
	dx1_.fijarRefTraspuesta(dx);
	dy1_.fijarRefTraspuesta(dy1);

	fem::arr_ref<double> dx2_(dx[0][0]), dy2_(dy2[0][0]);

	L_MatrizFortranTr::daxpy(10, 0.5, dx1_, 1, dy1_, 1);

	lapack_fem::daxpy(10, 0.5, dx2_, 1, dy2_, 1);

	dy1.print("dy1");
	dy2.print("dy2");

	printf("Parecido: %f\n", dy1.cosineDistanceTo(dy2));
#else
	printf("Templates variadicas no soportadas por turbo c...\n");
        return;
#endif
}



#ifdef L_FORTRAN_INTERFAZ
#define IF(X) if (X)
#define THEN {
#define ELSEIF(X) } else if (X)
#define ELSE } else
#define ENDIF }
#define RETURN return
#define LSAME(X,Y) (X==Y)
#define CALL
#define END
#define SUBROUTINE void
#define CONTINUE
#define GOTO goto
#define GO goto
#define TO
#define CHARACTER char &
#define INTEGER int
#define DOUBLE
#define PRECISION double
#define PRECISION_X L_MatrizFortranTr &
#define LOGICAL bool
#define PARAMETER const double

// Aca pueden ir funciones de Fortran en lenguaje nativo

#undef IF
#undef THEN
#undef ELSEIF
#undef ELSE
#undef ENDIF
#undef RETURN
#undef LSAME
#undef CALL
#undef END
#undef CONTINUE
#undef GOTO
#undef GO
#undef TO
#undef CHARACTER
#undef INTEGER
#undef DOUBLE
#undef PRECISION
#undef PRECISION_X
#undef LOGICAL
#undef PARAMETER
#endif //L_FORTRAN_INTERFAZ





bool L_EssentialMatrix::descompone_Rt(double tanDer1, double tanAba1, double tanDer2, double tanAba2, L_Matrix &RFinalDAA, L_Matrix &tFinalDAA, double *c1x, double *c2x)
{
	L_Matrix U, D, V, VT, Dcruz, Ra, t, RT, Pa, q, Q, PQ, Pb, H, HQ, PHQ, Pc, p1, v1, p2, v2tmp, v2, UD;
	double c1, c2;
	bool UPos, VPos;
	U = *this;
	U.svd_ordenado(D,V);
	UPos = (U.det() > 0);
	VPos = (V.det() > 0);
	// Si E es matriz esencial, -E tambien => podria ser semidefinida negativa en vez de positiva
	if (UPos == false)
		U.OP_mult(-1);
	if (VPos == false)
		V.OP_mult(-1);

	VT.transpOf(V);
	Dcruz.antiSymmetricElemental(); // es [0 1 0 ; -1 0 0 ; 0 0 1]
	// "Ra=U*Dcruz*VT;"
	{
		UD.OP_mult(U,Dcruz);
		Ra.OP_mult(UD,VT);
	}
	t.reallocate(3,1); // Debiera dar unitario por salir de U
	t(0,0) = U(0,2); // izq
	t(1,0) = U(1,2); // arr
	t(2,0) = U(2,2); // ade
	Pa.reallocate(4,4); // Una humilde matriz homogenea [Ra|t;0001)
	Pa.identity();
	Pa.setRotationComponents(Ra);
	Pa.setTranslationComponents(t);

	// Triangular punto suponiendo Pa = [Ra|t). Definir p2 y v2 (p1 y v1 ya definidos)
	// Punto es la interseccion de las lineas (0;0;0) + alfa*(tanDer1;tanAba1;1), y -RT*t + alfa*RT*(tanDer2;tanAba2;1)
	// Todo el codeMapping de la triangulacion se puede hacer muy eficiente si se trae la funcion de interseccion de lineas aca

	Q.reallocate(4,1); // Punto triangulado en el sistema global = sistema de la primera camara
	Q(3,0)=1;

	if (triangula(Ra, t, tanDer1, tanAba1, tanDer2, tanAba2, Q(0,0), Q(1,0), Q(2,0)) == false)
		return false;

	PQ.OP_mult(Pa,Q); // Punto llevado al sistema de la segunda camara

	c1 = Q(2,0) * Q(3,0); // indicaria punto atras de la primera camara
	c2 = PQ(2,0) * PQ(3,0); // indicaria punto atras de la segunda camara

	if (c1x!=NULL)
		*c1x = c1;
	if (c2x!=NULL)
		*c2x = c2;

	if (c1*c2 > 0)
	{
		if (c1 > 0)
		{
			// [Ra|t) es la configuracion correcta
			RFinalDAA = Ra;
			tFinalDAA = t;
		}
		else
		{
			// [Ra|-t) es probablemente la configuracion correcta.
			// Es lo mismo que fijar nueva camara [Ra|t)Hr y nuevo punto triangulado HrQ => se invierte adelante/atras en ambas
			RFinalDAA = Ra;
			// "tFinalDAA = -t";
			tFinalDAA.OP_mult(t,-1);
		}
	}
	else
	{
		// { [I|0), Pa = [Ra|t), Q } => { [I|0)Ht=[I|0), Pb=Pa*(-Ht), (-Ht)Q }.
		// Esto afecta el adelante/atras de la camara principal y no de la secundaria
		H.reallocate(4,4);
		H.setZero();
		H(0,0)=1;
		H(1,1)=1;
		H(2,2)=1;
		H(3,0)=-2*V(0,2);
		H(3,1)=-2*V(1,2);
		H(3,2)=-2*V(2,2);
		H(3,3)=-1;
		H.OP_mult(-1); // Se cambiaron los signos respecto al paper porque si no, se altera el signo de la rotacion y traslacion
		Pc.OP_mult(Pa,H); // Nueva matriz de camara secundaria
		HQ.OP_mult(H,Q); // Nuevo punto triangulado = punto al sistema de la camara principal (esa no hace nada)
		Pc.getRotationComponents(RFinalDAA);
		Pc.getTranslationComponents(tFinalDAA);
		// Hay que revisar que el nuevo punto triangulado este adelante de la camara principal
		c1 = HQ(2,0)*HQ(3,0);
		if (c1 < 0)
			tFinalDAA.OP_mult(-1); // Se invierte el adelante/atras de ambas
	}
	return true;
}

bool L_EssentialMatrix::triangula(L_Matrix &RaDAA, L_Matrix &tDAA, double tanDer1, double tanAba1, double tanDer2, double tanAba2, double &xDer, double &yAba, double &zAde)
{
	double v1x, v1y, v1z, p2x, p2y, p2z, v2x, v2y, v2z, mod1, mod2, dis1, v1pv2, den;
	v1x = tanDer1;
	v1y = tanAba1;
	v1z = 1;
	p2x = -RaDAA(0,0)*tDAA(0,0) -RaDAA(1,0)*tDAA(1,0) -RaDAA(2,0)*tDAA(2,0); // -RaT * tDAA
	p2y = -RaDAA(0,1)*tDAA(0,0) -RaDAA(1,1)*tDAA(1,0) -RaDAA(2,1)*tDAA(2,0);
	p2z = -RaDAA(0,2)*tDAA(0,0) -RaDAA(1,2)*tDAA(1,0) -RaDAA(2,2)*tDAA(2,0);
	v2x = RaDAA(0,0)*tanDer2 +RaDAA(1,0)*tanAba2 +RaDAA(2,0); // RaT * (tanizq2;tanAba2;1)
	v2y = RaDAA(0,1)*tanDer2 +RaDAA(1,1)*tanAba2 +RaDAA(2,1);
	v2z = RaDAA(0,2)*tanDer2 +RaDAA(1,2)*tanAba2 +RaDAA(2,2);

	mod1 = v1x*v1x + v1y*v1y + v1z*v1z;
	mod2 = v2x*v2x + v2y*v2y + v2z*v2z;
	if (mod1 > 1.01 || mod1 < 0.99)
	{
		mod1=sqrt(mod1);
		v1x/=mod1;
		v1y/=mod1;
		v1z/=mod1;
	}
	if (mod2 > 1.01 || mod2 < 0.99)
	{
		mod2=sqrt(mod2);
		v2x/=mod2;
		v2y/=mod2;
		v2z/=mod2;
	}
	v1pv2 = v1x*v2x + v1y*v2y + v1z*v2z;
	den = 1 - v1pv2*v1pv2;
	if (den == 0) // Paralelos, no intersectan en un punto del espacio euclidiano
		return false;
	dis1 = ( p2x*v1x + p2y*v1y + p2z*v1z  - ( p2x*v2x + p2y*v2y + p2z*v2z )*v1pv2
		) / den;
	xDer=v1x*dis1;
	yAba=v1y*dis1;
	zAde=v1z*dis1;
	return true;
}


namespace Sturm
{
	//
	// structure type for representing a polynomial
	//
	struct	poly {
				 int	ord;
				 double	coef[STURM_MAX_ORDER+1];
	};

	bool imprimeWasSturm = false;

	static int modp(poly *u, poly *v, poly *r);
	int modrf(int ord, double *coef, double a, double b, double *val);
	int numroots(int np, poly *sseq, int *atneg, int *atpos);
	int numchanges(int np, poly *sseq, double a);
	int buildsturm(int ord, poly *sseq);
	void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots);
	double evalpoly (int ord, double *coef, double x);



	// main.c
	//
	//	a sample driver program.

#ifdef DEBUG_STURM
	bool main_sturm(int ord, double *coefs, int &nroots, L_ptr<double> &roots)
#else
	bool main_sturm(int ord, double *coefs, int &nroots, double *roots)
#endif
	{
		poly	sseq[STURM_MAX_ORDER];
		double 	min, max; //, roots[STURM_MAX_ORDER];
		int		i, order, nchanges, np, atmin, atmax; //, nroots

		//
		// get the details...
		//

		order = ord;
		for (i = 0; i <= order; i++)
			sseq[0].coef[i]=coefs[i]; // coefs[0] es el coeficiente de menor orden

		//
		// build the Sturm sequence
		//

		np = buildsturm(order, sseq);

		//
		// get the number of real roots
		//
		nroots = numroots(np, sseq, &atmin, &atmax);

		if (nroots <= 0) // Solve: no real roots   // Error: decia == 0, pero puede ser negativo
			return false;

		//
		// calculate the bracket that the roots live in
		//
		min = -1.0;
		nchanges = numchanges(np, sseq, min);
		for (i = 0; nchanges != atmin && i != STURM_MAXPOW; i++) {
				min *= 10.0;
				nchanges = numchanges(np, sseq, min);
		}

		if (nchanges != atmin) { //solve: unable to bracket all negative roots
				printf("main_sturm(): No se pudieron atrapar todas las raices negativas\n");
				atmin = nchanges;
		}

		max = 1.0;
		nchanges = numchanges(np, sseq, max);
		for (i = 0; nchanges != atmax && i != STURM_MAXPOW; i++) {
				max *= 10.0;
				nchanges = numchanges(np, sseq, max);
		}

		if (nchanges != atmax) { // solve: unable to bracket all positive roots
				printf("main_sturm(): No se pudieron atrapar todas las raices positivas\n");
				atmax = nchanges;
		}

		nroots = atmin - atmax;

		//
		// perform the bisection.
		//
		sbisect(np, sseq, min, max, atmin, atmax, roots);
		return true;
	}

	 // util.c
	 //
	 // some utlity functions for root polishing and evaluating
	 // polynomials.

	//
	// evalpoly
	//
	//	evaluate polynomial defined in coef returning its value.
	//
	double	evalpoly (int ord, double *coef, double x)
	{
		double	*fp, f;

		fp = &coef[ord];
		f = *fp;

		for (fp--; fp >= coef; fp--)
		f = x * f + *fp;

		return(f);
	}


	//
	// modrf
	//
	//	uses the modified regula-falsi method to evaluate the root
	// in interval [a,b] of the polynomial described in coef. The
	// root is returned is returned in *val. The routine returns zero
	// if it can't converge.
	//
	int modrf(int ord, double *coef, double a, double b, double *val)
	{
		int		its;
		double	fa, fb, x, fx=0, lfx;
		double	*fp, *scoef, *ecoef;

		scoef = coef;
		ecoef = &coef[ord];

		fb = fa = *ecoef;
		for (fp = ecoef - 1; fp >= scoef; fp--) {
			fa = a * fa + *fp;
			fb = b * fb + *fp;
		}

		//
		// if there is no sign difference the method won't work
		//
		if (fa * fb > 0.0)
			return(0);

		if (fabs(fa) < STURM_RELERROR) {
			*val = a;
			return(1);
		}

		if (fabs(fb) < STURM_RELERROR) {
			*val = b;
			return(1);
		}

		lfx = fa;


		for (its = 0; its < STURM_MAXIT; its++) {

			x = (fb * a - fa * b) / (fb - fa);

			fx = *ecoef;
			for (fp = ecoef - 1; fp >= scoef; fp--)
					fx = x * fx + *fp;

			if (fabs(x) > STURM_RELERROR) {
					if (fabs(fx / x) < STURM_RELERROR) {
						*val = x;
						return(1);
					}
			} else if (fabs(fx) < STURM_RELERROR) {
					*val = x;
					return(1);
			}

			if ((fa * fx) < 0) {
					b = x;
					fb = fx;
					if ((lfx * fx) > 0)
						fa /= 2;
			} else {
					a = x;
					fa = fx;
					if ((lfx * fx) > 0)
						fb /= 2;
			}

			lfx = fx;
		}

		if (imprimeWasSturm) fprintf(stderr, "modrf overflow %f %f %f\n", a, b, fx);

		return(0);
	}

	//
	// sturm.c
	//
	//	the functions to build and evaluate the Sturm sequence
	//

	//
	// modp
	//
	//	calculates the modulus of u(x) / v(x) leaving it in r, it
	//  returns 0 if r(x) is a constant.
	//  note: this function assumes the leading coefficient of v
	//	is 1 or -1
	//
	static int modp(poly *u, poly *v, poly *r)
	{
		int		k, j;
		double	*nr, *end, *uc;

		nr = r->coef;
		end = &u->coef[u->ord];

		uc = u->coef;
		while (uc <= end)
				*nr++ = *uc++;

		if (v->coef[v->ord] < 0.0) {


				for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
					r->coef[k] = -r->coef[k];

				for (k = u->ord - v->ord; k >= 0; k--)
					for (j = v->ord + k - 1; j >= k; j--)
						r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
						* v->coef[j - k];
		} else {
				for (k = u->ord - v->ord; k >= 0; k--)
					for (j = v->ord + k - 1; j >= k; j--)
					r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
		}

		k = v->ord - 1;
		while (k >= 0 && fabs(r->coef[k]) < STURM_SMALL_ENOUGH) {
			r->coef[k] = 0.0;
			k--;
		}

		r->ord = (k < 0) ? 0 : k;

		return(r->ord);
	}

	//
	// buildsturm
	//
	//	build up a sturm sequence for a polynomial in smat, returning
	// the number of polynomials in the sequence
	//
	int buildsturm(int ord, poly *sseq)
	{
		int		i;
		double	f, *fp, *fc;
		poly	*sp;

		sseq[0].ord = ord;
		sseq[1].ord = ord - 1;


		//
		// calculate the derivative and normalise the leading
		// coefficient.
		//
		f = fabs(sseq[0].coef[ord] * ord);
		fp = sseq[1].coef;
		fc = sseq[0].coef + 1;
		for (i = 1; i <= ord; i++)
				*fp++ = *fc++ * i / f;

		//
		// construct the rest of the Sturm sequence
		//
		for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {

			//
			// reverse the sign and normalise
			//
			f = -fabs(sp->coef[sp->ord]);
			for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
					*fp /= f;
		}

		sp->coef[0] = -sp->coef[0];	// reverse the sign//

		return((int)(sp - sseq)); // Esta distancia supuestamente es pequena y se puede castear de size_t a int con seguridad
	}

	//
	// numroots
	//
	//	return the number of distinct real roots of the polynomial
	// described in sseq.
	//
	int numroots(int np, poly *sseq, int *atneg, int *atpos)
	{
			int		atposinf, atneginf;
			poly	*s;
			double	f, lf;

			atposinf = atneginf = 0;


		//
		// changes at positive infinity
		//
		lf = sseq[0].coef[sseq[0].ord];

		for (s = sseq + 1; s <= sseq + np; s++) {
				f = s->coef[s->ord];
				if (lf == 0.0 || lf * f < 0)
					atposinf++;
			lf = f;
		}

		//
		// changes at negative infinity
		//
		if (sseq[0].ord & 1)
				lf = -sseq[0].coef[sseq[0].ord];
		else
				lf = sseq[0].coef[sseq[0].ord];

		for (s = sseq + 1; s <= sseq + np; s++) {
				if (s->ord & 1)
					f = -s->coef[s->ord];
				else
					f = s->coef[s->ord];
				if (lf == 0.0 || lf * f < 0)
					atneginf++;
				lf = f;
		}

		*atneg = atneginf;
		*atpos = atposinf;

		return(atneginf - atposinf);
	}

	//
	// numchanges
	//
	//	return the number of sign changes in the Sturm sequence in
	// sseq at the value a.
	//
	int numchanges(int np, poly *sseq, double a)
	{
		int		changes;
		double	f, lf;
		poly	*s;

		changes = 0;

		lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

		for (s = sseq + 1; s <= sseq + np; s++) {
				f = evalpoly(s->ord, s->coef, a);
				if (lf == 0.0 || lf * f < 0)
					changes++;
				lf = f;
		}

		return(changes);
	}

	//
	// sbisect
	//
	//	uses a bisection based on the sturm sequence for the polynomial
	// described in sseq to isolate intervals in which roots occur,
	// the roots are returned in the roots array in order of magnitude.
	//
	void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots)
	{
		double	mid=0;
		int		n1 = 0, n2 = 0, its, atmid, nroot;

		if ((nroot = atmin - atmax) == 1) {

			//
			// first try a less expensive technique.
			//
			if (modrf(sseq->ord, sseq->coef, min, max, &roots[0]))
				return;


			//
			// if we get here we have to evaluate the root the hard
			// way by using the Sturm sequence.
			//
			for (its = 0; its < STURM_MAXIT; its++) {
					mid = (min + max) / 2;

					atmid = numchanges(np, sseq, mid);

					if (fabs(mid) > STURM_RELERROR) {
						if (fabs((max - min) / mid) < STURM_RELERROR) {
							roots[0] = mid;
							return;
						}
					} else if (fabs(max - min) < STURM_RELERROR) {
						roots[0] = mid;
						return;
					}

					if ((atmin - atmid) == 0)
						min = mid;
					else
						max = mid;
				}

			if (its == STURM_MAXIT) {
					if (imprimeWasSturm) fprintf(stderr, "sbisect: overflow min %f max %f\
						diff %e nroot %d n1 %d n2 %d\n",
						min, max, max - min, nroot, n1, n2);
				roots[0] = mid;
			}

			return;
		}

		//
		// more than one root in the interval, we have to bisect...
		//
		for (its = 0; its < STURM_MAXIT; its++) {

				mid = (min + max) / 2;

				atmid = numchanges(np, sseq, mid);

				n1 = atmin - atmid;
				n2 = atmid - atmax;


				if (n1 != 0 && n2 != 0) {
					sbisect(np, sseq, min, mid, atmin, atmid, roots);
					sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
					break;
				}

				if (n1 == 0)
					min = mid;
				else
					max = mid;
		}

		if (its == STURM_MAXIT) {
				if (imprimeWasSturm) fprintf(stderr, "sbisect: roots too close together\n");
				if (imprimeWasSturm) fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
					nroot %d n1 %d n2 %d\n",
					min, max, max - min, nroot, n1, n2);
				for (n1 = atmax; n1 < atmin; n1++)
				roots[n1 - atmax] = mid;
		}
	}
};


L_PolinGr11 operator*(const L_PolinGr11 &p1, const L_PolinGr11 &p2)
{
	L_PolinGr11 p3;
	int i, j;
	for (i=0; i<=p1.gr; i++) //Error: decia i<
	{
		for (j=0; j<=p2.gr; j++) // Error: decia i<
		{
			p3.c[i+j]+=p1.c[i]*p2.c[j];
		}
	}
	p3.gr=p1.gr+p2.gr;
	return p3;
}

L_PolinGr11 operator+(const L_PolinGr11 &p1, const L_PolinGr11 &p2)
{
	L_PolinGr11 p3;
	int i;
	if (p1.gr > p2.gr)
		p3.gr=p1.gr;
	else
		p3.gr=p2.gr;
	for (i=0; i<=p3.gr; i++)
		p3.c[i]=p1.c[i]+p2.c[i];
	return p3;
}

L_PolinGr11 operator-(const L_PolinGr11 &p1, const L_PolinGr11 &p2)
{
	L_PolinGr11 p3;
	int i;
	if (p1.gr > p2.gr)
		p3.gr=p1.gr;
	else
		p3.gr=p2.gr;
	for (i=0; i<=p3.gr; i++)
		p3.c[i]=p1.c[i]-p2.c[i];
	return p3;
}


L_PolinGr11 L_PolinGr11_matriz::det()
{
	if (li!=lj)
		throw L_ArgException();
	if (li==1)
		return elem[0][0];
	else if (li==2)
		return elem[0][0]*elem[1][1]-elem[0][1]*elem[1][0];
	L_PolinGr11 d;
	int i;
	L_PolinGr11_matriz sub;
	int signo=1;
	for (i=0; i<li; i++)
	{
		subMatrElimFilCol(sub,i,0);
		d+=elem[i][0]*sub.det()*signo;
		signo*=-1;
	}
	return d;
}

void L_PolinGr11_matriz::subMatrElimFilCol(L_PolinGr11_matriz &ret, int i0, int j0)
{
	int i, j;
	ret.reallocate(li-1, lj-1);
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
		{
			ret.elem[i][j]=elem[i][j];
		}
		for (j=j0; j<lj-1; j++)
		{
			ret.elem[i][j]=elem[i][j+1];
		}
	}
	for (i=i0; i<li-1; i++)
	{
		for (j=0; j<j0; j++)
		{
			ret.elem[i][j]=elem[i+1][j];
		}
		for (j=j0; j<lj-1; j++)
		{
			ret.elem[i][j]=elem[i+1][j+1];
		}
	}
}

void L_PolinGr11_matriz::evaluaEn(L_Matrix &m, double variable)
{
	int i, j;
	m.reallocate(li, lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			m(i,j) = elem[i][j].evaluate(variable);
}


void L_Polin_xyz::multX(L_Polin_xyz &fil)
{
	int ex, ey, ez, i, k;
	for (i=0; i<L_POLIN_XYZ_NCOL; i++)
		c[i]=0;
	for (int i=0; i<L_POLIN_XYZ_NCOL; i++)
	{
		if (fil.c[i] == 0)
			continue;
		exps(i, ex, ey, ez);
		ex++;
		k = col(ex, ey, ez);
		if (k == -1)
			continue; // Se sale del arreglo... el coeficiente deberia ser zero
		c[k]=fil.c[i];
	}
}

void L_Polin_xyz::multY(L_Polin_xyz &fil)
{
	int ex, ey, ez, i, k;
	for (i=0; i<L_POLIN_XYZ_NCOL; i++)
		c[i]=0;
	for (int i=0; i<L_POLIN_XYZ_NCOL; i++)
	{
		if (fil.c[i] == 0)
			continue;
		exps(i, ex, ey, ez);
		ey++;
		k = col(ex, ey, ez);
		if (k == -1)
			continue; // Se sale del arreglo... el coeficiente deberia ser zero
		c[k]=fil.c[i];
	}
}

void L_Polin_xyz::multZ(L_Polin_xyz &fil)
{
	int ex, ey, ez, i, k;
	for (i=0; i<L_POLIN_XYZ_NCOL; i++)
		c[i]=0;
	for (int i=0; i<L_POLIN_XYZ_NCOL; i++)
	{
		if (fil.c[i] == 0)
			continue;
		exps(i, ex, ey, ez);
		ez++;
		k = col(ex, ey, ez);
		if (k == -1)
			continue; // Se sale del arreglo... el coeficiente deberia ser zero
		c[k]=fil.c[i];
	}
}

void L_Polin_xyz::factoriza_xy_x_y_1(L_PolinGr11 &p1, L_PolinGr11 &p2, L_PolinGr11 &p3, L_PolinGr11 &p4, int exp1, int exp2, int exp3, int exp4)
{
	int i;
	for (i=0; i<=exp1; i++) // Corresponde al polinomio de z con factor xy //ERROR: decia < en vez de <=
		p1.c[i]=c[col(1,1,i)];
	p1.gr=exp1;
	for (i=0; i<=exp2; i++) // Corresponde al polinomio de z con factor x
		p2.c[i]=c[col(1,0,i)];
	p2.gr=exp2;
	for (i=0; i<=exp3; i++) // Corresponde al polinomio de z con factor y
		p3.c[i]=c[col(0,1,i)];
	p3.gr=exp3;
	for (i=0; i<=exp4; i++) // Corresponde al polinomio de z con factor 1
		p4.c[i]=c[col(0,0,i)];
	p4.gr=exp4;
}

void L_Polin_xyz::factorizaInv_xy_x_y_1(L_PolinGr11 &p1, L_PolinGr11 &p2, L_PolinGr11 &p3, L_PolinGr11 &p4, int exp1, int exp2, int exp3, int exp4)
{
	int i;
	for (i=0; i<L_POLIN_XYZ_NCOL; i++)
		c[i]=0;
	for (i=0; i<=exp1; i++) // Corresponde al polinomio de z con factor xy //ERROR: decia < en vez de <=
		c[col(1,1,i)]=p1.c[i];
	for (i=0; i<=exp2; i++) // Corresponde al polinomio de z con factor x
		c[col(1,0,i)]=p2.c[i];
	for (i=0; i<=exp3; i++) // Corresponde al polinomio de z con factor y
		c[col(0,1,i)]=p3.c[i];
	for (i=0; i<=exp4; i++) // Corresponde al polinomio de z con factor 1
		c[col(0,0,i)]=p4.c[i];
}

bool L_Algoritmo5Puntos::algoritmo5Puntos(L_Matrix E[10], int &nSols, const L_Matrix &u1v1_u2v2)
{
	// Matrices de las camaras: K1[I|0] (asociada a q) , K2[R|t] (asociada a q'), con x=izq, y=arr, z=ade
	// Entonces q'T E q = 0, con E = [t]x R
	L_Matrix QagrT;
	L_Matrix X, Y, Z, W;
	L_Matrix Res;
	const L_Matrix &m = u1v1_u2v2;
	int i;

	if (m.li < 5 || m.lj != 4 || m.begin() == NULL)
	{
		nSols = 0;
		return false;
	}

	QagrT.reallocate(m.li, 9);

	for (i=0; i<m.li; i++)
	{
		QagrT(i,0) = m(i,0)*m(i,2);
		QagrT(i,1) = m(i,1)*m(i,2);
		QagrT(i,2) =      1      *m(i,2);
		QagrT(i,3) = m(i,0)*m(i,3);
		QagrT(i,4) = m(i,1)*m(i,3);
		QagrT(i,5) =      1      *m(i,3);
		QagrT(i,6) = m(i,0)*      1;
		QagrT(i,7) = m(i,1)*      1;
		QagrT(i,8) =      1      *      1;
	}

	if (m.li == 5) // Despejar espacio nulo [X Y Z W] de la matriz esencial
	{
		if (despejar_XYZW_5p(X, Y, Z, W, QagrT) == false)
		{
			nSols = 0;
			return false;
		}
	}
	else // Revisar este...
	{
		if (despejar_XYZW_np(X, Y, Z, W, QagrT) == false)
		{
			nSols = 0;
			return false;
		}
	}
	if (calcMatrizCubicas9x20(Res, X, Y, Z, W) == false) // Calcular la matriz de restricciones cubicas de 9x20
		return false;
	// Ahora hay que pivotear la matriz de Res hasta dejarla como la tabla que está entre las ecuaciones (9) y (10)... :S
	if (Res.diagonalize_me() == false)
		return false;
	// Ahora, a calcular las matrices B y C
	L_PolinGr11_matriz B, C;
	int d_a = 100;
	L_PolinGr11 pol; // Hasta aca es 100% seguro que esta bueno :). De aca para abajo...
	int d_b = 100;
	calcMatricesBC(B, C, Res);
	calcPolinGr10(pol, B, C);
	// Ahora, calcular las "10" raices del polinomio
	int nRaices;
	int d_c = 100;
	double raices[11];
	int d_d = 100;
	pol.calculaRaicesSturm(nRaices, raices); // Deben ser a lo mas 10...

	if (d_c != 100 || d_d != 100)
		printf("L_Algoritmo5Puntos : Corrupcion del stack\n");

	if (nRaices == 0 || nRaices < 0)  // nRaices < 0 si las raices estan al infinito
		return false;

	if (nRaices > 10) // Aca aun no se ha usado E[]
	{
		printf("algoritmo5Puntos(): Error: %d matrices esenciales obtenidas!!!\n", nSols);
		return false;
	}

	// Ahora, se tienen hasta 10 raices para la variable z
	// Cada root de la variable z permite calcular un par (x,y) => se obtienen hasta 10 (x,y,z) que dan hasta 10 matrices E=xX+yY+zZ+W
	L_PolinGr11_matriz Bpeq = B.subMatrix(0,2 , 0,3); // pasa de 4x4 a 3x4, tiene solo 3 filas linealmente independientes
	int d_e = 100;
	L_Matrix Bnum, nulo;
	double x, y, z;
	
	for (i=0; i<nRaices; i++)
	{
		// Para calcular las variables (x,y) hay que evaluar la matriz (B) del paper y calcular su vector nulo
		// Como tiene que tener un vector nulo, se puede descartar una fila y queda una matriz de 3x4
		z = raices[i];
		Bpeq.evaluaEn(Bnum, raices[i]);
		nulo.espacioNuloPivoteo(Bnum); // Bnum.elem
		x = nulo(1,0) / nulo(3,0);
		y = nulo(2,0) / nulo(3,0);
		// "E[i] = x*X + y*Y + z*Z + W;"
		{
			E[i].reallocate(3,3);
#define E_xX_yY_zZ_W(u,v) ( E[i](u,v) = x*X(u,v) + y*Y(u,v) + z*Z(u,v) + W(u,v) )
			E_xX_yY_zZ_W(0,0);
			E_xX_yY_zZ_W(0,1);
			E_xX_yY_zZ_W(0,2);
			E_xX_yY_zZ_W(1,0);
			E_xX_yY_zZ_W(1,1);
			E_xX_yY_zZ_W(1,2);
			E_xX_yY_zZ_W(2,0);
			E_xX_yY_zZ_W(2,1);
			E_xX_yY_zZ_W(2,2);
#undef E_xX_yY_zZ_W
		}
		// La matriz esencial tiene valores propios {x,x,0}
		// La idea es dejarla con valores propios {1,1,0}
		double tr = E[i](0,0) + E[i](1,1) + E[i](2,2);
		double fac = 1/(2*tr);
		for (int u=0; u<3; u++)
			for (int v=0; v<3; v++)
				E[i](u,v) *= fac;
	}
	nSols = nRaices;
	return true;
}


bool L_Algoritmo5Puntos::calcMatrizCubicas9x20(L_Matrix &Rc, L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W)
{
	int grx[5], col;
	int i, j, a, b, c;
	double trm;
	L_Matrix U, Vabc;
	L_Matrix M[4], MT[4];
	M[0]=X;
	M[1]=Y;
	M[2]=Z;
	M[3]=W;
	MT[0].transpOf(X);
	MT[1].transpOf(Y);
	MT[2].transpOf(Z);
	MT[3].transpOf(W);

	L_Polin_xyz::inicializa(); // Para que funcione L_Polin_xyz::col( )

	Rc.reallocate(9, 20); // Matriz con 9 polinomios cubicos con 20 monomios c/u
	Rc.setZero();

	for (a=0; a<4; a++)
	{
		for (b=0; b<4; b++)
		{
			U.OP_mult(M[a],MT[b]);
			trm=0.5*(U(0,0)+U(1,1)+U(2,2));
			U(0,0)-=trm;
			U(1,1)-=trm;
			U(2,2)-=trm;
			for (c=0; c<4; c++)
			{
				// Se construye Vabc
				Vabc.OP_mult(U,M[c]); // Decia M[b]
				// Se tiene que: E*ET*E - 1/2 * Tr(E*ET)*E = sum(a,b,c=0..3 ; Vabc*var[a]*var[b]*var[c])
				// Se buildFrom el termino var[a]*var[b]*var[c] donde var[0]="x", var[1]="y", var[2]="z", var[3]="1"
				grx[0]=0; // exponente en x para este monomio
				grx[1]=0; // exponente en y para este monomio
				grx[2]=0; // exponente en z para este monomio
				grx[3]=0; // exponente en w para este monomio (ojo, w=1)
				grx[a]++;
				grx[b]++;
				grx[c]++; // El monomio es: x^grx[0]*y^grx[1]*z^grx[2]
				col = L_Polin_xyz::col(grx[0], grx[1], grx[2]); // Columna que corresponde a la tripleta de exponentes x,y,z
				if (col == -1 || col > 19)
				{
					printf("Error de indices en L_Algoritmo5Puntos::calcMatrizCubicas()\n");
					return false;
				}
				// Se tiene la matriz "Vabc", luego hay que agregar esa informacion
				// a la columna "pos" que corresponde a var[a]*var[b]*var[c]
				for (i=0; i<3; i++)
					for (j=0; j<3; j++)
						Rc(i+j*3,col)+=Vabc(i,j);
			}
		}
	}
	return true;
}

bool L_Algoritmo5Puntos::calcMatrizCubicas10x20(L_Matrix &Rc, L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W)
{
	double e00,e01,e02,e03,e04,e05,e06,e07,e08;
	double e10,e11,e12,e13,e14,e15,e16,e17,e18;
	double e20,e21,e22,e23,e24,e25,e26,e27,e28;
	double e30,e31,e32,e33,e34,e35,e36,e37,e38;

	double e002,e012,e022,e032,e042,e052,e062,e072,e082;
	double e102,e112,e122,e132,e142,e152,e162,e172,e182;
	double e202,e212,e222,e232,e242,e252,e262,e272,e282;
	double e302,e312,e322,e332,e342,e352,e362,e372,e382;

	double e003,e013,e023,e033,e043,e053,e063,e073,e083;
	double e103,e113,e123,e133,e143,e153,e163,e173,e183;
	double e203,e213,e223,e233,e243,e253,e263,e273,e283;
	double e303,e313,e323,e333,e343,e353,e363,e373,e383;

	e00 = X(0,0);
	e10 = Y(0,0);
	e20 = Z(0,0);
	e30 = W(0,0);
	e01 = X(0,1);
	e11 = Y(0,1);
	e21 = Z(0,1);
	e31 = W(0,1);
	e02 = X(0,2);
	e12 = Y(0,2);
	e22 = Z(0,2);
	e32 = W(0,2);
	e03 = X(1,0);
	e13 = Y(1,0);
	e23 = Z(1,0);
	e33 = W(1,0);
	e04 = X(1,1);
	e14 = Y(1,1);
	e24 = Z(1,1);
	e34 = W(1,1);
	e05 = X(1,2);
	e15 = Y(1,2);
	e25 = Z(1,2);
	e35 = W(1,2);
	e06 = X(2,0);
	e16 = Y(2,0);
	e26 = Z(2,0);
	e36 = W(2,0);
	e07 = X(2,1);
	e17 = Y(2,1);
	e27 = Z(2,1);
	e37 = W(2,1);
	e08 = X(2,2);
	e18 = Y(2,2);
	e28 = Z(2,2);
	e38 = W(2,2);


	e002 =e00*e00;
	e102 =e10*e10;
	e202 =e20*e20;
	e302 =e30*e30;
	e012 =e01*e01;
	e112 =e11*e11;
	e212 =e21*e21;
	e312 =e31*e31;
	e022 =e02*e02;
	e122 =e12*e12;
	e222 =e22*e22;
	e322 =e32*e32;
	e032 =e03*e03;
	e132 =e13*e13;
	e232 =e23*e23;
	e332 =e33*e33;
	e042 =e04*e04;
	e142 =e14*e14;
	e242 =e24*e24;
	e342 =e34*e34;
	e052 =e05*e05;
	e152 =e15*e15;
	e252 =e25*e25;
	e352 =e35*e35;
	e062 =e06*e06;
	e162 =e16*e16;
	e262 =e26*e26;
	e362 =e36*e36;
	e072 =e07*e07;
	e172 =e17*e17;
	e272 =e27*e27;
	e372 =e37*e37;
	e082 =e08*e08;
	e182 =e18*e18;
	e282 =e28*e28;
	e382 =e38*e38;

	e003 =e00*e00*e00;
	e103 =e10*e10*e10;
	e203 =e20*e20*e20;
	e303 =e30*e30*e30;
	e013 =e01*e01*e01;
	e113 =e11*e11*e11;
	e213 =e21*e21*e21;
	e313 =e31*e31*e31;
	e023 =e02*e02*e02;
	e123 =e12*e12*e12;
	e223 =e22*e22*e22;
	e323 =e32*e32*e32;
	e033 =e03*e03*e03;
	e133 =e13*e13*e13;
	e233 =e23*e23*e23;
	e333 =e33*e33*e33;
	e043 =e04*e04*e04;
	e143 =e14*e14*e14;
	e243 =e24*e24*e24;
	e343 =e34*e34*e34;
	e053 =e05*e05*e05;
	e153 =e15*e15*e15;
	e253 =e25*e25*e25;
	e353 =e35*e35*e35;
	e063 =e06*e06*e06;
	e163 =e16*e16*e16;
	e263 =e26*e26*e26;
	e363 =e36*e36*e36;
	e073 =e07*e07*e07;
	e173 =e17*e17*e17;
	e273 =e27*e27*e27;
	e373 =e37*e37*e37;
	e083 =e08*e08*e08;
	e183 =e18*e18*e18;
	e283 =e28*e28*e28;
	e383 =e38*e38*e38;

	Rc.reallocate(10,20);

	//         1 3 4 2 5 9 6 11 14 17 7 10 8 12 15 18 13 16 19 20 -> cambios de columna en MATLAB
	int ord[]={0,2,3,1,4,8,5,10,13,16,6,9,7,11,14,17,12,15,18,19};
	Rc(0,ord[0])=0.5*e003+0.5*e00*e012+0.5*e00*e022+0.5*e00*e032+e03*e01*e04+e03*e02*e05+0.5*e00*e062+e06*e01*e07+e06*e02*e08-0.5*e00*e042-0.5*e00*e052-0.5*e00*e072-0.5*e00*e082;
	Rc(0,ord[1])=e00*e11*e01+e00*e12*e02+e03*e00*e13+e03*e11*e04+e03*e01*e14+e03*e12*e05+e03*e02*e15+e13*e01*e04+e13*e02*e05+e06*e00*e16+1.5*e10*e002+0.5*e10*e012+0.5*e10*e022+0.5*e10*e062-0.5*e10*e042-0.5*e10*e052-0.5*e10*e072+0.5*e10*e032+e06*e11*e07+e06*e01*e17+e06*e12*e08+e06*e02*e18+e16*e01*e07+e16*e02*e08-e00*e14*e04-e00*e17*e07-e00*e15*e05-e00*e18*e08-0.5*e10*e082;
	Rc(0,ord[2])=e16*e02*e18+e03*e12*e15+e10*e11*e01+e10*e12*e02+e03*e10*e13+e03*e11*e14+e13*e11*e04+e13*e01*e14+e13*e12*e05+e13*e02*e15+e06*e10*e16+e06*e12*e18+e06*e11*e17+e16*e11*e07+e16*e01*e17+e16*e12*e08-e10*e14*e04-e10*e17*e07-e10*e15*e05-e10*e18*e08+1.5*e00*e102+0.5*e00*e122+0.5*e00*e112+0.5*e00*e132+0.5*e00*e162-0.5*e00*e152-0.5*e00*e172-0.5*e00*e182-0.5*e00*e142;
	Rc(0,ord[3])=0.5*e103+0.5*e10*e122+0.5*e10*e112+0.5*e10*e132+e13*e12*e15+e13*e11*e14+0.5*e10*e162+e16*e12*e18+e16*e11*e17-0.5*e10*e152-0.5*e10*e172-0.5*e10*e182-0.5*e10*e142;
	Rc(0,ord[4])=-e00*e28*e08-e00*e25*e05-e00*e27*e07-e00*e24*e04+e26*e02*e08+e26*e01*e07+e06*e02*e28+e06*e22*e08+e06*e01*e27+e06*e21*e07+e23*e02*e05+e23*e01*e04+e03*e02*e25+e03*e22*e05+e03*e01*e24+e03*e21*e04+e00*e22*e02+e00*e21*e01-0.5*e20*e082-0.5*e20*e052-0.5*e20*e072-0.5*e20*e042+e06*e00*e26+0.5*e20*e062+e03*e00*e23+0.5*e20*e022+1.5*e20*e002+0.5*e20*e032+0.5*e20*e012;
	Rc(0,ord[5])=-e10*e24*e04-e10*e27*e07-e10*e25*e05-e10*e28*e08-e20*e14*e04-e20*e17*e07-e20*e15*e05-e20*e18*e08-e00*e24*e14-e00*e25*e15-e00*e27*e17-e00*e28*e18+e06*e21*e17+e06*e22*e18+e06*e12*e28+e16*e00*e26+e16*e21*e07+e16*e01*e27+e16*e22*e08+e16*e02*e28+e26*e11*e07+e26*e01*e17+e26*e12*e08+e26*e02*e18+e06*e11*e27+e23*e11*e04+e23*e01*e14+e23*e12*e05+e23*e02*e15+e06*e20*e16+e06*e10*e26+e03*e21*e14+e03*e22*e15+e03*e12*e25+e13*e00*e23+e13*e21*e04+e13*e01*e24+e13*e22*e05+e13*e02*e25+e03*e11*e24+e03*e20*e13+e03*e10*e23+e00*e21*e11+3*e00*e20*e10+e00*e22*e12+e20*e12*e02+e20*e11*e01+e10*e22*e02+e10*e21*e01;
	Rc(0,ord[6])=-0.5*e20*e152+e26*e11*e17-e10*e24*e14-e10*e25*e15-e10*e27*e17-e10*e28*e18+0.5*e20*e162+e13*e10*e23+e13*e22*e15+e23*e12*e15+e23*e11*e14+e16*e10*e26+e16*e21*e17+e16*e11*e27+e16*e22*e18+e16*e12*e28+e26*e12*e18+e13*e12*e25+0.5*e20*e132+1.5*e20*e102+0.5*e20*e122+0.5*e20*e112+e10*e21*e11+e10*e22*e12+e13*e11*e24-0.5*e20*e172-0.5*e20*e182-0.5*e20*e142+e13*e21*e14;
	Rc(0,ord[7])=-e20*e25*e05-e20*e28*e08-0.5*e00*e272-0.5*e00*e282-0.5*e00*e242+0.5*e00*e262-0.5*e00*e252+e06*e20*e26+0.5*e00*e232+e06*e22*e28+e06*e21*e27+e26*e21*e07+e26*e01*e27+e26*e22*e08+e26*e02*e28-e20*e24*e04-e20*e27*e07+e03*e20*e23+e03*e22*e25+e03*e21*e24+e23*e21*e04+e23*e01*e24+e23*e22*e05+e23*e02*e25+e20*e21*e01+e20*e22*e02+1.5*e00*e202+0.5*e00*e222+0.5*e00*e212;
	Rc(0,ord[8])=e23*e21*e14+e23*e11*e24+e23*e22*e15+e23*e12*e25+e16*e20*e26+e16*e22*e28+e16*e21*e27+e26*e21*e17+e26*e11*e27+e26*e22*e18+e26*e12*e28+1.5*e10*e202+0.5*e10*e222+0.5*e10*e212+0.5*e10*e232+e20*e21*e11+e20*e22*e12+e13*e20*e23+e13*e22*e25+e13*e21*e24-e20*e24*e14-e20*e25*e15-e20*e27*e17-e20*e28*e18-0.5*e10*e272-0.5*e10*e282-0.5*e10*e242-0.5*e10*e252+0.5*e10*e262;
	Rc(0,ord[9])=0.5*e203+0.5*e20*e222+0.5*e20*e212+0.5*e20*e232+e23*e22*e25+e23*e21*e24+0.5*e20*e262+e26*e22*e28+e26*e21*e27-0.5*e20*e252-0.5*e20*e272-0.5*e20*e282-0.5*e20*e242;
	Rc(0,ord[10])=e06*e32*e08-0.5*e30*e082-0.5*e30*e042-0.5*e30*e052-0.5*e30*e072+0.5*e30*e012+0.5*e30*e022+0.5*e30*e032+0.5*e30*e062+1.5*e30*e002+e00*e31*e01+e00*e32*e02+e03*e31*e04+e03*e01*e34+e03*e32*e05+e03*e02*e35+e33*e01*e04+e33*e02*e05+e06*e00*e36+e06*e31*e07+e06*e01*e37+e06*e02*e38+e36*e01*e07+e36*e02*e08-e00*e34*e04-e00*e37*e07-e00*e35*e05-e00*e38*e08+e03*e00*e33;
	Rc(0,ord[11])=e06*e30*e16+e03*e30*e13+e16*e31*e07+e06*e10*e36-e10*e37*e07+3*e00*e30*e10+e00*e32*e12-e00*e38*e18-e10*e34*e04-e10*e35*e05-e10*e38*e08-e30*e14*e04-e30*e17*e07-e30*e15*e05-e30*e18*e08+e00*e31*e11+e10*e31*e01+e10*e32*e02+e30*e11*e01+e30*e12*e02+e03*e10*e33-e00*e34*e14-e00*e35*e15-e00*e37*e17+e03*e31*e14+e03*e11*e34+e03*e32*e15+e03*e12*e35+e13*e00*e33+e13*e31*e04+e13*e01*e34+e13*e32*e05+e13*e02*e35+e33*e11*e04+e33*e01*e14+e33*e12*e05+e33*e02*e15+e06*e31*e17+e06*e11*e37+e06*e32*e18+e06*e12*e38+e16*e00*e36+e16*e01*e37+e16*e32*e08+e16*e02*e38+e36*e11*e07+e36*e01*e17+e36*e12*e08+e36*e02*e18;
	Rc(0,ord[12])=e13*e10*e33+e33*e11*e14+e16*e10*e36+e16*e31*e17+e16*e11*e37+e16*e32*e18+e16*e12*e38+e36*e12*e18+e36*e11*e17-e10*e34*e14-e10*e35*e15-e10*e37*e17-e10*e38*e18+e10*e31*e11+e10*e32*e12+e13*e31*e14+e13*e11*e34+e13*e32*e15+e13*e12*e35+e33*e12*e15+1.5*e30*e102+0.5*e30*e122+0.5*e30*e112+0.5*e30*e132+0.5*e30*e162-0.5*e30*e152-0.5*e30*e172-0.5*e30*e182-0.5*e30*e142;
	Rc(0,ord[13])=e00*e32*e22+3*e00*e30*e20+e00*e31*e21+e20*e31*e01+e20*e32*e02+e30*e21*e01+e30*e22*e02+e03*e20*e33+e03*e32*e25+e03*e22*e35+e03*e31*e24+e03*e21*e34+e23*e00*e33+e23*e31*e04+e23*e01*e34+e23*e32*e05+e23*e02*e35+e33*e21*e04+e33*e01*e24+e33*e22*e05+e33*e02*e25+e06*e30*e26+e06*e20*e36+e06*e32*e28+e06*e22*e38+e06*e31*e27+e06*e21*e37+e26*e00*e36+e26*e31*e07+e03*e30*e23+e26*e01*e37+e26*e32*e08+e26*e02*e38+e36*e21*e07+e36*e01*e27+e36*e22*e08+e36*e02*e28-e00*e35*e25-e00*e37*e27-e00*e38*e28-e00*e34*e24-e20*e34*e04-e20*e37*e07-e20*e35*e05-e20*e38*e08-e30*e24*e04-e30*e27*e07-e30*e25*e05-e30*e28*e08;
	Rc(0,ord[14])=e16*e30*e26+e13*e21*e34+3*e10*e30*e20+e10*e32*e22+e10*e31*e21+e20*e31*e11+e20*e32*e12+e30*e21*e11+e30*e22*e12+e13*e30*e23+e13*e20*e33+e13*e32*e25+e13*e22*e35+e13*e31*e24+e23*e10*e33+e23*e31*e14+e23*e11*e34+e23*e32*e15+e23*e12*e35+e33*e21*e14+e33*e11*e24+e33*e22*e15+e33*e12*e25+e16*e20*e36+e16*e32*e28+e16*e22*e38+e16*e31*e27+e16*e21*e37+e26*e10*e36+e26*e31*e17+e26*e11*e37+e26*e32*e18+e26*e12*e38+e36*e21*e17+e36*e11*e27+e36*e22*e18+e36*e12*e28-e10*e35*e25-e10*e37*e27-e10*e38*e28-e10*e34*e24-e20*e34*e14-e20*e35*e15-e20*e37*e17-e20*e38*e18-e30*e24*e14-e30*e25*e15-e30*e27*e17-e30*e28*e18;
	Rc(0,ord[15])=-e20*e34*e24+0.5*e30*e262-0.5*e30*e252-0.5*e30*e272-0.5*e30*e282-0.5*e30*e242+1.5*e30*e202+0.5*e30*e222+0.5*e30*e212+0.5*e30*e232+e20*e32*e22+e20*e31*e21+e23*e20*e33+e23*e32*e25+e23*e22*e35+e23*e31*e24+e23*e21*e34+e33*e22*e25+e33*e21*e24+e26*e20*e36+e26*e32*e28+e26*e22*e38+e26*e31*e27+e26*e21*e37+e36*e22*e28+e36*e21*e27-e20*e35*e25-e20*e37*e27-e20*e38*e28;
	Rc(0,ord[16])=0.5*e00*e322+e30*e32*e02+e30*e31*e01+1.5*e00*e302+0.5*e00*e312+e03*e32*e35+e33*e31*e04+e33*e01*e34+e33*e32*e05+e33*e02*e35+e06*e30*e36+e06*e31*e37+e06*e32*e38+e36*e31*e07+e36*e01*e37+e36*e32*e08+e36*e02*e38-e30*e34*e04-e30*e37*e07-e30*e35*e05-e30*e38*e08+0.5*e00*e332+0.5*e00*e362-0.5*e00*e382-0.5*e00*e352-0.5*e00*e342-0.5*e00*e372+e03*e30*e33+e03*e31*e34;
	Rc(0,ord[17])=0.5*e10*e362-0.5*e10*e382-0.5*e10*e352-0.5*e10*e342-0.5*e10*e372+e36*e31*e17+e36*e11*e37+e36*e32*e18+e36*e12*e38-e30*e34*e14-e30*e35*e15-e30*e37*e17-e30*e38*e18+1.5*e10*e302+0.5*e10*e312+0.5*e10*e322+0.5*e10*e332+e30*e31*e11+e30*e32*e12+e13*e30*e33+e13*e31*e34+e13*e32*e35+e33*e31*e14+e33*e11*e34+e33*e32*e15+e33*e12*e35+e16*e30*e36+e16*e31*e37+e16*e32*e38;
	Rc(0,ord[18])=e33*e31*e24+e33*e21*e34+e26*e30*e36+e26*e31*e37+e26*e32*e38+e36*e32*e28+e36*e22*e38+e36*e31*e27+e36*e21*e37-e30*e35*e25-e30*e37*e27-e30*e38*e28-e30*e34*e24+e33*e22*e35+1.5*e20*e302+0.5*e20*e312+0.5*e20*e322+0.5*e20*e332+0.5*e20*e362-0.5*e20*e382-0.5*e20*e352-0.5*e20*e342-0.5*e20*e372+e30*e32*e22+e30*e31*e21+e23*e30*e33+e23*e31*e34+e23*e32*e35+e33*e32*e25;
	Rc(0,ord[19])=0.5*e303+0.5*e30*e312+0.5*e30*e322+0.5*e30*e332+e33*e31*e34+e33*e32*e35+0.5*e30*e362+e36*e31*e37+e36*e32*e38-0.5*e30*e382-0.5*e30*e352-0.5*e30*e342-0.5*e30*e372;
	Rc(1,ord[0])=e00*e01*e04+0.5*e002*e03+e00*e02*e05+0.5*e033+0.5*e03*e042+0.5*e03*e052+0.5*e03*e062+e06*e04*e07+e06*e05*e08-0.5*e03*e012-0.5*e03*e072-0.5*e03*e022-0.5*e03*e082;
	Rc(1,ord[1])=e03*e14*e04+e10*e01*e04+e16*e05*e08+e00*e10*e03+e00*e11*e04+e00*e01*e14+e00*e12*e05+e00*e02*e15+e10*e02*e05+e03*e15*e05+e06*e03*e16+e06*e14*e07+e06*e04*e17+e06*e15*e08+e06*e05*e18+0.5*e002*e13+1.5*e13*e032+0.5*e13*e042+0.5*e13*e052+0.5*e13*e062-0.5*e13*e012-0.5*e13*e072-0.5*e13*e022-0.5*e13*e082+e16*e04*e07-e03*e12*e02-e03*e11*e01-e03*e17*e07-e03*e18*e08;
	Rc(1,ord[2])=-e13*e11*e01+e00*e10*e13+e00*e12*e15+e00*e11*e14+e10*e11*e04+e10*e01*e14+e10*e12*e05+e10*e02*e15+e13*e14*e04+e13*e15*e05+e06*e13*e16+e06*e15*e18+e06*e14*e17+e16*e14*e07+e16*e04*e17+e16*e15*e08+e16*e05*e18-e13*e12*e02-e13*e17*e07-e13*e18*e08+0.5*e102*e03+1.5*e03*e132+0.5*e03*e152+0.5*e03*e142+0.5*e03*e162-0.5*e03*e112-0.5*e03*e172-0.5*e03*e122-0.5*e03*e182;
	Rc(1,ord[3])=0.5*e102*e13+e10*e11*e14+e10*e12*e15+0.5*e133+0.5*e13*e152+0.5*e13*e142+0.5*e13*e162+e16*e15*e18+e16*e14*e17-0.5*e13*e112-0.5*e13*e122-0.5*e13*e172-0.5*e13*e182;
	Rc(1,ord[4])=-e03*e28*e08-e03*e27*e07-e03*e21*e01-e03*e22*e02+e26*e05*e08+e26*e04*e07+e06*e05*e28+e06*e25*e08+e06*e04*e27+e06*e24*e07+e03*e25*e05+e03*e24*e04+e20*e02*e05+e20*e01*e04+e00*e02*e25+e00*e22*e05+e00*e01*e24+e00*e21*e04+e00*e20*e03-0.5*e23*e072-0.5*e23*e082-0.5*e23*e022-0.5*e23*e012+e06*e03*e26+0.5*e23*e052+0.5*e23*e062+1.5*e23*e032+0.5*e23*e042+0.5*e002*e23;
	Rc(1,ord[5])=e00*e21*e14+e00*e11*e24+e00*e10*e23+e00*e22*e15+e00*e12*e25+e20*e12*e05+e20*e01*e14+e20*e11*e04+e00*e20*e13+e10*e02*e25+e10*e22*e05+e10*e01*e24+e10*e21*e04+e10*e20*e03+e23*e15*e05+e23*e14*e04+e13*e25*e05+e13*e24*e04+e03*e24*e14+e03*e25*e15+3*e03*e23*e13+e20*e02*e15+e16*e03*e26+e06*e14*e27-e23*e18*e08+e06*e24*e17+e06*e15*e28+e06*e25*e18+e06*e13*e26+e06*e23*e16+e26*e04*e17+e26*e14*e07+e16*e05*e28+e16*e25*e08+e16*e04*e27+e16*e24*e07-e03*e22*e12-e03*e21*e11+e26*e05*e18+e26*e15*e08-e03*e27*e17-e03*e28*e18-e13*e22*e02-e13*e28*e08-e13*e27*e07-e13*e21*e01-e23*e17*e07-e23*e11*e01-e23*e12*e02;
	Rc(1,ord[6])=-0.5*e23*e182-0.5*e23*e172-0.5*e23*e112-0.5*e23*e122-e13*e22*e12-e13*e27*e17-e13*e28*e18+e26*e15*e18+e26*e14*e17-e13*e21*e11+e20*e12*e15+e13*e25*e15+e13*e24*e14+e16*e13*e26+e16*e25*e18+e16*e15*e28+e16*e24*e17+e16*e14*e27+1.5*e23*e132+0.5*e23*e152+0.5*e23*e142+0.5*e23*e162+e10*e20*e13+e10*e21*e14+e10*e11*e24+e10*e22*e15+e10*e12*e25+e20*e11*e14+0.5*e102*e23;
	Rc(1,ord[7])=e26*e04*e27+e00*e22*e25-e23*e28*e08+0.5*e03*e262-0.5*e03*e212-0.5*e03*e272-0.5*e03*e222-0.5*e03*e282+e23*e24*e04+e23*e25*e05+0.5*e202*e03+e06*e23*e26+e06*e24*e27+e06*e25*e28+e26*e24*e07+e26*e25*e08+e26*e05*e28-e23*e22*e02-e23*e21*e01-e23*e27*e07+e00*e20*e23+e00*e21*e24+e20*e21*e04+e20*e01*e24+e20*e22*e05+e20*e02*e25+1.5*e03*e232+0.5*e03*e242+0.5*e03*e252;
	Rc(1,ord[8])=e20*e11*e24-0.5*e13*e212-0.5*e13*e272-0.5*e13*e222-0.5*e13*e282-e23*e27*e17-e23*e28*e18+e26*e25*e18+e26*e24*e17+e26*e14*e27-e23*e21*e11-e23*e22*e12+e26*e15*e28+e23*e25*e15+e23*e24*e14+e16*e23*e26+e16*e24*e27+e16*e25*e28+0.5*e13*e262+e20*e21*e14+e20*e22*e15+e20*e12*e25+0.5*e13*e242+0.5*e13*e252+0.5*e202*e13+1.5*e13*e232+e10*e20*e23+e10*e22*e25+e10*e21*e24;
	Rc(1,ord[9])=0.5*e202*e23+e20*e22*e25+e20*e21*e24+0.5*e233+0.5*e23*e242+0.5*e23*e252+0.5*e23*e262+e26*e24*e27+e26*e25*e28-0.5*e23*e212-0.5*e23*e272-0.5*e23*e222-0.5*e23*e282;
	Rc(1,ord[10])=e00*e30*e03+0.5*e33*e062-0.5*e33*e012-0.5*e33*e022-0.5*e33*e072+e03*e35*e05+e06*e03*e36+e06*e34*e07+e06*e04*e37+e06*e35*e08+e06*e05*e38+e36*e04*e07+e36*e05*e08-e03*e32*e02-e03*e31*e01-e03*e37*e07+e00*e31*e04+e00*e01*e34+e00*e32*e05+e00*e02*e35+e30*e01*e04+e30*e02*e05+e03*e34*e04-e03*e38*e08+0.5*e002*e33+1.5*e33*e032+0.5*e33*e042+0.5*e33*e052-0.5*e33*e082;
	Rc(1,ord[11])=e06*e35*e18+e06*e33*e16+e00*e30*e13+e00*e10*e33+e00*e31*e14+e00*e11*e34+e00*e32*e15+e00*e12*e35+e10*e30*e03-e33*e17*e07-e33*e18*e08+e10*e31*e04+e10*e01*e34+e10*e32*e05+e10*e02*e35+e30*e11*e04+e30*e01*e14+e30*e12*e05+e30*e02*e15+3*e03*e33*e13+e03*e35*e15+e03*e34*e14+e13*e34*e04+e13*e35*e05+e33*e14*e04+e33*e15*e05+e06*e13*e36+e06*e15*e38+e06*e34*e17+e06*e14*e37+e16*e03*e36+e16*e34*e07+e16*e04*e37+e16*e35*e08+e16*e05*e38+e36*e14*e07+e36*e04*e17+e36*e15*e08+e36*e05*e18-e03*e31*e11-e03*e32*e12-e03*e37*e17-e03*e38*e18-e13*e32*e02-e13*e31*e01-e13*e37*e07-e13*e38*e08-e33*e12*e02-e33*e11*e01;
	Rc(1,ord[12])=e16*e13*e36+e10*e11*e34+0.5*e33*e152+0.5*e33*e142+0.5*e33*e162-0.5*e33*e112-0.5*e33*e122-0.5*e33*e172-0.5*e33*e182+0.5*e102*e33+1.5*e33*e132+e10*e30*e13+e10*e31*e14+e10*e32*e15+e10*e12*e35+e30*e11*e14+e30*e12*e15+e13*e35*e15+e13*e34*e14+e16*e35*e18+e16*e15*e38+e16*e34*e17+e16*e14*e37+e36*e15*e18+e36*e14*e17-e13*e31*e11-e13*e32*e12-e13*e37*e17-e13*e38*e18;
	Rc(1,ord[13])=e06*e35*e28+e36*e04*e27+e00*e20*e33+e00*e30*e23+3*e03*e33*e23+e03*e34*e24+e03*e35*e25+e23*e34*e04+e23*e35*e05+e33*e24*e04+e33*e25*e05+e06*e33*e26+e06*e23*e36+e06*e34*e27+e06*e24*e37+e06*e25*e38+e26*e03*e36+e26*e34*e07+e26*e04*e37+e26*e35*e08+e26*e05*e38+e36*e24*e07+e36*e25*e08+e36*e05*e28-e03*e31*e21-e03*e37*e27-e03*e32*e22-e03*e38*e28-e23*e32*e02-e23*e31*e01-e23*e37*e07-e23*e38*e08-e33*e22*e02-e33*e21*e01-e33*e27*e07-e33*e28*e08+e00*e32*e25+e00*e22*e35+e00*e31*e24+e00*e21*e34+e20*e30*e03+e20*e31*e04+e20*e01*e34+e20*e32*e05+e20*e02*e35+e30*e21*e04+e30*e01*e24+e30*e22*e05+e30*e02*e25;
	Rc(1,ord[14])=e10*e30*e23+e10*e20*e33+e10*e22*e35+e10*e32*e25+e10*e31*e24+e10*e21*e34+e20*e30*e13+e20*e31*e14+e20*e11*e34+e20*e32*e15+e20*e12*e35+e30*e21*e14+e30*e11*e24+e30*e22*e15+e30*e12*e25+3*e13*e33*e23+e13*e34*e24+e13*e35*e25+e23*e35*e15+e23*e34*e14+e33*e25*e15+e33*e24*e14+e16*e33*e26+e16*e23*e36+e16*e34*e27+e16*e24*e37+e16*e35*e28+e16*e25*e38+e26*e13*e36+e26*e35*e18+e26*e15*e38+e26*e34*e17+e26*e14*e37+e36*e25*e18+e36*e15*e28+e36*e24*e17+e36*e14*e27-e13*e31*e21-e13*e37*e27-e13*e32*e22-e13*e38*e28-e23*e31*e11-e23*e32*e12-e23*e37*e17-e23*e38*e18-e33*e21*e11-e33*e22*e12-e33*e27*e17-e33*e28*e18;
	Rc(1,ord[15])=-0.5*e33*e212-0.5*e33*e272-0.5*e33*e222-0.5*e33*e282+e26*e23*e36+e20*e30*e23+e20*e32*e25+e20*e22*e35+e20*e31*e24+e20*e21*e34+e30*e22*e25+e30*e21*e24+e23*e34*e24+e23*e35*e25+e26*e34*e27+e26*e24*e37+e26*e35*e28+e26*e25*e38+e36*e24*e27+e36*e25*e28-e23*e31*e21-e23*e37*e27-e23*e32*e22-e23*e38*e28+0.5*e202*e33+1.5*e33*e232+0.5*e33*e242+0.5*e33*e252+0.5*e33*e262;
	Rc(1,ord[16])=e33*e35*e05+e30*e32*e05+0.5*e03*e362+0.5*e302*e03+1.5*e03*e332+0.5*e03*e352+0.5*e03*e342+e00*e30*e33+e00*e31*e34+e00*e32*e35+e30*e31*e04+e30*e01*e34+e30*e02*e35+e33*e34*e04+e06*e33*e36+e06*e35*e38+e06*e34*e37+e36*e34*e07+e36*e04*e37+e36*e35*e08+e36*e05*e38-e33*e32*e02-e33*e31*e01-e33*e37*e07-e33*e38*e08-0.5*e03*e322-0.5*e03*e382-0.5*e03*e312-0.5*e03*e372;
	Rc(1,ord[17])=-e33*e31*e11-e33*e32*e12-e33*e38*e18+e30*e11*e34+e30*e32*e15+e30*e12*e35+e33*e35*e15+e33*e34*e14+e16*e33*e36+e16*e35*e38+e16*e34*e37+e36*e35*e18+e36*e15*e38+e36*e34*e17+e36*e14*e37-e33*e37*e17+0.5*e302*e13+1.5*e13*e332+0.5*e13*e352+0.5*e13*e342+0.5*e13*e362-0.5*e13*e322-0.5*e13*e382-0.5*e13*e312-0.5*e13*e372+e10*e30*e33+e10*e31*e34+e10*e32*e35+e30*e31*e14;
	Rc(1,ord[18])=e36*e25*e38+0.5*e302*e23+1.5*e23*e332+0.5*e23*e352+0.5*e23*e342+0.5*e23*e362-0.5*e23*e322-0.5*e23*e382-0.5*e23*e312-0.5*e23*e372+e20*e30*e33+e20*e31*e34+e20*e32*e35+e30*e32*e25+e30*e22*e35+e30*e31*e24+e30*e21*e34+e33*e34*e24+e33*e35*e25+e26*e33*e36+e26*e35*e38+e26*e34*e37+e36*e34*e27+e36*e24*e37+e36*e35*e28-e33*e31*e21-e33*e37*e27-e33*e32*e22-e33*e38*e28;
	Rc(1,ord[19])=0.5*e302*e33+e30*e31*e34+e30*e32*e35+0.5*e333+0.5*e33*e352+0.5*e33*e342+0.5*e33*e362+e36*e35*e38+e36*e34*e37-0.5*e33*e322-0.5*e33*e382-0.5*e33*e312-0.5*e33*e372;
	Rc(2,ord[0])=0.5*e002*e06+e00*e01*e07+e00*e02*e08+0.5*e032*e06+e03*e04*e07+e03*e05*e08+0.5*e063+0.5*e06*e072+0.5*e06*e082-0.5*e06*e012-0.5*e06*e022-0.5*e06*e042-0.5*e06*e052;
	Rc(2,ord[1])=e00*e10*e06+0.5*e002*e16+0.5*e032*e16+1.5*e16*e062+0.5*e16*e072+0.5*e16*e082-0.5*e16*e012-0.5*e16*e022-0.5*e16*e042-0.5*e16*e052+e00*e11*e07+e00*e01*e17+e00*e12*e08+e00*e02*e18+e10*e01*e07+e10*e02*e08+e03*e13*e06+e03*e14*e07+e03*e04*e17+e03*e15*e08+e03*e05*e18+e13*e04*e07+e13*e05*e08+e06*e17*e07+e06*e18*e08-e06*e12*e02-e06*e11*e01-e06*e14*e04-e06*e15*e05;
	Rc(2,ord[2])=e13*e14*e07+0.5*e102*e06+e00*e10*e16+e00*e12*e18+e00*e11*e17+e10*e11*e07+e10*e01*e17+e10*e12*e08+e10*e02*e18+e03*e13*e16+e03*e15*e18+e03*e14*e17+e13*e04*e17+e13*e15*e08+e13*e05*e18+e16*e17*e07+e16*e18*e08-e16*e12*e02-e16*e11*e01-e16*e14*e04-e16*e15*e05+0.5*e132*e06+1.5*e06*e162+0.5*e06*e182+0.5*e06*e172-0.5*e06*e112-0.5*e06*e122-0.5*e06*e142-0.5*e06*e152;
	Rc(2,ord[3])=0.5*e102*e16+e10*e12*e18+e10*e11*e17+0.5*e132*e16+e13*e15*e18+e13*e14*e17+0.5*e163+0.5*e16*e182+0.5*e16*e172-0.5*e16*e112-0.5*e16*e122-0.5*e16*e142-0.5*e16*e152;
	Rc(2,ord[4])=e06*e27*e07+e23*e05*e08+e23*e04*e07+e03*e05*e28+e03*e25*e08+e03*e04*e27+e03*e24*e07+e20*e02*e08+e20*e01*e07+e00*e02*e28+e00*e22*e08+e00*e01*e27+e00*e21*e07+e00*e20*e06-e06*e25*e05-e06*e24*e04-e06*e21*e01-e06*e22*e02+e06*e28*e08-0.5*e26*e042-0.5*e26*e052-0.5*e26*e012-0.5*e26*e022+0.5*e26*e082+0.5*e26*e072+1.5*e26*e062+0.5*e002*e26+e03*e23*e06+0.5*e032*e26;
	Rc(2,ord[5])=e13*e05*e28+e00*e12*e28+e13*e25*e08+e13*e04*e27+e13*e24*e07+e13*e23*e06+e03*e14*e27+e03*e24*e17+e03*e15*e28+e03*e25*e18+e03*e13*e26+e03*e23*e16+e20*e02*e18+e20*e12*e08+e20*e01*e17+e20*e11*e07+e00*e21*e17+e10*e02*e28+e10*e22*e08+e10*e01*e27+e10*e21*e07+e10*e20*e06+e00*e11*e27-e26*e15*e05-e26*e14*e04-e26*e11*e01-e26*e12*e02-e16*e25*e05-e16*e24*e04-e16*e21*e01-e16*e22*e02-e06*e24*e14-e06*e22*e12-e06*e21*e11-e06*e25*e15+e00*e20*e16+e00*e22*e18+e00*e10*e26+e26*e18*e08+e26*e17*e07+e16*e28*e08+e16*e27*e07+e06*e27*e17+e06*e28*e18+3*e06*e26*e16+e23*e05*e18+e23*e15*e08+e23*e04*e17+e23*e14*e07;
	Rc(2,ord[6])=e10*e22*e18+0.5*e26*e182+0.5*e26*e172+e16*e28*e18+e16*e27*e17-e16*e25*e15-e16*e21*e11-e16*e22*e12+1.5*e26*e162+e13*e15*e28+e13*e24*e17+e13*e14*e27+e23*e15*e18+e23*e14*e17+e10*e12*e28+e10*e21*e17+e10*e11*e27+e20*e12*e18+e20*e11*e17+e13*e23*e16+e13*e25*e18+e10*e20*e16+0.5*e102*e26-0.5*e26*e122-0.5*e26*e142-0.5*e26*e152-e16*e24*e14-0.5*e26*e112+0.5*e132*e26;
	Rc(2,ord[7])=-0.5*e06*e212-0.5*e06*e252-0.5*e06*e242+0.5*e06*e272+0.5*e06*e282-0.5*e06*e222+e20*e02*e28+e03*e23*e26+e03*e24*e27+e03*e25*e28+e23*e24*e07+e23*e04*e27+e23*e25*e08+e23*e05*e28+e26*e28*e08-e26*e22*e02-e26*e21*e01-e26*e24*e04-e26*e25*e05+e26*e27*e07+e00*e20*e26+e00*e21*e27+e00*e22*e28+e20*e21*e07+e20*e01*e27+e20*e22*e08+0.5*e202*e06+0.5*e232*e06+1.5*e06*e262;
	Rc(2,ord[8])=-e26*e24*e14-0.5*e16*e212-0.5*e16*e252-0.5*e16*e242-e26*e25*e15-0.5*e16*e222-e26*e21*e11+e26*e28*e18+e26*e27*e17-e26*e22*e12+e23*e15*e28+e23*e24*e17+e23*e14*e27+0.5*e232*e16+1.5*e16*e262+0.5*e16*e272+0.5*e16*e282+e10*e20*e26+e10*e21*e27+e10*e22*e28+e20*e22*e18+e20*e12*e28+e20*e21*e17+e20*e11*e27+e13*e23*e26+e13*e24*e27+e13*e25*e28+e23*e25*e18+0.5*e202*e16;
	Rc(2,ord[9])=0.5*e202*e26+e20*e21*e27+e20*e22*e28+0.5*e232*e26+e23*e24*e27+e23*e25*e28+0.5*e263+0.5*e26*e272+0.5*e26*e282-0.5*e26*e222-0.5*e26*e212-0.5*e26*e252-0.5*e26*e242;
	Rc(2,ord[10])=e03*e34*e07+0.5*e032*e36+1.5*e36*e062+e03*e33*e06+e00*e31*e07+e00*e01*e37+e00*e32*e08+e00*e02*e38+e30*e01*e07+e30*e02*e08+e03*e04*e37+e03*e35*e08+e03*e05*e38+0.5*e002*e36-0.5*e36*e022-0.5*e36*e042-0.5*e36*e052+0.5*e36*e072+0.5*e36*e082-0.5*e36*e012+e33*e04*e07+e33*e05*e08+e06*e37*e07+e06*e38*e08-e06*e32*e02-e06*e31*e01-e06*e34*e04-e06*e35*e05+e00*e30*e06;
	Rc(2,ord[11])=e13*e33*e06+e13*e34*e07+e13*e04*e37+e13*e35*e08+e13*e05*e38+e33*e14*e07+e33*e04*e17+e33*e15*e08+e33*e05*e18+3*e06*e36*e16+e06*e38*e18+e06*e37*e17+e16*e37*e07+e16*e38*e08+e36*e17*e07+e36*e18*e08-e06*e35*e15-e06*e31*e11-e06*e32*e12+e00*e31*e17+e00*e11*e37+e10*e30*e06+e10*e31*e07+e10*e01*e37+e10*e32*e08+e10*e02*e38+e30*e11*e07+e30*e01*e17+e30*e12*e08+e30*e02*e18+e03*e33*e16+e03*e13*e36+e03*e35*e18+e03*e15*e38+e03*e34*e17+e03*e14*e37+e00*e30*e16+e00*e12*e38-e06*e34*e14-e16*e32*e02-e16*e31*e01-e16*e34*e04-e16*e35*e05-e36*e12*e02-e36*e11*e01-e36*e14*e04-e36*e15*e05+e00*e10*e36+e00*e32*e18;
	Rc(2,ord[12])=0.5*e36*e182+0.5*e36*e172-0.5*e36*e112-0.5*e36*e122-0.5*e36*e142-0.5*e36*e152+0.5*e102*e36+0.5*e132*e36+1.5*e36*e162+e10*e30*e16+e10*e32*e18+e10*e12*e38+e10*e31*e17+e10*e11*e37+e30*e12*e18+e30*e11*e17+e13*e33*e16+e13*e35*e18+e13*e15*e38+e13*e34*e17+e13*e14*e37+e33*e15*e18+e33*e14*e17+e16*e38*e18+e16*e37*e17-e16*e35*e15-e16*e31*e11-e16*e32*e12-e16*e34*e14;
	Rc(2,ord[13])=e00*e20*e36+e00*e31*e27+e00*e21*e37+e00*e32*e28+e00*e22*e38+e20*e30*e06+e20*e31*e07+e20*e01*e37+e20*e32*e08+e20*e02*e38+e30*e21*e07+e30*e01*e27+e30*e22*e08+e30*e02*e28+e03*e33*e26+e03*e23*e36+e03*e34*e27+e03*e24*e37+e03*e35*e28-e26*e31*e01-e26*e35*e05-e36*e22*e02-e36*e21*e01-e36*e24*e04-e36*e25*e05-e26*e34*e04+e03*e25*e38+e23*e34*e07+e23*e04*e37+e23*e35*e08+e23*e05*e38+e33*e24*e07+e33*e04*e27+e33*e25*e08+e33*e05*e28+3*e06*e36*e26+e06*e37*e27+e06*e38*e28+e26*e37*e07+e26*e38*e08+e36*e27*e07+e36*e28*e08-e06*e32*e22-e06*e31*e21-e06*e35*e25-e06*e34*e24-e26*e32*e02+e00*e30*e26+e23*e33*e06;
	Rc(2,ord[14])=e10*e30*e26+e10*e20*e36+e10*e31*e27+e10*e21*e37+e10*e32*e28+e10*e22*e38+e20*e30*e16+e20*e32*e18+e20*e12*e38+e20*e31*e17+e20*e11*e37+e30*e22*e18+e30*e12*e28+e30*e21*e17+e30*e11*e27+e13*e33*e26+e13*e23*e36+e13*e34*e27+e13*e24*e37+e13*e35*e28+e13*e25*e38+e23*e33*e16+e23*e35*e18+e23*e15*e38+e23*e34*e17+e23*e14*e37+e33*e25*e18+e33*e15*e28+e33*e24*e17+e33*e14*e27+3*e16*e36*e26+e16*e37*e27+e16*e38*e28+e26*e38*e18+e26*e37*e17+e36*e28*e18+e36*e27*e17-e16*e32*e22-e16*e31*e21-e16*e35*e25-e16*e34*e24-e26*e35*e15-e26*e31*e11-e26*e32*e12-e26*e34*e14-e36*e25*e15-e36*e21*e11-e36*e22*e12-e36*e24*e14;
	Rc(2,ord[15])=e33*e25*e28+e20*e30*e26+e20*e32*e28+e20*e31*e27+e20*e21*e37+e20*e22*e38+e30*e21*e27+e30*e22*e28+e23*e33*e26+e23*e34*e27+e23*e24*e37+e23*e35*e28+e23*e25*e38+e33*e24*e27+e26*e37*e27+e26*e38*e28-e26*e32*e22-e26*e31*e21-e26*e35*e25-e26*e34*e24+0.5*e202*e36+0.5*e232*e36+1.5*e36*e262+0.5*e36*e272+0.5*e36*e282-0.5*e36*e222-0.5*e36*e212-0.5*e36*e252-0.5*e36*e242;
	Rc(2,ord[16])=e00*e30*e36+e00*e32*e38+e00*e31*e37+e30*e31*e07+e30*e01*e37+e30*e32*e08+e30*e02*e38+e03*e33*e36-0.5*e06*e342+e03*e35*e38+e33*e34*e07+e33*e04*e37+e33*e35*e08+e33*e05*e38+e36*e37*e07+e36*e38*e08-e36*e32*e02-e36*e31*e01-e36*e34*e04-e36*e35*e05+e03*e34*e37+0.5*e302*e06+0.5*e332*e06+1.5*e06*e362+0.5*e06*e382+0.5*e06*e372-0.5*e06*e352-0.5*e06*e312-0.5*e06*e322;
	Rc(2,ord[17])=-e36*e35*e15+e10*e30*e36+0.5*e302*e16+0.5*e332*e16+1.5*e16*e362+0.5*e16*e382+0.5*e16*e372-0.5*e16*e352-0.5*e16*e312-0.5*e16*e322-0.5*e16*e342+e10*e32*e38+e10*e31*e37+e30*e32*e18+e30*e12*e38+e30*e31*e17+e30*e11*e37+e13*e33*e36+e13*e35*e38+e13*e34*e37+e33*e35*e18+e33*e15*e38+e33*e34*e17+e33*e14*e37+e36*e38*e18+e36*e37*e17-e36*e31*e11-e36*e32*e12-e36*e34*e14;
	Rc(2,ord[18])=-e36*e35*e25+e30*e32*e28+0.5*e302*e26+0.5*e332*e26+1.5*e26*e362+0.5*e26*e382+0.5*e26*e372-0.5*e26*e352-0.5*e26*e312-0.5*e26*e322-0.5*e26*e342+e20*e30*e36+e20*e32*e38+e20*e31*e37+e30*e31*e27+e30*e21*e37+e30*e22*e38+e23*e33*e36+e23*e35*e38+e23*e34*e37+e33*e34*e27+e33*e24*e37+e33*e35*e28+e33*e25*e38+e36*e37*e27+e36*e38*e28-e36*e32*e22-e36*e31*e21-e36*e34*e24;
	Rc(2,ord[19])=0.5*e302*e36+e30*e32*e38+e30*e31*e37+0.5*e332*e36+e33*e35*e38+e33*e34*e37+0.5*e363+0.5*e36*e382+0.5*e36*e372-0.5*e36*e352-0.5*e36*e312-0.5*e36*e322-0.5*e36*e342;
	Rc(3,ord[0])=0.5*e01*e002+0.5*e013+0.5*e01*e022+e04*e00*e03+0.5*e01*e042+e04*e02*e05+e07*e00*e06+0.5*e01*e072+e07*e02*e08-0.5*e01*e032-0.5*e01*e052-0.5*e01*e062-0.5*e01*e082;
	Rc(3,ord[1])=1.5*e11*e012+0.5*e11*e002+0.5*e11*e022+0.5*e11*e042+0.5*e11*e072-0.5*e11*e032-0.5*e11*e052-0.5*e11*e062-0.5*e11*e082+e01*e10*e00+e01*e12*e02+e04*e10*e03+e04*e00*e13+e04*e01*e14+e04*e12*e05+e04*e02*e15+e14*e00*e03+e14*e02*e05+e07*e10*e06+e07*e00*e16+e07*e01*e17+e07*e12*e08+e07*e02*e18+e17*e00*e06+e17*e02*e08-e01*e13*e03-e01*e16*e06-e01*e15*e05-e01*e18*e08;
	Rc(3,ord[2])=e17*e02*e18+e14*e10*e03+e11*e12*e02-e11*e18*e08+0.5*e01*e102+0.5*e01*e122+1.5*e01*e112+0.5*e01*e142+0.5*e01*e172-0.5*e01*e132-0.5*e01*e152-0.5*e01*e162-0.5*e01*e182+e11*e10*e00+e04*e10*e13+e04*e12*e15+e04*e11*e14+e14*e00*e13+e14*e12*e05+e14*e02*e15+e07*e10*e16+e07*e12*e18+e07*e11*e17+e17*e10*e06+e17*e00*e16+e17*e12*e08-e11*e13*e03-e11*e16*e06-e11*e15*e05;
	Rc(3,ord[3])=0.5*e11*e102+0.5*e11*e122+0.5*e113+e14*e10*e13+e14*e12*e15+0.5*e11*e142+e17*e10*e16+e17*e12*e18+0.5*e11*e172-0.5*e11*e132-0.5*e11*e152-0.5*e11*e162-0.5*e11*e182;
	Rc(3,ord[4])=-e01*e25*e05-e01*e26*e06-e01*e23*e03+e27*e02*e08+e27*e00*e06+e07*e02*e28+e07*e22*e08+e07*e01*e27+e07*e00*e26+e24*e02*e05+e24*e00*e03+e04*e02*e25+e04*e22*e05+e04*e01*e24+e04*e00*e23+e04*e20*e03+e01*e22*e02+e01*e20*e00-e01*e28*e08+e07*e20*e06+0.5*e21*e072+0.5*e21*e042+0.5*e21*e022+0.5*e21*e002+1.5*e21*e012-0.5*e21*e082-0.5*e21*e052-0.5*e21*e062-0.5*e21*e032;
	Rc(3,ord[5])=e11*e20*e00+e07*e20*e16+3*e01*e21*e11+e01*e22*e12-e21*e18*e08-e21*e15*e05-e21*e16*e06-e21*e13*e03-e11*e28*e08-e11*e25*e05-e11*e26*e06-e11*e23*e03-e01*e28*e18-e01*e23*e13-e01*e25*e15-e01*e26*e16+e27*e02*e18+e27*e12*e08+e27*e00*e16+e27*e10*e06+e17*e02*e28+e17*e22*e08+e17*e01*e27+e17*e00*e26+e17*e20*e06+e07*e11*e27+e07*e21*e17+e07*e12*e28+e07*e22*e18+e07*e10*e26+e24*e02*e15+e24*e12*e05+e24*e00*e13+e24*e10*e03+e14*e02*e25+e14*e22*e05+e14*e01*e24+e14*e00*e23+e14*e20*e03+e04*e11*e24+e04*e21*e14+e04*e12*e25+e04*e22*e15+e21*e12*e02+e04*e20*e13+e01*e20*e10+e11*e22*e02+e21*e10*e00+e04*e10*e23;
	Rc(3,ord[6])=1.5*e21*e112+0.5*e21*e102+0.5*e21*e122+e11*e20*e10+e11*e22*e12+e14*e10*e23+e14*e22*e15+e14*e12*e25-0.5*e21*e162-0.5*e21*e152-0.5*e21*e132-0.5*e21*e182+e27*e12*e18-e11*e26*e16-e11*e25*e15-e11*e23*e13-e11*e28*e18+e17*e20*e16+e17*e10*e26+e17*e22*e18+e17*e12*e28+e17*e11*e27+e27*e10*e16+0.5*e21*e172+e14*e11*e24+e24*e10*e13+e24*e12*e15+0.5*e21*e142+e14*e20*e13;
	Rc(3,ord[7])=-0.5*e01*e262-0.5*e01*e282-0.5*e01*e252-0.5*e01*e232+0.5*e01*e272+e27*e22*e08+e27*e02*e28-e21*e23*e03-e21*e26*e06-e21*e25*e05-e21*e28*e08+e04*e22*e25+e24*e20*e03+e24*e00*e23+e24*e22*e05+e24*e02*e25+e07*e20*e26+e07*e21*e27+e07*e22*e28+e27*e20*e06+e27*e00*e26+e21*e20*e00+e21*e22*e02+e04*e20*e23+e04*e21*e24+0.5*e01*e222+0.5*e01*e242+1.5*e01*e212+0.5*e01*e202;
	Rc(3,ord[8])=-0.5*e11*e282-0.5*e11*e252-e21*e26*e16+e27*e12*e28-e21*e25*e15-e21*e23*e13-e21*e28*e18+e17*e20*e26+e17*e21*e27+e17*e22*e28+e27*e20*e16+e27*e10*e26+e27*e22*e18+0.5*e11*e242+0.5*e11*e272-0.5*e11*e232-0.5*e11*e262+0.5*e11*e202+1.5*e11*e212+0.5*e11*e222+e21*e20*e10+e14*e20*e23+e14*e21*e24+e14*e22*e25+e24*e20*e13+e24*e10*e23+e24*e22*e15+e24*e12*e25+e21*e22*e12;
	Rc(3,ord[9])=0.5*e21*e202+0.5*e213+0.5*e21*e222+e24*e20*e23+0.5*e21*e242+e24*e22*e25+e27*e20*e26+0.5*e21*e272+e27*e22*e28-0.5*e21*e232-0.5*e21*e262-0.5*e21*e282-0.5*e21*e252;
	Rc(3,ord[10])=-0.5*e31*e032-0.5*e31*e052-0.5*e31*e062-0.5*e31*e082+e07*e30*e06+e07*e00*e36+e07*e01*e37+e07*e32*e08+e07*e02*e38+e37*e00*e06+e37*e02*e08-e01*e33*e03-e01*e36*e06-e01*e35*e05-e01*e38*e08+0.5*e31*e072+e04*e30*e03+e04*e00*e33+e04*e01*e34+e04*e32*e05+e04*e02*e35+e34*e00*e03+e34*e02*e05+0.5*e31*e002+0.5*e31*e022+0.5*e31*e042+e01*e30*e00+e01*e32*e02+1.5*e31*e012;
	Rc(3,ord[11])=e34*e12*e05+e34*e02*e15+e07*e10*e36+e07*e32*e18+e07*e12*e38+e07*e31*e17+e07*e11*e37+e17*e30*e06+e17*e00*e36+e17*e01*e37+e17*e32*e08+e17*e02*e38+e37*e10*e06+e37*e00*e16+e37*e12*e08+e37*e02*e18-e01*e36*e16-e01*e35*e15-e01*e33*e13-e01*e38*e18-e11*e33*e03-e11*e36*e06-e11*e35*e05+e01*e30*e10+e01*e32*e12+3*e01*e31*e11+e11*e30*e00+e11*e32*e02+e31*e10*e00+e31*e12*e02+e04*e30*e13+e04*e10*e33+e04*e32*e15+e04*e12*e35+e04*e31*e14+e04*e11*e34+e14*e30*e03+e14*e00*e33+e14*e01*e34+e14*e32*e05+e14*e02*e35+e34*e10*e03+e34*e00*e13+e07*e30*e16-e11*e38*e08-e31*e13*e03-e31*e16*e06-e31*e15*e05-e31*e18*e08;
	Rc(3,ord[12])=-e11*e33*e13-e11*e38*e18+0.5*e31*e142+0.5*e31*e172-0.5*e31*e162-0.5*e31*e152-0.5*e31*e132-0.5*e31*e182+0.5*e31*e122+0.5*e31*e102+e11*e30*e10+e11*e32*e12+e14*e30*e13+e14*e10*e33+e14*e32*e15+e14*e12*e35+e14*e11*e34+e34*e10*e13+e34*e12*e15+e17*e30*e16+e17*e10*e36+e17*e32*e18+e17*e12*e38+e17*e11*e37+e37*e10*e16+e37*e12*e18-e11*e36*e16-e11*e35*e15+1.5*e31*e112;
	Rc(3,ord[13])=-e21*e35*e05+e07*e32*e28+e01*e30*e20-e21*e33*e03-e21*e36*e06-e21*e38*e08-e31*e23*e03-e31*e26*e06-e31*e25*e05-e31*e28*e08+3*e01*e31*e21+e01*e32*e22+e21*e30*e00+e21*e32*e02+e31*e20*e00+e31*e22*e02+e04*e30*e23+e04*e20*e33+e04*e31*e24+e04*e21*e34+e04*e32*e25+e04*e22*e35+e24*e30*e03+e24*e00*e33+e24*e01*e34+e24*e32*e05+e24*e02*e35+e34*e20*e03+e34*e00*e23+e34*e22*e05+e34*e02*e25+e07*e30*e26+e07*e20*e36+e07*e31*e27+e07*e21*e37+e07*e22*e38+e27*e30*e06+e27*e00*e36+e27*e01*e37+e27*e32*e08+e27*e02*e38+e37*e00*e26+e37*e22*e08+e37*e02*e28-e01*e33*e23-e01*e36*e26-e01*e38*e28-e01*e35*e25+e37*e20*e06;
	Rc(3,ord[14])=e11*e32*e22+e34*e12*e25+e11*e30*e20+3*e11*e31*e21+e21*e30*e10+e21*e32*e12+e34*e10*e23+e34*e22*e15+e17*e30*e26+e17*e20*e36+e17*e31*e27+e17*e21*e37+e17*e32*e28+e17*e22*e38+e27*e30*e16+e27*e10*e36+e27*e32*e18+e27*e12*e38+e27*e11*e37+e37*e20*e16+e37*e10*e26+e37*e22*e18+e37*e12*e28-e11*e33*e23-e11*e36*e26-e11*e38*e28-e11*e35*e25-e21*e36*e16-e21*e35*e15-e21*e33*e13-e21*e38*e18-e31*e26*e16-e31*e25*e15-e31*e23*e13-e31*e28*e18+e31*e20*e10+e31*e22*e12+e14*e30*e23+e14*e20*e33+e14*e31*e24+e14*e21*e34+e14*e32*e25+e14*e22*e35+e24*e30*e13+e24*e10*e33+e24*e32*e15+e24*e12*e35+e24*e11*e34+e34*e20*e13;
	Rc(3,ord[15])=-e21*e36*e26+e37*e22*e28-e21*e33*e23-e21*e38*e28-e21*e35*e25+0.5*e31*e222+0.5*e31*e242+0.5*e31*e272-0.5*e31*e232-0.5*e31*e262-0.5*e31*e282-0.5*e31*e252+e21*e30*e20+e21*e32*e22+e24*e30*e23+e24*e20*e33+e24*e21*e34+e24*e32*e25+e24*e22*e35+e34*e20*e23+e34*e22*e25+e27*e30*e26+e27*e20*e36+e27*e21*e37+e27*e32*e28+e27*e22*e38+e37*e20*e26+1.5*e31*e212+0.5*e31*e202;
	Rc(3,ord[16])=e04*e32*e35+0.5*e01*e372-0.5*e01*e352-0.5*e01*e362-0.5*e01*e332-0.5*e01*e382+e04*e31*e34+e34*e30*e03+e34*e00*e33+e34*e32*e05+e34*e02*e35+e07*e30*e36+e07*e32*e38+e07*e31*e37+e37*e30*e06+e37*e00*e36+e37*e32*e08+e04*e30*e33+e37*e02*e38-e31*e33*e03-e31*e36*e06-e31*e35*e05-e31*e38*e08+0.5*e01*e302+0.5*e01*e322+1.5*e01*e312+0.5*e01*e342+e31*e30*e00+e31*e32*e02;
	Rc(3,ord[17])=e31*e32*e12+e14*e30*e33+e14*e32*e35+e14*e31*e34+e34*e30*e13+e34*e10*e33+e34*e32*e15+e34*e12*e35+e17*e30*e36+e17*e32*e38+e17*e31*e37+e37*e30*e16+e37*e10*e36+e37*e32*e18+e37*e12*e38-e31*e36*e16-e31*e35*e15+0.5*e11*e302+0.5*e11*e322+1.5*e11*e312-e31*e33*e13-e31*e38*e18+0.5*e11*e342+0.5*e11*e372+e31*e30*e10-0.5*e11*e352-0.5*e11*e362-0.5*e11*e332-0.5*e11*e382;
	Rc(3,ord[18])=e34*e32*e25+0.5*e21*e342+0.5*e21*e372-0.5*e21*e352-0.5*e21*e362-0.5*e21*e332-0.5*e21*e382+0.5*e21*e302+0.5*e21*e322+1.5*e21*e312+e31*e30*e20+e31*e32*e22+e24*e30*e33+e24*e32*e35+e24*e31*e34+e34*e30*e23+e34*e20*e33+e34*e22*e35+e27*e30*e36+e27*e32*e38+e27*e31*e37+e37*e30*e26+e37*e20*e36+e37*e32*e28+e37*e22*e38-e31*e33*e23-e31*e36*e26-e31*e38*e28-e31*e35*e25;
	Rc(3,ord[19])=0.5*e31*e302+0.5*e31*e322+0.5*e313+e34*e30*e33+e34*e32*e35+0.5*e31*e342+e37*e30*e36+e37*e32*e38+0.5*e31*e372-0.5*e31*e352-0.5*e31*e362-0.5*e31*e332-0.5*e31*e382;
	Rc(4,ord[0])=e01*e00*e03+0.5*e012*e04+e01*e02*e05+0.5*e04*e032+0.5*e043+0.5*e04*e052+e07*e03*e06+0.5*e04*e072+e07*e05*e08-0.5*e04*e002-0.5*e04*e022-0.5*e04*e062-0.5*e04*e082;
	Rc(4,ord[1])=e07*e13*e06+e01*e10*e03-0.5*e14*e002-0.5*e14*e022-0.5*e14*e062-0.5*e14*e082+e01*e00*e13+e01*e11*e04+e01*e12*e05+e01*e02*e15+e11*e00*e03+e11*e02*e05+e04*e13*e03+e04*e15*e05+e07*e03*e16+e07*e04*e17+e07*e15*e08+e07*e05*e18+e17*e03*e06+e17*e05*e08-e04*e10*e00-e04*e12*e02-e04*e16*e06-e04*e18*e08+0.5*e012*e14+1.5*e14*e042+0.5*e14*e032+0.5*e14*e052+0.5*e14*e072;
	Rc(4,ord[2])=e11*e10*e03+0.5*e112*e04+0.5*e04*e132+0.5*e04*e152+1.5*e04*e142+0.5*e04*e172-0.5*e04*e102-0.5*e04*e162-0.5*e04*e122-0.5*e04*e182+e01*e10*e13+e01*e12*e15+e01*e11*e14+e11*e00*e13+e11*e12*e05+e11*e02*e15+e14*e13*e03+e14*e15*e05+e07*e13*e16+e07*e15*e18+e07*e14*e17+e17*e13*e06+e17*e03*e16+e17*e15*e08+e17*e05*e18-e14*e10*e00-e14*e12*e02-e14*e16*e06-e14*e18*e08;
	Rc(4,ord[3])=e11*e10*e13+e11*e12*e15+0.5*e112*e14+0.5*e14*e132+0.5*e14*e152+0.5*e143+e17*e13*e16+e17*e15*e18+0.5*e14*e172-0.5*e14*e102-0.5*e14*e162-0.5*e14*e122-0.5*e14*e182;
	Rc(4,ord[4])=-e04*e28*e08-e04*e26*e06-e04*e22*e02-e04*e20*e00+e27*e05*e08+e27*e03*e06+e07*e05*e28+e07*e25*e08+e07*e04*e27+e07*e03*e26+e07*e23*e06+e04*e25*e05+e04*e23*e03+e21*e02*e05+e21*e00*e03+e01*e02*e25+e01*e22*e05+e01*e21*e04+e01*e00*e23+e01*e20*e03+0.5*e012*e24+0.5*e24*e072+0.5*e24*e052+0.5*e24*e032+1.5*e24*e042-0.5*e24*e022-0.5*e24*e002-0.5*e24*e082-0.5*e24*e062;
	Rc(4,ord[5])=e11*e02*e25+e11*e22*e05+e11*e21*e04-e24*e18*e08-e24*e16*e06-e24*e12*e02-e14*e28*e08-e14*e26*e06-e14*e22*e02-e14*e20*e00-e04*e28*e18-e04*e22*e12-e04*e26*e16-e04*e20*e10+e27*e05*e18+e27*e15*e08+e27*e03*e16+e27*e13*e06+e17*e05*e28+e17*e25*e08+e17*e04*e27+e17*e03*e26+e17*e23*e06+e07*e14*e27+e07*e24*e17+e07*e15*e28+e07*e25*e18+e07*e13*e26+e07*e23*e16+e24*e15*e05+e24*e13*e03+e14*e25*e05+e14*e23*e03+3*e04*e24*e14+e04*e25*e15+e04*e23*e13+e21*e02*e15+e21*e12*e05+e21*e00*e13+e21*e10*e03-e24*e10*e00+e01*e20*e13+e01*e10*e23+e01*e22*e15+e01*e12*e25+e01*e11*e24+e11*e20*e03+e11*e00*e23+e01*e21*e14;
	Rc(4,ord[6])=e11*e12*e25-0.5*e24*e182-0.5*e24*e102-0.5*e24*e162-0.5*e24*e122+e27*e13*e16+e27*e15*e18-e14*e20*e10-e14*e26*e16-e14*e22*e12-e14*e28*e18+e17*e15*e28+e17*e14*e27+1.5*e24*e142+0.5*e24*e132+0.5*e24*e152+0.5*e112*e24+0.5*e24*e172+e11*e21*e14+e11*e22*e15+e11*e10*e23+e11*e20*e13+e21*e10*e13+e21*e12*e15+e17*e13*e26+e17*e23*e16+e14*e25*e15+e14*e23*e13+e17*e25*e18;
	Rc(4,ord[7])=e27*e03*e26+e27*e25*e08+e27*e05*e28-e24*e20*e00-e24*e22*e02-e24*e26*e06-e24*e28*e08+0.5*e04*e232+1.5*e04*e242+0.5*e04*e252+0.5*e04*e272-0.5*e04*e202-0.5*e04*e222-0.5*e04*e262-0.5*e04*e282+e24*e23*e03+e24*e25*e05+e07*e23*e26+e07*e24*e27+e07*e25*e28+e27*e23*e06+e21*e20*e03+e21*e00*e23+e21*e22*e05+e21*e02*e25+0.5*e212*e04+e01*e20*e23+e01*e21*e24+e01*e22*e25;
	Rc(4,ord[8])=-e24*e22*e12-e24*e28*e18+0.5*e14*e272-0.5*e14*e202-0.5*e14*e222-0.5*e14*e262-0.5*e14*e282+e17*e23*e26+e17*e24*e27+e17*e25*e28+e27*e23*e16+e27*e13*e26+e27*e25*e18+e27*e15*e28-e24*e20*e10-e24*e26*e16+0.5*e14*e232+1.5*e14*e242+0.5*e14*e252+e21*e10*e23+e21*e22*e15+e21*e12*e25+e24*e23*e13+e24*e25*e15+e21*e20*e13+0.5*e212*e14+e11*e20*e23+e11*e22*e25+e11*e21*e24;
	Rc(4,ord[9])=e21*e20*e23+0.5*e212*e24+e21*e22*e25+0.5*e24*e232+0.5*e243+0.5*e24*e252+e27*e23*e26+0.5*e24*e272+e27*e25*e28-0.5*e24*e202-0.5*e24*e222-0.5*e24*e262-0.5*e24*e282;
	Rc(4,ord[10])=-e04*e38*e08-e04*e32*e02-e04*e36*e06-0.5*e34*e002-0.5*e34*e022-0.5*e34*e062-0.5*e34*e082+e37*e03*e06+e37*e05*e08-e04*e30*e00+0.5*e34*e032+0.5*e34*e052+0.5*e34*e072+1.5*e34*e042+e01*e30*e03+e01*e00*e33+e01*e31*e04+e01*e32*e05+e01*e02*e35+e31*e00*e03+e31*e02*e05+e04*e33*e03+e04*e35*e05+e07*e03*e36+e07*e04*e37+e07*e35*e08+e07*e05*e38+0.5*e012*e34+e07*e33*e06;
	Rc(4,ord[11])=e07*e13*e36+e01*e12*e35-e04*e30*e10+e17*e04*e37+e17*e35*e08+e17*e05*e38+e37*e13*e06+e37*e03*e16+e37*e15*e08+e37*e05*e18-e04*e36*e16+e17*e33*e06+e04*e33*e13+e04*e35*e15+3*e04*e34*e14+e14*e33*e03+e14*e35*e05+e34*e13*e03+e34*e15*e05+e07*e33*e16+e07*e35*e18+e07*e15*e38+e07*e34*e17+e07*e14*e37+e17*e03*e36+e31*e10*e03+e01*e30*e13+e01*e10*e33+e01*e32*e15+e01*e31*e14+e01*e11*e34+e11*e30*e03+e11*e00*e33+e11*e31*e04+e11*e32*e05+e11*e02*e35+e31*e00*e13+e31*e12*e05+e31*e02*e15-e34*e12*e02-e34*e16*e06-e34*e18*e08-e14*e32*e02-e14*e36*e06-e14*e38*e08-e34*e10*e00-e04*e32*e12-e04*e38*e18-e14*e30*e00;
	Rc(4,ord[12])=e11*e32*e15-0.5*e34*e102-0.5*e34*e162-0.5*e34*e122-0.5*e34*e182+e37*e13*e16+0.5*e112*e34+1.5*e34*e142+0.5*e34*e132+0.5*e34*e152+0.5*e34*e172+e11*e30*e13+e11*e10*e33+e11*e12*e35+e11*e31*e14+e31*e10*e13+e31*e12*e15+e14*e33*e13+e14*e35*e15+e17*e33*e16+e17*e13*e36+e17*e35*e18+e17*e15*e38+e17*e14*e37+e37*e15*e18-e14*e30*e10-e14*e36*e16-e14*e32*e12-e14*e38*e18;
	Rc(4,ord[13])=e01*e22*e35-e04*e30*e20-e04*e32*e22+e01*e31*e24+e01*e21*e34+e01*e32*e25+e21*e30*e03+e21*e00*e33+e21*e31*e04+e21*e32*e05+e21*e02*e35+e31*e20*e03+e31*e00*e23+e31*e22*e05+e31*e02*e25+e04*e33*e23+3*e04*e34*e24+e04*e35*e25+e24*e33*e03+e37*e05*e28-e04*e36*e26-e04*e38*e28-e24*e30*e00-e24*e32*e02-e24*e36*e06-e24*e38*e08+e24*e35*e05+e34*e23*e03+e34*e25*e05+e07*e33*e26+e07*e23*e36+e07*e34*e27+e07*e24*e37+e07*e35*e28+e07*e25*e38+e27*e33*e06+e27*e03*e36+e27*e04*e37+e27*e35*e08+e27*e05*e38+e37*e23*e06+e37*e03*e26+e37*e25*e08-e34*e20*e00-e34*e22*e02-e34*e26*e06-e34*e28*e08+e01*e30*e23+e01*e20*e33;
	Rc(4,ord[14])=e21*e10*e33+e11*e30*e23+e11*e20*e33+e11*e31*e24+e11*e21*e34+e11*e32*e25+e11*e22*e35+e21*e30*e13+e21*e32*e15+e21*e12*e35+e21*e31*e14+e31*e20*e13+e31*e10*e23+e31*e22*e15+e31*e12*e25+e14*e33*e23+3*e14*e34*e24+e14*e35*e25+e24*e33*e13+e24*e35*e15+e34*e23*e13+e34*e25*e15+e17*e33*e26+e17*e23*e36+e17*e34*e27+e17*e24*e37+e17*e35*e28+e17*e25*e38+e27*e33*e16+e27*e13*e36+e27*e35*e18+e27*e15*e38+e27*e14*e37+e37*e23*e16+e37*e13*e26+e37*e25*e18+e37*e15*e28-e34*e28*e18-e34*e22*e12-e14*e32*e22-e14*e36*e26-e14*e38*e28-e24*e30*e10-e24*e36*e16-e24*e32*e12-e24*e38*e18-e34*e20*e10-e34*e26*e16-e14*e30*e20;
	Rc(4,ord[15])=-0.5*e34*e202-0.5*e34*e222-0.5*e34*e262-0.5*e34*e282+e37*e25*e28-e24*e32*e22-e24*e36*e26-e24*e38*e28-e24*e30*e20+0.5*e212*e34+1.5*e34*e242+0.5*e34*e232+0.5*e34*e252+0.5*e34*e272+e21*e30*e23+e21*e20*e33+e21*e31*e24+e21*e32*e25+e21*e22*e35+e31*e20*e23+e31*e22*e25+e24*e33*e23+e24*e35*e25+e27*e33*e26+e27*e23*e36+e27*e24*e37+e27*e35*e28+e27*e25*e38+e37*e23*e26;
	Rc(4,ord[16])=e37*e33*e06+e01*e30*e33+e01*e31*e34+e31*e30*e03+e31*e02*e35+e34*e33*e03+e34*e35*e05+e07*e33*e36+e07*e35*e38+e07*e34*e37+e37*e03*e36+e37*e35*e08+e37*e05*e38-e34*e32*e02-e34*e36*e06-e34*e38*e08+e31*e32*e05+e31*e00*e33+0.5*e312*e04+0.5*e04*e332+0.5*e04*e352+1.5*e04*e342+0.5*e04*e372+e01*e32*e35-0.5*e04*e302-0.5*e04*e322-0.5*e04*e362-0.5*e04*e382-e34*e30*e00;
	Rc(4,ord[17])=0.5*e14*e372-0.5*e14*e302-0.5*e14*e322-0.5*e14*e362-0.5*e14*e382+0.5*e312*e14+0.5*e14*e332+0.5*e14*e352+1.5*e14*e342+e11*e30*e33+e11*e32*e35+e11*e31*e34+e31*e30*e13+e31*e10*e33+e31*e32*e15+e31*e12*e35+e34*e33*e13+e34*e35*e15+e17*e33*e36+e17*e35*e38+e17*e34*e37+e37*e33*e16+e37*e13*e36+e37*e35*e18+e37*e15*e38-e34*e30*e10-e34*e36*e16-e34*e32*e12-e34*e38*e18;
	Rc(4,ord[18])=-e34*e32*e22-e34*e36*e26-e34*e38*e28+0.5*e24*e332+0.5*e24*e352+1.5*e24*e342+0.5*e24*e372-0.5*e24*e302-0.5*e24*e322-0.5*e24*e362-0.5*e24*e382+e21*e30*e33+0.5*e312*e24+e21*e32*e35+e21*e31*e34+e31*e30*e23+e31*e20*e33+e31*e32*e25+e31*e22*e35+e34*e33*e23+e34*e35*e25+e27*e33*e36+e27*e35*e38+e27*e34*e37+e37*e33*e26+e37*e23*e36+e37*e35*e28+e37*e25*e38-e34*e30*e20;
	Rc(4,ord[19])=e31*e30*e33+e31*e32*e35+0.5*e312*e34+0.5*e34*e332+0.5*e34*e352+0.5*e343+e37*e33*e36+e37*e35*e38+0.5*e34*e372-0.5*e34*e302-0.5*e34*e322-0.5*e34*e362-0.5*e34*e382;
	Rc(5,ord[0])=e01*e00*e06+0.5*e012*e07+e01*e02*e08+e04*e03*e06+0.5*e042*e07+e04*e05*e08+0.5*e07*e062+0.5*e073+0.5*e07*e082-0.5*e07*e002-0.5*e07*e022-0.5*e07*e032-0.5*e07*e052;
	Rc(5,ord[1])=e04*e13*e06+0.5*e042*e17+1.5*e17*e072+0.5*e17*e062+0.5*e17*e082-0.5*e17*e002-0.5*e17*e022-0.5*e17*e032-0.5*e17*e052+e01*e10*e06+e07*e16*e06+e07*e18*e08-e07*e10*e00-e07*e12*e02-e07*e13*e03-e07*e15*e05+e01*e00*e16+e01*e11*e07+e01*e12*e08+e01*e02*e18+e11*e00*e06+e11*e02*e08+e04*e03*e16+e04*e14*e07+e04*e15*e08+e04*e05*e18+e14*e03*e06+e14*e05*e08+0.5*e012*e17;
	Rc(5,ord[2])=-e17*e10*e00+0.5*e112*e07+0.5*e142*e07+0.5*e07*e162+0.5*e07*e182+1.5*e07*e172-0.5*e07*e102-0.5*e07*e152-0.5*e07*e132-0.5*e07*e122+e01*e10*e16+e01*e12*e18+e01*e11*e17+e11*e10*e06+e11*e00*e16+e11*e12*e08+e11*e02*e18+e04*e13*e16+e04*e15*e18+e04*e14*e17+e14*e13*e06+e14*e03*e16+e14*e15*e08+e14*e05*e18+e17*e16*e06+e17*e18*e08-e17*e12*e02-e17*e13*e03-e17*e15*e05;
	Rc(5,ord[3])=e11*e10*e16+e11*e12*e18+0.5*e112*e17+e14*e13*e16+e14*e15*e18+0.5*e142*e17+0.5*e17*e162+0.5*e17*e182+0.5*e173-0.5*e17*e102-0.5*e17*e152-0.5*e17*e132-0.5*e17*e122;
	Rc(5,ord[4])=e01*e22*e08+e07*e28*e08-e07*e20*e00-e07*e23*e03-e07*e22*e02-e07*e25*e05+0.5*e012*e27+0.5*e042*e27+1.5*e27*e072+0.5*e27*e062+0.5*e27*e082-0.5*e27*e002-0.5*e27*e022-0.5*e27*e032-0.5*e27*e052+e07*e26*e06+e01*e20*e06+e01*e00*e26+e01*e21*e07+e01*e02*e28+e21*e00*e06+e21*e02*e08+e04*e23*e06+e04*e03*e26+e04*e24*e07+e04*e25*e08+e04*e05*e28+e24*e03*e06+e24*e05*e08;
	Rc(5,ord[5])=e14*e24*e07+e14*e03*e26+e14*e23*e06+e04*e14*e27+e04*e24*e17+e04*e15*e28+e04*e25*e18+e04*e13*e26+e04*e23*e16+e21*e02*e18+e21*e12*e08-e27*e15*e05-e27*e13*e03-e27*e12*e02-e27*e10*e00-e17*e25*e05-e17*e23*e03-e17*e22*e02-e17*e20*e00-e07*e22*e12-e07*e23*e13-e07*e25*e15-e07*e20*e10+e27*e18*e08+e27*e16*e06+e17*e28*e08+e17*e26*e06+3*e07*e27*e17+e07*e28*e18+e07*e26*e16+e24*e05*e18+e24*e15*e08+e24*e03*e16+e24*e13*e06+e14*e05*e28+e14*e25*e08+e01*e12*e28+e01*e20*e16+e01*e10*e26+e01*e22*e18+e01*e21*e17+e11*e20*e06+e01*e11*e27+e21*e00*e16+e21*e10*e06+e11*e21*e07+e11*e22*e08+e11*e02*e28+e11*e00*e26;
	Rc(5,ord[6])=-0.5*e27*e102-0.5*e27*e152-0.5*e27*e132-0.5*e27*e122+0.5*e142*e27+1.5*e27*e172+0.5*e27*e162+0.5*e27*e182+0.5*e112*e27+e11*e22*e18+e11*e10*e26+e11*e20*e16-e17*e22*e12-e17*e23*e13-e17*e25*e15-e17*e20*e10+e17*e28*e18+e17*e26*e16+e24*e15*e18+e24*e13*e16+e14*e24*e17+e14*e15*e28+e14*e25*e18+e14*e13*e26+e14*e23*e16+e21*e12*e18+e21*e10*e16+e11*e21*e17+e11*e12*e28;
	Rc(5,ord[7])=-0.5*e07*e252+e27*e26*e06+e27*e28*e08-e27*e20*e00-e27*e22*e02-e27*e23*e03-e27*e25*e05+1.5*e07*e272+0.5*e07*e282+e01*e22*e28+e21*e20*e06+e21*e00*e26+e21*e22*e08+e21*e02*e28+e04*e23*e26+e04*e24*e27+e04*e25*e28+e24*e23*e06+e24*e03*e26+e24*e25*e08+e24*e05*e28+0.5*e212*e07+0.5*e242*e07+0.5*e07*e262+e01*e20*e26+e01*e21*e27-0.5*e07*e202-0.5*e07*e232-0.5*e07*e222;
	Rc(5,ord[8])=-e27*e25*e15-e27*e23*e13-e27*e22*e12-0.5*e17*e252-0.5*e17*e202-0.5*e17*e222-0.5*e17*e232+0.5*e17*e262+1.5*e17*e272+0.5*e17*e282+e24*e23*e16+e24*e13*e26+e24*e25*e18+e24*e15*e28+e27*e26*e16+e27*e28*e18-e27*e20*e10+e14*e24*e27+e14*e25*e28+0.5*e212*e17+0.5*e242*e17+e11*e20*e26+e11*e21*e27+e11*e22*e28+e21*e20*e16+e21*e10*e26+e21*e22*e18+e21*e12*e28+e14*e23*e26;
	Rc(5,ord[9])=e21*e20*e26+0.5*e212*e27+e21*e22*e28+e24*e23*e26+0.5*e242*e27+e24*e25*e28+0.5*e27*e262+0.5*e273+0.5*e27*e282-0.5*e27*e202-0.5*e27*e222-0.5*e27*e232-0.5*e27*e252;
	Rc(5,ord[10])=e04*e05*e38+e01*e30*e06-0.5*e37*e002-0.5*e37*e022-0.5*e37*e032-0.5*e37*e052-e07*e32*e02-e07*e35*e05-e07*e33*e03+e07*e36*e06+e07*e38*e08-e07*e30*e00+1.5*e37*e072+0.5*e37*e062+0.5*e37*e082+e01*e02*e38+e31*e00*e06+e31*e02*e08+e04*e33*e06+e04*e03*e36+e04*e34*e07+e04*e35*e08+e34*e03*e06+e34*e05*e08+0.5*e012*e37+0.5*e042*e37+e01*e00*e36+e01*e31*e07+e01*e32*e08;
	Rc(5,ord[11])=e14*e33*e06+e11*e30*e06+e11*e00*e36+e11*e31*e07+e31*e10*e06+e11*e32*e08+e11*e02*e38+e31*e00*e16+e31*e12*e08+e31*e02*e18+e04*e33*e16+e04*e13*e36+e04*e35*e18+e04*e15*e38+e01*e10*e36+e01*e32*e18+e01*e12*e38+e01*e31*e17+e01*e11*e37+e01*e30*e16-e17*e35*e05-e37*e10*e00-e37*e12*e02-e37*e13*e03-e37*e15*e05+e37*e18*e08-e07*e30*e10-e07*e35*e15-e07*e33*e13-e07*e32*e12-e17*e30*e00-e17*e32*e02-e17*e33*e03+e07*e38*e18+3*e07*e37*e17+e17*e36*e06+e17*e38*e08+e37*e16*e06+e04*e34*e17+e04*e14*e37+e14*e03*e36+e14*e34*e07+e14*e35*e08+e14*e05*e38+e34*e13*e06+e34*e03*e16+e34*e15*e08+e34*e05*e18+e07*e36*e16;
	Rc(5,ord[12])=e11*e32*e18-0.5*e37*e102-0.5*e37*e152-0.5*e37*e132-0.5*e37*e122+0.5*e112*e37+0.5*e142*e37+1.5*e37*e172+0.5*e37*e162+0.5*e37*e182+e11*e10*e36+e11*e12*e38+e11*e31*e17+e31*e10*e16+e31*e12*e18+e14*e33*e16+e14*e13*e36+e14*e35*e18+e14*e15*e38+e14*e34*e17+e34*e13*e16+e34*e15*e18+e17*e36*e16+e17*e38*e18-e17*e30*e10-e17*e35*e15-e17*e33*e13-e17*e32*e12+e11*e30*e16;
	Rc(5,ord[13])=e01*e20*e36+e01*e31*e27+e01*e21*e37+e01*e32*e28+e01*e22*e38+e21*e30*e06+e21*e00*e36+e21*e31*e07+e21*e32*e08+e21*e02*e38+e01*e30*e26+e31*e20*e06+e31*e00*e26+e31*e22*e08+e31*e02*e28+e04*e33*e26+e04*e23*e36+e04*e34*e27+e04*e24*e37+e04*e35*e28+e04*e25*e38+e24*e33*e06+e24*e03*e36+e24*e34*e07+e24*e35*e08+e24*e05*e38+e34*e23*e06+e34*e03*e26+e34*e25*e08+e34*e05*e28+e07*e36*e26+3*e07*e37*e27+e07*e38*e28+e27*e36*e06+e27*e38*e08+e37*e26*e06+e37*e28*e08-e07*e30*e20-e07*e32*e22-e07*e33*e23-e07*e35*e25-e27*e30*e00-e27*e32*e02-e27*e33*e03-e27*e35*e05-e37*e20*e00-e37*e22*e02-e37*e23*e03-e37*e25*e05;
	Rc(5,ord[14])=e11*e30*e26+e11*e20*e36+e11*e31*e27+e11*e21*e37+e11*e32*e28+e11*e22*e38+e21*e10*e36+e21*e32*e18+e21*e12*e38+e21*e31*e17+e31*e20*e16+e31*e10*e26+e31*e22*e18+e31*e12*e28+e14*e33*e26+e14*e23*e36+e14*e34*e27+e14*e24*e37+e14*e35*e28+e14*e25*e38+e24*e33*e16+e24*e13*e36+e24*e35*e18+e24*e15*e38+e24*e34*e17+e34*e23*e16+e34*e13*e26+e34*e25*e18+e34*e15*e28+e17*e36*e26+3*e17*e37*e27+e17*e38*e28+e27*e36*e16+e27*e38*e18+e37*e26*e16+e37*e28*e18-e17*e30*e20-e17*e32*e22-e17*e33*e23-e17*e35*e25-e27*e30*e10-e27*e35*e15-e27*e33*e13-e27*e32*e12-e37*e20*e10-e37*e25*e15-e37*e23*e13-e37*e22*e12+e21*e30*e16;
	Rc(5,ord[15])=e21*e20*e36+e21*e31*e27+e21*e32*e28+e21*e22*e38+e31*e22*e28+e24*e33*e26+e24*e23*e36+e24*e34*e27+e24*e35*e28+e24*e25*e38+e34*e23*e26+e34*e25*e28+e27*e36*e26+e27*e38*e28-e27*e30*e20-e27*e32*e22-e27*e33*e23-e27*e35*e25+0.5*e242*e37+1.5*e37*e272+0.5*e37*e262+0.5*e37*e282+e31*e20*e26+e21*e30*e26+0.5*e212*e37-0.5*e37*e202-0.5*e37*e222-0.5*e37*e232-0.5*e37*e252;
	Rc(5,ord[16])=e01*e30*e36+e01*e32*e38+e01*e31*e37+e31*e30*e06+e31*e00*e36+e31*e32*e08+e31*e02*e38+e04*e33*e36+e04*e35*e38+e04*e34*e37+e34*e33*e06+e34*e03*e36+e34*e35*e08+e34*e05*e38+e37*e36*e06+e37*e38*e08-e37*e30*e00-e37*e32*e02-e37*e33*e03-e37*e35*e05+0.5*e312*e07+0.5*e342*e07+0.5*e07*e362+0.5*e07*e382+1.5*e07*e372-0.5*e07*e302-0.5*e07*e352-0.5*e07*e322-0.5*e07*e332;
	Rc(5,ord[17])=0.5*e312*e17+0.5*e342*e17+0.5*e17*e362+0.5*e17*e382+1.5*e17*e372-0.5*e17*e302-0.5*e17*e352-0.5*e17*e322-0.5*e17*e332-e37*e32*e12-e37*e33*e13+e11*e30*e36+e11*e32*e38+e11*e31*e37+e31*e30*e16+e31*e10*e36+e31*e32*e18+e31*e12*e38+e14*e33*e36+e14*e35*e38+e14*e34*e37+e34*e33*e16+e34*e13*e36+e34*e35*e18+e34*e15*e38+e37*e36*e16+e37*e38*e18-e37*e30*e10-e37*e35*e15;
	Rc(5,ord[18])=e21*e31*e37-0.5*e27*e332+e21*e30*e36+e21*e32*e38+e31*e30*e26+e31*e20*e36+e31*e32*e28+e31*e22*e38+e24*e33*e36+e24*e35*e38+e24*e34*e37+e34*e33*e26+e34*e23*e36+e34*e35*e28+e34*e25*e38+e37*e36*e26+e37*e38*e28-e37*e30*e20-e37*e32*e22-e37*e33*e23-e37*e35*e25+0.5*e312*e27+0.5*e342*e27+0.5*e27*e362+0.5*e27*e382+1.5*e27*e372-0.5*e27*e302-0.5*e27*e352-0.5*e27*e322;
	Rc(5,ord[19])=e31*e30*e36+e31*e32*e38+0.5*e312*e37+e34*e33*e36+e34*e35*e38+0.5*e342*e37+0.5*e37*e362+0.5*e37*e382+0.5*e373-0.5*e37*e302-0.5*e37*e352-0.5*e37*e322-0.5*e37*e332;
	Rc(6,ord[0])=0.5*e02*e002+0.5*e02*e012+0.5*e023+e05*e00*e03+e05*e01*e04+0.5*e02*e052+e08*e00*e06+e08*e01*e07+0.5*e02*e082-0.5*e02*e032-0.5*e02*e042-0.5*e02*e062-0.5*e02*e072;
	Rc(6,ord[1])=-0.5*e12*e042-0.5*e12*e062-0.5*e12*e072+0.5*e12*e082-0.5*e12*e032+1.5*e12*e022+0.5*e12*e002+0.5*e12*e012+0.5*e12*e052+e02*e10*e00+e02*e11*e01+e05*e10*e03+e05*e00*e13+e05*e11*e04+e05*e01*e14+e05*e02*e15+e15*e00*e03+e15*e01*e04+e08*e10*e06+e08*e00*e16+e08*e11*e07+e08*e01*e17+e08*e02*e18+e18*e00*e06+e18*e01*e07-e02*e13*e03-e02*e14*e04-e02*e16*e06-e02*e17*e07;
	Rc(6,ord[2])=0.5*e02*e102+1.5*e02*e122+0.5*e02*e112+0.5*e02*e152+0.5*e02*e182-0.5*e02*e162-0.5*e02*e172-0.5*e02*e132-0.5*e02*e142+e12*e10*e00+e12*e11*e01+e05*e10*e13+e05*e12*e15+e05*e11*e14+e15*e10*e03+e15*e00*e13+e15*e11*e04+e15*e01*e14+e08*e10*e16+e08*e12*e18+e08*e11*e17+e18*e10*e06+e18*e00*e16+e18*e11*e07+e18*e01*e17-e12*e13*e03-e12*e14*e04-e12*e16*e06-e12*e17*e07;
	Rc(6,ord[3])=0.5*e12*e102+0.5*e123+0.5*e12*e112+e15*e10*e13+0.5*e12*e152+e15*e11*e14+e18*e10*e16+0.5*e12*e182+e18*e11*e17-0.5*e12*e162-0.5*e12*e172-0.5*e12*e132-0.5*e12*e142;
	Rc(6,ord[4])=-0.5*e22*e032-0.5*e22*e042-0.5*e22*e062-0.5*e22*e072+0.5*e22*e082+1.5*e22*e022+0.5*e22*e002+0.5*e22*e012+0.5*e22*e052+e02*e20*e00+e02*e21*e01+e05*e20*e03+e05*e00*e23+e05*e21*e04+e05*e01*e24+e05*e02*e25+e25*e00*e03+e25*e01*e04+e08*e20*e06+e08*e00*e26+e08*e21*e07+e08*e01*e27+e08*e02*e28+e28*e00*e06+e28*e01*e07-e02*e27*e07-e02*e23*e03-e02*e24*e04-e02*e26*e06;
	Rc(6,ord[5])=-e22*e17*e07-e22*e16*e06-e22*e14*e04-e22*e13*e03-e12*e26*e06-e12*e24*e04-e12*e23*e03-e12*e27*e07-e02*e24*e14-e02*e23*e13-e02*e27*e17-e02*e26*e16+e28*e01*e17+e28*e11*e07+e28*e00*e16+e28*e10*e06+e18*e02*e28+e18*e01*e27+e18*e21*e07+e18*e00*e26+e18*e20*e06+e08*e11*e27+e08*e21*e17+e08*e12*e28+e08*e22*e18+e08*e10*e26+e25*e01*e14+e25*e11*e04+e25*e00*e13+e25*e10*e03+e15*e01*e24+e02*e21*e11+e12*e21*e01+e15*e02*e25+e15*e21*e04+e05*e22*e15+e05*e11*e24+e15*e20*e03+e15*e00*e23+e05*e10*e23+e05*e12*e25+e05*e21*e14+e22*e10*e00+e22*e11*e01+e02*e20*e10+3*e02*e22*e12+e12*e20*e00+e08*e20*e16+e05*e20*e13;
	Rc(6,ord[6])=-e12*e24*e14-e12*e23*e13-e12*e27*e17-e12*e26*e16+e28*e11*e17+e28*e10*e16+e18*e11*e27+e18*e21*e17+e18*e12*e28+e18*e10*e26+e18*e20*e16+e25*e11*e14+e25*e10*e13+e15*e11*e24+e15*e21*e14+e15*e12*e25+e15*e10*e23+e15*e20*e13+e12*e21*e11+0.5*e22*e182+0.5*e22*e152+1.5*e22*e122+0.5*e22*e102+e12*e20*e10+0.5*e22*e112-0.5*e22*e172-0.5*e22*e132-0.5*e22*e142-0.5*e22*e162;
	Rc(6,ord[7])=0.5*e02*e282+e28*e01*e27-e22*e27*e07-e22*e23*e03-e22*e24*e04-e22*e26*e06+0.5*e02*e252+e05*e20*e23+e05*e22*e25+e25*e20*e03+e25*e00*e23+e25*e21*e04+e25*e01*e24+e08*e20*e26+e08*e21*e27+e08*e22*e28+e28*e20*e06+e28*e00*e26+e28*e21*e07+e05*e21*e24+0.5*e02*e202+0.5*e02*e212+1.5*e02*e222+e22*e20*e00+e22*e21*e01-0.5*e02*e272-0.5*e02*e242-0.5*e02*e232-0.5*e02*e262;
	Rc(6,ord[8])=-e22*e27*e17-e22*e23*e13-e22*e24*e14-0.5*e12*e232-0.5*e12*e262-0.5*e12*e242-0.5*e12*e272+0.5*e12*e282+e18*e21*e27+e28*e20*e16+e28*e10*e26+e28*e21*e17+e28*e11*e27-e22*e26*e16+e18*e22*e28+0.5*e12*e252+0.5*e12*e202+0.5*e12*e212+1.5*e12*e222+e22*e20*e10+e15*e20*e23+e15*e21*e24+e15*e22*e25+e25*e20*e13+e25*e10*e23+e25*e21*e14+e25*e11*e24+e18*e20*e26+e22*e21*e11;
	Rc(6,ord[9])=0.5*e22*e202+0.5*e22*e212+0.5*e223+e25*e20*e23+e25*e21*e24+0.5*e22*e252+e28*e20*e26+e28*e21*e27+0.5*e22*e282-0.5*e22*e232-0.5*e22*e262-0.5*e22*e242-0.5*e22*e272;
	Rc(6,ord[10])=e08*e31*e07-0.5*e32*e032-e02*e33*e03-e02*e34*e04-e02*e36*e06-0.5*e32*e042-0.5*e32*e062-0.5*e32*e072+e38*e01*e07+e38*e00*e06-e02*e37*e07+e05*e31*e04+e05*e01*e34+e05*e02*e35+e35*e01*e04+e35*e00*e03+e08*e30*e06+e08*e00*e36+e08*e01*e37+e08*e02*e38+0.5*e32*e052+e02*e30*e00+e02*e31*e01+e05*e30*e03+e05*e00*e33+1.5*e32*e022+0.5*e32*e012+0.5*e32*e002+0.5*e32*e082;
	Rc(6,ord[11])=e05*e32*e15+e32*e11*e01+e38*e10*e06+e08*e12*e38-e32*e14*e04-e32*e16*e06-e32*e17*e07-e12*e36*e06-e32*e13*e03-e02*e34*e14-e12*e37*e07-e12*e33*e03-e12*e34*e04-e02*e37*e17-e02*e33*e13+e38*e01*e17-e02*e36*e16+e18*e01*e37+e18*e02*e38+e38*e00*e16+e38*e11*e07+e08*e30*e16+e08*e10*e36+e08*e32*e18+e08*e31*e17+e08*e11*e37+e18*e30*e06+e18*e00*e36+e18*e31*e07+e35*e10*e03+e35*e00*e13+e35*e11*e04+e35*e01*e14+e15*e02*e35+e05*e10*e33+e05*e12*e35+e05*e31*e14+e05*e11*e34+e15*e30*e03+e15*e00*e33+e15*e31*e04+e15*e01*e34+e05*e30*e13+e02*e30*e10+e02*e31*e11+3*e02*e32*e12+e12*e30*e00+e12*e31*e01+e32*e10*e00;
	Rc(6,ord[12])=0.5*e32*e102+0.5*e32*e112+e12*e30*e10+1.5*e32*e122+e12*e31*e11+e15*e30*e13+e15*e10*e33+e15*e12*e35+e15*e31*e14+e15*e11*e34+e35*e10*e13-0.5*e32*e162-0.5*e32*e172-0.5*e32*e132-0.5*e32*e142-e12*e37*e17-e12*e33*e13-e12*e34*e14+0.5*e32*e182+0.5*e32*e152+e35*e11*e14+e18*e30*e16+e18*e10*e36+e18*e12*e38+e18*e31*e17+e18*e11*e37+e38*e10*e16+e38*e11*e17-e12*e36*e16;
	Rc(6,ord[13])=3*e02*e32*e22+e05*e31*e24+e08*e22*e38+e02*e31*e21+e22*e30*e00+e22*e31*e01+e32*e20*e00+e32*e21*e01+e05*e30*e23+e05*e20*e33+e05*e21*e34-e22*e37*e07-e22*e33*e03-e22*e34*e04-e22*e36*e06-e32*e27*e07-e32*e23*e03-e32*e24*e04-e32*e26*e06+e05*e32*e25+e25*e30*e03+e25*e00*e33+e25*e31*e04+e25*e01*e34+e25*e02*e35+e35*e20*e03+e35*e00*e23+e35*e21*e04+e35*e01*e24+e08*e30*e26+e08*e20*e36+e08*e31*e27+e08*e21*e37+e08*e32*e28+e28*e30*e06+e28*e00*e36+e28*e31*e07+e28*e01*e37+e28*e02*e38+e38*e20*e06+e38*e00*e26+e38*e21*e07+e38*e01*e27-e02*e33*e23-e02*e36*e26-e02*e34*e24-e02*e37*e27+e05*e22*e35+e02*e30*e20;
	Rc(6,ord[14])=e18*e22*e38+e12*e31*e21+3*e12*e32*e22+e22*e30*e10+e22*e31*e11+e32*e20*e10+e32*e21*e11+e15*e30*e23+e15*e20*e33+e15*e31*e24+e15*e21*e34+e15*e32*e25+e15*e22*e35+e25*e30*e13+e25*e10*e33+e25*e12*e35+e25*e31*e14+e25*e11*e34+e35*e20*e13+e35*e10*e23+e35*e21*e14+e35*e11*e24+e18*e30*e26+e18*e20*e36+e18*e31*e27+e18*e21*e37+e18*e32*e28+e28*e30*e16+e28*e10*e36+e28*e12*e38+e28*e31*e17+e28*e11*e37+e38*e20*e16+e38*e10*e26+e12*e30*e20-e22*e37*e17-e22*e33*e13-e22*e34*e14-e32*e26*e16-e32*e27*e17-e32*e23*e13-e32*e24*e14-e22*e36*e16+e38*e21*e17+e38*e11*e27-e12*e33*e23-e12*e36*e26-e12*e34*e24-e12*e37*e27;
	Rc(6,ord[15])=e25*e30*e23+e22*e30*e20+e22*e31*e21+e25*e20*e33+e25*e31*e24+e25*e21*e34+e25*e22*e35+e35*e20*e23+e35*e21*e24+e28*e30*e26+e28*e20*e36+e28*e31*e27+e28*e21*e37+e28*e22*e38+e38*e20*e26+e38*e21*e27-e22*e33*e23-e22*e36*e26-e22*e34*e24-e22*e37*e27+0.5*e32*e212+0.5*e32*e252+1.5*e32*e222+0.5*e32*e202+0.5*e32*e282-0.5*e32*e232-0.5*e32*e262-0.5*e32*e242-0.5*e32*e272;
	Rc(6,ord[16])=0.5*e02*e302+1.5*e02*e322+0.5*e02*e312+0.5*e02*e352+0.5*e02*e382-0.5*e02*e342-0.5*e02*e362-0.5*e02*e332-0.5*e02*e372+e38*e30*e06+e32*e30*e00+e32*e31*e01+e05*e30*e33+e05*e32*e35+e05*e31*e34+e35*e30*e03+e35*e00*e33+e35*e31*e04+e35*e01*e34+e08*e30*e36+e08*e32*e38+e08*e31*e37+e38*e00*e36+e38*e31*e07+e38*e01*e37-e32*e37*e07-e32*e33*e03-e32*e34*e04-e32*e36*e06;
	Rc(6,ord[17])=e32*e30*e10+e32*e31*e11+e15*e30*e33+e15*e32*e35+e15*e31*e34+e35*e30*e13+e35*e10*e33+e35*e31*e14+e35*e11*e34+e18*e30*e36+e18*e32*e38+e18*e31*e37+e38*e30*e16+e38*e10*e36+e38*e31*e17+e38*e11*e37-e32*e36*e16-e32*e37*e17-e32*e33*e13-e32*e34*e14+0.5*e12*e382-0.5*e12*e342-0.5*e12*e362-0.5*e12*e332-0.5*e12*e372+0.5*e12*e352+1.5*e12*e322+0.5*e12*e312+0.5*e12*e302;
	Rc(6,ord[18])=0.5*e22*e302+0.5*e22*e312+0.5*e22*e352+0.5*e22*e382-0.5*e22*e342-0.5*e22*e362-0.5*e22*e332-0.5*e22*e372+1.5*e22*e322+e32*e30*e20+e32*e31*e21+e25*e30*e33+e25*e32*e35+e25*e31*e34+e35*e30*e23+e35*e20*e33+e35*e31*e24+e35*e21*e34+e28*e30*e36+e28*e32*e38+e28*e31*e37+e38*e30*e26+e38*e20*e36+e38*e31*e27+e38*e21*e37-e32*e33*e23-e32*e36*e26-e32*e34*e24-e32*e37*e27;
	Rc(6,ord[19])=0.5*e32*e302+0.5*e323+0.5*e32*e312+e35*e30*e33+0.5*e32*e352+e35*e31*e34+e38*e30*e36+0.5*e32*e382+e38*e31*e37-0.5*e32*e342-0.5*e32*e362-0.5*e32*e332-0.5*e32*e372;
	Rc(7,ord[0])=e02*e01*e04+e02*e00*e03+0.5*e022*e05+0.5*e05*e032+0.5*e05*e042+0.5*e053+e08*e03*e06+e08*e04*e07+0.5*e05*e082-0.5*e05*e002-0.5*e05*e062-0.5*e05*e012-0.5*e05*e072;
	Rc(7,ord[1])=e08*e13*e06+e02*e10*e03+e02*e00*e13+e02*e11*e04+e02*e01*e14+e02*e12*e05+e12*e01*e04+e12*e00*e03+e05*e13*e03+e05*e14*e04+e08*e03*e16+e08*e14*e07+e08*e04*e17+e08*e05*e18+e18*e03*e06+e18*e04*e07-e05*e10*e00-e05*e11*e01-e05*e16*e06-e05*e17*e07+0.5*e022*e15+1.5*e15*e052+0.5*e15*e032+0.5*e15*e042+0.5*e15*e082-0.5*e15*e002-0.5*e15*e062-0.5*e15*e012-0.5*e15*e072;
	Rc(7,ord[2])=0.5*e122*e05+0.5*e05*e132+1.5*e05*e152+0.5*e05*e142+0.5*e05*e182-0.5*e05*e102-0.5*e05*e162-0.5*e05*e112-0.5*e05*e172+e02*e10*e13+e02*e12*e15+e02*e11*e14+e12*e10*e03+e12*e00*e13+e12*e11*e04+e12*e01*e14+e15*e13*e03+e15*e14*e04+e08*e13*e16+e08*e15*e18+e08*e14*e17+e18*e13*e06+e18*e03*e16+e18*e14*e07+e18*e04*e17-e15*e11*e01-e15*e16*e06-e15*e17*e07-e15*e10*e00;
	Rc(7,ord[3])=e12*e10*e13+0.5*e122*e15+e12*e11*e14+0.5*e15*e132+0.5*e153+0.5*e15*e142+e18*e13*e16+0.5*e15*e182+e18*e14*e17-0.5*e15*e102-0.5*e15*e162-0.5*e15*e112-0.5*e15*e172;
	Rc(7,ord[4])=0.5*e25*e082-0.5*e25*e002-0.5*e25*e062-0.5*e25*e012-0.5*e25*e072+e02*e20*e03+e02*e00*e23+e02*e21*e04+e02*e01*e24+e02*e22*e05+e22*e01*e04+e22*e00*e03+e05*e23*e03+e05*e24*e04+e08*e23*e06+e08*e03*e26+e08*e24*e07+e08*e04*e27+e08*e05*e28+e28*e03*e06+e28*e04*e07-e05*e20*e00-e05*e27*e07-e05*e21*e01-e05*e26*e06+0.5*e022*e25+1.5*e25*e052+0.5*e25*e032+0.5*e25*e042;
	Rc(7,ord[5])=-e25*e17*e07-e25*e16*e06-e25*e11*e01-e25*e10*e00-e15*e26*e06-e15*e21*e01-e15*e27*e07-e15*e20*e00-e05*e27*e17-e05*e21*e11-e05*e26*e16-e05*e20*e10+e28*e04*e17+e28*e14*e07+e28*e03*e16+e28*e13*e06+e18*e05*e28+e18*e04*e27+e18*e24*e07+e18*e03*e26+e18*e23*e06+e08*e14*e27+e08*e24*e17+e08*e15*e28+e08*e25*e18+e08*e13*e26+e08*e23*e16+e25*e14*e04+e25*e13*e03+e15*e24*e04+e15*e23*e03+e05*e24*e14+3*e05*e25*e15+e05*e23*e13+e22*e01*e14+e22*e11*e04+e22*e00*e13+e22*e10*e03+e12*e22*e05+e12*e01*e24+e12*e21*e04+e12*e00*e23+e12*e20*e03+e02*e11*e24+e02*e21*e14+e02*e12*e25+e02*e22*e15+e02*e10*e23+e02*e20*e13;
	Rc(7,ord[6])=-e15*e27*e17-e15*e21*e11-e15*e26*e16+e28*e14*e17+e28*e13*e16+e18*e14*e27+e18*e24*e17+e18*e15*e28+e18*e13*e26+e15*e24*e14+e15*e23*e13+e22*e11*e14+e22*e10*e13+e12*e11*e24+e12*e21*e14+e12*e22*e15+e12*e10*e23+e18*e23*e16+0.5*e25*e142+0.5*e25*e182+1.5*e25*e152+0.5*e25*e132+0.5*e122*e25+e12*e20*e13-0.5*e25*e172-0.5*e25*e162-0.5*e25*e112-0.5*e25*e102-e15*e20*e10;
	Rc(7,ord[7])=e28*e24*e07-0.5*e05*e272-0.5*e05*e262-0.5*e05*e212+0.5*e05*e282-0.5*e05*e202+e28*e23*e06+e08*e23*e26+e08*e25*e28+e08*e24*e27+e28*e03*e26+e28*e04*e27-e25*e27*e07-e25*e21*e01-e25*e26*e06+e02*e20*e23+e02*e22*e25+e02*e21*e24+e22*e20*e03+e22*e00*e23+e22*e21*e04+e22*e01*e24+e25*e23*e03+e25*e24*e04+0.5*e222*e05+0.5*e05*e232+1.5*e05*e252+0.5*e05*e242-e25*e20*e00;
	Rc(7,ord[8])=-0.5*e15*e202-0.5*e15*e262-0.5*e15*e212-0.5*e15*e272+e18*e23*e26+e18*e25*e28+e18*e24*e27+e28*e23*e16+e28*e13*e26+e28*e24*e17+e28*e14*e27-e25*e20*e10-e25*e26*e16-e25*e21*e11-e25*e27*e17+0.5*e15*e282+0.5*e15*e232+1.5*e15*e252+0.5*e15*e242+0.5*e222*e15+e12*e21*e24+e22*e20*e13+e22*e10*e23+e22*e21*e14+e22*e11*e24+e25*e23*e13+e25*e24*e14+e12*e20*e23+e12*e22*e25;
	Rc(7,ord[9])=e22*e20*e23+0.5*e222*e25+e22*e21*e24+0.5*e25*e232+0.5*e253+0.5*e25*e242+e28*e23*e26+0.5*e25*e282+e28*e24*e27-0.5*e25*e202-0.5*e25*e262-0.5*e25*e212-0.5*e25*e272;
	Rc(7,ord[10])=-0.5*e35*e062-0.5*e35*e012-0.5*e35*e072-e05*e30*e00-e05*e31*e01-e05*e36*e06-e05*e37*e07-0.5*e35*e002+0.5*e35*e082+e05*e34*e04+e08*e33*e06+e08*e03*e36+e08*e34*e07+e08*e04*e37+e08*e05*e38+e38*e04*e07+e38*e03*e06+0.5*e022*e35+1.5*e35*e052+0.5*e35*e042+0.5*e35*e032+e02*e30*e03+e02*e00*e33+e02*e31*e04+e02*e01*e34+e02*e32*e05+e32*e01*e04+e32*e00*e03+e05*e33*e03;
	Rc(7,ord[11])=e08*e33*e16-e35*e16*e06-e35*e17*e07-e15*e30*e00-e15*e37*e07-e15*e31*e01-e15*e36*e06-e35*e10*e00-e35*e11*e01-e05*e37*e17-e05*e31*e11+e38*e04*e17-e05*e30*e10-e05*e36*e16+e18*e33*e06+e18*e03*e36+e18*e34*e07+e18*e04*e37+e18*e05*e38+e38*e13*e06+e38*e03*e16+e38*e14*e07+e35*e14*e04+e08*e13*e36+e08*e35*e18+e08*e15*e38+e08*e34*e17+e08*e14*e37+e35*e13*e03+e05*e33*e13+3*e05*e35*e15+e05*e34*e14+e15*e33*e03+e15*e34*e04+e12*e01*e34+e12*e32*e05+e32*e10*e03+e32*e00*e13+e32*e11*e04+e32*e01*e14+e12*e30*e03+e02*e30*e13+e02*e32*e15+e02*e10*e33+e02*e12*e35+e12*e00*e33+e02*e31*e14+e02*e11*e34+e12*e31*e04;
	Rc(7,ord[12])=-0.5*e35*e162-0.5*e35*e172-e15*e36*e16-e15*e31*e11-e15*e37*e17-0.5*e35*e102-0.5*e35*e112-e15*e30*e10+e18*e13*e36+e18*e15*e38+e18*e34*e17+e18*e14*e37+e38*e13*e16+e38*e14*e17+e18*e33*e16+1.5*e35*e152+0.5*e35*e132+0.5*e35*e142+0.5*e35*e182+0.5*e122*e35+e32*e10*e13+e32*e11*e14+e15*e33*e13+e15*e34*e14+e12*e10*e33+e12*e32*e15+e12*e31*e14+e12*e11*e34+e12*e30*e13;
	Rc(7,ord[13])=e05*e33*e23+3*e05*e35*e25+e05*e34*e24+e25*e33*e03+e25*e34*e04+e35*e23*e03+e35*e24*e04+e08*e33*e26+e08*e23*e36+e08*e35*e28+e02*e20*e33+e02*e32*e25+e02*e22*e35+e02*e31*e24+e02*e21*e34+e22*e30*e03+e22*e00*e33+e22*e31*e04+e22*e01*e34+e22*e32*e05+e32*e20*e03+e32*e00*e23+e32*e21*e04+e32*e01*e24+e02*e30*e23-e35*e27*e07-e35*e21*e01-e35*e26*e06+e08*e25*e38+e08*e34*e27+e08*e24*e37+e28*e33*e06+e28*e03*e36+e28*e34*e07+e28*e04*e37+e28*e05*e38+e38*e23*e06+e38*e03*e26+e38*e24*e07+e38*e04*e27-e05*e30*e20-e05*e36*e26-e05*e31*e21-e05*e37*e27-e25*e30*e00-e25*e37*e07-e25*e31*e01-e25*e36*e06-e35*e20*e00;
	Rc(7,ord[14])=e12*e21*e34+e18*e25*e38+e12*e30*e23+e12*e20*e33+e12*e32*e25+e12*e22*e35+e12*e31*e24+e22*e30*e13+e22*e10*e33+e22*e32*e15+e22*e31*e14+e22*e11*e34+e32*e20*e13+e32*e10*e23+e32*e21*e14-e25*e30*e10-e25*e36*e16-e25*e31*e11-e25*e37*e17-e35*e20*e10-e35*e26*e16-e35*e21*e11-e35*e27*e17+e15*e33*e23+3*e15*e35*e25+e15*e34*e24+e25*e33*e13+e25*e34*e14+e35*e23*e13+e35*e24*e14+e18*e33*e26+e18*e23*e36+e18*e35*e28+e18*e34*e27+e18*e24*e37+e28*e33*e16+e28*e13*e36+e28*e15*e38+e28*e34*e17+e28*e14*e37+e38*e23*e16+e38*e13*e26+e38*e24*e17+e38*e14*e27-e15*e30*e20-e15*e36*e26-e15*e31*e21-e15*e37*e27+e32*e11*e24;
	Rc(7,ord[15])=-0.5*e35*e202-0.5*e35*e262-0.5*e35*e212-0.5*e35*e272+e25*e34*e24+e28*e23*e36+e28*e25*e38+e28*e34*e27+e28*e24*e37+e38*e23*e26+e38*e24*e27-e25*e30*e20-e25*e36*e26-e25*e31*e21-e25*e37*e27+e25*e33*e23+0.5*e222*e35+1.5*e35*e252+0.5*e35*e232+0.5*e35*e242+0.5*e35*e282+e22*e30*e23+e22*e20*e33+e22*e32*e25+e22*e31*e24+e22*e21*e34+e32*e20*e23+e32*e21*e24+e28*e33*e26;
	Rc(7,ord[16])=-e35*e30*e00-e35*e31*e01-e35*e36*e06-e35*e37*e07+0.5*e322*e05+0.5*e05*e332+0.5*e05*e342+1.5*e05*e352+0.5*e05*e382-0.5*e05*e302-0.5*e05*e362-0.5*e05*e312-0.5*e05*e372+e02*e30*e33+e02*e31*e34+e02*e32*e35+e32*e30*e03+e32*e00*e33+e32*e31*e04+e32*e01*e34+e35*e33*e03+e35*e34*e04+e08*e33*e36+e08*e34*e37+e08*e35*e38+e38*e33*e06+e38*e03*e36+e38*e34*e07+e38*e04*e37;
	Rc(7,ord[17])=-e35*e30*e10+e12*e32*e35-0.5*e15*e362-0.5*e15*e312-0.5*e15*e372-e35*e36*e16+0.5*e322*e15+0.5*e15*e332+0.5*e15*e342+1.5*e15*e352+0.5*e15*e382-0.5*e15*e302+e12*e30*e33+e12*e31*e34+e32*e30*e13+e32*e10*e33+e32*e31*e14+e32*e11*e34+e35*e33*e13+e35*e34*e14+e18*e33*e36+e18*e34*e37+e18*e35*e38+e38*e33*e16+e38*e13*e36+e38*e34*e17+e38*e14*e37-e35*e31*e11-e35*e37*e17;
	Rc(7,ord[18])=-0.5*e25*e302-0.5*e25*e362-0.5*e25*e312-0.5*e25*e372+0.5*e322*e25+0.5*e25*e332+0.5*e25*e342+1.5*e25*e352+0.5*e25*e382+e22*e30*e33+e22*e31*e34+e22*e32*e35+e32*e30*e23+e32*e20*e33+e32*e31*e24+e32*e21*e34+e35*e33*e23+e35*e34*e24+e28*e33*e36+e28*e34*e37+e28*e35*e38+e38*e33*e26+e38*e23*e36+e38*e34*e27+e38*e24*e37-e35*e30*e20-e35*e36*e26-e35*e31*e21-e35*e37*e27;
	Rc(7,ord[19])=e32*e30*e33+e32*e31*e34+0.5*e322*e35+0.5*e35*e332+0.5*e35*e342+0.5*e353+e38*e33*e36+e38*e34*e37+0.5*e35*e382-0.5*e35*e302-0.5*e35*e362-0.5*e35*e312-0.5*e35*e372;
	Rc(8,ord[0])=e02*e00*e06+e02*e01*e07+0.5*e022*e08+e05*e04*e07+e05*e03*e06+0.5*e052*e08+0.5*e08*e062+0.5*e08*e072+0.5*e083-0.5*e08*e042-0.5*e08*e002-0.5*e08*e012-0.5*e08*e032;
	Rc(8,ord[1])=e02*e10*e06+e02*e00*e16+e02*e11*e07+e02*e01*e17+e02*e12*e08+e12*e00*e06+e12*e01*e07+e05*e13*e06+e05*e03*e16+e05*e14*e07+e05*e04*e17+e05*e15*e08+e15*e04*e07+e15*e03*e06+e08*e16*e06+e08*e17*e07-e08*e10*e00-e08*e11*e01-e08*e13*e03-e08*e14*e04+0.5*e022*e18+0.5*e052*e18+1.5*e18*e082+0.5*e18*e062+0.5*e18*e072-0.5*e18*e042-0.5*e18*e002-0.5*e18*e012-0.5*e18*e032;
	Rc(8,ord[2])=e12*e01*e17+0.5*e152*e08+0.5*e08*e162+1.5*e08*e182+0.5*e08*e172-0.5*e08*e102-0.5*e08*e112-0.5*e08*e132-0.5*e08*e142+e05*e13*e16+e05*e14*e17+e05*e15*e18+e15*e13*e06+e15*e03*e16+e15*e14*e07+e15*e04*e17+e18*e16*e06+e18*e17*e07-e18*e10*e00-e18*e11*e01-e18*e13*e03-e18*e14*e04+0.5*e122*e08+e02*e10*e16+e02*e12*e18+e02*e11*e17+e12*e10*e06+e12*e00*e16+e12*e11*e07;
	Rc(8,ord[3])=e12*e10*e16+0.5*e122*e18+e12*e11*e17+e15*e13*e16+e15*e14*e17+0.5*e152*e18+0.5*e18*e162+0.5*e183+0.5*e18*e172-0.5*e18*e102-0.5*e18*e112-0.5*e18*e132-0.5*e18*e142;
	Rc(8,ord[4])=-e08*e20*e00+e08*e27*e07-e08*e21*e01-e08*e23*e03-e08*e24*e04+e02*e20*e06+e02*e00*e26+e02*e21*e07+e02*e01*e27+e02*e22*e08+e22*e00*e06+e22*e01*e07+e05*e23*e06+e05*e03*e26+e05*e24*e07+e05*e04*e27+e05*e25*e08+e25*e04*e07+e25*e03*e06+e08*e26*e06+0.5*e022*e28+0.5*e052*e28+1.5*e28*e082+0.5*e28*e062+0.5*e28*e072-0.5*e28*e042-0.5*e28*e002-0.5*e28*e012-0.5*e28*e032;
	Rc(8,ord[5])=e22*e10*e06+e22*e11*e07+e22*e01*e17+e05*e23*e16+e05*e13*e26+e05*e25*e18+e05*e15*e28+e05*e24*e17+e05*e14*e27+e15*e23*e06+e15*e03*e26+e15*e24*e07+e15*e04*e27+e15*e25*e08+e25*e13*e06+e25*e03*e16+e25*e14*e07+e25*e04*e17+e08*e26*e16+3*e08*e28*e18+e08*e27*e17+e18*e26*e06+e18*e27*e07+e22*e00*e16+e28*e16*e06+e28*e17*e07-e08*e20*e10-e08*e21*e11-e08*e23*e13-e08*e24*e14-e18*e20*e00-e18*e21*e01-e18*e23*e03-e18*e24*e04-e28*e10*e00-e28*e11*e01-e28*e13*e03-e28*e14*e04+e02*e20*e16+e02*e10*e26+e02*e22*e18+e02*e12*e28+e02*e21*e17+e02*e11*e27+e12*e20*e06+e12*e00*e26+e12*e21*e07+e12*e01*e27+e12*e22*e08;
	Rc(8,ord[6])=-e18*e24*e14-e18*e21*e11-e18*e23*e13-e18*e20*e10+e18*e27*e17+e18*e26*e16+e25*e14*e17+e25*e13*e16+e15*e25*e18+e15*e14*e27+e15*e24*e17+e15*e13*e26+e15*e23*e16+e22*e11*e17+e22*e10*e16+e12*e11*e27+e12*e21*e17+e12*e22*e18+e12*e10*e26+e12*e20*e16+0.5*e28*e162+0.5*e28*e172+1.5*e28*e182+0.5*e152*e28-0.5*e28*e142-0.5*e28*e112-0.5*e28*e132-0.5*e28*e102+0.5*e122*e28;
	Rc(8,ord[7])=-e28*e24*e04-e28*e21*e01-e28*e23*e03-e28*e20*e00+e28*e27*e07+e28*e26*e06+e25*e04*e27+e25*e24*e07+e25*e03*e26+e05*e24*e27+e05*e25*e28+e05*e23*e26+e22*e01*e27+e22*e21*e07+e22*e00*e26+e22*e20*e06+e02*e22*e28+e02*e20*e26+e02*e21*e27+0.5*e222*e08-0.5*e08*e242-0.5*e08*e212-0.5*e08*e232-0.5*e08*e202+0.5*e08*e262+0.5*e08*e272+1.5*e08*e282+0.5*e252*e08+e25*e23*e06;
	Rc(8,ord[8])=e25*e24*e17+e25*e14*e27+e28*e26*e16+e28*e27*e17-e28*e21*e11-e28*e24*e14+e12*e22*e28+e22*e10*e26+e22*e21*e17+e22*e11*e27+e15*e23*e26+e15*e25*e28+e15*e24*e27+e25*e23*e16+e25*e13*e26+e22*e20*e16+0.5*e222*e18+0.5*e252*e18+0.5*e18*e262+0.5*e18*e272+e12*e20*e26+e12*e21*e27-e28*e20*e10-0.5*e18*e232-0.5*e18*e242-e28*e23*e13-0.5*e18*e212+1.5*e18*e282-0.5*e18*e202;
	Rc(8,ord[9])=e22*e20*e26+e22*e21*e27+0.5*e222*e28+e25*e23*e26+0.5*e252*e28+e25*e24*e27+0.5*e28*e262+0.5*e28*e272+0.5*e283-0.5*e28*e202-0.5*e28*e212-0.5*e28*e232-0.5*e28*e242;
	Rc(8,ord[10])=-e08*e30*e00-0.5*e38*e042-0.5*e38*e002-0.5*e38*e012-0.5*e38*e032+1.5*e38*e082+0.5*e38*e062+0.5*e38*e072+e32*e01*e07+e05*e33*e06+e05*e03*e36+e05*e34*e07+e05*e04*e37+e05*e35*e08+e35*e04*e07+e35*e03*e06+e08*e36*e06+e08*e37*e07+0.5*e052*e38+e32*e00*e06+e02*e30*e06+e02*e00*e36+e02*e31*e07+e02*e01*e37+e02*e32*e08+0.5*e022*e38-e08*e33*e03-e08*e31*e01-e08*e34*e04;
	Rc(8,ord[11])=-e38*e11*e01-e38*e14*e04-e38*e10*e00-e38*e13*e03-e18*e30*e00-e18*e33*e03-e18*e31*e01-e18*e34*e04-e08*e30*e10-e08*e33*e13-e08*e31*e11-e08*e34*e14+3*e08*e38*e18+e08*e37*e17+e18*e36*e06+e18*e37*e07+e38*e16*e06+e38*e17*e07+e15*e35*e08+e35*e13*e06+e35*e03*e16+e35*e14*e07+e35*e04*e17+e08*e36*e16+e05*e35*e18+e05*e15*e38+e15*e33*e06+e15*e03*e36+e15*e34*e07+e15*e04*e37+e05*e14*e37+e12*e30*e06+e12*e31*e07+e12*e01*e37+e12*e00*e36+e12*e32*e08+e32*e10*e06+e32*e00*e16+e32*e11*e07+e32*e01*e17+e05*e33*e16+e05*e13*e36+e05*e34*e17+e02*e30*e16+e02*e10*e36+e02*e32*e18+e02*e12*e38+e02*e31*e17+e02*e11*e37;
	Rc(8,ord[12])=e12*e30*e16+e12*e10*e36+e12*e32*e18+e12*e31*e17+e12*e11*e37+e32*e10*e16+e32*e11*e17+e15*e33*e16+e15*e13*e36-0.5*e38*e102-0.5*e38*e112-0.5*e38*e132-0.5*e38*e142+0.5*e38*e162+0.5*e38*e172+e15*e34*e17+e15*e14*e37+e15*e35*e18+e35*e13*e16+e35*e14*e17+e18*e36*e16+e18*e37*e17-e18*e30*e10-e18*e33*e13-e18*e31*e11-e18*e34*e14+0.5*e122*e38+0.5*e152*e38+1.5*e38*e182;
	Rc(8,ord[13])=e22*e30*e06-e28*e34*e04+e05*e35*e28+e02*e22*e38+e22*e00*e36+e22*e31*e07+e22*e01*e37+e02*e32*e28+e02*e21*e37-e38*e20*e00-e28*e31*e01-e38*e23*e03-e38*e21*e01-e38*e24*e04-e28*e30*e00-e08*e30*e20-e08*e31*e21-e08*e33*e23-e08*e34*e24-e28*e33*e03+e35*e24*e07+e35*e04*e27+e08*e36*e26+e08*e37*e27+3*e08*e38*e28+e28*e36*e06+e28*e37*e07+e38*e26*e06+e38*e27*e07+e25*e04*e37+e25*e35*e08+e35*e23*e06+e35*e03*e26+e05*e23*e36+e05*e25*e38+e05*e34*e27+e05*e24*e37+e25*e33*e06+e25*e03*e36+e25*e34*e07+e05*e33*e26+e32*e21*e07+e32*e01*e27+e22*e32*e08+e32*e20*e06+e32*e00*e26+e02*e30*e26+e02*e20*e36+e02*e31*e27;
	Rc(8,ord[14])=e35*e13*e26-e38*e21*e11-e38*e24*e14+e35*e24*e17+e35*e14*e27+e18*e36*e26+e18*e37*e27+3*e18*e38*e28+e28*e36*e16+e28*e37*e17+e38*e26*e16+e38*e27*e17-e18*e30*e20-e18*e31*e21-e18*e33*e23-e18*e34*e24-e28*e30*e10-e28*e33*e13-e28*e31*e11-e28*e34*e14-e38*e20*e10-e38*e23*e13+e35*e23*e16+e12*e20*e36+e12*e30*e26+e12*e31*e27+e12*e21*e37+e12*e32*e28+e12*e22*e38+e22*e30*e16+e22*e10*e36+e22*e32*e18+e22*e31*e17+e22*e11*e37+e32*e20*e16+e32*e10*e26+e32*e21*e17+e32*e11*e27+e15*e33*e26+e15*e23*e36+e15*e35*e28+e15*e25*e38+e15*e34*e27+e15*e24*e37+e25*e33*e16+e25*e13*e36+e25*e34*e17+e25*e14*e37+e25*e35*e18;
	Rc(8,ord[15])=-e28*e30*e20+e22*e30*e26+e22*e20*e36+e22*e31*e27+e22*e21*e37+e22*e32*e28+e32*e20*e26+e32*e21*e27+e25*e33*e26+e25*e23*e36+e25*e35*e28+e25*e34*e27+e25*e24*e37+e35*e23*e26+e35*e24*e27+e28*e36*e26+e28*e37*e27-e28*e31*e21-e28*e33*e23-e28*e34*e24-0.5*e38*e242+0.5*e252*e38+1.5*e38*e282+0.5*e38*e262+0.5*e38*e272-0.5*e38*e202-0.5*e38*e212-0.5*e38*e232+0.5*e222*e38;
	Rc(8,ord[16])=-0.5*e08*e312-0.5*e08*e342+0.5*e352*e08+0.5*e08*e362+1.5*e08*e382+0.5*e08*e372-0.5*e08*e302-0.5*e08*e332+e02*e30*e36+e02*e32*e38+e02*e31*e37+e32*e30*e06+e32*e00*e36+e32*e31*e07+e32*e01*e37+e05*e33*e36+e05*e34*e37+e05*e35*e38+e35*e33*e06+e35*e03*e36+e35*e34*e07+e35*e04*e37+0.5*e322*e08+e38*e36*e06+e38*e37*e07-e38*e30*e00-e38*e33*e03-e38*e31*e01-e38*e34*e04;
	Rc(8,ord[17])=-e38*e30*e10+e38*e36*e16+e38*e37*e17-e38*e33*e13-e38*e31*e11-e38*e34*e14+0.5*e18*e362+e12*e30*e36+e12*e32*e38+e12*e31*e37+e32*e30*e16+e32*e10*e36+e32*e31*e17+e32*e11*e37+e15*e33*e36+e15*e34*e37+e15*e35*e38+e35*e33*e16+e35*e13*e36+e35*e34*e17+e35*e14*e37+0.5*e322*e18+0.5*e352*e18+1.5*e18*e382+0.5*e18*e372-0.5*e18*e302-0.5*e18*e332-0.5*e18*e312-0.5*e18*e342;
	Rc(8,ord[18])=-e38*e30*e20+e25*e35*e38+e22*e30*e36+e22*e32*e38+e22*e31*e37+e32*e30*e26+e32*e20*e36+e32*e31*e27+e32*e21*e37+e25*e33*e36+e25*e34*e37+e35*e33*e26+e35*e23*e36+e35*e34*e27+e35*e24*e37+e38*e36*e26+e38*e37*e27-e38*e31*e21-e38*e33*e23-e38*e34*e24-0.5*e28*e332-0.5*e28*e312-0.5*e28*e342+0.5*e322*e28+0.5*e352*e28+0.5*e28*e362+1.5*e28*e382+0.5*e28*e372-0.5*e28*e302;
	Rc(8,ord[19])=e32*e30*e36+0.5*e322*e38+e32*e31*e37+e35*e33*e36+e35*e34*e37+0.5*e352*e38+0.5*e38*e362+0.5*e383+0.5*e38*e372-0.5*e38*e302-0.5*e38*e332-0.5*e38*e312-0.5*e38*e342;
	Rc(9,ord[0])=e00*e04*e08-e00*e05*e07+e03*e02*e07-e03*e01*e08-e06*e02*e04+e06*e01*e05;
	Rc(9,ord[1])=e06*e01*e15-e16*e02*e04+e16*e01*e05+e03*e02*e17-e13*e01*e08+e06*e11*e05+e13*e02*e07+e00*e04*e18+e00*e14*e08-e00*e05*e17-e10*e05*e07-e00*e15*e07-e06*e12*e04-e06*e02*e14-e03*e01*e18-e03*e11*e08+e10*e04*e08+e03*e12*e07;
	Rc(9,ord[2])=-e13*e01*e18-e13*e11*e08+e13*e12*e07+e13*e02*e17+e03*e12*e17-e10*e15*e07+e10*e04*e18+e10*e14*e08-e10*e05*e17-e00*e15*e17+e00*e14*e18+e16*e01*e15+e06*e11*e15-e06*e12*e14-e16*e12*e04-e16*e02*e14+e16*e11*e05-e03*e11*e18;
	Rc(9,ord[3])=e10*e14*e18-e10*e15*e17-e13*e11*e18+e13*e12*e17+e16*e11*e15-e16*e12*e14;
	Rc(9,ord[4])=-e20*e05*e07+e03*e22*e07+e06*e21*e05+e06*e01*e25-e23*e01*e08+e23*e02*e07+e00*e24*e08-e00*e25*e07-e00*e05*e27+e00*e04*e28-e06*e22*e04-e06*e02*e24-e03*e21*e08-e03*e01*e28-e26*e02*e04+e26*e01*e05+e03*e02*e27+e20*e04*e08;
	Rc(9,ord[5])=e23*e12*e07-e26*e02*e14+e16*e21*e05-e23*e11*e08+e10*e24*e08-e20*e05*e17+e26*e11*e05+e26*e01*e15+e10*e04*e28+e00*e24*e18-e00*e15*e27+e03*e22*e17-e13*e01*e28+e23*e02*e17+e16*e01*e25+e20*e04*e18+e06*e11*e25+e13*e02*e27-e23*e01*e18-e20*e15*e07-e10*e25*e07+e13*e22*e07-e06*e22*e14-e26*e12*e04-e03*e11*e28-e03*e21*e18-e16*e22*e04-e16*e02*e24-e06*e12*e24+e06*e21*e15+e00*e14*e28-e00*e25*e17+e20*e14*e08-e13*e21*e08-e10*e05*e27+e03*e12*e27;
	Rc(9,ord[6])=-e13*e11*e28+e13*e12*e27+e13*e22*e17+e16*e11*e25+e10*e14*e28-e13*e21*e18-e23*e11*e18+e23*e12*e17+e20*e14*e18-e20*e15*e17+e26*e11*e15-e10*e15*e27-e10*e25*e17-e16*e22*e14-e16*e12*e24+e16*e21*e15-e26*e12*e14+e10*e24*e18;
	Rc(9,ord[7])=e26*e21*e05+e26*e01*e25+e20*e04*e28+e20*e24*e08-e20*e25*e07+e23*e22*e07+e03*e22*e27-e03*e21*e28-e26*e22*e04-e20*e05*e27-e00*e25*e27+e06*e21*e25-e06*e22*e24+e00*e24*e28-e26*e02*e24-e23*e21*e08-e23*e01*e28+e23*e02*e27;
	Rc(9,ord[8])=-e10*e25*e27+e10*e24*e28-e20*e15*e27-e20*e25*e17+e20*e14*e28+e20*e24*e18+e26*e11*e25+e23*e22*e17-e23*e11*e28+e23*e12*e27-e23*e21*e18-e13*e21*e28+e13*e22*e27-e26*e12*e24+e26*e21*e15-e16*e22*e24+e16*e21*e25-e26*e22*e14;
	Rc(9,ord[9])=-e20*e25*e27+e20*e24*e28-e23*e21*e28+e23*e22*e27-e26*e22*e24+e26*e21*e25;
	Rc(9,ord[10])=e03*e02*e37-e03*e31*e08-e03*e01*e38+e03*e32*e07-e00*e35*e07+e30*e04*e08+e06*e31*e05-e36*e02*e04+e36*e01*e05-e06*e32*e04-e06*e02*e34+e06*e01*e35+e00*e04*e38-e00*e05*e37+e33*e02*e07-e33*e01*e08-e30*e05*e07+e00*e34*e08;
	Rc(9,ord[11])=-e36*e12*e04+e30*e04*e18-e30*e15*e07-e36*e02*e14-e30*e05*e17+e30*e14*e08-e00*e35*e17-e00*e15*e37+e33*e02*e17-e06*e32*e14-e06*e12*e34-e16*e32*e04+e06*e31*e15+e06*e11*e35+e00*e34*e18-e10*e35*e07-e33*e11*e08-e33*e01*e18+e16*e01*e35-e16*e02*e34+e16*e31*e05-e03*e31*e18-e03*e11*e38+e03*e32*e17+e13*e02*e37-e13*e31*e08-e13*e01*e38+e10*e34*e08+e00*e14*e38+e36*e11*e05+e36*e01*e15+e03*e12*e37-e10*e05*e37+e10*e04*e38+e33*e12*e07+e13*e32*e07;
	Rc(9,ord[12])=-e36*e12*e14-e30*e15*e17+e13*e32*e17-e13*e31*e18-e33*e11*e18+e33*e12*e17+e10*e14*e38+e30*e14*e18-e13*e11*e38+e13*e12*e37-e10*e35*e17+e10*e34*e18-e16*e12*e34-e16*e32*e14+e16*e11*e35+e16*e31*e15+e36*e11*e15-e10*e15*e37;
	Rc(9,ord[13])=-e06*e22*e34-e06*e32*e24-e00*e25*e37-e00*e35*e27+e23*e02*e37+e00*e24*e38-e23*e01*e38-e03*e31*e28-e33*e01*e28+e03*e22*e37+e03*e32*e27+e33*e02*e27-e03*e21*e38-e26*e32*e04-e33*e21*e08+e36*e01*e25+e36*e21*e05-e20*e05*e37+e20*e04*e38+e30*e04*e28-e20*e35*e07+e33*e22*e07+e30*e24*e08-e30*e25*e07-e23*e31*e08+e23*e32*e07+e00*e34*e28+e06*e21*e35+e06*e31*e25-e36*e02*e24+e26*e01*e35-e36*e22*e04+e26*e31*e05-e26*e02*e34+e20*e34*e08-e30*e05*e27;
	Rc(9,ord[14])=e33*e22*e17+e33*e12*e27+e16*e21*e35-e16*e22*e34-e16*e32*e24+e23*e32*e17-e23*e11*e38-e23*e31*e18+e23*e12*e37-e13*e21*e38-e13*e31*e28+e13*e22*e37+e36*e21*e15-e36*e12*e24+e36*e11*e25-e26*e12*e34-e20*e35*e17+e20*e14*e38+e20*e34*e18+e30*e24*e18-e30*e15*e27-e30*e25*e17+e30*e14*e28-e33*e21*e18+e10*e34*e28+e10*e24*e38-e10*e35*e27-e10*e25*e37-e20*e15*e37-e26*e32*e14+e26*e11*e35+e26*e31*e15-e36*e22*e14+e13*e32*e27+e16*e31*e25-e33*e11*e28;
	Rc(9,ord[15])=-e20*e35*e27-e20*e25*e37+e20*e34*e28+e20*e24*e38+e30*e24*e28-e30*e25*e27+e23*e32*e27+e23*e22*e37-e23*e31*e28-e23*e21*e38+e33*e22*e27-e26*e22*e34-e26*e32*e24+e26*e21*e35+e26*e31*e25-e36*e22*e24+e36*e21*e25-e33*e21*e28;
	Rc(9,ord[16])=-e33*e01*e38-e03*e31*e38+e00*e34*e38+e33*e32*e07+e03*e32*e37+e06*e31*e35-e00*e35*e37-e36*e32*e04-e06*e32*e34-e36*e02*e34+e36*e01*e35+e36*e31*e05+e30*e04*e38+e30*e34*e08-e33*e31*e08+e33*e02*e37-e30*e05*e37-e30*e35*e07;
	Rc(9,ord[17])=-e33*e31*e18-e33*e11*e38+e10*e34*e38+e30*e14*e38-e10*e35*e37-e30*e15*e37-e13*e31*e38+e13*e32*e37-e30*e35*e17+e33*e12*e37+e30*e34*e18+e33*e32*e17+e16*e31*e35-e16*e32*e34-e36*e12*e34-e36*e32*e14+e36*e11*e35+e36*e31*e15;
	Rc(9,ord[18])=-e20*e35*e37+e20*e34*e38+e30*e24*e38-e30*e35*e27-e30*e25*e37+e30*e34*e28+e23*e32*e37-e23*e31*e38-e33*e21*e38-e33*e31*e28+e33*e22*e37+e33*e32*e27+e26*e31*e35-e26*e32*e34-e36*e22*e34-e36*e32*e24+e36*e21*e35+e36*e31*e25;
	Rc(9,ord[19])=-e33*e31*e38-e30*e35*e37+e36*e31*e35+e33*e32*e37+e30*e34*e38-e36*e32*e34;
	return true;
}

bool L_Algoritmo5Puntos::calcMatricesBC(L_PolinGr11_matriz &B, L_PolinGr11_matriz &C, L_Matrix &Rc)
{
	B.reallocate(4,4);
	C.reallocate(4,4);
	L_Polin_xyz c, d, e, f, g, h, i;
	L_Polin_xyz j, k, l, m;
	L_Polin_xyz zg, zh, xh, ze, yg, zf;
	double L, M, N, O, P, Q, R, S;
	int u, v;
	c.sacaFila(Rc,2);
	d.sacaFila(Rc,3);
	e.sacaFila(Rc,4);
	f.sacaFila(Rc,5);
	g.sacaFila(Rc,6);
	h.sacaFila(Rc,7);
	i.sacaFila(Rc,8);
	L=Rc(6,9);
	M=Rc(6,13);
	N=Rc(6,14);
	O=Rc(6,15);
	P=Rc(7,9);
	Q=Rc(7,10);
	R=Rc(7,11);
	S=Rc(7,12);
	zg.multZ(g);
	zh.multZ(h);
	xh.multX(h);
	ze.multZ(e); // Error: decia multX
	yg.multY(g);
	zf.multZ(f);
	j = e - zg; // ERROR: decia i =
	k = f - zh;
	l = d - xh + P*c + Q*ze + R*e + S*g;
	m = c - yg + L*d + M*zf + N*f + O*h;
	// Ahora tenemos las expresiones para i, j, k, l, m
	// Ahora hay que factorizarlos como xy[) + x[) + y[) + [)

	i.factoriza_xy_x_y_1(B(0,0), B(0,1), B(0,2), B(0,3), 1, 2, 2, 3);
	j.factoriza_xy_x_y_1(B(1,0), B(1,1), B(1,2), B(1,3), 1, 3, 3, 4);
	k.factoriza_xy_x_y_1(B(2,0), B(2,1), B(2,2), B(2,3), 1, 3, 3, 4);
	l.factoriza_xy_x_y_1(B(3,0), B(3,1), B(3,2), B(3,3), 2, 3, 3, 4);

	for (u=0; u<3; u++)
		for (v=0; v<4; v++)
			C(u,v)=B(u,v);
	m.factoriza_xy_x_y_1(C(3,0), C(3,1), C(3,2), C(3,3), 2, 3, 3, 4);

	return true;
}

bool L_SegmentadorColoresLUTdif::pintaImagen(L_ImageRGBUchar &output)
{
	int i, j;
	L_uchar c;
	output.reallocate(imSegmentada.lx, imSegmentada.ly);
	for (j=0; j<output.ly; j++)
	{
		for (i=0; i<output.lx; i++)
		{
			c=imSegmentada.pix(i,j);
			output.pix(i,j,0)=(L_uchar)colorIndice[c].colorCentralR;
			output.pix(i,j,1)=(L_uchar)colorIndice[c].colorCentralG;
			output.pix(i,j,2)=(L_uchar)colorIndice[c].colorCentralB;
		}
	}
	return true;
}

bool L_SegmentadorColoresLUTdif::segmentaDif(const L_ImageRGBUchar &original)
{
	int i, j;
	int u, v;
	L_uchar R, G, B;
	bool enBlob = false;
	if (original.lx < 1 || original.ly < 1 || original.data()==NULL)
		return false;
	imSegmentada.reallocate(original.lx, original.ly);
	imSegmentada.setZero();
	// Informacion contenida en cada recursion.pix[i][j]
	// 0 = sin procesar, 1 = avance derecha, 2 = avance arriba, 4 = izq, 8 = abajo
	recursion.reallocate(original.lx, original.ly);
	recursion.setZero();
	u=0;
	v=0;
	blobs.resize(0);

	//int uStop=100; // Para debug
	//int vStop=100;

	for (j=0; j<original.ly; j++)
	{
		for (i=0; i<original.lx; i++)
		{
			if (enBlob == false)
			{
				switch(recursion.pix(i,j))
				{
				case 0:
					R=original.pix(i,j,0);
					G=original.pix(i,j,1);
					B=original.pix(i,j,2);
					blobs.resize(blobs.size()+1);
					enBlob=true;
					blobs[blobs.size()-1].resetea();
					blobs[blobs.size()-1].indiceColor=LUT[R/dR][G/dG][B/dB];
					u=i;
					v=j;
					break;					
				default:
					break;
				}
			}
			else
			{
				while(enBlob)
				{
					//if (u==uStop && v==vStop)
					//	printf("Aca\n");
					R=original.pix(u,v,0);
					G=original.pix(u,v,1);
					B=original.pix(u,v,2);
					if (LUT[R/dR][G/dG][B/dB] != blobs[blobs.size()-1].indiceColor)
					{
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							printf("Error interno de difusion en L_SegmentadorColoresLUTdif::segmentaDif()\n");
							break;
						case 16:
							recursion.pix(u,v)=0;
							u--;
							break;
						case 32:
							recursion.pix(u,v)=0;
							v++;
							break;
						case 64:
							recursion.pix(u,v)=0;
							u++;
							break;
						case 128:
							recursion.pix(u,v)=0;
							v--;
							break;
						}
					}
					switch(recursion.pix(u,v)&0x0F)
					{
					case 0:
						imSegmentada.pix(u,v)=(L_uchar)(blobs[blobs.size()-1].indiceColor);
						blobs[blobs.size()-1].agregaPixel(R,G,B,u,v);
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 1;
						if (u<recursion.lx-1 && (recursion.pix(u+1,v)&0x0F)==0)
						{
							u++;
							recursion.pix(u,v) = 16;
						}
						break;
					case 1:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 2;
						if (v>0 && (recursion.pix(u,v-1)&0x0F)==0)
						{
							v--;
							recursion.pix(u,v) = 32;
						}
						break;
					case 2:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 4;
						if (u>0 && (recursion.pix(u-1,v)&0x0F)==0)
						{
							u--;
							recursion.pix(u,v) = 64;
						}
						break;
					case 4:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 8;
						if (v<recursion.ly-1 && (recursion.pix(u,v+1)&0x0F)==0)
						{
							v++;
							recursion.pix(u,v) = 128;
						}
						break;
					case 8:
						// Aca supuestamente se vuelve a u=i y v=j
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							enBlob=false;
							if (blobs[blobs.size()-1].area < areaMin)
								blobs.resize(blobs.size()-1); // que pena, pero se clean :)
							break;
						case 16:
							u--;
							break;
						case 32:
							v++;
							break;
						case 64:
							u++;
							break;
						case 128:
							v--;
							break;
						}
						break;
					}
				}
			}
		}
	}
	return true;
}

bool L_SegmentadorColoresLUTdif::segmentaDif(const L_ImageGrayUchar &original)
{
	int i, j;
	int u, v;
	bool enBlob = false;
	if (original.lx < 1 || original.ly < 1 || original.data()==NULL)
		return false;
	imSegmentada.reallocate(original.lx, original.ly);
	imSegmentada.setZero();
	// Informacion contenida en cada recursion.pix[i][j]
	// 0 = sin procesar, 1 = avance derecha, 2 = avance arriba, 4 = izq, 8 = abajo
	recursion.reallocate(original.lx, original.ly);
	recursion.setZero();
	u=0;
	v=0;
	blobs.resize(0);

	//int uStop=100; // Para debug
	//int vStop=100;

	for (j=0; j<original.ly; j++)
	{
		for (i=0; i<original.lx; i++)
		{
			if (enBlob == false && original.pix(i,j) != 0)
			{
				switch(recursion.pix(i,j))
				{
				case 0:
					blobs.resize(blobs.size()+1);
					enBlob=true;
					blobs[blobs.size()-1].resetea();
					blobs[blobs.size()-1].indiceColor=original.pix(i,j);
					u=i;
					v=j;
					break;					
				default:
					break;
				}
			}
			else
			{
				while(enBlob)
				{
					//if (u==uStop && v==vStop)
					//	printf("Aca\n");
					if (original.pix(u,v) != blobs[blobs.size()-1].indiceColor)
					{
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							printf("Error interno de difusion en L_SegmentadorColoresLUTdif::segmentaDif()\n");
							break;
						case 16:
							recursion.pix(u,v)=0;
							u--;
							break;
						case 32:
							recursion.pix(u,v)=0;
							v++;
							break;
						case 64:
							recursion.pix(u,v)=0;
							u++;
							break;
						case 128:
							recursion.pix(u,v)=0;
							v--;
							break;
						}
					}
					blobs[blobs.size()-1].agregaKeypoint(u,v);
					switch(recursion.pix(u,v)&0x0F)
					{
					case 0:
						imSegmentada.pix(u,v)=(L_uchar)(blobs[blobs.size()-1].indiceColor);
						blobs[blobs.size()-1].agregaPixel(u,v);
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 1;
						if (u<recursion.lx-1 && (recursion.pix(u+1,v)&0x0F)==0)
						{
							u++;
							recursion.pix(u,v) = 16;
						}
						break;
					case 1:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 2;
						if (v>0 && (recursion.pix(u,v-1)&0x0F)==0)
						{
							v--;
							recursion.pix(u,v) = 32;
						}
						break;
					case 2:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 4;
						if (u>0 && (recursion.pix(u-1,v)&0x0F)==0)
						{
							u--;
							recursion.pix(u,v) = 64;
						}
						break;
					case 4:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 8;
						if (v<recursion.ly-1 && (recursion.pix(u,v+1)&0x0F)==0)
						{
							v++;
							recursion.pix(u,v) = 128;
						}
						break;
					case 8:
						// Aca supuestamente se vuelve a u=i y v=j
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							enBlob=false;
							if (blobs[blobs.size()-1].area < areaMin)
								blobs.resize(blobs.size()-1); // que pena, pero se clean :)
							break;
						case 16:
							u--;
							break;
						case 32:
							v++;
							break;
						case 64:
							u++;
							break;
						case 128:
							v--;
							break;
						}
						break;
					}
				}
			}
		}
	}
	return true;
}

bool L_Segmentador2ColoresAdaptivo::segmentaDif(const L_ImageRGBUchar &original)
{
	int i, j;
	int u, v;
	L_uchar R, G, B;
	bool enBlob = false;
	int indice;
	if (original.lx < 1 || original.ly < 1 || original.data()==NULL)
		return false;
	imSegmentada.reallocate(original.lx, original.ly);
	recursion.reallocate(original.lx, original.ly);
	// Informacion contenida en cada recursion.pix[i][j]
	// 0 = sin procesar, 1 = avance derecha, 2 = avance arriba, 4 = izq, 8 = abajo
	for (j=0; j<imSegmentada.ly; j++)
		for (i=0; i<imSegmentada.lx; i++)
			recursion.pix(i,j)=0;
	u=0;
	v=0;
	blobs.resize(0);

	//int uStop=100; // Para debug
	//int vStop=100;

	for (j=0; j<original.ly; j++)
	{
		for (i=0; i<original.lx; i++)
		{
			if (enBlob == false)
			{
				switch(recursion.pix(i,j))
				{
				case 0:
					R=original.pix(i,j,0);
					G=original.pix(i,j,1);
					B=original.pix(i,j,2);
					blobs.resize(blobs.size()+1);
					enBlob=true;
					blobs[blobs.size()-1].resetea();
					blobs[blobs.size()-1].indiceColor=(facR*R+facG*G+facB*B+fac1 > 0) ? 1 : 0;
					u=i;
					v=j;
					break;					
				default:
					break;
				}
			}
			else
			{
				while(enBlob)
				{
					//if (u==uStop && v==vStop)
					//	printf("Aca\n");
					R=original.pix(u,v,0);
					G=original.pix(u,v,1);
					B=original.pix(u,v,2);
					indice = (facR*R+facG*G+facB*B+fac1 > 0) ? 1 : 0;
					if (indice != blobs[blobs.size()-1].indiceColor)
					{
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							printf("Error interno de difusion en L_SegmentadorColoresLUTdif::segmentaDif()\n");
							break;
						case 16:
							recursion.pix(u,v)=0;
							u--;
							break;
						case 32:
							recursion.pix(u,v)=0;
							v++;
							break;
						case 64:
							recursion.pix(u,v)=0;
							u++;
							break;
						case 128:
							recursion.pix(u,v)=0;
							v--;
							break;
						}
					}
					switch(recursion.pix(u,v)&0x0F)
					{
					case 0:
						imSegmentada.pix(u,v)=(L_uchar)(blobs[blobs.size()-1].indiceColor);
						blobs[blobs.size()-1].agregaPixel(R,G,B,u,v);
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 1;
						if (u<recursion.lx-1 && (recursion.pix(u+1,v)&0x0F)==0)
						{
							u++;
							recursion.pix(u,v) = 16;
						}
						break;
					case 1:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 2;
						if (v>0 && (recursion.pix(u,v-1)&0x0F)==0)
						{
							v--;
							recursion.pix(u,v) = 32;
						}
						break;
					case 2:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 4;
						if (u>0 && (recursion.pix(u-1,v)&0x0F)==0)
						{
							u--;
							recursion.pix(u,v) = 64;
						}
						break;
					case 4:
						recursion.pix(u,v) = (recursion.pix(u,v)&0xF0) | 8;
						if (v<recursion.ly-1 && (recursion.pix(u,v+1)&0x0F)==0)
						{
							v++;
							recursion.pix(u,v) = 128;
						}
						break;
					case 8:
						// Aca supuestamente se vuelve a u=i y v=j
						switch(recursion.pix(u,v)&0xF0)
						{
						case 0:
							enBlob=false;
							if (blobs[blobs.size()-1].area < areaMin)
								blobs.resize(blobs.size()-1); // que pena, pero se clean :)
							break;
						case 16:
							u--;
							break;
						case 32:
							v++;
							break;
						case 64:
							u++;
							break;
						case 128:
							v--;
							break;
						}
						break;
					}
				}
			}
		}
	}
	return true;
}



bool L_CaracterizadorBlobs::ajustaCirculoAlgebraico(double *x, double *y, bool *acept, int n, double &r, double &cx, double &cy)
{
	// De los nubots...
	double sumXY,sumXX,sumYY,sumXZ,sumYZ;
	double meanXY,meanXX,meanYY,meanXZ,meanYZ;
	double B,C,G11,G12,G22,D1,D2;
	double meanX,meanY;
	double meanD = 0;
	int i, nPoints;
	double pound;
	double sumPounds = 0;

	meanX = meanY = 0;
	nPoints = n;

	for (i=0; i < nPoints; i++)
	{
		pound = 1;
		if (acept != NULL && acept[i] == false)
			continue; // pound zero
		sumPounds += pound;
		meanX += x[i]*pound;
		meanY += y[i]*pound;
	}

	meanX = (int) meanX/sumPounds;
	meanY = (int) meanY/sumPounds;

	double Xi,Yi,Zi;
	sumXY = sumXX = sumYY = sumXZ = sumYZ = 0;
	for (i=0; i < nPoints; i++)
	{
		pound = 1;
		if (acept != NULL && acept[i] == false)
			continue; // pound zero

		Xi = x[i] - meanX;
		Yi = y[i] - meanY;
		Zi = Xi*Xi + Yi*Yi;
		sumXY = sumXY + Xi*Yi*pound;
		sumXX = sumXX + Xi*Xi*pound;
		sumYY = sumYY + Yi*Yi*pound;
		sumXZ = sumXZ + Xi*Zi*pound;
		sumYZ = sumYZ + Yi*Zi*pound;
	}

	meanXX = sumXX/sumPounds;
	meanYY = sumYY/sumPounds;
	meanXY = sumXY/sumPounds;
	meanXZ = sumXZ/sumPounds;
	meanYZ = sumYZ/sumPounds;

	G11 = sqrt(meanXX);
	if (G11 == 0)
		return false;

	G12 = meanXY/G11;

	if ((meanYY - G12*G12) < 0)
		return false;

	G22 = sqrt(meanYY - G12*G12);

	if (G22==0)
		return false;

	D1 = meanXZ/G11;
	D2 = (meanYZ - D1*G12)/G22;
	C = D2/G22;
	B = (D1 - G12*C)/G11;

	cx = B/2;
	cy = C/2;
	r = sqrt(cx*cx + cy*cy + meanXX + meanYY);
	cx += meanX;
	cy += meanY;
	return true;
}

bool L_CaracterizadorBlobs::ajustaCirculoRansac(double *x, double *y, int n, int nIntentos, bool *acept, int nMin, double maxErr, double &r, double &cx, double &cy)
{
	int in, k, cand, temp, nacept;
	std::vector<int> nums;
	double t1, t3, t5, t7, t9, t15, t19, t27, t41;
	double x0, y0, x1, y1, x2, y2, err2, maxErr2;
	nums.resize(n);
	maxErr2=maxErr*maxErr;

	if (n < nMin)
		return false;

	for (k=0; k<n; k++)
		nums[k]=k;

	for (in=0; in < nIntentos; in++)
	{
		// Elegir 3 puntos distintos
		cand = rand()%n;
		temp = nums[cand];
		nums[cand] = nums[0];
		nums[0] = temp;
		//
		cand = rand()%(n-1)+1;
		temp = nums[cand];
		nums[cand] = nums[1];
		nums[1] = temp;
		//
		cand = rand()%(n-2)+2;
		temp = nums[cand];
		nums[cand] = nums[2];
		nums[2] = temp;
		//
		x0 = x[nums[0]];
		y0 = y[nums[0]];
		x1 = x[nums[1]];
		y1 = y[nums[1]];
		x2 = x[nums[2]];
		y2 = y[nums[2]];
		// Calcular circulo usando los 3 puntos
		t1 = x1 * x1;
		t3 = y0 * y0;
		t5 = y1 * y1;
		t7 = x0 * x0;
		t9 = y2 * y2;
		t15 = x2 * x2;
		t19 = -y2 * t1 - y1 * t3 - y2 * t5 + y2 * t7 - y0 * t9 + y2 * t3 + y0 * t1 + t9 * y1 + y0 * t5 + t15 * y1 - y0 * t15 - y1 * t7;
		t27 = 0.1e1 / (-x1 * y2 + y2 * x0 + y0 * x1 + x2 * y1 - y1 * x0 - y0 * x2);
		cx = t19 * t27 / 0.2e1;
		t41 = x1 * t9 - t9 * x0 + x0 * t5 + x2 * t7 - x2 * t1 - x1 * t7 + x2 * t3 + x1 * t15 + x0 * t1 - x2 * t5 - t3 * x1 - t15 * x0;
		cy = -t41 * t27 / 0.2e1;
		r = sqrt( (x1-cx)*(x1-cx) + (y1-cy)*(y1-cy));
		// test de consenso
		nacept = 0;
		acept[nums[0]] = true;
		acept[nums[1]] = true;
		acept[nums[2]] = true;
		for (k=3; k<n; k++)
			acept[nums[k]] = false;
		//extern L_ShapeArray wa;
		for (k=0; k<n; k++)
		{
			err2 = (x[nums[k]] - cx)*(x[nums[k]] - cx) + (y[nums[k]] - cy)*(y[nums[k]] - cy) - r*r;
			if ( err2 > -maxErr2 && err2 < maxErr2 )
			{
				nacept++;
				acept[nums[k]] = true;
			}
		}
		if (nacept >= nMin)
		{
			/*
			for (k=0; k<n; k++)
			{
				if (acept[k])
					wa.drawLine((int)cx, (int)cy, (int)x[k], (int)y[k], 255, 0, 0);
				else
					wa.drawLine((int)cx, (int)cy, (int)x[k], (int)y[k], 0, 255, 255);
			}
			*/
			return ajustaCirculoAlgebraico(x, y, acept, n, r, cx, cy);
		}
		/*
		else
			for (k=0; k<n; k++)
				wa.drawLine((int)cx, (int)cy, (int)x[k], (int)y[k], 0, 0, 0);
				*/
	}
	return false;
}



bool L_CaracterizadorBlobs::scanLineFromImage(L_ImageGrayUchar &imSegm, L_uchar c, double x, double y, double dx, double dy, int toler, int &xf, int &yf)
{
	int xi, yi;
	int nFuera = 0;
	xi = (int)x;
	yi = (int)y;

	while (xi>=0 && yi>=0 && xi<imSegm.lx && yi<imSegm.ly)
	{
		if (imSegm.pix(xi,yi) == c)
			nFuera = 0;
		else
			nFuera++;
		if (nFuera > toler)
		{
			xf = (int)(x - dx*nFuera);
			yf = (int)(y - dy*nFuera);
			return true;
		}
		x=x+dx;
		y=y+dy;
		xi = (int)x;
		yi = (int)y;
	}
	return false;
}

int L_RelGeom2D::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_RelGeom2D::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err;
	int i;
	if (c < 0)
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));
			corr += (err <= errMax);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));
			if (err <= errMax)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}

void L_RelEpipolarFundamental::calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni, int nRev)
{
	int nMax = x1y1x2y2.li;
	if (nRev >= 0 && iIni + nRev < nMax)
		nMax = iIni + nRev;
	err.resize(x1y1x2y2.li);
	for (int i=0; i<nMax-iIni; i++)
		err[i] = errorRelGeom(x1y1x2y2(i+iIni,0), x1y1x2y2(i+iIni,1), x1y1x2y2(i+iIni,2), x1y1x2y2(i+iIni,3));
}

int L_RelEpipolarFundamental::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_RelEpipolarFundamental::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err;
	int i;

	if (c < 0) // Revisarlos todos rapido
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));
			corr += (err <= errMax);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));

			if (err <= errMax)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}


int L_RelEpipolarFundamental::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_MatrizFundamental Fs[3];
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_RelEpipolarFundamental *ptr;
	int k, ihip=0;
	int rango;
	int nSols;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip;
		if (rango >= x1y1x2y2.li)
			rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (algoritmo7puntos(Fs, nSols, xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_RelEpipolarFundamental();
			ptr->F = Fs[k];
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_RelEpipolarFundamental::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_MatrizFundamental Fs[3];
	std::vector<int> perm(gradosLibertad());
	L_Matrix xyxy;
	L_RelEpipolarFundamental *ptr;
	int k, ihip=0;
	int rango;
	int nSols;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = x1y1x2y2.li-1;
		L_randomSelection(perm, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (algoritmo7puntos(Fs, nSols, xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_RelEpipolarFundamental();
			ptr->F = Fs[k];
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

bool L_RelEpipolarFundamental::calcHartley(L_Matrix &F, const L_Matrix &x1y1x2y2, bool escalamiento)
{
	// "In Defense of the Eight-Point Algorithm", Richard Hartley
	L_Matrix a, aT, aTa, U, d, V, VT;
	L_Matrix tmp, uv, T1, T2;
	int i, j;

	// 1: Normalizar
	uv = x1y1x2y2;
	if (escalamiento)
	{
		if (escalamientoIsotherpico(uv, T1, T2) == false) // Puntos repetidos
			return false;
	}

	// 2: Minimos cuadrados. Se busca que (uR vR) F (uL vL)T = 0
	//      [u2 u3 1]*[e00 e01 e02;e10 e11 e12;e20 e21 e22]*[u0;u1;1] = 0
	//      Maple  ->  u0*u2*e00 + u1*u2*e01 + u2*e02 + u0*u3*e10 + u1*u3*e11 + u3*e12 + u0*e20 + u1*e21 + e22 = 0

	a.reallocate(uv.li, 9);
	for (i=0; i<uv.li; i++)
	{
		a(i,0)=uv(i,0)*uv(i,2); // E00
		a(i,1)=uv(i,1)*uv(i,2); // E01
		a(i,2)=uv(i,2); // E02
		a(i,3)=uv(i,0)*uv(i,3); //E10
		a(i,4)=uv(i,1)*uv(i,3); // E11
		a(i,5)=uv(i,3); // E12
		a(i,6)=uv(i,0); // E20
		a(i,7)=uv(i,1); // E21
		a(i,8)=1; // E22
	}
	// svd_ordenado: F = U d V,  si F[n x a) =>  U[n x a) , d[n x n) , V[n x n) 
	aT.transpOf(a);
	aTa.OP_mult(aT,a);
	if (aTa.svd_ordenado(d, V) == false)
		return false;

	// 3: Menor value propio es V[0..8,8), redimensionarlo a U[012,012)
	U.reallocate(3,3);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			U(i,j)=V(i*3+j,8); // v(:,8) es el vector propio asociado al menor value propio

	// 4: Reparacion matriz esencial para que tenga rango 2
	U.svd_ordenado(d, V); // De mayor a menor value propio
	VT.transpOf(V);
	d(2,0)=0.0;
	tmp.OP_multByDiagonal(U,d);
	F.OP_mult(tmp,VT);

	// Normalizacion inversa
	if (escalamiento)
		escalamientoIsotherpicoInv(F, T1, T2);

	return true;
}

bool L_RelEpipolarFundamental::escalamientoIsotherpico(L_Matrix &uv, L_Matrix &T1, L_Matrix &T2)
{
	// "In Defense of the Eight-Point Algorithm", Richard Hartley
	double cx, cy, icx, icy, r2, ir;
	int i;

	T1.reallocate(3,3);
	T2.reallocate(3,3);
	T1.identity();
	T2.identity();

	// uv_nuevo[01] = T1*uv[01]
	// uv_nuevo[23] = T2*uv[23]

	// Dejar centrados los datos en zero
	cx = 0;
	cy = 0;
	for (i=0; i<uv.li; i++)
	{
		cx += uv(i,0);
		cy += uv(i,1);
	}
	icx = -cx / uv.li;
	icy = -cy / uv.li;
	for (i=0; i<uv.li; i++)
	{
		uv(i,0) += icx;
		uv(i,1) += icy;
	}
	T1(0,2) = icx;
	T1(1,2) = icy;

	//
	cx = 0;
	cy = 0;
	for (i=0; i<uv.li; i++)
	{
		cx += uv(i,2);
		cy += uv(i,3);
	}
	icx = -cx / uv.li;
	icy = -cy / uv.li;
	for (i=0; i<uv.li; i++)
	{
		uv(i,2) += icx;
		uv(i,3) += icy;
	}
	T2(0,2) = icx;
	T2(1,2) = icy;

	// Ahora el escalamiento isotherpico
	// Dejar con variance unitaria
	r2 = 0;
	for (i=0; i<uv.li; i++)
		r2 += uv(i,0)*uv(i,0) + uv(i,1)*uv(i,1);
	if (r2 == 0 || r2!=r2) // 0 o NaN
		return false;
	ir = SQRT2 / sqrt(r2/uv.li); // sqrt(2) / sqrt(mean de los r*r)
	for (i=0; i<uv.li; i++)
	{
		uv(i,0) *= ir;
		uv(i,1) *= ir;
	}
	T1(0,0) = ir;
	T1(1,1) = ir;
	T1(0,2) = ir*T1(0,2); // Esto es medio raro, pero es para que cuadre la transformacion
	T1(1,2) = ir*T1(1,2); // Esto es medio raro, pero es para que cuadre la transformacion

	//
	r2 = 0;
	for (i=0; i<uv.li; i++)
		r2 += uv(i,2)*uv(i,2) + uv(i,3)*uv(i,3);
	if (r2 == 0 || r2!=r2) // 0 o NaN
		return false;
	ir = SQRT2 / sqrt(r2/uv.li); // sqrt(2) / sqrt(mean de los r*r)
	for (i=0; i<uv.li; i++)
	{
		uv(i,2) *= ir;
		uv(i,3) *= ir;
	}
	T2(0,0) = ir;
	T2(1,1) = ir;
	T2(0,2) = ir*T2(0,2); // Esto es medio raro, pero es para que cuadre la transformacion
	T2(1,2) = ir*T2(1,2); // Esto es medio raro, pero es para que cuadre la transformacion

	return true;
}


void L_RelEpipolarFundamental::escalamientoIsotherpicoInv(L_Matrix &F, L_Matrix &T1, L_Matrix &T2)
{
	// "In Defense of the Eight-Point Algorithm", Richard Hartley
	// F_nuevo = T2T*F*T1
	L_Matrix T2T, tmp;
	T2T.transpOf(T2);

	tmp.OP_mult(T2T, F);
	F.OP_mult(tmp,T1);
	return;
}

bool L_RelEpipolarFundamental::despejar_XY_7p(L_Matrix &X, L_Matrix &Y, const L_Matrix &QagrT)
{
	L_Matrix nu;
	int i;

	nu.espacioNuloPivoteo(QagrT); // QagrT ya es de 7x9 supuestamente
	// Ortonormalizar los vectores; nu es de 9x4
	nu.gramSchmidtColumnas();
	// Copiar los vectores finales
	X.reallocate(3,3);
	Y.reallocate(3,3);
	for (i=0; i<9; i++)
	{
		X(i/3,i%3) = nu(i,0);
		Y(i/3,i%3) = nu(i,1);
	}
	return true;
}


bool L_RelEpipolarFundamental::algoritmo7puntos(L_Matrix F[3], int &nSols, const L_Matrix &x1y1x2y2)
{
	// "In Defense of the Eight-Point Algorithm", Richard Hartley
	L_Matrix a, aT, aTa, U1, U2;
	L_Matrix tmp, uv, T1, T2;
	int i;
	bool escalamiento = false;

	// 1:
	uv = x1y1x2y2; // 7 supuestamente

	// 2: Minimos cuadrados. Se busca que (uR vR) F (uL vL)T = 0
	//      [u2 u3 1]*[e00 e01 e02;e10 e11 e12;e20 e21 e22]*[u0;u1;1] = 0
	//      Maple  ->  u0*u2*e00 + u1*u2*e01 + u2*e02 + u0*u3*e10 + u1*u3*e11 + u3*e12 + u0*e20 + u1*e21 + e22 = 0

	a.reallocate(7, 9); // Lo dejo en 7 por si son mas
	for (i=0; i<7; i++)
	{
/*
		a(i,0)=uv(i,0)*uv(i,2); // E00
		a(i,1)=uv(i,1)*uv(i,2); // E01
		a(i,2)=      1      *uv(i,2); // E02
		a(i,3)=uv(i,0)*uv(i,3); // E10
		a(i,4)=uv(i,1)*uv(i,3); // E11
		a(i,5)=      1      *uv(i,3); // E12
		a(i,6)=uv(i,0)*      1      ; // E20
		a(i,7)=uv(i,1)*      1      ; // E21
		a(i,8)=      1      *      1      ; // E22
*/
		a(i,0)=uv(i,2)*uv(i,0); // E00
		a(i,1)=uv(i,3)*uv(i,0); // E01
		a(i,2)=      1      *uv(i,0); // E02
		a(i,3)=uv(i,2)*uv(i,1); // E10
		a(i,4)=uv(i,3)*uv(i,1); // E11
		a(i,5)=      1      *uv(i,1); // E12
		a(i,6)=uv(i,2)*      1      ; // E20
		a(i,7)=uv(i,3)*      1      ; // E21
		a(i,8)=      1      *      1      ; // E22

	}

	despejar_XY_7p(U1, U2, a);

	// Si U1 tiene rango 2, se devuelve de inmediato
	if (U1.det() == 0)
	{
		F[0] = U1;
		nSols = 1;
		return true;
	}

	// Elegir lambda tal que el determinante de U1+lambda*U2 = 0

	// queda: det(U1+lambda*U2) = a*lambda^3 + b*lambda^2 + c*lambda + d = 0

	// res := collect(det(U1+lambda*U2), lambda);
	// A[0] := subs(lambda=0, res);
	// A[1] := subs(lambda=0, diff(res,lambda));
	// A[2] := subs(lambda=0, diff(res,lambda,lambda) / 2);
	// A[3] := diff(res,lambda,lambda,lambda) / 6;
	double fa = U1(2,0) * U1(0,1) * U1(1,2) - U1(2,0) * U1(0,2) * U1(1,1) + U1(0,0) * U1(1,1) * U1(2,2) - U1(1,0) * U1(0,1) * U1(2,2) - U1(0,0) * U1(1,2) * U1(2,1) + U1(1,0) * U1(0,2) * U1(2,1);
	double fb = -U1(1,0) * U2(0,1) * U1(2,2) + U2(1,0) * U1(0,2) * U1(2,1) + U1(1,0) * U1(0,2) * U2(2,1) + U1(0,0) * U1(1,1) * U2(2,2) + U2(0,0) * U1(1,1) * U1(2,2) - U1(2,0) * U1(0,2) * U2(1,1) - U2(2,0) * U1(0,2) * U1(1,1) + U1(2,0) * U2(0,1) * U1(1,2) - U2(0,0) * U1(1,2) * U1(2,1) + U2(2,0) * U1(0,1) * U1(1,2) - U1(2,0) * U2(0,2) * U1(1,1) + U1(2,0) * U1(0,1) * U2(1,2) - U1(0,0) * U2(1,2) * U1(2,1) - U1(1,0) * U1(0,1) * U2(2,2) - U2(1,0) * U1(0,1) * U1(2,2) - U1(0,0) * U1(1,2) * U2(2,1) + U1(1,0) * U2(0,2) * U1(2,1) + U1(0,0) * U2(1,1) * U1(2,2);
	double fc = -U2(1,0) * U1(0,1) * U2(2,2) + U2(0,0) * U2(1,1) * U1(2,2) + U2(2,0) * U2(0,1) * U1(1,2) - U2(0,0) * U1(1,2) * U2(2,1) - U2(0,0) * U2(1,2) * U1(2,1) + U1(0,0) * U2(1,1) * U2(2,2) - U2(1,0) * U2(0,1) * U1(2,2) + U1(2,0) * U2(0,1) * U2(1,2) + U2(0,0) * U1(1,1) * U2(2,2) + U2(1,0) * U1(0,2) * U2(2,1) + U2(1,0) * U2(0,2) * U1(2,1) - U1(0,0) * U2(1,2) * U2(2,1) + U2(2,0) * U1(0,1) * U2(1,2) - U1(1,0) * U2(0,1) * U2(2,2) - U1(2,0) * U2(0,2) * U2(1,1) + U1(1,0) * U2(0,2) * U2(2,1) - U2(2,0) * U2(0,2) * U1(1,1) - U2(2,0) * U1(0,2) * U2(1,1);
	double fd = -U2(1,0) * U2(0,1) * U2(2,2) + U2(1,0) * U2(0,2) * U2(2,1) + U2(2,0) * U2(0,1) * U2(1,2) - U2(0,0) * U2(1,2) * U2(2,1) + U2(0,0) * U2(1,1) * U2(2,2) - U2(2,0) * U2(0,2) * U2(1,1);

	// Ahora encontrar lambda
	L_ComplexDouble x1, x2, x3;
	L_ComplexDouble::solveCubicEquation(fa, fb, fc, fd, x1, x2, x3);

	// Ahora ver si hay 1 o 3 reales... se admite una pequenisima parte imaginaria
	double lambda[3];
	nSols = 0;
	if (x1.isReal(0.01)) // real...
		lambda[nSols++] = x1.re;
	if (x2.isReal(0.01)) // real...
		lambda[nSols++] = x2.re;
	if (x3.isReal(0.01)) // real...
		lambda[nSols++] = x3.re;

	// Ahora calcular las soluciones
	for (i=0; i<nSols; i++)
	{
		F[i].reallocate(3,3);
		for (int j=0; j<9; j++)
			F[i](j/3,j%3) = U1(j/3,j%3) + lambda[i]*U2(j/3,j%3);
	}
	return true;
}

bool L_RelEpipolarFundamental::algoritmo7_mas_1_puntos(L_Matrix &F, const L_Matrix &x1y1x2y2)
{
	L_Matrix Fs[3];
	int nSols;
	double eMin = 1.0e60;
	double x1, y1, x2, y2;
	double ux, uy, uz, d, e1, e2, e;
	int i, iMin=-1;
	if (algoritmo7puntos(Fs, nSols, x1y1x2y2) == false)
		return false;

	for (i=0; i<nSols; i++)
	{
		x1 = x1y1x2y2(7,0);
		y1 = x1y1x2y2(7,1);
		x2 = x1y1x2y2(7,2);
		y2 = x1y1x2y2(7,3);

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
		e = sqrt(e1*e1+e2*e2);

		if (e < eMin)
		{
			eMin = e;
			iMin = i;
		}
#undef x1
#undef y1
#undef x2
#undef y2
	}
	F = Fs[iMin];
	return true;
}


void L_RelEpipolarEsencial::calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni, int nRev)
{
	int nMax = x1y1x2y2.li;
	if (nRev >= 0 && iIni + nRev < nMax)
		nMax = iIni + nRev;
	err.resize(x1y1x2y2.li);
	for (int i=0; i<nMax-iIni; i++)
		err[i] = errorRelGeom(x1y1x2y2(i+iIni,0), x1y1x2y2(i+iIni,1), x1y1x2y2(i+iIni,2), x1y1x2y2(i+iIni,3));
}

int L_RelEpipolarEsencial::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_RelEpipolarEsencial::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err;
	int i;

	if (c < 0) // Revisarlos todos rapido
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li; // iFin es una output (para saber hasta donde se llego)
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));
			corr += (err <= errMax);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li; // iFin es una output (para saber hasta donde se llego)
			err = errorRelGeom(x1y1x2y2(iUlt,0), x1y1x2y2(iUlt,1), x1y1x2y2(iUlt,2), x1y1x2y2(iUlt,3));

			if (err <= errMax)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}



int L_RelEpipolarEsencial::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_EssentialMatrix Es[10];
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_RelEpipolarEsencial *ptr;
	int k, ihip=0;
	int rango;
	int nSols;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip;
		if (rango >= x1y1x2y2.li)
			rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (L_Algoritmo5Puntos::algoritmo5Puntos(Es, nSols, xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_RelEpipolarEsencial();
			ptr->E = Es[k];
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_RelEpipolarEsencial::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_EssentialMatrix Es[10];
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_RelEpipolarEsencial *ptr;
	int k, ihip=0;
	int rango;
	int nSols;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = permIni[k];
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (L_Algoritmo5Puntos::algoritmo5Puntos(Es, nSols, xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_RelEpipolarEsencial();
			ptr->E = Es[k];
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

int L_TransfAfin2D::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_TransfAfin2D::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err2, errMax2 = errMax*errMax;
	double ux, uy;
	int i;
	if (c<0)
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li; // Para evitar overflow
			ux = m11*x1y1x2y2(iUlt,0)+m12*x1y1x2y2(iUlt,1)+tx - x1y1x2y2(iUlt,2);
			uy = m21*x1y1x2y2(iUlt,0)+m22*x1y1x2y2(iUlt,1)+ty - x1y1x2y2(iUlt,3);
			err2 = ux*ux+uy*uy;
			corr += (err2 <= errMax2);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li; // Para evitar overflow
			ux = m11*x1y1x2y2(iUlt,0)+m12*x1y1x2y2(iUlt,1)+tx - x1y1x2y2(iUlt,2);
			uy = m21*x1y1x2y2(iUlt,0)+m22*x1y1x2y2(iUlt,1)+ty - x1y1x2y2(iUlt,3);
			err2 = ux*ux+uy*uy;
			if (err2 <= errMax2)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}


int L_TransfAfin2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfAfin2D tr;
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_TransfAfin2D *ptr;
	int k, ihip=0;
	int rango;
	int nSols = 1;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip;
		if (rango >= x1y1x2y2.li)
			rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (tr.calcExacto(xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfAfin2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_TransfAfin2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfAfin2D tr;
	std::vector<int> perm(gradosLibertad());
	L_Matrix xyxy;
	L_TransfAfin2D *ptr;
	int k, ihip=0;
	int rango;
	int nSols = 1;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = x1y1x2y2.li-1;
		L_randomSelection(perm, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (tr.calcExacto(xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfAfin2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

bool L_TransfAfin2D::calcExacto(const L_Matrix &x1y1x2y2)
{
	L_Matrix A,x,b;
	int i;
	A.reallocate(6,6);
	b.reallocate(6,1);
	for (i=0; i<3; i++)
	{
		A(2*i+0,0)=x1y1x2y2(i,0);
		A(2*i+0,1)=x1y1x2y2(i,1);
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=x1y1x2y2(i,0);
		A(2*i+1,4)=x1y1x2y2(i,1);
		A(2*i+1,5)=1;

		b(2*i+0,0)=x1y1x2y2(i,2);
		b(2*i+1,0)=x1y1x2y2(i,3);
	}
	b.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(A) == false) // Ahora b = solucion de (Ax=b)
		return false;

	m11=x(0,0);
	m12=x(1,0);
	tx=x(2,0);
	m21=x(3,0);
	m22=x(4,0);
	ty=x(5,0);
	return true;
}



bool L_TransfAfin2D::calcMinCuad(const L_Matrix &x1y1x2y2)
{
	L_Matrix A,AT,ATA,ATb,x,b;
	int i;

	A.reallocate(2*x1y1x2y2.li,6);
	b.reallocate(2*x1y1x2y2.li,1);
	for (i=0; i<x1y1x2y2.li; i++)
	{
		A(2*i+0,0)=x1y1x2y2(i,0);
		A(2*i+0,1)=x1y1x2y2(i,1);
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=x1y1x2y2(i,0);
		A(2*i+1,4)=x1y1x2y2(i,1);
		A(2*i+1,5)=1;

		b(2*i+0,0)=x1y1x2y2(i,2);
		b(2*i+1,0)=x1y1x2y2(i,3);
	}
	AT.transpOf(A);
	ATA.OP_mult(AT,A);
	ATb.OP_mult(AT,b);

	ATb.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(ATA) == false)
		return false;

	m11=x(0,0);
	m12=x(1,0);
	tx=x(2,0);
	m21=x(3,0);
	m22=x(4,0);
	ty=x(5,0);
	return true;
}


void L_TransfAfin2D::calcMatricesRayos(const L_CamaraPinhole &cR, const L_CamaraPinhole &cP, L_Matrix &M, L_Matrix &t)
{
	M.reallocate(2,2);
	t.reallocate(2,1);
	// (xP) = (m11 m12)(xR) + (tX)
	// (yP)   (m21 m22)(yR)   (tY)
	
	// xP = a11 + b11*uP, xR = a12 + b12*uR;
	// yP = a21 + b21*vP, yR = a22 + b22*vR

	// (a11 + b11*uP) = (m11 m12)(a12 + b12*uR) + (tX)
	// (a21 + b21*vP)   (m21 m22)(a22 + b22*vR)   (tY)

	// uP = m11*b12/b11*uR + m12*b22/b11*vR + (-a11+m11*a12+m12*a22+tX)/b11
	// vP = m21*b12/b21*uR + m22*b22/b21*vR + (-a21+m21*a12+m22*a22+tY)/b21
	double a11, a12, a21, a22, b11, b12, b21, b22;

	a11 = cP.cenX; a12 = cR.cenX;
	a21 = cP.cenY; a22 = cR.cenY;
	b11 = -cP.dfX; b12 = -cR.dfX;
	b21 = -cP.dfY; b22 = -cR.dfY;

	M(0,0) = m11*b12/b11; M(0,1) = m12*b22/b11; t(0,0) = (-a11+m11*a12+m12*a22+tx)/b11;
	M(1,0) = m21*b12/b21; M(1,1) = m22*b22/b21; t(1,0) = (-a21+m21*a12+m22*a22+ty)/b21;
}


void L_TransfAfin2D::descomponeSVD(double &escMax, double &escMin, double &rot, const L_CamaraPinhole *cR, const L_CamaraPinhole *cP)
{
	L_Matrix M(2,2), t(2,1);
	if (cR != NULL && cP!=NULL)
		calcMatricesRayos(*cR, *cP, M, t);
	else
	{
		M(0,0) = m11;
		M(0,1) = m12;
		M(1,0) = m21;
		M(1,1) = m22;
	}
	M.svd2x2(rot, escMax, escMin);
	rot = -rot; // Por el sistema de ejes que uso
}

bool L_TransfAfin2D::invertMe()
{
	L_Matrix A(3,3);
	A(0,0) = m11;
	A(0,1) = m12;
	A(0,2) = tx;
	A(1,0) = m21;
	A(1,1) = m22;
	A(1,2) = ty;
	A(2,0) = 0;
	A(2,1) = 0;
	A(2,2) = 1;
	if (A.invertMe() == false)
		return false;
	m11 = A(0,0);
	m12 = A(0,1);
	tx  = A(0,2);
	m21 = A(1,0);
	m22 = A(1,1);
	ty  = A(1,2);
	return true;
}

void L_TransfSimil2D::calcVectorError(std::vector<double> &err, const L_Matrix &x1y1x2y2, int iIni, int nRev)
{
	int nMax = x1y1x2y2.li;
	if (nRev >= 0 && iIni + nRev < nMax)
		nMax = iIni + nRev;
	double resx, resy, dx, dy;
	err.resize(x1y1x2y2.li);
	for (int i=iIni; i<nMax; i++)
	{
		resx = m11*x1y1x2y2(i,0) + m12*x1y1x2y2(i,1) + tx;
		resy = m21*x1y1x2y2(i,0) + m22*x1y1x2y2(i,1) + ty;
		dx = resx - x1y1x2y2(i,2);
		dy = resy - x1y1x2y2(i,3);
		err[i] = sqrt(dx*dx + dy*dy);
	}
}

int L_TransfSimil2D::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_TransfSimil2D::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err2;
	double ux, uy;
	int i;
	double errMax2 = errMax*errMax;
	if (c < 0)
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			ux = m11*x1y1x2y2(iUlt,0)+m12*x1y1x2y2(iUlt,1)+tx - x1y1x2y2(iUlt,2);
			uy = m21*x1y1x2y2(iUlt,0)+m22*x1y1x2y2(iUlt,1)+ty - x1y1x2y2(iUlt,3);
			err2 = ux*ux+uy*uy;
			corr += (err2 <= errMax2);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			ux = m11*x1y1x2y2(iUlt,0)+m12*x1y1x2y2(iUlt,1)+tx - x1y1x2y2(iUlt,2);
			uy = m21*x1y1x2y2(iUlt,0)+m22*x1y1x2y2(iUlt,1)+ty - x1y1x2y2(iUlt,3);
			err2 = ux*ux+uy*uy;
			if (err2 <= errMax2)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}


int L_TransfSimil2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfSimil2D tr;
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_TransfSimil2D *ptr;
	int k, ihip=0;
	int rango;
	int iter = 0;
	int nSols = 1;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip + iter;
		if (rango >= x1y1x2y2.li)
			rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		iter++;
		if (tr.calcExacto(xyxy)==true)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfSimil2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_TransfSimil2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfSimil2D tr;
	std::vector<int> perm(gradosLibertad());
	L_Matrix xyxy;
	L_TransfSimil2D *ptr;
	int k, ihip=0;
	int rango;
	int nSols = 1;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = x1y1x2y2.li-1;
		L_randomSelection(perm, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (tr.calcExacto(xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfSimil2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

bool L_TransfSimil2D::calcExacto(const L_Matrix &x1y1x2y2)
{
	double xc, yc, dx, dy;
	double xcT, ycT, dxT, dyT;
	double rot, rotT;
	double r, rT;
	// Se considera [0][0],[0][1] ; [1][0],[1][1] la primera referencia
	// Se considera [0][2],[0][3] ; [1][2],[1][3] la primera referencia

	xc = x1y1x2y2(0,0);
	yc = x1y1x2y2(0,1);
	dx = x1y1x2y2(1,0) - xc;
	dy = x1y1x2y2(1,1) - yc;

	xcT = x1y1x2y2(0,2);
	ycT = x1y1x2y2(0,3);
	dxT = x1y1x2y2(1,2) - xcT;
	dyT = x1y1x2y2(1,3) - ycT;

	rot = -atan2(dy, dx); // Se mide invertido
	rotT = -atan2(dyT, dxT);
	ang = rotT-rot;

	r = sqrt(dx*dx + dy*dy);
	rT = sqrt(dxT*dxT + dyT*dyT);
	if (r == 0)
		return false; // Escala indefinida
	esc = rT / r;
	m11 = esc * cos(-ang); // Se mide invertido
	m21 = esc * sin(-ang);
	m12 = -m21;
	m22 = m11;

	tx = xcT - (m11*xc + m12*yc);
	ty = ycT - (m21*xc + m22*yc);
	return true;
}


bool L_TransfSimil2D::invertMe()
{
	double c, s;
	double tx_, ty_;
	if (esc == 0)
		return false;
	esc = 1/esc;
	ang = -ang;
	c = cos(ang);
	s = sin(ang);
	tx_ = -esc*c*tx + esc*s*ty; // REVISAR ESTA FUNCION Y SVD2x2
	ty_ = -esc*s*tx - esc*c*ty;
	tx = tx_;
	ty = ty_;
	m11 = esc*c;
	m21 = esc*s;
	m12 = -m21;
	m22 = m11;
	return true;
}

bool L_TransfSimil2D::calcMinCuad(const L_Matrix &x1y1x2y2)
{
	L_Matrix A,AT,ATA,ATb,x,b, Matr(2,2);
	int i;
	double e1, e2, ang1, ang2;

	A.reallocate(2*x1y1x2y2.li,6);
	b.reallocate(2*x1y1x2y2.li,1);
	for (i=0; i<x1y1x2y2.li; i++)
	{
		A(2*i+0,0)=x1y1x2y2(i,0);
		A(2*i+0,1)=x1y1x2y2(i,1);
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=x1y1x2y2(i,0);
		A(2*i+1,4)=x1y1x2y2(i,1);
		A(2*i+1,5)=1;

		b(2*i+0,0)=x1y1x2y2(i,2);
		b(2*i+1,0)=x1y1x2y2(i,3);
	}
	AT.transpOf(A);
	ATA.OP_mult(AT,A);
	ATb.OP_mult(AT,b);

	ATb.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(ATA) == false)
		return false;

	m11 = x(0,0);
	m12 = x(1,0);
	tx  = x(2,0);
	m21 = x(3,0);
	m22 = x(4,0);
	ty  = x(5,0);
	Matr(0,0) = m11;
	Matr(0,1) = m12;
	Matr(1,0) = m21;
	Matr(1,1) = m22;
	Matr.svd2x2(ang1,ang2,e1,e2);
	ang = ang1 + ang2;
	esc = sqrt(e1*e2); // Es lo mejor que se me ocurre
	m11 = esc*cos(ang);
	m21 = esc*sin(ang);
	m12 = -m21;
	m22 = m11;
	return true;
}


int L_TransfProyectiva2D::test_c_de_d(const L_Matrix &x1y1x2y2, double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > x1y1x2y2.li, "L_TransfProyectiva2D::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	double err2, errMax2 = errMax*errMax;
	double resx, resy, resz, dx, dy;
	int i;
	if (c<0)
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			resx=m11*x1y1x2y2(iUlt,0)+m12*m11*x1y1x2y2(iUlt,1)+tx;
			resy=m21*m11*x1y1x2y2(iUlt,0)+m22*m11*x1y1x2y2(iUlt,1)+ty;
			resz = m31*m11*x1y1x2y2(iUlt,0)+m32*m11*x1y1x2y2(iUlt,1) + 1;
			dx = (resx-x1y1x2y2(iUlt,2))/resz;
			dy = (resy-x1y1x2y2(iUlt,3))/resz;
			err2 = dx*dx + dy*dy;

			corr += (err2 <= errMax2);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (i+iIni)%x1y1x2y2.li;
			resx=m11*x1y1x2y2(iUlt,0)+m12*m11*x1y1x2y2(iUlt,1)+tx;
			resy=m21*m11*x1y1x2y2(iUlt,0)+m22*m11*x1y1x2y2(iUlt,1)+ty;
			resz = m31*m11*x1y1x2y2(iUlt,0)+m32*m11*x1y1x2y2(iUlt,1) + 1;
			dx = (resx-x1y1x2y2(iUlt,2))/resz;
			dy = (resy-x1y1x2y2(iUlt,3))/resz;
			err2 = dx*dx + dy*dy;

			if (err2 <= errMax2)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (incorr > d-c)
					return corr;
			}
		}
	}
	return corr;
}


int L_TransfProyectiva2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfProyectiva2D tr;
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_TransfProyectiva2D *ptr;
	int k, ihip=0;
	int rango;
	int nSols = 1;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip;
		if (rango >= x1y1x2y2.li)
			rango = x1y1x2y2.li-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (tr.calcExacto(xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfProyectiva2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_TransfProyectiva2D::generarHipotesis(L_Array<L_RelGeom2D*> &hip, const L_Matrix &x1y1x2y2)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfProyectiva2D tr;
	std::vector<int> perm(gradosLibertad());
	L_Matrix xyxy;
	L_TransfProyectiva2D *ptr;
	int k, ihip=0;
	int rango;
	int nSols = 1;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = x1y1x2y2.li-1;
		L_randomSelection(perm, rango);
		for (k=0; k<(int)perm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (tr.calcExacto(xyxy) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfProyectiva2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

int L_TransfPoseProy2D::test_c_de_d(double errMax, int c, int d, int iIni, int &iUlt)
{
	throw_L_ArgException_if(iIni+d > (int)pun.size(), "L_TransfPoseProy2D::test_c_de_d() d fuera del arreglo");
	int corr=0, incorr=0;
	L_RayoCamara rProy;
	double xs, errMax2, errx, erry, err2;
	int i;
	errMax2 = errMax*errMax;
	if (c<0)
	{
		for (i=0; i<d; i++)
		{
			iUlt = (int)( (i+iIni)%pun.size() );
			pose.prmovi_v(rProy, pun[i], xs);
			errx = rProy.tanIzq - ray[iUlt].tanIzq;
			erry = rProy.tanArr - ray[iUlt].tanArr;
			err2 = errx*errx + erry*erry;

			corr += (xs>0 && err2 <= errMax2);
		}
	}
	else
	{
		for (i=0; i<d; i++)
		{
			iUlt = (int)( (i+iIni)%pun.size() );
			pose.prmovi_v(rProy, pun[i], xs);
			errx = rProy.tanIzq - ray[iUlt].tanIzq;
			erry = rProy.tanArr - ray[iUlt].tanArr;
			err2 = errx*errx + erry*erry;
			if (xs>0 && err2 <= errMax2)
			{
				// Para devolver true, tienen que haber c+ correctos de d posibles
				corr++;
				if (c>0 && corr >=c)
					return corr;
			}
			else
			{
				// Para devolver false, no tienen que haber c+ correctos en el rango 0,...,d-1
				incorr++;
				if (c>0 && incorr > d-c)
					return corr;
			}
		}
	}

	return corr;
}



int L_TransfPoseProy2D::generarHipotesis(L_Array<L_RelGeomMuestra2D*> &hip, const std::vector<double> &puntajes)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfPoseProy2D tr;
	std::vector<int> indices;
	std::vector<int> perm(gradosLibertad()), permIni(gradosLibertad());
	L_Matrix xyxy;
	L_TransfPoseProy2D *ptr;
	int k, ihip=0;
	int rango, tamD;
	int nSols = 1;

	L_getOrderedIndexesMajorMinor(indices, puntajes);
	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	tamD = tamDatos();
	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = gradosLibertad() + ihip;
		if (rango >= tamD)
			rango = tamD-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
		if (calcExacto(perm) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfPoseProy2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}


int L_TransfPoseProy2D::generarHipotesis(L_Array<L_RelGeomMuestra2D*> &hip)
{
	// hip.size() ya debe existir, e indica la cantidad maxima
	L_TransfPoseProy2D tr;
	std::vector<int> perm(gradosLibertad());
	L_Matrix xyxy;
	L_TransfPoseProy2D *ptr;
	int k, ihip=0;
	int rango, tamD;
	int nSols = 1;

	xyxy.reallocate(gradosLibertad(), 4);

	for (ihip=0; ihip<(int)hip.size(); ihip++)
		hip[ihip] = NULL;

	tamD = tamDatos();
	for (ihip=0; ihip<(int)hip.size();  )
	{
		rango = tamD-1;
		L_randomSelection(perm, rango);
		if (calcExacto(perm) == false)
			continue;
		for (k=0; k<nSols; k++)
		{
			if (ihip >= (int)hip.size())
				break;
			ptr = new L_TransfPoseProy2D();
			*ptr = tr;
			hip[ihip] = ptr;
			ihip++;
		}
	}
	return ihip;
}

bool L_TransfProyectiva2D::calcExacto(const L_Matrix &x1y1x2y2)
{
	L_Matrix A,x,b;
	int i;
	A.reallocate(8,8);
	b.reallocate(8,1);
	for (i=0; i<4; i++)
	{
		A(2*i+0,0)=x1y1x2y2(i,0);
		A(2*i+0,1)=x1y1x2y2(i,1);
		A(2*i+0,2)=1;
		A(2*i+0,3)=0;
		A(2*i+0,4)=0;
		A(2*i+0,5)=0;

		A(2*i+1,0)=0;
		A(2*i+1,1)=0;
		A(2*i+1,2)=0;
		A(2*i+1,3)=x1y1x2y2(i,0);
		A(2*i+1,4)=x1y1x2y2(i,1);
		A(2*i+1,5)=1;

		b(2*i+0,0)=x1y1x2y2(i,2);
		b(2*i+1,0)=x1y1x2y2(i,3);
	}
	b.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(A) == false) // Ahora b = solucion de (Ax=b)
		return false;

	m11=x(0,0);
	m12=x(1,0);
	tx=x(2,0);
	m21=x(3,0);
	m22=x(4,0);
	ty=x(5,0);
	return true;
}

bool L_TransfProyectiva2D::calcMinCuad(const L_Matrix &x1y1x2y2)
{
	L_Matrix A,x,b;
	int i;
	A.reallocate(2*x1y1x2y2.li,8);
	b.reallocate(2*x1y1x2y2.li,1);
	for (i=0; i<x1y1x2y2.li; i++)
	{
		A(2*i,0) = x1y1x2y2(i,0);
		A(2*i,1) = x1y1x2y2(i,1);
		A(2*i,2) = 1;
		A(2*i,3) = 0;
		A(2*i,4) = 0;
		A(2*i,5) = 0;
		A(2*i,6) = -x1y1x2y2(i,0)*x1y1x2y2(i,2);
		A(2*i,7) = -x1y1x2y2(i,1)*x1y1x2y2(i,2);
		b(2*i,0) = x1y1x2y2(i,2);
		//
		A(2*i+1,0) = 0;
		A(2*i+1,1) = 0;
		A(2*i+1,2) = 0;
		A(2*i+1,3) = x1y1x2y2(i,0);
		A(2*i+1,4) = x1y1x2y2(i,1);
		A(2*i+1,5) = 1;
		A(2*i+1,6) = -x1y1x2y2(i,0)*x1y1x2y2(i,3);
		A(2*i+1,7) = -x1y1x2y2(i,1)*x1y1x2y2(i,3);
		b(2*i+1,0) = x1y1x2y2(i,3);
	}
	b.swap(x);
	if (x.solveLinearSystem_x_Ab_result_in_b(A) == false) // Ahora b = solucion de (Ax=b)
		return false;

	for (i=0; i<8; i++)
		el(i)=x(i,0);
	return true;
}

bool L_TransfProyectiva2D::invertMe()
{
	L_Matrix A(3,3);
	A(0,0) = m11;
	A(0,1) = m12;
	A(0,2) = tx;
	A(1,0) = m21;
	A(1,1) = m22;
	A(1,2) = ty;
	A(2,0) = 0;
	A(2,1) = 0;
	A(2,2) = 1;
	if (A.invertMe() == false)
		return false;
	m11 = A(0,0);
	m12 = A(0,1);
	tx  = A(0,2);
	m21 = A(1,0);
	m22 = A(1,1);
	ty  = A(1,2);
	return true;
}



bool L_TransfPoseProy2D::calcExacto(const std::vector<int> &selecc)
{
	L_Pose3D_cuat poseInv;
	int i;
	for (i=0; i<4; i++)
	{
		punSel[i] = pun[selecc[i]];
		raySel[i] = ray[selecc[i]];
	}
	poseInv.inverseOf(pose);
	if (poseInv.calcPose_movido_a_fijo_4rayos(punSel,raySel) == false) // Entrega la inversa de la que se desea
		return false;
	pose.inverseOf(poseInv);
	return true;
}

bool L_TransfPoseProy2D::calcMinCuad()
{
	L_Pose3D_cuat poseInv;
	L_LevenbergMarquardt opt;
	opt.epsilon = 1e-6;
	opt.errorToStop = 1e-10;
	opt.nIterationsMax = 10;
	opt.lenVector = 7;
	poseInv.inverseOf(pose);
	poseInv.calcPose_movido_a_fijo_Nrayos(pun,ray,opt); // Entrega la inversa de la que se desea
	pose.inverseOf(poseInv);
	return true;
}


//! Esta funcion erase_preserving_order la escala y la distorsion afin de las componentes de rotacion de la matriz homogenea
double L_Rot3DMatrix::normalizaSVD(double *eMax, double *eMed, double *eMin)
{
	//void L_HomogeneousMatrix::svdcmp(double **a, int m, int n, double w[), double **v)
	// Given a matrix a[1..m,1__n), this routine computes its singular value decomposition, A = U·W·VT.
	// The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1__n).
	// The matrix V (not the transpose VT ) is output as v[1__n,1__n).
	L_Matrix u, d, v, vT, rot;
	int i, j;
	u.reallocate(3,3);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			u(i,j)=operator()(i,j);
	u.svd(d, v);
	vT.transpOf(v);
	rot.OP_mult(u,vT);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			operator()(i,j)=rot(i,j);

	if (eMax != NULL)
		*eMax = d(0,0);
	if (eMed != NULL)
		*eMed = d(1,0);
	if (eMin != NULL)
		*eMin = d(2,0);
	double e1, e2, e3, error;
	e1=log(L_ABS(d(0,0)));
	e2=log(L_ABS(d(1,0)));
	e3=log(L_ABS(d(2,0)));
	error=e1*e1+e2*e2+e3*e3;
	return error;
}

double L_Rot3DMatrix::normalizaSVD_corrigeJ(L_Matrix &J9xN)
{
	//void L_HomogeneousMatrix::svdcmp(double **a, int m, int n, double w[), double **v)
	// Given a matrix a[1..m,1__n), this routine computes its singular value decomposition, A = U·W·VT.
	// The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[1__n).
	// The matrix V (not the transpose VT ) is output as v[1__n,1__n).
	L_Matrix u, d, v, vT, rot;
	int i, j, k;
	u.reallocate(3,3);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			u(i,j)=operator()(i,j);
	u.svd(d, v);
	vT.transpOf(v);
	rot.OP_mult(u,vT);
	// Intento de arreglar la matriz J12xN
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			for (k=0; k<J9xN.lj; k++)
				J9xN(3*i+j,k) = J9xN(3*i+j,k) * (operator()(i,j)!=0 ? (rot(i,j) / operator()(i,j)) : 0);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			operator()(i,j)=rot(i,j);

	double e1, e2, e3, error;
	e1=log(L_ABS(d(0,0)));
	e2=log(L_ABS(d(1,0)));
	e3=log(L_ABS(d(2,0)));
	error=e1*e1+e2*e2+e3*e3;
	return error;
}


bool L_CalculadorPoseTriangulo::metodoDirecto()
{
	double a2, b2, c2;
	double cosAlfa, cosBeta, cosGamma;
	double cos2Alfa, cos2Beta, cos2Gamma;
	double sin2Alfa, sin2Beta, sin2Gamma;

	if (b==0)
		return false;

	j1 = q1.pideUnitario();
	j2 = q2.pideUnitario();
	j3 = q3.pideUnitario();

	a2 = a*a;
	b2 = b*b;
	c2 = c*c;

	cosAlfa = j2.punto(j3);
	cosBeta = j1.punto(j3);
	cosGamma = j1.punto(j2);
	cos2Alfa = cosAlfa*cosAlfa;
	cos2Beta = cosBeta*cosBeta;
	cos2Gamma = cosGamma*cosGamma;
	sin2Alfa = 1 - cos2Alfa;
	sin2Beta = 1 - cosBeta*cosBeta;
	sin2Gamma = 1 - cos2Gamma;
	double A = (a2-c2)/b2;
	double A4 = (A - 1)*(A - 1) - 4*c2/b2*cos2Alfa;
	double A3 = 4 * ( A*(1-A)*cosBeta - (1 - (a2+c2)/b2)*cosAlfa*cosGamma + 2*c2/b2*cos2Alfa*cosBeta );
	double A2 = 2 * ( A*A - 1 + 2*A*A*cos2Beta + 2*(b2-c2)/b2*cos2Alfa
		- 4*(a2+c2)/b2*cosAlfa*cosBeta*cosGamma + 2*(b2-a2)/b2*cos2Gamma );
	double A1 = 4 * ( -A*(1+A)*cosBeta + 2*a2/b2*cos2Gamma*cosBeta - (1-(a2+c2)/b2)*cosAlfa*cosGamma );
	double A0 = (1+A)*(1+A) - 4*a2/b2*cos2Gamma;

	L_ComplexDouble v[4];
	double u[4];
	int i;
	if (L_ComplexDouble::solveQuarticEquation(A4, A3, A2, A1, A0, v[0], v[1], v[2], v[3]) == -1)
		return false;
	for (i=0; i<4; i++)
	{
		if (-1.0e-4 < v[i].im && v[i].im < 1.0e-4 && (cosGamma - v[i].re*cosAlfa)!=0 && (u[i]*u[i] + v[i].re*v[i].re - 2*u[i]*v[i].re*cosAlfa)!=0)
		{
			valida[i] = true;
			u[i] = ( (-1+A )*v[i].re*v[i].re - 2*A*cosBeta*v[i].re + 1 + A ) / (2*(cosGamma - v[i].re*cosAlfa));
			s1[i] = a / sqrt(u[i]*u[i] + v[i].re*v[i].re - 2*u[i]*v[i].re*cosAlfa);
			s2[i] = u[i] * s1[i];
			s3[i] = v[i].re * s1[i];
		}
		else
		{
			valida[i] = false;
			s1[i]=s2[i]=s3[i]=sqrt(-1.0);
		}

	}
	return true;
}


void L_Quaternion::jacob_izq_rotarVector(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0, int j0) const
{
	double t1 = c*c;
	double t2 = a*a;
	double t3 = b*b;
	double t4 = d*d;
	double t5 = t2 + t3 + t1 + t4;
	double t6 = t5 * t5;
	double t7 = 0.1e1 / t6;
	double t8 = t1 * t7;
	double t9 = t8 * a;
	double t10 = t4 * t7;
	double t11 = t10 * a;
	double t15 = t7 * a;
	double t17 = 0.2e1 * b * c * t15;
	double t18 = 0.1e1 / t5;
	double t19 = d * t18;
	double t22 = 0.2e1 * t2 * d * t7;
	double t25 = c * t18;
	double t28 = 0.2e1 * t2 * c * t7;
	double t29 = b * d;
	double t31 = 0.2e1 * t29 * t15;
	double t35 = t8 * b;
	double t36 = t10 * b;
	double t40 = t3 * c * t7;
	double t41 = 0.2e1 * t40;
	double t45 = t3 * d * t7;
	double t46 = 0.2e1 * t45;
	double t51 = t1 * c * t7;
	double t52 = t10 * c;
	double t55 = b * t18;
	double t56 = 0.2e1 * t35;
	double t58 = t7 * c;
	double t60 = 0.2e1 * a * d * t58;
	double t63 = a * t18;
	double t64 = 0.2e1 * t9;
	double t66 = 0.2e1 * t29 * t58;
	double t70 = t8 * d;
	double t72 = t4 * d * t7;
	double t75 = 0.2e1 * t11;
	double t78 = 0.2e1 * t36;
	double t85 = t3 * t7 * a;
	double t90 = 0.2e1 * t2 * b * t7;
	double t97 = t3 * b * t7;
	double t100 = 0.2e1 * t85;
	double t108 = 0.2e1 * t70;
	double t116 = 0.2e1 * t52;

	J3x4(0+i0,0+j0) =0.4e1 * (t9 + t11) * inicial.x + 0.2e1 * (-t17 - t19 + t22) * inicial.y + 0.2e1 * (t25 - t28 - t31) * inicial.z;
	J3x4(0+i0,1+j0) =0.4e1 * (t35 + t36) * inicial.x + 0.2e1 * (t25 - t41 + t31) * inicial.y + 0.2e1 * (-t17 + t19 - t46) * inicial.z;
	J3x4(0+i0,2+j0) =0.4e1 * (-t25 + t51 + t52) * inicial.x + 0.2e1 * (t55 - t56 + t60) * inicial.y + 0.2e1 * (t63 - t64 - t66) * inicial.z;
	J3x4(0+i0,3+j0) =0.4e1 * (t70 - t19 + t72) * inicial.x + 0.2e1 * (-t66 - t63 + t75) * inicial.y + 0.2e1 * (-t60 + t55 - t78) * inicial.z;
	J3x4(1+i0,0+j0) =0.2e1 * (t19 - t22 - t17) * inicial.x + 0.4e1 * (t85 + t11) * inicial.y + 0.2e1 * (-t60 - t55 + t90) * inicial.z;
	J3x4(1+i0,1+j0) =0.2e1 * (-t31 + t25 - t41) * inicial.x + 0.4e1 * (-t55 + t97 + t36) * inicial.y + 0.2e1 * (-t66 - t63 + t100) * inicial.z;
	J3x4(1+i0,2+j0) =0.2e1 * (-t60 + t55 - t56) * inicial.x + 0.4e1 * (t40 + t52) * inicial.y + 0.2e1 * (t19 - t108 + t17) * inicial.z;
	J3x4(1+i0,3+j0) =0.2e1 * (t63 - t75 - t66) * inicial.x + 0.4e1 * (t45 - t19 + t72) * inicial.y + 0.2e1 * (t25 - t116 + t31) * inicial.z;
	J3x4(2+i0,0+j0) =0.2e1 * (-t31 - t25 + t28) * inicial.x + 0.2e1 * (t55 - t90 - t60) * inicial.y + 0.4e1 * (t85 + t9) * inicial.z;
	J3x4(2+i0,1+j0) =0.2e1 * (t19 - t46 + t17) * inicial.x + 0.2e1 * (t63 - t100 - t66) * inicial.y + 0.4e1 * (-t55 + t97 + t35) * inicial.z;
	J3x4(2+i0,2+j0) =0.2e1 * (-t66 - t63 + t64) * inicial.x + 0.2e1 * (-t17 + t19 - t108) * inicial.y + 0.4e1 * (t40 - t25 + t51) * inicial.z;
	J3x4(2+i0,3+j0) =0.2e1 * (t55 - t78 + t60) * inicial.x + 0.2e1 * (-t31 + t25 - t116) * inicial.y + 0.4e1 * (t45 + t70) * inicial.z;
}



void L_Quaternion::jacob_izq_rotarVector_pre(double *pre) const
{
	double t1 = c*c;
	double t2 = a*a;
	double t3 = b*b;
	double t4 = d*d;
	double t5 = t2 + t3 + t1 + t4;
	double t6 = t5 * t5;
	double t7 = 0.1e1 / t6;
	double t8 = t1 * t7;
	double t9 = t8 * a;
	double t10 = t4 * t7;
	double t11 = t10 * a;
	double t15 = t7 * a;
	double t17 = 0.2e1 * b * c * t15;
	double t18 = 0.1e1 / t5;
	double t19 = d * t18;
	double t22 = 0.2e1 * t2 * d * t7;
	double t25 = c * t18;
	double t28 = 0.2e1 * t2 * c * t7;
	double t29 = b * d;
	double t31 = 0.2e1 * t29 * t15;
	double t35 = t8 * b;
	double t36 = t10 * b;
	double t40 = t3 * c * t7;
	double t41 = 0.2e1 * t40;
	double t45 = t3 * d * t7;
	double t46 = 0.2e1 * t45;
	double t51 = t1 * c * t7;
	double t52 = t10 * c;
	double t55 = b * t18;
	double t56 = 0.2e1 * t35;
	double t58 = t7 * c;
	double t60 = 0.2e1 * a * d * t58;
	double t63 = a * t18;
	double t64 = 0.2e1 * t9;
	double t66 = 0.2e1 * t29 * t58;
	double t70 = t8 * d;
	double t72 = t4 * d * t7;
	double t75 = 0.2e1 * t11;
	double t78 = 0.2e1 * t36;
	double t85 = t3 * t7 * a;
	double t90 = 0.2e1 * t2 * b * t7;
	double t97 = t3 * b * t7;
	double t100 = 0.2e1 * t85;
	double t108 = 0.2e1 * t70;
	double t116 = 0.2e1 * t52;


	pre[0] = 0.4e1 * (t9 + t11);
	pre[1] = 0.2e1 * (-t17 - t19 + t22);
	pre[2] = 0.2e1 * (t25 - t28 - t31);
	pre[3] = 0.4e1 * (t35 + t36);
	pre[4] = 0.2e1 * (t25 - t41 + t31);
	pre[5] = 0.2e1 * (-t17 + t19 - t46);
	pre[6] = 0.4e1 * (-t25 + t51 + t52);
	pre[7] = 0.2e1 * (t55 - t56 + t60);
	pre[8] = 0.2e1 * (t63 - t64 - t66);
	pre[9] = 0.4e1 * (t70 - t19 + t72);
	pre[10] = 0.2e1 * (-t66 - t63 + t75);
	pre[11] = 0.2e1 * (-t60 + t55 - t78);
	pre[12] = 0.2e1 * (t19 - t22 - t17);
	pre[13] = 0.4e1 * (t85 + t11);
	pre[14] = 0.2e1 * (-t60 - t55 + t90);
	pre[15] = 0.2e1 * (-t31 + t25 - t41);
	pre[16] = 0.4e1 * (-t55 + t97 + t36);
	pre[17] = 0.2e1 * (-t66 - t63 + t100);
	pre[18] = 0.2e1 * (-t60 + t55 - t56);
	pre[19] = 0.4e1 * (t40 + t52);
	pre[20] = 0.2e1 * (t19 - t108 + t17);
	pre[21] = 0.2e1 * (t63 - t75 - t66);
	pre[22] = 0.4e1 * (t45 - t19 + t72);
	pre[23] = 0.2e1 * (t25 - t116 + t31);
	pre[24] = 0.2e1 * (-t31 - t25 + t28);
	pre[25] = 0.2e1 * (t55 - t90 - t60);
	pre[26] = 0.4e1 * (t85 + t9);
	pre[27] = 0.2e1 * (t19 - t46 + t17);
	pre[28] = 0.2e1 * (t63 - t100 - t66);
	pre[29] = 0.4e1 * (-t55 + t97 + t35);
	pre[30] = 0.2e1 * (-t66 - t63 + t64);
	pre[31] = 0.2e1 * (-t17 + t19 - t108);
	pre[32] = 0.4e1 * (t40 - t25 + t51);
	pre[33] = 0.2e1 * (t55 - t78 + t60);
	pre[34] = 0.2e1 * (-t31 + t25 - t116);
	pre[35] = 0.4e1 * (t45 + t70);
}


void L_Quaternion::jacob_izq_rotarVectorInv(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0, int j0) const
{
	double t1 = c*c;
	double t2 = a*a;
	double t3 = b*b;
	double t4 = d*d;
	double t5 = t2 + t3 + t1 + t4;
	double t6 = t5 * t5;
	double t7 = 0.1e1 / t6;
	double t8 = t1 * t7;
	double t9 = t8 * a;
	double t10 = t4 * t7;
	double t11 = t10 * a;
	double t15 = t7 * a;
	double t17 = 0.2e1 * b * c * t15;
	double t18 = 0.1e1 / t5;
	double t19 = d * t18;
	double t22 = 0.2e1 * t2 * d * t7;
	double t25 = c * t18;
	double t28 = 0.2e1 * t2 * c * t7;
	double t29 = b * d;
	double t31 = 0.2e1 * t29 * t15;
	double t35 = t8 * b;
	double t36 = t10 * b;
	double t40 = t3 * c * t7;
	double t41 = 0.2e1 * t40;
	double t45 = t3 * d * t7;
	double t46 = 0.2e1 * t45;
	double t51 = t1 * c * t7;
	double t52 = t10 * c;
	double t55 = b * t18;
	double t56 = 0.2e1 * t35;
	double t58 = t7 * c;
	double t60 = 0.2e1 * a * d * t58;
	double t63 = a * t18;
	double t64 = 0.2e1 * t9;
	double t66 = 0.2e1 * t29 * t58;
	double t70 = t8 * d;
	double t72 = t4 * d * t7;
	double t75 = 0.2e1 * t11;
	double t78 = 0.2e1 * t36;
	double t85 = t3 * t7 * a;
	double t90 = 0.2e1 * t2 * b * t7;
	double t97 = t3 * b * t7;
	double t100 = 0.2e1 * t85;
	double t108 = 0.2e1 * t70;
	double t116 = 0.2e1 * t52;
	J3x4(0+i0,0+j0) = 0.4e1 * (t9 + t11) * inicial.x + 0.2e1 * (-t17 + t19 - t22) * inicial.y + 0.2e1 * (-t25 + t28 - t31) * inicial.z;
	J3x4(0+i0,1+j0) = 0.4e1 * (t35 + t36) * inicial.x + 0.2e1 * (t25 - t41 - t31) * inicial.y + 0.2e1 * (t17 + t19 - t46) * inicial.z;
	J3x4(0+i0,2+j0) = 0.4e1 * (-t25 + t51 + t52) * inicial.x + 0.2e1 * (t55 - t56 - t60) * inicial.y + 0.2e1 * (-t63 + t64 - t66) * inicial.z;
	J3x4(0+i0,3+j0) = 0.4e1 * (t70 - t19 + t72) * inicial.x + 0.2e1 * (-t66 + t63 - t75) * inicial.y + 0.2e1 * (t60 + t55 - t78) * inicial.z;
	J3x4(1+i0,0+j0) = 0.2e1 * (-t19 + t22 - t17) * inicial.x + 0.4e1 * (t85 + t11) * inicial.y + 0.2e1 * (-t60 + t55 - t90) * inicial.z;
	J3x4(1+i0,1+j0) = 0.2e1 * (t31 + t25 - t41) * inicial.x + 0.4e1 * (-t55 + t97 + t36) * inicial.y + 0.2e1 * (-t66 + t63 - t100) * inicial.z;
	J3x4(1+i0,2+j0) = 0.2e1 * (t60 + t55 - t56) * inicial.x + 0.4e1 * (t40 + t52) * inicial.y + 0.2e1 * (t19 - t108 - t17) * inicial.z;
	J3x4(1+i0,3+j0) = 0.2e1 * (-t63 + t75 - t66) * inicial.x + 0.4e1 * (t45 - t19 + t72) * inicial.y + 0.2e1 * (t25 - t116 - t31) * inicial.z;
	J3x4(2+i0,0+j0) = 0.2e1 * (-t31 + t25 - t28) * inicial.x + 0.2e1 * (-t55 + t90 - t60) * inicial.y + 0.4e1 * (t85 + t9) * inicial.z;
	J3x4(2+i0,1+j0) = 0.2e1 * (t19 - t46 - t17) * inicial.x + 0.2e1 * (-t63 + t100 - t66) * inicial.y + 0.4e1 * (-t55 + t97 + t35) * inicial.z;
	J3x4(2+i0,2+j0) = 0.2e1 * (-t66 + t63 - t64) * inicial.x + 0.2e1 * (t17 + t19 - t108) * inicial.y + 0.4e1 * (t40 - t25 + t51) * inicial.z;
	J3x4(2+i0,3+j0) = 0.2e1 * (t55 - t78 - t60) * inicial.x + 0.2e1 * (t31 + t25 - t116) * inicial.y + 0.4e1 * (t45 + t70) * inicial.z;
}

void L_Quaternion::jacob_izq_rotarVectorInv_pre(double *pre) const
{
	double t1 = c*c;
	double t2 = a*a;
	double t3 = b*b;
	double t4 = d*d;
	double t5 = t2 + t3 + t1 + t4;
	double t6 = t5 * t5;
	double t7 = 0.1e1 / t6;
	double t8 = t1 * t7;
	double t9 = t8 * a;
	double t10 = t4 * t7;
	double t11 = t10 * a;
	double t15 = t7 * a;
	double t17 = 0.2e1 * b * c * t15;
	double t18 = 0.1e1 / t5;
	double t19 = d * t18;
	double t22 = 0.2e1 * t2 * d * t7;
	double t25 = c * t18;
	double t28 = 0.2e1 * t2 * c * t7;
	double t29 = b * d;
	double t31 = 0.2e1 * t29 * t15;
	double t35 = t8 * b;
	double t36 = t10 * b;
	double t40 = t3 * c * t7;
	double t41 = 0.2e1 * t40;
	double t45 = t3 * d * t7;
	double t46 = 0.2e1 * t45;
	double t51 = t1 * c * t7;
	double t52 = t10 * c;
	double t55 = b * t18;
	double t56 = 0.2e1 * t35;
	double t58 = t7 * c;
	double t60 = 0.2e1 * a * d * t58;
	double t63 = a * t18;
	double t64 = 0.2e1 * t9;
	double t66 = 0.2e1 * t29 * t58;
	double t70 = t8 * d;
	double t72 = t4 * d * t7;
	double t75 = 0.2e1 * t11;
	double t78 = 0.2e1 * t36;
	double t85 = t3 * t7 * a;
	double t90 = 0.2e1 * t2 * b * t7;
	double t97 = t3 * b * t7;
	double t100 = 0.2e1 * t85;
	double t108 = 0.2e1 * t70;
	double t116 = 0.2e1 * t52;

	pre[0] = 0.4e1 * (t9 + t11);
	pre[1] = 0.2e1 * (-t17 + t19 - t22);
	pre[2] = 0.2e1 * (-t25 + t28 - t31);
	pre[3] = 0.4e1 * (t35 + t36);
	pre[4] = 0.2e1 * (t25 - t41 - t31);
	pre[5] = 0.2e1 * (t17 + t19 - t46);
	pre[6] = 0.4e1 * (-t25 + t51 + t52);
	pre[7] = 0.2e1 * (t55 - t56 - t60);
	pre[8] = 0.2e1 * (-t63 + t64 - t66);
	pre[9] = 0.4e1 * (t70 - t19 + t72);
	pre[10] = 0.2e1 * (-t66 + t63 - t75);
	pre[11] = 0.2e1 * (t60 + t55 - t78);
	pre[12] = 0.2e1 * (-t19 + t22 - t17);
	pre[13] = 0.4e1 * (t85 + t11);
	pre[14] = 0.2e1 * (-t60 + t55 - t90);
	pre[15] = 0.2e1 * (t31 + t25 - t41);
	pre[16] = 0.4e1 * (-t55 + t97 + t36);
	pre[17] = 0.2e1 * (-t66 + t63 - t100);
	pre[18] =  0.2e1 * (t60 + t55 - t56);
	pre[19] = 0.4e1 * (t40 + t52);
	pre[20] = 0.2e1 * (t19 - t108 - t17);
	pre[21] = 0.2e1 * (-t63 + t75 - t66);
	pre[22] = 0.4e1 * (t45 - t19 + t72);
	pre[23] = 0.2e1 * (t25 - t116 - t31);
	pre[24] = 0.2e1 * (-t31 + t25 - t28);
	pre[25] = 0.2e1 * (-t55 + t90 - t60);
	pre[26] = 0.4e1 * (t85 + t9);
	pre[27] =  0.2e1 * (t19 - t46 - t17);
	pre[28] = 0.2e1 * (-t63 + t100 - t66);
	pre[29] = 0.4e1 * (-t55 + t97 + t35);
	pre[30] = 0.2e1 * (-t66 + t63 - t64); // Corregido, decia 20
	pre[31] = 0.2e1 * (t17 + t19 - t108);
	pre[32] = 0.4e1 * (t40 - t25 + t51);
	pre[33] = 0.2e1 * (t55 - t78 - t60);
	pre[34] = 0.2e1 * (t31 + t25 - t116);
	pre[35] = 0.4e1 * (t45 + t70);
}


void L_Pose3D_cuat::err_mf_vect(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_nube_puntos *nube = (const L_Pose3D_cuat_nube_puntos *)pose3D_cuat_nube;
	L_Pose3D_cuat p;
	const L_CoordsCart3D *m, *f;
	L_CoordsCart3D pun, despl;
	int i;
	errores.reallocate((int)( nube->fijo->size()*3 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);
	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(3*i+0,0) = 1.0e60; // ERROR: decia [3*i+0,1)
			errores(3*i+1,0) = 1.0e60;
			errores(3*i+2,0) = 1.0e60;
		}
	}

	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		m = &(*nube->movido)[i];
		f = &(*nube->fijo)[i];
		pun = p.movido_a_fijo(*m); // Aca da lo mismo cual de las 2 formas usar, esta es la mas lineal
		despl = pun - (*f);
		errores(3*i+0,0) = despl.x; // ERROR: decia [3*i+0,1)
		errores(3*i+1,0) = despl.y;
		errores(3*i+2,0) = despl.z;
	}
}

void L_Pose3D_cuat::err_mf_vect_rapido(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_nube_puntos *nube = (const L_Pose3D_cuat_nube_puntos *)pose3D_cuat_nube;
	L_Pose3D_cuat p;
	L_CoordsCart3D pun;
	L_HomogeneousMatrix H;
	int i;
	errores.reallocate((int)( nube->fijo->size()*3 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);
	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(3*i+0,0) = 1.0e60; // ERROR: decia [3*i+0,1)
			errores(3*i+1,0) = 1.0e60;
			errores(3*i+2,0) = 1.0e60;
		}
	}
	H.fijaPose3D_cuat(p);

	// Este ciclo debe ser mas rapido
	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		const L_CoordsCart3D &m = nube->movido->operator[](i);
		const L_CoordsCart3D &f = nube->fijo->operator[](i);
		
		L_MatrizHomogenea_OP_mult_P(pun,H,m);
		errores(3*i+0,0) = pun.x - f.x; // ERROR: decia [3*i+0,1)
		errores(3*i+1,0) = pun.y - f.y;
		errores(3*i+2,0) = pun.z - f.z;
	}
}

void L_Pose3D_cuat::err_mf_vect_rapido_sqrt(const void *pose3D_cuat_nube, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_nube_puntos *nube = (const L_Pose3D_cuat_nube_puntos *)pose3D_cuat_nube;
	L_Pose3D_cuat p;
	L_CoordsCart3D pun;
	L_HomogeneousMatrix H;
	int i;
	errores.reallocate((int)( nube->fijo->size()*3 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);
	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(3*i+0,0) = 1.0e60; // ERROR: decia [3*i+0,1)
			errores(3*i+1,0) = 1.0e60;
			errores(3*i+2,0) = 1.0e60;
		}
	}
	H.fijaPose3D_cuat(p);

	// Este ciclo debe ser mas rapido
	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		const L_CoordsCart3D &m = nube->movido->operator[](i);
		const L_CoordsCart3D &f = nube->fijo->operator[](i);
		
		L_MatrizHomogenea_OP_mult_P(pun,H,m);
		errores(3*i+0,0) = pun.x - f.x >= 0 ? sqrt(pun.x - f.x) : -sqrt(pun.x - f.x);
		errores(3*i+1,0) = pun.y - f.y >= 0 ? sqrt(pun.y - f.y) : -sqrt(pun.y - f.y);
		errores(3*i+2,0) = pun.z - f.z >= 0 ? sqrt(pun.z - f.z) : -sqrt(pun.z - f.z);
	}
}

void L_Pose3D_cuat::err_proy_fm_vect(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_rayos_rmpf *nube = (const L_Pose3D_cuat_rayos_rmpf *)pose3D_cuat_nube_ray;
	L_Pose3D_cuat p;
	const L_RayoCamara *m;
	L_RayoCamara mov;
	const L_CoordsCart3D *f;
	L_CoordsCart3D pR, pC;
	int i;
	errores.reallocate((int)( nube->fijo->size()*2 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);
	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(2*i+0,0) = 1.0e60;
			errores(2*i+1,0) = 1.0e60;
		}
	}

	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		m = &nube->movido->operator[](i);  // Rayo
		f = &nube->fijo->operator[](i);  // Punto

		// nube->pCam es un arreglo de punto-cuaterniones correspondientes a la pose relativa de las sub-camaras
		pR = p.fijo_a_movido(*f); // Del sistema fijo al del robot
		if (nube->pCam != NULL)
			pC = (*nube->pCam)[i]->fijo_a_movido(pR); // Del sistema del robot al de la camara puesta en el robot
		else
			pC = pR;

		mov.define(pC.y/pC.x, pC.z/pC.x);
		errores(2*i+0,0) = mov.tanIzq - m->tanIzq;
		errores(2*i+1,0) = mov.tanArr - m->tanArr;
	}
}


void L_Pose3D_cuat::err_proy_mf_vect(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_rayos_pmrf *nube = (const L_Pose3D_cuat_rayos_pmrf *)pose3D_cuat_nube_ray;
	L_Pose3D_cuat p;
	const L_CoordsCart3D *m;
	const L_RayoCamara *f;
	L_RayoCamara fijo;
	L_CoordsCart3D pR, pC;
	int i;
	errores.reallocate((int)( nube->fijo->size()*2 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);
	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(2*i+0,0) = 1.0e60;
			errores(2*i+1,0) = 1.0e60;
		}
	}

	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		m = &nube->movido->operator[](i);  // Punto
		f = &nube->fijo->operator[](i);  // Rayo

		pR = p.movido_a_fijo(*m); // Del sistema movil (del landmark) al fijo (camara1)

		// Ver si el punto #i esta siendo visto desde una camara secundaria cam2
		if (nube->pCam21->size() != 0 && (*nube->pCam21)[i] != NULL)
			pC = (*nube->pCam21)[i]->fijo_a_movido(pR); // Del sistema de la camara1 al de la camara2
		else
			pC = pR;

		fijo.define(pC.y/pC.x, pC.z/pC.x);
		errores(2*i+0,0) = fijo.tanIzq - f->tanIzq;
		errores(2*i+1,0) = fijo.tanArr - f->tanArr;
	}
}

void L_Pose3D_cuat::err_proy_mf_vect_rapido(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_rayos_pmrf *nube = (const L_Pose3D_cuat_rayos_pmrf *)pose3D_cuat_nube_ray;
	L_Pose3D_cuat p;
	L_CoordsCart3D pR, pC;
	L_HomogeneousMatrix H;
	int i;
	errores.reallocate((int)( nube->fijo->size()*2 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);

	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(2*i+0,0) = 1.0e60;
			errores(2*i+1,0) = 1.0e60;
		}
	}
	
	if (nube->pCam21 == NULL) // Caso con una sola camara
	{
		H.fijaPose3D_cuat(p);
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			const L_CoordsCart3D &m = nube->movido->operator[](i);  // Punto
			const L_RayoCamara &f = nube->fijo->operator[](i);  // Rayo

			L_MatrizHomogenea_OP_mult_P(pC, H, m);
			errores(2*i+0,0) = pC.y/pC.x - f.tanIzq;
			errores(2*i+1,0) = pC.z/pC.x - f.tanArr;
		}
	}
	else // dos camaras
	{

		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			const L_CoordsCart3D &m = nube->movido->operator[](i);  // Punto
			const L_RayoCamara &f = nube->fijo->operator[](i);  // Rayo

			pR = p.movido_a_fijo(m); // Del sistema movil (del landmark) al fijo (camara1)

			// Ver si el punto #i esta siendo visto desde una camara secundaria cam2
			pC = (*nube->pCam21)[i]->fijo_a_movido(pR); // Del sistema de la camara1 al de la camara2

			errores(2*i+0,0) = pC.y/pC.x - f.tanIzq;
			errores(2*i+1,0) = pC.z/pC.x - f.tanArr;
		}
	}
}

void L_Pose3D_cuat::err_proy_mf_vect_rapido_sqrt(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_rayos_pmrf *nube = (const L_Pose3D_cuat_rayos_pmrf *)pose3D_cuat_nube_ray;
	L_Pose3D_cuat p;
	L_RayoCamara fijo;
	L_CoordsCart3D pR, pC;
	L_HomogeneousMatrix H;
	int i;
	errores.reallocate((int)( nube->fijo->size()*2 ), 1);
	p.pos.x = pcuat(0,0); p.pos.y = pcuat(1,0); p.pos.z = pcuat(2,0);
	p.ori.a = pcuat(3,0); p.ori.b = pcuat(4,0); p.ori.c = pcuat(5,0); p.ori.d = pcuat(6,0);

	if (p.ori.normaliza_ret() == false)
	{
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			errores(2*i+0,0) = 1.0e60;
			errores(2*i+1,0) = 1.0e60;
		}
	}
	
	if (nube->pCam21 == NULL) // Caso con una sola camara
	{
		H.fijaPose3D_cuat(p);
		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			const L_CoordsCart3D &m = nube->movido->operator[](i);  // Punto
			const L_RayoCamara &f = nube->fijo->operator[](i);  // Rayo

			L_MatrizHomogenea_OP_mult_P(pC, H, m);
			errores(2*i+0,0) = pC.y/pC.x - f.tanIzq;
			errores(2*i+1,0) = pC.z/pC.x - f.tanArr;
			if (errores(2*i+0,0)>0)
				errores(2*i+0,0) = sqrt(errores(2*i+0,0));
			else
				errores(2*i+0,0) = -sqrt(-errores(2*i+0,0));

			if (errores(2*i+1,0)>0)
				errores(2*i+1,0) = sqrt(errores(2*i+1,0));
			else
				errores(2*i+1,0) = -sqrt(-errores(2*i+1,0));

		}
	}
	else // dos camaras
	{

		for (i=0; i<(int)nube->fijo->size(); i++)
		{
			const L_CoordsCart3D &m = nube->movido->operator[](i);  // Punto
			const L_RayoCamara &f = nube->fijo->operator[](i);  // Rayo

			pR = p.movido_a_fijo(m); // Del sistema movil (del landmark) al fijo (camara1)

			// Ver si el punto #i esta siendo visto desde una camara secundaria cam2
			pC = (*nube->pCam21)[i]->fijo_a_movido(pR); // Del sistema de la camara1 al de la camara2

			fijo.define(pC.y/pC.x, pC.z/pC.x);
			errores(2*i+0,0) = fijo.tanIzq - f.tanIzq;
			errores(2*i+1,0) = fijo.tanArr - f.tanArr;
		}
	}
}



// Funciones que transforman rayos movidos en puntos fijos
bool L_Pose3D_cuat::_calcPose_movido_a_fijo_4rayos(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, int i0)
{
	L_CalculadorPoseTriangulo cpt;
	std::vector<L_CoordsCart3D> cm;
	std::vector<L_CoordsCart3D> cf;
	L_Pose3D_cuat t[4];
	L_RayoCamara ray;
	double e[4];
	int i, iMin;
	cm.resize(4);
	cf.resize(4);

	cf[0] = pf[i0+0];
	cf[1] = pf[i0+1];
	cf[2] = pf[i0+2];
	cf[3] = pf[i0+3];

	// Lados del triangulo, usar "cf" que es el que existe
	cpt.a = (cf[1]-cf[2]).pideR();
	cpt.b = (cf[0]-cf[2]).pideR();
	cpt.c = (cf[0]-cf[1]).pideR();
	cpt.q1 = rm[0+i0];
	cpt.q2 = rm[1+i0];
	cpt.q3 = rm[2+i0];
	
	if (cpt.metodoDirecto() == false) // Usando los primeros 3 puntos da 4 posibles soluciones
		return false;

	// Usar el cuarto punto para seleccionar la mejor solucion de las 4 posibles
	for (i=0; i<4; i++)
	{
		e[i] = 1e100;
		if (cpt.valida[i] == false)
			continue;
		// Calcular el triangulo movido numero i
		cm[0] = cpt.j1*cpt.s1[i];
		cm[1] = cpt.j2*cpt.s2[i];
		cm[2] = cpt.j3*cpt.s3[i];
		// calcular la pose que lleva el triangulo movido al fijo
		if (t[i].calcpose_movido_a_fijo_3puntos(cm, cf) == false)  // Lo deja con error 1e100
			continue;
		// calcular el cuarto punto movido
		cm[3] = t[i].fijo_a_movido(cf[3]);
		// errores
		if (cm[3].x == 0)  // Lo deja con error 1e100
			continue;
		ray.tanIzq = cm[3].y / cm[3].x;
		ray.tanArr = cm[3].z / cm[3].x;
		e[i] = (ray.tanIzq - rm[3+i0].tanIzq)*(ray.tanIzq - rm[3+i0].tanIzq) +
			(ray.tanArr - rm[3+i0].tanArr)*(ray.tanArr - rm[3+i0].tanArr);
	}
	iMin = 0;
	for (i=1; i<4; i++)
		if (e[i] < e[iMin])
			iMin = i;
	*this = t[iMin];

	if (pos.x > 1.0e5 || pos.x < -1.0e5 || pos.y > 1.0e5 || pos.y < -1.0e5 || pos.z > 1.0e5 || pos.z < -1.0e5)
		return false;

	return true;
}


// Funciones que transforman rayos movidos en puntos fijos

void L_Pose3D_cuat::_jacob_rm_calcPose_movido_a_fijo_Nrayos_matriz_malo(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, L_Matrix &J7x2n, int i0, int j0)
{
	// Al llamar esta funcion, la pose del punto-cuaternion ya debe estar estimada
	// Se pide transformar de rm a pf y sacar el jacobiano respecto a pf
	// Solo se puede transformar de pf a rm (al reves) lo cual es algo inoportuno
	L_CoordsCart3D f, m;
	L_RayoCamara r;
	L_Matrix A, b, AT, ATA, ATAinv, J12x3n, J7x3n, J12x2n, J7x12(7,12), J7x7i(7,7), J7x12i(7,12), ATAATb;
	L_HomogeneousMatrix H;
	L_Pose3D_cuat qfm, q;
	int i, j;
	double delta = 1.0e-8;
	A.reallocate((int)( 3*rm.size() ), 12);
	b.reallocate((int)( 3*rm.size() ), 1);
	A.setZero();

	for (i=0; i<(int)rm.size(); i++)
	{
		// pm = mov(q,pf)
		//
		// Si se asume xm constante durante la variacion de los r se puede despejar el jacobiano
		// (xm, ym, zm, 1) = H*(xf , yf, zf, 1)
		// (xm/xm, ym/xm, zm/xm) = H*(xf/xm , yf/xm, zf/xm, 1/xm)
		// (1, um, vm) = H*(xf/xm, yf/xm, zf/xm, 1/xm)
		// (1, Um, Vm) = (Xf/xm, Yf/xm, Zf/xm, 1/xm)*h
		// b = A*h
		// h = (ATA)^-1AT*(1, Um, Vm)  ->  dh/d(Um,Vm) = (ATA)^-1AT = J12x3n
		// qfm = cuat(H) = cuat(J12x3n*1uv);
		// q = inv(cuat(H)) == dinv(qfm)/dq * dcuat(H)/dH * J12x3n
				
		f = pf[i];
		m = fijo_a_movido(f);
		r = rm[i];

		A(3*i+0,0) = f.x/m.x;  A(3*i+0,1) = f.y/m.x; A(3*i+0,2) = f.z/m.x;  A(3*i+0,3) = 1.0/m.x; 
		A(3*i+1,4) = f.x/m.x;  A(3*i+1,5) = f.y/m.x; A(3*i+1,6) = f.z/m.x;  A(3*i+1,7) = 1.0/m.x; 
		A(3*i+2,8) = f.x/m.x;  A(3*i+2,9) = f.y/m.x; A(3*i+2,10) = f.z/m.x;  A(3*i+2,11) = 1.0/m.x; 
		b(3*i+0,0) = 1;
		b(3*i+1,0) = r.tanIzq;
		b(3*i+2,0) = r.tanArr;
	}
	AT.transpOf(A); // A:(3*n)x12, AT:12x(3*n)
	ATA.OP_mult(AT, A); // ATA: 12x12
	ATAinv.inverseOf(ATA); // 12x12

	J12x3n.OP_mult(ATAinv,AT); //12x3n -> dh/d(1,u,v)

	// Extraer las componentes que se refieren al rayo [0,1,2) -> [1,2)
	J12x2n.reallocate(12, (int)( 2*rm.size() ));
	for (i=0; i<12; i++)
		for (j=0; j<2*(int)rm.size(); j++)
			J12x2n(i,j) = J12x3n(i,3*(j/2) + (j%2)+1);

	// movido = H*fijo
	// -> H esta al reves

	// Ahora a calcular el jacobiano 7x12
	ATAATb.OP_mult(J12x3n, b);
	for (i=0; i<12; i++)
		H(i/4,i%4) = ATAATb(i,0);
	H.normalizaSVD();
	qfm = H.jacob_calcPose3D_cuat(J7x12);
	q.jacob_inverseOf(qfm, J7x7i);
	J7x12i.OP_mult(J7x7i, J7x12);
	if (qfm.ori.OP_mult_elementwise(q.ori) * q.ori.OP_mult_elementwise(ori) < 0) // La orientacion quedo al reves, cambiar el signo del jacobiano
	{
		for (i=3; i<7; i++)
			for (j=0; j<12; j++)
				J7x12i(i,j) = -J7x12i(i,j);
	}

	J7x12.checkRangeValues(-1e10, 1e10, true, "J7x12");
	J12x2n.checkRangeValues(-1e10, 1e10, true, "J7x2n");
	J7x2n.reallocate(7, (int)( 2*rm.size() ));
	J7x2n.OP_mult_sub(i0, j0, J7x12i, J12x2n, 0, 0, J7x12.li, J7x12.lj, 0, 0, J12x2n.li, J12x2n.lj);
}






void L_CamaraDistorsionRadial::distorsionaPixel_k1_Tsai(double xPix, double yPix, double &xPixDist, double &yPixDist) const
{
	//(tanIzqD,tanArrD)*( 1 + k1*rD^2 + k2*rD^4 + ...) = (tanIzq,tanArr)

#define CBRT(x) ( ((x)>=0) ? pow((x), 1.0 / 3.0) : -pow(-(x), 1.0 / 3.0) )

    double    rTan,
              rTanDist,
              factRad,
              c,
              d,
              Q,
              R,
              D,
              S,
              T,
              sinT,
              cosT;

	L_RayoCamara rayo;

	if (k1==0)
	{
		xPixDist=xPix;
		yPixDist=yPix;
		return;
	}

	rayo.tanIzq=(cenX-xPix)/dfX;
	rayo.tanArr=(cenY-yPix)/dfY;

    rTan = sqrt(rayo.tanIzq*rayo.tanIzq+rayo.tanArr*rayo.tanArr);

    c = 1 / k1;
    d = -c * rTan;

    Q = c / 3;
    R = -d / 2;
    D = Q*Q*Q + R*R;

    if (D >= 0)
	{		/* one real root */
		D = sqrt (D);
		S = CBRT (R + D);
		T = CBRT (R - D);
		rTanDist = S + T;

		if (rTanDist < 0)
			rTanDist = sqrt (-1 / (3 * k1)); // Distorsion excesiva, "rebote" de la imagen
	}
	else
	{			/* three real roots */
		D = sqrt (-D);
		S = CBRT (sqrt(R*R + D*D));
		T = atan2 (D, R) / 3;
		sinT=sin(T);
		cosT=cos(T);

		/* the larger positive root is    2*S*cos(T)                   */
		/* the smaller positive root is   -S*cos(T) + SQRT(3)*S*sin(T) */
		/* the negative root is           -S*cos(T) - SQRT(3)*S*sin(T) */

		rTanDist = -S * cosT + SQRT3 * S * sinT;	/* use the smaller positive root */
    }

    factRad = rTanDist / rTan;

	xPixDist = ( -dfX*(factRad*rayo.tanIzq) + cenX );
	yPixDist = ( -dfY*(factRad*rayo.tanArr) + cenY );
#undef CBRT
}


void L_CamaraDistorsionRadial::distorsionaPixel_k1k2k3k4_Heikkil(double xPix, double yPix, double &xPixDist, double &yPixDist) const
{
	double r2, du, dv, ud, vd;
	double x, y;
	// Pixel ideal a rayo
	x = -(xPix-cenX)/dfX;
	y = -(xPix-cenY)/dfY;
	// Rayo a pixel distorsionado
	r2 = x*x + y*y;
	du = 2*k3*x*y + k4*(r2+2*x*x);
	dv = k3*(r2+2*y*y) + 2*k4*x*y;
	ud = (1 + k1*r2 + k2*r2*r2) * x + du;
	vd = (1 + k1*r2 + k2*r2*r2) * y + dv;
	xPixDist = -dfX * ud + cenX;
	yPixDist = -dfY * vd + cenY;
}

// Para estimar el rayo usando Levenberg-Marquardt (poco recomendable)
void L_CamaraDistorsionRadial::error_Heikkil(const void *modelo, const L_Matrix &x, L_Matrix &e)
{
	e.reallocate(2,1);
	double r2, du, dv, u, v, ud, vd, uim, vim;
	double k1, k2, k3, k4;
	const L_CamaraDistorsionRadial *cam = (const L_CamaraDistorsionRadial *)modelo;
	k1 = cam->k1;
	k2 = cam->k2;
	k3 = cam->k3;
	k4 = cam->k4;
	u = x(0,0); // El rayo en el espacio (incognita)
	v = x(1,0); // El rayo en el espacio (incognita)
	r2 = u*u+v*v;
	du = 2*k3*u*v + k4*(r2+2*u*u);
	dv = k3*(r2+2*v*v) + 2*k4*u*v;
	ud = (1 + k1*r2 + k2*r2*r2) * u + du;
	vd = (1 + k1*r2 + k2*r2*r2) * v + dv;
	uim = -cam->dfX * ud + cam->cenX;
	vim = -cam->dfY * vd + cam->cenY;
	e(0,0) = cam->xim - uim;
	e(1,0) = cam->yim - vim;
}

#define L_BBFPRIORNODO3D_DISTINF (1.0e9)  // Cantidad que puede considerarse infinita respecto a las distancias entre descriptores
class L_BBFPriorNodo3D
{
public:
	L_KdNodo3D *dir;
	double dist;
	L_BBFPriorNodo3D *sig;
	L_BBFPriorNodo3D() {dir=NULL; dist=L_BBFPRIORNODO3D_DISTINF; sig=NULL;}
	//! No debe destruir dir
	void destroyRec() {L_BBFPriorNodo3D *tmp; while (sig!=NULL) {tmp=sig; sig=sig->sig; delete tmp;}}
	~L_BBFPriorNodo3D() {}
};

#define new_L_BBFPriorNodo3D(ret,L_BBFPriorLista_li) {if (L_BBFPriorLista_li.memo!=NULL) {ret = L_BBFPriorLista_li.memo; L_BBFPriorLista_li.memo=L_BBFPriorLista_li.memo->sig; ret->sig=NULL;} else ret = new L_BBFPriorNodo3D();}
#define delete_L_BBFPriorNodo3D(ret,L_BBFPriorLista_li) {ret->sig = L_BBFPriorLista_li.memo; L_BBFPriorLista_li.memo = ret;}
//#define new_L_BBFPriorNodo3D(ret,L_BBFPriorLista_li) {ret = new L_BBFPriorNodo3D();}
//#define delete_L_BBFPriorNodo3D(ret,L_BBFPriorLista_li) {delete ret;}

class L_BBFPriorLista3D
{
public:
	L_BBFPriorNodo3D *root;
	L_BBFPriorNodo3D *memo; // Aca se guardan los nodos desechados para no hacer new, delete todo el rato
	double cotaDist;
	int n;
	L_BBFPriorLista3D() {root=NULL; memo = NULL; cotaDist=100000.0; n=0;} //value muy grande
	void addFrom(L_KdNodo3D &nododiv, L_KdNodo3D &nododir, const L_CoordsCart3D &orig, double dist=-1.0);
	L_KdNodo3D *popNodoMasCercano();
	~L_BBFPriorLista3D() {if (root!=NULL) {root->destroyRec(); delete root;} if (memo!=NULL) {memo->destroyRec(); delete memo;}}
};



// addFrom los nodos en primer lugar por distancia al objetivo, y en segundo lugar por tiempo de llegada
void L_BBFPriorLista3D::addFrom(L_KdNodo3D &nododiv, L_KdNodo3D &nododir, const L_CoordsCart3D &orig, double dist) // addFrom nodo en orden de mas cercania a orig
{
	L_BBFPriorNodo3D *nuevo;
	L_BBFPriorNodo3D **pptr;
	double val, origval;
	int d;

	if (dist<-0.5)
	{
		d = nododiv.d;
		val = (d==0)*nododiv.val[0] + (d==1)*nododiv.val[1] + (d==2)*nododiv.val[2];
		origval = (d==0)*orig.x + (d==1)*orig.y + (d==2)*orig.z;
		L_SETABS(dist, origval - val);
	}
	if (dist>cotaDist)
		return;

	new_L_BBFPriorNodo3D(nuevo,(*this));
	nuevo->dist=dist;
	nuevo->dir=&nododir;
	pptr=&root;
	while (*pptr!=NULL && (*pptr)->dist<=dist)
		pptr=&(*pptr)->sig;
	nuevo->sig=*pptr;
	*pptr=nuevo;
	n++;
}

L_KdNodo3D *L_BBFPriorLista3D::popNodoMasCercano() // devuelve nodo sin memoria propia
{
	L_BBFPriorNodo3D *pri;
	L_KdNodo3D *ret;

	if (root==NULL)
		return NULL;

	if (root->dist>cotaDist)
	{
		root->destroyRec();
		delete root;
		root=NULL;
		return NULL;
	}

	pri=root;
	root=root->sig;
	ret=pri->dir;
	delete_L_BBFPriorNodo3D(pri, (*this));
	n--;
	return ret;
}


///

void L_KdNodo3D::destroyRec()
{
	if (izq != NULL)
	{
		izq->destroyRec();
		delete izq;
	}
	if (der != NULL)
	{
		der->destroyRec();
		delete der;
	}
}

void L_KdNodo3D::destruyeNoRec()
{
	int si[100];
	L_KdNodo3D *sp[100];
	int n=0;

	sp[n] = this;
	si[n] = 0;
	do
	{
		switch(si[n])
		{
		case 0:
			si[n] = 1;
			if (sp[n]->izq != NULL)
				sp[n+1] = sp[n]->izq;
			else
				break;
			si[n+1] = 0;
			n++;
			break;
		case 1:
			si[n] = 2;
			if (sp[n]->izq != NULL)
				delete sp[n]->izq;
			if (sp[n]->der != NULL)
				sp[n+1] = sp[n]->der;
			else
				break;
			si[n+1] = 0;
			n++;
			break;
		case 2:
			si[n] = 3;
			if (sp[n]->der != NULL)
				delete sp[n]->der;
			n--;
		}
	} while (n>=0);
}

int L_KdNodo3D::cmpRepetidos(const void *kdNodo3D_1, const void *kdNodo3D_2)
{
	L_KdNodo3D *p1 = (L_KdNodo3D *)kdNodo3D_1;
	L_KdNodo3D *p2 = (L_KdNodo3D *)kdNodo3D_2;
	if      (p1->val[0] > p2->val[0])
		return 1;
	else if (p1->val[0] < p2->val[0])
		return -1;
	else if (p1->val[1] > p2->val[1])
		return 1;
	else if (p1->val[1] < p2->val[1])
		return -1;
	else if (p1->val[2] > p2->val[2])
		return 1;
	else if (p1->val[2] < p2->val[2])
		return -1;
	return 0;
}

#define L_KD3_dist(pun,vec) (sqrt((orig[0]-vec[0])*(orig[0]-vec[0]) + (orig[1]-vec[1])*(orig[1]-vec[1]) + (orig[2]-vec[2])*(orig[2]-vec[2])))

int L_KdNodo3D::cmpDimension(const void *kdNodo3D_1, const void *kdNodo3D_2)
{
	double v1 = ((L_KdNodo3D *)kdNodo3D_1)->val[((L_KdNodo3D *)kdNodo3D_1)->d];
	double v2 = ((L_KdNodo3D *)kdNodo3D_2)->val[((L_KdNodo3D *)kdNodo3D_2)->d];
	return (v1 > v2) - (v1 < v2);
}

// Esta funcion se deja aca para evitar que ensucie la class_name en el .h
int L_KdTree3D_calcDirMayorVarianza(L_KdNodo3D *inic, int n)
{
	int i;
	double sx[3]={0,0,0}, sxx[3]={0,0,0};
	for (i=0; i<n; i++)
	{
		sx[0] += inic[i].val[0];
		sxx[0] += inic[i].val[0] * inic[i].val[0];
		sx[1] += inic[i].val[1];
		sxx[1] += inic[i].val[1] * inic[i].val[1];
		sx[2] += inic[i].val[2];
		sxx[2] += inic[i].val[2] * inic[i].val[2];
	}
	sxx[0] = sxx[0] / n  - (sx[0]/n * sx[0]/n);
	sxx[1] = sxx[1] / n  - (sx[1]/n * sx[1]/n);
	sxx[2] = sxx[2] / n  - (sx[2]/n * sx[2]/n);
	if (sxx[0] == 0 && sxx[1] ==0 && sxx[2] == 0)
		return -1;
	else if (sxx[0] >= sxx[1] && sxx[0] >= sxx[2])
		return 0;
	else if (sxx[1] >= sxx[0] && sxx[1] >= sxx[2])
		return 1;
	return 2;
	
}

L_KdNodo3D *L_KdNodo3D::creaRec(L_KdNodo3D *arr, int n)
{
	int i, d;
	int nmedio;
	if (n==1)
	{
		L_KdNodo3D *ret = new L_KdNodo3D;
		*ret = arr[0]; // Se copia .i
		ret->d = -1; // Es hoja
		return ret;
	}
	d = L_KdTree3D_calcDirMayorVarianza(arr, n);
	if (d<0 || d>2)
		printf("d no valido\n");
	for (i=0; i<n; i++)
		arr[i].d = d;
	qsort(arr, n, sizeof(L_KdNodo3D), &L_KdNodo3D::cmpDimension);
	// Tratar de buscar un punto de divisionl la idea es que arr[nmedio-1] != arr[nmedio]
	// Primero hacia la izquierda
	nmedio = n/2;
	while (nmedio > 1 && arr[nmedio-1].val[d] == arr[nmedio].val[d])
		nmedio--;
	// Luego hacia la derecha
	if (nmedio == 0)
	{
		nmedio = n/2;
		while (nmedio < n && arr[nmedio-1].val[d] == arr[nmedio].val[d])
			nmedio++;
	}
	throw_L_ArgException_if(arr[nmedio-1].val[d] == arr[nmedio].val[d], "L_KdTree3D_creaRec() : vector repetido");


	// Dividir arr[0.__n-1] en arr[0...nmedio-1] y arr[nmedio.__n-1]
	L_KdNodo3D *ret = new L_KdNodo3D;
	*ret = arr[nmedio];
	ret->d = d; // Es interior
	ret->i = -1;
	ret->umbral = arr[nmedio-1].val[d];
	ret->izq = L_KdNodo3D::creaRec(arr, nmedio); // arr[0,nmedio-1] son <=umbral
	ret->der = L_KdNodo3D::creaRec(arr+nmedio, n-nmedio); // arr[nmedio,n-1] son > umbral
	return ret; // Decia arr, que mongo...
}

struct L_KdNodo3D_stack_data
{
	L_KdNodo3D *inic;
	L_KdNodo3D *ret;
	int nmedio;
	int n;
	int dir; //0 = izq; 1 = der
};

L_KdNodo3D *L_KdNodo3D::creaNoRec(L_KdNodo3D *arr, int nTot)
{
	int nit = 0;
	int i, d;
	int nmedio;
	L_KdNodo3D_stack_data st[100]; // Para eliminar recursividad (hasta donde se puede)
	int ist = 0;
	if (nTot==1)
	{
		L_KdNodo3D *ret = new L_KdNodo3D;
		*ret = arr[0];
		ret->d = -1; // Es hoja
		return ret;
	}
	st[ist].inic = arr;
	st[ist].n = nTot;
	st[ist].ret = NULL;
	st[ist].nmedio = 0;
	st[ist].dir = 0;

	do
	{
		nit++;
		if (st[ist].dir == 0)
		{
			if (st[ist].n == 1)
			{
				st[ist].ret = new L_KdNodo3D;
				*st[ist].ret = st[ist].inic[0];
				st[ist].ret->d = -1; // Es hoja
				st[ist].dir = 2;
				ist--;
				continue;
			}
			d = L_KdTree3D_calcDirMayorVarianza(st[ist].inic, st[ist].n);
			throw_L_ArgException_if(d==-1, "L_KdTree3D_creaNoRec() : nodos repetidos");
			for (i=0; i<st[ist].n; i++)
				st[ist].inic[i].d = d;
			qsort(st[ist].inic, st[ist].n, sizeof(L_KdNodo3D), &L_KdNodo3D::cmpDimension);
			// Tratar de buscar un punto de divisionl la idea es que inic[nmedio-1] != inic[nmedio]
			// Primero hacia la izquierda
			nmedio = st[ist].n/2;
			while (nmedio > 0 && st[ist].inic[nmedio-1].val[d] == st[ist].inic[nmedio].val[d])
				nmedio--;
			// Luego hacia la derecha
			if (nmedio == 0)
			{
				nmedio = st[ist].n/2;
				while (nmedio < st[ist].n && st[ist].inic[nmedio-1].val[d] == st[ist].inic[nmedio].val[d])
					nmedio++;
			}
			throw_L_ArgException_if(st[ist].inic[nmedio-1].val[d] == st[ist].inic[nmedio].val[d], "L_KdTree3D_creaRec() : vector repetido");
			// Dividir inic[0.__n-1] en inic[0...nmedio-1] y inic[nmedio.__n-1]
			L_KdNodo3D *ret = new L_KdNodo3D;
			*ret = st[ist].inic[nmedio];
			ret->d = d; // Es interior
			ret->i = -1;
			ret->umbral = st[ist].inic[nmedio-1].val[d];

			// Empujar lso datos al stack
			st[ist].ret = ret;
			st[ist].nmedio = nmedio;
			st[ist].dir = 1;

			st[ist+1].inic = st[ist].inic;
			st[ist+1].n = nmedio;
			st[ist+1].ret = NULL;
			st[ist+1].nmedio = 0;
			st[ist+1].dir = 0;
			ist++;
			continue;
		}
		else if (st[ist].dir == 1)
		{
			st[ist].ret->izq = st[ist+1].ret;
			st[ist].dir = 2;
			st[ist+1].inic = st[ist].inic+st[ist].nmedio;
			st[ist+1].n = st[ist].n-st[ist].nmedio;
			st[ist+1].ret = NULL;
			st[ist+1].nmedio = 0;
			st[ist+1].dir = 0;
			ist++;
			continue;
		}
		else if (st[ist].dir == 2)
		{
			st[ist].ret->der = st[ist+1].ret;
			st[ist].dir = 3;
			ist--;
		}
		throw_L_ArgException_if(st[ist].dir > 2, "L_KdTree3D_creaNoRec() : simulacion de recursion incorrecta");
	} while (ist >= 0);
	return st[0].ret; // Decia inic, que mongo...
}


bool L_KdTree3D::createFrom(const std::vector<L_CoordsCart3D> &arr)
{
	int i;
	if (root != NULL)
	{
		//root->destroyRec();
		root->destruyeNoRec();
		delete root;
		root = NULL;
	}
	L_Array<L_KdNodo3D> pun((int)arr.size());
	if (pun.size() == 0)
		return false;
	for (i=0; i<(int)pun.size(); i++)
	{
		pun[i].i = i;
		pun[i].val[0] = arr[i].x;
		pun[i].val[1] = arr[i].y;
		pun[i].val[2] = arr[i].z;
	}
	pun.resize(pun.size() - pun.sortIsolatingRepeatedAtEnd(&L_KdNodo3D::cmpRepetidos)); // Mueve los repetidos al final y luego disminuye n
	if (pun.size() == 0)
		return false;
	//root = L_KdNodo3D::creaRec(pun.elem, pun.size());
	root = L_KdNodo3D::creaNoRec(pun.data(), pun.size());
	nHojas = pun.size();
	nElim = (int)( arr.size() - pun.size() );
	return true;
}

//! Devuelve los n nodos del árbol que son más parecidos a pun en {dest[0],dest[1],...}.
void L_KdTree3D::buscaMasCercanos(L_CoordsCart3D &pun, int n, int *indices, double *dist)
{
	int nComp=0;
	L_BBFPriorLista3D lisPrio;
	L_KdNodo3D *ptr;
	int des, desTemp;
	double dif, difTemp;
	int i;
	int nBus;
	double orig[3];

	orig[0] = pun.x;
	orig[1] = pun.y;
	orig[2] = pun.z;

	if (usePercentage)
		nBus=(int)(1+floor(porcRev*nHojas));
	else
		nBus=nCompMax;

	if (nBus>nHojas)
		nBus=nHojas;

	for (i=0; i<n; i++)
	{
		indices[i]=-1;
		dist[i]=L_BBFPRIORNODO3D_DISTINF; // número grande
	}

	if (root==NULL)
		return;

	else if (root->d < 0)
	{
		dist[0]=L_KD3_dist(pun,root->val);
		indices[0]=root->i;
		return;
	}

	dif = orig[root->d] - root->umbral;

	if (dif>0)
	{
		lisPrio.addFrom(*root, *root->der, pun, 0); // ir primero por la derecha, distancia puesta a 0
		lisPrio.addFrom(*root, *root->izq, pun);
	}
	else
	{
		lisPrio.addFrom(*root, *root->izq, pun, 0); // ir primero por la izquierda, distancia puesta a 0
		lisPrio.addFrom(*root, *root->der, pun);
	}

	while (true)
	{
		ptr=lisPrio.popNodoMasCercano();
		if (ptr==NULL) // se acabó la búsqueda
			break;
		if (nComp>=nBus) // se acabó la búsqueda
			break;
		// Bajar por el árbol hasta llegar a una hoja => des
		while (true)
		{
			if (ptr->d < 0) // ptr es nodo hoja
				break;
			else // ptr es nodo de division
			{
				if (orig[ptr->d] > ptr->umbral) // bajar por la derecha
				{
					lisPrio.addFrom(*ptr, *ptr->izq, pun);
					ptr=ptr->der;
				}
				else // bajar por la izquierda
				{
					lisPrio.addFrom(*ptr, *ptr->der, pun);
					ptr=ptr->izq;
				}
			}
		}
		// Aca hay un nodo hoja ptr, insertarlo en algun lugar de la lista
		des=ptr->i;
		dif=L_KD3_dist(orig, ptr->val);
		nComp++;
		for (i=0; i<n; i++)
		{
			if (indices[i]==-1)
			{
				indices[i]=des;
				dist[i]=dif;
				break;
			}
			if (dif<dist[i])
			{
				difTemp=dif;
				desTemp=des;
				dif=dist[i];
				des=indices[i];
				dist[i]=difTemp;
				indices[i]=desTemp;
			}
		}
	}
}

int L_KdTree3D::buscaMasCercano(L_CoordsCart3D &pun)
{
	L_KdNodo3D *ptr;
	double val[3];
	val[0] = pun.x;
	val[1] = pun.y;
	val[2] = pun.z;
	ptr = root;
	while (ptr!=NULL && ptr->d >= 0) // Hasta llegar al nodo hoja
	{
		if (val[ptr->d] <= ptr->umbral)
			ptr = ptr->izq;
		else
			ptr = ptr->der;
	}
	return (ptr==NULL) ? -1 : ptr->i;
}



//! Devuelve los n nodos del árbol que son más parecidos a pun en {dest[0],dest[1],...}.
void L_ArregloP3D::buscaMasCercanos(L_CoordsCart3D &pun, int n, int *indices, double *dist)
{
	int i, j, k;
	double d;

	for (i=0; i<n; i++)
	{
		indices[i] = 0;
		dist[i] = 1e60*(1+i);
	}

	for (i=0; i<(int)arr.size(); i++)
	{
		d = arr[i].distanciaA(pun);
		for (j=0; j<n; j++)
		{
			if (d < dist[i])
			{
				for (k=n-1; k>j; k--)
				{
					indices[k] = indices[k-1];
					dist[k] = dist[k-1];
				}
				dist[j] = d;
				indices[j] = i;
				break;
			}
		}
	}
}

int L_ArregloP3D::buscaMasCercano(L_CoordsCart3D &pun)
{
	int i, imin=-1;
	double d, d2=1e60;
	for (i=0; i<(int)arr.size(); i++)
	{
		d = arr[i].distancia2A(pun);
		if (d < d2)
		{
			d2 = d;
			imin = i;
		}
	}
	return imin;
}


void L_CurvaPuntos::suavizar(L_CurvaPuntos &resultado, double sigma)
{
	int i, j;
	L_CoordsCart3D sum;
	double anchoGaus2, dist2, den, factor, exponente, valgaus;

	resultado.resize(size());

	anchoGaus2 = (3*sigma)*(3*sigma);

	factor = -0.5/(sigma*sigma);

	for (i=0; i<(int)size(); i++)
	{
		sum.x = 0;
		sum.y = 0;
		sum.z = 0;
		den = 0;

		// sum de los valores dentro de la gaussian hacia la mitad positiva
		for (j=0; j<(int)size(); j++)
		{
			dist2 = L_CoordsCart3D_distancia2A(operator[](i), operator[](j));
			//if (dist2 > anchoGaus2)
			//	break;
			exponente = factor * dist2;
			valgaus = exp(exponente);
			sum.x += valgaus * operator[](j).x;
			sum.y += valgaus * operator[](j).y;
			sum.z += valgaus * operator[](j).z;
			den += valgaus;
		}

		resultado[i].x = sum.x/den;
		resultado[i].y = sum.y/den;
		resultado[i].z = sum.z/den;
	}
}

void L_CurvaPuntos::calcularNormales2D(L_CurvaPuntos &normales, L_KdTree3D &arb)
{
	L_Statistics2D estad;
	int i, j, inds[10];
	double dists[10], ang, c, s, cPr=0, sPr=0;
	normales.resize(size());
	for (i=0; i<(int)size(); i++)
	{
		arb.buscaMasCercanos(operator[](i), 10, inds, dists);
		for (j=0; j<10 && inds[j] >= 0; j++)
			L_Statistics2D_push(estad,operator[](inds[j]).x, operator[](inds[j]).y); // Decia i...
		ang = estad.getLineOrientationDeming(1);
		ang += M_PI/2; // Angulo perpendicular a los puntos
		L_Statistics2D_clear(estad);

		c = cos(ang);
		s = sin(ang);

		if (c*cPr + s*sPr < 0) // Para que no cambie de "signo" la normal
		{
			c = -c;
			s = -s;
		}
		cPr = 0.8*cPr + 0.2*c;
		sPr = 0.8*sPr + 0.2*s;

		normales[i].x = c;
		normales[i].y = s;
		normales[i].z = 0;
	}
}

#ifdef __COMPAT_OPENNI__
bool L_CapturadorImagenKinectNI::capturarImagen()
{
	XnStatus nRetVal = XN_STATUS_OK;

	if (!activo && crear() == false)
		return false;

	// Wait for new data to be available
	//nRetVal = context.WaitOneUpdateAll(depth);
	nRetVal = context.WaitAndUpdateAll();
	if (nRetVal != XN_STATUS_OK)
	{
		printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		return false;
	}

    // Take current depth map
	const XnUInt8 *pImageMap = imageGen.GetImageMap();

	if (pImageMap == NULL)
		return false;

	XnMapOutputMode imMapMode, depthMapMode;
	nRetVal = imageGen.GetMapOutputMode(imMapMode);
	if (nRetVal != XN_STATUS_OK)
	{
		printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		return false;
	}

	im.reallocate(imMapMode.nXRes, imMapMode.nYRes);
	im.copyFromBuffer(pImageMap);

	// Resize image if needed
	while(redimensionar.im_activo && (im.lx > redimensionar.im_lx || im.ly > redimensionar.im_ly))
		im.halfResolution();

	if (prof.capturarProf == false && prof.calcularProf == true)
		printf("Error in L_CapturadorImagenKinectNI: prof.calcularProf == false && prof.calcularProf == true\n");
	if (prof.calcularProf == false && prof.calcularProfRGB == true)
		printf("Error in L_CapturadorImagenKinectNI: prof.calcularProf == false && prof.calcularProfRGB == true\n");

	if (prof.capturarProf)
	{
		throw_L_ArgException_if(sizeof(XnDepthPixel)!=sizeof(L_uint16), "L_CapturadorImagenKinectNI::capturarImagen() : profundidad no tiene 16 bits");
		const XnDepthPixel* pDepthMap = depthGen.GetDepthMap();

		if (pDepthMap == NULL)
			return false;
		nRetVal = depthGen.GetMapOutputMode(depthMapMode);
		if (nRetVal != XN_STATUS_OK)
		{
			printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
			return false;
		}
		// Save depth image (unsigned short)
		imKin.reallocate_fc(depthMapMode.nYRes, depthMapMode.nXRes);
		imKin.copyFromBuffer(pDepthMap);
		
		// Resize image if needed
		while(redimensionar.imProf_activo && (imKin.lx > redimensionar.imProf_lx || imKin.ly > redimensionar.imProf_ly))
			imKin.halfResolution_noResample();

		if (prof.calcularProf)
		{
			imDepth.reallocate(imKin.lx, imKin.ly);
			for (int i=0; i<imKin.lx*imKin.ly; i++)
				imDepth.data()[i] = imKin.data()[i]; // It is yet in [mm]
			/*
			// Calcular imagen de profundidad (double)
			imDepth.reallocate(depthMapMode.nXRes, depthMapMode.nYRes);
			double *ip = &imDepth.pix(0,0);
			const XnDepthPixel *pd, *pdend = pDepthMap + imDepth.lx*imDepth.ly;

			// Data is yet in [mm]
			for (pd = pDepthMap; pd < pdend; ip++, pd++)
				*ip = *pd;
				*/
		}

		if (prof.calcularProf && prof.calcularProfRGB)
		{
			L_RayoCamara r;
			int ivis, jvis;
			L_CamaraPinhole camKinV, camKinD;
			L_CoordsCart3D p;
			camKinV.fijaParametrosKinectVisible(im.lj, im.li);
			camKinD.fijaParametrosKinectIR(imDepth.lj, imDepth.li);
			imDepthRGB.reallocate(imDepth.lj, imDepth.li);
			memset(imDepthRGB.data(), 0, imDepthRGB.ljStep*imDepthRGB.li);
			for (int i=0; i<imDepth.li; i++)
			{
				// We suppose parallel cameras => ivis don't depend on j
				double depth = 4000;
				camKinD.calcRayo(r, 0, i);
				camKinV.calcPixel(r, jvis, ivis);
				if (ivis < 1 || ivis > im.lj)
					continue;

				for (int j=0; j<imDepth.lj; j++)
				{
					depth = imDepth(i,j);
					if (depth == 0)
						continue;

					// Camera are assumed to be parallel; in other case, uncomment the code
					r.tanIzq=(camKinD.cenX-j)/camKinD.dfX; //camKinD.calcRayo(r, j, i);

					//p.x = depth; // depth in [mm]
					p.y = r.tanIzq * depth;
					//p.z = r.tanArr * depth;

					// Kinect component order : [  proyector   RGB IR   "XBOX360"  ]
					p.y += 25; // kinect baseline in [mm]

					r.tanIzq = p.y/depth;
					//r.tanArr = p.z/depth;
					
					jvis = (int)(camKinV.cenX-r.tanIzq*camKinV.dfX+0.5); // camKinV.calcPixel(r, jvis, ivis);

					if (jvis < 1 || jvis >= im.lj)
						continue;

					if (imDepthRGB(ivis-1, jvis-1) == 0)
						imDepthRGB(ivis-1, jvis-1) = depth;
					if (imDepthRGB(ivis-1, jvis) == 0)
						imDepthRGB(ivis-1, jvis) = depth;
					if (imDepthRGB(ivis, jvis-1) == 0)
						imDepthRGB(ivis, jvis-1) = depth;
					if (imDepthRGB(ivis, jvis) == 0)
						imDepthRGB(ivis, jvis) = depth;
				}
			}
		}
	}

	return true;
}
#endif // __COMPAT_OPENNI__

long L_calcNumVerif(const char *str, int len0)
{
	long res=0;
	int i;
	long len=(long)strlen(str);
	if (len!=len0)
		return 127;
	for (i=0; i<len; i++)
		res=(long)(str[i] + str[i]*len + len*len);
	return res;
}

#ifdef VERIF_INT
char strIntegr_2[]="\x2d\x20\x28\x63\x29\x20\x50\x61\x74\x72\x69\x63\x69\x6f\x20\x4c\x6f\x6e\x63\x6f\x6d\x69\x6c\x6c\x61\x20\x2d\x20\x32\x30\x30\x36\x20\x2d";

bool L_VerificaIntegridadCodigo()
{
	if (L_calcNumVerif(strIntegr_1, 17)!=1099)
	{
		printf("Error: codeMapping en memoria corrupto\n");
		return false;
	}
	if (L_calcNumVerif(strIntegr_2, 34)!=2731)
	{
		printf("Error: codeMapping en memoria corrupto\n");
		return false;
	}
	return true;
}

bool L_VerificaIntegridadCodigoRand(double prob)
{
	if (rand()%10000L < prob*10000L)
		return L_VerificaIntegridadCodigo();
	return true;
}

#endif

void L_ImprimeDatosIntegridad(const char *ch)
{
	printf("%ld %ld\n", long(strlen(ch)), L_calcNumVerif(ch, (int)strlen(ch)));
}














































