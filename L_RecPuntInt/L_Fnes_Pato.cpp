#include "L_Fnes_Pato.h"
#include <cfloat>
//#include "fnes_difgaus.h"

// Incluye manipulacion de: numeros complejos, fracciones, listas, arboles, compresion Huffman,
// arreglos, strings, parametros (servidor de parametros), imagenes, matrices y modelos de camaras
// Incluye interfaces a ATL, OpenCV, SVS si se usan "defines" adecuados
// Nota: Las imagenes son representadas como [columna][fila], es decir, se almacenan como columnas consecutivas
// Compatibilidad: Ansi C++
// Ojo: Poner acentos en los comentarios es anti ANSI
// Autor: Patricio Loncomilla (ploncomi@gmail.com)


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


#ifdef __COMPAT_FLTK__
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Menu_Bar.H>
#endif // __COMPAT_FLTK__

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#ifdef L_DEFINE_PRIME_NUMBERS
size_t L_Primos[100] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29 // Para que tener el arreglo si los primos salen aca...
, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71 
, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113 
, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173 
, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229 
, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281 
, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349 
, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409 
, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463 
, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541};
#endif



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
	int res = lu_factorize(A, pm);
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



















int throw_L_ArgException_fn(const char *s)
{
	printf("L_ArgException en %s\n", s);
	throw L_ArgException();
	return -1;
}

int printing_L_ArgException_fn(const char *s)
{
	printf("\n");
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf("         L_ArgException en %s\n", s);
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf("\n");
	return -1;
}

int throw_L_ArgException_fn_NULL(const char *s)
{
	char *puntero = NULL;
	printf("L_ArgException en %s\n", s);
	*puntero = 3; // Aca se deberia caer limpiamente
	return -1;
}

#ifdef L_USE_SMOOTH_SIGMA_STORING_TREE
const L_ImagenBNxfloatFnConvArreglo L_imBNxFnArr(false); // El unico object de tipo L_ImagenBNxfloatFnConvArreglo
#endif //L_USE_SMOOTH_SIGMA_STORING_TREE

bool L_VerificaIntegridadCodigo();
bool L_VerificaIntegridadCodigoRand(double prob);

void L_hard_shutdown_fn(const char *file, long line, const char *x, ...) //! Función para lanzar excepción grave en el sistema
{
	char *u=NULL;
	va_list ppp;
	va_start(ppp, x);
	vprintf(x, ppp);
	*u = 5;
	//throw L_HardErrorException(file, line, x, ppp); // ppp se cierra acá dentro
}

bool L_printError(const char *formato, ...) //! Función para lanzar excepción grave en el sistema
{
	va_list ppp;
	va_start(ppp, formato);
	vprintf(formato, ppp);
	va_end(ppp);
	return true;
}


void L_pruebaRand()
{
	int i;
	L_Rand randF;
	randF.init(-5);
	printf("Probando L_RANDOM()\n");
	for (i=0; i<30; i++)
		printf("%f entre 3 y 6\n", L_RANDOM(3,6));
	printf("\n");
	for (i=0; i<30; i++)
		printf("%f entre -3 y -1\n", L_RANDOM(-3,-1));
	printf("\n");
	printf("Probando L_Rand::random()\n");
	for (i=0; i<30; i++)
		printf("%f entre 3 y 6\n", randF.random(3,6));
	printf("\n");
	for (i=0; i<30; i++)
		printf("%f entre -3 y -1\n", randF.random(-3,-1));
	printf("Presione ENTER para salir\n");
	getchar();

	L_Rand azar, azar2;

	azar.init(-5);

	for (i=0; i<10; i++)
		printf("%f\n", azar.random(-10, 10));

	azar2.iterate_until_state(azar.seed0, azar.numCalls, azar.numCalls2);
	printf("\nLas siguientes dos secuencias deben ser iguales\n");

	for (i=0; i<10; i++)
		printf("%f\n", azar.random(-10, 10));

	printf("\n");

	for (i=0; i<10; i++)
		printf("%f\n", azar2.random(-10, 10));

	printf("Presione ENTER para salir\n");
	getchar();
}


L_FunctionCallList llamFnListaPato; // Dificil que el name se repita por mala suerte (si ocurre esto, nada funciona)
L_FunctionCall llamFnNodoPato;
int llamFnMaxNivel = 2;

void L_DebugPushEnterFuncion(const char *s)
{
	int i;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
	{
		L_hard_shutdown2("Error de integridad al ingresar a la funcion \"%s\"\n", s);
	}
#endif
	for (i=0; s[i]!=0 && i<L_FUNCTION_CALL_SIZE-1; i++)
		llamFnNodoPato.name[i]=s[i];
	llamFnNodoPato.name[i]=0;

	llamFnListaPato.push_front(llamFnNodoPato);
	if ((int)llamFnListaPato.size() <= llamFnMaxNivel)
	{
		printf("*** ");
		for (i=0; i<(int)llamFnListaPato.size(); i++)
			printf(".");
		printf("%s( )\n",s);
	}
}

void L_DebugPopExitFuncion(const char *s)
{
	int i;
	for (i=0; s[i]!=0 && i<L_FUNCTION_CALL_SIZE-1; i++)
		llamFnNodoPato.name[i]=s[i];
	llamFnNodoPato.name[i]=0;

	if (llamFnListaPato.size()==0)
	{
		printf("XXX Retorno desde %s( ) sin llamados anteriores registrados\n", s);
	}
	else if ( strcmp(llamFnNodoPato.name, llamFnListaPato.root->c.name) != 0 )
	{
		printf("XXX Retorno desde %s( ) y se esperaba %s( )\n", s, llamFnListaPato.root->c.name);
	}
	else if ((int)llamFnListaPato.size() <= llamFnMaxNivel)
	{
		printf("*** ");
		for (i=0; i<(int)llamFnListaPato.size(); i++)
			printf(".");
		printf("~%s\n", s);
	}
	if (llamFnListaPato.root!=NULL)
	{
		llamFnListaPato.root=llamFnListaPato.root->sig;
		if (llamFnListaPato.root==NULL)
			llamFnListaPato.pult=&llamFnListaPato.root;
		llamFnListaPato.n--;
	}
}

char *L_gets(char *buf, int sizebuf)
{
	char *buffn = buf+sizebuf-1;
	char *ret = fgets(buf, sizebuf-1, stdin);
	if (ret == NULL)
		return NULL;
	while(buf < buffn && *buf!=0 && *buf!='\n')
		buf++;
	*buf=0;
	return ret;
}

void L_verifTypeSizes()
{
#ifdef __TURBOC__
	#pragma warn -8008 // Para que el compilador no reclame que las siguientes condiciones son obvias...
	#pragma warn -8066
#endif
	if (sizeof(L_int8)!=1 || sizeof(L_int16)!=2 || sizeof(L_int32) !=4)
	{
		fprintf(stderr, "Tipos de datos de tamano incorrecto. Modifique L_int8, L_int16 y L_int32 en el codeMapping\n");
		L_hard_shutdown("Tipos de datos de tamano incorrecto. Modifique L_int8, L_int16 y L_int32 en el codeMapping");
	}
	if (sizeof(L_uint8)!=1 || sizeof(L_uint16)!=2 || sizeof(L_uint32) !=4)
	{
		fprintf(stderr, "Tipos de datos de tamano incorrecto. Modifique L_uint8, L_uint16 y L_uint32 en el codeMapping\n");
		L_hard_shutdown("Tipos de datos de tamano incorrecto. Modifique L_uint8, L_uint16 y L_uint32 en el codeMapping");
	}
#ifdef __TURBOC__
	#pragma warn +8008
	#pragma warn +8066
#endif
}

#ifdef VERIF_INT
char strIntegr_1[]="\x2d\x20\x4c\x5f\x46\x6e\x65\x73\x5f\x50\x61\x74\x6f\x2e\x68\x20\x2d";
#endif


#ifdef __FORCE_ALIGN_FLOAT_DOUBLE__
void* L_malloc_aligned(size_t bytes, size_t alin)
{
	if(alin == 0 || alin == 1)
		return malloc(bytes);

	size_t ptr = (size_t)malloc(bytes + alin + 2*sizeof(size_t));



	if(ptr == NULL)
		return NULL;

	size_t aligned_ptr = (ptr + alin + sizeof(size_t)) / alin * alin;

	*(size_t*)(aligned_ptr - sizeof(size_t)) = ptr;

	return (void*)aligned_ptr;
}
 
void L_free_aligned(void* p)
{
	size_t mem = ( *((size_t*)((size_t)p - sizeof(size_t))) );
	free((void*)mem);
}


double *L_malloc_double(int n)
{
	return (double *)L_malloc_aligned(sizeof(double)*n, sizeof(double[4]));
}

void L_free_double(double *mem)
{
	L_free_aligned(mem);
}

float *L_malloc_float(int n)
{
	return (float *)L_malloc_aligned(sizeof(float)*n, sizeof(float[8]));
}

void L_free_float(float *mem)
{
	L_free_aligned(mem);
}


template <> double** L_new2d<double>(int lx0, int ly0) //! Pide memoria bidimensional alineada para double[4] -> __m128d
{
	int i;
    double ** im=NULL;
	try{
#ifdef L_BASIC_DEBUG
		if (lx0<=0 || ly0<=0)
			L_hard_shutdown("L_New2d<double> : non-positive size allocation");
#endif
		im=new double*[lx0]; // lanza std::bad_alloc
		if (im==NULL) throw std::bad_alloc();
		try
		{
			im[0]=L_malloc_double( lx0*ly0 );
			if (im[0]==NULL) throw std::bad_alloc();
			for (i=1; i<lx0; i++)
			{
				im[i]=im[0]+ly0*i;
			}
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

template <> void L_delete2d<double>(double **im) //! destroy memoria bidimensional alineada para double[4] -> __m128d
{
	L_free_double( im[0] );
	delete[] im;
}



template <> float** L_new2d<float>(int lx0, int ly0) //! Pide memoria bidimensional alineada para float[4] -> __m128d
{
	int i;
    float ** im=NULL;
	try{
#ifdef L_BASIC_DEBUG
		if (lx0<=0 || ly0<=0)
			L_hard_shutdown("L_New2d<float> : non-positive size allocation");
#endif
		im=new float*[lx0]; // lanza std::bad_alloc
		if (im==NULL) throw std::bad_alloc();
		try
		{
			im[0]=L_malloc_float( lx0*ly0 );
			if (im[0]==NULL) throw std::bad_alloc();
			for (i=1; i<lx0; i++)
			{
				im[i]=im[0]+ly0*i;
			}
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

template <> void L_delete2d<float>(float **im) //! destroy memoria bidimensional alineada para float[4] -> __m128d
{
	L_free_float( im[0] );
	delete[] im;
}

void L_prueba_alineamiento_L_malloc(size_t blo)
{
	double *buf[100];
	int i;
	size_t delta;
	bool alineado = true;

	for (i=0; i<100; i++)
	{
		buf[i] = L_malloc_double(rand()%100);
		delta = (size_t)(buf[i]) % blo;
		if (delta != 0)
			alineado = false;
	}
	for (i=0; i<100; i++)
		L_free_double(buf[i]);
	if (alineado)
		printf("L_malloc_double() devuelve bloques alineados\n");
	else
		printf("L_malloc_double() no devuelve bloques alineados (blo=%d delta = %d)\n", blo, delta);
}

#endif // __FORCE_ALIGN_FLOAT_DOUBLE__


void L_test_for_aligned_malloc(size_t blo)
{
	double *buf[100];
	int i;
	size_t delta;
	bool alineado = true;

	for (i=0; i<100; i++)
	{
		buf[i] = new double[rand()%100];
		delta = (size_t)(buf[i]) % blo;
		if (delta != 0)
			alineado = false;
	}
	for (i=0; i<100; i++)
		delete[] buf[i];
	if (alineado)
		printf("new[] devuelve bloques alineados\n");
	else
		printf("new[] no devuelve bloques alineados (blo=%d delta = %d)\n", blo, delta);
}

void L_ArgcArgv::define(const char *cmdList) // Los ° permiten especificar espacios en un mismo argumento
{
	L_PUSH_EXECUTING_FN("L_ArgcArgv::define");
	int largoMax=0;
	int length;
	int i;
	int len;
	int nelem=0;
	char *cmdTemp;
	bool primerCaracter;

	len=(int)strlen(cmdList);
	cmdTemp=new char[len+1];
	strcpy(cmdTemp,cmdList);

	primerCaracter=true;
	length=0;
	for (i=0; i<=len; i++)
	{
		if (primerCaracter==true && cmdTemp[i]!=0 && cmdTemp[i]!=' ')
		{
			nelem++;
			length=0;
			primerCaracter=false;
		}
		if (cmdTemp[i]==' '||cmdTemp[i]==0)
		{
			cmdTemp[i]=0;
			if (length>largoMax)
				largoMax=length;
			primerCaracter=true;
		}
		else
		{
			length++;
			if (cmdTemp[i]=='°')
				cmdTemp[i]=' ';
		}
	}
	if (argv!=NULL)
		L_delete2d<char>(argv);
	argv=L_new2d<char>(nelem,largoMax+1);
	primerCaracter=true;
	argc=0;
	for (i=0; i<=len; i++)
	{
		if (primerCaracter==true && cmdTemp[i]!=0)
		{
			strcpy(argv[argc++],&cmdTemp[i]);
			primerCaracter=false;
		}
		if (cmdTemp[i]==0)
			primerCaracter=true;
	}
	delete[] cmdTemp;
	L_POP_EXECUTING_FN("L_ArgcArgv::define");
}

void L_calcGausH(std::vector<double> &buf, double s, int width, int facpaso) // compute kernel gaussiano
{
	double cen;
	int i;
	double facnorm;
	double facexp;

	buf.resize(width);

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif
	cen=(width-1 - 0)/2.0;
	facnorm=1/(sqrt(2*M_PI)*s);
	facexp=-0.5/(s*s);
	facnorm*=facpaso;
	for (i=0; i<width; i++)
		buf[i]=exp(facexp*(i-cen)*(i-cen))*facnorm;
	return;
}

void L_calcGausHNorm(std::vector<double> &buf, double s, int width, int facpaso) // compute kernel gaussiano
{
	double cen;
	int i;
	double facexp;
	double sum=0;

	buf.resize(width);

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif
	cen=(width-1 - 0)/2.0;
	facexp=-0.5/(s*s);
	for (i=0; i<width; i++)
		sum+=buf[i]=exp(facexp*(i-cen)*(i-cen));
	sum=facpaso/sum;
	for (i=0; i<width; i++)
		buf[i]*=sum;
	return;
}

void L_calcGausHInt(std::vector<double> &buf, double s, int width, int facpaso) // compute kernel gaussiano
{
	int i, j;
	double facnorm;
	double facexp;
	double x;

	buf.resize(width);

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif
	x=-width/2.0 + 0.5/10;
	facexp=-0.5/(s*s);
	facnorm=1/(sqrt(2*M_PI)*s)/10;
	for (i=0; i<width; i++)
	{
		buf[i]=0;
		for (j=0; j<10; j++)
		{
			buf[i]+=exp(facexp*x*x);
			x+=1.0/10;
		}
		buf[i]*=facnorm*facpaso;
	}
	return;
}

void L_calcGaborSeparable_cos_cos(std::vector<double> &h, double s, double f, double theta, int width, int facpaso)
{
	double cen;
	int i;
	double facexp;
	double factrig;
	double sum=0;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif
	h.resize(width);

	cen=(width-1 - 0)/2.0;
	facexp=-0.5/(s*s);
	factrig=2*M_PI*f*cos(theta);
	for (i=0; i<width; i++)
	{
		sum+=h[i]=exp(facexp*(i-cen)*(i-cen));
		h[i]*=cos(factrig*(i-cen));
	}
	sum=facpaso/sum;
	for (i=0; i<width; i++)
		h[i]*=sum;
	return;
}

void L_calcGaborSeparable_cos_sin(std::vector<double> &h, double s, double f, double theta, int width, int facpaso)
{
	double cen;
	int i;
	double facexp;
	double factrig;
	double sum=0;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif

	h.resize(width);

	cen=(width-1 - 0)/2.0;
	facexp=-0.5/(s*s);
	factrig=2*M_PI*f*sin(theta);
	for (i=0; i<width; i++)
	{
		sum+=h[i]=exp(facexp*(i-cen)*(i-cen));
		h[i]*=cos(factrig*(i-cen));
	}
	sum=facpaso/sum;
	for (i=0; i<width; i++)
		h[i]*=sum;
	return;
}

void L_calcGaborSeparable_sin_cos(std::vector<double> &h, double s, double f, double theta, int width, int facpaso)
{
	double cen;
	int i;
	double facexp;
	double factrig;
	double sum=0;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif

	h.resize(width);

	cen=(width-1 - 0)/2.0;
	facexp=-0.5/(s*s);
	factrig=2*M_PI*f*cos(theta);
	for (i=0; i<width; i++)
	{
		sum+=h[i]=exp(facexp*(i-cen)*(i-cen));
		h[i]*=sin(factrig*(i-cen));
	}
	sum=facpaso/sum;
	for (i=0; i<width; i++)
		h[i]*=sum;
	return;
}

void L_calcGaborSeparable_sin_sin(std::vector<double> &h, double s, double f, double theta, int width, int facpaso)
{
	double cen;
	int i;
	double facexp;
	double factrig;
	double sum=0;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

#ifdef L_BASIC_DEBUG
	if (width<=0)
		L_hard_shutdown("L_calcGausH : non-positive size allocation");
#endif

	h.resize(width);

	cen=(width-1 - 0)/2.0;
	facexp=-0.5/(s*s);
	factrig=2*M_PI*f*sin(theta);
	for (i=0; i<width; i++)
	{
		sum+=h[i]=exp(facexp*(i-cen)*(i-cen));
		h[i]*=sin(factrig*(i-cen));
	}
	sum=facpaso/sum;
	for (i=0; i<width; i++)
		h[i]*=sum;
	return;
}

void L_intBits(int *bits, int n, int log2n) // rius_yong
{
	throw_L_ArgException_if(bits == NULL, "L_intBits()");

	int k, l, r, i;
	bits[0]=0; bits[1]=n/2;
	for(i=2,l=2,r=n/4; r; l*=2,r/=2)
	{
		for (k=0; k<l; k++, i++)
		{
			bits[i] =bits[k] +r;
		}
	}
}

int L_ComplexDouble::solveQuadraticEquation(double a, double b, double c, L_ComplexDouble &x1, L_ComplexDouble &x2)
{
	L_ComplexDouble Vb2_4ac;
	Vb2_4ac.calcSqrt(b*b-4*a*c);
	x1 = ( -b - Vb2_4ac ) / (2*a);
	x2 = ( -b + Vb2_4ac ) / (2*a);
	return 2;
}

bool L_Random_plus_LM_optimizer::readyToWork(bool print)
{
	LM.readyToWork(print);
	if (!listoParaTrabajarBase())
	{
		if (print)
			printf("L_OptimizaAzar::listoParaTrabajarBase() no listo\n");
		return false;
	}
	if (vectorLimInf.size() == 0)
	{
		if (print)
			printf("L_OptimizaAzar::vectorLimInf no definido\n");
		return false;
	}
	if (vectorLimSup.size() == 0)
	{
		if (print)
			printf("L_OptimizaAzar::vectorLimSup no definido\n");
		return false;
	}
	if (nIterationsMax==16777215L)
	{
		if (print)
			printf("L_OptimizaAzar::nIterationsMax no definido\n");
		return false;
	}
	return true;
}


double L_Random_plus_LM_optimizer::minimize(const void *object, std::vector<double> &vect, double (*fnToMinimize)(const void *, double *))
{
	int i, j;
	L_Matrix errVect;
	double errMin = 1.0e60;
	for (i=0; i<nIterationsMax; i++)
	{
		paths.resize(paths.size()+1);
		paths[i].puntos.resize(2);
		paths[i].errores.resize(2);

		for (j=0; j<(int)vect.size(); j++)
			vect[i] = L_RANDOM(vectorLimInf[i],  vectorLimSup[i]);

		paths[i].puntos[0].copyVectorFrom(vect);
		paths[i].errores[0] = (*fnToMinimize)(object, &(vect[0]));
		if (paths[i].errores[0] < errMin)
			errMin = paths[i].errores[0];

		if (usarLM)
		{
			LM.minimize(object, vect, fnToMinimize);

			paths[i].puntos[1].copyVectorFrom(vect);
			paths[i].errores[1] = (*fnToMinimize)(object, &(vect[0]));
			if (paths[i].errores[0] < errMin)
				errMin = paths[i].errores[0];
		}
	}
	return errMin;
}


#if defined __TURBOC__
const L_ImageGrayDouble ** L_newArreglo(int nelem, const L_ImageGrayDouble *elem1, ...)
{
	va_list ppp;
	const L_ImageGrayDouble ** ptr;
	va_start(ppp, elem1);
	ptr=L_vNewArray<const L_ImageGrayDouble *>(nelem, elem1, ppp);
	va_end(ppp);
	return ptr;
}
#endif

void L_Int::print(int d)
{
	printf("%d\n", d);
}

void L_Double::print(double d)
{
	printf("%f\n", d);
}

double L_GaussianNoise()
{
	static bool iset=false;
	static double gset;
	double fac,rsq,v1,v2;
	if (iset == false)
	{
		do
		{
			v1=2.0*rand()/(RAND_MAX+1.0)-1.0;
			v2=2.0*rand()/(RAND_MAX+1.0)-1.0;
			rsq=v1*v1+v2*v2;// Ver si (v1,v2) esta en el circulo unitario
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		//Transformacion de Box-Muller entrega 2 valores gaussianos. Uno se guarda para despues
		gset=v1*fac;
		iset=true;
		return v2*fac;
	}
	else
	{
		iset=false;
		return gset;
	}
}


////// Clases propias



void L_HuffmanEncoder::test(bool soloCompr)
{
	size_type i, j, L1=320*200, L2=640*480;
	long nint = 10000;
	std::vector<char> arr1, arr2, arrc;
	double entrop=0, ratio = 0;
	double t1, t2;

	srand((unsigned int)time(NULL));

	t1 = L_TIME();
	for (i=0; i<nint; i++)
	{
		arr1.resize((rand()%(L2-L1))+L1);
		for (j=0; j<(size_type)arr1.size(); j++)
			arr1[j] = (char)(127 + (j + j*j + j*j*j) % 30 + (j*j*j*j) % 30 + (j*j*j*j*j) % 30);
		// entropia
		entrop += entropy_in_bits(arr1);
		arr2 = arr1;
		encodeAll(arr1, arrc);
		ratio += arrc.size() / (double)arr1.size();
		if (soloCompr)
			continue;
		if ( decodeAll(arrc, arr2) == false )
		{
			printf("\n");
			printf("Error en L_HuffmanTree : error al descomprimir\n");
			return;
		}
		for (j=0; j<(size_type)arr1.size(); j++)
		{
			if (arr1[j] != arr2[j])
			{
				printf("\n");
				printf("Error en L_HuffmanTree : decodificacion incorrecta\n");
				return;
			}
		}
		if ( i % (nint / 100) == 0)
			printf("\rAvance test: %ld%% (%ld compresiones/descompresiones)   ", 100*i/nint, i);
	}
	printf("\n");
	entrop /= nint;
	t2 = L_TIME();
	printf("L_ArbolHuffman paso la test con %ld intentos, %.2f [fps]\n",nint, nint/(t2-t1));
	printf("Compresion %.3f%% (length (%ld-%ld))\n", 100*ratio/nint, L1, L2);
	printf("Entropia real: %.2f[bits]\n", entrop);
	printf("Entropia equivalente compresion: %.2f[bits]\n", 8*ratio/nint);
}

void L_HuffmanEncoder::testEncoderRLE()
{
	size_type i, j, k, enne;
	std::vector<char> input, RLE, entrada2;
	int nint = 1000;

	srand((unsigned int)time(NULL));

	FILE *fp = fopen("sale.txt", "w");

	for (i=0; i<nint; i++)
	{
		input.resize(3000);
		for (j=0; j<(size_type)input.size(); j++)
		{
			input[j] = rand()%2;
			if (j > 3 && input[j] == 0 && input[j-1] == 0 && input[j-2] == 0)
			{
				enne = rand()%1024;
				for (k=0; k<enne && k+j < (size_type)input.size(); k++)
					input[k+j] = 0;
				j+=k;
				if (j<(size_type)input.size())
					input[j] = 1;
			}
		}
		encodeRLE0(input, RLE);
		decodeRLE0(RLE, entrada2);
		for (j=0; j<(size_type)input.size(); j++)
		{
			k = (unsigned char)input[j];
			if (k < 10)
				fprintf(fp, "%d", k);
			else
				fprintf(fp, "(%d)", k);
		}
		fprintf(fp, "\n");
		for (j=0; j<(size_type)RLE.size(); j++)
		{
			k = (unsigned char)RLE[j];
			if (k < 10)
				fprintf(fp, "%d", k);
			else
				fprintf(fp, "(%d)", k);
		}
		fprintf(fp, "\n");
		for (j=0; j<(size_type)entrada2.size(); j++)
		{
			k = (unsigned char)entrada2[j];
			if (k < 10)
				fprintf(fp, "%d", k);
			else
				fprintf(fp, "(%d)", k);
		}
		fprintf(fp, "\n");
		if (entrada2.size() != input.size())
		{
			fprintf(fp,	"Error en L_HuffmanEncoder::testEncoderRLE()  (1)\n");
			printf("Error en L_HuffmanEncoder::testEncoderRLE()  (1)\n");
			fclose(fp);
			return;
		}
		for (j=0; j<(size_type)input.size(); j++)
		{
			if (input[j] != entrada2[j])
			{
				fprintf(fp, "Error en L_HuffmanEncoder::testEncoderRLE()  (2)\n");
				printf("Error en L_HuffmanEncoder::testEncoderRLE()  (2)\n");
				fclose(fp);
				return;
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "L_HuffmanEncoder::testEncoderRLE() paso la test\n");
	printf("L_HuffmanEncoder::testEncoderRLE() paso la test\n");
	fclose(fp);
}

void L_HuffmanEncoder::testTree()
{
	size_type i, j;
	int nint = 100;
	std::vector<char> arr1;
	L_HuffmanEncoder cod1;

	srand((unsigned int)time(NULL));

	for (i=0; i<nint; i++)
	{
		arr1.resize((rand()%1000)+10);
		for (j=0; j<(size_type)arr1.size(); j++)
			arr1[j] = 127 + (char)(rand()%60);
		cod1.clear();

		cod1.build(&(arr1[0]), (size_type)arr1.size());
		if (cod1.tree.search_for_incoherences_in_tree() == false)
		{
			printf("L_ArbolHuffman::testTree() fallo la test (1)\n");
			return;
		}
		if (cod1.check_for_incoherences_beetween_tree_and_paths() == false)
		{
			printf("L_ArbolHuffman::testTree() fallo la test (2)\n");
			return;
		}
		cod1.tree.clear();
		if (cod1.buildTreeFromPaths() == false)
		{
			printf("L_ArbolHuffman::testTree() fallo la test (3)\n");
			return;
		}
		if (cod1.tree.search_for_incoherences_in_tree() == false)
		{
			printf("L_ArbolHuffman::testTree() fallo la test (4)\n");
			return;
		}
		if (cod1.check_for_incoherences_beetween_tree_and_paths() == false)
		{
			printf("L_ArbolHuffman::testTree() fallo la test (5)\n");
			return;
		}
	}
	printf("L_ArbolHuffman::testTree() paso la test con %d intentos\n", nint);
}

void L_HuffmanEncoder::testReadWriteHuffman()
{
	size_type i, j;
	int nint = 100;
	std::vector<char> arr1, huf;
	L_HuffmanEncoder cod1, cod2;
	double t1, t2;

	srand((unsigned int)time(NULL));

	t1 = L_TIME();
	for (i=0; i<nint; i++)
	{
		arr1.resize((rand()%1000)+10);
		for (j=0; j<(size_type)arr1.size(); j++)
			arr1[j] = 127 + (char)(rand()%60);
		cod1.clear();
		cod2.clear();

		cod1.build(&(arr1[0]), (size_type)arr1.size());
		cod1.tree.clear();
		if (cod1.buildTreeFromPaths() == false)
		{
			printf("L_ArbolHuffman::testReadWriteHuffman() fallo la test (1)\n");
			return;
		}

		cod1.writeHuffman(huf);
		cod2.readHuffman(huf);

		for (j=0; j<256; j++)
		{
			if (cod1.codeMapping[j].trainOfBits != cod2.codeMapping[j].trainOfBits || cod1.codeMapping[j].numBits != cod2.codeMapping[j].numBits)
			{
				printf("L_ArbolHuffman::testReadWriteHuffman() fallo la test (2)\n");
				return;
			}
		}
	}
	t2 = L_TIME();
	printf("L_ArbolHuffman::testReadWriteHuffman() paso la test con %d intentos, %.3f[fps]\n", nint, nint/(t2-t1));
}

void L_HuffmanEncoder::testEncoding()
{
	size_type i, j;
	size_type lenFinal;
	int nint = 100;
	std::vector<char> arr1, huf, arrc;
	L_HuffmanEncoder cod1;
	double t1, t2;

	srand((unsigned int)time(NULL));

	t1 = L_TIME();
	for (i=0; i<nint; i++)
	{
		arr1.resize((rand()%1000)+10);
		for (j=0; j<(size_type)arr1.size(); j++)
			arr1[j] = 127 + (char)(rand()%60);
		cod1.clear();
		cod1.build(&(arr1[0]), (size_type)arr1.size());
		cod1.buildCanonicalOrdering();
		cod1.tree.destroyTree();
		cod1.buildTreeFromPaths();
		cod1.check_for_incoherences_beetween_tree_and_paths();
		cod1.encode(&(arr1[0]), (size_type)arr1.size());
		lenFinal = (size_type)cod1.res.size();
		arrc = cod1.res;
		cod1.decode(&(arrc[0]), (size_type)arrc.size(), (size_type)arr1.size());
		for (j=0; j<(size_type)arr1.size(); j++)
		{
			if (arr1[j] != cod1.res[j])
			{
				printf("L_HuffmanEncoder::testEncoding() fallo la test\n");
			}
		}
	}
	t2 = L_TIME();
	printf("L_HuffmanEncoder::testEncoding() paso la test con %d intentos\n", nint);
	printf("%lf fps\n", nint/(t2-t1));
}

void L_HuffmanEncoder::testBits()
{
	size_type i, k;
	size_type bit, bitVal, nBits;
	int nint = 1000;
	std::vector<char> arr1, arr2;

	srand((unsigned int)time(NULL));

	for (i=0; i<nint; i++)
	{
		arr1.resize((rand()%1000)+10);
		for (k=0; k<(size_type)arr1.size(); k++)
			arr1[k] = rand()%256;
		arr2.resize(arr1.size());

		nBits = (size_type)(arr1.size()*8);

		for (bit=0; bit < nBits; bit++)
		{
			bitVal = L_CH_getBit(&(arr1[0]), bit);
			bitVal = bitVal * (1+rand()%100); // El bit no tiene por que valer 1...
			L_CH_setBit(&(arr2[0]), bit, bitVal);
		}

		for (k=0; k<(size_type)arr2.size(); k++)
		{
			if (arr1[k] != arr2[k])
			{
				printf("L_ArbolHuffman::testBits() fallo la test\n");
				return;
			}
		}
	}
	printf("L_ArbolHuffman::testBits() paso la test con %d intentos\n", nint);
}


L_HardErrorException::L_HardErrorException(const char *file, long line, const char *message, va_list ppp) nothrows()
{
	char mens[1024];
	vsprintf(mens, message, ppp);
	va_end(ppp);  // cerramos aqui ppp

	this->message = NULL;
	this->messageLen = 0;
	this->fileName = NULL;
	this->fileNameLen = 0;

	if (message != NULL)
	{
		messageLen=(int)(strlen(message)+1);
		this->message=new char[messageLen];
		if (this->message != NULL)
			strcpy(this->message, message);
	}

	if (file != NULL)
	{
		fileNameLen=(int)(strlen(file)+1);
		this->fileName=new char[fileNameLen];
		if (this->fileName != NULL)
			strcpy(this->fileName, file);
	}
	lineOfCode=line;
}

L_HardErrorException::L_HardErrorException(const L_HardErrorException& other) nothrows()
{
	message = NULL;
	messageLen = 0;
	fileName = NULL;
	fileNameLen = 0;

	if (other.message!=NULL)
	{
		message=new char[other.messageLen];
    	if (message!=NULL)
			strcpy(message, other.message);
	}
	if (other.fileName!=NULL)
	{
		fileName=new char[other.fileNameLen];
    	if (fileName!=NULL)
			strcpy(fileName, other.fileName);
	}
	lineOfCode=other.lineOfCode;
}

L_HardErrorException& L_HardErrorException::operator =(L_HardErrorException& other) nothrows()
{
	if (message != NULL)
		delete[] message;
	message = NULL;
	messageLen = 0;
	if (fileName != NULL)
		delete[] fileName;
	fileName = NULL;
	fileNameLen = 0;

	if (other.message!=NULL)
	{
		message=new char[other.messageLen];
    	if (message!=NULL)
			strcpy(message, other.message);
	}
	if (other.fileName!=NULL)
	{
		fileName=new char[other.fileNameLen];
    	if (fileName!=NULL)
			strcpy(fileName, other.fileName);
	}
	lineOfCode=other.lineOfCode;
	return *this;
}

#ifdef L_USE_NEW_PARAMETER_SERVER
void L_RefData::copia(const L_RefData &other)
{
	if (other.tipo!=L_ty_undefined && tipo!=L_ty_undefined)
	{
		throw_L_ArgException_if(other.tipo==tipo,"L_RefData::copia")
		other.copyValueOn(ptr);
	}
	else if (other.hasNum())
	{
		if (tipo==L_ty_undefined)
			setType(other.tipo);
		setValue(other.getValueDouble());
	}
	else
	{
		if (tipo==L_ty_undefined)
			setType(other.tipo);
		setValue(other.getValueString());
	}
}

bool L_RefData::readValuesFromLine(const char *buf)
{
	double num;

	while((*buf==' ' || *buf=='\t') && *buf!=0)
		buf++;
	if (*buf=='(')
	{
		while ( *buf!='0' && *buf!=')')
			buf++;
		if (*buf==')')
			buf++;
	}
	if (*buf==0)
		return false;
	while((*buf==' ' || *buf=='\t') && *buf!=0)
		buf++;

	if (buf[0]==L_REFDATA_SEPARATOR) // Comilla de apertura, es un string
	{
		L_String s(&buf[1]), s2;
		bool cerrado = false;
		int i;
		for (i=0; i<s.length; i++)
		{
			if (s.str[i] == L_REFDATA_SEPARATOR) // Comilla de cierre
			{
				s.str[i]=0;
				s2=s.str;
				cerrado=true;
				break;
			}
		}
		if (cerrado)
			setValue(s2);
		else
		{
			setValue(L_String(""));
			printf("Archivo de parametros con string incompleto\n");
		}
	}
	else if (buf[0]=='t' && buf[1]=='r' && buf[2]=='u' && buf[3]=='e')
	{
		setValue(1);
	}
	else if (buf[0]=='f' && buf[1]=='a' && buf[2]=='l' && buf[3]=='s' && buf[4]=='e')
	{
		setValue(0);
	}
	else
	{
		if (sscanf(buf,"%lf",&num) == 0)
			return false;
		setValue(num);
	}
	return true;
}

void L_ParamLabel::setName(const char *name)
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	if (this->name != NULL)
		delete[] this->name;
	this->name=new char[strlen(name)+1];
	strcpy(this->name,name);
}

void L_ParamLabel::setClass(const char *class_name)
{
	if (class_name==NULL)
		return;
	if (this->class_name != NULL)
		delete[] this->class_name;
	this->class_name=new char[strlen(class_name)+1];
	strcpy(this->class_name,class_name);
}

void L_ParamLabel::setFromString(const char *value)
{
	data.readValuesFromLine(value);
}

int L_ParamLabel::cmpInv(const void *a, const void *b)
{
	const char *s1,*s2;
	int i, j;
	s1=((L_ParamLabel *)a)->name;
	s2=((L_ParamLabel *)b)->name;
	i=(int)strlen(s1);
	j=(int)strlen(s2);
	if (i>j)
		return 1;
	else if(i<j)
		return -1;
	while (true)
	{
		if (s1[i]>s1[j])
			return 1;
		else if (s1[i]<s1[j])
			return -1;
		if (i==0)
		{
			if (j==0)
				return 0;
			return -1;
		}
		else if (j==0)
			return 1;
		i--;
		j--;
	}
	#ifndef __TURBOC__
	return 0;
	#endif
}
int L_ParamLabel::cmpNorm(const void *a, const void *b)
{
	return strcmp(((L_ParamLabel *)a)->name,((L_ParamLabel *)b)->name);
}

void L_ParamBlock::revMemo()
{
	if (nPaTot >= nPaMem)
	{
		L_ParamLabel *pa_other;
		L_ParamLabel *basura;
		int i;
		pa_other=new L_ParamLabel[nPaMem + 5];
		for (i=0; i<nPaTot; i++)
		{
			pa[i].swap(pa_other[i]);
		}
		nPaMem+=5;
		basura=pa;
		pa=pa_other;
		if (basura!=NULL)
			delete[] basura;
	}
}

L_RefData *L_ParamBlock::findVariable(const char *name)
{
	L_ParamLabel temp;
	L_ParamLabel *ptr=NULL;
	temp.name=const_cast<char *>(name); // Copia de referencia temporal, se le quita el const, no se debe modificar!!!
	if (!ordered)
		sort();
	if (nPaTot>0)
		ptr=(L_ParamLabel *)bsearch(&temp, pa, nPaTot, sizeof(L_ParamLabel), &L_ParamLabel::cmpNorm);
	if (ptr==NULL)
	{
		temp.name=NULL; // Destruccion de referencia tenporal
		return NULL;
	}
	temp.name=NULL; // Destruccion de referencia temporal
	return &(ptr->data);
}

bool L_ParamLabelList::readFile(FILE *fp, bool cerrarArchivo)
{
	L_ParamLabel param;
	char lin[200];
	char cla[100], name[100]; //, val[40];
	int i;
	if (fp==NULL)
		return false;
	while(fgets(lin, 199, fp)!=NULL)
	{
		//if (sscanf(lin,"%s%s%s",cla,name,val)<3) // linea incompleta
		//	continue;
		if (sscanf(lin,"%s%s",cla,name)<2)
			continue;
		if (*cla=='#') // comentario
			continue;
		i=0;
		while (lin[i]==' ')
			i++;
		while (lin[i]!=' ')
			i++;
		while (lin[i]==' ')
			i++;
		while (lin[i]!=' ')
			i++;
		while (lin[i]==' ')
			i++;
		param.setClass(cla);
		param.setName(name);
		//param.setFromString(val);
		param.setFromString(&lin[i]);
		push_back_swapping(param);
	}
	if (cerrarArchivo)
		fclose(fp);
	return true;
}

bool L_ParamLabelList::saveFile(FILE *fp, bool cerrarArchivo)
{
	L_ParamLabelNode *param;
	L_ParamLabelNode *paramInicioClase;
	L_ParamLabelNode *paramFinalClase;

	char lin[200];
	char name[200];
	char val[40];
	int len;
	int maxAnchoNombre=0;

	if (fp==NULL)
		return false;

	paramInicioClase=root;
	for(param=root; param!=NULL; param=param->sig)
	{
		if (strcmp(param->class_name,paramInicioClase->class_name)!=0) // Procesar una class_name
		{
			paramFinalClase=param;
			for (param=paramInicioClase; param!=paramFinalClase; param=param->sig)
			{
				strcpy(name,param->name);
				len=(int)strlen(name);
				while(len<maxAnchoNombre)
					name[len++]=' ';
				name[len]=0;
				sprintf(lin,"%s  %s %s\n", param->class_name, name, param->data.allToText(val));
				fputs(lin,fp);
			}
			paramInicioClase=paramFinalClase;
			fputs("\n",fp);
			maxAnchoNombre=0;
		}
		len=(int)strlen(param->name);
		if (len>maxAnchoNombre)
			maxAnchoNombre=len;
	}
	if (paramInicioClase!=NULL) // Procesar la última class_name
	{
		paramFinalClase=param;
		for (param=paramInicioClase; param!=paramFinalClase; param=param->sig)
		{
			strcpy(name,param->name);
			len=(int)strlen(name);
			while(len<maxAnchoNombre)
				name[len++]=' ';
			name[len]=0;
			sprintf(lin,"%s  %s %s\n", param->class_name, name, param->data.allToText(val));
			fputs(lin,fp);
		}
		paramInicioClase=paramFinalClase;
		fputs("\n",fp);
	}
	if (cerrarArchivo)
		fclose(fp);
	return true;
}

bool L_ParamBlock::updateValue(const L_ParamLabel &other)
{
	L_RefData *var;
	if (!belongTo(other))
		return false;
	if (!ordered)
		sort();

	var=findVariable(other.name);
	if (var==NULL)
		return false;
	*var=other.data;
	return true;
}

void L_ParamBlock::updateAssociatedValues(L_ParamLabelList &lista, bool borrarValoresUsados)
{
	L_ParamLabelNode **pptr;
	L_ParamLabelNode *basura;
	for (pptr=&lista.root; *pptr!=NULL;)
	{
		if (updateValue(**pptr)==true && borrarValoresUsados)
		{
			basura=*pptr;
			*pptr=(*pptr)->sig;
			basura->sig=NULL;
			delete basura;
			lista.size()--;
			if (lista.root==NULL)
				pptr=&lista.root;
		}
		else
			pptr=&(*pptr)->sig;
	}
	lista.pult=pptr;
}

void L_ParamBlock::buildListFromAssociatedValues(L_ParamLabelList &lista)
{
	int i;
	if (!ordered)
		sort();
	for (i=0; i<nPaTot; i++)
		lista.push_back(pa[i]);
}

void L_ParamManagerLocal::swap(L_ParamManagerLocal &other)
{
	// Permutar la componente params
	params.swap(other.params);
	// Permutar las demas componentes
	typedef L_ParamManagerLocal_POD sub;
	sub tmp;
	tmp=*this;
	(*this).sub::operator=(other);
	other.sub::operator=(tmp);
}

bool L_ParamManagerLocal::addChildren(L_ParamManagerLocal *nuevoHijo)
{
	if (nuevoHijo==NULL)
		return false;
	if (nHijosTot==nChildrenMem)
	{
		printf("L_ParamManagerLocal: Inicializado con memoria insuficiente\n");
		return false;
	}
	children[nHijosTot++]=nuevoHijo;
	return true;
}

void L_ParamManagerLocal::updateValues_me_and_children(L_ParamLabelList &lista)
{
	int i;
	params.updateAssociatedValues(lista, modifiedFlowOfInfo);
	for (i=0; i<nHijosTot; i++)
		children[i]->updateValues_me_and_children(lista);
}
void L_ParamManagerLocal::buildListOfValues_me_and_children(L_ParamLabelList &lista)
{
	int i;
	if (modifiedFlowOfInfo)
		params.buildListFromAssociatedValues(lista);
	for (i=0; i<nHijosTot; i++)
		children[i]->buildListOfValues_me_and_children(lista);
}
#else
bool L_RefData::readValuesFromLine(const char *buf)
{
	double num;

	while((*buf==' ' || *buf=='\t') && *buf!=0)
		buf++;
	if (*buf=='(')
	{
		while ( *buf!='0' && *buf!=')')
			buf++;
		if (*buf==')')
			buf++;
	}
	if (*buf==0)
		return false;
	while((*buf==' ' || *buf=='\t') && *buf!=0)
		buf++;

	if (buf[0]=='t' && buf[1]=='r' && buf[2]=='u' && buf[3]=='e')
	{
		setValue(1);
	}
	else if (buf[0]=='f' && buf[1]=='a' && buf[2]=='l' && buf[3]=='s' && buf[4]=='e')
	{
		setValue(0);
	}
	else
	{
		if (sscanf(buf,"%lf",&num) == 0)
			return false;
		setValue(num);
	}
	return true;
}

void L_ParamLabel::setName(const char *name)
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	if (this->name != NULL)
		delete[] this->name;
	this->name=new char[strlen(name)+1];
	strcpy(this->name,name);
}

void L_ParamLabel::setClass(const char *class_name)
{
	if (class_name==NULL)
		return;
	if (this->class_name != NULL)
		delete[] this->class_name;
	this->class_name=new char[strlen(class_name)+1];
	strcpy(this->class_name,class_name);
}

void L_ParamLabel::setFromString(const char *value)
{
	data.readValuesFromLine(value);
}

int L_ParamLabel::cmpInv(const void *a, const void *b)
{
	const char *s1,*s2;
	int i, j;
	s1=((L_ParamLabel *)a)->name;
	s2=((L_ParamLabel *)b)->name;
	i=(int)strlen(s1);
	j=(int)strlen(s2);
	if (i>j)
		return 1;
	else if(i<j)
		return -1;
	while (true)
	{
		if (s1[i]>s1[j])
			return 1;
		else if (s1[i]<s1[j])
			return -1;
		if (i==0)
		{
			if (j==0)
				return 0;
			return -1;
		}
		else if (j==0)
			return 1;
		i--;
		j--;
	}
	#ifndef __TURBOC__
	return 0;
	#endif
}
int L_ParamLabel::cmpNorm(const void *a, const void *b)
{
	return strcmp(((L_ParamLabel *)a)->name,((L_ParamLabel *)b)->name);
}

void L_ParamBlock::revMemo()
{
	if (nPaTot >= nPaMem)
	{
		L_ParamLabel *pa_other;
		L_ParamLabel *basura;
		int i;
		pa_other=new L_ParamLabel[nPaMem + 5];
		for (i=0; i<nPaTot; i++)
		{
			pa[i].swap(pa_other[i]);
		}
		nPaMem+=5;
		basura=pa;
		pa=pa_other;
		if (basura!=NULL)
			delete[] basura;
	}
}

L_RefData *L_ParamBlock::findVariable(const char *name)
{
	L_ParamLabel temp;
	L_ParamLabel *ptr=NULL;
	temp.name=const_cast<char *>(name); // Copia de referencia temporal, se le quita el const, no se debe modificar!!!
	if (!ordered)
		sort();
	if (nPaTot>0)
		ptr=(L_ParamLabel *)bsearch(&temp, pa, nPaTot, sizeof(L_ParamLabel), &L_ParamLabel::cmpNorm);
	if (ptr==NULL)
	{
		temp.name=NULL; // Destruccion de referencia tenporal
		return NULL;
	}
	temp.name=NULL; // Destruccion de referencia temporal
	return &(ptr->data);
}

bool L_ParamLabelList::readFile(FILE *fp, bool cerrarArchivo)
{
	L_ParamLabel param;
	char lin[200];
	char cla[100], name[100], val[40];
	if (fp==NULL)
		return false;
	while(fgets(lin, 199, fp)!=NULL)
	{
		if (sscanf(lin,"%s%s%s",cla,name,val)<3) // linea incompleta
			continue;
		if (*cla=='#') // comentario
			continue;
		param.setClass(cla);
		param.setName(name);
		param.setFromString(val);
		push_back_swapping(param);
	}
	if (cerrarArchivo)
		fclose(fp);
	return true;
}

bool L_ParamLabelList::saveFile(FILE *fp, bool cerrarArchivo)
{
	L_ParamLabelNode *param;
	L_ParamLabelNode *paramInicioClase;
	L_ParamLabelNode *paramFinalClase;

	char lin[200];
	char name[200];
	char val[40];
	int len;
	int maxAnchoNombre=0;

	if (fp==NULL)
		return false;

	paramInicioClase=root;
	for(param=root; param!=NULL; param=param->sig)
	{
		if (strcmp(param->c.class_name,paramInicioClase->c.class_name)!=0) // Procesar una class_name
		{
			paramFinalClase=param;
			for (param=paramInicioClase; param!=paramFinalClase; param=param->sig)
			{
				strcpy(name,param->c.name);
				len=(int)strlen(name);
				while(len<maxAnchoNombre)
					name[len++]=' ';
				name[len]=0;
				sprintf(lin,"%s  %s %s\n", param->c.class_name, name, param->c.data.allToText(val));
				fputs(lin,fp);
			}
			paramInicioClase=paramFinalClase;
			fputs("\n",fp);
			maxAnchoNombre=0;
		}
		len=(int)strlen(param->c.name);
		if (len>maxAnchoNombre)
			maxAnchoNombre=len;
	}
	if (paramInicioClase!=NULL) // Procesar la última class_name
	{
		paramFinalClase=param;
		for (param=paramInicioClase; param!=paramFinalClase; param=param->sig)
		{
			strcpy(name,param->c.name);
			len=(int)strlen(name);
			while(len<maxAnchoNombre)
				name[len++]=' ';
			name[len]=0;
			sprintf(lin,"%s  %s %s\n", param->c.class_name, name, param->c.data.allToText(val));
			fputs(lin,fp);
		}
		paramInicioClase=paramFinalClase;
		fputs("\n",fp);
	}
	if (cerrarArchivo)
		fclose(fp);
	return true;
}

bool L_ParamBlock::updateValue(const L_ParamLabel &other)
{
	L_RefData *var;
	if (!belongTo(other))
		return false;
	if (!ordered)
		sort();

	var=findVariable(other.name);
	if (var==NULL || var->ptr==NULL)
		return false;
	if (other.data.tipo == var->tipo)
	{
		other.data.copyValueOn(var->ptr);
	}
	else if (other.data.tipo==L_ty_undefined)
	{
		var->setValue(*(double*)other.data.ptr);
	}
	else
		return false;
	return true;
}

void L_ParamBlock::updateAssociatedValues(L_ParamLabelList &lista, bool borrarValoresUsados)
{
	L_ParamLabelNode **pptr;
	L_ParamLabelNode *basura;
	for (pptr=&lista.root; *pptr!=NULL;)
	{
		if (updateValue((**pptr).c)==true && borrarValoresUsados)
		{
			basura=*pptr;
			*pptr=(*pptr)->sig;
			basura->sig=NULL;
			delete basura;
			lista.n--;
			if (lista.root==NULL)
				pptr=&lista.root;
		}
		else
			pptr=&(*pptr)->sig;
	}
	lista.pult=pptr;
}

void L_ParamBlock::buildListFromAssociatedValues(L_ParamLabelList &lista)
{
	int i;
	if (!ordered)
		sort();
	for (i=0; i<nPaTot; i++)
		lista.push_back(pa[i]);
}


void L_ParamManagerLocal::swap(L_ParamManagerLocal &other)
{
	L_ParamManagerLocal **hijos_t = other.children;
	int nHijosTot_t = other.nHijosTot;
	int nHijosMem_t = other.nChildrenMem;
	bool alterarFlujoInfo_t = other.modifiedFlowOfInfo;

	params.swap(other.params);

	other.children=children;
	other.nHijosTot = nHijosTot;
	other.nChildrenMem = nChildrenMem;
	other.modifiedFlowOfInfo = modifiedFlowOfInfo;

	children=hijos_t;
	nHijosTot = nHijosTot_t;
	nChildrenMem = nHijosMem_t;
	modifiedFlowOfInfo = alterarFlujoInfo_t;
}

bool L_ParamManagerLocal::addChildren(L_ParamManagerLocal *nuevoHijo)
{
	if (nuevoHijo==NULL)
		return false;
	if (nHijosTot==nChildrenMem)
	{
		printf("L_ParamManagerLocal: Inicializado con memoria insuficiente\n");
		return false;
	}
	children[nHijosTot++]=nuevoHijo;
	return true;
}

void L_ParamManagerLocal::olvidaHijos()
{
	nHijosTot = 0;
}

void L_ParamManagerLocal::updateValues_me_and_children(L_ParamLabelList &lista)
{
	int i;
	params.updateAssociatedValues(lista, modifiedFlowOfInfo);
	for (i=0; i<nHijosTot; i++)
		children[i]->updateValues_me_and_children(lista);
}

void L_ParamManagerLocal::buildListOfValues_me_and_children(L_ParamLabelList &lista)
{
	int i;
	if (modifiedFlowOfInfo)
		params.buildListFromAssociatedValues(lista);
	for (i=0; i<nHijosTot; i++)
		children[i]->buildListOfValues_me_and_children(lista);
}
#endif
bool L_LevenbergMarquardt::readyToWork(bool print)
{
	if (!listoParaTrabajarBase())
	{
		if (print)
			printf("L_LevenbergMarquardt::listoParaTrabajarBase() no listo\n");
		return false;
	}
	if (nIterationsMax==16777215L)
	{
		if (print)
			printf("L_LevenbergMarquardt::nIterationsMax no definido\n");
		return false;
	}
	return true;
}

long L_Frac::mcd(long num, long den)
{
	long r;
	if (den > num) return mcd(den,num);
	r = num % den;
	if (r == 0)
	return den;
	else
		return mcd(den, r);
}

void L_Frac::simplify()
{
	long m;
	if (num>0)
	{
		m=mcd(num, den);
		num/=m;
		den/=m;
	}
	else
	{
		m=mcd(-num, den);
		num/=m;
		den/=m;
	}
}

int L_searchInSortedArray(const std::vector<double> valOrd, double val)
{
	int k1=0;
	int k2=(int)valOrd.size()-1;
	int k;
	while (k2-k1 > 1)
	{
		k=(k2+k1) >> 1;
		if (valOrd[k] > val)
			k2=k;
		else
			k1=k;
	}
	return k1;
}

void L_randomSelection(std::vector<int> &seleccion, int nCota)
{
	int i, j, act;
	throw_L_ArgException_if(seleccion.size() == 0 || (int)seleccion.size() > nCota, "L_randomSelection() : arreglo inadecuado");
	for (i=0; i<(int)seleccion.size(); i++)
	{
		act = rand() % (nCota - i) + i;
		for (j=0; j<i; j++)
		{
			if (act == seleccion[j])
				act = j;
		}
		seleccion[i] = act;
	}
}

// Algoritmo "University of Exeter", adaptado de http://www.bearcave.com/random_hacks/permute.html, se elimino la recursividad con paciencia...
void L_shiftPermutation_right(std::vector<int> &v, std::vector<int> &stackTmp)
{
	int n = (int)v.size();
	if (stackTmp.size() == 0)
	{
		stackTmp.resize(3*n+1);
		stackTmp[0] = 0; // zona del if
		stackTmp[1] = 0; // tmp
		stackTmp[2] = 0; // i
		stackTmp[stackTmp.size()-1] = 0; // profundidad
	}
	while (stackTmp[stackTmp.size()-1] < n-1)
	{
		if (stackTmp[3*stackTmp[stackTmp.size()-1]+2] < n)
		{
			if (stackTmp[3*stackTmp[stackTmp.size()-1]+0] == 0)
			{
				stackTmp[3*stackTmp[stackTmp.size()-1]+1] = v[stackTmp[3*stackTmp[stackTmp.size()-1]+2]];
				v[stackTmp[3*stackTmp[stackTmp.size()-1]+2]] = v[stackTmp[stackTmp.size()-1]];
				v[stackTmp[stackTmp.size()-1]] = stackTmp[3*stackTmp[stackTmp.size()-1]+1];
				stackTmp[3*stackTmp[stackTmp.size()-1]+0]=1;
				stackTmp[stackTmp.size()-1]++; // llamado
				stackTmp[3*stackTmp[stackTmp.size()-1]+0] = 0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+1] = 0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+2] = stackTmp[stackTmp.size()-1];
			}
			else
			{
				v[stackTmp[stackTmp.size()-1]] = v[stackTmp[3*stackTmp[stackTmp.size()-1]+2]];
				v[stackTmp[3*stackTmp[stackTmp.size()-1]+2]] = stackTmp[3*stackTmp[stackTmp.size()-1]+1];
				stackTmp[3*stackTmp[stackTmp.size()-1]+0]=0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+2]++;
			}
		}
		else
			stackTmp[stackTmp.size()-1]--;
		if (stackTmp[stackTmp.size()-1] < 0) // Cerrando el ciclo
		{
			// Simular condiciones iniciales
			stackTmp[stackTmp.size()-1] = 0;
			stackTmp[0] = 0;
			stackTmp[1] = 0;
			stackTmp[2] = 0;
		}
	}
	stackTmp[stackTmp.size()-1]--;
}

// El mismo de arriba, pero se reflejaron las permutaciones
void L_shiftPermutation_left(std::vector<int> &v, std::vector<int> &stackTmp)
{
	int n = (int)v.size();
	if (stackTmp.size() == 0)
	{
		stackTmp.resize(3*n+1);
		stackTmp[0] = 0; // zona del if
		stackTmp[1] = 0; // tmp
		stackTmp[2] = 0; // i
		stackTmp[stackTmp.size()-1] = 0; // profundidad
	}
	while (stackTmp[stackTmp.size()-1] < n-1)
	{
		if (stackTmp[3*stackTmp[stackTmp.size()-1]+2] < n)
		{
			if (stackTmp[3*stackTmp[stackTmp.size()-1]+0] == 0)
			{
				stackTmp[3*stackTmp[stackTmp.size()-1]+1] = v[v.size()-1-stackTmp[3*stackTmp[stackTmp.size()-1]+2]];
				v[v.size()-1-stackTmp[3*stackTmp[stackTmp.size()-1]+2]] = v[v.size()-1-stackTmp[stackTmp.size()-1]];
				v[v.size()-1-stackTmp[stackTmp.size()-1]] = stackTmp[3*stackTmp[stackTmp.size()-1]+1];
				stackTmp[3*stackTmp[stackTmp.size()-1]+0]=1;
				stackTmp[stackTmp.size()-1]++; // llamado
				stackTmp[3*stackTmp[stackTmp.size()-1]+0] = 0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+1] = 0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+2] = stackTmp[stackTmp.size()-1];
			}
			else
			{
				v[v.size()-1-stackTmp[stackTmp.size()-1]] = v[v.size()-1-stackTmp[3*stackTmp[stackTmp.size()-1]+2]];
				v[v.size()-1-stackTmp[3*stackTmp[stackTmp.size()-1]+2]] = stackTmp[3*stackTmp[stackTmp.size()-1]+1];
				stackTmp[3*stackTmp[stackTmp.size()-1]+0]=0;
				stackTmp[3*stackTmp[stackTmp.size()-1]+2]++;
			}
		}
		else
			stackTmp[stackTmp.size()-1]--;
		if (stackTmp[stackTmp.size()-1] < 0) // Cerrando el ciclo
		{
			// Simular condiciones iniciales
			stackTmp[stackTmp.size()-1] = 0;
			stackTmp[0] = 0;
			stackTmp[1] = 0;
			stackTmp[2] = 0;
		}
	}
	stackTmp[stackTmp.size()-1]--;
}

void L_getOrderedIndexesMajorMinor(std::vector<int> &inds, const std::vector<double> &puntaje)
{
	L_Array<L_IndOrd> ord((L_Array<L_IndOrd>::size_type)puntaje.size());
	int i;
	inds.resize(puntaje.size());
	for (i=0; i<(int)inds.size(); i++)
	{
		ord[i].i = i;
		ord[i].d = puntaje[i];
	}
	ord.sort(L_IndOrd::cmp);
	for (i=0; i<(int)inds.size(); i++)
		inds[i] = ord[(int)inds.size()-1-i].i;
}

void L_getOrderedIndexesMinorMajor(std::vector<int> &inds, const std::vector<double> &puntaje)
{
	L_Array<L_IndOrd> ord((L_Array<L_IndOrd>::size_type)puntaje.size());
	int i;
	inds.resize(puntaje.size());
	for (i=0; i<(int)inds.size(); i++)
	{
		ord[i].i = i;
		ord[i].d = puntaje[i];
	}
	ord.sort(L_IndOrd::cmp);
	for (i=0; i<(int)inds.size(); i++)
		inds[i] = ord[i].i;
}

void L_String::erase_spacing_characters_at_beginning_and_at_ending()
{
	size_type pos1=0; // El caracter en la cual parte el texto
	size_type pos2=size()-2; // El caracter en el cual termina el texto
	size_type i;

	if (size()==0 || c_str()==NULL)
		return;

	for (i=0; i<(size_type)size(); i++)
	{
		pos1 = i;
		if ((*this)[i] != ' ' && (*this)[i] != '\t' && (*this)[i] != '\n')
			break;
	}

	for (i=(size_type)size()-2; i>0; i--)
	{
		pos2 = i;
		if ((*this)[i] != ' ' && (*this)[i] != '\t' && (*this)[i] != '\n')
			break;
	}

	resize(pos2 - pos1 + 1 + 1); // Se agrega el 0 final
	if (size()<0)
		resize(1); // En este caso raro el string va a ser igual a ""
	for (i=0; i<(size_type)size()-1; i++)
		(*this)[i] = (*this)[i+pos1];
	(*this)[size()-1] = 0;
}

int L_String::sprintf_cpy(size_type lenmax, const char *format, ...)
{
	std::vector<char> buf(lenmax+1);
	int ret;
	va_list ppp;
#ifdef L_BASIC_DEBUG
	if (lenmax<=0)
		L_hard_shutdown("L_String::sprintf_cpy : non-positive size allocation");
#endif
	va_start(ppp, format);
	ret = vsprintf(&(buf[0]),format,ppp);
	va_end(ppp);
	resize((size_type)strlen(&(buf[0]))+1);
	strcpy(data(),&(buf[0]));
	return ret;
}

int L_String::sprintf_cat(size_type lenmax, const char *format, ...)
{
	std::vector<char> buf;
	int ret;
	va_list ppp;
#ifdef L_BASIC_DEBUG
	if (lenmax<=0)
		L_hard_shutdown("L_String::sprintf_cat : non-positive size allocation");
#endif
	buf.resize(lenmax + 1);
	va_start(ppp, format);
	ret = vsprintf(&(buf[0]),format,ppp);
	va_end(ppp);
	resize(size() + (size_type)strlen(&(buf[0])));
	strcat(data(), &(buf[0]));
	return ret;
}

void L_String::encript(unsigned long semilla)
{
	size_type i;
	char c;
	char mask;
	std::vector<char> copia0;

	resize((size_type)strlen(data())+2); // Cuidado: debe ser +2
	copia0.resize(size());
	strcpy(&(copia0[0]), data());

	for (i=0; i<(size_type)size()-1; i++)
	{
		c=operator[](i);
		operator[](i)&=0xF0;
		operator[](i)|=((c>>3)&1)<<0;
		operator[](i)|=((c>>2)&1)<<1;
		operator[](i)|=((c>>1)&1)<<2;
		operator[](i)|=((c>>0)&1)<<3;
	}

	for (i=0; i<(size_type)size()-1; i++)
		semilla+=((long)i%5+2)*(operator[](i)&(1<<(i%4)));
	for (i=0; i<4; i++)
		semilla=(semilla*7+11)%1715;
	for (i=0; i<(size_type)size()-1; i++)
	{
		semilla=(semilla*16+12)%1715;
		while (   ( semilla%16 & ~(1<<(i%4)) ) == 0   )
			semilla=(semilla*16+12)%1715;
		mask=(char)(semilla%16);
		mask&=~(1<<(i%4));
		operator[](i)=L_XOR(operator[](i),mask);
	}

	for (i=0; i<(size_type)size()-1; i++)
	{
		c=operator[](i);
		operator[](i)&=0xF0;
		operator[](i)|=((c>>3)&1)<<0;
		operator[](i)|=((c>>2)&1)<<1;
		operator[](i)|=((c>>1)&1)<<2;
		operator[](i)|=((c>>0)&1)<<3;
		if (operator[](i)==(char)127)
			operator[](i)=copia0[i];
	}
}

void L_String::codifBase16(const char *txt)
{
	size_type len;
	size_type i;
	len=(size_type)strlen(txt);
	resize(2*len+1);
	for (i=0; i<len; i++)
	{
		operator[](2*i)=cHexDigit(txt[i]>>4);
		operator[](2*i+1)=cHexDigit(txt[i]&0x0F);
	}
	operator[](2*i)=0;
}

void L_String::decodifBase16(const char *txt)
{
	size_type len;
	size_type i;
	len=(size_type)strlen(txt);
	resize(len/2+1);
	for (i=0; i<(size_type)size()-1; i++)
		operator[](i)=(char)(   ( cHexDigitInv(txt[2*i])<<4 ) + cHexDigitInv(txt[2*i+1])   );
	operator[](i)=0;
}

L_String& L_String::operator =(const L_String &s)
{
	if (s.size() == 0)
	{
		resize(1);
		operator[](0) = 0;
	}
	else
	{
		resize(s.size());
		throw_L_ArgException_if(data()==NULL, "L_String::operator =()");
		strcpy(data(), s.c_str());
	}
	return *this;
}

L_String operator+ (const L_String &s1, const L_String &s2)
{
	L_String ret;
	ret.resize(s1.size()+s2.size()-1);
	L_snprintf(&(ret[0]),ret.size(),"%s%s",(char *)s1.c_str(),(char *)s2.c_str());
	return ret; // Devolver ret sin destruirlo
}

L_String& L_String::operator +=(const L_String &s)
{
	L_String buf;
	buf.resize(size() + s.size());
	L_snprintf(buf.data(),buf.size(),"%s%s",c_str(),s.c_str());
	buf.swap(*this);
	return *this;
}

void L_String::setInteger(long num, int numberWidthZeros) // ej: (2,3) = 007
{
	size_type i;
	int numDig;
	long nTmp; // Numero de caracteres
	int signo=0;

	if (num<0) // Dejar "num" positivo y el signo en "signo"
	{
		signo=1;
		num=-num;
	}
	// Contar el numero de digitos
	nTmp=num;
	numDig = 0;
	while (nTmp>0)
	{
		nTmp/=10;
		numDig++;
	}
	// Adaptar el numero de digitos para que sea >= numberWidthZeros
	if (numberWidthZeros>numDig)
		numDig=numberWidthZeros;
	resize(numDig+signo+1); // Cantidad de bytes a pedir
	(*this)[numDig+signo]=0; // Fin del string
	i = numDig-1; // Indice de la ultima cifra
	for (; i>=0; i--)
	{
		(*this)[i+signo]=( (char)(num%10) )+'0'; // Escribir ultima cifra actual del numero
		num/=10; // Eliminar ultima cifra actual del numero
	}
	if (signo)
		(*this)[0]='-'; // Escribir el '-' si corresponde
}

void L_String::setDouble(double num, int nDecimals)
{
	char buf[100];
	char s[100];
	sprintf(s,"%%.%df",nDecimals);
	sprintf(buf,s,num);
	resize((size_type)strlen(buf) + 1);
	strcpy(data(),buf);
}

bool L_String::find(const char *subString, int &pos)
{
	char *r;
	r=strstr(data()+pos, subString);
	if (r==NULL)
		return false;
	pos = (int)(r-subString);
	return true;
}

void L_String::tokens(L_Array<L_String> &arrStr, const char *delimitadores)
{
	L_String s2(c_str());
	char *ptr;
	arrStr.resize(0);
	ptr = strtok(s2.data(), delimitadores);
	while (ptr!=NULL)
	{
		arrStr.resize(arrStr.size()+1);
		arrStr[arrStr.size()-1] = ptr;
		ptr = strtok(NULL, delimitadores);
	}
}


// Elige un color RGB al azar que se aleje del gris
void L_ShapeArray::genRandomColor(L_uchar& R, L_uchar& G, L_uchar& B)
{
	G=(L_uchar)(rand()*255.9999/RAND_MAX);
	if (G>127)
		B=(L_uchar)(rand()*99.9999/RAND_MAX);
	else
		B=(L_uchar)(255-rand()*99.9999/RAND_MAX);
	if (G+B<100)
		R=(L_uchar)(127L+rand()*126.9999/RAND_MAX);
	else
		R=0;
}

void L_ShapeArray::_drawShape(int im1, int im2, int x1, int y1, int x2, int y2, L_uchar R, L_uchar G, L_uchar B)
{
	v.push_back(L_Shape(x1, y1, im1, x2, y2, im2, 0.0,0.0,0.0, R, G, B, L_Shape_line));
}

void L_ShapeArray::_drawRectangle(int im1, int im2, int x1, int y1, int x2, int y2, L_uchar R, L_uchar G, L_uchar B)
{
	this->_drawShape(im1, im2, x1, y1, x2, y1, R, G, B);
	this->_drawShape(im1, im2, x2, y1, x2, y2, R, G, B);
	this->_drawShape(im1, im2, x2, y2, x1, y2, R, G, B);
	this->_drawShape(im1, im2, x1, y2, x1, y1, R, G, B);
}

void L_ShapeArray::_drawRotatedSquare(int im1, int im2, int xc, int yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B)
{
	int c, s;
	c=(int)(r*cos(ang));
	s=(int)(r*sin(ang));
	this->_drawShape(im1, im2, xc+s, yc-c, xc+c, yc+s, R, G, B);
	this->_drawShape(im1, im2, xc+c, yc+s, xc-s, yc+c, R, G, B);
	this->_drawShape(im1, im2, xc-s, yc+c, xc-c, yc-s, R, G, B);
	this->_drawShape(im1, im2, xc-c, yc-s, xc+s, yc-c, R, G, B);
}

void L_ShapeArray::_dibFlecha(int im1, int im2, int xc, int yc, double r, double ang, L_uchar R, L_uchar G, L_uchar B)
{
	double c=r*cos(ang);
	double s=r*sin(ang);
	this->_drawShape(im1, im2, xc, yc, (int)(xc+c), (int)(yc-s), R, G, B);
	this->_drawShape(im1, im2, (int)(xc+c), (int)(yc-s), (int)(xc+c-0.3*c+0.3*s), (int)(yc-s+0.3*s+0.3*c), R, G, B);
	this->_drawShape(im1, im2, (int)(xc+c), (int)(yc-s), (int)(xc+c-0.3*c-0.3*s), (int)(yc-s+0.3*s-0.3*c), R, G, B);
}

void L_ShapeArray::_drawCircle(int im1, int im2, int xc, int yc, double r, L_uchar R, L_uchar G, L_uchar B)
{
	v.push_back(L_Shape(xc, yc, im1, xc, yc, im2, r,r,0, R, G, B, L_Shape_circle));
}

void L_ShapeArray::_drawEllipse(int im1, int im2, int xc, int yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B)
{
	v.push_back(L_Shape(xc, yc, im1, xc, yc, im2, rx,ry,ang, R, G, B, L_Shape_ellipse));
}

void L_ShapeArray::_drawDashedEllipse(int im1, int im2, int xc, int yc, double rx, double ry, double ang, L_uchar R, L_uchar G, L_uchar B)
{
	v.push_back(L_Shape(xc, yc, im1, xc, yc, im2, rx,ry,ang, R, G, B, L_Shape_dashed_ellipse));
}

void L_ShapeArray::drawNumber(int xc, int yc, int num, int radio, L_uchar R, L_uchar G, L_uchar B)
{
	int cifra, iter = 0;
	double avance;
	int a=-1, b=-1, c=-1, d=-1, e=-1, f=-1, g=-1;
	avance = (log10((double)num) + (num<0))* radio * 1.3;
	//  aaa
	//  b c
	//  ddd
	//  e f
	//  ggg
	while (num != 0 || iter == 0)
	{
		if (num < 0)
		{
			drawLine(int(xc+avance + 0),  int(yc + radio/2),  int(xc+avance + radio),  int(yc + radio/2), R, G, B);
			avance -= radio*1.3;
			num = num * -1;
		}
		cifra = num % 10;
		num = num/10;
		switch(cifra)
		{

		case 0:
			a=1; b=1; c=1; d=0; e=1; f=1; g=1;
			break;
		case 1:
			a=0; b=0; c=1; d=0; e=0; f=1; g=0;
			drawLine(int(xc+avance + 0.75*radio),  int(yc + 0.25*radio),  int(xc+avance + radio),  int(yc + 0), R, G, B);
			break;
		case 2:
			a=1; b=0; c=1; d=1; e=1; f=0; g=1;
			break;
		case 3:
			a=1; b=0; c=1; d=1; e=0; f=1; g=1;
			break;
		case 4:
			a=0; b=1; c=1; d=1; e=0; f=1; g=0;
			break;
		case 5:
			a=1; b=1; c=0; d=1; e=0; f=1; g=1;
			break;
	//  aaa
	//  b c
	//  ddd
	//  e f
	//  ggg
		case 6:
			a=0; b=1; c=0; d=1; e=1; f=1; g=1;
			break;
		case 7:
			a=1; b=0; c=1; d=0; e=0; f=1; g=0;
			break;
		case 8:
			a=1; b=1; c=1; d=1; e=1; f=1; g=1;
			break;
		case 9:
			a=1; b=1; c=1; d=1; e=0; f=1; g=0;
			break;
		}
		if (a>0)
			drawLine(int(xc+avance + 0),  int(yc + 0),  int(xc+avance + radio),  int(yc + 0), R, G, B);
		if (b>0)
			drawLine(int(xc+avance + 0),  int(yc + 0),  int(xc+avance + 0),  int(yc + 0.5*radio), R, G, B);
		if (c>0)
			drawLine(int(xc+avance + radio),  int(yc + 0),  int(xc+avance + radio),  int(yc + 0.5*radio), R, G, B);
		if (d>0)
			drawLine(int(xc+avance + 0),  int(yc + 0.5*radio),  int(xc+avance + radio),  int(yc + 0.5*radio), R, G, B);
		if (e>0)
			drawLine(int(xc+avance + 0),  int(yc + 0.5*radio),  int(xc+avance + 0),  int(yc + radio), R, G, B);
		if (f>0)
			drawLine(int(xc+avance + radio),  int(yc + 0.5*radio),  int(xc+avance + radio),  int(yc + radio), R, G, B);
		if (g>0)
			drawLine(int(xc+avance + 0),  int(yc + radio),  int(xc+avance + radio),  int(yc + radio), R, G, B);
		avance -= 1.3*radio;
		iter++;
	}
}


void L_ShapeArray::changeSizeFactor(double sizeFactor)
{
	for (size_type i=0; i<(size_type)size(); i++)
	{
		v[i].x1=(size_type)(sizeFactor*v[i].x1);
		v[i].y1=(size_type)(sizeFactor*v[i].y1);
		v[i].x2=(size_type)(sizeFactor*v[i].x2);
		v[i].y2=(size_type)(sizeFactor*v[i].y2);
	}
}

void L_ShapeArray::drawTracks(L_Matrix &m, L_Array<bool> *v)
{
	size_type i, j, a,b,c,d;
	if (m.li < 1 || m.lj < 4)
		return;
	for (i=0; i<m.li; i++)
	{
		if (v != NULL && v->operator[](i) == false)
			continue;
		for (j=1; j<m.lj/2; j++)
		{
			a = (size_type)m(i,2*j-2);
			b = (size_type)m(i,2*j-1);
			c = (size_type)m(i,2*j+0);
			d = (size_type)m(i,2*j+1);
			drawLine(a, b, c, d, 255,255,255);
			drawLine(a+1, b+1, c+1, d+1, 255,255,255);
			drawLine(a+1, b-1, c+1, d-1, 255,255,255);
			drawLine(a-1, b-1, c-1, d-1, 255,255,255);
			drawLine(a-1, b+1, c-1, d+1, 255,255,255);
		}
	}
}

void L_ShapeArray::drawLineByParameters(double a, double b, double c, int lx, int ly, L_uchar R, L_uchar G, L_uchar B)
{
	double xb[4], yb[4];
	double x0=-1, y0=-1, x1=-1, y1=-1;

	// Recta: a*x+b*y+c=0
	// Se compute la interseccion con las 4 rectas que definen los bordes de la pantalla
	// Luego de eso, se eligen los puntos apropiados

	// Interseccion con x = 0 (vertical izquierda)
	xb[0] = 0;
	yb[0] = -c/b;

	// Interseccion con x=lx-1 (vertical derecha)
	xb[1] = lx-1;
	yb[1] = -a*(lx-1)/b -c/b;

	// Interseccion con y = 0 (horizontal superior)
	xb[2] = -c/a;
	yb[2] = 0;

	// Interseccion con y=ly-1 (horizontal inferior)
	xb[3] = -b*(ly-1)/a -c/a;
	yb[3] = ly-1;

	// Ver cuales de los cuatro puntos son validos -> no validos con x=-1
	if (yb[0] < 0 || yb[1] >= ly)
		xb[0] = -1;
	if (yb[1] < 0 || yb[1] >= ly)
		xb[1] = -1;
	if (xb[2] < 0 || xb[2] >= lx)
		xb[2] = -1;
	if (xb[3] < 0 || xb[3] >= lx)
		xb[3] = -1;

	// Elegir los dos puntos que definen la linea
	int i;
	for (i=0; i<4; i++)
	{
		if (xb[i] != -1)
		{
			if (x0 == -1)
			{
				x0 = xb[i];
				y0 = yb[i];
			}
			else
			{
				x1 = xb[i];
				y1 = yb[i];
			}
		}
	}

	// Finalmente agregar los vertices de la recta

	drawLine((int)x0, (int)y0, (int)x1, (int)y1, R, G, B);
}

void L_ShapeArray::write(FILE *fp)
{
	for (int i=0; i<(int)size(); i++)
	{
		fprintf(fp,"%d %d %d %d %d %d ", v[i].x1,v[i].y1,v[i].n1,v[i].x2,v[i].y2,v[i].n2);
		fprintf(fp,"%lg %lg %lg %d %d %d\n", v[i].rx, v[i].ry, v[i].ang, (int)v[i].R, (int)v[i].G, (int)v[i].B );
	}
	fprintf(fp,"---\n");
}

void L_ShapeArray::read(FILE *fp)
{
	L_Shape linea;
	char lin[1000];
	int R, G, B;
	while (fgets(lin, 999, fp)!=NULL)
	{
		if (strlen(lin)>=4 && strcmp(lin,"---\n")==0)
			break;
		sscanf(lin,"%d%d%d%d%d%d"
		"%lg%lg%lg"
		"%d%d%d", &linea.x1,&linea.y1,&linea.n1,&linea.x2,&linea.y2,&linea.n2,
			&linea.rx,&linea.ry,&linea.ang,&R,&G,&B);
		linea.R=(L_uchar)R;
		linea.G=(L_uchar)G;
		linea.B=(L_uchar)B;
		push_back(linea);
	}
	return;
}

/*
void L_ShapeArray::testShapes()
{
	L_ShapeArray lins;
	L_ImageRGBUchar imRGB;
	int x[5], y[5];
	int nClick = 0;
	int tipo = 0;
	L_VentanaImagen vent(0, 0, 320, 240, "Click para dibujar");
	L_VentanaImagen ventSel(400, 0, 320, 240, "Click para cambiar tipo");
	imRGB.reallocate(320,240);
	while (true)
	{
		if (vent.mouseE() != L_MouseLeido)
		{
			if (vent.mouseE() == L_MousePresiona)
			{
				x[nClick] = vent.mouseX();
				y[nClick] = vent.mouseY();
				vent.mouseE() = L_MouseLeido;
				nClick++;
			}
			vent.mouseE() = L_MouseLeido;
		}
		switch(tipo)
		{
		case 0: // Linea
			if (nClick == 2)
			{
				lins.drawLine(x[0],y[0],x[1],y[1]);
				imRGB.setZero();
				imRGB.genDrawing(lins);
				nClick = 0;
			}
			break;
		case 1: // Elipse
			if (nClick == 3)
			{
				double rx, ry, ang;
				double dx1, dy1, dx2, dy2;
				dx1 = x[1] - x[0];
				dy1 = y[1] - y[0];
				dx2 = x[2] - x[0];
				dy2 = y[2] - y[0];
				rx = sqrt(dx1*dx1 + dy1*dy1);
				ry = fabs(dx1*-dy2 + dy1*dx2)/rx; // d1 cruz d2
				ang = -atan2(dy1,dx1); // El angulo se mide hacia ese lado
				lins.drawDashedEllipse(x[0],y[0],rx,ry,ang,255,255,255);
				imRGB.setZero();
				imRGB.genDrawing(lins);
				nClick = 0;
			}
			break;
		}

		if (ventSel.mouseE() != L_MouseLeido)
		{
			if (ventSel.mouseE() == L_MousePresiona)
			{
				ventSel.mouseE() = L_MouseLeido;
				tipo++;
				if (tipo >= 2)
					tipo = 0;
			}
			ventSel.mouseE() = L_MouseLeido;
		}
		vent.dibuja(imRGB);
		vent.check();
	}
}
*/

void L_MatParser::genCanonicalExpr(const char *expression, L_String &ultExpr)
{
	const char *c1;
	char *c2;
	int uesize;
	ultExpr.resize((int)strlen(expression)+3);
	uesize = (int)ultExpr.size(); // It cannot be so long...

	c1=expression;
	c2=ultExpr.data();
	// Eliminacion de ' ' y '\t'
	*(c2++)='(';
	while (*c1!=0)
	{
		if (*c1 == '*' && *(c1-1) == '*') // Revisar ** -> ^
		{
			*(c2-1) = '^';
			uesize--;
		}
		else if (*c1!=' ' && *c1!='\t') // Revisar caracter no blanco
			*(c2++)=*c1;
		else
			uesize--;
		c1++;
	}
	ultExpr.resize(uesize);
	*(c2++)=')';
	*(c2++)=0;
}

void L_MatParser::buildFrom(const char *expression, bool debug)
{
	L_Array<L_MatParsParenNode *> mapaParen;
	typedef L_Array<L_MatParsParenNode *>::size_type arr_size_type;
	L_MatParsParenNode *parenPtr, *parenPtr2;
	L_MatParsParenNode *parenRaiz;
	L_String ultExpr;
	int i;

	genCanonicalExpr(expression, ultExpr);

	mapaParen.resize((arr_size_type)ultExpr.size()-1);
	for (i=0; i<(arr_size_type)ultExpr.size()-1; i++)
		mapaParen[i]=NULL;
	parenRaiz=new L_MatParsParenNode;
	parenRaiz->i1=0;
	mapaParen[0]=parenRaiz;
	parenRaiz->i2=(int)( ultExpr.size()-2 );
	parenPtr=parenRaiz;
	if (root!=NULL)
		root->destroyRec();
	else
		root=new L_MatParsNode;
	root->pre_i1=1;
	root->pos_i2=(int)( ultExpr.size()-3 );
	for (i=(int)(ultExpr.size()-2); i>=0; i--)
	{
		if (ultExpr[i]==')' || ultExpr[i]=='}' || ultExpr[i]==']')
		{
			parenPtr2=new L_MatParsParenNode;
			parenPtr2->i2=i;
			mapaParen[i]=parenPtr2;
			parenPtr2->father=parenPtr;
			parenPtr->children.push_back(parenPtr2);
			parenPtr=parenPtr2;
		}
		if (ultExpr[i]=='(' || ultExpr[i]=='{' || ultExpr[i]=='[')
		{
			parenPtr->i1=i;
			parenPtr=parenPtr->father;
		}
	}
	// Ahora, en mapaParen se tiene info sobre los grupos con parentesis
	analyzeParenthesisRecursive(root, mapaParen, ultExpr.c_str(), 0);
	if (debug)
		showStructureDebug();
	buildNodeRecursive(root, ultExpr);
	parenRaiz->destroyRec();
	delete parenRaiz;
}

void L_MatParser::analyzeParenthesisRecursive(L_MatParsNode *raizPars, L_Array<L_MatParsParenNode *> &mapaParen, const char *expr, int level)
{
	int primSuma=-1;
	int primExp=-1;
	int primMult=-1;
	int div=-1;
	int i;
	L_MatParsNode *nodo;

	// Transformar (expres) en expres

	raizPars->level=level;
	// Debug
	L_String strDeb;
	strDeb.copySubString(expr, raizPars->pre_i1, raizPars->pos_i2);

	while (mapaParen[raizPars->pos_i2] != NULL && mapaParen[raizPars->pos_i2]->i1==raizPars->pre_i1)
	{
		raizPars->pre_i1++;
		raizPars->pos_i2--;
	}

	if ( raizPars->pre_i1 >= raizPars->pos_i2 )
		return;

	for ( i=raizPars->pos_i2 ; i>=raizPars->pre_i1 ; i-- )
	{
		if (mapaParen[i]!=NULL) // Saltarse subexpresiones
		{
			i=mapaParen[i]->i1;
			continue;
		}
		if (expr[i]==',') // El de menor precedencia
		{
			div=i;
			raizPars->preced=3; // Para no poner zero, por si aca
			break;
		}
		else if ( (primSuma==-1 || primSuma==i+1) && (expr[i]=='+' || expr[i]=='-') )
		{
			primSuma=i;
			raizPars->preced=6;
		}
		else if ( primMult==-1 && (expr[i]=='*' || expr[i]=='/' || expr[i]==':') )
		{
			if (primSuma==i+1)
				primSuma=-1;
			primMult=i;
			raizPars->preced=9;
		}
		else if (primExp==-1 && expr[i]=='^') // '**' = '^'
		{
			if (primSuma==i+1)
				primSuma=-1;
			primExp=i;
			raizPars->preced=12;
		}
	}
	if (div==-1 && primSuma!=-1)
		div=primSuma;
	if (div==-1 && primMult!=-1)
		div=primMult;
	if (div==-1 && primExp!=-1)
		div=primExp;

	if (div!=-1)
	{
		raizPars->pre_i2=div-1;
		raizPars->pos_i1=div+1;
		if (raizPars->pre_i1<=raizPars->pre_i2)
		{
			nodo=new L_MatParsNode;
			nodo->pre_i1=raizPars->pre_i1;
			nodo->pos_i2=div-1;
			raizPars->preArg=nodo;
			analyzeParenthesisRecursive(raizPars->preArg, mapaParen, expr, level+1);
		}
		if (raizPars->pos_i1<=raizPars->pos_i2)
		{
			nodo=new L_MatParsNode;
			nodo->pre_i1=div+1;
			nodo->pos_i2=raizPars->pos_i2;
			raizPars->posArg=nodo;
			analyzeParenthesisRecursive(raizPars->posArg, mapaParen, expr, level+1);
		}
	}
	return;
}

void L_MatParser::showStructureDebug(FILE *fp)
{
	if (root!=NULL)
		showStructureDebugRecursive(fp, root, 0);
}

void L_MatParser::showStructureDebugRecursive(FILE *fp, L_MatParsNode *nodo, int level)
{
	int i;
	for (i=0; i<level; i++) fprintf(fp, "|");
	if (nodo->preArg!=NULL || nodo->posArg!=NULL)
	{
		fprintf(fp, "%c (op)\n", nodo->str.c_str()[0]); //ultExpr[nodo->pos_i1-1]);
		if (nodo->preArg!=NULL)
			showStructureDebugRecursive(fp, nodo->preArg, level+1);
		if (nodo->posArg!=NULL)
			showStructureDebugRecursive(fp, nodo->posArg, level+1);
	}
	else
	{
		fprintf(fp, "%s\n",nodo->str.c_str());
	}
}

bool L_MatParser::buildNodeRecursive(L_MatParsNode *raizPars, L_String &ultExpr)
{
	bool izqUtil=false, derUtil=false;
	//bool esFuncion=false, esOperador=false;
	int i;

	if (raizPars->preArg!=NULL && raizPars->posArg!=NULL)
	{
		//if (raizPars->pre_i1+1==raizPars->pos_i1) // No hay caracter entre izq y der
		//	esFuncion=true; // Hasta ahora, es funcion
		//else
		//	esOperador=true; // Hasta ahora, es operador
	}

	if (raizPars->preArg!=NULL)
	{
		izqUtil=buildNodeRecursive(raizPars->preArg, ultExpr);
		if (!izqUtil)
		{
			raizPars->preArg->destroyRec();
			delete raizPars->preArg;
			raizPars->preArg = NULL;
		}
	}
	if (raizPars->posArg!=NULL)
	{
		derUtil=buildNodeRecursive(raizPars->posArg, ultExpr);
		if (!derUtil)
		{
			raizPars->posArg->destroyRec();
			delete raizPars->posArg;
			raizPars->posArg = NULL;
		}
	}

	//if (!izqUtil || !derUtil)
	//{
		//esFuncion = false; // No puede ser función en este caso
		//esOperador = false; // No puede ser operador en este caso
	//}

	// Ahora hay que ver si el nodo actual es useful
	raizPars->useful=false; // Hasta ahora

	// Destruir los espacios en blanco a la derecha e izquierda de la expression
	while (raizPars->pre_i1+1<(int)ultExpr.size()-1 && (ultExpr[raizPars->pre_i1+1] == ' ' || ultExpr[raizPars->pre_i1+1] == '\t'))
		raizPars->pre_i1++;
	while (raizPars->pos_i2-1>=0 && (ultExpr[raizPars->pos_i2-1] == ' ' || ultExpr[raizPars->pos_i2-1] == '\t') )
		raizPars->pos_i2--;

	if (raizPars->pre_i1>=raizPars->pos_i2)
		raizPars->useful=false;
	if (raizPars->preArg==NULL)
	{
		if (raizPars->posArg==NULL) // Constante o variable
		{
			raizPars->mytype=LMP_constant;
			for (i=raizPars->pre_i1; i<=raizPars->pos_i2; i++)
			{
				if ((ultExpr[i]<'0' || ultExpr[i]>'9') && ultExpr[i]!='.')
				{
					raizPars->mytype=LMP_variable;
					break;
				}
			}
			raizPars->useful=true;
			raizPars->str.copySubString(ultExpr.c_str(), raizPars->pre_i1, raizPars->pos_i2);
		}
		else // pre-operador (ejemplo: signo -)
		{
			raizPars->mytype=LMP_preoperator;
			while (raizPars->pre_i2+1<(int)ultExpr.size()-1 && (ultExpr[raizPars->pre_i2+1] == ' ' || ultExpr[raizPars->pre_i2+1] == '\t') )
				raizPars->pre_i2++;
			while (raizPars->pos_i1-1>=0 && (ultExpr[raizPars->pos_i1-1] == ' ' || ultExpr[raizPars->pos_i1-1] == '\t'))
				raizPars->pos_i1--;
			if (raizPars->pre_i2 < raizPars->pos_i1)
			{
				raizPars->useful=true;
				raizPars->str.copySubString(ultExpr.c_str(), raizPars->pre_i2+1, raizPars->pos_i1-1);
			}
		}
	}
	else
	{
		while (raizPars->pre_i2+1<(int)ultExpr.size()-1 && (ultExpr[raizPars->pre_i2+1] == ' ' || ultExpr[raizPars->pre_i2+1] == '\t'))
			raizPars->pre_i2++;
		while (raizPars->pos_i1-1>=0 && (ultExpr[raizPars->pos_i1-1] == ' ' || ultExpr[raizPars->pos_i1-1] == '\t'))
			raizPars->pos_i1--;
		if (raizPars->pre_i2 < raizPars->pos_i1)
		{
			raizPars->useful=true;
			raizPars->str.copySubString(ultExpr.c_str(), raizPars->pre_i2+1, raizPars->pos_i1-1);
		}
		if (raizPars->posArg==NULL) // Post-operador
			raizPars->mytype=LMP_postoperator;
		else // Operador binario o funcion (izq: name, der:argumentos)
		{
			if (raizPars->useful)
				raizPars->mytype=LMP_operator;
			else
				raizPars->mytype=LMP_function;
		}
	}
	return (raizPars->useful || izqUtil || derUtil);
}

const char *L_MatParsNode::strUndef="LMP_undef";
const char *L_MatParsNode::strVariable="LMP_variable";
const char *L_MatParsNode::strConstant="LMP_constant";
const char *L_MatParsNode::strOperator="LMP_operator";
const char *L_MatParsNode::strFunction="LMP_function";
const char *L_MatParsNode::strPreoperator="LMP_preoperator";
const char *L_MatParsNode::strPostoperator="LMP_postoperator";

const char *L_MatParsNode::tipoStr()
{
	switch(mytype)
	{
	case LMP_variable: return strVariable;
	case LMP_constant: return strConstant;
	case LMP_operator: return strOperator;
	case LMP_function: return strFunction;
	case LMP_preoperator: return strPreoperator;
	case LMP_postoperator: return strPostoperator;
	default: break;
	}
	return strUndef;
}

double L_MatParsNode::evaluateRealRec(const L_MatParser &mat)
{
	double d=0;
	switch(mytype)
	{
	case LMP_variable:
		if (mat.searchVariableValue(str.c_str(), d))
			return d;
		else
			throw_L_ArgException_if(true, "L_MatParsNode::evaluateRealRec()");
	case LMP_constant:
		sscanf(str.c_str(),"%lf",&d);
		break;
	case LMP_operator:
		switch(str[0])
		{
		case'+':
			d=preArg->evaluateRealRec(mat) + posArg->evaluateRealRec(mat);
			break;
		case'-':
			d=preArg->evaluateRealRec(mat) - posArg->evaluateRealRec(mat);
			break;
		case'*':
			d=preArg->evaluateRealRec(mat) * posArg->evaluateRealRec(mat);
			break;
		case'/':
		case':':
			d=preArg->evaluateRealRec(mat) / posArg->evaluateRealRec(mat);
			break;
		case'^':   // '\94' = '^' en dos
			d=pow(preArg->evaluateRealRec(mat), posArg->evaluateRealRec(mat));
			break;
		default:
			throw_L_ArgException_if(true, "L_MatParsNode::evaluateRealRec()");
			#ifndef __TURBOC__
			break;
			#endif
		}
		break;
	case LMP_preoperator:
		if (str[0]=='-')
		{
			d=-posArg->evaluateRealRec(mat);
			break;
		}
	default:
		throw_L_ArgException_if(true, "L_MatParsNode::evaluateRealRec()");
		#ifndef __TURBOC__
		break;
		#endif
	}
	return d;

}

double L_MatParser::evaluateReal()
{
	if (root==NULL)
		return 0;
	return root->evaluateRealRec(*this);
}

L_ComplexDouble L_MatParsNode::evaluateComplexRec(const L_MatParser &mat)
{
	L_ComplexDouble d(0,0);
	double num;
	switch(mytype)
	{
	case LMP_variable:
		if (str[0]=='i' || str[0]=='j' || str[0]=='I' || str[0]=='J')
			d.im=1;
		else if (mat.searchVariableValue(str.c_str(), d))
			return d;
		else
			throw_L_ArgException_if(true,"L_MatParsNode::evaluateComplexRec()");
		break;
	case LMP_constant:
		sscanf(str.c_str(),"%lf",&num);
		d.re=num;
		break;
	case LMP_operator:
		switch(str[0])
		{
		case'+':
			d=preArg->evaluateComplexRec(mat);
			d+=posArg->evaluateComplexRec(mat);
			break;
		case'-':
			d=preArg->evaluateComplexRec(mat);
			d-=posArg->evaluateComplexRec(mat);
			break;
		case'*':
			d=preArg->evaluateComplexRec(mat);
			d*=posArg->evaluateComplexRec(mat);
			break;
		case'/':
		case':':
			d=preArg->evaluateComplexRec(mat);
			d/=posArg->evaluateComplexRec(mat);
			break;
		case'^': //'\94' = '^' en dos
			d=preArg->evaluateComplexRec(mat);
			d.calcPow(posArg->evaluateComplexRec(mat));
			break;
		default:
			throw_L_ArgException_if(true, "L_MatParsNode::evaluateComplexRec()");
			#ifndef __TURBOC__
			break;
			#endif
		}
		break;
	case LMP_preoperator:
		if (str[0]=='-')
		{
			d-=posArg->evaluateComplexRec(mat);
			break;
		}
	default:
		throw_L_ArgException_if(true, "L_MatParsNode::evaluateComplexRec()");
		#ifndef __TURBOC__
		break;
		#endif
	}
	return d;
}

L_ComplexDouble L_MatParser::evaluateComplex()
{
	if (root==NULL)
		return L_ComplexDouble(0,0);
	return root->evaluateComplexRec(*this);
}

long L_MatParsNode::evaluateIntegerRec(const L_MatParser &mat)
{
	long d=0;
	switch(mytype)
	{
	case LMP_variable:
		if (mat.searchVariableValue(str.c_str(), d))
			return d;
		else
			throw_L_ArgException_if(true, "L_MatParsNode::evaluateIntegerRec()");
	case LMP_constant:
		sscanf(str.c_str(),"%ld",&d);
		break;
	case LMP_operator:
		switch(str[0])
		{
		case'+':
			d=preArg->evaluateIntegerRec(mat) + posArg->evaluateIntegerRec(mat);
			break;
		case'-':
			d=preArg->evaluateIntegerRec(mat) - posArg->evaluateIntegerRec(mat);
			break;
		case'*':
			d=preArg->evaluateIntegerRec(mat) * posArg->evaluateIntegerRec(mat);
			break;
		case'/':
		case':':
			d=preArg->evaluateIntegerRec(mat) / posArg->evaluateIntegerRec(mat);
			break;
		case'^': // '\94' = '^' en dos
			d=(long)pow(preArg->evaluateIntegerRec(mat)*1.0, posArg->evaluateIntegerRec(mat)*1.0);
			break;
		default:
			throw_L_ArgException_if(true, "L_MatParsNode::evaluateIntegerRec()");
			#ifndef __TURBOC__
			break;
			#endif
		}
		break;
	case LMP_preoperator:
		if (str[0]=='-')
		{
			d=-posArg->evaluateIntegerRec(mat);
			break;
		}
	default:
		throw_L_ArgException_if(true, "L_MatParsNode::evaluateIntegerRec()");
		#ifndef __TURBOC__
		break;
		#endif
	}
	return d;

}

long L_MatParser::evaluateInteger()
{
	if (root==NULL)
		return 0;
	return root->evaluateIntegerRec(*this);
}


void L_MatParser::showStructure(FILE *fp)
{
	if (root!=NULL)
		showStructureRecursive(fp, root, 0);
}

void L_MatParser::showStructureRecursive(FILE *fp, L_MatParsNode *nodo, int nivel)
{
	int i;
	for (i=0; i<nivel; i++) fprintf(fp, "|");
	if (nodo->preArg!=NULL || nodo->posArg!=NULL)
	{
		fprintf(fp, "%s (%s)\n", nodo->str.c_str(), nodo->tipoStr());
		if (nodo->preArg!=NULL)
			showStructureRecursive(fp, nodo->preArg, nivel+1);
		if (nodo->posArg!=NULL)
			showStructureRecursive(fp, nodo->posArg, nivel+1);
	}
	else
	{
		L_String s;
		fprintf(fp, "%s (%s)\n", nodo->str.c_str(), nodo->tipoStr());
	}
}

int L_MatParser::searchIndexOfPlaceholder(const char *str) const
{
	int i1=0, i2, imed, cmp;
	throw_L_ArgException_if(reemplazos == NULL, "L_MatParser::searchIndexOfPlaceholder()");
	i2=(int)(reemplazos->size()-1);
	while (i1<i2)
	{
		imed=(i1+i2)/2;
		cmp=strcmp(str, (*reemplazos)[imed].varNom.c_str());
		if (cmp>0)
			i1=imed+1;
		else if (cmp<0)
			i2=imed-1;
		else
			break;
	}
	if (strcmp(str, (*reemplazos)[i1].varNom.c_str())==0)
		return i1;
	return -1;
}

const char *L_StringWithNumber::generateString(int i)
{
	L_String formatoPrintf;

	if (numDigits < 0)
	{
		printf("L_StringWithNumber::generateString() : num de digitos no especificado\n");
		*this = prefix + suffix;
		return NULL;
	}

	if (i < 0)
	{
		printf("L_StringWithNumber::generateString() : indice negativo\n");
		*this = prefix + suffix;
		return NULL;
	}

	formatoPrintf = prefix + "%0" + numDigits + "d" + suffix;
	resize(formatoPrintf.size() + numDigits + 10); // faltaba el numDigits, obvio
	sprintf(data(), formatoPrintf.c_str(), i);
	resize((int)strlen(c_str()) + 2); // CUIDADO ACA

	if (i >= numMax)
	{
		printf("L_StringWithNumber::generateString() : indice fuera de rango\n");
		return NULL;
	}

	return c_str();
}

void L_StringWithNumber::write(FILE *fp)
{
	//           p % 0n d
	if (numDigits >= 0)
		fprintf(fp, "%s%%0%dd%s\n", prefix.c_str(), numDigits, suffix.c_str());
	else
		fprintf(fp, "%s\n", c_str());
	fprintf(fp, "%d\n", numMax);
}

bool L_StringWithNumber::read(FILE *fp)
{
	char buf[400];
	int i, len, nleidos;
	if (fgets(buf, 390, fp) == NULL)
		{printf("L_StringWithNumber::read() : formato de archivo incorrecto\n"); return false;}
	len = (int)strlen(buf);

	for (i=0; i<len; i++)
	{
		if (buf[i] == '%')
			break;
	}
	if (i == len) // No es parametrico
	{
		*this = buf;
		erase_spacing_characters_at_beginning_and_at_ending();
		numDigits = -1;
		return false;
	}
	nleidos = sscanf(&buf[i+1], "%d", &numDigits); // Leer digitos
	if (nleidos == 0) // El archivo terminaba en '%'
	{
		*this = buf;
		erase_spacing_characters_at_beginning_and_at_ending();
		numDigits = -1;
		return false;

	}
	buf[i] = 0;  // % = fin del primer string
	i++; // Este es 0
	for (; i<len; i++)
	{
		if (buf[i] == 'd')
			break;
		if (buf[i] < '0' || buf[i] > '9')
		{
			*this = buf;
			erase_spacing_characters_at_beginning_and_at_ending();
			numDigits = -1;
			return false;
		}
	}
	if (i == len)
	{
		prefix = buf;
		suffix = "";
		return false;
	}
	prefix = buf;
	suffix = &buf[i+1]; // despues de la d
	prefix.erase_spacing_characters_at_beginning_and_at_ending();
	suffix.erase_spacing_characters_at_beginning_and_at_ending();
	if (fscanf(fp, "%d", &numMax) == 0)
	{
		printf("L_StringWithNumber::read() : formato de archivo incorrecto\n");
		return false;
	}
	if (fgets(buf, 390, fp) != NULL) // \n
		{printf("L_StringWithNumber::read() : falta \\n final en archivo\n"); return false;}
	return true;
}

void L_StringWithNumber::test()
{
	L_StringWithNumber param;
	char pref[400], suf[400];
	int i;

	printf("Ingrese prefix:");
	if (fgets(pref, 390, stdin) == NULL)
		{printf("L_StringWithNumber::test() : error en stdin\n"); return;}
	param.prefix = pref;
	param.prefix.erase_spacing_characters_at_beginning_and_at_ending();
	printf("Ingrese numero de cifras: ");
	while (scanf("%d", &param.numDigits) == 0)
	{
		if (fgets(suf,390,stdin) == NULL)
			{printf("L_StringWithNumber::test() : error en stdin\n"); return;}
		printf("Ingrese numero de cifras: ");
	}
	if (fgets(suf, 390, stdin) == NULL) // \n
		{printf("L_StringWithNumber::test() : error en stdin\n"); return;}
	printf("Ingrese suffix:");
	if (fgets(suf, 390, stdin) == NULL)
		{printf("L_StringWithNumber::test() : error en stdin\n"); return;}
	param.suffix = suf;
	param.suffix.erase_spacing_characters_at_beginning_and_at_ending();
	param.numMax = 10;

	printf("Probando 0-9:\n");
	for (i=0; i<10; i++)
		puts(param.generateString(i));
	printf("\n");
	param.write(stdout);
	printf("\n");

	printf("Ingrese un string parametrico en el formato mostrado (string, nummax):");
	while (param.read(stdin) == false)
		printf("String parametrico no valido. Ingrese other por favor: ");

	for (i=0; i<param.numMax; i++)
		puts(param.generateString(i));

	printf("Presione ENTER para terminar\n");
	getchar();
}




void L_Spline1D::test()
{
	L_Spline1D spl;
	std::vector<int> x,y;
	int i, nMax=6;

	spl.resize(nMax);
	x.resize(nMax);
	y.resize(nMax);

	for (i=0; i<nMax; i++)
	{
		x[i] = 3*i;
		y[i] = rand()%8;
		spl.x[i] = x[i];
		spl.y[i] = y[i];
	}

	spl.xEntero = false;
	spl.compute();

	for (i=0; i<nMax; i++)
		printf("%.1f\t?\t?\t", (double)y[i]); // value correspondiente a x=2*i

	printf("\n");

	for (i=0; i<3*nMax; i++)
		printf("%.1f\t", spl.evaluate(i));
	printf("\n");

	printf("Presione ENTER para salir\n");
	getchar();
}

void L_Spline1D::pruebaGen()
{
	L_Spline1DGen<double> spl;
	std::vector<int> x,y;
	int i, nMax=6;

	x.resize(nMax);
	y.resize(nMax);

	for (i=0; i<nMax; i++)
	{
		x[i] = 3*i;
		y[i] = rand()%8;
		spl.addFrom(x[i],y[i]);
	}

	spl.compute();

	for (i=0; i<nMax; i++)
		printf("%d ? ", (int)y[i]); // value correspondiente a x=2*i

	printf("\n");

	for (i=0; i<3*nMax; i++)
		printf("%.1f ", spl.evaluate(i));

	printf("Presione ENTER para salir\n");
	getchar();
}

L_BMPheaderInfo::L_BMPheaderInfo()
{
	L_verifTypeSizes();
	_bmp_B='B';
	_bmp_M='M';
	_sizeOfFile=0;
	_reserv1=0;
	_reserv2=0;
	_offsetData=54;
	headerSize=40;
	width=0;
	height=0;
	planes=1;
	bits=24;
	compr=0;
	sizeOfImage=0;
	pixByEachMeterX=196;
	pixByEachMeterY=196;
	nColors=0;
	nColorsImp=0;
	whiteSpace=0;
}

bool L_BMPheaderInfo::setSize(L_uint32 width, L_uint32 Alto, L_uint16 Bits)
{
	long tamMemo;
	long tamLinea;
	bits=Bits;
	width=width;
	height=Alto;

	if (bits!=16 && bits!=24 && bits!=32)
		return false;
	tamLinea=(width*bits+7)/8;
	tamLinea=(tamLinea+3)/4*4;
	sizeOfImage=tamLinea*height;

	tamMemo=width*(bits/8)*height;
	whiteSpace=(sizeOfImage-tamMemo)/height;
	_sizeOfFile=_offsetData+sizeOfImage;
	return true;
}

void L_BMPheaderInfo::write(FILE *fp)
{
	fwrite(&_bmp_B,sizeof(_bmp_B),1,fp);
	fwrite(&_bmp_M,sizeof(_bmp_M),1,fp);
	fwrite(&_sizeOfFile,sizeof(_sizeOfFile),1,fp);
	fwrite(&_reserv1,sizeof(_reserv1),1,fp);
	fwrite(&_reserv2,sizeof(_reserv2),1,fp);
	fwrite(&_offsetData,sizeof(_offsetData),1,fp);
	//
	fwrite(&headerSize,sizeof(headerSize),1,fp);
	fwrite(&width,sizeof(width),1,fp);
	fwrite(&height,sizeof(height),1,fp);
	fwrite(&planes,sizeof(planes),1,fp);
	fwrite(&bits,sizeof(bits),1,fp);
	fwrite(&compr,sizeof(compr),1,fp);
	fwrite(&sizeOfImage,sizeof(sizeOfImage),1,fp);
	fwrite(&pixByEachMeterX,sizeof(pixByEachMeterX),1,fp);
	fwrite(&pixByEachMeterY,sizeof(pixByEachMeterY),1,fp);
	fwrite(&nColors,sizeof(nColors),1,fp);
	fwrite(&nColorsImp,sizeof(nColorsImp),1,fp);
}

bool L_BMPheaderInfo::read(FILE *fp)
{
	long tamMemo;
	if (fread(&_bmp_B,sizeof(_bmp_B),1,fp)==0)
		return false;
	if (_bmp_B!='B')
		return false;
	if (fread(&_bmp_M,sizeof(_bmp_M),1,fp)==0)
		return false;
	if (_bmp_M!='M')
		return false;
	size_t nleidos = 0;
	nleidos += fread(&_sizeOfFile,sizeof(_sizeOfFile),1,fp);
	nleidos += fread(&_reserv1,sizeof(_reserv1),1,fp);
	nleidos += fread(&_reserv2,sizeof(_reserv2),1,fp);
	nleidos += fread(&_offsetData,sizeof(_offsetData),1,fp);
	//
	nleidos += fread(&headerSize,sizeof(headerSize),1,fp);
	nleidos += fread(&width,sizeof(width),1,fp);
	nleidos += fread(&height,sizeof(height),1,fp);
	nleidos += fread(&planes,sizeof(planes),1,fp);
	nleidos += fread(&bits,sizeof(bits),1,fp);
	nleidos += fread(&compr,sizeof(compr),1,fp);
	nleidos += fread(&sizeOfImage,sizeof(sizeOfImage),1,fp);
	nleidos += fread(&pixByEachMeterX,sizeof(pixByEachMeterX),1,fp);
	nleidos += fread(&pixByEachMeterY,sizeof(pixByEachMeterY),1,fp);
	nleidos += fread(&nColors,sizeof(nColors),1,fp);
	nleidos += fread(&nColorsImp,sizeof(nColorsImp),1,fp);
	if (nleidos != 15)
		{printf("L_BMPheaderInfo::read() : formato de archivo incorrecto\n"); return false;}
	if (bits!=16 && bits!=24 && bits!=32)
		return false;
	tamMemo=width*(bits/8)*height;
	if ((L_uint32)sizeOfImage<(L_uint32)tamMemo) // Algunos archivos BMP no definen el campo sizeOfImage...
		whiteSpace=(4-width%4)%4;
	else
		whiteSpace=(sizeOfImage-tamMemo)/height;
	return true;
}


int L_SigmaConvGaus::cmp(const void *a, const void *b)
{
	if ( ((L_SigmaConvGaus*)a)->codeMapping > ((L_SigmaConvGaus*)b)->codeMapping )
		return 1;
	else if ( ((L_SigmaConvGaus*)a)->codeMapping < ((L_SigmaConvGaus*)b)->codeMapping )
		return -1;
	return 0;
}

void L_SigmasConvGausCounter::agregaSigma(double s)
{
	int i;
	long cod=L_ImagenBNxfloatFnConvNodo::calcN(s);
	if (nSigmas>=100)
		return;
	for (i=0; i<nSigmas; i++)
	{
		if (arr[i].codeMapping==cod)
			return;
	}
	arr[nSigmas].sigma=s;
	arr[nSigmas].codeMapping=cod;
    nSigmas++;
}

void L_SigmasConvGausCounter::imprCodigo(FILE *fp)
{
	int i;
	if (fp==NULL)
		return;

	qsort(arr, nSigmas, sizeof(L_SigmaConvGaus), L_SigmaConvGaus::cmp);

	fprintf(fp,"/////////// codeMapping para L_Fnes_pato.h -> L_ImageGrayDouble ////////// \n\n");

	fprintf(fp, "\t#define L_defFnBNxfloat(S) L_ImageGrayDouble& L_ImageGrayDouble::z_difumX_##S##"
		"(const L_ImageGrayDouble &im1, int paso);"
		"L_ImageGrayDouble& L_ImageGrayDouble::z_difumY_##S##"
		"(const L_ImageGrayDouble &im1, int paso)\n");

	for (i=0; i<nSigmas; i++)
		fprintf(fp,"\tL_defFnBNxfloat(%ld);\n",arr[i].codeMapping);

	fprintf(fp,"\n\n/////////// codeMapping para L_Fnes_pato.h -> L_ImagenBNxfloatFnConvArreglo::L_ImagenBNxfloatFnConvArreglo() ////////// \n\n");

	fprintf(fp,"\t\t#define L_defAgregaFnBNxfloat(S) ( addFrom(&L_ImageGrayDouble::z_difumX_##S, &L_ImageGrayDouble::z_difumY_##S, S) )\n");

	for (i=0; i<nSigmas; i++)
		fprintf(fp,"\t\tL_defAgregaFnBNxfloat(%ld);\n",arr[i].codeMapping);

	fprintf(fp,"\n\n/////////// codeMapping para L_Fnes_pato_Precalc.cpp ////////// \n\n");

	for (i=0; i<nSigmas; i++)
	{
		L_ImageGrayDouble::z_difumXGenCodigo(arr[i].sigma, fp);
		L_ImageGrayDouble::z_difumYGenCodigo(arr[i].sigma, fp);
	}
}

#ifdef L_USE_SMOOTH_SIGMA_STORING_TREE
L_SigmasConvGausCounter L_ImageGrayDouble::cuentaSigmas; // cuentaSigmas vive en este archivo
#endif

// Senales

bool L_SignalDouble::readWAV(const char *name)
{
	L_Arr_SignalDouble arr;
	if (arr.readWAV(name) == false)
		return false;
	arr[0].swap(*this);
	return true;
}


bool L_SignalDouble::writeWAV(const char *name, int bitsPerSample)
{
	L_Arr_SignalDouble arr;
	bool ret;
	arr.resize(1);
	arr[0].swap(*this);
	ret = arr.writeWAV(name, bitsPerSample);
	arr[0].swap(*this);
	return ret;
}

class L_RiffHeader
{
public:
	char riff[4]; // "RIFF" para wav. "RIFX" para bytes invertidos.
	L_uint32 chunkSize; // El tamano del elemento que sigue hacia abajo
	char wave[4]; // "WAVE" para wav
	char fmt_[4]; // "fmt "
	L_uint32 tam1_16; // 16
	L_uint16 pcm; // 1 para PCM
	L_uint16 numCanales; // 1 o 2
	L_uint32 fs; // 8000, 44100, etc
	L_uint32 byteRate; // = fs * numCanales * bitsPerSample/8
	L_uint16 blockAlign; // = numCanales * bitsPerSample/8
	L_uint16 bitsPerSample; // 8 bits, o 16 bits
	char data[4]; // "data" para PCM
	L_uint32 bytes; // numMuestras * numCanales * bitsPerSample/8 = cantidad que hay que read a continuacion

	bool esCorrecto();
	long numMuestras() {return bytes/numCanales/(bitsPerSample/8);}
	void set(L_uint32 fs, L_uint16 bitsPerSample, L_uint16 numCanales, long numMuestras); // bitsPerSample = 8 o 16
};

bool L_RiffHeader::esCorrecto()
{
	return (
		riff[0] == 'R' &&
		riff[1] == 'I' &&
		riff[2] == 'F' &&
		riff[3] == 'F' &&
		chunkSize >= bytes + 36 && // Deberia ser igual pero puede haber basura
		fmt_[0] == 'f' &&
		fmt_[1] == 'm' &&
		fmt_[2] == 't' &&
		fmt_[3] == ' ' &&
		wave[0] == 'W' &&
		wave[1] == 'A' &&
		wave[2] == 'V' &&
		wave[3] == 'E' &&
		tam1_16 == 16 &&
		pcm == 1 &&
		byteRate == fs * numCanales * bitsPerSample/8 &&
		blockAlign == numCanales * bitsPerSample/8 &&
		(bitsPerSample == 8 || bitsPerSample == 16 || bitsPerSample == 32) &&
		data[0] == 'd' &&
		data[1] == 'a' &&
		data[2] == 't' &&
		data[3] == 'a' &&
		bytes == numMuestras() * numCanales * bitsPerSample/8
		);
}

void L_RiffHeader::set(L_uint32 fs, L_uint16 bitsPerSample, L_uint16 numCanales, long numMuestras)
{
	this->fs = fs;
	this->bitsPerSample = bitsPerSample;
	this->numCanales = numCanales;

	riff[0] = 'R';
	riff[1] = 'I';
	riff[2] = 'F';
	riff[3] = 'F';
	chunkSize = bytes + 36;
	wave[0] = 'W';
	wave[1] = 'A';
	wave[2] = 'V';
	wave[3] = 'E';
	fmt_[0] = 'f';
	fmt_[1] = 'm';
	fmt_[2] = 't';
	fmt_[3] = ' ';
	tam1_16 = 16;
	pcm = 1;
	byteRate = fs * numCanales * bitsPerSample/8;
	blockAlign = numCanales * bitsPerSample/8;
	data[0] = 'd';
	data[1] = 'a';
	data[2] = 't';
	data[3] = 'a';
	bytes = numMuestras * numCanales * bitsPerSample/8;
}

bool L_Arr_SignalDouble::readWAV(const char *name)
{
	L_RiffHeader header;
	L_Array<L_uint8> arr8;
	L_Array<L_int16> arr16;
	L_Array<L_int32> arr32;
	FILE *fp;
	long i, j;
	double max;

	fp = fopen(name, "rb");
	if (fp == NULL)
		return false;
	if (fread(&header, sizeof(header), 1, fp) != 1)
		{fclose(fp); return false;}		
	if (header.esCorrecto() == false)
	{
		fclose(fp);
		return false;
	}
	resize(header.numCanales);

	for (i=0; i<(int)size(); i++)
	{
		(*this)[i].fs = header.fs; // Faltaba esto
		(*this)[i].resize(header.numMuestras());
	}

	if (header.bitsPerSample == 8)
	{
		arr8.resize(header.numMuestras() * header.numCanales);
		if (fread(arr8.data(), 1, header.numMuestras() * header.numCanales, fp) != header.numMuestras() * header.numCanales)
			{fclose(fp); return false;}
		max = pow(2.0, 8);
		for (i=0; i<(int)arr8.size(); i++)
			for (j=0; j<(int)size(); j++)
				(*this)[j][i] = arr8[size()*i+j]*2.0/(max - 1) - 0.5;
	}
	else if (header.bitsPerSample == 16)
	{
		arr16.resize(header.numMuestras() * header.numCanales);
		if (fread(arr16.data(), 2, header.numMuestras() * header.numCanales, fp) != header.numMuestras() * header.numCanales)
			{fclose(fp); return false;}
		max = pow(2.0, 16);
		for (i=0; i<(int)arr16.size(); i++)
			for (j=0; j<(int)size(); j++)
				(*this)[j][i] = arr16[size()*i+j]*2.0/(max - 1);
	}
	else if (header.bitsPerSample == 32) // Poco estandar
	{
		arr32.resize(header.numMuestras() * header.numCanales);
		if (fread(arr32.data(), 4, header.numMuestras() * header.numCanales, fp) != header.numMuestras() * header.numCanales)
			{fclose(fp); return false;}
		max = pow(2.0, 32);
		for (i=0; i<(int)arr32.size(); i++)
			for (j=0; j<(int)size(); j++)
				(*this)[j][i] = arr32[size()*i+j]*2.0/(max - 1);
	}
	else
	{
		fclose(fp);
		return false;
	}
	fclose(fp);
	return true;
}


bool L_Arr_SignalDouble::writeWAV(const char *name, int bitsPerSample)
{
	L_RiffHeader header;
	L_Array<L_uint8> arr8;
	L_Array<L_int16> arr16;
	L_Array<L_int32> arr32;
	FILE *fp;
	long i, j;
	long numM;
	double max, v;

	throw_L_ArgException_if(size()==0, "L_Arr_SignalDouble::writeWAV()");

	fp = fopen(name, "wb");
	if (fp == NULL)
		return false;

	numM = (*this)[0].size();
	for (i=1; i<(int)size(); i++)
		throw_L_ArgException_if((*this)[i].size() != numM, "L_Arr_SignalDouble::writeWAV()");

	header.set((*this)[0].fs, bitsPerSample, size(), (*this)[0].size());
	throw_L_ArgException_if(!header.esCorrecto(), "L_Arr_SignalDouble::writeWAV()");

	fwrite(&header, sizeof(header), 1, fp);

	if (header.bitsPerSample == 8)
	{
		arr8.resize(header.numMuestras() * header.numCanales);
		max = pow(2.0, 8);
		for (i=0; i<(int)arr8.size(); i++)
			for (j=0; j<(int)size(); j++)
			{
				v = (*this)[j][i]; // La senal debe estar entre -1 y 1
				if (v < -1)
					v = -1; // clipping
				if (v > 1)
					v = 1;
				arr8[size()*i+j] = (L_uint8)(0.5*v*(max-1)+max/2);
			}
		fwrite(arr8.data(), 1, header.numMuestras() * header.numCanales, fp);
	}
	else if (header.bitsPerSample == 16)
	{
		arr16.resize(header.numMuestras() * header.numCanales);
		max = pow(2.0, 16);
		for (i=0; i<(int)arr16.size(); i++)
			for (j=0; j<(int)size(); j++)
			{
				v = (*this)[j][i]; // La senal debe estar entre -1 y 1
				if (v < -1)
					v = -1; // clipping
				if (v > 1)
					v = 1;
				arr16[size()*i+j] = (L_uint16)(0.5*v*(max - 1));
			}
		fwrite(arr16.data(), 2, header.numMuestras() * header.numCanales, fp);
	}
	else if (header.bitsPerSample == 32) // Poco estandar
	{
		arr32.resize(header.numMuestras() * header.numCanales);
		max = pow(2.0, 32);
		for (i=0; i<(int)arr32.size(); i++)
			for (j=0; j<(int)size(); j++)
			{
				v = (*this)[j][i]; // La senal debe estar entre -1 y 1
				if (v < -1)
					v = -1; // clipping
				if (v > 1)
					v = 1;
				arr32[size()*i+j] = (L_uint32)(0.5*v*(max - 1));
			}
		fwrite(arr32.data(), 4, header.numMuestras() * header.numCanales, fp);
	}
	else
	{
		fclose(fp);
		return false;
	}
	fclose(fp);
	return true;
}

void L_SignalDouble::resample_simple(const L_SignalDouble &other, long fsNuevo)
{
	int i;
	double dt, odt;
	fs = fsNuevo;
	dt = 1.0/fs;
	odt = 1.0/other.fs;
	resize((long)(other.size()*odt/dt+1)); //other.size()*odt = n*dt
	for (i=0; i<(int)size(); i++)
		(*this)[i] = other[(long)(i*dt/odt)];
}

void L_SignalDouble::fourier(int dir, L_SignalDouble &imag)
{
	int log2x=0, lxx, i;

	// Revisar que el tamano es potencia de 2
	lxx = size();
	while (lxx > 1)
	{
		lxx /= 2;
		log2x++;
	}
	lxx = 1;
	for (i=0; i<log2x; i++)
		lxx *= 2;

	throw_L_ArgException_if(lxx != size(), "L_SignalDouble::fourier()");

	imag.resize(size());

	L_FFT1D(dir, &(s[0]), &(imag.s[0]), log2x);
}

void L_SignalDouble::fourierByWindowing(int dir, int width, L_SignalDouble &imag)
{
	int log2x=0, lxx, i;
	int nt;
	double *este, *other;
	nt = (int)size() / width;
	throw_L_ArgException_if(nt*width != (int)size(), "L_SignalDouble::fourierByWindows()");

	// Revisar que el tamano es potencia de 2
	lxx = width;
	while (lxx > 1)
	{
		lxx /= 2;
		log2x++;
	}
	lxx = 1;
	for (i=0; i<log2x; i++)
		lxx *= 2;

	throw_L_ArgException_if(lxx != width, "L_SignalDouble::fourierByWindows()");

	imag.resize(size());
	
	for (i=0; i<nt; i++)
	{
		este = &(*this)[i*width];
		other = &imag[i*width];
		L_FFT1D(dir, este, other, log2x);
	}
}


bool L_SignalDouble::genGraphic(L_ImageRGBUchar &im, long i1, long i2)
{
#ifdef L_normGGR
#error Error: L_normGGR ya definido
#endif
#define L_normGGR(X) (int)( (ly-1) - (ly-1)*((X)-min)/(max-min) )
	double sum=0, min, max;
	int i, ly=100;
	L_ShapeArray li;
	if (i2==-1)
		i2=size();
	if (i1>=i2)
		return false;
	sum=min=max=(*this)[i1];
	for (i=i1+1; i<i2; i++)
	{
		sum+=(*this)[i];
		if ((*this)[i] < min)
			min=(*this)[i];
		if ((*this)[i] > max)
			max=(*this)[i];
	}
	if (min==max) // senal plana
	{
		min-=1;
		max+=1;
	}
	im.reallocate(i2-i1+1, ly);
	for (i=i1+1; i<i2; i++)
		li.drawLine(i-1, L_normGGR((*this)[i-1]), i, L_normGGR((*this)[i]));
	im.genDrawing(li);
	return true;
#undef L_normGGR
}

bool L_SignalDouble::scanLineFromImage(const L_ImageGrayDouble &im, double i1, double j1, double i2, double j2, bool bilineal)
{
	double length, dx, dy, x, y;
	long u;
	dx=i2-i1;
	dy=j2-j1;
	length=sqrt((i1-i2)*(i1-i2)+(j1-j2)*(j1-j2));
	dx/=length;
	dy/=length;
	resize((int)(length+1));
	x=i1;
	y=j1;
	if (bilineal)
	{
		double sumaArr, sumaAba, distX, distY;
		int iIzq, jArr;
		for (u=0; u<length; u++)
		{
			iIzq=(int)x;
			jArr=(int)y;
			if (iIzq<0 || jArr<0 || iIzq+1>=im.lx || jArr+1>=im.ly)
				return false;
			distX=x-iIzq;
			distY=y-jArr;
			sumaArr=im.pix(iIzq,jArr)*(1-distX)+im.pix(iIzq+1,jArr)*distX;
			sumaAba=im.pix(iIzq,jArr+1)*(1-distX)+im.pix(iIzq+1,jArr+1)*distX;
			(*this)[u]=(double)(sumaArr*(1-distY)+sumaAba*distY);
			x+=dx;
			y+=dy;
		}
	}
	else
	{
		int i, j;
		for (u=0; u<length; u++)
		{
			i=(int)(x+0.5);
			j=(int)(y+0.5);
			if (i<0 || j<0 || i>=im.lx || i>=im.ly)
				return false;
			(*this)[u]=im.pix(i,j);
			x+=dx;
			y+=dy;
		}

	}
	return true;
}

bool L_SignalDouble::countMaxMin(int &nMax, int &nMin, double umbralMax, double umbralMin, int semiAnchoVentana)
{
	int i, j, j1, j2;
	bool esMax, esMin;
	nMax=0;
	nMin=0;
	for (i=0; i<(int)size(); i++)
	{
		j1=L_MAX(i-semiAnchoVentana,0);
		j2=L_MIN(i+semiAnchoVentana, (int)size()-1);
		esMax=true;
		esMin=true;
		for (j=j1; j<i; j++)
		{
			if ((*this)[i] <= (*this)[j])
			{
				esMax=false;
				if (!esMin)
					break;
			}
			if ((*this)[i] >= (*this)[j])
			{
				esMin=false;
				if (!esMax)
					break;
			}
		}
		for (j=i+1; j<j2; j++)
		{
			if ((*this)[i] <= (*this)[j])
			{
				esMax=false;
				if (!esMin)
					break;
			}
			if ((*this)[i] >= (*this)[j])
			{
				esMin=false;
				if (!esMax)
					break;
			}
		}
		if (esMax && (*this)[i]>umbralMax)
			nMax++;
		if (esMin && (*this)[i]<umbralMin)
			nMin++;
	}
	return true;
}

L_SignalDouble& L_SignalDouble::convSim(const L_SignalDouble &other, L_SignalDouble &h, int paso)
{
	long i, k, u, v, s3;
	double sum;

	resize(other.size());
	fs = other.fs;

	s3=h.size()/2;

	for (i=0; i<(int)size(); i++)
	{
		u=i-s3;
		v=i+s3;
		if (u<0)
			u=0;
		if (v>=(int)size())
			v=(int)size()-1;
		sum=0;
		for (k=u; k<=v; k+=paso)
			sum+=other[k]*h[i-k+s3];
		(*this)[i]=(double)sum;
	}
	return *this;
}

bool L_SignalDouble::normalizeHistogram(int borde, double valMin, double valMax)
{
	int i;
	double min, max;
	if (borde >= (int)size())
		return false;
	min=(*this)[borde];
	max=(*this)[borde];
	for (i=borde+1; i<(int)size()-borde; i++)
	{
		if ((*this)[i]>max)
			max=(*this)[i];
		if ((*this)[i]<min)
			min=(*this)[i];
	}
	if (max == min)
	{
		for (i=0; i<(int)size(); i++)
			(*this)[i]=0;
	}
	else
	{
		for (i=0; i<(int)size(); i++)
			(*this)[i]=(double)( ((*this)[i]-min)/(max-min) * (valMax-valMin) + valMin );
	}
	return true;
}

void L_SignalDouble::filter_IR_right(double a)
{
	int i;
	for (i=1; i<(int)size(); i++)
		(*this)[i]=(double)a*(*this)[i]+(double)(1-a)*(*this)[i-1];
}

void L_SignalDouble::filter_IR_left(double a)
{
	int i;
	for (i=size()-2; i>=0; i--)
		(*this)[i]=(double)a*(*this)[i]+(double)(1-a)*(*this)[i+1];
}


int L_ImageGrayDouble::readBMP(const char *name)
{
	int nb;
	L_ImageRGBUchar imRGB;
	L_ImageGrayDouble imBN;
	nb=imRGB.readBMP(name);
	if (nb==0) return 0;
	imBN=imRGB;
	imBN.swap(*this);
	return nb;
}
bool L_ImageGrayDouble::writeBMP(const char *name, short nb)
{
	L_ImageRGBUchar imRGB;
	imRGB=*this;
	return imRGB.writeBMP(name, nb);
}
bool L_ImageGrayDouble::writeBMP(const char *name) {return writeBMP(name,24);}

bool L_ImageGrayDouble::writeBMP_RB(const char *name, int marg)
{
	L_ImageRGBUchar imRGB;
	int i, j;
	double min, max;
	double val;
	imRGB.reallocate(lx, ly);
	min=max=pix(0,0);
	for (j=marg; j<ly-marg; j++)
	{
		for (i=marg; i<lx-marg; i++)
		{
			if (pix(i,j)<min)
				min=pix(i,j);
			if (pix(i,j)>max)
				max=pix(i,j);
		}
	}
	if (max<0)
		max=-max;
	if (min<0)
		min=-min;
	if (min>max)
		max=min;
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			if (pix(i,j)>=0)
			{
				val=pix(i,j)/max*255;
				if (val>255)
					val=255;
				imRGB.pix(i,j,0)=(L_uchar)val;
				imRGB.pix(i,j,1)=0;
				imRGB.pix(i,j,2)=0;
			}
			else
			{
				val=-pix(i,j)/max*255;
				if (val>255)
					val=255;
				imRGB.pix(i,j,0)=0;
				imRGB.pix(i,j,1)=0;
				imRGB.pix(i,j,2)=(L_uchar)val;
			}
		}
	}
	return imRGB.writeBMP(name);
}

bool L_ImageGrayDouble::readImage(const char *name)
{
	bool b;
	L_ImageRGBUchar im;
	b=im.readImage(name);
	if (b==true)
		*this=im;
	return b;
}
bool L_ImageGrayDouble::saveImage(const char *name)
{
	L_ImageRGBUchar im;
	im=*this;
	return im.saveImage(name);
}

#ifdef __COMPAT_ATLIMAGE__
int L_ImageGrayDouble::readImageATL(const char *name)
{
	int func;
	L_ImageRGBUchar imRGB;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	func=imRGB.readImageATL(name);
	if (func==0) return false;
	reallocate(imRGB.lx,imRGB.ly);
	*this=imRGB;
	return true;
}
bool L_ImageGrayDouble::writeImageATL(const char *name)
{
	L_ImageRGBUchar imRGB;
	imRGB=*this;
	return imRGB.writeImageATL(name);
}
#endif


L_ImageGrayDouble& L_ImageGrayDouble::smooth(const L_ImageGrayDouble &im1, double s, int paso) // La funcion que hay que adaptar
{
#define L_SIGMA_CASE(S) case S:	imt.z_difumX_##S(im1,paso); z_difumY_##S(imt,paso);	break;

#ifdef L_USE_SMOOTH_SIGMA_STORING_TREE
	long n;
	const L_ImagenBNxfloatFnConvNodo *fnNodo;
	if (cuentaSigmas.activo)
		cuentaSigmas.agregaSigma(s);
	n=L_ImagenBNxfloatFnConvNodo::calcN(s);
	fnNodo=L_imBNxFnArr.buscaFn(n); // Intenta buscar fn de convolucion precalculada para el s dado
#else
	const L_ImagenBNxfloatFnConvNodo *fnNodo = NULL;
#endif
	if (fnNodo==NULL)
	{
		if (paso==1)
			z_smooth0Tr(im1,s);
		else
			z_smooth0Tr_step(im1,s,paso);
	}
	else
	{
		try // Intentar usar convoluciones precalculadas
		{
			L_ImageGrayDouble imt;
			(imt.*(fnNodo->fX))(im1,paso);
			((*this).*(fnNodo->fY))(imt,paso);
		}
		catch(L_ArgException e) // Gaussiana precalculada muy grande para el tamano de la imagen
		{
			if (paso==1)
				z_smooth0Tr(im1,s);
			else
				z_smooth0Tr_step(im1,s,paso);
		}
	}

	return *this;
#undef L_SIGMA_CASE
}

L_ImageGrayDouble &L_ImageGrayDouble::smooth_sub(L_ImageGrayDouble &imTmp, const L_ImageGrayDouble &im1, double s, int x0, int y0, int rx1, int ry1, int rx2, int ry2)
{
	z_smooth0Tr_sub(imTmp, im1, s, x0, y0, rx1, ry1, rx2, ry2);
	return *this;
}

// Difumina la imagen im1 desde un sigma sInic hasta un sigma sFinal. El paso se compute en funcion de los sigmas.
L_ImageGrayDouble& L_ImageGrayDouble::smooth_changeSigma(const L_ImageGrayDouble &im1, double sInic, double sFinal, double sMinSubm, FILE *fpo)
{
	int paso=1;
	double sPaso;
	double sPru1, sPru2;

	sPaso=sqrt(sFinal*sFinal-sInic*sInic);
	sPru1=sPaso;
	sPru2=sInic;

	while(true)
	{
		sPru1/=2;
		sPru2/=2;
		if (sPru1<=sMinSubm || sPru2<=sMinSubm)
			break;
		paso++;
	}
	if (fpo!=NULL)
	{
		z_difumXGenCodigo(sPaso, -1, fpo);
		z_difumYGenCodigo(sPaso, -1, fpo);
	}

	return smooth(im1,sPaso,paso);
}

/////////
// Convoluciones

//! build codeMapping en C++ para hacer una convolucion gaussian en direccion X con coeficientes explicitos.
void L_ImageGrayDouble::z_difumXGenCodigo(double s, int s3, FILE *fp)
{
	const double sMin=0.5; // minElement value aceptable de sigma para una gaussian
	char pasotxt[100]="{";
	char minitxt[30];
	char lin1[400];
	char lin2[400];
	std::vector<double> gaus;
	double st;
	int npas=1, paso;
	int i,j,l,v;
	if (s3==-1)
		s3=L_GaussianHalfWidth(s);
	throw_L_ArgException_if(fp==NULL, "L_ImageGrayDouble::z_difumXGenCodigo()");

	//Calcular número de pasos a implementar en la función de convolución
	st=s;
	while (true)
	{
		sprintf(minitxt,"%d,",npas);
		strcat(pasotxt,minitxt);

		st/=2;
		if (st<sMin)
			break;
		npas++;
	}
	l=(int)strlen(pasotxt);
	pasotxt[l-1]='}';

	///// Generar funcion para hacer convolucion en eje X
	// Encabezado y declaracion de variables
	fprintf(fp, "//Funcion de convolucion gaussian en eje X precalculada para sigma %.10lf\n",s);
	fprintf(fp, "L_ImageGrayDouble& L_ImageGrayDouble::z_difumX_%ld(const L_ImageGrayDouble &im1, int paso) // paso: %s\n",L_ImagenBNxfloatFnConvNodo::calcN(s),pasotxt);
	fprintf(fp, "{\n");  // ABRE CORCH 1
	fprintf(fp, "\tregister double **p;\n");
	fprintf(fp, "\tregister double **pim;\n");
	fprintf(fp, "\tregister double **pimfin;\n");
	fprintf(fp, "\tregister int j;\n");
	fprintf(fp, "\throw_L_ArgException_if(im1.lx<=%d, \"L_ImageGrayDouble::z_difumX_%ld()\");\n" , 2*s3+1, L_ImagenBNxfloatFnConvNodo::calcN(s));
	fprintf(fp, "\treconstruye(im1.lx,im1.ly);\n");
	fprintf(fp, "\tpimfin=pix+im1.lx-%d;\n",s3+1);

	for (paso=1; paso<=npas; paso++)
	{
		// Calcular gaussian a usar
		L_calcGausHNorm(gaus,s,2*s3+1,paso);

		if (paso==1)
			fprintf(fp,"\tif (paso==%d)\n",paso);
		else
			fprintf(fp,"\telse if (paso==%d)\n",paso);

		fprintf(fp,"\t{\n");  // ABRE CORCH 2
		fprintf(fp,"\t\tfor (j=0; j<im1.ly; j++)\n");
		fprintf(fp,"\t\t{\n");  // ABRE CORCH 3
		fprintf(fp,"\t\t\tp=im1.pix;\n");
		fprintf(fp,"\t\t\tpim=pix;\n");

		for (i=0; i<s3; i++)
		{
			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<i; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t// i=%d\n",i);
			fprintf(fp,"\t\t\t(*(pim++))[j]=%.6f*p[0][j]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\t\tp++;\n");
		}

		// while( 1 vez )
		{
			fprintf(fp,"\t\t\t// ies\n");
			fprintf(fp, "\t\t\tfor (; pim<=pimfin;)\n"); // , s3);
			fprintf(fp,"\t\t\t{\n"); // ABRE CORCH 4

			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t\t(*(pim++))[j]=%.6f*p[0][j]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\t\tp++;\n");
			fprintf(fp,"\t\t\t}\n"); // CIERRA CORCH 4
		}

		for (i=s3-1; i>=0; i--)
		{
			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<i; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d][j]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t// i=lx-%d\n",i+1);
			fprintf(fp,"\t\t\t(*(pim++))[j]=%.6f*p[0][j]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\tp++;\n");
		}
		fprintf(fp,"\t\t}\n"); // CIERRA CORCH 3 for
		fprintf(fp,"\t}\n"); // CIERRA CORCH 2 else if
	}
	fprintf(fp,"\telse\n");
	fprintf(fp,"\t{\n"); // ABRE CORCH 2
	fprintf(fp,"\t\throw_L_ArgException_if(true,\"L_ImageGrayDouble::z_difumX_%ld()\");\n", L_ImagenBNxfloatFnConvNodo::calcN(s));
	fprintf(fp,"\t}\n"); // CIERRA CORCH 2
	fprintf(fp,"\treturn *this;\n");
	fprintf(fp,"}\n"); // CIERRA CORCH 1
	return;
}

//! build codeMapping en C++ para hacer una convolucion gaussian en direccion Y con coeficientes explicitos.
void L_ImageGrayDouble::z_difumYGenCodigo(double s, int s3, FILE *fp)
{
	const double sMin=0.5; // minElement value aceptable de sigma para una gaussian
	char pasotxt[100]="{";
	char minitxt[30];
	char lin1[400];
	char lin2[400];
	std::vector<double> gaus;
	double st;
	int npas=1, paso;
	int i,j,l,v;
	if (s3==-1)
		s3=L_GaussianHalfWidth(s);
	throw_L_ArgException_if(fp==NULL, "L_ImageGrayDouble::z_difumYGenCodigo()");

	//Calcular número de pasos a implementar en la función de convolución
	st=s;
	while (true)
	{
		sprintf(minitxt,"%d,",npas);
		strcat(pasotxt,minitxt);
		st/=2;
		if (st<sMin)
			break;
		npas++;
	}
	l=(int)strlen(pasotxt);
	pasotxt[l-1]='}';

	///// Generar funcion para hacer convolucion en eje X
	// Encabezado y declaracion de variables
	fprintf(fp, "//Funcion de convolucion gaussian en eje Y precalculada para sigma %.10lf\n",s);
	fprintf(fp, "L_ImageGrayDouble& L_ImageGrayDouble::z_difumY_%ld(const L_ImageGrayDouble &im1, int paso) // paso: %s\n",L_ImagenBNxfloatFnConvNodo::calcN(s),pasotxt);
	fprintf(fp, "{\n");  // ABRE CORCH 1
	fprintf(fp, "\tregister double *p;\n");
	fprintf(fp, "\tregister double *pim;\n");
	fprintf(fp, "\tint i;\n");
	fprintf(fp,"\tregister double *pimfin;\n");
	fprintf(fp, "\throw_L_ArgException_if(im1.ly<=%d, \"L_ImageGrayDouble::z_difumY_%ld()\");\n",2*s3+1, L_ImagenBNxfloatFnConvNodo::calcN(s));
	fprintf(fp,"\treconstruye(im1.lx,im1.ly);\n");

	for (paso=1; paso<=npas; paso++)
	{
		// Calcular gaussian a usar
		L_calcGausHNorm(gaus,s,2*s3+1,paso);

		if (paso==1)
			fprintf(fp,"\tif (paso==%d)\n",paso);
		else
			fprintf(fp,"\telse if (paso==%d)\n",paso);

		fprintf(fp,"\t{\n");  //   ABRE CORCH 2
		fprintf(fp,"\t\tfor (i=0; i<im1.lx; i++)\n");
		fprintf(fp,"\t\t{\n"); // ABRE CORCH 3
		fprintf(fp,"\t\t\tp=im1.pix[i];\n");
		fprintf(fp,"\t\t\tpim=pix[i];\n");
		fprintf(fp,"\t\t\tpimfin=pim+im1.ly-%d;\n",s3+1);

		for (i=0; i<s3; i++)
		{
			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<i; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t// j=%d\n",i);
			fprintf(fp,"\t\t\t*(pim++)=%.6f*p[0]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\tp++;\n");
		}

		// while( 1 vez )
		{
			fprintf(fp,"\t\t\t// js\n");
			fprintf(fp, "\t\t\tfor (; pim<=pimfin;)\n"); //  , s3);
			fprintf(fp,"\t\t\t{\n"); //   ABRE CORCH 4

			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t\t*(pim++)=%.6f*p[0]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\t\tp++;\n");
			fprintf(fp,"\t\t\t}\n"); //   CIERRA CORCH 4
		}

		for (i=s3-1; i>=0; i--)
		{
			lin1[0]=0;
			lin2[0]=0;
			for (j=0; j<s3; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3-v],-v);
					strcat(lin1,minitxt);
				}
			}
			for (j=0; j<i; j++)
			{
				v=j+1;
				if (v%paso==0)
				{
					sprintf(minitxt, "%.6f*p[%+d]+",gaus[s3+v],v);
					strcat(lin2,minitxt);
				}
			}
			l=(int)strlen(lin2);
			if (l==0)
			{
				lin2[0]='0';
				lin2[1]=';';
				lin2[2]=0;
			}
			else
				lin2[l-1]=';';
			fprintf(fp,"\t\t\t// j=ly-%d\n",i+1);
			fprintf(fp,"\t\t\t*(pim++)=%.6f*p[0]+\n",gaus[s3]);
			if (*lin1!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin1);
			if (*lin2!=0)
				fprintf(fp,"\t\t\t\t%s\n",lin2);
			fprintf(fp,"\t\t\tp++;\n");
		}
		fprintf(fp,"\t\t}\n"); // CIERRA CORCH 3 for
		fprintf(fp,"\t}\n"); // CIERRA CORCH 2 else if
	}
	fprintf(fp,"\telse\n");
	fprintf(fp,"\t{\n"); // ABRE CORCH 2
	fprintf(fp,"\t\throw_L_ArgException_if(true,\"L_ImageGrayDouble::z_difumY_%ld()\");\n", L_ImagenBNxfloatFnConvNodo::calcN(s));
	fprintf(fp,"\t}\n"); // CIERRA CORCH 2
	fprintf(fp,"\treturn *this;\n");
	fprintf(fp,"}\n"); // CIERRA CORCH 1
	return;
}

L_ImageGrayDouble& L_ImageGrayDouble::z_smooth0Tr(const L_ImageGrayDouble &im1, double s)
{
	std::vector<double> h;
	L_ImageGrayDouble im2;

	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width);
	im2.convXSim(im1, h);
	convYSim(im2, h);
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::z_smooth0Tr_sub(L_ImageGrayDouble &imTmp, const L_ImageGrayDouble &im1, double s, int x0, int y0, int rx1, int ry1, int rx2, int ry2)
{
	std::vector<double> h;

	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width);
	imTmp.convXSim_sub(im1, h, x0, y0, rx1, ry1);
	convYSim_sub(imTmp, h, x0, y0, rx2, ry2);
	return *this;
}


L_ImageGrayDouble& L_ImageGrayDouble::z_smooth1Tr(const L_ImageGrayDouble &im1, double s)
{
	std::vector<double> h;
	L_ImageGrayDouble im2;
	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width);
	L_ImageGrayDouble im0, im3;
	im0.transpOf(im1);
	im3.convYSim(im0, h);
	im2.transpOf(im3);
	convYSim(im2, h);
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::z_smooth2Tr(const L_ImageGrayDouble &im1, double s)
{
	std::vector<double> h;
	L_ImageGrayDouble im0;
	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width);
	im0.convYSimTr(im1, h);
	convYSimTr(im0, h);
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::filterGabor(const L_ImageGrayDouble &other, L_ImageGrayDouble *gaborIm, double s, double lambda, double ang, double diametroEnSigmas)
{
	///// Ejemplo de uso, filtros Gabor:
	// L_ImageGrayDouble im;
	// L_ImageGrayDouble GabRe, GabIm;
	// im.readImageATL("1_1.tif");
	// GabRe.filterGabor(im, &GabIm, 10.0, 1/7.0, 45.0*M_PI/180);
	// GabRe.normalizeHistogram(10);
	// GabIm.normalizeHistogram(10);
	// GabRe.writeImageATL("GabRe.tif");
	// GabIm.writeImageATL("GabIm.tif");
	/////

	std::vector<double> gabor_cos_cos, gabor_cos_sin, gabor_sin_cos, gabor_sin_sin;
	int width;
	double ang_int;
	double f;
	f=1/lambda;
	L_ImageGrayDouble other_cos_cos;
	L_ImageGrayDouble other_sin_cos;
	L_ImageGrayDouble other_cos_cos__cos_sin;
	L_ImageGrayDouble other_sin_cos__sin_sin;
	L_ImageGrayDouble other_sin_cos__cos_sin;
	L_ImageGrayDouble other_cos_cos__sin_sin;

	if (diametroEnSigmas<0)
		diametroEnSigmas=6.0;
	width=2*(int)(diametroEnSigmas*s/2+1)+1;
	ang_int=-ang; // Debido a que uso un sistema de referencia medio raro en mi formato de imagen

	L_calcGaborSeparable_cos_cos(gabor_cos_cos, s, f, ang_int, width);
	L_calcGaborSeparable_cos_sin(gabor_cos_sin, s, f, ang_int, width);
	L_calcGaborSeparable_sin_cos(gabor_sin_cos, s, f, ang_int, width);
	L_calcGaborSeparable_sin_sin(gabor_sin_sin, s, f, ang_int, width);

// #define L_IMBN_Gab_usaTr  //No conviene usar traspuestas

	other_cos_cos.convXSim(other, gabor_cos_cos);
	other_sin_cos.convXSim(other, gabor_sin_cos);
#ifdef L_IMBN_Gab_usaTr
	other_cos_cos__cos_sin.convYSimTr(other_cos_cos, gabor_cos_sin);
	other_sin_cos__sin_sin.convYSimTr(other_sin_cos, gabor_sin_sin);
#else
	other_cos_cos__cos_sin.convYSim(other_cos_cos, gabor_cos_sin);
	other_sin_cos__sin_sin.convYSim(other_sin_cos, gabor_sin_sin);
#endif

	if (gaborIm!=NULL)
	{
#ifdef L_IMBN_Gab_usaTr
		other_sin_cos__cos_sin.convYSimTr(other_sin_cos, gabor_cos_sin);
		other_cos_cos__sin_sin.convYSimTr(other_cos_cos, gabor_sin_sin);
#else
		other_sin_cos__cos_sin.convYSim(other_sin_cos, gabor_cos_sin);
		other_cos_cos__sin_sin.convYSim(other_cos_cos, gabor_sin_sin);
#endif
	}
	this->OP_subtract(other_cos_cos__cos_sin, other_sin_cos__sin_sin);
	if (gaborIm!=NULL)
		gaborIm->OP_add(other_sin_cos__cos_sin, other_cos_cos__sin_sin);

	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::monogenicDecomposition_slow(int semiancho, L_ImageGrayDouble &I_Cos, L_ImageGrayDouble &I_SinCos, L_ImageGrayDouble &I_SinSin)
{
	L_ImageGrayDouble masc_RieszX;
	L_ImageGrayDouble masc_RieszY;
	double x, y, val_abs3, factor;
	int i, j;
	int u, v;

	masc_RieszX.reallocate(semiancho*2+1, semiancho*2+1);
	masc_RieszY.reallocate(semiancho*2+1, semiancho*2+1);

	factor=1.0;

	for (j=0; j<masc_RieszX.ly; j++)
	{
		for (i=0; i<masc_RieszX.lx; i++)
		{
			masc_RieszX.pix(i,j)=0;
			masc_RieszY.pix(i,j)=0;
			for (u=-5; u<=5; u+=2)
			{
				for (v=-5; v<=5; v+=2)
				{
					x=-(i-semiancho + v/10.0)*factor;
					y=+(j-semiancho + v/10.0)*factor;
					val_abs3=sqrt(x*x+y*y);
					val_abs3=val_abs3*val_abs3*val_abs3;
					masc_RieszX.pix(i,j)=(double)( x/(2*M_PI*val_abs3) / 36.0 *factor );
					masc_RieszY.pix(i,j)=(double)( y/(2*M_PI*val_abs3) / 36.0 *factor );
				}
			}
		}
	}

	I_SinCos.reallocate(lx,ly);
	I_SinSin.reallocate(lx,ly);

	if (&I_Cos != this) // Por si acaso
		I_Cos=*this;
	I_SinCos.conv2d(*this, masc_RieszX, semiancho, semiancho);
	I_SinSin.conv2d(*this, masc_RieszY, semiancho, semiancho);

	return *this;
}



L_ImageGrayDouble& L_ImageGrayDouble::convYSimTr(const L_ImageGrayDouble &other, const std::vector<double> &h)
{
	L_ImageGrayDouble tr1, tr2;
	int n = (int)h.size();
	tr1.transpOf(other);
	tr2.convXSim(tr1, h);
	transpOf(tr2);
	return *this;
}
/////////
// Convoluciones con paso > 1
L_ImageGrayDouble& L_ImageGrayDouble::z_smooth0Tr_step(const L_ImageGrayDouble &im1, double s, int paso)
{
	std::vector<double> h;
	L_ImageGrayDouble im2;
	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width, paso);
	im2.convXSim_step(im1, h, paso);
	convYSim_step(im2, h, paso);
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::z_smoothTr_step(const L_ImageGrayDouble &im1, double s, int paso)
{
	std::vector<double> h;
	L_ImageGrayDouble im2;
	int width=L_GaussianWidth(s);

	L_calcGausHNorm(h, s, width);
	L_ImageGrayDouble im0, im3;
	im0.transpOf(im1);
	im3.convYSim_step(im0, h, paso);
	im2.transpOf(im3);
	convYSim_step(im2, h, paso);
	return *this;
}
L_ImageGrayDouble& L_ImageGrayDouble::z_smooth2Tr_step(const L_ImageGrayDouble &im1, double s, int paso)
{
	std::vector<double> h;
	L_ImageGrayDouble im0;
	int width=L_GaussianWidth(s);
	L_calcGausHNorm(h, s, width);
	im0.convYSimTr_step(im1, h, paso);
	convYSimTr_step(im0, h, paso);
	return *this;
}

// Convolución sin submuestreo ultra rápida (?) por P. Loncomilla
L_ImageGrayDouble& L_ImageGrayDouble::convYSimTr_step(const L_ImageGrayDouble &other, const std::vector<double> &h, int paso)
{
	L_ImageGrayDouble tr1, tr2;
	tr1.transpOf(other);
	tr2.convXSim_step(tr1, h, paso);
	transpOf(tr2);
	return *this;
}

L_ImageGrayDouble& L_ImageGrayDouble::harris(const L_ImageGrayDouble &other, double sd, double si, double alfa)
{
	int c, f;
	double det, trace;
	L_ImageGrayDouble im1;
	L_ImageGrayDouble gradx;
	L_ImageGrayDouble grady;
	L_ImageGrayDouble gradxx;
	L_ImageGrayDouble gradxy;
	L_ImageGrayDouble gradyy;
	L_ImageGrayDouble mxx;
	L_ImageGrayDouble mxy;
	L_ImageGrayDouble myy;

	reallocate(other.lx, other.ly);

	if (sd > 0)
		im1.smooth(other,sd);
	else
		im1 = other;
	gradx.gradRight_3(im1);
	grady.gradUp_3(im1);
	gradxx.OP_mult_elementwise(gradx, gradx);
	gradxy.OP_mult_elementwise(gradx, grady);
	gradyy.OP_mult_elementwise(grady, grady);

	mxx.smooth(gradxx, si);
	mxy.smooth(gradxy, si);
	myy.smooth(gradyy, si);

	for (f=0; f<ly; f++)
	{
		for (c=0; c<lx; c++)
		{
			trace=mxx.pix(c,f)+myy.pix(c,f);
			det=mxx.pix(c,f)*myy.pix(c,f)-mxy.pix(c,f)*mxy.pix(c,f);
			pix(c,f)=(double)(det-alfa*trace*trace);
		}
	}
	return *this;
}

void L_ImageGrayDouble::harrisEliminaLineas(double cornernessMin, const L_ImageGrayDouble &other, double *x, double *y, double n, bool *eliminar, double sd, double si, double alfa)
{
	int c, f;
	int iPun;
	int r;
	double det, trace;
	double cornerness = 1.0;
	if (n==0)
		return;
	reallocate(other.lx, other.ly);
	L_ImageGrayDouble im1(lx,ly), gradx(lx,ly), grady(lx,ly);
	L_ImageGrayDouble gradxx(lx,ly), gradxy(lx,ly), gradyy(lx,ly), mxx(lx,ly), mxy(lx,ly), myy(lx,ly), imTmp(lx,ly);

	r = L_GaussianHalfWidth(si);
	for (iPun=0; iPun<n; iPun++)
	{
		c = (int)(x[iPun]+0.5);
		f = (int)(y[iPun]+0.5);

		if (c<r || c>lx-r || f<r || f>ly-r)
		{
			eliminar[iPun] = true;
			continue;
		}

		// Estas funciones requieren un width previo de 6*si+1
		im1.smooth_sub(imTmp, other, sd, c, f, r-1, r-1, 2*(r+1), 2*(r+1));

		gradx.gradRight_3_sub(im1, c, f, r, r);
		grady.gradUp_3_sub(im1, c, f, r, r);

		gradxx.OP_mult_elementwise_sub(gradx, gradx, c, f, r, r);
		gradxy.OP_mult_elementwise_sub(gradx, grady, c, f, r, r);
		gradyy.OP_mult_elementwise_sub(grady, grady, c, f, r, r);

		mxx.smooth_sub(imTmp, gradxx, si, c, f, 0, r, 0, 0);
		mxy.smooth_sub(imTmp, gradxy, si, c, f, 0, r, 0, 0);
		myy.smooth_sub(imTmp, gradyy, si, c, f, 0, r, 0, 0);

		trace=mxx.pix(c,f)+myy.pix(c,f);
		det=mxx.pix(c,f)*myy.pix(c,f)-mxy.pix(c,f)*mxy.pix(c,f);
		cornerness = det-alfa*trace*trace;
		eliminar[iPun] = cornerness < cornernessMin;
	}
}

double L_ImageGrayDouble::harrisIndividual(L_Array<L_ImageGrayDouble> &imArrTmp, const L_ImageGrayDouble &other, double x, double y, double sd, double si, double alfa)
{
	int c, f, r;
	if (imArrTmp.size() == 0)
	{
		imArrTmp.resize(11);
		for (c=0; c<11; c++)
		{
			imArrTmp[c].reallocate(other.lx, other.ly);
			imArrTmp[c].setZero();
		}
	}
	L_ImageGrayDouble &im1=imArrTmp[0], &gradx=imArrTmp[1], &grady=imArrTmp[2];
	L_ImageGrayDouble &gradxx=imArrTmp[3], &gradxy=imArrTmp[4], &gradyy=imArrTmp[5];
	L_ImageGrayDouble &mxx=imArrTmp[6], &mxy=imArrTmp[7], &myy=imArrTmp[8];
	L_ImageGrayDouble &imTmp=imArrTmp[9], &imHarris=imArrTmp[10];
	double trace, det, cornerness;
	c = (int)(x+0.5);
	f = (int)(y+0.5);
	r = L_GaussianHalfWidth(si);

	if (c<r || c>other.lx-r || f<r || f>other.ly-r)
	{
		return -1.0e30;
	}

	// Estas funciones requieren un width previo de 6*si+1
	im1.smooth_sub(imTmp, other, sd, c, f, r-1, r-1, 2*(r+1), 2*(r+1));

	gradx.gradRight_3_sub(im1, c, f, r, r);
	grady.gradUp_3_sub(im1, c, f, r, r);

	gradxx.OP_mult_elementwise_sub(gradx, gradx, c, f, r, r);
	gradxy.OP_mult_elementwise_sub(gradx, grady, c, f, r, r);
	gradyy.OP_mult_elementwise_sub(grady, grady, c, f, r, r);

	mxx.smooth_sub(imTmp, gradxx, si, c, f, 0, r, 0, 0);
	mxy.smooth_sub(imTmp, gradxy, si, c, f, 0, r, 0, 0);
	myy.smooth_sub(imTmp, gradyy, si, c, f, 0, r, 0, 0);

	trace=mxx.pix(c,f)+myy.pix(c,f);
	det=mxx.pix(c,f)*myy.pix(c,f)-mxy.pix(c,f)*mxy.pix(c,f);
	cornerness = det-alfa*trace*trace;
	return cornerness;
}

void L_ImageGrayDouble::halfResolution(const L_ImageGrayDouble &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v)=(L_uchar)(0.25*(1.0*other.pix(i,j)+other.pix(i+1,j)+other.pix(i,j+1)+other.pix(i+1,j+1)));
		}
}

void L_ImageGrayDouble::halfResolution_noResample(const L_ImageGrayDouble &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v)=other.pix(i,j);
		}
}

void L_ImageGrayDouble::calcSpline2D(L_Spline2DSlow &spline) const
{
	int c, f;
	spline.resize(lx,ly);
	spline.fillIndexes();
	for (f=0; f<ly; f++)
	{
		spline.splArr[f].xEntero=true;
		for (c=0; c<lx; c++)
		{
			spline.splArr[f].x[c]=c;
			spline.splArr[f].y[c]=pix(c,f);
		}
	}
	spline.compute();
}

void L_ImageGrayDouble::cambiaTamanoSpline(const L_ImageGrayDouble &im, double factor)
{
	L_Spline2DSlow spline2D;
	int c, f;
	reallocate((int)(im.lx*factor +0.5), (int)(im.ly*factor +0.5));
	im.calcSpline2D(spline2D);
	for (f=0; f<ly; f++)
	{
		for (c=0; c<lx; c++)
		{
			pix(c,f)=(double)spline2D.evaluate(c/factor, f/factor);
		}
	}
}

void L_ImageGrayDouble::cambiaTamanoPixelado(const L_ImageGrayDouble &im, int nuevolx, int nuevoly)
{
	int c, f;
	double fx, fy;
	reallocate(nuevolx, nuevoly);
	fx = im.lx / (double) nuevolx;
	fy = im.ly / (double) nuevoly;
	for (f=0; f<nuevoly; f++)
		for (c=0; c<nuevolx; c++)
			pix(c,f)=im.pix((int)(c*fx),(int)(f*fy));
}

bool L_ImageGrayDouble::AplicaTransfProyectivaPixelado(const L_ImageGrayDouble &im, const double *coef_8)
{
	int c, f;
	double iOrig, jOrig;
	L_Matrix M;
	double coef[8], den;
	M.reallocate(3,3);
	for (c=0; c<8; c++)
		M(c/3,c%3) = coef_8[c];
	M(2,2) = 1;
	if (M.invertMe() == false)
		return false;
	for (c=0; c<8; c++)
		coef[c]=M(c/3,c%3)/M(2,2);
	for (f=0; f<ly; f++)
	{
		for (c=0; c<lx; c++)
		{
			den = coef[6]*c + coef[7]*f +1;
			iOrig=(coef[0]*c + coef[1]*f + coef[2])/den;
			jOrig=(coef[3]*c + coef[4]*f + coef[5])/den;
			if (iOrig>=0 && iOrig<im.lx-1 && jOrig>=0 && jOrig<im.ly-1)
			{
				pix(c,f)=im.pix(L_ROUND(iOrig),L_ROUND(jOrig));
			}
			else
				pix(c,f)=0;
		}
	}
	return true;
}

void L_ImageGrayDouble::AplicaTransfAfinPixelado(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano)
{
	double detI, mxxI, mxyI, myxI, myyI, txI, tyI;
	//L_Spline2DSlow spline2D;
	int c, f;
	double iOrig, jOrig;
	// Invertiendo la transformacion...
	detI=1/(mxx*myy-mxy*myx);
	mxxI=myy*detI;
	mxyI=-mxy*detI;
	myxI=-myx*detI;
	myyI=mxx*detI;
	txI=(mxy*ty-tx*myy)*detI;
    tyI=-(mxx*ty-tx*myx)*detI;

	if (calcTamano)
	{
		lx=0; ly=0;
		c=(int)(0*mxx+0*mxy+tx+0.5);
		f=(int)(0*myx+0*myy+ty+0.5);
		if (c>lx)
			lx=c;
		if (f>ly)
			ly=f;
		c=(int)(im.lx*mxx+0*mxy+tx+0.5);
		f=(int)(im.lx*myx+0*myy+ty+0.5);
		if (c>lx)
			lx=c;
		if (f>ly)
			ly=f;
		c=(int)(0*mxx+im.ly*mxy+tx+0.5);
		f=(int)(0*myx+im.ly*myy+ty+0.5);
		if (c>lx)
			lx=c;
		if (f>ly)
			ly=f;
		c=(int)(im.lx*mxx+im.ly*mxy+tx+0.5);
		f=(int)(im.ly*myx+im.ly*myy+ty+0.5);
		if (c>lx)
			lx=c;
		if (f>ly)
			ly=f;
		reallocate(lx, ly);
	}
	for (f=0; f<ly; f++)
	{
		for (c=0; c<lx; c++)
		{
			iOrig=mxxI*c+mxyI*f+txI;
			jOrig=myxI*c+myyI*f+tyI;
			if (iOrig>=0 && iOrig<im.lx-1 && jOrig>=0 && jOrig<im.ly-1)
			{
				pix(c,f)=im.pix(L_ROUND(iOrig),L_ROUND(jOrig));
			}
			else
				pix(c,f)=0;
		}
	}
}

void L_ImageGrayDouble::AplicaTransfSemejanzaPixelado(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano)
{
	return AplicaTransfAfinPixelado(im, e*cos(rot), e*sin(rot), -e*sin(rot), e*cos(rot),
		xFijoSal-(e*cos(rot)*xFijoOrig+e*sin(rot)*yFijoOrig), yFijoSal-(-e*sin(rot)*xFijoOrig+e*cos(rot)*yFijoOrig), calcTamano);
}

void L_ImageGrayDouble::AplicaTransfAfinBilineal(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano)
{
	double detI, mxxI, mxyI, myxI, myyI, txI, tyI, a, b, c;
	L_Spline2DSlow spline2D;
	int i, j, u, v;
	double iOrig, jOrig;
	detI=1/(mxx*myy-mxy*myx);
	mxxI=myy*detI;
	mxyI=-mxy*detI;
	myxI=-myx*detI;
	myyI=mxx*detI;
	txI=(mxy*ty-tx*myy)*detI;
    tyI=-(mxx*ty-tx*myx)*detI;

	if (calcTamano)
	{
		lx=0; ly=0;
		i=(int)(0*mxx+0*mxy+tx+0.5);
		j=(int)(0*myx+0*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(im.lx*mxx+0*mxy+tx+0.5);
		j=(int)(im.lx*myx+0*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(0*mxx+im.ly*mxy+tx+0.5);
		j=(int)(0*myx+im.ly*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(im.lx*mxx+im.ly*mxy+tx+0.5);
		j=(int)(im.ly*myx+im.ly*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		reallocate(lx, ly);
	}
	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			iOrig=mxxI*i+mxyI*j+txI;
			jOrig=myxI*i+myyI*j+tyI;
			if (iOrig>=0 && iOrig<im.lx-1 && jOrig>=0 && jOrig<im.ly-1)
			{
				u=(int)iOrig;
				v=(int)jOrig;
				a=iOrig-u;
				b=(1-a)*im.pix(u,v)+a*im.pix(u+1,v);
				c=(1-a)*im.pix(u,v+1)+a*im.pix(u+1,v+1);
				a=jOrig-v;
				pix(i,j)=(double)( (1-a)*b+a*c );
			}
			else
				pix(i,j)=0;
		}
	}
}

void L_ImageGrayDouble::AplicaTransfSemejanzaBilineal(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano)
{
	return AplicaTransfAfinBilineal(im, e*cos(rot), e*sin(rot), -e*sin(rot), e*cos(rot),
		xFijoSal-(e*cos(rot)*xFijoOrig+e*sin(rot)*yFijoOrig), yFijoSal-(-e*sin(rot)*xFijoOrig+e*cos(rot)*yFijoOrig), calcTamano);
}

void L_ImageGrayDouble::AplicaTransfAfinSpline(const L_ImageGrayDouble &im, double mxx, double mxy, double myx, double myy, double tx, double ty, bool calcTamano)
{
	double detI, mxxI, mxyI, myxI, myyI, txI, tyI;
	L_Spline2DSlow spline2D;
	int i, j;
	double iOrig, jOrig;
	detI=1/(mxx*myy-mxy*myx);
	mxxI=myy*detI;
	mxyI=-mxy*detI;
	myxI=-myx*detI;
	myyI=mxx*detI;
	txI=(mxy*ty-tx*myy)*detI;
    tyI=-(mxx*ty-tx*myx)*detI;

	if (calcTamano)
	{
		lx=0; ly=0;
		i=(int)(0*mxx+0*mxy+tx+0.5);
		j=(int)(0*myx+0*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(im.lx*mxx+0*mxy+tx+0.5);
		j=(int)(im.lx*myx+0*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(0*mxx+im.ly*mxy+tx+0.5);
		j=(int)(0*myx+im.ly*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		i=(int)(im.lx*mxx+im.ly*mxy+tx+0.5);
		j=(int)(im.ly*myx+im.ly*myy+ty+0.5);
		if (i>lx)
			lx=i;
		if (j>ly)
			ly=j;
		reallocate(lx, ly);
	}
	im.calcSpline2D(spline2D);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			iOrig=mxxI*i+mxyI*j+txI;
			jOrig=myxI*i+myyI*j+tyI;
			//if (iOrig>=0 && iOrig<=im.lx-1 && jOrig>=0 && jOrig<=im.ly-1)
				pix(i,j)=(double)spline2D.evaluate(iOrig, jOrig);
			//else
				//pix[i][j]=0;
		}
	}
}

void L_ImageGrayDouble::AplicaTransfSemejanzaSpline(const L_ImageGrayDouble &im, double e, double rot, double xFijoOrig, double yFijoOrig, double xFijoSal, double yFijoSal, bool calcTamano)
{
	AplicaTransfAfinSpline(im, e*cos(rot), e*sin(rot), -e*sin(rot), e*cos(rot),
		xFijoSal-(e*cos(rot)*xFijoOrig+e*sin(rot)*yFijoOrig), yFijoSal-(-e*sin(rot)*xFijoOrig+e*cos(rot)*yFijoOrig), calcTamano);
}


#ifdef __COMPAT_ATLIMAGE__
void L_ImageGrayDouble::operator=(const ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=other;
	*this=im;
}
void L_ImageGrayDouble::copyTo(ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=*this;
	im.copyTo(other);
}
#endif

#ifdef __COMPAT_IPLIMAGE__
void L_ImageGrayDouble::operator=(const IplImage &imagen1)
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
					pix(j,i_)= (pixel[0]+(double)pixel[1]+pixel[2]) / (3.0 * 255.0); // 0
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
					pix(j,i)=(pixel[0]+(double)pixel[1]+pixel[2]) / (3.0 * 255.0); // 0
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
					pix(j,i_)=pixel[0] / 255.0;
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
					pix(j,i)=pixel[0] / 255.0;
				}
			}
		}
		return;
	}
	#undef i_
	return;
}

void L_ImageGrayDouble::copyTo(IplImage **imagen2)
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	unsigned char *pixel;
	int i, j;

	if ( (*imagen2) == NULL || (*imagen2)->width!=lx || (*imagen2)->height!=ly )
		(*imagen2) = cvCreateImage(cvSize(lx,ly),IPL_DEPTH_8U,1);

	int origin = (*imagen2)->origin;

	#if defined(i_)
		#error #define conflicts
	#endif
	#define i_ (ly-i-1)
	if((*imagen2)->nChannels == 3)  // RGB
	{
		if (origin==IPL_ORIGIN_BL) // Eje y "invertido" respecto a L_ImageGrayDouble
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j*3];
					pixel[0]=(unsigned char)(pix(j,i_)*255); // 0
					pixel[1]=(unsigned char)(pix(j,i_)*255); // 1
					pixel[2]=(unsigned char)(pix(j,i_)*255); // 2
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
					pixel[0]=(unsigned char)(pix(j,i)*255); // 0
					pixel[1]=(unsigned char)(pix(j,i)*255); // 1
					pixel[2]=(unsigned char)(pix(j,i)*255); // 2
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
					pixel[0]=(unsigned char)(pix(j,i_)*255);
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
					pixel[0]=(unsigned char)(pix(j,i)*255);
				}
			}
		}
		return;
	}
	#undef i_
	return;
}
#endif


#ifdef __COMPAT_CVMAT__
void L_ImageGrayDouble::operator=(const cv::Mat& other)
{
	int i, j;

	reallocate(other.cols, other.rows);

	if (other.type() == CV_32FC1)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = other.at<float>(i, j) / 255.0;
			}
		}
	}
	else if (other.type() == CV_32FC3)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = (other.at<cv::Vec3f>(i, j)[0] + other.at<cv::Vec3f>(i, j)[1] + other.at<cv::Vec3f>(i, j)[2]) / 3 / 255.0;
			}
		}
	}
	else if (other.type() == CV_8UC1)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = other.at<unsigned char>(i, j) / 255.0;
			}
		}
	}
	else if (other.type() == CV_8UC3)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = (other.at<cv::Vec3b>(i, j)[0]/3 + other.at<cv::Vec3b>(i, j)[1]/3 + other.at<cv::Vec3b>(i, j)[2]/3) / 255.0;
			}
		}
	}
	else
	{
		printf("L_ImageGrayDouble::operator=(const cv::Mat& other): bad image format");
	}

	return;
}

void L_ImageGrayDouble::copyTo(cv::Mat& other)
{
	other = cv::Mat::zeros(rows(), cols(), CV_32FC1);
	for (int i = 0; i < other.rows; i++)
	{
		for (int j = 0; j < other.cols; j++)
		{
			other.at<float>(i,j) = pix(j,i) * 255.0;
		}
	}
}

#endif // __COMPAT_CVMAT__


#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageGrayDouble::crearEigen() const
{
	if (lj*sizeof(double) == ljStep)
	{
		size_t lI=(size_t)li, lJ=(size_t)lj;
		Eigen::MatrixXd mat(lI, lJ);
		memcpy(mat.data(), data(), sizeof(double)*li*lj);
		return mat;
	}
	else
		return crearEigenLenta();
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageGrayDouble::crearEigenLenta() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			mat(i,j) = operator()(i,j);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageGrayDouble::crearEigenRef() const
{
	if (lj*sizeof(double) == ljStep)
		throw L_ArgException();
	size_t lI=(size_t)li, lJ=(size_t)lj;
#if EIGEN_VERSION_AT_LEAST(2,90,0)
	Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor> > mat(const_cast<double*>(data()), lI, lJ); // 1 = row major, revisar esto
#else
	Eigen::Map<Eigen::Matrix<double,10000,10000,Eigen::RowMajor> > mat(&elem[0][0], lI, lJ); // 1 = row major
#endif
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
void L_ImageGrayDouble::operator=(const Eigen::MatrixXd &other)
{
	Eigen::MatrixXd m;
	reallocate_fc((int)other.rows(), (int)other.cols());
	//memcpy(&elem[0][0],  Eigen::Transpose<Eigen::MatrixXd>(other).data(), sizeof(double)*li*lj);
	m = Eigen::Transpose<Eigen::MatrixXd>(const_cast<Eigen::MatrixXd&>(other));
	if (lj*sizeof(double) == ljStep)
		memcpy(data(),  m.data(), sizeof(double)*li*lj);
	else
	{
		int i, j;
		for (i=0; i<li; i++)
			for (j=0; j<lj; j++)
				operator()(i,j) = other(i,j);
	}
}
#endif // __COMPAT_EIGEN__



L_ImageGrayUchar& L_ImageGrayUchar::genDibujo_CenGrav(const L_ShapeArray &lins, double *cx, double *cy)
{
	int a;
	int u, u1, u2, v1, uc;
	int dx, dy;
	int dxA, dyA;
	double m;
	long sumx_c=0, sumy_c=0, n_c=0;

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

	for (a=0; a<(int)lins.size(); a++)
	{
		dx=lins[a].x2-lins[a].x1;
		dy=lins[a].y2-lins[a].y1;
		sumx_c+=lins[a].x1;
		sumy_c+=lins[a].y1;
		n_c++;

		if (lins[a].rx!=0 || lins[a].ry!=0)
		{
			// Dibujar círculo o elipse ... wa!!!
			// La sum de las distancias de los ptos de la elipse a los focos es 2R
			if (lins[a].rx == lins[a].ry) // circulo
			{
				double x, y;
				int xi, yi;
				u1=(int)(lins[a].x1-lins[a].rx)-1;
				u2=(int)(lins[a].x1+lins[a].rx)+1;
				uc=lins[a].x1;
				v1=lins[a].y1;
				for (u=u1; u<=u2; u++)
				{
					if (u<0)
						u=0;
					if (u>=lx)
						break;
					y=lins[a].rx*lins[a].rx-(u-uc)*(u-uc);
					if (y>0)
					{
						yi=(int)sqrt(y);
						if (v1+yi>=0 && v1+yi<ly)
						{
							pix(u,v1+yi)=lins[a].R;
						}
						if (v1-yi>=0 && v1-yi<ly)
						{
							pix(u,v1-yi)=lins[a].R;
						}
					}
				}
				u1=(int)(lins[a].y1-lins[a].rx)-1;
				u2=(int)(lins[a].y1+lins[a].rx)+1;
				uc=lins[a].y1;
				v1=lins[a].x1;
				for (u=u1; u<=u2; u++)
				{
					if (u<0)
						u=0;
					if (u>=ly)
						break;
					x=lins[a].rx*lins[a].rx-(u-uc)*(u-uc);
					if (x>0)
					{
						xi=(int)sqrt(x);
						if (v1+xi>=0 && v1+xi<lx)
						{
							pix(v1+xi,u)=lins[a].R;
						}
						if (v1-xi>=0 && v1-xi<lx)
						{
							pix(v1-xi,u)=lins[a].R;
						}
					}
				}
			}
			else // elipse
			{
				double maxr = (lins[a].rx > lins[a].ry) ? lins[a].rx : lins[a].ry;
				double minr = (lins[a].rx < lins[a].ry) ? lins[a].rx : lins[a].ry;
				int numPuntos = (int)(2.0*M_PI*maxr + 5);
				 
				// Ver si es necesario dibujar la elipse
				if (minr*minr > lx*lx + ly*ly || minr < 0 || maxr < 0)
					continue;

				double beta = -lins[a].ang;
				double sinbeta = sin(beta);
				double cosbeta = cos(beta);

				for (double i = 0; i < 360; i += 360.0 / numPuntos) 
				{
					double alfa = i * (M_PI / 180) ;
					double sinalfa = sin(alfa);
					double cosalfa = cos(alfa);

					int X = (int)( lins[a].x1 + (lins[a].rx * cosalfa * cosbeta - lins[a].ry * sinalfa * sinbeta) );
					int Y = (int)( lins[a].y1 + (lins[a].rx * cosalfa * sinbeta + lins[a].ry * sinalfa * cosbeta) );

					if (X>=0 && Y>=0 && X<lx && Y<ly)
					{
						pix(X,Y) = lins[a].R;
					}
				}
			}
		}
		else
		{
			dx=lins[a].x2-lins[a].x1;
			dy=lins[a].y2-lins[a].y1;
			if (dx==0 && dy==0)
			{
				pix(lins[a].x1,lins[a].y1)=lins[a].R;
				continue;
			}
			if (dx>=0)
				dxA=dx;
			else
				dxA=-dx;
			if (dy>=0)
				dyA=dy;
			else
				dyA=-dy;
			if (dxA>=dyA)
			{
				int y;
				m=dy/(double)dx;
				if (dx>0)
				{
					u1=lins[a].x1;
					u2=lins[a].x2;
					v1=lins[a].y1;
				}
				else
				{
					u1=lins[a].x2;
					u2=lins[a].x1;
					v1=lins[a].y2;
				}
				for (u=u1; u<=u2; u++)
				{
					if (u<0)
						u=0;
					if (u>=lx)
						break;
					y=(int)(v1+(u-u1)*m);
					if (y>=0 && y<ly)
					{
						pix(u,y)=lins[a].R;
					}
				}
			}
			else
			{
				int x;
				m=dx/(double)dy;
				if (dy>0)
				{
					u1=lins[a].y1;
					u2=lins[a].y2;
					v1=lins[a].x1;
				}
				else
				{
					u1=lins[a].y2;
					u2=lins[a].y1;
					v1=lins[a].x2;
				}
				for (u=u1; u<=u2; u++)
				{
					if (u<0)
						u=0;
					if (u>=ly)
						break;
					x=(int)(v1+(u-u1)*m);
					if (x>=0 && x<lx)
					{
						pix(x,u)=lins[a].R;
					}
				}
			}
		}
	}
	if (cx!=NULL && cy!=NULL)
	{
		*cx = sumx_c / (double) n_c;
		*cy = sumy_c / (double) n_c;
	}
	return *this;
}


void L_ImageGrayUchar::lbp(L_ImageGrayUchar &orig)
{
	int i, j;
	reallocate(orig.lx, orig.ly);
	for (i=0; i<lx; i++)
		pix(i,0) = pix(i,ly-1) = 0;
	for (j=0; j<ly; j++)
		pix(0,j) = pix(lx-1,j) = 0;
	for (j=1; j<ly-1; j++)
		for (i=1; i<lx-1; i++)
			pix(i,j) =
			((pix(i-1,j-1) > pix(i,j)))     + ((pix(i,j-1) > pix(i,j))<<1)  + ((pix(i+1,j-1) > pix(i,j))<<2) +
			((pix(i-1,j)   > pix(i,j))<<3)  +                                 + ((pix(i+1,j)   > pix(i,j))<<4)+
			((pix(i-1,j+1) > pix(i,j))<<5)  + ((pix(i,j+1) > pix(i,j))<<6)  + ((pix(i+1,j+1) > pix(i,j))<<7);
}


void L_ImageGrayUchar::lbp(L_ImageGrayDouble &orig)
{
	int i, j;
	reallocate(orig.lx, orig.ly);
	for (i=0; i<lx; i++)
		pix(i,0) = pix(i,ly-1) = 0;
	for (j=0; j<ly; j++)
		pix(0,j) = pix(lx-1,j) = 0;
	for (j=1; j<ly-1; j++)
		for (i=1; i<lx-1; i++)
			pix(i,j) =
			((pix(i-1,j-1) > pix(i,j)))     + ((pix(i,j-1) > pix(i,j))<<1)  + ((pix(i+1,j-1) > pix(i,j))<<2) +
			((pix(i-1,j)   > pix(i,j))<<3)  +                                 + ((pix(i+1,j)   > pix(i,j))<<4)+
			((pix(i-1,j+1) > pix(i,j))<<5)  + ((pix(i,j+1) > pix(i,j))<<6)  + ((pix(i+1,j+1) > pix(i,j))<<7);
}




int L_ImageGrayUchar::readBMP(const char *name)
{
	int func;
	L_ImageRGBUchar imRGB;
	func=imRGB.readBMP(name);
	if (func==0) return false;
	reallocate(imRGB.lx,imRGB.ly);
	*this=imRGB;
	return func;
}
bool L_ImageGrayUchar::writeBMP(const char *name, short nb)
{
	L_ImageRGBUchar imRGB;
	imRGB=*this;
	return imRGB.writeBMP(name, nb);
}
bool L_ImageGrayUchar::writeBMP(const char *name) {return writeBMP(name,24);}

bool L_ImageGrayUchar::readImage(const char *name)
{
	bool b;
	L_ImageRGBUchar im;
	if ((b=im.readImage(name))==true)
		*this=im;
	return b;
}
bool L_ImageGrayUchar::saveImage(const char *name)
{
	L_ImageRGBUchar im;
	im=*this;
	return im.saveImage(name);
}


#ifdef __COMPAT_ATLIMAGE__
int L_ImageGrayUchar::readImageATL(const char *name)
{
	int func;
	L_ImageRGBUchar imRGB;
	func=imRGB.readImageATL(name);
	if (func==0) return false;
	reallocate(imRGB.lx,imRGB.ly);
	*this=imRGB;
	return true;
}
bool L_ImageGrayUchar::writeImageATL(const char *name)
{
	L_ImageRGBUchar imRGB;
	imRGB=*this;
	return imRGB.writeImageATL(name);
}
#endif

bool L_ImageGrayUchar::leeB8(const char *name)
{
	FILE *fp;
	int i, j;
	int c;
	fp=fopen(name,"rb");
	if (fp==NULL)
		return false;
	lx=(unsigned char)getc(fp);
	lx=lx*256+(unsigned char)getc(fp);
	ly=(unsigned char)getc(fp);
	ly=ly*256+(unsigned char)getc(fp);
	reallocate(lx, ly);
	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			if ((c=getc(fp))==EOF) {fclose(fp);return false;}
			pix(i,j)=(L_uchar)c;
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUchar::grabaB8(const char *name)
{
	FILE *fp;
	int i, j;
	fp=fopen(name,"wb");
	if (fp==NULL)
		return false;
	putc(lx>>8, fp);
	putc(lx&0xFF, fp);
	putc(ly>>8, fp);
	putc(ly&0xFF, fp);
	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			putc(pix(i,j),fp);
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUchar::grabaPGM(const char *name)
{
	FILE *fp;
	int i, j;
	fp=fopen(name,"wb");
	if (fp==NULL)
		return false;
	fprintf(fp, "P5 ");
	fprintf(fp, "%d ", lx);
	fprintf(fp, "%d ", ly);
	fprintf(fp, "255 ");
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			putc(pix(i,j),fp);
		}
	}
	fclose(fp);
	return true;
}

int L_ImageGrayUchar::grabaComprGris(FILE *fp, double *le, int forzMetodo, double t, long ncam) const
{
	std::vector<char> buf; // No tiene memoria propia, solo se usa como referencia
	std::vector<char> bufD;
	std::vector<char> bufDD;
	std::vector<char> bufRLE;
	std::vector<char> bufT;
	double largoEst[5]; // 0 = sin compr;  1 = compr directo;  2 = delta y compr
	int i, iMin, sz;
	long lxs, lys;
	if (le == NULL)
		le = largoEst;

	// 0: sin comprimir, 1: Huffman, 2:resta-huffman, 3:resta-RLE_huffman
	// Calculos compresion Hough simple
	switch(forzMetodo)
	{
	case -1: // Todos
		//buf.setData((char *)&pix(0,0), lx*ly);  // Better slow programming than bad one...
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0), lx*ly);
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		L_HuffmanEncoder::encodeRLE0(bufDD, bufRLE);
		break;
	case 0:
		//buf.setData((char *)&pix(0,0), lx*ly);
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0), lx*ly);
		break;
	case 1:
		//buf.setData( (char *)&pix(0,0), lx*ly);
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0), lx*ly);
		break;
	case 2:
		encodeDelta(bufD);
		break;
	case 3:
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		break;
	case 4:
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		L_HuffmanEncoder::encodeRLE0(bufDD, bufRLE);
		break;
	default:
		//if (buf.size()!=0)
		//	buf.setData( NULL, 0 );
		printf("Modo de compresion no disponible\n");
		return -1;
	}

	// Evaluar el mas conveniente
	if (forzMetodo == -1)
	{
		le[0] = buf.size();  // Largo sin comprimir
		le[1] = 165 + L_HuffmanEncoder::entropy_in_bits(buf)/8*buf.size(); // Largo aproximado del encabezado + comprimido
		le[2] = 165 + L_HuffmanEncoder::entropy_in_bits(bufD)/8*bufD.size(); // Largo aproximado del encabezado + comprimido con delta
		le[3] = 165 + L_HuffmanEncoder::entropy_in_bits(bufDD)/8*bufDD.size(); // Largo aproximado del encabezado + comprimido con delta
		le[4] = 165 + L_HuffmanEncoder::entropy_in_bits(bufRLE)/8*bufRLE.size(); // Largo aproximado del encabezado + comprimido con delta

		// El calculo de bufDD es tonto (es restar saltandose 1), mejor probar con algo que elimine ruido como:
		//   buf[i] - buf[i-2]/2 - buf[i-1]/2
		//   buf[i] - (buf[i-2]+2*(int)buf[i-1])/3

		iMin = 0;
		for (i=1; i<5; i++)
			if (le[i] < le[iMin])
				iMin = i;
	}
	else
		iMin = forzMetodo;

	//iMin = 2; // Forzar las cosas

	lxs = lx;
	lys = ly;

	// Escritura 1
	fputc('B',fp);
	fputc('M',fp);
	fputc('L',fp);
	fputc('1',fp);

	fwrite(&lxs, sizeof(lxs), 1, fp);
	fwrite(&lys, sizeof(lys), 1, fp);
	fwrite(&t, sizeof(t), 1, fp);
	fwrite(&ncam, sizeof(ncam), 1, fp);
	fputc(iMin,fp);

	switch(iMin)  // 0 = sin compr;  1 = compr directo;  2 = delta y compr
	{
	case 0:
		sz = (int)buf.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(buf[0]), 1, buf.size(), fp);
		break;
	case 1:
		L_HuffmanEncoder::encodeAll(buf, bufT);
		sz = (size_type)bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 2:
		L_HuffmanEncoder::encodeAll(bufD, bufT);
		sz = (size_type)bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 3:
		L_HuffmanEncoder::encodeAll(bufDD, bufT);
		sz = (size_type)bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 4:
		L_HuffmanEncoder::encodeAll(bufRLE, bufT);
		sz = (size_type)bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	}
	//buf.setData( NULL, 0 );
	return iMin;
}

int L_ImageGrayUchar::leeComprGris(FILE *fp, double *t, long *ncam)
{
	std::vector<char> buf;
	std::vector<char> bufD;
	std::vector<char> bufDD;
	std::vector<char> bufRLE;
	std::vector<char> bufT;
	int iMin, sz;
	long lxs, lys, nc;
	double ti;
	int head[4];

	head[0] = fgetc(fp);
	head[1] = fgetc(fp);
	head[2] = fgetc(fp);
	head[3] = fgetc(fp);

	// Lectura 1
	if ( (head[0] != 'B' || head[1] != 'M' || head[2] != 'Z' || head[3] != '1') &&
		(head[0] != 'B' || head[1] != 'M' || head[2] != 'L' || head[3] != '1'))
		return -1;
	size_t nleidos = 0;	
	nleidos += fread(&lxs, sizeof(lxs), 1, fp);
	nleidos += fread(&lys, sizeof(lys), 1, fp);
	nleidos += fread(&ti, sizeof(ti), 1, fp);
	nleidos += fread(&nc, sizeof(nc), 1, fp);
	if (nleidos != 4)
		return -1;
	iMin = fgetc(fp);
	
	if (iMin == EOF)
		return -1;

	if (t!=NULL)
		*t = ti;
	if (ncam != NULL)
		*ncam = nc;

	reallocate(lxs, lys);
	buf.resize(lx*ly);

	switch(iMin)  // 0 = sin compr;  1 = compr directo;  2 = delta y compr  LECTURA
	{
	case 0:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		buf.resize(sz);
		if (fread(&(buf[0]), 1, buf.size(), fp) != buf.size())
			return -1;
		decodeFromBuffer(buf);
		break;
	case 1:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, buf);
		decodeFromBuffer(buf);
		break;
	case 2:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufD);
		decodeDelta(bufD);
		break;
	case 3:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufDD);
		L_HuffmanEncoder::decodeDelta(bufDD, bufD);
		decodeDelta(bufD);
		break;
	case 4:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufRLE);
		L_HuffmanEncoder::decodeRLE0(bufRLE, bufDD);
		L_HuffmanEncoder::decodeDelta(bufDD, bufD);
		decodeDelta(bufD);
	}
	if (head[2] = 'Z')  // BMZ
		transposeMe();
	return iMin;
}

int L_ImageGrayUchar::leeComprGris_contar(FILE *fp, int *histMetodo)
{
	long pos0, pos, posFin;
	int iMin, num=0;
	long lxs, lys, nc;
	int N;
	size_t nleidos;
	double ti;
	int head[4];

	pos0 = ftell(fp);
	fseek(fp, 0, SEEK_END);
	posFin = ftell(fp);
	fseek(fp, pos0, SEEK_SET);

	while(!feof(fp))
	{
		pos = ftell(fp);
		if (pos == posFin)
			break;

		head[0] = fgetc(fp);
		head[1] = fgetc(fp);
		head[2] = fgetc(fp);
		head[3] = fgetc(fp);

		// Lectura 1
		if ( (head[0] != 'B' || head[1] != 'M' || head[2] != 'Z' || head[3] != '1') &&
			(head[0] != 'B' || head[1] != 'M' || head[2] != 'L' || head[3] != '1'))
			return -1;
		nleidos = 0;
		nleidos += fread(&lxs, sizeof(lxs), 1, fp);
		nleidos += fread(&lys, sizeof(lys), 1, fp);
		nleidos += fread(&ti, sizeof(ti), 1, fp);
		nleidos += fread(&nc, sizeof(nc), 1, fp);
		if (nleidos != 4)
			return -1;
		iMin = getc(fp);
		if (iMin == EOF)
			return -1;
		if (histMetodo != NULL)
			(histMetodo)[iMin] ++;
		if (fread(&N, sizeof(N), 1, fp) != 1) // bufT
			return -1;
		fseek(fp, N, SEEK_CUR);
		num++;
	}
	fseek(fp, pos0, SEEK_SET);
	return num;
}


void L_ImageGrayUchar::pruebaComprGris()
{
	FILE *fp = NULL;
	int u, i, j, val=0;
	L_ImageGrayUchar im1, im2;
	const char name[] = "pruebitaaa.bmz";
	srand((unsigned int)time(NULL));

	for (u=0; u<60; u++)
	{
		im1.reallocate(10+rand()%310, 10+rand()%190);
		//im1.fillRandomValues(10,30);
		for (j=0; j<im1.ly; j++)
		{
			for (i=0; i<im1.lx; i++)
			{
				// metodo 1
				val = rand()%10;
				// metodo 2
				//val = val + rand()%10 - 5;
				im1.pix(i,j) = val;
			}
		}
		fp = fopen(name, "wb");
		if (fp == NULL)
		{
			printf("Fallo la escritura de %s\n", name);
			return;
		}
		printf("metodo %d\n", im1.grabaComprGris(fp));
		fclose(fp);
		fp = fopen(name, "rb");
		if (fp == NULL)
		{
			printf("Fallo la lectura de %s\n", name);
			return;
		}
		im2.leeComprGris(fp);
		fclose(fp);
		for (j=0; j<im1.ly; j++)
			for (i=0; i<im1.lx; i++)
				if (im1.pix(i,j) != im2.pix(i,j))
				{
					printf("pruebaComprGris() fallo\n");
					return;
				}
	}
	printf("pruebaComprGris() funciono adecuadamente\n");
}

#ifdef __COMPAT_ATLIMAGE__
void L_ImageGrayUchar::operator=(const ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=other;
	*this=im;
}
void L_ImageGrayUchar::copyTo(ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=*this;
	im.copyTo(other);
}
#endif

#ifdef __COMPAT_IPLIMAGE__
void L_ImageGrayUchar::operator=(const IplImage &imagen1)
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
					pix(j,i_)=(L_uchar)((pixel[0]+(double)pixel[1]+pixel[2]) / 3.0); // 0
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
					pix(j,i)=(L_uchar)((pixel[0]+(double)pixel[1]+pixel[2]) / 3.0); // 0
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
					pix(j,i_)=(L_uchar)pixel[0];
				}
			}
		}
		else
		{
			if (imagen1.width == imagen1.widthStep)
			{
				memcpy(&pix(0,0), imagen1.imageData, sizeof(char)*imagen1.width*imagen1.height);
				return;
			}
			for (i=0; i<imagen1.height; i++)
			{
				for (j=0; j<imagen1.width; j++)
				{
					pixel   = &((unsigned char*)(imagen1.imageData + imagen1.widthStep*i))[j];
					pix(j,i)=(L_uchar)pixel[0];
				}
			}
		}
		return;
	}
	#undef i_
	return;
}
void L_ImageGrayUchar::copyTo(IplImage **imagen2)
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	unsigned char *pixel;
	int i, j;

	if ( (*imagen2) == NULL || (*imagen2)->width!=lx || (*imagen2)->height!=ly )
		(*imagen2) = cvCreateImage(cvSize(lx,ly),IPL_DEPTH_8U,1);

	int origin = (*imagen2)->origin;

	#if defined(i_)
		#error #define conflicts
	#endif
	#define i_ (ly-i-1)
	if((*imagen2)->nChannels == 3)  // RGB
	{
		if (origin==IPL_ORIGIN_BL) // Eje y "invertido" respecto a L_ImageGrayDouble
		{
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j*3];
					pixel[0]=(unsigned char)pix(j,i_); // 0
					pixel[1]=(unsigned char)pix(j,i_); // 1
					pixel[2]=(unsigned char)pix(j,i_); // 2
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
					pixel[0]=(unsigned char)pix(j,i); // 0
					pixel[1]=(unsigned char)pix(j,i); // 1
					pixel[2]=(unsigned char)pix(j,i); // 2
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
					pixel[0]=(unsigned char)pix(j,i_);
				}
			}
		}
		else
		{
			if ((*imagen2)->widthStep == (*imagen2)->width)
			{
				memcpy((*imagen2)->imageData, &pix(0,0), sizeof(L_uchar)*(*imagen2)->width*(*imagen2)->height);
				return;
			}
			for (i=0; i<(*imagen2)->height; i++)
			{
				for (j=0; j<(*imagen2)->width; j++)
				{
					pixel   = &((unsigned char*)((*imagen2)->imageData + (*imagen2)->widthStep*i))[j];
					pixel[0]=(unsigned char)pix(j,i);
				}
			}
		}
		return;
	}
	#undef i_
	return;
}
#endif


#ifdef __COMPAT_CVMAT__
void L_ImageGrayUchar::operator=(const cv::Mat& other)
{
	int i, j;

	reallocate(other.cols, other.rows);

	if (other.type() == CV_32FC1)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = other.at<float>(i, j);
			}
		}
	}
	else if (other.type() == CV_32FC3)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = (other.at<cv::Vec3f>(i, j)[0] + other.at<cv::Vec3f>(i, j)[1] + other.at<cv::Vec3f>(i, j)[2]) / 3;
			}
		}
	}
	else if (other.type() == CV_8UC1)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = other.at<unsigned char>(i, j);
			}
		}
	}
	else if (other.type() == CV_8UC3)
	{
		for (i = 0; i < other.rows; i++)
		{
			for (j = 0; j < other.cols; j++)
			{
				pix(j, i) = other.at<cv::Vec3b>(i, j)[0] / 3 + other.at<cv::Vec3b>(i, j)[1] / 3 + other.at<cv::Vec3b>(i, j)[2] / 3;
			}
		}
	}
	else
	{
		printf("L_ImageGrayUchar::operator=(const cv::Mat& other): bad image format");
	}

	return;
}
void L_ImageGrayUchar::copyTo(cv::Mat& other)
{
	other = cv::Mat::zeros(rows(), cols(), CV_8UC1);
	for (int i = 0; i < other.rows; i++)
	{
		for (int j = 0; j < other.cols; j++)
		{
			other.at<unsigned char>(i, j) = pix(j, i);
		}
	}
}

#endif // __COMPAT_CVMAT__

void L_ImageGrayUint16::doubleResolution(const L_ImageGrayUint16 &other)
{
	int i, j;
	int u, v;
	reallocate(2*other.lx, 2*other.ly);
	for (j=0; j<other.ly; j++)
	{
		for (i=0; i<other.lx; i++)
		{
			u=2*i;
			v=2*j;
			pix(u,v)=pix(u+1,v)=pix(u,v+1)=pix(u+1,v+1)=other.pix(i,j);
		}
	}
}

void L_ImageGrayUint16::halfResolution(const L_ImageGrayUint16 &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
	{
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v)=(L_uchar)(0.25*(1.0*other.pix(i,j)+other.pix(i+1,j)+other.pix(i,j+1)+other.pix(i+1,j+1)));
		}
	}
}

void L_ImageGrayUint16::halfResolution_noResample(const L_ImageGrayUint16 &other)
{
	int i, j;
	int u, v;
	reallocate(other.lx/2, other.ly/2);
	for (j=0; j<other.ly-1; j+=2)
	{
		for (i=0; i<other.lx-1; i+=2)
		{
			u=(int)(0.5*i);
			v=(int)(0.5*j);
			pix(u,v)=other.pix(i,j);
		}
	}
}

#ifdef __COMPAT_IPLIMAGE__
void L_ImageGrayUint16::copyOn_16(IplImage **imagen2)
{
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	(*imagen2) = cvCreateImage(cvSize(lx,ly),IPL_DEPTH_16U,1);
	copyToBuffer((L_uint16*)(*imagen2)->imageData, (*imagen2)->widthStep);
	return;
}
#endif // __COMPAT_IPLIMAGE__



bool L_ImageGrayUint16::saveImage_16(const char *name)
{
	L_FileName na;
	na.asign(name);
	if (na.hasExtension("txt"))
		return saveImage_16_txt(name);
#ifdef __COMPAT_IPLIMAGE__
	return saveImageIPL_16(name);
#else
	return false;
#endif
}

bool L_ImageGrayUint16::saveImage_16_txt(const char *name)
{
	FILE *fp = fopen(name, "w");
	if (fp == NULL)
		return false;
	// lx   ly   nbits   ncanales    es_punto_flotante
	fprintf(fp, "%d %d %d %d %d\n", lx, ly, 16, 1, 0);
	for (int i=0; i<li; i++)
	{
		for (int j=0; j<lj; j++)
			fprintf(fp, "%d ", operator()(i,j));
		fprintf(fp, "\n");
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUint16::saveImage_16_pgm(const char *name)
{
	FILE *fp = fopen(name, "wb");
	if (fp == NULL)
		return false;
	// lx   ly   val_max  bytes....
	fprintf(fp, "P5\n%d %d\n%d\n", lx, ly, 65535);
	for (int i=0; i<li; i++)
	{
		for (int j=0; j<lj; j++)
		{
			unsigned int v = operator()(i,j);
			fputc((v>>8)&0xFF, fp);
			fputc(v&0xFF, fp);
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUint16::savePoints_xyzrgb(const char *name, const L_ImageRGBUchar &im, double dfx, double dfy, double factor)
{
	double cx, cy, d;
	FILE *fp = fopen(name, "w");
	if (fp == NULL)
		return false;
	cx = (lj+1) / 2.0;
	cy = (li+1) / 2.0;
	for (int i=0; i<li; i++)
	{
		for (int j=0; j<lj; j++)
		{
			d = operator()(i,j) * factor;
			fprintf(fp, "%f %f %f %d %d %d\n", (j-cx)/dfx*d, -(i-cy)/dfy*d, -d, im(i,j,0), im(i,j,1), im(i,j,2)); // right, up, back
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUint16::savePoints_xyz(const char *name, double dfx, double dfy, double factor)
{
	double cx, cy, d;
	FILE *fp = fopen(name, "w");
	if (fp == NULL)
		return false;
	cx = (lj+1) / 2.0;
	cy = (li+1) / 2.0;
	for (int i=0; i<li; i++)
	{
		for (int j=0; j<lj; j++)
		{
			d = operator()(i,j) * factor;
			fprintf(fp, "%lf %lf %lf\n", (j-cx)/dfx*d, -(i-cy)/dfy*d, -d); // right, up, back
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageGrayUint16::readImage_16_txt(const char *name)
{
	int nb = 0, nc=1, es_flotante=false;
	int nread, val;
	char lin[100];
	FILE *fp = fopen(name, "r");
	if (fp == NULL)
		return false;
	if (fgets(lin, 99, fp) == NULL)
	{
		printf("L_ImageGrayUint16::readImage_16_txt() : formato de archivo incorrecto (1)\n");
		fclose(fp);
		return false;
	}
	nread = sscanf(lin, "%d%d%d%d%d", &lx, &ly, &nb, &nc, &es_flotante);
	if (nread < 3 || nb!=16 || nc!=1 || es_flotante != 0)
	{
		printf("L_ImageGrayUint16::readImage_16_txt() : formato de archivo incorrecto (1)\n");
		fclose(fp);
		return false;
	}
	reallocate(lx, ly);
	for (int i=0; i<li; i++)
	{
		for (int j=0; j<lj; j++)
		{
			if (fscanf(fp, "%d", &val) != 1)
			{
				printf("L_ImageGrayUint16::readImage_16_txt() : formato de archivo incorrecto (2)\n");
				fclose(fp);
				return false;
			}
			operator()(i,j) = (unsigned short)val;
		}
	}
	fclose(fp);
	return true;
}

int L_ImageRGBDouble::readBMP(const char *name)
{
	int a;
	L_ImageRGBUchar im;
	a=im.readBMP(name);
	if (a==0)
		return 0;
	*this=im;
	return a;
}
bool L_ImageRGBDouble::writeBMP(const char *name, short nb)
{
	bool a;
	L_ImageRGBUchar im;
	im=*this;
	a=im.writeBMP(name);
	return a;
}

bool L_ImageRGBDouble::readImage(const char *name)
{
	bool b;
	L_ImageRGBUchar im;
	if ((b=im.readImage(name))==true)
		*this=im;
	return b;
}
bool L_ImageRGBDouble::saveImage(const char *name)
{
	L_ImageRGBUchar im;
	im=*this;
	return im.saveImage(name);
}

L_ImageRGBDouble &L_ImageRGBDouble::conv2d(const L_ImageRGBDouble &other, const L_ImageRGBDouble &mask, int xCenMasc, int yCenMasc)
{
	int i, j, u, v;
	reallocate(other.lx, other.ly);
	for (j=0; j<ly; j++)
	{
		for (i=0; i<lx; i++)
		{
			pix(i,j,0)=0;
			pix(i,j,1)=0;
			pix(i,j,2)=0;
		}
	}

	//  0  1  ... xCenMasc ... ... lx-2  lx-1
	//  <--- xC+1 ---> <-----  lx-xC  ---->

	for (i=xCenMasc; i<lx-(mask.lx-xCenMasc); i++)
	{
		for (j=yCenMasc; j<ly-(mask.ly-yCenMasc); j++)
		{
			//pix(i,j)=0;
			for (v=0; v<mask.ly; v++)
				for (u=0; u<mask.lx; u++)
				{
					pix(i,j,0)+=other.pix(i+u-xCenMasc,j+v-yCenMasc,0)*mask.pix(u,v,0);
					pix(i,j,1)+=other.pix(i+u-xCenMasc,j+v-yCenMasc,1)*mask.pix(u,v,1);
					pix(i,j,2)+=other.pix(i+u-xCenMasc,j+v-yCenMasc,2)*mask.pix(u,v,2);
				}
		}
	}
	// De ahi veo los bordes...
	return *this;
}


#ifdef __COMPAT_ATLIMAGE__
int L_ImageRGBDouble::readImageATL(const char *name)
{
	int a;
	L_ImageRGBUchar im;
	a=im.readImageATL(name);
	if (a==0)
		return 0;
	*this=im;
	return a;
}
bool L_ImageRGBDouble::writeImageATL(const char *name)
{
	bool a;
	L_ImageRGBUchar im;
	im=*this;
	a=im.writeImageATL(name);
	return a;
}
#endif


#ifdef __COMPAT_ATLIMAGE__
void L_ImageRGBDouble::operator=(const ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=other;
	*this=im;
}
void L_ImageRGBDouble::copyTo(ATL_CImage &other)
{
	L_ImageRGBUchar im;
	im=*this;
	im.copyTo(other);
}
#endif

#ifdef __COMPAT_IPLIMAGE__
void L_ImageRGBDouble::operator=(const IplImage &other)
{
	L_ImageRGBUchar im;
	im=other;
	*this=im;
}
void L_ImageRGBDouble::copyTo(IplImage **other)
{
	L_ImageRGBUchar im;
	im=*this;
	im.copyTo(other);
}
#endif


#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageRGBDouble::crearEigenLentaR() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			mat(i,j) = operator()(i,j,0);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageRGBDouble::crearEigenLentaG() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			mat(i,j) = operator()(i,j,1);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_ImageRGBDouble::crearEigenLentaB() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			mat(i,j) = operator()(i,j,2);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
void L_ImageRGBDouble::operator=(const Eigen::MatrixXd &other)
{
	Eigen::MatrixXd m;
	int i, j;
	reallocate_fc((int)other.rows(), (int)other.cols());
	//memcpy(&elem[0][0],  Eigen::Transpose<Eigen::MatrixXd>(other).data(), sizeof(double)*li*lj);
	m = Eigen::Transpose<Eigen::MatrixXd>(const_cast<Eigen::MatrixXd&>(other));
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			operator()(i,j,0) = other(i,j);
			operator()(i,j,1) = other(i,j);
			operator()(i,j,2) = other(i,j);
		}
	}
}
#endif // __COMPAT_EIGEN__

int L_ImageRGBUchar::readBMP(const char *name)
{
	FILE *fp;
	L_BMPheaderInfo inf;
	struct { unsigned int b:5;unsigned int g:5;unsigned int r:5;unsigned int z:1; } rgb16;
	int i, j;

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

	fp=fopen(name,"rb");
	if (fp==NULL)
		return 0;
	if (inf.read(fp)==false)
	{fclose(fp); printf("Error: Archivo BMP debe ser de 16, 24 o 32 bpp\n"); return 0;}
	reallocate(inf.width, inf.height);
	fseek(fp, inf._offsetData, SEEK_SET);

	if (inf.bits==16)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				if (fread(&rgb16,sizeof(rgb16),1,fp) != 1)
					{fclose(fp); return 0;}
				pix(i,j,0)=((L_uchar)rgb16.r);
				pix(i,j,1)=((L_uchar)rgb16.g);
				pix(i,j,2)=((L_uchar)rgb16.b);
			}
			fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}

	else if (inf.bits==24)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				int a, b, c;
				a = getc(fp);
				b = getc(fp);
				c = getc(fp);
				if (a == EOF || b == EOF || c == EOF)
					{fclose(fp); return 0;}
				pix(i,j,2)=(L_uchar)a;
				pix(i,j,1)=(L_uchar)b;
				pix(i,j,0)=(L_uchar)c;
			}
			fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}

	else if (inf.bits==32)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				int a, b, c, d;
				a = getc(fp);
				b = getc(fp);
				c = getc(fp);
				d = getc(fp);
				if (a == EOF || b == EOF || c == EOF || d == EOF)
					{fclose(fp); return 0;}
				pix(i,j,2)=(L_uchar)a;
				pix(i,j,1)=(L_uchar)b;
				pix(i,j,0)=(L_uchar)c;
			}
			fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}
	else
	{
		fprintf(stderr," > No puedo read BMP con paletas de colores\n\r");
		destroy();
		fclose(fp);
		return 0;
	}
	fclose(fp);
	return inf.bits;
}

bool L_ImageRGBUchar::writeBMP(const char *name, short nb)
{
	FILE *fp;
	L_BMPheaderInfo inf;
	int i, j;
	struct {unsigned int b:5;unsigned int g:5;unsigned int r:5;unsigned int z:1;} rgb16;

	fp=fopen(name,"wb");
	if (fp==NULL)
		return false;
	inf.setSize(lx, ly, nb);
	inf.write(fp);
	fseek(fp, inf._offsetData, SEEK_SET);
	if (inf.bits==16)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				rgb16.r=pix(i,j,0);
				rgb16.g=pix(i,j,1);
				rgb16.b=pix(i,j,2);
				rgb16.z=0;
				fwrite(&rgb16,sizeof(rgb16),1,fp);
			}
			for (i=0; i<inf.whiteSpace; i++)
				putc(0, fp);
			//fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}
	else if (inf.bits==24)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				putc(pix(i,j,2),fp);
				putc(pix(i,j,1),fp);
				putc(pix(i,j,0),fp);
			}
			for (i=0; i<inf.whiteSpace; i++)
				putc(0, fp);
			//fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}
	else if (inf.bits==32)
	{
		for (j=ly-1; j>=0; j--)
		{
			for (i=0; i<lx; i++)
			{
				putc(pix(i,j,2),fp);
				putc(pix(i,j,1),fp);
				putc(pix(i,j,0),fp);
				putc(0,fp);
			}
			for (i=0; i<inf.whiteSpace; i++)
				putc(0, fp);
			//fseek(fp,inf.whiteSpace,SEEK_CUR);
		}
	}
	else
	{
		fprintf(stderr,"Problemas con nb al write BMP\n");
		return false;
	}
	fclose(fp);
	return true;
}


// Lector / grabador de imagenes generico, usa lo que haya

bool L_ImageRGBUchar::readImage(const char *name_str)
{
#if defined __ATLIMAGE_H__
	return readImageATL(name);
#elif defined __COMPAT_IPLIMAGE__
	return leeImagenIPL(name);
#elif defined __COMPAT_TBITMAP__
	return leeImagenTBitmap(name);
#else
	L_String name;
	name = name_str;
	name.lowercase();
	if (name.size() < 5 || name[name.size()-5] != '.' || name[name.size()-4] != 'b' || name[name.size()-3] != 'm' || name[name.size()-2] != 'p')
		printf(" (leyendo %s como .BMP)\n", name.data());
	return readBMP(name.data())!=0;
#endif
}

bool L_ImageRGBUchar::saveImage(const char *name_str)
{
#if defined __ATLIMAGE_H__
	return writeImageATL(name);
#elif defined __COMPAT_IPLIMAGE__
	return grabaImagenIPL(name);
#elif defined __COMPAT_TBITMAP__
	return grabaImagenTBitmap(name);
#else
	L_String name;
	name = name_str;
	name.lowercase();
	if (name.size() < 5 || name[name.size()-5] != '.' || name[name.size()-4] != 'b' || name[name.size()-3] != 'm' || name[name.size()-2] != 'p')
		printf(" (grabando %s como .BMP)\n", name.data());
	return writeBMP(name.data());
#endif
}




#ifdef __ATLIMAGE_H__
bool L_ImageRGBUchar::readImageATL(const char *name)
{
	ATL_CImage imATL;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
#ifdef UNICODE
#error Usando unicode
    size_t origsize = strlen(name) + 1;
    const size_t newsize = 100;
    size_t convertedChars = 0;
    wchar_t wcstring[newsize];
    mbstowcs_s(&convertedChars, wcstring, origsize, name, _TRUNCATE);
	if (FAILED(imATL.Load(wcstring))) // Si reclama elegir Project -> X Properties -> General -> Character Set = Not Set
		return false;
#else
	if (FAILED(imATL.Load(name))) // Si reclama elegir Project -> X Properties -> General -> Character Set = Not Set
		return false;
#endif
	*this=imATL;
	return true;
}

bool L_ImageRGBUchar::writeImageATL(const char *name)
{
	ATL_CImage imATL;
	copyTo(imATL);
#ifdef UNICODE
    size_t origsize = strlen(name) + 1;
    const size_t newsize = 100;
    size_t convertedChars = 0;
    wchar_t wcstring[newsize];
    mbstowcs_s(&convertedChars, wcstring, origsize, name, _TRUNCATE);
	if (FAILED(imATL.Save(wcstring)))
		return false;
#else
	if (FAILED(imATL.Save(name)))
		return false;
#endif
	return true;
}
#endif //__ATLIMAGE_H__


#ifdef __COMPAT_IPLIMAGE__
bool L_ImageRGBUchar::leeImagenIPL(const char *name)
{
	IplImage *imIPL;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	imIPL=cvLoadImage(name,1); // 1 : read varios colores
	if (imIPL == NULL)
		return false;
	*this=*imIPL;
	cvReleaseImage(&imIPL);
	return true;
}

bool L_ImageRGBUchar::grabaImagenIPL(const char *name)
{
	IplImage *imIPL;
	copyTo(&imIPL);
	cvSaveImage(name, imIPL);
	cvReleaseImage(&imIPL);
	return true;
}
#endif

#ifdef __COMPAT_CVMAT__
void L_ImageRGBUchar::copyTo(cv::Mat& other)
{
	other = cv::Mat::zeros(rows(), cols(), CV_8UC3);
	for (int i = 0; i < other.rows; i++)
	{
		for (int j = 0; j < other.cols; j++)
		{
			other.at<cv::Vec3b>(i, j)[0] = pix(j, i, 0);
			other.at<cv::Vec3b>(i, j)[1] = pix(j, i, 1);
			other.at<cv::Vec3b>(i, j)[2] = pix(j, i, 2);
		}
	}
}
#endif

#ifdef GraphicsHPP
bool L_ImageRGBUchar::leeImagenTBitmap(const char *name)
{
	Graphics::TBitmap *pBitmap = new Graphics::TBitmap();
	PRGBTriple srcRow;
	long len;
	try
	{
		len=strlen(name);
		if (len>4 && (name[len-2]=='p' || name[len-2]=='P') && (name[len-1]=='g' || name[len-1]=='G')) // Un jpg ?
		{
			TJPEGImage *jp = new TJPEGImage();
			jp->LoadFromFile(name);
			pBitmap->Assign(jp);
		}
		else
			pBitmap->LoadFromFile(name);  // Un bmp ?
		reallocate(pBitmap->Width, pBitmap->Height);
		for (int y = 0; y < pBitmap->Height; y++)
		{
			srcRow=(PRGBTriple)pBitmap->ScanLine[y];
			for (int x = 0; x < pBitmap->Width; x++)
			{
				pix[x][y][0]=srcRow[x].rgbtRed;
				pix[x][y][1]=srcRow[x].rgbtGreen;
				pix[x][y][2]=srcRow[x].rgbtBlue;
			}
		}
	}
	catch (...)
	{
		return false;
	}

	delete pBitmap;
	return true;
}

bool L_ImageRGBUchar::grabaImagenTBitmap(const char *name)
{
	Graphics::TBitmap *pBitmap = new Graphics::TBitmap();
	pBitmap->Width=lx;
	pBitmap->Height=ly;
	PRGBTriple srcRow;
	long len;
	reallocate(pBitmap->Width, pBitmap->Height);
	try
	{
		for (int y = 0; y < pBitmap->Height; y++)
		{
			srcRow=(PRGBTriple)pBitmap->ScanLine[y];
			for (int x = 0; x < pBitmap->Width; x++)
			{
				srcRow[x].rgbtRed=pix[x][y][0];
				srcRow[x].rgbtGreen=pix[x][y][1];
				srcRow[x].rgbtBlue=pix[x][y][2];
			}
		}
		len=strlen(name);
		if (len>4 && (name[len-2]=='p' || name[len-2]=='P') && (name[len-1]=='g' || name[len-1]=='G')) // Un jpg ?
		{
			TJPEGImage *jp = new TJPEGImage();
			jp->Assign(pBitmap);
			jp->SaveToFile(name);
		}
		else
			pBitmap->SaveToFile(name); // Un bmp ?
	}
	catch (...)
	{
		return false;
	}

	delete pBitmap;
	return true;
}
#endif


bool L_ImageRGBUchar::leeB24(const char *name)
{
	FILE *fp;
	int i, j;
	int c;
#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif
	fp=fopen(name,"rb");
	if (fp==NULL)
		return false;
	lx=(unsigned char)getc(fp);
	lx=lx*256+(unsigned char)getc(fp);
	ly=(unsigned char)getc(fp);
	ly=ly*256+(unsigned char)getc(fp);
	reallocate(lx, ly);
	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			if ((c=getc(fp))==EOF) {fclose(fp);return false;}
			pix(i,j,0)=(L_uchar)c;
			if ((c=getc(fp))==EOF) {fclose(fp);return false;}
			pix(i,j,1)=(L_uchar)c;
			if ((c=getc(fp))==EOF) {fclose(fp);return false;}
			pix(i,j,2)=(L_uchar)c;
		}
	}
	fclose(fp);
	return true;
}

bool L_ImageRGBUchar::grabaB24(const char *name)
{
	FILE *fp;
	int i, j;
	if (lx<0 || ((long)lx)>0xFFFFL || ly<0 || ((long)ly)>=0xFFFFL)
		return false;
	fp=fopen(name,"wb");
	if (fp==NULL)
		return false;
	putc(lx>>8, fp);
	putc(lx&0xFF, fp);
	putc(ly>>8, fp);
	putc(ly&0xFF, fp);
	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			putc(pix(i,j,0),fp);
			putc(pix(i,j,1),fp);
			putc(pix(i,j,2),fp);
		}
	}
	fclose(fp);
	return true;
}


int L_ImageRGBUchar::writeComprRGB(FILE *fp, double *le, int forzMetodo, double t, long ncam) const
{
	std::vector<char> buf;
	std::vector<char> bufD;
	std::vector<char> bufDD;
	std::vector<char> bufRLE;
	std::vector<char> bufT;
	double largoEst[5]; // 0 = sin compr;  1 = compr directo;  2 = delta y compr
	int i, iMin;
	std::vector<char>::size_type sz, lxs, lys;
	if (le == NULL)
		le = largoEst;

	// 0: sin comprimir, 1: Huffman, 2:resta-huffman, 3:resta-RLE_huffman
	// Calculos compresion Hough simple
	switch(forzMetodo)
	{
	case -1: // Todos
		//buf.setData( (char *)&pix(0,0,0), lx*ly*3);  // Better slow than bad coding
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0,0), lx*ly*3);
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		L_HuffmanEncoder::encodeRLE0(bufDD, bufRLE);
		break;
	case 0:
		//buf.setData( (char *)&pix(0,0,0), lx*ly*3);
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0,0), lx*ly*3);
		break;
	case 1:
		//buf.setData( (char *)&pix(0,0,0), lx*ly*3);
		buf.resize(lx*ly);
		memcpy(&(buf[0]), (char *)&pix(0,0,0), lx*ly*3);
		break;
	case 2:
		encodeDelta(bufD);
		break;
	case 3:
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		break;
	case 4:
		encodeDelta(bufD);
		L_HuffmanEncoder::encodeDelta(bufD, bufDD);
		L_HuffmanEncoder::encodeRLE0(bufDD, bufRLE);
		break;
	default:
		//if (buf.size()!=0)
		//	buf.setData(NULL,0);
		printf("Modo de compresion no disponible\n");
		return -1;
	}

	// Evaluar el mas conveniente
	if (forzMetodo == -1)
	{
		le[0] = buf.size();  // Largo sin comprimir
		le[1] = 165 + L_HuffmanEncoder::entropy_in_bits(buf)/8*buf.size(); // Largo aproximado del encabezado + comprimido
		le[2] = 165 + L_HuffmanEncoder::entropy_in_bits(bufD)/8*bufD.size(); // Largo aproximado del encabezado + comprimido con delta
		le[3] = 165 + L_HuffmanEncoder::entropy_in_bits(bufDD)/8*bufDD.size(); // Largo aproximado del encabezado + comprimido con delta
		le[4] = 165 + L_HuffmanEncoder::entropy_in_bits(bufRLE)/8*bufRLE.size(); // Largo aproximado del encabezado + comprimido con delta

		// El calculo de bufDD es tonto (es restar saltandose 1), mejor probar con algo que elimine ruido como:
		//   buf[i] - buf[i-2]/2 - buf[i-1]/2
		//   buf[i] - (buf[i-2]+2*(int)buf[i-1])/3

		iMin = 0;
		for (i=1; i<5; i++)
			if (le[i] < le[iMin])
				iMin = i;
	}
	else
		iMin = forzMetodo;

	//iMin = 2; // Forzar las cosas

	lxs = lx;
	lys = ly;

	// Escritura 1
	fputc('B',fp);
	fputc('M',fp);
	fputc('L',fp);  // ex BMZ...
	fputc('3',fp);
	fwrite(&lxs, sizeof(lxs), 1, fp);
	fwrite(&lys, sizeof(lys), 1, fp);
	fwrite(&t, sizeof(t), 1, fp);
	fwrite(&ncam, sizeof(ncam), 1, fp);
	fputc((char)iMin,fp);

	switch(iMin)  // 0 = sin compr;  1 = compr directo;  2 = delta y compr
	{
	case 0:
		sz = buf.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(buf[0]), 1, buf.size(), fp);
		break;
	case 1:
		L_HuffmanEncoder::encodeAll(buf, bufT);
		sz = bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 2:
		L_HuffmanEncoder::encodeAll(bufD, bufT);
		sz = bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 3:
		L_HuffmanEncoder::encodeAll(bufDD, bufT);
		sz = bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	case 4:
		L_HuffmanEncoder::encodeAll(bufRLE, bufT);
		sz = bufT.size();
		fwrite(&sz, sizeof(sz), 1, fp);
		fwrite(&(bufT[0]), sizeof(char), bufT.size(), fp);
		break;
	}
	//if (buf.size()!=0)
	//	buf.setData(NULL,0);
	return iMin;
}

int L_ImageRGBUchar::readComprRGB(FILE *fp, double *t, long *ncam)
{
	std::vector<char> buf;
	std::vector<char> bufD;
	std::vector<char> bufDD;
	std::vector<char> bufRLE;
	std::vector<char> bufT;
	int iMin, sz;
	long lxs, lys, nc;
	double ti;

	// Lectura 1
	int head[4];

	head[0] = fgetc(fp);
	head[1] = fgetc(fp);
	head[2] = fgetc(fp);
	head[3] = fgetc(fp);

	// Lectura 1
	if ( (head[0] != 'B' || head[1] != 'M' || head[2] != 'Z' || head[3] != '3') &&
		(head[0] != 'B' || head[1] != 'M' || head[2] != 'L' || head[3] != '3'))
		return -1;
	size_t nleidos = 0;
	nleidos += fread(&lxs, sizeof(lxs), 1, fp);
	nleidos += fread(&lys, sizeof(lys), 1, fp);
	nleidos += fread(&ti, sizeof(ti), 1, fp);
	nleidos += fread(&nc, sizeof(nc), 1, fp);
	if (nleidos != 4)
		return -1;		
	iMin = fgetc(fp);
	if (iMin == EOF)
		return -1;

	if (t!=NULL)
		*t = ti;
	if (ncam != NULL)
		*ncam = nc;

	reallocate(lxs, lys);
	buf.resize(3*lx*ly);

	switch(iMin)  // 0 = sin compr;  1 = compr directo;  2 = delta y compr  LECTURA
	{
	case 0:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		buf.resize(sz);
		if (fread(&(buf[0]), 1, buf.size(), fp) != buf.size())
			return -1;
		decodeFromBuffer(buf);
		break;
	case 1:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, buf);
		decodeFromBuffer(buf);
		break;
	case 2:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufD);
		decodeDelta(bufD);
		break;
	case 3:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufDD);
		L_HuffmanEncoder::decodeDelta(bufDD, bufD);
		decodeDelta(bufD);
		break;
	case 4:
		if (fread(&sz, sizeof(sz), 1, fp) != 1)
			return -1;
		bufT.resize(sz);
		if (fread(&(bufT[0]), 1, bufT.size(), fp) != bufT.size())
			return -1;
		L_HuffmanEncoder::decodeAll(bufT, bufRLE);
		L_HuffmanEncoder::decodeRLE0(bufRLE, bufDD);
		L_HuffmanEncoder::decodeDelta(bufDD, bufD);
		decodeDelta(bufD);
	}
	if (head[2] == 'Z')
		transposeMe();
	return iMin;
}

int L_ImageRGBUchar::leeComprRGB_contar(FILE *fp, int *histMetodo)
{
	long pos0, pos, posFin;
	int iMin, num=0;
	long lxs, lys, nc;
	int N;
	size_t nleidos;
	double ti;
	int head[4];

	pos0 = ftell(fp);
	fseek(fp, 0, SEEK_END);
	posFin = ftell(fp);
	fseek(fp, pos0, SEEK_SET);

	while(!feof(fp))
	{
		pos = ftell(fp);
		if (pos == posFin)
			break;

		head[0] = fgetc(fp);
		head[1] = fgetc(fp);
		head[2] = fgetc(fp);
		head[3] = fgetc(fp);

		if ( (head[0] != 'B' || head[1] != 'M' || head[2] != 'Z' || head[3] != '3') &&
			(head[0] != 'B' || head[1] != 'M' || head[2] != 'L' || head[3] != '3'))
			return -1;
		nleidos = 0;
		nleidos += fread(&lxs, sizeof(lxs), 1, fp);
		nleidos += fread(&lys, sizeof(lys), 1, fp);
		nleidos += fread(&ti, sizeof(ti), 1, fp);
		nleidos += fread(&nc, sizeof(nc), 1, fp);
		if (nleidos != 4)
			return -1;
		iMin = getc(fp);
		if (iMin == EOF)
			return -1;
		if (histMetodo != NULL)
			(histMetodo)[iMin] ++;
		if (fread(&N, sizeof(N), 1, fp) != 1) // bufT
			return -1;
		fseek(fp, N, SEEK_CUR);
		num++;
	}
	fseek(fp, pos0, SEEK_SET);
	return num;
}


void L_ImageRGBUchar::testComprRGB()
{
	FILE *fp = NULL;
	int u, i, j, val=0;
	L_ImageRGBUchar im1, im2;
	const char name[] = "pruebitaaa.bmz";
	srand((unsigned int)time(NULL));

	for (u=0; u<60; u++)
	{
		im1.reallocate(10+rand()%310, 10+rand()%190);
		//im1.fillRandomValues(10,30);
		for (j=0; j<im1.ly; j++)
		{
			for (i=0; i<im1.lx; i++)
			{
				// metodo 1
				val = rand()%10;
				// metodo 2
				//val = val + rand()%10 - 5;
				im1.pix(i,j,0) = val; im1.pix(i,j,1) = val+7; im1.pix(i,j,2) = val+12;
			}
		}
		fp = fopen(name, "wb");
		if (fp == NULL)
		{
			printf("Fallo la escritura de %s\n", name);
			return;
		}
		printf("metodo %d\n", im1.writeComprRGB(fp));
		fclose(fp);
		fp = fopen(name, "rb");
		if (fp == NULL)
		{
			printf("Fallo la lectura de %s\n", name);
			return;
		}
		im2.readComprRGB(fp);
		fclose(fp);
		for (j=0; j<im1.ly; j++)
			for (i=0; i<im1.lx; i++)
				if (im1.pix(i,j,0) != im2.pix(i,j,0) || im1.pix(i,j,1) != im2.pix(i,j,1) || im1.pix(i,j,2) != im2.pix(i,j,2))
				{
					printf("testComprRGB() fallo\n");
					return;
				}
	}
	printf("testComprRGB() funciono adecuadamente\n");
}


void L_ImageRGBUchar::pruebaVelocComprRGB()
{
	FILE *fp = NULL;
	int u, i, j, val=0, vald=0, vald2=0; //, M;
	L_ImageRGBUchar im1, im2;
	double t1, t2;
	const char name[] = "pruebitaaa.bmz";

	srand((unsigned int)time(NULL));
	int v0 = rand();

	fp = fopen(name, "wb");
	t1 = L_TIME();
	for (u=0; u<120; u++)
	{
		im1.reallocate(320, 200);
		for (j=0; j<im1.ly; j++)
		{
			for (i=0; i<im1.lx; i++)
			{
				v0 = v0*13828 + v0*v0*2728 + 21892;
				vald2 = ((unsigned char)v0)%10;
				vald = vald + vald2;
				val = val + vald;
				im1.pix(i,j,0) = val; im1.pix(i,j,1) = val+20; im1.pix(i,j,2) = val+40;
			}
		}
		if (fp == NULL)
		{
			printf("Fallo la escritura de %s\n", name);
			return;
		}
		im1.writeComprRGB(fp, NULL, 4);
	}
	t2 = L_TIME();
	fclose(fp);
	printf("pruebaVelocComprRGB() funciono adecuadamente a %d[fps]\n", (int)(120/(t2-t1)));
}

bool L_ImageRGBUchar::normalizeHistogram()
{
	int i, j, c;
	L_uchar max, min;
	L_uchar dif;
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
				pix(i,j,c)-=(L_uchar)min;
				pix(i,j,c)/=(L_uchar)dif;
			}
		}
	}
	return ret;
}


L_ImageRGBUchar& L_ImageRGBUchar::transformaDaltonicoLonco(bool transfInversa)
{
	L_ImageRGBDouble imfloat;
	int i, j;
	imfloat=*this;
	imfloat.RGB_to_YPbPr();
	if (transfInversa)
		for (j=0; j<ly; j++)
			for (i=0; i<lx; i++)
				imfloat.pix(i,j,2)*=2;
	else
		for (j=0; j<ly; j++)
			for (i=0; i<lx; i++)
				imfloat.pix(i,j,2)*=0.5;
	imfloat.YPbPr_to_RGB();
	*this=imfloat;
	return *this;
}

L_ImageRGBUchar& L_ImageRGBUchar::transformaDaltonicoLonco2(bool transfInversa)
{
	int i, j;
	if (transfInversa)
		for (j=0; j<ly; j++)
			for (i=0; i<lx; i++)
			{
				pix(i,j,1)/=2;pix(i,j,2)/=2;
			}
	else
		for (j=0; j<ly; j++)
			for (i=0; i<lx; i++)
				pix(i,j,0)/=2;
	return *this;
}



void L_ImageRGBUchar::computePositionsForComposingImages(int n, int *imLx, int *imLy, int *imPosX, int *imPosY, int &tamTotalX, int &tamTotalY)
{
	int i;
	int iGrupo[4]; // izq, arr, aba, der
	long sumaX=0, sumaY=0, s;
	int anchoGrupo[4]={0,0,0,0}, altoGrupo[4]={0,0,0,0};
	int alto01, ancho012;

	// Calcular anchos y altos
	for (i=1; i<n; i++)
	{
		sumaX+=imLx[i];
		sumaY+=imLy[i];
	}
	iGrupo[0]=1;
	for (i=iGrupo[0]; i<n; i++)
	{
		if (altoGrupo[0]>=sumaY/4)
			break;
		anchoGrupo[0]+=imLx[i];
		altoGrupo[0]+=imLy[i];
	}
	iGrupo[1]=i;
	for (i=iGrupo[1]; i<n; i++)
	{
		if (anchoGrupo[1]>=sumaX/4)
			break;
		anchoGrupo[1]+=imLx[i];
		altoGrupo[1]+=imLy[i];
	}
	iGrupo[2]=i;
	for (i=iGrupo[2]; i<n; i++)
	{
		if (anchoGrupo[2]>=sumaX/4)
			break;
		anchoGrupo[2]+=imLx[i];
		altoGrupo[2]+=imLy[i];
	}
	iGrupo[3]=i;
	for (i=iGrupo[3]; i<n; i++)
	{
		anchoGrupo[3]+=imLx[i];
		altoGrupo[3]+=imLy[i];
	}
	// Determinar posiciones
	s=0;
	for (i=iGrupo[0]; i<iGrupo[1]; i++)
	{
		imPosX[i]=0;
		imPosY[i]=s;
		s+=imLy[i];
	}
	s=0;
	for (i=iGrupo[1]; i<iGrupo[2]; i++)
	{
		imPosX[i]=s+anchoGrupo[0];
		imPosY[i]=0;
		s+=imLy[i];
	}
	s=0;
	alto01=L_MAX(imLy[0],L_MAX(altoGrupo[0],altoGrupo[3])-altoGrupo[1]);
	for (i=iGrupo[2]; i<iGrupo[3]; i++)
	{
		imPosX[i]=s+anchoGrupo[0];
		imPosY[i]=altoGrupo[1]+alto01;
		s+=imLy[i];
	}
	s=0;
	ancho012=L_MAX(L_MAX(imLx[0],anchoGrupo[1]),anchoGrupo[2]);
	for (i=iGrupo[3]; i<n; i++)
	{
		imPosX[i]=ancho012+anchoGrupo[0];
		imPosY[i]=altoGrupo[1]+imLy[0];
		s+=imLy[i];
	}
	imPosX[0]=anchoGrupo[0]+ancho012/2-imLx[0]/2;
	imPosY[0]=altoGrupo[1]+alto01/2-imLy[0]/2;
	tamTotalX=anchoGrupo[0]+ancho012+anchoGrupo[3];
	tamTotalY=altoGrupo[1]+alto01+altoGrupo[2];
}


L_ImageRGBUchar& L_ImageRGBUchar::genDrawing(const L_ShapeArray &lins)
{
	int a;
	int u, u1, u2, v1, uc;
	int dx, dy;
	int dxA, dyA;
	double m;
	double x, y;
	int xi, yi, ui;
	// Para elipse
	double maxr, minr;
	int numPuntos;
	double alfa, sinalfa, cosalfa;
	double beta, sinbeta, cosbeta;

#ifdef VERIF_INT
	if (L_VerificaIntegridadCodigoRand(0.005)==false)
		L_hard_shutdown("Error de integridad\n");
#endif

	for (a=0; a<(int)lins.size(); a++)
	{
		switch(lins[a].tipo)
		{
		case L_Shape_undefined:
			throw_L_ArgException_if(true, "L_ImageRGBUchar::genDrawing() : tipo de linea indefinido");
			break;
		case L_Shape_line:
			if (lins[a].x1!=lins[a].x1 || lins[a].x2!=lins[a].x2 || lins[a].y1!=lins[a].y1 || lins[a].y2!=lins[a].y2) // NaN
				continue;
			dx=lins[a].x2-lins[a].x1;
			dy=lins[a].y2-lins[a].y1;
			if (dx==0 && dy==0 && lins[a].x1>=0 && lins[a].x1<lx && lins[a].y1>=0 && lins[a].y1<ly)
			{
				pix(lins[a].x1,lins[a].y1,0)=lins[a].R;
				pix(lins[a].x1,lins[a].y1,1)=lins[a].G;
				pix(lins[a].x1,lins[a].y1,2)=lins[a].B;
				continue;
			}
			if (dx>=0)
				dxA=dx;
			else
				dxA=-dx;
			if (dy>=0)
				dyA=dy;
			else
				dyA=-dy;
			if (dxA>=dyA)
			{
				int y;
				m=dy/(double)dx;
				if (dx>0)
				{
					u1=lins[a].x1;
					u2=lins[a].x2;
					v1=lins[a].y1;
				}
				else
				{
					u1=lins[a].x2;
					u2=lins[a].x1;
					v1=lins[a].y2;
				}
				u = u1;
				if (u<0)
					u=0;
				if (u > lx-1)
					u=lx-1;
				if (u2 > lx-1)
					u2 = lx-1;
				for ( ; u<=u2 ; u++)
				{
					y=(int)(v1+(u-u1)*m);
					if (y>=0 && y<ly)
					{
						pix(u,y,0)=lins[a].R;
						pix(u,y,1)=lins[a].G;
						pix(u,y,2)=lins[a].B;
					}
				}
			}
			else
			{
				int x;
				m=dx/(double)dy;
				if (dy>0)
				{
					u1=lins[a].y1;
					u2=lins[a].y2;
					v1=lins[a].x1;
				}
				else
				{
					u1=lins[a].y2;
					u2=lins[a].y1;
					v1=lins[a].x2;
				}

				u = u1;
				if (u < 0)
					u = 0;
				if (u > ly-1)
					u = ly-1;
				if (u2 > ly-1)
					u2 = ly-1;
				for ( ; u<=u2 ; u++)
				{
					x=(int)(v1+(u-u1)*m);
					if (x>=0 && x<lx)
					{
						pix(x,u,0)=lins[a].R;
						pix(x,u,1)=lins[a].G;
						pix(x,u,2)=lins[a].B;
					}
				}
			}
			break;
		case L_Shape_circle:
			u1=(int)(lins[a].x1-lins[a].rx)-1;
			u2=(int)(lins[a].x1+lins[a].rx)+1;
			uc=lins[a].x1;
			v1=lins[a].y1;
			if (lins[a].x1!=lins[a].x1 || lins[a].rx!=lins[a].rx) // NaN
				continue;

			for (u=u1; u<=u2; u++)
			{
				if (u<0)
					u=0;
				if (u>=lx)
					break;
				y=lins[a].rx*lins[a].rx-(u-uc)*(u-uc);
				if (y>0)
				{
					yi=(int)sqrt(y);
					if (v1+yi>=0 && v1+yi<ly)
					{
						pix(u,v1+yi,0)=lins[a].R;
						pix(u,v1+yi,1)=lins[a].G;
						pix(u,v1+yi,2)=lins[a].B;
					}
					if (v1-yi>=0 && v1-yi<ly)
					{
						pix(u,v1-yi,0)=lins[a].R;
						pix(u,v1-yi,1)=lins[a].G;
						pix(u,v1-yi,2)=lins[a].B;
					}
				}
			}
			u1=(int)(lins[a].y1-lins[a].rx)-1;
			u2=(int)(lins[a].y1+lins[a].rx)+1;
			uc=lins[a].y1;
			v1=lins[a].x1;
			for (u=u1; u<=u2; u++)
			{
				if (u<0)
					u=0;
				if (u>=ly)
					break;
				x=lins[a].rx*lins[a].rx-(u-uc)*(u-uc);
				if (x>0)
				{
					xi=(int)sqrt(x);
					if (v1+xi>=0 && v1+xi<lx)
					{
						pix(v1+xi,u,0)=lins[a].R;
						pix(v1+xi,u,1)=lins[a].G;
						pix(v1+xi,u,2)=lins[a].B;
					}
					if (v1-xi>=0 && v1-xi<lx)
					{
						pix(v1-xi,u,0)=lins[a].R;
						pix(v1-xi,u,1)=lins[a].G;
						pix(v1-xi,u,2)=lins[a].B;
					}
				}
			}
			break;
		case L_Shape_ellipse:
		case L_Shape_dashed_ellipse:
			maxr = (lins[a].rx > lins[a].ry) ? lins[a].rx : lins[a].ry;
			minr = (lins[a].rx < lins[a].ry) ? lins[a].rx : lins[a].ry;
			numPuntos = (int)(2.0*M_PI*maxr + 5);

			// Ver si es necesario dibujar la elipse
			if (minr*minr > lx*lx + ly*ly || minr < 0 || maxr < 0)
				continue;
			if (maxr!=maxr || minr!=minr || numPuntos <= 0)
				continue;

			if (numPuntos > 10000)
				continue;

			beta = -lins[a].ang;
			sinbeta = sin(beta);
			cosbeta = cos(beta);

			for (double i = 0; i < 360; i += 360.0 / numPuntos) 
			{
				alfa = i * (M_PI / 180) ;
				sinalfa = sin(alfa);
				cosalfa = cos(alfa);
				ui = (int)(i*numPuntos/200);

				if (lins[a].tipo == L_Shape_dashed_ellipse && ui%2 == 0) // Tachado
					continue;

				int X = (int)( lins[a].x1 + (lins[a].rx * cosalfa * cosbeta - lins[a].ry * sinalfa * sinbeta) );
				int Y = (int)( lins[a].y1 + (lins[a].rx * cosalfa * sinbeta + lins[a].ry * sinalfa * cosbeta) );

				if (X>=0 && Y>=0 && X<lx && Y<ly)
				{
					pix(X,Y,0) = lins[a].R;
					pix(X,Y,1) = lins[a].G;
					pix(X,Y,2) = lins[a].B;
				}
			}
			break;
		}
	}
	return *this;
}

L_ImageRGBUchar & L_ImageRGBUchar::genDrawingMatches(L_Array<const L_ImageGrayDouble *> &imarr, L_ShapeArray &lins) // imarr[0] es la imagen de test
{
	std::vector<int> imLx(imarr.size()), imLy(imarr.size()), imPosX(imarr.size()), imPosY(imarr.size());
	int tamTotalX, tamTotalY, i;

#ifdef L_BASIC_DEBUG
	if (imarr.size()<=0)
		L_hard_shutdown("L_ImageRGBUchar::genDrawingMatches : non-positive size allocation");
#endif
	for (i=0; i<(int)imarr.size(); i++)
	{
		imLx[i]=imarr[i]->lx;
		imLy[i]=imarr[i]->ly;
	}

	computePositionsForComposingImages(imarr.size(), &(imLx[0]), &(imLy[0]), &(imPosX[0]), &(imPosY[0]), tamTotalX, tamTotalY);
	reallocate(tamTotalX, tamTotalY);
	fill(255,255,255);
	for (i=0; i<(int)imarr.size(); i++)
		copyImageWithOffset(*imarr[i],imPosX[i],imPosY[i]);

	for (int a=0; a<(int)lins.size(); a++)
	{
		lins[a].x1+=imPosX[lins[a].n1];
		lins[a].y1+=imPosY[lins[a].n1];
		lins[a].x2+=imPosX[lins[a].n2];
		lins[a].y2+=imPosY[lins[a].n2];
		lins[a].n1=0;
		lins[a].n2=0;
	}

	genDrawing(lins);
	return *this;
}

double L_ImageRGBUchar::computeMeanShiftTracking_private(double &xCen, double &yCen, std::vector<double> &hist, double radio_h_x, double radio_h_y, int nBinsPorCanal)
{
	std::vector<double> histZ;
	double parec, parecZ;
	double x = xCen, y = yCen;
	int i, niter = 0, u = 0;
	bool inicial = true;

	// Histograma del modelo = hist
	// Histograma actual inicial
	computeMeanShiftHistogram(x, y, histZ, radio_h_x, radio_h_y, nBinsPorCanal);

	//Error inicial
	parecZ = 0;
	for (i=0; i<(int)hist.size(); i++)
		parecZ += sqrt(hist[i]*histZ[i]);

	while (inicial || fabs(x-xCen) > 1 || fabs(y-yCen) > 1)
	{
		inicial = false;
		xCen = x; // Actualizar valores de la iteracion anterior
		yCen = y;
		parec = parecZ;

		computeMeanShiftTracking_private(x, y, hist, histZ, radio_h_x, radio_h_y, nBinsPorCanal);

		// Calcular el nuevo histograma actual y el nuevo error
		computeMeanShiftHistogram(x, y, histZ, radio_h_x, radio_h_y, nBinsPorCanal);
		parecZ = 0;
		for (i=0; i<(int)hist.size(); i++)
			parecZ += sqrt(hist[i]*histZ[i]);

		// Ver si el parecido a la distribucion original mejoro; de lo contrario devolverse a la mitad del camino
		for (u=0; u < 5 && parecZ < parec; u++)
		{
			// Hacer que los valores desplazados se acerquen a los antiguos
			x = (x + xCen) * 0.5;
			y = (y + yCen) * 0.5;
			// Reevaluar el error
			computeMeanShiftHistogram(x, y, histZ, radio_h_x, radio_h_y, nBinsPorCanal);
			parecZ = 0;
			for (i=0; i<(int)hist.size(); i++)
				parecZ += sqrt(hist[i]*histZ[i]);
		}
		niter++;
	}
	xCen = x; // Actualizar valores de la iteracion anterior
	yCen = y;
	return parecZ;
}

double L_ImageRGBUchar::computeMeanShiftTracking(double &xCen, double &yCen, std::vector<double> &hist, double &radio_h_x, double &radio_h_y, bool adaptTamano, int nBinsPorCanal)
{
	double p1, p2, p3, p;
	double x1, y1, x2, y2, x3, y3;
	double rx1, rx2, rx3, ry1, ry2, ry3;

	if (adaptTamano)
	{
		rx1 = radio_h_x * 0.9;
		rx2 = radio_h_x * 1.0;
		rx3 = radio_h_x * 1.1;
		ry1 = radio_h_y * 0.9;
		ry2 = radio_h_y * 1.0;
		ry3 = radio_h_y * 1.1;
		x1 = x2 = x3 = xCen;
		y1 = y2 = y3 = yCen;
		
		// Hacer mean-shift en tres escalas y evaluar los parecidos de los histogramas obtenidos
		p1 = computeMeanShiftTracking_private(x1, y1, hist, rx1, ry1, nBinsPorCanal); 
		p2 = computeMeanShiftTracking_private(x2, y2, hist, rx2, ry2, nBinsPorCanal); 
		p3 = computeMeanShiftTracking_private(x3, y3, hist, rx3, ry3, nBinsPorCanal); 

		// La escala que logra el histograma mas parecido es la elegida
		if (p1 > p2 && p1 > p3)
		{
			xCen = x1;
			yCen = y1;
			radio_h_x = rx1;
			radio_h_y = ry1;
			p = p1;
		}
		else if (p2 > p1 && p2 > p3)
		{
			xCen = x2;
			yCen = y2;
			radio_h_x = rx2;
			radio_h_y = ry2;
			p = p2;
		}
		else
		{
			xCen = x3;
			yCen = y3;
			radio_h_x = rx3;
			radio_h_y = ry3;
			p = p3;
		}
	}
	else
		p = computeMeanShiftTracking_private(xCen, yCen, hist, radio_h_x, radio_h_y, nBinsPorCanal);
	return p;
}

#ifdef __COMPAT_ATLIMAGE__ // Opcional, se debe define para activar compatibilidad con ATL (Windows)
bool L_ImageRGBUchar::copyFromsde_lento(const ATL_CImage &other)
{
	if (other.IsDIBSection()) // El caso facil
	{
		*this = other;
		return true;
	}

	BITMAP bmpScreen;
	HBITMAP hbmScreen = other.operator HBITMAP();
    GetObject(hbmScreen,sizeof(BITMAP),&bmpScreen);

    BITMAPINFOHEADER   bi;
     
    bi.biSize = sizeof(BITMAPINFOHEADER);    
    bi.biWidth = bmpScreen.bmWidth;    
    bi.biHeight = bmpScreen.bmHeight;  
    bi.biPlanes = 1;    
    bi.biBitCount = 32;    
    bi.biCompression = BI_RGB;    
    bi.biSizeImage = 0;  
    bi.biXPelsPerMeter = 0;    
    bi.biYPelsPerMeter = 0;    
    bi.biClrUsed = 0;    
    bi.biClrImportant = 0;

    DWORD dwBmpSize = ((bmpScreen.bmWidth * bi.biBitCount + 31) / 32) * 4 * bmpScreen.bmHeight;
	HDC hdcWindow = other.GetDC();

    // Starting with 32-bit Windows, GlobalAlloc and LocalAlloc are implemented as wrapper functions that 
    // call HeapAlloc using a handle to the process's default heap. Therefore, GlobalAlloc and LocalAlloc 
    // have greater overhead than HeapAlloc.
    char *lpbitmap = new char[dwBmpSize];    

    // Gets the "bits" from the bitmap and copies them into a buffer 
    // which is pointed to by lpbitmap.
    GetDIBits(hdcWindow, hbmScreen, 0,
        (UINT)bmpScreen.bmHeight,
        lpbitmap,
        (BITMAPINFO *)&bi, DIB_RGB_COLORS);

	int i, j;

	if (bmpScreen.bmWidth <= 0 || bmpScreen.bmHeight <= 0 || bmpScreen.bmWidth > 10000 || bmpScreen.bmHeight > 10000)
		return false;

	reallocate(bmpScreen.bmWidth, bmpScreen.bmHeight);

	for (i=0; i<lx; i++)
	{
		for (j=0; j<ly; j++)
		{
			pix(i,ly-j-1,0) = lpbitmap[4*j*lx+4*i+2];
			pix(i,ly-j-1,1) = lpbitmap[4*j*lx+4*i+1];
			pix(i,ly-j-1,2) = lpbitmap[4*j*lx+4*i+0];
		}
	}
	delete[] lpbitmap;
	other.ReleaseDC();
    DeleteObject(hbmScreen);
	return true;
}
#endif // __COMPAT_ATLIMAGE__


//////////


void L_Matrix::OP_mult(const L_Matrix &other)
{
	L_Matrix aux;
	aux.OP_mult(*this, other);
	aux.swap(*this);
}


void L_Matrix::OP_mult(const L_Matrix &other, double val)
{
	int i, j;
	throw_L_ArgException_if(other.begin()==NULL, "L_Matrix::OP_mult"); // Caso en que la other matriz no esta inicializada, esto podria tapar algunos errores...
	reallocate(other.li, other.lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j)=other(i,j)*val;
}

void L_Matrix::repairPositiveDefinite_valPr()
{
	L_Matrix Pr, D, PrT;
	swap(Pr);
	if (li == 2 && lj == 2)
		Pr.vectPr_2x2(D);
	else
		Pr.vectPr_sim(D);  // -> Pr * D * Pr^-1  =  Pr * D * PrT // Requiere matriz simetrica
	D.repairRangeValues(1e-10, 1e100);
	PrT.transpOf(Pr);
	Pr.OP_multByDiagonal(D);
	Pr.OP_mult(PrT);
	Pr.swap(*this);
}


//bool L_Matrix::repararDefinidaPositivaSimetrica_LLT()
//{
//	L_Matrix L, LT;
//	if (L.cholesky_L_De(*this) == false)
//		return false;
//	LT.transpOf(L);
//	OP_mult(L,LT);
//	return true;
//}

//void L_Matrix::repararDefinidaPositivaSimetrica_LLT_forzado()
//{
//	L_Matrix L, LT;
//	L.cholesky_L_forzado_De(*this);
//	LT.transpOf(L);
//	OP_mult(L,LT);
//	return;
//}


#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<double> L_Matrix::crearUblas() const
{
	boost::numeric::ublas::matrix<double> mb(li,lj);
	memcpy(&(mb.data()[0]), data(), sizeof(double)*li*lj);
	return mb;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
boost::numeric::ublas::matrix<double> L_Matrix::crearUblasRef() const
{
	return boost::numeric::ublas::L_make_matrix_from_pointer(li, lj, data());
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
L_Matrix &L_Matrix::operator=(const boost::numeric::ublas::matrix<double> &other)
{
	reallocate((int)other.size1(), (int)other.size2());
	memcpy(data(), &(other.data()[0]), sizeof(double)*li*lj);
	return *this;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_UBLAS__
bool L_Matrix::invertMe_ublas()
{
	boost::numeric::ublas::matrix<double> input = crearUblasRef(), output(li, lj);
	bool ret;
	ret = boost::numeric::ublas::L_InvertMatrix(input, output);
	if (ret == false)
		return false;
	*this = output;
	return ret;
}
#endif // __COMPAT_UBLAS__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_Matrix::crearEigen() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	memcpy(mat.data(), data(), sizeof(double)*li*lj);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_Matrix::crearEigenLenta() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
	Eigen::MatrixXd mat(lI, lJ);
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			mat(i,j) = operator()(i,j);
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
Eigen::MatrixXd L_Matrix::crearEigenRef() const
{
	size_t lI=(size_t)li, lJ=(size_t)lj;
#if EIGEN_VERSION_AT_LEAST(2,90,0)
	Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor> > mat(const_cast<double*>(begin()), lI, lJ); // 1 = row major, revisar esto
#else
	Eigen::Map<Eigen::Matrix<double,10000,10000,Eigen::RowMajor> > mat(&elem[0][0], lI, lJ); // 1 = row major
#endif
	return mat;
}
#endif // __COMPAT_EIGEN__

#ifdef __COMPAT_EIGEN__
L_Matrix &L_Matrix::operator=(const Eigen::MatrixXd &other)
{
	Eigen::MatrixXd m;
	reallocate((int)other.rows(), (int)other.cols());
	//memcpy(&elem[0][0],  Eigen::Transpose<Eigen::MatrixXd>(other).data(), sizeof(double)*li*lj);
	m = Eigen::Transpose<Eigen::MatrixXd>(const_cast<Eigen::MatrixXd&>(other));
	memcpy(data(),  m.data(), sizeof(double)*li*lj);
	return *this;
}
#endif // __COMPAT_EIGEN__


void L_Matrix::fill(const double *lista)
{
	memcpy(begin(), lista, sizeof(double)*li*lj);
}

bool L_Matrix::readFile(FILE *fp)
{
	bool enEspacio=true;
	int i, j;
	L_Array<std::vector<double> > m;
	int maxnum = 0;

	L_String lin;  // El constructor de lin no permite reserve espacio
	int nMax = 300000;
	lin.resize(nMax);
	L_Array<L_String> nums;

	while(true)
	{
		if (fgets(lin.data(), nMax-2, fp) == NULL)
			break;
		lin.resize((int)( strlen(lin.c_str())+1 ));
		lin[lin.size()-1] = 0;
		lin.tokens(nums, " ,\t;\r");
		if (nums.size() == 0)
			continue;
		if (maxnum < (int)nums.size())
			maxnum = (int)nums.size();
		m.resize_swapping(m.size()+1);
		m[m.size()-1].resize(nums.size());
		for (i=0; i<(int)nums.size(); i++)
			sscanf(nums[i].c_str(), "%lf", &m[m.size()-1][i]);
	}

	reallocate(m.size(), maxnum);
	setZero();

	for (i=0; i<(int)m.size(); i++)
		for (j=0; j<(int)m[i].size(); j++)
			operator()(i,j) = m[i][j];
	return true;
}

bool L_Matrix::saveFile(FILE *fp)
{
	int i, j;

	for (i=0; i<li; i++)
	{
		for (j=0; j<lj-1; j++)
			fprintf(fp, "%lf ", operator()(i,j));
		if (j == lj-1)
			fprintf(fp, "%lf", operator()(i,j));
		fprintf(fp, "\n");
	}

	return true;
}

void L_Matrix::integralDebugOf(const L_Matrix &other)
{
	int i, j, u, v;
	double res;
	reallocate(other.li,other.lj);
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			res = 0;
			for (u=0; u<=i; u++)
				for (v=0; v<=j; v++)
					res += other(u,v);
			operator()(i,j) = res;
		}
	}
}

void L_Matrix::integralCuadradosDe_prueba(const L_Matrix &other)
{
	int i, j, u, v;
	reallocate(other.li,other.lj);
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			operator()(i,j) = 0;
			for (u=0; u<=i; u++)
				for (v=0; v<=j; v++)
					operator()(i,j) += other(u,v)*other(u,v);
		}
	}
}


void L_Matrix::repairSymmetry()
{
	int i, j;
	throw_L_ArgException_if (li != lj, "L_Matrix::repairSymmetry() : no cuadrada");
	for (i=0; i<li; i++)
		for (j=i+1; j<lj; j++)
			operator()(j,i) = operator()(i,j) = 0.5*(operator()(i,j)+operator()(j,i));
}

void L_Matrix::repairRangeValues(double limInf, double limMax)
{
	int i, j;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
		{
			if (operator()(i,j) < limInf)
				operator()(i,j) = limInf;
			if (operator()(i,j) > limMax)
				operator()(i,j) = limMax;
		}
}


double L_Matrix::covarianceThroughJacobianTest(const std::vector<int> &a, const std::vector<int> &b, const L_Matrix &J, const L_Matrix *Ji, double umbral, bool *aceptado, std::vector<L_StaticMatrix<3,3> > *cov3x3, L_Matrix *rot)
{
	typedef L_Matrix::size_type m_szt;
	L_Matrix Paa((m_szt)a.size(),(m_szt)a.size()), PaM((m_szt)a.size(),li-(m_szt)a.size()), PMa(li-(m_szt)a.size(),(m_szt)a.size()), PMM(li-(m_szt)a.size(),li-(m_szt)a.size()), PD;
	L_Matrix Pbb, PbN, PNb; // PNN = PMM
	L_Matrix JT, JPaa;
	L_Matrix M3x3(3,3);
	std::vector<int> M, N;
	m_szt i, j;
	double perdida = 1e60;

	L_complementOf(a, M, (int)li);
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

	// Paa: 3n x 3n
	// Pbb: 7 x 7

	// test, calculo de la perdida de covarianza (+- = informacion) ocurrida durante la transformacion
	if (Ji != NULL)
	{
		// double PAAn, Paan;
		L_Matrix PAA, PAM;    // Reconstruccion de Paa, PaM
		L_Matrix JiT, PbbJiT; // Temporales
		//

		L_Matrix JiJ;
		JiJ.OP_mult(*Ji,J);

		// PAA = Ji*Pbb*JiT
		JiT.transpOf(*Ji);
		PbbJiT.OP_mult(Pbb,JiT);
		PAA.OP_mult(*Ji,PbbJiT);

		// PAAn = PAA.normL2();
		// Paan = Paa.normL2();

		// Se supone que el landmark 3D contiene la covarianza conjunta de los landmarks
		// Falta extraer la covarianza individual de cada landmark y guardarla aparte
		// Es una opcion...

		if (cov3x3 != NULL)
		{
			// Sacar covarianzas individuales
			// ESTA PARTE HAY QUE REVISARLA MAS
			int u, v;
			cov3x3->resize(a.size()/3);
			M3x3.reallocate(3,3);
			for (i=0; i < (int)cov3x3->size(); i++)
			{
				for (u=0; u<3; u++)
					for (v=0; v<3; v++)
						M3x3(u,v) = Paa(3*i+u,3*i+v) - PAA(3*i+u,3*i+v);
				if (M3x3(0,0) > 0 || M3x3(1,1) > 0 || M3x3(2,2) > 0)
				{
					for (u=0; u<3; u++)
						M3x3(u,u) = M3x3(u,u) > 0 ? M3x3(u,u) : 0;  // decia [u,v)
					M3x3.repairSymmetry();
					M3x3.repairPositiveDefinite_valPr();
				}
				else
					M3x3.identity(1e-10);
				L_StaticMatrix_OP_assign((*cov3x3)[i], M3x3);
			}
			// Sumar las covarianzas individuales a la covarianza aproximada
			for (i=0; i < (int)cov3x3->size(); i++)
			{
				for (u=0; u<3; u++)
					for (v=0; v<3; v++)
						PAA(3*i+u,3*i+v) += (*cov3x3)[i](u,v);
			}
		}
		PD.OP_subtract(PAA, Paa);
		perdida = PD.normL2();
		if (perdida > umbral && umbral != 0)
		{
			if (aceptado != NULL)
				*aceptado = false;
			return perdida; // Hasta aca llegamos, no se transformUsing la covarianza
		}
	}

	reallocate((m_szt)(b.size() + N.size()), (m_szt)(b.size() + N.size()));

	for (i=0; i<(m_szt)b.size(); i++)
		for (j=0; j<(m_szt)b.size(); j++)
			operator()(b[i],b[j]) = Pbb(i,j);
	for (i=0; i<(m_szt)b.size(); i++)
		for (j=0; j<(m_szt)N.size(); j++)
			operator()(b[i],N[j]) = PbN(i,j);
	for (i=0; i<(m_szt)N.size(); i++)
		for (j=0; j<(m_szt)b.size(); j++)
			operator()(N[i],b[j]) = PNb(i,j);
	for (i=0; i<(m_szt)N.size(); i++)
		for (j=0; j<(m_szt)N.size(); j++)
			operator()(N[i],N[j]) = PMM(i,j);

	if (aceptado != NULL)
		*aceptado = true;

	return perdida;
}


void L_Matrix::covarianceThroughJacobianTestNew(const std::vector<int> &a, const std::vector<int> &b, const L_Matrix &J, std::vector<L_StaticMatrix<3,3> > &cov3x3, bool forzar_maxima_cov)
{
	typedef L_Matrix::size_type m_szt;
	L_Matrix Paa((m_szt)a.size(),(m_szt)a.size()), PaM((m_szt)a.size(),(m_szt)(li-a.size())), PMa((m_szt)(li-a.size()),(m_szt)a.size()), PMM((m_szt)(li-a.size()),(m_szt)(li-a.size())), PD;
	L_Matrix Pbb, PbN, PNb; // PNN = PMM
	L_Matrix JT, JPaa;
	L_Matrix M3x3(3,3);
	std::vector<int> M, N;
	m_szt i, j;
	// double perdida = 1e60;


	L_complementOf(a, M, (int)li);
	L_complementOf(b, N, (int)(li+b.size()-a.size()));

	for (i=0; i<(m_szt)a.size(); i++)
		for (j=0; j<(m_szt)a.size(); j++)
			Paa(i,j) = operator()(a[i],a[j]);
	for (i=0; i<(m_szt)a.size(); i++)
		for (j=0; j<(m_szt)M.size(); j++)
			PaM(i,j) = operator()(a[i],M[j]);
	for (i=0; i<(m_szt)M.size(); i++)
		for (j=0; j<(m_szt)a.size(); j++)
			PMa(i,j) = operator()(M[i],a[j]);
	for (i=0; i<(m_szt)M.size(); i++)
		for (j=0; j<(m_szt)N.size(); j++)
			PMM(i,j) = operator()(M[i],M[j]);
	
	if (Paa.isGoodCovarianceMatrix() == false)
	{
		printf("covarianceThroughJacobianTestNew() : submatriz original no simetrica definida positiva (1)\n");
		Paa.repairPositiveDefinite_valPr();
	}

	L_TransfCov_tresMatrices tres;
	L_Matrix Pbodypoints;
	int u;

	tres.encontrarCovarianzasIndividuales_vect(Paa, Pbodypoints);

	L_Matrix miniCov(3,3);
	cov3x3.resize(Paa.li/3);
	for (u=0; u<Paa.li/3; u++)
	{
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				miniCov(i,j) = Pbodypoints(3*u+i,3*u+j);
		miniCov *= 0.95;
		miniCov.repairPositiveDefinite_valPr();
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				Pbodypoints(3*u+i,3*u+j) = miniCov(i,j);
		L_StaticMatrix_OP_assign(cov3x3[u], miniCov);
	}

	Paa.OP_subtract(Paa, Pbodypoints); // Reduccion de Paa para sacarle la cov de los body points
	if (Paa.isGoodCovarianceMatrix() == false)
		Paa.repairPositiveDefinite_valPr();

	// Truco: dejar la maxima covarianza individual como value final
	if (forzar_maxima_cov)
		for (u=0; u<Paa.li/3; u++)
			for (i=0; i<3; i++)
				for (j=0; j<3; j++)
					Pbodypoints(3*u+i,3*u+j) = Paa(i,j);

	// Supuestamente, la reduccion de Paa no afecta las covarianzas cruzadas PaM
	JT.transpOf(J);

	JPaa.OP_mult(J,Paa);
	Pbb.OP_mult(JPaa,JT);
	PbN.OP_mult(J,PaM);
	PNb.OP_mult(PMa,JT);

	reallocate((m_szt)(b.size() + N.size()), (m_szt)(b.size() + N.size()));

	for (i=0; i<(m_szt)b.size(); i++)
		for (j=0; j<(m_szt)b.size(); j++)
			operator()(b[i],b[j]) = Pbb(i,j);
	for (i=0; i<(m_szt)b.size(); i++)
		for (j=0; j<(m_szt)N.size(); j++)
			operator()(b[i],N[j]) = PbN(i,j);
	for (i=0; i<(m_szt)N.size(); i++)
		for (j=0; j<(m_szt)b.size(); j++)
			operator()(N[i],b[j]) = PNb(i,j);
	for (i=0; i<(m_szt)N.size(); i++)
		for (j=0; j<(m_szt)N.size(); j++)
			operator()(N[i],N[j]) = PMM(i,j);

	if (isGoodCovarianceMatrix() == false)
	{
		printf("Se ensucio con el cambio de representacion (1)\n");
		L_Matrix L, JL;
		if (L.cholesky_L_De(Paa) == true)
		{
			JL.OP_mult(J, L);
			Pbb.OP_mult_ABT(JL, JL);
			for (i=0; i<(int)b.size(); i++)
				for (j=0; j<(int)b.size(); j++)
					operator()(b[i],b[j]) = Pbb(i,j);
		}
		else
		{
			printf("covarianceThroughJacobianTestNew() : submatriz original no simetrica definida positiva (2)\n");
			Paa.repairPositiveDefinite_valPr();
			if (L.cholesky_L_De(Paa) == true)
			{
				JL.OP_mult(J, L);
				Pbb.OP_mult_ABT(JL, JL);
				for (i=0; i<(int)b.size(); i++)
					for (j=0; j<(int)b.size(); j++)
						operator()(b[i],b[j]) = Pbb(i,j);
			}
		}
	}
	if (isGoodCovarianceMatrix() == false)
	{
		printf("Se ensucio con el cambio de representacion (2)\n");
		this->repairPositiveDefinite_valPr();
	}
}


double L_Matrix::trace()
{
	throw_L_ArgException_if(li!=lj, "L_Matrix::trace() : no cuadrada");
	int i;
	double trace=0;
	for (i=0; i<li; i++)
		trace+=operator()(i,i);
	return trace;
}

double L_Matrix::traceRange(int i1, int i2)
{
	throw_L_ArgException_if(li!=lj, "L_Matrix::trace() : no cuadrada");
	int i;
	double trace=0;
	throw_L_ArgException_if(i1 < 0 || i2 < 0 || i1 >= li || i2 >= li || i1 >= lj || i2 >= lj || i1 > i2, "L_Matrix::traceRange()");
	for (i=i1; i<=i2; i++)
		trace+=operator()(i,i);
	return trace;
}

double L_Matrix::det_v0()
{
	throw_L_ArgException_if(li!=lj, "L_Matrix::det_v0() : no cuadrada");
	if (li==1)
		return operator()(0,0);
	else if (li==2)
		return operator()(0,0)*operator()(1,1)-operator()(0,1)*operator()(1,0);
	double d=0;
	int i;
	L_Matrix sub;
	int signo=1;
	for (i=0; i<li; i++)
	{
		subMatrDeleteRowCol_at(sub,i,0);
		d+=operator()(i,0)*sub.det()*signo;
		signo*=-1;
	}
	return d;
}


void L_Matrix::genDrawing(L_ImageRGBUchar &im)
{
	int i, j;
	double mx, prm=0, a, b;
	throw_L_ArgException_if(li==0 || lj==0 || begin()==NULL, "L_Matrix::genDrawing() : imagen nula");
	im.reallocate(lj, li);
	mx = operator()(0,0);

	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			prm += L_ABS(operator()(i,j));
			if (L_ABS(operator()(i,j)) > mx)
				mx = L_ABS(operator()(i,j));
		}
	}
	prm/=(li*1.0*lj);

	// 0 -> 0
	// prm -> 128
	// max -> 255
	// a*x^2 + b*x = y
	// a*128^2 + b*128 = prm  =>  a = prm/(128^2) - b/128
	// a*255^2 + b*255 = max => 255^2*prm/(128^2) - 255^2*b/128 + b*255 = max  => b*(255 - 255^2/128)  = max - 255^2*prm/(128^2)
	// => b = (max - 255*255*prm/(128*128)) / (255 - 255*255 / 128) , a = prm/(128*128) - b/128

	b = (mx - 255.0*255.0*prm/(128.0*128.0)) / (255.0 - 255.0*255.0 / 128.0);
	a = prm/(128.0*128.0) - b/128.0;


	for (j=0; j<lj; j++)
	{
		for (i=0; i<li; i++)
		{
			if (operator()(i,j) > 0)
			{
				im.pix(j,i,0) = 0;
				im.pix(j,i,1) = 0;
				im.pix(j,i,2) = (L_uchar)( a*operator()(i,j)*operator()(i,j) + b*operator()(i,j) );
			}
			else
			{
				im.pix(j,i,0) = 0;
				im.pix(j,i,1) = (L_uchar)( a*(-operator()(i,j))*(-operator()(i,j)) + b*(-operator()(i,j)) );
				im.pix(j,i,2) = 0;
			}
		}
	}
}


void L_Matrix::OP_mult_Abstr(const L_Matrix *m1, const L_Matrix *m2) // Para objetos polimorficos...
{
	int i,j,u;
	double res;

	throw_L_ArgException_if(m1->lj!=m2->li, "L_Matrix::OP_mult_Abstr() : tamanos inconcordantes");
	if (begin()==NULL || li!=m1->li || lj!=m2->lj)
		reallocate(m1->li, m2->lj);

	for (i=0; i<m1->li; i++) for (j=0; j<m2->lj; j++)
	{
		res=0;
		for (u=0; u<m2->li; u++)
			res+=m1->operator()(i,u)*m2->operator()(u,j);
		operator()(i,j) = res;
	}
	return;
}

void L_Matrix::OP_div(const L_Matrix &m1, const L_Matrix &m2) // m1*(m2^-1)
{
	L_Matrix divisor;
	divisor=m2;
	if (m1.lj==1)
	{
		*this = m1;
		if (!solveLinearSystem_x_Ab_result_in_b(divisor))
			throw L_ZeroDivisionException();
	}
	else
	{
		if (!divisor.invertMe())
			throw L_ZeroDivisionException();
		(*this).OP_mult(m1, divisor);
	}
	return;
}




void L_Matrix::OP_div_NULL(const L_Matrix &other)
{
	L_Matrix aux;
	aux=other;

	throw_L_ArgException_if(other.begin() == NULL, "L_Matrix::OP_div_NULL() : division por nula"); // Division por nulo es error siempre
	if (begin() == NULL)
		return; // zero dividido por algo es zero

	if (lj == 1)
	{
		if (!solveLinearSystem_x_Ab_result_in_b(aux))
			throw L_ZeroDivisionException(); // No es L_ArgException
		return;
	}
	if (!aux.invertMe())
		throw L_ZeroDivisionException();
	(*this)*=aux;
	return;
}



void L_Matrix::OP_mult_Abstr_NULL(const L_Matrix *m1, const L_Matrix *m2)
{
	int i_m1,j_m2,u;
	double res;

	throw_L_ArgException_if(m1->lj!=m2->li, "L_Matrix::OP_mult_Abstr_NULL() : tamanos incorcondantes");
	if (m1->begin() == NULL || m2->begin() == NULL)
	{
		destroy();
		li = m1->li;
		lj = m2->lj;
		return;
	}
	if (begin()==NULL || li!=m1->li || lj!=m2->lj)
		reallocate(m1->li, m2->lj);

	for (i_m1=0; i_m1<m1->li; i_m1++) for (j_m2=0; j_m2<m2->lj; j_m2++)
	{
		res = 0;
		for (u=0; u<m2->li; u++)
			res+=m1->operator()(i_m1,u) * m2->operator()(u,j_m2);
		operator()(i_m1,j_m2) = res;
	}
	return;
}

void L_Matrix::OP_div_NULL(const L_Matrix &m1, const L_Matrix &m2)
{
	L_Matrix divisor;
	divisor=m2;
	throw_L_ArgException_if(m2.begin() == NULL, "L_Matrix::OP_div_NULL() : division por nula");
	if (m1.begin() == NULL)
	{
		destroy();
		li = m1.li;
		lj = m2.lj;
		return;
	}
	if (m1.lj==1)
	{
		*this = m1;
		if (!solveLinearSystem_x_Ab_result_in_b(divisor))
			throw L_ZeroDivisionException();
	}
	else
	{
		if (!divisor.invertMe())
			throw L_ZeroDivisionException();
		(*this).OP_mult(m1, divisor);
	}
	return;
}

bool L_Matrix::invertMegj()
{
	bool v;
	L_Matrix_intercambia a, b;
	a.OP_assign(*this);
	b.reallocate(li, lj);
	b.identity();
	v = a.gaussj(b);
	if (v)
		b.copyTo(*this);
	return v;
}

bool L_Matrix::resolverSisLinealgj(const L_Matrix &A, const L_Matrix &b)
{
	bool v;
	L_Matrix_intercambia aI, bI;
	aI.OP_assign(A);
	bI.OP_assign(b);
	v = aI.gaussj(bI);
	if (v)
		bI.copyTo(*this);
	return v;
}

void L_Matrix::prueba_inversas_veloc()
{
	L_Matrix A(20,20), AT, H, I, J;
	A.fillRandomValues(-100, 100);
	double t1, t2, t3, t4, t5;
	int i;

	AT.transpOf(A);  // Crear matriz simetrica definida positiva
	H.OP_mult(A, AT);

	t1 = L_TIME();
	I.OP_assign(H);
	for (i=0; i<300; i++)
		I.invertMe();
	t2 = L_TIME();
	I.OP_assign(H);
	for (i=0; i<300; i++)
		I.invertMegj();
	t3 = L_TIME();
	I.OP_assign(H);
	for (i=0; i<300; i++)
		I.invertMeDefPos();
	t4 = L_TIME();
	I.OP_assign(H);
	for (i=0; i<300; i++)
		I.invertMeSimetricaDefPos();
	t5 = L_TIME();
	printf("mia: %g    recipes: %g   posdef: %g   posdefsim: %g\n" , t2-t1, t3-t2, t4-t3, t5-t4);
	getchar();
}

void L_Matrix::espacioNuloPivoteo(const L_Matrix &original_ancha)
{
	L_Matrix m;
	m = original_ancha; // La input es de h x (h+d) ; la output es de (h+d) x d
	int i, j, d = m.lj - m.li, h = m.li;
	m.diagonalize_me(); // Queda [ Ihxh | v1 ... vd ]
	reallocate(h+d, d); // Debe quedar [-v1 ... -vd ; Idxd]
	for (i=0; i<h; i++)
		for (j=0; j<d; j++)
			operator()(i,j) = -m(i,h+j); // Copiar [-v1 ... -vd ;
	for (i=0; i<d; i++)
		for (j=0; j<d; j++)
			operator()(h+i,j) = (i==j) ? 1 : 0; // Copiar  ; Idxd)
}

void L_Matrix::gramSchmidtColumnas()
{
	int i, j, w;
	double p;
	for (j=0; j<lj; j++)
	{
		// Restar proyeccion sobre columnas anteriores
		for (w=0; w<j; w++)
		{
			p=0;
			for (i=0; i<li; i++)
				p+=operator()(i,j)*operator()(i,w);
			for (i=0; i<li; i++)
				operator()(i,j)-=operator()(i,w)*p;

		}
		// Normalizar la columna actual
		p=0;
		for (i=0; i<li; i++)
			p+=operator()(i,j)*operator()(i,j);
		p=sqrt(p);
		for (i=0; i<li; i++)
			operator()(i,j)/=p;

	}
}

void L_Matrix::prueba_inversas()
{
	L_Matrix A, Ainv, b, xa, xb;
	int i, j, vez;
	bool ret1, ret2;

	A.reallocate(10,10);
	b.reallocate(10,1);

	// test inversas
	for (vez = 0; vez < 200; vez ++)
	{
		printf("%d\t", vez);
		for (i=0; i<10; i++)
		{
			for (j=0; j<10; j++)
			{
				A(i,j)=(rand()%1000-500)*0.1;
				if (vez > 100 && i==j)
					A(i,j)=0;
			}
		}
		for (i=0; i<10; i++)
			b(i,0)=(rand()%1000-500)*0.1;
		Ainv=A;
		ret1 = Ainv.invertMe();
		if (ret1)
			xa.OP_mult(Ainv,b);
		ret2 = xb.solveLinearSystem(A, b);
		if (ret1 != ret2)
		{
			printf(" true/false incompatibles\n");
			continue;
		}
		if (ret1 == false)
			continue;
		for (i=0; i<10; i++)
		{
			if (fabs(xa(i,0) - xb(i,0)) > 0.00001)
				printf(" Error en calculos de despeje de matrices\n");
			break;
		}
	}
	printf("\nPresione ENTER para salir\n");
	getchar();
}


void L_Matrix::print(const char *name, FILE *fp) const
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

	for (i=0; i<li; i++)
	{
		if (str.size() != 0)
		{
			if (imprNom==false && i>=li/2)
			{
				fprintf(fp, "%s", str.c_str());
				imprNom=true;
			}
			else
				fprintf(fp, "%s", str2.c_str());
		}
		fprintf(fp,"[");
		for (j=0; j<lj; j++)
			fprintf(fp, "%+.3g ", operator()(i,j));
		fprintf(fp,"]\n");
	}
	fprintf(fp,"\n");
}

void L_Matrix::imprime_arreglo(double *v, int n, const char *name, FILE *fp)
{
	int i;
	fprintf(fp, "%s = {", name);
	for (i=0; i<n; i++)
		fprintf(fp, "%.2f ", v[i]);
	fprintf(fp,"}\n");
}
void L_Matrix::printMatlabFormat(FILE *fp, const char *name) const
{
	int i, j;
	fprintf(fp, "%s = [", name);
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
		{
			if (operator()(i,j) == operator()(i,j))
				fprintf(fp, "%.10g, ", operator()(i,j));
			else
				fprintf(fp, "NaN, ");
		}
		fprintf(fp,";");
	}
	fprintf(fp,");\n");
}

double L_Matrix::debugValue(int i, int j)
{
	return operator()(i,j);
}

void L_Matrix::prueba_det()
{
	L_Matrix m;
	long vez, nMax;
	int li, i;
	double det1, det2;

	nMax = 100000;
	for (vez=0; vez<nMax; vez++)
	{
		if (vez*100/nMax % 10 == 0)
			printf("%ld%%\n", vez*100/nMax);
		li = (int)L_RANDOM(5, 20);
		m.reallocate((int)li, (int)li);
		m.fillRandomValues(-100, 100);
		for (i=0; i<li*li; i++)
			if (fabs(m(i/100,i%100)) < 1 && L_RANDOM(-1,1) > 0)
				m(i/100,i%100) = 0;
		det1 = m.det_v0();
		det2 = m.det();

	}
}

void L_Matrix::prueba_svd(double error)
{
	L_Matrix u;
	u=*this;
	L_Matrix d, v, vT, F, G;
	int i, j, k;
	double sum;
	bool converge;
	converge=u.svd(d, v);
	if (converge==false)
	{
		printf("SVD no converge\n");
		return;
	}
	F.OP_multByDiagonal(u,d);
	vT.transpOf(v);
	G.OP_mult(F,vT);
	print("Original");
	G.print("Reconstruida");
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			sum=0;
			for (k=0; k<3; k++)
			{
				sum+=u(i,k)*u(j,k);
			}
			if (i==j)
			{
				if (fabs(sum-1.0) > error)
					printf("matriz u de SVD no es ortonormal: error %f\n", fabs(sum-1.0));
			}
			else
			{
				if (fabs(sum) > error)
					printf("matriz u de SVD no es ortonormal: error %f\n", fabs(sum));
			}
		}
	}
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			sum=0;
			for (k=0; k<3; k++)
			{
				sum+=v(i,k)*v(j,k);
			}
			if (i==j)
			{
				if (fabs(sum-1.0) > error)
					printf("matriz v de SVD no es ortonormal: error %f\n", fabs(sum-1.0));
			}
			else
			{
				if (fabs(sum) > error)
					printf("matriz v de SVD no es ortonormal, error %f\n", fabs(sum));
			}
		}
	}

}

void L_Matrix::prueba_svd_ordenado(double error)
{
	L_Matrix u;
	u=*this;
	L_Matrix d, v, vT, F, G;
	int i, j, k;
	double sum;
	bool converge;
	converge=u.svd_ordenado(d, v);
	if (converge==false)
	{
		printf("SVD no converge\n");
		return;
	}
	for (i=1; i<d.li; i++)
		if (d(i,0) < d(i-1,0))
			printf("Error en el orden de la diagonal\n");
	for (i=0; i<d.li; i++)
		if (d(i,0) < 0)
			printf("Matriz con elemento diagonal negativo\n");

	F.OP_multByDiagonal(u,d);
	vT.transpOf(v);
	G.OP_mult(F,vT);
	print("Original");
	G.print("Reconstruida");
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			sum=0;
			for (k=0; k<3; k++)
			{
				sum+=u(i,k)*u(j,k);
			}
			if (i==j)
			{
				if (fabs(sum-1.0) > error)
					printf("matriz u de SVD no es ortonormal, error %f\n", fabs(sum-1.0));
			}
			else
			{
				if (fabs(sum) > error)
					printf("matriz u de SVD no es ortonormal, error %f\n", fabs(sum));
			}
		}
	}
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			sum=0;
			for (k=0; k<3; k++)
			{
				sum+=v(i,k)*v(j,k);
			}
			if (i==j)
			{
				if (fabs(sum-1.0) > error)
					printf("matriz v de SVD no es ortonormal, error %f\n", fabs(sum-1.0));
			}
			else
			{
				if (fabs(sum) > error)
					printf("matriz v de SVD no es ortonormal, error %f\n", fabs(sum));
			}
		}
	}

}


void L_Matrix::pruebasvd_cholesky_vectpr()
{
	L_Matrix A, AT, ATA, AL, AU, ALAU, D;
	L_Matrix uns;
	A.reallocate(11,11);
	A.fillRandomValues(-0.5, 0.5);

	AT.transpOf(A);
	ATA.OP_mult(AT, A);

	AL.cholesky_L_De(ATA);
	AU.transpOf(AL);

	A.print("A");
	AL.print("AL");
	ATA.print("ATA");
	ALAU.OP_mult(AL, AU);
	ALAU.print("AL*AU");

	L_Matrix U, S, V, P, PT, S2;

	U = AT;
	U.svd_ordenado(S, V);
	V.transposeMe();

	P = ATA;
	P.vectPr_sim(D);
	PT = P;
	PT.transposeMe();

	U.ordenaDiagonal(S,V);
	P.ordenaDiagonal(D,PT);

	printf("ATA->valPr\n");
	P.print("P");
	D.print("D");
	printf("AT->svd\n");
	U.print("U");
	S.print("S");
	S2.OP_productElementwise(S,S);
	S2.print("S2");
	getchar();

	printf("Presione enter\n");
	getchar();
}

void L_Matrix::vectPr_sim(L_Matrix &valPr)
{
#if defined(__COMPAT_EIGEN__) && defined(OJO______USE_OPTIMIZED_EIGENDECOMPOSITION)
	return vectPr_Eigen_sim(valPr); // No logro hacer que funcione
#else
	return vectPr_jacobi_sim(valPr);
#endif
}


void L_Matrix::vectPr_jacobi_sim(L_Matrix &valPr)
{
	throw_L_ArgException_if(li == 0 ||  li != lj || begin() == NULL, "L_Matrix::vectPr_sim() : tamano inadecuado");

	L_Matrix A;
	std::vector<double> val;
	bool columna = true;
	val.resize(li);
	A=*this;
	Jacobi_Cyclic_Method(&(val[0]), begin(), A.begin(), li);

	if (columna)
		valPr.diagonalVectorOf(&(val[0]), li);
	else
		valPr.diagonalMatrixOf(&(val[0]), li);
}

bool L_Matrix::interseccionRayos(L_Matrix &p1, L_Matrix &v1, L_Matrix &p2, L_Matrix &v2, double &dis1, double &dis2)
{
	// Calcula la interseccion de los rayos: r1 = p1 + dis1*v1  y  r2 = p2 + dis2*v2
	double v1x, v1y, v1z, v2x, v2y, v2z, p1x, p1y, p1z, p2x, p2y, p2z;
	double mod1, mod2, v1pv2, den, p2_p1x, p2_p1y, p2_p1z;
	v1x = v1(0,0);
	v1y = v1(1,0);
	v1z = v1(2,0);
	v2x = v2(0,0);
	v2y = v2(1,0);
	v2z = v2(2,0);
	p1x = p1(0,0);
	p1y = p1(1,0);
	p1z = p1(2,0);
	p2x = p2(0,0);
	p2y = p2(1,0);
	p2z = p2(2,0);

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
	p2_p1x = p2x - p1x;
	p2_p1y = p2y - p1y;
	p2_p1z = p2z - p1z;
	dis1 = ( p2_p1x*v1x + p2_p1y*v1y + p2_p1z*v1z  - ( p2_p1x*v2x + p2_p1y*v2y + p2_p1z*v2z )*v1pv2
		) / den;
	dis2 = (-p2_p1x*v2x +-p2_p1y*v2y +-p2_p1z*v2z  - (-p2_p1x*v1x +-p2_p1y*v1y +-p2_p1z*v1z )*v1pv2
		) / den;
	reallocate(3,1);
	operator()(0,0) = p1x + v1x*dis1;
	operator()(1,0) = p1y + v1y*dis1;
	operator()(2,0) = p1z + v1z*dis1;
	return true;
}





// http://www.netlib.org/blas/zrotg.f
void L_Matrix::blas_zrotg(double &CA, double &CB, double &C, double &S) // SUBROUTINE ZROTG(CA,CB,C,S)
{
	double ALPHA; //      DOUBLE COMPLEX ALPHA
	double NORM, SCALE; //      DOUBLE PRECISION NORM,SCALE

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

// Se apply por columnas, igual que en fortran
// http://www.netlib.org/blas/dnrm2.f
double L_Matrix::blas_dnrm2(int N, const L_Matrix &X, int i0, int j0, int INCX)
{
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
		for (IX=0; IX<=(N-1)*INCX; IX+=INCX)//DO 10 IX+1 = 1,1 + (N-1)*INCX,INCX
		{
			if (X(IX+i0,j0) != 0) //IF (DBLE(X(IX+1)).NE.ZERO) THEN
			{
				TEMP = fabs(X(IX+i0,j0));//TEMP = ABS(DBLE(X(IX+1)))
				if (SCALE < TEMP) //IF (SCALE.LT.TEMP) THEN
				{
					SSQ = 1 + SSQ* (SCALE/TEMP)*(SCALE/TEMP);
					SCALE = TEMP;
				}
				else //ELSE
					SSQ = SSQ + (TEMP/SCALE)*(TEMP/SCALE);
				//END IF
			}// END IF
			// IF (DIMAG(X(IX+1)).NE.ZERO) THEN ... ENDIF -> aca no hay imaginarios
		} // 10     CONTINUE
		NORM = SCALE*sqrt(SSQ);
	}//END IF

	return NORM;
}

// Se apply por columnas, igual que en fortran
// http://www.netlib.org/blas/zdotc.f
double L_Matrix::blas_zdotc(int N, const L_Matrix &ZX, int i0, int j0, int INCX, const L_Matrix &ZY, int i1, int j1, int INCY)
{
	int I, IX, IY;
	double ZTEMP = 0;
	if (N == 0)  // IF (N.LE.0) RETURN
		return 0;
	if (INCX < 0)
		IX = (-N+1)*INCX;
	if (INCY < 0)
		IY = (-N+1)*INCY;

	if (INCX != 1 || INCY != 1) // IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
	{
		for (I=0; I<N; I++) // DO 10 I+1 = 1,N
		{
			ZTEMP += ZX(IX+i0,j0) * ZY(IY+i1,j1);  // ZTEMP = ZTEMP + DCONJG(ZX(IX+1))*ZY(IY+1)
			IX += INCX;
			IY += INCY;
		} // 10 CONTINUE
	}
	else
	{
		for (I=0; I<N; I++)  // 20   DO 30 I+1 = 1,N
			ZTEMP += ZX(I+i0,j0)*ZY(I+i1,j1);  // ZTEMP = ZTEMP + DCONJG(ZX(I+1))*ZY(I+1)
	} // 30 CONTINUE
	return ZTEMP;  // ZDOTC = ZTEMP    RETURN      END
}



// http://www.netlib.org/linpack/zchud.f
void L_Matrix::linpack_zchud(L_Matrix &r, int ldr, int p, const L_Matrix &x,L_Matrix &z, int ldz, int nz, L_Matrix &y, L_Matrix &rho, L_Matrix &c, L_Matrix &s)
{
	int i, j, jm1;
	double azeta, scale;
	double t, xj, zeta;  // complex*16

	// c   update r.

	for (j=0; j<p; j++)//do 30 j+1 = 1, p
	{
		xj = x(j,0);// xj = x(j+1)

		// c      apply the previous rotations.

		jm1 = j;
		if (jm1 >= 1)  //   if (jm1 .lt. 1) go to 20
		{
			for (i=0; i<jm1; i++)// do 10 i+1 = 1, jm1
			{
				t = c(i,0)*r(i,j) + s(i,0)*xj; // t = c(i+1)*r(i+1,j+1) + s(i+1)*xj
				xj = c(i,0)*xj - s(i,0)*r(i,j); // xj = c(i+1)*xj - dconjg(s(i+1))*r(i+1,j+1), todo es real aca
				r(i,j) = t; // r(i+1,j+1) = t
			} //	10    continue
		} //20    continue

	// c      compute the next rotation.

		blas_zrotg(r(j,j),xj,c(j,0),s(j,0));   // call zrotg(r(j+1,j+1),xj,c(j+1),s(j+1))
	} // 30 continue

	// c   if required, update z and rho.

	if (nz >= 1) //if (nz .lt. 1) go to 70
	{
		for (j=0; j<nz; j++)//do 60 j+1 = 1, nz
		{
			zeta = y(j,0); // zeta = y(j+1)
			for (i=0; i<p; i++)// do 40 i+1 = 1, p
			{
				t = c(i,0)*z(i,j) + s(i,0)*zeta; // t = c(i+1)*z(i+1,j+1) + s(i+1)*zeta  // Error: decia c(i,j)
				zeta = c(i,0)*zeta - s(i,0)*z(i,j); // zeta = c(i+1)*zeta - dconjg(s(i+1))*z(i+1,j+1)
				z(i,j) = t;  // z(i+1,j+1) = t
			} // 40    continue
			azeta = fabs(zeta);  // azeta = cdabs(zeta)
			if (azeta != 0.0 && rho(j,0) >= 0)  // //if (azeta .eq. 0.0d0 .or. rho(j+1) .lt. 0.0d0) go to 50
			{
				scale = azeta + rho(j,0); // scale = azeta + rho(j+1)
				double azs = azeta/scale;
				double rhs = rho(j,0)/scale;
				rho(j,0) = scale*sqrt(azs*azs+rhs*rhs); // rho(j+1) = scale*dsqrt((azeta/scale)**2+(rho(j+1)/scale)**2)
			} // 50    continue
		}//	60 continue
	} // 70 continue
	return;
} // end


void L_Matrix::linpack_zchdd(L_Matrix &r, int ldr, int p, const L_Matrix &x, L_Matrix &z, int ldz, int nz, L_Matrix &y, L_Matrix &rho, L_Matrix &c, L_Matrix &s, int &info)
{
	int i,ii,j; //, iter = 0;
	double a,alpha,azeta,norm;
	double t,zeta,b,xx; // complex*16
	double scale; // No se por que no estaba declarada
	//complex*16 zdumr,zdumi   //complex*16
	//dreal(zdumr) = zdumr
	//dimag(zdumi) = (0.0d0,-1.0d0)*zdumi

// c   solve the system ctrans(r)*a = x, placing the result
// c   in the array s.

	info = 0;
	s(0,0) = x(0,0) / r(0,0); //s(1) = dconjg(x(1))/dconjg(r(1,1))
	if (p >= 2) //if (p .lt. 2) go to 20
	{
		for (j=1; j<p; j++)//do 10 j+1 = 2, p
		{
			s(j,0) = x(j,0) - blas_zdotc(j, r, 0,j, 1, s, 0,0, 1); //s(j+1) = dconjg(x(j+1)) - zdotc(j,r(1,j+1),1,s,1)
			s(j,0) = s(j,0)/r(j,j); // s(j+1) = s(j+1)/dconjg(r(j+1,j+1))
		} //10 continue
	} //20 continue
	norm = blas_dnrm2(p, s, 0,0, 1); // dznrm2(p,s,1)

	if (norm >= 1.0e0)   // if (norm .lt. 1.0d0) go to 30
	{
		info = -1;
		return;  // go to 120
	}

	// 30 continue
	alpha = sqrt(1.0e0 - norm*norm);// dsqrt(1.0d0-norm**2)

	// c      determine the transformations.

	for (ii=0; ii<p ; ii++)//do 40 ii+1 = 1, p
	{
		i = p - ii - 1;
		scale = alpha + fabs(s(i,0)); // scale = alpha + cdabs(s(i+1))
		a = alpha/scale;
		b = s(i,0)/scale;// b = s(i+1)/scale
		norm = sqrt(a*a+b*b); // norm = dsqrt(a**2+dreal(b)**2+dimag(b)**2)
		c(i,0) = a/norm; // c(i+1) = a/norm
		s(i,0) = b/norm; //s(i+1) = dconjg(b)/norm
		alpha = scale*norm;
	} // 40    continue

	// c      apply the transformations to r.

	for (j=0; j<p; j++)//do 60 j+1 = 1, p
	{
		xx = 0.0; //(0.0d0,0.0d0)
		for (ii=0; ii<j+1; ii++)  //do 50 ii+1 = 1, j+1
		{
			i = j - ii-1 + 1;
			t = c(i,0)*xx + s(i,0)*r(i,j); // t = c(i+1)*xx + s(i+1)*r(i+1,j+1)
			r(i,j) = c(i,0) * r(i,j) - s(i,0)*xx; // r(i+1,j+1) = c(i+1)*r(i+1,j+1) - dconjg(s(i+1))*xx   // Error: decia * en vez de -
			xx = t;
		}  // 50       continue
	} //60    continue

	// c      if required, downdate z and rho.

	if (nz >= 1) // if (nz .lt. 1) go to 110
	{
		for (j=0; j<nz; j++) //do 100 j+1 = 1, nz
		{
			zeta = y(j,0); //zeta = y(j+1)
			for (i=0; i<p; i++)  // do 70 i+1 = 1, p
			{
				z(i,j) = (z(i,j) - s(i,0)*zeta) / c(i,0); //z(i+1,j+1) = (z(i+1,j+1) - dconjg(s(i+1))*zeta)/c(i+1)  // Error: faltaba *zeta  // Error: decia c(i+1,0)
				zeta = c(i,0)*zeta - s(i,0)*z(i,j); // zeta = c(i+1)*zeta - s(i+1)*z(i+1,j+1)
			} // 70       continue
			azeta = fabs(zeta); //azeta = cdabs(zeta)
			if (azeta >= rho(j,0))// if (azeta .le. rho(j+1)) go to 80    // Error: puse <
			{
				// OJO ACA
				double adentro = 1.0e0-(azeta/rho(j,0))*(azeta/rho(j,0));
				if (adentro < 0)
					adentro = 0;
				rho(j,0) = rho(j,0)*sqrt(adentro); // rho(j+1) = rho(j+1)*dsqrt(1.0d0-(azeta/rho(j+1))**2)
			}
			else // go to 90
			{// 80       continue
				info = 1;
				rho(j,0) = -1.0e0; // rho(j+1) = -1.0d0
			}  // 90       continue
		} //100    continue
	}// 110  continue
	//  120 continue
	return;
	//      end
}



void L_Matrix::linpack_zchex(L_Matrix &r, int ldr, int p, int k, int l, L_Matrix &z, int ldz, int nz, L_Matrix &c, L_Matrix &s, int job)
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

         for (i=0; i<l; i++)  //do 20 i+1 = 1, l
		 {
            ii = l - i + 1;
			s(i,0) = r(ii,l-1); //s(i+1) = r(ii+1,l)
		 }//20    continue
         for (jj=k-1; jj<lm1; jj++) //do 40 jj+1 = k, lm1
		 {
            j = lm1 - jj + k;
            for (i=0; i<j+1; i++)  //do 30 i+1 = 1, j+1
				r(i,j+1) = r(i,j); //r(i+1,j+1+1) = r(i+1,j+1)
            //30       continue
			r(j+1+1,j+1+1) = 0; //r(j+1+1,j+1+1) = (0.0d0,0.0d0)
		 }//40    continue
         if (k != 1) //if (k .eq. 1) go to 60
		 {
            for (i=0; i<km1; i++) //do 50 i+1 = 1, km1
			{
               ii = l - i - 1;
			   r(i,k-1) = s(ii,0); //r(i+1,k) = s(ii+1)
		    } // 50       continue
		 } //60    continue

//c      calculate the rotations.

		 t = s(0,0); //t = s(1)
         for (i=0; i<lmk; i++) //do 70 i+1 = 1, lmk
		 {
			 blas_zrotg(s(i+1,0), t, c(i,0), s(i,0)); //call zrotg(s(i+1+1),t,c(i+1),s(i+1))
			 t = s(i+1,0); //t = s(i+1+1)
		 } // 70    continue
		 r(k-1,k-1) = t;  //r(k,k) = t
         for (j=kp1-1; j<p; j++) //do 90 j+1 = kp1, p
		 {
            il = L_MAX(1,l-j-1+1); // max0(1,l-j+1+1)
            for (ii=il-1; ii<lmk; ii++) //do 80 ii+1 = il, lmk
			{
               i = l - ii;
			   t = c(ii,0)*r(i,j) + s(ii,0)*r(i+1,j); // t = c(ii+1)*r(i+1,j+1) + s(ii+1)*r(i+1+1,j+1)
			   r(i+1,j) = c(ii,0)*r(i+1,j) - s(ii,0)*r(i,j); // r(i+1+1,j+1) = c(ii+1)*r(i+1+1,j+1) - dconjg(s(ii+1))*r(i+1,j+1)
			   r(i,j) = t; // r(i+1,j+1) = t
			} //80       continue
		 } //90    continue

//c      if required, apply the transformations to z.

         if (nz >= 1) //if (nz .lt. 1) go to 120
		 {
         for (j=0; j<nz; j++) //do 110 j+1 = 1, nz
		 {
            for (ii=0; ii<lmk; ii++) //do 100 ii+1 = 1, lmk
			{
               i = l - ii-2;
			   t = c(ii,0)*z(i,j) + s(ii,0)*z(i+1,j); // t = c(ii+1)*z(i+1,j+1) + s(ii+1)*z(i+1+1,j+1)
			   z(i+1,j) = c(ii,0)*z(i+1,j) - s(ii,0)*z(i,j); //  z(i+1+1,j+1) = c(ii+1)*z(i+1+1,j+1) - dconjg(s(ii+1))*z(i+1,j+1)
			   z(i,j) = t; // z(i+1,j+1) = t
			}// 100       continue
		 } //110    continue
		 } // 120    continue
      //go to 260
  	  } // 130 continue
	  else  // job
	  {
//c   left circular shift
//c      reorder the columns

         for (i=0; i<k; i++) // do 140 i+1 = 1, k
		 {
            ii = lmk + i - 2;
			s(ii,0) = r(i,k-1); // s(ii+1) = r(i+1,k)
  
		 } // 140    continue
         for (j=k-1; j<lm1; j++) // do 160 j+1 = k, lm1
		 {
            for (i=0; i<j+1; i++) // do 150 i+1 = 1, j+1
				r(i,j) = r(i,j+1); // r(i+1,j+1) = r(i+1,j+1+1)
			// 150       continue
            jj = j - km1 - 2;
			s(jj,0) = r(j+1,j+1);// s(jj+1) = r(j+1+1,j+1+1)
		 } // 160    continue
         for (i=0; i<k; i++) // do 170 i+1 = 1, k
		 {
            ii = lmk + i;
			r(i,l-1) = s(ii,0); // r(i+1,l) = s(ii+1)
		 } // 170    continue
         for (i=kp1-1; i<l; i++) //do 180 i+1 = kp1, l
		 {
			 r(i,l-1) = 0; //r(i+1,l) = (0.0d0,0.0d0)
		 } // 180    continue

//c      reduction loop.

         for (j=k-1; j<p; j++) // do 220 j+1 = k, p
		 {
            if (j+1 != k) // if (j+1 .eq. k) go to 200
			{

//c            apply the rotations.

               iu = L_MIN(j,l-2); // iu = min0(j,l-1)
               for (i=k-1; i<iu; i++) // do 190 i+1 = k, iu
			   {
                  ii = i - k - 1;
				  t = c(ii,0)*r(i,j) + s(ii,0)*r(i+1,j); // t = c(ii+1)*r(i+1,j+1) + s(ii+1)*r(i+1+1,j+1)
				  r(i+1,j) = c(ii,0)*r(i+1,j) - s(ii,0)*r(i,j); // r(i+1+1,j+1) = c(ii+1)*r(i+1+1,j+1) - dconjg(s(ii+1))*r(i+1,j+1)  // Error: decir *r(i+1,j+1)
				  r(i,j) = t; // r(i+1,j+1) = t
			   } // 190          continue
			} // 200       continue
            if (j+1 < l) // if (j+1 .ge. l) go to 210
			{
			   jj = j - k - 1;
			   t = s(jj,0); // t = s(jj+1)
			   blas_zrotg(r(j,j),t,c(jj,0),s(jj,0)); // call zrotg(r(j+1,j+1),t,c(jj+1),s(jj+1))
			} // 210       continue
		 } // 220    continue

  //c      apply the rotations to z.

         if (nz > 1) // if (nz .lt. 1) go to 250
		 {
         for (j=0; j<nz; j++) // do 240 j+1 = 1, nz
		 {
            for (i=k-1; i<lm1; i++) // do 230 i+1 = k, lm1
			{
               ii = i - km1 - 2;
			   t = c(ii,0)*z(i,j) + s(ii,0)*z(i+1,j); // t = c(ii+1)*z(i+1,j+1) + s(ii+1)*z(i+1+1,j+1)
			   z(i+1,j) = c(ii,0)*z(i+1,j) - s(ii,0)*z(i,j); // z(i+1+1,j+1) = c(ii+1)*z(i+1+1,j+1) - dconjg(s(ii+1))*z(i+1,j+1)
			   z(i,j) = t; // z(i+1,j+1) = t
			} // 230       continue
		 } // 240    continue
		 } // 250    continue
      } // 260 continue, job
      return;
      //end
}



bool L_Matrix::cholupdate_U(const L_Matrix &z, char op)
{
	// La matriz aa no se usa, representa columnas extra que se someten al mismo proceso que *this
	bool ret = true;
	throw_L_ArgException_if((op!='+' && op!='-') || begin() == NULL || z.begin() == NULL || z.li!=li || li!=lj, "L_Matrix::cholupdate_U() : parametros incorrectos");
	if (op == '+')
	{
		L_Matrix aa, c(z.li, 1), s(z.li, 1);
		linpack_zchud(*this, li, li, z, aa, z.li, 0, aa, aa, c, s);  // Ese 0 indicaba columnas extra
	}
	else if (op == '-')
	{
		L_Matrix aa, c(z.li, 1), s(z.li, 1);
		int info;
		linpack_zchdd(*this, li, li, z, aa, z.li, 0, aa, aa, c, s, info);
		ret = (info==0);
	}
	return ret;
}

bool L_Matrix::cholupdate_L(const L_Matrix &z, char op)
{
	// La matriz aa no se usa, representa columnas extra que se someten al mismo proceso que *this
	L_Matrix ce(z.li,1), se(z.li,1);
	L_MatrizFortranTr U, c, s, aa, x;
	int info = 0;

	U.fijarRefTraspuesta(*this); // 
	c.fijarRefTraspuesta(ce);
	s.fijarRefTraspuesta(se);
	x.fijarRefTraspuesta(z);
	aa.fijarRefNula();

	throw_L_ArgException_if((op!='+' && op!='-') || begin() == NULL || z.begin() == NULL || z.li!=li || li!=lj, "L_Matrix::cholupdate_L_Fortran() : parametros incorrectos");
	if (op == '+')
		L_MatrizFortranTr::zchud(U, li, li, x, aa, li, 0, aa, aa, c, s);  // Ese 0 indicaba columnas extra
	else if (op == '-')
		L_MatrizFortranTr::zchdd(U, li, li, x, aa, li, 0, aa, aa, c, s, info);
	return info == 0;
}


int L_Matrix::restarCovarianzaRecomponer(L_Matrix &J, L_Matrix &C, double factor)
{
	// P := P - J*C*JT
	// P := P - J*W*WT*JT
	// L*LT := L*LT - Z*ZT
	//                           [z1T]
	// L*LT := L*LT - [z1|z2|..]*[z2T]
	//                           [...]
	//
	// L*LT := L*LT - z1*z1T - z2*z2T - ...     
	L_Matrix L, W, Z;
	int ret;
	if (L.cholesky_L_De(*this) == false) // nx * nx * nx
		return false;
	if (W.cholesky_L_De(C) == false) // nc * nc * nc
		return false;
	Z.OP_mult(J,W);
	ret = L.restarCovarianzaCholesky_L(Z, factor);
	OP_mult_ABT(L,L); // LTT * LT = L * LT
	return ret;
}


int L_Matrix::restarCovarianzaRecomponer(L_Matrix &C, double factor)
{
	// P := P - C
	// L*LT := L*LT - Z*ZT
	//                           [z1T]
	// L*LT := L*LT - [z1|z2|..]*[z2T]
	//                           [...]
	//
	// L*LT := L*LT - z1*z1T - z2*z2T - ...     
	L_Matrix L, Z, D, V;
	int ret;
	if (L.cholesky_L_De(*this) == false) // nx * nx * nx
		return false;
	if (Z.cholesky_L_De(C) == false) // nx * nx * nx
		return false;
	ret = L.restarCovarianzaCholesky_L(Z, factor);
	OP_mult_ABT(L,L); // LTT * LT = L * LT
	return ret;
}

int L_Matrix::restarCovarianzaCholesky_U(L_Matrix &Z, double factor)
{
	// UT*U := UT*U - Z*ZT
	//                           [z1T]
	// UT*U := UT*U - [z1|z2|..]*[z2T]
	//                           [...]
	//
	// UT*U := UT*U - z1*z1T - z2*z2T - ...     
	L_Matrix z;
	int u, i;
	int numExitos = 0;
	bool pudo;
	for (u=0; u<Z.lj; u++)  // 1 vez por cada observacion
	{
		z.reallocate(Z.li, 1);
		for (i=0; i<z.li; i++)
			z(i,0) = Z(i,u);
		pudo = cholupdate_U(z, '-'); // nx * nx
		numExitos += pudo;
		if (factor < 1.0 && pudo == false)  // Probar por segunda vez
		{
			z *= factor;
			cholupdate_U(z, '-'); // nx * nx
		}
	}
	return numExitos;
}


int L_Matrix::restarCovarianzaCholesky_L(L_Matrix &Z, double factor)
{
	// L*LT := L*LT - Z*ZT
	//                           [z1T)
	// L*LT := L*LT - [z1|z2|..)*[z2T)
	//                           [...)
	//
	// L*LT := L*LT - z1*z1T - z2*z2T - ...     
	L_Matrix z;
	int u, i;
	int numExitos = 0;
	bool pudo;
	for (u=0; u<Z.lj; u++)  // 1 vez por cada observacion
	{
		z.reallocate(Z.li, 1);
		for (i=0; i<z.li; i++)
			z(i,0) = Z(i,u);
		pudo = cholupdate_L(z, '-'); // nx * nx
		numExitos += pudo;
		if (factor < 1.0 && pudo == false)  // Probar por segunda vez
		{
			z *= factor;
			cholupdate_L(z, '-'); // nx * nx
		}
	}
	return numExitos;
}



void L_Matrix::PCA_ordenado(L_Matrix &P, L_Matrix &d_nx1) const
{
	// Espacio 1:  x, cov(x) = D  (sin covarianzas cruzadas)
	// Espacio 2:  X = P*x,   cov(X)=P*cov(x)*PT  (P = rotacion)
	// Entonces, cov(X) = P*D*PT
	// Entonces, x = PT*X => se obtiene el espacio sin covarianzas cruzadas
	// P lleva del espacio sin covarianzas cruzadas al entregado

	// datos = concatenacion de vectores = [v1 v2 v3 ... vn).

	throw_L_ArgException_if (li > lj, "PCA_ordenado: datos insuficientes para calcular covarianzas");
	L_Matrix cen(li,1), x_cen(li, lj);
	int i, j, u;
	// Centrar los datos
	cen.setZero();
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			cen(i,0) += operator()(i,j);
	for (i=0; i<li; i++)
		cen(i,0) /= lj;
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			x_cen(i,j) = operator()(i,j) - cen(i,0);
	// Calcular covarianza
	P.reallocate(li, li);
	P.setZero();
	for (i=0; i<li; i++)
		for (j=0; j<li; j++)
			for (u=0; u<lj; u++)
				P(i,j) += operator()(i,u)*operator()(j,u);
	P /= lj-1; // Por el tema de que sea insesgado y todo eso....
	// Descomponer
	d_nx1.reallocate(li, 1);
	P.vectPr_sim(d_nx1);
	// Definicion wikipedia: matriz rotacion = ortogonal con det 1
	double deter = P.det();
	if (deter < 0) // En este caso, cambiar el signo de alguna fila; por ejemplo, la ultima de todas
	{
		// Los vectores propios se pueden amplificar
		// Si se multiplican por -1 no cambia la ortonormalidad de la matriz
		for (i=0; i<P.li; i++)
			P(i,P.lj-1) *= -1;
	}
}

void L_Matrix::PCA_proyectar(L_Matrix &P, int ndims)
{
	// Espacio 1:  x, cov(x) = D  (sin covarianzas cruzadas)
	// Espacio 2:  X = P*x,   cov(X)=P*cov(x)*PT  (P = rotacion)
	// Entonces, cov(X) = P*D*PT
	// Entonces, x = PT*X => se obtiene el espacio sin covarianzas cruzadas
	// Entonces, x[:,1:ndims) = (PT*X)[:,1:ndims) => se obtiene el espacio sin covarianzas cruzadas

	int i, j;
	L_Matrix PT, PTX;
	throw_L_ArgException_if(ndims > li, "PCA_proyectar() : dimensionalidad del espacio resultante mayor que la original");
	if (ndims < 0)
		ndims = li;
	PT.transpOf(P);
	PTX.OP_mult(PT,*this);
	reallocate(ndims, lj);
	for (i=0; i<li; i++)
		for (j=0; j<lj; j++)
			operator()(i,j) = PTX(i,j);
}


void L_Matrix::pruebaDefinesMatrices()
{
	L_Matrix A, B, C, BT;
	L_Matrix_reallocate(A,5,5);
	A.fillRandomValues(-10, 10);

	B = A;
	B.print("B");
	L_Matrix_OP_assign(B, A);
	B.print("B");

	printf("---------\n");

	B.fillRandomValues(-10, 10);

	// Se usa * para probar consistencia de operadores
#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	C = A*B;
	C.print("A*B");
#endif
	L_Matrix_OP_mult(C,A,B);
	C.print("A*B");
	C.OP_mult_private(A,B);
	C.print("A*B");
	printf("presione ENTER para continuar\n");
	getchar();

	// Se usa + para probar consistencia de operadores
#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	C = A+B;
	C.print("A+B");
#endif
	L_Matrix_OP_add(C,A,B);
	C.print("A+B");

	printf("---------\n");

	// Se usa - para probar consistencia de operadores
#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	C = A-B;
	C.print("A-B");
#endif
	L_Matrix_OP_subtract(C,A,B);
	C.print("A-B");

	printf("presione ENTER para continuar\n");
	getchar();

	C = B; C.OP_amplify(3);
	C.print("B*3");
	L_Matrix_OP_assign(C,B); L_Matrix_OP_amplify(C,3);
	C.print("B*3");

	printf("---------\n");

	C.transpOf(B);
	C.print("BT");
	L_Matrix_transpOf(C,B);
	C.print("BT");
#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	C = ~B;
	C.print("BT");
#endif

	printf("---------\n");

#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
	C = A*~B;
	C.print("A*BT");
#endif
	BT.transpOf(B);
	L_Matrix_OP_mult(C,A,BT);
	C.print("A*BT");
	C.OP_mult_ABT(A,B);
	C.print("A*BT");
	printf("presione ENTER para continuar\n");
	getchar();

	C.setZero();
	C.OP_mult_sub(0,0, A, B, 0,0,2,2, 0,0,2,2);
	C.print("A*B(sub)");
	C.setZero();
	L_Matrix_OP_mult_sub(C, 0,0, A, B, 0,0,2,2, 0,0,2,2);
	C.print("A*B(sub)");

	printf("presione ENTER para continuar\n");
	getchar();
}

void L_Matrix::pruebaVectPr_2x2()
{
	L_Matrix A(2,2), P(2,2), D(2,1), Pinv(2,2), tmp(2,2), A2(2,2);
	int i, nMax = 40;
	
	for (i=0; i<nMax; i++)
	{
		A.fillRandomValues(-100,100);
		P = A;
		if (P.vectPr_2x2(D) == false)
			continue;
		Pinv.inverseOf(P);
		tmp.OP_multByDiagonal(P,D);
		A2.OP_mult(tmp,Pinv);
		A.print("A");
		A2.print("A2");
	}


}

void L_Matrix::benchmarks()
{
	L_Matrix m(100,100);
	L_Matrix a(100,100), aT;
	int frames;
	double t1, t, fps;

	m.fillRandomValues(-10, 10);

	frames = 0;
	t1 = L_TIME();
	t = L_TIME();
	while(t - t1 < 5)
	{
		m.invertMe();
		frames++;
		t = L_TIME();
	}
	fps = frames / (t-t1);
	printf("Inversa de matriz 100x100: %f[fps]\n", fps);

	a.fillRandomValues(-10, 10);
	aT.transpOf(a);
	m.OP_mult(a,aT);

	frames = 0;
	t1 = L_TIME();
	t = L_TIME();
	while(t - t1 < 5)
	{
		m.invertMeDefPos();
		frames++;
		t = L_TIME();
	}
	fps = frames / (t-t1);
	printf("Inversa de matriz 100x100 def + usando LU: %f[fps]\n", fps);
	printf("presione ENTER para salir\n");
	getchar();
}

void L_Matrix::pruebaSvd2x2()
{
	L_Matrix M1(2,2), M2(2,2), D(2,2), M(2,2), MD;
	double ang1, ang2, sang, esang;
	double e1, e2, eang1, eang2;
	int nmax = 100;
	bool retbool;

	printf("Probando %d iteraciones de svd2x2()\n", nmax);
	D.setZero();
	int nerr = 0, nerrDet = 0;
	for (int i=0; i<nmax; i++)
	{
		ang1 = L_RANDOM(-6*M_PI, 6*M_PI);
		ang2 = L_RANDOM(-6*M_PI, 6*M_PI);
		sang = fmod(200*M_PI+ang1+ang2,2*M_PI);
		D(0,0) = L_RANDOM(0, 10000);
		D(1,1) = D(0,0) + L_RANDOM(0, 10000);

		M1(0,0) = cos(ang1);
		M1(0,1) = -sin(ang1);
		M1(1,0) = sin(ang1);
		M1(1,1) = cos(ang1);

		M2(0,0) = cos(ang2);
		M2(0,1) = -sin(ang2);
		M2(1,0) = sin(ang2);
		M2(1,1) = cos(ang2);

		MD.OP_mult(M1,D);
		M.OP_mult(MD,M2);

		retbool = M.svd2x2(eang1, eang2, e1, e2);
		esang = fmod(200*M_PI+eang1+eang2,2*M_PI);
		while (ang1 > M_PI)
			ang1 -= 2*M_PI;
		while (ang1 < -M_PI)
			ang1 += M_PI;
		while (ang2 > M_PI)
			ang2 -= 2*M_PI;
		while (ang2 < -M_PI)
			ang2 += M_PI;

		double d1 = D(0,0);
		double d2 = D(1,1);
		if (retbool == false)
		{
			printf("Error numerico detectado en iteracion %i\n", i);
			printf("orig: ang = %.3f   %.3f    esc = %.3g   %.3g    + = %.3g\n", ang1, ang2, d1, d2, sang);
			printf("esti: ang = %.3f   %.3f    esc = %.3g   %.3g    + = %.3g\n", eang1, eang2, e1, e2, esang);
			nerrDet++;
			continue;
		}
		if (fabs(D(0,0) - e1)/e1 > 1.0e-4 || fabs(D(1,1) - e2)/e2 > 1.0e-4 || fabs(sang - esang) > 1.0e-4)
		{
			printf("Waaa... error no detectado en iteracion %d\n", i);
			printf("orig: ang = %.3f   %.3f    esc = %.3g   %.3g    + = %.3g\n", ang1, ang2, d1, d2, sang);
			printf("esti: ang = %.3f   %.3f    esc = %.3g   %.3g    + = %.3g\n", eang1, eang2, e1, e2, esang);
			nerr++;
		}
	}
	printf("Yaaaa\n");
	printf("Errores detectados: %d de %d\n", nerrDet, nmax);
	printf("Porcentaje: %g\n", nerrDet*100.0/nmax);
	printf("Errores no detectados: %d de %d\n", nerr, nmax);
	printf("Porcentaje: %g\n", nerr*100.0/nmax);
	getchar();
}


bool L_Matrix::ordenaDiagonal(L_Matrix &Dvector, L_Matrix &V)
{
	int i, j;
	std::vector<int> arr;
	arr.resize(Dvector.li);
	int temp;
	L_Matrix u2, d2, v2;
	u2 = *this;
	d2 = Dvector;
	v2=  V;
	for (i=0; i<d2.li; i++)
		arr[i]=i;
	// ordenamiento por burbuja de mayor a menor
	for (i=0; i<d2.li-1; i++)
	{
		for (j=i+1; j<d2.li; j++)
		{
			if (d2(arr[j-1],0) < d2(arr[j],0))
			{
				temp=arr[j-1];
				arr[j-1]=arr[j];
				arr[j]=temp;
			}
		}
	}
	for (j=0; j<d2.li; j++)
	{
		Dvector(j,0)=d2(arr[j],0);
		for (i=0; i<u2.li; i++)
			operator()(i,j)=u2(i,arr[j]);
		for (i=0; i<v2.li; i++)
			V(i,j)=v2(i,arr[j]);
	}
	return true;
}


bool L_Matrix::svd_reparaNegativos()
{
	L_Matrix D, V, VT, UD;
	if (svd_ordenado(D, V) == false)
		return false;
	int i;
	for (i=0; i<D.li; i++)
		if (D(i,0) < 0)
			D(i,0) = 0;
	UD.OP_multByDiagonal(*this, D);
	VT.transpOf(V);
	OP_mult(UD, VT);
	return true;
}

bool L_Matrix::svd_reparaFundamental()
{
	L_Matrix D, V, VT, UD;
	if (svd_ordenado(D, V) == false)
		return false;
	D(D.li-1,0) = 0;
	UD.OP_multByDiagonal(*this, D);
	VT.transpOf(V);
	OP_mult(UD, VT);
	return true;
}

bool L_Matrix::svd_reparaEsencial()
{
	L_Matrix D, V, VT, UD;
	int i;
	if (svd_ordenado(D, V) == false)
		return false;
	for (i=0; i<D.li-1; i++)
		D(i,0) = 1;
	D(D.li-1,0) = 0;
	UD.OP_multByDiagonal(*this, D);
	VT.transpOf(V);
	OP_mult(UD, VT);
	return true;
}

bool L_Matrix::svd_reparaOrtonormal()
{
	L_Matrix D, V, VT, UD;
	int i;
	if (svd_ordenado(D, V) == false)
		return false;
	for (i=0; i<D.li; i++)
		D(i,0) = 1;
	UD.OP_multByDiagonal(*this, D);
	VT.transpOf(V);
	OP_mult(UD, VT);
	return true;
}

void L_MegaMatriz::destroy()
{
	if (m!=NULL)
	{
		L_delete2d<L_Matrix>(m); // Usa delete, asi que destroy los subobjetos... es un cacho eso
		m = NULL;
	}
	ni=0;
	nj=0;
	tamanosFijados = false;
}

void L_MegaMatriz::reallocate(int ni, int nj)
{
	if (this->ni != ni || this->nj != nj)
	{
		destroy();
		this->ni = ni;
		this->nj = nj;
		m = L_new2d<L_Matrix>(ni, nj);
	}
}

bool L_MegaMatriz::coherent(std::vector<int> *lis, std::vector<int> *ljs, int *numTotI, int *numTotJ) const
{
	int i, j;
	bool noHayCero = true;
	std::vector<int> liss;
	std::vector<int> ljss;
	if (lis == NULL)
		lis = &liss;
	if (ljs == NULL)
		ljs = &ljss;
	lis->resize(ni);
	ljs->resize(nj);

	if (numTotI != NULL)
		*numTotI = 0;
	if (numTotJ != NULL)
		*numTotJ = 0;

	for (i=0; i<ni; i++)
	{
		if (i==0)
		{
			(*lis)[i] = m[i][0].li;
		if (numTotI != NULL)
			*numTotI += (*lis)[i];
			if ((*lis)[i] == 0)
				noHayCero = false;
		}
		for (j=0; j<nj; j++)
		{
			if (j==0)
			{
				(*ljs)[j] = m[0][j].li;
				if (numTotJ != NULL)
					*numTotJ = (*ljs)[j];
				if ((*ljs)[j] == 0)
					noHayCero = false;
			}
			if (m[i][j].li != (*lis)[i] || m[i][j].lj != (*ljs)[j])
			{
				throw_L_ArgException_if(true, "L_MegaMatriz::coherent()");
				noHayCero = false;
				return false; // Por si estan deshabilitadas las excepciones
			}
		}
	}
	return noHayCero;
}

bool L_MegaMatriz::pideTamanos(std::vector<int> &lis, std::vector<int> &ljs, int *numTotI, int *numTotJ)
{
	int i, j;
	bool noHayCero = true;
	lis.resize(ni);
	ljs.resize(nj);
	if (numTotI != NULL)
		*numTotI = 0;
	if (numTotJ != NULL)
		*numTotJ = 0;

	for (i=0; i<ni; i++)
	{
		lis[i] = m[i][0].li;
		if (numTotI != NULL)
			*numTotI += lis[i];
		if (lis[i] == 0)
			noHayCero = false;
	}
	for (j=0; j<nj; j++)
	{
		ljs[j] = m[0][j].li;
		if (numTotJ != NULL)
			*numTotJ = ljs[j];
		if (ljs[j] == 0)
			noHayCero = false;
	}
	return noHayCero; // Esto solo revisa si no hay ceros en la primera fila y columna
}

void L_MegaMatriz::reconstruyeFijaTamanos(const std::vector<int> &lis, const std::vector<int> &ljs)
{
	reallocate((size_type)lis.size(), (size_type)ljs.size());
	size_type i, j;
	for (i=0; i<ni; i++)
	{
		for (j=0; j<nj; j++)
		{
			if (m[i][j].begin() == NULL)
			{
				m[i][j].li = lis[i];
				m[i][j].lj = ljs[j];
			}
			else
				m[i][j].reallocate(lis[i], ljs[j]);
		}
	}
	tamanosFijados = true;
}

void L_MegaMatriz::pideMatriz(L_Matrix &resultado, std::vector<int> *lis, std::vector<int> *ljs) const
{
	std::vector<int> liss;
	std::vector<int> ljss;
	int i, j, u, v, a, b;
	int numTotI;
	int numTotJ;

	if (lis == NULL)
		lis = &liss;
	if (ljs == NULL)
		ljs = &ljss;
	coherent(lis, ljs, &numTotI, &numTotJ); // Si no es coherent esto se va a caer estripitosamente
	resultado.reallocate(numTotI, numTotJ);

	u=0;
	for (i=0; i<ni; i++)
	{
		v=0;
		for (j=0; j<nj; j++)
		{
			if (m[i][j].begin() == NULL) // Una matriz sin constriurse esta fill de ceros
			{
				for (a=0; a<m[i][j].li; a++)
					for (b=0; b<m[i][j].lj; b++)
						resultado(u+a,v+b)=0;
			}
			else
			{
				for (a=0; a<m[i][j].li; a++)
					for (b=0; b<m[i][j].lj; b++)
						resultado(u+a,v+b)=m[i][j](a,b);
			}
			v+=m[i][j].lj;
		}
		u+=m[i][0].li;
	}
}

void L_MegaMatriz::fijaMatriz(const L_Matrix &origen, const std::vector<int> &lis, const std::vector<int> &ljs, double zero)
{
	size_type i, j, u, v, a, b;
	double d;
	reallocate((size_type)lis.size(), (size_type)ljs.size());
	u=0;
	for (i=0; i<ni; i++)
	{
		v=0;
		for (j=0; j<nj; j++)
		{
			d=0;
			m[i][j].reallocate(lis[i], ljs[j]);
			for (a=0; a<m[i][j].li; a++)
				for (b=0; b<m[i][j].lj; b++)
					d+= m[i][j](a,b)=origen(u+a,v+b);
			if (d <= zero && d >= -zero) // Esto puede variar dependiendo del zero
			{
				m[i][j].destroy();
				m[i][j].li = lis[i];
				m[i][j].lj = ljs[j];
			}
			v+=ljs[j];
		}
		u+=lis[i];
	}
}

void L_MegaMatriz::fijaMatrizDeltas(const L_Matrix &origen, int deltai, int deltaj)
{
	std::vector<int> lis;
	std::vector<int> ljs;
	int i, j;
	lis.resize(origen.li / deltai);
	if (origen.li % deltai > 0)
		lis.resize(lis.size()+1);
	ljs.resize(origen.lj / deltaj);
	if (origen.lj % deltaj > 0)
		ljs.resize(ljs.size()+1);
	for (i=0; i<origen.li / deltai; i++)
		lis[i] = deltai;
	if (origen.li % deltai > 0)
		lis[lis.size()-1] = origen.li % deltai;
	for (j=0; i<origen.lj / deltaj; j++)
		ljs[j] = deltaj;
	if (origen.lj % deltaj > 0)
		ljs[ljs.size()-1] = origen.lj % deltaj;
	fijaMatriz(origen, lis, ljs);
}

void L_MegaMatriz::OP_add(const L_MegaMatriz &a, const L_MegaMatriz &b)
{
	int i, j;
	throw_L_ArgException_if(a.ni != b.ni || a.nj != b.nj, "L_MegaMatriz::OP_add()");
	reallocate(a.ni, b.nj);
	for (i=0; i<ni; i++)
		for (j=0; j<nj; j++)
			m[i][j].OP_add_NULL(a.m[i][j], b.m[i][j]);
}


void L_MegaMatriz::OP_subtract(const L_MegaMatriz &a, const L_MegaMatriz &b)
{
	int i, j;
	throw_L_ArgException_if(a.ni != b.ni || a.nj != b.nj, "L_MegaMatriz::OP_subtract()");
	reallocate(a.ni, b.nj);
	for (i=0; i<ni; i++)
		for (j=0; j<nj; j++)
			m[i][j].OP_subtract_NULL(a.m[i][j], b.m[i][j]);
}


void L_MegaMatriz::OP_mult(const L_MegaMatriz &a, const L_MegaMatriz &b)
{
	int i, j, u;
	L_Matrix sumita, res;
	bool algo = false;
	if (a.nj!=b.ni)
		throw L_ArgException();
	reallocate(a.ni, b.nj);
	for (i=0; i<ni; i++)
	{
		for (j=0; j<nj; j++)
		{
			sumita.destroy();
			sumita.li = a.m[i][j].li;
			sumita.lj = a.m[i][j].lj;
			for (u=0; u<a.ni; u++)
			{
				res.OP_mult_NULL(a.m[i][u], b.m[u][j]);
				sumita.OP_add_NULL(res);
			}
		}
	}
}

void L_MegaMatriz::inverseOf(const L_MegaMatriz &other)
{
	// Es muy complicado invertir una megamatriz
	// mejor hacerlo en el espacio normal
	L_Matrix inv;
	std::vector<int> lis;
	std::vector<int> ljs;
	other.pideMatriz(inv, &lis, &ljs);
	inv.invertMe();
	fijaMatriz(inv, ljs, lis);
}

void L_MegaMatriz::invertMe()
{
	// Es muy complicado invertir una megamatriz
	// mejor hacerlo en el espacio normal
	L_Matrix inv;
	std::vector<int> lis;
	std::vector<int> ljs;
	pideMatriz(inv, &lis, &ljs);
	inv.invertMe();
	fijaMatriz(inv, ljs, lis);
}

void L_MegaMatriz::traspuestaDe(const L_MegaMatriz &other)
{
	int i, j;
	reallocate(other.nj, other.ni);
	for (i=0; i<ni; i++)
		for (j=0; j<nj; j++)
			if (other.m[j][i].begin() == NULL)
			{
				m[i][j].destroy();
				m[i][j].li = other.m[j][i].lj;
				m[i][j].lj = other.m[j][i].li;
			}
			else
				m[i][j].transpOf(other.m[j][i]);
}

void L_MegaMatriz::insertaFilaColumna_regala(int i0, int li, int j0, int lj, const L_MegaMatriz &other)
{
	int i, j;
	reallocate(other.ni+1, other.nj+1);
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
			other.m[i][j].moveObjectTo_NULL(m[i][j]);
		m[i][j0].destroy();
		m[i][j0].li = li;
		m[i][j0].lj = lj;
		for (j=j0+1; j<nj; j++)
			other.m[i][j-1].moveObjectTo_NULL(m[i][j]);
	}
	for (j=0; j<nj; j++)
	{
		m[i0][j].destroy();
		m[i0][j].li = li;
		m[i0][j].lj = lj;
	}
	for (i=i0+1; i<li; i++)
	{
		for (j=0; j<j0; j++)
			other.m[i-1][j].moveObjectTo_NULL(m[i][j]);
		m[i][j0].destroy();
		m[i][j0].li = li;
		m[i][j0].lj = lj;
		for (j=j0+nj; j<lj; j++)
			other.m[i-1][j-1].moveObjectTo_NULL(m[i][j]);
	}
}

void L_MegaMatriz::eliminaFilasColumnas_regala(int *indi, int numi, int *indj, int numj, const L_MegaMatriz &other)
{
	reallocate(other.ni-numi, other.nj-numj);

	int i, u=0, j, v;
	for (i=0; i<ni; i++)
	{
		while (u < numi && indi[u] == i+u)
			u++;
		v = 0;
		for (j=0; j<nj; j++)
		{
			while (v < numj && indj[v] == j+v)
				v++;
			other.m[i+u][j+v].moveObjectTo_NULL(m[i][j]);
		}
	}
}

void L_MegaMatriz::eliminaFilasColumnas_regala(int i0, int numi, int j0, int numj, const L_MegaMatriz &other)
{
	reallocate(other.ni-numi, other.nj-numj);

	int i, j;
	for (i=0; i<i0; i++)
	{
		for (j=0; j<j0; j++)
			other.m[i][j].moveObjectTo_NULL(m[i][j]);
		for (j=j0; j<other.nj-numj; j++)
			other.m[i][j+numj].moveObjectTo_NULL(m[i][j]);
	}
	for (i=i0; i<other.ni-numi; i++)
	{
		for (j=0; j<j0; j++)
			other.m[i+numi][j].moveObjectTo_NULL(m[i][j]);
		for (j=j0; j<other.nj-numj; j++)
			other.m[i+numi][j+numj].moveObjectTo_NULL(m[i][j]);
	}
}

void L_MatrizFundamental::createFrom(L_Matrix &R, L_Matrix &t)
{
	L_Matrix tX;
	double length;
	length = sqrt(t(0,0)*t(0,0) + t(1,0)*t(1,0) + t(2,0)*t(2,0));
	tX.crossProductMatrix(t(0,0)/length, t(1,0)/length, t(2,0)/length);
	OP_mult(tX, R);
}

void L_MatrizFundamental::createFrom(L_Matrix &H)
{
	L_HomogeneousMatrix HDAA;
	L_Matrix RDAA(3,3), tDAA(3,1);
	L_StaticMatrix_OP_assign(HDAA, H);
	HDAA.cambiaEjes_a_DerAbaAde();
	HDAA.getRotationComponents(RDAA);
	HDAA.getTranslationComponents(tDAA);
	createFrom(RDAA,tDAA);
}

bool L_EssentialMatrix::descompone_Rt(double tanIzq1, double tanArr1, double tanIzq2, double tanArr2, L_HomogeneousMatrix &H)
{
	L_Matrix RDAA(3,3), tDAA(3,1), R, t;
	if (descompone_Rt(tanIzq1, tanArr1, tanIzq2, tanArr2, RDAA, tDAA) == false)
		return false;
	H.setRotationComponents(RDAA);
	H.fijaTraslacion(tDAA);
	H.cambiaEjes_a_AdeIzqArr();
	return true;
}

double L_MatrizFundamental::calcError(double der1, double aba1, double der2, double aba2)
{
	double ux, uy, uz, fac, e1, e2;
	ux = operator()(0,0)*der1+operator()(0,1)*aba1+operator()(0,2);
	uy = operator()(1,0)*der1+operator()(1,1)*aba1+operator()(1,2);
	uz = operator()(2,0)*der1+operator()(2,1)*aba1+operator()(2,2);
	fac = sqrt(ux*ux + uy*uy);
	ux /= fac;
	uy /= fac;
	uz /= fac;
	e1 = der2*ux+aba2*uy+uz;

	ux = operator()(0,0)*der2+operator()(1,0)*aba2+operator()(2,0);
	uy = operator()(0,1)*der2+operator()(1,1)*aba2+operator()(2,1);
	uz = operator()(0,2)*der2+operator()(1,2)*aba2+operator()(2,2);
	fac = sqrt(ux*ux + uy*uy);
	ux /= fac;
	uy /= fac;
	uz /= fac;
	e2 = der1*ux+aba1*uy+uz;

	return sqrt(e1*e1+e2*e2);
} // err = q2' E q1

double L_MatrizFundamental::calcErrorProm_qTEq(L_Matrix &u1v1_u2v2)
{
	double error=0;
	int i;
	for (i=0; i<u1v1_u2v2.li; i++)
		error += calcError(u1v1_u2v2(i,0), u1v1_u2v2(i,1), u1v1_u2v2(i,2), u1v1_u2v2(i,3));
	return error/u1v1_u2v2.li;
}

void L_MatrizFundamental::dibLinsF(L_ShapeArray &lins, int lx, int ly)
{

}

void L_EssentialMatrix::dibLinsE(L_ShapeArray &lins, double fx, double fy, int lx, int ly)
{
	// Matriz esencial: (izq2 arr2 1)' *E * (izq1 arr1 1) = 0
	// Se van a seleccionar algunos puntos (izq1 arr1) para poder dibujar las lineas

}

double L_PolinGr11::pruebaSturm()
{
	double raices[11];
	int nRaices;
	int i;
	double error = 0;
	gr = 11;
	for (i=0; i<=11; i++)
		c[i]=(5000-(rand()%10000))*0.01;
	calculaRaicesSturm(nRaices, raices);
	for (i=0; i<nRaices; i++)
		error += L_ABS(evaluate(raices[i])) / nRaices;
	return error;
}


void L_PolinGr11_matriz::print(FILE *fp, const char *name) const
{
	char *str=NULL;
	char *str2=NULL;
	int lar=0;
	int i, j;
	bool imprNom=false;

	if (name!=NULL)
	{
		lar=(int)( strlen(name)+3 );
		str=new char[lar+1];
		str2=new char[lar+1];
		strcpy(str, name);
		strcat(str, " = ");
		for (i=0; i<lar; i++)
			str2[i]=' ';
		str2[lar]=0;
	}

	for (i=0; i<li; i++)
	{
		if (str!=NULL)
		{
			if (imprNom==false && i>=li/2)
			{
				fprintf(fp, "%s", str);
				imprNom=true;
			}
			else
				fprintf(fp, "%s", str2);
		}
		fprintf(fp,"[");
		for (j=0; j<lj; j++)
		{
			elem[i][j].print(fp);
			fprintf(fp, " , ");
		}
		fprintf(fp,"]\n");
	}
	fprintf(fp,"\n");
	if (str!=NULL)
	{
		delete[] str;
		delete[] str2;
	}
}

void L_PolinGr11_matriz::printMatlabFormat(FILE *fp, const char *name) const
{
	int i, j;
	fprintf(fp, "%s = [", name);
	for (i=0; i<li; i++)
	{
		for (j=0; j<lj; j++)
			elem[i][j].print(fp);
		fprintf(fp,";");
	}
	fprintf(fp,"]\n");
}

#ifdef L_P_POS
#error L_P_POS no deberia estar definido aca
#endif
#define L_P_POS(x,y,z) ( (x)+(y)*4+(z)*16 )

// overhead aca
bool L_Polin_xyz::inicializado = false;
std::vector<int> L_Polin_xyz::grInd; // LUT para saber la relacion entre numero de columna y tripleta de exponentes
std::vector<int> L_Polin_xyz::grIndInv;

//                  0 1   2   3   4   5   6  7  8   9  10 11 12 13 14 15 16 17 18 19
// Orden columnas: x3 y3 x2y xy2 x2z y2z x2 y2 xyz xy xz2 xz x yz2 yz y  z3 z2  z  1

void L_Polin_xyz::inicializa()
{
	if (inicializado)
		return;
	inicializado = true;

	int i = 0;
	grIndInv.resize(L_POLIN_XYZ_NCOL);
	grIndInv[i++] = L_P_POS(3,0,0); grIndInv[i++] = L_P_POS(0,3,0); grIndInv[i++] = L_P_POS(2,1,0); grIndInv[i++] = L_P_POS(1,2,0);
	grIndInv[i++] = L_P_POS(2,0,1); grIndInv[i++] = L_P_POS(0,2,1); grIndInv[i++] = L_P_POS(2,0,0); grIndInv[i++] = L_P_POS(0,2,0);
	grIndInv[i++] = L_P_POS(1,1,1); grIndInv[i++] = L_P_POS(1,1,0); grIndInv[i++] = L_P_POS(1,0,2); grIndInv[i++] = L_P_POS(1,0,1);
	grIndInv[i++] = L_P_POS(1,0,0); grIndInv[i++] = L_P_POS(0,1,2); grIndInv[i++] = L_P_POS(0,1,1); grIndInv[i++] = L_P_POS(0,1,0);
	grIndInv[i++] = L_P_POS(0,0,3); grIndInv[i++] = L_P_POS(0,0,2); grIndInv[i++] = L_P_POS(0,0,1); grIndInv[i++] = L_P_POS(0,0,0);
	// Para poder mult por z
	grIndInv[i++] = L_P_POS(3,0,1); grIndInv[i++] = L_P_POS(2,0,2); grIndInv[i++] = L_P_POS(0,2,2);
	grIndInv[i++] = L_P_POS(1,0,3); grIndInv[i++] = L_P_POS(0,1,3); grIndInv[i++] = L_P_POS(0,0,4);
	// Para poder mult por x
	grIndInv[i++] = L_P_POS(1,1,2); grIndInv[i++] = L_P_POS(1,0,3);
	// Para poder mult por y
	grIndInv[i++] = L_P_POS(0,1,3);

	grInd.resize(80);
	for (i=0; i<80; i++)
		grInd[i]=-1;
	for (i=0; i<L_POLIN_XYZ_NCOL; i++)
		grInd[grIndInv[i]]=i; // grInd[pos(a,b,c)] indica el numero de columna asociado al monomio x^a*y^b*z^c
};

#undef L_P_POS



bool L_Algoritmo5Puntos::algoritmo5_mas_1_puntos(L_Matrix &E, const L_Matrix &u1v1_u2v2)
{
	L_EssentialMatrix Es[10];
	int i, j, nSols;
	double err, errMin;
	int iMin;
	L_Matrix uv(5,4);

	for (i=0; i<5; i++) // decia < 4
		for (j=0; j<4; j++) // decia < 5
			uv(i,j) = u1v1_u2v2(i,j);

	if (algoritmo5Puntos(Es, nSols, uv) == false || nSols == 0)
		return false;

	errMin = Es[0].calcError(u1v1_u2v2(5,0),u1v1_u2v2(5,1),u1v1_u2v2(5,2),u1v1_u2v2(5,3));
	errMin = errMin*errMin;
	iMin = 0;

	for (i=1; i<nSols; i++)
	{
		err = L_ABS( Es[i].calcError(u1v1_u2v2(5,0),u1v1_u2v2(5,1),u1v1_u2v2(5,2),u1v1_u2v2(5,3)) );
		err = err*err;
		if (err < errMin)  // Era afectado por el signo de err...
		{
			errMin = err;
			iMin = i;
		}
	}
	E = Es[iMin];
	return true;
}

bool L_Algoritmo5Puntos::despejar_XYZW_5p(L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W, const L_Matrix &QagrT)
{
	L_Matrix nu;
	int i;

	nu.espacioNuloPivoteo(QagrT); // QagrT ya es de 5x9 supuestamente
	// Ortonormalizar los vectores; nu es de 9x4
	nu.gramSchmidtColumnas();
	// Copiar los vectores finales
	X.reallocate(3,3);
	Y.reallocate(3,3);
	Z.reallocate(3,3);
	W.reallocate(3,3);
	for (i=0; i<9; i++)
	{
		X(i/3,i%3) = nu(i,0);
		Y(i/3,i%3) = nu(i,1);
		Z(i/3,i%3) = nu(i,2);
		W(i/3,i%3) = nu(i,3);
	}
	return true;
}

bool L_Algoritmo5Puntos::despejar_XYZW_np(L_Matrix &X, L_Matrix &Y, L_Matrix &Z, L_Matrix &W, const L_Matrix &QagrT)
{
	L_Matrix m, Qagr; // QagrT de Nx9, Qagr de 9xN
	L_Matrix d, v;
	int i, j;
	Qagr.transpOf(QagrT);
	m.OP_mult(Qagr, QagrT);
	m.svd_ordenado(d, v);
	X.reallocate(3,3);
	Y.reallocate(3,3);
	Z.reallocate(3,3);
	W.reallocate(3,3);
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			X(i,j) = m(i+j*3,5); // Ojala que sea asi... revisar!!!
			Y(i,j) = m(i+j*3,6);
			Z(i,j) = m(i+j*3,7);
			W(i,j) = m(i+j*3,8);
		}
	}
	return true;
}


bool L_Algoritmo5Puntos::calcPolinGr10(L_PolinGr11 &p, L_PolinGr11_matriz &B, L_PolinGr11_matriz &C)
{
	L_PolinGr11 n, o;
	// Copiado textual del paper
	n = B.det();
	o = C.det();
	p = n*o.c[11] - o*n.c[11]; // Con esto se erase_preserving_order la componente 11. So o.c[11]==0 y n.c[11]==0, el resultado se va a zero
	p.gr = 10;
	return true;
}



bool L_Algoritmo5Puntos::ejemploAlgoritmo5Puntos()
{
	L_CoordsCart3D p[50], pt, pi;
	L_Pose3D poseCam;
	L_Matrix Rvc, tvc;
	L_HomogeneousMatrix H, He[10];
	L_Matrix uv;
	L_EssentialMatrix E[10];
	int nSols, i, j, k;
	double e;
	bool wena = false;

	uv.reallocate(5,4);

	// La camara inicialmente esta en el origen y toma la primera foto
	// Luego se mueve a poseCam y toma la segunda foto
	poseCam.pos = L_CoordsCart3D(L_RANDOM(-100, -50), L_RANDOM(-20, 20), L_RANDOM(-20, 20));
	poseCam.ori = L_PanTiltRoll(L_RANDOM(-0.4, 0.4), L_RANDOM(-0.4, 0.4), L_RANDOM(-M_PI, M_PI)); 

	H.fijaPose3D(poseCam); // transformUsing coords de la segunda(poseCam) a la primera(I) vista
	H.invertMeTrasponiendo(); // transformUsing coords de la primera(I) vista a la segunda(poseCam) vista

	for (i=0; i<50; i++)
	{
		p[i].x = L_RANDOM(1, 100); // Coords en la primera vista
		p[i].y = p[i].x/100 * L_RANDOM(-100, 100);
		p[i].z = p[i].x/100 * L_RANDOM(-100, 100);
	}
	for (i=0; i<5; i++)
	{
		pi = p[i];  // primera vista
		pt = H*p[i]; // 
		pi.cambiaEjes_a_DerAbaAde();  /////////////// OJO CON ESTAS 2 LINEAS ///////////////////////
		pt.cambiaEjes_a_DerAbaAde();
		uv(i,0) = pi.x/pi.z; // tanDer, pos en la primera foto
		uv(i,1) = pi.y/pi.z; // tanAba
		uv(i,2) = pt.x/pt.z; // tanDer, pos en la segunda foto
		uv(i,3) = pt.y/pt.z; // tanAba
	}
	if (algoritmo5Puntos(E, nSols, uv) == false)
		printf("Fallo\n");

	for (i=0; i<nSols; i++)
	{
		if (E[i].descompone_Rt(uv(0,0), uv(0,1), uv(0,2), uv(0,3), Rvc, tvc) == false)
		{
			printf("%d (no descompone en R,t)\n", i);
			continue;
		}
		He[i].setRotationComponents(Rvc);
		He[i].fijaTraslacion(tvc);
		He[i].cambiaEjes_a_AdeIzqArr();  ///////////  Ojo con esta linea  ////////////////
		He[i].normalizaTraslacion();
	}

	H.normalizaTraslacion(); // Para que la traslacion quede con modulo 1 y sea comparable a las He[]

	for (i=0; i<nSols; i++)
	{
		e = 0;
		for (j=0; j<4; j++)
			for (k=0; k<4; k++)
				e+=(H(j,k)-He[i](j,k))*(H(j,k)-He[i](j,k));
		if (e < 1e-5)
		{
			wena = true;
			printf("Sol %d: error %.4g  :)\n", i, e);
		}
		else
		{
			printf("Sol %d: error %.4g\n", i, e);
		}
		He[i].print("H");
	}
	return wena;
}


bool L_Algoritmo5Puntos::ejemploAlgoritmo5Puntos_xAde()
{
	L_CoordsCart3D p[50], pt, pi;
	L_Pose3D poseCam;
	L_Matrix Rvc, tvc;
	L_HomogeneousMatrix H, He[10];
	L_Matrix uv;
	L_EssentialMatrix E[10];
	int nSols, i, j, k;
	double e;
	bool wena = false;

	uv.reallocate(5,4);

	// La camara inicialmente esta en el origen y toma la primera foto
	// Luego se mueve a poseCam y toma la segunda foto
	poseCam.pos = L_CoordsCart3D(L_RANDOM(-100, -50), L_RANDOM(-20, 20), L_RANDOM(-20, 20));
	poseCam.ori = L_PanTiltRoll(L_RANDOM(-0.4, 0.4), L_RANDOM(-0.4, 0.4), L_RANDOM(-M_PI, M_PI)); 

	H.fijaPose3D(poseCam); // transformUsing coords de la segunda(poseCam) a la primera(I) vista
	H.invertMeTrasponiendo(); // transformUsing coords de la primera(I) vista a la segunda(poseCam) vista

	for (i=0; i<50; i++)
	{
		p[i].x = L_RANDOM(1, 100); // Coords en la primera vista
		p[i].y = p[i].x/100 * L_RANDOM(-100, 100);
		p[i].z = p[i].x/100 * L_RANDOM(-100, 100);
	}
	for (i=0; i<5; i++)
	{
		pi = p[i];  // primera vista
		pt = H*p[i]; // 
		uv(i,0) = pi.y/pi.x; // tanIzq en la primera foto
		uv(i,1) = pi.z/pi.x; // tanArr
		uv(i,2) = pt.y/pt.x; // tanIzq en la segunda foto
		uv(i,3) = pt.z/pt.x; // tanArr
	}
	if (algoritmo5Puntos(E, nSols, uv) == false)
		printf("Fallo\n");

	for (i=0; i<nSols; i++)
	{
		E[i].print("E");
		if (E[i].descompone_Rt(uv(0,0), uv(0,1), uv(0,2), uv(0,3), He[i]) == false)
		{
			printf("%d (no descompone en R,t)\n", i);
			continue;
		}
		He[i].normalizaTraslacion();
	}

	H.normalizaTraslacion(); // Para que la traslacion quede con modulo 1 y sea comparable a las He[]

	for (i=0; i<nSols; i++)
	{
		e = 0;
		for (j=0; j<4; j++)
			for (k=0; k<4; k++)
				e+=(H(j,k)-He[i](j,k))*(H(j,k)-He[i](j,k));
		if (e < 1e-5)
		{
			wena = true;
			printf("Sol %d: error %.4g  :)\n", i, e);
		}
		else
		{
			printf("Sol %d: error %.4g\n", i, e);
		}
		He[i].print("H");
	}
	return wena;
}

bool L_Algoritmo5Puntos::ejemploAlgoritmo6Puntos()
{
	L_CoordsCart3D p[50], pt, pi;
	L_Pose3D poseCam;
	L_Matrix Rvc, tvc;
	L_HomogeneousMatrix H, He;
	L_Matrix uv;
	L_EssentialMatrix E;
	int i, j, k;
	double e;
	bool wena = false;

	uv.reallocate(6,4);

	// La camara inicialmente esta en el origen y toma la primera foto
	// Luego se mueve a poseCam y toma la segunda foto
	poseCam.pos = L_CoordsCart3D(L_RANDOM(-100, -50), L_RANDOM(-20, 20), L_RANDOM(-20, 20));
	poseCam.ori = L_PanTiltRoll(L_RANDOM(-0.4, 0.4), L_RANDOM(-0.4, 0.4), L_RANDOM(-M_PI, M_PI)); 

	H.fijaPose3D(poseCam); // transformUsing coords de la segunda vista a la primera vista
	H.invertMeTrasponiendo(); // transformUsing coords de la primera vista a la segunda vista

	for (i=0; i<50; i++)
	{
		p[i].x = L_RANDOM(1, 100); // Coords en la primera vista
		p[i].y = p[i].x/100 * L_RANDOM(1, 100);
		p[i].z = p[i].x/100 * L_RANDOM(1, 100);
	}
	for (i=0; i<6; i++)
	{
		pi = p[i];
		pt = H*p[i]; // H, al multiplicar, transformUsing coords en la primera vista a la segunda vista
		pi.cambiaEjes_a_DerAbaAde();
		pt.cambiaEjes_a_DerAbaAde();
		uv(i,0) = pi.x/pi.z; // tanDer, pos en la primera foto
		uv(i,1) = pi.y/pi.z; // tanAba
		uv(i,2) = pt.x/pt.z; // tanDer, pos en la segunda foto
		uv(i,3) = pt.y/pt.z; // tanAba
	}
	algoritmo5_mas_1_puntos(E, uv); // uv = [uv H*uv)

	if (E.descompone_Rt(uv(0,0), uv(0,1), uv(0,2), uv(0,3), Rvc, tvc) == false)
	{
		printf("E no descompone en R,t\n"); // , i);
		return false;
	}
	He.setRotationComponents(Rvc);
	He.fijaTraslacion(tvc);
	He.cambiaEjes_a_AdeIzqArr();

	He.normalizaTraslacion();
	H.normalizaTraslacion(); // Para que la traslacion quede con modulo 1 y sea comparable a las He[)

	e = 0;
	for (j=0; j<4; j++)
		for (k=0; k<4; k++)
			e+=(H(j,k)-He(j,k))*(H(j,k)-He(j,k));
	if (e < 1e-5)
	{
		wena = true;
		printf("Sol: error %.4g  :)\n", e);

	}
	else
	{
		printf("Sol: error %.4g\n", e);
	}
	He.print("H");
	return wena;
}

bool L_Algoritmo5Puntos::ejemploAlgoritmo6Puntos_xAde()
{
	L_CoordsCart3D p[50], pt, pi;
	L_Pose3D poseCam;
	L_Matrix Rvc, tvc;
	L_HomogeneousMatrix H, He;
	L_Matrix uv;
	L_EssentialMatrix E;
	int i, j, k;
	double e;
	bool wena = false;

	uv.reallocate(6,4);

	// La camara inicialmente esta en el origen y toma la primera foto
	// Luego se mueve a poseCam y toma la segunda foto
	poseCam.pos = L_CoordsCart3D(L_RANDOM(-100, -50), L_RANDOM(-20, 20), L_RANDOM(-20, 20));
	poseCam.ori = L_PanTiltRoll(L_RANDOM(-0.4, 0.4), L_RANDOM(-0.4, 0.4), L_RANDOM(-M_PI, M_PI)); 

	H.fijaPose3D(poseCam); // transformUsing coords de la segunda vista a la primera vista
	H.invertMeTrasponiendo(); // transformUsing coords de la primera vista a la segunda vista

	for (i=0; i<50; i++)
	{
		p[i].x = L_RANDOM(1, 100); // Coords en la primera vista
		p[i].y = p[i].x/100 * L_RANDOM(1, 100);
		p[i].z = p[i].x/100 * L_RANDOM(1, 100);
	}
	for (i=0; i<6; i++)
	{
		pi = p[i];
		pt = H*p[i]; // H, al multiplicar, transformUsing coords en la primera vista a la segunda vista
		uv(i,0) = pi.y/pi.x; // tanIzq en la primera foto
		uv(i,1) = pi.z/pi.x; // tanArr
		uv(i,2) = pt.y/pt.x; // tanIzq en la segunda foto
		uv(i,3) = pt.z/pt.x; // tanArr
	}
	algoritmo5_mas_1_puntos(E, uv); // uv = [uv H*uv)

	if (E.descompone_Rt(uv(0,0), uv(0,1), uv(0,2), uv(0,3), He) == false)
	{
		printf("E no descompone en R,t\n"); //, i);
		return false;
	}

	He.normalizaTraslacion();
	H.normalizaTraslacion(); // Para que la traslacion quede con modulo 1 y sea comparable a las He[)

	e = 0;
	for (j=0; j<4; j++)
		for (k=0; k<4; k++)
			e+=(H(j,k)-He(j,k))*(H(j,k)-He(j,k));
	if (e < 1e-5)
	{
		wena = true;
		printf("Sol: error %.4g  :)\n", e);
	}
	else
	{
		printf("Sol: error %.4g\n", e);
	}

	printf("Error qTEq: %f\n", E.calcErrorProm_qTEq(uv));

	He.print("H");
	return wena;
}

// Funciones para segmentar imagenes usando colores

void L_BlobColor::calcSaleImUsaVar(int lx, int ly, bool modeloCuadrado0_circulo1)
{
	double x, y, rx, ry;
	x=pos.getMeanX1();
	y=pos.getMeanX2();
	if (modeloCuadrado0_circulo1 == 0) // Modelo de cuadrado
	{
		rx = pos.getSideOfSquareFromVariance_X() / 2;
		ry = pos.getSideOfSquareFromVariance_Y() / 2;
	}
	else // Modelo de circulo
	{
		rx = pos.getRadioFromVariance_X() / 2;
		ry = pos.getRadioFromVariance_Y() / 2;
	}
	saleIzq|=((x-rx) < L_BLOB_BORDE);
	saleDer|=((x+rx) > lx-L_BLOB_BORDE-1);
	saleArr|=((y-ry) < L_BLOB_BORDE);
	saleAba|=((y+ry) > ly-L_BLOB_BORDE-1);
}


void L_SegmentadorColoresLUTdif::pideMemoHist()
{
	hist=new int***[nIndices];
	for (int a=0; a<nIndices; a++)
	{
		hist[a] = L_new3d<int>(256/dR, 256/dG, 256/dB);
		for (int i=0; i<256/dR; i++)
			for (int j=0; j<256/dG; j++)
				for (int k=0; k<256/dB; k++)
					hist[a][i][j][k]=0;
	}
}

void L_SegmentadorColoresLUTdif::entrenaHistograma(const L_ImageRGBUchar &im, const L_ImageGrayUchar &indices)
{
	int i, j;
	int u, v, w;
	if (hist==NULL)
		pideMemoHist();
	for (j=0; j<im.ly; j++)
	{
		for (i=0; i<im.lx; i++)
		{
			if (indices.pix(i,j)!=255) // indice 255 = saltar el pixel
			{
				u=im.pix(i,j,0)/dR;
				v=im.pix(i,j,1)/dG;
				w=im.pix(i,j,2)/dB;
				hist[indices.pix(i,j)][u][v][w]++;
			}
		}
	}
	nPixEjemplos+=im.lx*im.ly;
}

void L_SegmentadorColoresLUTdif::entrenaHistograma(const L_ImageRGBUchar &im, const L_ShapeArray &polig)
{
	L_ImageGrayUchar imBN;
	double xCen=0, yCen=0;
	L_uchar class_name;

	if (polig.size()==0)
		return;
	class_name = polig[0].R;
	imBN.reallocate(im.lx, im.ly);
	imBN.setConstant(255);
	imBN.genDibujo_CenGrav(polig, &xCen, &yCen);
	imBN.difusion_pinta((int)xCen, (int)yCen, (L_uchar)class_name);
	entrenaHistograma(im, imBN);
}

void L_SegmentadorColoresLUTdif::fusionaHistograma(const int ***histExt, int nPixEjemplosExt)
{
	if (hist==NULL)
		pideMemoHist();
	int nR, nG, nB;
	nR=256/dR;
	nG=256/dG;
	nB=256/dB;
	for (int a=0; a<nIndices; a++)
		for (int i=0; i<nR; i++)
			for (int j=0; j<nG; j++)
				for (int k=0; k<nB; k++)
					hist[a][i][j][k]+=histExt[i][j][k];
	nPixEjemplos+=nPixEjemplosExt;
}

void L_SegmentadorColoresLUTdif::borraHistograma(bool borrarLUT)
{
	if (hist==NULL)
		pideMemoHist();
	else
		for (int a=0; a<nIndices; a++)
			memset(**hist[a], 0, sizeof(hist[0][0][0][0])*(256/dR)*(256/dG)*(256/dB));
	nPixEjemplos=0;
	if (borrarLUT)
		memset(LUT[0][0], 0, (256/dR)*(256/dG)*(256/dB));
}

bool L_SegmentadorColoresLUTdif::generaLUT_mayorVotacion()
{
	int mayor;
	int vot;
	int votMayor;
	int numR, numG, numB;
	numR=256/dR;
	numG=256/dG;
	numB=256/dB;
	if (hist == NULL)
		return false;
	memset(LUT[0][0],0,numR*numG*numB);
	for (int i=0; i<nIndices; i++)
		colorIndice[i].resetea();
	for (int i=0; i<numR; i++)
	{
		for (int j=0; j<numG; j++)
		{
			for (int k=0; k<numB; k++)
			{
				mayor=0;
				votMayor=0;
				for (int a=0; a<nIndices; a++)
				{
					vot=hist[a][i][j][k];
					if (vot > votMayor)
					{
						votMayor = vot;
						mayor = a;
					}
				}
				LUT[i][j][k]=(char)mayor; // No tiene sentido segmentar mas de 255 niveles
				colorIndice[mayor].addFrom((L_uchar)((i+0.5)*dR),(L_uchar)((j+0.5)*dG),(L_uchar)((k+0.5)*dB), votMayor);
			}
		}
	}
	return true;
}

bool L_SegmentadorColoresLUTdif::segmentaImagen(const L_ImageRGBUchar &original)
{
	if (usarRLE)
		return segmentaRLE(original);
	return segmentaDif(original);
}


void L_SegmentadorColoresLUTdif::grabarHist(FILE *fp)
{
	int zero=0;
	fwrite(&nIndices, sizeof(int), 1, fp);
	fwrite(&dR, sizeof(int), 1, fp);
	fwrite(&dG, sizeof(int), 1, fp);
	fwrite(&dB, sizeof(int), 1, fp);
	fwrite(&zero, sizeof(int), 1, fp);
	if (hist==NULL)
		pideMemoHist();

	for (int a=0; a<nIndices; a++)
	{
		for (int i=0; i<256/dR; i++)
		{
			for (int j=0; j<256/dG; j++)
			{
				fwrite(hist[a][i][j], sizeof(long), 256/dB, fp);
			}
		}
	}
}

void L_SegmentadorColoresLUTdif::leerHist(FILE *fp)
{
	int niXD, drXD, dgXD, dbXD, ceroXD;
	size_t nleidos = 0;
	nleidos += fread(&niXD, sizeof(int), 1, fp);
	nleidos += fread(&drXD, sizeof(int), 1, fp);
	nleidos += fread(&dgXD, sizeof(int), 1, fp);
	nleidos += fread(&dbXD, sizeof(int), 1, fp);
	nleidos += fread(&ceroXD, sizeof(int), 1, fp);

	if (nleidos != 5)
	{
		printf("L_SegmentadorColoresLUTdif::leerHist() : formato de archivo incorrecto\n");
		return;
	}

	if (ceroXD!=0)
	{
		printf("L_SegmentadorColoresLUTdif::leerHist() : formato de archivo incorrecto\n");
		return;
	}
	if (niXD!=nIndices || drXD!=dR || dgXD!=dG || dbXD!=dB)
		reallocate(256/drXD, 256/dgXD, 256/dbXD, niXD);

	if (hist == NULL)
		pideMemoHist();

	for (int a=0; a<nIndices; a++)
	{
		for (int i=0; i<256/dR; i++)
		{
			for (int j=0; j<256/dG; j++)
			{
				if (fread(hist[a][i][j], sizeof(long), 256/dB, fp) != 256/dB)
				{
					printf("L_SegmentadorColoresLUTdif::leerHist() : formato de archivo incorrecto\n");
					return;
				}
			}
		}
	}
}

void L_SegmentadorColoresLUTdif::reallocate(int numR, int numG, int numB, int nIndicesTot)
{
	dR=256/numR;
	dG=256/numG;
	dB=256/numB;
	nIndices=nIndicesTot;
	nPixEjemplos=0;
	usarRLE=false;
	areaMin=16;
	if (LUT!=NULL)
		L_delete3d<char>(LUT);
	if (hist!=NULL)
	{
		for (int i=0; i<nIndices; i++)
			L_delete3d<int>(hist[i]);
		delete[] hist;
	}
	delete[] colorIndice;
	LUT = L_new3d<char>(numR, numG, numB);
	memset(LUT[0][0],0,numR*numG*numB);
	hist=NULL;
	colorIndice=new L_ColorIndice[nIndices];
}

bool L_SegmentadorColoresLUTdif::segmentaRLE(const L_ImageRGBUchar &original)
{
	return false;
}


void L_Segmentador2ColoresAdaptivo::init(L_uchar r1, L_uchar g1, L_uchar b1, L_uchar r0, L_uchar g0, L_uchar b0)
{
	double modulo;
	facR = (r1-r0);
	facG = (g1-g0);
	facB = (b1-b0);
	modulo=sqrt(facR*facR+facG*facG+facB*facB);
	facR/=modulo;
	facG/=modulo;
	facB/=modulo;
	fac1 = -facR * (r1+r0)/2 - facG * (g1+g0)/2 - facB * (b1+b0)/2;
	// idea: facR * rMed + facG * gMed + facB * bMed + fac1 = 0
}

void L_Segmentador2ColoresAdaptivo::adaptar(const L_ImageRGBUchar &original)
{
	int uMin, uMax, u, u1, u2;
	int hMin, uHMin;
	int i, j;
	u = (int)floor( original.pix(0,0,0)*facR + original.pix(0,0,1)*facG + original.pix(0,0,2)*facB + fac1 + 0.5);
	uMin = u;
	uMax = u;
	for (i=0; i<600; i++)
		histogramaProy[i] = 0;
	for (j=0; j<original.ly; j++)
	{
		for (i=0; i<original.lx; i++)
		{
			u = (int)floor( original.pix(i,j,0)*facR + original.pix(i,j,1)*facG + original.pix(i,j,2)*facB + fac1 + 0.5);
			if (u < uMin)
				uMin = u;
			if (u > uMax)
				uMax = u;
			histogramaProy[u+300]++;
		}
	}
	u1 = -30;
	if (u1 < uMin)
		u1 = uMin;
	u2 = 30;
	if (u2 > uMax)
		u2 = uMax;
	uHMin = u1;
	hMin = histogramaProy[uHMin+300];
	for (u=uMin; u<=uMax; u++)
	{
		if (histogramaProy[uHMin+300] < hMin)
		{
			uHMin = u;
			hMin = histogramaProy[u+300];
		}
	}
	// Regla: class_name = (facR*R+facG*G+facB*B+fac1 > 0) ? 1 : 0;
	// Valores extremos indican cambio no gradual, probablemente error (como ausencia de 1 de las clases)
	if (uHMin > uMin && uHMin < uMax)
		fac1 = fac1 - uHMin;  // La idea es mantener invarianza en la clasificacion
}

#ifdef NDIR
#error NDIR definido
#endif
#define NDIR 16
bool L_CaracterizadorBlobs::calcCirculo(L_ImageGrayUchar &imSegm, L_uchar c, int x, int y, double maxErr, double &r, double &cx, double &cy)
{
	int i, xf, yf;
	double dx[NDIR], dy[NDIR];
	double xs[NDIR], ys[NDIR];
	int n=0;
	bool acept[NDIR];
	for (i=0; i<NDIR; i++)
	{
		dx[i] = cos(2*M_PI*i / NDIR);
		dy[i] = sin(2*M_PI*i / NDIR);
	}
	for (i=0; i<NDIR; i++)
	{
		if (scanLineFromImage(imSegm, c, x, y, dx[i], dy[i], 8, xf, yf) == true)
		{
			xs[n]=xf;
			ys[n]=yf;
			n++;
		}
	}
	return ajustaCirculoRansac(xs, ys, n, 16, acept, 8, maxErr, r, cx, cy);
}
#undef NDIR

//! Esta funcion transformUsing "dir\x1.algo"  en  "dir", "x1" y "algo"
void L_div_name_file2(const char *dirnameext, char *dir, char *name, char *ext)
{
	long i,j,l;
	char *nomext;

	l=(int)strlen(dirnameext);

	nomext=new char[l+1];

	sprintf(nomext, "%s", dirnameext);

	for (i=l; i>=0; i--)
	{
		if (nomext[i]=='.' || nomext[i]=='\\' || nomext[i]=='/') break;
	}

	if (i==-1)
	{
		sprintf(name,"%s",nomext);
		*ext=0;
		*dir=0;
	}
	else
	{
		if (nomext[i]=='.')
		{
			for (j=i-1; j>=0; j--)
			{
				if (nomext[j]=='\\' || nomext[j]=='/') break;
			}
			if (j==-1)
			{
				nomext[i]=0;
				*dir=0;
				sprintf(name, "%s", nomext);
				sprintf(ext, "%s", &nomext[i+1]);
			}
			else
			{
				nomext[i]=0;
				nomext[j]=0;
				sprintf(dir, "%s", nomext);
				sprintf(name, "%s", &nomext[j+1]);
				sprintf(ext, "%s", &nomext[i+1]);
			}
		}
		else // hay un '/' sin '.'
		{
			nomext[i]=0;
			sprintf(dir, "%s", nomext);
			sprintf(name, "%s", &nomext[i+1]);
			*ext=0;
		}
	}
	delete[] nomext;
	return;
}

class L_ContInfo
{
public:
	char str[2];
	int n1, n2, n, dn;
	int idx;
	int width;
	L_ContInfo() {str[0]=0; str[1]=0; dn=1;}
};
class L_TokenInfo
{
public:
	L_String str;
	L_MatParser mat;
	bool esTexto;
	int width;
	L_TokenInfo() { esTexto=true; }
	L_TokenInfo &operator =(const L_TokenInfo &other) { str=other.str; mat=other.mat; esTexto=other.esTexto; return *this; }
};


int L_RelGeom2D::calcRansac(const L_Matrix &x1y1x2y2, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c, int d)
{
	throw_L_ArgException_if(x1y1x2y2.li == 0 || x1y1x2y2.lj != 4, "L_RelGeom2D::calcRansac() : matriz de datos incorrecta");
	if (x1y1x2y2.li < gradosLibertadU()+3) // Datos insuficientes
		return 0;
	std::vector<int> perm(gradosLibertadU());
	std::vector<int> mejorPerm;
	std::vector<double> errores;
	int consenso, mejorConsenso=0;
	int in, k, iIni, iFin;
	L_Matrix xyxy(gradosLibertadU(), 4);

	if (c == -1 || c > x1y1x2y2.li)
		c = L_MIN(x1y1x2y2.li, consensoMin);
	if (d == -1 || d > x1y1x2y2.li)
		d = x1y1x2y2.li;

	if (consensoAceptAutom == -1)
		consensoAceptAutom = x1y1x2y2.li;

	for (in=0; in<nIntentosMax; in++)
	{
		L_randomSelection(perm, x1y1x2y2.li);
		for (k=0; k<(int)perm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(perm[k],0);
			xyxy(k,1) = x1y1x2y2(perm[k],1);
			xyxy(k,2) = x1y1x2y2(perm[k],2);
			xyxy(k,3) = x1y1x2y2(perm[k],3);
		}
		if (calcExacto(xyxy) == false)
			continue;
		consenso = test_c_de_d(x1y1x2y2, umbralError, c, d, 0, iFin);

		if (consenso >= c && iFin+1 < x1y1x2y2.li)
		{
			iIni = iFin+1;
			consenso += test_c_de_d(x1y1x2y2, umbralError, -1, x1y1x2y2.li-iIni, iIni, iFin);
		}

		if (consenso > mejorConsenso)
		{
			mejorConsenso = consenso;
			mejorPerm = perm;
		}

		if (consenso > consensoAceptAutom)
			break;
	}

	if (mejorConsenso > 0)
	{
		for (k=0; k<(int)mejorPerm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(mejorPerm[k],0);
			xyxy(k,1) = x1y1x2y2(mejorPerm[k],1);
			xyxy(k,2) = x1y1x2y2(mejorPerm[k],2);
			xyxy(k,3) = x1y1x2y2(mejorPerm[k],3);
		}
		if (calcExacto(xyxy) == false)
			return 0;
		seleccionados.resize(x1y1x2y2.li);
		calcVectorError(errores, x1y1x2y2);

		for (k=0; k<x1y1x2y2.li; k++)
		{
			if (errores[k] < umbralError)
				seleccionados[k] = true;
			else
				seleccionados[k] = false;
		}
	}
	return mejorConsenso;
}


int L_RelGeom2D::calcProsacSimple(const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c, int d)
{
	throw_L_ArgException_if(x1y1x2y2.li == 0 || x1y1x2y2.lj != 4, "L_RelGeom2D::calcProsacSimple() : matriz de datos incorrecta");
	if (x1y1x2y2.li < gradosLibertadU()+3) // Datos insuficientes
		return 0;
	std::vector<int> indices;
	std::vector<double> errores;
	std::vector<int> perm(gradosLibertadU()), permIni(gradosLibertadU());
	std::vector<int> mejorPerm;
	std::vector<int> stack;
	int consenso, mejorConsenso=0;
	int in, k, rango, iIni, iFin;
	L_Matrix xyxy(gradosLibertadU(), 4);
	L_getOrderedIndexesMajorMinor(indices, puntajes);

	if (c == -1 || c > x1y1x2y2.li)
		c = L_MIN(x1y1x2y2.li, consensoMin);
	if (d == -1 || d > x1y1x2y2.li)
		d = x1y1x2y2.li;

	if (consensoAceptAutom == -1)
		consensoAceptAutom = x1y1x2y2.li;

	for (in=0; in<nIntentosMax; in++)
	{
		rango = gradosLibertadU() + in;
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
		if (calcExacto(xyxy) == false)
			continue;
		consenso = test_c_de_d(x1y1x2y2, umbralError, c, d, 0, iFin);

		if (consenso >= c && iFin+1 < x1y1x2y2.li)
		{
			iIni = iFin+1;
			consenso += test_c_de_d(x1y1x2y2, umbralError, -1, x1y1x2y2.li-iIni, iIni, iFin);
		}

		if (consenso > mejorConsenso)
		{
			mejorConsenso = consenso;
			mejorPerm = perm;
		}

		if (consenso > consensoAceptAutom)
			break;
	}

	if (mejorConsenso > 0)
	{
		for (k=0; k<(int)mejorPerm.size(); k++)
		{
			xyxy(k,0) = x1y1x2y2(mejorPerm[k],0);
			xyxy(k,1) = x1y1x2y2(mejorPerm[k],1);
			xyxy(k,2) = x1y1x2y2(mejorPerm[k],2);
			xyxy(k,3) = x1y1x2y2(mejorPerm[k],3);
		}
		if (calcExacto(xyxy) == false)
			return 0;
		seleccionados.resize(x1y1x2y2.li);
		calcVectorError(errores, x1y1x2y2);

		for (k=0; k<x1y1x2y2.li; k++)
		{
			if (errores[k] < umbralError)
				seleccionados[k] = true;
			else
				seleccionados[k] = false;
		}
	}
	return mejorConsenso;
}


int L_RelGeom2D::calcRansacPreemptivo(const L_Matrix &x1y1x2y2, L_Array<bool> &seleccionados, double umbralError, int nHipotesis, int nScorings, int c, int d)
{
	throw_L_ArgException_if(x1y1x2y2.li == 0 || x1y1x2y2.lj != 4, "L_RelGeom2D::calcProsacPreemptivo() : matriz de datos incorrecta");
	if (x1y1x2y2.li < gradosLibertad()+3) // Datos insuficientes
		return 0;
	std::vector<double> errores;
	L_Array<L_RelGeom2D *> hipotesis; // CUIDADO CON ESTE object
	std::vector<int> consensos;
	int consenso;
	int ihip, isco, iIni, iFin, iEval = 0, k;
	bool salir = false;

	hipotesis.resize(nHipotesis);
	generarHipotesis(hipotesis, x1y1x2y2);
	consensos.resize(nHipotesis);
	consensos.assign(consensos.size(), 0);

	// Ir matando de a poco las hipotesis
	for (ihip=0, isco=0; isco<nScorings && hipotesis.size()>1; isco++,ihip++)
	{
		iEval = iEval + d;
		if (iEval + d >= x1y1x2y2.li)
			iEval = iEval % (x1y1x2y2.li-d);

		if (ihip >= (int)hipotesis.size())
			ihip = 0;

		iIni = iEval;
		consenso = hipotesis[ihip]->test_c_de_d(x1y1x2y2, umbralError, c, d, 0, iFin);
		consensos[ihip] += consenso;

		if (consenso < c) // No paso el test
		{
			// chao con esta hipotesis
			delete hipotesis[ihip];
			hipotesis[ihip] = hipotesis[hipotesis.size()-1];
			hipotesis[hipotesis.size()-1] = NULL; // Por condicion autoimpuesta
			consensos[ihip] = consensos[hipotesis.size()-1];
			hipotesis.resize(hipotesis.size()-1);
			consensos.resize(consensos.size()-1);
		}
		ihip++;
		isco++;
	}

	// Elegir la mejor hipotesis
	int hipMax = -1;
	int consMax = 0;
	for (ihip = 0; ihip < (int)consensos.size(); ihip++)
	{
		if (consensos[ihip] > consMax)
		{
			consMax = consensos[ihip];
			hipMax = ihip;
		}
	}
	if (consMax == 0)
		return 0;
	OP_assign(hipotesis[hipMax]);
	seleccionados.resize(x1y1x2y2.li);
	calcVectorError(errores, x1y1x2y2);

	consenso = 0;
	for (k=0; k<x1y1x2y2.li; k++)
	{
		if (errores[k] < umbralError)
		{
			seleccionados[k] = true;
			consenso++;
		}
		else
			seleccionados[k] = false;
	}

	// clean las hipotesis del arreglo
	for (ihip=0; ihip < (int)hipotesis.size(); ihip++)
	{
		delete hipotesis[ihip];
		hipotesis[ihip] = NULL; // Por restriccion autoimpuesta
	}
	hipotesis.resize(0);

	return consenso;
}


int L_RelGeom2D::calcProsacPreemptivo(const L_Matrix &x1y1x2y2, const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nHipotesis, int nScorings, int c, int d)
{
	throw_L_ArgException_if(x1y1x2y2.li == 0 || x1y1x2y2.lj != 4, "L_RelGeom2D::calcProsacPreemptivo() : matriz de datos incorrecta");
	if (x1y1x2y2.li < gradosLibertad()+3) // Datos insuficientes
		return 0;
	std::vector<double> errores;
	L_Array<L_RelGeom2D *> hipotesis; // CUIDADO CON ESTE object
	std::vector<int> consensos;
	int consenso;
	int ihip, isco, iIni, iFin, iEval = 0, k;
	bool salir = false;

	hipotesis.resize(nHipotesis);
	generarHipotesis(hipotesis, x1y1x2y2, puntajes);
	consensos.resize(nHipotesis);
	consensos.assign(consensos.size(), 0);

	// Ir matando de a poco las hipotesis
	for (ihip=0, isco=0; isco<nScorings && hipotesis.size()>1; isco++,ihip++)
	{
		iEval = iEval + d;
		if (iEval + d >= x1y1x2y2.li)
			iEval = iEval % (x1y1x2y2.li-d);

		if (ihip >= (int)hipotesis.size())
			ihip = 0;

		iIni = iEval;
		consenso = hipotesis[ihip]->test_c_de_d(x1y1x2y2, umbralError, c, d, 0, iFin);
		consensos[ihip] += consenso;

		if (consenso < c) // No paso el test
		{
			// chao con esta hipotesis
			delete hipotesis[ihip];
			hipotesis[ihip] = hipotesis[hipotesis.size()-1];
			hipotesis[hipotesis.size()-1] = NULL; // Por condicion autoimpuesta
			consensos[ihip] = consensos[hipotesis.size()-1];
			hipotesis.resize(hipotesis.size()-1);
			consensos.resize(consensos.size()-1);
		}
		ihip++;
		isco++;
	}

	// Elegir la mejor hipotesis
	int hipMax = -1;
	int consMax = 0;
	for (ihip = 0; ihip < (int)consensos.size(); ihip++)
	{
		if (consensos[ihip] > consMax)
		{
			consMax = consensos[ihip];
			hipMax = ihip;
		}
	}

	if (consMax == 0)
		return 0;

	OP_assign(hipotesis[hipMax]);
	seleccionados.resize(x1y1x2y2.li);
	calcVectorError(errores, x1y1x2y2);

	consenso = 0;
	for (k=0; k<x1y1x2y2.li; k++)
	{
		if (errores[k] < umbralError)
		{
			seleccionados[k] = true;
			consenso++;
		}
		else
			seleccionados[k] = false;
	}

	// clean las hipotesis del arreglo
	for (ihip=0; ihip < (int)hipotesis.size(); ihip++)
	{
		delete hipotesis[ihip];
		hipotesis[ihip] = NULL; // Por restriccion autoimpuesta
	}
	hipotesis.resize(0);

	return consenso;
}

int L_RelGeomMuestra2D::calcRansac(L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c, int d)
{
	int tamD;
	tamD = tamDatos();
	throw_L_ArgException_if(tamD == 0, "L_RelGeomMuestra2D::calcRansac() : matriz de datos incorrecta");
	std::vector<int> perm(gradosLibertadU());
	std::vector<int> mejorPerm;
	std::vector<double> errores;
	int consenso, mejorConsenso=0;
	int in, k, iIni, iFin;

	if (c == -1 || c > tamD)
		c = L_MIN(tamD, consensoMin);
	if (d == -1 || d > tamD)
		d = tamD;

	if (consensoAceptAutom == -1)
		consensoAceptAutom = tamD;

	for (in=0; in<nIntentosMax; in++)
	{
		L_randomSelection(perm, tamD);
		if (calcExacto(perm) == false)
			continue;
		consenso = test_c_de_d(umbralError, c, d, 0, iFin);

		if (consenso >= c && iFin+1 < tamD)
		{
			iIni = iFin+1;
			consenso += test_c_de_d(umbralError, -1, tamD-iIni, iIni, iFin);
		}

		if (consenso > mejorConsenso)
		{
			mejorConsenso = consenso;
			mejorPerm = perm;
		}

		if (consenso > consensoAceptAutom)
			break;
	}

	seleccionados.resize(tamD);

	if (mejorConsenso > 0)
	{
		if (calcExacto(mejorPerm) == false)
			return 0;
		seleccionados.resize(tamD);
		calcVectorError(errores);

		for (k=0; k<tamD; k++)
		{
			if (errores[k] < umbralError)
				seleccionados[k] = true;
			else
				seleccionados[k] = false;
		}
	}
	return mejorConsenso;
}


int L_RelGeomMuestra2D::calcProsacSimple(const std::vector<double> &puntajes, L_Array<bool> &seleccionados, double umbralError, int nIntentosMax, int consensoMin, int consensoAceptAutom, int c, int d)
{
	int tamD;
	tamD = tamDatos();
	throw_L_ArgException_if(tamD == 0, "L_RelGeomMuestra2D::calcProsacSimple() : matriz de datos incorrecta");
	if (tamD < gradosLibertadU()+2 || tamD < gradosLibertad()+3) // Datos insuficientes
		return false;
	std::vector<int> indices;
	std::vector<double> errores;
	std::vector<int> perm(gradosLibertadU()), permIni(gradosLibertadU());
	std::vector<int> mejorPerm;
	int consenso, mejorConsenso=0;
	int in, k, rango, iIni, iFin;
	L_getOrderedIndexesMajorMinor(indices, puntajes);

	if (c == -1 || c > tamD)
		c = L_MIN(tamD, consensoMin);
	if (d == -1 || d > tamD)
		d = tamD;

	if (consensoAceptAutom == -1 || consensoAceptAutom > tamD)
		consensoAceptAutom = tamD;

	for (in=0; in<nIntentosMax; in++)
	{
		rango = gradosLibertadU() + in;
		if (rango >= tamD)
			rango = tamD-1;
		L_randomSelection(permIni, rango);
		for (k=0; k<(int)perm.size(); k++)
			perm[k] = indices[permIni[k]]; // permIni contiene inicialmente las mejores observaciones
		if (calcExacto(perm) == false)
			continue;
		consenso = test_c_de_d(umbralError, c, d, 0, iFin);

		if (consenso >= c && iFin+1 < tamD)
		{
			iIni = iFin+1;
			consenso += test_c_de_d(umbralError, -1, tamD-iIni, iIni, iFin);
		}

		if (consenso > mejorConsenso)
		{
			mejorConsenso = consenso;
			mejorPerm = perm;
		}

		if (consenso > consensoAceptAutom)
			break;
	}

	if (mejorConsenso > 0)
	{
		if (calcExacto(mejorPerm) == false)
			return 0;
		seleccionados.resize(tamD);
		calcVectorError(errores);

		for (k=0; k<tamD; k++)
		{
			if (errores[k] < umbralError)
				seleccionados[k] = true;
			else
				seleccionados[k] = false;
		}
	}
	return mejorConsenso;
}

void L_RelEpipolarFundamental::probando_probando()
{
	L_RelEpipolarFundamental relF;
	L_Matrix uv, F;
	L_Pose3D pose;
	L_HomogeneousMatrix H;
	L_CamaraPinhole c;
	L_RayoCamara r;
	L_CoordsCart3D p;
	L_Matrix D, V;
	int i;

	uv.reallocate(30, 4);
	pose.fijaAzar();
	H.fijaPose3D(pose);

	c.fijaParametros(320, 240, 60*M_PI/180, 45*M_PI/180);

	for (i=0; i<uv.li; i++)
	{
		p.x = L_RANDOM(1,100);
		p.y = p.x * L_RANDOM(-0.5,0.5);
		p.z = p.x * L_RANDOM(-0.5,0.5);
		r.define(p.y/p.x, p.z/p.x);
		c.calcPixel(r, uv(i,0), uv(i,1));
		p = H*p;
		r.define(p.y/p.x, p.z/p.x);
		c.calcPixel(r, uv(i,2), uv(i,3));
	}

	printf("----- Probando con el algoritmo de los 8 puntos ---\n");

	relF.calcMinCuad(uv);

	relF.F.print("F");

	printf("Error qT*F*q: %f\n", relF.F.calcErrorProm_qTEq(uv));
	relF.F.svd_ordenado(D,V);
	printf("Valores singulares: %f %f %f\n", D(0,0), D(1,0), D(2,0));


	printf("\n----- Probando con el algoritmo de los 7 puntos ---\n");

	L_RelEpipolarFundamental::algoritmo7_mas_1_puntos(relF.F, uv);

	relF.F.print("F");

	printf("Error qT*F*q: %f\n", relF.F.calcErrorProm_qTEq(uv));
	relF.F.svd_ordenado(D,V);
	printf("Valores singulares: %f %f %f\n", D(0,0), D(1,0), D(2,0));
}


void L_RelEpipolarEsencial::probando_probando()
{
	L_RelEpipolarEsencial relE;
	L_Matrix uv, F;
	L_Pose3D pose;
	L_HomogeneousMatrix H;
	L_Matrix D, V;
	L_CoordsCart3D p;
	L_Array<bool> sel;
	int i;

	uv.reallocate(30, 4);
	pose.fijaAzar();
	H.fijaPose3D(pose);

	printf("Generando 30 puntos y sus proyecciones\n");

	for (i=0; i<uv.li; i++)
	{
		p.x = L_RANDOM(1,100);
		p.y = p.x * L_RANDOM(-0.5,0.5);
		p.z = p.x * L_RANDOM(-0.5,0.5);
		uv(i,0) = p.y/p.x; // tanIzq
		uv(i,1) = p.z/p.x; // tanArr
		p = H*p;
		uv(i,2) = p.y/p.x; // tanIzq
		uv(i,3) = p.z/p.x; // tanArr
	}
	// uv tiene convencion de ejes (der,aba) por ser en pixeles

	printf("\n--- Calculando con el algoritmo de los 5 puntos\n");
	relE.calcExacto(uv);

	relE.E.print("E");

	printf("Error qT*E*q: %f\n", relE.E.calcErrorProm_qTEq(uv));
	printf("Consenso: %d\n", 
		relE.calcRansac(uv, sel, 0.0001, 5, uv.li-4, uv.li));
	relE.E.svd_ordenado(D,V);
	printf("Valores singulares: %f %f %f\n", D(0,0), D(1,0), D(2,0));


	printf("\n--- Calculando con el algoritmo de los 8 puntos\n");
	relE.calcMinCuad(uv);

	relE.E.print("E");

	printf("Error qT*E*q: %f\n", relE.E.calcErrorProm_qTEq(uv));
	printf("Consenso: %d\n", 
		relE.calcRansac(uv, sel, 0.0001, 5, uv.li-4, uv.li));
	relE.E.svd_ordenado(D,V);
	printf("Valores singulares: %f %f %f\n", D(0,0), D(1,0), D(2,0));
}


void L_TransfPoseProy2D::calcVectorError(std::vector<double> &err, int iIni, int nRev)
{
	int corr=0, incorr=0;
	L_RayoCamara rProy;
	double xs, errx, erry, err2;
	int i;
	err.resize(pun.size());
	for (i=0; i<(int)pun.size(); i++)
	{
		pose.prmovi_v(rProy, pun[i], xs);
		errx = rProy.tanIzq - ray[i+iIni].tanIzq;
		erry = rProy.tanArr - ray[i+iIni].tanArr;
		err2 = errx*errx + erry*erry;

		err[i+iIni] = (xs>0) ? sqrt(err2) : err2*1e10*xs;
	}
	return;
}


L_CoordsCart3D L_CoordsCart3D::componerVectorRot(L_CoordsCart3D &w2)
 {
	 L_Quaternion q1, q2, q3;
	 q1.fijaVectorRotacion(*this);
	 q2.fijaVectorRotacion(w2);
	 q3.OP_mult(q1,q2);
	 return q3.pideVectorRotacion();
 }

bool L_CoordsCart3D::normaliza_ret()
{
	double d=sqrt(x*x+y*y+z*z);
	if (d!=0 && d==d)
	{
		x/=d;
		y/=d;
		z/=d;
	}
	else
	{
		x=1;
		y=0;
		z=0;
		return false;
	}
	return true;
}

void L_CoordsCart3D::cambiaEjes(const L_CoordsCart3D &orig, L_SistemaEjes3D inicial, L_SistemaEjes3D final)
{
	double izq, arr, ade;
	switch(inicial)
	{
	case L_AdeIzqArr:
		ade = orig.x; izq = orig.y; arr = orig.z;
		break;
	case L_DerAbaAde:
		izq = -orig.x; arr = -orig.y; ade = orig.z;
		break;
	case L_DerArrAtr:
		izq = -orig.x; arr = orig.y; ade = -orig.z;
		break;
	case L_IzqArrAde:
		izq = orig.x; arr = orig.y; ade = orig.z;
		break;
	}
	switch(final)
	{
	case L_AdeIzqArr:
		x = ade; y = izq; z = arr;
		break;
	case L_DerAbaAde:
		x = -izq; y = -arr; z = ade;
		break;
	case L_DerArrAtr:
		x = -izq; y = arr; z = -ade;
		break;
	case L_IzqArrAde:
		x = izq; y = arr; z = ade;
		break;
	}
}

void L_CoordsCart3D::test()
{
	L_CoordsCart3D p1, p2, p3, p4, w;
	p1.fijaAzar();
	p2.fijaAzar();
	
	p1.normaliza();
	p2 = p2 - p1*(p1.punto(p2));
	p2.normaliza();

	p3 = p1.cruz(p2);

	p1.print("p1");
	p2.print("p2");
	p3.print("p3");
	printf("p2 perp p1; p3 = p1 x p2\n");
	printf("zero = %f %f %f\n", p1.punto(p2), p1.punto(p3), p2.punto(p3));
	printf("uno = %f %f %f\n", p1.pideR(), p2.pideR(), p3.pideR());
	printf("\n");

	w.define(0,0,M_PI);
	p4 = w.rotarVectRodrigues(p1);
	(w*180/M_PI).print("w(grad)");
	p1.print("p1");
	p4.print("w.rot(p1)");
	printf("\n");

	p1.print("p1");
	p2.print("p2");
	(p1+p2).print("p1+p2");
	printf("\n");

	p1.print("p1");
	p2.print("p2");
	(p1-p2).print("p1-p2");
	printf("\n");

	p1.print("p1");
	(p1*3).print("3*p1");

	printf("Presione ENTER para continuar\n");
}


L_CoordsEsfInv3D operator+(const L_CoordsEsfInv3D &a, const L_CoordsEsfInv3D &b)
{
	L_CoordsCart3D p1,p2;
	p1.fijaRThetaPhi(1/a.invR,a.theta,a.phi);
	p2.fijaRThetaPhi(1/b.invR, b.theta, b.phi);
	p1+=p2;
	L_CoordsEsfInv3D x;
	return x.fijaXYZ(p1.x, p1.y, p1.z);
}

L_CoordsEsfInv3D operator-(const L_CoordsEsfInv3D &a, const L_CoordsEsfInv3D &b)
{
	L_CoordsCart3D p1,p2;
	p1.fijaRThetaPhi(1/a.invR,a.theta,a.phi);
	p2.fijaRThetaPhi(1/b.invR, b.theta, b.phi);
	p1-=p2;
	L_CoordsEsfInv3D x;
	return x.fijaXYZ(p1.x, p1.y, p1.z);
}

void L_CoordsEsfInv3D::test()
{
	L_CoordsEsfInv3D ci, ci2;
	L_CoordsCart3D p, p2, pM, pm;
	L_Matrix J1(3,3), J2(3,3), J3(3,3);
	int i, j;
	double delta = 1.0e-4;
	
	p.fijaAzar();
	ci.fijaCoordsCart3D(p);
	p = ci.pideCoordsCart3D();
	L_CoordsEsfInv3D_pideCoordsCart3D(p2,ci);
	p.print("p");
	p2.print("p");
	printf("Parecido: %f\n", p.cosineDistanceTo(p2));

	ci2.fijaCoordsCart3D(p);
	ci.print("c");
	ci2.print("c");
	printf("Parecido: %f\n", ci.cosineDistanceTo(ci2));

	p2 = ci2.pideCoordsCart3D();
	p.print("p");
	p2.print("p");
	printf("Parecido: %f\n", p.cosineDistanceTo(p2));

	ci.jacob_pideCoordsCart3D(J1);
	L_CoordsEsfInv3D_jacob_pideCoordsCart3D(J2, ci, 0, 0);

	for (j=0; j<3; j++)
	{
		ci.el(j) += delta;
		pM = ci.pideCoordsCart3D();
		ci.el(j) -= 2*delta;
		pm = ci.pideCoordsCart3D();
		ci.el(j) += delta;
		for (i=0; i<3; i++)
			J3(i,j) = (pM.el(i) - pm.el(i)) / (2*delta);
	}

	J1.print("J(func)");
	J2.print("J(#def)");
	J3.print("J(dffn)");
	printf("Parecido: %f %f %f\n", J1.cosineDistanceTo(J2), J1.cosineDistanceTo(J3), J2.cosineDistanceTo(J3));

	printf("Presione ENTER para continuar\n");
	getchar();
}

L_CoordsCart3D L_Rayo3D::intersectionWith(const L_Rayo3D &other, L_CoordsCart3D *errorInt) // Requiere que ambos rayos estén normalizados
{
	double d1, d2;
	L_CoordsCart3D p2_p1;
	double  p2_p1v1, p2_p1v2, v1v2;
	L_CoordsCart3D r1, r2, inter;
	p2_p1=other.origen-origen; // = p2-p1
	v1v2=direcc.punto(other.direcc); // v1.v2
	p2_p1v1=p2_p1.punto(direcc); // = (p2-p1).v1
	p2_p1v2=p2_p1.punto(other.direcc); // = (p2-p1).v2
	d1=( p2_p1v1 - p2_p1v2*v1v2 )/( 1-v1v2*v1v2);
	d2=( -p2_p1v2 + p2_p1v1*v1v2 )/( 1-v1v2*v1v2); // p1-p2 = -(p2-p1)
	r1=origen+direcc*d1;
	r2=other.origen+other.direcc*d2;
	inter=(r1+r2)*0.5;
	if (errorInt!=NULL)
		*errorInt=r2-r1;
	return inter;
}


void L_CoordsInv3D::jacob_pideCoords3D(L_Matrix &J3x6, int i0, int j0) const
{
	// calcular dx3d/dx6d, dado x3d = x6d.puni.pos + x6d.(invR,theta,phi).pideCoordsCart3D()
	// Parte cartesiana de las coordenadas 6d: dx3d / dx6d.puni.pos
	J3x6(i0+0,j0+0) = 1; J3x6(i0+0,j0+1) = 0; J3x6(i0+0,j0+2) = 0; // Decia 110
	J3x6(i0+1,j0+0) = 0; J3x6(i0+1,j0+1) = 1; J3x6(i0+1,j0+2) = 0; // Decia 000
	J3x6(i0+2,j0+0) = 0; J3x6(i0+2,j0+1) = 0; J3x6(i0+2,j0+2) = 1;
	// Parte (rho, pan, tilt)
	double rho = 1/dir.invR;
	double rho2 = rho*rho;
	double cos_theta = cos(dir.theta);
	double sin_theta = sin(dir.theta);
	double cos_phi = cos(dir.phi);
	double sin_phi = sin(dir.phi);

	// Parte esferica-inversa de las coordenadas 6d: dx3d.(invR,theta,phi) / dx6d.puni.pos
	J3x6(i0+0,j0+3) =  -cos_theta*cos_phi*rho2;   // dx/dinvR
	J3x6(i0+0,j0+4) =  -sin_theta*cos_phi*rho;    // dx/dtheta
	J3x6(i0+0,j0+5) =  -cos_theta*sin_phi*rho;    // dx/dphi
	J3x6(i0+1,j0+3) =  -sin_theta*cos_phi*rho2;   // dy/dinvR
	J3x6(i0+1,j0+4) =   cos_theta*cos_phi*rho;    // dy/dtheta
	J3x6(i0+1,j0+5) =  -sin_theta*sin_phi*rho;    // dy/dphi
	J3x6(i0+2,j0+3)  = -sin_phi*rho2;             //dz/dinvR
	J3x6(i0+2,j0+4)  =  0;                        //dz/dtheta
	J3x6(i0+2,j0+5)  =  cos_phi*rho;              //dz/dphi
}

void L_CoordsInv3D::test()
{
	L_CoordsInv3D Vi, ViM, Vim;
	L_CoordsCart3D P, PM, Pm;
	double delta = 1.0e-5;
	L_Matrix j3x6(3,6), J3x6(3,6), JJ3x6(3,6);
	int i, j;

	Vi.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));
	L_CoordsInv3D_pideCoords3D(P,Vi);
	Vi.print("Pi");
	P.print("P");

	L_CoordsInv3D_jacob_pideCoords3D(Vi,j3x6,0,0);
	Vi.jacob_pideCoords3D(J3x6);

	for (j=0; j<6; j++)
	{
		ViM = Vi; ViM.el(j) += delta;
		Vim = Vi; Vim.el(j) -= delta;
		L_CoordsInv3D_pideCoords3D(PM,ViM);
		L_CoordsInv3D_pideCoords3D(Pm,Vim);
		for (i=0; i<3; i++)
			JJ3x6(i,j) = (PM.el(i)-Pm.el(i)) / (2*delta);
	}

	j3x6.print("j3x6(#def)");
	J3x6.print("j3x6(func)");
	JJ3x6.print("j3x6(dffn)");

	printf("Parecido: %f %f %f\n", j3x6.cosineDistanceTo(J3x6), j3x6.cosineDistanceTo(JJ3x6), J3x6.cosineDistanceTo(JJ3x6));

	printf("Presione ENTER para continuar\n");
	getchar();
}

double L_PanTiltRoll::errorRespectoA(L_PanTiltRoll &other)
{
	L_Rot3DMatrix M1, M2, M3, M4;
	L_PanTiltRoll xOri;
	M1.fijarPanTiltRoll(*this);
	M2.fijarPanTiltRoll(other);
	M3.transpOf(M2); // La inversa de una rotacion es la transpuesta
	M4.OP_mult(M1, M3);
	xOri=M4.pidePanTiltRoll();
	// Dejar la orientacion expresada usando los menores angulos posibles
	if (xOri.tilt > M_PI/2)
	{
		xOri.tilt = M_PI - xOri.tilt;
		xOri.pan += M_PI;
		xOri.roll = -xOri.roll;
	}
	if (xOri.tilt < -M_PI/2)
	{
		xOri.tilt = M_PI + xOri.tilt;
		xOri.pan += M_PI;
		xOri.roll = -xOri.roll;
	}
	while (xOri.pan > M_PI)
		xOri.pan-=2*M_PI;
	while (xOri.pan < -M_PI)
		xOri.pan+=2*M_PI;
	while (xOri.roll > 2*M_PI)
		xOri.roll-=2*M_PI;
	while (xOri.roll < -2*M_PI)
		xOri.roll+=2*M_PI;
	return fabs(xOri.pan)+fabs(xOri.tilt)+fabs(xOri.roll);
}

L_CoordsCart3D L_Quaternion::rotarVector_rapido(const L_CoordsCart3D &inicial) const
{
	double ab = a*b;
	double ac = a*c;
	double ad = a*d;
	double bb = b*b;
	double bc = b*c;
	double bd = b*d;
	double cc = c*c;
	double cd = c*d;
	double dd = d*d;
	L_CoordsCart3D final;
	final.x = 2*( (-cc -dd)*inicial.x + (bc - ad)*inicial.y + (ac + bd)*inicial.z ) + inicial.x;
	final.y = 2*( (ad + bc)*inicial.x + (-bb -dd)*inicial.y + (cd - ab)*inicial.z ) + inicial.y;
	final.z = 2*( (bd - ac)*inicial.x + (ab + cd)*inicial.y + (-bb -cc)*inicial.z ) + inicial.z;
	return final;
}

L_CoordsCart3D L_Quaternion::rotarVector(const L_CoordsCart3D &inicial) const
{
	double f2 = 1.0/(a*a+b*b+c*c+d*d);
	double ab = a*b*f2;
	double ac = a*c*f2;
	double ad = a*d*f2;
	double bb = b*b*f2;
	double bc = b*c*f2;
	double bd = b*d*f2;
	double cc = c*c*f2;
	double cd = c*d*f2;
	double dd = d*d*f2;
	L_CoordsCart3D final;
	final.x = 2*( (-cc -dd)*inicial.x + (bc - ad)*inicial.y + (ac + bd)*inicial.z ) + inicial.x;
	final.y = 2*( (ad + bc)*inicial.x + (-bb -dd)*inicial.y + (cd - ab)*inicial.z ) + inicial.y;
	final.z = 2*( (bd - ac)*inicial.x + (ab + cd)*inicial.y + (-bb -cc)*inicial.z ) + inicial.z;
	return final;
}

void L_Quaternion::jacob_izq_rotarVector(const L_CoordsCart3D &inicial, const double *pre, L_StaticMatrix<3,4> &J3x4)
{
	L_Cuaternion_jacob_izq_rotarVector(inicial,pre,J3x4,0,0);
}

void L_Quaternion::jacob_der_rotarVector(L_Matrix &J3x3, int i0, int j0) const
{
	double f2 = 1.0/(a*a+b*b+c*c+d*d);
	double ab = a*b*f2;
	double ac = a*c*f2;
	double ad = a*d*f2;
	double bb = b*b*f2;
	double bc = b*c*f2;
	double bd = b*d*f2;
	double cc = c*c*f2;
	double cd = c*d*f2;
	double dd = d*d*f2;

	J3x3(0+i0,0+j0) = 2*(-cc -dd) + 1; //dfinal.x / dinicial.x
	J3x3(0+i0,1+j0) = 2*(bc - ad); //dfinal.x / dinicial.y
	J3x3(0+i0,2+j0) = 2*(ac + bd); //dfinal.x / dinicial.z
	J3x3(1+i0,0+j0) = 2*(ad + bc); //dfinal.y / dinicial.x
	J3x3(1+i0,1+j0) = 2*(-bb -dd) + 1; //dfinal.y / dinicial.y
	J3x3(1+i0,2+j0) = 2*(cd - ab); //dfinal.y / dinicial.z
	J3x3(2+i0,0+j0) = 2*(bd - ac); //dfinal.z / dinicial.x
	J3x3(2+i0,1+j0) = 2*(ab + cd); //dfinal.z / dinicial.y
	J3x3(2+i0,2+j0) = 2*(-bb -cc) + 1; //dfinal.z / dinicial.z
}

void L_Quaternion::jacob_izq_rotarVector_rapido(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0, int j0) const
{
	J3x4(i0+0,j0+0) = -2*(d*inicial.y)+2*(c*inicial.z);  // dfinal.x / da
	J3x4(i0+0,j0+1) = 2*(c*inicial.y)+2*(d*inicial.z);  // dfinal.x / db
	J3x4(i0+0,j0+2) = -4*(c*inicial.x)+2*(b*inicial.y)+2*(a*inicial.z);  // dfinal.x / dc
	J3x4(i0+0,j0+3) = -4*(d*inicial.x)-2*(a*inicial.y)+2*(b*inicial.z);  // dfinal.x / dd

	J3x4(i0+1,j0+0) = 2*(d*inicial.x)-2*(b*inicial.z);  // dfinal.y / da
	J3x4(i0+1,j0+1) = 2*(c*inicial.x)-4*(b*inicial.y)-2*(a*inicial.z);  // dfinal.y / db
	J3x4(i0+1,j0+2) = 2*(b*inicial.x)+2*(d*inicial.z);  // dfinal.y / dc
	J3x4(i0+1,j0+3) = 2*(a*inicial.x)-4*(d*inicial.y)+2*(c*inicial.z);  // dfinal.y / dd

	J3x4(i0+2,j0+0) = -2*(c*inicial.x)+2*(b*inicial.y);  // dfinal.z / da
	J3x4(i0+2,j0+1) = 2*(d*inicial.x)+2*(a*inicial.y)-4*(b*inicial.z);  // dfinal.z / db
	J3x4(i0+2,j0+2) = -2*(a*inicial.x)+2*(d*inicial.y)-4*(c*inicial.z);  // dfinal.z / dc
	J3x4(i0+2,j0+3) = 2*(b*inicial.x)+2*(c*inicial.y);  // dfinal.z / dd
}

void L_Quaternion::jacob_der_rotarVector_rapido(L_Matrix &J3x3, int i0, int j0) const
{
	double ab = a*b;
	double ac = a*c;
	double ad = a*d;
	double bb = b*b;
	double bc = b*c;
	double bd = b*d;
	double cc = c*c;
	double cd = c*d;
	double dd = d*d;

	J3x3(0+i0,0+j0) = 2*(-cc -dd) + 1; //dfinal.x / dinicial.x
	J3x3(0+i0,1+j0) = 2*(bc - ad); //dfinal.x / dinicial.y
	J3x3(0+i0,2+j0) = 2*(ac + bd); //dfinal.x / dinicial.z
	J3x3(1+i0,0+j0) = 2*(ad + bc); //dfinal.y / dinicial.x
	J3x3(1+i0,1+j0) = 2*(-bb -dd) + 1; //dfinal.y / dinicial.y
	J3x3(1+i0,2+j0) = 2*(cd - ab); //dfinal.y / dinicial.z
	J3x3(2+i0,0+j0) = 2*(bd - ac); //dfinal.z / dinicial.x
	J3x3(2+i0,1+j0) = 2*(ab + cd); //dfinal.z / dinicial.y
	J3x3(2+i0,2+j0) = 2*(-bb -cc) + 1; //dfinal.z / dinicial.z
}


L_CoordsCart3D L_Quaternion::rotarVectorInv(const L_CoordsCart3D &inicial) const
{
	double f2 = 1.0/(a*a+b*b+c*c+d*d);
	double ab = a*b*f2;
	double ac = a*c*f2;
	double ad = a*d*f2;
	double bb = b*b*f2;
	double bc = b*c*f2;
	double bd = b*d*f2;
	double cc = c*c*f2;
	double cd = c*d*f2;
	double dd = d*d*f2;
	L_CoordsCart3D final;
	final.x = 2*( (-cc -dd)*inicial.x + (bc +ad)*inicial.y + (-ac + bd)*inicial.z ) + inicial.x;
	final.y = 2*( (-ad + bc)*inicial.x + (-bb -dd)*inicial.y + (cd +ab)*inicial.z ) + inicial.y;
	final.z = 2*( (bd +ac)*inicial.x + (-ab + cd)*inicial.y + (-bb -cc)*inicial.z ) + inicial.z;
	return final;
}

L_CoordsCart3D L_Quaternion::rotarVectorInv_rapido(const L_CoordsCart3D &inicial) const
{
	double ab = a*b; // Simplemente se usa el conjugado -> se invierte el signo de (-b,-c,-d)
	double ac = a*c;
	double ad = a*d;
	double bb = b*b;
	double bc = b*c;
	double bd = b*d;
	double cc = c*c;
	double cd = c*d;
	double dd = d*d;
	L_CoordsCart3D final;
	final.x = 2*( (-cc -dd)*inicial.x + (bc +ad)*inicial.y + (-ac + bd)*inicial.z ) + inicial.x;
	final.y = 2*( (-ad + bc)*inicial.x + (-bb -dd)*inicial.y + (cd +ab)*inicial.z ) + inicial.y;
	final.z = 2*( (bd +ac)*inicial.x + (-ab + cd)*inicial.y + (-bb -cc)*inicial.z ) + inicial.z;
	return final;
}



void L_Quaternion::jacob_izq_rotarVectorInv(const L_CoordsCart3D &inicial, const double *pre, L_StaticMatrix<3,4> &J3x4)
{
	J3x4(0,0) = pre[0] * inicial.x + pre[1] * inicial.y + pre[2] * inicial.z;
	J3x4(0,1) = pre[3] * inicial.x + pre[4] * inicial.y + pre[5] * inicial.z;
	J3x4(0,2) = pre[6] * inicial.x + pre[7] * inicial.y + pre[8] * inicial.z;
	J3x4(0,3) = pre[9] * inicial.x + pre[10] * inicial.y + pre[11] * inicial.z;
	J3x4(1,0) = pre[12] * inicial.x + pre[13] * inicial.y + pre[14] * inicial.z;
	J3x4(1,1) = pre[15] * inicial.x + pre[16] * inicial.y + pre[17] * inicial.z;
	J3x4(1,2) = pre[18] * inicial.x + pre[19] * inicial.y + pre[20] * inicial.z; // Malo
	J3x4(1,3) = pre[21] * inicial.x + pre[22] * inicial.y + pre[23] * inicial.z;
	J3x4(2,0) = pre[24] * inicial.x + pre[25] * inicial.y + pre[26] * inicial.z;
	J3x4(2,1) = pre[27] * inicial.x + pre[28] * inicial.y + pre[29] * inicial.z;
	J3x4(2,2) = pre[30] * inicial.x + pre[31] * inicial.y + pre[32] * inicial.z; // Malo
	J3x4(2,3) = pre[33] * inicial.x + pre[34] * inicial.y + pre[35] * inicial.z;
}

void L_Quaternion::jacob_der_rotarVectorInv(L_Matrix &J3x3, int i0, int j0) const
{
	double f2 = 1.0/(a*a+b*b+c*c+d*d);
	double ab = a*b*f2;
	double ac = a*c*f2;
	double ad = a*d*f2;
	double bb = b*b*f2;
	double bc = b*c*f2;
	double bd = b*d*f2;
	double cc = c*c*f2;
	double cd = c*d*f2;
	double dd = d*d*f2;

	J3x3(i0+0,j0+0) = -2*cc-2*dd+1;
	J3x3(i0+0,j0+1) = 2*bc+2*ad;
	J3x3(i0+0,j0+2) = -2*ac+2*bd;

	J3x3(i0+1,j0+0) = -2*ad+2*bc;
	J3x3(i0+1,j0+1) = -2*bb-2*dd+1;
	J3x3(i0+1,j0+2) = 2*cd+2*ab;

	J3x3(i0+2,j0+0) = 2*bd+2*ac;
	J3x3(i0+2,j0+1) = -2*ab+2*cd;
	J3x3(i0+2,j0+2) = -2*bb-2*cc+1;
}


void L_Quaternion::jacob_izq_rotarVectorInv_rapido(const L_CoordsCart3D &inicial, L_Matrix &J3x4, int i0, int j0) const
{
	J3x4(i0+0,j0+0) = 2*d*inicial.y-2*c*inicial.z;
	J3x4(i0+0,j0+1) = 2*c*inicial.y+2*d*inicial.z;
	J3x4(i0+0,j0+2) = -4*c*inicial.x+2*b*inicial.y-2*a*inicial.z;
	J3x4(i0+0,j0+3) = -4*d*inicial.x+2*a*inicial.y+2*b*inicial.z;

	J3x4(i0+1,j0+0) = -2*d*inicial.x+2*b*inicial.z;
	J3x4(i0+1,j0+1) = 2*c*inicial.x-4*b*inicial.y+2*a*inicial.z;
	J3x4(i0+1,j0+2) = 2*b*inicial.x+2*d*inicial.z;
	J3x4(i0+1,j0+3) = -2*a*inicial.x-4*d*inicial.y+2*c*inicial.z;

	J3x4(i0+2,j0+0) = 2*c*inicial.x-2*b*inicial.y;
	J3x4(i0+2,j0+1) = 2*d*inicial.x-2*a*inicial.y-4*b*inicial.z;
	J3x4(i0+2,j0+2) = 2*a*inicial.x+2*d*inicial.y-4*c*inicial.z;
	J3x4(i0+2,j0+3) = 2*b*inicial.x+2*c*inicial.y;
}

void L_Quaternion::jacob_der_rotarVectorInv_rapido(L_Matrix &J3x3, int i0, int j0) const
{
	double ab = a*b;
	double ac = a*c;
	double ad = a*d;
	double bb = b*b;
	double bc = b*c;
	double bd = b*d;
	double cc = c*c;
	double cd = c*d;
	double dd = d*d;

	J3x3(i0+0,j0+0) = -2*cc-2*dd+1;
	J3x3(i0+0,j0+1) = 2*bc+2*ad;
	J3x3(i0+0,j0+2) = -2*ac+2*bd;

	J3x3(i0+1,j0+0) = -2*ad+2*bc;
	J3x3(i0+1,j0+1) = -2*bb-2*dd+1;
	J3x3(i0+1,j0+2) = 2*cd+2*ab;

	J3x3(i0+2,j0+0) = 2*bd+2*ac;
	J3x3(i0+2,j0+1) = -2*ab+2*cd;
	J3x3(i0+2,j0+2) = -2*bb-2*cc+1;
}


bool L_Quaternion::fijaDosPuntos(const L_CoordsCart3D &fijo, L_CoordsCart3D &rotado)
{
	double m, p, c2, s2, ang, c, s;
	L_CoordsCart3D w;

	m = fijo.pideR()*rotado.pideR(); // producto de los modulos
	p = fijo.punto(rotado); // producto punto
	w = fijo.cruz(rotado); // producto cruz
	if (m==0) // un vector de length zero
		return false;
	if (p==m)
		{a=1; b=0; c=0; d=0; return true;}  // Rotacion 0
	if (p==-m)
		{a=0; b=0; c=0; d=1; return true;}  // Rotacion 180° en z
	s2 = w.pideR() / m; // seno
	c2 = p / m;
	if (w.normaliza_ret() == false)
		return false;
	ang = atan2(s2,c2);
	c = cos(ang/2);
	s = sin(ang/2);
	// cuaternion
	this->a = c;
	this->b = s*w.x;
	this->c = s*w.y;
	this->d = s*w.z;
	return true;
}

double L_Quaternion::fijaCuatroPuntos(const L_CoordsCart3D &fijo1, const L_CoordsCart3D &fijo2, const L_CoordsCart3D &rotado1, const L_CoordsCart3D &rotado2)
{
	double peso;
	L_CoordsCart3D f1(fijo1), f2(fijo2), r1(rotado1), r2(rotado2), d1, d2, w, a1, a2;
	double cs1, sn1, cs2, sn2, ang1, ang2;
	if (f1.normaliza_ret() == false || f2.normaliza_ret() == false || r1.normaliza_ret() == false || r2.normaliza_ret() == false)
		return 1.0e90;
	d1 = r1 - f1;
	d2 = r2 - f2;
	w = d1.cruz(d2); // Se tiene el eje de rotacion, pero sin signo y sin angulo
	peso = w.pideR();
	if (w.normaliza_ret()==false)   // Hasta aca vamos bien
		return 0;
	if (peso == 0) // No sirve tanto
		return 0;
	// Llevar los puntos al plano de rotacion usando gram-schmidt
	f1 = f1 - w*w.punto(f1);
	f2 = f2 - w*w.punto(f2);
	r1 = r1 - w*w.punto(r1);
	r2 = r2 - w*w.punto(r2);
	if (f1.normaliza_ret() == false || f2.normaliza_ret() == false || r1.normaliza_ret() == false || r2.normaliza_ret() == false)
		return 1.0e90;
	// Ahora estan todos los puntos en el plano de rotacion
	// Es trivial calcular la rotacion usando producto punto y cruz

	cs1 = f1.punto(r1);
	sn1 = f1.cruz(r1).punto(w); // Para ver si el signo de w representa bien la rotacion
	cs2 = f2.punto(r2);
	sn2 = f2.cruz(r2).punto(w); // Para ver el signo

	ang1 = atan2(sn1, cs1); // Lo que se me ocurre aca
	ang2 = atan2(sn2, cs2); // Lo que se me ocurre aca
	ang1 = (ang1 + ang2) / 2;  // Nos quedamos con el mean de los 2 angulos
	cs1 = cos(ang1 / 2); // Para el cuaternion
	sn1 = sin(ang1 / 2); // Para el cuaternion
	a = cs1;
	b = sn1 * w.x;
	c = sn1 * w.y;
	d = sn1 * w.z;

	return peso;
}

double L_Quaternion::fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, int niter)
{
	L_CoordsCart3D w, dw, wsuma;
	std::vector<L_CoordsCart3D> m(fijo.size()), f(fijo.size()), dif(fijo.size());
	std::vector<double> difR(fijo.size());
	double ang, cs, sn, sumamods, tmp;
	int i, j;
	if (fijo.size() < 2 || fijo.size() != movido.size())
		throw L_ArgException();
	if (fijo.size() == 2)
		return fijaCuatroPuntos(fijo[0], fijo[1], movido[0], movido[1]); // Pa que mas
	for (j=0; j<(int)f.size(); j++)
	{
		m[j]=movido[j];
		if (m[j].normaliza_ret() == false)
			return 1.0e90;
		f[j]=fijo[j];
		if (f[j].normaliza_ret() == false)
			return 1.0e90;
		dif[j] = m[j] - f[j];
		difR[j] = dif[j].pideR();
		if (difR[j] != 0)
			dif[j] /= difR[j];
		else
			dif[j].fijaCero();
	}
	// Ya se tiene una estimacion inicial de w, ahora a mejorarla
	// Para eso se va a hacer un gram-schmidt iterativo sobre las restas
	// Se tiene:
	//    f[] = fijos normalizados
	//    m[] = movidos normalizados
	//    dif[] = vectores que apuntan de los fijos a los movidos (dentro del plano de rotacion), normalizados
	//    f x m es paralelo a w
	//    dif es perpendicular a w
	w.fijaCero();
	for (j=0; j<(int)f.size(); j++)
		w = w + f[j].cruz(m[j]); // Alguna estimacion inicial
	if (w.normaliza_ret() == false)
		return 1.0e90;
	// Ya se tiene una estimacion inicial de w, ahora a mejorarla
	// Para eso se va a hacer un gram-schmidt iterativo sobre las restas
	// Se tiene:
	//    f[] = fijos normalizados
	//    m[] = movidos normalizados
	//    dif[] = vectores que apuntan de los fijos a los movidos (dentro del plano de rotacion), normalizados
	//    f x m es paralelo a w
	//    dif es perpendicular a w
	for (j=0; j<(int)f.size(); j++) // Para entrar en regimen permanente, no deberia quedar w en zero...
	{
		w = w - dif[j]*(dif[j].punto(w));
		tmp = w.pideR();
		if (tmp!=0)
			w*=1/tmp;
		else
			return 1e90;
	}
	sumamods = 0;
	wsuma.fijaCero();
	for (i=0; i<niter; i++)
	{
		for (j=0; j<(int)f.size(); j++) // Ya estamos en regimen permanente
		{
			dw = dif[j]*(dif[j].punto(w));
			sumamods += dw.pideR();
			w = w - dw;
			tmp = w.pideR();
			if (tmp!=0)
				w*=1/tmp;
			else
				return 1e90;
			wsuma += w*difR[j];
		}
	}
	sumamods /= niter*f.size();
	// Ahora wsuma deberia tener una muy buena estimacion de la rotacion ya que es
	// "casi" perpendicular a todos los desplazamientos
	if (wsuma.normaliza_ret() == false)
		return 1.0e90;
	// Llevar todos los puntos al plano
	for (j=0; j<(int)f.size(); j++)
	{
		m[j] = m[j] - wsuma*(wsuma.punto(m[j]));
		f[j] = f[j] - wsuma*(wsuma.punto(f[j]));
		dif[j] = m[j] - f[j];
	}
	// Se tiene:
	//    f[] = fijos normalizados dentro del plano de rotacion
	//    m[] = movidos normalizados dentro del plano de rotacion
	//    dif[] = vectores que apuntan de los fijos a los movidos (dentro del plano de rotacion), normalizados
	//    f x m es paralelo a w
	//    dif es perpendicular a w	// Ahora calcular los angulos
	ang = 0;
	for (j=0; j<(int)f.size(); j++)
	{
		cs = f[j].punto(m[j]);
		sn = f[j].cruz(m[j]).punto(w); // Para ver si el signo de w representa bien la rotacion
		ang += atan2(sn, cs);
	}
	ang /= f.size();
	cs = cos(ang / 2); // Para el cuaternion
	sn = sin(ang / 2); // Para el cuaternion
	a = cs;
	b = sn * wsuma.x;
	c = sn * wsuma.y;
	d = sn * wsuma.z;
	return sumamods; // A medida que crece el resultado es peor
}

double L_Quaternion::fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido)
{
	L_Matrix A, AT, ATA_, b, ATb, x, U, D, V, VT;
	L_Rot3DMatrix H;
	L_CoordsCart3D f, m;
	double mini;
	int i;
	A.reallocate(3*(int)fijo.size(), 9);  // It is impossible that the cast overflow...
	A.setZero();
	b.reallocate(3*(int)fijo.size(),1);
	if (fijo.size() < 3)
		throw L_ArgException(); // Esto sirve cuando hay una cantidad decente de datos
	for (i=0; i<(int)fijo.size(); i++)
	{
		f = fijo[i];
		if (f.normaliza_ret() == false)
			return 1.0e90;
		m = movido[i];
		if (m.normaliza_ret() == false)
			return 1.0e90;
		A(3*i+0,0) = f.x;  A(3*i+0,1) = f.y;  A(3*i+0,2) = f.z;
		A(3*i+1,3) = f.x;  A(3*i+1,4) = f.y;  A(3*i+1,5) = f.z;
		A(3*i+2,6) = f.x;  A(3*i+2,7) = f.y;  A(3*i+2,8) = f.z;
		b(3*i,0) = m.x;
		b(3*i,1) = m.y;
		b(3*i,2) = m.z;
	}
	AT.transpOf(A);
	ATA_.OP_mult(AT, A);
	if (ATA_.invertMe() == false)
		return 0;
	ATb.OP_mult(AT, b);
	x.OP_mult(ATA_,ATb);
	for (i=0; i<9; i++)
		H(i/3,i%3) = x(i,0);
	mini = H.normalizaSVD();
	*this = H.pideCuaternion();
	if (normaliza_ret() == false)
		return 1.0e60;
	return mini; // A medida que se acerca de uno a zero, el resultado es peor
}

double L_Quaternion::jacob_fijo_fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int niter, int i0, int j0) // Mas ambicioso
{
	int i, j;
	L_Quaternion q0, q1;
	double delta = 1e-4, ret;
	std::vector<L_CoordsCart3D> fijo2 = fijo;
	ret = fijaNpuntos_promediaw(fijo, movido);

	for (j=0; j<3*(int)fijo.size(); j++)
	{
		fijo2[j/3].el(j%3) += delta;
		q0.fijaNpuntos_promediaw(fijo2, movido, niter);
		if (q0.OP_mult_elementwise(*this) < 0)
			q0.invSigno();
		fijo2[j/3].el(j%3) -= 2*delta;
		q1.fijaNpuntos_promediaw(fijo2, movido, niter);
		if (q1.OP_mult_elementwise(*this) < 0)
			q1.invSigno();
		fijo2[j/3].el(j%3) += delta;

		for (i=0; i<4; i++)
			J4x3n(i+i0,j+j0) = (q0.el(i)-q1.el(i)) / (2*delta);
	}
	return ret;
}

double L_Quaternion::jacob_fijo_fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int i0, int j0) // Requiere invertir matriz
{
	double mini;
	L_Quaternion q;
	L_Matrix J4x4;
	// q = inv(fijaNp(F=movido, M=fijo));
	// dq/dfijo = dinv/dfijaNp * dfijaNp/dM
	// dq/dfijo = (dinv/dq)|fijaNp * (dfijaNp/dM)|(movido,fijo)
	mini = q.jacob_movido_fijaNpuntos_matrRot(movido, fijo, J4x3n, i0, j0);
	J4x4.reallocate(4,4);
	// Invertir el cuaternion
	jacob_inverseOf(q, J4x4);
	// Calcular el jacobiano
	J4x3n.OP_mult(J4x4, J4x3n);

	return mini; // A medida que crece el resultado es peor
}

double L_Quaternion::jacob_movido_fijaNpuntos_promediaw(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int niter, int i0, int j0) // Mas ambicioso
{
	int i, j;
	L_Quaternion q0, q1;
	double delta = 1e-4, ret;
	std::vector<L_CoordsCart3D> movido2 = movido;

	ret = fijaNpuntos_promediaw(fijo, movido);
	for (j=0; j<3*(int)movido2.size(); j++)
	{
		movido2[j/3].el(j%3) += delta;
		q0.fijaNpuntos_promediaw(fijo, movido2, niter);
		if (q0.OP_mult_elementwise(*this) < 0)
			q0.invSigno();
		movido2[j/3].el(j%3) -= 2*delta;
		q1.fijaNpuntos_promediaw(fijo, movido2, niter);
		if (q1.OP_mult_elementwise(*this) < 0)
			q1.invSigno();
		movido2[j/3].el(j%3) += delta;

		for (i=0; i<4; i++)
			J4x3n(i+i0,j+j0) = (q0.el(i)-q1.el(i)) / (2*delta);
	}
	return ret;
}


double L_Quaternion::jacob_movido_fijaNpuntos_matrRot(const std::vector<L_CoordsCart3D> &fijo, const std::vector<L_CoordsCart3D> &movido, L_Matrix &J4x3n, int i0, int j0) // Requiere invertir matriz
{
	L_Matrix A, AT, ATA_, b, ATb, x, ATA_AT, J4x9;
	L_Rot3DMatrix H;
	L_CoordsCart3D f, m;
	L_Quaternion q1, q2;
	double mini;
	int i;
	A.reallocate(3*(int)fijo.size(), 9);
	A.setZero();
	b.reallocate(3*(int)fijo.size(),1);
	if (fijo.size() < 3)
		throw L_ArgException(); // Esto sirve cuando hay una cantidad decente de datos
	for (i=0; i<(int)fijo.size(); i++)
	{
		f = fijo[i];
		if (f.normaliza_ret() == false)
			return 1.0e90;
		m = movido[i];
		if (m.normaliza_ret() == false)
			return 1.0e90;
		A(3*i+0,0) = f.x;  A(3*i+0,1) = f.y;  A(3*i+0,2) = f.z;
		A(3*i+1,3) = f.x;  A(3*i+1,4) = f.y;  A(3*i+1,5) = f.z;
		A(3*i+2,6) = f.x;  A(3*i+2,7) = f.y;  A(3*i+2,8) = f.z;
		b(3*i,0) = m.x;
		b(3*i,1) = m.y;
		b(3*i,2) = m.z;
	}
	AT.transpOf(A); // A:(3*n)x9, AT:9x(3*n)
	ATA_.OP_mult(AT, A); // 9x9
	if (ATA_.invertMe() == false)
		return 0;
	ATb.OP_mult(AT, b); // 9x1
	x.OP_mult(ATA_,ATb); // 9x1
	for (i=0; i<9; i++)
		H(i/3,i%3) = x(i,0);
	ATA_AT.OP_mult(ATA_, AT); // 9x(3*n)
	mini = H.normalizaSVD_corrigeJ(ATA_AT);
	*this = H.pideCuaternion();
	if (normaliza_ret() == false)
		return 1.0e90;

	J4x9.reallocate(4,9);
	*this = H.jacob_pideCuaternion(J4x9); // Esta funcion ya fue probada en pruebaCuaternion()

	J4x3n.OP_mult_sub(i0, j0, J4x9, ATA_AT, 0, 0, J4x9.li, J4x9.lj, 0, 0, ATA_AT.li, ATA_AT.lj);
	return mini; // A medida que crece el resultado es peor
}

void L_Quaternion::OP_mult(const L_Quaternion &izq, const L_Quaternion &der)
{
	a=izq.a*der.a - izq.b*der.b - izq.c*der.c - izq.d*der.d;
	b=izq.a*der.b + izq.b*der.a + izq.c*der.d - izq.d*der.c;
	c=izq.a*der.c - izq.b*der.d + izq.c*der.a + izq.d*der.b;
	d=izq.a*der.d + izq.b*der.c - izq.c*der.b + izq.d*der.a;
}

void L_Quaternion::OP_div(const L_Quaternion &izq, const L_Quaternion &der)
{
	L_Quaternion invDer;
	invDer.inverseOf(der);
	OP_mult(izq,invDer);
}

void L_Quaternion::jacob_izq_OP_mult(const L_Quaternion &izq, const L_Quaternion &der, L_Matrix &J4x4, int i0, int j0)
{
	J4x4(0+i0,0+j0) = der.a; // da / dder.a
	J4x4(0+i0,1+j0) = -der.b; // da / dder.b
	J4x4(0+i0,2+j0) = -der.c; // da / dder.c
	J4x4(0+i0,3+j0) = -der.d; // da / dder.d
	J4x4(1+i0,0+j0) = der.b; // db / dder.a
	J4x4(1+i0,1+j0) = der.a; // db / dder.b
	J4x4(1+i0,2+j0) = der.d; // db / dder.c
	J4x4(1+i0,3+j0) = -der.c; // db / dder.d
	J4x4(2+i0,0+j0) = der.c; // dc / dder.a
	J4x4(2+i0,1+j0) = -der.d; // dc / dder.b
	J4x4(2+i0,2+j0) = der.a; // dc / dder.c
	J4x4(2+i0,3+j0) = der.b; // dc / dder.d
	J4x4(3+i0,0+j0) = der.d; // dd / dder.a
	J4x4(3+i0,1+j0) = der.c; // dd / dder.b
	J4x4(3+i0,2+j0) = -der.b; // dd / dder.c
	J4x4(3+i0,3+j0) = der.a; // dd / dder.d
}

void L_Quaternion::jacob_der_OP_mult(const L_Quaternion &izq, const L_Quaternion &der, L_Matrix &J4x4, int i0, int j0)
{
	// La multiplicacion es funcion del primer argumento y se pide el jacobiano
	J4x4(0+i0,0+j0) = izq.a; // da / dizq.a
	J4x4(0+i0,1+j0) = -izq.b; // da / dizq.b
	J4x4(0+i0,2+j0) = -izq.c; // da / dizq.c
	J4x4(0+i0,3+j0) = -izq.d; // da / dizq.d
	J4x4(1+i0,0+j0) = izq.b; // db / dizq.a
	J4x4(1+i0,1+j0) = izq.a; // db / dizq.b
	J4x4(1+i0,2+j0) = -izq.d; // db / dizq.c
	J4x4(1+i0,3+j0) = izq.c; // db / dizq.d
	J4x4(2+i0,0+j0) = izq.c; // dc / dizq.a
	J4x4(2+i0,1+j0) = izq.d; // dc / dizq.b
	J4x4(2+i0,2+j0) = izq.a; // dc / dizq.c
	J4x4(2+i0,3+j0) = -izq.b; // dc / dizq.d
	J4x4(3+i0,0+j0) = izq.d; // dd / dizq.a
	J4x4(3+i0,1+j0) = -izq.c; // dd / dizq.b
	J4x4(3+i0,2+j0) = izq.b; // dd / dizq.c
	J4x4(3+i0,3+j0) = izq.a; // dd / dizq.d
}


void L_Quaternion::jacob_normaliza(L_Matrix &J4x4)
{
	// A(a) = a/pow(a^2+b^2+c^2+d^2);
	// df/da = b^2+c^2+d^2 / (a^2+b^2+c^2+d^2)^(3/2)
	// df / db = -ab / (a^2+b^2+c^2+d^2)^(3/2)
	double suma_3_2 = std::pow(a*a+b*b+c*c+d*d, 3.0/2);
	//
	J4x4(0,0) = (b*b+c*c+d*d)/suma_3_2;
	J4x4(0,1) = -a*b / suma_3_2;
	J4x4(0,2) = -a*c / suma_3_2;
	J4x4(0,3) = -a*d / suma_3_2;
	//
	J4x4(1,0) = -b*a / suma_3_2;
	J4x4(1,1) = (a*a+c*c+d*d)/suma_3_2;
	J4x4(1,2) = -b*c / suma_3_2;
	J4x4(1,3) = -b*d / suma_3_2;
	//
	J4x4(2,0) = -c*a / suma_3_2;
	J4x4(2,1) = -c*b / suma_3_2;
	J4x4(2,2) = (a*a+b*c+d*d)/suma_3_2;
	J4x4(2,3) = -c*d / suma_3_2;
	//
	J4x4(3,0) = -d*a / suma_3_2;
	J4x4(3,1) = -d*b / suma_3_2;
	J4x4(3,2) = -d*c / suma_3_2;
	J4x4(3,3) = (b*b+c*c+d*d)/suma_3_2;
}


void L_Quaternion::jacob_inverseOf(const L_Quaternion &other, L_Matrix &J4x4, int i0, int j0)
{
	double ir=1/(other.a*other.a+other.b*other.b+other.c*other.c+other.d*other.d);
	double ir2=ir*ir;
	double q[4];
	int i, j;
	a=ir*other.a;
	b=-ir*other.b;
	c=-ir*other.c;
	d=-ir*other.d;
	q[0] = other.a; q[1]=other.b; q[2]=other.c; q[3]=other.d;
	for (i=0; i<1; i++)
		for (j=0; j<4; j++)
			J4x4(i+i0,j+j0) = (i==j)*ir - 2*q[i]*q[j]*ir2;
	for (i=1; i<4; i++)
		for (j=0; j<4; j++)
			J4x4(i+i0,j+j0) = -( (i==j)*ir - 2*q[i]*q[j]*ir2 );
}

void L_Quaternion::pruebaCuaternion()
{
	int nMax = 10, u=1, i, j, k, signos, NTOT = 30;
	L_Quaternion q1, q2, q3, q4, q5;
	L_Rot3DMatrix R1, R2, R3;
	L_CoordsCart3D p1, p2, p3, p4;
	std::vector<L_CoordsCart3D> fijo, movido;
	L_Matrix J, J2, J4x3n, J4x3n_v2, J4x3n_v3, J3x4(3,4), J3x3(3,3), J3x4_(3,4), J3x3_(3,3);
	double pre36[36];
	double delta = 1.0e-6, e;	

	printf("Esta funcion test las clases L_Quaternion, L_Rot3DMatrix y L_CoordsCart3D juntas\n");

	q1.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q1.normaliza();
	q2.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q2.normaliza();
	p1.define(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10));

	printf("Probando q1.rot(p) = R(q)*p\n");
	p2 = q1.rotarVector(p1);
	p2.print("q.rot(p)");
	R1.fijarCuaternion(q1);
	p3 = R1.OP_mult(p1);
	p3.print("R(q)*p");
	// Multiplicacion directa con cuaterniones
	q3.define(0, p1.x, p1.y, p1.z);
	q4.inverseOf(q1);
	q5 = q1*q3*q4;
	p4.define(q5.b, q5.c, q5.d);
	p4.print("q*p*q^-1");
	printf("Parecidos: %f  %f\n\n", p2.cosineDistanceTo(p3), p2.cosineDistanceTo(p4));


	q1.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q1.normaliza();
	q2.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q2.normaliza();
	p1.define(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10));

	printf("Probando q1.roti(p) = R(q^-1)*p\n");
	p2 = q1.rotarVectorInv(p1);
	p2.print("q.roti(p)");
	q2.inverseOf(q1);
	R1.fijarCuaternion(q2);
	p3 = R1.OP_mult(p1);
	p3.print("R(q^-1)*p");
	// Multiplicacion directa con cuaterniones
	q3.define(0, p1.x, p1.y, p1.z);
	q4.inverseOf(q1);
	q5 = q4*q3*q1;
	p4.define(q5.b, q5.c, q5.d);
	p4.print("q^-1*p*q");
	printf("Parecidos: %f  %f\n", p2.cosineDistanceTo(p3), p2.cosineDistanceTo(p4));


	q1.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q1.normaliza();
	q2.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
	q2.normaliza();

	printf("Probando multiplicacion de cuaterniones\n");

	q3.OP_mult(q1, q2);
	q3.print("q1*q2");
	
	R1.fijarCuaternion(q1);
	R2.fijarCuaternion(q2);
	R3.OP_mult(R1,R2);
	q4 = R3.pideCuaternion();

	q4.print("c(R(q1)*R(q2))");

	printf("Parecido: %f\n", q3.cosineDistanceTo(q4));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando inversa de cuaternion: q3 = inv(q1)\n");
	q3.inverseOf(q1);
	(q3*q1).print("q3*q1");
	(q1*q3).print("q1*q3");

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando jacobiano de inversa de cuaternion: q3 = inv(q1)\n");
	J.reallocate(4,4);
	q3.jacob_inverseOf(q1, J);
	J.print("J");
	J2.reallocate(J.li, J.lj);
	q4 = q1;
	for (j=0; j<4; j++)
	{
		q4.el(j)+=1e-7;
		q5.inverseOf(q4);
		for (i=0; i<4; i++)
			J2(i,j) = ( q5.el(i) - q3.el(i) ) / 1e-7;
		q4.el(j)-=1e-7;
	}
	J2.print("J(dif.fin)");
	printf("Parecido %f\n", J.cosineDistanceTo(J2));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();


	printf("Probando rotarVector y rotarVectorInv\n");
	q3.inverseOf(q1);
	p2 = q1.rotarVectorInv(p1);
	p2.print("q1.rotarVectorInv(p1)");
	p3 = q3.rotarVector(p1);
	p3.print("qinv.rotarVector(p1)");
	R3.fijarCuaternion(q3);
	p4 = R3.OP_mult(p1);
	p4.print("R(qinv) * p1");

	printf("Parecido %f, %f\n", p2.cosineDistanceTo(p3), p2.cosineDistanceTo(p4));
	printf("\n");

	printf("Probando jacobianos de rotarVector y rotarVectorInv\n");

	q1.jacob_izq_rotarVector(p1, J3x4);
	q1.jacob_der_rotarVector(J3x3);
	J3x4.print("J3x4");
	J3x3.print("J3x3");

	L_StaticMatrix<3,4> M3x4;
	q1.jacob_izq_rotarVector_pre(pre36);
	q1.jacob_izq_rotarVector(p1, pre36, M3x4);
	L_Matrix_OP_assign(J3x4_,M3x4);


	J3x4_.print("J3x4(pre)");

	printf("Parecidos: %f\n", J3x4_.cosineDistanceTo(J3x4));
	// p1 y q1 no deben ser alterados

	q2 = q1;
	p2 = p1;
	// jacob der rotarVector
	for (j=0; j<3; j++)
	{
		p2.el(j)+=delta;
		p3 = q2.rotarVector(p2);
		p2.el(j)-=2*delta;
		p4 = q2.rotarVector(p2);
		p2.el(j)+=delta;
		for (i=0; i<3; i++)
			J3x3_(i,j) = (p3.el(i)-p4.el(i))/(2*delta);
	}
	// jacob izq rotarVector
	for (j=0; j<4; j++)
	{
		q2.el(j)+=delta;
		p3 = q2.rotarVector(p2);
		q2.el(j)-=2*delta;
		p4 = q2.rotarVector(p2);
		q2.el(j)+=delta;
		for (i=0; i<3; i++)
			J3x4_(i,j) = (p3.el(i)-p4.el(i))/(2*delta);
	}

	J3x3_.print("J3x3(dfn)");
	J3x4_.print("J3x4(dfn)");
	printf("Parecidos: %f %f\n", J3x3_.cosineDistanceTo(J3x3), J3x4_.cosineDistanceTo(J3x4));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();


	p2 = q1.rotarVectorInv(p1);
	p2.print("q1.rotarVectorInv(p1)");
	q2.inverseOf(q1);
	p3 = q2.rotarVector(p1);
	p3.print("qinv.rotarVector(p1)");
	R3.fijarCuaternion(q2);
	p4 = R3.OP_mult(p1);
	p4.print("R(qinv) * p1");

	printf("Parecidos: %f %f\n", p2.cosineDistanceTo(p3), p2.cosineDistanceTo(p4));


	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();


	q1.jacob_izq_rotarVectorInv(p1, J3x4);
	q1.jacob_der_rotarVectorInv(J3x3);
	J3x4.print("J3x4i");
	J3x3.print("J3x3i");

	//L_Cuaternion_calc_pre(q1,pre36);
	//L_Cuaternion_jacob_izq_rotarVectorInv(p1,pre36,J3x4_,0,0);
	q1.jacob_izq_rotarVectorInv_pre(pre36);
	q1.jacob_izq_rotarVectorInv(p1, pre36, M3x4);
	L_Matrix_OP_assign(J3x4_,M3x4);
	J3x4_.print("J3x4(pre)");

	printf("Parecidos: %f\n", J3x4_.cosineDistanceTo(J3x4));

	q2 = q1;
	p2 = p1;

	for (j=0; j<4; j++)
	{
		q2.el(j)+=delta;
		p3 = q2.rotarVectorInv(p2);
		q2.el(j)-=2*delta;
		p4 = q2.rotarVectorInv(p2);
		q2.el(j)+=delta;
		for (i=0; i<3; i++)
			J3x4_(i,j) = (p3.el(i)-p4.el(i))/(2*delta);
	}

	for (j=0; j<3; j++)
	{
		p2.el(j)+=delta;
		p3 = q2.rotarVectorInv(p2);
		p2.el(j)-=2*delta;
		p4 = q2.rotarVectorInv(p2);
		p2.el(j)+=delta;
		for (i=0; i<3; i++)
			J3x3_(i,j) = (p3.el(i)-p4.el(i))/(2*delta);
	}

	J3x4_.print("J3x4i(dfn)");
	J3x3_.print("J3x3i(dfn)");

	printf("Parecidos: %f %f\n", J3x4_.cosineDistanceTo(J3x4), J3x3_.cosineDistanceTo(J3x3));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando jacobiano de pideCuaternion()\n");
	J.reallocate(4,9);
	q1.normaliza();
	R1.fijarCuaternion(q1);
	q1 = R1.pideCuaternion(); // Para que quede con el mismo signo que q3
	q3 = R1.jacob_pideCuaternion(J, 0, 0);
	q1.print("q1");
	q3.print("q3=q1");
	J.print("J");
	J2 = J;
	J2.setZero();
	printf("Ahora por diferencias finitas:\n");
	for (j=0; j<9; j++)
	{
		R1(j/3,j%3) += 1.0e-8;
		q3 = R1.pideCuaternion();
		for (i=0; i<4; i++)
			J2(i,j) = (q3.el(i)-q1.el(i))/1.0e-8;
		R1(j/3,j%3) -= 1.0e-8;
	}
	J2.print("J2");

	printf("Parecido: %f\n", J.cosineDistanceTo(J2));

	printf("Probando para varios valores al azar\n");

	e=0;
	signos = 0;
	for (k=0; k<NTOT; k++)
	{
		double delta = 1.0e-6;
		q1.define(L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1), L_RANDOM(-1,1));
		q1.normaliza();
		R1.fijarCuaternion(q1);
		q3 = R1.jacob_pideCuaternion(J, 0, 0);
		if (q1.OP_mult_elementwise(q3) < 0)
		{
			q3.invSigno();
			J*=-1;
			signos++;
		}
		J2 = J;
		J2.setZero();
		q1 = q3;
		for (j=0; j<9; j++)
		{
			R1(j/3,j%3) += delta;
			q3 = R1.pideCuaternion();
			if (q1.OP_mult_elementwise(q3) < 0)
				q3.invSigno();
			for (i=0; i<4; i++)
				J2(i,j) = (q3.el(i)-q1.el(i))/delta;
			R1(j/3,j%3) -= delta;
			for (i=0; i<4; i++)
				e += fabs(J2(i,j) - J(i,j));
		}
	}
	printf("Error acumulado: %g\n", e);
	printf("Cambios de signo en cuat(R(q)): %d/%d\n", signos, NTOT);

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando encontrar rotacion a partir de varios puntos rotados\n");

	fijo.resize(NTOT);
	movido.resize(NTOT);
	for (i=0; i<NTOT; i++)
		fijo[i].define(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10));
	for (i=0; i<NTOT; i++)
		movido[i] = q1.rotarVector(fijo[i]);
	printf("original\n");
	q1.print("q1");
	//
	q2.fijaRotacionCero();
	q2.fijaNpuntos_promediaw(fijo, movido);
	if (q2.OP_mult_elementwise(q1) < 0)	{q2.invSigno();}
	printf("fijaNPuntos_promediaw()\n");
	q2.print("q2");
	//
	q3.fijaRotacionCero();
	q3.fijaNpuntos_matrRot(fijo, movido);
	if (q3.OP_mult_elementwise(q1) < 0)	{q3.invSigno();}
	printf("fijaNpuntos_matrRot()\n");
	q3.print("q3");

	printf("Parecidos: %f %f\n", q1.cosineDistanceTo(q2), q1.cosineDistanceTo(q3));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando encontrar rotacion a partir de varios puntos rotados y jacobianos\n");

	q1.print("q1");
	q3.fijaRotacionCero();
	J4x3n.reallocate(4, NTOT*3);
	J4x3n.setZero();
	J4x3n_v2.reallocate(4, NTOT*3);
	J4x3n.setZero();
	J4x3n_v3.reallocate(4, NTOT*3);
	J4x3n.setZero();
	q2.jacob_movido_fijaNpuntos_promediaw(fijo, movido, J4x3n, 20); // J4x3n
	if (q2.OP_mult_elementwise(q1) < 0)	{q2.invSigno();	J4x3n*=-1;}
	q2.print("jacob_movido_fijaNpuntos_promediaw");
	J4x3n.print("J4x3n");
	q3.fijaRotacionCero();
	q3.jacob_movido_fijaNpuntos_matrRot(fijo, movido, J4x3n_v2); // J4x3n_v2
	if (q3.OP_mult_elementwise(q1) < 0)	{q3.invSigno();	J4x3n_v2*=-1;}
	q3.print("jacob_movido_fijaNpuntos_matrRot");
	J4x3n_v2.print("J4x3n");
	for (j=0; j<3*NTOT; j++)
	{
		movido[j/3].el(j%3) += delta;
		q4.fijaNpuntos_matrRot(fijo, movido);
		if (q4.OP_mult_elementwise(q1) < 0)	{q4.invSigno();}
		movido[j/3].el(j%3) -= 2*delta;
		q5.fijaNpuntos_matrRot(fijo, movido);
		if (q5.OP_mult_elementwise(q1) < 0)	{q5.invSigno();}
		movido[j/3].el(j%3) += delta;
		for (i=0; i<4; i++)
			J4x3n_v3(i,j) = (q4.el(i)-q5.el(i))/(2*delta);
	}
	printf("jacob_movido_fijaNpuntos_matrRot(dif.fin)\n");
	J4x3n_v3.print("J4x3n");
	
	printf("Parecido q: %f %f %f\n", q1.cosineDistanceTo(q2), q1.cosineDistanceTo(q3), q1.cosineDistanceTo(q4));

	printf("Parecido 1(jpromw) 2(matr) = %g\n", J4x3n.normL2between(J4x3n_v2) / sqrt(J4x3n.normL2() * J4x3n_v2.normL2()));
	printf("Parecido 1(jpromw) 3(matr.itr) = %g\n", J4x3n.normL2between(J4x3n_v3) / sqrt(J4x3n.normL2() * J4x3n_v3.normL2()));
	printf("Parecido 2(matr) 3(matr.itr) = %g\n", J4x3n_v2.normL2between(J4x3n_v3) / sqrt(J4x3n_v2.normL2() * J4x3n_v3.normL2()));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();

	printf("Probando encontrar rotacion a partir de varios puntos rotados y jacobianos\n");

	q1.print("q1");
	q3.fijaRotacionCero();
	J4x3n.reallocate(4, NTOT*3);
	J4x3n.setZero();
	J4x3n_v2.reallocate(4, NTOT*3);
	J4x3n.setZero();
	J4x3n_v3.reallocate(4, NTOT*3);
	J4x3n.setZero();
	q3.jacob_fijo_fijaNpuntos_promediaw(fijo, movido, J4x3n, 20);
	if (q3.OP_mult_elementwise(q1) < 0)	{q3.invSigno();	J4x3n*=-1;}
	q3.print("jacob_fijo_fijaNpuntos_promediaw");
	J4x3n.print("J4x3n");
	q3.fijaRotacionCero();
	q3.jacob_fijo_fijaNpuntos_matrRot(fijo, movido, J4x3n_v2);
	if (q3.OP_mult_elementwise(q1) < 0)	{q3.invSigno();	J4x3n_v2*=-1;}
	q3.print("jacob_fijo_fijaNpuntos_matrRot");
	J4x3n_v2.print("J4x3n");
	for (j=0; j<3*NTOT; j++)
	{
		fijo[j/3].el(j%3) += delta;
		q3.fijaNpuntos_matrRot(fijo, movido);
		if (q3.OP_mult_elementwise(q1) < 0)	{q3.invSigno();}
		fijo[j/3].el(j%3) -= 2*delta;
		q4.fijaNpuntos_matrRot(fijo, movido);
		if (q4.OP_mult_elementwise(q1) < 0)	{q4.invSigno();}
		fijo[j/3].el(j%3) += delta;
		for (i=0; i<4; i++)
			J4x3n_v3(i,j) = (q3.el(i)-q4.el(i))/(2*delta);
	}
	J4x3n_v3.print("J4x3n(iter)");

	printf("Parecido 1(jpromw) 2(matr) = %g\n", J4x3n.normL2between(J4x3n_v2) / sqrt(J4x3n.normL2() * J4x3n_v2.normL2()));
	printf("Parecido 1(jpromw) 3(matr.itr) = %g\n", J4x3n.normL2between(J4x3n_v3) / sqrt(J4x3n.normL2() * J4x3n_v3.normL2()));
	printf("Parecido 2(matr) 3(matr.itr) = %g\n", J4x3n_v2.normL2between(J4x3n_v3) / sqrt(J4x3n_v2.normL2() * J4x3n_v3.normL2()));

	printf("\nImprima ENTER para continuar (%d de %d)\n\n", u++, nMax);
	getchar();
}


void L_Quaternion::pruebaCuaternion_def()
{
	L_Quaternion q1, q2, q3, q4;
	L_CoordsCart3D p1, p2, p3;
	L_Matrix J4x4(4,4), J4x4_(4,4), J3x4(3,4), J3x4_(3,4), J3x4_2(3,4);
	L_StaticMatrix<3,4> J3x4_e;
	double pre[36], prem[36];

	q1.fijaAzar();
	q2.fijaAzar();
	q3.fijaAzar();
	q4.fijaAzar();
	p1.fijaAzar();
	p2.fijaAzar();

	printf("Probando funciones versus #defines\n");

	L_Cuaternion_rotarVector(p2, q1, p1);
	p3 = q1.rotarVector(p1);
	p2.print("rotarVector(#def)");
	p3.print("rotarVector(func)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));
	printf("\n");

	L_Cuaternion_rotarVectorInv(p2, q1, p1);
	p3 = q1.rotarVectorInv(p1);
	p2.print("rotarVectorInv(#def)");
	p3.print("rotarVectorInv(func)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));
	printf("\n");

	q1.fijaAzar();
	q2.fijaAzar();
	q3.fijaAzar();
	q4.fijaAzar();
	p1.fijaAzar();
	p2.fijaAzar();

	L_Cuaternion_OP_mult(q3,q1,q2);
	q4.OP_mult(q1,q2);
	q3.print("OP_mult(#def)");
	q4.print("OP_mult(func)");
	printf("Parecido: %f\n", q3.cosineDistanceTo(q4));
	printf("\n");

	L_Cuaternion_inverseOf(q2, q1);
	q3.inverseOf(q1);
	q2.print("inverseOf(#def)");
	q3.print("inverseOf(func)");
	printf("Parecido: %f\n", q2.cosineDistanceTo(q3));
	printf("\n");

	printf("Presione ENTER para continuar\n");
	getchar();

	L_Cuaternion_jacob_izq_OP_mult(q1, q2, J4x4, 0, 0);
	q3.jacob_izq_OP_mult(q1, q2, J4x4_, 0, 0);
	J4x4.print("dOP_MULT/dq1(#def)");
	J4x4_.print("dOP_MULT/dq1(func)");
	printf("Parecido: %f\n", J4x4.cosineDistanceTo(J4x4_));
	printf("\n");


	printf("Presione ENTER para continuar\n");
	getchar();

	L_Cuaternion_jacob_der_OP_mult(q1, q2, J4x4, 0, 0);
	q3.jacob_der_OP_mult(q1, q2, J4x4_, 0, 0);
	J4x4.print("dOP_MULT/dq2(#def)");
	J4x4_.print("dOP_MULT/dq2(func)");
	printf("Parecido: %f\n", J4x4.cosineDistanceTo(J4x4_));
	printf("\n");

	printf("Presione ENTER para continuar\n");
	getchar();

	L_Cuaternion_calc_pre(q1, pre);

	L_Cuaternion_jacob_izq_rotarVector(p1, pre, J3x4, 0, 0);
	L_Quaternion::jacob_izq_rotarVector(p1, pre, J3x4_e);	
	L_Matrix_OP_assign(J3x4_, J3x4_e);
	q1.jacob_izq_rotarVector(p1, J3x4_2);	
	J3x4.print("dRotVect/dq(#def.pre)");
	J3x4_.print("dRotVect/dq(func.pre)");
	J3x4_2.print("dRotVect/dq(func.cal)");
	printf("Parecido: 1 2   %f\n", J3x4.cosineDistanceTo(J3x4_));
	printf("Parecido: 1 3   %f\n", J3x4.cosineDistanceTo(J3x4_2));
	printf("\n");
	
	L_Cuaternion_calc_prem(q1, prem);

	L_Cuaternion_jacob_izq_rotarVectorInv(p1, prem, J3x4, 0, 0);
	L_Quaternion::jacob_izq_rotarVectorInv(p1, prem, J3x4_e);	
	L_Matrix_OP_assign(J3x4_, J3x4_e);
	q1.jacob_izq_rotarVectorInv(p1, J3x4_2);	
	J3x4.print("dRotVectInv/dq(#def.prem)");
	J3x4_.print("dRotVectInv/dq(func.prem)");
	J3x4_2.print("dRotVectInv/dq(func.calc)");
	printf("Parecido: 1 2   %f\n", J3x4.cosineDistanceTo(J3x4_));
	printf("Parecido: 1 3   %f\n", J3x4.cosineDistanceTo(J3x4_2));
	printf("\n");

	printf("Presione ENTER para continuar\n");
	getchar();

}


L_Quaternion L_Rot3DMatrix::pideCuaternion()
{
	L_Quaternion q;
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*q.a 
		q.a = 0.25 * S;
		q.b = (operator()(2,1) - operator()(1,2)) / S;
		q.c = (operator()(0,2) - operator()(2,0)) / S; 
		q.d = (operator()(1,0) - operator()(0,1)) / S; 
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*q.b 
		q.a = (operator()(2,1) - operator()(1,2)) / S;
		q.b = 0.25 * S;
		q.c = (operator()(0,1) + operator()(1,0)) / S; 
		q.d = (operator()(0,2) + operator()(2,0)) / S; 
	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*q.c
		q.a = (operator()(0,2) - operator()(2,0)) / S;
		q.b = (operator()(0,1) + operator()(1,0)) / S; 
		q.c = 0.25 * S;
		q.d = (operator()(1,2) + operator()(2,1)) / S; 
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*q.d
		q.a = (operator()(1,0) - operator()(0,1)) / S;
		q.b = (operator()(0,2) + operator()(2,0)) / S;
		q.c = (operator()(1,2) + operator()(2,1)) / S;
		q.d = 0.25 * S;
	}
	return q;
}


L_Quaternion L_Rot3DMatrix::jacob_pideCuaternion(L_Matrix &J4x9, int i0, int j0)
{
	L_Quaternion q;
	int i, j;
	for (i=0; i<4; i++)
		for (j=0; j<9; j++)
			J4x9(i+i0,j+j0) = 0;
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	int e00=j0+0, e01=j0+1, e02=j0+2, e10=j0+3, e11=j0+4, e12=j0+5, e20=j0+6, e21=j0+7, e22=j0+8;
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*q.a 
		q.a = 0.25 * S;
		q.b = (operator()(2,1) - operator()(1,2)) / S;
		q.c = (operator()(0,2) - operator()(2,0)) / S;
		q.d = (operator()(1,0) - operator()(0,1)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = A; J4x9(i0+0,e11) = A; J4x9(i0+0,e22) = A;
		J4x9(i0+1,e00) = -q.b*C; J4x9(i0+1,e11) = -q.b*C; J4x9(i0+1,e22) = -q.b*C;
		J4x9(i0+1,e21) = B; J4x9(i0+1,e12) = -B;
		J4x9(i0+2,e00) = -q.c*C; J4x9(i0+2,e11) = -q.c*C; J4x9(i0+2,e22) = -q.c*C;
		J4x9(i0+2,e02) = B; J4x9(i0+2,e20) = -B;
		J4x9(i0+3,e00) = -q.d*C; J4x9(i0+3,e11) = -q.d*C; J4x9(i0+3,e22) = -q.d*C;
		J4x9(i0+3,e10) = B; J4x9(i0+3,e01) = -B;
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*q.b 
		q.a = (operator()(2,1) - operator()(1,2)) / S;
		q.b = 0.25 * S;
		q.c = (operator()(0,1) + operator()(1,0)) / S;
		q.d = (operator()(0,2) + operator()(2,0)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = -q.a*C; J4x9(i0+0,e11) = q.a*C; J4x9(i0+0,e22) = q.a*C;
		J4x9(i0+0,e21) = B; J4x9(i0+0,e12) = -B; // ERROR
		J4x9(i0+1,e00) = A; J4x9(i0+1,e11) = J4x9(i0+1,e22) = -A;
		J4x9(i0+2,e00) = -q.c*C; J4x9(i0+2,e11) = q.c*C;J4x9(i0+2,e22) = q.c*C;
		J4x9(i0+2,e01) = B; J4x9(i0+2,e10) = B;
		J4x9(i0+3,e00) = -q.d*C; J4x9(i0+3,e11) = q.d*C; J4x9(i0+3,e22) = q.d*C;
		J4x9(i0+3,e02) = B; J4x9(i0+3,e20) = B;

	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*q.c
		q.a = (operator()(0,2) - operator()(2,0)) / S;
		q.b = (operator()(0,1) + operator()(1,0)) / S;
		q.c = 0.25 * S;
		q.d = (operator()(1,2) + operator()(2,1)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = q.a*C; J4x9(i0+0,e11) = -q.a*C; J4x9(i0+0,e22) = q.a*C;
		J4x9(i0+0,e02) = B; J4x9(i0+0,e20) = -B;
		J4x9(i0+1,e00) = q.b*C; J4x9(i0+1,e11) = -q.b*C; J4x9(i0+1,e22) = q.b*C;
		J4x9(i0+1,e01) = B; J4x9(i0+1,e10) = B;
		J4x9(i0+2,e00) = -A; J4x9(i0+2,e11) = A; J4x9(i0+2,e22) = -A;
		J4x9(i0+3,e00) = q.d*C; J4x9(i0+3,e11) = -q.d*C; J4x9(i0+3,e22) = q.d*C;
		J4x9(i0+3,e12) = B; J4x9(i0+3,e21) = B;
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*q.d
		q.a = (operator()(1,0) - operator()(0,1)) / S;
		q.b = (operator()(0,2) + operator()(2,0)) / S;
		q.c = (operator()(1,2) + operator()(2,1)) / S;
		q.d = 0.25 * S;
		double rc = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1));
		double r32 = (1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1))*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = q.a*C; J4x9(i0+0,e11) = q.a*C; J4x9(i0+0,e22) = -q.a*C;
		J4x9(i0+0,e10) = B; J4x9(i0+0,e01) = -B;
		J4x9(i0+1,e00) = q.b*C; J4x9(i0+1,e11) = q.b*C; J4x9(i0+1,e22) = -q.b*C;
		J4x9(i0+1,e02) = B; J4x9(i0+1,e20) = B;
		J4x9(i0+2,e00) = q.c*C; J4x9(i0+2,e11) = q.c*C;J4x9(i0+2,e22) = -q.c*C;
		J4x9(i0+2,e12) = B; J4x9(i0+2,e21) = B;
		J4x9(i0+3,e00) = -A; J4x9(i0+3,e11) = -A; J4x9(i0+3,e22) = A;
	}
	return q;
}

L_Quaternion L_Rot3DMatrix::pideCuaternionInv()
{
	L_Quaternion q;
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*q.a 
		q.a = 0.25 * S;
		q.b = (operator()(1,2) - operator()(2,1)) / S;
		q.c = (operator()(2,0) - operator()(0,2)) / S; 
		q.d = (operator()(0,1) - operator()(1,0)) / S; 
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*q.b 
		q.a = (operator()(1,2) - operator()(2,1)) / S;
		q.b = 0.25 * S;
		q.c = (operator()(1,0) + operator()(0,1)) / S; 
		q.d = (operator()(2,0) + operator()(0,2)) / S; 
	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*q.c
		q.a = (operator()(2,0) - operator()(0,2)) / S;
		q.b = (operator()(1,0) + operator()(0,1)) / S; 
		q.c = 0.25 * S;
		q.d = (operator()(2,1) + operator()(1,2)) / S; 
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*q.d
		q.a = (operator()(0,1) - operator()(1,0)) / S;
		q.b = (operator()(2,0) + operator()(0,2)) / S;
		q.c = (operator()(2,1) + operator()(1,2)) / S;
		q.d = 0.25 * S;
	}
	return q;
}


L_Quaternion L_Rot3DMatrix::jacob_pideCuaternionInv(L_Matrix &J4x9, int i0, int j0)
{
	L_Quaternion q;
	int i, j;
	for (i=0; i<4; i++)
		for (j=0; j<9; j++)
			J4x9(i+i0,j+j0) = 0;
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	int e00=j0+0, e01=j0+1, e02=j0+2, e10=j0+3, e11=j0+4, e12=j0+5, e20=j0+6, e21=j0+7, e22=j0+8;

	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*q.a 
		q.a = 0.25 * S;
		q.b = (operator()(1,2) - operator()(2,1)) / S;
		q.c = (operator()(2,0) - operator()(0,2)) / S;
		q.d = (operator()(0,1) - operator()(1,0)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = A; J4x9(i0+0,e11) = A; J4x9(i0+0,e22) = A;
		J4x9(i0+1,e00) = -q.b*C; J4x9(i0+1,e11) = -q.b*C; J4x9(i0+1,e22) = -q.b*C;
		J4x9(i0+1,e12) = B; J4x9(i0+1,e21) = -B;
		J4x9(i0+2,e00) = -q.c*C; J4x9(i0+2,e11) = -q.c*C; J4x9(i0+2,e22) = -q.c*C;
		J4x9(i0+2,e20) = B; J4x9(i0+2,e02) = -B;
		J4x9(i0+3,e00) = -q.d*C; J4x9(i0+3,e11) = -q.d*C; J4x9(i0+3,e22) = -q.d*C;
		J4x9(i0+3,e01) = B; J4x9(i0+3,e10) = -B;
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*q.b 
		q.a = (operator()(1,2) - operator()(2,1)) / S;
		q.b = 0.25 * S;
		q.c = (operator()(1,0) + operator()(0,1)) / S;
		q.d = (operator()(2,0) + operator()(0,2)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = -q.a*C; J4x9(i0+0,e11) = q.a*C; J4x9(i0+0,e22) = q.a*C;
		J4x9(i0+0,e12) = B; J4x9(i0+0,e21) = -B; // ERROR
		J4x9(i0+1,e00) = A; J4x9(i0+1,e11) = J4x9(i0+1,e22) = -A;
		J4x9(i0+2,e00) = -q.c*C; J4x9(i0+2,e11) = q.c*C;J4x9(i0+2,e22) = q.c*C;
		J4x9(i0+2,e10) = B; J4x9(i0+2,e01) = B;
		J4x9(i0+3,e00) = -q.d*C; J4x9(i0+3,e11) = q.d*C; J4x9(i0+3,e22) = q.d*C;
		J4x9(i0+3,e20) = B; J4x9(i0+3,e02) = B;

	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*q.c
		q.a = (operator()(2,0) - operator()(0,2)) / S;
		q.b = (operator()(1,0) + operator()(0,1)) / S;
		q.c = 0.25 * S;
		q.d = (operator()(2,1) + operator()(1,2)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = q.a*C; J4x9(i0+0,e11) = -q.a*C; J4x9(i0+0,e22) = q.a*C;
		J4x9(i0+0,e20) = B; J4x9(i0+0,e02) = -B;
		J4x9(i0+1,e00) = q.b*C; J4x9(i0+1,e11) = -q.b*C; J4x9(i0+1,e22) = q.b*C;
		J4x9(i0+1,e10) = B; J4x9(i0+1,e01) = B;
		J4x9(i0+2,e00) = -A; J4x9(i0+2,e11) = A; J4x9(i0+2,e22) = -A;
		J4x9(i0+3,e00) = q.d*C; J4x9(i0+3,e11) = -q.d*C; J4x9(i0+3,e22) = q.d*C;
		J4x9(i0+3,e21) = B; J4x9(i0+3,e12) = B;
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*q.d
		q.a = (operator()(0,1) - operator()(1,0)) / S;
		q.b = (operator()(2,0) + operator()(0,2)) / S;
		q.c = (operator()(2,1) + operator()(1,2)) / S;
		q.d = 0.25 * S;
		double rc = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1));
		double r32 = (1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1))*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J4x9(i0+0,e00) = q.a*C; J4x9(i0+0,e11) = q.a*C; J4x9(i0+0,e22) = -q.a*C;
		J4x9(i0+0,e01) = B; J4x9(i0+0,e10) = -B;
		J4x9(i0+1,e00) = q.b*C; J4x9(i0+1,e11) = q.b*C; J4x9(i0+1,e22) = -q.b*C;
		J4x9(i0+1,e20) = B; J4x9(i0+1,e02) = B;
		J4x9(i0+2,e00) = q.c*C; J4x9(i0+2,e11) = q.c*C;J4x9(i0+2,e22) = -q.c*C;
		J4x9(i0+2,e21) = B; J4x9(i0+2,e12) = B;
		J4x9(i0+3,e00) = -A; J4x9(i0+3,e11) = -A; J4x9(i0+3,e22) = A;
	}
	return q;
}
L_CoordsCart3D L_Rot3DMatrix::pideVectorRotacion()
{
	L_Quaternion q;
	L_CoordsCart3D eje;
	q = pideCuaternion();
	return q.pideVectorRotacion();
}

L_PanTiltRoll L_Rot3DMatrix::pidePanTiltRoll()
{
	L_PanTiltRoll ori;
	double ct=sqrt(operator()(0,0)*operator()(0,0)+operator()(1,0)*operator()(1,0));
	if (ct > 0)
	{
		ori.pan=atan2(operator()(1,0),operator()(0,0));
		ori.tilt=atan2(operator()(2,0),ct);
		ori.roll=atan2(operator()(2,1),operator()(2,2));
	}
	else
	{
		ori.roll=0;
		if (operator()(2,0) > 0)
		{
			ori.tilt=M_PI;
			ori.pan = atan2(operator()(2,1),operator()(1,1));
		}
		else
		{
			ori.tilt=-M_PI;
			ori.pan = atan2(operator()(2,1),operator()(1,1));
		}
	}
	return ori;
}

void L_Rot3DMatrix::transpOf(const L_Rot3DMatrix &other)
{
	operator()(0,0) = other(0,0);
	operator()(0,1) = other(1,0);
	operator()(0,2) = other(2,0);
	operator()(1,0) = other(0,1);
	operator()(1,1) = other(1,1);
	operator()(1,2) = other(2,1);
	operator()(2,0) = other(0,2);
	operator()(2,1) = other(1,2);
	operator()(2,2) = other(2,2);
}


void L_Rot3DMatrix::OP_mult(const L_Rot3DMatrix &m1, const L_Rot3DMatrix &m2)
{
#ifdef _L_MatrizRot3D_ELEM_MULT
#error _L_MatrizRot3D_ELEM_MULT(A,B) ya definido
#endif
#define _L_MatrizRot3D_ELEM_MULT(A,B) operator()(A,B) = m1(A,0)*m2(0,B) + m1(A,1)*m2(1,B) + m1(A,2) * m2(2,B)

	_L_MatrizRot3D_ELEM_MULT(0,0);
	_L_MatrizRot3D_ELEM_MULT(0,1);
	_L_MatrizRot3D_ELEM_MULT(0,2);
	_L_MatrizRot3D_ELEM_MULT(1,0);
	_L_MatrizRot3D_ELEM_MULT(1,1);
	_L_MatrizRot3D_ELEM_MULT(1,2);
	_L_MatrizRot3D_ELEM_MULT(2,0);
	_L_MatrizRot3D_ELEM_MULT(2,1);
	_L_MatrizRot3D_ELEM_MULT(2,2);

#undef _L_MatrizRot3D_ELEM_MULT
}

void L_Rot3DMatrix::OP_mult(const L_Matrix &m1, const L_Matrix &m2)
{
#ifdef _L_MatrizRot3D_ELEM_MULT
#error _L_MatrizRot3D_ELEM_MULT(A,B) ya definido
#endif
#define _L_MatrizRot3D_ELEM_MULT(A,B) operator()(A,B) = m1(A,0)*m2(0,B) + m1(A,1)*m2(1,B) + m1(A,2) * m2(2,B)

	if (m1.li != 3 || m1.lj != 3 || m2.li != 3 || m2.lj != 3)
		throw L_ArgException();
	_L_MatrizRot3D_ELEM_MULT(0,0);
	_L_MatrizRot3D_ELEM_MULT(0,1);
	_L_MatrizRot3D_ELEM_MULT(0,2);
	_L_MatrizRot3D_ELEM_MULT(1,0);
	_L_MatrizRot3D_ELEM_MULT(1,1);
	_L_MatrizRot3D_ELEM_MULT(1,2);
	_L_MatrizRot3D_ELEM_MULT(2,0);
	_L_MatrizRot3D_ELEM_MULT(2,1);
	_L_MatrizRot3D_ELEM_MULT(2,2);

#undef _L_MatrizRot3D_ELEM_MULT
}

L_Pose3D &L_Pose3D::fijaPose_abs(double x, double y, double z, double xAde, double yAde, double zAde, double xArr, double yArr, double zArr)
{
	double r, p, yIzq, zIzq;
	xAde-=x; yAde-=y; zAde-=z; // Sacar offset de Ade y Arr
	xArr-=x; yArr-=y; zArr-=z;
	r = 1/sqrt(xAde*xAde+yAde*yAde+zAde*zAde);
	xAde*=r; yAde*=r; zAde*=r;// Normalizar Ade
	p = xArr*xAde + yArr*yAde + zArr*zAde;
	xArr-=p; yArr-=p; zArr-=p; // Gram-schmidt para que Arr quede perpendicular a Ade, Ade manda...
	r = 1/sqrt(xArr*xArr+yArr*yArr+zArr*zArr);
	xArr*=r; yArr*=r; zArr*=r; // Normalizar Arr
	yIzq = zArr*xAde; // Izq = Arr cruz Ade, los 3 son unitarios. Matriz rot = [Ade;Izq;Arr)
	zIzq = xArr*yAde;

	pos.x = x;
	pos.y = y;
	pos.z = z;

	r = sqrt(xAde*xAde+yAde*yAde);
	// Esto de aca copiado de L_HomogeneousMatrix::calcPose3D_fijo_a_movido()
	if (r > 0)
	{
		ori.pan = atan2(yAde,xAde);
		ori.tilt = atan2(zAde,r); // Entre -90° y 90°
		ori.roll = atan2(zIzq,zArr);
	}
	else
	{
		ori.roll = 0;
		if (xArr > 0)
		{
			ori.tilt = M_PI;
			ori.pan = -atan2(yArr,yIzq);
		}
		else
		{
			ori.tilt=-M_PI;
			ori.pan = atan2(yArr,yIzq);
		}
	}
	return *this;
}

L_Pose3D &L_Pose3D::fijaPose_rel(double x, double y, double z, double xDirAde, double yDirAde, double zDirAde, double xDirArr, double yDirArr, double zDirArr)
{
	double r, p, yIzq, zIzq;
	r = 1/sqrt(xDirAde*xDirAde+yDirAde*yDirAde+zDirAde*zDirAde);
	xDirAde*=r; yDirAde*=r; zDirAde*=r;// Normalizar Ade
	p = xDirArr*xDirAde + yDirArr*yDirAde + zDirArr*zDirAde;
	xDirArr-=p; yDirArr-=p; zDirArr-=p; // Gram-schmidt para que Arr quede perpendicular a Ade, Ade manda...
	r = 1/sqrt(xDirArr*xDirArr+yDirArr*yDirArr+zDirArr*zDirArr);
	xDirArr*=r; yDirArr*=r; zDirArr*=r; // Normalizar Arr
	yIzq = zDirArr*xDirAde; // Izq = Arr cruz Ade, los 3 son unitarios. Matriz rot = [Ade;Izq;Arr)
	zIzq = xDirArr*yDirAde;

	pos.x = x;
	pos.y = y;
	pos.z = z;

	r = sqrt(xDirAde*xDirAde+yDirAde*yDirAde);
	// Esto de aca copiado de L_HomogeneousMatrix::calcPose3D_fijo_a_movido()
	if (r > 0)
	{
		ori.pan = atan2(yDirAde,xDirAde);
		ori.tilt = atan2(zDirAde,r); // Entre -90° y 90°
		ori.roll = atan2(zIzq,zDirArr);
	}
	else
	{
		ori.roll = 0;
		if (xDirArr > 0)
		{
			ori.tilt = M_PI;
			ori.pan = -atan2(yDirArr,yIzq);
		}
		else
		{
			ori.tilt=-M_PI;
			ori.pan = atan2(yDirArr,yIzq);
		}
	}
	return *this;
}

void L_Pose3D::agregaDeltaPose(const L_Pose3D &deltaPose)
{
	L_HomogeneousMatrix H20, H21, H10;
	H10.fijaPose3D(*this);
	H21.fijaPose3D(deltaPose);
	H20.OP_mult(H21,H10); // nueva_acum = delta * antigua_acum
	*this = H20.calcPose3D();
}

void L_Pose3D::pruebaPose3D()
{
	L_Pose3D p1, p2, p3;
	L_HomogeneousMatrix H1, H2, H3;
	L_CoordsCart3D c1, c2;

	p1.fijaPose(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI));
	p2.fijaPose(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI));
	p3.fijaPose(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI), L_RANDOM(-M_PI, M_PI));

	printf("fijaPose3D v/s fijaPose3D_lento\n");

	H2.fijaPose3D(p2);
	H2.print("H2");
	H2.fijaPose3D_lento(p2);
	H2.print("H2");

	printf("\n");
	printf("fijaPose3D_fijo_a_movido v/s fijaPose3D_fijo_a_movido_lento\n");

	H3.fijaPose3D_fijo_a_movido(p3);
	H3.print("H3");
	H3.fijaPose3D_fijo_a_movido_lento(p3);
	H3.print("H3");

	printf("\n");
	printf("L_HomogeneousMatrix::* v/s agregaDeltaPose\n");

	H1.fijaPose3D(p1);
	H2.fijaPose3D(p2);

	H3 = H2*H1;
	p3 = H3.calcPose3D();
	p3.print("p3 (homog)");

	p3 = p1;
	p3.agregaDeltaPose(p2); // p3 = p1 -> p2
	p3.print("p3 (delta)");

	printf("\n");

	p1.print("p1");
	p1 = H1.calcPose3D();
	p1.print("p1");

	H1.fijaPose3D_fijo_a_movido(p1);
	p1 = H1.calcPose3D_fijo_a_movido();
	p1.print("p1");


	printf("Presione ENTER para continuar\n");
	getchar();
}

void L_Pose3D_cuat::imprimeGrados(const char *name, FILE *fp)
{
	L_PanTiltRoll ptt;
	L_Rot3DMatrix rot;
	rot.fijarCuaternion(ori);
	ptt=rot.pidePanTiltRoll();
	fprintf(fp, "%s = [%.4g %.4g %.4g) pan=%.4fgr tilt=%.4ggr roll=%.4ggr\n\n", name, pos.x, pos.y, pos.z, ptt.pan*180/M_PI, ptt.tilt*180/M_PI, ptt.roll*180/M_PI);
}

void L_Pose3D_cuat::OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	// eta1*eta2 es como : pc( H(eta2)*H(eta1) ) = pc( [R2 T2)*[R1 T1) ) = pc( [R2R1 R2T1+T2)) -> q1q2 ; q2't1q2 + t2
	// VA PRIMERO POS Y DESPUES ORI
	//pos = eta2.ori.rotarVector(eta1.pos) + eta2.pos;
	//ori = eta1.ori * eta2.ori;
	//
	// Ahora  pc( H(eta1)*H(eta2) ) = pc( [R1 T1)*[R2 T2) ) = pc( [R1R2 R1T2+T1))
	pos = eta1.ori.rotarVector(eta2.pos) + eta1.pos;
	ori = eta1.ori * eta2.ori;
}

void L_Pose3D_cuat::OP_div(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	L_Pose3D_cuat eta2Inv;
	eta2Inv.inverseOf(eta2);
	pos = eta1.ori.rotarVector(eta2Inv.pos) + eta1.pos;
	ori = eta1.ori * eta2Inv.ori;
}

void L_Pose3D_cuat::jacob_izq_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
	//pos = eta1.ori.rotarVector(eta2.pos) + eta1.pos;
	//ori = eta1.ori * eta2.ori;
	int i, j;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_izq_OP_mult() : matriz de tamano insuficiente");
	for (i=0; i<7; i++)
		for (j=0; j<7; j++)
			J7x7(i+i0,j+j0) = 0;
	// dpos / dp1.pos
	J7x7(i0+0,j0+0) = 1;
	J7x7(i0+1,j0+1) = 1;
	J7x7(i0+2,j0+2) = 1;
	// dpos / dp1.ori
	eta1.ori.jacob_izq_rotarVector(eta2.pos, J7x7, i0, j0+3);//  eta2.ori.jacob_der_rotarVector(eta1.pos, J7x7, 0, 0);// dpos / dp1.pos
	// dori / dp1.pos
	;
	// dori / dp1.ori
	L_Quaternion::jacob_izq_OP_mult(eta1.ori, eta2.ori, J7x7, i0+3, j0+3); // dori / dp1.ori 4x4 llenado
}

void L_Pose3D_cuat::jacob_der_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
	//pos = eta1.ori.rotarVector(eta2.pos) + eta1.pos;
	//ori = eta1.ori * eta2.ori;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_der_OP_mult() : matriz de tamano insuficiente");
	int i, j;
	for (i=0; i<7; i++)
		for (j=0; j<7; j++)
			J7x7(i+i0,j+j0) = 0;

	// dpos / dp2.pos
	eta1.ori.jacob_der_rotarVector(J7x7, i0, j0);//  eta2.ori.jacob_der_rotarVector(eta1.pos, J7x7, 0, 0);// dpos / dp1.pos
	// dpos / dp2.ori
	;
	// dori / dp2.pos
	;
	// dori / dp2.ori
	L_Quaternion::jacob_der_OP_mult(eta1.ori, eta2.ori, J7x7, i0+3, j0+3);
}

void L_Pose3D_cuat::jacob_izq_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, const double *pre_p2_36, L_StaticMatrix<7,7> &J7x7)
{
	int i, j;
	for (i=0; i<7; i++)
		for (j=0; j<7; j++)
			J7x7(i,j) = 0;

	// dpos / dp2.pos
	J7x7(0,0) = 1;
	J7x7(1,1) = 1;
	J7x7(2,2) = 1;
	// dpos / dp2.ori
	L_Cuaternion_jacob_izq_rotarVector(eta2.pos, pre_p2_36, J7x7, 0, 3); // Error: decia eta1.pos
	// dori / dp2.pos
	;
	// dori / dp2.ori
	L_Cuaternion_jacob_izq_OP_mult(eta1.ori, eta2.ori, J7x7, 3, 3);
}

void L_Pose3D_cuat::jacob_der_OP_mult(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_StaticMatrix<7,7> &J7x7)
{
	int i, j;
	for (i=0; i<7; i++)
		for (j=0; j<7; j++)
			J7x7(i,j) = 0;
	// dpos / dp1.pos
	L_Cuaternion_jacob_der_rotarVector(eta2.ori, J7x7, 0, 0);//  eta2.ori.jacob_der_rotarVector(eta1.pos, J7x7, 0, 0);// dpos / dp1.pos
	// dpos / dp1.ori
	;
	// dori / dp1.pos
	;
	// dori / dp1.ori
	L_Cuaternion_jacob_der_OP_mult(eta1.ori, eta2.ori, J7x7, 3, 3); // dori / dp1.ori 4x4 llenado
}
void L_Pose3D_cuat::OP_suma_punto_a_punto(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	pos.x = eta1.pos.x + eta2.pos.x;
	pos.y = eta1.pos.y + eta2.pos.y;
	pos.z = eta1.pos.z + eta2.pos.z;
	ori.a = eta1.ori.a + eta2.ori.a;
	ori.b = eta1.ori.b + eta2.ori.b;
	ori.c = eta1.ori.c + eta2.ori.c;
	ori.d = eta1.ori.d + eta2.ori.d;
}

void L_Pose3D_cuat::OP_resta_punto_a_punto(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	pos.x = eta1.pos.x - eta2.pos.x;
	pos.y = eta1.pos.y - eta2.pos.y;
	pos.z = eta1.pos.z - eta2.pos.z;
	ori.a = eta1.ori.a - eta2.ori.a;
	ori.b = eta1.ori.b - eta2.ori.b;
	ori.c = eta1.ori.c - eta2.ori.c;
	ori.d = eta1.ori.d - eta2.ori.d;
}

void L_Pose3D_cuat::mirarA(double xMirado, double yMirado, double zMirado, double roll)
{
	L_HomogeneousMatrix H;
	L_Pose3D p3d;
	L_CoordsCart3D ade(xMirado, yMirado, zMirado), arr(pos.x, pos.y, pos.z + 1000);
	if (xMirado == pos.x && yMirado == pos.y)
	{
		if (zMirado > pos.z)
			definePosPanTiltRoll(pos.x, pos.y, pos.z, 0, 90*M_PI/180, 0);
		else if (zMirado < pos.z)
			definePosPanTiltRoll(pos.x, pos.y, pos.z, 0, -90*M_PI/180, 0);
		else
			throw_L_ArgException_if(true, "L_Pose3D_cuat::mirarA() : punto mirado = punto base");
		return;
	}
	H.fijaOrigenAdeArr(pos, ade-pos, arr-pos);
	p3d = H.calcPose3D();
	p3d.ori.roll = 0;
	H.fijaPose3D(p3d);
	*this = H.calcPose3D_cuat();
}

void L_Pose3D_cuat::fijaPoseMirarA(double x, double y, double z, double xMirado, double yMirado, double zMirado, double roll)
{
	pos.x = x;
	pos.y = y;
	pos.z = z;
	mirarA(xMirado, yMirado, zMirado, roll);
}

void L_Pose3D_cuat::definePosPanTiltRoll(double x, double y, double z, double pan, double tilt, double roll)
{
	L_Pose3D p;
	L_Rot3DMatrix rot;
	p.fijaPose(x,y,z,pan,tilt,roll);
	pos=p.pos;
	ori = rot.fijarPanTiltRoll(p.ori).pideCuaternion();
}

void L_Pose3D_cuat::inverseOf(const L_Pose3D_cuat &eta)
{
	ori.inverseOf(eta.ori);
	pos = -1 * ori.rotarVector(eta.pos); // ori ya es la rotacion inversa
}

void L_Pose3D_cuat::jacob_inverseOf(const L_Pose3D_cuat &eta, L_Matrix &J7x7)
{
	int i, j;
	J7x7.reallocate(7,7);

	for (i=0; i<7; i++)
		for (j=0; j<7; j++)
			J7x7(i,j) = 0;

	L_Quaternion qinv;
	L_Matrix mJrotI, mJrotD, Jinv;

	// pos = -rot(I=inv(eta.ori),D=eta.pos)
	// ori = inv(eta.ori)

	// dpos/dp.pos = -(drot/dD)|eta.pos
	// dpos/dp.ori = -(dRot/dI)|inv(eta.ori) * (dinv/dq)|eta.ori
	// dori/dp.pos = 0
	// dori/dp.ori = (dinv/dq)|eta.ori

	Jinv.reallocate(4,4);
	mJrotI.reallocate(3,4);
	mJrotD.reallocate(3,3);
	qinv.inverseOf(eta.ori);
	qinv.jacob_inverseOf(eta.ori, Jinv);
	pos = -qinv.rotarVector(eta.pos);
	ori = qinv;
	qinv.jacob_izq_rotarVector(eta.pos, mJrotI);
	qinv.jacob_der_rotarVector(mJrotD);
	mJrotI*=-1;
	mJrotD*=-1;

	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			J7x7(i,j) = mJrotD(i,j);
	J7x7.OP_mult_sub(0, 3, mJrotI, Jinv, 0, 0, 3, 4, 0,  0, 4, 4);
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			J7x7(i+3,j+3) = Jinv(i,j);
}

void L_Pose3D_cuat::OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	L_Pose3D_cuat p2inv;
	p2inv.inverseOf(eta2);
	OP_mult(eta1, p2inv);
}

void L_Pose3D_cuat::jacob_izq_OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
	L_Pose3D_cuat p2inv;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_izq_OP_mult_p1_por_p2inv() : matriz de tamano insuficiente");
	p2inv.inverseOf(eta2);
	jacob_izq_OP_mult(eta1, p2inv, J7x7, i0, j0); // FALTABA ",i0, j0"
}

void L_Pose3D_cuat::jacob_der_OP_mult_p1_por_p2inv(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
	L_Pose3D_cuat p2inv;
	L_Matrix J7x7i, J7x7m;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_der_OP_mult_p1_por_p2inv() : matriz de tamano insuficiente");
	J7x7i.reallocate(7,7);
	J7x7m.reallocate(7,7);
	p2inv.jacob_inverseOf(eta2, J7x7i);
	jacob_der_OP_mult(eta1, p2inv, J7x7m);
	J7x7.OP_mult_sub(i0, j0, J7x7m, J7x7i, 0,0,7,7, 0,0,7,7);
}

void L_Pose3D_cuat::OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2)
{
	L_Pose3D_cuat p1inv;
	p1inv.inverseOf(eta1);
	OP_mult(p1inv, eta2);
}

void L_Pose3D_cuat::jacob_izq_OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
	L_Pose3D_cuat p1inv;
	L_Matrix J7x7i, J7x7m;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_izq_OP_mult_p1inv_por_p2() : matriz de tamano insuficiente");
	J7x7i.reallocate(7,7);
	J7x7m.reallocate(7,7);
	p1inv.jacob_inverseOf(eta1, J7x7i);
	jacob_izq_OP_mult(p1inv, eta2, J7x7m);
	J7x7.OP_mult_sub(i0, j0, J7x7m, J7x7i, 0,0,7,7, 0,0,7,7);
}

void L_Pose3D_cuat::jacob_der_OP_mult_p1inv_por_p2(const L_Pose3D_cuat &eta1, const L_Pose3D_cuat &eta2, L_Matrix &J7x7, int i0, int j0)
{
// Error: no esta i0, j0
	L_Pose3D_cuat p1inv;
	throw_L_ArgException_if(J7x7.li < i0+7 || J7x7.lj < j0+7, "L_Pose3D_cuat::jacob_der_OP_mult_p1inv_por_p2() : matriz de tamano insuficiente");
	p1inv.inverseOf(eta1);
	jacob_der_OP_mult(p1inv, eta2, J7x7, i0, j0); // Error: faltaba i0, j0
}

void L_Pose3D_cuat::pow(const L_Pose3D_cuat &pose, double exp)
{
	if (pose.ori.a != 0 || pose.ori.b != 0 || pose.ori.c != 0 || pose.ori.d != 0)
	{
		// Calcular potencia
		L_Rot3DMatrix rot;
		L_Matrix Id(3,3), m(3,3), t(3,1), tf(3,1);
		L_Pose3D_cuat p, d, pinv;
		rot.fijarCuaternion(pose.ori);
		L_Matrix_OP_subtract(m, Id, rot);
		if (m.invertMe() == false)
			{throw_L_ArgException_if (true, "L_Pose3D_cuat::pow() : I-R(q) no invertible");}
		t(0,0) = pose.pos.x;
		t(1,0) = pose.pos.y;
		t(2,0) = pose.pos.z;
		tf.OP_mult(m,t);
		p.pos.x = tf(0,0);
		p.pos.y = tf(1,0);
		p.pos.z = tf(2,0);
		p.ori.fijaRotacionCero();
		pinv.inverseOf(p);
		d.pos.fijaCero();
		d.ori.pow(ori, exp);
		*this = p*d*pinv;
	}
	else
	{
		ori.fijaRotacionCero();
		pos = pose.pos*exp;
	}
}

void L_Pose3D_cuat::movCov(L_StaticMatrix<3,3> &res, const L_StaticMatrix<3,3> &orig)
{
	L_Rot3DMatrix rot, rotT, tmp;
	rot.fijarCuaternion(ori);
	rotT.transpOf(rot);
	L_StaticMatrix_OP_mult(tmp, rot, orig);
	L_StaticMatrix_OP_mult(res, tmp, rotT);
}

void L_Pose3D_cuat::moviCov(L_StaticMatrix<3,3> &res, const L_StaticMatrix<3,3> &orig)
{
	L_Rot3DMatrix rot, rotT, tmp;
	rot.fijarCuaternion(ori);
	rotT.transpOf(rot);
	L_StaticMatrix_OP_mult(tmp, rotT, orig);
	L_StaticMatrix_OP_mult(res, tmp, rot);
}

bool L_Pose3D_cuat::promedioIR(const L_Pose3D_cuat &other, double factor)
{
	pos = (1-factor)*pos + factor*other.pos;
	ori.a = (1-factor)*ori.a + factor*other.ori.a;
	ori.b = (1-factor)*ori.b + factor*other.ori.b;
	ori.c = (1-factor)*ori.c + factor*other.ori.c;
	ori.d = (1-factor)*ori.d + factor*other.ori.d;
	return normaliza_ret();
}


L_Pose3D L_Pose3D_cuat::calcPose3D() const
{
	L_Pose3D p;
	L_HomogeneousMatrix h;
	L_MatrizHomogenea_fijaPose3D_cuat(h,*this);
	p = h.calcPose3D();
	return p;
}


void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2_ant(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, const double *pre36_pr, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const
{
	//f(pCM) = proy(fijo_mov(mult(pCM,pC2_C),fijo))
	//df/dpr = dproy/dfm * dfijo_mov/dmult * dmult/dpr
	//                         f(pre36)      (pre36pr)
	//df/fijo = dprfm(pCM*pC2_C)/dfijo
	//                J3x3s
	L_StaticMatrix<2,3> J2x3t; //dproy/dfm
	L_StaticMatrix<3,7> J3x7t; //dfijo_mov/dmult
	L_StaticMatrix<7,7> J7x7t; //dmult/dpr
	L_StaticMatrix<3,3> J3x3t;
	L_StaticMatrix<2,7> J2x3tJ3x7t;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x3.li < i0+2 || J2x3.lj < j0+3, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2_ant() : matriz de tamano insuficiente");

	L_CoordsCart3D movido;
	L_Pose3D_cuat_fijo_a_movido(movido,*this,fijo);
	L_RayoCamara_jacob_proyeccionDe(movido,J2x3t,0,0); // dproy/dfm
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(*this,fijo,pre36m,J3x3s,J3x7t,0,0); // dfm/dmult

	// OJO
	//
	// (pCM*pC2_C) -> jfm_pre -> J3x3s, pre36m
	// pCM  ->  jimult(pC2_C) -> pre36_pr

	//                               dmult/dpr
	jacob_izq_OP_mult(pCM,pC2_C,pre36_pr,J7x7t);
	//                              dproy/dfm dfm/dmult
	L_StaticMatrix_OP_mult(J2x3tJ3x7t,J2x3t,J3x7t);
	//                                        dproy/dmult dmult/dpr
	L_StaticMatrix_OP_mult_sub(J2x7, i0, j0, J2x3tJ3x7t, J7x7t, 0, 0, 2, 7, 0, 0, 7, 7);

	L_StaticMatrix_OP_mult_sub(J2x3, i1, j1, J2x3t,J3x3s,0,0,2,3,0,0,3,3);
}


void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2_ant(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &J3x3s, const double *pre36m, const double *pre36_pr, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const
{
	//f(pCM) = proy(fijo_mov(mult(pCM,pC2_C),c3d(fijo)))
	//df/dpr = dproy/dfm * dfijo_mov/dmult * dmult/dpr
	//                         f(pre36)     (pre36_pr)
	//df/fijo = dprfm(pCM*pC2_C)/dpfijo * dc3d/dfijo
	//                J3x3s              J3x6t
	L_StaticMatrix<2,3> J2x3t; //dproy/dfm
	L_StaticMatrix<3,7> J3x7t; //dfijo_mov/dmult
	L_StaticMatrix<7,7> J7x7t; //dmult/dpr
	L_StaticMatrix<2,7> J2x3tJ3x7t;
	L_StaticMatrix<3,6> J3x6t;
	L_StaticMatrix<2,3> J2x3tJ3x3t;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x6.li < i0+2 || J2x6.lj < j0+6, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2_ant() : matriz de tamano insuficiente");

	L_CoordsCart3D fijov, movido;
	fijov = fijo.pideCoords3D();
	L_Pose3D_cuat_fijo_a_movido(movido,*this,fijov);
	L_RayoCamara_jacob_proyeccionDe(movido,J2x3t,0,0);
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(*this,fijov,pre36m,J3x3s,J3x7t,0,0); //

	// OJO
	//
	// (pCM*pC2_C) -> jfm_pre -> J3x3s, pre36m
	// pCM  ->  jimult(pC2_C) -> pre36_pr

	jacob_izq_OP_mult(pCM,pC2_C,pre36_pr,J7x7t);
	L_StaticMatrix_OP_mult(J2x3tJ3x7t,J2x3t,J3x7t);
	L_StaticMatrix_OP_mult_sub(J2x7, i0, j0, J2x3tJ3x7t, J7x7t, 0, 0, 2, 7, 0, 0, 7, 7);

	L_CoordsInv3D_jacob_pideCoords3D(fijo,J3x6t,0,0);
	L_StaticMatrix_OP_mult(J2x3tJ3x3t,J2x3t,J3x3s);
	L_StaticMatrix_OP_mult_sub(J2x6, i1, j1, J2x3tJ3x3t,J3x6t,0,0,2,3,0,0,3,6);
}


void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &prJ3x3s, const L_StaticMatrix<3,3> &pc2J3x3s, const L_StaticMatrix<3,3> &prpc2J3x3s, const double *prpre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1) const
{
	//f(pCM) = proy(fijo_mov(pC2_C,fijo_mov(pCM,fijo)))
	//df/dpr = dproy/dp * dfm(pC2_C)/dder * dfm(pCM)/dizq

	//f(pCM) = proy(fijo_mov(pCM*pC2_C,fijo))
	//df/fijo = dproy/dp * dfm(pCM*pC2_C)/dder

	L_StaticMatrix<2,3> J2x3t; //dproy/dp
	L_StaticMatrix<3,7> J3x7t; // dfm(pCM)/dizq
	L_StaticMatrix<3,7> J3x3J3x7t; // dfm(pC2_C)/dder * dfm(pCM)/dizq
	L_StaticMatrix<3,6> J3x6; // dc3d/dfijo
	L_CoordsCart3D movido;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x3.li < i0+2 || J2x3.lj < j0+3, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2() : matriz de tamano insuficiente");

	movido = pCM.fijo_a_movido(pC2_C.fijo_a_movido(fijo));
	L_RayoCamara_jacob_proyeccionDe(movido,J2x3t,0,0);
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(pCM,fijo,prpre36m,prJ3x3s,J3x7t,0,0);

	L_StaticMatrix_OP_mult(J3x3J3x7t, pc2J3x3s, J3x7t);
	L_StaticMatrix_OP_mult_sub(J2x7, i0, j0, J2x3t, J3x3J3x7t, 0,0,2,3, 0,0,3,7);
	L_StaticMatrix_OP_mult_sub(J2x3, i1, j1, J2x3t, prpc2J3x3s, 0,0,2,3, 0,0,3,3);
}

void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, const L_StaticMatrix<3,3> &prJ3x3s, const L_StaticMatrix<3,3> &pc2J3x3s, const L_StaticMatrix<3,3> &prpc2J3x3s, const double *prpre36m, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1) const
{
	//f(pCM) = proy(fijo_mov(pC2_C,fijo_mov(pCM,c3d(fijo))))
	//df/dpr = dproy/dp * dfm(pC2_C)/dder * dfm(pCM)/dizq

	//f(pCM) = proy(fijo_mov(pCM*pC2_C,c3d(fijo)))
	//df/fijo = dproy/dp * dfm(pCM*pC2_C)/dder * dc3d/dfijo

	L_StaticMatrix<2,3> J2x3t; //dproy/dp
	L_StaticMatrix<3,7> J3x7t; // dfm(pCM)/dizq
	L_StaticMatrix<3,7> J3x3J3x7t; // dfm(pC2_C)/dder * dfm(pCM)/dizq
	L_StaticMatrix<3,6> J3x6t; // dc3d/dfijo
	L_StaticMatrix<3,6> J2x3; // dproy/dp * dfm(pCM*pC2_C)/dder
	L_CoordsCart3D movido, pfijo;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x6.li < i0+2 || J2x6.lj < j0+6, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2() : matriz de tamano insuficiente");

	pfijo = fijo.pideCoords3D();
	movido = pCM.fijo_a_movido(pC2_C.fijo_a_movido(pfijo));

	L_RayoCamara_jacob_proyeccionDe(movido,J2x3t,0,0);
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(pCM,pfijo,prpre36m,prJ3x3s,J3x7t,0,0); // dfm(pCM)/dizq
	L_CoordsInv3D_jacob_pideCoords3D(fijo,J3x6t,0,0);

	L_StaticMatrix_OP_mult(J3x3J3x7t, pc2J3x3s, J3x7t);
	L_StaticMatrix_OP_mult_sub(J2x7, i0, j0, J2x3t, J3x3J3x7t, 0,0,2,3, 0,0,3,7);
	L_StaticMatrix_OP_mult(J2x3, J2x3t, prpc2J3x3s);
	L_StaticMatrix_OP_mult_sub(J2x6, i1, j1, J2x3, J3x6t, 0,0,2,3, 0,0,3,6);
}

void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2_(const L_CoordsCart3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x3, int i1, int j1)
{
	int j;
	L_CoordsCart3D pun = fijo;
	L_Pose3D_cuat p1 = pCM;
	L_Pose3D_cuat p2 = pC2_C;
	L_RayoCamara ray1, ray2;
	double delta = 1.0e-3;
	L_StaticMatrix<3,3> J3x3s;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x3.li < i0+2 || J2x3.lj < j0+3, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_v_c2_() : matriz de tamano insuficiente");

	for (j=0; j<7; j++)
	{
		p1.el(j) += delta;
		ray1 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		p1.el(j) -= 2*delta;
		ray2 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		p1.el(j) += delta;
		J2x7(0+i0,j+j0) = (ray1.tanIzq - ray2.tanIzq) / (2*delta);
		J2x7(1+i0,j+j0) = (ray1.tanArr - ray2.tanArr) / (2*delta);
	}
	for (j=0; j<3; j++)
	{
		pun.el(j) += delta;
		ray1 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		pun.el(j) -= 2*delta;
		ray2 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		pun.el(j) += delta;
		J2x3(0+i1,j+j1) = (ray1.tanIzq - ray2.tanIzq) / (2*delta);
		J2x3(1+i1,j+j1) = (ray1.tanArr - ray2.tanArr) / (2*delta);
	}
}


void L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2_(const L_CoordsInv3D &fijo, const L_Pose3D_cuat &pCM, const L_Pose3D_cuat &pC2_C, L_Matrix &J2x7, int i0, int j0, L_Matrix &J2x6, int i1, int j1)
{
	int j;
	L_CoordsInv3D puni = fijo;
	L_CoordsCart3D pun;
	L_Pose3D_cuat p1 = pCM;
	L_Pose3D_cuat p2 = pC2_C;
	L_RayoCamara ray1, ray2;
	double delta = 1.0e-3;
	L_StaticMatrix<3,3> J3x3s;

	throw_L_ArgException_if(J2x7.li < i0+2 || J2x7.lj < j0+7 || J2x6.li < i0+2 || J2x6.lj < j0+6, "L_Pose3D_cuat::jacob_proyectar_fijo_a_movido_vi_c2_() : matriz de tamano insuficiente");

	pun = puni.pideCoords3D();
	for (j=0; j<7; j++)
	{
		p1.el(j) += delta;
		ray1 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		p1.el(j) -= 2*delta;
		ray2 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(pun));
		p1.el(j) += delta;
		J2x7(0+i0,j+j0) = (ray1.tanIzq - ray2.tanIzq) / (2*delta);
		J2x7(1+i0,j+j0) = (ray1.tanArr - ray2.tanArr) / (2*delta);
	}
	for (j=0; j<6; j++)
	{
		puni.el(j) += delta;
		ray1 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(puni.pideCoords3D()));
		puni.el(j) -= 2*delta;
		ray2 = p2.proyectar_fijo_a_movido_v(p1.fijo_a_movido(puni.pideCoords3D()));
		puni.el(j) += delta;
		J2x6(0+i1,j+j1) = (ray1.tanIzq - ray2.tanIzq) / (2*delta);
		J2x6(1+i1,j+j1) = (ray1.tanArr - ray2.tanArr) / (2*delta);
	}
}


bool L_Pose3D_cuat::calcpose_movido_a_fijo_3puntos(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo)
{
	L_HomogeneousMatrix A, B, Ai, C;
	L_CoordsCart3D p1, p2, p3, p4, p5, p6;

	p1 = movido[1]-movido[0];
	if (p1.normaliza_ret() == false)
		return false;
	p2 = movido[2]-movido[0];
	p2 = p2 - p1*(p1.punto(p2));
	if (p2.normaliza_ret() == false)
		return false;
	p3 = p1.cruz(p2);

	A(0,0) = p1.x; A(0,1) = p2.x; A(0,2) = p3.x; A(0,3) = movido[0].x;
	A(1,0) = p1.y; A(1,1) = p2.y; A(1,2) = p3.y; A(1,3) = movido[0].y;
	A(2,0) = p1.z; A(2,1) = p2.z; A(2,2) = p3.z; A(2,3) = movido[0].z;

	p4 = fijo[1]-fijo[0];
	if (p4.normaliza_ret() == false)
		return false;
	p5 = fijo[2]-fijo[1];
	p5 = p5 - p4*(p4.punto(p5));
	if (p5.normaliza_ret() == false)
		return false;
	p6 = p4.cruz(p5);

	B(0,0) = p4.x; B(0,1) = p5.x; B(0,2) = p6.x; B(0,3) = fijo[0].x;
	B(1,0) = p4.y; B(1,1) = p5.y; B(1,2) = p6.y; B(1,3) = fijo[0].y;
	B(2,0) = p4.z; B(2,1) = p5.z; B(2,2) = p6.z; B(2,3) = fijo[0].z;

	// A -> rayos, B -> puntos
	// A = (m1-m0 m2-m0 m3-m0 m0),  B=(f1-f0 f2-f0 f3-f0 f0)
	// A = (m1-m0 m2-m0 m3-m0 m0),  B=(Hm1-Hm0 Hm2-Hm0 Hm3-Hm0 Hm0)
	// B = B
	// (Hm1-Hm0 Hm2-Hm0 Hm3-Hm0 Hm0) = (f1-f0 f2-f0 f3-f0 f0)
	// H(m1-m0 m2-m0 m3-m0 m0)=(f1-f0 f2-f0 f3-f0 f0)
	// HA = B
	// H = B*A^-1
	// H = M(mov) * M(fij)^-1

	Ai.inversaTrasponiendoDe(A); // A*sisOrigen = sisMovido,  Ai*sisMovido = sisOrigen
	C.OP_mult(B, Ai); // B*sisOrigen = sisFijo;    B*AI*sisMovido = sisFijo

	*this = C.calcPose3D_cuat();
	return true;
}


// Funciones que transforman rayos movidos en puntos fijos
bool L_Pose3D_cuat::calcPose_movido_a_fijo_4rayos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, int i0)
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

	// Puntos movidos conocidos
	cm[0] = pm[i0+0];
	cm[1] = pm[i0+1];
	cm[2] = pm[i0+2];
	cm[3] = pm[i0+3];

	// Lados del triangulo, usar "cm" que es el que existe
	cpt.a = (cm[1]-cm[2]).pideR();
	cpt.b = (cm[0]-cm[2]).pideR();
	cpt.c = (cm[0]-cm[1]).pideR();

	// Rayos que apuntan hacia los vertices del triangulo
	cpt.q1 = rf[0+i0];
	cpt.q2 = rf[1+i0];
	cpt.q3 = rf[2+i0];
	
	if (cpt.metodoDirecto() == false) // Usando los primeros 3 puntos da 4 posibles soluciones
		return false;

	// Usar el cuarto punto para seleccionar la mejor solucion de las 4 posibles
	for (i=0; i<4; i++)
	{
		if (cpt.valida[i] == false)
		{
			e[i] = 1e100;
			continue;
		}
		// Calcular el triangulo fijo numero i
		cf[0] = cpt.j1*cpt.s1[i];
		cf[1] = cpt.j2*cpt.s2[i];
		cf[2] = cpt.j3*cpt.s3[i];
		// calcular la pose que lleva el triangulo movido (conocido) al fijo (calculado usando alg3p)
		t[i].calcpose_movido_a_fijo_3puntos(cm, cf); // cf = t*rm  -> t = H(cm)*H(cf)^-1
		// calcular el cuarto punto fijo
		cf[3] = t[i].movido_a_fijo(cm[3]);
		// errores
		ray.tanIzq = cf[3].y / cf[3].x;
		ray.tanArr = cf[3].z / cf[3].x;
		e[i] = (ray.tanIzq - rf[3+i0].tanIzq)*(ray.tanIzq - rf[3+i0].tanIzq) +
			(ray.tanArr - rf[3+i0].tanArr)*(ray.tanArr - rf[3+i0].tanArr);
	}
	iMin = 0;
	for (i=1; i<4; i++)
		if (e[i] < e[iMin])
			iMin = i;
	*this = t[iMin];

	if (pos.x > 1.0e5 || pos.x < -1.0e5 || pos.y > 1.0e5 || pos.y < -1.0e5 || pos.z > 1.0e5 || pos.z < -1.0e5 || ori.abs2() < 0.9 || ori.abs2() > 1.1)
		return false;

	return true;
}
void L_Pose3D_cuat::calcPose_movido_a_fijo_rapidoMuyImpreciso(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo)
{
	std::vector<L_CoordsCart3D> f, m;
	L_CoordsCart3D pf, pm;
	pf = L_VMC(fijo).mean();
	pm = L_VMC(movido).mean();
	pos = pf - pm;
	L_VM(f).subtractElement(fijo, pf);
	L_VM(m).subtractElement(movido, pm);
	ori.fijaNpuntos_promediaw(fijo, movido);
}

void L_Pose3D_cuat::calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo)
{
	L_Matrix A, b, AT, ATA, ATAinv, ATb, X;
	L_HomogeneousMatrix H;
	int i;
	A.reallocate(3*(int)movido.size(), 12);
	b.reallocate(3*(int)movido.size(), 1);
	A.setZero();

	for (i=0; i<(int)movido.size(); i++)
	{
		A(3*i+0,0)=movido[i].x;
		A(3*i+0,1)=movido[i].y;
		A(3*i+0,2)=movido[i].z;
		A(3*i+0,3)=1;

		A(3*i+1,4)=movido[i].x;
		A(3*i+1,5)=movido[i].y;
		A(3*i+1,6)=movido[i].z;
		A(3*i+1,7)=1;

		A(3*i+2,8)=movido[i].x;
		A(3*i+2,9)=movido[i].y;
		A(3*i+2,10)=movido[i].z;
		A(3*i+2,11)=1;

		b(3*i+0,0) = fijo[i].x;
		b(3*i+1,0) = fijo[i].y;
		b(3*i+2,0) = fijo[i].z;
	}
	AT.transpOf(A);
	ATA.OP_mult(AT, A);
	ATAinv.inverseOf(ATA);
	ATb.OP_mult(AT, b);
	X.OP_mult(ATAinv, ATb);
	for (i=0; i<12; i++)
		H(i/4,i%4) = X(i,0);
	H.normalizaSVD();
	*this = H.calcPose3D_cuat();
}

void L_Pose3D_cuat::calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Pose3D_cuat *inic, L_LevenbergMarquardt *opt)
{
	L_Matrix pose;
	L_LevenbergMarquardt optim;
	L_Pose3D_cuat_nube_puntos nube;

	nube.fijo = &fijo;
	nube.movido = &movido;
	if (inic == NULL)
		calcPose_movido_a_fijo_matriz(movido, fijo); // Estimacion inicial
	else
		*this = *inic;

	// Parametros del optimizador
	optim.epsilon = 1e-7;
	optim.nIterationsMax = niter;
	optim.lenVector = 7;
	
	if (opt == NULL)
		opt = &optim;

	//opt.useHessian = true; // Solo se usa en la minimizacion con error escalar
	pose.reallocate(7,1);
	pose(0,0) = pos.x;
	pose(1,0) = pos.y;
	pose(2,0) = pos.z;
	pose(3,0) = ori.a;
	pose(4,0) = ori.b;
	pose(5,0) = ori.c;
	pose(6,0) = ori.d;
	opt->minimize_vect((const void *)&nube, pose, &err_mf_vect);
	pos.x = pose(0,0);
	pos.y = pose(1,0);
	pos.z = pose(2,0);
	ori.a = pose(3,0);
	ori.b = pose(4,0);
	ori.c = pose(5,0);
	ori.d = pose(6,0);
	ori.normaliza(); // No tiene por que quedar normalizado, no es una restriccion de la optimizacion, solo que no afecta el resultado
}

void L_Pose3D_cuat::calcPose_fijo_a_movido_matriz_v2(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo)
{
	L_Matrix A, b, AT, ATA, ATAinv, ATb, X;
	L_HomogeneousMatrix H;
	int i;
	A.reallocate(3*(int)movido.size(), 12);
	b.reallocate(3*(int)movido.size(), 1);
	A.setZero();

	for (i=0; i<(int)movido.size(); i++)
	{
		A(3*i+0,0)=A(3*i+1,4)=A(3*i+2,8)=fijo[i].x;
		A(3*i+0,1)=A(3*i+1,5)=A(3*i+2,9)=fijo[i].y;
		A(3*i+0,2)=A(3*i+1,6)=A(3*i+2,10)=fijo[i].z;
		A(3*i+0,3)=A(3*i+1,7)=A(3*i+2,11)=1;
		b(3*i+0,0) = movido[i].x;
		b(3*i+1,0) = movido[i].y;
		b(3*i+2,0) = movido[i].z;
	}
	AT.transpOf(A);
	ATA.OP_mult(AT, A);
	ATAinv.inverseOf(ATA);
	ATb.OP_mult(AT, b);
	X.OP_mult(ATAinv, ATb);
	for (i=0; i<12; i++)
		H(i/4,i%4) = X(i,0);
	H.normalizaSVD();
	*this = H.calcPose3D_cuat();
}

void L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0, int j0)
{
	L_Matrix A, b, AT, ATA, ATAinv, ATb, X, J7x12, J12x3n;
	L_HomogeneousMatrix H;
	int i;
	double delta = 1.0e-8;
	A.reallocate(3*(int)movido.size(), 12);
	b.reallocate(3*(int)movido.size(), 1);
	A.setZero();

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz() : matriz de tamano insuficiente");

	for (i=0; i<(int)movido.size(); i++)
	{
		A(3*i+0,0)=A(3*i+1,4)=A(3*i+2,8)=movido[i].x;
		A(3*i+0,1)=A(3*i+1,5)=A(3*i+2,9)=movido[i].y;
		A(3*i+0,2)=A(3*i+1,6)=A(3*i+2,10)=movido[i].z;
		A(3*i+0,3)=A(3*i+1,7)=A(3*i+2,11)=1;
		b(3*i+0,0) = fijo[i].x;
		b(3*i+1,0) = fijo[i].y;
		b(3*i+2,0) = fijo[i].z;
	}
	AT.transpOf(A); // A:(3*n)x12, AT:12x(3*n)
	ATA.OP_mult(AT, A); // ATA: 12x12
	if (ATAinv.inverseOf(ATA) == false) // 12x12
	{
		printf("L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz() : Matriz ATA no invertible\n");
	}
	ATb.OP_mult(AT, b); // 12x1
	X.OP_mult(ATAinv, ATb); // 12x1
	for (i=0; i<12; i++)
		H(i/4,i%4) = X(i,0);

	//H.normalizaSVD();
	J12x3n.OP_mult(ATAinv,AT); //12x3n
	H.normalizaSVD_corrigeJ(J12x3n);

	J7x12.reallocate(7,12); // Esta chiquita la calculo por diferencias finitas para evitar problemas de signo del cuaternion
	*this = H.jacob_calcPose3D_cuat(J7x12);

	J7x12.checkRangeValues(-1e10, 1e10, true, "J7x12");
	J12x3n.checkRangeValues(-1e10, 1e10, true, "J7x3n");
	J7x3n.OP_mult_sub(i0, j0, J7x12, J12x3n, 0, 0, J7x12.li, J7x12.lj, 0, 0, J12x3n.li, J12x3n.lj);
}

void L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz_diffin(const std::vector<L_CoordsCart3D> &movido, std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0, int j0)
{
	// Funka decente
	L_Pose3D_cuat p1, p2;
	int i, j;
	double delta = 1.0e-5;

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz_diffin() : matriz de tamano insuficiente");

	calcPose_movido_a_fijo_matriz(movido, fijo);
	for (j=0; j<3*(int)movido.size(); j++)
	{
		fijo[j/3].el(j%3) += delta;
		p1.calcPose_movido_a_fijo_matriz(movido, fijo);
		if (p1.ori.OP_mult_elementwise(ori) < 0)
			p1.ori.invSigno();
		fijo[j/3].el(j%3) -= 2*delta;
		p2.calcPose_movido_a_fijo_matriz(movido, fijo);
		if (p2.ori.OP_mult_elementwise(ori) < 0)
			p2.ori.invSigno();
		fijo[j/3].el(j%3) += delta;
		for (i=0; i<7; i++)
			J7x3n(i+i0,j+j0) = (p1.el(i) - p2.el(i)) / (2*delta);  // Error: faltaba i0, j0
	}
}


void L_Pose3D_cuat::jacob_movido_calcPose_movido_a_fijo_matriz_diffin(std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0, int j0)
{
	L_Pose3D_cuat p1, p2;
	int i, j;
	double delta = 1.0e-5;

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_matriz_diffin() : matriz de tamano insuficiente");

	calcPose_movido_a_fijo_matriz(movido, fijo);
	for (j=0; j<3*(int)movido.size(); j++)
	{
		movido[j/3].el(j%3) += delta;
		p1.calcPose_movido_a_fijo_matriz(movido, fijo);
		if (p1.ori.OP_mult_elementwise(ori) < 0)
			p1.ori.invSigno();
		movido[j/3].el(j%3) -= 2*delta;
		p2.calcPose_movido_a_fijo_matriz(movido, fijo);
		if (p2.ori.OP_mult_elementwise(ori) < 0)
			p2.ori.invSigno();
		movido[j/3].el(j%3) += delta;
		for (i=0; i<7; i++)
			J7x3n(i+i0,j+j0) = (p1.el(i) - p2.el(i)) / (2*delta);  // Error: faltaba i0, j0
	}
}

void L_Pose3D_cuat::jacob_movido_calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0, int j0)
{
	// cmf(m=M, f=F) = inv(q=cmf(m=F, f=M))
	//	dcmf / dM = dinv/dq|cmf(m=F, f=M) * dcmf/df|(m=F,f=M)
	L_Matrix J, Ji;
	L_Pose3D_cuat p;

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_movido_calcPose_movido_a_fijo_matriz() : matriz de tamano insuficiente");

	J.reallocate(7, 3*(int)movido.size());
	p.jacob_fijo_calcPose_movido_a_fijo_matriz(fijo, movido, J);
	Ji.reallocate(7,7);
	jacob_inverseOf(p, Ji);
	J7x3n.OP_mult_sub(i0, j0, Ji, J, 0,0, 7,7, 0,0, 7,3*(int)movido.size());
	// cmf(m=M, f=F) = cfm(m=F, f=M)
	// dcmf / dM = dcfm/df|(m=F, f=M) ... ...
}


void L_Pose3D_cuat::jacob_movido_calcPose_fijo_a_movido_matriz_v2(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, L_Matrix &J7x3n, int i0, int j0)
{
	L_Matrix A, b, AT, ATA, ATAinv, ATb, X, J7x12, J12x3n;
	L_HomogeneousMatrix H;
	int i;
	double delta = 1.0e-8;
	A.reallocate(3*(int)movido.size(), 12);
	b.reallocate(3*(int)movido.size(), 1);
	A.setZero();

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_movido_calcPose_fijo_a_movido_matriz_v2() : matriz de tamano insuficiente");

	for (i=0; i<(int)movido.size(); i++)
	{
		A(3*i+0,0)=A(3*i+1,4)=A(3*i+2,8)=fijo[i].x;
		A(3*i+0,1)=A(3*i+1,5)=A(3*i+2,9)=fijo[i].y;
		A(3*i+0,2)=A(3*i+1,6)=A(3*i+2,10)=fijo[i].z;
		A(3*i+0,3)=A(3*i+1,7)=A(3*i+2,11)=1;
		b(3*i+0,0) = movido[i].x;
		b(3*i+1,0) = movido[i].y;
		b(3*i+2,0) = movido[i].z;
	}
	AT.transpOf(A); // A:(3*n)x12, AT:12x(3*n)
	ATA.OP_mult(AT, A); // ATA: 12x12
	ATAinv.inverseOf(ATA); // 12x12
	ATb.OP_mult(AT, b); // 12x1
	X.OP_mult(ATAinv, ATb); // 12x1
	for (i=0; i<12; i++)
		H(i/4,i%4) = X(i,0);
	H.normalizaSVD();

	J7x12.reallocate(7,12); // Esta chiquita la calculo por diferencias finitas para evitar problemas de signo del cuaternion
	*this = H.jacob_calcPose3D_cuat(J7x12);
	J12x3n.OP_mult(ATAinv,AT); //12x3n

	J7x12.checkRangeValues(-1e10, 1e10, true, "J7x12");
	J12x3n.checkRangeValues(-1e10, 1e10, true, "J7x3n");
	J7x3n.OP_mult_sub(i0, j0, J7x12, J12x3n, 0, 0, J7x12.li, J7x12.lj, 0, 0, J12x3n.li, J12x3n.lj);
}

void L_Pose3D_cuat::jacob_movido_calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Matrix &J7x3n, int i0, int j0, double delta)
{
	std::vector<L_CoordsCart3D> movido2;
	L_Pose3D_cuat ppp;
	int i, j;

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_iter() : matriz de tamano insuficiente");

	L_LevenbergMarquardt optim;

	// Parametros del optimizador
	optim.epsilon = 1e-7;
	optim.nIterationsMax = niter;
	optim.lenVector = 7;

	movido2 = movido;

	// Estimacion inicial para las iteraciones
	calcPose_movido_a_fijo_iter(movido, fijo, niter, NULL, &optim); // Para hacer converger el optimizador
	ppp = *this;

	// Parametros del optimizador
	for (j=0; j<3*(int)movido2.size(); j++)
	{
		movido2[j/3].el(j%3) += delta;
		ppp.calcPose_movido_a_fijo_iter(movido2, fijo, 3, this, &optim); // Ya no lleva niter
		if (ppp.ori.OP_mult_elementwise(ori) < 0)
			ppp.ori.invSigno();
		movido2[j/3].el(j%3) -= delta;

		for (i=0; i<7; i++)
			J7x3n(i+i0,j+j0) = (ppp.el(i)-el(i)) / delta;
	}
}

void L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_iter(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo, int niter, L_Matrix &J7x3n, int i0, int j0, double delta)
{
	std::vector<L_CoordsCart3D> fijo2;
	L_Pose3D_cuat ppp;
	int i, j;

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)movido.size(), "L_Pose3D_cuat::jacob_fijo_calcPose_movido_a_fijo_iter() : matriz de tamano insuficiente");

	L_LevenbergMarquardt optim;

	// Parametros del optimizador
	optim.epsilon = 1e-7;
	optim.nIterationsMax = niter;
	optim.lenVector = 7;

	fijo2 = fijo;

	// Estimacion inicial para las iteraciones
	calcPose_movido_a_fijo_iter(movido, fijo, niter, NULL, &optim); // Para hacer converger el optimizador
	ppp = *this;

	// Parametros del optimizador
	for (j=0; j<3*(int)fijo2.size(); j++)
	{
		fijo2[j/3].el(j%3) += delta;
		ppp.calcPose_movido_a_fijo_iter(movido, fijo2, 3, this, &optim); // Ya no lleva niter
		if (ppp.ori.OP_mult_elementwise(ori) < 0)
			ppp.ori.invSigno();
		fijo2[j/3].el(j%3) -= delta;

		for (i=0; i<7; i++)
			J7x3n(i+i0,j+j0) = (ppp.el(i)-el(i)) / delta;
	}
}


bool L_Pose3D_cuat::_calcPose_movido_a_fijo_4rayos3D(const L_Array<L_Rayo3D> &rm, const std::vector<L_CoordsCart3D> &pf, int i0)
{
	// Hay que implementar "A Minimal Solution to the Generalised 3-Point Pose Problem", David Nister
	// Da 8 soluciones posibles con 3 puntos, el cuarto se puede usar para seleccionar la mejor de ellas
	// No tengo tiempo para implementar esto ahora
	throw L_NonImplementedException();
}

void L_Pose3D_cuat::_calcPose_movido_a_fijo_Nrayos(const std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, int niter, L_Pose3D_cuat *inicial)
{
	L_LevenbergMarquardt opt;
	L_Pose3D_cuat_rayos_rmpf nube;
	L_Pose3D_cuat inic, ptmp;
	L_Matrix vect;
	int i;

	// Optimizacion iterativa
	nube.movido = &rm;
	nube.fijo = &pf;
	nube.pCam = NULL;

	opt.epsilon = 1.0e-8; // epsilon para las variaciones
	opt.errorToStop = 1.0e-6; // si la variacion del error entre iteraciones es menor a esto, se para
	opt.nIterationsMax = niter;
	opt.lenVector = 7; // pose-cuaternion


	// Estimacion inicial del punto-cuaternion

	if (inicial == NULL)
	{
		L_Matrix errv, p(7,1);
		double err, errMejor;
		int numIn;
		// Probar si la pose (10,0,0) es un buen punto inicial
		L_Pose3D_cuat mejor(L_CoordsCart3D(10,0,0),L_Quaternion(1,0,0,0));
		mejor.puntocuat_a_matriz(p);
		err_proy_fm_vect(&nube, p, errv);
		errMejor = errv.normL2();
		
		// Probar si la pose *this es un buen punto inicial
		ptmp = *this;
		ptmp.puntocuat_a_matriz(p);
		err_proy_fm_vect(&nube, p, errv);
		err = errv.normL2();
		if (err < errMejor)
		{
			errMejor = err;
			mejor = ptmp;
		}

		numIn = L_MIN((int)rm.size()-3, 10);
		for (i=0; i<numIn; i++)
		{
			if (ptmp._calcPose_movido_a_fijo_4rayos(rm, pf, i) == true)
			{
				ptmp.puntocuat_a_matriz(p);
				err_proy_fm_vect(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}
		inic = mejor;
	}
	else
		inic = *inicial;


	vect.reallocate(7, 1);
	inic.puntocuat_a_matriz(vect);

	opt.minimize_vect(&nube, vect, &err_proy_fm_vect);

	matriz_a_puntocuat(vect);

	normaliza();

	//printf("dist opt: %f\n", pos.distanciaA(inic.pos));
}

void L_Pose3D_cuat::calcPose_movido_a_fijo_Nrayos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial, L_Array<L_Pose3D_cuat> *candidatos, int nIntentos3pt)
{
	L_Pose3D_cuat_rayos_pmrf nube;
	L_Pose3D_cuat inic, ptmp;
	L_Matrix vect;
	void (*err_mf)(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores);
	int i;

	// Optimizacion iterativa
	nube.movido = &pm;
	nube.fijo = &rf;
	nube.pCam21 = NULL;

	//err_mf = &err_proy_mf_vect;
	err_mf = &err_proy_mf_vect_rapido;

	// Se le puede pasar el optimizador ya creado para aprovechar el punto de convergencia anterior
	opt.lenVector = 7; // pose-cuaternion

	// Punto inicial para la optimizacion
	if (inicial == NULL)
	{
		L_Matrix errv, p(7,1);
		double err, errMejor;
		int numIn;

		// Probar si la pose (10,0,0) es un buen punto inicial
		L_Pose3D_cuat mejor(L_CoordsCart3D(10,0,0),L_Quaternion(1,0,0,0));
		mejor.puntocuat_a_matriz(p);
		(*err_mf)(&nube, p, errv);
		errMejor = errv.normL2();

		// Probar si la pose *this es un buen punto inicial
		if (ori.a*ori.a + ori.b*ori.b + ori.c*ori.c + ori.d*ori.d > 0) // Podria ser zero
		{
			ptmp = *this;
			ptmp.puntocuat_a_matriz(p);
			(*err_mf)(&nube, p, errv);
			err = errv.normL2();
			if (err < errMejor)
			{
				errMejor = err;
				mejor = ptmp;
			}
		}

		// Probar si los candidatos[] son buenos puntos iniciales
		if (candidatos != NULL)
		{
			for (i=0; i<(int)candidatos->size(); i++)
			{
				ptmp = (*candidatos)[i];
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}

		// Probar si el algoritmo de los 3 puntos da un buen punto inicial
		numIn = L_MIN((int)rf.size()-3, nIntentos3pt);
		for (i=0; i<numIn; i++)
		{
			if (ptmp.calcPose_movido_a_fijo_4rayos(pm, rf, i) == true)
			{
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}
		inic = mejor;
	}
	else
		inic = *inicial;


	vect.reallocate(7, 1);
	inic.puntocuat_a_matriz(vect);

	opt.minimize_vect(&nube, vect, err_mf);

	matriz_a_puntocuat(vect);

	normaliza();
}

void L_Pose3D_cuat::calcPose_movido_a_fijo_Nrayos_sqrt(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial, L_Array<L_Pose3D_cuat> *candidatos, int nIntentos3pt)
{
	L_Pose3D_cuat_rayos_pmrf nube;
	L_Pose3D_cuat inic, ptmp;
	L_Matrix vect;
	void (*err_mf)(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores);
	int i;

	// Optimizacion iterativa
	nube.movido = &pm;
	nube.fijo = &rf;
	nube.pCam21 = NULL;

	//err_mf = &err_proy_mf_vect;
	err_mf = &err_proy_mf_vect_rapido_sqrt;

	// Se le puede pasar el optimizador ya creado para aprovechar el punto de convergencia anterior
	opt.lenVector = 7; // pose-cuaternion

	// Punto inicial para la optimizacion
	if (inicial == NULL)
	{
		L_Matrix errv, p(7,1);
		double err, errMejor;
		int numIn;

		// Probar si la pose (10,0,0) es un buen punto inicial
		L_Pose3D_cuat mejor(L_CoordsCart3D(10,0,0),L_Quaternion(1,0,0,0));
		mejor.puntocuat_a_matriz(p);
		(*err_mf)(&nube, p, errv);
		errMejor = errv.normL2();

		// Probar si la pose *this es un buen punto inicial
		if (ori.a*ori.a + ori.b*ori.b + ori.c*ori.c + ori.d*ori.d > 0) // Podria ser zero
		{
			ptmp = *this;
			ptmp.puntocuat_a_matriz(p);
			(*err_mf)(&nube, p, errv);
			err = errv.normL2();
			if (err < errMejor)
			{
				errMejor = err;
				mejor = ptmp;
			}
		}

		// Probar si los candidatos[] son buenos puntos iniciales
		if (candidatos != NULL)
		{
			for (i=0; i<(int)candidatos->size(); i++)
			{
				ptmp = (*candidatos)[i];
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}

		// Probar si el algoritmo de los 3 puntos da un buen punto inicial
		numIn = L_MIN((int)rf.size()-3, nIntentos3pt);
		for (i=0; i<numIn; i++)
		{
			if (ptmp.calcPose_movido_a_fijo_4rayos(pm, rf, i) == true)
			{
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}
		inic = mejor;
	}
	else
		inic = *inicial;


	vect.reallocate(7, 1);
	inic.puntocuat_a_matriz(vect);

	opt.minimize_vect(&nube, vect, err_mf);

	matriz_a_puntocuat(vect);

	normaliza();
}


void L_Pose3D_cuat::calcPose_movido_a_fijo_Npuntos(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_CoordsCart3D> &pf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial, L_Array<L_Pose3D_cuat> *candidatos, int nIntentos3pt)
{
	L_Pose3D_cuat_nube_puntos nube;
	L_Pose3D_cuat inic, ptmp;
	L_Matrix vect;
	void (*err_mf)(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores);
	int i;

	// Optimizacion iterativa
	nube.movido = &pm;
	nube.fijo = &pf;

	//err_mf = &err_proy_mf_vect;
	err_mf = &err_mf_vect;

	// Se le puede pasar el optimizador ya creado para aprovechar el punto de convergencia anterior
	opt.lenVector = 7; // pose-cuaternion

	// Punto inicial para la optimizacion
	if (inicial == NULL)
	{
		L_Matrix errv, p(7,1);
		double err, errMejor;

		// Probar si la pose (10,0,0) es un buen punto inicial
		L_Pose3D_cuat mejor(L_CoordsCart3D(10,0,0),L_Quaternion(1,0,0,0));
		mejor.puntocuat_a_matriz(p);
		(*err_mf)(&nube, p, errv);
		errMejor = errv.normL2();

		// Probar si la pose *this es un buen punto inicial
		if (ori.a*ori.a + ori.b*ori.b + ori.c*ori.c + ori.d*ori.d > 0) // Podria ser zero
		{
			ptmp = *this;
			ptmp.puntocuat_a_matriz(p);
			(*err_mf)(&nube, p, errv);
			err = errv.normL2();
			if (err < errMejor)
			{
				errMejor = err;
				mejor = ptmp;
			}
		}

		// Probar si los candidatos[] son buenos puntos iniciales
		if (candidatos != NULL)
		{
			for (i=0; i<(int)candidatos->size(); i++)
			{
				ptmp = (*candidatos)[i];
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}
		inic = mejor;
	}
	else
		inic = *inicial;


	vect.reallocate(7, 1);
	inic.puntocuat_a_matriz(vect);

	opt.minimize_vect(&nube, vect, err_mf);

	matriz_a_puntocuat(vect);

	normaliza();
}

void L_Pose3D_cuat::calcPose_movido_a_fijo_Npuntos_sqrt(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_CoordsCart3D> &pf, L_LevenbergMarquardt &opt, L_Pose3D_cuat *inicial, L_Array<L_Pose3D_cuat> *candidatos, int nIntentos3pt)
{
	L_Pose3D_cuat_nube_puntos nube;
	L_Pose3D_cuat inic, ptmp;
	L_Matrix vect;
	void (*err_mf)(const void *pose3D_cuat_nube_ray, const L_Matrix &pcuat, L_Matrix &errores);
	int i;

	// Optimizacion iterativa
	nube.movido = &pm;
	nube.fijo = &pf;

	//err_mf = &err_proy_mf_vect;
	err_mf = &err_mf_vect;

	// Se le puede pasar el optimizador ya creado para aprovechar el punto de convergencia anterior
	opt.lenVector = 7; // pose-cuaternion

	// Punto inicial para la optimizacion
	if (inicial == NULL)
	{
		L_Matrix errv, p(7,1);
		double err, errMejor;

		// Probar si la pose (10,0,0) es un buen punto inicial
		L_Pose3D_cuat mejor(L_CoordsCart3D(10,0,0),L_Quaternion(1,0,0,0));
		mejor.puntocuat_a_matriz(p);
		(*err_mf)(&nube, p, errv);
		errMejor = errv.normL2();

		// Probar si la pose *this es un buen punto inicial
		if (ori.a*ori.a + ori.b*ori.b + ori.c*ori.c + ori.d*ori.d > 0) // Podria ser zero
		{
			ptmp = *this;
			ptmp.puntocuat_a_matriz(p);
			(*err_mf)(&nube, p, errv);
			err = errv.normL2();
			if (err < errMejor)
			{
				errMejor = err;
				mejor = ptmp;
			}
		}

		// Probar si los candidatos[] son buenos puntos iniciales
		if (candidatos != NULL)
		{
			for (i=0; i<(int)candidatos->size(); i++)
			{
				ptmp = (*candidatos)[i];
				ptmp.puntocuat_a_matriz(p);
				(*err_mf)(&nube, p, errv);
				err = errv.normL2();
				if (err < errMejor)
				{
					errMejor = err;
					mejor = ptmp;
				}
			}
		}
		inic = mejor;
	}
	else
		inic = *inicial;


	vect.reallocate(7, 1);
	inic.puntocuat_a_matriz(vect);

	opt.minimize_vect(&nube, vect, err_mf);

	matriz_a_puntocuat(vect);

	normaliza();
}



void L_Pose3D_cuat::calcPose_movido_a_fijo_Nrayos3D(const std::vector<L_CoordsCart3D> &pf, const L_Array<L_Rayo3D> &rf, int niter, L_Pose3D_cuat *inicial)
{
	// Hay que implementar "A Minimal Solution to the Generalised 3-Point Pose Problem", David Nister
	// Hay que tomar subconjuntos de 4 puntos y llamar a calcPose_movido_a_fijo_4rayos3D() para el punto inicial
	// No tengo tiempo para implementar esto ahora
	throw L_NonImplementedException();
}

void L_Pose3D_cuat::_calcPose_movido_a_fijo_Nrayos3D(const L_Array<L_Rayo3D> &rm, const std::vector<L_CoordsCart3D> &pf, int niter, L_Pose3D_cuat *inicial)
{
	// Hay que implementar "A Minimal Solution to the Generalised 3-Point Pose Problem", David Nister
	// Hay que tomar subconjuntos de 4 puntos y llamar a calcPose_movido_a_fijo_4rayos3D() para el punto inicial
	// No tengo tiempo para implementar esto ahora
	throw L_NonImplementedException();
}


void L_Pose3D_cuat::_jacob_rm_calcPose_movido_a_fijo_Nrayos_diffin(std::vector<L_RayoCamara> &rm, const std::vector<L_CoordsCart3D> &pf, L_Matrix &J7x2n, int i0, int j0)
{
	L_Pose3D_cuat p1, p2;
	int i, j;
	double delta = 1.0e-2;
	L_Pose3D_cuat zero;
	zero.fijaCero();
	_calcPose_movido_a_fijo_Nrayos(rm, pf, 10); // Para poder evaluar el signo del cuaternion
	//calcPose_movido_a_fijo_Nrayos(rm, pf, 10, &zero); // Para poder evaluar el signo del cuaternion
	for (j=0; j<2*(int)rm.size(); j++)
	{
		rm[j/2].el(j%2) += delta;
		p1._calcPose_movido_a_fijo_Nrayos(rm, pf, 3, this);
		if (p1.ori.OP_mult_elementwise(ori) < 0)
			p1.ori.invSigno();
		rm[j/2].el(j%2) -= 2*delta;
		p2._calcPose_movido_a_fijo_Nrayos(rm, pf, 3, this);
		if (p2.ori.OP_mult_elementwise(ori) < 0)
			p2.ori.invSigno();
		rm[j/2].el(j%2) += delta;
		for (i=0; i<7; i++)
			J7x2n(i+i0,j+j0) = (p1.el(i) - p2.el(i)) / (2*delta);
	}
}




// Funciones que transforman puntos movidos en rayos fijos


void L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_matriz_malo(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x2n, int i0, int j0)
{
	L_Matrix JJ7x2n(7,2*(int)pm.size()), J7x7i(7,7);

	throw_L_ArgException_if(J7x2n.li < i0+7 || J7x2n.lj < j0+2*(int)pm.size(), "L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_matriz_malo() : matriz de tamano insuficiente");

	_jacob_rm_calcPose_movido_a_fijo_Nrayos_matriz_malo(rf, pm, JJ7x2n, 0, 0);
	jacob_inverseOf(*this, J7x7i);
	J7x2n.OP_mult_sub(i0, j0, J7x7i, JJ7x2n, 0, 0, 7, 7, 0, 0, 7, 2*(int)pm.size());
}

void L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x2n, int i0, int j0, L_Pose3D_cuat *inic, L_Array<L_Pose3D_cuat> *candidatos, double delta, double epsilon, double errorToStop, int nIterationsMax, int nIteracParar2)
{
	L_Pose3D_cuat p1, p2;
	int i, j;
	L_LevenbergMarquardt opt;
	std::vector<L_RayoCamara> r;
	r = rf;

	throw_L_ArgException_if(J7x2n.li < i0+7 || J7x2n.lj < j0+2*(int)r.size(), "L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin() : matriz de tamano insuficiente");

	opt.epsilon = epsilon; // 1.0e-8; // epsilon para las variaciones
	opt.errorToStop = errorToStop; // 1.0e-6; // si la variacion del error entre iteraciones es menor a esto, se para
	opt.nIterationsMax = nIterationsMax; // 10;  // 10 iteraciones como maxElement

	if (inic == NULL)
		calcPose_movido_a_fijo_Nrayos(pm, r, opt, NULL, candidatos);
	else
		*this = *inic;
	for (j=0; j<2*(int)r.size(); j++)
	{
		opt.nIterationsMax = nIteracParar2; // 4;
		r[j/2].el(j%2) += delta;
		p1.calcPose_movido_a_fijo_Nrayos(pm, r, opt, this);
		if (p1.ori.OP_mult_elementwise(ori) < 0)
			p1.ori.invSigno();
		r[j/2].el(j%2) -= 2*delta;
		p2.calcPose_movido_a_fijo_Nrayos(pm, r, opt, this);
		if (p2.ori.OP_mult_elementwise(ori) < 0)
			p2.ori.invSigno();
		r[j/2].el(j%2) += delta;
		for (i=0; i<7; i++)
			J7x2n(i+i0,j+j0) = (p1.el(i) - p2.el(i)) / (2*delta);
	}
}


bool L_Pose3D_cuat::jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &Hess7x7inv, double delta)
{
	//   Aca hay que calcular el hessiano del error de proyeccion respecto al puntocuaternion
	// y luego invertirlo
	//   Este calculo no se puede hacer muy rapido usando ningun truco, aunque el hessiano es simetrico y
	// eso si se puede aprovechar
	int i, j, k; // 7x7
	Hess7x7inv.reallocate(7,7);
	L_Pose3D_cuat pose;
	L_CoordsCart3D pun;
	L_HomogeneousMatrix H;
	double dI, dA;
	double err[4];
	// Aca el error que se considera es el correspondiente a err_proy_mf_vect_rapido()

	for (i=0; i<7; i++)
	{
		// Calcular los hessianos cruzados
		for (j=i; j<7; j++)
		{
			err[0] = 0;
			err[1] = 0;
			err[2] = 0;
			err[3] = 0;
			pose = *this; // Para recuperar el original
			pose.el(i) += delta;
			pose.el(j) += delta;
			H.fijaPose3D_cuat(pose); // movido a fijo
			for (k=0; k<(int)pm.size(); k++)
			{
				L_MatrizHomogenea_OP_mult_P(pun, H, pm[k]);
				dI = pun.y/pun.x - rf[k].tanIzq;
				dA = pun.z/pun.x - rf[k].tanArr;
				err[0] += dI*dI + dA*dA; 
			}
			pose = *this; // Para recuperar el original
			pose.el(i) += delta;
			pose.el(j) -= delta;
			H.fijaPose3D_cuat(pose); // movido a fijo
			for (k=0; k<(int)pm.size(); k++)
			{
				L_MatrizHomogenea_OP_mult_P(pun, H, pm[k]);
				dI = pun.y/pun.x - rf[k].tanIzq;
				dA = pun.z/pun.x - rf[k].tanArr;
				err[1] += dI*dI + dA*dA; 
			}
			pose = *this; // Para recuperar el original
			pose.el(i) -= delta;
			pose.el(j) += delta;
			H.fijaPose3D_cuat(pose); // movido a fijo
			for (k=0; k<(int)pm.size(); k++)
			{
				L_MatrizHomogenea_OP_mult_P(pun, H, pm[k]);
				dI = pun.y/pun.x - rf[k].tanIzq;
				dA = pun.z/pun.x - rf[k].tanArr;
				err[2] += dI*dI + dA*dA; 
			}
			pose = *this; // Para recuperar el original
			pose.el(i) -= delta;
			pose.el(j) -= delta;
			H.fijaPose3D_cuat(pose); // movido a fijo
			for (k=0; k<(int)pm.size(); k++)
			{
				L_MatrizHomogenea_OP_mult_P(pun, H, pm[k]);
				dI = pun.y/pun.x - rf[k].tanIzq;
				dA = pun.z/pun.x - rf[k].tanArr;
				err[3] += dI*dI + dA*dA; 
			}
			Hess7x7inv(i,j) = ( (err[0]-err[2])/(2*delta) - (err[1]-err[3])/(2*delta) ) / (2*delta);
			Hess7x7inv(j,i) = Hess7x7inv(i,j);
		}
	}
	// Hess7x7inv es singular ya que hay una direccion que no aporta informacion
	// En consecuencia, ese value propio puede ser transformado de zero a un value arbitrario
	Hess7x7inv.invertMe_invertirValoresPropios(1.0e-5);
	return true;
}

void L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, const L_Matrix &Hess7x7inv, L_Matrix &J7x2n, int i0, int j0, double delta)
{
	throw_L_ArgException_if(J7x2n.li < i0+7 || J7x2n.lj < j0+2*(int)rf.size(), "L_Pose3D_cuat::jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin_analitica() : matriz de tamano insuficiente");
	// Se sabe que la pose actual es la que lleva los puntos
	// pm a los rayos rf en forma optima
	// error = suma_i ( (eta*pm[i])[2]/(eta*pm[i])[1] - rf[i][1] )^2 + ( (eta*pm[i])[3]/(eta*pm[i])[1] - rf[i][2] )^2
	//
	// El error es una sum de terminos, cada elemento de la sum depende de solo algunos datos
	// En consecuencia, hay que evaluar solo la variacion de un termino cada vez y no de la sum completa

	int i, j; // 7 x rf.size()
	L_Pose3D_cuat pose;
	L_CoordsCart3D pun;
	L_HomogeneousMatrix H;
	double dI, dA;
	L_Matrix mFxa;
	L_Matrix J(4,2*(int)rf.size());
	
	// Aca se hacen variar los rayos
	// Hay que usar la formula: dx*/da = Fxx^-1 * Fxa
	// donde: Fxx(i,j) = d2f /( dx(i) dx(j) )
	//        Fxa(i,j) = d2f /( dx(i) da(j) )
	//       mFxa(i,j) = -Fxa(i,j)

	mFxa.reallocate(7, 2*(int)rf.size());

	for (i=0; i<7; i++)
	{
		pose = *this;
		pose.el(i) += delta;
		H.fijaPose3D_cuat(pose);
		// Calcular los hessianos cruzados
		for (j=0; j<(int)rf.size(); j++)
		{
			// El calculo del error es el siguiente:
			L_MatrizHomogenea_OP_mult_P(pun, H, pm[j]);
			dI = pun.y/pun.x - rf[j].tanIzq;
			dA = pun.z/pun.x - rf[j].tanArr;
			// err += dI^2 + dA^2; 
			// En consecuencia:
			//    d(err)/d(rf[i].tanIzq) = -2*dI
			//    d(err)/d(rf[i].tanArr) = -2*dA
			J(0,j) = -2*dI;
			J(1,j) = -2*dA;
		}
		pose = *this;
		pose.el(i) -= delta;
		H.fijaPose3D_cuat(pose);
		// Calcular los hessianos cruzados
		for (j=0; j<(int)rf.size(); j++)
		{
			// El calculo del error es el siguiente:
			L_MatrizHomogenea_OP_mult_P(pun, H, pm[j]);
			dI = pun.y/pun.x - rf[j].tanIzq;
			dA = pun.z/pun.x - rf[j].tanArr;
			// err += dI^2 + dA^2; 
			// En consecuencia:
			//    d(err)/d(rf[i].tanIzq) = -2*dI
			//    d(err)/d(rf[i].tanArr) = -2*dA
			J(2,j) = -2*dI;
			J(3,j) = -2*dA;
		}
		for (j=0; j<(int)rf.size(); j++)  // Aca se define -Fxa
		{
			mFxa(i,2*j+0) = - (J(0,j)-J(2,j)) / (2*delta);
			mFxa(i,2*j+1) = - (J(1,j)-J(3,j)) / (2*delta);
		}
	}
	// Calcular el resultado final
	//J7x2n.OP_mult(Hess7x7inv, J);
	L_Matrix_OP_mult_sub(J7x2n, i0, j0, Hess7x7inv, mFxa, 0, 0, Hess7x7inv.li, Hess7x7inv.lj, 0, 0, mFxa.li, mFxa.lj)
}

void L_Pose3D_cuat::jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin_analitica(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, const L_Matrix &Hess7x7inv, L_Matrix &J7x3n, int i0, int j0, double delta)
{
	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)pm.size(), "L_Pose3D_cuat::jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin_analitica() : matriz de tamano insuficiente");
	// Se sabe que la pose actual es la que lleva los puntos
	// pm a los rayos rf en forma optima
	// error = suma_i ( (eta*pm[i])[2]/(eta*pm[i])[1] - rf[i][1] )^2 + ( (eta*pm[i])[2]/(eta*pm[i])[1] - rf[i][1] )^2 + ( (eta*pm[i])[3]/(eta*pm[i])[1] - rf[i][2] )^2
	//
	// El error es una sum de terminos, cada elemento de la sum depende de solo algunos datos
	// En consecuencia, hay que evaluar solo la ariacion de un termino cada vez y no de la sum completa

	int i, j; // 7 x rf.size()
	L_Pose3D_cuat pose;
	L_CoordsCart3D pun;
	L_HomogeneousMatrix H;
	double dI, dA;
	L_Matrix mFxa;
	L_Matrix J(6,(int)pm.size());
	
	// Aca se hacen variar los rayos
	// Hay que usar la formula: dx*/da = Fxx^-1 * Fxa
	// donde: Fxx(i,j) = d2f /( dx(i) dx(j) )
	//        Fxa(i,j) = d2f /( dx(i) da(j) )
	//       mFxa(i,j) = -Fxa(i,j)

	mFxa.reallocate(7, 3*(int)pm.size());

	for (i=0; i<7; i++)
	{
		pose = *this;
		pose.el(i) += delta;
		H.fijaPose3D_cuat(pose);
		// Calcular los hessianos cruzados
		for (j=0; j<(int)pm.size(); j++)
		{
			// El calculo del error es el siguiente:
			L_MatrizHomogenea_OP_mult_P(pun, H, pm[j]);
			dI = pun.y/pun.x - rf[j].tanIzq;
			dA = pun.z/pun.x - rf[j].tanArr;
			// err += dI^2 + dA^2; 
			// En consecuencia:
			//    d(err)/d(pm[i].x) = 2*dI*-pun.y/(pun.x*pun.x) + 2*dA*-pun.z/(pun.x*pun.x)
			//    d(err)/d(pm[i].y) = 2*dI * 1/pun.x
			//    d(err)/d(pm[i].z) = 2*dA * 1/pun.x
			J(0,j) = 2*(-dI*pun.y - dA*pun.z)/(pun.x*pun.x);
			J(1,j) = 2*dI * 1/pun.x;
			J(2,j) = 2*dA * 1/pun.x;
		}
		pose = *this;
		pose.el(i) -= delta;
		H.fijaPose3D_cuat(pose);
		// Calcular los hessianos cruzados
		for (j=0; j<(int)pm.size(); j++)
		{
			// El calculo del error es el siguiente:
			L_MatrizHomogenea_OP_mult_P(pun, H, pm[j]);
			dI = pun.y/pun.x - rf[j].tanIzq;
			dA = pun.z/pun.x - rf[j].tanArr;
			// err += dI^2 + dA^2; 
			// En consecuencia:
			//    d(err)/d(pm[i].x) = dI*-pun.y/(pun.x*pun.x) + dA*-pun.z/(pun.x*pun.x)
			//    d(err)/d(pm[i].y) = dI * 1/pun.x
			//    d(err)/d(pm[i].z) = dA * 1/pun.x
			J(3,j) = 2*(-dI*pun.y - dA*pun.z)/(pun.x*pun.x);
			J(4,j) = 2*dI * 1/pun.x;
			J(5,j) = 2*dA * 1/pun.x;
		}
		for (j=0; j<(int)pm.size(); j++) // Aca se compute -Fxa
		{
			mFxa(i,3*j+0) = - (J(0,j)-J(3,j)) / (2*delta);
			mFxa(i,3*j+1) = - (J(1,j)-J(4,j)) / (2*delta);
			mFxa(i,3*j+2) = - (J(2,j)-J(5,j)) / (2*delta);
		}
	}
	// Calcular el resultado final
	//J7x3n.OP_mult(Hess7x7inv, J);
	L_Matrix_OP_mult_sub(J7x3n, i0, j0, Hess7x7inv, mFxa, 0, 0, Hess7x7inv.li, Hess7x7inv.lj, 0, 0, mFxa.li, mFxa.lj)
}


void L_Pose3D_cuat::jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin(const std::vector<L_CoordsCart3D> &pm, const std::vector<L_RayoCamara> &rf, L_Matrix &J7x3n, int i0, int j0, L_Pose3D_cuat *inic, double delta, double epsilon, double errorToStop, int nIterationsMax, int nIteracParar2)
{
	L_Pose3D_cuat p1, p2;
	std::vector<L_CoordsCart3D> p;
	int i, j;
	L_LevenbergMarquardt opt;

	opt.epsilon = epsilon; // 1.0e-8; // epsilon para las variaciones
	opt.errorToStop = errorToStop; // 1.0e-6; // si la variacion del error entre iteraciones es menor a esto, se para
	opt.nIterationsMax = nIterationsMax; // 10;  // 10 iteraciones como maxElement

	throw_L_ArgException_if(J7x3n.li < i0+7 || J7x3n.lj < j0+3*(int)p.size(), "L_Pose3D_cuat::jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin() : matriz de tamano insuficiente");

	p = pm;

	if (inic == NULL)
		calcPose_movido_a_fijo_Nrayos(p, rf, opt);
	else
		*this = *inic;
	for (j=0; j<3*(int)p.size(); j++)
	{
		opt.nIterationsMax = nIteracParar2; // 4;
		p[j/3].el(j%3) += delta;   // Error: decia j/2
		p1.calcPose_movido_a_fijo_Nrayos(p, rf, opt, this);
		if (p1.ori.OP_mult_elementwise(ori) < 0)
			p1.ori.invSigno();
		p[j/3].el(j%3) -= 2*delta;
		p2.calcPose_movido_a_fijo_Nrayos(p, rf, opt, this);
		if (p2.ori.OP_mult_elementwise(ori) < 0)
			p2.ori.invSigno();
		p[j/3].el(j%3) += delta;
		for (i=0; i<7; i++)
			J7x3n(i+i0,j+j0) = (p1.el(i) - p2.el(i)) / (2*delta);
	}
}


L_CoordsCart3D L_HomogeneousMatrix::operator *(const L_CoordsCart3D &other)
{
	L_CoordsCart3D ret;
	ret.x=other.x*operator()(0,0)+other.y*operator()(0,1)+other.z*operator()(0,2)+operator()(0,3);
	ret.y=other.x*operator()(1,0)+other.y*operator()(1,1)+other.z*operator()(1,2)+operator()(1,3);
	ret.z=other.x*operator()(2,0)+other.y*operator()(2,1)+other.z*operator()(2,2)+operator()(2,3);
	return ret;
}

L_CoordsCart3D L_HomogeneousMatrix::movido_a_fijo(const L_CoordsCart3D &other)
{
	L_CoordsCart3D ret;
	ret.x=other.x*operator()(0,0)+other.y*operator()(0,1)+other.z*operator()(0,2)+operator()(0,3);
	ret.y=other.x*operator()(1,0)+other.y*operator()(1,1)+other.z*operator()(1,2)+operator()(1,3);
	ret.z=other.x*operator()(2,0)+other.y*operator()(2,1)+other.z*operator()(2,2)+operator()(2,3);
	return ret;
}

L_CoordsCart3D L_HomogeneousMatrix::fijo_a_movido(const L_CoordsCart3D &other)
{
	L_CoordsCart3D ret; // H^-1*p = [RT -RT*t) * p = (RT*p - RT*t) + RT(p-t)
	double dx = other.x-operator()(0,3);
	double dy = other.y-operator()(1,3);
	double dz = other.z-operator()(2,3);
	ret.x=dx*operator()(0,0)+dy*operator()(1,0)+dz*operator()(2,0);
	ret.y=dx*operator()(0,1)+dy*operator()(1,1)+dz*operator()(2,1);
	ret.z=dx*operator()(0,2)+dy*operator()(1,2)+dz*operator()(2,2);
	return ret;
}

L_CoordsCart3D L_HomogeneousMatrix::OP_rota(const L_CoordsCart3D &other)
{
	L_CoordsCart3D ret;
	ret.x=other.x*operator()(0,0)+other.y*operator()(0,1)+other.z*operator()(0,2);
	ret.y=other.x*operator()(1,0)+other.y*operator()(1,1)+other.z*operator()(1,2);
	ret.z=other.x*operator()(2,0)+other.y*operator()(2,1)+other.z*operator()(2,2);
	return ret;
}

void L_CoordsEsfInv3D::jacob_pideCoordsCart3D(L_Matrix &J3x3) const
{
	// ret.x=1/invR*cos(theta)*cos(phi); ret.y=1/invR*sin(theta)*cos(phi); ret.z=1/invR*sin(phi);  pos -> invR, theta, phi
	double rho = 1/invR;
	double rho2 = rho*rho;
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);

	// Parte esferica-inversa de las coordenadas 6d: dx3d.(invR,theta,phi) / dx6d.puni.pos
	J3x3(0,0) =  -cos_theta*cos_phi*rho2;   // dx/dinvR
	J3x3(0,1) =  -sin_theta*cos_phi*rho;    // dx/dtheta
	J3x3(0,2) =  -cos_theta*sin_phi*rho;    // dx/dphi
	J3x3(1,0) =  -sin_theta*cos_phi*rho2;   // dy/dinvR
	J3x3(1,1) =   cos_theta*cos_phi*rho;    // dy/dtheta
	J3x3(1,2) =  -sin_theta*sin_phi*rho;    // dy/dphi
	J3x3(2,0)  = -sin_phi*rho2;             //dz/dinvR
	J3x3(2,1)  =  0;                        //dz/dtheta
	J3x3(2,2)  =  cos_phi*rho;              //dz/dphi}
}

void L_Pose3D_cuat::pruebaPoseCuaternion(int nPrueba)
{
	L_Pose3D po1, po2, po3;
	L_HomogeneousMatrix H1, H2, H3;
	L_Pose3D_cuat pc1, pc2, pc3, pctres, pc4, ptmp;
	L_CoordsCart3D pun1, pun2, pun3;
	L_Matrix J, J2, J3;
	L_StaticMatrix<7,7> Jv2;
	std::vector<L_CoordsCart3D> pa1(30), pa2(30);
	std::vector<L_RayoCamara> rM(30), rF(30);
	L_Pose3D_cuat p3dc, p3dd;
	//double coefs[7], coefs3[7], coefstres[7];
	double pre36[36];
	double delta = 1.0e-8;
	int i, j;
	int npr = 1, nprmax = 21;
	int nIterTiempos = 5;
	double t1, t2;

	po1.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),
		L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	po2.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),
		L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	H1.fijaPose3D(po1);
	H2.fijaPose3D(po2);
	pc1=H1.calcPose3D_cuat();
	pc2=H2.calcPose3D_cuat();
	pc3 = pc1;
	pctres = pc2;
	ptmp = pc1;
	p3dc = pc1;
	p3dd = pc2;
	pun1.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	pun2.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	pun3.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	for (i=0; i<(int)pa1.size(); i++)
		pa1[i].define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	for (i=0; i<(int)pa2.size(); i++)
		pa2[i].define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	for (i=0; i<(int)rM.size(); i++)
		rM[i].define(L_RANDOM(-1,1), L_RANDOM(-1,1));
	for (i=0; i<(int)rF.size(); i++)
		rF[i].define(L_RANDOM(-1,1), L_RANDOM(-1,1));

	if (nPrueba == -1 || nPrueba == 0)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Generando poses pan-tilt-roll y homogeneas\n");

		po1.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),
			L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
		po2.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),
			L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
		po1.imprime_gr("po1");
		po2.imprime_gr("po2");
		printf("## probando pose3D a homogenea y viceversa\n");
		H1.fijaPose3D(po1);
		H2.fijaPose3D(po2);
		po1=H1.calcPose3D();
		po2=H2.calcPose3D();
		po1.imprime_gr("po1");
		po2.imprime_gr("po2");

		printf(" ## Probando homogenea a pose-cuaternion: calcPose3D_cuat()\n");
		pc1=H1.calcPose3D_cuat();
		pc2=H2.calcPose3D_cuat();
		pc1.print("pc1");
		pc2.print("pc2");
		printf(" ## Probando pose-cuaternion a homogenea: fijaPose3D_cuat()\n");
		H1.fijaPose3D_cuat(pc1);
		H2.fijaPose3D_cuat(pc2);
		printf("normas homogeneas(=3): %.3g %.3g\n", H1.normaL2rot(), H2.normaL2rot());
		printf(" ## Probando homogenea a pose-cuaternion: calcPose3D_cuat()\n");
		pc1=H1.calcPose3D_cuat();
		pc2=H2.calcPose3D_cuat();
		pc1.print("pc1");
		pc2.print("pc2");
		printf(" ## Probando pose-cuaternion a homogenea: fijaPose3D_cuat_fijo_a_movido()\n");
		H1.fijaPose3D_cuat_fijo_a_movido(pc1);
		H2.fijaPose3D_cuat_fijo_a_movido(pc2);
		printf("normas homogeneas(=3): %.3g %.3g\n", H1.normaL2rot(), H2.normaL2rot());
		printf(" ## Probando homogenea a pose-cuaternion: calcPose3D_cuat_fijo_a_movido()\n");
		pc1=H1.calcPose3D_cuat_fijo_a_movido();
		pc2=H2.calcPose3D_cuat_fijo_a_movido();
		pc1.print("pc1");
		pc2.print("pc2");
		printf(" ## Probando matrices homogeneas a poses pan-tilt-roll: fijaPose3D_cuat(), calcPose3D()\n");
		po1.imprime_gr("po1");
		po2.imprime_gr("po2");
		H1.fijaPose3D_cuat(pc1);
		H2.fijaPose3D_cuat(pc2);
		po1=H1.calcPose3D();
		po2=H2.calcPose3D();
		po1.imprime_gr("po1");
		po2.imprime_gr("po2");

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 1)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando jacobiano de calcPose3D_cuat()\n");
		J.reallocate(7,12);
		pc1 = H1.jacob_calcPose3D_cuat(J);
		J2.reallocate(7,12);
		for (j=0; j<12; j++)
		{
			H1(j/4,j%4) += delta;
			pc3 = H1.calcPose3D_cuat();
			H1(j/4,j%4) -= 2*delta;
			pctres = H1.calcPose3D_cuat();
			H1(j/4,j%4) += delta;
			for (i=0; i<7; i++)
				J2(i,j) = (pc3.el(i)-pctres.el(i)) / (2*delta);
		}
		J.print("J7x12");
		J2.print("J7x12");

		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();

		printf(" ## Probando pq1*pq2 = p3dcuat( H(pq1)*H(pq2) )\n");
		pc3 = pc1*pc2;
		H1.fijaPose3D_cuat(pc1);
		H2.fijaPose3D_cuat(pc2);
		H3 = H1*H2;

		pctres = H3.calcPose3D_cuat();

		pc3.print("q1*q2");
		pctres.print("cuat( H(q1)*H(q2) )");

		printf("Parecido: %f\n", pc3.cosineDistanceTo(pctres));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 2)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Rotacion vector con pose-cuaternion\n");
		pun1.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
		pun2 = pc1.movido_a_fijo(pun1);
		pun1.print("pun1");
		pun2.print("pun2");

		printf(" ## Rotacion vector con homogenea\n");
		H1.fijaPose3D_cuat(pc1);
		pun3 = H1*pun1;
		pun1.print("pun1");
		pun3.print("pun2");
		printf("Parecido: %f\n", pun2.cosineDistanceTo(pun3));
		printf("\n");

		printf(" ## Rotacion vector inversa con pose-cuaternion\n");
		pun1.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
		pun2 = pc1.fijo_a_movido(pun1);
		pun1.print("pun1");
		pun2.print("pun2");

		printf(" ## Rotacion vector con homogenea inversa\n");
		H1.fijaPose3D_cuat_fijo_a_movido(pc1);
		pun3 = H1*pun1;
		pun1.print("pun1");
		pun3.print("pun2");
		printf("Parecido: %f\n", pun2.cosineDistanceTo(pun3));

		printf("\n");
		printf(" ## Rotacion vector con pose-cuaternion usando inversa de pose-cuaternion\n");
		pc2.inverseOf(pc1);
		pun3 = pc2.movido_a_fijo(pun1);
		pun1.print("pun1");
		pun3.print("pun2");
		printf("Parecido: %f\n", pun2.cosineDistanceTo(pun3));

		printf(" ## Rotacion vector con homogenea usando inversa de homogenea\n");
		H1.fijaPose3D_cuat(pc1);
		H2.inversaTrasponiendoDe(H1);
		pun3 = H2*pun1;
		pun1.print("pun1");
		pun3.print("pun2");
		printf("Parecido: %f\n", pun2.cosineDistanceTo(pun3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 3)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Jacobiano de rotacion vector con pose-cuaternion\n");
		J.reallocate(3,7);
		pc1.jacob_izq_movido_a_fijo(pun1, J);
		J.print("J3x7");
		J2.reallocate(3,7);
		for (j=0; j<7; j++)
		{
			pc2 = pc1; pc2.el(j)+=delta;
			pc3 = pc1; pc3.el(j)-=delta;
			pun2 = pc2.movido_a_fijo(pun1);
			pun3 = pc3.movido_a_fijo(pun1);
			for (i=0; i<3; i++)
				J2(i,j) = (pun2.el(i)-pun3.el(i)) / (2*delta);
		}
		J2.print("J3x7(dfn)");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		J.reallocate(3,3);
		pc1.jacob_der_movido_a_fijo(J);
		J.print("J3x3");

		J2.reallocate(3,3);
		for (j=0; j<3; j++)
		{
			pun2 = pun1; pun2.el(j)+=delta;
			pun3 = pun1; pun3.el(j)-=delta;
			pun2 = pc1.movido_a_fijo(pun2);
			pun3 = pc1.movido_a_fijo(pun3);
			for (i=0; i<3; i++)
				J2(i,j) = (pun2.el(i)-pun3.el(i)) / (2*delta);
		}
		J2.print("J3x3(dfn)");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 4)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando inversa de pose-cuaternion\n");
		pctres.inverseOf(pc1);
		ptmp.OP_mult(pc1,pctres);
		ptmp.print("pc1*pc1inv");
		ptmp.OP_mult(pctres, pc1);
		ptmp.print("pc1inv*pc1");

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 5)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando jacobiano de la inversa de pose-cuaternion\n");
		J.reallocate(7,7);
		pctres.jacob_inverseOf(pc1, J);
		J.print("Jinv");
		pctres.inverseOf(pc1);
		J2.reallocate(7,7);
		for (j=0; j<7; j++)
		{
			pc1.el(j)+=delta;
			ptmp.inverseOf(pc1);
			for (i=0; i<7; i++)
				J2(i,j) = (ptmp.el(i) - pctres.el(i)) / delta;
			pc1.el(j)-=delta;
		}
		J2.print("Jinv(dif.fin)");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 6)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano izquierdo de la multiplicacion (directo)\n");
		J.reallocate(7,7);
		L_Pose3D_cuat::jacob_izq_OP_mult(pc1, pc2, J);
		J.print("J");

		printf(" ## Probando jacobiano izquierdo de la multiplicacion (precalc)\n");
		// OJO
		L_Pose3D_cuat::jacob_izq_OP_mult_pre(pc1,pre36);  // Decia pc2
		L_Pose3D_cuat::jacob_izq_OP_mult(pc1, pc2, pre36, Jv2); // Jv2 es L_StaticMatrix<7,7>
		Jv2.print("J(pre)");
		L_Matrix_OP_assign(J2,Jv2);

		printf(" ## Probando jacobiano izquierdo de la multiplicacion (dif finitas)\n");
		pc3 = pc1*pc2;
		J3.reallocate(7,7);
		for (j=0; j<7; j++)
		{
			// Variacion a pc1
			pc1.el(j) += delta;
			pctres = pc1*pc2;
			pc1.el(j) -= delta;
			// derivadas
			for (i=0; i<7; i++)
				J3(i,j) = (pctres.el(i) - pc3.el(i)) / delta;
		}
		J3.print("J");
		printf("Parecido: %f %f\n", J.cosineDistanceTo(J2), J.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);

		getchar();
	}

	if (nPrueba == -1 || nPrueba == 7)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano derecho de la multiplicacion (directo)\n");
		J.reallocate(7,7);
		J3.reallocate(7,7);
		L_Pose3D_cuat::jacob_der_OP_mult(pc1, pc2, J);
		J.print("J");

		printf(" ## Probando jacobiano derecho de la multiplicacion (dif finitas)\n");
		J3.reallocate(7,7);
		pc3 = pc1*pc2;
		for (j=0; j<7; j++)
		{
			// Variacion a pc1
			pc2.el(j) += delta;
			pctres = pc1*pc2;
			pc2.el(j) -= delta;
			// derivadas
			for (i=0; i<7; i++)
				J3(i,j) = (pctres.el(i) - pc3.el(i)) / delta;
		}
		J3.print("J");
		printf("Parecido: %f\n", J.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 8)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		pc1.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));
		pc1.normaliza();
		pc2.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));
		pc2.normaliza();
		printf(" ## Probando OP_mult_p1inv_por_p2()\n");
		pc3.OP_mult_p1inv_por_p2(pc1, pc2);
		pc3.print("OP_mult_p1inv_por_p2");
		ptmp.inverseOf(pc1);
		pc4 = ptmp*pc2;
		pc4.print("inverseOf(pc1)*pc2");
		printf("Parecido: %f\n", pc3.cosineDistanceTo(pc4));

		printf(" ## Probando OP_mult_p1_por_p2inv()\n");
		pc3.OP_mult_p1_por_p2inv(pc1, pc2);
		pc3.print("OP_mult_p1_por_p2inv");
		ptmp.inverseOf(pc2);
		pc4 = pc1*ptmp;
		pc4.print("pc1*inverseOf(pc2)");
		printf("Parecido: %f\n", pc3.cosineDistanceTo(pc4));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 9)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano izq de OP_mult_p1inv_por_p2() (directo)\n");
		J.reallocate(7,7);
		J2.reallocate(7,7);
		L_Pose3D_cuat::jacob_izq_OP_mult_p1inv_por_p2(pc1, pc2, J);
		J.print("J");
		printf("\n");
		printf(" ## Probando jacobiano izq de OP_mult_p1inv_por_p2() (dif finitas)\n");
		pc3.OP_mult_p1inv_por_p2(pc1, pc2);
		for (j=0; j<7; j++)
		{
			pc1.el(j) += delta;
			pc3.OP_mult_p1inv_por_p2(pc1, pc2);
			pc1.el(j) -= 2*delta;
			pc4.OP_mult_p1inv_por_p2(pc1, pc2);
			pc1.el(j) += delta;
			
			// derivadas
			for (i=0; i<7; i++)
				J2(i,j) = (pc3.el(i) - pc4.el(i)) / (2*delta);
		}
		J2.print("J");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 10)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano der de OP_mult_p1inv_por_p2() (directo)\n");
		J.reallocate(7,7);
		J2.reallocate(7,7);
		L_Pose3D_cuat::jacob_der_OP_mult_p1inv_por_p2(pc1, pc2, J);
		J.print("J");
		printf("\n");
		printf(" ## Probando jacobiano der de OP_mult_p1inv_por_p2() (dif finitas)\n");
		for (j=0; j<7; j++)
		{
			pc2.el(j) += delta;
			pc3.OP_mult_p1inv_por_p2(pc1, pc2);
			pc2.el(j) -= 2*delta;
			pc4.OP_mult_p1inv_por_p2(pc1, pc2);
			pc2.el(j) += delta;
			
			// derivadas
			for (i=0; i<7; i++)
				J2(i,j) = (pc3.el(i) - pc4.el(i)) / (2*delta);
		}
		J2.print("J");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}


	if (nPrueba == -1 || nPrueba == 11)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano izq de OP_mult_p1_por_p2inv() (directo)\n");
		J.reallocate(7,7);
		J2.reallocate(7,7);
		L_Pose3D_cuat::jacob_izq_OP_mult_p1_por_p2inv(pc1, pc2, J);
		J.print("J");
		printf("\n");
		printf(" ## Probando jacobiano izq de OP_mult_p1_por_p2inv() (dif finitas)\n");
		pc3.OP_mult_p1inv_por_p2(pc1, pc2);
		for (j=0; j<7; j++)
		{
			pc1.el(j) += delta;
			pc3.OP_mult_p1_por_p2inv(pc1, pc2);
			pc1.el(j) -= 2*delta;
			pc4.OP_mult_p1_por_p2inv(pc1, pc2);
			pc1.el(j) += delta;
			
			// derivadas
			for (i=0; i<7; i++)
				J2(i,j) = (pc3.el(i) - pc4.el(i)) / (2*delta);
		}
		J2.print("J");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 12)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf(" ## Probando jacobiano der de OP_mult_p1_por_p2inv() (directo)\n");
		J.reallocate(7,7);
		J2.reallocate(7,7);
		L_Pose3D_cuat::jacob_der_OP_mult_p1_por_p2inv(pc1, pc2, J);
		J.print("J");
		printf("\n");
		printf(" ## Probando jacobiano der de OP_mult_p1_por_p2inv() (dif finitas)\n");
		for (j=0; j<7; j++)
		{
			pc2.el(j) += delta;
			pc3.OP_mult_p1_por_p2inv(pc1, pc2);
			pc2.el(j) -= 2*delta;
			pc4.OP_mult_p1_por_p2inv(pc1, pc2);
			pc2.el(j) += delta;
			
			// derivadas
			for (i=0; i<7; i++)
				J2(i,j) = (pc3.el(i) - pc4.el(i)) / (2*delta);
		}
		J2.print("J");
		printf("Parecido: %f\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 13)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando calc_movido_a_fijo con varios puntos (%ld puntos)\n", long(pa1.size()));
		pc1.print("Pose correcta");

		pa1.resize(30);
		pa2.resize(pa1.size());

		printf("Probando calc_movido_a_fijo con varios puntos (%ld puntos)\n", long(pa1.size()));

		for (i=0; i<(int)pa1.size(); i++)
			pa1[i].define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]);

		pc1.print("pc1");

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_rapidoMuyImpreciso(pa1, pa2);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_rapidoMuyImpreciso() sin ruido %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_matriz(pa1, pa2);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_matriz() sin ruido %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_iter(pa1, pa2, 10);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_iter() (10 iter) sin ruido %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_iter(pa1, pa2, 30);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_iter() (30 iter) sin ruido %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcpose_movido_a_fijo_3puntos(pa1, pa2);
		t2 = L_TIME();
		printf("calcpose_movido_a_fijo_3puntos() %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 14)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		pc1.print("Pose correcta");

		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]) + L_CoordsCart3D(L_RANDOM(-0.5,0.5), L_RANDOM(-0.5,0.5), L_RANDOM(-0.5,0.5));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_rapidoMuyImpreciso(pa1, pa2);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_rapidoMuyImpreciso() con ruido 5%% %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_matriz(pa1, pa2);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_matriz() con ruido 5%% %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_iter(pa1, pa2, 10);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_iter() (10 iter) con ruido 5%% %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcPose_movido_a_fijo_iter(pa1, pa2, 30);
		t2 = L_TIME();
		printf("calcPose_movido_a_fijo_iter() (30 iter) con ruido 5%% %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		t1 = L_TIME();
		for (i=0; i<nIterTiempos; i++)
			p3dc.calcpose_movido_a_fijo_3puntos(pa1, pa2);
		t2 = L_TIME();
		printf("calcpose_movido_a_fijo_3puntos() con ruido 5%% %g seg:\n", (t2-t1)/50);
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		printf("Probando calc_fijo_a_movido con varios puntos (%ld puntos, roles invertidos)\n", long(pa1.size()));

		p3dc.calcPose_fijo_a_movido_matriz(pa2, pa1);
		printf("calcPose_fijo_a_movido_matriz() con ruido 5%%\n");
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		p3dc.calcPose_fijo_a_movido_matriz_v2(pa2, pa1);
		printf("calcPose_fijo_a_movido_matriz_v2() con ruido 5%%\n");
		p3dc.print("p3dc");
		printf("Parecido: %f\n", pc1.cosineDistanceTo(p3dc));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 15)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		double t1, t2;
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]);

		J.reallocate(7, 3*(int)pa1.size());
		J2.reallocate(7, 3*(int)pa1.size());
		J3.reallocate(7, 3*(int)pa1.size());

		printf("Probando jacob_fijo_calcPose_movido_a_fijo_matriz\n");
		t1 = L_TIME();
		pc1.jacob_fijo_calcPose_movido_a_fijo_matriz(pa1, pa2, J);
		t2 = L_TIME();
		J.print("J(matriz)");
		printf("fps: %f\n", 1/(t2-t1));
		t1 = L_TIME();
		pc1.jacob_fijo_calcPose_movido_a_fijo_matriz_diffin(pa1, pa2, J2);
		t2 = L_TIME();
		J2.print("J(matr.df)"); // PROBANDO ACA
		printf("fps: %f\n", 1/(t2-t1));
		t1 = L_TIME();
		pc1.jacob_fijo_calcPose_movido_a_fijo_iter(pa1, pa2, 10, J3);
		t2 = L_TIME();
		J3.print("J(itr.df)"); // PROBANDO ACA
		printf("fps: %f\n", 1/(t2-t1));
		printf("Parecido 1(matr) 2(matr.df): %g\n", J.cosineDistanceTo(J2));
		printf("Parecido 1(matr) 3(itr.df): %g\n", J.cosineDistanceTo(J3));
		printf("Parecido 2(matr.df) 3(itr.df): %g\n", J2.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 16)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		double t1, t2;
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]);

		J.reallocate(7, 3*(int)pa1.size());
		J2.reallocate(7, 3*(int)pa1.size());
		J3.reallocate(7, 3*(int)pa1.size());

		printf("Probando jacob_movido_calcPose_movido_a_fijo_matriz\n");

		t1 = L_TIME();
		pc1.jacob_movido_calcPose_movido_a_fijo_matriz(pa1, pa2, J);
		t2 = L_TIME();
		J.print("J(matriz)");
		printf("fps: %f\n", 1/(t2-t1));
		t1 = L_TIME();
		pc2.fijaAzar();
		pc2.jacob_movido_calcPose_movido_a_fijo_matriz_diffin(pa1, pa2, J2);
		t2 = L_TIME();
		J2.print("J(matr.df)"); // PROBANDO ACA
		printf("fps: %f\n", 1/(t2-t1));
		t1 = L_TIME();
		pc3.fijaAzar();
		pc3.jacob_movido_calcPose_movido_a_fijo_iter(pa1, pa2, 10, J3);
		t2 = L_TIME();
		J3.print("J(dif.itr)"); // PROBANDO ACA
		printf("fps: %f\n", 1/(t2-t1));
		printf("Concordancia orientaciones: %f %f %f\n", pc1.ori.OP_mult_elementwise(pc2.ori), pc1.ori.OP_mult_elementwise(pc3.ori), pc2.ori.OP_mult_elementwise(pc3.ori));
		printf("Parecido res: %f %f %f\n", pc1.cosineDistanceTo(pc2), pc1.cosineDistanceTo(pc3), pc2.cosineDistanceTo(pc3));	
		printf("Parecido 1(matr) 2(matr.df): %g (mula)\n", J.cosineDistanceTo(J2));
		printf("Parecido 1(matr) 3(itr.df): %g (mula)\n", J.cosineDistanceTo(J3));
		printf("Parecido 2(matr.df) 3(itr.df): %g\n", J2.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 17)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		pc1.print("Pose correcta");
		printf("Probando calcPose_movido_a_fijo_4rayos(rm,pf)\n");

		pa1.resize(8);
		pa2.resize(pa1.size());

		for (i=0; i<(int)pa1.size(); i++)
		{
			pa1[i].x = L_RANDOM(1, 10);
			pa1[i].y = pa1[i].x * L_RANDOM(-2, 2);
			pa1[i].z = pa1[i].x * L_RANDOM(-2, 2);
			rM[i].define(pa1[i].y/pa1[i].x, pa1[i].z/pa1[i].x);
		}
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]);

		pc2.fijaAzar();
		pc2._calcPose_movido_a_fijo_4rayos(rM, pa2);
		pc2.print("p3dc");
		printf("Probando calcPose_movido_a_fijo_Nrayos(rm,pf) (%ld rayos)\n", long(rM.size()));
		pc3.fijaAzar();
		pc3._calcPose_movido_a_fijo_Nrayos(rM, pa2, 10);
		pc3.print("p3dc");

		printf("Parecido: %f %f %f\n", pc1.cosineDistanceTo(pc2), pc1.cosineDistanceTo(pc3), pc2.cosineDistanceTo(pc3));

		printf("Probando calcPose_movido_a_fijo_4rayos(rm,pf) ruido 1 pixel\n");
		for (i=0 ; i<(int)rM.size(); i++)
			rM[i].define(rM[i].tanIzq + L_RANDOM(-1.0/320, 1.0/320), rM[i].tanArr + L_RANDOM(-1.0/320, 1.0/320));

		pc2._calcPose_movido_a_fijo_4rayos(rM, pa2);
		pc2.print("p3dc");
		printf("Probando calcPose_movido_a_fijo_Nrayos(rm,pf) (%ld rayos) 10 iter. ruido 1 pixel\n", long(rM.size()));
		pc3._calcPose_movido_a_fijo_Nrayos(rM, pa2, 10);
		pc3.print("p3dc");
		printf("Probando calcPose_movido_a_fijo_Nrayos(rm,pf) (%ld rayos) 30 iter. ruido 1 pixel\n", long(rM.size()));
		pc4._calcPose_movido_a_fijo_Nrayos(rM, pa2, 30);
		pc4.print("p3dc");

		printf("Parecido: %f %f %f %f\n", pc1.cosineDistanceTo(pc2), pc1.cosineDistanceTo(pc3), pc1.cosineDistanceTo(pc4), pc2.cosineDistanceTo(pc4));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 18)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando jacobianos de calcPose_movido_a_fijo_Nrayos(rm,pf) - el metodo matricial es penca\n");

		pa1.resize(8);
		pa2.resize(pa1.size());
		rM.resize(pa1.size());

		for (i=0; i<(int)pa1.size(); i++)
		{
			pa1[i].x = L_RANDOM(1, 10);
			pa1[i].y = pa1[i].x * L_RANDOM(-2, 2);
			pa1[i].z = pa1[i].x * L_RANDOM(-2, 2);
			rM[i].define(pa1[i].y/pa1[i].x, pa1[i].z/pa1[i].x);
		}
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.movido_a_fijo(pa1[i]);

		// pc1*rM -> pa2

		J.reallocate(7, 2*(int)rM.size());
		p3dc = pc1;
		p3dc._jacob_rm_calcPose_movido_a_fijo_Nrayos_matriz_malo(rM, pa2, J); // Requiere la pose ya calculada
		J.print("J(matriz)");
		printf("\n");
		J2.reallocate(7, 2*(int)rM.size());
		pc2.fijaAzar();
		p3dc = pc2; // Para que no sea predecible
		p3dc._jacob_rm_calcPose_movido_a_fijo_Nrayos_diffin(rM, pa2, J2);
		J2.print("J(itr.df)");

		printf("Parecido 1(matriz) 2(itr.df): %g (mula)\n", J.cosineDistanceTo(J2));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}



	if (nPrueba == -1 || nPrueba == 19)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		pc1.print("Pose correcta");
		printf("Probando calcPose_movido_a_fijo_4rayos(pm,rf)\n");
		L_LevenbergMarquardt opt;

		opt.epsilon = 1.0e-6;
		opt.errorToStop = 1e-10;
		opt.nIterationsMax = 10;

		pa1.resize(8);
		pa2.resize(pa1.size());
		rF.resize(pa1.size());

		for (i=0; i<(int)pa1.size(); i++)
		{
			pa1[i].x = L_RANDOM(1, 10);
			pa1[i].y = pa1[i].x * L_RANDOM(-2, 2);
			pa1[i].z = pa1[i].x * L_RANDOM(-2, 2);
			rF[i].define(pa1[i].y/pa1[i].x, pa1[i].z/pa1[i].x);
		}
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.fijo_a_movido(pa1[i]);

		p3dc.calcPose_movido_a_fijo_4rayos(pa2, rF);
		p3dc.print("p3dc");
		printf("Probando calcPose_movido_a_fijo_Nrayos(pm,rf) (%ld rayos)\n", long(rF.size()));
		p3dd.calcPose_movido_a_fijo_Nrayos(pa2, rF, opt);
		p3dd.print("p3dc");
		printf("Parecido: %f\n", p3dc.cosineDistanceTo(p3dd));

		printf("Probando calcPose_movido_a_fijo_4rayos(pm,rf) ruido 1 pixel\n");
		for (i=0 ; i<(int)rM.size(); i++)
			rM[i].define(rM[i].tanIzq + L_RANDOM(-1.0/320, 1.0/320), rM[i].tanArr + L_RANDOM(-1.0/320, 1.0/320));

		p3dc.calcPose_movido_a_fijo_4rayos(pa2, rF);
		p3dc.print("p3dc");
		printf("Probando calcPose_movido_a_fijo_Nrayos(pm,rf) (%ld rayos) 10 iter. ruido 1 pixel\n", long(rF.size()));
		pc4.calcPose_movido_a_fijo_Nrayos(pa2, rF, opt);
		pc4.print("p3dc");
		printf("Parecido: %f\n", p3dc.cosineDistanceTo(pc4));
		printf("Probando calcPose_movido_a_fijo_Nrayos(pm,rf) (%ld rayos) 30 iter. ruido 1 pixel\n", long(rF.size()));
		opt.nIterationsMax = 30;
		p3dd.calcPose_movido_a_fijo_Nrayos(pa2, rF, opt);
		p3dd.print("p3dc");
		printf("Parecido: %f\n", p3dc.cosineDistanceTo(p3dd));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 20)
	{
		if (nPrueba != -1)
			npr = nPrueba;
		printf("Probando jacobianos de calcPose_movido_a_fijo_Nrayos(pm,rf) - el metodo matricial es \"incorrecto\"\n");

		pa1.resize(10);
		pa2.resize(pa1.size());
		rF.resize(pa1.size());

		for (i=0; i<(int)pa1.size(); i++)
		{
			pa1[i].x = L_RANDOM(1, 10);
			pa1[i].y = pa1[i].x * L_RANDOM(-2, 2);
			pa1[i].z = pa1[i].x * L_RANDOM(-2, 2);
			rF[i].define(pa1[i].y/pa1[i].x, pa1[i].z/pa1[i].x);
		}
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.fijo_a_movido(pa1[i]);

		// rF fijo, pa2 movido
		// Metodo matricial (incorrecto)
		J.reallocate(7, 2*(int)rF.size());
		p3dc = pc1; // el resultado correcto es pc1
		p3dc.jacob_rf_calcPose_movido_a_fijo_Nrayos_matriz_malo(pa2, rF, J);
		J.print("J(matriz)");
		printf("\n");
		// Metodo usando diferencias finitas sobre metodo iterativo
		J2.reallocate(7, 2*(int)rF.size());
		//p3dc = pc2; // Para que no sea predecible
		p3dc = pc1; // Para que no sea predecible
		p3dc.jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin(pa2, rF, J2, 0, 0, NULL, NULL, 1.0e-3, 1.0e-6, 0.0, 30, 10);
		J2.print("J(itr.df)");
		// metodo analitico
		J3.reallocate(7, 2*(int)rF.size());
		p3dc = pc1;
		L_Matrix Hi;
		if (p3dc.jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica(pa2, rF, Hi) == false)
			printf("Fallo en jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica()\n");
		p3dc.jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin_analitica(pa2, rF, Hi, J3);
		J3.print("J(analitica)");

		printf("Relaciones entre las matrices J(dif.fin)./J(analitica)\n");
		for (int i=0; i<J2.li; i++)
		{
			for (int j=0; j<J2.lj; j++)
				printf("%.2g ", J2(i,j) / J3(i,j));
			printf("\n");
		}

		printf("Parecido: 1(matriz)  2(itr.df) %g\n", J.cosineDistanceTo(J2));
		printf("Parecido: 1(matriz)  3(analitica) %g\n", J.cosineDistanceTo(J3));
		printf("Parecido: 2(itr.df)  3(analitica) %g\n", J2.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}

	if (nPrueba == -1 || nPrueba == 21)
	{
		if (nPrueba != -1)
			npr = nPrueba;

		pa1.resize(10);
		pa2.resize(pa1.size());
		rF.resize(pa1.size());

		for (i=0; i<(int)pa1.size(); i++)
		{
			pa1[i].x = L_RANDOM(1, 10);
			pa1[i].y = pa1[i].x * L_RANDOM(-2, 2);
			pa1[i].z = pa1[i].x * L_RANDOM(-2, 2);
			rF[i].define(pa1[i].y/pa1[i].x, pa1[i].z/pa1[i].x);
		}
		for (i=0; i<(int)pa2.size(); i++)
			pa2[i] = pc1.fijo_a_movido(pa1[i]);

		// Metodo usando diferencias finitas sobre metodo iterativo
		J2.reallocate(7, 3*(int)rF.size());
		//p3dc = pc2; // Para que no sea predecible
		p3dc = pc1; // Para que no sea predecible
		p3dc.jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin(pa2, rF, J2, 0, 0, NULL, 1.0e-3, 1.0e-6, 0.0, 30, 10);
		J2.print("J(itr.df)");
		// metodo analitico
		J3.reallocate(7, 3*(int)rF.size());
		p3dc = pc1;
		L_Matrix Hi;
		if (p3dc.jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica(pa2, rF, Hi) == false)
			printf("Fallo en jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica()\n");
		p3dc.jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin_analitica(pa2, rF, Hi, J3);
		J3.print("J(analitica)");

		printf("Relaciones entre las matrices J(dif.fin)./J(analitica)\n");
		for (int i=0; i<J2.li; i++)
		{
			for (int j=0; j<J2.lj; j++)
				printf("%.2g ", J2(i,j) / J3(i,j));
			printf("\n");
		}

		printf("Parecido: J(itr.df)  J(analitica) %g\n", J2.cosineDistanceTo(J3));

		printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
		getchar();
	}
}




void L_Pose3D_cuat::pruebaPoseCuaternion_2()
{
	L_Pose3D po;
	L_HomogeneousMatrix H;
	L_Pose3D_cuat P3d, P3dM, P3dm;
	L_CoordsCart3D V, VM, Vm;
	L_CoordsCart3D P, PM, Pm;
	L_CoordsInv3D Vi, ViM, Vim;
	L_RayoCamara R, RM, Rm;
	double Pre36m[36];
	L_StaticMatrix<3,3> JJ3x3s;
	L_StaticMatrix<2,7> JJ2x7;
	L_StaticMatrix<2,3> JJ2x3;
	L_StaticMatrix<2,6> JJ2x6;
	L_StaticMatrix<3,4> JJ3x4;
	L_StaticMatrix<3,7> JJ3x7;
	L_Matrix j3x7(3,7), j2x3(2,3), j2x7(2,7), j2x6(2,6), j3x4(3,4), j3x3s(3,3);
	L_Matrix J3x7(3,7), J2x3(2,3), J2x7(2,7), J2x6(2,6), J3x4(3,4), J3x3s(3,3);
	int i, j;
	double delta = 1.0e-3;
	int npr = 1, nprmax = 10;

	// Proyeccion L_CoordsCart3D

	po.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	H.fijaPose3D(po); // movido a fijo, es lo normal
	P3d = H.calcPose3D_cuat(); // movido a fijo, es lo normal

	V.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	//
	printf("\nProbando L_Pose3D_cuat_proyectar_movido_a_fijo_v()\n");
	L_Pose3D_cuat_proyectar_movido_a_fijo_v(R,P3d,V);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*V;
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);
	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));

	//
	po.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	H.fijaPose3D(po); // movido a fijo, es lo normal
	P3d = H.calcPose3D_cuat(); // movido a fijo, es lo normal

	V.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	printf("\nProbando L_Pose3D_cuat_proyectar_movido_a_fijo_v()\n");
	L_Pose3D_cuat_proyectar_movido_a_fijo_v(R,P3d,V);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*V;
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);
	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));

	//
	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_pre()\n");
	L_Pose3D_cuat_jacob_movido_a_fijo_pre(P3d,JJ3x3s,Pre36m);
	JJ3x3s.print("JJ3x3s (P3d)");
	H.print("H");

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();

	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_v()\n");
	L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_v(P3d,V,JJ3x3s,Pre36m,JJ2x7,0,0,JJ2x3,0,0);
	JJ2x7.print("JJ2x7(P3d)");
	JJ2x3.print("JJ2x3(P3d)");
	L_StaticMatrix_OP_assign(j2x7, JJ2x7);
	L_StaticMatrix_OP_assign(j2x3, JJ2x3);

	for (j=0; j<7; j++)
	{
		P3d.el(j) += delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_v(RM,P3d,V);
		P3d.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_v(Rm,P3d,V);
		P3d.el(j) += delta;
		for (i=0; i<2; i++)
			J2x7(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}
	for (j=0; j<3; j++)
	{
		V.el(j) += delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_v(RM,P3d,V);
		V.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_v(Rm,P3d,V);
		V.el(j) += delta;
		for (i=0; i<2; i++)
			J2x3(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}

	J2x7.print("JJ2x7(diffin)");
	J2x3.print("JJ2x3(diffin)");

	printf("Parecido: %f %f\n", j2x3.cosineDistanceTo(J2x3), j2x7.cosineDistanceTo(J2x7));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();

	// Proyeccion L_CoordsInv3D

	printf("\nProbando L_Pose3D_cuat_proyectar_movido_a_fijo_vi()\n");

	Vi.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	L_Pose3D_cuat_proyectar_movido_a_fijo_vi(R,P3d,Vi);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*Vi.pideCoords3D();
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);

	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));

	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_vi()\n");
	L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_vi(P3d,Vi,JJ3x3s,Pre36m,JJ2x7,0,0,JJ2x6,0,0);
	JJ2x7.print("JJ2x7(P3d)");
	JJ2x6.print("JJ2x6(P3d)");
	L_StaticMatrix_OP_assign(j2x7, JJ2x7);
	L_StaticMatrix_OP_assign(j2x6, JJ2x6);

	for (j=0; j<7; j++)
	{
		P3dM = P3d;
		P3dm = P3d;
		P3dM.el(j) += delta;
		P3dm.el(j) -= delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(R,P3d,Vi);
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(RM,P3dM,Vi);
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(Rm,P3dm,Vi);
		for (i=0; i<2; i++)
			J2x7(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}
	for (j=0; j<6; j++)
	{
		ViM = Vi;
		Vim = Vi;
		ViM.el(j) += delta;
		Vim.el(j) -= delta;
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(R,P3d,Vi);
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(RM,P3d,ViM);
		L_Pose3D_cuat_proyectar_movido_a_fijo_vi(Rm,P3d,Vim);
		for (i=0; i<2; i++)
			J2x6(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}

	JJ2x7.print("JJ2x7(diffin)");
	JJ2x6.print("JJ2x6(diffin)");

	printf("Parecido: %f %f\n", j2x7.cosineDistanceTo(J2x7), j2x6.cosineDistanceTo(J2x6));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();



	printf("\n////////////////////////\n\n");

	printf("\nProbando L_Cuaternion_jacob_izq_rotarVectorInv()\n");

	P3d.ori.jacob_izq_rotarVectorInv_pre(Pre36m);
	L_Cuaternion_jacob_izq_rotarVectorInv(V,Pre36m,JJ3x4,0,0);

	JJ3x4.print("J3x4(pre)");
	L_StaticMatrix_OP_assign(J3x4, JJ3x4);

	P3d.ori.jacob_izq_rotarVectorInv(V,j3x4,0,0);

	j3x4.print("J3x4(dir)");
	
	printf("Parecido: %f\n", J3x4.cosineDistanceTo(j3x4));

	for (j=0; j<4; j++)
	{
		P3dM.ori = P3d.ori; P3dM.ori.el(j) += delta;
		P3dm.ori = P3d.ori; P3dm.ori.el(j) -= delta;
		L_Cuaternion_rotarVectorInv(VM,P3dM.ori,V);
		L_Cuaternion_rotarVectorInv(Vm,P3dm.ori,V);
		for (i=0; i<3; i++)
			J3x4(i,j) = (VM.el(i) - Vm.el(i)) / (2*delta);
	}

	J3x4.print("J3x4(dfn)");

	printf("Parecido: %f\n", J3x4.cosineDistanceTo(j3x4));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();

	printf("\nProbando L_Cuaternion_jacob_der_rotarVectorInv()\n");

	L_Cuaternion_jacob_der_rotarVectorInv(P3d.ori,JJ3x3s,0,0);

	JJ3x3s.print("J3x3s(def)");
	L_StaticMatrix_OP_assign(J3x3s, JJ3x3s);

	P3d.ori.jacob_der_rotarVectorInv(j3x3s,0,0);
	j3x3s.print("j3x3s(fn)");
	printf("Parecido: %f\n", J3x3s.cosineDistanceTo(j3x3s));

	for (j=0; j<3; j++)
	{
		VM = V; VM.el(j) += delta;
		Vm = V; Vm.el(j) -= delta;
		L_Cuaternion_rotarVectorInv(PM,P3d.ori,VM);
		L_Cuaternion_rotarVectorInv(Pm,P3d.ori,Vm);
		for (i=0; i<3; i++)
			J3x3s(i,j) = (PM.el(i) - Pm.el(i)) / (2*delta);
	}

	J3x3s.print("J3x3s(dfn)");
	printf("Parecido: %f\n", J3x3s.cosineDistanceTo(j3x3s));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();

	printf("Probando L_Pose3D_cuat_jacob_izq_fijo_a_movido()\n");
	
	L_Pose3D_cuat_jacob_fijo_a_movido_pre(P3d,JJ3x3s,Pre36m);
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(P3d,V,Pre36m,JJ3x3s,JJ3x7,0,0);

	JJ3x7.print("J3x7(pre)");
	L_StaticMatrix_OP_assign(j3x7, JJ3x7);

	for (j=0; j<7; j++)
	{
		P3dM = P3d; P3dM.el(j) += delta;
		P3dm = P3d; P3dm.el(j) -= delta;
		L_Pose3D_cuat_fijo_a_movido(PM,P3dM,V);
		L_Pose3D_cuat_fijo_a_movido(Pm,P3dm,V);
		for (i=0; i<3; i++)
			J3x7(i,j) = (PM.el(i) - Pm.el(i)) / (2*delta);
	}

	J3x7.print("J3x7(dfn)");

	printf("Parecido: %f\n", J3x7.cosineDistanceTo(j3x7));

	printf("Probando L_Pose3D_cuat_jacob_der_fijo_a_movido()\n");
	
	//L_Pose3D_cuat_jacob_fijo_a_movido_pre(P3d,JJ3x3s,Pre36m);
	L_Pose3D_cuat_jacob_der_fijo_a_movido(P3d,J3x3s,0,0);

	J3x3s.print("J3x3(#def)");

	for (j=0; j<3; j++)
	{
		V.el(j) += delta;
		L_Pose3D_cuat_fijo_a_movido(PM,P3d,V);
		V.el(j) -= 2*delta;
		L_Pose3D_cuat_fijo_a_movido(Pm,P3d,V);
		V.el(j) += delta;
		for (i=0; i<3; i++)
			j3x3s(i,j) = (PM.el(i) - Pm.el(i)) / (2*delta);
	}

	j3x3s.print("J3x3(dfn)");

	printf("Parecido: %f\n", J3x3s.cosineDistanceTo(j3x3s));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();

	printf("\n////////////////////////\n\n");

	// Proyeccion L_CoordsCart3D movido a fijo

	po.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	H.fijaPose3D_fijo_a_movido(po); // fijo a movido, es el uso inverso
	P3d = H.calcPose3D_cuat_fijo_a_movido(); // fijo a movido, es el uso inverso

	V.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	printf("\nProbando L_Pose3D_cuat_proyectar_fijo_a_movido_v()\n");
	L_Pose3D_cuat_proyectar_fijo_a_movido_v(R,P3d,V);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*V;
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);
	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));


	po.fijaPose(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-180, 180)*M_PI/180, L_RANDOM(-90, 90)*M_PI/180, L_RANDOM(-180, 180)*M_PI/180);
	H.fijaPose3D_fijo_a_movido(po); // fijo a movido, es el uso inverso
	P3d = H.calcPose3D_cuat_fijo_a_movido(); // fijo a movido, es el uso inverso

	V.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));

	printf("\nProbando L_Pose3D_cuat_proyectar_fijo_a_movido_v()\n");
	L_Pose3D_cuat_proyectar_fijo_a_movido_v(R,P3d,V);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*V;
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);
	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));


	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_pre()\n");
	L_Pose3D_cuat_jacob_fijo_a_movido_pre(P3d,JJ3x3s,Pre36m);
	JJ3x3s.print("JJ3x3s (P3d)");
	H.print("H");

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();


	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v()\n");

	L_Pose3D_cuat_jacob_fijo_a_movido_pre(P3d,J3x3s,Pre36m);
	L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v(P3d,V,J3x3s,Pre36m,j2x7,0,0,j2x3,0,0);

	j2x7.print("JJ2x7(#def)");
	j2x3.print("JJ2x3(#def)");

	for (j=0; j<7; j++)
	{
		P3d.el(j) += delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_v(RM,P3d,V);
		P3d.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_v(Rm,P3d,V);
		P3d.el(j) += delta;
		for (i=0; i<2; i++)
			J2x7(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}
	for (j=0; j<3; j++)
	{
		V.el(j) += delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_v(RM,P3d,V);
		V.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_v(Rm,P3d,V);
		V.el(j) += delta;
		for (i=0; i<2; i++)
			J2x3(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}

	J2x7.print("JJ2x7(diffin)");
	J2x3.print("JJ2x3(diffin)");

	printf("Parecido: %f %f\n", J2x7.cosineDistanceTo(j2x7), J2x3.cosineDistanceTo(j2x3));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();


	printf("\nProbando L_Pose3D_cuat_proyectar_fijo_a_movido_vi()\n");

	Vi.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10),L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	L_Pose3D_cuat_proyectar_fijo_a_movido_vi(R,P3d,Vi);
	printf(" >> R=(%.3f %.3f)   (P3d)\n", R.tanIzq, R.tanArr);

	VM = H*Vi.pideCoords3D();
	Rm.define(VM.y/VM.x,VM.z/VM.x);
	printf(" >> R=(%.3f %.3f)   (H)\n", Rm.tanIzq, Rm.tanArr);
	printf("Parecido: %f\n", R.cosineDistanceTo(Rm));

	printf("\nProbando L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_vi()\n");

	L_Pose3D_cuat_jacob_fijo_a_movido_pre(P3d,j3x3s,Pre36m);
	L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_vi(P3d,Vi,j3x3s,Pre36m,j2x7,0,0,j2x6,0,0);

	j2x7.print("JJ2x7(#def)");
	j2x6.print("JJ2x6(#def)");

	L_Pose3D_cuat_proyectar_fijo_a_movido_vi(R,P3d,Vi);
	for (j=0; j<7; j++)
	{
		P3d.el(j) += delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_vi(RM,P3d,Vi);
		P3d.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_vi(Rm,P3d,Vi);
		P3d.el(j) += delta;
		for (i=0; i<2; i++)
			J2x7(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}
	for (j=0; j<6; j++)
	{
		Vi.el(j) += delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_vi(RM,P3d,Vi);
		Vi.el(j) -= 2*delta;
		L_Pose3D_cuat_proyectar_fijo_a_movido_vi(Rm,P3d,Vi);
		Vi.el(j) += delta;
		for (i=0; i<2; i++)
			J2x6(i,j) = (RM.el(i)-Rm.el(i)) / (2*delta);
	}

	J2x7.print("JJ2x7(diffin)");
	J2x6.print("JJ2x6(diffin)");

	printf("Parecido: %f %f\n", J2x7.cosineDistanceTo(j2x7), J2x6.cosineDistanceTo(j2x6));

	printf("  >> presione ENTER para continuar (%d de %d) fin :)\n", npr++, nprmax);
	getchar();
}




void L_Pose3D_cuat::pruebaPoseCuaternion_3()
{
	L_Pose3D_cuat pCM, pC2_C, p_prc2, prM, prm, p3d_2, p3d_3;
	L_CoordsCart3D pun1, pun1M, pun1m, pun2, pun3;
	L_CoordsInv3D puni, puniM, punim;
	L_Matrix J3x3(3,3), j3x3(3,3), J3x7t(3,7), J7x7t(7,7), J3x7(3,7);
	L_Matrix J2x7(2,7), j2x7(2,7), JJ2x7(2,7), J2x3(2,3), j2x3(2,3), JJ2x3(2,3), J2x6(2,6), j2x6(2,6), JJ2x6(2,6);
	L_Matrix J3x7df(3,7);
	L_StaticMatrix<3,3> J3x3s, prJ3x3s, pc2J3x3s;
	L_RayoCamara r, rM, rm;
	double delta = 1.0e-5;
	double pre36m[36], pre36_pr[36], pre36m_pr[36];
	int i, j, npr=1, nprmax = 3;

	pCM.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));
	pCM.normaliza();
	pC2_C.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));
	pun1.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));

	printf("Probando funciones que involucran una camara secundaria\n\n");

	// Primero del mapa al robot, luego del robot a la 2a camara
	// (a*b)fm(c) = b.fm(a.fm(c))
	pun2 = pC2_C.fijo_a_movido(pCM.fijo_a_movido(pun1));
	pun2.print("pun2");

	p_prc2 = pCM*pC2_C;
	pun3 = p_prc2.fijo_a_movido(pun1);
	pun3.print("pun3");

	printf("Parecido: %f\n", pun2.cosineDistanceTo(pun3));

	// Derivadas respecto a pun1
	p_prc2.jacob_der_fijo_a_movido(J3x3);
	J3x3.print("jacob_der(fn)");

	for (j=0; j<3; j++)
	{
		pun1.el(j) += delta;
		pun2 = p_prc2.fijo_a_movido(pun1);
		pun1.el(j) -= 2*delta;
		pun3 = p_prc2.fijo_a_movido(pun1);
		pun1.el(j) += delta;

		for (i=0; i<3; i++)
			j3x3(i,j) = (pun2.el(i) - pun3.el(i)) / (2*delta);
	}
	j3x3.print("jacob_der(dif.fin)");

	printf("Parecido: %f\n", J3x3.cosineDistanceTo(j3x3));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();



	// Derivadas respecto a pCM
	// f(pCM) = (pCM*pC2_C).fijo_a_movido(pun1)
	//f(pCM) = fijo_movido(mult(pCM,pC2_C),pun1);
	//f(pCM) = fijo_movido(p_prc2,pun1);
	// df/dpr = fijo_movido(p_prc2,pun1)/dp_prc2 * mult(pCM,pC2_C)/dpr
	// J3x7   =      J3x7t                 *  J7x7t

	printf("Probando dfm(pCM*pc2,pun)/dpr\n");

	L_Pose3D_cuat_jacob_fijo_a_movido_pre(p_prc2,J3x3s,pre36m);
	L_Pose3D_cuat_jacob_izq_fijo_a_movido(p_prc2,pun1,pre36m,J3x3s,J3x7t,0,0); // Depende del punto
	L_Pose3D_cuat::jacob_izq_OP_mult(pCM,pC2_C,J7x7t); // Se reutiliza  // pC2_C,pCM
	J3x7.OP_mult(J3x7t,J7x7t);

	J3x7t.print("J3x7t");
	J7x7t.print("J7x7t");
	printf("J3x7 = J3x7t*J7x7t\n");
	J3x7.print("jacob_izq(deriv)");

	// Repetir todos

	for (j=0; j<7; j++)
	{
		pCM;
		pCM.el(j) += delta;
		pun2 = (pCM*pC2_C).fijo_a_movido(pun1);
		pCM.el(j) -= 2*delta;
		pun3 = (pCM*pC2_C).fijo_a_movido(pun1);
		pCM.el(j) += delta;
		for (i=0; i<3; i++)
			J3x7df(i,j) = (pun2.el(i) - pun3.el(i)) / (2*delta);
	}

	J3x7df.print("jacob_izq(diffn)");

	printf("Parecido: %f\n", J3x7.cosineDistanceTo(J3x7df));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();
	


	// Ahora calculos para la proyeccion

	// Derivadas respecto a pCM
	// f(pCM) = (pCM*pC2_C).pr_fijo_a_movido(pun1)
	//f(pCM) = pr_fijo_movido(mult(pCM,pC2_C),pun1);
	//f(pCM) = pr_fijo_movido(p_prc2,pun1);
	// df/dpCM = pr_fijo_movido(p_prc2,pun1)/dp_prc2 * mult(pCM,pC2_C)/dpCM
	// J2x7   =      J2x7t                 *  J7x7t

	printf("Probando dprfm(pCM*pc2,pun)/dpCM\n");
	p_prc2 = pCM*pC2_C;
	r = p_prc2.proyectar_fijo_a_movido_v(pun1);
	r.print("proy");

	p_prc2.jacob_proyectar_fijo_a_movido_pre(J3x3s, pre36m);
	L_Pose3D_cuat::jacob_izq_OP_mult_pre(pC2_C,pre36_pr); // Se compute UNA SOLA VEZ
	L_Cuaternion_calc_prem(pCM.ori,pre36m_pr);
	L_Cuaternion_jacob_der_rotarVectorInv(pCM.ori,prJ3x3s,0,0);
	L_Cuaternion_jacob_der_rotarVectorInv(pC2_C.ori,pc2J3x3s,0,0);

	p_prc2.jacob_proyectar_fijo_a_movido_v_c2(pun1, pCM, pC2_C, prJ3x3s, pc2J3x3s, J3x3s, pre36m_pr, j2x7, 0, 0, j2x3, 0, 0);
	p_prc2.jacob_proyectar_fijo_a_movido_v_c2_ant(pun1, pCM, pC2_C, J3x3s, pre36m, pre36_pr, J2x7, 0, 0, J2x3, 0, 0);
	p_prc2.jacob_proyectar_fijo_a_movido_v_c2_(pun1, pCM, pC2_C, JJ2x7, 0, 0, JJ2x3, 0, 0);

	j2x7.print("J2x7(prec)");
	J2x7.print("J2x7(prec)");
	JJ2x7.print("J2x7(dffn)");

	printf("Parecido: %f %f %f\n", j2x7.cosineDistanceTo(J2x7), j2x7.cosineDistanceTo(JJ2x7), J2x7.cosineDistanceTo(JJ2x7));

	j2x3.print("J2x3(prec)");
	J2x3.print("J2x3(prec)");
	JJ2x3.print("J2x3(dffn)");

	printf("Parecido: %f %f %f\n", j2x3.cosineDistanceTo(J2x3), j2x3.cosineDistanceTo(J2x3), J2x3.cosineDistanceTo(JJ2x3));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();



	puni.define(L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10),L_RANDOM(-10,10));

	printf("Probando dprfm(pCM*pc2,puni)/dpCM\n");
	p_prc2 = pCM*pC2_C;
	r = p_prc2.proyectar_fijo_a_movido_vi(puni);
	r.print("proy");

	// r = proy(movi(mu(pCM,pC2_C),puni))
	// dr/dpCM = dprmovi/dmu * dmu/dpCM
	p_prc2.jacob_proyectar_fijo_a_movido_pre(J3x3s, pre36m); // dprmovi/dmu (precalc)
	L_Pose3D_cuat::jacob_izq_OP_mult_pre(pC2_C,pre36_pr); // dmu/dpCM (precalc)

	// r = proy(movi(pC2_C,movi(pCM,puni)))
	// dr/dpCM = dproy/dmovi * dmovi(pC2_C)/dmovi(pCM) * dmovi(pCM)/dpCM
	L_Cuaternion_calc_prem(pCM.ori,pre36m_pr); //
	L_Cuaternion_jacob_der_rotarVectorInv(pCM.ori,prJ3x3s,0,0);
	L_Cuaternion_jacob_der_rotarVectorInv(pC2_C.ori,pc2J3x3s,0,0);

	p_prc2.jacob_proyectar_fijo_a_movido_vi_c2(puni, pCM, pC2_C, prJ3x3s, pc2J3x3s, J3x3s, pre36m_pr, j2x7, 0, 0, j2x6, 0, 0);
	p_prc2.jacob_proyectar_fijo_a_movido_vi_c2_ant(puni, pCM, pC2_C, J3x3s, pre36m, pre36_pr, J2x7, 0, 0, J2x6, 0, 0);
	// jprmovi_c2 es esta de abajo
	p_prc2.jacob_proyectar_fijo_a_movido_vi_c2_(puni, pCM, pC2_C, JJ2x7, 0, 0, JJ2x6, 0, 0);

	j2x7.print("J2x7(prec)");
	J2x7.print("J2x7(prec)");
	JJ2x7.print("J2x7(dffn)");

	printf("Parecido: %f %f %f\n", j2x7.cosineDistanceTo(J2x7), j2x7.cosineDistanceTo(JJ2x7), J2x7.cosineDistanceTo(JJ2x7));

	j2x6.print("J2x6(prec)");
	J2x6.print("J2x6(prec)");
	JJ2x6.print("J2x6(dffn)");

	printf("Parecido: %f %f %f\n", j2x6.cosineDistanceTo(J2x6), j2x6.cosineDistanceTo(JJ2x6), J2x6.cosineDistanceTo(JJ2x6));

	printf("  >> presione ENTER para continuar (%d de %d)\n", npr++, nprmax);
	getchar();
}

void L_Pose3D_cuat::pruebaPoseCuaternion_4()
{
	L_Pose3D_cuat p, pi, pim;
	L_Matrix J7x7, Ji, J7x7Ji;

	p.define(L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10), L_RANDOM(-10, 10));
	p.normaliza();

	printf("pi.inverseOf(p); pim.cambiaSigno(pi)\n");
	pi.inverseOf(p);
	pim.cambiaSigno(pi);

	p.print("p");
	pi.print("pi");
	pim.print("pim");

	getchar();

	pi.jacob_inverseOf(p, Ji);
	pim.jacob_cambiaSigno(pi, J7x7);
	
	printf("pi.jacob_inverseOf(p, Ji); pim.jacob_cambiaSigno(pi, J7x7)\n");

	Ji.print("Ji");
	J7x7.print("J7x7");
	J7x7Ji.OP_mult(J7x7, Ji); // Error: no la calculaba, usaba el producto directamente
	J7x7Ji.print("J7x7*Ji(mult)");

	pim.jacob_cambiaSigno_porM(pi, J7x7, Ji, J7x7Ji);
	J7x7Ji.print("J7x7*Ji(jacob)");

	getchar();
}


void L_Pose3D_cuat::pruebaPoseCuaternion_d()
{
	L_HomogeneousMatrix H1, H2, H3;
	L_Pose3D_cuat q1, q2, q3, q4;
	L_CoordsCart3D p1, p2, p3;
	L_CoordsInv3D e1;
	double pre36[36], pre36m[36], pre36mb[36];
	L_StaticMatrix<3,3> JJ3x3s;
	L_Matrix J3x3s(3,3), j3x3s(3,3), J3x7(3,7), j3x7(3,7), J3x3(3,3), j3x3(3,3);
	L_Matrix J2x7(2,7), j2x7(2,7), J2x3(2,3), j2x3(2,3), J2x6(2,6), j2x6(2,6);
	L_RayoCamara r1, r2;

	q1.fijaAzar();
	q2.fijaAzar();
	q3.fijaAzar();
	q4.fijaAzar();
	p1.fijaAzar();
	p2.fijaAzar();
	p3.fijaAzar();

	L_Pose3D_cuat_movido_a_fijo(p2,q1,p1);
	p3 = q1.movido_a_fijo(p1);
	p2.print("mov_a_fijo(#def)");
	p3.print("mov_a_fijo(func)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));

	L_Pose3D_cuat_movido_a_fijo(p2,q1,p1);
	q1.mov(p3,p1);
	p2.print("mov_a_fijo(#def)");
	p3.print("mov_a_fijo(mov)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));

	L_Pose3D_cuat_fijo_a_movido(p2,q1,p1);
	p3 = q1.fijo_a_movido(p1);
	p2.print("mov_a_fijo(#def)");
	p3.print("mov_a_fijo(func)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));

	L_Pose3D_cuat_fijo_a_movido(p2,q1,p1);
	q1.movi(p3,p1);
	p2.print("mov_a_fijo(#def)");
	p3.print("mov_a_fijo(movi)");
	printf("Parecido: %f\n", p2.cosineDistanceTo(p3));

	printf("Presione ENTER para continuar\n");
	getchar();

	L_Pose3D_cuat_inverseOf(q2,q1);
	q3.inverseOf(q1);
	q2.print("invDe(#def)");
	q3.print("invDe(func)");
	printf("Parecido: %f\n", q2.cosineDistanceTo(q3));

	q1.fijaAzar();
	q2.fijaAzar();

	L_Pose3D_cuat_OP_mult(q3,q1,q2);
	q4.OP_mult(q1,q2);
	q3.print("OP_mult(#def)");
	q4.print("OP_mult(func)");
	printf("Parecido: %f\n", q3.cosineDistanceTo(q4));
	
	printf("Presione ENTER para continuar\n");
	getchar();

	L_Pose3D_cuat_OP_mult_p1_por_p2inv(q3,q1,q2);
	q4.OP_mult_p1_por_p2inv(q1,q2);
	q3.print("q1*q2^-1(#def)");
	q4.print("q1*q2^-1(func)");
	printf("Parecido: %f\n", q3.cosineDistanceTo(q4));

	L_Pose3D_cuat_jacob_fijo_a_movido_pre(q1,J3x3s,pre36m);
	q1.jacob_proyectar_fijo_a_movido_pre(JJ3x3s, pre36mb);
	L_StaticMatrix_OP_assign(j3x3s, JJ3x3s);
	L_Matrix::imprime_arreglo(pre36m, 36, "pre36m(#def)");
	L_Matrix::imprime_arreglo(pre36mb, 36, "pre36m(func)");
	printf("Parecido: %f\n", L_Matrix::cosineDistanceBetweenArrays(pre36m, pre36mb, 36));
	J3x3s.print("J3x3s(#def)");
	j3x3s.print("J3x3s(#def)");
	printf("Parecido: %f\n", J3x3s.cosineDistanceTo(j3x3s));

	L_Pose3D_cuat_jacob_movido_a_fijo_pre(q1,J3x3,pre36);
	L_Matrix::imprime_arreglo(pre36, 36, "pre36(#def)");
	J3x3.print("J3x3d(#def)");

	printf("Presione ENTER para continuar\n");
	getchar();

	//  Funciones inutiles para el SLAM
	//L_Pose3D_cuat_jacob_izq_fijo_a_movido(...); -> no tiene equivalente
	L_Pose3D_cuat_jacob_izq_movido_a_fijo(q1,p1,pre36,J3x7);
	q1.jacob_izq_movido_a_fijo(p1,j3x7);
	J3x7.print("ji_mf(#def)");
	j3x7.print("ji_mf(func)");
	printf("Parecido: %f\n", J3x7.cosineDistanceTo(j3x7));
	
	L_Pose3D_cuat_jacob_der_movido_a_fijo(q1, J3x3);
	q1.jacob_der_movido_a_fijo(j3x3);
	J3x3.print("jd_mf(#def)");
	j3x3.print("jd_mf(func)");
	printf("Parecido: %f\n", J3x3.cosineDistanceTo(j3x3));

	printf("Presione ENTER para continuar\n");
	getchar();

	double xs;
	L_Pose3D_cuat_proyectar_fijo_a_movido_v(r1,q1,p1);
	r2 = q1.proyectar_fijo_a_movido_v(p1);
	r1.print("proyfm(#def)");
	r2.print("proyfm(func)");
	printf("Parecido: %f\n", r1.cosineDistanceTo(r2));
	q1.prmovi_v(r2,p1,xs);
	r1.print("proyfm(#def)");
	r2.print("proyfm(prmovi)");
	printf("Parecido: %f\n", r1.cosineDistanceTo(r2));

	e1.fijaAzar();
	L_Pose3D_cuat_proyectar_fijo_a_movido_vi(r1,q1,e1);
	r2 = q1.proyectar_fijo_a_movido_vi(e1);
	r1.print("proyfm6(#def)");
	r2.print("proyfm6(func)");
	printf("Parecido: %f\n", r1.cosineDistanceTo(r2));
	q1.prmovi_vi(r2,e1,xs);
	r1.print("proyfm6(#def)");
	r2.print("proyfm6(prmovi)");
	printf("Parecido: %f\n", r1.cosineDistanceTo(r2));

	printf("Presione ENTER para continuar\n");
	getchar();

	L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_v(q1,p1,JJ3x3s,pre36m,J2x7,0,0,J2x3,0,0);
	q1.jacob_proyectar_fijo_a_movido_v(p1,JJ3x3s,pre36m,j2x7,0,0,j2x3,0,0);

	J2x7.print("j1_pfm(#def)");
	j2x7.print("j1_pfm(func)");
	printf("Parecido:%f\n", J2x7.cosineDistanceTo(j2x7));
	J2x3.print("j2_pfm(#def)");
	j2x3.print("j2_pfm(func)");
	printf("Parecido:%f\n", J2x3.cosineDistanceTo(j2x3));

	q1.jprmovi_v(p1,JJ3x3s,pre36m,j2x7,0,0,j2x3,0,0);

	J2x7.print("j1_pfm(#def)");
	j2x7.print("j1_pfm(jprmovi)");
	printf("Parecido:%f\n", J2x7.cosineDistanceTo(j2x7));
	J2x3.print("j2_pfm(#def)");
	j2x3.print("j2_pfm(jprmovi)");
	printf("Parecido:%f\n", J2x3.cosineDistanceTo(j2x3));

	printf("Presione ENTER para continuar\n");
	getchar();

	L_Pose3D_cuat_jacob_proyectar_fijo_a_movido_vi(q1,e1,JJ3x3s,pre36m,J2x7,0,0,J2x6,0,0);
	q1.jacob_proyectar_fijo_a_movido_vi(e1,JJ3x3s,pre36m,j2x7,0,0,j2x6,0,0);
	J2x7.print("j1_pfm(#def)");
	j2x7.print("j1_pfm(func)");
	printf("Parecido:%f\n", J2x7.cosineDistanceTo(j2x7));
	J2x6.print("j2_pfm(#def)");
	j2x6.print("j2_pfm(func)");
	printf("Parecido:%f\n", J2x6.cosineDistanceTo(j2x6));

	q1.jprmovi_vi(e1,JJ3x3s,pre36m,j2x7,0,0,j2x6,0,0);
	J2x7.print("j1_pfm(#def)");
	j2x7.print("j1_pfm(jprmovi)");
	printf("Parecido:%f\n", J2x7.cosineDistanceTo(j2x7));
	J2x6.print("j2_pfm(#def)");
	j2x6.print("j2_pfm(jprmovi)");
	printf("Parecido:%f\n", J2x6.cosineDistanceTo(j2x6));

	printf("Presione ENTER para continuar\n");
	getchar();

	// Ahora vienen funciones inutiles para el SLAM
	L_Pose3D_cuat_proyectar_movido_a_fijo_v(r1,q1,p1);
	r2 = q1.proyectar_movido_a_fijo(p1);
	r1.print("proymf(#def)");
	r2.print("proymf(func)");
	printf("Parecido: %f (inutiles)\n", r1.cosineDistanceTo(r2));

	// Ahora vienen funciones inutiles para el SLAM
	e1.fijaAzar();
	L_Pose3D_cuat_proyectar_movido_a_fijo_vi(r1,q1,e1);
	r2 = q1.proyectar_movido_a_fijo(e1);
	r1.print("proymf6(#def)");
	r2.print("proymf6(func)");
	printf("Parecido: %f (inutiles)\n", r1.cosineDistanceTo(r2));

	// Ahora vienen funciones inutiles para el SLAM
	L_Pose3D_cuat_jacob_movido_a_fijo_pre(q1,J3x3,pre36);
	L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_v(q1,p1,J3x3,pre36,J2x7,0,0,J2x3,0,0);
	q1.jacob_izq_proyectar_movido_a_fijo(p1,j2x7, 0, 0);
	q1.jacob_der_proyectar_movido_a_fijo(p1,j2x3,0,0);

	J2x7.print("j1_pmf(#def)");
	j2x7.print("j1_pmf(func)");
	printf("Parecido:%f   (inutiles)\n", J2x7.cosineDistanceTo(j2x7));
	J2x3.print("j2_pmf(#def)");
	j2x3.print("j2_pmf(func)");
	printf("Parecido:%f   (inutiles)\n", J2x3.cosineDistanceTo(j2x3));

	printf("Presione ENTER para continuar\n");

	getchar();

	// Ahora vienen funciones inutiles para el SLAM
	L_Pose3D_cuat_jacob_movido_a_fijo_pre(q1,J3x3,pre36);
	L_Pose3D_cuat_jacob_proyectar_movido_a_fijo_vi(q1,e1,J3x3,pre36,J2x7,0,0,J2x6,0,0);
	q1.jacob_izq_proyectar_movido_a_fijo(e1,j2x7,0,0);
	q1.jacob_der_proyectar_movido_a_fijo(e1,j2x6,0,0);
	J2x7.print("j1_pmf(#def)");
	j2x7.print("j1_pmf(func)");
	printf("Parecido:%f   (inutiles)\n", J2x7.cosineDistanceTo(j2x7));
	J2x6.print("j2_pmf(#def)");
	j2x6.print("j2_pmf(func)");
	printf("Parecido:%f   (inutiles)\n", J2x6.cosineDistanceTo(j2x6));


	printf("Presione ENTER para continuar\n");
	getchar();

}

// Funcion chanta
void L_Pose3D_cuat::pruebaPoseCuaternion_converg()
{
	std::vector<L_CoordsCart3D> pos;
	std::vector<L_RayoCamara> r;
	L_RayoCamara rtmp;
	L_Pose3D_cuat pEstim, p3d, p3dM, p3dm, pReal;
	double xs, delta = 1.0e-6;
	L_Matrix dx, z, h, hM, hm, H, HT, A, K, e;
	int i, j, iter;

	printf("test de convergencia de puntocuaterniones\n\n");

	pReal.pos.fijaCero();
	pReal.ori.setRotationComponents(160*M_PI/180, 0, 0, 1);
	//pReal.ori.setRotationComponents(1*M_PI/180, 0, 0, 1);

	pEstim.pos.fijaCero();
	pEstim.ori.setRotationComponents(180*M_PI/180, 0, 0, 1);
	//pEstim.ori.setRotationComponents(1.1*M_PI/180, 0, 0, 1);

	int metodo = 1; // 0 o 1

	//
	printf("Calculo de pose por minimizacion LM de proyecciones\n");

	int nMax = 20;
	pos.resize(nMax);
	r.resize(nMax);

	for (i=0; i<nMax; i++)
	{
		pos[i].x = L_RANDOM(-100, -0.1);
		pos[i].y = L_RANDOM(-1,1)*pos[i].x;
		pos[i].z = L_RANDOM(-1,1)*pos[i].x;
		pReal.prmovi_v(r[i], pos[i], xs);
	}

	p3d = pEstim;
	p3d._calcPose_movido_a_fijo_Nrayos(r, pos, 10, &p3d);
	pReal.print("pReal");
	p3d.print("pEsti");

	printf("\n");

	//
	printf("Calculo de pose por metodo del gradiente\n");

	p3d = pEstim;

	dx.reallocate(7, 1);
	z.reallocate(nMax*2, 1);
	h.reallocate(nMax*2, 1);
	hM.reallocate(nMax*2, 1);
	hm.reallocate(nMax*2, 1);
	H.reallocate(nMax*2, 7);

	// z es invariable
	for (i=0; i<nMax; i++)
	{
		z(2*i+0,0) = r[i].tanIzq;
		z(2*i+1,0) = r[i].tanArr;
	}

	for (iter=0; iter<50; iter++)
	{
		// h varia de iteracion en iteracion
		for (i=0; i<nMax; i++)
		{
			p3d.prmovi_v(rtmp, pos[i], xs);
			h(2*i+0,0) = rtmp.tanIzq;
			h(2*i+1,0) = rtmp.tanArr;
		}
		for (j=0; j<7; j++)
		{
			p3dM = p3d;
			p3dM.el(j) += delta;
			for (i=0; i<nMax; i++)
			{
				p3dM.prmovi_v(rtmp, pos[i], xs);
				hM(2*i+0,0) = rtmp.tanIzq;
				hM(2*i+1,0) = rtmp.tanArr;
			}
			p3dm = p3d;
			p3dm.el(j) -= delta;
			for (i=0; i<nMax; i++)
			{
				p3dm.prmovi_v(rtmp, pos[i], xs);
				hm(2*i+0,0) = rtmp.tanIzq;
				hm(2*i+1,0) = rtmp.tanArr;
			}
			for (i=0; i<2*nMax; i++)
				H(i,j) = (hM(i,j) - hm(i,j)) / (2*delta);
		}

		// z = f(xRe)
		// h = f(xEs)
		// f(xRe + err) = f(xRe) + dfdx*err
		// h = z + H*err
		// err = inv(HT*H)*HT*(h-z);

		e.OP_subtract(z, h);

		if (metodo == 0) // Tipo gradiente
		{
			// xEs := xEs + inv(HT*H)*HT*(z-h);
			HT.transpOf(H);
			A.OP_mult(HT,H);
			L_Matrix id(A.li,A.lj); id.identity(1e-5); A.OP_add(id); // Normalizacion: K no es invertible
			if (A.invertMe() == false)
				printf("HTH no invertible\n");
			K.OP_mult(A,HT);
			printf("Tr(A)=%g  Det(A)=%g  Det(K)=%g\n", A.trace(), A.det(), K.det());
		}
		else // Tipo kalman
		{
			// xEs := xEs + P*HT*inv(H*P*HT+R)*(z-h);
			HT.transpOf(H);
			A.OP_mult(H,HT);
			L_Matrix id(A.li,A.lj); id.identity(1e-5); A.OP_add(id);
			if (A.invertMe() == false)
				printf("HHT no invertible\n");
			K.OP_mult(HT,A);
			printf("Tr(A)=%g  Det(A)=%g  Det(K)=%g\n", A.trace(), A.det(), K.det());
		}

		// De todos modos K es parecido a HT

		dx.OP_mult(K, e);
		for (i=0; i<7; i++)
			p3d.el(i) += dx(i,0);
		p3d.normaliza();

	}
	pReal.print("pReal");
	p3d.print("pEsti");

	getchar();
}


void L_Pose3D_cuat::pruebaPoseCuaternion_veloc()
{
	L_Pose3D_cuat pose;
	int i, n=10, nf;
	double t1;
	std::vector<L_CoordsCart3D> arrP(n);
	std::vector<L_RayoCamara> arrR(n);
	L_Matrix Hi, J7x2n(7, 2*n), J7x3n(7, 3*n);

	pose.fijaAzar();
	for (i=0; i<n; i++)
		arrP[i].fijaAzar();
	for (i=0; i<n; i++)
	{
		L_CoordsCart3D pun;
		pose.movi(pun, arrP[i]);
		arrR[i].tanIzq = pun.y/pun.x + L_RANDOM(-0.05, 0.05);
		arrR[i].tanArr = pun.z/pun.x + L_RANDOM(-0.05, 0.05);
	}

	nf = 0;
	t1 = L_TIME();
	while (L_TIME() - t1 < 0.5)
	{
		if (pose.jacob_pre_calcPose_movido_a_fijo_Nrayos_diffin_analitica(arrP, arrR, Hi, 1.0e-4) == false)
			printf("false\n");
		nf++;
	}
	printf("Hi: %f fps\n", nf/0.5);

	nf = 0;
	t1 = L_TIME();
	while (L_TIME() - t1 < 0.5)
	{
		pose.jacob_rf_calcPose_movido_a_fijo_Nrayos_diffin_analitica(arrP, arrR, Hi, J7x2n);
		nf++;
	}
	printf("J7x2n: %f fps\n", nf/0.5);

	nf = 0;
	t1 = L_TIME();
	while (L_TIME() - t1 < 0.5)
	{
		pose.jacob_pm_calcPose_movido_a_fijo_Nrayos_diffin_analitica(arrP, arrR, Hi, J7x3n);
		nf++;
	}
	printf("J7x3n: %f fps\n", nf/0.5);

	getchar();

}

void L_Pose3DSpline::agregaPose(double t, const L_Pose3D_cuat &offset, double x, double y, double z, double pan, double tilt, double roll)
{
	L_Pose3D_cuat poseRes, poseFin;
	L_Pose3D_cuat_sumable pSum;
	poseRes.definePosPanTiltRoll(x, y, z, pan, tilt, roll);
	poseFin = offset * poseRes;  // pOM * PCO
	if (this->y.size() > 0)
	{
		double distCos = this->y[this->y.size()-1].ori.OP_mult_elementwise(poseFin.ori);
		if (distCos < 0)
			poseFin.ori = -poseFin.ori;
	}
	pSum.L_Pose3D_cuat::operator=(poseFin);
	addFrom(t, pSum);
}


void L_Pose3DSpline::agregaPoseMirarA(double t, const L_Pose3D_cuat &offset, double x, double y, double z, double xMirado, double yMirado, double zMirado, double roll)
{
	L_Pose3D_cuat poseRes, poseFin;
	L_Pose3D_cuat_sumable pSum;
	poseRes.fijaPoseMirarA(x,y,z,xMirado,yMirado,zMirado,roll);
	poseFin = offset * poseRes;  // pOM * PCO
	if (this->y.size() > 0)
	{
		double distCos = this->y[this->y.size()-1].ori.OP_mult_elementwise(poseFin.ori);
		if (distCos < 0)
			poseFin.ori = -poseFin.ori;
	}
	pSum.L_Pose3D_cuat::operator=(poseFin);
	addFrom(t, pSum);
}

bool L_Pose3DSpline::leerPuntosControl(const char *nomArch)
{
	FILE *fp=NULL;
	L_Array<L_Pose3D> p3d;
	std::vector<double> t;
	L_Pose3D_cuat offset;
	int i, j, ene;
	int retr;

	offset.fijaCero();

	if (nomArch == NULL)
	{
		printf("L_Pose3DSpline::leerPuntosControl() : name de archivo no especificado\n");
		return false;
	}

	fp = fopen(nomArch, "r");
	if (fp == NULL)
	{
		printf("No se pudo abrir el archivo de puntos de control %s\n", nomArch);
		return false;
	}
	int nleidos = 0;
	nleidos += fscanf(fp, "%d", &ene);
	nleidos += fscanf(fp, "%d", &retr);
	nleidos += fscanf(fp, "%lg", &t1);
	nleidos += fscanf(fp, "%lg", &t2);
	if (nleidos != 4)
	{
		printf("L_Pose3DSpline::leerPuntosControl() : formato de archivo incorrecto\n");
		fclose(fp);
		return false;
	}

	retroceder_t = (retr != 0);

	x.reserve(ene);
	y.reserve(ene);
	y2.reserve(ene);
	x.resize(0);
	y.resize(0);
	y2.resize(0);
	p3d.resize(ene);
	t.resize(ene);

	for (i=0; i<6; i++)
	{
		for (j=0; j<ene; j++)
		{
			if (fscanf(fp, "%lg", &p3d[j].el(i)) != 1)
			{
				printf("Error en formato de archivo %s\n", nomArch);
				fclose(fp);
				return false;
			}
		}
	}
	for (j=0; j<ene; j++)
	{
		if (fscanf(fp, "%lg", &t[j]) != 1)
		{
			printf("Error en formato de archivo %s\n", nomArch);
			fclose(fp);
			return false;
		}
	}
	fclose(fp);
	fp = NULL;

	for (j=0; j<ene; j++)
		agregaPose(t[j], offset, p3d[j].pos.x, p3d[j].pos.y, p3d[j].pos.z, p3d[j].ori.pan, p3d[j].ori.tilt, p3d[j].ori.roll);
	compute();
	return true;
}

bool L_Pose3DSpline::grabarPuntosControl(const char *nomArch)
{
	FILE *fp=NULL;
	std::vector<L_Pose3D> p3d;
	int i, j;

	if (nomArch == NULL)
	{
		printf("L_Pose3DSpline::leerPuntosControl() : name de archivo no especificado\n");
		return false;
	}

	fp = fopen(nomArch, "w");
	if (fp == NULL)
	{
		printf("No se pudo abrir el archivo de puntos de control %s\n", nomArch);
		return false;
	}
	if (retroceder_t == false)
	{
		t1 = 0; // Para evitar NaN
		t2 = 0;
	}
	fprintf(fp, "%d %d %f %f\n", n, retroceder_t?1:0, t1, t2);
	p3d.resize(y.size());

	for (j=0; j<(int)y.size(); j++)
		p3d[j] = y[j].calcPose3D();

	for (i=0; i<6; i++)
	{
		for (j=0; j<(int)y.size(); j++)
		{
			fprintf(fp, "%.5g  ", p3d[j].el(i));
		}
		fprintf(fp, "\n");
	}
	for (j=0; j<(int)x.size(); j++)
		fprintf(fp, "%g  ", x[j]);
	fprintf(fp, "\n");
	fclose(fp);
	return true;
}


void L_Pose3DSpline::generarRutaEjemplo(int numRuta, L_Rand &randF)
{
	L_Pose3D_cuat poseOffset;
	L_CoordsCart3D puntoMirado;
	L_CoordsCart3D ultPos, deltaPos;
	double tFactor;
	double rad = 40;
	int i;
	// dt = 1.0/15 en prueba_v6_simulada()
	// solve({1.0*i = frames*1.0/15, frames = 1200}, {frames, i}); -> i = 80

	// 0 : forma luna
	// 1 : forma ese
	// 2 : ese + teleport cuadrado

	x.resize(0);
	y.resize(0);
	y2.resize(0);
	n = 0;
	retroceder_t = false;
	t1 = -1;
	t2 = -1;
	switch (numRuta)
	{
	case 0:
		// Ruta repetitiva en forma de U 
		tFactor = 3.0;
		poseOffset.definePosPanTiltRoll(-60, 0, 0, 0, 0, 0);
		agregaPose(0*tFactor, poseOffset,  0,  0, 0,  0, 0, 0);
		agregaPose(0.05*tFactor, poseOffset,  0,  0, 0,  0, 0, 0);
		agregaPose(0.1*tFactor, poseOffset,  0,  0, 0,  0, 0, 0);
		agregaPose(2*tFactor, poseOffset, 30, 60, 0, -30*M_PI/180,0,0);
		agregaPose(3*tFactor, poseOffset, 70, 90, 0, -60*M_PI/180,0,0);
		agregaPose(4*tFactor, poseOffset, 30, 60, 0, -30*M_PI/180,0,0);
		agregaPose(5*tFactor, poseOffset, 00,  0, 0,   0,0,0);
		agregaPose(6*tFactor, poseOffset, 30,-60, 0,  30*M_PI/180,0,0);
		agregaPose(7*tFactor, poseOffset, 70,-90, 0,  60*M_PI/180,0,0);
		agregaPose(8*tFactor, poseOffset, 30,-60, 0,  30*M_PI/180,0,0);
		agregaPose(9*tFactor, poseOffset, 00,  0, 0,  0,0,0);
		agregaPose(10*tFactor, poseOffset, 30, 60, 0, -30*M_PI/180,0,0); // El 10 -> 2
		agregaPose(11*tFactor, poseOffset, 70, 90, 0, -60*M_PI/180,0,0);
		fijarSaltoTiempoPeriodico(10*tFactor, 2*tFactor); // Ahora el 10 vuelve al 2
		compute(); // build trayectoria usando spline
		break;
	case 1:
		// Ruta con forma de S, se aleja progresivamente del punto de mira
		tFactor = 1.0;
		poseOffset.definePosPanTiltRoll(0, 0, 0, 0, 0, 0);
		puntoMirado.define(30, 0, 0);
		i=0;
		agregaPoseMirarA(0*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.05*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.1*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		for (i=1; i<300; i++) // Hasta el tiempo t*300
			agregaPoseMirarA(i*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		compute(); // build trayectoria usando spline
			break;
	case 2:
		// Ruta con forma de S, luego teleportaciones entre las 4 esquinas de un rectangulo
		tFactor = 1.0;
		poseOffset.definePosPanTiltRoll(0, 0, 0, 0, 0, 0);
		puntoMirado.define(30, 0, 0);
		i=0;
		agregaPoseMirarA(0*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.05*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.1*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		for (i=1; i<80; i++)
			agregaPoseMirarA(i*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		for (; i<140; i++)
		{
			switch((i%16)/4) // Teleportaciones
			{
			case 0:
				agregaPoseMirarA(i*tFactor, poseOffset, -90, -40, 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
				break;
			case 1:
				agregaPoseMirarA(i*tFactor, poseOffset, -90, 40, 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
				break;
			case 2:
				agregaPoseMirarA(i*tFactor, poseOffset, -30, 40, 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
				break;
			case 3:
				agregaPoseMirarA(i*tFactor, poseOffset, -30, -40, 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
				break;
			}
		}
		fijarSaltoTiempoPeriodico(tFactor*128, tFactor*112);
		compute(); // build trayectoria usando spline
			break;
	case 3:
		// Ruta en forma de S, luego poses al azar
		tFactor = 1.0;
		poseOffset.definePosPanTiltRoll(0, 0, 0, 0, 0, 0);
		puntoMirado.define(30, 0, 0);
		i=0;
		agregaPoseMirarA(0*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.05*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		agregaPoseMirarA(0.1*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		for (i=1; i<80; i++)
			agregaPoseMirarA(i*tFactor, poseOffset, -rad*(1+i/54.0)*cos(cos(1.0/3*i)), rad*(1+i/54.0)*sin(cos(1.0/3*i)), 0,   puntoMirado.x,puntoMirado.y,puntoMirado.z);
		compute();
		ultPos.define(y[y.size()-1].pos.x, y[y.size()-1].pos.y, y[y.size()-1].pos.z);
		for (; i<300; i++)
		{
			deltaPos.x = randF.random(-10, 10);
			deltaPos.y = randF.random(-10, 10);
			deltaPos.z = randF.random(-10, 10);
			ultPos += deltaPos;
			agregaPoseMirarA(i*tFactor, poseOffset, ultPos.x, ultPos.y, ultPos.z, puntoMirado.x, puntoMirado.y, puntoMirado.z);
		}
		break;
	default:
		throw_L_ArgException_if(true, "L_Pose3DSpline::generarRutaEjemplo() : numero de ruta inexistente");
	}
}

bool L_Pose3DSpline::generarGrabarRutasAzar(L_StringWithNumber &nomArch, int nRutas, L_Rand &randF)
{
	int i;
	for (i=0; i<nRutas; i++)
	{
		generarRutaEjemplo(3, randF);
		if (grabarPuntosControl(nomArch.generateString(i)) == false)
			return false;
	}
	return true;
}

void L_Pose3DSpline::genDibujoFlaite(L_ImageRGBUchar &im)
{
	int i;
	double minIzq, minArr, maxIzq, maxArr;
	double u, v, a, b, xVar;
	double deltamin = 1e10;
	L_Pose3D p3d;
	L_Array<L_Pose3D_cuat> other;
	L_Pose3D_cuat aca;
	L_ShapeArray lins;
	if (im.data() == NULL)
		im.reallocate(640, 480);
	
	im.setZero();

	if (x.size() > 0)
	{
		minIzq = maxIzq = y[0].pos.x;
		minArr = maxArr = y[0].pos.y;
	}
	else
	{
		minIzq = 0;
		minArr = 0;
		maxIzq = 1;
		maxArr = 1;
	}
	for (i=1; i<(int)x.size(); i++)
	{
		if (y[i].pos.y < minIzq)
			minIzq = y[i].pos.y;
		if (y[i].pos.y > maxIzq)
			maxIzq = y[i].pos.y;
		if (y[i].pos.x < minArr)
			minArr = y[i].pos.x;
		if (y[i].pos.x > maxArr)
			maxArr = y[i].pos.x;
	}

	if (minIzq > 0)
		minIzq = 0;
	if (maxIzq < 0)
		maxIzq = 0;
	if (minArr > 0)
		minArr = 0;
	if (maxArr < 0)
		maxArr = 0;

	double mx, nx, my, ny;

	// Los sis ref estan invertidos
	// Los pixeles se definen centrados en zero
	// u = minIzq -> pix = 0.45*L
	// u = maxIzq -> pix = -0.45*L
	// pix = m * u + n
	// solve({0.45*L = m * minIzq + n, -0.45*L = m * maxIzq + n},{m,n});
	
	mx = .9000000000*im.lx/(minIzq-maxIzq);
	nx = -.4500000000*im.lx*(minIzq+maxIzq)/(minIzq-maxIzq);

	my = .9000000000*im.ly/(minArr-maxArr);
	ny = -.4500000000*im.ly*(minArr+maxArr)/(minArr-maxArr);

	// Seleccionar el menor m
	mx = L_MAX(mx, my);
	my = mx;

	// Agregar tx/2 a nx y a ny
	nx += im.lx/2;
	ny += im.ly/2;

	// Graficar los puntos

	// Lineas de demarcacion del sistema de referencia

	u = -100*mx + nx;
	v = 0*my + ny;
	a = 100*mx + nx;
	b = 0*my + ny;
	lins.drawLine((int)u, (int)v, (int)a, (int)b, 255, 0, 0);

	u = 0*mx + nx;
	v = -100*my + ny;
	a = 0*mx + nx;
	b = 100*my + ny;
	lins.drawLine((int)u, (int)v, (int)a, (int)b, 255, 0, 0);

	// Dibujar puntos control
	for (i=1; i<(int)x.size(); i++)
	{
		u = y[i-1].pos.y*mx + nx;
		v = y[i-1].pos.x*my + ny;
		a = y[i].pos.y*mx + nx;
		b = y[i].pos.x*my + ny;

		lins.drawLine((int)u, (int)v, (int)a, (int)b, 255, 255, 255);
		if (i%4 == 0)
		{
			p3d = y[i].calcPose3D();
			lins.drawArrow((int)a, (int)b, 10, p3d.ori.pan+M_PI/2, 0, 0, 255);
		}
	}
	if (retroceder_t == true)
	{
		// Dibujar los dos extremos con un circulo rojo
		aca = evaluate(t1);
		u = aca.pos.y*mx + nx;
		v = aca.pos.x*my + ny;

		aca = evaluate(t2);
		a = aca.pos.y*mx + nx;
		b = aca.pos.x*my + ny;

		lins.drawCircle((int)u, (int) v, 10, 255, 0, 0);
		lins.drawCircle((int)u, (int) v, 8, 255, 0, 0);
		lins.drawCircle((int)a, (int) b, 10, 255, 0, 0);
		lins.drawCircle((int)a, (int) b, 8, 255, 0, 0);
	}

	// Dibujar spline
	deltamin = x[x.size()-1] / 100000;
	aca = evaluate(0);
	for (xVar = 0; xVar < x[x.size()-1]; xVar += deltamin)
	{
		u = aca.pos.y*mx + nx;
		v = aca.pos.x*my + ny;
		aca = evaluate(xVar);
		a = aca.pos.y*mx + nx;
		b = aca.pos.x*my + ny;
		lins.drawLine((int)u, (int)v, (int)a, (int)b, 255, 0, 0);

		//if (a<0 || a>im.lx || b<0 || b>im.ly)
		//	printf("Aca es\n");
	}

	lins.drawNumber(20, 10, (int)(x[x.size()-1]), 10, 255, 255, 255);
	lins.drawNumber(20, 30, (int)x.size(), 10, 255, 255, 255);
	lins.drawArrow(10, 21, 12, 45, 0, 0, 255);
	lins.drawCircle(10, 35, 8, 0, 0, 255);
	lins.drawLine(10, 35, 10+8, 35, 255, 0, 0);
	lins.drawLine(10, 35, 10, 35-8, 255, 0, 0);
	im.genDrawing(lins);
}


L_Rot3DMatrix &L_Rot3DMatrix::fijarCuaternionInv(const L_Quaternion &q)
{
	double a2, b2, c2, d2, ab, ac, ad, bc, bd, cd;
	double fac = q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d;
	if (fac == 0)
		printf("L_Rot3DMatrix::fijarCuaternionInv() : cuaternion con modulo zero\n");
	fac = 1/fac;
	a2 = q.a*q.a*fac;
	b2 = q.b*q.b*fac;
	c2 = q.c*q.c*fac;
	d2 = q.d*q.d*fac;
	ab = q.a*q.b*fac;
	ac = q.a*q.c*fac;
	ad = q.a*q.d*fac;
	bc = q.b*q.c*fac;
	bd = q.b*q.d*fac;
	cd = q.c*q.d*fac;

	operator()(0,0) = a2+b2-c2-d2;
	operator()(0,1) = 2*bc+2*ad;
	operator()(0,2) = 2*bd-2*ac;

	operator()(1,0) = 2*bc-2*ad;
	operator()(1,1) = a2-b2+c2-d2;
	operator()(1,2) = 2*cd+2*ab;

	operator()(2,0) = 2*bd+2*ac;
	operator()(2,1) = 2*cd-2*ab;
	operator()(2,2) = a2-b2-c2+d2;
	return *this;
}


void L_Rot3DMatrix::jacob_fijarCuaternionInv(const L_Quaternion &q, L_Matrix &J4x9, int i0, int j0)
{
	double a, b, c, d;
	a = q.a;
	b = q.b;
	c = q.c;
	d = q.d;

	// operator()(0,0) = a2+b2-c2-d2;
	J4x9(0,0) = 2*a;
	J4x9(1,0) = 2*b;
	J4x9(2,0) = -2*c;
	J4x9(3,0) = -2*d;

	//1
	// operator()(0,1) = 2*bc+2*ad;
	J4x9(0,1) = 2*d;
	J4x9(1,1) = 2*c;
	J4x9(2,1) = 2*b;
	J4x9(3,1) = 2*a;

	//2
	// operator()(0,2) = 2*bd - 2*ac;
	J4x9(0,2) = 2*c;
	J4x9(1,2) = -2*d;
	J4x9(2,2) = 2*a;
	J4x9(3,2) = -2*b;

	//3
	// operator()(1,0) = 2*bc-2*ad;
	J4x9(0,3) = -2*d;
	J4x9(1,3) = 2*c;
	J4x9(2,3) = 2*b;
	J4x9(3,3) = -2*a;

	//4
	// operator()(1,1) = a2-b2+c2-d2;
	J4x9(0,4) = 2*a;
	J4x9(1,4) = -2*b;
	J4x9(2,4) = 2*c;
	J4x9(3,4) = -2*d;

	//5
	// operator()(1,2) = 2*cd + 2*ab;
	J4x9(0,5) = 2*b;
	J4x9(1,5) = 2*a;
	J4x9(2,5) = 2*d;
	J4x9(3,5) = 2*c;

	//6
	//operator()(2,0) = 2*bd+2*ac;
	J4x9(0,6) = 2*c;
	J4x9(1,6) = 2*d;
	J4x9(2,6) = 2*a;
	J4x9(3,6) = 2*b;

	//7
	// operator()(2,1) = 2*cd - 2*ab;
	J4x9(0,7) = -2*b;
	J4x9(1,7) = -2*a;
	J4x9(2,7) = 2*d;
	J4x9(3,7) = 2*c;

	//8
	// operator()(2,2) = a2-b2-c2+d2;
	J4x9(0,8) = 2*a;
	J4x9(1,8) = -2*b;
	J4x9(2,8) = -2*c;
	J4x9(3,8) = -2*d;
}


L_Rot3DMatrix &L_Rot3DMatrix::fijarCuaternion(const L_Quaternion &q)
{
	double a2, b2, c2, d2, ab, ac, ad, bc, bd, cd, fac2;
	fac2 = q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d;
	if (fac2 == 0)
		printf("L_Rot3DMatrix::fijarCuaternion() : cuaternion con modulo zero\n");
	fac2 = 1/fac2;
	a2 = fac2*q.a*q.a;
	b2 = fac2*q.b*q.b;
	c2 = fac2*q.c*q.c;
	d2 = fac2*q.d*q.d;
	ab = fac2*q.a*q.b;
	ac = fac2*q.a*q.c;
	ad = fac2*q.a*q.d;
	bc = fac2*q.b*q.c;
	bd = fac2*q.b*q.d;
	cd = fac2*q.c*q.d;
	operator()(0,0) = a2+b2-c2-d2;
	operator()(0,1) = 2*bc-2*ad; // -
	operator()(0,2) = 2*bd+2*ac; // +

	operator()(1,0) = 2*bc+2*ad; // +
	operator()(1,1) = a2-b2+c2-d2;
	operator()(1,2) = 2*cd - 2*ab;  // -

	operator()(2,0) = 2*bd - 2*ac; // -
	operator()(2,1) = 2*cd + 2*ab; // +
	operator()(2,2) = a2-b2-c2+d2;
	return *this;
}


void L_Rot3DMatrix::jacob_fijarCuaternion(const L_Quaternion &q, L_Matrix &J4x9, int i0, int j0)
{
	double a, b, c, d;
	a = q.a;
	b = q.b;
	c = q.c;
	d = q.d;
	//operator()(0,0) = a2+b2-c2-d2;
	J4x9(0,0) = 2*a;
	J4x9(1,0) = 2*b;
	J4x9(2,0) = -2*c;
	J4x9(3,0) = -2*d;

	//operator()(0,1) = 2*bc-2*ad; // -
	J4x9(0,1) = -2*d;
	J4x9(1,1) = 2*c;
	J4x9(2,1) = 2*b;
	J4x9(3,1) = -2*a;

	//operator()(0,2) = 2*bd+2*ac; // +
	J4x9(0,2) = 2*c;
	J4x9(1,2) = 2*d;
	J4x9(2,2) = 2*a;
	J4x9(3,2) = 2*b;

	//operator()(1,0) = 2*bc+2*ad; // +
	J4x9(0,3) = 2*d;
	J4x9(1,3) = 2*c;
	J4x9(2,3) = 2*b;
	J4x9(3,3) = 2*a;

	//operator()(1,1) = a2-b2+c2-d2;
	J4x9(0,4) = 2*a;
	J4x9(1,4) = -2*b;
	J4x9(2,4) = 2*c;
	J4x9(3,4) = -2*d;

	//operator()(1,2) = 2*cd - 2*ab;  // -
	J4x9(0,5) = -2*b;
	J4x9(1,5) = -2*a;
	J4x9(2,5) = 2*d;
	J4x9(3,5) = 2*c;

	//operator()(2,0) = 2*bd - 2*ac; // -
	J4x9(0,6) = -2*c;
	J4x9(1,6) = 2*d;
	J4x9(2,6) = -2*a;
	J4x9(3,6) = 2*b;

	//operator()(2,1) = 2*cd + 2*ab; // +
	J4x9(0,7) = 2*b;
	J4x9(1,7) = 2*a;
	J4x9(2,7) = 2*d;
	J4x9(3,7) = 2*c;

	//operator()(2,2) = a2-b2-c2+d2;
	J4x9(0,8) = 2*a;
	J4x9(1,8) = -2*b;
	J4x9(2,8) = -2*c;
	J4x9(3,8) = -2*d;
}

L_Rot3DMatrix &L_Rot3DMatrix::fijarVectorRotacion(const L_CoordsCart3D &eje)
{
	L_Quaternion q;
	q.fijaVectorRotacion(eje);
	fijarCuaternion(q);
	return *this;
}

L_Rot3DMatrix &L_Rot3DMatrix::fijarPanTiltRoll(const L_PanTiltRoll &ori)
{
	double cp=cos(ori.pan), sp=sin(ori.pan);
	double ct=cos(ori.tilt), st=sin(ori.tilt);
	double cr=cos(ori.roll), sr=sin(ori.roll);

	// Ahora, la matriz de rotacion compuesta
	operator()(0,0)=ct*cp;     operator()(0,1)=-sr*st*cp-cr*sp;     operator()(0,2)=-cr*st*cp+sr*sp;
	operator()(1,0)=ct*sp;     operator()(1,1)=-sr*st*sp+cr*cp;     operator()(1,2)=-cr*st*sp-sr*cp;
	operator()(2,0)=st;		  operator()(2,1)=sr*ct;               operator()(2,2)=cr*ct;
	return *this;
}




void L_HomogeneousMatrix::OP_mult(const L_HomogeneousMatrix &m1, const L_HomogeneousMatrix &m2)
{
	int i,j,u;

	for (i=0; i<3; i++)
	{
		for (j=0; j<4; j++)
		{
			operator()(i,j)=0;
			for (u=0; u<4; u++)
				operator()(i,j)+=m1(i,u)*m2(u,j);
		}
	}
	operator()(3,0) = operator()(3,1) = operator()(3,2) = 0.0; operator()(3,3) = 1.0;
	return;
}


L_HomogeneousMatrix &L_HomogeneousMatrix::operator*=(const L_HomogeneousMatrix &other)
{
	L_HomogeneousMatrix temp;
	L_StaticMatrix_OP_assign(temp,*this);
	OP_mult(temp, other);
	return *this;
}

void L_HomogeneousMatrix::izqArrAde_2_adeIzqArr(L_HomogeneousMatrix &otraDAA)
{
	operator()(0,0) = otraDAA(2,2);
	operator()(0,1) = otraDAA(2,0);
	operator()(0,2) = otraDAA(2,1);
	operator()(1,0) = otraDAA(0,2);
	operator()(1,1) = otraDAA(0,0);
	operator()(1,2) = otraDAA(0,1);
	operator()(2,0) = otraDAA(1,2);
	operator()(2,1) = otraDAA(1,0);
	operator()(2,2) = otraDAA(1,1);
	operator()(0,3)=otraDAA(2,3);
	operator()(1,3)=otraDAA(0,3);
	operator()(2,3)=otraDAA(1,3);
}

void L_HomogeneousMatrix::adeIzqArr_2_izqArrAde(L_HomogeneousMatrix &otraAIA)
{
	operator()(0,0) = otraAIA(1,1);
	operator()(0,1) = otraAIA(1,2);
	operator()(0,2) = otraAIA(1,0);
	operator()(1,0) = otraAIA(2,1);
	operator()(1,1) = otraAIA(2,2);
	operator()(1,2) = otraAIA(2,0);
	operator()(2,0) = otraAIA(0,1);
	operator()(2,1) = otraAIA(0,2);
	operator()(2,2) = otraAIA(0,0);
	operator()(0,3) = otraAIA(1,3);
	operator()(1,3) = otraAIA(2,3);
	operator()(2,3) = otraAIA(0,3);
}

void L_HomogeneousMatrix::inversaTrasponiendoDe(const L_HomogeneousMatrix& other)
{
	// inv([R,t;0,1)) = [RT,-RTt;0,1)
	operator()(0,0) = other(0,0);
	operator()(0,1) = other(1,0);
	operator()(0,2) = other(2,0);
	
	operator()(1,0) = other(0,1);
	operator()(1,1) = other(1,1);
	operator()(1,2) = other(2,1);

	operator()(2,0) = other(0,2);
	operator()(2,1) = other(1,2);
	operator()(2,2) = other(2,2);

	operator()(0,3) = - other(0,0)*other(0,3) - other(1,0)*other(1,3) - other(2,0)*other(2,3);
	operator()(1,3) = - other(0,1)*other(0,3) - other(1,1)*other(1,3) - other(2,1)*other(2,3);
	operator()(2,3) = - other(0,2)*other(0,3) - other(1,2)*other(1,3) - other(2,2)*other(2,3);
}

void L_HomogeneousMatrix::invertMeTrasponiendo()
{
	L_HomogeneousMatrix other(*this);
	// inv([R,t;0,1)) = [RT,-RTt;0,1)
	operator()(0,0) = other(0,0);
	operator()(0,1) = other(1,0);
	operator()(0,2) = other(2,0);
	
	operator()(1,0) = other(0,1);
	operator()(1,1) = other(1,1);
	operator()(1,2) = other(2,1);

	operator()(2,0) = other(0,2);
	operator()(2,1) = other(1,2);
	operator()(2,2) = other(2,2);

	operator()(0,3) = - other(0,0)*other(0,3) - other(1,0)*other(1,3) - other(2,0)*other(2,3);
	operator()(1,3) = - other(0,1)*other(0,3) - other(1,1)*other(1,3) - other(2,1)*other(2,3);
	operator()(2,3) = - other(0,2)*other(0,3) - other(1,2)*other(1,3) - other(2,2)*other(2,3);}

void L_HomogeneousMatrix::jacob_inversaTrasponiendoDe(const L_HomogeneousMatrix& other, L_Matrix &J12x12)
{
	J12x12.reallocate(12, 12);
	J12x12.setZero();

	J12x12(0,0) = 1;
	J12x12(1,4) = 1;
	J12x12(2,8) = 1;
	
	J12x12(4,1) = 1;
	J12x12(5,5) = 1;
	J12x12(6,9) = 1;

	J12x12(8,2) = 1;
	J12x12(9,6) = 1;
	J12x12(10,10) = 1;

	J12x12(3,0) = -other(0,3);
	J12x12(3,3) = -other(0,0);

	J12x12(3,4) = -other(1,3);
	J12x12(3,7) = -other(1,0);

	J12x12(3,8) = -other(2,3);
	J12x12(3,11) = -other(2,0);

	J12x12(7,1) = -other(0,3);
	J12x12(7,3) = -other(0,1);
	
	J12x12(7,5) = -other(1,3);
	J12x12(7,7) = -other(1,1);

	J12x12(7,9) = -other(2,3);
	J12x12(7,11) = -other(2,1);

	J12x12(11,2) = -other(0,3);
	J12x12(11,3) = -other(0,2);

	J12x12(11,6) = -other(1,3);
	J12x12(11,7) = -other(1,2);

	J12x12(11,10) = -other(2,3);
	J12x12(11,11) = -other(2,2);
}

void L_HomogeneousMatrix::jacob_inversaAproxDe(const L_HomogeneousMatrix& other, L_Matrix &J12x12)
{
	J12x12.reallocate(12, 12);
	J12x12.setZero();
	double det =
		+ other(0,0)*(other(1,1)*other(2,2)-other(1,2)*other(2,1))
		- other(0,1)*(other(1,0)*other(2,2)-other(1,2)*other(2,0))
		+ other(0,2)*(other(1,0)*other(2,1)-other(2,2)*other(3,0));
	double fac2 = -1/(det*det); // diff(1/det, det)
	double fac3 = -2/(det*det*det); // diff(1/det^2, det)

	J12x12(0,0) = fac3*1;
	J12x12(1,4) = fac3*1;
	J12x12(2,8) = fac3*1;
	
	J12x12(4,1) = fac3*1;
	J12x12(5,5) = fac3*1;
	J12x12(6,9) = fac3*1;

	J12x12(8,2) = fac3*1;
	J12x12(9,6) = fac3*1;
	J12x12(10,10) = fac3*1;

	J12x12(3,0) = -fac2*other(0,3);
	J12x12(3,3) = -other(0,0);

	J12x12(3,4) = -fac2*other(1,3);
	J12x12(3,7) = -other(1,0);

	J12x12(3,8) = -fac2*other(2,3);
	J12x12(3,11) = -other(2,0);

	J12x12(7,1) = -fac2*other(0,3);
	J12x12(7,3) = -other(0,1);
	
	J12x12(7,5) = -fac2*other(1,3);
	J12x12(7,7) = -other(1,1);

	J12x12(7,9) = -fac2*other(2,3);
	J12x12(7,11) = -other(2,1);

	J12x12(11,2) = -fac2*other(0,3);
	J12x12(11,3) = -other(0,2);

	J12x12(11,6) = -fac2*other(1,3);
	J12x12(11,7) = -other(1,2);

	J12x12(11,10) = -fac2*other(2,3);
	J12x12(11,11) = -other(2,2);
}

void L_HomogeneousMatrix::errorRot_prepara(void *matriz, std::vector<double> &panTiltRoll, bool prenormalizar)
{
	int i, j;
	double a=0;
	L_Pose3D pose;
	if (prenormalizar)
	{
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				a+=operator()(i,j)*operator()(i,j);
		a = sqrt(a/3);
		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				operator()(i,j)/=a;
	}
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			((double *)matriz)[i*3+j]=operator()(i,j);
	pose=calcPose3D_fijo_a_movido();
	panTiltRoll[0]=pose.ori.pan;
	panTiltRoll[1]=pose.ori.tilt;
	panTiltRoll[2]=pose.ori.roll;
}

double L_HomogeneousMatrix::errorRot(const void *matriz, double *panTiltRoll)
{
	L_HomogeneousMatrix m;
	double *matDato=(double *)matriz;
	double err = 0;
	int i, j;
	m.fijaPose3D_fijo_a_movido(L_Pose3D(0,0,0,panTiltRoll[0],panTiltRoll[1],panTiltRoll[2]));
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			err+=(m(i,j)-matDato[i*3+j])*(m(i,j)-matDato[i*3+j]);
	return err;
}


L_Rayo3D L_HomogeneousMatrix::operator *(const L_Rayo3D &other)
{
	L_Rayo3D ret;
	ret.origen=(*this)*other.origen;
	ret.direcc=OP_rota(other.direcc);
	ret.normaliza(); // Por si acaso
	return ret;
};




L_Pose3D L_HomogeneousMatrix::calcPose3D() const
{
	L_Pose3D pose;

	if (L_MatrizHomogenea_revisa(*this) == false)
		throw L_ArgException();

	pose.pos.x = operator()(0,3);
	pose.pos.y = operator()(1,3);
	pose.pos.z = operator()(2,3);
	double ct=sqrt(operator()(0,0)*operator()(0,0)+operator()(1,0)*operator()(1,0));
	if (ct > 0)
	{
		pose.ori.pan=atan2(operator()(1,0),operator()(0,0));
		pose.ori.tilt=atan2(operator()(2,0),ct);
		pose.ori.roll=atan2(operator()(2,1),operator()(2,2));
	}
	else
	{
		pose.ori.roll=0;
		if (operator()(2,0) > 0)
		{
			pose.ori.tilt=M_PI;
			pose.ori.pan = atan2(operator()(2,1),operator()(1,1));
		}
		else
		{
			pose.ori.tilt=-M_PI;
			pose.ori.pan = atan2(operator()(2,1),operator()(1,1));
		}
	}
	return pose;
}

L_Pose3D L_HomogeneousMatrix::calcPose3D_fijo_a_movido() const
{
	L_Pose3D pose;

	if (L_MatrizHomogenea_revisa(*this) == false)
		throw L_ArgException();

	pose.pos.x = -operator()(0,3)*operator()(0,0)-operator()(1,3)*operator()(1,0)-operator()(2,3)*operator()(2,0);
	pose.pos.y = -operator()(0,3)*operator()(0,1)-operator()(1,3)*operator()(1,1)-operator()(2,3)*operator()(2,1);
	pose.pos.z = -operator()(0,3)*operator()(0,2)-operator()(1,3)*operator()(1,2)-operator()(2,3)*operator()(2,2);
	double ct=sqrt(operator()(0,0)*operator()(0,0)+operator()(0,1)*operator()(0,1));
	if (ct > 0)
	{
		pose.ori.pan=atan2(operator()(0,1),operator()(0,0));
		pose.ori.tilt=atan2(operator()(0,2),ct);
		pose.ori.roll=atan2(operator()(1,2),operator()(2,2));
	}
	else
	{
		pose.ori.roll=0;
		if (operator()(0,2) > 0)
		{
			pose.ori.tilt=M_PI;
			pose.ori.pan = atan2(operator()(1,2),operator()(1,1));
		}
		else
		{
			pose.ori.tilt=-M_PI;
			pose.ori.pan = atan2(operator()(1,2),operator()(1,1));
		}
	}
	return pose;
}

void L_HomogeneousMatrix::calcEsencial(L_Matrix &E) const
{
	E.reallocate(3,3);

}

void L_HomogeneousMatrix::fijaPose3D(const L_Pose3D &poseMovido)
{
	double cp=cos(poseMovido.ori.pan), sp=sin(poseMovido.ori.pan);
	double ct=cos(poseMovido.ori.tilt), st=sin(poseMovido.ori.tilt);
	double cr=cos(poseMovido.ori.roll), sr=sin(poseMovido.ori.roll);
	// Ahora, la matriz de rotacion compuesta

	operator()(0,0)=ct*cp;     operator()(0,1)=-sr*st*cp-cr*sp;     operator()(0,2)=-cr*st*cp+sr*sp;
	operator()(1,0)=ct*sp;     operator()(1,1)=-sr*st*sp+cr*cp;     operator()(1,2)=-cr*st*sp-sr*cp;
	operator()(2,0)=st;		  operator()(2,1)=sr*ct;               operator()(2,2)=cr*ct;

	operator()(0,3)=poseMovido.pos.x;
	operator()(1,3)=poseMovido.pos.y;
	operator()(2,3)=poseMovido.pos.z;
}

void L_HomogeneousMatrix::fijaPose3D_lento(const L_Pose3D &poseMovido)
{
	// Supongamos que tenemos los ejes {i,j,k}(fijos),  {i',j',k'}(movidos)
	// Sea R la rotacion que lleva los ejes {i,j,k} a {i',j',k'}
	// Sea t la traslacion que lleva el origen de los ejes {i,j,k} al origen de los ejes {i',j',k'}
	// entonces si:
	//     P = xi + yj + zk = x'i' + y'j' + z'k'
	//    (x,y,z)=H(x',y',z')
	//
	//     => H=[R t)
	L_HomogeneousMatrix pa, ti, ro, t;
	
	pa.fijaRotEjeZ(poseMovido.ori.pan);
	ti.fijaRotEjeY(-poseMovido.ori.tilt);
	ro.fijaRotEjeX(poseMovido.ori.roll);
	t.fijaTr(poseMovido.pos.x,poseMovido.pos.y,poseMovido.pos.z);

	*this= t * pa * ti * ro; // Primero roll, luego tilt, luego pan, luego la traslacion
}

void L_HomogeneousMatrix::fijaPose3D_fijo_a_movido(const L_Pose3D &poseMovido)
{
	double cp=cos(poseMovido.ori.pan), sp=sin(poseMovido.ori.pan);
	double ct=cos(poseMovido.ori.tilt), st=sin(poseMovido.ori.tilt);
	double cr=cos(poseMovido.ori.roll), sr=sin(poseMovido.ori.roll);
	// Ahora, la matriz de rotacion compuesta
	operator()(0,0)=ct*cp;			operator()(0,1)=ct*sp;			operator()(0,2)=st;		
	operator()(1,0)=-sr*st*cp-cr*sp;	operator()(1,1)=-sr*st*sp+cr*cp;	operator()(1,2)=sr*ct;	
	operator()(2,0)=-cr*st*cp+sr*sp;	operator()(2,1)=-cr*st*sp-sr*cp;	operator()(2,2)=cr*ct;	
	operator()(0,3) = -operator()(0,0)*poseMovido.pos.x-operator()(0,1)*poseMovido.pos.y-operator()(0,2)*poseMovido.pos.z;
	operator()(1,3) = -operator()(1,0)*poseMovido.pos.x-operator()(1,1)*poseMovido.pos.y-operator()(1,2)*poseMovido.pos.z;
	operator()(2,3) = -operator()(2,0)*poseMovido.pos.x-operator()(2,1)*poseMovido.pos.y-operator()(2,2)*poseMovido.pos.z;
}

void L_HomogeneousMatrix::fijaPose3D_fijo_a_movido_lento(const L_Pose3D &poseMovido)
{
	L_HomogeneousMatrix other;
	other.fijaPose3D_lento(poseMovido);
	inversaTrasponiendoDe(other);
}

L_Pose3D_cuat L_HomogeneousMatrix::calcPose3D_cuat() const
{
	L_Pose3D_cuat poseMovido;

	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*poseMovido.ori.a 
		poseMovido.ori.a = 0.25 * S;
		poseMovido.ori.b = (operator()(2,1) - operator()(1,2)) / S;
		poseMovido.ori.c = (operator()(0,2) - operator()(2,0)) / S; 
		poseMovido.ori.d = (operator()(1,0) - operator()(0,1)) / S; 
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*poseMovido.ori.b 
		poseMovido.ori.a = (operator()(2,1) - operator()(1,2)) / S;
		poseMovido.ori.b = 0.25 * S;
		poseMovido.ori.c = (operator()(0,1) + operator()(1,0)) / S; 
		poseMovido.ori.d = (operator()(0,2) + operator()(2,0)) / S; 
	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*poseMovido.ori.c
		poseMovido.ori.a = (operator()(0,2) - operator()(2,0)) / S;
		poseMovido.ori.b = (operator()(0,1) + operator()(1,0)) / S; 
		poseMovido.ori.c = 0.25 * S;
		poseMovido.ori.d = (operator()(1,2) + operator()(2,1)) / S; 
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*poseMovido.ori.d
		poseMovido.ori.a = (operator()(1,0) - operator()(0,1)) / S;
		poseMovido.ori.b = (operator()(0,2) + operator()(2,0)) / S;
		poseMovido.ori.c = (operator()(1,2) + operator()(2,1)) / S;
		poseMovido.ori.d = 0.25 * S;
	}
	poseMovido.pos.x=operator()(0,3);
	poseMovido.pos.y=operator()(1,3);
	poseMovido.pos.z=operator()(2,3);
	return poseMovido;
}


bool L_HomogeneousMatrix::fijaPose3D_cuat(const L_Pose3D_cuat &poseMovido)
{
	L_Pose3D_cuat pm;
	double a2, b2, c2, d2, ab, ac, ad, bc, bd, cd;
	pm = poseMovido;
	if (pm.normaliza_ret() == false)
		return false;

	a2 = pm.ori.a*pm.ori.a;
	b2 = pm.ori.b*pm.ori.b;
	c2 = pm.ori.c*pm.ori.c;
	d2 = pm.ori.d*pm.ori.d;
	ab = pm.ori.a*pm.ori.b;
	ac = pm.ori.a*pm.ori.c;
	ad = pm.ori.a*pm.ori.d;
	bc = pm.ori.b*pm.ori.c;
	bd = pm.ori.b*pm.ori.d;
	cd = pm.ori.c*pm.ori.d;
	operator()(0,0) = a2+b2-c2-d2;
	operator()(0,1) = 2*bc-2*ad; // -
	operator()(0,2) = 2*bd+2*ac; // +

	operator()(1,0) = 2*bc+2*ad; // +
	operator()(1,1) = a2-b2+c2-d2;
	operator()(1,2) = 2*cd - 2*ab;  // -

	operator()(2,0) = 2*bd - 2*ac; // -
	operator()(2,1) = 2*cd + 2*ab; // +
	operator()(2,2) = a2-b2-c2+d2;

	operator()(0,3)=pm.pos.x;
	operator()(1,3)=pm.pos.y;
	operator()(2,3)=pm.pos.z;
	return true;
}

void L_HomogeneousMatrix::fijaPose3D_cuat_fijo_a_movido(const L_Pose3D_cuat &poseMovido)
{
	L_Pose3D_cuat pm;
	double a2, b2, c2, d2, ab, ac, ad, bc, bd, cd;
	pm = poseMovido;
	pm.normaliza();

	a2 = pm.ori.a*pm.ori.a;
	b2 = pm.ori.b*pm.ori.b;
	c2 = pm.ori.c*pm.ori.c;
	d2 = pm.ori.d*pm.ori.d;
	ab = pm.ori.a*pm.ori.b;
	ac = pm.ori.a*pm.ori.c;
	ad = pm.ori.a*pm.ori.d;
	bc = pm.ori.b*pm.ori.c;
	bd = pm.ori.b*pm.ori.d;
	cd = pm.ori.c*pm.ori.d;

	operator()(0,0) = a2+b2-c2-d2;
	operator()(0,1) = 2*bc+2*ad; // +
	operator()(0,2) = 2*bd - 2*ac; // -

	operator()(1,0) = 2*bc-2*ad; // -
	operator()(1,1) = a2-b2+c2-d2;
	operator()(1,2) = 2*cd + 2*ab; // +

	operator()(2,0) = 2*bd+2*ac; // +
	operator()(2,1) = 2*cd - 2*ab;  // -
	operator()(2,2) = a2-b2-c2+d2;

	operator()(0,3) = -operator()(0,0)*poseMovido.pos.x - operator()(0,1)*poseMovido.pos.y - operator()(0,2)*poseMovido.pos.z;
	operator()(1,3) = -operator()(1,0)*poseMovido.pos.x - operator()(1,1)*poseMovido.pos.y - operator()(1,2)*poseMovido.pos.z;
	operator()(2,3) = -operator()(2,0)*poseMovido.pos.x - operator()(2,1)*poseMovido.pos.y - operator()(2,2)*poseMovido.pos.z;
}


L_Pose3D_cuat L_HomogeneousMatrix::calcPose3D_cuat_fijo_a_movido() const
{
	L_Pose3D_cuat poseMovido;
	// Calculo de rotacion igual a la version normal pero con la matriz traspuesta
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*poseMovido.ori.a 
		poseMovido.ori.a = 0.25 * S;
		poseMovido.ori.b = (operator()(1,2) - operator()(2,1)) / S;
		poseMovido.ori.c = (operator()(2,0) - operator()(0,2)) / S; 
		poseMovido.ori.d = (operator()(0,1) - operator()(1,0)) / S; 
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*poseMovido.ori.b 
		poseMovido.ori.a = (operator()(1,2) - operator()(2,1)) / S;
		poseMovido.ori.b = 0.25 * S;
		poseMovido.ori.c = (operator()(1,0) + operator()(0,1)) / S; 
		poseMovido.ori.d = (operator()(2,0) + operator()(0,2)) / S; 
	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*poseMovido.ori.c
		poseMovido.ori.a = (operator()(2,0) - operator()(0,2)) / S;
		poseMovido.ori.b = (operator()(1,0) + operator()(0,1)) / S; 
		poseMovido.ori.c = 0.25 * S;
		poseMovido.ori.d = (operator()(2,1) + operator()(1,2)) / S; 
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*poseMovido.ori.d
		poseMovido.ori.a = (operator()(0,1) - operator()(1,0)) / S;
		poseMovido.ori.b = (operator()(2,0) + operator()(0,2)) / S;
		poseMovido.ori.c = (operator()(2,1) + operator()(1,2)) / S;
		poseMovido.ori.d = 0.25 * S;
	}
	poseMovido.pos.x = -operator()(0,3)*operator()(0,0)-operator()(1,3)*operator()(1,0)-operator()(2,3)*operator()(2,0);
	poseMovido.pos.y = -operator()(0,3)*operator()(0,1)-operator()(1,3)*operator()(1,1)-operator()(2,3)*operator()(2,1);
	poseMovido.pos.z = -operator()(0,3)*operator()(0,2)-operator()(1,3)*operator()(1,2)-operator()(2,3)*operator()(2,2);
	return poseMovido;
}

L_Pose3D_cuat L_HomogeneousMatrix::jacob_calcPose3D_cuat(L_Matrix &J7x12, int i0, int j0)
{
	L_Pose3D_cuat q;
	int i, j;
	for (i=0; i<7; i++)
		for (j=0; j<12; j++)
			J7x12(i+i0,j+j0) = 0;
	double tr = operator()(0,0) + operator()(1,1) + operator()(2,2);
	int e00=j0+0, e01=j0+1, e02=j0+2,  e10=j0+4, e11=j0+5, e12=j0+6,  e20=j0+8, e21=j0+9, e22=j0+10;
	if (tr > 0) { 
		double S = sqrt(tr+1.0) * 2; // S=4*q.ori.a 
		q.ori.a = 0.25 * S;
		q.ori.b = (operator()(2,1) - operator()(1,2)) / S;
		q.ori.c = (operator()(0,2) - operator()(2,0)) / S;
		q.ori.d = (operator()(1,0) - operator()(0,1)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J7x12(i0+3,e00) = A; J7x12(i0+3,e11) = A; J7x12(i0+3,e22) = A;
		J7x12(i0+4,e00) = -q.ori.b*C; J7x12(i0+4,e11) = -q.ori.b*C; J7x12(i0+4,e22) = -q.ori.b*C;
		J7x12(i0+4,e21) = B; J7x12(i0+4,e12) = -B;
		J7x12(i0+5,e00) = -q.ori.c*C; J7x12(i0+5,e11) = -q.ori.c*C; J7x12(i0+5,e22) = -q.ori.c*C;
		J7x12(i0+5,e02) = B; J7x12(i0+5,e20) = -B;
		J7x12(i0+6,e00) = -q.ori.d*C; J7x12(i0+6,e11) = -q.ori.d*C; J7x12(i0+6,e22) = -q.ori.d*C;
		J7x12(i0+6,e10) = B; J7x12(i0+6,e01) = -B;
	} else if ((operator()(0,0) > operator()(1,1))&&(operator()(0,0) > operator()(2,2))) { 
		double S = sqrt(1.0 + operator()(0,0) - operator()(1,1) - operator()(2,2)) * 2; // S=4*q.ori.b 
		q.ori.a = (operator()(2,1) - operator()(1,2)) / S;
		q.ori.b = 0.25 * S;
		q.ori.c = (operator()(0,1) + operator()(1,0)) / S;
		q.ori.d = (operator()(0,2) + operator()(2,0)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J7x12(i0+3,e00) = -q.ori.a*C; J7x12(i0+3,e11) = q.ori.a*C; J7x12(i0+3,e22) = q.ori.a*C;
		J7x12(i0+3,e21) = B; J7x12(i0+3,e12) = -B; // ERROR
		J7x12(i0+4,e00) = A; J7x12(i0+4,e11) = J7x12(i0+4,e22) = -A;
		J7x12(i0+5,e00) = -q.ori.c*C; J7x12(i0+5,e11) = q.ori.c*C;J7x12(i0+5,e22) = q.ori.c*C;
		J7x12(i0+5,e01) = B; J7x12(i0+5,e10) = B;
		J7x12(i0+6,e00) = -q.ori.d*C; J7x12(i0+6,e11) = q.ori.d*C; J7x12(i0+6,e22) = q.ori.d*C;
		J7x12(i0+6,e02) = B; J7x12(i0+6,e20) = B;

	} else if (operator()(1,1) > operator()(2,2)) { 
		double S = sqrt(1.0 + operator()(1,1) - operator()(0,0) - operator()(2,2)) * 2; // S=4*q.ori.c
		q.ori.a = (operator()(0,2) - operator()(2,0)) / S;
		q.ori.b = (operator()(0,1) + operator()(1,0)) / S;
		q.ori.c = 0.25 * S;
		q.ori.d = (operator()(1,2) + operator()(2,1)) / S;
		double rc = S/2;
		double r32 = rc*rc*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J7x12(i0+3,e00) = q.ori.a*C; J7x12(i0+3,e11) = -q.ori.a*C; J7x12(i0+3,e22) = q.ori.a*C;
		J7x12(i0+3,e02) = B; J7x12(i0+3,e20) = -B;
		J7x12(i0+4,e00) = q.ori.b*C; J7x12(i0+4,e11) = -q.ori.b*C; J7x12(i0+4,e22) = q.ori.b*C;
		J7x12(i0+4,e01) = B; J7x12(i0+4,e10) = B;
		J7x12(i0+5,e00) = -A; J7x12(i0+5,e11) = A; J7x12(i0+5,e22) = -A;
		J7x12(i0+6,e00) = q.ori.d*C; J7x12(i0+6,e11) = -q.ori.d*C; J7x12(i0+6,e22) = q.ori.d*C;
		J7x12(i0+6,e12) = B; J7x12(i0+6,e21) = B;
	} else { 
		double S = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1)) * 2; // S=4*q.ori.d
		q.ori.a = (operator()(1,0) - operator()(0,1)) / S;
		q.ori.b = (operator()(0,2) + operator()(2,0)) / S;
		q.ori.c = (operator()(1,2) + operator()(2,1)) / S;
		q.ori.d = 0.25 * S;
		double rc = sqrt(1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1));
		double r32 = (1.0 + operator()(2,2) - operator()(0,0) - operator()(1,1))*rc;
		double A=1/(4*rc), B=1/(2*rc), C=S/(4*r32);
		J7x12(i0+3,e00) = q.ori.a*C; J7x12(i0+3,e11) = q.ori.a*C; J7x12(i0+3,e22) = -q.ori.a*C;
		J7x12(i0+3,e10) = B; J7x12(i0+3,e01) = -B;
		J7x12(i0+4,e00) = q.ori.b*C; J7x12(i0+4,e11) = q.ori.b*C; J7x12(i0+4,e22) = -q.ori.b*C;
		J7x12(i0+4,e02) = B; J7x12(i0+4,e20) = B;
		J7x12(i0+5,e00) = q.ori.c*C; J7x12(i0+5,e11) = q.ori.c*C;J7x12(i0+5,e22) = -q.ori.c*C;
		J7x12(i0+5,e12) = B; J7x12(i0+5,e21) = B;
		J7x12(i0+6,e00) = -A; J7x12(i0+6,e11) = -A; J7x12(i0+3,e22) = A;
	}
	q.pos.x = operator()(0,3);
	q.pos.y = operator()(0,7);
	q.pos.z = operator()(0,11);
	J7x12(i0+0,j0+3)=1;
	J7x12(i0+1,j0+7)=1;
	J7x12(i0+2,j0+11)=1;
	return q;
}

void L_HomogeneousMatrix::jacob_der_fijaPose3D_cuat_fijo_a_movido_lento(L_Pose3D_cuat &poseMovido, L_Array<L_Matrix> &J)
{
	L_HomogeneousMatrix mFMD[7], mFMD2[7];
	const double delta = 1.0e-7;
	int i, j, k;
	
	fijaPose3D_cuat_fijo_a_movido(poseMovido);
	// Calcular valores modificados
	poseMovido.pos.x += delta;
	mFMD[0].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.x -= 2*delta;
	mFMD2[0].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.x += delta;

	poseMovido.pos.y += delta;
	mFMD[1].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.y -= 2*delta;
	mFMD2[1].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.y += delta;

	poseMovido.pos.z += delta;
	mFMD[2].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.z -= 2*delta;
	mFMD2[2].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.pos.z += delta;

	poseMovido.ori.a += delta;
	mFMD[3].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.a -= 2*delta;
	mFMD2[3].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.a += delta;

	poseMovido.ori.b += delta;
	mFMD[4].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.b -= 2*delta;
	mFMD2[4].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.b += delta;

	poseMovido.ori.c += delta;
	mFMD[5].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.c -= 2*delta;
	mFMD2[5].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.c += delta;

	poseMovido.ori.d += delta;
	mFMD[6].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.d -= 2*delta;
	mFMD2[6].fijaPose3D_cuat_fijo_a_movido(poseMovido);
	poseMovido.ori.d += delta;

	// restar media a los valores modificados para obtener derivadas
	for (k=0; k<7; k++)
	{
		J[k].reallocate(4,4);
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				J[k](i,j) = (mFMD[k](i,j) - mFMD2[k](i,j))/(2*delta);
	}
}

void L_HomogeneousMatrix::guardaInfoEn(L_Matrix &vector)
{
	vector.reallocate(12,1);
	vector(0,0)=operator()(0,0);
	vector(1,0)=operator()(1,0);
	vector(2,0)=operator()(2,0);
	vector(3,0)=operator()(0,1);
	vector(4,0)=operator()(1,1);
	vector(5,0)=operator()(2,1);
	vector(6,0)=operator()(0,2);
	vector(7,0)=operator()(1,2);
	vector(8,0)=operator()(2,2);
	vector(9,0)=operator()(0,3);
	vector(10,0)=operator()(1,3);
	vector(11,0)=operator()(2,3);
}

void L_HomogeneousMatrix::leeInfoDe(const L_Matrix &vector)
{
	operator()(0,0)=vector(0,0);
	operator()(1,0)=vector(1,0);
	operator()(2,0)=vector(2,0);
	operator()(0,1)=vector(3,0);
	operator()(1,1)=vector(4,0);
	operator()(2,1)=vector(5,0);
	operator()(0,2)=vector(6,0);
	operator()(1,2)=vector(7,0);
	operator()(2,2)=vector(8,0);
	operator()(0,3)=vector(9,0);
	operator()(1,3)=vector(10,0);
	operator()(2,3)=vector(11,0);
}

//! Esta funcion erase_preserving_order la escala y la distorsion afin de las componentes de rotacion de la matriz homogenea
double L_HomogeneousMatrix::normalizaSVD(double *eMax, double *eMed, double *eMin)
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
	u.svd_ordenado(d, v);
	vT.transpOf(v);
	rot.OP_mult(u,vT);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			operator()(i,j)=rot(i,j);
	operator()(3,0) = operator()(3,1) = operator()(3,2) = 0; operator()(3,3) = 1;
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

double L_HomogeneousMatrix::normalizaSVD_corrigeJ(L_Matrix &J12xN)
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
			for (k=0; k<J12xN.lj; k++)
				J12xN(4*i+j,k) = J12xN(3*i+j,k) * (operator()(i,j)!=0 ? (rot(i,j) / operator()(i,j)) : 0);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			operator()(i,j)=rot(i,j);
	operator()(3,0) = operator()(3,1) = operator()(3,2) = 0; operator()(3,3) = 1;

	double e1, e2, e3, error;
	e1=log(L_ABS(d(0,0)));
	e2=log(L_ABS(d(1,0)));
	e3=log(L_ABS(d(2,0)));
	error=e1*e1+e2*e2+e3*e3;
	return error;
}

double L_HomogeneousMatrix::normaliza_mincuad()
{
	double matriz[9];
	std::vector<double> panTiltRoll(3);
	L_CoordsCart3D v1, v2, v3;
	L_LevenbergMarquardt optimizador;
	L_Pose3D pose;
	double err=0;

	v1.x=operator()(0,0);
	v1.y=operator()(1,0);
	v1.z=operator()(2,0);
	v2.x=operator()(0,1);
	v2.y=operator()(1,1);
	v2.z=operator()(2,1);
	v3.x=operator()(0,2);
	v3.y=operator()(1,2);
	v3.z=operator()(2,2);

	err+=fabs(log(v1.punto(v1))); // QUE CUATICO, DECIA ABS()
	err+=fabs(log(v2.punto(v2)));
	err+=fabs(log(v3.punto(v3)));
	err+=v1.punto(v2);
	err+=v1.punto(v3);
	err+=v2.punto(v3);

	// Calcular la matriz de rotacion mas cercana a la dada
	errorRot_prepara(matriz, panTiltRoll, true);
	optimizador.nIterationsMax=10;
	optimizador.errorToStop=1e-24;
	optimizador.lenVector=3;
	optimizador.epsilon=1.0e-7;
	optimizador.minimize(matriz, panTiltRoll, &errorRot);  // OJO REVISAR
	// Ahora, a corregir la pose
	pose=calcPose3D_fijo_a_movido();
	pose.ori.pan=panTiltRoll[0];
	pose.ori.tilt=panTiltRoll[1];
	pose.ori.roll=panTiltRoll[2];
	fijaPose3D_fijo_a_movido(pose);
	return err;
}

bool L_HomogeneousMatrix::fijaOrigenAdeArr(const L_CoordsCart3D &ori, const L_CoordsCart3D &dirAde, const L_CoordsCart3D &dirArr)
{
	L_CoordsCart3D ade, izq, arr;

	ade = dirAde;
	if (ade.normaliza_ret() == false)
		return false;
	arr = dirArr;
	if (arr.normaliza_ret() == false)
		return false;
	arr = arr - ade*arr.punto(ade);
	if (arr.normaliza_ret() == false)
		return false;
	izq = arr.cruz(ade);
		
	operator()(0,0)=ade.x;  operator()(0,1) = izq.x;  operator()(0,2)=arr.x;  operator()(0,3)=ori.x;
	operator()(1,0)=ade.y;  operator()(1,1) = izq.y;  operator()(1,2)=arr.y;  operator()(1,3)=ori.y;
	operator()(2,0)=ade.z;  operator()(2,1) = izq.z;  operator()(2,2)=arr.z;  operator()(2,3)=ori.z;
	operator()(3,0)=0.0;    operator()(3,1) = 0.0;    operator()(3,2)=0.0;    operator()(3,3)=1.0;
	return true;
}

void L_HomogeneousMatrix::cambiaEjes_a_AdeIzqArr()
{
	// Estan inicialmente en el sistema DerAbaAde
	L_HomogeneousMatrix H;

	H(0,0) = operator()(2,2);
	H(0,1) = -operator()(2,0);
	H(0,2) = -operator()(2,1);
	H(0,3) = operator()(2,3);
	H(1,0) = -operator()(0,2);
	H(1,1) = operator()(0,0);
	H(1,2) = operator()(0,1);
	H(1,3) = -operator()(0,3);
	H(2,0) = -operator()(1,2);
	H(2,1) = operator()(1,0);
	H(2,2) = operator()(1,1);
	H(2,3) = -operator()(1,3);
	H(3,0) = 0;
	H(3,1) = 0;
	H(3,2) = 0;
	H(3,3) = 1;

	*this = H;
}


void L_HomogeneousMatrix::cambiaEjes_a_DerAbaAde()
{
	L_HomogeneousMatrix H;

	H(0,0) = operator()(1,1);
	H(0,1) = operator()(1,2);
	H(0,2) = -operator()(1,0);
	H(0,3) = -operator()(1,3);
	H(1,0) = operator()(2,1);
	H(1,1) = operator()(2,2);
	H(1,2) = -operator()(2,0);
	H(1,3) = -operator()(2,3);
	H(2,0) = -operator()(0,1);
	H(2,1) = -operator()(0,2);
	H(2,2) = operator()(0,0);
	H(2,3) = operator()(0,3);
	H(3,0) = 0;
	H(3,1) = 0;
	H(3,2) = 0;
	H(3,3) = 1;

	*this = H;
}

void L_HomogeneousMatrix::cambiaEjes_transf(const L_HomogeneousMatrix &H, L_SistemaEjes3D inicial, L_SistemaEjes3D final)
{
	double m[3][3];
	double izq, arr, ade;
	// C:\Documents and Settings\Administrador\Mis documentos\Misdoc\Ramos\Tesis SLAM\Calculos\Cambios sistemas ejes v2.mws
	// (xyz)=R*(uvw), (XYZ)=R*(UVW)
	// (UVW)=M*(uvw) -> R*(UVW)=R*M*(uvw) -> R*(UVW)=R*M*R^-1*R*(uvw) -> (XYZ) = R*M*R^-1*(xyz)

	// Primera parte: (ade,izq,arr)=R*(x,y,z) -> Maia = R*Mxyz*R^-1
	switch(inicial)
	{
	case L_AdeIzqArr:
		m[0][0] = H(0,0); m[0][1] = H(0,1); m[0][2] = H(0,2);
		m[1][0] = H(1,0); m[1][1] = H(1,1); m[1][2] = H(1,2);
		m[2][0] = H(2,0); m[2][1] = H(2,1); m[2][2] = H(2,2);
		ade = H(0,3); izq = H(1,3); arr = H(2,3);
		break;
	case L_DerAbaAde:
		m[0][0] = H(2,2); m[0][1] = -H(2,0); m[0][2] = -H(2,1);
		m[1][0] = -H(0,2); m[1][1] = H(0,0); m[1][2] = H(0,1);
		m[2][0] = -H(1,2); m[2][1] = H(1,0); m[2][2] = H(1,1);
		ade = H(2,3); izq = -H(0,3); arr = -H(1,3);
		break;
	case L_DerArrAtr:
		m[0][0] = H(2,2); m[0][1] = H(2,0); m[0][2] = -H(2,1);
		m[1][0] = H(0,2); m[1][1] = H(0,0); m[1][2] = -H(0,1);
		m[2][0] = -H(1,2); m[2][1] = -H(1,0); m[2][2] = H(1,1);
		ade = -H(2,3); izq = -H(0,3); arr = H(1,3);
		break;
	case L_IzqArrAde:
		m[0][0] = H(2,2); m[0][1] = H(2,0); m[0][2] = H(2,1);
		m[1][0] = H(0,2); m[1][1] = H(0,0); m[1][2] = H(0,1);
		m[2][0] = H(1,2); m[2][1] = H(1,0); m[2][2] = H(1,1);
		ade = H(2,3); izq = H(0,3); arr = H(1,3);
		break;
	}

	// Segunda parte: (x,y,z)=R*(ade,izq,arr) -> Mxyz = R*Maia*R^-1
	switch(final)
	{
	case L_AdeIzqArr:
		operator()(0,0) = m[0][0]; operator()(0,1) = m[0][1]; operator()(0,2) = m[0][2];
		operator()(1,0) = m[1][0]; operator()(1,1) = m[1][1]; operator()(1,2) = m[1][2];
		operator()(2,0) = m[2][0]; operator()(2,1) = m[2][1]; operator()(2,2) = m[2][2];  // tenia error...
		operator()(0,3) = ade; operator()(1,3) = izq; operator()(2,3) = arr;
		break;
	case L_DerAbaAde:
		operator()(0,0) = m[1][1]; operator()(0,1) = m[1][2]; operator()(0,2) = -m[1][0];
		operator()(1,0) = m[2][1]; operator()(1,1) = m[2][2]; operator()(1,2) = -m[2][0];
		operator()(2,0) = -m[0][1]; operator()(2,1) = -m[0][2]; operator()(2,2) = m[0][0];
		operator()(0,3) = -izq; operator()(1,3) = -arr; operator()(2,3) = ade;
		break;
	case L_DerArrAtr:
		operator()(0,0) = m[1][1]; operator()(0,1) = -m[1][2]; operator()(0,2) = m[1][0];
		operator()(1,0) = -m[2][1]; operator()(1,1) = m[2][2]; operator()(1,2) = -m[2][0];
		operator()(2,0) = m[0][1]; operator()(2,1) = -m[0][2]; operator()(2,2) = m[0][0];
		operator()(0,3) = -izq; operator()(1,3) = arr; operator()(2,3) = -ade;
		break;
	case L_IzqArrAde:
		operator()(0,0) = m[1][1]; operator()(0,1) = m[1][2]; operator()(0,2) = m[1][0];
		operator()(1,0) = m[2][1]; operator()(1,1) = m[2][2]; operator()(1,2) = m[2][0];
		operator()(2,0) = m[0][1]; operator()(2,1) = m[0][2]; operator()(2,2) = m[0][0];
		operator()(0,3) = izq; operator()(1,3) = arr; operator()(2,3) = ade;
		break;
	}
}

void L_HomogeneousMatrix::cambiaEjes_pose(const L_HomogeneousMatrix &H, L_SistemaEjes3D inicial, L_SistemaEjes3D final)
{
	double m[3][3];
	double izq, arr, ade;
	// C:\Documents and Settings\Administrador\Mis documentos\Misdoc\Ramos\Tesis SLAM\Calculos\Cambios sistemas ejes v2.mws
	// (xyz)=R*(uvw), (XYZ)=R*(UVW)
	// (UVW)=M*(uvw) -> R*(UVW)=R*M*(uvw) -> R*(UVW)=R*M*R^-1*R*(uvw) -> (XYZ) = R*M*R^-1*(xyz)

	// Primera parte: (ade,izq,arr)=R*(x,y,z) -> Maia = R*Mxyz*R^-1
	switch(inicial)
	{
	case L_AdeIzqArr:
		m[0][0] = H(0,0); m[0][1] = H(0,1); m[0][2] = H(0,2);
		m[1][0] = H(1,0); m[1][1] = H(1,1); m[1][2] = H(1,2);
		m[2][0] = H(2,0); m[2][1] = H(2,1); m[2][2] = H(2,2);
		ade = H(0,3); izq = H(1,3); arr = H(2,3);
		break;
	case L_DerAbaAde:
		m[0][0] = H(2,0); m[0][1] = H(2,1); m[0][2] = H(2,2);
		m[1][0] = -H(0,0); m[1][1] = -H(0,1); m[1][2] = -H(0,2);
		m[2][0] = -H(1,0); m[2][1] = -H(1,1); m[2][2] = -H(1,2);
		ade = H(2,3); izq = -H(0,3); arr = -H(1,3);
		break;
	case L_DerArrAtr:
		m[0][0] = -H(2,0); m[0][1] = -H(2,1); m[0][2] = -H(2,2);
		m[1][0] = -H(0,0); m[1][1] = -H(2,1); m[1][2] = -H(2,2);
		m[2][0] = H(1,0); m[2][1] = H(1,1); m[2][2] = H(1,2);
		ade = -H(2,3); izq = -H(0,3); arr = H(1,3);
		break;
	case L_IzqArrAde:
		m[0][0] = H(2,0); m[0][1] = H(2,1); m[0][2] = H(2,2);
		m[1][0] = H(0,0); m[1][1] = H(0,1); m[1][2] = H(0,2);
		m[2][0] = H(1,0); m[2][1] = H(1,1); m[2][2] = H(1,2);
		ade = H(2,3); izq = H(0,3); arr = H(1,3);
		break;
	}

	// Segunda parte: (x,y,z)=R*(ade,izq,arr) -> Mxyz = R*Maia*R^-1
	switch(final)
	{
	case L_AdeIzqArr:
		operator()(0,0) = m[0][0]; operator()(0,1) = m[0][1]; operator()(0,2) = m[0][2];
		operator()(1,0) = m[1][0]; operator()(1,1) = m[1][1]; operator()(1,2) = m[1][2];
		operator()(2,0) = m[2][0]; operator()(2,1) = m[2][1]; operator()(2,2) = m[2][2];  // tenia error...
		operator()(0,3) = ade; operator()(1,3) = izq; operator()(2,3) = arr;
		break;
	case L_DerAbaAde:
		operator()(0,0) = -m[1][0]; operator()(0,1) = -m[1][1]; operator()(0,2) = -m[1][2];
		operator()(1,0) = m[2][0]; operator()(1,1) = m[2][1]; operator()(1,2) = m[2][2];
		operator()(2,0) = -m[0][0]; operator()(2,1) = -m[0][1]; operator()(2,2) = -m[0][2];
		operator()(0,3) = -izq; operator()(1,3) = -arr; operator()(2,3) = ade;
		break;
	case L_DerArrAtr:
		operator()(0,0) = -m[1][0]; operator()(0,1) = -m[1][1]; operator()(0,2) = -m[1][2];
		operator()(1,0) = m[2][0]; operator()(1,1) = m[2][1]; operator()(1,2) = m[2][2];
		operator()(2,0) = -m[0][0]; operator()(2,1) = -m[0][1]; operator()(2,2) = -m[0][2];
		operator()(0,3) = -izq; operator()(1,3) = arr; operator()(2,3) = -ade;
		break;
	case L_IzqArrAde:
		operator()(0,0) = m[1][0]; operator()(0,1) = m[1][1]; operator()(0,2) = m[1][2];
		operator()(1,0) = m[0][0]; operator()(1,1) = m[0][1]; operator()(1,2) = m[0][2];
		operator()(2,0) = m[1][0]; operator()(2,1) = m[1][1]; operator()(2,2) = m[1][2];
		operator()(0,3) = izq; operator()(1,3) = arr; operator()(2,3) = ade;
		break;
	}
}

void L_HomogeneousMatrix::test()
{
	L_Pose3D_cuat p;
	L_HomogeneousMatrix H1, H2, HM, Hm;
	L_Matrix J(12,12), J2(12,12);
	int i, j;
	double delta = 1.0e-3;
	p.define(L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10), L_RANDOM(-10,10));
	p.normaliza();
	H1.fijaPose3D_cuat(p);

	jacob_inversaTrasponiendoDe(H1, J);
	J.print("J");

	for (j=0; j<12; j++)
	{
		H1(j/4,j%4) += delta;
		HM.inversaTrasponiendoDe(H1);
		H1(j/4,j%4) -= 2*delta;
		Hm.inversaTrasponiendoDe(H1);
		H1(j/4,j%4) += delta;

		for (i=0; i<12; i++)
			J2(i,j) = (HM(i/4,i%4) - Hm(i/4,i%4)) / (2*delta);
	}
	J2.print("J2");
	printf("Parecido: %f\n", J.cosineDistanceTo(J2));
	getchar();

	jacob_inversaAproxDe(H1, J);
	J.print("J");

	for (j=0; j<12; j++)
	{
		H1(j/4,j%4) += delta;
		HM.inverseOf(H1);
		H1(j/4,j%4) -= 2*delta;
		Hm.inverseOf(H1);
		H1(j/4,j%4) += delta;

		for (i=0; i<12; i++)
			J2(i,j) = (HM(i/4,i%4) - Hm(i/4,i%4)) / (2*delta);
	}
	J2.print("J2");
	printf("Parecido: %f\n", J.cosineDistanceTo(J2));
	getchar();
}


void L_Ellipsoid3D::setCovarianceAndTranslation(const double *cov3x3, const double *t3x1, double chi2)
{
	// (x-t)T * C * (x-t) = 1
	// (x-t)T * RT * ST * S * R  * (x-t) = 1
	// (x-t)T   U  *    D   * VT * (x-t) = 1


	// u = R*x - R*t = H*x  // del sistema global (x) al de la elipse(u)
	// v = S*u;   // del sistema de la elipse(u) al sistema del circulo (v)

	// vT*v = 1;  // ecuacion del circulo

	// u = S^-1 * v;  // del sistema del circulo al de la elipse
	// x = RTu + t;   // del sistema de la elipse al global

	L_Matrix u(3,3), d(3,1), v(3,3);
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			u(i,j) = cov3x3[i*3+j];
	u.svd(d,v);
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			H(i,j) = v(j,i);  // R = VT
	H(0,3) = -H(0,0)*t3x1[0] - H(0,1)*t3x1[1] - H(0,2)*t3x1[2];
	H(1,3) = -H(1,0)*t3x1[0] - H(1,1)*t3x1[1] - H(1,2)*t3x1[2];
	H(2,3) = -H(2,0)*t3x1[0] - H(2,1)*t3x1[1] - H(2,2)*t3x1[2];
	sx = sqrt(d(0,0)*chi2);
	sy = sqrt(d(1,0)*chi2);
	sz = sqrt(d(2,0)*chi2);
}


void L_ConjuntoRectas3D::agregarVertices(double x1, double y1, double x2, double y2)
{
	vertices.resize(vertices.size()+2);
	vertices[vertices.size()-2].define(x1,y1,0);
	vertices[vertices.size()-1].define(x2,y2,0);
	a.resize(a.size()+1);
	b.resize(b.size()+1);
	c.resize(c.size()+1);
	a[a.size()-1] = y2 - y1;
	b[b.size()-1] = x1 - x2;
	c[c.size()-1] = x2*y1 - x1*y2;
	double r = sqrt(a[a.size()-1]*a[a.size()-1] + b[b.size()-1]*b[b.size()-1]);
	a[a.size()-1] /= r;
	b[a.size()-1] /= r;
	c[a.size()-1] /= r;
}

void L_ConjuntoRectas3D::agregarElipseVertices(double lx, double ly, int nVertices, double theta0, double sentido)
{
	int i;
	double theta1, theta2;
	double x1, y1, x2, y2;
	double c1, s1, c2, s2;
	for (i=0; i<nVertices; i++)
	{
		if (sentido > 0)
		{
			theta1 = i * 2*M_PI / nVertices + theta0;
			theta2 = (i+1) * 2*M_PI / nVertices + theta0;
		}
		else
		{
			theta1 = -i * 2*M_PI / nVertices + theta0;
			theta2 = -(i+1) * 2*M_PI / nVertices + theta0;
		}
		c1 = cos(theta1);
		s1 = sin(theta1);
		c2 = cos(theta2);
		s2 = sin(theta2);
		x1 = lx*c1/4;
		y1 = ly*s1/4;
		x2 = lx*c2/4;
		y2 = ly*s2/4;
		agregarVertices(x1, y1, x2, y2);
	}
}



void L_ConjuntoRectas3D::generarAsociacionesSegmentos()
{
	int k=0, nsegmt;
	pSegm.resize(p.size());
	nSegm.resize(p.size());
	nsegmt = (int)segmentos.size()/2;

	// Revisar para cada punto a cual segmento belongTo
	// Los puntos que no estan asociados a ningun segmento se ignoran
	for (int i=0; i<(int)p.size(); i++)
	{
		for (int j=0; j < nsegmt; j++)
		{
			if (i >= segmentos[2*j+0] && i <= segmentos[2*j+1])
			{
				pSegm[k] = p[i];
				nSegm[k] = j;
				k++;
			}
		}
	}
	nSegm.resize(k);
	pSegm.resize(k);
}

double L_ConjuntoRectas3D::calcularMejorRectaPorPunto(std::vector<L_CoordsCart3D> &puntosMov, std::vector<int> &mejorRecta)
{
	double ex, ey, e2, exMin, eyMin, e2Min;
	double errAcum = 0;
	int i, k, kMin;
	mejorRecta.resize(puntosMov.size());

	for (i=0; i<(int)puntosMov.size(); i++)
	{
		exMin = 0;
		eyMin = 0;
		e2Min = 1.0e60;
		kMin = 0;
		// Probar cual de las rectas es mas cercana al punto
		for (k=0; k<(int)a.size(); k++)
		{
			ex = a[k]*puntosMov[i].x + b[k]*puntosMov[i].y + c[k];
			ey = puntosMov[i].z;
			e2 = ex*ex + ey*ey;
			if (e2 < e2Min)
			{
				kMin = k;
				exMin = ex;
				eyMin = ey;
				e2Min = e2;
			}
		}
		mejorRecta[i] = kMin;
		//errAcum += sqrt(e2Min);
		errAcum += fabs(exMin);
	}
	return errAcum;
}

void L_ConjuntoRectas3D::listaVideos(int numeroVideo)
{
	if (numeroVideo >= 60 && numeroVideo < 80)
	{
		agregarVertices(2,0,  3,0);   // * * *
		agregarVertices(3,0,  3,-2);  // *   *
		agregarVertices(3,-2, 1,-2);  // *   *
		agregarVertices(1,-2, 0,-1);  // * *
		if (numeroVideo == 70)
		{
			agregarSegmento(0, 200);
			agregarSegmento(258, 409);
			agregarSegmento(465, 690);
			agregarSegmento(749, 867);

		}
	}
	if (numeroVideo >= 110)
	{
		agregarVertices(0,0,  3,0);   //  Figura de 3x2 pentagonal
		agregarVertices(3,0,  3,-2);  // * * 
		agregarVertices(3,-2, 0,-2);  // * * 
		agregarVertices(0,-2, 0,0);   // * 
		if (numeroVideo == 101)
		{
			agregarSegmento(0, 200);
			agregarSegmento(258, 409);
			agregarSegmento(465, 690);
			agregarSegmento(749, 867);
		}
	}
	if (numeroVideo >= 110 && numeroVideo < 300)
	{
		agregarVertices(0,0,  3,0);   //  Rectangulo de 3x2
		agregarVertices(3,0,  3,-2);  // * *
		agregarVertices(3,-2, 0,-2);  // * *
		agregarVertices(0,-2, 0,0);   // * *
		if (numeroVideo == 101)
		{
			agregarSegmento(0, 200);
			agregarSegmento(258, 409);
			agregarSegmento(465, 690);
			agregarSegmento(749, 867);
		}
		if (numeroVideo == 165)
		{
			agregarSegmento(0, 200);
			agregarSegmento(280, 430);
			agregarSegmento(490, 700);
			agregarSegmento(790, 952);
		}
		if (numeroVideo == 166)
		{
			agregarSegmento(0, 190);  // 390 = 1/3 de arriba
			agregarSegmento(290, 455);
			agregarSegmento(505, 735); // 505 inventado
			agregarSegmento(843, 1007);
		}
		if (numeroVideo == 167)
		{
			// Erroneo
		}
		if (numeroVideo == 168)
		{
			//agregarSegmento(0, 225);
			//agregarSegmento(, 495); ... se interrumpio, iba malo
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
		if (numeroVideo == 169)
		{
			agregarSegmento(0, 220); // 0 220 conocidos
			agregarSegmento(280, 529); // 559 conocido
			agregarSegmento(560, 720);
			agregarSegmento(840, 2000);
		}
		if (numeroVideo == 170)
		{
			agregarSegmento(0, 204);
			agregarSegmento(285, 445);
			agregarSegmento(525, 730);
			agregarSegmento(820, 985);
		}
		if (numeroVideo == 171)
		{
			//agregarSegmento(0, );
			//agregarSegmento(300, 440);
			//agregarSegmento(, 750);
			//agregarSegmento(, 1058);
		}
		if (numeroVideo == 172)
		{
		}
		//
		//
		if (numeroVideo == 251) // Se pierde, apunta mal
		{
			//agregarSegmento(0, 240);
			//agregarSegmento(, 280);
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
		if (numeroVideo == 252) // Jugoso
		{
			//agregarSegmento(0, 240);
			//agregarSegmento(, 280);
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
		if (numeroVideo == 253) // Jugoso
		{
			//agregarSegmento(0, );
			//agregarSegmento(, );
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
		if (numeroVideo == 261) // Jugoso
		{
			//agregarSegmento(0, );
			//agregarSegmento(, );
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
		//
		//
		if (numeroVideo == 273) // Jugoso
		{
			//agregarSegmento(0, );
			//agregarSegmento(, );
			//agregarSegmento(, );
			//agregarSegmento(, );
		}
	}

	if (numeroVideo >= 310 && numeroVideo < 400)
	{
		double lx = 2;
		double ly = 3;
		agregarElipseVertices(lx, ly, 8, 90*M_PI/180, -1);
		if (numeroVideo == 310) // 1
		{
			agregarSegmento(0, 70);
			agregarSegmento(70, 145);
			agregarSegmento(145, 315);
			agregarSegmento(315, 388);
			agregarSegmento(388, 464);
			agregarSegmento(464, 570);
			agregarSegmento(570, 660);
			agregarSegmento(660, 696);
		}
		if (numeroVideo == 311) // 2
		{
			agregarSegmento(0, 80);
			agregarSegmento(80, 160);
			agregarSegmento(160, 245);
			agregarSegmento(245, 320);
			agregarSegmento(320, 420);
			agregarSegmento(420, 500);
			agregarSegmento(500, 564);
			agregarSegmento(564, 610);
		}
		if (numeroVideo == 312) // 3
		{
			agregarSegmento(0, 150);
			agregarSegmento(150, 180);
			agregarSegmento(180, 286);
			agregarSegmento(286, 350);
			agregarSegmento(350, 499);
			agregarSegmento(499, 540);
			agregarSegmento(540, 650);
			agregarSegmento(650, 690);
		}
		if (numeroVideo == 313)  // Este es bueno y no lo medi  // 4
		{
			agregarSegmento(0, 116);
			agregarSegmento(116, 175);
			agregarSegmento(175, 276);
			agregarSegmento(276, 345);
			agregarSegmento(345, 440);
			agregarSegmento(440, 500);
			agregarSegmento(500, 580); // ?
			agregarSegmento(580, 650); // ?
		}
		if (numeroVideo == 315) // 5
		{
			agregarSegmento(0, 130);
			agregarSegmento(130, 190);
			agregarSegmento(190, 315);
			agregarSegmento(315, 390); // ?
			agregarSegmento(390, 460);
			agregarSegmento(460, 510);
			agregarSegmento(510, 595);
			agregarSegmento(595, 650); // ?
		}
		if (numeroVideo == 316) // 6
		{
			agregarSegmento(0, 90);
			agregarSegmento(90, 170);
			agregarSegmento(170, 320);
			agregarSegmento(320, 340);
			agregarSegmento(340, 510);
			agregarSegmento(510, 550);
			agregarSegmento(550, 617);
			agregarSegmento(617, 850); // ?
		}
		if (numeroVideo == 324) // 7
		{
			agregarSegmento(0, 130);
			agregarSegmento(130, 190);
			agregarSegmento(190, 274);
			agregarSegmento(274, 349);
			agregarSegmento(349, 506);
			agregarSegmento(506, 564);
			agregarSegmento(564, 660);
			agregarSegmento(660, 850); // ?
		}
	}

	// Para nuevo trozo
			//agregarSegmento(0, );
			//agregarSegmento(, );
			//agregarSegmento(, );
			//agregarSegmento(, );
	//
}

void L_ConjuntoRectas3D::calcularErrorGroundTruth(const char *nomArchM)
{
	std::vector<L_CoordsCart3D> ruta;
	std::vector<L_CoordsCart3D> rutaAlin;
	std::vector<int> recta;
	L_RotEscTrasl3D tr;
	int i;
	double errAcum;

	if (a.size() == 0)
	{
		printf("L_ConjuntoRectas3D::calcularErrorGroundTruth() : poligono no definido\nPresione ENTER\n");
		getchar();
		return;
	}

	// Leer archivo con los puntos (son puntos+cuaternion)
	FILE *fp;
	fp = fopen(nomArchM, "r");
	if (fp == NULL)
		{printf("No se puede abrir %s.\nPresione ENTER para continuar\n", nomArchM); getchar(); return;}

	int rs;
	if (fscanf(fp, "%d", &rs) != 1)
		{printf("L_ConjuntoRectas3D::calcularErrorGroundTruth() : formato de archivo incorrecto.\nPresione ENTER para continuar\n"); getchar(); return;}

	ruta.resize(rs);

	for (i=0; i<(int)ruta.size(); i++)
		if (fscanf(fp, "%lf%lf%lf%*f%*f%*f%*f", &ruta[i].x, &ruta[i].y, &ruta[i].z) != 3)
			{printf("L_ConjuntoRectas3D::calcularErrorGroundTruth() : formato de archivo incorrecto.\nPresione ENTER para continuar\n"); getchar(); return;}

	fclose(fp);
	
	// Definir las rectas del modelo poligonal
	p = ruta;

	generarAsociacionesSegmentos();

	tr.calcTransfPoligono_1(*this, 60);  // Mejor recta dados los puntos con asociacion conocida
	//tr.calcTransfPoligono_2(*this, 60);  // Mejor recta dados todos los puntos

	rutaAlin.resize(ruta.size());
	for (i=0; i<(int)ruta.size(); i++)
		rutaAlin[i] = tr.transformar(ruta[i]);

	errAcum = calcularMejorRectaPorPunto(rutaAlin, recta);

	FILE *fp2 = fopen("psal.txt", "w");
	L_CoordsCart3D vmov;
	vmov = tr.transformarInv(vertices[0]);
	fprintf(fp2, "%g %g %g\n", vmov.x, vmov.y, vmov.z);
	vmov = tr.transformarInv(vertices[2]);
	fprintf(fp2, "%g %g %g\n", vmov.x, vmov.y, vmov.z);
	vmov = tr.transformarInv(vertices[4]);
	fprintf(fp2, "%g %g %g\n", vmov.x, vmov.y, vmov.z);
	vmov = tr.transformarInv(vertices[6]);
	fprintf(fp2, "%g %g %g\n", vmov.x, vmov.y, vmov.z);
	vmov = tr.transformarInv(vertices[7]);
	fprintf(fp2, "%g %g %g\n", vmov.x, vmov.y, vmov.z);
	for (i=0; i<(int)rutaAlin.size(); i++)
	{
		double A = a[recta[i]], B = b[recta[i]], C = c[recta[i]];
		double X = rutaAlin[i].x, Y = rutaAlin[i].y, Z = rutaAlin[i].z;
		double exy = fabs(A*X + B*Y + C);
		fprintf(fp2, "%g    %g %g %g    %g %g\n", exy, X, Y, Z, (B!=0) ? X : -(B*Y+C)/A, (B!=0) ? -(A*C+C)/B : Y);
	}
	fclose(fp2);

	L_CoordsCart3D despl = rutaAlin[0] - rutaAlin[rutaAlin.size()-1];
	double errIniFin = sqrt(despl.x*despl.x + despl.y*despl.y);
	printf("Error acumulado: %f\n", errAcum);
	printf("Error mean: %f\n", errAcum / ruta.size());
	printf("Distancia entre punto inicial y final: %f\n", errIniFin);
	printf("Presione ENTER para continuar\n");
	getchar();
}

// Calcula un error de los puntos respecto a un conjunto de rectas
// Sirve para encontrar la mejor transformacion de un conjunto de puntos a un poligono
void L_ConjuntoRectas3D::calcError1(const void *obj, const L_Matrix &x, L_Matrix &err)
{
	const L_ConjuntoRectas3D *dat = (const L_ConjuntoRectas3D *)obj;
	L_CoordsCart3D res;
	L_RotEscTrasl3D tr;
	int i, ns;
	tr.pcuat.pos.x = x(0,0);
	tr.pcuat.pos.y = x(0,1);
	tr.pcuat.pos.z = x(0,2);
	tr.pcuat.ori.a = x(0,3);
	tr.pcuat.ori.b = x(0,4);
	tr.pcuat.ori.c = x(0,5);
	tr.pcuat.ori.d = x(0,6);
	tr.escala = x(0,7);

	// Solamente algunos puntos se consideran en esta minimizacion
	err.reallocate(2*(int)dat->pSegm.size(), 1);
	for (i=0; i<(int)dat->pSegm.size(); i++)
	{
		res = tr.transformar(dat->pSegm[i]);
		ns = dat->nSegm[i];
		err(2*i+0,0) = dat->a[ns]*res.x + dat->b[ns]*res.y + dat->c[ns];
		err(2*i+1,0) = res.z;
	}
}


// Calcula un error de los puntos respecto a un conjunto de rectas
// Sirve para encontrar la mejor transformacion de un conjunto de puntos a un poligono
// Esta funcion de error es considerablemente no lineal y puede requierir muchas iteraciones
void L_ConjuntoRectas3D::calcError2(const void *obj, const L_Matrix &x, L_Matrix &err)
{
	const L_ConjuntoRectas3D *dat = (const L_ConjuntoRectas3D *)obj;
	L_CoordsCart3D res;
	L_RotEscTrasl3D tr;
	int i;
	tr.pcuat.pos.x = x(0,0);
	tr.pcuat.pos.y = x(0,1);
	tr.pcuat.pos.z = x(0,2);
	tr.pcuat.ori.a = x(0,3);
	tr.pcuat.ori.b = x(0,4);
	tr.pcuat.ori.c = x(0,5);
	tr.pcuat.ori.d = x(0,6);
	tr.escala = x(0,7);

	// Se consideran todos los puntos pero sin asociacion conocida
	double ex, ey, e2, exMin, eyMin, e2Min;
	int k, kMin;

	err.reallocate(2*(int)dat->p.size(), 1);
	for (i=0; i<(int)dat->p.size(); i++)
	{
		exMin = 0;
		eyMin = 0;
		e2Min = 1.0e60;
		kMin = 0;
		res = tr.transformar(dat->p[i]);
		// Probar cual de las rectas es mas cercana al punto
		for (k=0; k<(int)dat->a.size(); k++)
		{
			ex = dat->a[k]*res.x + dat->b[k]*res.y + dat->c[k];
			ey = res.z;
			e2 = ex*ex + ey*ey;
			if (e2 < e2Min)
			{
				kMin = k;
				exMin = ex;
				eyMin = ey;
				e2Min = e2;
			}
		}
		// Calcular vector de error para Levenberg-Marquardt; usar root cuadrada con signo para MAE
		if (exMin > 0)
			err(2*i+0,0) = sqrt(exMin);
		else
			err(2*i+0,0) = -sqrt(-exMin);
		if (eyMin > 0)
			err(2*i+1,0) = sqrt(eyMin);
		else
			err(2*i+1,0) = -sqrt(-eyMin);
	}
}


void L_RotEscTrasl3D::err_mf_vect(const void *pose3D_cuat_nube_pesos, const L_Matrix &pcuatesc, L_Matrix &errores)
{
	// Supuestamente debe ser rapida
	const L_Pose3D_cuat_nube_puntos_pesos *nube = (const L_Pose3D_cuat_nube_puntos_pesos *)pose3D_cuat_nube_pesos;
	L_RotEscTrasl3D p;
	const L_CoordsCart3D *m, *f;
	L_CoordsCart3D pun, despl;
	int i;
	errores.reallocate((int)nube->fijo->size()*3, 1);
	p.pcuat.pos.x = pcuatesc(0,0); p.pcuat.pos.y = pcuatesc(1,0); p.pcuat.pos.z = pcuatesc(2,0);
	p.pcuat.ori.a = pcuatesc(3,0); p.pcuat.ori.b = pcuatesc(4,0); p.pcuat.ori.c = pcuatesc(5,0); p.pcuat.ori.d = pcuatesc(6,0);
	p.escala = pcuatesc(7,0);
	if (p.pcuat.ori.normaliza_ret() == false)
		p.pcuat.ori.fijaRotacionCero();

	for (i=0; i<(int)nube->fijo->size(); i++)
	{
		m = &nube->movido->operator[](i);
		f = &nube->fijo->operator[](i);
		pun = p.pcuat.ori.rotarVector(*m);
		pun = pun * p.escala;
		pun = pun + p.pcuat.pos;
		despl = pun - (*f);
		errores(3*i+0,0) = despl.x * (*nube->pesos)[i]; // ERROR: decia [3*i+0][1]
		errores(3*i+1,0) = despl.y * (*nube->pesos)[i];
		errores(3*i+2,0) = despl.z * (*nube->pesos)[i];
	}
}

void L_RotEscTrasl3D::encontrarTransfMinCuad(const std::vector<L_CoordsCart3D> &orig, const std::vector<L_CoordsCart3D> &res)
{
	L_Matrix pose;
	L_LevenbergMarquardt optim, *opt=NULL;
	L_Pose3D_cuat_nube_puntos_pesos nube;
	std::vector<double> pesos(orig.size());
	int i;

	for (i=0; i<(int)orig.size(); i++)
		pesos[i] = 1;

	pcuat.pos.fijaCero();
	pcuat.ori.fijaRotacionCero();
	escala = 1;

	nube.movido = &orig;
	nube.fijo = &res;

	// Parametros del optimizador
	optim.epsilon = 1e-7;
	optim.nIterationsMax = 8;
	optim.errorToStop = 0;
	optim.lenVector = 8;
	
	if (opt == NULL)
		opt = &optim;

	//opt.useHessian = true; // Solo se usa en la minimizacion con error escalar
	pose.reallocate(8,1);
	pose(0,0) = pcuat.pos.x;
	pose(1,0) = pcuat.pos.y;
	pose(2,0) = pcuat.pos.z;
	pose(3,0) = pcuat.ori.a;
	pose(4,0) = pcuat.ori.b;
	pose(5,0) = pcuat.ori.c;
	pose(6,0) = pcuat.ori.d;
	pose(7,0) = escala;
	opt->minimize_vect((const void *)&nube, pose, &err_mf_vect);
	pcuat.pos.x = pose(0,0);
	pcuat.pos.y = pose(1,0);
	pcuat.pos.z = pose(2,0);
	pcuat.ori.a = pose(3,0);
	pcuat.ori.b = pose(4,0);
	pcuat.ori.c = pose(5,0);
	pcuat.ori.d = pose(6,0);
	escala      = pose(7,0);
	if (pcuat.ori.normaliza_ret() == false)
		pcuat.ori.fijaRotacionCero();		
}


double L_RotEscTrasl3D::encontrarTransfMinCuadTrayectorias(const std::vector<L_CoordsCart3D> &orig, const std::vector<L_CoordsCart3D> &res)
{
	L_LevenbergMarquardt optim, *opt=NULL;
	L_Pose3D_cuat_nube_puntos_pesos nube;
	L_RotEscTrasl3D v2;
	std::vector<double> pesos(orig.size());
	L_Matrix pose;
	double errMin;
	int i;

	pesos[0] = 0; // Desconocido...
	for (i=1; i<(int)orig.size(); i++)
		pesos[i-1] = res[i-1].distanciaA(res[i]);

	pcuat.pos.fijaCero();
	pcuat.ori.fijaRotacionCero();
	escala = 1;

	v2.pcuat.pos.fijaCero();
	v2.pcuat.ori.fijaRotacionCero();
	v2.escala = 1;

	// Encontrar transf inicial para asegurar convergencia... tiene escala incorrecta pero es mejor que nada
	v2.calcPose_movido_a_fijo_matriz(orig, res);

	nube.movido = &orig;
	nube.fijo = &res;
	nube.pesos = &pesos;

	// Parametros del optimizador
	optim.epsilon = 1e-7;
	optim.nIterationsMax = 15;  // Esta parte no demora tanto en realidad...
	optim.errorToStop = 1e-10;
	optim.lenVector = 8;
	
	if (opt == NULL)
		opt = &optim;

	// Ver si v2 es mejor que "zero"
	{
		L_Matrix v(8,1), err(3*(int)orig.size(), 1);
		double err0, err1, err2;
		for (i=0; i<8; i++)
			v(i,0) = (i==3 || i==8) ? 1 : 0;
		err_mf_vect((void *)&nube, v, err);
		err0 = err.normL2();
		for (i=0; i<8; i++)
			v(i,0) = el(i);
		err_mf_vect((void *)&nube, v, err);
		err1 = err.normL2();
		for (i=0; i<8; i++)
			v(i,0) = v2.el(i);
		err_mf_vect((void *)&nube, v, err);
		err2 = err.normL2();
		if (err0 < err1 && err0 < err2)
			fijaCero();
		else if (err1 < err2)
			; // Ya esta bien
		else
			*this = v2;
	}


	//opt.useHessian = true; // Solo se usa en la minimizacion con error escalar
	pose.reallocate(8,1);
	pose(0,0) = pcuat.pos.x;
	pose(1,0) = pcuat.pos.y;
	pose(2,0) = pcuat.pos.z;
	pose(3,0) = pcuat.ori.a;
	pose(4,0) = pcuat.ori.b;
	pose(5,0) = pcuat.ori.c;
	pose(6,0) = pcuat.ori.d;
	pose(7,0) = escala;
	errMin = opt->minimize_vect((const void *)&nube, pose, &err_mf_vect);
	pcuat.pos.x = pose(0,0);
	pcuat.pos.y = pose(1,0);
	pcuat.pos.z = pose(2,0);
	pcuat.ori.a = pose(3,0);
	pcuat.ori.b = pose(4,0);
	pcuat.ori.c = pose(5,0);
	pcuat.ori.d = pose(6,0);
	escala      = pose(7,0);
	if (pcuat.ori.normaliza_ret() == false)
		pcuat.ori.fijaRotacionCero();
	return errMin;
}


void L_RotEscTrasl3D::calcPose_movido_a_fijo_matriz(const std::vector<L_CoordsCart3D> &movido, const std::vector<L_CoordsCart3D> &fijo)
{
	L_Matrix A, b, AT, ATA, ATAinv, ATb, X;
	L_HomogeneousMatrix H;
	double emax, emed, emin;
	int i;
	A.reallocate(3*(int)movido.size(), 12);
	b.reallocate(3*(int)movido.size(), 1);
	A.setZero();

	for (i=0; i<(int)movido.size(); i++)
	{
		A(3*i+0,0)=movido[i].x;
		A(3*i+0,1)=movido[i].y;
		A(3*i+0,2)=movido[i].z;
		A(3*i+0,3)=1;

		A(3*i+1,4)=movido[i].x;
		A(3*i+1,5)=movido[i].y;
		A(3*i+1,6)=movido[i].z;
		A(3*i+1,7)=1;

		A(3*i+2,8)=movido[i].x;
		A(3*i+2,9)=movido[i].y;
		A(3*i+2,10)=movido[i].z;
		A(3*i+2,11)=1;

		b(3*i+0,0) = fijo[i].x;
		b(3*i+1,0) = fijo[i].y;
		b(3*i+2,0) = fijo[i].z;
	}
	AT.transpOf(A);
	ATA.OP_mult(AT, A);
	ATAinv.inverseOf(ATA);
	ATb.OP_mult(AT, b);
	X.OP_mult(ATAinv,ATb);
	for (i=0; i<12; i++)
		H(i/4,i%4) = X(i,0);
	H.normalizaSVD(&emax, &emed, &emin); // Puede que un val propio de 0 pq no hay info en un eje...
	pcuat = H.calcPose3D_cuat();
	escala = emed;
}


void L_RotEscTrasl3D::pideMatriz4x4(L_Matrix &m)
{
	L_HomogeneousMatrix H;
	m.reallocate(4,4);
	H.fijaPose3D_cuat(pcuat);
	m.identity();
	m(0,0) = escala*H(0,0);
	m(0,1) = escala*H(0,1);
	m(0,2) = escala*H(0,2);
	m(1,0) = escala*H(1,0);
	m(1,1) = escala*H(1,1);
	m(1,2) = escala*H(1,2);
	m(2,0) = escala*H(2,0);
	m(2,1) = escala*H(2,1);
	m(2,2) = escala*H(2,2);
	m(0,3) = H(0,3);
	m(1,3) = H(1,3);
	m(2,3) = H(2,3);
}

void L_RotEscTrasl3D::calcTransfPoligono_1(const L_ConjuntoRectas3D &rectas, int niter)
{
	escala = 1.0;
	L_LevenbergMarquardt opt;
	L_Matrix x(8,1);
	x(0,0) = 0;
	x(0,1) = 0;
	x(0,2) = 0;
	x(0,3) = 1;
	x(0,4) = 0;
	x(0,5) = 0;
	x(0,6) = 0;
	x(0,7) = 1;
	opt.nIterationsMax = niter;
	opt.lenVector = 8;
	opt.epsilon = 1.0e-5;
	opt.minimize_vect(&rectas, x, &L_ConjuntoRectas3D::calcError1);
	pcuat.pos.x = x(0,0);
	pcuat.pos.y = x(1,0);
	pcuat.pos.z = x(2,0);
	pcuat.ori.a = x(3,0);
	pcuat.ori.b = x(4,0);
	pcuat.ori.c = x(5,0);
	pcuat.ori.d = x(6,0);
	escala = x(7,0);
}

void L_RotEscTrasl3D::calcTransfPoligono_2(const L_ConjuntoRectas3D &rectas, int niter)
{
	escala = 1.0;
	L_LevenbergMarquardt opt;
	L_Matrix x(8,1);
	x(0,0) = pcuat.pos.x;
	x(0,1) = pcuat.pos.y;
	x(0,2) = pcuat.pos.z;
	x(0,3) = pcuat.ori.a;
	x(0,4) = pcuat.ori.b;
	x(0,5) = pcuat.ori.c;
	x(0,6) = pcuat.ori.d;
	x(0,7) = escala;
	opt.nIterationsMax = niter;
	opt.lenVector = 8;
	opt.epsilon = 1.0e-5;
	opt.minimize_vect(&rectas, x, &L_ConjuntoRectas3D::calcError2);
	pcuat.pos.x = x(0,0);
	pcuat.pos.y = x(1,0);
	pcuat.pos.z = x(2,0);
	pcuat.ori.a = x(3,0);
	pcuat.ori.b = x(4,0);
	pcuat.ori.c = x(5,0);
	pcuat.ori.d = x(6,0);
	escala = x(7,0);
}

bool L_CalculadorPoseTriangulo::metodoFinsterwalder()
{
	// Entradas: q1, q2, q3, a, b, c
	// Salidas: s1, s2, s3
	double a2, b2, c2;
	double cosAlfa, cosBeta, cosGamma;
	double cos2Alfa;
	double sin2Alfa, sin2Beta, sin2Gamma;
	double G, H, I, J;
	L_ComplexDouble lambda0, lambda1, lambda2;

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
	sin2Alfa = 1 - cos2Alfa;
	sin2Beta = 1 - cosBeta*cosBeta;
	sin2Gamma = 1 - cosGamma*cosGamma;

	G = c2*(c2*sin2Beta - b2*sin2Gamma);
	H = b2*(b2 - a2)*sin2Gamma + c2*(c2 + 2*a2)*sin2Beta + 2*b2*c2*(-1 + cosAlfa*cosBeta*cosGamma);
	I = b2*(b2 - c2)*sin2Alfa + a2*(a2 + 2*c2)*sin2Beta + 2*a2*b2*(-1 + cosAlfa*cosBeta*cosGamma);
	J = a2*(a2*sin2Beta - b2*sin2Alfa);

	L_ComplexDouble::solveCubicEquation(G, H, I, J, lambda0, lambda1, lambda2);

	double lambda;

	if (lambda0.isReal(0.01))
		lambda = lambda0.re;
	else if (lambda1.isReal(0.01))
		lambda = lambda1.re;
	else if (lambda2.isReal(0.01))
		lambda = lambda2.re;
	else
	{
		printf("No hay soluciones reales puras en la ecuacion cubica de Finsterwalder\n");
		printf("Opciones: %.3f +i*%.3f; %.3f+i*%.3f, %.3f+i*%.3f\n", lambda0.re, lambda0.im, lambda1.re, lambda1.im, lambda2.re, lambda2.im);
		return false;
	}

	double p, q;
	double A, B, C, D, E, F, eme[2], ene[2], AA[2], BB[2], CC[2], u[4], v[4];
	A = 1 + lambda;
	B = -cosAlfa;
	C = (b2-a2)/b2 - lambda*c2/b2;
	D = -lambda*cosGamma;
	E = (a2/b2 + lambda*c2/b2) * cosBeta;
	F = -a2/b2 + lambda*(b2-c2)/b2;

	p = sqrt(B*B-A*C);
	q = L_sgn(B*E-C*D) * sqrt(E*E - C*F); // v = - (B*u + E) +- (p*u + q)  =   u*m + n

	eme[0] = (-B + p)/C;
	ene[0] = ( -(E - q) )/C;
	eme[1] = (-B - p)/C;
	ene[1] = ( -(E + q) )/C;

	int z;
	for (z=0; z<2; z++) // z considera los dos signos posibles del termino variable de v que generan dos v, dos m y dos n
	{
		AA[z] = b*b - eme[z]*c*c;
		BB[z] = c2*(cosBeta - ene[z])*eme[z] - b2*cosGamma;
		CC[z] = -c*ene[z]*ene[z] + 2*c2*ene[z]*cosBeta + b2 - c2;

		// Eleccion del signo del termino variable de u
		if (BB[z]*BB[z]-AA[z]*CC[z] >= 0)
		{
			valida[2*z] = true;
			valida[2*z+1] = true;
			u[2*z] = -L_sgn(BB[z])/AA[z] * ( L_ABS(BB[z]) + sqrt(BB[z]*BB[z]-AA[z]*CC[z]) );
			u[2*z+1] = CC[z]/(AA[z]*u[2*z]);

			v[2*z] = u[2*z]*eme[z] + ene[z];
			v[2*z+1] = u[2*z+1]*eme[z] + ene[z];

			s1[2*z] = sqrt( a2 / (u[2*z]*u[2*z] + v[2*z]*v[2*z] - 2*u[2*z]*v[2*z]*cosAlfa) );
			s2[2*z] = u[2*z]*s1[2*z];
			s3[2*z] = v[2*z]*s1[2*z];

			s1[2*z+1] = sqrt( a2 / (u[2*z+1]*u[2*z+1] + v[2*z+1]*v[2*z+1] - 2*u[2*z+1]*v[2*z+1]*cosAlfa) );
			s2[2*z+1] = u[2*z+1]*s1[2*z+1];
			s3[2*z+1] = v[2*z+1]*s1[2*z+1];
		}
		else
		{
			valida[2*z] = false;
			valida[2*z+1] = false;
			s1[2*z]=s2[2*z]=s3[2*z]=s1[2*z+1]=s2[2*z+1]=s3[2*z+1]=sqrt(-1.0);
		}
	}

	return true;
}



bool L_CalculadorPoseTriangulo::pruebame()
{
	L_CoordsCart3D p[3];
	//double error = 0;
	int i, n;
	// double s1Real, s2Real, s3Real;

	for (n=0; n<4; n++)
	{
		for (i=0; i<3; i++)
		{
			//p[i].x = L_RANDOM(0.001, 1);
			//p[i].y = L_RANDOM(-1, 1) * p[i].x; // Para que quede en un cono
			//p[i].z = L_RANDOM(-1, 1) * p[i].x; // Para que quede en un cono
			p[i].x = 0.1 + (i+1)*(n+1)/15.0;
			p[i].y = (fmod((i+1)*(n+1)*0.3, 2) - 1) * p[i].x; // Para que quede en un cono
			p[i].z = (fmod(0.4 + (i+1)*(n+1)*0.5, 2) - 1) * p[i].x; // Para que quede en un cono
			printf("p[%d]: %.5f %.5f %.5f\n", i, p[i].x, p[i].y, p[i].z);
		}
		q1.tanIzq = p[1-1].y / p[1-1].x;
		q1.tanArr = p[1-1].z / p[1-1].x;
		q2.tanIzq = p[2-1].y / p[2-1].x;
		q2.tanArr = p[2-1].z / p[2-1].x;
		q3.tanIzq = p[3-1].y / p[3-1].x;
		q3.tanArr = p[3-1].z / p[3-1].x;

		a = p[2-1].distanciaA(p[3-1]);
		b = p[1-1].distanciaA(p[3-1]);
		c = p[1-1].distanciaA(p[2-1]);

		// para debug
		//s1Real = p[1-1].pideR();
		//s2Real = p[2-1].pideR();
		//s3Real = p[3-1].pideR();

		printf("v: %.5f %.5f %.5f = %.5f %.5f %.5f\n", p[0].x, p[0].y, p[0].z, q1.pideUnitario().x*p[0].pideR(), q1.pideUnitario().y*p[0].pideR(), q1.pideUnitario().z*p[0].pideR());
		printf("v: %.5f %.5f %.5f = %.5f %.5f %.5f\n", p[1].x, p[1].y, p[1].z, q2.pideUnitario().x*p[1].pideR(), q2.pideUnitario().y*p[1].pideR(), q2.pideUnitario().z*p[1].pideR());
		printf("v: %.5f %.5f %.5f = %.5f %.5f %.5f\n", p[2].x, p[2].y, p[2].z, q3.pideUnitario().x*p[2].pideR(), q3.pideUnitario().y*p[2].pideR(), q3.pideUnitario().z*p[2].pideR());

		//if ( metodoFinsterwalder() == false )
		//	return false;
		if ( metodoDirecto() == false )
			return false;

		printf("Iteracion %d:\n", n+1);
		printf("  Original: %.3f %.3f %.3f\n", p[0].pideR(), p[1].pideR(), p[2].pideR());
		printf("  Aprox  1: %.3f %.3f %.3f\n", s1[0], s2[0], s3[0]);
		printf("  Aprox  2: %.3f %.3f %.3f\n", s1[1], s2[1], s3[1]);
		printf("  Aprox  3: %.3f %.3f %.3f\n", s1[2], s2[2], s3[2]);
		printf("  Aprox  4: %.3f %.3f %.3f\n", s1[3], s2[3], s3[3]);

		//printf("  P1 Original: %.3f %.3f %.3f\n", p[0].x, p[0].y, p[0].z);
		//printf("  P1 Aprox  1: %.3f %.3f %.3f\n", (j1*s1[0]).x, (j1*s1[0]).y, (j1*s1[0]).z);
		//printf("  P1 Aprox  1: %.3f %.3f %.3f\n", (j1*s1[1]).x, (j1*s1[1]).y, (j1*s1[1]).z);
		//printf("  P1 Aprox  1: %.3f %.3f %.3f\n", (j1*s1[2]).x, (j1*s1[2]).y, (j1*s1[2]).z);
		//printf("  P1 Aprox  1: %.3f %.3f %.3f\n", (j1*s1[3]).x, (j1*s1[3]).y, (j1*s1[3]).z);
	}
	return true;
}

bool L_CalculadorPoseTrianguloGeneralizado::metodoNister()
{
	// David Nister, "A Minimal Solution to the Generalised 3-Point Pose Problem"
	L_HomogeneousMatrix H1, H2, H3, H4,  H21, H22, H22t;
	L_CoordsCart3D d1i, d2i, d3i, d4i, d5i, zero, t;
	L_CoordsCart3D p1i, p2i, p3i, p4i, R1p4i;
	double alfa, s;

	throw L_NonImplementedException(); // Le falta, dejo esto por ahora

	zero.fijaCero();
	d1i = r1.direcc;
	d2i = r2.direcc;
	d3i = r3.direcc;
	d4i = d1i.cruz(d2i);
	if (d4i.normaliza_ret() == false) // Podrian ser colineales
		return false;
	d5i = d1i.cruz(d4i);
	
	H1.fija4columnas(d5i, d1i, d4i, zero);

	p1i = r1.origen;
	p2i = r2.origen;
	p3i = r3.origen;
	s = d1i.punto(d2i) / d5i.punto(d2i);
	alfa = (d1i-s*d5i).punto(p2i-p1i);
	p4i = p1i + alfa*d1i;
	R1p4i = H1*p4i;
	H1(0,3) = -R1p4i.x;
	H1(1,3) = -R1p4i.y;
	H1(2,3) = -R1p4i.z;

	L_CoordsCart3D d1, d2, d3, d6, d7, d8, d9;
	L_CoordsCart3D p1, p2, p3, R2q1;
	double q2_q1;
	// Aplicar H1 a los rayos para dejarlos en una pose predeterminada
	p1 = H1*p1i;
	p2 = H1*p2i;
	p3 = H1*p3i;
	d1 = H1*(p1i+d1i) - p1;
	d2 = H1*(p2i+d2i) - p1;
	d3 = H1*(p1i+d1i) - p1;

	q2_q1 = q2.distanciaA(q1);
	d6.x = sqrt(q2_q1*q2_q1-p2.z*p2.z) / q2_q1;
	d6.y = 0;
	d6.z = p2.z / q2_q1;
	d7 = (q2-q1) / q2_q1;
	d8.define(0,1,0);
	d9 = d7.cruz(q3-q1);
	if (d9.normaliza_ret() == false)
		return false;
	
	t = d6.cruz(d8);
	H21.fija4columnas(d6, d8, t, zero);
	t = d7.cruz(d9);
	H22.fija4columnas(d7, d9, t, zero);
	L_StaticMatrix_transpDe(H22t,H22);
	H2.OP_mult(H21, H22t);

	R2q1 = H2*q1;
	H2(0,3) = -R2q1.x;
	H2(1,3) = -R2q1.y;
	H2(2,3) = -R2q1.z;

	q1 = H2*q1;
	q2 = H2*q2;
	q3 = H2*q3;

	L_CoordsCart3D uoo(1,0,0), ouo(0,1,0), ncand1, ncand2, n, nP;
	
	ncand1 = d3.cruz(uoo);
	ncand2 = d3.cruz(ouo);

	if (ncand1.pideR2() > ncand2.pideR2())
		n = ncand1;
	else
		n = ncand2;

	nP = d3.cruz(n);
	
	double L[5], LP[5]; // No voy a usar L[0] ni LP[0], en este caso voy a partir del 1 para copiar la notacion del paper
	L[1] = n.x;
	L[2] = n.y;
	L[3] = n.z;
	L[4] = -n.punto(p3);
	LP[1] = nP.x;
	LP[2] = nP.y;
	LP[3] = nP.z;
	LP[4] = -nP.punto(p3);

	throw L_NonImplementedException(); // Le falta, dejo esto por ahora

	return false;
}

void L_Camara::calcRayos(L_Matrix &uv_nx2) const
{
	L_RayoCamara r;
	int i;
	for (i=0; i<uv_nx2.li; i++)
	{
		calcRayo(r, uv_nx2(i,0), uv_nx2(i,1));
		uv_nx2(i,0)=r.tanIzq;
		uv_nx2(i,1)=r.tanArr;
	}
}

void L_Camara::calcPixeles(L_Matrix &uv_nx2) const
{
	L_RayoCamara r;
	int i;
	for (i=0; i<uv_nx2.li; i++)
	{
		r.tanIzq=uv_nx2(i,0);
		r.tanArr=uv_nx2(i,1);
		calcPixel(r, uv_nx2(i,0), uv_nx2(i,1));
	}
}

void L_CamaraPinhole::calcRayos(L_Matrix &uv_nx2) const
{
	L_RayoCamara r;
	int i;
	for (i=0; i<uv_nx2.li; i++)
	{
		r.tanIzq=(cenX-uv_nx2(i,0))/dfX;
		r.tanArr=(cenY-uv_nx2(i,1))/dfY;
		uv_nx2(i,0)=r.tanIzq;
		uv_nx2(i,1)=r.tanArr;
	}
}

void L_CamaraPinhole::calcPixeles(L_Matrix &uv_nx2) const
{
	L_RayoCamara r;
	int i;
	for (i=0; i<uv_nx2.li; i++)
	{
		r.tanIzq=uv_nx2(i,0);
		r.tanArr=uv_nx2(i,1);
		uv_nx2(i,0) = cenX-r.tanIzq*dfX;
		uv_nx2(i,1) = cenY-r.tanArr*dfY;
	}
}

void L_CamaraPinhole::fijaMatrizCalibracion(const L_Matrix &K3x3)
{
	// xPixDer*w = K3x3[0,0)*tanDer + K3x3[0,1)*tanAba + K3x3[0,2)*1
	// yPixAba*w = K3x3[1,0)*tanDer + K3x3[1,1)*tanAba + K3x3[1,2)*1
	// w    =      K3x3[2,0)*tanDer + K3x3[2,1)*tanAba + K3x3[2,2)*1

	// xPixDer*w = -K3x3[0,0)*tanIzq                      + K3x3[0,2)*1
	// yPixAba*w =                      -K3x3[1,1)*tanArr + K3x3[1,2)*1
	// w    =                                                K3x3[2,2)

	// xPixDer = K3x3[0,2)/K3x3[2,2) - K3x3[0,0)/K3x3[2,2)*tanIzq
	// xPixArr = K3x3[1,2)/K3x3[2,2) - K3x3[1,1)/K3x3[2,2)*tanArr

	// xPixDer = cenX-rayo.tanIzq*dfX
	// yPixAba = cenY-rayo.tanArr*dfY;

	cenX = K3x3(0,2) / K3x3(2,2);
	cenY = K3x3(1,2) / K3x3(2,2);

	dfX = K3x3(0,0) / K3x3(2,2);
	dfY = K3x3(1,1) / K3x3(2,2);

	// Los parametros que siguen son opcionales, se obtienen aproximados
	resX = (int)(2*cenX+1); // Una aproximacion
	resY = (int)(2*cenY+1); // Una aproximacion
	// dfX = (resX/2)/tan(fovX/2)
	// tan(fovX/2) = (resX/2) / (dfX)
	fovX = 2*atan2(resX/2, dfX);
	fovY = 2*atan2(resY/2, dfY);
}

void L_CamaraPinhole::pedirMatrizCalibracion(L_Matrix &K3x3) const
{
	K3x3.reallocate(3,3);
	K3x3(0,0) = dfX;
	K3x3(0,1) = 0;
	K3x3(0,2) = cenX;
	//
	K3x3(1,0) = 0;
	K3x3(1,1) = dfY;
	K3x3(1,2) = cenY;
	//
	K3x3(2,0) = 0;
	K3x3(2,1) = 0;
	K3x3(2,2) = 1;
}

void L_CamaraDistorsionRadial::calcRayo(L_RayoCamara &rayo, double xPix, double yPix) const
{
	if (tipoModelo == 0)
	{
		rayo.tanIzq = (cenX-xPix)/dfX;
		rayo.tanArr = (cenY-yPix)/dfY;
		double r2 = rayo.tanIzq*rayo.tanIzq + rayo.tanArr*rayo.tanArr;
		double factRad = 1 + k1*r2 + k2*r2*r2;
		rayo.tanIzq*=factRad;
		rayo.tanArr*=factRad;
	}
	else
	{
		throw_L_ArgException_if(true, "L_CamaraDistorsionRadial::calcRayo() : intento de calcular rayo en modelo de Heikkil");
	}
}

void L_CamaraDistorsionRadial::calcPixel(const L_RayoCamara &rayo, double &xPix, double &yPix) const
{
	if (tipoModelo == 0)
	{
		double xPixPinh, yPixPinh;
		xPixPinh = -dfX*rayo.tanIzq + cenX; // La pendiente del rayo cambia al pasar por el lente en "factRad%"
		yPixPinh = -dfY*rayo.tanArr + cenY; // La pendiente del rayo cambia al pasar por el lente en "factRad%"
		distorsionaPixel_k1_Tsai(xPixPinh, yPixPinh, xPix, yPix);
	}
	else
		distorsionaPixel_k1k2k3k4_Heikkil(rayo.tanIzq, rayo.tanArr, xPix, yPix);
}

void L_CamaraDistorsionRadial::calcRayos(L_Matrix &uv_nx2) const
{
	double r2, factRad;
	int i;
	if (tipoModelo == 0)
	{
		for (i=0; i<uv_nx2.li; i++)
		{
			uv_nx2(i,0)=(cenX-uv_nx2(i,0))/dfX;
			uv_nx2(i,1)=(cenY-uv_nx2(i,1))/dfY;
			r2 = uv_nx2(i,0)*uv_nx2(i,0) + uv_nx2(i,1)*uv_nx2(i,1);
			factRad = 1 + k1*r2 + k2*r2*r2;
			uv_nx2(i,0)*=factRad;
			uv_nx2(i,1)*=factRad;
		}
	}
	else
	{
		throw_L_ArgException_if(true, "L_CamaraDistorsionRadial::calcRayos() : intento de calcular rayo en modelo de Heikkil");
	}
}

void L_CamaraDistorsionRadial::calcPixeles(L_Matrix &uv_nx2) const
{
	double x, y;
	int i;
	if (tipoModelo == 0)
	{
		for (i=0; i<uv_nx2.li; i++)
		{
			x = -dfX*uv_nx2(i,0) + cenX;
			y = -dfY*uv_nx2(i,1) + cenY;
			distorsionaPixel_k1_Tsai(x, y, uv_nx2(i,0), uv_nx2(i,1));
		}
	}
	else
	{
		for (i=0; i<uv_nx2.li; i++)
		{
			x = -dfX*uv_nx2(i,0) + cenX;
			y = -dfY*uv_nx2(i,1) + cenY;
			distorsionaPixel_k1k2k3k4_Heikkil(x, y, uv_nx2(i,0), uv_nx2(i,1));
		}
	}
}

void L_CamaraDistorsionRadial::corrigeImagenPixelado(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr)
{
	int i, j, u, v;
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	corr.reallocate(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			distorsionaPixel_LUT(i, j, u, v);
			if (u>=0 && v>=0 && u<distor.lx && v<distor.ly)
			{
				corr.pix(i,j,0) = distor.pix(u,v,0);
				corr.pix(i,j,1) = distor.pix(u,v,1);
				corr.pix(i,j,2) = distor.pix(u,v,2);
			}
			else
			{
				corr.pix(i,j,0) = 0;
				corr.pix(i,j,1) = 0;
				corr.pix(i,j,2) = 0;
			}
		}
	}
}

void L_CamaraDistorsionRadial::corrigeImagenBilineal(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr)
{
	int i, j, u, v;
	L_RayoCamara rayo;
	double alfa, beta;
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	corr.reallocate(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			distorsionaPixel_LUT(i, j, alfa, beta);
			u = (int)alfa;
			v = (int)beta;
			alfa = alfa-u; // alfa > u -> alfa entre 0 y 1
			beta = beta-v; // beta > v -> beta entre 0 y 1
			if (u>=0 && v>=0 && u<distor.lx-1 && v<distor.ly-1)
			{
				// Mmm... voy alfa tratar de dejarlo lo mas claro "posible"
#if defined(a) || defined(b) || defined(c) || defined(b)
#error A quien se le ocurre define a, b, c o d...
#endif

#define			a   (1-alfa)*distor.pix(u,v,0)+alfa*distor.pix(u+1,v,0)
#define			b   (1-alfa)*distor.pix(u,v+1,0)+alfa*distor.pix(u+1,v+1,0)
				corr.pix(i,j,0) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b

#define			a   (1-alfa)*distor.pix(u,v,1)+alfa*distor.pix(u+1,v,1)
#define			b   (1-alfa)*distor.pix(u,v+1,1)+alfa*distor.pix(u+1,v+1,1)
				corr.pix(i,j,1) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b

#define			a   (1-alfa)*distor.pix(u,v,2)+alfa*distor.pix(u+1,v,2)
#define			b   (1-alfa)*distor.pix(u,v+1,2)+alfa*distor.pix(u+1,v+1,2)
				corr.pix(i,j,2) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b
			}
			else
			{
				corr.pix(i,j,0) = 0;
				corr.pix(i,j,1) = 0;
				corr.pix(i,j,2) = 0;
			}
		}
	}
}

void L_CamaraDistorsionRadial::corrigeImagenBilineal_k1(const L_ImageRGBUchar &distor, L_ImageRGBUchar &corr)
{
	int i, j, u, v;
	L_RayoCamara rayo;
	double alfa, beta;
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	corr.reallocate(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			distorsionaPixel_LUT(i, j, alfa, beta);
			u = (int)alfa;
			v = (int)beta;
			alfa = alfa-u; // alfa > u -> alfa entre 0 y 1
			beta = beta-v; // beta > v -> beta entre 0 y 1
			if (u>=0 && v>=0 && u<distor.lx-1 && v<distor.ly-1)
			{
				// Mmm... voy alfa tratar de dejarlo lo mas claro "posible"
#if defined(a) || defined(b) || defined(c) || defined(b)
#error A quien se le ocurre define a, b, c o d...
#endif

#define			a   (1-alfa)*distor.pix(u,v,0)+alfa*distor.pix(u+1,v,0)
#define			b   (1-alfa)*distor.pix(u,v+1,0)+alfa*distor.pix(u+1,v+1,0)
				corr.pix(i,j,0) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b

#define			a   (1-alfa)*distor.pix(u,v,1)+alfa*distor.pix(u+1,v,1)
#define			b   (1-alfa)*distor.pix(u,v+1,1)+alfa*distor.pix(u+1,v+1,1)
				corr.pix(i,j,1) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b

#define			a   (1-alfa)*distor.pix(u,v,2)+alfa*distor.pix(u+1,v,2)
#define			b   (1-alfa)*distor.pix(u,v+1,2)+alfa*distor.pix(u+1,v+1,2)
				corr.pix(i,j,2) = (L_uchar)( (1-beta)* (a) + beta * (b) );
#undef a
#undef b
			}
			else
			{
				corr.pix(i,j,0) = 0;
				corr.pix(i,j,1) = 0;
				corr.pix(i,j,2) = 0;
			}
		}
	}
}

void L_CamaraDistorsionRadial::corrigeImagenPixelado(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr)
{
	int i, j, u, v;
	L_RayoCamara rayo;
	//double r2, factRad;
	corr.reallocate(distor.lx, distor.ly);
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			distorsionaPixel_LUT(i, j, u, v);
			if (u>=0 && v>=0 && u<distor.lx && v<distor.ly)
				corr.pix(i,j) = distor.pix(u,v);
			else
				corr.pix(i,j) = 0;
		}
	}
}

void L_CamaraDistorsionRadial::corrigeImagenBilineal(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr)
{
	int i, j, u, v;
	L_RayoCamara rayo;
//	double r2, factRad;
	double alfa, beta;
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	corr.reallocate(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			distorsionaPixel_LUT(i, j, alfa, beta);
			u = (int)alfa;
			v = (int)beta;
			alfa = alfa-u; // alfa > u -> alfa entre 0 y 1
			beta = beta-v; // beta > v -> beta entre 0 y 1
			if (u>=0 && v>=0 && u<distor.lx-1 && v<distor.ly-1)
			{
				// Mmm... voy alfa tratar de dejarlo lo mas claro "posible"
#if defined(a) || defined(b) || defined(c) || defined(b)
#error A quien se le ocurre define a, b, c o d...
#endif

#define			a   (1-alfa)*distor.pix(u,v)+alfa*distor.pix(u+1,v)
#define			b   (1-alfa)*distor.pix(u,v+1)+alfa*distor.pix(u+1,v+1)
				corr.pix(i,j) = (1-beta)* (a) + beta * (b);
#undef a
#undef b
			}
			else
			{
				corr.pix(i,j) = 0;
			}
		}
	}
}


void L_CamaraDistorsionRadial::corrigeImagenBilineal_k1(const L_ImageGrayDouble &distor, L_ImageGrayDouble &corr)
{
	int i, j, u, v;
	L_RayoCamara rayo;
	double alfa, beta;
	if (tipoModelo == 0)
		calcLUT_k1_Tsai(distor.lx, distor.ly);
	else
		calcLUT_k1k2k3k4_Heikkil(distor.lx, distor.ly);
	corr.reallocate(distor.lx, distor.ly);
	for (j=0; j<corr.ly; j++)
	{
		for (i=0; i<corr.lx; i++)
		{
			// Transformar (i,j) sin deformar en rayo
			distorsionaPixel_LUT(i, j, alfa, beta);
			u = (int)alfa;
			v = (int)beta;
			alfa = alfa-u; // alfa > u -> alfa entre 0 y 1
			beta = beta-v; // beta > v -> beta entre 0 y 1
			if (u>=0 && v>=0 && u<distor.lx-1 && v<distor.ly-1)
			{
				// Mmm... voy alfa tratar de dejarlo lo mas claro "posible"
#if defined(a) || defined(b) || defined(c) || defined(b)
#error A quien se le ocurre define a, b, c o d...
#endif

#define			a   (1-alfa)*distor.pix(u,v)+alfa*distor.pix(u+1,v)
#define			b   (1-alfa)*distor.pix(u,v+1)+alfa*distor.pix(u+1,v+1)
				corr.pix(i,j) = (1-beta)* (a) + beta * (b);
#undef a
#undef b
			}
			else
			{
				corr.pix(i,j) = 0;
			}
		}
	}
}

void L_CamaraDistorsionRadial::muestraImagenError(const L_ImageRGBUchar &orig, L_ImageRGBUchar &final, int niter)
{
	int i, j, k, u, v;
	double x, y;
	L_RayoCamara r;
	final.reallocate(orig.lx, orig.ly);
	for (j=0; j<final.ly; j++)
	{
		for (i=0; i<final.lx; i++)
		{
			x = i;
			y = j;
			for (k=0; k<niter; k++)
			{
				calcRayo(r, x, y);
				calcPixel(r, x, y);
			}
			u = (int) (x + 0.5);
			v = (int) (y + 0.5);
			if (u>=0 && v>=0 && u<orig.lx && v<orig.ly)
			{
				final.pix(i,j,0)=orig.pix(u,v,0);
				final.pix(i,j,1)=orig.pix(u,v,1);
				final.pix(i,j,2)=orig.pix(u,v,2);
			}
			else
			{
				final.pix(i,j,0)=0;
				final.pix(i,j,1)=0;
				final.pix(i,j,2)=0;
			}
		}
	}
}

void L_CamaraDistorsionRadial::muestraImagenError(const L_ImageGrayDouble &orig, L_ImageGrayDouble &final, int niter)
{
	int i, j, k, u, v;
	double x, y;
	L_RayoCamara r;
	final.reallocate(orig.lx, orig.ly);
	for (j=0; j<final.ly; j++)
	{
		for (i=0; i<final.lx; i++)
		{
			x = i;
			y = j;
			for (k=0; k<niter; k++)
			{
				calcRayo(r, x, y);
				calcPixel(r, x, y);
			}
			u = (int) (x+0.5);
			v = (int) (y+0.5);
			if (u>=0 && v>=0 && u<orig.lx && v<orig.ly)
				final.pix(i,j)=orig.pix(u,v);
			else
				final.pix(i,j)=0;
		}
	}
}

void L_CamaraDistorsionRadial::distorsionaPixel(double xPix, double yPix, double &xPixDist, double &yPixDist)
{
	if (tipoModelo == 0)
		distorsionaPixel_k1_Tsai(xPix, yPix, xPixDist, yPixDist);
	else
		distorsionaPixel_k1k2k3k4_Heikkil(xPix, yPix, xPixDist, yPixDist);
}

void L_CamaraDistorsionRadial::guardarParametros(FILE *fp)
{
	fprintf(fp,
		"resX=%d\n"
		"resY=%d\n"
		"cenX=%lf\n"
		"cenY=%lf\n"
		"k1=%lf\n"
		"k2=%lf\n"
		"dfX=%lf\n"
		"dfY=%lf\n", resX, resY, cenX, cenY, k1, k2, dfX, dfY);
}
bool L_CamaraDistorsionRadial::leerParametros(FILE *fp)
{
	char str[100];
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "resX=%d", &resX) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "resY=%d", &resY) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "cenX=%lf", &cenX) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "cenY=%lf", &cenY) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "k1=%lf", &k1) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "k2=%lf", &k2) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "dfX=%lf", &dfX) < 1)
		return false;
	if (fgets(str, 99, fp) == NULL)
		return false;
	if (sscanf(str, "dfY=%lf", &dfY) < 1)
		return false;
	return true;
}


bool L_CuadroCalibradorCamara::guardar(FILE *fp)
{
	int i;
	fprintf(fp, "npuntos=%d\n", v.size());
	for (i=0; i<(int)v.size(); i++)
		fprintf(fp, "%lf %lf %lf %lf %lf\n", v[i].derMM, v[i].abaMM, v[i].adeMM, v[i].xPix, v[i].yPix);
	return true;
}

bool L_CuadroCalibradorCamara::read(FILE *fp)
{
	int i, n;
	char str[200];
	if (fgets(str, 200, fp) == NULL)
		return false;
	if ( sscanf(str, "npuntos=%d\n", &n) < 1)
		return false;
	v.resize(n);
	for (i=0; i<(int)v.size(); i++)
	{
		if (fgets(str, 200, fp) == NULL)
			return false;
		if (sscanf(str, "%lf%lf%lf%lf%lf", &v[i].derMM, &v[i].abaMM, &v[i].adeMM, &v[i].xPix, &v[i].yPix) < 5)
			return false;
	}
	return true;
}

bool L_CuadroCalibradorCamara::guardar_TsaiCM(const char *nomarch)
{
	FILE *fp;
	int i;
	fp = fopen(nomarch, "w");
	if (fp == NULL)
		return false;
	for (i=0; i<(int)v.size(); i++)
		fprintf(fp, "%lf %lf %lf %lf %lf\n", v[i].derMM, v[i].abaMM, v[i].adeMM, v[i].xPix, v[i].yPix);
	fclose(fp);
	return true;
}

bool L_CuadroCalibradorCamaraArr::guardar(FILE *fp)
{
	int i;
	fprintf(fp, "ncuadros=%ld\n", long(size()));
	for (i=0; i<(int)v.size(); i++)
		if ( v[i].guardar(fp) == false)
			return false;
	return true;
}

bool L_CuadroCalibradorCamaraArr::read(FILE *fp)
{
	int i;
	int nExtra;
	char str[200];
	if (fgets(str, 199, fp) == NULL)
		return false;
	if ( sscanf(str, "ncuadros=%d", &nExtra) < 1)
		return false;
	v.resize(v.size()+nExtra);
	for (i=(int)v.size()-nExtra; i<(int)v.size(); i++)
	{
		if ( v[i].read(fp) == false)
			return false;
		nPuntosEstad.push(v[i].v.size());
	}
	return true;
}

void L_TransfCov_tresMatrices::init(const L_Matrix &Pppe)
{
	int i, u, v;
	Ppp = Pppe;
	Ptridiag.reallocate(Ppp.li, Ppp.lj);
	Ptridiag.setZero();
	for (i=0; i<Ppp.li/3; i++)
	{
		for (u=0; u<3; u++)
			for (v=0; v<3; v++)
				Ptridiag(3*i+u,3*i+v) = Ppp(3*i+u,3*i+v);
	}
	Pdiag.reallocate(Ppp.li, Ppp.lj);
	Pdiag.setZero();
	for (i=0; i<Ppp.li; i++)
		Pdiag(i,i) = Ppp(i,i); // OJO: faltaba una [i) !!!
}

double L_TransfCov_tresMatrices::error(double lambda1, double lambda2) const
{
	L_Matrix PDIF, Ptmp1, Ptmp2, valPr;
	PDIF=Ppp;
	Ptmp1.OP_mult(Ptridiag,lambda1);
	Ptmp2.OP_mult(Pdiag,lambda2);
	PDIF-=Ptmp1;
	PDIF-=Ptmp2;
	Ptmp1 = PDIF;
	Ptmp1.svd_ordenado(valPr, Ptmp2);
	return valPr(valPr.li-1,0); // Para que de positivo
}

double L_TransfCov_tresMatrices::errorGen(const void *object, double *vect)
{
	const L_TransfCov_tresMatrices *obj;
	double err;
	obj = (const L_TransfCov_tresMatrices *)object;
	err = obj->error(vect[0], vect[1]);
	return (err-1.0e-30)*(err-1.0e-30);
}

void L_TransfCov_tresMatrices::errorGen_vect(const void *object, const L_Matrix &vect, L_Matrix &err)
{
	const L_TransfCov_tresMatrices *obj;
	obj = (const L_TransfCov_tresMatrices *)object;
	err.reallocate(1,1);
	err(0,0) = obj->error(vect(0,0), vect(1,0));
}

void L_TransfCov_tresMatrices::encontrarCovarianzasIndividuales(const L_Matrix &Ppp, L_Matrix &P_arr_3x3)
{
	L_LevenbergMarquardt opt;
	std::vector<double> lambdas(2);
	init(Ppp);
	opt.epsilon = 1.0e-6;
	opt.lenVector = 2;
	opt.nIterationsMax = 10;
	lambdas[0] = 0;
	lambdas[1] = 0;
	opt.minimize((void *)this, lambdas, &errorGen);
	L_Matrix Ptmp1, Ptmp2;
	Ptmp1.OP_mult(Ptridiag, lambdas[0]);
	Ptmp2.OP_mult(Pdiag, lambdas[1]);
	P_arr_3x3.OP_add(Ptmp1, Ptmp2);
}

void L_TransfCov_tresMatrices::encontrarCovarianzasIndividuales_vect(const L_Matrix &Ppp, L_Matrix &P_arr_3x3)
{
	L_LevenbergMarquardt opt;
	L_Matrix lambdas(2,1);
	init(Ppp);
	opt.epsilon = 1.0e-6;
	opt.lenVector = 2;
	opt.nIterationsMax = 10;
	lambdas(0,0) = 0;
	lambdas(1,0) = 0;
	opt.minimize_vect((void *)this, lambdas, &errorGen_vect);
	L_Matrix Ptmp1, Ptmp2;
	Ptmp1.OP_mult(Ptridiag, lambdas(0,0));
	Ptmp2.OP_mult(Pdiag, lambdas(1,0));
	P_arr_3x3.OP_add(Ptmp1, Ptmp2);
}

void L_TransfCov_tresMatrices::pruebame()
{
	L_TransfCov_tresMatrices tres;
	L_Matrix Pe(10,10), PeT, M, I, MD, D2, M2, err;
	double fac;
	srand(40);
	I = Pe;
	I.identity();
	I(I.li-1,I.li-1) = 0;
	Pe.fillRandomValues(-10, 10);
	PeT.transpOf(Pe);
	err.OP_mult(Pe, I);
	M.OP_mult(err, PeT);
	MD=M;
	fac = L_RANDOM(1, 3);
	for (int i=0; i<MD.li; i++)
			MD(i,i) *= fac;
	M.print("M");
	MD.print("MD");
	tres.encontrarCovarianzasIndividuales_vect(MD, D2);
	D2.print("D2");
	M2.OP_subtract(MD, D2);
	M2.print("M2");
	err.OP_subtract(M, M2);
	err.print("err");
	if (M2.isGoodCovarianceMatrix() == true)
		printf("La matriz M2 es valida\n");
	getchar();
}

#ifdef __COMPAT_IPLIMAGE__
bool L_CapturadorImagenCV::crear()
{
	activo=true;
	if (name.size() != 0)
	{
		capture = cvCaptureFromAVI(name.c_str());
		dt = 1/cvGetCaptureProperty(capture, CV_CAP_PROP_FPS);
	}
	else
	{
		capture = cvCaptureFromCAM(-1);
		dt = 0;
	}
	activo = (capture != NULL);
	return activo;
}
#endif // __COMPAT_SVS__


#ifdef __COMPAT_IPLIMAGE__
bool L_CapturadorImagenCV::capturarImagen()
{
	if (!activo && crear()== false)
		return false;
	imCam = cvQueryFrame( capture );
	if (imCam==NULL)
		return false;

	estereo.esEstereo = false;
	im=*imCam;
	// Si la imagen es mas ancha que 2 veces su height, se supone estereo
	if (esArchivo() == true && imCam->width > 2*imCam->height)
	{
		estereo.esEstereo = true;
		L_ImageRGBUchar imC;
		im.swap(imC);
		imC.divideHorz(im, imDer);
		if (estereo.usarIzq0Der1Ambas2 == 1)
			imDer.swap(im);
	}

	nFrame++;
	return true;
}
#endif // __COMPAT_SVS__





#ifdef __COMPAT_SVS__
#define L_USA_SETJMP_VIDERE
#endif


#ifdef L_USA_SETJMP_VIDERE
jmp_buf buf_eseveese;   // Estado del stack, por alguna razon esto no funciona en algunos lados aunque es ANSI C
#endif

void L_Catchhh_eseveese(int i)
{
#ifdef L_USA_SETJMP_VIDERE
	// Intento de que el programa supere el error del driver Videre
	longjmp(buf_eseveese, 1); // Se salta al estado del stack en buf_eseveese
#else
	printf("* Error en el driver de la camara Videre, segmentation fault en svs::getVideoObject()\n");
	printf("* La camara Videre ya estaba transmitiendo\n");
	printf("* Presione ENTER para terminar el programa\n");
	getchar();
	exit(1);
#endif

	return;
}

#ifdef __COMPAT_SVS__
bool L_CapturadorImagenSVS::crearRec(int n)
{
	int i = 0, valjmp=0;
	bool bul;
	if (n<0)
		return false;

#ifdef L_USA_SETJMP_VIDERE
	valjmp = setjmp(buf_eseveese); // Se guarda el estado del stack en buf_eseveese
#endif
	if (valjmp == 0)  // "TRY"
	{
		closeVideoObject(); // Esta wea no sirve

		signal(SIGSEGV, L_Catchhh_eseveese); // Funcion de manejo de error que se llama si el driver Videre se cae
		videoObject = getVideoObject(); // Esta wa se cae cuando la camara ya esta trabajando
		signal(SIGSEGV, SIG_DFL); // Eliminar la funcion de manejo de error, si llegamos aca no se cayo el driver
	}
	else   // "CATCH SIGSEGV"
	{
		signal(SIGSEGV, SIG_DFL); // Se cayo el driver Videre, eliminar la funcion de manejo de error
		activo = false;
		return crearRec(n-1);
	}
	
	if (videoObject->Open() == false)
	{
		activo = false;
		bul = crearRec(n-1);
		if (bul == false)
			closeVideoObject();
		return bul;
	}

	videoObject->SetBrightness(true, 1000);

	videoObject->SetSize(lx, ly);
	if (videoObject->Start() == false)
	{
		videoObject->Close(); // Puede que la vez anterior haya quedado abierta
		return crearRec(n-1);
	}
	videoObject->SetColor(true, true);
	activo=true;
	return true;
}
#endif // COMPAT_SVS


#ifdef __COMPAT_SVS__
bool L_CapturadorImagenSVS::capturarImagen()
{
	if (!activo && crear() == false)
		return false;
	imageObject = videoObject->GetImage(400);
	
	if (estereo.usarIzq0Der1Ambas2 == 0)
	{
		return im.SVSPideIzq(imageObject);
	}
	else if (estereo.usarIzq0Der1Ambas2 == 1)
	{
		return im.SVSPideDer(imageObject);
	}
	return im.SVSPideIzq(imageObject) && imDer.SVSPideDer(imageObject);
}
#endif // COMPAT_SVS

bool L_CapturadorImagenBMZ::capturarImagen()
{
	long ncamaras = 0;
	if (!activo || fp==NULL)
		return false;
	if (im.readComprRGB(fp, &t, &ncamaras) == -1)
		return false;
	estereo.esEstereo = false;
	if (ncamaras == 2) // Estereo supuestamente
	{
		estereo.esEstereo = true;
		L_ImageRGBUchar imC;
		t = -t;
		im.swap(imC);
		imC.divideHorz(im, imDer);
		if (estereo.usarIzq0Der1Ambas2 == 1)
			imDer.swap(im);
	}
	return true;
}

bool L_CapturadorImagenesCarpeta::capturarImagen()
{
	if (i > name.numMax)
		return false;
	if (im.readImage(name.generateString(i)) == 0)
		return false;
	if (nombreDer.size() != 0 && imDer.readImage(nombreDer.generateString(i)) == 0)
		return false;
	if (nombreKin.size() != 0 && imKin.readImage_16_txt(nombreKin.generateString(i)) == 0)
		return false;
	if (nombreProf.size() != 0 && imDepth.readImage(nombreProf.generateString(i)) == 0)
		return false;
	i++;
	return true;
}



#ifdef __COMPAT_OPENNI__
bool L_CapturadorImagenKinectNI::crear()
{
	XnStatus nRetVal = XN_STATUS_OK;
	activo = false;

	nRetVal = context.Init();
	if (nRetVal != XN_STATUS_OK)
	{
		//printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}

	nRetVal = depthGen.Create(context);
	if (nRetVal != XN_STATUS_OK)
	{
		//printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}
	
	nRetVal = imageGen.Create(context);
	if (nRetVal != XN_STATUS_OK)
	{
		//printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}

	XnMapOutputMode mapMode;

	mapMode.nXRes = XN_VGA_X_RES;
	mapMode.nYRes = XN_VGA_Y_RES;
	mapMode.nFPS = 30;

	nRetVal = depthGen.SetMapOutputMode(mapMode);
	if (nRetVal != XN_STATUS_OK)
	{
		printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}

	nRetVal = imageGen.SetMapOutputMode(mapMode);
	if (nRetVal != XN_STATUS_OK)
	{
		printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}

	nRetVal = context.StartGeneratingAll();
	if (nRetVal != XN_STATUS_OK)
	{
		printf("Error en Open NI: %s\n", xnGetStatusString(nRetVal));
		//context.Shutdown();
		context.Release();
		return false;
	}
	
	activo = true;
	return true;
}
#endif // __COMPAT_OPENNI__


#ifdef __COMPAT_OPENNI__
bool L_CapturadorImagenKinectNI::cerrar()
{
	//context.Shutdown();
	context.Release();
	activo = false;
	return true;
}
#endif // __COMPAT_OPENNI__


#ifdef __COMPAT_OPENNI__
void L_CapturadorImagenKinectNI::test_me()
{
	L_CapturadorImagenKinectNI capt;
	L_ImageRGBUchar imp;
	L_VentanaImagen vent1(640,480,"RGB"), vent2(640,480,"Depth");

	while(true)
	{
		capt.capturarImagen();
		vent1.dibuja(capt.im);
		capt.imDepth.multiply_to_each_element(1.0/4096);
		imp = capt.imDepth;
		vent2.dibuja(imp);
		vent1.check();
	}
}
#endif // __COMPAT_OPENNI__


#ifdef __COMPAT_KINECTSDK__
bool L_CapturadorImagenKinectSDK::crear()
{
	HRESULT ret;
	
	ret = NuiInitialize( NUI_INITIALIZE_FLAG_USES_COLOR | NUI_INITIALIZE_FLAG_USES_DEPTH );
	if ( ret != S_OK )
	{
		return false;
	}

	NUI_IMAGE_RESOLUTION resolution = NUI_IMAGE_RESOLUTION_640x480;

	imageEvent = ::CreateEvent( 0, TRUE, FALSE, 0 );
	imageHandle = 0;
	ret = NuiImageStreamOpen( NUI_IMAGE_TYPE_COLOR, resolution, 0, 2, imageEvent, &imageHandle );
	if ( ret != S_OK )
	{
		return false;
	}

	depthEvent = ::CreateEvent( 0, TRUE, FALSE, 0 );
	depthHandle = 0;
	ret = NuiImageStreamOpen( NUI_IMAGE_TYPE_DEPTH, resolution, 0, 2, depthEvent, &depthHandle );
	if ( ret != S_OK )
	{
		return false;
	}
	return true;
}
#endif // __COMPAT_KINECTSDK__


#ifdef __COMPAT_KINECTSDK__
bool L_CapturadorImagenKinectSDK::capturarImagen()
{
	long i;
	DWORD width, height;
	HRESULT ret;

	WaitForSingleObject( imageEvent, INFINITE );

	CONST NUI_IMAGE_FRAME *imageFrame = 0;
	CONST NUI_IMAGE_FRAME *depthFrame = 0;
	ret = NuiImageStreamGetNextFrame( imageHandle, 0, &imageFrame );
	if (ret != S_OK)
		return false;

	KINECT_LOCKED_RECT rect;
	ret = imageFrame->pFrameTexture->LockRect( 0, &rect, 0, 0 );
	if (ret != S_OK)
		return false;

	NuiImageResolutionToSize(imageFrame->eResolution, width, height);
	im.reallocate(width, height);
	for (i=0; i<(long)(width*height); i++)
	{
		im.data()[i][0] = ((unsigned char *)rect.pBits)[4*i+0];
		im.data()[i][1] = ((unsigned char *)rect.pBits)[4*i+1];
		im.data()[i][2] = ((unsigned char *)rect.pBits)[4*i+2];
	}

	ret = NuiImageStreamReleaseFrame( imageHandle, imageFrame );
	if (ret != S_OK)
		return false;

	while(redimensionar.im_activo && (im.lx > redimensionar.im_lx || im.ly > redimensionar.im_ly))
		im.halfResolution();

	if (prof.capturarProf)
	{
		WaitForSingleObject( depthEvent, INFINITE );
		depthFrame = 0;
		ret = NuiImageStreamGetNextFrame( depthHandle, 0, &depthFrame );
		if (ret != S_OK)
			return false;
		ret = depthFrame->pFrameTexture->LockRect( 0, &rect, 0, 0 );
		if (ret != S_OK)
			return false;

		NuiImageResolutionToSize(depthFrame->eResolution, width, height);
		imDepth.reallocate(width, height);
		imKin.reallocate(width, height);
		
		//memcpy(imKin.data(), rect.pBits, width*height*sizeof(unsigned short));
		for (i=0; i<(long)(width*height); i++)
		{
			imKin.data()[i] = ((unsigned short *)rect.pBits)[i+0];
		}
		ret = NuiImageStreamReleaseFrame( depthHandle, depthFrame );
		if (ret != S_OK)
			return false;
		// Resize image if needed
		while(redimensionar.imProf_activo && (imKin.lx > redimensionar.imProf_lx || imKin.ly > redimensionar.imProf_ly))
			imKin.halfResolution_noResample();

		if (prof.calcularProf)
		{
			imDepth.reallocate(imKin.lx, imKin.ly);
			for (int i=0; i<imKin.lx*imKin.ly; i++)
				imDepth.data()[i] = imKin.data()[i]; // It is yet in [mm]
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
#endif // __COMPAT_KINECTSDK__

#ifdef __COMPAT_KINECTSDK__
bool L_CapturadorImagenKinectSDK::capturarImagenOld()
{
	long i;
	DWORD width, height;
	HRESULT ret;

	WaitForSingleObject( imageEvent, INFINITE );

	CONST NUI_IMAGE_FRAME *imageFrame = 0;
	ret = NuiImageStreamGetNextFrame( imageHandle, 0, &imageFrame );
	if (ret != S_OK)
		return false;

	KINECT_LOCKED_RECT rect;
	ret = imageFrame->pFrameTexture->LockRect( 0, &rect, 0, 0 );
	if (ret != S_OK)
		return false;

	NuiImageResolutionToSize(imageFrame->eResolution, width, height);
	im.reallocate(width, height);
	for (i=0; i<(long)(width*height); i++)
	{
		im.data()[i][0] = ((unsigned char *)rect.pBits)[4*i+0];
		im.data()[i][1] = ((unsigned char *)rect.pBits)[4*i+1];
		im.data()[i][2] = ((unsigned char *)rect.pBits)[4*i+2];
	}

	ret = NuiImageStreamReleaseFrame( imageHandle, imageFrame );
	if (ret != S_OK)
		return false;

	WaitForSingleObject( depthEvent, INFINITE );
	imageFrame = 0;
	CONST NUI_IMAGE_FRAME *depthFrame = 0;
	ret = NuiImageStreamGetNextFrame( depthHandle, 0, &depthFrame );
	if (ret != S_OK)
		return false;
	ret = depthFrame->pFrameTexture->LockRect( 0, &rect, 0, 0 );
	if (ret != S_OK)
		return false;

	NuiImageResolutionToSize(depthFrame->eResolution, width, height);
	imDepth.reallocate(width, height);
	imKin.reallocate(width, height);
	memcpy(imKin.data(), rect.pBits, width*height*sizeof(unsigned short));

	for (i=0; i<(long)(width*height); i++)
	{
		imDepth.data()[i] = imKin.data()[i];//0x0fff & imKin.data()[i];
	}
	ret = NuiImageStreamReleaseFrame( depthHandle, depthFrame );
	if (ret != S_OK)
		return false;
	return true;
}
#endif // __COMPAT_KINECTSDK__


#ifdef __COMPAT_KINECTSDK__
bool L_CapturadorImagenKinectSDK::cerrar()
{
	NuiShutdown();
	return true;
}
#endif // __COMPAT_KINECTSDK__





L_CapturadorImagen::L_CapturadorImagen()
{
#ifdef __COMPAT_KINECTSDK__
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenKinectSDK();
#endif
#ifdef __COMPAT_OPENNI__
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenKinectNI();
#endif
#ifdef __COMPAT_IPLIMAGE__
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenCV();
#endif
#ifdef __COMPAT_SVS__
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenSVS();
#endif
	redimensionar.im_activo = false;
	redimensionar.im_resample = true;
	redimensionar.imDer_activo = false;
	redimensionar.imDer_resample = false;
	redimensionar.imProf_activo = false;
	redimensionar.imProf_resample = false;
}

L_CapturadorImagen::L_CapturadorImagen(const char *name)
{
	L_FileName ar;
	bool es_bmz;
	ar.asign(name);
	if (ar.hasExtension("bmz") == true)
		es_bmz = true;
	else
		es_bmz = false;

	if (es_bmz == true)
	{
		capt.push_back(new L_CapturadorImagenBMZ(name));
	}
#ifdef __COMPAT_IPLIMAGE__
	if (es_bmz == false)
	{
		capt.push_back(new L_CapturadorImagenCV(name));
	}
#endif
	if (capt.size() == 0)
		printf("L_CapturadorImagen() : No se puede abrir el archivo %s\n", name);
	redimensionar.im_activo = false;
	redimensionar.im_resample = true;
	redimensionar.imDer_activo = false;
	redimensionar.imDer_resample = false;
	redimensionar.imProf_activo = false;
	redimensionar.imProf_resample = false;
}


L_CapturadorImagen::L_CapturadorImagen(L_StringWithNumber name, int iIni)
{
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenesCarpeta(name, iIni);
	redimensionar.im_activo = false;
	redimensionar.im_resample = true;
	redimensionar.imDer_activo = false;
	redimensionar.imDer_resample = false;
	redimensionar.imProf_activo = false;
	redimensionar.imProf_resample = false;
}

L_CapturadorImagen::L_CapturadorImagen(L_StringWithNumber name, L_StringWithNumber nomKin, int iIni)
{
	capt.resize(capt.size()+1);
	capt[capt.size()-1]=new L_CapturadorImagenesCarpeta(name, nomKin, iIni);
	redimensionar.im_activo = false;
	redimensionar.im_resample = true;
	redimensionar.imDer_activo = false;
	redimensionar.imDer_resample = false;
	redimensionar.imProf_activo = false;
	redimensionar.imProf_resample = false;
}


bool L_CapturadorImagen::crear()
{
	// Caso lectura de imagen
	if (capt.size() == 0)
	{
		printf("L_CapturadorImagen::crear() : No hay capturadores disponibles\n");
		return false;
	}
	// Si es la primera vez, buscar cual de todos los capturadores entrega imagenes
	if (activo == true)
		return true;

	for (numActivo = 0; numActivo < (int)capt.size(); numActivo++)
	{
		if (capt[numActivo]->crear() == true)
		{
			activo = true;
			estereo.esEstereo = capt[numActivo]->estereo.esEstereo;
			prof.tieneProf = capt[numActivo]->prof.tieneProf;
			return true;
		}
	}
	// Ninguno se abre
	activo = false;
	numActivo = 0;
	capt.resize(0);
	return false;
}

bool L_CapturadorImagen::capturarImagen()
{
	// Caso lectura de imagen
	if (capt.size() == 0)
	{
		printf("L_CapturadorImagen::capturarImagen() : No hay capturadores disponibles\n");
		return false;
	}
	// Si es la primera vez, buscar cual de todos los capturadores entrega imagenes
	if (activo == false)
		crear();
	capt[numActivo]->estereo = estereo;
	capt[numActivo]->prof = prof;
	capt[numActivo]->redimensionar = redimensionar;
	if (capt[numActivo]->capturarImagen() == true)
	{
		im.swap(capt[numActivo]->im);
		imDer.swap(capt[numActivo]->imDer);
		imKin.swap(capt[numActivo]->imKin);
		imDepth.swap(capt[numActivo]->imDepth);
		imDepthRGB.swap(capt[numActivo]->imDepthRGB);
		redimensionarImagenes();
		return true;
	}
	return false;
}


bool L_CapturadorImagen::redimensionarImagenes()
{
	if (redimensionar.im_activo == true && im.data() != NULL && (im.lx != redimensionar.im_lx || im.ly != redimensionar.im_ly))
	{
		while (im.lx > redimensionar.im_lx && im.ly > redimensionar.im_ly)
			if (redimensionar.im_resample == true)
				im.halfResolution();
			else
				im.halfResolution_noResample();
		while (im.lx < redimensionar.im_lx && im.ly < redimensionar.im_ly)
			im.doubleResolution();
	}
	if (redimensionar.imDer_activo == true && imDer.data() != NULL && (imDer.lx != redimensionar.imDer_lx || im.ly != redimensionar.imDer_ly))
	{
		while (imDer.lx > redimensionar.imDer_lx && imDer.ly > redimensionar.imDer_ly)
			if (redimensionar.imDer_resample == true)
				imDer.halfResolution();
			else
				imDer.halfResolution_noResample();
		while (imDer.lx < redimensionar.imDer_lx && imDer.ly < redimensionar.imDer_ly)
			imDer.doubleResolution();
	}
	if (redimensionar.imProf_activo == true && imKin.data() != NULL && (imKin.lx != redimensionar.imProf_lx || imKin.ly != redimensionar.imProf_ly))
	{
		while (imKin.lx > redimensionar.imProf_lx && imKin.ly > redimensionar.imProf_ly)
			if (redimensionar.imDer_resample == true)
				imKin.halfResolution();
			else
				imKin.halfResolution_noResample();
		while (imKin.lx < redimensionar.imProf_lx && imKin.ly < redimensionar.imProf_ly)
			imKin.doubleResolution();
	}
	if (redimensionar.imProf_activo == true && imDepth.data() != NULL && (imDepth.lx != redimensionar.imProf_lx || imDepth.ly != redimensionar.imProf_ly))
	{
		while (imDepth.lx > redimensionar.imProf_lx && imDepth.ly > redimensionar.imProf_ly)
			if (redimensionar.imDer_resample == true)
				imDepth.halfResolution();
			else
				imDepth.halfResolution_noResample();
		while (imDepth.lx < redimensionar.imProf_lx && imDepth.ly < redimensionar.imProf_ly)
			imDepth.doubleResolution();
	}
	if (redimensionar.imProf_activo == true && imDepthRGB.data() != NULL && (imDepthRGB.lx != redimensionar.imProf_lx || imDepthRGB.ly != redimensionar.imProf_ly))
	{
		while (imDepthRGB.lx > redimensionar.imProf_lx && imDepthRGB.ly > redimensionar.imProf_ly)
			if (redimensionar.imDer_resample == true)
				imDepthRGB.halfResolution();
			else
				imDepthRGB.halfResolution_noResample();
		while (imDepthRGB.lx < redimensionar.imProf_lx && imDepthRGB.ly < redimensionar.imProf_ly)
			imDepthRGB.doubleResolution();
	}
	return true;
}

double L_CapturadorImagen::tiempoSeg()
{
	if (activo == true)
		return capt[numActivo]->tiempoSeg();
	return L_TIME();
}

const char *L_CapturadorImagen::nombreVideo()
{
	L_CapturadorImagenBMZ *ptrBMZ;
#ifdef __COMPAT_IPLIMAGE__
	L_CapturadorImagenCV *ptrCV;
#endif // __COMPAT_IPLIMAGE__

	if (activo == true)
	{
		if (capt[numActivo]->esArchivo())
		{
			ptrBMZ = dynamic_cast<L_CapturadorImagenBMZ *>(capt[numActivo]);
			if (ptrBMZ != NULL)
				return ptrBMZ->name.c_str();
#ifdef __COMPAT_IPLIMAGE__
			ptrCV = dynamic_cast<L_CapturadorImagenCV *>(capt[numActivo]);
			if (ptrCV != NULL)
				return ptrCV->name.c_str();
#endif // __COMPAT_IPLIMAGE__
		}
	}
	return NULL;
}

bool L_CapturadorImagen::fijarResolucion(int width, int height)
{
	int u;
	for (u = 0; u < (int)capt.size(); u++)
		capt[u]->fijarResolucion(width,height);
	return false;
}

bool L_CapturadorImagen::cerrar()
{
	int u;
	for (u=0; u<(int)capt.size(); u++)
		capt[u]->cerrar();
	numActivo = -1;
	activo = false;
	return true;
}

L_CapturadorImagen::~L_CapturadorImagen()
{
	int i;
	cerrar();
	for (i=0; i<(int)capt.size(); i++)
		delete capt[i];
}

#ifdef __COMPAT_FLTK__
int L_Fl_Double_Window::handle(int ev)
{
	if (ev==FL_PUSH || ev==FL_DRAG || ev==FL_RELEASE)
	{
		if (_mouseE==L_MouseLeido)
		{
			switch(ev)
			{
			case FL_PUSH:
				_mouseE = L_MousePresiona;
				break;
			case FL_RELEASE:
				_mouseE = L_MouseSuelta;
				break;
			case FL_DRAG:
				_mouseE = L_MouseArrastra;
				break;
			}
			_mouseX=Fl::event_x();
			_mouseY=Fl::event_y();
		}
		return 1;
	}
	return Fl_Double_Window::handle(ev);
}
#endif //__COMPAT_FLTK__



#ifdef __COMPAT_FLTK__
double L_PideNumFLTK(const char *mensaje, double porDefecto)
{
	const char *resp;
	resp=fl_input(mensaje, "%lf", porDefecto);
	if (resp==NULL || *resp == 0)
		return porDefecto;
	return atof(resp);
}
#endif

#ifdef __COMPAT_IPLIMAGE__
void L_cvCallback(int event, int x, int y, int flags, void* param)
{
	L_VentanaImagenCV *vent = (L_VentanaImagenCV *)param;
	switch(event)
	{
	case CV_EVENT_LBUTTONDOWN:
		vent->_mouseE = L_MousePresiona;
		break;
	case CV_EVENT_LBUTTONUP:
		vent->_mouseE = L_MouseSuelta;
		break;
	case CV_EVENT_MOUSEMOVE:
		vent->_mouseE = (vent->_mouseE == L_MousePresiona || vent->_mouseE == L_MouseArrastra) ? L_MouseArrastra : vent->_mouseE;
		break;
	}
	vent->_mouseX = x;
	vent->_mouseY = y;
}
#endif



#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaPunto3D_cubo(const L_CoordsCart3D &punto, long color, double radio)
{
	puntos.resize(puntos.size()+1);
	colores.resize(colores.size()+1);
	tipos.resize(tipos.size()+1);
	radios.resize(radios.size()+1);
	puntos[puntos.size()-1].pos=punto;
	puntos[puntos.size()-1].ori.fijaRotacionCero();
	colores[colores.size()-1]=color;
	tipos[tipos.size()-1]=L_cubo; // cubito
	radios[radios.size()-1] = radio;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaPunto3D_esfera(const L_CoordsCart3D &punto, long color, double radio)
{
	puntos.resize(puntos.size()+1);
	colores.resize(colores.size()+1);
	tipos.resize(tipos.size()+1);
	radios.resize(radios.size()+1);
	puntos[puntos.size()-1].pos=punto;
	puntos[puntos.size()-1].ori.fijaRotacionCero();
	colores[colores.size()-1]=color;
	tipos[tipos.size()-1]=L_esfera; // esfera
	radios[radios.size()-1] = radio;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaPunto3D_punto(const L_CoordsCart3D &punto, long color)
{
	puntos.resize(puntos.size()+1);
	colores.resize(colores.size()+1);
	tipos.resize(tipos.size()+1);
	radios.resize(radios.size()+1);
	puntos[puntos.size()-1].pos=punto;
	puntos[puntos.size()-1].ori.fijaRotacionCero();
	colores[colores.size()-1]=color;
	tipos[tipos.size()-1]=L_punto; // esfera
	radios[radios.size()-1] = 1;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaElipsoide3D(const double *cov3x3, const double *tr3x1, long color)
{
	L_Ellipsoid3D elip;
	elip.setCovarianceAndTranslation(cov3x3, tr3x1);
	elipsoides.push_back(elip);
	coloresElipsoides.push_back(color);
}
#endif // __COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaPose3D(const L_Pose3D_cuat &punto, long color, double radio)
{
	puntos.resize(puntos.size()+1);
	colores.resize(colores.size()+1);
	tipos.resize(tipos.size()+1);
	radios.resize(radios.size()+1);
	puntos[puntos.size()-1] = punto;
	colores[colores.size()-1]=color;
	tipos[tipos.size()-1] = L_sis_coords; // sis de coordenadas
	radios[radios.size()-1] = radio;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaLinea(const L_CoordsCart3D &ini, const L_CoordsCart3D &fin)
{
	lineas.resize(lineas.size() + 2);
	lineas[lineas.size()-2] = ini;
	lineas[lineas.size()-1] = fin;
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::agregaLineas(const std::vector<L_CoordsCart3D> &lins)
{
	int i, i0, nMax;
	nMax = (int)lins.size()-1;
	i0 = (int)lineas.size();
	lineas.resize(lineas.size()+2*nMax);  // Se debe agregar el inicio y fin
	for (i=0; i<nMax; i++)
	{
		lineas[i0 + 2*i+0] = lins[i+0];
		lineas[i0 + 2*i+1] = lins[i+1];  // no se pasa del ultimo supuestamente
	}
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw()
{
	int i;

	GLUquadricObj *quadratic;
	quadratic=gluNewQuadric();			// Create A Pointer To The Quadric Object ( NEW )
	gluQuadricNormals(quadratic, GLU_SMOOTH);	// Create Smooth Normals ( NEW )

	if (!valid())
	{
		// Fijar parametros de escena 3d
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-1, 1, -1, 1, 1, maxDistFrustrum); // izq, der, aba, arr, cerca, lejos
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glShadeModel(GL_SMOOTH);
		GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
		GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
		//GLfloat position[] = { -1.5f, 1.0f, -4.0f, 1.0f };
		//GLfloat position[] = { -10.0f, 0.0f, 30.0f, 1.0f };
		GLfloat position[4];
		position[0] = (GLfloat)posLuz.x;
		position[1] = (GLfloat)posLuz.y;
		position[2] = (GLfloat)posLuz.z;
		position[3] = (GLfloat)0.0;  // Luz al infinito

		// Assign created components to GL_LIGHT0
		glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
		glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
		draw_fijaPoseGL(L_Pose3D(0,0,0,0,0,0));
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		glEnable(GL_LIGHT0);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_COLOR_MATERIAL);
	glColor4f(0.8f, 0.8f, 0.0f,1.0f);
	draw_fijaPoseGL(poseRobot);
	draw_dibRobot(radioRobot);

	for (i=0; i<(int)puntos.size(); i++)
	{
		if (tipos[i] == L_cubo)  // Cubo
		{
			glColor4f(((colores[i])&0xFF)/255.0f, ((colores[i]>>8)&0xFF)/255.0f, ((colores[i]>>16)&0xFF)/255.0f,1.0f);
			draw_fijaPoseGL(puntos[i]);
			draw_dibCubo(radios[i]);
		}
		if (tipos[i] == L_esfera)  // Esfera
		{
			glColor4f(((colores[i])&0xFF)/255.0f, ((colores[i]>>8)&0xFF)/255.0f, ((colores[i]>>16)&0xFF)/255.0f,1.0f);
			draw_fijaPoseGL(puntos[i]);
			draw_dibEsfera(radios[i], quadratic);
		}
		else if (tipos[i] == L_sis_coords)  // Sis coordenadas
		{
			glColor4f(((colores[i])&0xFF)/255.0f, ((colores[i]>>8)&0xFF)/255.0f, ((colores[i]>>16)&0xFF)/255.0f,1.0f);
			draw_fijaPoseGL(puntos[i]);
			draw_dibSisCoords(radios[i], radios[i] / 15);
		}
		else if (tipos[i] == L_punto)  // Punto
		{
			glColor4f(((colores[i])&0xFF)/255.0f, ((colores[i]>>8)&0xFF)/255.0f, ((colores[i]>>16)&0xFF)/255.0f,1.0f);
			draw_fijaPoseGL(puntos[i]);
			glBegin(GL_POINTS);
			glVertex3f(0,0,0);
			glEnd();
		}
	}

	for (i=0; i<(int)elipsoides.size(); i++)
	{
		glColor4f(((coloresElipsoides[i])&0xFF)/255.0f, ((coloresElipsoides[i]>>8)&0xFF)/255.0f, ((coloresElipsoides[i]>>16)&0xFF)/255.0f,((coloresElipsoides[i]>>24)&0xFF)/255.0f);
		L_HomogeneousMatrix HT;
		HT.inversaTrasponiendoDe(elipsoides[i].H);
		draw_fijaPoseGL(HT);
		draw_dibElipsoide(elipsoides[i].rx(), elipsoides[i].ry(), elipsoides[i].rz());
	}

	L_Pose3D_cuat ident;
	ident.fijaCero();
	draw_fijaPoseGL(ident);
	draw_dibLineas();
	gluDeleteQuadric(quadratic);
	glDisable(GL_COLOR_MATERIAL);

	glFlush();
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_fijaPoseGL(const L_Pose3D &pose)
{
	L_HomogeneousMatrix m1, m2;
	L_StaticMatrix<4,4> m3;
	// mGL es la matriz que lleva del sistema de la camara al sistema central de OpenGL
	m1.fijaPose3D(pose); // Del sis local al de la camara
	L_StaticMatrix_OP_mult4x4(m2,mGL,m1); // Del sis. local al de la camara, y del de la camara al de OpenGL.
	L_StaticMatrix_transpOf4x4(m3,m2); // Transponer para cumplir el formato de OpenGL
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(&m3(0,0));
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_fijaPoseGL(const L_Pose3D_cuat &pose)
{
	L_HomogeneousMatrix m1, m2;
	L_StaticMatrix<4,4> m3;
	L_MatrizHomogenea_fijaPose3D_cuat(m1,pose); // Del sis local al de la camara
	L_StaticMatrix_OP_mult4x4(m2,mGL,m1); // Del sis. local al de la camara, y del de la camara al de OpenGL.
	L_StaticMatrix_transpOf4x4(m3,m2);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(&m3(0,0));
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_fijaPoseGL(const L_HomogeneousMatrix &pose)
{
	L_HomogeneousMatrix m2;
	L_StaticMatrix<4,4> m3;
	L_StaticMatrix_OP_mult4x4(m2,mGL,pose); // Del sis. local al de la camara, y del de la camara al de OpenGL.
	L_StaticMatrix_transpOf4x4(m3,m2);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(&m3(0,0));
}
#endif // __COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibCubo(double radio)
{
	GLfloat rad=(GLfloat)radio;
	glBegin(GL_QUADS);

	// cara adelante
	glNormal3f(1,0,0);
	glVertex3f(rad,-rad,rad);
	glNormal3f(1,0,0);
	glVertex3f(rad,-rad,-rad);
	glNormal3f(1,0,0);
	glVertex3f(rad,rad,-rad);
	glNormal3f(1,0,0);
	glVertex3f(rad,rad,rad);
	// cara atras
	glNormal3f(-1,0,0);
	glVertex3f(-rad,-rad,rad);
	glNormal3f(-1,0,0);
	glVertex3f(-rad,rad,rad);
	glNormal3f(-1,0,0);
	glVertex3f(-rad,rad,-rad);
	glNormal3f(-1,0,0);
	glVertex3f(-rad,-rad,-rad);
	// cara izq
	glNormal3f(0,1,0);
	glVertex3f(rad,rad,rad);
	glNormal3f(0,1,0);
	glVertex3f(rad,rad,-rad);
	glNormal3f(0,1,0);
	glVertex3f(-rad,rad,-rad);
	glNormal3f(0,1,0);
	glVertex3f(-rad,rad,rad);
	// cara der
	glNormal3f(0,-1,0);
	glVertex3f(rad,-rad,rad);
	glNormal3f(0,-1,0);
	glVertex3f(-rad,-rad,rad);
	glNormal3f(0,-1,0);
	glVertex3f(-rad,-rad,-rad);
	glNormal3f(0,-1,0);
	glVertex3f(rad,-rad,-rad);
	// cara arr
	glNormal3f(0,0,1);
	glVertex3f(rad,rad,rad);
	glNormal3f(0,0,1);
	glVertex3f(-rad,rad,rad);
	glNormal3f(0,0,1);
	glVertex3f(-rad,-rad,rad);
	glNormal3f(0,0,1);
	glVertex3f(rad,-rad,rad);
	// cara aba
	glNormal3f(0,0,-1);
	glVertex3f(rad,rad,-rad);
	glNormal3f(0,0,-1);
	glVertex3f(rad,-rad,-rad);
	glNormal3f(0,0,-1);
	glVertex3f(-rad,-rad,-rad);
	glNormal3f(0,0,-1);
	glVertex3f(-rad,rad,-rad);
	glEnd();
}
#endif // __COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibParalelepipedo(double x0, double y0, double z0, double x1, double y1, double z1)
{
	// x0 < x1, y0 < y1, z0 < z1
	// Las cantidades mayores apuntan hacia los ejes
	GLfloat ade = (GLfloat) x1;
	GLfloat izq = (GLfloat) y1;
	GLfloat arr = (GLfloat) z1;

	GLfloat atr = (GLfloat) x0;
	GLfloat der = (GLfloat) y0;
	GLfloat aba = (GLfloat) z0;

	glBegin(GL_QUADS);
	// cara adelante:
	glNormal3f(1,0,0);
	glVertex3f(ade, izq, arr);
	glNormal3f(1,0,0);
	glVertex3f(ade, der, arr);
	glNormal3f(1,0,0);
	glVertex3f(ade, der, aba);
	glNormal3f(1,0,0);
	glVertex3f(ade, izq, aba);
	// cara atras
	glNormal3f(-1,0,0);
	glVertex3f(atr, izq, arr);
	glNormal3f(-1,0,0);
	glVertex3f(atr, izq, aba);
	glNormal3f(-1,0,0);
	glVertex3f(atr, der, aba);
	glNormal3f(-1,0,0);
	glVertex3f(atr, der, arr);
	// cara izq
	glNormal3f(0,1,0);
	glVertex3f(ade, izq, arr);
	glNormal3f(0,1,0);
	glVertex3f(atr, izq, arr);
	glNormal3f(0,1,0);
	glVertex3f(atr, izq, aba);
	glNormal3f(0,1,0);
	glVertex3f(ade, izq, aba);
	// cara der
	glNormal3f(0,-1,0);
	glVertex3f(ade, der, arr);
	glNormal3f(0,-1,0);
	glVertex3f(atr, der, arr);
	glNormal3f(0,-1,0);
	glVertex3f(atr, der, aba);
	glNormal3f(0,-1,0);
	glVertex3f(ade, der, aba);
	// cara arr
	glNormal3f(0,0,1);
	glVertex3f(ade, izq, arr);
	glNormal3f(0,0,1);
	glVertex3f(atr, izq, arr);
	glNormal3f(0,0,1);
	glVertex3f(atr, der, arr);
	glNormal3f(0,0,1);
	glVertex3f(ade, der, arr);
	// cara aba
	glNormal3f(0,0,-1);
	glVertex3f(ade, izq, aba);
	glNormal3f(0,0,-1);
	glVertex3f(ade, der, aba);
	glNormal3f(0,0,-1);
	glVertex3f(atr, der, aba);
	glNormal3f(0,0,-1);
	glVertex3f(atr, izq, aba);
	glEnd();
}
#endif //__COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibEsfera(double radio, GLUquadricObj *quadratic)
{
	if (quadratic == NULL)
	{
		quadratic=gluNewQuadric();			// Create A Pointer To The Quadric Object ( NEW )
		gluQuadricNormals(quadratic, GLU_SMOOTH);	// Create Smooth Normals ( NEW )
		gluSphere(quadratic,radio,16,16); // decia 32,32
		gluDeleteQuadric(quadratic);
	}
	else
		gluSphere(quadratic,radio,16,16); // decia 32,32
}
#endif //__COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibElipsoide(double rx, double ry, double rz, GLUquadricObj *quadratic)
{
	//printf("            %f  %f  %f\n", rx, ry, rz);
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	if (quadratic == NULL)
	{
		quadratic=gluNewQuadric();			// Create A Pointer To The Quadric Object ( NEW )
		gluQuadricNormals(quadratic, GLU_SMOOTH);	// Create Smooth Normals ( NEW )
		glScaled(1/rx, 1/ry, 1/rz);
		gluSphere(quadratic,1,16,16); // decia 32,32
		gluDeleteQuadric(quadratic);
	}
	else
	{
		glScaled(1/rx, 1/ry, 1/rz);
		gluSphere(quadratic,1,16,16); // decia 32,32
	}
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPopMatrix();
}
#endif //__COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibLineas()
{
	int i;
	glBegin(GL_LINES);
	glColor3d(1.0, 1.0, 1.0);
	for (i=0; i<(int)lineas.size()/2; i++)
	{
		glVertex3d(lineas[2*i+0].x, lineas[2*i+0].y, lineas[2*i+0].z);
		glVertex3d(lineas[2*i+1].x, lineas[2*i+1].y, lineas[2*i+1].z);
	}
	glEnd();
}
#endif //__COMPAT_FLTK__

#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibRobot(double radio)
{
	double m[]={
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0, // tAde, tIzq, tArr, 1.0
	};
	draw_dibCubo(radio);
	m[12]=radio;
	glMultMatrixd(m);
	draw_dibCubo(radio/2);
	m[12]=-radio;
	glMultMatrixd(m);
}
#endif //__COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::draw_dibSisCoords(double radio, double grosor)
{
	draw_dibParalelepipedo(0, -grosor, -grosor, radio, grosor, grosor); // El eje x
	draw_dibParalelepipedo(-grosor, 0, -grosor, grosor, radio, grosor); // El eje y
	draw_dibParalelepipedo(-grosor, -grosor, 0, grosor, grosor, radio); // El eje z
	draw_dibParalelepipedo(radio-2*grosor, -2*grosor, -2*grosor, radio+2*grosor, +2*grosor, +2*grosor); // Cubito para el eje x
}
#endif //__COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
void L_VentMapa3D::test()
{
	L_Pose3D puntitos[100];
	long coloritos[100];
	L_Pose3D pose;
	long i, j;

	show();

	for (i=0; i<100; i++)
	{
		puntitos[i].pos.x=(rand()%1000 - 500)*0.2;
		puntitos[i].pos.y=(rand()%1000 - 500)*0.2;
		puntitos[i].pos.z=(rand()%1000 - 500)*0.2;
		puntitos[i].ori.fijaRotacionCero();
		coloritos[i] = ((rand()%255)) +  ((rand()%255)<<8) + ((rand()%255)<<16);
	}
	
	pose.fijaPose(10,0,0,0,0,0);
	fijaPoseCamara(L_Pose3D().fijaPose_abs(-50,0,20, 0,0,0, -50,0,50)); // La pose de la camara, recalcula la transformacion camara-OpenGL
	fijaPosLuz(L_CoordsCart3D(-300.0, 0.0, 100.0));

	while(1)
	{
		clear();
		pose.ori.pan+=0.01;
		fijaPoseRobot(pose);
		puntos.resize(100);
		colores.resize(100);
		for (j=0; j<100; j++)
		{
			puntitos[j].pos.x+=(rand()%1000 - 500)*0.002;
			puntitos[j].pos.y+=(rand()%1000 - 500)*0.002;
			puntitos[j].pos.z+=(rand()%1000 - 500)*0.002;
			puntitos[j].ori.pan+=0.01;
			puntitos[j].ori.tilt+=0.01;
			puntos[j].pos=puntitos[j].pos;
			colores[j]=coloritos[j];
		}
		redraw();
		Fl::check();
	}
}
#endif //__COMPAT_FLTK__



bool L_CurvaPuntos::leerBarridoTxt(FILE *fp, int nPuntos)
{
	int i;
	int nleidos = 0;
	resize(nPuntos);

	for (i=0; i<nPuntos; i++)
	{
		nleidos += fscanf(fp, "%lf%lf", &v[i].x, &v[i].y);
		v[i].z = 0;
	}

	return nleidos == nPuntos * 2;
}

bool L_CurvaPuntos::guardarBarridoTxt(FILE *fp)
{
	int i;
	for (i=0; i<(int)size(); i++)
		fprintf(fp, "%.1lf %.1lf ", v[i].x, v[i].y);
	fprintf(fp, "\n");
	return true;
}

bool L_CurvaPuntos::guardarBarridoTxt2col(FILE *fp)
{
	int i;
	for (i=0; i<(int)size(); i++)
		fprintf(fp, "%.1f %.1f\n", v[i].x, v[i].y);
	return true;
}

void L_CurvaPuntos::OP_mult(L_HomogeneousMatrix &H, L_CurvaPuntos &curva)
{
	int i;
	resize(curva.size());
	for (i=0; i<(int)size(); i++)
		L_MatrizHomogenea_OP_mult_P(v[i],H,curva[i]);
}



bool L_ArregloCurvaPuntos::leerTxt(FILE *fp, int nPuntos)
{
	resize_swapping(1);

	while ((*this)[size()-1].leerBarridoTxt(fp, nPuntos))
	{
		resize_swapping(size()+1);
	}
	resize_swapping(size()-1);
	if (size()==0)
		return false;
	return true;
}

bool L_ArregloCurvaPuntos::guardarTxt(FILE *fp)
{
	int i;
	for (i=0; i<(int)size(); i++)
		(*this)[i].guardarBarridoTxt(fp);
	return true;
}

bool L_ArregloCurvaPuntos::guardarTxt2col(FILE *fp)
{
	int i;
	for (i=0; i<(int)size(); i++)
		(*this)[i].guardarBarridoTxt2col(fp);
	return true;
}

bool L_ArregloCurvaPuntos::transformar_dat_a_txt(const char *nombreDat)
{
	L_FileName ar;
	FILE *fpin, *fpout;
	L_ArregloLaserDato laser;
	L_ArregloCurvaPuntos puntos;
	ar.asign(nombreDat);
	ar.ext = "txt";
	fpin = fopen(nombreDat, "rb");
	if (fpin == NULL)
		return false;
	laser.leerDat(fpin);
	fclose(fpin);

	laser.copiarPuntosEn(puntos);

	fpout = fopen(ar.dirnameext(), "w");
	if (fpout == NULL)
		return false;
	puntos.guardarTxt2col(fpout);
	fclose(fpout);
	return true;
}


void L_ArregloLaserDato::leerDat(FILE *fp)
{
	resize(0);
	while (!feof(fp))
	{
		resize(size()+1);
		if (fread(&(*this)[size()-1], sizeof(L_LaserDato), 1, fp) != 1)
			{printf("L_ArregloLaserDato::leerDat() : formato de archivo incorrecto"); return;}
	}
}
void L_ArregloLaserDato::grabarDat(FILE *fp)
{
	int i;
	for (i=0; i<(int)size(); i++)
		fwrite(&(*this)[i], sizeof(L_LaserDato), 1, fp);
}
void L_ArregloLaserDato::guardarPuntosEn(L_CurvaPuntos &arr, int nBarrido)
{
	int i;
	arr.resize(MAXLASERPOINTS);
	double theta = -M_PI/2;
	for (i=0; i<MAXLASERPOINTS; i++)
	{
		arr[i].x = (*this)[nBarrido].rayos[i] * cos(theta);
		arr[i].y = (*this)[nBarrido].rayos[i] * sin(theta);
		arr[i].z = 0;
		theta += 0.5 * M_PI / 180;
	}
}

void L_ArregloLaserDato::copiarPuntosEn(L_ArregloCurvaPuntos &arr)
{
	int i;
	arr.resize_swapping(size());
	for (i=0; i<(int)size(); i++)
		guardarPuntosEn(arr[i], i);
}


void L_ICP::fnObjetivo_calcAsoc(const void *datos, L_Matrix &xytheta)
{
	int i, i1;
	double c, s, tx, ty;
	L_CoordsCart3D p1, p1t; // Transformacion de coordenadas y variables
	double di, di1, dim, diM;

	L_ICP_Dato *obj = (L_ICP_Dato *)datos;

	// Recordar que si *R,+t llevan los ejes {i,j} a {i',j'} => (x,y) = (RT t) (x',y')
	tx = xytheta(0,0);
	ty = xytheta(1,0);
	c = cos(xytheta(2,0));
	s = sin(xytheta(2,0));

	for (i=0; i<(int)obj->arr1->size(); i++)
	{
		p1 = (*obj->arr1)[i];
		p1t.x = c*p1.x - s*p1.y + tx;
		p1t.y = s*p1.x + c*p1.y + ty;
		p1t.z = 0;
		i1 = obj->arb2->buscaMasCercano(p1t);
		di = p1t.distancia2A((*obj->arr2)[(*obj->asoc)[i]]);
		// Truco: (a,b) == b
		// Pequeno test para no empeorar las asociaciones
		if (i1 != -1 && (di1 = p1t.distancia2A((*obj->arr2)[i1]), di1 < di))
		{
			(*obj->asoc)[i] = i1;
			di = di1;
		}
		// Buscar en los puntos anterior y posterior por si hay un calce mejor
		if ((*obj->asoc)[i]>0 &&  (dim = p1t.distancia2A((*obj->arr2)[(*obj->asoc)[i]-1]), dim < di))
		{
			(*obj->asoc)[i]--;
			di = dim;
		}
		if ((*obj->asoc)[i]<(int)obj->arr1->size()-1 && (diM = p1t.distancia2A((*obj->arr2)[(*obj->asoc)[i]+1]), diM < di))
		{
			(*obj->asoc)[i]++;
			di = diM;
		}
	}
}

void L_ICP::fnObjetivo_gradiente(const void *datos, const L_Matrix &xytheta, L_Matrix &errVect)
{
	int i, j, signo;
	double c, s, tx, ty;
	L_CoordsCart3D p1, p1t, p2, resta; // Transformacion de coordenadas y variables

	L_ICP_Dato *obj = (L_ICP_Dato *)datos;
	errVect.reallocate((int)obj->arr1->size(), 1);
	// Recordar que si *R,+t llevan los ejes {i,j} a {i',j'} => (x,y) = (RT t) (x',y')
	tx = xytheta(0,0);
	ty = xytheta(1,0);
	c = cos(xytheta(2,0));
	s = sin(xytheta(2,0));

	for (i=0; i<(int)obj->arr1->size(); i++)
	{
		p1 = (*obj->arr1)[i];
		p1t.x = c*p1.x - s*p1.y + tx;
		p1t.y = s*p1.x + c*p1.y + ty;
		p1t.z = 0;
		j = (*obj->asoc)[i];
		p2 = (*obj->arr2)[j];
		L_CoordsCart3D_OP_resta(resta,p2,p1t);
		errVect(i,0) = L_CoordsCart3D_punto(resta, (*obj->per2)[j]); // Proyeccion sobre la perpendicular
		signo = 2*(errVect(i,0)>=0)-1;
		if (errVect(i,0)*signo > obj->distMax)
			errVect(i,0) = obj->distMax*signo;
	}
}


void L_ICP::fnObjetivo_puntos(const void *datos, const L_Matrix &xytheta, L_Matrix &errVect)
{
	int i, j;
	double c, s, tx, ty, signo1, signo2;
	L_CoordsCart3D p1, p1t, p2, resta; // Transformacion de coordenadas y variables

	L_ICP_Dato *obj = (L_ICP_Dato *)datos;
	errVect.reallocate(2*(int)obj->arr1->size(), 1);
	// Recordar que si *R,+t llevan los ejes {i,j} a {i',j'} => (x,y) = (RT t) (x',y')
	tx = xytheta(0,0);
	ty = xytheta(1,0);
	c = cos(xytheta(2,0));
	s = sin(xytheta(2,0));

	for (i=0; i<(int)obj->arr1->size(); i++)
	{
		p1 = (*obj->arr1)[i];
		p1t.x = c*p1.x - s*p1.y + tx;
		p1t.y = s*p1.x + c*p1.y + ty;
		p1t.z = 0;
		j = (*obj->asoc)[i];
		p2 = (*obj->arr2)[j];
		L_CoordsCart3D_OP_resta(resta,p2,p1t);
		errVect(2*i+0,0) = resta.x; // Proyeccion sobre la perpendicular
		errVect(2*i+1,0) = resta.y; // Proyeccion sobre la perpendicular

		signo1 = (errVect(2*i+0,0)>=0)*2-1;
		signo2 = (errVect(2*i+1,0)>=0)*2-1;
		if (errVect(2*i+0,0)*signo1 > obj->distMax)
			errVect(2*i+0,0) = obj->distMax*signo1;
		if (errVect(2*i+1,0)*signo2 > obj->distMax)
			errVect(2*i+1,0) = obj->distMax*signo2;
	}
}


double L_ICP::ICP_lineal_norm_2d(L_ICP_Dato *datos, L_Matrix &xytheta)
{
	// Kok-Lim Low, "Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration"
	int it;
	int i;
	int ene = (int)datos->arr1->size();
	L_CoordsCart3D s, d, n;
	L_Pose3D pose;
	double alfa=0, beta=0, gamma=0, tx=0, ty=0, tz=0;
	double e, e1, e2;
	L_Matrix A(ene,3),x(3,1),b(ene,1),AT(3,ene), ATA(3,3), ATAi(3,3), ATb(3,1);
	L_HomogeneousMatrix H, Htot;
	std::vector<L_CoordsCart3D> arrMov(ene);
	arrMov = *datos->arr1;
	xytheta(0,0) = 0;
	xytheta(1,0) = 0;
	xytheta(2,0) = 0;
	Htot.fijaIdent();

	for (it=0; it<20; it++)
	{
		fnObjetivo_calcAsoc(datos, xytheta);
		e1 = 0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			n = (*datos->per2)[(*datos->asoc)[i]];
			A(i,0) = n.y*s.x - n.x*s.y;
			A(i,1) = n.x;
			A(i,2) = n.y;
			b(i,0) = n.x*d.x + n.y*d.y - n.x*s.x - n.y*s.y;
			e = n.punto(s-d);
			e1 +=e*e;
		}
		AT.transpOf(A);
		ATA.OP_mult(AT,A);
		ATAi.inverseOf(ATA);
		ATb.OP_mult(AT,b);
		x.OP_mult(ATAi,ATb);
		// [alfa beta gamma tx ty tz)
		gamma = x(0,0);
		tx = x(1,0);
		ty = x(2,0);
		// Recuperacion de la rotacion
		H(0,0) = 1;
		H(0,1) = -gamma;
		H(0,2) = beta;
		H(1,0) = gamma;
		H(1,1) = 1;
		H(1,2) = -alfa;
		H(2,0) = -beta;
		H(2,1) = alfa;
		H(2,2) = 1;
		H(0,3) = tx;
		H(1,3) = ty;
		H(2,3) = tz;
		H.normalizaSVD(); // Ahora la transfromacion afin se aproxima a rotacion
		Htot = H*Htot;
		// Mover los puntos iniciales
		for (i=0; i<ene; i++)
			arrMov[i] = Htot * (*datos->arr1)[i];
		pose = Htot.calcPose3D();
		e2=0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			n = (*datos->per2)[(*datos->asoc)[i]];
			e = n.punto(s-d);
			e2 +=e*e;
		}
		if (e2 > e1) // Para que empeorar los datos...
			break;
		xytheta(0,0) = pose.pos.x;
		xytheta(1,0) = pose.pos.y;
		xytheta(2,0) = pose.ori.pan;
	}
	return e1;
}


double L_ICP::ICP_lineal_norm_3d(L_ICP_Dato *datos, L_Matrix &xytheta)
{
	// Kok-Lim Low, "Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration"
	int it;
	int i;
	int ene = (int)datos->arr1->size();
	L_CoordsCart3D s, d, n;
	L_Pose3D pose;
	double alfa, beta, gamma, tx, ty, tz;
	double e, e1, e2;
	L_Matrix A(ene,6),x(6,1),b(ene,1),AT(6,ene), ATA(6,6), ATAi(6,6), ATb(6,1);
	L_HomogeneousMatrix H, Htot;
	std::vector<L_CoordsCart3D> arrMov(ene);
	arrMov = *datos->arr1;
	xytheta(0,0) = 0;
	xytheta(1,0) = 0;
	xytheta(2,0) = 0;
	Htot.fijaIdent();

	for (it=0; it<10; it++)
	{
		fnObjetivo_calcAsoc(datos, xytheta);
		e1 = 0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			n = (*datos->per2)[(*datos->asoc)[i]];
			A(i,0) = n.z*s.y - n.y*s.z;
			A(i,1) = n.x*s.z - n.z*s.x;
			A(i,2) = n.y*s.x - n.x*s.y;
			A(i,3) = n.x;
			A(i,4) = n.y;
			A(i,5) = n.z;
			b(i,0) = n.x*d.x + n.y*d.y + + n.z*d.z - n.x*s.x - n.y*s.y - n.z*s.z;
			e = n.punto(s-d);
			e1 +=e*e;
		}
		AT.transpOf(A);
		ATA.OP_mult(AT,A);
		ATAi.inverseOf(ATA);
		ATb.OP_mult(AT,b);
		x.OP_mult(ATAi,ATb);
		// [alfa beta gamma tx ty tz]
		alfa = x(0,0);
		beta = x(1,0);
		gamma = x(2,0);
		tx = x(3,0);
		ty = x(4,0);
		tz = x(5,0);
		// Recuperacion de la rotacion
		H(0,0) = 1;
		H(0,1) = -gamma;
		H(0,2) = beta;
		H(1,0) = gamma;
		H(1,1) = 1;
		H(1,2) = -alfa;
		H(2,0) = -beta;
		H(2,1) = alfa;
		H(2,2) = 1;
		H(0,3) = tx;
		H(1,3) = ty;
		H(2,3) = tz;
		H.normalizaSVD(); // Ahora la transformacion afin se aproxima a rotacion
		Htot = H*Htot;
		// Mover los puntos iniciales
		for (i=0; i<ene; i++)
			arrMov[i] = H*arrMov[i];
		pose = Htot.calcPose3D();
		e2=0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			n = (*datos->per2)[(*datos->asoc)[i]];
			e = n.punto(s-d);
			e2 +=e*e;
		}
		if (e2 > e1) // Para que empeorar los datos...
			break;
		xytheta(0,0) = pose.pos.x;
		xytheta(1,0) = pose.pos.y;
		xytheta(2,0) = pose.ori.pan;
	}
	return e1;
}


double L_ICP::ICP_lineal_pun_2d(L_ICP_Dato *datos, L_Matrix &xytheta)
{
	// Kok-Lim Low, "Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration"
	int it;
	int i;
	int ene = (int)datos->arr1->size();
	L_CoordsCart3D s, d, n;
	L_Pose3D pose;
	double alfa=0, beta=0, gamma=0, tx=0, ty=0, tz=0;
	double e1, e2;
	L_Matrix A(2*ene,3),x(3,1),b(2*ene,1),AT(3,2*ene), ATA(3,3), ATAi(3,3), ATb(3,1);
	L_HomogeneousMatrix H, Htot;
	std::vector<L_CoordsCart3D> arrMov(ene);
	arrMov = *datos->arr1;
	xytheta(0,0) = 0;
	xytheta(1,0) = 0;
	xytheta(2,0) = 0;
	Htot.fijaIdent();

	for (it=0; it<20; it++)
	{
		fnObjetivo_calcAsoc(datos, xytheta);
		e1 = 0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			A(2*i+0,0) = 1;
			A(2*i+0,1) = 0;
			A(2*i+0,2) = s.y;
			A(2*i+1,0) = 0;
			A(2*i+1,1) = 1;
			A(2*i+1,2) = -s.x;
			b(2*i+0,0) = d.x-s.x;
			b(2*i+1,0) = d.y-s.y;
			e1 +=L_CoordsCart3D_distancia2A(s,d);
		}
		AT.transpOf(A);
		ATA.OP_mult(AT,A);
		ATAi.inverseOf(ATA);
		ATb.OP_mult(AT,b);
		x.OP_mult(ATAi,ATb);
		// [alfa beta gamma tx ty tz)
		gamma = x(0,0);
		tx = x(1,0);
		ty = x(2,0);
		// Recuperacion de la rotacion
		H(0,0) = 1;
		H(0,1) = -gamma;
		H(0,2) = beta;
		H(1,0) = gamma;
		H(1,1) = 1;
		H(1,2) = -alfa;
		H(2,0) = -beta;
		H(2,1) = alfa;
		H(2,2) = 1;
		H(0,3) = tx;
		H(1,3) = ty;
		H(2,3) = tz;
		H.normalizaSVD(); // Ahora la transformacion afin se aproxima a rotacion
		Htot = H*Htot;
		// Mover los puntos iniciales
		for (i=0; i<ene; i++)
			arrMov[i] = Htot * (*datos->arr1)[i];
		pose = Htot.calcPose3D();
		e2=0;
		for (i=0; i<ene; i++)
		{
			s = arrMov[i];
			d = (*datos->arr2)[(*datos->asoc)[i]];
			e2 +=L_CoordsCart3D_distancia2A(s,d);
		}
		if (e2 > e1) // Para que empeorar los datos...
			break;
		xytheta(0,0) = pose.pos.x;
		xytheta(1,0) = pose.pos.y;
		xytheta(2,0) = pose.ori.pan;
	}
	return e1;
}

double L_ICP::calcICP(L_CurvaPuntos &inic, L_CurvaPuntos &final, double dMax, L_Pose3D &pose, int nMetodo, L_CurvaPuntos *normales)
{
	L_Matrix v(3,1);
	L_CurvaPuntos gradPerFin;
	L_KdTree3D arbFin;
	std::vector<int> asoc(inic.size());
	L_LevenbergMarquardt opt;
	L_ICP_Dato data;
	int i;

	if (normales == NULL)
		normales = &gradPerFin;

	data.distMax = dMax;
	data.arr1 = &inic.v;
	data.arr2 = &final.v;
	data.arb2 = &arbFin;
	data.per2 = &normales->v;
	data.asoc = &asoc;

	opt.epsilon = 1e-5;
	opt.errorToStop = 1e-3;
	opt.lenVector = 3;
	opt.nIterationsMax = 20;

	v(0,0) = 0;
	v(1,0) = 0;
	v(2,0) = 0;

	arbFin.createFrom(final.v);

	for (i=0; i<(int)asoc.size(); i++) // Asociaciones iniciales, mejoran usando el kd-tree
		asoc[i] = i;

	switch(nMetodo)
	{
	case 0:
		opt.minimize_vect(&data, v, &L_ICP::fnObjetivo_puntos, NULL, &L_ICP::fnObjetivo_calcAsoc);
		break;
	case 1:
		final.calcularNormales2D(*normales, arbFin);
		opt.minimize_vect(&data, v, &L_ICP::fnObjetivo_gradiente, NULL, &L_ICP::fnObjetivo_calcAsoc);
		break;
	case 2:
		ICP_lineal_pun_2d(&data, v);
		break;
	case 3:
		final.calcularNormales2D(*normales, arbFin);
		ICP_lineal_norm_2d(&data, v);
		break;
	};

	pose.pos.x = v(0,0);
	pose.pos.y = v(1,0);
	pose.pos.z = 0;
	pose.ori.pan = v(2,0);
	pose.ori.tilt = 0;
	pose.ori.roll = 0;
	return 0;
}

int L_ICP::main_icp(int argc, char *argv[])
{
	if (argc != 3)
	{
		printf(" Hay 2 opciones de uso:\n");
		printf("LICP puntos361.txt puntosCorregidos361.txt\n");
		printf("LICP puntos361.dat puntosCorregidos361.txt\n");
		printf("  Formato del texto: 2 columnas de 361*N puntos\n");
		return 1;
	}
	L_FileName ar;
	L_ArregloCurvaPuntos arr1, arr2;
	L_HomogeneousMatrix H, Hac, Haci;
	int i;
	L_Pose3D pose;

	ar.asign(argv[1]);
	if (ar.hasExtension("dat"))
	{
		FILE *fp = fopen(argv[1], "rb");
		if (fp == NULL)
			{printf("No esta el archivo %s\n", argv[1]); return 0;}
		L_ArregloLaserDato laser;
		laser.leerDat(fp);
		fclose(fp);
		laser.copiarPuntosEn(arr1);
	}
	else
	{
		FILE *fp = fopen(argv[1], "r");
		if (fp == NULL)
			{printf("No esta el archivo %s\n", argv[1]); return 0;}
		arr1.leerTxt(fp, 361);
		fclose(fp);
	}

	arr2.resize(arr1.size());
	
	H.fijaIdent();
	Hac.fijaIdent();
	
	L_ArregloCurvaPuntos arrNorm(arr1.size());

	// p2 = H*p1
	// p2 = (H1*H2)*p1

	for (i=0; i<(int)arr1.size()-1; i++)
	{
		L_ICP::calcICP(arr1[i], arr1[i+1], 100, pose, 1);
		H.fijaPose3D(pose);
		Hac = H*Hac; // ptos_actuales = Hac*ptos_iniciales
		Haci.inversaTrasponiendoDe(Hac);
		arr2[i].OP_mult(Haci, arr1[i]);
		if (i%(arr1.size()/100) == 0)
		{
			L_Pose3D p3d = Hac.calcPose3D();
			printf("\r%d %%  pose actual: %.2f[cm]  %.2f[cm]  %.2f[grados]      ", i*100/arr1.size(), p3d.pos.x, p3d.pos.y, p3d.ori.pan*180/M_PI);
		}
	}
	FILE *output = fopen(argv[2], "w");
	if (output == NULL)
		{printf("No esta el archivo %s\n", argv[2]); return 0;}
	arr2.guardarTxt2col(output);
	fclose(output);
	return 0;
}





int L_ProgramasFnesPato::main_L_for(int argc, char *argv[], bool debug)
{
	L_Array<L_ContInfo> rang;
	L_Array<L_TokenInfo> tok;
	L_ContInfo tmp;
	L_TokenInfo tmpTok;
	std::vector<L_MatParsPlaceholder> reemplazos;
	L_MatParsPlaceholder tmpReem;
	char c;
	char c2[2]={0,0};
	int i, n1, dn, n2, l;
	long numLong;
	char expr[400]="";
	int enExpr=0;

	if (argc<3)
	{
		printf("\n\n  L_for /i=n1:dn:n2 /j=n3:n4 ... comando\n\n"
			"Ejecuta \"comando\", reemplazando los $(i) en \"comando\" por su value.\n"
			"Se pueden usar expresiones matematicas dentro de los parentesis: $(f(i))\n"
			"Los rangos se definen como en el comando for de matlab\n\n");
		return 1;
	}

	// Almacenar los rangos en "rang" y las variables en "reemplazos"
	// Cada rango tiene un (c,n1,n2,width), donde "c" es la letra de la variable, y "width" es el width que ocupa el numero reemplazado
	for ( i=1 ; i<argc && *argv[i]=='/' ; i++ )
	{
		if ( sscanf(argv[i],"/%c=%d:%d:%d", &c, &n1, &dn, &n2) != 4)
		{
			dn=1;
			if ( sscanf(argv[i],"/%c=%d:%d", &c, &n1, &n2) != 3)
			{
				printf("Error en formato en %s\n", argv[i]);
				return 2;
			}
		}
		l=3;
		tmp.width=0;
		if (argv[i][l]=='-' || argv[i][l]=='+')
			l++;
		while(argv[i][l]>='0' && argv[i][l]<='9')
		{
			l++;
			tmp.width++;
		}
		tmp.str[0]=c;
		tmp.n=n1;
		tmp.n1=n1;
		tmp.n2=n2;
		tmp.dn=dn;

		tmpReem.varNom=tmp.str;
		tmpReem.varRe=n1;

		reemplazos.push_back(tmpReem);
		tmp.idx=(int)reemplazos.size()-1;
		rang.push_back(tmp);
	}

	// Si no hay rangos, salir
	if (rang.size()==0)
		return 0;

	// Copiar el comando a evaluar en expr
	for ( ; i<argc ; i++)
	{
		if (strlen(expr)+strlen(argv[i])+1>390)
		{
			printf("comando muy length\n");
			return 2;
		}
		strcat(expr, " ");
		strcat(expr, argv[i]);
	}

	// Descomponer el comando "expr" en tokens (string-texto o string-expresion), y guardarlos en "tok"
	l=(int)strlen(expr);
	tmpTok.esTexto=true;
	n1=0;
	enExpr=0;
	for (i=0; i<l; i++)
	{
		n2=i;
		if (expr[i]=='$')
		{
			// Guardar lo que se llevaba hasta ahora (tmpTok) como string no matematico (texto)
			tmpTok.esTexto=true;
			tmpTok.str.copySubString(expr, n1, n2-1); // Hasta antes del $
			tok.push_back(tmpTok);
			// Comenzar el parseo de la expresion
			if (i<l-1 && expr[i+1]=='(')
			{
				n1=i;
				enExpr=1;
				i++;
			}
			else if (i<l-1 && expr[i+1]=='$') // Caso $$: dejarlo como "$"
			{
				n1=i+1;
				i++;
			}
			else // caso $i: value directo de la variable, guardar como string matematico
			{
				tmpTok.esTexto=false;
				tmpTok.mat.destroy();
				c2[0]=expr[i+1];
				tmpTok.str = c2;
				tmpTok.mat.buildFrom(tmpTok.str.c_str());
				//tmpTok.mat.showStructure();
				tmpTok.mat.reemplazos=&reemplazos;
				tok.push_back(tmpTok);
				n1=i+2;
				n2=i+2;
				enExpr=0;
				i++;
			}
		}
		else if (enExpr==1 && expr[i]==')') // Se cierra el ultimo parentesis, guardar lo que se lleva como string matematico
		{
			tmpTok.esTexto=false;
			tmpTok.str.copySubString(expr, n1+2, n2-1); // Sin el parentesis
			tmpTok.mat.destroy();
			tmpTok.mat.buildFrom(tmpTok.str.c_str());
			//tmpTok.mat.showStructure();
			tmpTok.mat.reemplazos=&reemplazos;
			tok.push_back(tmpTok);
			n1=i+1;
			n2=i+1;
			enExpr=0;
		}
		else if (enExpr>0 && expr[i]=='(') // Se abre un parentesis
			enExpr++;
		else if (enExpr>0 && expr[i]==')') // Se cierra un parentesis (no es el ultimo)
			enExpr--;
	}
	if (n1<l && n2<l && n1<=n2) // Guardar la parte final del string
	{
		tmpTok.esTexto=true;
		tmpTok.str.copySubString(expr, n1, n2);
		tok.push_back(tmpTok);
	}

	// Definir anchos minimos de tokens numericos "tok" a partir de los anchos contenidos en los rangos "rang"
	for (i=0; i<(int)tok.size(); i++)
	{
		tok[i].width=0; // Inicialmente el width minElement del token evaluado es zero
		for (n1=0; n1<(int)rang.size(); n1++)
		{
			if (rang[n1].width<=tok[i].width)
				continue;
			int j;
			for (j=0; j<(int)tok[i].str.size()-1; j++) // Ver si la variable rang[n1] esta contenida en el token
			{
				if (tok[i].esTexto == false && tok[i].str[j]==rang[n1].str[0])
				{
					tok[i].width=rang[n1].width; // Copiar el width minElement de la variable rang[n1] al token
					break;
				}
			}
		}
	}

	// Ahora se evaluate el comando en cada uno de los valores dentro de los rangos
	// Cuando la variable del primer rango llegue a su value maxElement se han evaluado todas las posibilidades
	while(rang[0].n<=rang[0].n2)
	{
		*expr=0;
		for (i=0; i<(int)tok.size(); i++) // Evaluacion de los tokens y copia a "expr"
		{
			if (tok[i].esTexto)
			{
				strcat(expr, tok[i].str.c_str()); // Copiar el texto
			}
			else
			{
				try
				{
					numLong=tok[i].mat.evaluateInteger(); // Evaluar el token i
				}
				catch (L_ArgException)
				{
					printf("Intento de usar una variable no definida\n");
					return 0;
				}
				tok[i].str=L_String(numLong, tok[i].width); // Transformar el resultado numerico en string
				strcat(expr, tok[i].str.c_str()); // Copiar el string
			}
		}
		if (debug)
			printf("%s\n",expr);
		else
			if( system(expr) ) {}  // Silenciar warning por no uso de value de retorno
		rang[rang.size()-1].n+=rang[rang.size()-1].dn;
		reemplazos[rang[rang.size()-1].idx].varRe+=rang[rang.size()-1].dn;
		i=rang.size()-1;
		while (i>0 && rang[i].n>rang[i].n2)
		{
			rang[i].n=rang[i].n1;
			reemplazos[rang[i].idx].varRe=rang[i].n1;
			rang[i-1].n+=rang[i-1].dn;
			reemplazos[rang[i-1].idx].varRe+=rang[i-1].dn;
		}
	}

	return 0;
}

int L_ProgramasFnesPato::main_L_encriptado(int argc, char *argv[], char *sal1, char *sal2)
{
	L_String comando, temp;
	char expr1[800]="";
	char expr2[800]="";
	bool echo=false;
	bool vez1=true;
	bool esHex=false;
	bool transfHex=false;
	int i;
	unsigned long semilla=514;


	if (argc<2)
	{
		printf("\n  L_Traductor comando\n\n"
			"Ejecuta \"comando\" con su value original. Usar $ para terminar traduccion.\n\n");
		return 1;
	}

	i=1;

	if (strlen(argv[i])>=2 && (argv[i][0]=='/' || argv[i][0]=='-') && (argv[i][1]=='e' || argv[i][1]=='E'))
	{
		echo=true;
		i++;
	}

	if (strlen(argv[i])>=2 && (argv[i][0]=='/' || argv[i][0]=='-') && (argv[i][1]=='h' || argv[i][1]=='H'))
	{
		transfHex=true;
		i++;
	}

	if (strlen(argv[i])>=2 && (argv[i][0]=='/' || argv[i][0]=='-') && (argv[i][1]=='e' || argv[i][1]=='E'))
	{
		echo=true;
		i++;
	}

	if (argv[i][0]=='~')
	{
		esHex=true;
		i++;
	}

	for ( ; i<argc ; i++)
	{
		if (argv[i][0]=='$')
		{
			if (argv[i][1]==0)
			{
				i++;
				break;
			}
			else if (argv[i][1]=='$')
			{
				strcat(expr1,"$");
				continue;
			}
		}
		if (vez1)
			vez1=false;
		else
			strcat(expr1, " ");
		if (strlen(expr1)+strlen(argv[i])+1>790)
		{
			printf("comando muy length\n");
			return 2;
		}
		strcat(expr1, argv[i]);
	}
	for ( ; i<argc ; i++)
	{
		if (vez1)
			vez1=false;
		else
			strcat(expr2, " ");
		if (strlen(expr2)+strlen(argv[i])+1>790)
		{
			printf("comando muy length\n");
			return 2;
		}
		strcat(expr2, argv[i]);
	}
	for (i=0; expr2[i]!=0; i++)
	    semilla+=expr2[i];

	if (esHex)
		comando.decodifBase16(expr1);
	else
		comando=expr1;
	comando.encript(semilla);
	if (transfHex)
	{
		temp.codifBase16(comando.c_str());
		temp.swap(comando);
	}
	if (sal1!=NULL && sal2!=NULL)
	{
		strcpy(sal1, comando.c_str());
		strcpy(sal2, expr2);
	}
	comando+=expr2;

	if (echo)
		puts(comando.c_str());
	else
		if( system(comando.c_str()) ) {} // Para silenciar warning por no usar value de retorno
	if (sal1!=NULL && sal2==NULL)
		strcpy(sal1, comando.c_str());
	return 0;
}

int L_ProgramasFnesPato::main_L_bmplz(int argc, char *argv[])
{
	L_ImageRGBUchar im;
	L_String nomBmp, nomBmz, nomTmp;
	L_FileName dneBmp, dneTmp;
	FILE *fp = NULL;
	double le[10];
	int hist[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	long i;
	bool argv1bmz = false;//, argv2bmz = false;

	if (argc == 1 || ((argv[1][0] == '/' || argv[1][0] == '-') && argv[1][1] == '?'))
	{
		printf("L_bmz im.bmp im.bmplz       (compresion)\n");
		printf("L_bmz im.bmplz              (info sobre el archivo)\n");
		printf("L_bmz im.bmplz im.bmp       (descompresion)\n");
		printf("\n");
		printf("Los archivos descomprimidos se llamaran im1.bmp ... imN.bmp\n");
		return 1;
	}
	nomTmp = argv[1];
	nomTmp.lowercase();

	argv1bmz = nomTmp.size()-1 > 6 && nomTmp[nomTmp.size()-1-6] == '.' && nomTmp[nomTmp.size()-1-5] == 'b' && nomTmp[nomTmp.size()-1-4] == 'm' && nomTmp[nomTmp.size()-1-3] == 'p'&& nomTmp[nomTmp.size()-1-2] == 'l'&& nomTmp[nomTmp.size()-1-1] == 'z';

	if (argv1bmz && argc == 2)
	{
		// Info sobre el archivo
		nomBmz = argv[1];
		fp = fopen(nomBmz.c_str(), "rb");
		if (fp == NULL)
			{printf("No se pudo read %s\n", nomBmz.c_str()); return 1;}
		printf("%d imagenes\n", L_ImageRGBUchar::leeComprRGB_contar(fp, hist));
		printf("Metodos: ");
		for (i=0; i<10; i++)
			printf("%d ", hist[i]);
		printf("\n");
		fclose(fp);

	}
	else if (argv1bmz && argc == 3)
	{
		// Descomprimir
		nomBmz = argv[1];
		nomBmp = argv[2];
		dneBmp.asign(nomBmp.c_str());
		dneTmp.asign(nomBmp.c_str());

		fp = fopen(nomBmz.c_str(), "rb");
		if (fp == NULL)
			{printf("No se pudo read %s\n", nomBmz.c_str()); return 1;}
		i = 0;

		long lo1 = ftell(fp);
		fseek(fp, 0, SEEK_END);
		long lo2 = ftell(fp);
		fseek(fp, lo1, SEEK_SET);
		long lo3;

		// Ver cuantos ceros hay que agregar a los nombres de archivo
		int nFr = im.leeComprRGB_contar(fp);
		int nFr2 = nFr, nceros = 0;
		while (nFr2 > 0)
		{
			nFr = nFr / 10;
			nceros++;
		}

		while (! feof(fp))
		{
			lo3 = ftell(fp);
			if (lo3 == lo2) // A veces feof no funka...
				break;
			if (im.readComprRGB(fp) == 0)
				{printf("Error en el formato de %s\n", nomBmz.c_str()); fclose(fp); return 2;}
			dneTmp.name = dneBmp.name + L_String(i,nceros);
			if (im.saveImage(dneTmp.c_str()) == 0)
				{printf("No se pudo write %s\n", dneTmp.c_str()); fclose(fp); return 3;}
			i++;
		}
		fclose(fp);
	}
	else
	{
		// Comprimir
		nomBmp = argv[1];
		nomBmz = argv[2];
		if (im.readImage(nomBmp.c_str()) == 0)
			{printf("Error al read imagen %s\n", argv[1]); return 2; }
		fp = fopen(nomBmz.c_str(), "wb");
		if (fp == NULL)
			{printf("Error al write imagen %s\n", argv[2]); return 3;}
		if (im.writeComprRGB(fp, le) == 0)
			{printf("Error interno\n"); fclose(fp); return 4;}
		printf("or:%ld  hu:%ld  d1:%ld\n", (long)le[0], (long)le[1], (long)le[2]);
		fclose(fp);
	}
	return 0;
}

inline double L_multttt(double a, double b) {return a*b;}
inline double L_multtttp(double *a, double *b) {return (*a)*(*b);}
extern "C" inline double L_multtttpC(double *a, double *b) {return (*a)*(*b);}
inline double L_multtttpNT(double *a, double *b) nothrows() {return (*a)*(*b);}


static FILE *fp_L_PFPato_cBmz = NULL;

void L_ProgramasFnesPato_main_capturaBmz_atexit()
{
	if (fp_L_PFPato_cBmz != NULL)
		fclose(fp_L_PFPato_cBmz); // Para que grabe la imagen
	fp_L_PFPato_cBmz = NULL;
	printf("Archivo cerrado\n");
}

int L_ProgramasFnesPato::main_capturaBmplz(int argc, char *argv[])
{
	L_ImageRGBUchar im;
	double t1, t2;
	int modo;
	long nframe;

	// Parametros de grabacion
	if (argc < 2 || argc > 4 || ((argv[1][0] == '/' || argv[1][0] == '-') && argv[1][1] == '?'))
	{
		printf("Uso:\n");
		printf("   L_Capturabmz arch.bmz [modo] [v]\n\n");
		printf("modo = -1 (auto), 0(no compr), 1(huff), 2(D-huff), 3(DD-huff), 4(DDRLE)\n");
		printf("v = mostrar imagen en ventana (puede afectar algo la velocidad)\n");
		return 1;
	}

	L_CapturadorImagen cap;

	if (cap.capturarImagen() == false)
	{
		printf("No se encuentra la camara\n");
		exit(0);
	}
	if (cap.estereo.esEstereo == true)
	{
		cap.estereo.usarIzq0Der1Ambas2 = 2;
	}

	fp_L_PFPato_cBmz = fopen(argv[1], "wb");
	if (fp_L_PFPato_cBmz == NULL)
	{
		printf("No se pudo write %s\n", argv[1]);
		exit(0);
	}

	if (argc >= 3 && argv[2][0] >= '0' && argv[2][0] <= '9')
		modo = argv[2][0] - '0';
	else
		modo = -1;
	if (modo < 0 || modo > 4)
		modo = -1;
	// Grabacion
	atexit(&L_ProgramasFnesPato_main_capturaBmz_atexit);

#if defined (__L_VENTANA__)
	bool mostrarVentana = false;
	if (argc >= 4 && (argv[3][0] == 'v' || argv[3][0] == 'V'))
		mostrarVentana = true;
	if (cap.capturarImagen() == false)
		{printf("No se encuentra la camara\n");exit(0);}
	L_VentanaImagen vent(0,0,cap.im.lx+cap.imDer.lx, cap.im.ly, "captura (click para comenzar)");
	t1 = L_TIME();
	nframe = 0;
	// Esperar hasta que se haga click en la imagen
	printf("Haga click en la ventana para comenzar a write\n");
	while(vent.mouseE() == L_MouseLeido)
	{
		if (cap.capturarImagen() == false)
		{
			printf("No se encuentra la camara\n");
			exit(0);
		}
		if (cap.estereo.esEstereo == false)
		{
			if (mostrarVentana)
			{
				vent.dibuja(cap.im);
				vent.check();
			}
		}
		else
		{
			im.concatenateHorz(cap.im, cap.imDer);
			if (mostrarVentana)
			{
				vent.dibuja(im);
				vent.check();
			}
		}
	}
	vent.mouseE() = L_MouseLeido;
#else
	printf("Presione ENTER para comenzar a write\n");
	getchar();
#endif
	printf("Grabando %s(modo %d). Presione ctrl-c para terminar\n", argv[1], modo);
	t1 = L_TIME();
	nframe = 0;
	while(true)
	{
		if (cap.capturarImagen() == false)
		{
			printf("No se encuentra la camara\n");
			exit(0);
		}
		t2 = L_TIME();
		if (cap.estereo.esEstereo == false)
		{
			cap.im.writeComprRGB(fp_L_PFPato_cBmz, NULL, modo, t2, 1); // Forzar metodo Huffman-resta, es mas rapido asi
			fflush(fp_L_PFPato_cBmz);
#if defined (__L_VENTANA__)
			if (mostrarVentana)
			{
				vent.dibuja(cap.im);
				vent.check();
			}
#endif
		}
		else
		{
			im.concatenateHorz(cap.im, cap.imDer);
			im.writeComprRGB(fp_L_PFPato_cBmz, NULL, modo, t2, 2);
			fflush(fp_L_PFPato_cBmz);
#if defined (__L_VENTANA__)
			if (mostrarVentana)
			{
				vent.dibuja(im);
				vent.check();
			}
#endif
		}
		nframe++;
		if (nframe%30 == 0) // 5 seg
		{
			printf("%d fps    \r", (int)(30 / (t2-t1)));
			t1 = t2;
		}
#if defined (__L_VENTANA__)
		if (nframe < 20)
		{
			vent.mouseE() = L_MouseLeido;
		}
		else
		{
			if (vent.mouseE() != L_MouseLeido)
				break;
		}
#endif // __COMPAT_FLTK__
	}
	fclose(fp_L_PFPato_cBmz);
	fp_L_PFPato_cBmz = NULL;

	printf("\n");
	return 0;
}

int L_ProgramasFnesPato::main_pruebaClases(int argc, char *argv[])
{
	int nEl = 20;
	L_Array<L_Quaternion> arr(nEl), arr2(nEl);
	L_Array<bool> elim(nEl);
	L_Node<L_Quaternion> *nodo;
	L_List<L_Quaternion> lis;
	int i;

	printf("Probando eliminacion en arreglos y listas con L_Array<bool>\n");

	for (i=0; i<nEl; i++)
	{
		arr[i].fijaAzar();
		arr[i].b = i;
	}

	for (i=0; i<nEl; i++)
		arr[i].print("q");
	printf("\n");

	for (i=0; i<nEl; i++)
		elim[i] = arr[i].a < 0;
	elim.resize(i);

	arr2 = arr;
	arr2.erase_preserving_order(elim);

	for (i=0; i<(int)arr2.size(); i++)
		arr2[i].print("q");
	printf("\n");

	for (i=0; i<(int)arr.size(); i++)
		lis.push_back(arr[i]);
	lis.erase_preserving_order(elim);

	for (nodo=lis.root; nodo!=NULL; nodo=nodo->sig)
		nodo->c.print("q");
	
	printf("\nPresione ENTER\n\n");
	getchar();
	return 0;
}

int L_ProgramasFnesPato::main_imprimirArgv(int argc, char *argv[])
{
	int i;
	for (i=0; i<argc; i++)
		printf("%s ", argv[i]);
	printf("\n");
	printf("presione ENTER para continuar\n");
	getchar();
	return 0;
}


int L_ProgramasFnesPato::main_interpolarVideo(int argc, char *argv[])
{
	L_ImageRGBUchar im;
	L_ImageGrayDouble imBN;
	FILE *fpIn = NULL, *fpOut = NULL;
	L_ImageRGBUchar imRGB;

	if (argc!=3)
	{
		printf("main_interpolarVideo vidEntrada.bmz vidSalida.bmz\n");
		printf("\n");
		return 1;
	}

	fpIn = fopen(argv[1], "rb");
	if (fpIn == NULL)
	{
		printf("No se puede read archivo %s\n", argv[1]);
		return 1;
	}
	fpOut = fopen(argv[2], "rb");
	if (fpIn == NULL)
	{
		fclose(fpIn);
		printf("No se puede write archivo %s\n", argv[2]);
		return 1;
	}

	throw L_NonImplementedException();

	fclose(fpIn);
	fclose(fpOut);
	return 0;
}

struct datosHomogr
{
	double xi[100];
	double xf[100];
	double yi[100];
	double yf[100];
	int np;
};

void errHomogr(const void *obj, const L_Matrix &h, L_Matrix &e)
{
	datosHomogr *dat = (datosHomogr *)obj;
	double xh, yh, zh;
	int i;
	e.reallocate(2*dat->np, 1);
	for (i=0; i<dat->np; i++)
	{
		xh = h(0,0)*dat->xi[i] + h(0,1)*dat->yi[i] + h(0,2);
		yh = h(0,3)*dat->xi[i] + h(0,4)*dat->yi[i] + h(0,5);
		zh = h(0,6)*dat->xi[i] + h(0,7)*dat->yi[i] + 1;
		e(2*i+0,0) = dat->xf[i] - xh/zh;
		e(2*i+1,0) = dat->yf[i] - yh/zh;
	}
}


int L_ProgramasFnesPato::main_rectificar(int argc, char *argv[])
{
	L_ImageRGBUchar im, im2;
	datosHomogr dat;
	L_Matrix h(8,1);
	int i;
	//double xi[4] = {84, 225, 226,  80}, xf[4]={0, 147, 147,   0};
	//double yi[4] = {39, 42,  209, 203}, yf[4]={0,   0, 167, 167};
	//double xi[4] = {105, 191, 184, 102}, xf[4]={0, 85,  85,   0};
	//double yi[4] = {40,  38,  203, 209}, yf[4]={0,  0, 167, 167};
	//double xi[4] = {113, 240, 238, 105}, xf[4]={0, 147, 147,   0};
	//double yi[4] = { 24,  22, 167, 169}, yf[4]={0,  0, 167, 167};
	//double xi[4] = {147, 237, 238, 150}, xf[4]={0, 85,  85,   0};
	//double yi[4] = {48,  48,  217, 217}, yf[4]={0,  0, 167, 167};

	//double xi[4] = {69, 238, 226,  80}, xf[4]={0, 147, 147,   0};
	//double yi[4] = {45, 46,  209, 203}, yf[4]={0,   0, 167, 167};
	//double xi[4] = {95, 191, 184, 96}, xf[4]={0, 85,  85,   0};
	//double yi[4] = {40,  38,  203, 205}, yf[4]={0,  0, 167, 167};
	//double xi[4] = {113, 240, 238, 105}, xf[4]={0, 147, 147,   0};
	//double yi[4] = { 24,  22, 167, 169}, yf[4]={0,  0, 167, 167};
	double xi[4] = {147, 251, 238, 150}, xf[4]={0, 85,  85,   0};
	double yi[4] = {48,  48,  217, 217}, yf[4]={0,  0, 167, 167};


	//im.readImage("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\screenshot1.png");
	//im.readImage("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\screenshot30.png");
	//im.readImage("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\screenshot60.png");
	im.readImage("C:\\Documents and Settings\\Administrador\\Mis documentos\\Mis vídeos\\Brasil\\Tracking3D\\screenshot90.png");

	L_LevenbergMarquardt opt;

	for (i=0; i<4; i++)
		dat.xi[i] = xi[i];
	for (i=0; i<4; i++)
		dat.yi[i] = yi[i];
	for (i=0; i<4; i++)
		dat.xf[i] = xf[i];
	for (i=0; i<4; i++)
		dat.yf[i] = yf[i];
	dat.np = 4;

	h.setZero();
	h(0,0) = 1;
	h(4,0) = 1;

	opt.nIterationsMax = 100;
	opt.lenVector = 8;
	opt.minimize_vect((const void *)&dat, h, &errHomogr);

	printf("Parametros de la homografia: ");
	for (i=0; i<8; i++)
		printf("%f ", h(i,0));
	printf("\n");

	im2.applyHomography(h, im);
	im2.saveImage("sale.png");
	return 0;	
}


/*
int L_ProgramasFnesPato::main_muestraBmz(int argc, char *argv[])
{
	if (argc != 2 || argv[1][0] == '/' && argv[1][1] == '?')
	{
		printf("\nL_muestrabmz archivo.bmz\n\n");
		return 1;
	}

	L_ImageRGBUchar im;
	FILE *fp;
	double t;
	long ncam;

	fp = fopen(argv[1],"rb");
	if (fp == NULL)
	{
		printf("No se puede abrir el archivo %s\n", argv[1]);
		return 2;
	}

	L_VentanaImagen vent(0,0,1200,800, "video");

	while( !feof(fp) && im.readComprRGB(fp, &t, &ncam) != -1 )
	{
		vent.dibujaRedimensiona(im);
		vent.check();
	}
	fclose(fp);
	return 0;
}
*/

#if defined (__COMPAT_FLTK__)
int L_ProgramasFnesPato::main_muestra_plpts(int argc, char *argv[])
{
	L_Pose3DSpline spl;
	L_ImageRGBUchar im;
	int lx = Fl::w(), ly = Fl::h();
	L_VentanaImagen vent(0,0,lx,ly, "Puntos control (click para salir)");

	if (argc != 2 || argv[1][0] == '/' && argv[1][1] == '?')
	{
		printf("main_muestra_plpts archivo.bmz\n\n");
		return 1;
	}

	im.reallocate(lx, ly);

	spl.leerPuntosControl(argv[1]);
	spl.genDibujoFlaite(im);

	vent.dibujaRedimensiona(im);
	while( vent.mouseE() == L_MouseLeido )
	{
		vent.check();
	}
	return 0;
}
#endif // __COMPAT_FLTK__


#ifdef __COMPAT_FLTK__
int L_ProgramasFnesPato::main_dib3d(int argc, char *argv[])
{
	L_VentMapa3D v3d(0,0,640,480, "Puntos");
	L_VentanaImagenFLTK vPanTilt(0, 500, 200, 200, "pan-tilt");
	L_VentanaImagenFLTK vXY(250, 500, 200, 200, "x-y");
	L_VentanaImagenFLTK vZ(500, 500, 200, 200, "z");
	L_VentanaImagenFLTK vSal(750, 500, 200, 200, "Click para salir");
	L_Pose3D poseCam;

	v3d.show();

	if (argc != 2 || argv[0][0] == '/' || argv[0][0] == '?')
	{
		printf(
			"dib3d arch.txt\n"
			"\n"
			"arch.txt:  puntos del tipo: x y z r g b tipo radio dx dy dz\n"
			"           solo es obligatorio especificar x y z\n"
			"           tipo: 0=punto  1=cubo  2=esfera\n"
			"           dx dy dz: para dibujar flecha \n");
		return 1;
	}

	std::vector<L_CoordsCart3D> pun(20);
	std::vector<long> color(20);
	
	for (unsigned int i=0; i<(int)pun.size(); i++)
	{
		pun[i].x = L_RANDOM(-10, 10);
		pun[i].y = L_RANDOM(-10, 10);
		pun[i].z = L_RANDOM(-10, 10);
		color[i] = rand() % 0x00FFFFFF;
	}

	poseCam.fijaCero();
	poseCam.pos.x = -50;
	while (vSal.mouseE() == L_MouseLeido)
	{
		v3d.clear();
		for (unsigned int i=0; i<(int)pun.size(); i++)
			v3d.agregaPunto3D_cubo(pun[i], color[i]);
		if (vPanTilt.mouseE() != L_MouseLeido)
		{
			poseCam.ori.pan = ((vPanTilt.mouseX() / (double)vPanTilt.vent.x()) - 0.5) * M_PI;
			poseCam.ori.tilt = ((vPanTilt.mouseY() / (double)vPanTilt.vent.y()) - 0.5) * M_PI;
			vPanTilt.mouseE() = L_MouseLeido;
		}
		v3d.fijaPoseCamara(poseCam);
		v3d.redraw();
		Fl::check();
	}
	return 0;
}
#endif // __COMPAT_FLTK__


#ifdef __COMPAT_IPLIMAGE__
int L_ProgramasFnesPato::main_L_GrabaFramesOpenCV(int argc, char *argv[])
{
	CvCapture *capture;
	IplImage *frame;
	int keypr;
	int numIm=0;
	char carpeta[100]="", nomarch[150]="";
	int n1=0, n2=0, n3=0, n4=0, n5=0;

	capture = cvCaptureFromCAM( CV_CAP_ANY ); //110+i );
	if( !capture )
	{
		fprintf( stderr, "ERROR: capture is NULL \n");
		getchar();
		return -1;
	}

	// Create a window in which the captured images will be presented
	cvNamedWindow( "Camara 1", CV_WINDOW_AUTOSIZE );

	printf("name de la carpeta:");
	L_gets(carpeta,99);
	sprintf(nomarch, "if not exist %s md %s\n", carpeta, carpeta);
	if( system(nomarch) ) {} // Silenciar warning por no uso de value de retorno

	printf("Botones: (espacio)=write imagen, (esc)=salir\n");
	printf("Seleccione la ventana con el maouse para comenzar\n");

	// Show the image captured from the camera in the window and repeat
	while( true )
	{
		if (!cvGrabFrame( capture ))
		{
			printf("No se pudo capturar\n");
			getchar();
			return -1;
		}
		frame = cvRetrieveFrame(capture);

		if (frame==NULL)
		{
			printf("No se pudo capturar\n");
			getchar();
			return -1;
		}

		cvShowImage( "Camara 1", frame );
		// Do not release the frame!

		//If ESC key pressed, Key=0x10001B under OpenCV 0.9.7(linux version),
		//remove higher bits using AND operator
		keypr = cvWaitKey(10);
		if ( (keypr & 255) == ' ' )
		{
			// Grabar imagen
			n5=numIm%10;
			n4=(numIm/10)%10;
			n3=(numIm/100)%10;
			n2=(numIm/1000)%10;
			n1=numIm/10000;

			sprintf(nomarch, "%s/im%d%d%d%d%d.bmp", carpeta, n1, n2, n3, n4, n5);
			cvSaveImage(nomarch, frame);

			printf("im: %d\n", numIm);

			numIm++;
		}
		if( (keypr & 255) == 27 )
			break;
	}

	// Release the capture device housekeeping
	cvReleaseCapture( &capture );
	cvDestroyWindow( "Camara 1" );
	return 0;
}
#endif

#ifdef __COMPAT_IPLIMAGE__
int L_ProgramasFnesPato::main_L_bmplzavi(int argc, char *argv[])
{
	L_ImageRGBUchar im;
	L_String nomAvi, nomBmz, nomTmp;
	IplImage *imCV = NULL;
	FILE *fp = NULL;
	long pos, posFin, nim;
	bool argv1bmz = false;
	CvCapture *captCV = NULL;
	CvVideoWriter *wrCV = NULL;

	if (argc==1 || argc == 2)
	{
		printf("L_bmplzavi vid.avi vid.bmplz  (avi a bmplz)\n");
		printf("L_bmplzavi vid.bmplz vid.avi CODEC(bmplz a avi)\n");
		printf("\n");
		printf("codec = PIM1, MJPG, MP42, DIV3, DIVX, U263, I263, FLV1, ... (4 letras)\n");
		printf("por defecto codec = DIB\n");
		return 1;
	}
	nomTmp = argv[1];
	nomTmp.lowercase();

	argv1bmz = nomTmp.size()-1 > 6 && nomTmp[nomTmp.size()-1-6] == '.' && nomTmp[nomTmp.size()-1-5] == 'b' && nomTmp[nomTmp.size()-1-4] == 'm' && nomTmp[nomTmp.size()-1-3] == 'p' && nomTmp[nomTmp.size()-1-2] == 'l' && nomTmp[nomTmp.size()-1-1] == 'z';

	if (argv1bmz)
	{
		// bmz a avi
		nomBmz = argv[1];
		nomAvi = argv[2];
		fp = fopen(nomBmz.c_str(), "rb");
		if (fp == NULL)
			{printf("No se pudo read %s\n", nomBmz.c_str()); return 1;}
		nim = L_ImageRGBUchar::leeComprRGB_contar(fp);
		if (nim == -1)
			{printf("Error al read %s\n", nomBmz.c_str()); return 2;}
		printf("%ld imagenes\n", nim);
		pos = ftell(fp);
		if (im.readComprRGB(fp) == -1)
			{printf("Error al read %s\n", nomBmz.c_str()); return 3;}
		fseek(fp, 0, SEEK_END);
		posFin = ftell(fp);
		fseek(fp, pos, SEEK_SET);
		if (argc == 4)
			wrCV = cvCreateVideoWriter(nomAvi.c_str(), CV_FOURCC(argv[3][0],argv[3][1],argv[3][2],argv[3][3]), 30, cvSize(im.lx, im.ly), 1);
		else
			wrCV = cvCreateVideoWriter(nomAvi.c_str(), CV_FOURCC('D','I','B',' '), 30, cvSize(im.lx, im.ly), 1);
		if (wrCV == NULL)
			{printf("OpenCV no pudo crear el archivo %s %c%c%c%c\n", nomAvi.c_str(), argv[3][0],argv[3][1],argv[3][2],argv[3][3]); return 4;}
		while (!feof(fp))
		{
			if (ftell(fp) == posFin)
				break;
			if (im.readComprRGB(fp) == -1)
				{printf("Error al read %s\n", nomBmz.c_str()); return 5;}
			im.copyTo(&imCV);
			cvWriteFrame(wrCV, imCV);
		}
		cvReleaseImage(&imCV);
		cvReleaseVideoWriter(&wrCV);
		fclose(fp);
	}
	else
	{
		// avi a bmz
		nomAvi = argv[1];
		nomBmz = argv[2];
		int nim = 0;
		fp = fopen(nomBmz.c_str(), "wb");
		if (fp == NULL)
			{printf("No se pudo read %s\n", nomBmz.c_str()); return 1;}
		captCV = cvCaptureFromAVI(nomAvi.c_str());
		if (captCV == NULL)
			{printf("No se pudo abrir %s\n", nomAvi.c_str()); return 2;}
		while (cvGrabFrame(captCV))
		{
			imCV = cvRetrieveFrame(captCV);
			im = *imCV;
			im.writeComprRGB(fp);
			nim ++;
		}
		cvReleaseCapture(&captCV);
		fclose(fp);
		printf("%d imagenes\n", nim);
	}
	return 0;

}
#endif //__COMPAT_IPLIMAGE__

#if defined(__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_L_capturaCV(int argc, char *argv[])
{
	L_ImageRGBUchar im;
	IplImage *imCV = NULL;
	FILE *fp = NULL;
	CvCapture *captCV = NULL;
	CvVideoWriter *wrCV = NULL;
	// Parametros de grabacion
	if (argc < 2 || argc > 3 || ((argv[1][0] == '/' || argv[1][0] == '-') && argv[1][1] == '?'))
	{
		printf("Uso:\n");
		printf("   L_CapturaCV arch.avi [DIB |MJPG|...]\n\n");
		return 1;
	}

	L_CapturadorImagen cap;

	if (cap.capturarImagen() == false)
	{
		printf("No se encuentra la camara\n");
		exit(0);
	}
	if (cap.estereo.esEstereo == true)
	{
		cap.estereo.usarIzq0Der1Ambas2 = 2;
	}

	if (cap.capturarImagen() == false)
		{printf("No se encuentra la camara\n");exit(0);}

	L_VentanaImagenCV vent(0,0,cap.im.lx+cap.imDer.lx, cap.im.ly, "captura (click para comenzar)");
	// Esperar hasta que se haga click en la imagen
	printf("Haga click en la ventana para comenzar a write\n");
	while(vent.mouseE() != L_MouseSuelta)
	{
		vent.mouseE() = L_MouseLeido;
		if (cap.capturarImagen() == false)
		{
			printf("No se encuentra la camara\n");
			exit(0);
		}
		if (cap.estereo.esEstereo == false)
			im = cap.im;
		else
			im.concatenateHorz(cap.im, cap.imDer);
		vent.dibuja(im);
		vent.check();
	}
	vent.mouseE() = L_MouseLeido;

	vent.check();

	if (vent.mouseE() != L_MouseLeido)
		vent.mouseE() = L_MouseLeido;

	printf("Grabando %s. Click para terminar\n", argv[1]);

	if (argc == 3)
	{
		wrCV = cvCreateVideoWriter(argv[1], CV_FOURCC(argv[2][0],argv[2][1],argv[2][2],argv[2][3]), 30, cvSize(im.lx, im.ly), 1);
		if (wrCV == NULL)
			{printf("OpenCV no pudo crear el archivo %s %c%c%c%c\n", argv[1], argv[2][0],argv[2][1],argv[2][2],argv[2][3]); return 4;}
	}
	else
	{
		wrCV = cvCreateVideoWriter(argv[1], CV_FOURCC('D','I','B',' '), 30, cvSize(im.lx, im.ly), 1);
		if (wrCV == NULL)
			{printf("OpenCV no pudo crear el archivo %s DIB\n", argv[2]); return 5;}
	}

	while (vent.mouseE() == L_MouseLeido)
	{
		if (cap.capturarImagen() == false)
		{
			printf("No se encuentra la camara\n");
			exit(0);
		}
		if (cap.estereo.esEstereo == false)
			im = cap.im;
		else
			im.concatenateHorz(cap.im, cap.imDer);

		vent.dibuja(im);
		vent.check();
		im.copyTo(&imCV);
		cvWriteFrame(wrCV, imCV);
	}

	cvReleaseImage(&imCV);
	cvReleaseVideoWriter(&wrCV);

	printf("\n");
	return 0;
}
#endif //__COMPAT_IPLIMAGE__

int L_ProgramasFnesPato::main_capturaPNG(int argc, char *argv[])
{
	L_CapturadorImagen capt;
	int i, nframes;
	char prefix[100], nom_im[100];
	if (capt.capturarImagen() == false)
	{
		printf("No se puede capturar. Presione ENTER para salir\n");
		return 1;
	}
	
	printf("Prefijo del archivo de imagenes a write (sin extension): ");
	L_gets(prefix, 100);

	printf("Numero de imagenes a write: ");
	while (scanf("%d", &nframes) != 1)
		printf("** Ingrese solamente el numero seguido de ENTER por favor... **\n");

	// Configurar captura
	if (capt.estereo.esEstereo)
		capt.estereo.usarIzq0Der1Ambas2 = 2;

	for (i=0; i<nframes; i++)
	{
		if (capt.capturarImagen() == false)
			break;
		L_snprintf(nom_im, 99, "%s_im_%04d.png", prefix, i+1);
		capt.im.saveImage(nom_im);
		if (capt.imDer.data() != NULL)
		{
			L_snprintf(nom_im, 99, "%s_der_%04d.png", prefix, i+1);
			capt.im.saveImage(nom_im);
		}
		if (capt.imKin.data() != NULL)
		{
			L_snprintf(nom_im, 99, "%s_prof_%04d.txt", prefix, i+1); // png
			capt.imKin.saveImage_16_txt(nom_im);
			L_snprintf(nom_im, 99, "%s_prof_%04d.png", prefix, i+1); // png
			capt.imKin.saveImage_8(nom_im, 255.0 / 4096.0); 
		}
	}
	return 0;
}




#if defined (__COMPAT_FLTK__)
int L_ProgramasFnesPato::main_pruebaFourier(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imx, imag, immod;
	L_ImageRGBUchar imFou;
	int i, j, lxlog2, lylog2, lx2, ly2, lx2min, ly2min;
	double med;
	double dt;
	double t[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	long frame = 0;

	//cap.fijarResolucion(320,200);  ... :(
	if (cap.capturarImagen()==false) {printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;	}

	lxlog2 = (int)(log((double)cap.im.lx-1)/L_LOG_2);
	lylog2 = (int)(log((double)cap.im.ly-1)/L_LOG_2);

	lx2 = (int)pow(2.0, lxlog2);
	ly2 = (int)pow(2.0, lylog2);
	imx.reallocate(lx2,ly2);

	imx.setZero();

	L_VentanaImagen vent1(0, 0, imx.lx, imx.ly,"Click para salir");
	L_VentanaImagen vent2(imx.lx + 25, 0, imx.lx, imx.ly,"Fourier");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		if (vent1.mouseE() != L_MouseLeido)
			break;

		frame++;

		lx2min = L_MIN(cap.im.lx, lx2);
		ly2min = L_MIN(cap.im.ly, ly2);

		med = 0;
		for (j=0; j<ly2min; j++)
			for (i=0; i<lx2min; i++)
				med += ((double)cap.im.pix(i,j,0) + (double)cap.im.pix(i,j,1) + (double)cap.im.pix(i,j,2)) / 3.0;
		med = med / lx2min / ly2min;

		imx.reallocate(lx2,ly2);

		imx.setConstant(med);

		for (j=0; j<ly2min; j++)
			for (i=0; i<lx2min; i++)
				imx.pix(i,j) = ((double)cap.im.pix(i,j,0) + (double)cap.im.pix(i,j,1) + (double)cap.im.pix(i,j,2)) / 3.0 - med;

		imag.reallocate(imx.lx, imx.ly);
		imag.setZero();

		imx.fourier(1, imag); // 165 [ms)

		immod.reallocate(imx.lx, imx.ly);
		for (j=0; j<immod.ly; j++)
			for (i=0; i<immod.lx; i++)
				immod.pix(i,j) = log(sqrt(imx.pix(i,j)*imx.pix(i,j) + imag.pix(i,j)*imag.pix(i,j)) + 1e-30);

		immod.normalizeHistogram(5);
		immod.fftshift();

		imFou = immod;

		t[frame%10] = L_TIME();
		dt = (t[frame%10] - t[(frame+1)%10]) / 10;
		printf("%f[ms]\n", dt*1000);

		vent1.dibuja(cap.im);
		vent2.dibuja(imFou);
		vent1.check();
	}

	return 0;
}
#endif

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_pruebaHough(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imHough, imx;
	L_ImageGrayUchar binariz;
	L_ImageRGBUchar imRGB;
	std::vector<double> px, py, ang;
	L_ShapeArray lins;
	double dt;
	double t[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	long frame = 0;
	double umbralBin = 0.1, umbralHough = 105;
	int nRho=150, nTheta=90;
	bool usaGrad = true, supresionNoMaximo = true;

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	L_VentanaImagen vent1(0, 0, cap.im.lx, cap.im.ly, "Click para salir");
	L_VentanaImagen vent2(0, cap.im.ly + 100, cap.im.lx, cap.im.ly, "Lineas");
	L_VentanaImagen vent3(cap.im.lx + 100, cap.im.ly + 100, nTheta, nRho, "Hough");
	L_VentanaImagen vent4(cap.im.lx + 100, 0, cap.im.lx, cap.im.ly, "Grad");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		vent1.dibuja(cap.im);

		imx = cap.im;
		imHough.houghDetPuntosBordes(imx, binariz, px, py, ang, umbralBin, supresionNoMaximo);
		imHough.houghDetLineas(imx, px, py, ang, nTheta, nRho, usaGrad); // Aca la imagen se transformUsing en espacio de Hough
		imHough.houghDibLineas(imx, lins, umbralHough);

		imRGB = cap.im;
		imRGB.genDrawing(lins);
		vent2.dibuja(imRGB);

		imHough.multiply_to_each_element(0.05);
		imRGB = imHough;
		vent3.dibuja(imRGB);

		imRGB = binariz;
		vent4.dibuja(imRGB);

		if (vent1.mouseE() != L_MouseLeido)
			break;
		
		if (vent4.mouseE() == L_MouseArrastra)
		{
			umbralBin = vent4.mouseX() *1.0 / cap.im.lx;
			umbralHough = vent4.mouseY() *250.0 / cap.im.ly;
		}
		vent4.mouseE() = L_MouseLeido;

		if (vent1.mouseE() == L_MouseSuelta)
			supresionNoMaximo = !supresionNoMaximo;
		vent1.mouseE() = L_MouseLeido;

		if (vent3.mouseE() == L_MouseSuelta)
			usaGrad = !usaGrad;
		vent3.mouseE() = L_MouseLeido;

		t[frame%10] = L_TIME();
		dt = (t[frame%10] - t[(frame+1)%10]) / 10;
		printf("%f[ms]  umbralBin:%f  umbralHough:%f\n", dt*1000, umbralBin, umbralHough);

		vent1.check();
	}

	return 0;
}
#endif

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_pruebaGabor(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imx, imGabRe, imGabIm, imGabAmpl;
	L_ImageRGBUchar imGabReIm;
	int i, j;
	double vald;
	double s, lambda = 5.0, theta = 0;
	double dt;
	double t[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double dx, dy;
	long frame = 0;

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	imx.reallocate(cap.im.lx, cap.im.ly);
	imx.setZero();
	imGabReIm.reallocate(cap.im.lx, cap.im.ly);

	L_VentanaImagen vent1(0, 0, imx.lx, imx.ly,"Click para cambiar parametros");
	L_VentanaImagen vent2(imx.lx + 25, 0, imx.lx, imx.ly,"Filtrado Gabor - click para salir");
	L_VentanaImagen vent3(0, imx.ly+100, imx.lx, imx.ly,"Amplitud Gabor");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		imx = cap.im;

		s = 1.6*lambda;
		imGabRe.filterGabor(imx, &imGabIm, s, lambda, theta);
		imGabReIm.reallocate(imGabRe.lx, imGabRe.ly);

		for (j=0; j<imGabReIm.ly; j++)
		{
			for (i=0; i<imGabReIm.lx; i++)
			{
				imGabReIm.pix(i,j,0) = 0;
				vald = 127 + 128*imGabRe.pix(i,j)/2 * 20;
				if (vald>255)
					vald=255;
				if (vald<0)
					vald = 0;
				imGabReIm.pix(i,j,1) = (L_uchar)vald;
				vald = 127 + 128*imGabIm.pix(i,j)/2;
				if (vald>255)
					vald=255;
				if (vald<0)
					vald = 0;
				imGabReIm.pix(i,j,2) = (L_uchar)vald;
			}
		}

		imGabAmpl.reallocate(imGabRe.lx, imGabRe.ly);

		for (j=0; j<imGabAmpl.ly; j++)
			for (i=0; i<imGabAmpl.lx; i++)
				imGabAmpl.pix(i,j) = 20 * sqrt(imGabRe.pix(i,j)*imGabRe.pix(i,j) + imGabIm.pix(i,j)*imGabIm.pix(i,j));

		if (vent1.mouseE() == L_MouseArrastra)
		{
			dx = vent1.mouseX() - imGabReIm.lx/2;
			dy = vent1.mouseY() - imGabReIm.ly/2;
			lambda = sqrt(dx*dx+dy*dy) / 10;
			theta = -atan2(dy, dx);
		}
		vent1.mouseE() = L_MouseLeido;

		if (vent2.mouseE() != L_MouseLeido)
			break;

		t[frame%10] = L_TIME();
		dt = (t[frame%10] - t[(frame+1)%10]) / 10;
		printf("%f[ms]\n", dt*1000);

		vent1.dibuja(cap.im);
		vent2.dibuja(imGabReIm);
		imGabReIm = imGabAmpl;
		vent3.dibuja(imGabReIm);
		vent1.check();
	}

	return 0;
}
#endif

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_pruebaGaborUmbral(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imx, imGabRe, imGabIm, imGabAmpl;
	L_ImageRGBUchar imGabReIm;
	int i, j;
	double vald;
	double s, lambda = 5.0, theta = 0;
	double dt;
	double t[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double dx, dy;
	long frame = 0;

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	imx.reallocate(cap.im.lx, cap.im.ly);
	imx.setZero();
	imGabReIm.reallocate(cap.im.lx, cap.im.ly);

	L_VentanaImagen vent1(0, 0, imx.lx, imx.ly,"Click para cambiar parametros");
	L_VentanaImagen vent2(imx.lx + 25, 0, imx.lx, imx.ly,"Filtrado Gabor - click para salir");
	L_VentanaImagen vent3(0, imx.ly+100, imx.lx, imx.ly,"Amplitud Gabor");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		imx = cap.im;

		s = 1.6*lambda;
		imGabRe.filterGabor(imx, &imGabIm, s, lambda, theta);
		imGabReIm.reallocate(imGabRe.lx, imGabRe.ly);

		for (j=0; j<imGabReIm.ly; j++)
		{
			for (i=0; i<imGabReIm.lx; i++)
			{
				vald = imGabRe.pix(i,j);
				imGabReIm.pix(i,j,0) = (vald > 0.001) ? 255 : 0;
				imGabReIm.pix(i,j,1) = (vald > 0.001) ? 255 : 0;
				imGabReIm.pix(i,j,2) = (vald > 0.001) ? 255 : 0;
			}
		}

		imGabAmpl.reallocate(imGabRe.lx, imGabRe.ly);

		for (j=0; j<imGabAmpl.ly; j++)
			for (i=0; i<imGabAmpl.lx; i++)
				imGabAmpl.pix(i,j) = 20 * sqrt(imGabRe.pix(i,j)*imGabRe.pix(i,j) + imGabIm.pix(i,j)*imGabIm.pix(i,j));

		if (vent1.mouseE() == L_MouseArrastra)
		{
			dx = vent1.mouseX() - imGabReIm.lx/2;
			dy = vent1.mouseY() - imGabReIm.ly/2;
			lambda = sqrt(dx*dx+dy*dy) / 10;
			theta = -atan2(dy, dx);
		}
		vent1.mouseE() = L_MouseLeido;

		if (vent2.mouseE() != L_MouseLeido)
			break;

		t[frame%10] = L_TIME();
		dt = (t[frame%10] - t[(frame+1)%10]) / 10;
		printf("%f[ms]\n", dt*1000);

		vent1.dibuja(cap.im);
		vent2.dibuja(imGabReIm);
		imGabReIm = imGabAmpl;
		vent3.dibuja(imGabReIm);
		vent1.check();
	}

	return 0;
}
#endif

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_pruebaMeanShift(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	std::vector<double> hist;
	L_ImageRGBUchar im;
	L_ShapeArray lins;
	double xCen, yCen, rx, ry;
	int x1, y1, x2, y2;
	long frame = 0;

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	L_VentanaImagen vent1(0, 0, cap.im.lx, cap.im.ly,"Captura");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		if (vent1.mouseE() == L_MousePresiona)
		{
			x1 = vent1.mouseX();
			y1 = vent1.mouseY();
		}
		if (vent1.mouseE() == L_MouseSuelta)
		{
			x2 = vent1.mouseX();
			y2 = vent1.mouseY();
			xCen = (x1+x2)/2;
			yCen = (y1+y2)/2;
			rx = fabs((double)x1-x2) / 2;
			ry = fabs((double)y1-y2) / 2;
			printf("Inicial: cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
			cap.im.computeMeanShiftHistogram(xCen, yCen, hist, rx, ry);
		}
		vent1.mouseE() = L_MouseLeido;

		if (hist.size() > 0) // Existe el histograma
		{
			cap.im.computeMeanShiftTracking(xCen, yCen, hist, rx, ry, true);
			lins.drawRectangle((int)(xCen-rx), (int)(yCen-ry), (int)(xCen+rx), (int)(yCen+ry));
			im = cap.im;
			im.genDrawing(lins);
			vent1.dibuja(im);
			lins.resize(0);
			if (frame % 10 == 0)
				printf("cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
		}
		else
			vent1.dibuja(cap.im);
		vent1.check();
	}

	return 0;
}
#endif


#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_pruebaPiel(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	const double M[16][3] = {{73.53,29.94,17.76},{249.71,233.94,217.49},{161.68,116.25,96.95},{186.07,136.62,114.40},{189.26,98.37,51.18},{247.00,152.20,90.84},{150.10,72.66,37.76},{206.85,171.09,156.34},{212.78,152.82,120.04},{234.87,175.43,138.94},{151.19,97.74,74.59},{120.52,77.55,59.82},{192.20,119.62,82.32},{214.29,136.08,87.24},{99.57,54.33,38.06},{238.88,203.08,176.91}};
	const double C[16][3] = {{765.40,121.44,112.80},{39.94,154.44,396.05},{291.03,60.48,162.85},{274.95,64.60,198.27},{633.18,222.40,250.69},{65.23,691.53,609.92},{408.63,200.77,257.57},{530.08,155.08,572.79},{160.57,84.52,243.90},{163.80,121.57,279.22},{425.40,73.56,175.11},{330.45,70.34,151.82},{152.76,92.14,259.15},{204.90,140.17,270.19},{448.13,90.18,151.29},{178.38,156.27,404.99}};
	const double W[16][1] = {{0.0294},{0.0331},{0.0654},{0.0756},{0.0554},{0.0314},{0.0454},{0.0469},{0.0956},{0.0763},{0.1100},{0.0676},{0.0755},{0.0500},{0.0667},{0.0749}};
	const double MO[16][3] = {{254.37,254.41,253.82},{9.39,8.09,8.52},{96.57,96.95,91.53},{160.44,162.49,159.06},{74.98,63.23,46.33},{121.83,60.88,18.31},{202.18,154.88,91.04},{193.06,201.93,206.55},{51.88,57.14,61.55},{30.88,26.84,25.32},{44.97,85.96,131.95},{236.02,236.27,230.70},{207.86,191.20,164.12},{99.83,148.11,188.17},{135.06,131.92,123.10},{135.96,103.89,66.88}};
	const double CO[16][3] = {{2.77,2.81,5.46},{46.84,33.59,32.48},{280.69,156.79,436.58},{355.98,115.89,591.24},{414.84,245.95,361.27},{2502.24,1383.53,237.18},{957.42,1766.94,1582.52},{562.88,190.23,447.28},{344.11,191.77,433.40},{222.07,118.65,182.41},{651.32,840.52,963.67},{225.03,117.29,331.95},{494.04,237.69,533.52},{955.88,654.95,916.70},{350.35,130.30,388.43},{806.44,642.20,350.36}};
	const double WO[16][1] = {{0.0637},{0.0516},{0.0864},{0.0636},{0.0747},{0.0365},{0.0349},{0.0649},{0.0656},{0.01189},{0.0362},{0.0849},{0.0368},{0.0389},{0.0943},{0.0477}};
	double fM[16], fMO[16];
	struct A
	{
		double lutPiel[64][64][64], lutNoPiel[64][64][64];
	};

	struct A *a;
	a = new A;  // Muy grande para el stack, mejor en el heap

	L_ImageRGBUchar imPiel;
	double umbral = 2.0;
	double R, G, B, dist;
	double exponPiel, exponNoPiel, probPiel, probNoPiel, varPiel, varNoPiel;
	long frame = 0;
	int i, j, k, u;

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	L_VentanaImagen vent1(0, 0, cap.im.lx, cap.im.ly,"Captura");

	for (u=0; u<16; u++)
		fM[u] = 1/sqrt(2*M_PI*2*M_PI*2*M_PI*C[u][0]*C[u][1]*C[u][2]) * W[u][0];
	for (u=0; u<16; u++)
		fMO[u] = 1/sqrt(2*M_PI*2*M_PI*2*M_PI*CO[u][0]*CO[u][1]*CO[u][2]) * WO[u][0];

	for (i=0; i<64; i++)
	{
		for (j=0; j<64; j++)
		{
			for (k=0; k<64; k++)
			{
				R = i*4+2;
				G = j*4+2;
				B = k*4+2;
				probPiel = 0;
				probNoPiel = 0;
				for (u=0; u<16; u++)
				{
					exponPiel = 0;
					exponNoPiel = 0;
					varPiel = 0;
					varNoPiel = 0;
					//
					dist =  R - M[u][0];
					exponPiel += dist*dist / C[u][0];
					dist =  R - MO[u][0];
					exponNoPiel += dist*dist / CO[u][0];
					//
					dist =  G - M[u][1];
					exponPiel += dist*dist / C[u][1];
					dist =  G - MO[u][1];
					exponNoPiel += dist*dist / CO[u][1];
					//
					dist =  B - M[u][2];
					exponPiel += dist*dist / C[u][2];
					dist =  B - MO[u][2];
					exponNoPiel += dist*dist / CO[u][2];
					//
					probPiel += fM[u] * exp(-0.5*exponPiel);
					probNoPiel += fMO[u] * exp(-0.5*exponNoPiel);
				}
				a->lutPiel[i][j][k] = probPiel;
				a->lutNoPiel[i][j][k] = probNoPiel;
			}
		}
	}

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		imPiel.reallocate(cap.im.lx, cap.im.ly);
		
		for (j=0; j<imPiel.ly; j++)
		{
			for (i=0; i<imPiel.lx; i++)
			{
				probPiel = a->lutPiel[cap.im.pix(i,j,0)/4][cap.im.pix(i,j,1)/4][cap.im.pix(i,j,2)/4];
				probNoPiel = a->lutNoPiel[cap.im.pix(i,j,0)/4][cap.im.pix(i,j,1)/4][cap.im.pix(i,j,2)/4];
				if (probPiel / (probNoPiel + 1.0e-30) > umbral)
				{
					imPiel.pix(i,j,0) = cap.im.pix(i,j,0);
					imPiel.pix(i,j,1) = cap.im.pix(i,j,1);
					imPiel.pix(i,j,2) = cap.im.pix(i,j,2);
				}
				else
				{
					imPiel.pix(i,j,0) = 0;
					imPiel.pix(i,j,1) = 0;
					imPiel.pix(i,j,2) = 0;
				}
			}
		}

		if (vent1.mouseE() != L_MouseLeido)
			umbral = vent1.mouseX() / 60.0;
		vent1.mouseE() = L_MouseLeido;

		vent1.dibuja(imPiel);
		vent1.check();
		frame++;
	}
	delete a;
	return 0;
}
#endif


#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato::main_fractal(int argc, char *argv[])
{
	bool dibujando = false;
	double escape;
	double xc=0, yc=0, xc0=0, yc0=0;
	double e=0.01, lum=1;
	int coloreado = 1;
	L_ComplexDouble z, c(0.0,0.0), ex_z;
	L_Quaternion q, ex_q;
	bool usarCuat;
	L_ImageRGBUchar im, im2;
	int i, j, k;
	int W=400, H=400;
	int w=300, h=250;
	int facGr = 12;
	int nmax=0;
	bool salir = false;
	L_VentanaImagen vent1(0, 0, W, H,"Conjunto de Julia: Traslacion");
	L_VentanaImagen vent2(W+100, 0, w, h,"value de c");
	L_VentanaImagen vent3(W+100, h+50, w, h,"Escala y luminosidad");
	L_VentanaImagen vent4(W+100, 2*h+100, w, h, "Metodo de coloreo");
	L_VentanaImagen vent5(W+60, 2*h+150, w, h, "Clic para write");
	L_VentanaImagen vent6(W+30, 2*h+200, w, h, "Ingresar parametros");

	srand((int)time(NULL));

	while(true)
	{
		if (vent1.mouseE() == L_MousePresiona)
		{
			xc0 = xc - vent1.mouseX()*e;
			yc0 = yc - vent1.mouseY()*e;
		}
		if (vent1.mouseE() == L_MouseArrastra)
		{
			xc = xc0 + vent1.mouseX()*e;
			yc = yc0 + vent1.mouseY()*e;
		}
		vent1.mouseE() = L_MouseLeido;

		if (vent2.mouseE() == L_MousePresiona || vent2.mouseE() == L_MouseArrastra)
		{
			c.re = (vent2.mouseX() - w/2.0) / w * 8;
			c.im = (vent2.mouseY() - h/2.0) / h * 8;
		}
		vent2.mouseE() = L_MouseLeido;

		if (vent3.mouseE() == L_MousePresiona || vent3.mouseE() == L_MouseArrastra)
		{
			e = log(1 + (vent3.mouseX()*1.0) / w * 0.01);
			lum = (vent3.mouseY()*1.0) / h * 10;
			printf("Lum: %f\n", lum);
		}
		vent3.mouseE() = L_MouseLeido;

		if (vent4.mouseE() == L_MousePresiona || vent4.mouseE() == L_MouseArrastra)
		{
			coloreado++;
			coloreado = coloreado % 2;
		}
		vent4.mouseE() = L_MouseLeido;

		if (vent5.mouseE() == L_MousePresiona || vent5.mouseE() == L_MouseArrastra)
		{
			printf("Dibujando\n");
			W = facGr*W;
			H = facGr*H;
			e = e/facGr;
			dibujando = true;
		}
		vent5.mouseE() = L_MouseLeido;

		if (vent6.mouseE() == L_MousePresiona || vent6.mouseE() == L_MouseArrastra)
		{
			printf("Ingresar parametros:\n");
			printf("c.re,c.im: ");
			while (scanf("%lf%lf", &c.re,&c.im) != 2)
				{fflush(stdin);printf("** Ingrese 2 numeros separados por espacio y luego ENTER **\n");}
			printf("cenx,ceny: ");
			while (scanf("%lf%lf", &xc,&yc) != 2)
				{fflush(stdin);printf("** Ingrese 2 numeros separados por espacio y luego ENTER **\n");}
			printf("esc: ");
			while (scanf("%lf", &e) != 1)
				{fflush(stdin);printf("Ingrese un numero seguido de ENTER\n");}
			printf("lum: ");
			while (scanf("%lf", &lum) != 1)
				{fflush(stdin); printf("Ingrese un numero seguido de ENTER\n");}
		}
		vent6.mouseE() = L_MouseLeido;

		im.reallocate(W,H);
		for (j=0; j<H; j++)
		{
			for (i=0; i<W; i++)
			{
				usarCuat = true;
				escape = 0;
				z.re = (i - W/2)*e - xc;
				z.im = (j - H/2)*e - yc;

				if (usarCuat)
				{
					q.define(z.re, z.im, 0, 0);
					for (k=0; k<nmax; k++)
					{
						q = q*q*q*q*q + L_Quaternion(c.re,0,0,c.im);
						if (q.abs2() > 16)
							break;
						escape = q.abs2() / ex_q.abs2();
						ex_q = q;
					}
				}
				else
				{
					ex_z = z;
					for (k=0; k<nmax; k++)
					{
						//z.calcCos(); z = z + c; // Figuras complicadas se repiten en eje x
						//z.calcSqrt(); z = z + c;  // Figuras simples
						z = z*z*z*z*z + c;
						if (z.abs2() > 16)
							break;
						escape = z.abs2() / ex_z.abs2();
						ex_z = z;
					}
				}
				im.pix(i,j,0) = 0;
				im.pix(i,j,1) = 0;
				im.pix(i,j,2) = 0;

				switch(coloreado)
				{
				case 0:
					im.pix(i,j,2) = (k >= 256) ? 255 : (L_uchar)k;
					nmax = 255;
					break;
				case 1:
					escape = lum*escape;
					im.pix(i,j,0) = (escape >= 256) ? 255 : (L_uchar)escape;
					escape = lum*log(escape);
					im.pix(i,j,1) = (escape >= 256) ? 255 : (L_uchar)escape;
					escape = lum*log(escape);
					im.pix(i,j,2) = (escape >= 256) ? 255 : (L_uchar)escape;
					nmax = 10;
					break;
				}
			}
		}
		if (dibujando)
		{
			im.saveImage("sale.jpg");
			W = W/facGr;
			H = H/facGr;
			e = e*facGr;
			dibujando = false;
			FILE *fp;
			fp = fopen("sale.txt", "w");
			if (fp == NULL)
				printf("Error de grabacion de texto\n");
			fprintf(fp, "c = %f %f  cen = %f %f  es=%f  lum=%f  coloreado=%d  cuat=%d\n", c.re, c.im, xc, yc, e, lum, coloreado, usarCuat ? 1 : 0);
			fclose(fp);
			printf("ok\n");
		}
		vent1.dibuja(im);
		vent1.check();
	}

}
#endif


#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato_prueba::main_infoMutua(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imRef, imBN;
	L_ImageRGBUchar im, imRefRGB;
	L_ShapeArray lins;
	double xCen, yCen, rx, ry;
	double x=0, y=0, vx=0, vy=0;
	int x1, y1, x2, y2;
	long frame = 0;
	//
	double h12[8][8], h1[8], h2[8];

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	L_VentanaImagen vent1(0, 0, cap.im.lx, cap.im.ly,"Captura");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		if (vent1.mouseE() == L_MousePresiona)
		{
			x1 = vent1.mouseX();
			y1 = vent1.mouseY();
		}
		if (vent1.mouseE() == L_MouseSuelta)
		{
			x2 = vent1.mouseX();
			y2 = vent1.mouseY();
			xCen = (x1+x2)/2;
			yCen = (y1+y2)/2;
			rx = fabs((double)x1-x2) / 2;
			ry = fabs((double)y1-y2) / 2;
			printf("Inicial: cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
			imRefRGB.subImageOf(cap.im, x1, y1, x2-x1, y2-y1);
			imRef = imRefRGB;
			x = xCen;
			y = yCen;
			vx = 0;
			vy = 0;
		}
		vent1.mouseE() = L_MouseLeido;

		if (imRef.data() != NULL) // Existe el histograma
		{
			int i, j, li = x2-x1, lj = y2-y1;
			int u1, u2, xo, yo;
			double sum[4];
			imBN = cap.im;

			x = x+vx;
			y = y+vy;

			int dir, vez;

			for (vez = 5; vez > 0; vez--)
			{
				for (dir = 0; dir < 4; dir++)
				{
					xo = (int)x - (int)rx + ( (dir == 0 || dir == 1) ? -vez: +vez );
					yo = (int)y - (int)ry + ( (dir == 0 || dir == 2) ? -vez: +vez );

					for (i=0; i<8; i++)
					{
						h1[i] = 0;
						h2[i] = 0;
					}

					for (i=0; i<8; i++)
						for (j=0; j<8; j++)
							h12[i][j] = 0;
					
					if (xo < 0 || xo+li > imBN.lx || yo < 0 || yo+lj > imBN.ly)
						return 0;

					for (j=0; j<lj; j++)
					{
						for (i=0; i<li; i++)
						{
							u1 = (int)(imBN.pix(xo+i,yo+j) * 7.999);
							u2 = (int)(imRef.pix(i,j) * 7.999);
							h1[u1]++;
							h2[u2]++;
							h12[u1][u2]++;
						}
					}

					sum[dir] = 0;

					for (i=0; i<8; i++)
						for (j=0; j<8; j++)
							sum[dir] += h12[i][j] > 0  ?  h12[i][j]*log(h12[i][j]/(h1[i]*h2[j])) : 0;
				}
				if (sum[0] + sum[1] > sum[2] + sum[3])
					x -= vez;
				else
					x += vez;

				if (sum[0] + sum[2] > sum[1] + sum[3])
					y -= vez;
				else
					y += vez;
			}

			vx = 0.5*vx + 0.5*(x - xCen);
			vy = 0.5*vy + 0.5*(y - yCen);

			xCen = x;
			yCen = y;

			lins.drawRectangle((int)(xCen-rx), (int)(yCen-ry), (int)(xCen+rx), (int)(yCen+ry));
			im = cap.im;
			im.genDrawing(lins);
			vent1.dibuja(im);
			lins.resize(0);
			if (frame % 10 == 0)
				printf("cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
		}
		else
			vent1.dibuja(cap.im);
		vent1.check();
	}

	return 0;
}
#endif

#if defined (__COMPAT_FLTK__) || defined (__COMPAT_IPLIMAGE__)
int L_ProgramasFnesPato_prueba::main_infoMutuaRotEsc(int argc, char *argv[])
{
	L_CapturadorImagen cap;
	L_ImageGrayDouble imRef, imBN;
	L_ImageRGBUchar im, imRefRGB;
	L_ShapeArray lins;
	double xCen, yCen, rx, ry;
	double x=0, y=0, vx=0, vy=0;
	int x1, y1, x2, y2;
	long frame = 0;
	//
	double h12[8][8], h1[8], h2[8];

	if (cap.capturarImagen()==false)
	{
		printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
	}

	L_VentanaImagen vent1(0, 0, cap.im.lx, cap.im.ly,"Captura");

	while(true)
	{
		if (cap.capturarImagen()==false)
		{
			printf("No se puede capturar. Presione ENTER\n"); getchar(); return 1;
		}

		frame++;

		if (vent1.mouseE() == L_MousePresiona)
		{
			x1 = vent1.mouseX();
			y1 = vent1.mouseY();
		}
		if (vent1.mouseE() == L_MouseSuelta)
		{
			x2 = vent1.mouseX();
			y2 = vent1.mouseY();
			xCen = (x1+x2)/2;
			yCen = (y1+y2)/2;
			rx = fabs((double)x1-x2) / 2;
			ry = fabs((double)y1-y2) / 2;
			printf("Inicial: cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
			imRefRGB.subImageOf(cap.im, x1, y1, x2-x1, y2-y1);
			imRef = imRefRGB;
			x = xCen;
			y = yCen;
			vx = 0;
			vy = 0;
		}
		vent1.mouseE() = L_MouseLeido;

		if (imRef.data() != NULL) // Existe el histograma
		{
			int i, j, li = x2-x1, lj = y2-y1;
			int u1, u2, xo, yo;
			double sum[4];
			imBN = cap.im;

			x = x+vx;
			y = y+vy;

			int dir, vez;

			for (vez = 5; vez > 0; vez--)
			{
				for (dir = 0; dir < 4; dir++)
				{
					xo = (int)x - (int)rx + ( (dir == 0 || dir == 1) ? -vez: +vez );
					yo = (int)y - (int)ry + ( (dir == 0 || dir == 2) ? -vez: +vez );

					for (i=0; i<8; i++)
					{
						h1[i] = 0;
						h2[i] = 0;
					}

					for (i=0; i<8; i++)
						for (j=0; j<8; j++)
							h12[i][j] = 0;
					
					for (j=0; j<lj; j++)
					{
						for (i=0; i<li; i++)
						{
							u1 = (int)(imBN.pix(xo+i,yo+j) * 7.999);
							u2 = (int)(imRef.pix(i,j) * 7.999);
							h1[u1]++;
							h2[u2]++;
							h12[u1][u2]++;
						}
					}

					sum[dir] = 0;

					for (i=0; i<8; i++)
						for (j=0; j<8; j++)
							sum[dir] += h12[i][j] > 0  ?  h12[i][j]*log(h12[i][j]/(h1[i]*h2[j])) : 0;
				}
				if (sum[0] + sum[1] > sum[2] + sum[3])
					x -= vez;
				else
					x += vez;

				if (sum[0] + sum[2] > sum[1] + sum[3])
					y -= vez;
				else
					y += vez;
			}

			vx = 0.5*vx + 0.5*(x - xCen);
			vy = 0.5*vy + 0.5*(y - yCen);

			xCen = x;
			yCen = y;

			lins.drawRectangle((int)(xCen-rx), (int)(yCen-ry), (int)(xCen+rx), (int)(yCen+ry));
			im = cap.im;
			im.genDrawing(lins);
			vent1.dibuja(im);
			lins.resize(0);
			if (frame % 10 == 0)
				printf("cen = %.1f %.1f   radios = %.2f %.2f\n", xCen, yCen, rx, ry);
		}
		else
			vent1.dibuja(cap.im);
		vent1.check();
	}

	return 0;
}
#endif





int L_ProgramasFnesPato::main_pruebaMitos(int ntest)
{
	double t1=0, t2=0;
	long vez;
	double N;
	bool hacer_todos;
	// precalentamiento
	for (vez=0; vez < 1000; vez++)
		N = vez;
	t2 = L_TIME();

	hacer_todos = (ntest == -1);

	if (hacer_todos)
		ntest = 0;

	// stack frame
	if (ntest == 0 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
#ifdef _COMPAT_SSE2
		printf("(%d) Probando mult normal versus unrolling x 4 (double) vs SSE2\n", ntest-1);
#else
		printf("(%d) Probando mult normal versus unrolling x 4 (double)\n", ntest-1);
#endif
		std::vector<double> a(120), b(120), c(120);
		a.assign(a.size(), 2);
		b.assign(b.size(), 4);

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				for (int i=0; i<120; i++)
					c[i] = a[i] * b[i];
			}
			t2 = L_TIME();
		}
		printf(" >> Normal(double): %.2g[fps]\n", N/(t2-t1));

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				for (int i=0; i<120; i+=4)
				{
					c[i+0] = a[i+0] * b[i+0];
					c[i+1] = a[i+1] * b[i+1];
					c[i+2] = a[i+2] * b[i+2];
					c[i+3] = a[i+3] * b[i+3];
				}
			}
			t2 = L_TIME();
		}
		printf(" >> Unrolling x 4(double): %.2g[fps]\n", N/(t2-t1));
		
		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				double *cp=&c[0], *cpf=&c[0]+40-1, *ap=&a[0], *bp=&b[0];
				*cp = *ap * *bp;
				while (cp<cpf)
					*(++cp) = *(++ap) * (*(++bp));
			}
			t2 = L_TIME();
		}
		printf(" >> Punteros(double): %.2g[fps]\n", N/(t2-t1));

	}

	// stack frame
	if (ntest == 1 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		printf("(%d) Probando mult normal versus unrolling x 4 (float)\n", ntest-1);
		std::vector<float> a(40), b(40), c(40);
		a.assign(a.size(), 2);
		b.assign(b.size(), 4);

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				for (int i=0; i<40; i++)
					c[i] = a[i] * b[i];
			}
			t2 = L_TIME();
		}
		printf(" >> Normal(float): %.2g[fps]\n", N/(t2-t1));

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				for (int i=0; i<40; i+=4)
				{
					c[i+0] = a[i+0] * b[i+0];
					c[i+1] = a[i+1] * b[i+1];
					c[i+2] = a[i+2] * b[i+2];
					c[i+3] = a[i+3] * b[i+3];
				}
			}
			t2 = L_TIME();
		}
		printf(" >> Unrolling x 4(float): %.2g[fps]\n", N/(t2-t1));

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10000;
			for (vez=0; vez < 10000; vez++)
			{
				float *cp=&c[0], *cpf=&c[0]+40-1, *ap=&a[0], *bp=&b[0];
				*cp = *ap * *bp;
				while (cp<cpf)
					*(++cp) = *(++ap) * (*(++bp));
			}
			t2 = L_TIME();
		}
		printf(" >> Punteros(float): %.2g[fps]\n", N/(t2-t1));
	}

	// stack frame
	if (ntest == 2 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		printf("(%d) Probando peticion/liberacion memoria estatica versus memoria dinamica\n", ntest-1);

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				{
					double a[1000];
					a[0] = 1;
				}
			}
			t2 = L_TIME();
		}
		printf(" >> Estatica: %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				double *a;
				a = new double[1000];
				a[0] = 1;
				delete[] a;
			}
			t2 = L_TIME();
		}
		printf(" >> Dinamica: %.2g[fps]\n", N/(t2-t1));
	}

	if (ntest == 3 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double est[100][100];
		L_Matrix m;
		m.reallocate(100, 100);
		printf("(%d) Probando calculo en matriz estatica versus matriz dinamica\n", ntest-1);

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<100; i++)
					for (int j=0; j<100; j++)
						est[i][j] = i*j;
			}
			t2 = L_TIME();
		}
		printf(" >> Estatica: %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<100; i++)
					for (int j=0; j<100; j++)
						m(i,j) = i*j;
			}
			t2 = L_TIME();
		}
		printf(" >> Dinamica: %.2g[fps]\n", N/(t2-t1));
	}
	// stack frame
	if (ntest == 4 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double res;
		int i = 0;
		double a[5] = {1.0, 3.0, 2.0, 5.0, 4.0};
		printf("(%d)Probando funcion inline (mult) versus operacion directa\n", ntest-1);

		res=0;
		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				res += vez*(double)a[i%5];
			}
			t2 = L_TIME();
		}
		printf(" >> Directa: %.2g[fps]\n", N/(t2-t1));

		res = 0;
		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				res += L_multttt(vez, a[i%5]);
			}
			t2 = L_TIME();
		}
		printf(" >> Inline: %.2g[fps]\n", N/(t2-t1));
	}


	// stack frame
	if (ntest == 5 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double res;
		int i = 0;
		printf("(%d) Probando versus  +(==0)*  vs  ?: \n", ntest-1);

		res = 0;
		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				res += ((i++)%2) * vez + 0.0;
			}
			t2 = L_TIME();
		}
		printf(" >> +(==0)*: %.2g[fps]\n", N/(t2-t1));


		res = 0;
		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				res += ((i++)%2>0) ? vez : 0.0;
			}
			t2 = L_TIME();
		}
		printf(" >> ?:  : %.2g[fps]\n", N/(t2-t1));
	}

	// stack frame
	if (ntest == 6 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double res=0;
		int i=0;

		printf("(%d) Probando sum versus if\n", ntest-1);

		res = 0;
		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				res += ((i++)%2) * vez + 0;
			}
			t2 = L_TIME();
		}
		printf(">> += ()* : %.2g[fps]\n", N/(t2-t1));

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				if ((i++)%2>0)
					res += vez;
				else
					res += 0;
			}
			t2 = L_TIME();
		}
		printf(">> if (): %.2g[fps]\n", N/(t2-t1));
	}

	//stack frame
	if (ntest == 7 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		std::vector<double> a(2000), b(2000);
		a.assign(a.size(), 10.0);
		printf("(%d) Comparando for versus memcpy\n", ntest-1);

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<(int)a.size(); i++)
					b[i] = a[i];
			}
			t2 = L_TIME();
		}
		printf(">> for()  : %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				memcpy(&(b[0]), &(a[0]), sizeof(double)*a.size());
			}
			t2 = L_TIME();
		}
		printf(">> memcpy()  : %.2g[fps]\n", N/(t2-t1));
	}

	// stack frame
	if (ntest == 8 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_Matrix a, b, c;
		a.reallocate(3,3);
		a.identity();
		b.reallocate(3,3);
		b.identity();
		c.reallocate(3,3);
		c.identity();
		printf("(%d) Probando for()  vs  #define  vs  funcion sum  vs  operator + matr 3x3\n", ntest-1);

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<c.li; i++)
					for (int j=0; j<c.lj; j++)
						c(i,j) = a(i,j) + b(i,j);
			}
			t2 = L_TIME();
		}
		printf(">> for()  : %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				L_Matrix_OP_add(c,a,b);
			}
			t2 = L_TIME();
		}
		printf(">> #def f(): %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				c.OP_add(a,b);
			}
			t2 = L_TIME();
		}
		printf(">> c.f(a,b) : %.2g[fps]\n", N/(t2-t1));

#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				// Se usa + para benchmark
				c = a+b;
			}
			t2 = L_TIME();
		}
		printf(">> c = a+b : %.2g[fps]\n", N/(t2-t1));
#endif
	}

	// stack frame
	if (ntest == 9 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_Matrix a, b, c;
		a.reallocate(10,10);
		a.identity();
		b.reallocate(10,10);
		b.identity();
		c.reallocate(10,10);
		c.identity();
		printf("(%d) Probando for()  vs  #define  vs  funcion sum  vs  operator + matr 10x10\n", ntest-1);

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<c.li; i++)
					for (int j=0; j<c.lj; j++)
						c(i,j) = a(i,j) + b(i,j);
			}
			t2 = L_TIME();
		}
		printf(">> for()  : %.2g[fps]\n", N);

		N=0;
		t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				L_Matrix_OP_add(c,a,b);
			}
			t2 = L_TIME();
		}
		printf(">> #def f(): %.2g[fps]\n", N/(t2-t1));

		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				c.OP_add(a,b);
			}
			t2 = L_TIME();
		}
		printf(">> c.f(a,b) : %.2g[fps]\n", N);

#ifdef DEFINE_NONEFFICIENT_MATRIX_OPERATORS
		t1 = L_TIME();
		N=0;
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				// Se usa + para benchmark
				c = a+b;
			}
			t2 = L_TIME();
		}
		printf(">> c = a+b : %.2g[fps]\n", N/(t2-t1));
#endif
	}

	// stack frame
	if (ntest == 10 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_Matrix m1;
		double *m2ptr; // Necesario para que no se coma los if
		printf("(%d) Probando if(M.size() == 0...){printf} versus #define L_Matrix_checkIfDefined() matr 10x10\n", ntest-1);
		m1.reallocate(10,10);
		for (int i=0; i<10; i++)
			for (int j=0; j<10; j++)
				m1(i,j)=rand()%10;
		m2ptr = m1.begin();

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m2ptr == NULL || m1.li != 10 || m1.lj != 10)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> if(M.size() == 0...): %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				L_Matrix_checkIfDefined(m1,10,10);
			}
			t2 = L_TIME();
		}
		printf(" >> revisaDefinida: %.2g[fps]\n", N/(t2-t1));
	}

	//stack frame
	if (ntest == 11 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_Matrix M1(60,60), M1T, M2;
		L_StaticMatrix<60,60> ME;
		L_MatrixFloat Mf;
		M1.fillRandomValues(-100,100);
		L_StaticMatrix_OP_assign(ME,M1);
		Mf.copyFrom(M1);
		printf("(%d) Probando inversion usando varias alternativas matr 60x60\n", ntest-1);
		t1 = L_TIME();
		M1T.transpOf(M1);
		M1*=M1T; // M1 = M1*M1T

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (M1.invertMe() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> invertMe(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (M1.invertMegj() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> invertMegj(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (M1.invertMeDefPos() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> invertMeDefPos(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (M1.invertMeSimetricaDefPos() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> invertMeSimetricaDefPos(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (ME.invertMe() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> Estatica::invertMe(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (Mf.invertMeDefPos() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> float::invertMeDefPos(): %.3g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				if (Mf.invertMeDefPosx4() == false)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(" >> float::invertMeDefPosx4(): %.3g[fps]\n", N/(t2-t1));
	}

	// stack frame
	if (ntest == 12 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		struct L_imagenfloat
		{
			float **elem;
			int lx;
			int ly;
			inline float &operator()(int i, int j) {return elem[i][j];}
		};
		L_imagenfloat a, b, c;
		a.lx = b.lx = c.lx = 100;
		a.ly = b.ly = c.ly = 100;
		a.elem = L_new2d<float>(100,100);
		b.elem = L_new2d<float>(100,100);
		c.elem = L_new2d<float>(100,100);

		L_ImageGrayDouble aI, bI, cI;
		aI.reallocate(100,100);
		bI.reallocate(100,100);
		cI.reallocate(100,100);

		printf("(%d) Probando float[][] vs double[][] (sum) matr 100x100\n", ntest-1);

		for (int i=0; i<a.lx; i++)
			for (int j=0; j<a.ly; j++)
				a(i,j) = b(i,j) = c(i,j) = (float)(rand() % 10);

		for (int i=0; i<a.lx; i++)
			for (int j=0; j<a.ly; j++)
				aI(i,j) = bI(i,j) = cI(i,j) = rand() % 10;

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int jj=0; jj<c.ly; jj++)
					for (int ii=0; ii<c.lx; ii++)
						c(jj,ii) = a(jj,ii) + b(jj,ii);
			}
			t2 = L_TIME();
		}
		printf(" >> float[][]: %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				L_ImageGrayFloat_OP_add(cI, aI, bI);
			}
			t2 = L_TIME();
		}
		printf(" >> double[][]: %.2g[fps]\n", N/(t2-t1));

		L_delete2d<float>(a.elem);
		L_delete2d<float>(b.elem);
		L_delete2d<float>(c.elem);
	}

	// stack frame
	if (ntest == 13 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_ImageGrayDouble im1, im2, im3;
		int w = 11, wm=w/2;
		im1.reallocate(100, 100);
		im2.reallocate(w, w);
		im3.reallocate(100, 100);

		printf("(%d) Probando [i][j] vs inline(x,y) (conv) imag 800x600 & imag 30x30\n", ntest-1);

		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im1.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<w; j++)
			for (int i=0; i<w; i++)
				im2.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im3.pix(i,j) = 0;

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 2;
			for (vez=0; vez<2; vez++)
			{
				for (int j=w; j<im1.ly-w; j++)
					for (int i=w; i<im1.lx-w; i++)
						for (int v=0; v<w; v++)
							for (int u=0; u<w; u++)
								im3(j+v-wm,i+u-wm) += im1(j,i) * im2(v,u);
			}
			t2 = L_TIME();
		}
		printf(" >> im[y][x]: %.3g[fps]\n", N/(t2-t1));

		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im1.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<w; j++)
			for (int i=0; i<w; i++)
				im2.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im3.pix(i,j) = 0;

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 2;
			for (vez=0; vez<2; vez++)
			{
				for (int j=w; j<im1.ly-w; j++)
					for (int i=w; i<im1.lx-w; i++)
						for (int v=0; v<w; v++)
							for (int u=0; u<w; u++)
								im3.pix(i+u-wm,j+v-wm) += im1.pix(i,j) * im2.pix(u,v);
			}
			t2 = L_TIME();
		}
		printf(" >> im.pix(x,y): %.3g[fps]\n", N/(t2-t1));


		printf("Probando [i][j] vs inline(x,y) (conv) traspuesto imag 800x600 & imag 30x30\n");

		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im1.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<w; j++)
			for (int i=0; i<w; i++)
				im2.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im3.pix(i,j) = 0;

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 2;
			for (vez=0; vez<2; vez++)
			{
				for (int i=w; i<im1.lx-w; i++)
					for (int j=w; j<im1.ly-w; j++)
						for (int u=0; u<w; u++)
							for (int v=0; v<w; v++)
								im3(j+v-wm,i+u-wm) += im1(j,i) * im2(v,u);
			}
			t2 = L_TIME();
		}
		printf(" >> im[y][x] tr: %.3g[fps]\n", N/(t2-t1));

		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im1.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<w; j++)
			for (int i=0; i<w; i++)
				im2.pix(i,j) = L_RANDOM(-1, 1);
		for (int j=0; j<im1.ly; j++)
			for (int i=0; i<im1.lx; i++)
				im3.pix(i,j) = 0;

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 2;
			for (vez=0; vez<2; vez++)
			{
				for (int i=w; i<im1.lx-w; i++)
					for (int j=w; j<im1.ly-w; j++)
						for (int u=0; u<w; u++)
							for (int v=0; v<w; v++)
								im3.pix(i+u-wm,j+v-wm) += im1.pix(i,j) * im2.pix(u,v);
			}
			t2 = L_TIME();
		}
		printf(" >> im.pix(x,y) tr: %.3g[fps]\n", N/(t2-t1));
	}


	// stack frame
	if (ntest == 14 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double *d1, *d2;
		double arr1[1]={100.0}, arr2[1]={200.0};
		d1 = arr1;
		d2 = arr2;
		printf("(%d) Probando f() {}   vs   extern C f() {}   vs   f() throw() {}   vs   __fastcall f() (todas inline)\n", ntest-1);


		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				arr1[0] = L_multtttp(d1, d2);
			}
			t2 = L_TIME();
		}
		printf(" >> f() {}: %.2g[fps]\n", N/(t2-t1));


		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				arr1[0] = L_multtttpC(d1, d2);
			}
			t2 = L_TIME();
		}
		printf(" >> extern C f() {}: %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				arr1[0] = L_multtttpNT(d1, d2);
			}
			t2 = L_TIME();
		}
		printf(" >> f() throw() {}: %.2g[fps]\n", N/(t2-t1));
	}


	// stack frame
	if (ntest == 15 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		L_Matrix m1;
		double *m2ptr; // Necesario para que no se coma los if
		printf("(%d) Probando if() throw;  vs  if() printf(...);  vs  if() printf();  vs   if() j=j*2;  vs  if() *NULL=0; matr 10x10\n", ntest-1);
		m1.reallocate(10,10);
		for (int i=0; i<10; i++)
			for (int j=0; j<10; j++)
				m1(i,j)=rand()%10;
		m2ptr = m1.begin(); // cuidado

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					throw L_ArgException();
			}
			t2 = L_TIME();
		}
		printf(">> if() throw; : %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					printf("%d %d %d %ld %ld\n", m1.li, m1.lj, m1.li*m1.lj, vez*m1.li, vez*m1.lj);
			}
			t2 = L_TIME();
		}
		printf(">> if() printf(...); : %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					printf("Mala\n");
			}
			t2 = L_TIME();
		}
		printf(">> if() printf(); : %.2g[fps]\n", N/(t2-t1));


		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			int j=0;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					j=j*2;
			}
			t2 = L_TIME();
		}
		printf(">> if() j=j*2;  : %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					*(int*)NULL=0;
			}
			t2 = L_TIME();
		}
		printf(">> if() *NULL=0;  : %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				m2ptr = m1.begin();
				if (m1.begin() == NULL || m1.li != 10 || m1.lj != 10)
					*(int*)128=0;
			}
			t2 = L_TIME();
		}
		printf(">> if() *128=0;  : %.2g[fps]\n", N/(t2-t1));
	}


	// stack frame
	if (ntest == 16 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		int nc = 640, nf = 480;
		printf("(%d) Probando velocidad new[][] + memcpy en imagenes\n", ntest-1);
		L_ImageGrayDouble im1d(nc,nf), im2d(nc,nf);
		L_ImageGrayFloat im1f(nc,nf), im2f(nc,nf);
		L_ImageRGBUchar im1RGB(nc,nf), im2RGB(nc,nf);

		for (int i=0; i<im1d.li; i++)
			for (int j=0; j<im1d.lj; j++)
				im1d(i,j) = L_RANDOM(0,1);

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2d.destroy();
				im2d = im1d;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 640x480 double : %.6g[fps]\n", N/(t2-t1));

		for (int i=0; i<im1f.li; i++)
			for (int j=0; j<im1f.lj; j++)
				im1f(i,j) = (float)L_RANDOM(0,1);

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2f.destroy();
				im2f = im1f;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 640x480 float : %.6g[fps]\n", N/(t2-t1));

		im1RGB.fillRandomValues();

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2RGB.destroy();
				im2RGB = im1RGB;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 640x480 char[3] : %.6g[fps]\n", N/(t2-t1));

		/////

		im1d.reallocate(1024,768);
		for (int i=0; i<im1d.li; i++)
			for (int j=0; j<im1d.lj; j++)
				im1d(i,j) = L_RANDOM(0,1);

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)

		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2d.destroy();
				im2d = im1d;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 1024x768 double : %.6g[fps]\n", N/(t2-t1));

		im1f.reallocate(1024, 768);
		for (int i=0; i<im1f.li; i++)
			for (int j=0; j<im1f.lj; j++)
				im1f(i,j) = (float)L_RANDOM(0,1);

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2f.destroy();
				im2f = im1f;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 1024x768 float : %.6g[fps]\n", N/(t2-t1));

		im1RGB.reallocate(1024, 768);
		im1RGB.fillRandomValues();

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 10;
			for (vez=0; vez<10; vez++)
			{
				im2RGB.destroy();
				im2RGB = im1RGB;
			}
			t2 = L_TIME();
		}
		printf(">> new[][] + memcpy 1024x768 char[3] : %.6g[fps]\n", N/(t2-t1));
	}


	if (ntest == 17 || hacer_todos == true)
	{
		ntest += hacer_todos ? 1 : 0;
		double **pptr;
		double *ptr;
		int li = 1200, lj = 800, ljStep = sizeof(double)*800;

		pptr = L_new2d<double>(li, ljStep);
		ptr = (double *)L_aligned_malloc_128(li*ljStep);

		printf("Tiempos de acceso (escritura m[i][j] = i*j+vez)\n");

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<li; i++)
					for (int j=0; j<lj; j++)
						pptr[i][j] = i*j+vez;
			}
			t2 = L_TIME();
		}
		printf(">> elem[i][j] : %.2g[fps]\n", N/(t2-t1));

		N=0; t1 = L_TIME();
		while(t2-t1 < 1.0)
		{
			N+= 100;
			for (vez=0; vez<100; vez++)
			{
				for (int i=0; i<li; i++)
					for (int j=0; j<lj; j++)
						((double *)((char *)ptr + i*ljStep))[j] = i*j+vez;
			}
			t2 = L_TIME();
		}
		printf(">> (elem+i*ljStep)[j]; : %.2g[fps]\n", N/(t2-t1));
		L_delete2d<double>(pptr);
		L_aligned_free_128(ptr);
	}
	printf("Presione ENTER para salir\n");
	getchar();
	return 0;
}
