#ifndef __L_GENPUNTINT_H__
#define __L_GENPUNTINT_H__
#include "L_Descriptores.h"

////////////////////////////////////
/////
///// Clases de nivel intermedio para define generadores de puntos de interés
/////
//////////////////////////////


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




/*!
  class_name abstracta, cuyos children son generadores de descriptores.
  Todos los children deben cumplir las siguientes condiciones:<br>
  - Deben tener una funcion calcDescrPuntInt( ) que llene la lista desL con descriptores. El llenado de punL y dirL es opcional.<br>
  - Deben colocar el numero de object descriptor->nobj igual a la variable nobj en cada uno de los descriptores creados<br>
  - Deben reserve memoria para descriptor->vector del mismo length de descriptor->nt para cada uno de los descriptores creados<br>
  - Deben agregar los parametros de esta class_name a su servidor de parametros, si tienen alguno<br>
  - Deben hacer que descriptor->refgen apunte a "this", es decir, que apunte al mismo generador de descriptores<br>
*/

// Para ver el codeMapping antigu, ver __DIFGAUS_H_ ... :)

class L_GenDescrPuntInt_POD
{
public:
	int nobj; // numero de object que se colocara a los descriptores creados
	// Parametros (deben agregarse al servidor de parametros de los children, si dicho servidor existe)
	int minX;
	int minY;
	double calces_maxDist;   // Peor distancia entre descriptores que es aceptable al generar calces
	double calces_maxRatio;  // Peor relacion dist_primero_mas_cercano / dist_segundo_mas_cercano que es aceptable al generar calces
	int calces_numCalcesPorDescr; // Numero de calces que se generaran para cada descriptor de la imagen de test
	double calces_probCorrecto;  // Probabilidad de que un calce generado usando estos descriptores sea correcto (?)
	double calces_maxDifIntens; // Cota para la proporción valMin/valMax de los valores "DoG" de los puntos de interés de un calce

	L_PuntIntTipo tipoP;  // El tipo de puntos de interes que build
	L_DescriptorTipo tipoD;  // El tipo de descriptores que build
};

class L_GenDescrPuntInt : public L_GenDescrPuntInt_POD // Cuidado con nobj: debe fijarse este value en todos los descriptores
{
public:
	L_PuntIntLista punL; // almacenamiento temporal
	L_DireccLista dirL;  // almacenamiento temporal
	L_DescriptorLista desL; // almacenamiento definitivo (lo util)

	L_GenDescrPuntInt()
	{
		calces_maxDist = 2.0;
		calces_maxRatio = 0.9;
		calces_numCalcesPorDescr = 1;
		calces_probCorrecto = 0.01;
		calces_maxDifIntens = 0.0;
	}
	virtual bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0)=0; // La funcion util externamente para los detectores
	virtual L_ParamManagerLocal *pideParams() {return NULL;}

	void dibFlechasDescriptores(L_ShapeArray &lins);
	void dibCuadradosDescriptores(L_ShapeArray &lins);

	void swap(L_GenDescrPuntInt &other)
	{
		L_GenDescrPuntInt_POD g_t;
		g_t.L_GenDescrPuntInt_POD::operator=(other);
		other.L_GenDescrPuntInt_POD::operator=(*this);
		this->L_GenDescrPuntInt_POD::operator=(g_t);
		punL.swap(other.punL);
		dirL.swap(other.dirL);
		desL.swap(other.desL);
	}

	virtual ~L_GenDescrPuntInt() {}
};

class L_SSS // scale space submuestreado, solo se submuestrea a la mitad del tamano
{
public:
	L_ParamManagerLocal paramSSS;
	L_Array<L_ImageGrayDouble> imarr; // Octava actual, nie imagenes
	L_ImageGrayDouble *imorig; // puntero a la imagen original, sin memoria propia
	std::vector<double> sigm; // sigmas deseados para las imagenes del scalespace submuestreado, sig[0]=s1
	std::vector<double> sigmrel; // sigm=s1*sigmrel
	int nsub; // nº de submuestreos hechos a la imagen
	L_Frac reltam; // relacion de tamano entre imorig y imarr[0]. En general es <1
//param
	int nio; // num imagenes por octava
	int nie; // num imagenes por escalon, nie>=nio
	double s0;
	double s1; // sigma deseado para la primera imagen del scale space submuestreado cuando se doubleResolution la imagen, s1>=s0
	double s1_ND; // sigma deseado para la primera imagen del scale space submuestreado cuando no se doubleResolution la imagen, s1_ND>=s0
	double sMin; // El minElement sigma que se puede usar para hacer una convolucion gaussian (sMin>100 => convoluciones no submuestreadas)
	bool dupl; // indica si se debe duplicar inicialmente la imagen
	bool submEfic;
	bool usasqrt2;

	L_SSS();
	L_SSS(int nio_, int nie_);
	void fijaImagenOriginal(L_ImageGrayDouble& im, double sigma=0.5);
	void borraOctava();
	bool nuevaOctava();
	void grabaArchivosBMP(const L_FileName &nomarch);
	void destroy();
	virtual ~L_SSS() {destroy();}
};

class L_PSS // scale space piramide
{
#if defined(L_PDS_NPIR)
	#error #define conflicts
#endif
#define L_PDS_NPIR 12 // Numero de pisos de las piramides maxElement a calcular
public:
	L_ParamManagerLocal paramPSS;
	L_ImageGrayDouble imact;  // Imagen actual
	L_ImageGrayDouble imactG; // Imagen actual grande
	L_ImageGrayDouble imantG; // Imagen anterior grande
	L_ImageGrayDouble *imorig; // puntero a la imagen original, sin memoria propia
	L_Frac porcSub; // porcentaje para submuestrear. Valores posibles: {1/2, 2/3}. En general es <1
	int nsub; // nº de submuestreos hechos a la imagen
	L_Frac reltam; // relacion de tamano entre imorig y imactG. En general es <1
//param
	double s0; // sigma de la imagen original
	double s1; // sigma deseado para la imagen actual cuando la imagen se doubleResolution inicialmente, s1>=s0
	double s1_ND; // sigma deseado para la imagen actual cuando la imagen no se doubleResolution inicialmente, s1_ND>=s0
	double sMin; // El minElement sigma que se puede usar para hacer una convolucion gaussian (sMin>100 => convoluciones no submuestreadas)
	bool dupl; // Indica si duplicar o no la imagen inicial
	bool usasqrt2; // Usar las funciones sqrt2 para hacer las convoluciones
	bool imprSigmas; // Imprimir los sigmas de las convoluciones hechas

private:
	L_PSS(L_PSS& other); // {}  No debe existir
public:

	L_PSS();
	L_PSS(L_Frac &f);
	void fijaImagenOriginal(L_ImageGrayDouble& im, double sigma=0.5);
	bool nuevaImagen();
	void grabaArchivosBMP(const L_FileName &nomarch);
	~L_PSS() {};
};

class L_STrioEstat
{
public:
	static void calc_M_BTB_BT(L_Matrix &M_BTB_BT);
	static double inline ajustaPar(double u, double v, double w) {return 0.5*(w-u)/(2*v-u-w);}
};

class L_PTrioEstat
{
public:
	static void calc_M_BTB_BT(L_Matrix &M_BTB_BT);
	static double inline ajustaPar(double u, double v, double w) {return 0.5*(w-u)/(2*v-u-w);}
};

class L_STrioIm
{
public:
	L_ImageGrayDouble im[3]; // usar * puede ser algo lento
	int n; // nº de capa scalespace = nº de trio

	L_Matrix M_BTB_BT;

	L_STrioIm() {n=0;}
	//L_STrioIm(L_SSS &p) {n=0;}
	inline bool esMax3D(int i, int j) {return im[1].esMax2D(i,j) && im[0].esMax2D(im[1].pix(i,j),i,j) && im[2].esMax2D(im[1].pix(i,j),i,j);}
	inline bool esMax3D(int i, int j, double porc) {return im[1].esMax2D(i,j,porc) && im[0].esMax2D(im[1].pix(i,j),i,j,porc) && im[2].esMax2D(im[1].pix(i,j),i,j,porc);}
	inline bool esMin3D(int i, int j) {return im[1].esMin2D(i,j) && im[0].esMin2D(im[1].pix(i,j),i,j) && im[2].esMin2D(im[1].pix(i,j),i,j);}
	inline bool esMin3D(int i, int j, double porc) {return im[1].esMin2D(i,j,porc) && im[0].esMin2D(im[1].pix(i,j),i,j,porc) && im[2].esMin2D(im[1].pix(i,j),i,j,porc);}
	inline bool esMax1D(int i, int j) {return im[1].pix(i,j)>im[0].pix(i,j) && im[1].pix(i,j)>im[2].pix(i,j);}
	inline bool esMax1D(int i, int j, double porc) {return im[1].pix(i,j)*porc>im[0].pix(i,j) && im[1].pix(i,j)*porc>im[2].pix(i,j);}
	inline bool esMin1D(int i, int j) {return im[1].pix(i,j)<im[0].pix(i,j) && im[1].pix(i,j)<im[2].pix(i,j);}
	inline bool esMin1D(int i, int j, double porc) {return im[1].pix(i,j)*porc<im[0].pix(i,j) && im[1].pix(i,j)*porc<im[2].pix(i,j);}
	void ajustaParabol(int i, int j, L_Matrix &grad, L_Matrix &hess);
	void ajustaTaylor1(int i, int j, L_Matrix &grad, L_Matrix &hess);
	double interpola(L_PuntInt &p,  L_Matrix &xmax, const L_Matrix &grad, const L_Matrix &gradtr, const L_Matrix& hess, const L_Matrix &hessinv);
	~L_STrioIm() {}
};


class L_STrioImRef
{
public:
	L_ImageGrayDouble *im[3]; // usar * puede ser algo lento
	int n; // nº de capa scalespace = nº de trio
	L_STrioImRef() {n=0; im[0]=NULL; im[1]=NULL; im[2]=NULL;}
	//L_STrioImRef(L_SSS &p) {n=0; im[0]=NULL; im[1]=NULL; im[2]=NULL;}
	inline bool esMax3D(int i, int j){return im[1]->esMax2D(i,j) && im[0]->esMax2D(im[1]->pix(i,j),i,j) && im[2]->esMax2D(im[1]->pix(i,j),i,j);}
	inline bool esMin3D(int i, int j){return im[1]->esMin2D(i,j) && im[0]->esMin2D(im[1]->pix(i,j),i,j) && im[2]->esMin2D(im[1]->pix(i,j),i,j);}
	inline bool esMax1D(int i, int j){return im[1]->pix(i,j)>im[0]->pix(i,j) && im[1]->pix(i,j)>im[2]->pix(i,j);}
	inline bool esMin1D(int i, int j){return im[1]->pix(i,j)<im[0]->pix(i,j) && im[1]->pix(i,j)<im[2]->pix(i,j);}
	~L_STrioImRef() {}
};

#define L_PTrioIm_esMax3D(T,i,j) (L_ImageGrayFloat_isMax2D(*(T).im[1],(i),(j)) && L_ImageGrayFloat_isMax2D_c(*(T).im[0],(T).im[1]->pix((i),(j)),((i)*(T).f01->den+1)/(T).f01->num,((j)*(T).f01->den+1)/(T).f01->num) && L_ImageGrayFloat_isMax2D_c(*(T).im[2],(T).im[1]->pix((i),(j)),(i)*(T).f12->num/(T).f12->den,(j)*(T).f12->num/(T).f12->den))
#define L_PTrioIm_esMin3D(T,i,j) (L_ImageGrayFloat_isMin2D(*(T).im[1],(i),(j)) && L_ImageGrayFloat_isMin2D_c(*(T).im[0],(T).im[1]->pix((i),(j)),((i)*(T).f01->den+1)/(T).f01->num,((j)*(T).f01->den+1)/(T).f01->num) && L_ImageGrayFloat_isMin2D_c(*(T).im[2],(T).im[1]->pix((i),(j)),(i)*(T).f12->num/(T).f12->den,(j)*(T).f12->num/(T).f12->den))

class L_PTrioIm_POD
{
public:
	L_ImageGrayDouble *im[3];
	L_Frac *f01; // Referencia sin memo propia, usualmente es < 1. SE NECESITAN 2 PARA EL SURF
	L_Frac *f12; // referencia sin memo propia, usualmente es < 1. SE NECESITAN 2 PARA EL SURF
	int n; // nº de capa scalespace = nº de trio

	void swap(L_PTrioIm_POD &other) {L_PTrioIm_POD t = other; other = *this; *this = t;}
};

class L_PTrioIm : public L_PTrioIm_POD
{
private:
	L_PTrioIm(L_PTrioIm &dontCallMe);
public:
	L_Matrix M_BTB_BT;

	L_PTrioIm() {f01=NULL; f12=NULL; n=0;}

	void swap(L_PTrioIm &other) {L_PTrioIm_POD::swap(other); M_BTB_BT.swap(other.M_BTB_BT);}

#if defined(L_PTI1)
	#error #define conflicts
#endif
#define L_PTI1 ((1))

	inline double pixTrio(int n, int i, int j) {if (n==1) return im[1]->pix(i,j); else if (n>1) return im[2]->pix(i*f12->num/f12->den,j*f12->num/f12->den); else return im[0]->pix((i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num);}
	inline bool esMax3D(int i, int j)
	{
		return im[1]->esMax2D(i,j) && im[0]->esMax2D(im[1]->pix(i,j),(i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[2]->esMax2D(im[1]->pix(i,j),i*f12->num/f12->den,j*f12->num/f12->den);
	}
	inline bool esMin3D(int i, int j)
	{
		return im[1]->esMin2D(i,j) && im[0]->esMin2D(im[1]->pix(i,j),(i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[2]->esMin2D(im[1]->pix(i,j),i*f12->num/f12->den,j*f12->num/f12->den);
	}
	inline bool esMax1D_1(int i, int j)
	{
		return im[1]->pix(i,j)>im[0]->pix((i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[1]->pix(i,j)>im[2]->pix(i*f12->num/f12->den,j*f12->num/f12->den);
	}
	inline bool esMin1D_1(int i, int j)
	{
		return im[1]->pix(i,j)<im[0]->pix((i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[1]->pix(i,j)<im[2]->pix(i*f12->num/f12->den,j*f12->num/f12->den);
	}
	inline bool esMax1D_1_abs(int i, int j)
	{
		return
			(
			im[1]->pix(i,j) > 0 &&
			im[1]->pix(i,j)>im[0]->pix((i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[1]->pix(i,j)>im[2]->pix(i*f12->num/f12->den,j*f12->num/f12->den)
		) || (
			im[1]->pix(i,j) < 0 &&
			im[1]->pix(i,j)<im[0]->pix((i*f01->den+L_PTI1)/f01->num,(j*f01->den+L_PTI1)/f01->num) &&
			im[1]->pix(i,j)<im[2]->pix(i*f12->num/f12->den,j*f12->num/f12->den)
		);
	}
	inline bool esMax1D_5(int i, int j)
	{
		return  esMax1D_1(i-1,j) 
			&& esMax1D_1(i+0,j-1) && esMax1D_1(i+0,j) && esMax1D_1(i+0,j+1)
			&& esMax1D_1(i+1,j);
	}
	inline bool esMin1D_5(int i, int j)
	{
		return esMin1D_1(i-1,j)
			&& esMin1D_1(i+0,j-1) && esMin1D_1(i+0,j) && esMin1D_1(i+0,j+1)
			&& esMin1D_1(i+1,j);
	}
	inline bool esMax1D_9(int i, int j)
	{
		return esMax1D_1(i-1,j-1) && esMax1D_1(i-1,j) && esMax1D_1(i-1,j+1)
			&& esMax1D_1(i+0,j-1) && esMax1D_1(i+0,j) && esMax1D_1(i+0,j+1)
			&& esMax1D_1(i+1,j-1) && esMax1D_1(i+1,j) && esMax1D_1(i+1,j+1);
	}
	inline bool esMin1D_9(int i, int j)
	{
		return esMin1D_1(i-1,j-1) && esMin1D_1(i-1,j) && esMin1D_1(i-1,j+1)
			&& esMin1D_1(i+0,j-1) && esMin1D_1(i+0,j) && esMin1D_1(i+0,j+1)
			&& esMin1D_1(i+1,j-1) && esMin1D_1(i+1,j) && esMin1D_1(i+1,j+1);
	}
	void ajustaParabol(int i, int j, L_Matrix &grad, L_Matrix &hess);
	void ajustaTaylor1(int i, int j, L_Matrix &grad, L_Matrix &hess);
	double interpola(L_PuntInt &p,  L_Matrix &xmax, const L_Matrix &grad, const L_Matrix &gradtr, const L_Matrix& hess, const L_Matrix &hessinv);
	~L_PTrioIm() {}
};
#undef L_PTI1

class L_PTrioImRef
{
public:
	L_ImageGrayDouble *im[3]; // Usar *`puede ser algo lento
	L_Frac *f; // Referencia sin memo propia
	int n; // nº de capa scalespace = nº de trio
	bool mempropia;

	L_PTrioImRef() {f=NULL; n=0; im[0]=NULL; im[1]=NULL; im[2]=NULL;}
	L_PTrioImRef(L_PSS &p) {f=&p.porcSub; n=0; im[0]=NULL; im[1]=NULL; im[2]=NULL;}
	inline bool esMax3D(int i, int j)
	{
		return im[1]->esMax2D(i,j) && im[0]->esMax2D(im[1]->pix(i,j),i*f->den/f->num,j*f->den/f->num) &&
			im[2]->esMax2D(im[1]->pix(i,j),i*f->num/f->den,j*f->num/f->den);
	}
	inline bool esMin3D(int i, int j)
	{
		return im[1]->esMin2D(i,j) && im[0]->esMin2D(im[1]->pix(i,j),i*f->den/f->num,j*f->den/f->num) &&
			im[2]->esMin2D(im[1]->pix(i,j),i*f->num/f->den,j*f->num/f->den);
	}
	inline bool esMax1D(int i, int j)
	{
		return im[1]->pix(i,j)>im[0]->pix(i*f->den/f->num,j*f->den/f->num) &&
			im[1]->pix(i,j)>im[2]->pix(i*f->num/f->den,j*f->num/f->den);
	}
	inline bool esMin1D(int i, int j)
	{
		return im[1]->pix(i,j)<im[0]->pix(i*f->den/f->num,j*f->den/f->num) &&
			im[1]->pix(i,j)<im[2]->pix(i*f->num/f->den,j*f->num/f->den);
	}
	~L_PTrioImRef() {}
};

// Usar mejor L_PDoG_SIFT en vez de este
class L_SDoG_SIFT:public L_SSS,public L_GenDescrPuntInt // detector DoG scale space submuestreado
{
public:
	L_ParamManagerLocal paramSDOGSIFT;
	L_STrioIm trioSDoG; // imagenes
	L_STrioIm trioGradR, trioGradA; // modulo y angulo del gradiente de trioIm
	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros
	int it; // numero de trio en la octava que se esta procesando
	long ultId;
// param
	double r_lin; // relacion (mayor value propio del hessiano)/(menor value propio del hessiano) para descartar lineas
	double val_minP_DoG; // minElement value del punto maxElement en el SDoG que es aceptable a priori
	double val_min_DoG;  // minElement value del punto maxElement en el SDoG que es aceptable tras aprox parabolica
	int usa_parabol; // Interpolacion del maxElement: 0 = taylor, 1 = parabola 3d, 2 = parabola 1d x 3
	double porcCurvMax;
	double s0_def;
	double s1_def;
	double s1_def_chico;

	L_SDoG_SIFT();

	bool calcPuntInt_act();
	bool calcDirecc_act();
	bool calcDescr_act();
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0);

	L_ParamManagerLocal *pideParams();

#ifdef __DIFGAUS_H_ // incluir al sistema antiguo
	static SIFT1 *calcsift(uchar ***im3D, long lx, long ly);
#endif

	virtual ~L_SDoG_SIFT() {}
};

class L_PDoG_SIFT:public L_GenDescrPuntInt // detector DoG piramide scale space
{
public:
	L_ParamManagerLocal paramPDOGSIFT;
	L_PSS genPiram;
	L_PTrioIm trioPDoG; // imagenes
	L_PTrioIm trioGradR, trioGradA; // modulo y angulo del gradiente de trioIm
	L_ImageGrayDouble pirGauss[L_PDS_NPIR];
	L_ImageGrayDouble pirLaplace[L_PDS_NPIR];
	L_ImageGrayDouble pirGradR[L_PDS_NPIR];
	L_ImageGrayDouble pirGradA[L_PDS_NPIR];
	double sigm[L_PDS_NPIR];
	L_Frac reltam[L_PDS_NPIR];
	int nPisos;
	long ultId;

	double r_lin; // relacion (mayor value propio del hessiano)/(menor value propio del hessiano) para descartar lineas
	double val_minP_DoG;
	int usa_parabol;
	double porcCurvMax;
	double s0_def;
	double s1_def;
	double s1_def_chico;

	L_PTrioIm trioDoG;
	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros

	L_PDoG_SIFT();

	bool calcPiramides(L_ImageGrayDouble &im, double s0);
	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcDescriptores();

	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams() {return &paramPDOGSIFT;}
	virtual ~L_PDoG_SIFT() {}
};

class L_DoG_SIFT:public L_GenDescrPuntInt // detector DoG no invariante a escala
{
public:
	L_SIFTCalculador gen;
	long ultId;
	L_DoG_SIFT() {tipoP = DOG; tipoD = SIFT; nobj=-1; ultId = 0;}
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0);
	virtual ~L_DoG_SIFT() {}
};

class L_SHarris_SIFT:public L_SSS,public L_GenDescrPuntInt // detector Harris sale space submuestreado
{
public:
	L_ParamManagerLocal paramSHARRISSIFT;
	L_STrioIm trioSDoG; // imagenes
	L_STrioIm trioGradR, trioGradA; // modulo y angulo del gradiente de trioIm
	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros
	int it; // numero de trio en la octava que se esta procesando
// param
	double factorS;
	double alfa;
	double val_minP_DoG; // minElement value del punto maxElement en DoG que es aceptable a priori
	double val_minP_Harris; // minElement value del punto maxElement en Harris que es aceptable a priori
	double porcCurvMax;
	double s0_def;
	double s1_def;
	double s1_def_chico;
	int usa_parabol;
	long ultId;

	L_SHarris_SIFT();

	bool calcPuntInt_act();
	bool calcDirecc_act();
	bool calcDescr_act();
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams();

#ifdef __DIFGAUS_H_ // incluir al sistema antiguo
	static SIFT1 *calcsift(uchar ***im3D, long lx, long ly);
#endif

	virtual ~L_SHarris_SIFT() {}
};

class L_PHarris_SIFT:public L_GenDescrPuntInt // detector Harris piramide scale space
{
public:
	L_ParamManagerLocal paramPHARSIFT;
	L_PSS genPiram;
	L_PTrioIm trioPHarris; // imagenes
	L_PTrioIm trioPDoG; // imagenes
	L_PTrioIm trioGradR, trioGradA; // modulo y angulo del gradiente de trioIm
	L_ImageGrayDouble pirGauss[L_PDS_NPIR];
	L_ImageGrayDouble pirLaplace[L_PDS_NPIR];
	L_ImageGrayDouble pirHarris[L_PDS_NPIR];
	L_ImageGrayDouble pirGradR[L_PDS_NPIR];
	L_ImageGrayDouble pirGradA[L_PDS_NPIR];
	double sigm[L_PDS_NPIR];
	L_Frac reltam[L_PDS_NPIR];
	int nPisos;

	double val_minP_DoG; // minElement value del punto maxElement en el SDoG que es aceptable a priori
	double val_minP_Harris; // minElement value del punto maxElement en Harris que es aceptable a priori
	double porcCurvMax;
	double factorS;
	double alfa;
	double s0_def;
	double s1_def;
	double s1_def_chico;
	int usa_parabol;
	long ultId;

	L_PTrioIm trioDoG;
	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros

	L_PHarris_SIFT();

	bool calcPiramides(L_ImageGrayDouble &im, double s0);
	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcDescriptores();

	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams() {return &paramPHARSIFT;}
	virtual ~L_PHarris_SIFT() {}
};


class L_PDoG_SURFG:public L_GenDescrPuntInt // detector DoG piramide scale space
{
public:
	L_ParamManagerLocal paramPDOGSURFG;
	L_PSS genPiram;
	L_PTrioIm trioPDoG; // imagenes
	L_PTrioIm trioGradR, trioGradA; // modulo y angulo del gradiente de trioIm
	L_ImageGrayDouble pirGauss[L_PDS_NPIR];
	L_ImageGrayDouble pirLaplace[L_PDS_NPIR];
	L_ImageGrayDouble pirGradR[L_PDS_NPIR];
	L_ImageGrayDouble pirGradA[L_PDS_NPIR];
	L_ImageGrayDouble pirGradDer[L_PDS_NPIR];
	L_ImageGrayDouble pirGradArr[L_PDS_NPIR];
	double sigm[L_PDS_NPIR];
	L_Frac reltam[L_PDS_NPIR];
	int nPisos;
	long ultId;

	double r_lin; // relacion (mayor value propio del hessiano)/(menor value propio del hessiano) para descartar lineas
	double val_minP_DoG;
	int usa_parabol;
	double porcCurvMax;
	double s0_def;
	double s1_def;
	double s1_def_chico;

	L_PTrioIm trioDoG;
	L_SURFGaussCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros

	L_PDoG_SURFG();

	bool calcPiramides(L_ImageGrayDouble &im, double s0);
	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcDescriptores();

	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams() {return &paramPDOGSURFG;}
	virtual ~L_PDoG_SURFG() {}
};


// CUIDADO: HAY QUE LLAMAR ESTAS MACROS PASANDOLES VARIABLES, no expresiones a menos que incluyan solo operaciones con constantes
#define L_SURF_macros
#ifdef L_SURF_macros
#define L_SURF_intRect(R,I,i1,j1,i2,j2) ((R) = (I).pix((i1)-1,(j1)-1)+(I).pix((i2),(j2))-(I).pix((i1)-1,(j2))-(I).pix((i2),(j1)-1))
#define L_SURF_intRectA(R,I,ies,i,j) (L_SURF_intRect(R,(I),(ies)[0]+(i),(ies)[1]+(j),(ies)[2]+(i),(ies)[3]+(j)))
#define L_SURF_intRectAT(R,I,ies,i,j) (L_SURF_intRect(R,(I),(ies)[1]+(i),(ies)[0]+(j),(ies)[3]+(i),(ies)[2]+(j)))
#define L_SURF_promRect(R,I,i1,j1,i2,j2) ((R) = (((I)).pix((i1),(j1))+((I)).pix((i2)+1,(j2)+1)-((I)).pix((i1),(j2)+1)-((I)).pix((i2)+1,(j1)))/((i2)+1-(i1))/((j2)+1-(j1)))
#define L_SURF_promRectA(R,I,ies,i,j) (L_SURF_promRect((R),(I),(ies)[0]+(i),(ies)[1]+(j),(ies)[2]+(i),(ies)[3]+(j)))
#define L_SURF_promRectAT(R,I,ies,i,j) (L_SURF_promRect((R),(I),(ies)[1]+(i),(ies)[0]+(j),(ies)[3]+(i),(ies)[2]+(j)))
#define L_SURF_areaRectA(R,I,ies) ((R) = ((ies)[2]+1-(ies)[0])*(long)((ies)[3]+1-(ies)[1]))
#define L_SURF_calcDxx(R,I,esc,i,j) {double uxx, vxx, wxx; L_SURF_intRectA(uxx,(I),Dxx_area1_p[(esc)],(i),(j)); L_SURF_intRectA(vxx,(I),Dxx_area2_n[(esc)],(i),(j)); L_SURF_intRectA(wxx,(I),Dxx_area3_p[(esc)],(i),(j)); (R)=uxx-2*vxx+wxx; }
#define L_SURF_calcDxy(R,I,esc,i,j) {double uxy, vxy, wxy, xxy; L_SURF_intRectA(uxy,(I),Dxy_area1_p[(esc)],(i),(j)); L_SURF_intRectA(vxy,(I),Dxy_area2_n[(esc)],(i),(j)); L_SURF_intRectA(wxy,(I),Dxy_area3_n[(esc)],(i),(j)); L_SURF_intRectA(xxy,(I),Dxy_area4_p[(esc)],(i),(j)); (R)=uxy-vxy-wxy+xxy; }
#define L_SURF_calcDyy(R,I,esc,i,j) {double uyy, vyy, wyy; L_SURF_intRectAT(uyy,(I),Dxx_area1_p[(esc)],(i),(j)); L_SURF_intRectAT(vyy,(I),Dxx_area2_n[(esc)],(i),(j)); L_SURF_intRectAT(wyy,(I),Dxx_area3_p[(esc)],(i),(j)); (R)=uyy-2*vyy+wyy; }
#define L_SURF_calcDetH(R,I,esc,i,j) { double Dxy38, dxx38, dyy38; L_SURF_calcDxy(Dxy38,(I), (esc),(i),(j)); Dxy38 = 0.9*Dxy38; L_SURF_calcDxx(dxx38,(I),(esc),(i),(j)); L_SURF_calcDyy(dyy38,(I),(esc),(i),(j)); (R) = (dxx38*dyy38-Dxy38*Dxy38)/tam4[(esc)];}
#define L_SURF_calcGrDer(R,I,esc,i,j) {double uD, vD; L_SURF_intRectA(uD,(I),GrX_area1_n[(esc)],(i),(j)); L_SURF_intRectA(vD,(I),GrX_area2_p[(esc)],(i),(j)); R = (-uD+vD)/tam2[(esc)];} // Gradiente en la direccion X
#define L_SURF_calcGrArr(R,I,esc,i,j) {double uA, vA; L_SURF_intRectAT(uA,(I),GrX_area1_n[(esc)],(i),(j)); L_SURF_intRectAT(vA,(I),GrX_area2_p[(esc)],(i),(j)); R = (uA-vA)/tam2[(esc)];} // Gradiente en la direccion -Y
#define L_SURF_calcGrDerV(R,I,lado,i,j) {double uDV, vDV; L_SURF_intRect(uDV,(I),(i)-(lado)/2, (j)-(lado)/2, (i), (j)+(lado)/2); L_SURF_intRect(vDV,(I),(i), (j)-(lado)/2, (i)+(lado)/2, (j)+(lado)/2); (R) = (-uDV+vDV)/(lado)/(lado);} // Gradiente en la direccion X
#define L_SURF_calcGrArrV(R,I,lado,i,j) {double uAV, vAV; L_SURF_intRect(uAV,(I),(i)-(lado)/2, (j)-(lado)/2, (i)+(lado)/2, (j)); L_SURF_intRect(vAV,(I),(i)-(lado)/2, (j), (i)+(lado)/2, (j)+(lado)/2); (R) = (uAV-vAV)/(lado)/(lado);} // Gradiente en la direccion -Y
#define L_SURF_calcGrDerV_VG(R,I,lado,i,j,VG) {double uDV, vDV; L_SURF_intRect(uDV,(I),(i)-(lado)/2, (j)-(lado)/2, (i), (j)+(lado)/2); L_SURF_intRect(vDV,(I),(i), (j)-(lado)/2, (i)+(lado)/2, (j)+(lado)/2); (R) = (VG)*(-uDV+vDV)/(lado)/(lado);} // Gradiente en la direccion X
#define L_SURF_calcGrArrV_VG(R,I,lado,i,j,VG) {double uAV, vAV; L_SURF_intRect(uAV,(I),(i)-(lado)/2, (j)-(lado)/2, (i)+(lado)/2, (j)); L_SURF_intRect(vAV,(I),(i)-(lado)/2, (j), (i)+(lado)/2, (j)+(lado)/2); (R) = (VG)*(uAV-vAV)/(lado)/(lado);} // Gradiente en la direccion -Y
#define L_SURF_calcGrDerV_VG(R,I,lado,i,j,VG) {double uDV, vDV; L_SURF_intRect(uDV,(I),(i)-(lado)/2, (j)-(lado)/2, (i), (j)+(lado)/2); L_SURF_intRect(vDV,(I),(i), (j)-(lado)/2, (i)+(lado)/2, (j)+(lado)/2); (R) = (VG)*(-uDV+vDV)/(lado)/(lado);} // Gradiente en la direccion X
#define L_SURF_calcGrArrV_VG(R,I,lado,i,j,VG) {double uAV, vAV; L_SURF_intRect(uAV,(I),(i)-(lado)/2, (j)-(lado)/2, (i)+(lado)/2, (j)); L_SURF_intRect(vAV,(I),(i)-(lado)/2, (j), (i)+(lado)/2, (j)+(lado)/2); (R) = (VG)*(uAV-vAV)/(lado)/(lado);} // Gradiente en la direccion -Y
#else
#define L_SURF_intRect(R,I,i1,j1,i2,j2) ((R) = intRect(i1,j1,i2,j2))
#define L_SURF_intRectA(R,I,ies,i,j) ((R) = intRectA(ies,i2,j2))
#define L_SURF_intRectAT(R,I,ies,i,j) ((R) = intRectAT(ies,i2,j2))
#define L_SURF_promRect(R,I,i1,j1,i2,j2) ((R) = promRect(i1,j1,i2,j2))
#define L_SURF_promRectA(R,I,ies,i,j) ((R) = promRectA(ies,i2,j2))
#define L_SURF_promRectAT(R,I,ies,i,j) ((R) = promRectAT(ies,i2,j2))
#define L_SURF_areaRectA(R,I,ies) ((R) = areaRectA(ies))
#define L_SURF_calcDxx(R,I,esc,i,j) ((R)=calcDxx(esc,i,j))
#define L_SURF_calcDxy(R,I,esc,i,j) ((R)=calcDxy(esc,i,j))
#define L_SURF_calcDyy(R,I,esc,i,j) ((R)=calcDyy(esc,i,j))
#define L_SURF_calcDetH(R,I,esc,i,j) ((R)=calcDetH(esc,i,j))
#define L_SURF_calcGrDer(R,I,esc,i,j) ((R)=calcGrDer(esc,i,j))
#define L_SURF_calcGrArr(R,I,esc,i,j) ((R)=calcGrArr(esc,i,j))
#define L_SURF_calcGrDerV(R,I,lado,i,j) ((R)=calcGrDerV(esc,i,j))
#define L_SURF_calcGrArrV(R,I,lado,i,j) ((R)=calcGrArrV(esc,i,j))
#endif

class L_SURF_POD
{
public:
	int tam[L_PDS_NPIR];
	double tam2[L_PDS_NPIR];
	double tam4[L_PDS_NPIR];
	double sigma[L_PDS_NPIR];
	int l_gaus2_5;
	L_Frac porcSubm01, porcSubm12;
	int tam_0;
	int deltaTam_0;
	bool Dcalculado;
	int nPisos;
	int Dxx_area1_p[L_PDS_NPIR][4];
	int Dxx_area2_n[L_PDS_NPIR][4];
	int Dxx_area3_p[L_PDS_NPIR][4];
	int Dxy_area1_p[L_PDS_NPIR][4];
	int Dxy_area2_n[L_PDS_NPIR][4];
	int Dxy_area3_n[L_PDS_NPIR][4];
	int Dxy_area4_p[L_PDS_NPIR][4];
	int GrX_area1_n[L_PDS_NPIR][4];
	int GrX_area2_p[L_PDS_NPIR][4];
	double calcOri_mFinal_cuadr1[6];
	double calcOri_mFinal_cuadr2[6];
	double calcOri_sumaX_celdas[24];
	double calcOri_sumaY_celdas[24];
	bool calcOri_precalc;
	// Parametros
	double resMaxMinSigmas;
	double val_minP_HoG;
	double porcCurvMax;
	bool dupl;
	int nrx;
	int lx;
	long ultId;

	void swap(L_SURF_POD &other) {L_SURF_POD t = other; other = *this; *this = t;}
};

// Los tamanos de los filtros cambian de forma distinta a la indicada en el paper
// Aca hay mas tamanos intermedios
class L_SURF:public L_GenDescrPuntInt, public L_SURF_POD
{
public:
	L_ParamManagerLocal paramSURF;
	L_ImageGrayDouble imIntegral;
	L_ImageGrayDouble original; // Se usa solo durante tests de debug
	L_ImageGrayDouble pirHess[L_PDS_NPIR];
	L_ImageGrayDouble pirGradR[L_PDS_NPIR];
	L_ImageGrayDouble pirGradA[L_PDS_NPIR];
	L_PTrioIm trio;
	L_Descriptor desSURF;
	std::vector<double> gaus2_5;
	std::vector<double> gaus3_3;

	L_SURF();

	void calcD_pn();
	double intRect(int i1, int j1, int i2, int j2) {return imIntegral.pix(i1-1,j1-1)+imIntegral.pix(i2,j2)-imIntegral.pix(i1-1,j2)-imIntegral.pix(i2,j1-1);}
	double intRectA(int ies[4], int i, int j) {return intRect(ies[0]+i,ies[1]+j,ies[2]+i,ies[3]+j);}
	double intRectAT(int ies[4], int i, int j) {return intRect(ies[1]+i,ies[0]+j,ies[3]+i,ies[2]+j);}
	double promRect(int i1, int j1, int i2, int j2) {return (imIntegral.pix(i1,j1)+imIntegral.pix(i2+1,j2+1)-imIntegral.pix(i1,j2+1)-imIntegral.pix(i2+1,j1))/(i2+1-i1)/(j2+1-j1);}
	double promRectA(int ies[4], int i, int j) {return promRect(ies[0]+i,ies[1]+j,ies[2]+i,ies[3]+j);}
	double promRectAT(int ies[4], int i, int j) {return promRect(ies[1]+i,ies[0]+j,ies[3]+i,ies[2]+j);}
	static long areaRectA(int ies[4]) {return (ies[2]+1-ies[0])*(long)(ies[3]+1-ies[1]);}
	static bool revisaAreasIguales(int n, ...) {bool b; va_list ppp; va_start(ppp,n); b=vrevisaAreasIguales(n, ppp); va_end(ppp); return b;}
	static bool vrevisaAreasIguales(int n, va_list ppp);

	long genDxx(int esc, double x1, double y1, double x2, double y2);
	long genDxy(int esc, double x1, double y1, double x2, double y2);
	long genGrX(int esc, double x1, double y1, double x2, double y2);

	double calcDxx(int esc, int i, int j) { return intRectA(Dxx_area1_p[esc],i,j)-2*intRectA(Dxx_area2_n[esc],i,j)+intRectA(Dxx_area3_p[esc],i,j); }
	double calcDxy(int esc, int i, int j) { return intRectA(Dxy_area1_p[esc],i,j)-intRectA(Dxy_area2_n[esc],i,j)-intRectA(Dxy_area3_n[esc],i,j)+intRectA(Dxy_area4_p[esc],i,j); }
	double calcDyy(int esc, int i, int j) { return intRectAT(Dxx_area1_p[esc],i,j)-2*intRectAT(Dxx_area2_n[esc],i,j)+intRectAT(Dxx_area3_p[esc],i,j); }
	double calcDetH(int esc, int i, int j) { double Dxy38=0.9*calcDxy(esc,i,j); return (calcDxx(esc,i,j)*calcDyy(esc,i,j)-Dxy38*Dxy38)/tam4[esc];}
	double calcGrDer(int esc, int i, int j) {return (-intRectA(GrX_area1_n[esc],i,j)+intRectA(GrX_area2_p[esc],i,j))/tam2[esc];} // Gradiente en la direccion X
	double calcGrArr(int esc, int i, int j) {return (intRectAT(GrX_area1_n[esc],i,j)-intRectAT(GrX_area2_p[esc],i,j))/tam2[esc];} // Gradiente en la direccion -Y
	double calcGrDerV(int lado, int i, int j) {return (-intRect(i-lado/2, j-lado/2, i, j+lado/2)+intRect(i, j-lado/2, i+lado/2, j+lado/2))/lado/lado;} // Gradiente en la direccion X
	double calcGrArrV(int lado, int i, int j) {return (intRect(i-lado/2, j-lado/2, i+lado/2, j)-intRect(i-lado/2, j, i+lado/2, j+lado/2))/lado/lado;} // Gradiente en la direccion -Y

	double calcOrientacion(int esc, double i, double j, double sigma0);

	bool calcPiramides();
	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcVectorDescriptor(double *vector, int esc, double i, double j, double sigma0, double ang, int &nCompVect);
	bool calcDescriptores();
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams() {return &paramSURF;}

	void swap(L_SURF &other);

	~L_SURF() {}
};


class L_Harris:public L_GenDescrPuntInt
{
public:
	double val_minP_Harris; // minElement value del punto maxElement en Harris que es aceptable a priori
	double factorS;
	double alfa;
	int nobj;
	long ultId;

	L_ImageGrayDouble imHarris;
	L_ImageGrayDouble gradR, gradAng;

	L_PuntIntLista punL;
	L_SIFTdireccLista sdirL;
	L_DescriptorLista desL;
	
	L_Harris() {tipoP = HARRIS; tipoD = DESNADA; val_minP_Harris = 1e-6; factorS=1.5; alfa=0.04; nobj=-1; ultId = 0;}

	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);

};

class L_GenDescrExterno : public L_GenDescrPuntInt
{
public:
	L_String _lineaComando;
	L_String _nomArchImagen;
	L_String _nomArchDescr;
	bool formMiko;
	bool _BMP;
	bool _PGM;
	long ultId;

private:
	L_GenDescrExterno() {tipoP = P_EXT; tipoD = D_EXT; _BMP=false; _PGM=false; ultId = 0;} // No debe llamarse

public:
	L_GenDescrExterno(const char *lineaComando, const char *nomArchImagen, const char *nomArchDescr);
	L_GenDescrExterno(const char *nomArchImagen);
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0);
	~L_GenDescrExterno()
	{
	}
};


class L_MinuciasExterno:public L_GenDescrPuntInt
{
public:
	L_ParamManagerLocal paramMinExt;
	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros
	bool agregaEspejo;
	L_String lineaComando;     // Linea de comando para llamar al calculador de minucias, si es "" no se llama
	L_String nomArchImagenTmp; // name de archivo de imagen temporal de transferencia para llamar a la linea de comandos
	L_String nomArchMinucias;  // name del archivo de texto de donde se van a read las minucias
	bool formMiko;  // No implementado
	bool usaBMP; // Usar bmp para archivo de imagen temporal de transferencia
	bool usaPGM; // Usar pgm para archivo de imagen temporal de transferencia
	bool usaPiramide;
	bool muestraProceso;  // Para debug
	bool manejarArchivos; // Indica si se debe preocupar de la lectura de los archivos
	long ultId;

	// Imagenes derivadas de la deteccion de huellas (amplitud y frecuencia, opcionales)
	L_ImageGrayDouble imFrec;
	L_ImageGrayDouble imOri;

private:	
	L_MinuciasExterno():paramMinExt("MINUC_EXT",2) { L_hard_shutdown("Constructor no valido"); } // No debe existir
public:
	L_MinuciasExterno(const char *linComando, const char *nomArchivoImagen, const char *nomArchivoMinucias); // Linea de comandos
	L_MinuciasExterno(const char *nomArchivoMinucias); // En este caso se read directamente el archivo de minucias de una base de datos

	void fijaArchivoMinucias(const char *nomArchivoMinucias);

	L_ParamManagerLocal *pideParams() {return &paramMinExt;}
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0);

	bool buscaCamposOrientacionFrecuencia(); // busca el archivonomArchImagenTmp + "_OF"
	void grabaCamposOrientacionFrecuencia(FILE *fp, bool compr = true); // Revisa por si mismo si los campos estan definidos
	void leeCamposOrientacionFrecuencia(FILE *fp); // Revisa por si mismo si los campos estan definidos
	bool comparaCamposOrientacionFrec(double &difOri, double &difFrec, double &corrFrec, L_TransfAfinPI2D &tr, L_MinuciasExterno &other, double frecMin, double frecMax, bool usa_v, int *vx, int *vy); // 4 vertices, re=corrAng, im=corrFrec, tr transformUsing (1) en (2)

	void grabaFlechasMinucias(const char *nomArch, L_ImageGrayDouble &im);
};

class L_Harris_SIFT:public L_GenDescrPuntInt // detector Harris no invariante a escala
{
public:
	L_ParamManagerLocal paramH1SIFT;
	double val_min_Harris; // minElement value del punto maxElement en Harris que es aceptable a priori
	double porcCurvMax;
	double alfa;
	double sigma_i;
	double sigma_d;
	double sigma_i_dupl;
	double sigma_d_dupl;
	bool dupl;
	double sMin;
	double sigma0;
	long ultId;

	L_SIFTCalculador sgen; // almacenamiento temporal + parametros
	L_SIFTdireccLista sdirL; // almacenamiento temporal + parametros
	L_ImageGrayDouble imInt;
	L_ImageGrayDouble imGradR, imGradA, imHarris;

	L_Harris_SIFT();

	bool calcRegiones(L_ImageGrayDouble &im);
	bool calcDirecciones();
	bool calcDescriptores();
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	L_ParamManagerLocal *pideParams() {return &paramH1SIFT;}
	virtual ~L_Harris_SIFT() {}
};

class L_MSER:public L_GenDescrPuntInt
{
public:
	double deltaLum;
	long ultId;
	L_MSER() {nobj=-1; deltaLum=0.033; ultId = 0;}
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5);
	void calcRegionesMser(double valInic, double delta);
	virtual ~L_MSER() {}
};

// Para poner varios generadores de descriptores en paralelo
class L_GenDescrArr:public L_GenDescrPuntInt
{
private:
	L_ParamManagerLocal paramUnion; // Solo contiene referencias a los parametros de los children
	L_GenDescrArr(); //:paramUnion("GENDESCRARR_NULO",2); { L_hard_shutdown("Constructor no valido"); } // No debe existir
	L_GenDescrArr(const L_GenDescrArr &other); //:paramUnion("GENDESCRARR_NULO",2) { L_hard_shutdown("Constructor no valido"); } // No debe existir
	L_GenDescrArr& operator=(const L_GenDescrArr &other); // {L_hard_shutdown("Constructor no valido"); return *this;} // No debe existir
public:
	L_GenDescrPuntInt **arr; // sub-generadores de descriptores, colocados en paralelo
	bool *activo;
	int nMax;
	int n;
	L_GenDescrArr(int nMax):paramUnion("GENDESCRARR_NULO",nMax+2) {tipoP = P_UNION; tipoD = D_UNION; this->nMax=nMax; n=0; arr=new L_GenDescrPuntInt*[nMax]; activo=new bool[nMax];}
	bool agregaReferencia(L_GenDescrPuntInt *gpi);
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5); //!< Pasa la tarea de generar los descriptores a sus sub-generadores
	L_ParamManagerLocal *pideParams() {return &paramUnion;}
	~L_GenDescrArr() {delete[] arr; delete activo;}
};



class L_EvalGPIInfo
{
public:
	double rot;
	double esc;
	double cotaRuido;
	double radioPixeles;
	double repetibilidadPuntos;
	double repetibilidadDescriptores;
};

typedef L_Node <L_EvalGPIInfo> L_EvalGPIInfoNodo;
typedef L_List <L_EvalGPIInfo> L_EvalGPIInfoLista;

class L_EvalGenDescrPuntInt
{
public:
	L_GenDescrPuntInt *gpi;
	L_ImageGrayDouble imOrig;
	L_EvalGPIInfoLista infL;
	L_ParamManagerLocal paramEvalGPI;
//parametros
	double rotMin;       // rotacion inicial de las pruebas (en radianes)
	double rotMax;       // rotacion final de las pruebas (en radianes)
	double rotDelta;     // cambio de rotacion entre pruebas (en radianes)
	double escMin;       // escala inicial (en octavas)
	double escMax;       // escala final (en octavas)
	double escDelta;     // cambio de escala entre pruebas (en octavas)
	double cotaRuido;    // Cota para ruido uniforme sumado a las imagenes, entre 0 y 1.
	double radioPixeles; // radio de los puntos de interes (en sigmas)

	L_EvalGenDescrPuntInt(L_GenDescrPuntInt *genPuntInt):paramEvalGPI("EVALGPI",1)
	{
		gpi=genPuntInt;
		rotMin=0; rotMax=2*M_PI; rotDelta=2*M_PI/12; escMin=-2; escMax=2; escDelta=0.25; cotaRuido=0; radioPixeles=5;
		paramEvalGPI.addFrom("rotMin",&rotMin);
		paramEvalGPI.addFrom("rotMax",&rotMax);
		paramEvalGPI.addFrom("rotDelta",&rotDelta);
		paramEvalGPI.addFrom("escMin",&escMin);
		paramEvalGPI.addFrom("escMax",&escMax);
		paramEvalGPI.addFrom("escDelta",&escDelta);
		paramEvalGPI.addFrom("cotaRuido",&cotaRuido);
		paramEvalGPI.addFrom("radioPixeles",&radioPixeles);
		paramEvalGPI.addChildren(genPuntInt->pideParams());
	}

	void evalGPI();
};


#endif // __L_GENPUNTINT_H__
