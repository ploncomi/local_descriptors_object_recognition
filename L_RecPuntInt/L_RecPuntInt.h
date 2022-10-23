#ifndef __L_RECPUNTINT_H__
#define __L_RECPUNTINT_H__
#include "L_GenPuntInt.h"
//#include "../L_LoweSIFT/L_Lowe.h"


//#include "../OpenSURF/L_OpenSURF.h"
#include <time.h>

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



////////////////////////////////////
/////
///// Clases de nivel height para define reconocedores que usen los detectores,
/////    verificadores, etc. como bloques.
/////
//////////////////////////////


class L_TransfSimil2DCelda_POD
{
public:
	int i,j,k,z,nobj;
	void swap(L_TransfSimil2DCelda_POD &other) {L_TransfSimil2DCelda_POD t = other; other = *this; *this = t;}
};

class L_TransfSimil2DCelda : public L_TransfSimil2DCelda_POD
{
public:
	L_CalceLista calL;
	L_TransfSimil2DCelda() {};
	~L_TransfSimil2DCelda() {};
	void swap(L_TransfSimil2DCelda &other);
};

typedef L_Node<L_TransfSimil2DCelda> L_TransfSimil2DCeldaNodo;
typedef L_List<L_TransfSimil2DCelda> L_TransfSimil2DCeldaLista;

class L_TransfSimil2DCeldaPtr
{
public:
	L_TransfSimil2DCelda *c;

	// Los ordna del menos penca al mas penca
	bool operator < (const L_TransfSimil2DCeldaPtr &other) const {return c->calL.size() > other.c->calL.size();}
};

class L_TransfSimil2DCeldaHash_POD
{
public:
	int nTabla;
	int nTabla_param;
	double dxBin;
	double dyBin;
	double dangBin;
	double doctBin;
#ifdef __TURBOC__
	__OP_ASSIGN(L_TransfSimil2DCeldaHash_POD);
#endif
	void swap(L_TransfSimil2DCeldaHash_POD &other) {L_TransfSimil2DCeldaHash_POD t = other; other = *this; *this = t;}
};

class L_TransfSimil2DCeldaHash : public L_TransfSimil2DCeldaHash_POD
{
public:
	L_ParamManagerLocal paramTSCH;
	L_TransfSimil2DCeldaLista *tabla;
	L_CountingTree votObjs;

	L_TransfSimil2DCeldaHash():paramTSCH("TSCH",5)
	{
		tabla=NULL;
		nTabla_param=8000;
		dxBin=0.25;  // Celda de 1/4 del width de la imagen original
		dyBin=0.25;  // Celda de 1/4 del height de la imagen original
		dangBin=30*M_PI/180;  // Celda de 30 grados de width angular
		doctBin=1;  // Celda de 1 octava de width
		paramTSCH.addFrom("nTabla_param",&nTabla_param);
		paramTSCH.addFrom("dxBin",&dxBin);
		paramTSCH.addFrom("dyBin",&dyBin);
		paramTSCH.addFrom("dangBin",&dangBin);
		paramTSCH.addFrom("doctBin",&doctBin);
	}
	void voto(L_Calce& cal, L_List<L_Calce> *cal_pool = NULL);
	void minivoto(int i, int j, int k, int z, int nobj, L_Calce& cal, L_List<L_Calce> *cal_pool = NULL);
	int getVotos(int i, int j, int k, int z, int nobj) {L_TransfSimil2DCeldaNodo *bin; buscaCelda(i,j,k,z,nobj,&bin); if (bin==NULL) return 0; else return (int)bin->c.calL.size();}
	L_TransfSimil2DCeldaLista *buscaCelda(int i, int j, int k, int z, int nobj, L_TransfSimil2DCeldaNodo **bin);
	void destroy() {if (tabla!=NULL) {delete[] tabla; tabla=NULL;}}
	~L_TransfSimil2DCeldaHash() {destroy();}

	void swap(L_TransfSimil2DCeldaHash &other);
};


#if defined(L_DEFT_NFIL) || defined(L_DEFT_NCAR) || defined(L_DEFT_NPUN)
	#error #define conflicts
#endif
#define L_DEFT_NFIL (10) // Numero de filtros soportados
#define L_DEFT_NCAR (15) // Numero de caracteres por name de filtro
#define L_DEFT_NPUN (10)  // Numero de mejores candidatos eliminados a los que se les guarda el puntaje
enum L_Debug_filtrosTransf_criterios {nVotos, probTransf, probCelda};

class L_Debug_filtrosTransf
{
public:
	double mejoresPunt[L_DEFT_NPUN];
	double mejoresCaract[L_DEFT_NPUN];
	int mejoresFiltros[L_DEFT_NPUN];
	double mejoresRazon[L_DEFT_NPUN];
	int nTransfTot;
	int nTransf;
	int nFiltros;
	bool activo;
	L_Debug_filtrosTransf_criterios selPunt;
	char filtros[L_DEFT_NFIL][L_DEFT_NCAR]; // L_DEFT_NFIL filtros de L_DEFT_NCAR caracteres

	L_Debug_filtrosTransf()
	{
		activo=false;
		nTransf=0;
		nTransfTot=0;
		nFiltros=0;
		selPunt=nVotos;
	}
	void resetear() {nTransf=0; nTransfTot=0;}
	void fijaCriterio(L_Debug_filtrosTransf_criterios crit) {selPunt=crit;}
	void agregaFiltro(const char *name);
	void marcaEliminacion(L_TransfPuntInt2D *candTr, L_TransfSimil2DCelda*candCel, const char *name, double razon) {if (activo) marcaEliminacion_(candTr, candCel, name, razon);}
	void fusiona(const L_Debug_filtrosTransf &d1, const L_Debug_filtrosTransf &d2);
	void imprimeReporte();
	void swap(L_Debug_filtrosTransf &other) {L_Debug_filtrosTransf t = other; other = *this; *this = t;}
private:
	void marcaEliminacion_(L_TransfPuntInt2D *candTr, L_TransfSimil2DCelda *candCel, const char *name, double razon);
	void insertaElimDebug(int filtro, double puntaje, double caract, double razon);
};



class L_ImagenBNxfloat_ext
{
public:
	L_ImageGrayDouble im;
	L_ImageGrayDouble *ext; // Si esta definido, se apunta a una imagen externa
	L_ImagenBNxfloat_ext() : ext(NULL) { }
	void swap(L_ImagenBNxfloat_ext &other) {L_ImageGrayDouble *ext_t = other.ext;  im.swap(other.im); other.ext = ext; ext = ext_t;}
	void destroy() {im.destroy(); ext=NULL;}
};

typedef L_Node<L_ImagenBNxfloat_ext> L_ImagenBNNodo;
class L_ImagenBNLista : public L_List<L_ImagenBNxfloat_ext>
{
public:
	L_ImageGrayDouble *insertaRobaImagen(L_ImageGrayDouble &im) {L_ImagenBNxfloat_ext n; im.swap(n.im); return &push_back_swapping_ret(n)->c.im;}
	L_ImageGrayDouble *insertaPtrImagen(L_ImageGrayDouble &im) {L_ImagenBNxfloat_ext n; n.ext = &im; return push_back_swapping_ret(n)->c.ext;}
};
typedef L_Node<L_ImageRGBDouble> L_ImagenRGBNodo;
typedef L_List<L_ImageRGBDouble> L_ImagenRGBLista;

class L_VerifPixelsTransfBN_POD // Componentes de L_VerifPixelsTransfBN que se pueden copiar directamente
{
public:
	L_Debug_filtrosTransf debFiltr;
	bool verif;
	double umbCorr;
	L_ImageGrayDouble *ultIm;
#ifdef __TURBOC__
	__OP_ASSIGN(L_VerifPixelsTransfBN_POD);
#endif
	void swap(L_VerifPixelsTransfBN_POD &other) {L_VerifPixelsTransfBN_POD t = other; other = *this; *this = t;}
};

class L_VerifPixelsTransfBN : public L_VerifPixelsTransfBN_POD
{
public:
	L_ParamManagerLocal paramCPBN;
	L_ImagenBNLista imL; // Almacenando las imagenes aca, uno se asegura que no se intente acceder a una imagen que fue eliminada
	L_VerifPixelsTransfBN():paramCPBN("CPBN",1)
	{
		verif=true; umbCorr=0.4; ultIm=NULL;
		paramCPBN.addFrom("umbCorr",&umbCorr);
		paramCPBN.addFrom("verif",&verif);
		paramCPBN.addFrom("debFiltr.activo",&debFiltr.activo);
		debFiltr.agregaFiltro("corrPix");
	}
	bool agregaImagenEntrenRobaObjetoDe(L_ImageGrayDouble &imRef); // Se queda con la imagen que se le paso
	bool agregaImagenEntrenCopiaObjetoDe(L_ImageGrayDouble &imRef); // Se queda con una copia de la imagen que se le paso
	bool agregaImagenEntrenSinCopia(L_ImageGrayDouble &imRef); // No se queda con copia de la imagen que se le paso
	double calcCorrPix(const L_ImageGrayDouble &imPru, L_TransfPuntInt2D &tr);
	bool verifTransf(const L_ImageGrayDouble &imPru, L_TransfPuntInt2D &tr, double *val = NULL);
	void verifListaTransf(L_TransfPuntInt2DListaPtr &trL);
	void liberaPixeles();

	void swap(L_VerifPixelsTransfBN &other) {other.L_VerifPixelsTransfBN_POD::swap(*this); paramCPBN.swap(other.paramCPBN); imL.swap(other.imL);}
};

class L_VerifPixelsTransfRGB
{
public:
	L_ParamManagerLocal paramCPRGB;
	bool verif;
	double umbCorr;
	L_ImagenRGBLista imL; // Almacenando las imagenes aca, uno se asegura que no se intente acceder a una imagen que fue eliminada
	L_ImageRGBDouble *ultIm;
	L_VerifPixelsTransfRGB():paramCPRGB("CPRGB",1)
	{
		verif=true; umbCorr=0.4; ultIm=NULL;
		paramCPRGB.addFrom("umbCorr",&umbCorr);
		paramCPRGB.addFrom("verif",&verif);
	}
	bool agregaImagenEntrenRobaObjetoDe(L_ImageRGBDouble &imRef);
	bool agregaImagenEntrenCopiaObjetoDe(L_ImageRGBDouble &imRef);
	bool verifTransf(L_ImageRGBDouble &imPru, L_TransfPuntInt2D &tr);
	void verifListaTransf(L_ImageRGBDouble &imPru, L_TransfPuntInt2DListaPtr &trL);
	void liberaPixeles();
};


// class_name para hacer debug




enum L_GenTransfCalces_tipo {gtc_indef, genTrAfinHough, genTrProyRansac_1tr, genTrCamParal_1tr, genTrAfinEpipRansac};

class L_GenTransfCalces_POD
{
public:
	L_GenTransfCalces_tipo tipo;
	L_Debug_filtrosTransf debFiltr;
	struct strCorrLinDatos {bool verif; double rMin;} corrLinDatos; // verifica correlacion lineal de los datos de input (debe ser baja)
	void swap(L_GenTransfCalces_POD &other) {L_GenTransfCalces_POD t = other; other = *this; *this = t;}
};

class L_GenTransfCalces : public L_GenTransfCalces_POD
{
private:
	L_GenTransfCalces() {tipo = gtc_indef; corrLinDatos.verif=true; corrLinDatos.rMin=0.9;}
public:
	L_GenTransfCalces(L_GenTransfCalces_tipo enumTipo) {tipo = enumTipo; corrLinDatos.verif=true; corrLinDatos.rMin=0.9;}
	virtual bool genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool = NULL)=0; // No debe modificar calL
	virtual L_ParamManagerLocal *pideParams() {return NULL;}

	bool verifCorrLinDatos(L_CalceLista &calL, double *corrPreg=NULL) {return calL.verifCorrLinDatos(corrLinDatos.rMin, corrPreg);}
};


class L_GenTransfAfinHough_POD
{
public:
	bool topDown; // Intentar agregar calces de calL a la nueva transformación
	bool maxVot4D; // La celda (i,j,k,z) debe tener votacion maxima localmente
	int nVotMinCelda; // = 4
	int nMaxCeldasRev; // = -1
	bool mostrarTiemposInterno;

	struct strProbCelda {bool verif; bool calce16votos; double pUmb;} probCelda; // probabilidad antes de calcular transformacion
	struct strProbTransf {bool verif; bool corrPz; double pUmb;} probTransf; // probabilidad despues de calcular transformacion
	struct strDistorRomb {bool verif; double senAng; double lados;} distorRomb; // verifica distorsion afin pequeña
	struct strFusion {bool activo; bool separarRansac; double superp;} fusion; // fusion de transformaciones parecidas
	struct strRansac {bool activo; int nInt; double dx; double dy; double votMin; double porcVotMin; double porcVotInf; double factorExp; double coef;} ransac;
#ifdef __TURBOC__
        __OP_ASSIGN(L_GenTransfAfinHough_POD);
#endif
	void swap(L_GenTransfAfinHough_POD &other) {L_GenTransfAfinHough_POD t = other; other = *this; *this = t;}
};

class L_GenTransfAfinHough:public L_GenTransfCalces, public L_GenTransfAfinHough_POD
{
public:
	L_ParamManagerLocal paramGTAH;
	L_TransfSimil2DCeldaHash hash;

	L_GenTransfAfinHough():L_GenTransfCalces(genTrAfinHough),paramGTAH("GTAH",5)
	{
		topDown=true;
		maxVot4D=true;
		nVotMinCelda = 4;
		nMaxCeldasRev = -1;
		mostrarTiemposInterno = false;
		probCelda.verif=true; probCelda.calce16votos=true; probCelda.pUmb=0.85;
		probTransf.verif=true; probTransf.corrPz=false; probTransf.pUmb=0.95;
		distorRomb.verif=true; distorRomb.senAng=sin(55*M_PI/180); distorRomb.lados=2;
		fusion.activo=true; fusion.separarRansac = true; fusion.superp=0.85; //fusion.difang=30*M_PI/180; 
		//ransac
		ransac.activo=true; ransac.dx=0.125; ransac.dy=0.125; ransac.votMin=10; ransac.nInt=8;
		ransac.porcVotMin=0.625; ransac.porcVotInf=0.25; ransac.factorExp=0.944;

		debFiltr.agregaFiltro("probCelda");
		debFiltr.agregaFiltro("probTransf");
		debFiltr.agregaFiltro("distorRomb");
		debFiltr.agregaFiltro("corrLinPuntos");
		debFiltr.agregaFiltro("fusion_RANSAC");

		paramGTAH.addFrom("topDown",&topDown);
		paramGTAH.addFrom("maxVot4D",&maxVot4D);
		paramGTAH.addFrom("probCelda.verif",&probCelda.verif);
		paramGTAH.addFrom("probCelda.pUmb",&probCelda.pUmb);
		paramGTAH.addFrom("probCelda.calce16votos",&probCelda.calce16votos);
		paramGTAH.addFrom("distorRomb.verif",&distorRomb.verif);
		paramGTAH.addFrom("distorRomb.senAng",&distorRomb.senAng);
		paramGTAH.addFrom("distorRomb.lados",&distorRomb.lados);
		paramGTAH.addFrom("corrLinDatos.verif",&corrLinDatos.verif);
		paramGTAH.addFrom("corrLinDatos.rMin",&corrLinDatos.rMin);
		paramGTAH.addFrom("probTransf.verif",&probTransf.verif);
		paramGTAH.addFrom("probTransf.corrPz",&probTransf.corrPz);
		paramGTAH.addFrom("probTransf.pUmb",&probTransf.pUmb);
		paramGTAH.addFrom("fusion.activo",&fusion.activo);
		paramGTAH.addFrom("fusion.separarRansac",&fusion.separarRansac);
		paramGTAH.addFrom("fusion.superp",&fusion.superp);
		paramGTAH.addFrom("ransac.activo",&ransac.activo);
		paramGTAH.addFrom("ransac.dx",&ransac.dx);
		paramGTAH.addFrom("ransac.dy",&ransac.dy);
		paramGTAH.addFrom("ransac.votMin",&ransac.votMin);
		paramGTAH.addFrom("ransac.nInt",&ransac.nInt);
		paramGTAH.addFrom("ransac.porcVotMin",&ransac.porcVotMin);
		paramGTAH.addFrom("ransac.porcVotInf",&ransac.porcVotInf);
		paramGTAH.addFrom("ransac.factorExp", &ransac.factorExp);
		paramGTAH.addFrom("debFiltr.activo", &debFiltr.activo);

		paramGTAH.addChildren(&hash.paramTSCH);

		//coef ransac, para acelerar calculos de la funcion de cantidad de votos
		ransacCoefRecalcular();
	}

	virtual ~L_GenTransfAfinHough() {}

	inline void ransacCoefRecalcular() {ransac.coef=(ransac.porcVotMin-ransac.porcVotInf)/pow(ransac.factorExp,ransac.votMin);}

	bool genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool);

	L_ParamManagerLocal *pideParams();

	bool verifDistorRomb(L_TransfAfinPI2D &tr, double *val=NULL);
	bool fusiona(L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool=NULL);
	bool verifFusion(L_TransfAfinPI2D &tr1, L_TransfAfinPI2D &tr2);
	int fnNumRansac(int nvotos);
	void agregaTopDown(L_TransfAfinPI2D &tr, L_CalceLista &calL, L_List<L_Calce> *cal_pool=NULL);
	double calcProbCelda(L_TransfSimil2DCelda &bin, L_CalceLista &calL, int nVotCelda=-1); // calculo de probabilidad de bin (mío)
	double calcProbTransf(L_TransfAfinPI2D &tr, L_CalceLista &calL); // calculo de probabilidad de transf (Lowe)

	void swap(L_GenTransfAfinHough &other);
};

class L_GenTransfProyectivaRansac_1tr:public L_GenTransfCalces // Solo busca 1 transformacion
{
public:
	L_ParamManagerLocal paramGTPR_1;
	L_CamaraDistorsionRadial cam;
	double radioError;
	bool usar_4calces;
	bool considerarCam;
	int minCalces;
	L_GenTransfProyectivaRansac_1tr();
	virtual bool genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool=NULL);
	virtual L_ParamManagerLocal *pideParams() {return &paramGTPR_1;}
};

class L_GenTransfCamParalelas_1tr:public L_GenTransfCalces
{
public:
	L_ParamManagerLocal paramGTCP_1;
	double errorVertAnt;
	double errorVertPost;

	L_GenTransfCamParalelas_1tr();
	virtual bool genTransf(L_CalceLista &calL, L_TransfPuntInt2DListaPtr &traL, L_List<L_Calce> *cal_pool=NULL);
	virtual L_ParamManagerLocal *pideParams() {return &paramGTCP_1;}
};


class L_RecTimeSeg
{
public:
	L_PrintingLanguage idioma; // L_espanol o L_english
	double t1, t2;
	L_Statistics1D gpiTime;
	L_Statistics1D gbdTime;
	L_Statistics1D calTime;
	L_Statistics1D gtcTime;
	L_RecTimeSeg() {idioma=L_compilationLanguage();}
	void escribeResumen(FILE *fp = stdout)
	{
		if (idioma==L_espanol)
		{
			fprintf(fp, "generador descrip : %.4f[s] +- %.4f[s] #%d\n", gpiTime.getMean(), sqrt(gpiTime.getVariance()), gpiTime.size());
			fprintf(fp, "generador base dat: %.4f[s] +- %.4f[s] #%d\n", gbdTime.getMean(), sqrt(gbdTime.getVariance()), gbdTime.size());
			fprintf(fp, "generador calces  : %.4f[s] +- %.4f[s] #%d\n", calTime.getMean(), sqrt(calTime.getVariance()), calTime.size());
			fprintf(fp, "generador transfor: %.4f[s] +- %.4f[s] #%d\n", gtcTime.getMean(), sqrt(gtcTime.getVariance()), gtcTime.size());
			fprintf(fp, "frame rate total  : %.1f[fps]\n", 1/(gpiTime.getMean()+gbdTime.getMean()+calTime.getMean()+gtcTime.getMean()));
		}
		else if (idioma==L_english)
		{
			fprintf(fp, "Descript generation : %.4f[s] +- %.4f[s] #%d\n", gpiTime.getMean(), sqrt(gpiTime.getVariance()), gpiTime.size());
			fprintf(fp, "Database generation : %.4f[s] +- %.4f[s] #%d\n", gbdTime.getMean(), sqrt(gbdTime.getVariance()), gbdTime.size());
			fprintf(fp, "Matches generation  : %.4f[s] +- %.4f[s] #%d\n", calTime.getMean(), sqrt(calTime.getVariance()), calTime.size());
			fprintf(fp, "Transf generation   : %.4f[s] +- %.4f[s] #%d\n", gtcTime.getMean(), sqrt(gtcTime.getVariance()), gtcTime.size());
			fprintf(fp, "Total frame rate    : %.1f[fps]\n", 1/(gpiTime.getMean()+gbdTime.getMean()+calTime.getMean()+gtcTime.getMean()));
		}
	}
};

class L_RecPuntInt_subObjCopiable // Componentes de L_RecPuntInt que pueden copiarse rápidamente y directamente usando operador=
{
public:
	bool aprendiendo;
	bool cuenta_calces_1_a_n;
	bool cuentaDistrSigmas;
	bool cuentaDistrS_imprDiv;
	L_RecTimeSeg ti;
	struct restrSemiLocales {bool apply; bool corrige; int nVecinos; int minCorresp; bool verifEscalaAngs; double varEscala; double varAng;} RSL;
	struct deltaAngCalces_ {bool apply; double deltaMax;} deltaAngCalces;
	struct cuentaCerros_ {bool apply; bool corrige; double minDistCerros; double porcCerrosCorrectosMin; 
		double maxDistCerros; int minCerrosUtil; double parecidoCerrosMin; int nLineasMax;} cuentaCerros;
	
	struct caractRec_ {int nCal; int nTr;} caractRec;

	const L_ImageGrayDouble *impru; // Punteros a imagen de test, es usado por dibujaTransformaciones()
	int nobj; // numero correlativo para diferenciar los objetos almacenados en la base de datos
	bool conoceImagenesEntrenamiento;

	L_GenDescrPuntInt *gpi; // Generador de puntos de interes y descriptores (ej: sDoG+SIFT)
	L_ConjDescrBus *gbd;    // Generador de base de datos para busqueda de descriptores (ej: kdtree+BBF)
	L_GenTransfCalces *gtc; // Generador de transformaciones a partir de listas de calces (ej: transf afin)
	L_CalcDistrSigmas distrS; // Para visualizar distribucion de los sigmas
	L_FiltroDescriptores *filtD; // Para no considerar descriptores que no cumplan con un criterio

	long cal_pool_max;

#ifdef __TURBOC__
	__OP_ASSIGN(L_RecPuntInt_subObjCopiable);
#endif
	void swap(L_RecPuntInt_subObjCopiable &other) {L_RecPuntInt_subObjCopiable t = other; other = *this; *this = t;}
};

class L_RecPuntInt : public L_RecPuntInt_subObjCopiable
{
private:
	L_RecPuntInt(const L_RecPuntInt &other); //:paramRPI("RPI",5) { L_hard_shutdown("Constructor no valido");} // No debe ser llamado
	L_RecPuntInt& operator =(const L_RecPuntInt &other); // { L_hard_shutdown("Constructor no valido"); return *this; } // No debe ser llamado
public:
	L_ParamManagerLocal paramRPI;
	L_ParamLabelList paramNoUsados;
	L_VerifPixelsTransfBN corrPixBN; // Puede trabajar con cualquier tipo de transformacion

	L_Array<const L_ImageGrayDouble *> imsref; // Punteros a imagenes de entrenamiento, es usada por dibujaTransformaciones()

	L_CalceLista calL;       // Todos los calces generados
	L_DescriptorLista desLRef;  // Lista de descriptores para uso temporal
	L_DescriptorLista desLPru;  // Lista de descriptores para uso temporal
	L_TransfPuntInt2DListaPtr traL; // Resultado (output)

	L_List<L_Calce> cal_pool;  // Para guardar nodos y no tener que hacer new todo el rato

	L_RecPuntInt():paramRPI("RPI",5)
	{
		conoceImagenesEntrenamiento = true;
		impru=NULL;
		aprendiendo=false;  gpi=NULL; gbd=NULL; gtc=NULL; nobj=0; //maxDist=1; maxRatio=0.8;
		cuentaDistrSigmas=false; filtD=NULL; cuenta_calces_1_a_n=false; //numCalcesMasCercanos=1;
		deltaAngCalces.apply=false; deltaAngCalces.deltaMax=40*M_PI/180;
		RSL.apply=false; RSL.corrige=false; RSL.minCorresp=3; RSL.nVecinos=4; RSL.verifEscalaAngs=true; RSL.varAng=30*M_PI/180; RSL.varEscala=1.7;
		cuentaCerros.apply=false; cuentaCerros.minDistCerros=5.0; cuentaCerros.maxDistCerros=20; cuentaCerros.porcCerrosCorrectosMin=0.6;
		cuentaCerros.minCerrosUtil=4; cuentaCerros.parecidoCerrosMin=0.8; cuentaCerros.nLineasMax = 100;
		cal_pool_max = 100000;
		paramRPI.addFrom("cuentaDistrSigmas",&cuentaDistrSigmas);
		paramRPI.addFrom("cuentaDistrS_imprDiv",&cuentaDistrS_imprDiv);
		paramRPI.addFrom("cuenta_calces_1_a_n",&cuenta_calces_1_a_n);
		paramRPI.addFrom("deltaAngCalces.apply",&deltaAngCalces.apply);
		paramRPI.addFrom("deltaAngCalces.deltaMax",&deltaAngCalces.deltaMax);
		paramRPI.addFrom("RSL.apply",&RSL.apply);
		paramRPI.addFrom("RSL.corrige",&RSL.corrige);
		paramRPI.addFrom("RSL.nVecinos",&RSL.nVecinos);
		paramRPI.addFrom("RSL.minCorresp",&RSL.minCorresp);
		paramRPI.addFrom("RSL.verifEscalaAngs",&RSL.verifEscalaAngs);
		paramRPI.addFrom("RSL.varEscala",&RSL.varEscala);
		paramRPI.addFrom("RSL.varAng",&RSL.varAng);
		paramRPI.addFrom("cuentaCerros.apply",&cuentaCerros.apply);
		paramRPI.addFrom("cuentaCerros.corrige",&cuentaCerros.corrige);
		paramRPI.addFrom("cuentaCerros.minDistCerros",&cuentaCerros.minDistCerros);
		paramRPI.addFrom("cuentaCerros.maxDistCerros",&cuentaCerros.maxDistCerros);
		paramRPI.addFrom("cuentaCerros.minCerrosUtil",&cuentaCerros.minCerrosUtil);
		paramRPI.addFrom("cuentaCerros.parecidoCerrosMin",&cuentaCerros.parecidoCerrosMin);
		paramRPI.addFrom("cuentaCerros.porcCerrosCorrectosMin",&cuentaCerros.porcCerrosCorrectosMin);
		paramRPI.addFrom("cuentaCerros.nLineasMax",&cuentaCerros.nLineasMax);
		paramRPI.addChildren(&corrPixBN.paramCPBN);
	}
	void swap(L_RecPuntInt &other);
	inline void fijaReconocedor(L_GenDescrPuntInt *gpi, L_ConjDescrBus *gbd, L_GenTransfCalces *gtc, L_FiltroDescriptores *filtD=NULL)
	{
		this->gpi=gpi;
		this->gbd=gbd;
		this->gtc=gtc;
		this->filtD=filtD;
		paramRPI.olvidaHijos();
		if (gpi->pideParams()!=NULL)
			paramRPI.addChildren(gpi->pideParams());
		if (gbd->pideParams()!=NULL)
			paramRPI.addChildren(gbd->pideParams());
		if (gtc->pideParams()!=NULL)
			paramRPI.addChildren(gtc->pideParams());
		if (filtD != NULL && filtD->pideParams()!=NULL)
			paramRPI.addChildren(filtD->pideParams());
	}
	bool reconocedorFijado() {return gpi!=NULL && gbd!=NULL && gtc!=NULL;}
private:
	bool agregaImagenEntren_interno(L_ImageGrayDouble &imRef, int modoAlmacenamiento, double s0, L_DescriptorLista *precalc);
public:
	bool agregaImagenEntrenRobaObjetoDe(L_ImageGrayDouble &imRef, double s0=-1, L_DescriptorLista *precalc=NULL) {return agregaImagenEntren_interno(imRef, 0, s0, precalc);}
	bool agregaImagenEntrenCopiaObjetoDe(L_ImageGrayDouble &imRef, double s0=-1, L_DescriptorLista *precalc=NULL) {return agregaImagenEntren_interno(imRef, 1, s0, precalc);}
	bool agregaImagenesEntren_usa_desL(L_DescriptorLista &desL, L_ImageGrayDouble *imRef = NULL, bool copiarImagen = false); // Roba el contenido de desL. En este caso no se puede llamar dibujaTransformaciones()

	bool finalizaEntren(); // transformUsing this->desLPru en el tree con los descriptores del entrenamiento

	bool procesaImagenesPrueba(L_ImageGrayDouble &imPru, double s0=-1, L_ImageRGBUchar *dibCalcesDebug=NULL);
	bool procesaImagenesPrueba_arch(L_ImageGrayDouble &imPru, double s0=-1, const char *dibCalcesDebug=NULL);
	bool procesaImagenesPrueba_usa_desL_arch(const L_ImageGrayDouble *imPru, const char *dibCalcesDebug=NULL); // Esta funcion no compute los descriptores de test, usa los de this->desL
	bool procesaImagenesPrueba_usa_desL(L_DescriptorLista &desL, const L_ImageGrayDouble *imPru=NULL, const char *dibCalcesDebug=NULL); // Si se llama con imPru=NULL no se puede usar dibujaTransformaciones()
	// La que hace el trabajo de verdad
	bool procesaImagenesPrueba_interno(const L_ImageGrayDouble *imPru, L_ImageRGBUchar *dibCalcesDebug=NULL); // Esta funcion no compute los descriptores de test, usa los de this->desL

	bool procesaReferenciaPruebaExt(L_ImageGrayDouble &imRef, L_ImageGrayDouble &imPru, double s0=-1, L_ImageRGBUchar *dibCalcesDebug=NULL);

	void reseteaPrueba() {traL.destroyList(&cal_pool); desLPru.clear(); calL.moveListTo(cal_pool); if ((long)cal_pool.size() > cal_pool_max) cal_pool.clear();}
	void reseteaEntren(); //! Esta funcion sirve para dejar el sistema listo para ser reentrenado

	// Funciones para dibujar transformaciones. Al usuario le corresponde asegurarse de que las imagenes necesarias no han sido eliminadas
	void generaArregloImagenesParaCalces(L_Array<const L_ImageGrayDouble *> &imarr) const; // se llama im.dibCalces(imarr.elem, imarr.size(), lins);
	bool dibujaCalcesIniciales(L_ImageRGBUchar &imSal);
	bool dibujaCalcesIniciales(const char *nomarch);
	bool dibujaTransformaciones(L_ImageRGBUchar &imSal, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica=-1, int transfEspecifica=-1);
	bool dibujaTransformaciones(const char *nomarch, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica=-1, int transfEspecifica=-1);
	void imprimeReporte(FILE *fp=NULL);

	void transformaImagenPruebaEnEntrenamiento(); // Funcion rara para que la image actual de test sea la futura de entrenamiento

	bool debugActivo(int desact0_act1=-1) {if (desact0_act1==0) {gtc->debFiltr.activo=corrPixBN.debFiltr.activo=false;} else if (desact0_act1==1) {gtc->debFiltr.activo=corrPixBN.debFiltr.activo=true;}; return gtc->debFiltr.activo && corrPixBN.debFiltr.activo;}

	bool leeParametrosPropaga(const char *arch);
	bool grabaParametrosPropaga(const char *arch);
};

class L_RecPuntIntConjuntoTracks
{
public:
	L_Array<L_Track> arr;
	L_RecPuntInt rec;
	int nMax;
	int frameMax;
	int ultId;
	bool actualizarEntrenamiento;

	void fijaReconocedor(L_GenDescrPuntInt *gpi, L_ConjDescrBus *gbd, L_GenTransfCalces *gtc) {rec.fijaReconocedor(gpi, gbd, gtc);}
	int procesaImagen(const L_ImageGrayDouble &im, const L_ImageGrayDouble &imc2, long frame); // imc2.pix puede ser NULL
	void pide_xy_n_frames(int nFrames, int &nTracks, L_Matrix &xy, long frame);
	void pide_xy_n_frames_c2(int nFrames, int &nTracks, L_Matrix &xy, long frame);
	void pide_uv_n_frames(int nFrames, int &nTracks, const L_Camara *c, L_Matrix &uv, long frame);
};

// class_name nueva, usa nf frames hacia atras, es la que se hizo para Brasil
class L_RecPuntIntTrackPuntos
{
public:
	L_Array<L_ImageGrayDouble> imAnt;
	L_Array<L_DescriptorLista> desAnt;
	int nAnt;
	int nf;
	L_Matrix tracks;
	int frame;
	int ultTrack; // = num descriptores

#ifdef __L_LOWE_H_
	L_LoweGenDescr sift;
#else
	L_PDoG_SIFT sift;
#endif
	L_KdTreeBBF tree;
	//L_Flann tree;
	L_GenTransfAfinHough hough;
	L_RecPuntInt rec;


	// Formato de tracks:  x der,  y aba
	//
	//  numFeatures numFrames
	//  x  y  1    x  y  1    x  y  1     -> num frames
	// -1 -1 -1    x  y  1    x  y  1
	//
	//  |
	//  V  num features

	L_RecPuntIntTrackPuntos();

	static void main_procesar();
	//
	bool procesar(L_CapturadorImagen &capt, int nMax = -1);
	bool procesar_img(L_ImageGrayDouble &imNueva);
	void preparar_matriz();
	void resetear();

	static void genDibujosTracks();
};




// Nueva
class L_CalibradorCamaraRPI
{
private:
	L_ImageGrayDouble imInter; // La guardo para poder referenciarla al graficar las transformaciones
public:
	// Reconocedor puntos de interes
#ifdef __L_LOWE_H_
	L_LoweGenDescr genPuntos;
#else
	L_PDoG_SIFT genPuntos;
#endif
	L_KdTreeBBF tree;
	L_GenTransfProyectivaRansac_1tr genProy;
	L_RecPuntInt rpi;
	bool rpiEntrenado;
	// Datos de interes externo
	L_ImageGrayDouble imRef;
	char nomImagen[1000];
	double lxMM, lyMM;
	// Datos para calibrar camara
	L_CuadroCalibradorCamaraArr puntos; // Puntos para calibrar
	L_CamaraDistorsionRadial camInter;  // Modelo de camara intermedia
	L_CalibradorCamaraRPI();
	bool fijaImagenReferencia(const char *nomarch) {if (imRef.readImage(nomarch)==0) return false; strcpy(nomImagen, nomarch); rpiEntrenado=false; return true;}
	bool procesaImagenReferencia();
	bool procesaNuevaImagen(L_ImageGrayDouble &im);
	bool leerArchivoConfiguracion(const char *nomarch);
	bool grabarArchivoConfiguracion(const char *nomarch);
	bool leerArchivoPuntos(const char *nomarch) {FILE *fp = fopen(nomarch, "r"); if (fp==NULL) return false; puntos.read(fp); fclose(fp); return true;}
	bool grabarArchivoPuntos(const char *nomarch) {FILE *fp = fopen(nomarch, "w"); if (fp==NULL) return false; puntos.guardar(fp); fclose(fp); return true;}
	bool leerParametrosCamara(const char *nomarch) {FILE *fp=fopen(nomarch, "r"); if (fp==NULL) return false; bool a; a=camInter.leerParametros(fp); fclose(fp); return a;}
	bool grabarParametrosCamara(const char *nomarch) {FILE *fp=fopen(nomarch, "w"); if (fp==NULL) return false; camInter.guardarParametros(fp); fclose(fp); return true;}
	void resetearCamara() {camInter.resetearParametros();}
	bool calibrar_Tsai();
};

#ifdef __COMPAT_FLTK__
class L_CalibradorCamaraFLTK
{
private:
	// No declaro estas variables como static para no usar memoria innecesariamente (por si no se usa esta class_name)
	L_CalibradorCamaraRPI calib;
	bool salir;
	int estado; // 0 = descansando / calibracion manual, 1 = calibracion automatica
	bool modificark2;
	bool pedirCapturaCalces;
	bool pedirCambioCamara;
	L_VentanaImagenFLTK *venCam;
	L_VentanaImagenFLTK *venRef;
	L_VentanaImagenFLTK *venCal;
	L_ImageRGBUchar imCorr;

	static void cb_leer_calibracion_camara(Fl_Widget *widget, void *param);
	static void cb_guardar_calibracion_camara(Fl_Widget *widget, void *param);
	static void cb_leer_puntos(Fl_Widget *widget, void *param);
	static void cb_guardar_puntos(Fl_Widget *widget, void *param);
	static void cb_salir(Fl_Widget *widget, void *param);
	static void cb_resetear_calibracion_camara(Fl_Widget *widget, void *param);
	static void cb_fijar_fovx_manual(Fl_Widget *widget, void *param);
	static void cb_fijar_fovy_manual(Fl_Widget *widget, void *param);
	static void cb_fijar_dfx_manual(Fl_Widget *widget, void *param);
	static void cb_fijar_dfy_manual(Fl_Widget *widget, void *param);
	static void cb_fijar_k1_manual(Fl_Widget *widget, void *param);
	static void cb_fijar_k2_manual(Fl_Widget *widget, void *param);
	static void cb_permitir_modif_k2(Fl_Widget *widget, void *param);
	static void cb_fijar_referencia(Fl_Widget *widget, void *param);
	static void cb_fijar_ancho_mm(Fl_Widget *widget, void *param);
	static void cb_fijar_alto_mm(Fl_Widget *widget, void *param);
	static void cb_cambiar_camara(Fl_Widget *widget, void *param);
	static void cb_leer_config_referencia(Fl_Widget *widget, void *param);
	static void cb_guardar_config_referencia(Fl_Widget *widget, void *param);
	static void cb_capturar_calces(Fl_Widget *widget, void *param);
	static void cb_borrar_calces(Fl_Widget *widget, void *param);
	static void cb_calibracion_camara_automatica(Fl_Widget *widget, void *param);
public:
	int main_calibrCam_obj(int argc, char *argv[]);
	static int main_calibrCam(int argc, char *argv[]) {L_CalibradorCamaraFLTK obj; return obj.main_calibrCam_obj(argc, argv);}
};
#endif // __COMPAT_FLTK__



class L_TemplateDescr:public L_GenDescrPuntInt
{
public:
	// Generadores de puntos
	L_MinuciasExterno gMi; // No build los puntos, es solo para los parametros
#ifdef __L_LOWE_H_
	L_LoweGenDescr gLo;
#else
	L_PDoG_SIFT gLo;
#endif
	L_PHarris_SIFT gHa;
	L_Harris_SIFT gHa1;
	L_GenDescrArr gUnion;

	L_ImageGrayDouble im;
	double dt;
	char nomIm[512];
	bool compr;

	L_ParamManagerLocal paramTempl;

	L_TemplateDescr();
	void fijaGeneradores(bool usaMinuc, bool usaLSift, bool usaPHarris, bool usaHarris);
	bool calcDescrPuntInt(L_ImageGrayDouble &im, double s0=0.5); //!< Pasa la tarea de generar los descriptores a sus sub-generadores
	bool procesaEntrada(const char *name); // Procesa los archivos xxx.bmp, xxx_min.bmp, xxx_min.txt
	bool leeTemplate(const char *name);
	bool grabaTemplate(const char *name);
	virtual L_ParamManagerLocal *pideParams() {return gUnion.pideParams();}
};




class L_ProgramasRecPuntInt
{
public:
	static int main_MuestraPuntos(int argc, char *argv[]);
	static int main_GrafPuntInt(int argc, char *argv[]);
	static int main_RecPuntInt(int argc, char *argv[]); // build archivo texto con transformaciones
	static int main_ProcesaImagenesGenTxtCalces(int argc, char *argv[]);
	static int analizaResultadosUch(const char *metodo, int n1, int n2, int n3, int n4, int width=3);
	static int main_VerificaResultados(int argc, char *argv[]);
	static int main_CuentaDetecciones(int argc, char *argv[]);
	static int main_FVC_enroll(int argc, char *argv[]); // build template
	static int main_FVC_match(int argc, char *argv[]); // Compara imagen y template
	static int main_FVC_matchT(int argc, char *argv[]); // Compara template y lista de template
	static int main_FVC_matchT2(int argc, char *argv[]); // Compara 2 listas de templates
	static int main_pruebaCruzada(int argc, char *argv[]);
	static int main_pruebaCruzada_1v1(int argc, char *argv[]);
	static int main_CalcTransf(int argc, char *argv[]); // Proyecto detector objetos doblados

#ifdef __COMPAT_FLTK__
	static int main_calibrCam(int argc, char *argv[]) {return L_CalibradorCamaraFLTK::main_calibrCam(argc, argv);}
#endif

#if defined(__COMPAT_FLTK__) || defined(__COMPAT_IPLIMAGE__)
	static int main_captStereo(int argc, char *argv[]);
	static int main_camRecPuntInt(int argc, char *argv[]);
	static int main_camRecPuntInt_video_0(int argc, char *argv[]);
	static int main_camRecPuntInt_video_1(int argc, char *argv[]);
	static int main_camRecPuntInt_pose_plana(int argc, char *argv[]);

	static int main_trackCaja(int argc, char *argv[]);
	static int main_dibModelo(const L_Descriptor3DArreglo &des3DArr, int xIni = 200);

		//
	static int main_camCalcPosePuntInt_1(int argc, char *argv[]);
	static int main_camCalcPosePuntInt_2(int argc, char *argv[]);
#endif // __COMPAT_FLTK__ || __COMPAT_IPLIMAGE__
};

#endif //__L_RECPUNTINT_H__

