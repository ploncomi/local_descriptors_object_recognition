#ifndef __L_DESCRIPTORES_H__
#define __L_DESCRIPTORES_H__

#include "L_Fnes_Pato.h"

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


//#include "../flann/flann.h"
//#include "../RobHessSift/imgfeatures.h"

////////////////////////////////////
/////
///// Clases de bajo nivel para define puntos de interes, calces, etc
/////
//////////////////////////////
// Clases generalizadas

//////////
// Tipos de puntos de interes y de descriptores

enum L_PuntIntTipo {PUNNADA, P_UNION, DOG, HARRIS, SDOG, SHARRIS, PDOG, PHARRIS, LOWE_DOG, P_SURF, BAY_P_SURF, P_OPENSURF, P_EXT, MINUCIA, P_CVSIFT};
enum L_DescriptorTipo {DESNADA, D_UNION, SIFT, LOWE_SIFT, D_SURF, BAY_D_SURF, D_OPENSURF, D_SURFGAUSS, D_EXT, D_CVSIFT};
enum L_LandTipo {LANDNADA, LANDPUN, LANDPUNINV, LANDCONJ, LANDLASER};

class L_GenDescrPuntInt;

#define L_MAX_NT 128 // SIFT = 128, SURF = 64... es dificil que sea mayor a 128, si no se cambia

class L_Descriptor3D;
struct feature;

class L_Descriptor //! Nodo que contiene un descriptor genérico
{
public:
	// Para debug, para "marcar" los descriptores, prefiero comentarlo a poner defines
	//char zona[16];

	//Variables para uso general
	double x0; // posic x respecto a la imagen original
	double y0; // posic y respecto a la imagen original
	double sigma0; // nivel de blur en el scalespace (imagen original tiene blur 0.5 generalmente)
	double ang;

	double vector[L_MAX_NT];
	int nt; //! Largo del arreglo vector

	double intens; // value de intensidad del maxElement o minElement; no se hacen calces cuando los intens tienen signos opuestos
	int signo;
	int nobj;
	int imorig_lx;
	int imorig_ly;

	const L_ImageGrayDouble *imorig; // Referencia a la imagen original, no tiene memoria propia
	L_GenDescrPuntInt *refGen; //! Referencia al generador de puntos de interes usado para crear el descriptor
	L_PuntIntTipo tipoP;
	L_DescriptorTipo tipoD; // Tipo de descriptor

	// Para SLAM
	L_LandTipo tLand;
	int indLand;    // Posicion en slam.arrL[]
	int indSubLand; // Posicion en pun[] de un landmark solido
	int cam; // 0 = principal, 1 = secundaria
	L_RayoCamara ray; // Tangente respecto a la camara
	double cov[3][3];
	double xAde;  // Para ver si esta delante de la camara
	char hasDepth; // Only values 0, 1 allowed. Using bool generate data misalignment
	long id;  // Numeros correlativos

	// Para structure from motion
	L_Descriptor3D *d3d;


	//Variables temporales para generacion de puntos de interes (posicion en las piramides)
	double pir_x; // posic x en la piramide (uso temporal)
	double pir_y; // posic y en la piramide (uso temporal)
	int pir_n; // numero de piso en la piramide (uso temporal)
	L_Frac pir_reltam; // relacion de tamaño en la piramide, x=x0*reltam (aproximadamente)

public:
	void creaVectorPromedioDe(const L_Descriptor &d1, const L_Descriptor &d2);

	double x0PuntaFlecha() const {return x0 + sigma0 * cos(ang);} // Ese maldito const siempre se me va...
	double y0PuntaFlecha() const {return y0 - sigma0 * sin(ang);}

	void swap(L_Descriptor &other) {L_Descriptor t = other; other = *this; *this = t;}

	double distEuclidCuadr(const L_Descriptor & other) const;
	double distEuclidCuadr_mismoTipo(const L_Descriptor & other) const; // Penaliza mucho los calces entre tipos distintos de descriptores
	bool normalizaHist();
	bool normalizaHistCota(double maxValorComponente);
	int cmpPosSigmaAng(const L_Descriptor &other) const;
	int cmpPtr(const L_Descriptor &other) const;
	int cmpMula(const L_Descriptor &other) const;
	static int qCmpMula_nodo(const void *pp1, const void *pp2);
	static int qCmpMulaPtr_nodo(const void *pp1, const void *pp2);
	static int qCmpIntens(const void *p1, const void *p2);
	static int qCmpIntens_nodo(const void *pp1, const void *pp2);
	static int qCmpIntensInv_nodo(const void *pp1, const void *pp2);
	double medidaHarris(L_Array<L_ImageGrayDouble> &tmp, double sd=1.3, double si=2.0, double alfa=0.04);

	void computeRayFromPixel(const L_CamaraPinhole &cam) {cam.calcRayo(ray, x0, y0);}
	void computeCovariance3DFromDisparity(double sigmaPixNormX, double sigmaPixNormY, double sigmaPixNormDisparity, double baselinemm);

	void writeElementTxt(FILE *fp);
	void readElementTxt(FILE *fp, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen);
	void writeElementBin(FILE *fp);
	void readElementBin(FILE *fp, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen);
	void readElementBin(FILE *fp, va_list ppp);
	int writeElementMem(std::vector<char> &buf);
	int readElementMem(const std::vector<char> &buf, int i0, const L_ImageGrayDouble *im, L_GenDescrPuntInt *gen);
	int readElementMem(const std::vector<char> &buf, int i0, va_list va);

	inline static bool estaAtras(const L_Descriptor &d) {return d.xAde < 0;}
	static L_Descriptor creaDescriptorAzar(int lx, int ly, int dim, L_GenDescrPuntInt *gen);
	static L_Descriptor creaDescriptorAzarSemilla(int lx, int ly, int dim, L_GenDescrPuntInt *gen, L_Rand &randF);

#ifdef IMGFEATURES_H
	void OP_assign(const feature &feat);
	void copyTo(feature &feat);
#endif // IMGFEATURES_H

	inline int idTipo() const {return tipoP*20 + tipoD;}
#ifdef __TURBOC__
	__OP_ASSIGN(L_Descriptor);
#endif
};


// Tipos antiguos, ahora se homologan a L_Descriptor

typedef L_Node<L_Descriptor> L_DescriptorNodo;

class L_DescriptorLista : public L_List<L_Descriptor>
{
public:
	void dibujaFlechas(L_ShapeArray &lins, int r=-1, int g=-1, int b=-1, int grosor=1);
	void dibujaCuadrados(L_ShapeArray &lins);
	bool leeArchivoPato(const char *nomarch, const L_ImageGrayDouble &im, L_GenDescrPuntInt *gen, int nobj);
	bool grabaArchivoPato(const char *nomarch);
	bool leeArchivoPato_minu(const char *nomarch, L_ImageGrayDouble &imDummy, L_DescriptorLista &punL, bool agregaEspejo); // HAY QUE ACTUALIZAR ESTA FUNCION
	bool grabaArchivoPato_minu(const char *nomarch);
	void transformaEnMatriz(L_Matrix &uv);
	bool revisaCoherencia();

	void fijaPunteroImagen(const L_ImageGrayDouble *im) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.imorig=im;}
	void fijaPunteroGenerador(L_GenDescrPuntInt *gen) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.refGen=gen;}
	void fijaNumObj(int nobj) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.nobj=nobj;}

	// Para debug, para "marcar" los descriptores, prefiero comentarlo a poner defines
	//void fijaZona(const char *zona) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) strcpy(p->c.zona,zona);}

	void fijaIdCorrelativos(int id0 = 0) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.id=id0++;}
	void fijaIds(int id = 0) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.id=id;}
	void fijaCam(int numcam = 0) {L_DescriptorNodo *p; for (p=root; p!=NULL; p=p->sig) p->c.cam=numcam;}

	void copia_a_arreglo(std::vector<L_Descriptor> &arr);
	void copia_a_arreglo_punteros(std::vector<L_Descriptor *> &arr);
	static int qCmpIntens(const void *des1, const void *des2) { if ( (*(L_DescriptorNodo**)des1)->c.intens < (*(L_DescriptorNodo**)des2)->c.intens ) return 1; return -1;}

	void pedirMatrizDatos(L_Matrix &xy);
	void pedirMatrizDatosRayos(L_Matrix &xy, const L_CamaraPinhole &cam);

	void agregarRuidoUniforme(double deltaPix);
	void agregarRuidoUniformeSemilla(double deltaPix, L_Rand &randF);
	void eliminaLineasHarris(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	void eliminaLineasHarris_indiv(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	void eliminaLineasHarris_ant(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	int cuentaVisibles(int lx, int ly, double borde = 0);
	void eliminaNoVisibles(int lx, int ly, double borde = 0);
	// Estos eliminan calces correctos usando el algoritmo de los 3 puntos
	int eliminaPose_Ransac3Puntos(const L_CamaraPinhole &cam, double errMax, int nIter, int nMin, int nAceptAutom, const std::vector<L_CoordsCart3D> &pun, const std::vector<double> &puntajes, L_Pose3D_cuat &pose);

	template <class T> void fillDepthGen(const T &imDepth) {L_DescriptorNodo *d; int u, v; for (d=root; d!=NULL; d=d->sig) {d->c.hasDepth = 0; u = (int)d->c.x0; v = (int)d->c.y0; if (u>0 && u<imDepth.cols()-1 && v>0 && v<imDepth.rows()-1 && imDepth.pix(u,v) > 0) {d->c.xAde = imDepth.pix(u,v); d->c.hasDepth = 1;}}}
	void fillDepth(L_ImageGrayDouble &imDepth);
	void fillRayTangents(const L_CamaraPinhole &cam) {L_DescriptorNodo *d; for (d=root; d!=NULL; d=d->sig) cam.calcRayo(d->c.ray, d->c.x0, d->c.y0);}
	void fillCovarianceInverseDistance(double sigmaPixNormX, double sigmaPixNormY, double sigmaPixNormDisparity, double baselinemm); // Requires previous call to fillRayTangents

	void fillCovarianceInverseDistance_v0(double sigmar0, double sigmaw); // Requires previous call to fillRayTangents

#ifdef IMGFEATURES_H
	void OP_assign(const feature *feat, int n);
	void copyTo(feature *feat);
#endif // IMGFEATURES_H

};

typedef L_Descriptor L_PuntInt;
typedef L_DescriptorNodo L_PuntIntNodo;
typedef L_DescriptorLista L_PuntIntLista;

typedef L_Descriptor L_Direcc;
typedef L_DescriptorNodo L_DireccNodo;
typedef L_DescriptorLista L_DireccLista;


class L_DescriptorPtr
{
public:
	L_DescriptorNodo *ptr;
	static int qCmpMula(const void *p1, const void *p2);
};

typedef L_Node<L_DescriptorPtr> L_DescriptorPtrNodo;
typedef L_List<L_DescriptorPtr> L_DescriptorPtrLista; // Al destruir esta lista, no se destruyen los descriptores



class L_Descriptor3D
{
public:
	L_Descriptor d;
	double sigma1;  // sigma a una distancia de "1"
	L_CamaraPinhole cam;
	// {nor,gra} = sis de referencia local centrado en cen
	L_CoordsCart3D cen;  // vector posicion
	L_CoordsCart3D nor;  // vector normal
	L_CoordsCart3D gra;  // vector gradiente
	
	// Crear en una posicion pos, se debe indicar la direccion normal, una direccion horizontal y la distancia de la camara al descriptor (para normalizar la escala)
	void crear(const L_Descriptor &d, const L_CoordsCart3D &pos, const L_CoordsCart3D &normal, const L_CoordsCart3D &horz, double distAde, const L_CamaraPinhole &cam);

	void movido_a_fijo(const L_Pose3D_cuat &p, const L_Descriptor3D &orig); // Incluye el descriptor proyectado
	void movido_a_fijo_normalImprecisa(const L_Pose3D_cuat &p, const L_Descriptor3D &orig, double proporcionCorreccion = 0.5); // Incluye el descriptor proyectado
	void fijo_a_movido(const L_Pose3D_cuat &p, const L_Descriptor3D &orig); // Incluye el descriptor proyectado

	void swap(L_Descriptor3D &other) {L_Descriptor3D t = other; other = *this; *this = t;}
};

class kbParams
{
public:
	L_SistemaEjes3D sis;
	bool invertirOri;
	bool softNormal;
	kbParams() : sis(L_DerAbaAde), invertirOri(false), softNormal(true) {}
};

class L_Descriptor3DArreglo : public L_Array<L_Descriptor3D>
{
public:
	void leerFormatoB(const char *name);
	void calcPose(L_Pose3D_cuat &pose, std::vector<L_Descriptor> &dIm);

	static void errProy(const void *data, const L_Matrix &pq, L_Matrix &err);

	void generarDescriptorLista(L_DescriptorLista &desL, L_Pose3D_cuat &pCam);
	void generarDescriptorLista_normalImprecisa(L_DescriptorLista &desL, L_Pose3D_cuat &pCam);
	void generarDescriptorListaOriginal(L_DescriptorLista &desL);
	void genDrawing(L_ShapeArray &lins, L_Pose3D_cuat &pCam);
	void leerArchivoKb(L_CamaraPinhole &cam, L_Array<L_HomogeneousMatrix> &arrH, std::vector<int> &poses, const char *name, const kbParams &params);
};

struct L_Descr3D_datos
{
	L_Descriptor3DArreglo *arreglo;
	L_Matrix *uv;
};




class L_xyst
{
public:
	double x;
	double y;
	double sigma; // Si existe se usa
	double theta; // Si existe se usa
	long frame;
	double ti;
	L_xyst() {frame=-1;}
	void define(const L_Descriptor &d, long frameNum, double tiem = -1) {x=d.x0; y=d.y0; sigma=d.sigma0; theta=d.ang; frame=frameNum; ti = (tiem>0) ? tiem : frameNum;}
};

class L_xyst_2cam
{
public:
	double x;
	double y;
	double sigma; // Si existe se usa
	double theta; // Si existe se usa
	double xc2;
	double yc2;
	double sigmac2; // Si existe se usa
	double thetac2; // Si existe se usa
	long frame;
	long framec2;
	double ti;
	L_xyst_2cam() {frame=-1;}
	void define(const L_Descriptor &d, long frameNum, double tiem = -1) {x=d.x0; y=d.y0; sigma=d.sigma0; theta=d.ang; frame=frameNum; ti = (tiem>0) ? tiem : frameNum;}
	void define_c2(const L_Descriptor &d, long frameNum, double tiem = -1) {xc2=d.x0; yc2=d.y0; sigmac2=d.sigma0; thetac2=d.ang; framec2=frameNum; ti = (tiem>0) ? tiem : frameNum;}
};

class L_Track_POD
{
public:
	L_uchar cR, cG, cB;
	double x, y;
	double xc2, yc2;
	long id;
#ifdef __TURBOC__
        __OP_ASSIGN(L_Track_POD)
#endif // __TURBOC__
};

class L_Track : public L_Track_POD
{
public:
	L_Array<L_xyst_2cam> t;
	L_Descriptor d;
	void createFrom(const L_Descriptor &des, long frame, long idNum, L_uchar R=0, L_uchar G=0, L_uchar B=0, double tiem = -1);
	void addFrom(const L_Descriptor &des, long frame, double tiem = -1);
	void agrega_c2(const L_Descriptor &des, long frame, double tiem = -1);
	void addFrom(const L_Descriptor &des, const L_Descriptor &desc2, long frame, double tiem = -1);
	L_Track() {}
	L_Track(const L_Track &tr) {t=tr.t; d=tr.d; L_Track_POD::operator=(tr);}
	L_Track& operator=(const L_Track &tr) {t=tr.t; d=tr.d; L_Track_POD::operator=(tr); return *this;} // Hay que hacerlo pq los subobjetos tienen punteros
};

class L_ConjuntoTracks
{
public:
	L_Array<L_Track> arr;
	double dx, dy, dsigmaRatio, dtheta; //dsigmaRatio entre 0 y 1
	long nMax, frameMax; // numero maxElement de tracks y tiempo maxElement de olvido de tracks
	double maxDist, maxRatio;
	long ultId;

	L_ConjuntoTracks() { dx=-1.0e10; dy=-1.0e10; dsigmaRatio=-1.0e10; dtheta=-1.0e10; nMax = -1; frameMax = -1; maxDist = 1.0; maxRatio = 0.8; ultId = 0;} // Parte no inicializado
	bool inicializado() {return dx > -1.0e9 && dy > -1.0e9 && dsigmaRatio > -1.0e9 && dtheta > -1.0e9 && nMax >= 1 && frameMax >= 0;}
	void procesa1(const L_DescriptorLista &desL, long frame); // Mantiene un conjunto de nMax tracks de descriptores
	void procesa2(const L_DescriptorLista &desL, long frame, L_ImageRGBUchar *im=NULL); // Mantiene un conjunto de nMax tracks de descriptores, other forma
	void pideDescriptores(std::vector<L_Descriptor> &ret, long frame);
	void procesaFalso(std::vector<L_CoordsCart3D> &p3d, const L_CamaraPinhole &c, L_Pose3D &pose, long frame, int *nVisibles); // Mantiene un conjunto de nMax tracks de descriptores, other forma
	void pide_xy_n_frames(int nFrames, int &nTracks, L_Matrix &xy, long frame);
	void pide_uv_n_frames(int nFrames, int &nTracks, const L_Camara *c, L_Matrix &uv, long frame);
};

////////
// Clases para define filtros para eliminar descriptores en listas de descriptores

class L_FiltroDescriptores
{
public:
	virtual bool criterioEliminacion(const L_DescriptorNodo *des)=0; // Devuelve true => eliminar des
	virtual void filtrarNodos(L_DescriptorLista &desL);
	virtual L_ParamManagerLocal *pideParams() {return NULL;}
	virtual ~L_FiltroDescriptores() {}
};

class L_FiltroDesPolig : public L_FiltroDescriptores
{
public:
	int nVertices;
	int nVerticesMax;
	int *xVer, *yVer; // Define un poligono convexo dentro del cual sobreviven los descriptores
	//int xCen, yCen;   // Puntos contenidos estrictamente dentro del poligono
private:
	L_FiltroDesPolig() {}
	L_FiltroDesPolig(const L_FiltroDesPolig &other); // {}
	L_FiltroDesPolig &operator =(const L_FiltroDesPolig &other); // {return *this;}
public:
	L_FiltroDesPolig(int nVerticesMax) {this->nVertices=0; this->nVerticesMax=nVerticesMax; xVer=new int[nVerticesMax]; yVer=new int[nVerticesMax];}
	void agregaVertice(int x0, int y0); // Los vertices deben ser agregados en el sentido contrario a las agujas del reloj
	bool criterioEliminacion(const L_DescriptorNodo *des);
	bool criterioEliminacion(double x, double y);
	virtual ~L_FiltroDesPolig() {delete[] xVer; delete[] yVer;}
};

class L_FiltroImagenMascara : public L_FiltroDescriptores
{
public:
	L_ImageGrayDouble *mask;
	L_FiltroImagenMascara() {mask=NULL;}
	bool criterioEliminacion(const L_DescriptorNodo *des) {if (mask==NULL) return false; return mask->pix((int)(des->c.x0),(int)(des->c.y0)) == 0;}
};

class L_FiltroDescriptoresGrupo : public L_FiltroDescriptores
{
private:
	L_FiltroDescriptores *filtros[20];
	int nFiltros;
public:
	L_FiltroDescriptoresGrupo() {nFiltros = 0;}
	void addFrom(L_FiltroDescriptores *filtro) {if (nFiltros < 20) {filtros[nFiltros]=filtro; nFiltros ++;} }
	void resetea() {nFiltros=0;}
	bool criterioEliminacion(const L_DescriptorNodo *des) {for (int i=0; i<nFiltros; i++) if (filtros[i]->criterioEliminacion(des)==true) return true; return false;}
};

//////////
// class_name para contar los sigmas usados para generar gaussianas en el sistema

class L_CalcDistrSigmas
{
public:
	long nPun[50];
	L_CalcDistrSigmas() {int i; for (i=0; i<50; i++) nPun[i]=0;}
	int calcNum(double sigma) {return 4 + (int)floor(0.5 + 2 * log(sigma)/L_LOG_2);}
	void AgregaLista(const L_DescriptorLista &desL);
	void print(bool imprDiv=false);
};


/////////
// Clases para almacenar descriptores de forma ordenada para implementar busqueda rapida

class L_ConjDescrBus //! Estructura abstracta para almacenar descriptores ordenados y buscarlos de modo aproximado.
{
public:
	L_ConjDescrBus() {delete_or_isolate_repeated=true;}
	bool delete_or_isolate_repeated; //!< Si hay grupos de nodos iguales, borrarlos todos. Si no, dejar 1 al azar.
	virtual bool createFrom(L_DescriptorLista &desL)=0; //!< createFrom el object, no debe modificar desL
	virtual bool addFrom(L_DescriptorLista &desL)=0; //!< createFrom al object, no debe modificar desL
	virtual bool optimize()=0; //!< optimize el object para hacer que el acceso sea más rápido
	virtual bool destroy()=0; //!< destroy el object, no debe destruir los descriptores usados para entrenar

	//! Busca los n descriptores más parecidos en la base de datos. Debe devolver el mismo tipo de descriptor que se le pasa.
	virtual void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist) =0;

	//!< Busca los n descriptores más parecidos en la base de datos. Debe devolver el mismo tipo de descriptor que se le pasa.
	virtual void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist) =0;

	virtual L_ParamManagerLocal *pideParams() {return NULL;}
	virtual ~L_ConjDescrBus() {}  // Las clases derivadas tienen memoria dinamica propia
};

class L_KdTreeBBF;

class L_KdNodo //! Nodo de un kd-tree. Si se destroy el nodo, no se destroy el árbol
{
public:
	L_KdNodo() : umbral(0), dim(-1), der(NULL), izq(NULL), descr(NULL) {}
	// (umbral,dim,der,izq) sirven cuando es nodo de division
	double umbral;
	int dim;
	L_KdNodo *der; // der->vector->elem[dim] > umbral
	L_KdNodo *izq; // der->vector->elem[dim] <= umbral
	L_DescriptorNodo *descr; // Contiene la memoria del descriptor, no es const *

	// Funciones recursivas
	bool creaRec(void *LDRF_arr, int n, L_KdTreeBBF *kdt);
	void destroyRec() {if (der!=NULL) {der->destroyRec();delete der;} if (izq!=NULL) {izq->destroyRec();delete izq;} if (descr!=NULL) delete descr;}
	// Funciones no recursivas
	static L_KdNodo *creaNoRec(L_KdNodo *inicial, int nTot);
	void destruyeNoRec();	// Faltaba eliminar el descriptor
	static int cmpRepetidos(const void *d1, const void *d2);
	~L_KdNodo() {der=NULL; izq=NULL; descr=NULL;};
};

class L_KdTreeBBF : public L_ConjDescrBus //! kd-tree con búsqueda BBF
{
private:
	L_KdTreeBBF(const L_KdTreeBBF &other); // No existe
	L_KdTreeBBF &operator =(const L_KdTreeBBF &other); // No existe
public:
	L_ParamManagerLocal paramKDTBBF;
	L_KdNodo *root;
	int nHojas; // nº de hojas en el tree
	int nElim; // nº de nodos que no pudieron ser agregados al kd-tree porque ya habia other igual antes

	bool usePercentage; // Usar porcRev o nCompMax como parametro principal
	int nCompMax; // Parametro: nº máximo de comparaciones a realizar
	double porcRev; // Parametro: porcentaje de la base de datos que se revisan
	bool elimRepeatedElements; // Eliminar los nodos repetidos, o dejarlos en la lista original

	L_KdTreeBBF();
	bool createFrom(L_DescriptorLista &desL);
	bool addFrom(L_DescriptorLista &desL);
	bool optimize(); //! optimize el object para hacer que el acceso sea más rápido
	bool destroy(); //! destroy el object
	void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos
	void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos

	bool crea_rec(L_DescriptorLista &desL);
	bool crea_noRec(L_DescriptorLista &desL);

	void swap(L_KdTreeBBF &other)
	{
		L_KdNodo *raiz_t = other.root;
		int nHojas_t = other.nHojas;
		int nElim_t = other.nElim;
		bool usePercentage_t = other.usePercentage;
		int nCompMax_t = other.nCompMax;
		double porcRev_t = other.porcRev;
		bool elimRepeatedElements_t = other.elimRepeatedElements;

		paramKDTBBF.swap(other.paramKDTBBF);

		other.root = root;
		other.nHojas = nHojas;
		other.nElim = nElim;
		other.usePercentage = usePercentage;
		other.nCompMax = nCompMax;
		other.porcRev = porcRev;
		other.elimRepeatedElements = elimRepeatedElements;

		root = raiz_t;
		nHojas = nHojas_t;
		nElim = nElim_t;
		usePercentage = usePercentage_t;
		nCompMax = nCompMax_t;
		porcRev = porcRev_t;
		elimRepeatedElements = elimRepeatedElements_t;
	}
	bool verifyCoherenceOfObject(L_KdNodo *nodoRaiz=NULL, L_KdNodo *nodoRec=NULL);
	L_ParamManagerLocal *pideParams();
	virtual ~L_KdTreeBBF() {destroy();}
};

#define L_BBFPRIORNODO_DISTINF (1.0e9)  // Cantidad que puede considerarse infinita respecto a las distancias entre descriptores
class L_BBFPriorNodo
{
public:
	L_KdNodo *dir;
	double dist;
	L_BBFPriorNodo *sig;
	L_BBFPriorNodo() {dir=NULL; dist=L_BBFPRIORNODO_DISTINF; sig=NULL;}
	//! No debe destruir dir
	void destroyRec() {L_BBFPriorNodo *tmp; while (sig!=NULL) {tmp=sig; sig=sig->sig; delete tmp;}}
	~L_BBFPriorNodo() {}
};

#define new_L_BBFPriorNodo(ret,L_BBFPriorLista_li) {if (L_BBFPriorLista_li.memo!=NULL) {ret = L_BBFPriorLista_li.memo; L_BBFPriorLista_li.memo=L_BBFPriorLista_li.memo->sig; ret->sig=NULL;} else ret = new L_BBFPriorNodo();}
#define delete_L_BBFPriorNodo(ret,L_BBFPriorLista_li) {ret->sig = L_BBFPriorLista_li.memo; L_BBFPriorLista_li.memo = ret;}
//#define new_L_BBFPriorNodo(ret,L_BBFPriorLista_li) {ret = new L_BBFPriorNodo();}
//#define delete_L_BBFPriorNodo(ret,L_BBFPriorLista_li) {delete ret;}

class L_BBFPriorLista
{
public:
	L_BBFPriorNodo *root;
	L_BBFPriorNodo *memo; // Aca se guardan los nodos desechados para no hacer new, delete todo el rato
	double cotaDist;
	int n;
	L_BBFPriorLista() {root=NULL; memo = NULL; cotaDist=100000.0; n=0;} //value muy grande
	void addFrom(L_KdNodo &nododiv, L_KdNodo &nododir, const L_Descriptor &orig, double dist=-1.0);
	L_KdNodo *popNodoMasCercano();
	~L_BBFPriorLista() {if (root!=NULL) {root->destroyRec(); delete root;} if (memo!=NULL) {memo->destroyRec(); delete memo;}}
};

// class_name derivada para almacenamiento secuencial con busqueda secuencial
class L_BusqSecDescr : public L_ConjDescrBus //! Busqueda secuencial
{
public:
	enum boolVBC {sin_error, error_primero, error_segundo, error_ambos};
	L_DescriptorPtrLista lista;
	L_BusqSecDescr() {}
	bool createFrom(L_DescriptorLista &desL);
	bool addFrom(L_DescriptorLista &desL);
	bool optimize(); //! optimize el object para hacer que el acceso sea más rápido
	bool destroy(); //! destroy el object
	void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos
	void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos
	boolVBC verifBusqCorrecta(const L_Descriptor &orig, const L_Descriptor &prim, const L_Descriptor &seg); // Verifica que (prim,seg) sean realmente el primero y segundo más cercanos
	virtual ~L_BusqSecDescr() {destroy();};
};

// Permutacion asociada a un elemento
class L_PermPiv
{
public:
	L_DescriptorNodo *des; // descriptor, referencia a memo externa
	std::vector<int> perm; // pemutacion
	int dist; // distancia para ordenarlos con qsort
	L_PermPiv()
	{
		des=NULL;
		dist=0;
	}
	~L_PermPiv()
	{
	}
};

// class_name derivada para busqueda usando vectores de distancias a pivotes (Algoritmo: Karina Figueroa)
class L_BusqPermPiv : public L_ConjDescrBus
{
	struct dist_str {int ind; double d;};
public:
	L_ParamManagerLocal paramBusPP;
	std::vector<L_DescriptorNodo *> piv; // pivotes
	int nPiv; // numero de pivotes
	std::vector<L_PermPiv> bd; // base de datos
	int nBd; // numero de descriptores en base de datos
	std::vector<dist_str> dist; // Estructura para almacenamiento temporal
	std::vector<int> temp; // temporal para uso en calcSpearman()

	bool usePercentage; // Usar porcRev o nCompMax como parametro principal
	int nCompMax; // Parametro: nº máximo de comparaciones a realizar
	double porcRev; // Parametro: porcentaje de la base de datos que se revisan
	L_BusqPermPiv():paramBusPP("paramBusPiv",0)
	{
		nPiv=20; // El numero de pivotes debe ser bastante menor a la dimension del vector para que el metodo sea rapido.
		nBd=0;

		usePercentage=true; // Indica cual de las proximas dos opciones se usara para limitar la busqueda
		nCompMax=200; // cantidad maxima de comparaciones
		porcRev=0.25; // % de la base de datos que se revisa

		paramBusPP.addFrom("usePercentage",&usePercentage);
		paramBusPP.addFrom("nCompMax",&nCompMax);
		paramBusPP.addFrom("porcRev",&porcRev);
		paramBusPP.addFrom("nPiv",&nPiv);
	};
	bool createFrom(L_DescriptorLista &desL); // Primero que todo, transformUsing desc en un L_DescriptorRefOrd[]
	bool addFrom(L_DescriptorLista &desL) {throw L_NonImplementedException(); return false;} // Hay que ver como hacer esto!!
	bool optimize(); //! optimize el object para hacer que el acceso sea más rápido
	void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos
	void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist); //! Busca los n descriptores más parecidos en la base de datos
	L_ParamManagerLocal *pideParams();

	static int cmpDist_str(const void *a, const void *b);
	static int cmpPermDist(const void *a, const void *b);
	void calcPerm(const L_Descriptor &des, L_PermPiv &permObj);
	int calcSpearman(const L_PermPiv &p1, const L_PermPiv &p2);
	void ordenaPiv(const L_Descriptor &orig);
};


#ifdef FLANN_H
class L_Flann : public L_ConjDescrBus
{
public:
	FLANN_INDEX index;
	std::vector<L_Descriptor> desArr;
	std::vector<float> data;
	L_ParamManagerLocal params;
	FLANNParameters flannParams;

	L_Flann() : params("FLANN", 2) {index = NULL; /*flannParams.target_precision = (float)0.95;*/ fijaParamB(); params.addFrom("target_precision", &flannParams.target_precision);}

	void fijaParamB();

	bool createFrom(L_DescriptorLista &desc);
	bool addFrom(L_DescriptorLista &desc) {throw L_NonImplementedException(); return false;}
	bool optimize() {throw L_NonImplementedException(); return false;}
	bool destroy();
	void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist);
	void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist);
	virtual void encuentraMasCercanos(std::vector<L_Descriptor> &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist);
	L_ParamManagerLocal *pideParams() {return &params;}

	~L_Flann() {destroy();}
};
#endif // FLANN_H

template <class T>
class L_BusArr : public L_ConjDescrBus
{
private:
	L_Array<T> arr;
	L_Array<L_DescriptorLista> desLarr;
	int tipos[600];
	int nMaxTipos;
	int nMax;
	int n;
	L_ParamManagerLocal paramsBusArr;
	L_BusArr():paramsBusArr("L_BusArr", 2) {} // No debe existir
	L_BusArr(const L_BusArr &other):paramsBusArr("L_BusArr", other.nMax+2) { }
	const L_BusArr &operator =(const L_BusArr &other) { return *this; }
public:
	L_BusArr(int nMax):paramsBusArr("L_BusArr", nMax+2), nMaxTipos(600) {this->nMax=nMax; n=0;}
	bool createFrom(L_DescriptorLista &desc)
	{
		int i;
		L_DescriptorNodo *ptr;
		n=0;
		if (arr.size() == 0)
		{
			arr.resize_swapping(nMax);
			// Solamente arr[nMax-1] debe tener un servidor de parametros que altere el flujo de informacion
			for (i=0; i<nMax-1; i++)
				if (arr[i].pideParams() != NULL)
					arr[i].pideParams()->setValue_modifiedFlowOfInfo(false);
			for (i=0; i<nMax; i++)
				if (arr[i].pideParams() != NULL)
					paramsBusArr.addChildren(arr[i].pideParams());
		}
		if (desLarr.size() == 0)
		{
			desLarr.resize(nMax);
		}
		for (i=0; i<nMaxTipos; i++)
			tipos[i]=-1;
		for (ptr=desc.root; ptr!=NULL; ptr=ptr->sig)
		{
			if (tipos[ptr->c.idTipo()]==-1)
			{
				tipos[ptr->c.idTipo()]=n;
				n++;
			}
			if (n>nMax)
			{
				printf("L_BusArr: nMax demasiado chico\n");
				return false;
			}
		}
		for (ptr=desc.root; ptr!=NULL; ptr=ptr->sig)
		{
			if (tipos[ptr->c.idTipo()]==-1)
			{
				printf("L_BusArr: Error raro\n");
				return false;
			}
			*desLarr[tipos[ptr->c.idTipo()]].pult=ptr;
			desLarr[tipos[ptr->c.idTipo()]].pult=&(*desLarr[tipos[ptr->c.idTipo()]].pult)->sig;
			desLarr[tipos[ptr->c.idTipo()]].n++;
		}
		desc.root=NULL;
		desc.pult = &desc.root;
		desc.resize(0);
		for (i=0; i<n; i++)
		{
			arr[i].createFrom(desLarr[i]);
			desLarr[i].clear(); // Por si acaso
		}
		return true;
	}
	bool addFrom(L_DescriptorLista &desc) {throw L_NonImplementedException(); return false;}
	bool optimize() {int i; bool a=true; if (arr.size() != 0) for (i=0; i<nMax; i++) a&=arr[i].optimize(); return a;}
	bool destroy() {int i; bool a=true; if (arr.size() != 0) for (i=0; i<nMax; i++) a&=arr[i].destroy(); return a;}
	void findNearest_anyObject(const L_Descriptor &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
	{
		int i;
		if (tipos[orig.idTipo()]==-1)
		{
			//printf("L_BusArr: Tipo de punto de interés no presente durante el entrenamiento\n");
			for (i=0; i<n; i++)
			{
				dest[i]=NULL;
				dist[i]=10000.0;
			}
			return;
		}
		return arr[tipos[orig.idTipo()]].findNearest_anyObject(orig, n, dest, dist);
	}
	void findNearest_otherObject(const L_Descriptor &orig, int n,  std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
	{
		int i;
		if (tipos[orig.idTipo()]==-1)
		{
			//printf("L_BusArr: Tipo de punto de interés no presente durante el entrenamiento\n");
			for (i=0; i<n; i++)
			{
				dest[i]=NULL;
				dist[i]=10000.0;
			}
			return;
		}
		return arr[tipos[orig.idTipo()]].findNearest_otherObject(orig, n, dest, dist);
	}
	void encuentraMasCercanos(std::vector<L_Descriptor> &orig, int n, std::vector<L_Descriptor *> &dest, std::vector<double> &dist)
	{
		std::vector<L_Descriptor *> dest_tmp;
		std::vector<double> dist_tmp;
		int i;
		dest_tmp.resize(2);
		dist_tmp.resize(2);
		for (i=0; i<orig.size(); i++)
		{
			dest_tmp.elem = &dest[n*i];
			dist_tmp.elem = &dist[n*i];
			findNearest_anyObject(orig[i], n, dest_tmp, dist_tmp);
		}
		dest_tmp.elem = NULL;
		dist_tmp.elem = NULL;
	}
	L_ParamManagerLocal *pideParams() {return &paramsBusArr;}
};
#undef L_MAX_PUNTINT_TIPOS


///////
// Clases para define calces y listas de calces

class L_Calce // Es un object POD (estatico), no necesita copias. LOS DESCRIPTORES DEBEN SOBREVIVIR A ESTE object
{
public:
	// Se dejan como const, si se necesita modificar alguno de los dos usar const_cast<L_DescriptorNodo *>
	const L_Descriptor *dRef; // No tiene memoria propia, debe tener nobj>=0 (numero de imagen de entrenamiento)
	const L_Descriptor *dPru; // No tiene memoria propia, puede tener cualquier nobj (en principio no se usa)
	double vectDist;
	double vectDivDist;
	L_Calce() : dRef(NULL), dPru(NULL) { }
	static int qcmpRef(const void *calceNodopptr1, const void *calceNodopptr2);
	static int qcmpPosSigmaAng(const void *calceNodopptr1, const void *calceNodopptr2);
};

class L_TransfAfinPI2D;

typedef L_Node<L_Calce> L_CalceNodo;
class L_CalceLista:public L_List<L_Calce>
{
public:
	void insertaNodoAlInicio_(L_DescriptorLista &desLRef, L_DescriptorLista &desLPru, int xRef, int yRef, int xPru, int yPru, L_ImageGrayDouble *imRef=NULL, L_ImageGrayDouble *imPru=NULL);
	void insertaNodoInverso(L_Calce &cal) {*pult=new L_CalceNodo; (*pult)->c.dPru=cal.dRef; (*pult)->c.dRef=cal.dPru; pult=&(*pult)->sig; *pult=NULL; n++;}

	void genCalces(L_ConjDescrBus &db, L_DescriptorLista &desL, bool restrDeltaAng, double deltaAng, L_List<L_Calce> *cal_pool=NULL);

	void genCalces_nMasCercanos(L_ConjDescrBus &db, L_DescriptorLista &desL, bool restrDeltaAng, double deltaAng, L_List<L_Calce> *cal_pool=NULL);
	double genCalcesVerif(L_BusqSecDescr &sec, L_ConjDescrBus &db, L_DescriptorLista &desL, double distmax, double ratio);
	void parear(L_DescriptorLista &desL1, L_DescriptorLista &desL2);
	void fijaPunteroImagenReferencia(const L_ImageGrayDouble *imRef) {for (L_CalceNodo *cal=root; cal!=NULL; cal=cal->sig) const_cast<L_Descriptor *>(cal->c.dRef)->imorig=imRef;}
	void fijaPunteroImagenPrueba(const L_ImageGrayDouble *imPru) {for (L_CalceNodo *cal=root; cal!=NULL; cal=cal->sig) const_cast<L_Descriptor *>(cal->c.dPru)->imorig=imPru;}
	void invierteRefPru() {const L_Descriptor *d; for (L_CalceNodo *cal=root; cal!=NULL; cal=cal->sig) {d=cal->c.dPru; cal->c.dPru=cal->c.dRef; cal->c.dRef=d;}}
	void genLineas(L_ShapeArray &lins, int nImRef=-1);
	bool dibuja(const char *nomarch); // Solo puede dibujar si hay calces
	bool dibuja(L_ImageRGBUchar &im); // Solo puede dibujar si hay calces
	void dibuja(L_ImageRGBUchar &im, L_ImageGrayDouble &imRef, L_ImageGrayDouble &imPru); // Solo puede dibujar si hay calces
	double calcAngProm();
	bool revisaCoherencia();
	bool restrSemiLocales(const L_Calce &central, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng); // 2 <= minCorresp <= nVecinos <= this->n-1 , varEscala > 1 , 0 < varAng < M_PI
	int aplicaRestrSemiLocales(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng); // Devuelve el numero de elementos en los cuales fallan las restricciones semi locales
	double aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax);
	bool cuentaCerros(L_Calce &c1, L_Calce &c2, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin);
	static int cmpCalceXY(const void *a, const void *b);
	static int cmpCalceXY_v2(const void *a, const void *b);
	void transformaEnMatrices(L_Matrix &uv_1, L_Matrix &uv_2);
	bool verifCorrLinDatos(double rMin, double *corrPreg=NULL);
	void ordena_id_pru();

	void eliminaRepetidos() {sort(1, &L_CalceLista::cmpCalceXY);}

	int reemplazarReferenciasDescriptores(L_DescriptorLista &desPru);

	void pedirMatrizDatos(L_Matrix &x1y1x2y2);
	void pedirMatrizDatosRayos(L_Matrix &x1y1x2y2, const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru);
	void pedirMatrizPuntajes(std::vector<double> &puntajes);
	void eliminaDescriptoresSinLandmark();
	void eliminaLineasHarrisPru(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	void eliminaLineasHarrisPru_indiv(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	void eliminaLineasHarrisPru_ant(double cornernessMin, double sd=1.3, double si=2.0, double alfa=0.04);
	// Estos eliminan calces incorrectos usando el algoritmo de los 5 puntos con RANSAC o el de los 8 puntos, depende
	int eliminaEpipolar_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E);
	int eliminaEpipolarPuntajes_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E);
	int eliminaEpipolarPuntajes_Ransac6Puntos(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, L_TransfAfinPI2D &tr, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_EssentialMatrix &E);
	void eliminaEpipolar(const L_CamaraPinhole &cRef, const L_CamaraPinhole &cPru, double errMax, const L_EssentialMatrix &E); // Aca la matriz esencial ya se conoce
	int eliminaEpipolar_Ransac8Puntos(double errMax, int nHip, int nMinAcept, int nMin, L_MatrizFundamental &F);
	int eliminaEpipolarPuntajes_Ransac8Puntos(double errMax, int nHip, int nMinAcept, int nAceptAutom, L_MatrizFundamental &F);
	int eliminaEpipolarPuntajes_Ransac8Puntos(L_TransfAfinPI2D &tr, double errMax, int nHip, int nMinAcept, int nAceptAutom, L_MatrizFundamental &F);
	void eliminaEpipolar(double errMax, const L_MatrizFundamental &F); // Aca la matriz fundamental ya se conoce
	// Estos eliminan calces correctos usando el algoritmo de los 3 puntos
	int eliminaPose_Ransac3Puntos(int proyRef0Pru1, const L_CamaraPinhole &cam, double errMax, int nIter, int nMinAcept, int nAceptAutom, const std::vector<L_CoordsCart3D> &pun, const std::vector<double> &puntajes, L_Pose3D_cuat &pose);
	// Transf de semejanza
	int eliminaSemejanza_Ransac2puntos(double errMax, int nHip, int nMinAcept, int nAceptAutom, double &rot, double &esc, double &tx, double &ty);
	int eliminaDelta(double dMax, double dangMax, double descMax); // erase_preserving_order siempre, cuidado...

	static void ejemploRansac(int npuntos);  // requiere ventana
	static void ajustarRansac(int npuntos);
};


class L_Calcec2 // Es un object estatico, no necesita copias
{
public:
	const L_Descriptor *dRef; // No tiene memoria propia
	const L_Descriptor *dPru; // No tiene memoria propia
	const L_Descriptor *dPruc2; // No tiene memoria propia
	L_Calcec2() {dRef=NULL; dPru=NULL; dPruc2=NULL;}
};

typedef L_Node<L_Calcec2> L_Calcec2Nodo;

class L_Calcec2Lista : public L_List<L_Calcec2>
{
public:
	void parear(const L_CalceLista &c1, const L_CalceLista &c2);
	void agregar(const L_CalceLista &c1);
};

//////
// Clases para define relaciones geometricas entre dos imagenes

//////
// Clases para define transformaciones geometricas, son un caso particular de relacion geometrica

enum L_Transf2DTipo {L_Tr_Indef,L_Tr_Afin,L_Tr_Proyectiva};

enum L_Transf2DDeformTipo {
	PIXELADO, // No se usa filtro anti-aliasing ni interpolacion al transformar la imagen
	BILIN_1NIV_ISOTR, // Filtro bilineal para interpolar, filtro antialiasing isotrópico de 1 nivel (usa peor caso)
	BILIN_1NIV_ANISOTR, // Filtro bilineal para interpolar, filtro antialiasing anisotrópico de 1 nivel (usa peor caso)
	BILIN_XNIV_ISOTR,  // Filtro bilineal para interpolar, filtro antialiasing anisotrópico adaptivo (usa peor caso)
	BILIN_XNIV_ANISOTR, // Filtro bilineal para interpolar, filtro antialiasing anisotrópico adaptivo
	SPLINE_1NIV_ISOTR, // Filtro bi-spline para interpolar, filtro antialiasing isotrópico de 1 nivel (usa peor caso)
	SPLINE_1NIV_ANISOTR, // Filtro bi-spline para interpolar, filtro antialiasing anisotrópico de 1 nivel (usa peor caso)
	SPLINE_XNIV_ISOTR,  // Filtro bi-spline para interpolar, filtro antialiasing anisotrópico adaptivo (usa peor caso)
	SPLINE_XNIV_ANISOTR // Filtro bi-spline para interpolar, filtro antialiasing anisotrópico adaptivo
};

class L_Transf2D_caract // Vector de caracteristicas para cada transformacion
{
public:
	int n; // Numero de caracteristicas

	double numVotosTransfFus; // Despues de la fusion
	double numVotosTransf; // Antes de la fusion
	double probCelda;
	double probTransf;
	double numVotosCelda; // Numero inicial de votos de la celda (Hough)
	double distRomb;
	double corrLin;
	double corrPix;
	double ransac;
	double calces_1_a_n;
	double RSLfallados;
	double porcCerrosCorrectos;
	double punt;  //  = probTransf + ransac

	double &el(int i);  // No necesita ser tan rapido

	static L_String imprimeCaract(); // Entrega texto con las caracteristicas ordenadas

	L_Transf2D_caract() : n(13),
		numVotosTransfFus(0),
		numVotosTransf(0),
		probCelda(0),
		probTransf(0),
		numVotosCelda(0),
		distRomb(0),
		corrLin(0),
		corrPix(0),
		ransac(0),
		calces_1_a_n(0),
		RSLfallados(0),
		porcCerrosCorrectos(0),
		punt(0) {}
};

class L_Transf2D_POD
{
public:
	L_Transf2DTipo tipo; // Tipo de transformacion

	// Caracteristicas
	L_Transf2D_caract car;
#ifdef __TURBOC__
	__OP_ASSIGN(L_Transf2D_POD)
#endif
};

class L_TransfPuntInt2D : public L_Transf2D_POD //! Transformación geométrica abstracta. Mapea la imagen ref en la imagen pru.
{
public:
	L_CalceLista calL; // Lista de calces que generan la transformacion

	L_TransfPuntInt2D() {tipo=L_Tr_Indef;}
	L_TransfPuntInt2D(L_Transf2DTipo tipo) {this->tipo=tipo;}

private:
	// Deshabilitadas por eficiencia, mejor usar las swap() o copiaObjetoEn() explicitamente
	L_TransfPuntInt2D(const L_TransfPuntInt2D &other); // {L_Transf2D_POD::operator=(other); other.calL.copyListOn(calL);}
	L_TransfPuntInt2D &operator=(const L_TransfPuntInt2D &other);//  {L_Transf2D_POD::operator=(other); calL.clear(); other.calL.copyListOn(calL);}
public:

	virtual void proyeccionDe(double &resx, double &resy, double x, double y) const = 0;
	virtual bool calcTransfExacto(const L_Calce *calArr)=0; // #calArr = nParam
	virtual bool calcTransfMinCuad()=0;
	virtual bool invertMe()=0; // Devuelve false si no existe la inversa
	virtual bool componeCon(L_TransfPuntInt2D& t)=0;
	virtual bool calcTransfRansac(double errX, double errY, int nMinAcept, int nIntentos)=0; // Devuelve false si falla o no existe. El error es en el espacio de la imagen de referencia.
	virtual L_TransfPuntInt2D* clona()=0;
	virtual L_TransfPuntInt2D* clonaRegalaCalces()=0;
	virtual bool deformaImagen(L_ImageGrayDouble &dest, const L_ImageGrayDouble &orig, L_Transf2DDeformTipo tipo=BILIN_1NIV_ISOTR) const = 0;

	double errorRelGeom(double x1, double y1, double x2, double y2) {double xP, yP, ex, ey; proyeccionDe(xP,yP,x1,y1); ex = x2-xP; ey = y2-yP; return sqrt(ex*ex+ey*ey);}

	void proyeccionDe(const L_PuntInt &p1, L_PuntInt &resp2) const {this->proyeccionDe(resp2.x0, resp2.y0, p1.x0, p1.y0);}

	int pideNumObj() {if (calL.root!=NULL && calL.root->c.dRef!=NULL) return calL.root->c.dRef->nobj; else return -1;}

	void cuenta_calces_1_a_n();
	bool aplicaRSL(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng);
	bool aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax);

	static int qcmp(const void *tr1, const void *tr2);

	virtual ~L_TransfPuntInt2D() {}
};

typedef L_NodePtr<L_TransfPuntInt2D> L_TransfPuntInt2DNodoPtr;

class L_TransfPuntInt2DListaPtr:public L_ListPtr<L_TransfPuntInt2D>
{
public:
	int calcNumTransf();
	int calcNumTransf(L_CountingTree &numTransf);
	void sort();
	void eliminaCalcesRepetidos() {	L_TransfPuntInt2DNodoPtr *tr; for (tr=root; tr!=NULL; tr=tr->sig) tr->c->calL.sort(1, &L_CalceLista::cmpCalceXY);}

	void escribeResumen(FILE *fp);
	void escribeResumen(const char *nomarch);
	void cuenta_calces_1_a_n();
	void aplicaRSL(bool corrigeLista, int nVecinos, int minCorresp, bool verifEscalaAngs, double varEscala, double varAng);
	void aplicaCuentaCerros(bool corrigeLista, double minDistCerros, double maxDistCerros, int minCerrosUtil, double parecidoCerrosMin, int nLineasMax, double porcCerrosCorrectosMin);
	void dibujaTransformaciones(L_ShapeArray &lins, bool calces1_o_caja0,  int maxImRef, int maxTransfPorImagen, int maxTransfTotales, int imRefEspecifica=-1, int transfEspecifica=-1); // Mejor llamar a L_RecPuntInt.genDrawing()
	L_TransfPuntInt2D *buscaTransf(int nobj, int indice);

    void destroyList(L_List<L_Calce> *cal_pool);
};

class L_TransfAfinPI2D: public L_TransfPuntInt2D, public L_TransfAfin2D //! Transformación afín
{
	// param[] = {m11,m12,tx,m21,m22,ty}
private: // Por eficiencia
	L_TransfAfinPI2D(const L_TransfAfinPI2D &other);
	L_TransfAfinPI2D operator=(const L_TransfAfinPI2D &other);
public:
	L_TransfAfinPI2D():L_TransfPuntInt2D(L_Tr_Afin) {};
	void proyeccionDe(double &resx, double &resy, double x, double y) const {resx = m11*x + m12*y + tx; resy = m21*x + m22*y + ty;}
	void proyeccionDe(const L_PuntInt &p1, L_PuntInt &resp2) const {this->proyeccionDe(resp2.x0, resp2.y0, p1.x0, p1.y0);}
	bool calcTransfExacto(const L_Calce *calArr); // #calArr = nParam
	bool calcTransfMinCuad();
	bool calcParametrosInversosDe(const L_TransfAfinPI2D &other);
	bool invertMe() {return calcParametrosInversosDe(*this);}
	bool componeCon(L_TransfPuntInt2D& t) {return false;}
	bool calcTransfRansac(double errX, double errY, int nMinAcept, int nIntentos); // El error es en el espacio de la imagen de referencia
	bool deformaImagen(L_ImageGrayDouble &dest, const L_ImageGrayDouble &orig, L_Transf2DDeformTipo tipo=PIXELADO) const; // Solo esta implementado para PIXELADO por ahora
	bool calcVarParametros(L_Matrix &var); // Calcula la transformación y el error asociado
	void linealiza(L_TransfPuntInt2D *tr, double x0, double y0);
	L_TransfAfinPI2D *clona();
	L_TransfAfinPI2D *clonaRegalaCalces();

	void swap(L_TransfAfinPI2D &other);

	~L_TransfAfinPI2D() {}
};

class L_TransfProyectivaPI2D: public L_TransfProyectiva2D, public L_TransfPuntInt2D //! Transformación geométrica abstracta. Mapea la imagen ref en la imagen pru.
{
	// param[] = {m11,m12,tx,m21,m22,ty,m31,m32}
public:
	L_CamaraDistorsionRadial camInter;
	bool ransac_2cal;
	bool distorsionarCoordImagenPrueba;
	bool distorsionarCoordImagenReferencia;
	L_TransfProyectivaPI2D():L_TransfPuntInt2D(L_Tr_Proyectiva) {ransac_2cal=true; distorsionarCoordImagenPrueba=false; distorsionarCoordImagenReferencia=false;}
	void proyeccionDe(double &resx, double &resy, double x, double y) const {double den = m31*x + m32*y + 1; resx = (m11*x + m12*y + tx)/den; resy = (m21*x + m22*y + ty)/den;}
	void proyeccionDe(const L_PuntInt &p1, L_PuntInt &resp2) const {L_TransfProyectiva2D::proyeccionDe(resp2.x0, resp2.y0, p1.x0, p1.y0);}
	bool calcTransfExacto(const L_Calce *calArr); // #calArr = nParam = 8
	bool calcTransfExacto(double *xR, double *yR, double *xP, double *yP); // #calArr = nParam = 8
	bool calcTransfMinCuad();  // No es exacto, no minimize en el espacio de los pixeles sino en other...
	bool calcTransfMinCuad(double *xR, double *yR, double *xP, double *yP, int n);  // No es exacto, no minimize en el espacio de los pixeles sino en other...
	bool calcParametrosInversosDe(const L_TransfProyectivaPI2D &other);
	bool invertMe() {return calcParametrosInversosDe(*this);}
	bool componeCon(L_TransfPuntInt2D& t) {return false;}
	bool calcTransfRansac_2cal_ant(double errX, double errY, int nMinAcept, int nIntentos);
	bool calcTransfRansac_4cal_ant(double errX, double errY, int nMinAcept, int nIntentos);
	bool calcTransfRansac_ant(double errX, double errY, int nMinAcept, int nIntentos); // El error es en el espacio de la imagen de referencia
	bool calcTransfRansac(double errX, double errY, int nMinAcept, int nIntentos); // El error es en el espacio de la imagen de referencia
	bool deformaImagen(L_ImageGrayDouble &dest, const L_ImageGrayDouble &orig, L_Transf2DDeformTipo tipo=PIXELADO) const; // Solo esta implementado para PIXELADO por ahora
	L_TransfProyectivaPI2D *clona();
	L_TransfProyectivaPI2D *clonaRegalaCalces();
	~L_TransfProyectivaPI2D() {}
};


//////////////////////////////
// Clases para calculadores de descriptores comunes a varios generadores de descriptores

class L_SIFTdireccLista //! Lista de direcciones para puntos SIFT
{
public:
	L_ParamManagerLocal paramSIFTDL;
	int nd; //! nº de orientaciones posibles a considerar para colocar sistema de referencia local
	double porc; //! porcentaje del máximo para generar dirección principal
	std::vector<double> gaus;
	int ngaus;
	double sgaus;
	double facS;
	bool orig;
	L_DireccLista lis;
	L_SIFTdireccLista():paramSIFTDL("SIFTDL",0)
	{
		nd=8; porc=0.8; ngaus=0; sgaus=0; facS=1.5;
		orig=false; // True: usar funcion propia, false: usar funcion pirateada
		paramSIFTDL.addFrom("nd",&nd);
		paramSIFTDL.addFrom("porc",&porc);
		paramSIFTDL.addFrom("facS",&facS);
		paramSIFTDL.addFrom("orig",&orig);
	}
	void calcDirecciones(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin); // gradPolares=true => grad1=abs(grad), grad2=atan2(gradAba,gradDer) ; gradPolares=false => grad1=gradDer ; grad2=gradAba
	void clear()
	{
		lis.clear();
	}
	~L_SIFTdireccLista() {}
private:
	void calcDirecciones_PLZ(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin);
	// Funciones Lowe
	void calcDirecciones_Lowe(int xc, int yc, L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares, double s, L_PuntIntNodo *pin);
	void SmoothHistogram(std::vector<double> &hist, int bins);
	double FindOriPeaks(std::vector<double> &hist, int bins);
};

class L_SIFTCalculador:public L_Descriptor //! Calculador de descriptor SIFT
{
public:
	L_ParamManagerLocal paramSIFT;
	int nx, ny, nd, lx, ly; // (nx, ny) deben ser pares
	std::vector<double> gdx_h, gdy_h;
	double sigma;
	bool orig;
	L_SIFTCalculador(int nx_=4, int ny_=4, int nd_=8, int lx_=4, int ly_=4, double sigma_=0.5) : L_Descriptor(),paramSIFT("SIFT",0)
	{
		nx=nx_;
		ny=ny_;
		nd=nd_;
		lx=lx_;
		ly=ly_;
		sigma=sigma_;
		tipoD=SIFT;
		nt=nx*ny*nd;
		orig=false;  // true: usar funcion propia; false: usar funcion pirateada
		paramSIFT.addFrom("nx",&nx);
		paramSIFT.addFrom("ny",&ny);
		paramSIFT.addFrom("nd",&nd);
		paramSIFT.addFrom("lx",&lx);
		paramSIFT.addFrom("ly",&ly);
		paramSIFT.addFrom("orig",&orig);
	}
	//Funcion central de calculo de descriptores: llama a una funcion particular
	bool calcDescriptor(L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares);

private:
	bool calcDescriptorGenerico_norot(L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares);

	///// Funciones Lowe
	bool calcDescriptorGenerico_Lowe(L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares);
	void NormalizeVec(double *vec, int len);
	void KeySample(double *fvec, L_ImageGrayDouble & grad, L_ImageGrayDouble & ori, int row, int col, bool gradPolares);
	void AddSample(double ***index, L_ImageGrayDouble & grad, L_ImageGrayDouble & orim, int r, int c, double rloc, double cloc, bool gradPolares);
	void PlaceInIndex(double ***index, double mag, double ori, double rloc, double cloc);
	void PlaceInIndexVec(double ***index, double mag, double ori, double rloc, double cloc);
public:

	~L_SIFTCalculador() {}
};


// Tipo de datos L_SIFTGenFnDescriptor : "puntero a fn de calculo de descriptor en L_SIFTCalculador"
typedef bool (L_SIFTCalculador::*L_SIFTGenFnDescriptor)(L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares);


class L_SIFTGenFnDescriptorNodo
{
public:
	L_SIFTGenFnDescriptor f;
	long id; // info de lx,ly,nd,nx,ny
	long blur; // info de sigma =(long)(s*10000);
	static int cmp(const void *nodo1, const void *nodo2)
	{
		if ( (*(L_SIFTGenFnDescriptorNodo *)nodo1).id>(*(L_SIFTGenFnDescriptorNodo *)nodo2).id )
			return 1;
		else if ( (*(L_SIFTGenFnDescriptorNodo *)nodo1).id<(*(L_SIFTGenFnDescriptorNodo *)nodo2).id )
			return -1;
		else if ( (*(L_SIFTGenFnDescriptorNodo *)nodo1).blur>(*(L_SIFTGenFnDescriptorNodo *)nodo2).blur )
			return 1;
		else if ( (*(L_SIFTGenFnDescriptorNodo *)nodo1).blur<(*(L_SIFTGenFnDescriptorNodo *)nodo2).blur )
			return -1;
		else
			return 0;
	}
	static long calcId(int lx, int ly, int nd, int nx, int ny)
	{
		#if defined(MAXlx) || defined(MAXly) || defined(MAXnd) || defined(MAXnx) || defined(MAXny)
			#error #define conflicts
		#endif
		#define MAXlx 8L
		#define MAXly 8L
		#define MAXnd 16L
		#define MAXnx 8L
		#define MAXny 8L
		long res;
		res = lx;
		res = res*MAXly + ly;
		res = res*MAXnd + nd;
		res = res*MAXnx + nx;
		res = res*MAXny + ny;
		return res;
		#undef MAXlx
		#undef MAXly
		#undef MAXnd
		#undef MAXnx
		#undef MAXny
	}
	static long calcBlur(double sigma) {return (long)(sigma*10000);}
};

class L_SIFTGenFnDescriptorArreglo
{
public:
	L_SIFTGenFnDescriptorNodo *arr;
	int n;
	int nMem;
	bool activo;
	L_SIFTGenFnDescriptorArreglo(bool activar)
	{
		arr=NULL;
		n=0;
		nMem=0;
		activo=activar;

		//addFrom(&L_SIFTCalculador::funcion, lx, ly, nd, nx, ny);

		sort();
	};
	void addFrom(L_SIFTGenFnDescriptor f, int lx, int ly, int nd, int nx, int ny, double sigma)
	{
		if (n==nMem)
		{
			L_SIFTGenFnDescriptorNodo *arr2;
			arr2=(L_SIFTGenFnDescriptorNodo *)realloc(arr, (nMem+10) * sizeof(L_SIFTGenFnDescriptorNodo));
			if (arr2==NULL)
				throw std::bad_alloc();
			else
				arr=arr2;
			nMem+=10;
		}
		arr[n].f=f;
		arr[n].id=L_SIFTGenFnDescriptorNodo::calcId(lx, ly, nd, nx, ny);
		arr[n].blur=L_SIFTGenFnDescriptorNodo::calcBlur(sigma);
		n++;
	}
	void sort() {qsort(arr,n,sizeof(L_SIFTGenFnDescriptorNodo),L_SIFTGenFnDescriptorNodo::cmp);}
	L_SIFTGenFnDescriptorNodo *buscaFn(int lx, int ly, int nd, int nx, int ny) const
	{
		L_SIFTGenFnDescriptorNodo f0;
		if (!activo || n==0)
			return NULL;
		f0.f=(L_SIFTGenFnDescriptor)NULL;
		f0.id=L_SIFTGenFnDescriptorNodo::calcId(lx, ly, nd, nx, ny);
		return (L_SIFTGenFnDescriptorNodo *)bsearch(&f0,arr,n,sizeof(L_SIFTGenFnDescriptorNodo),L_SIFTGenFnDescriptorNodo::cmp);
	}
	~L_SIFTGenFnDescriptorArreglo()	{ if (arr!=NULL) free(arr); }
};






class L_SURFGaussCalculador:public L_Descriptor //! Calculador de descriptor SIFT
{
public:
	L_ParamManagerLocal paramSURFG;
	int nx, ny, lx, ly; // (nx, ny) deben ser pares
	std::vector<double> gdx_h, gdy_h;
	double sigma;

	L_SURFGaussCalculador(int nx_=4, int ny_=4, int lx_=5, int ly_=5, double sigma_=0.5) : L_Descriptor(),paramSURFG("SURFG",0)
	{
		nx=nx_;
		ny=ny_;
		lx=lx_;
		ly=ly_;
		sigma=sigma_;
		tipoD=D_SURFGAUSS;
		nt=nx*ny*4; // grder, grarr, |grder|, |grarr|
		paramSURFG.addFrom("nx",&nx);
		paramSURFG.addFrom("ny",&ny);
		paramSURFG.addFrom("lx",&lx);
		paramSURFG.addFrom("ly",&ly);
	}
	//Funcion central de calculo de descriptores: llama a una funcion particular
	bool calcDescriptor(L_ImageGrayDouble & grad1, L_ImageGrayDouble & grad2, bool gradPolares);

};


#endif //__L_DESCRIPTORES_H__

